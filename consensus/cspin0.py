import os
import time
import ctypes

import cProfile
import multiprocessing as mp

from itertools import combinations, product
from collections import defaultdict
from multiprocessing import Pool
from Queue import Empty

from .engraph import RNAStructureGraph, LabelsList, Label, Domain, DFSCode,\
                     DFSCodeFragment, RNAStructureEngrapher, DomainChange,\
                     DFSGenerator as DFSGen


class GraphList(list):

    IGNORE = (RNAStructureGraph.EDGELABELS[RNAStructureGraph.ISPARTOF.name],)

    def __init__(self, graphs):
        super(GraphList, self).__init__(graphs)

    def get_domain(self):
        nodelabels, edgelabels = set(), set()
        for graph in self:
            nodelabels.update(graph.nodes.values())
            edgelabels.update({edge[2] for edge in graph.edges})
        nodelabels = RNAStructureGraph.LABELS.node.create_subset(nodelabels)
        edgelabels = RNAStructureGraph.LABELS.edge.create_subset(edgelabels)
        domain = Domain(nodelabels, edgelabels)
        return domain

    def embed_code(self, base):
        dfslist = sum([graph.embed_code(base) for graph in self], list())
        return dfslist

    def remove_nodes_labelled(self, labels):
        for graph in self:
            nodeids = graph.nodes_with_label_in(labels)
            for nodeid in nodeids:
                graph.delete_node(nodeid)

    def remove_edges_labelled(self, edgelabels, ignore=list()):
        """ Removes edges with given labels *and* 
            deletes isolated nodes afterwards. 
        """
        for graph in self:
            edges = graph.edges_with_label_in(edgelabels)
            graph.edges -= edges
            incident_nodes = [[nodeid, edonid] for nodeid, edonid, _ in edges]
            incident_nodes = set(sum(incident_nodes), list())
            for nodeid in incident_nodes:
                incident = graph.incident_edges(nodeid)
                ignorable = graph.edges_with_label_in(cls.IGNORE, incident)
                if len(ignorable) == len(incident):
                    graph.delete_node(nodeid)


class GraphEmbedsEgg(object):
    def __init__(self, run, code, locus=None, legacy=None):
        self.run = run
        self.code = DFSCode.remulted(code, self.run.domain)
        self.locus = list() if locus is None else locus
        self.legacy = (list(), list(), list()) if legacy is None else legacy
        self.embeds = None

    def warm(self, run=None):
        if not self.run:
            if not run:
                raise ValueError("'run' attribute not initiated")
            self.run = run
        self.embeds = self.run.graphs.embed_code(self.code)

    def hatch(self):
        if self.specificity_early_termination():
            return GraphEmbeds(self, dict(), dict(), list())
        treepres, backpres, hintpres = self.legacy
        treeexts = self.extensions(treemode=True, oldfrags=treepres)
        hintexts = self.external_eqoccce(hintpres, treeexts)
        backexts = self.extensions(treemode=False, oldfrags=backpres)
        return GraphEmbeds(self, treeexts, backexts, hintexts)

    def external_eqoccce(self, hintpres, treeexts):
        exteqoccs = list()
        selfbeds = len(self.embeds)
        if len(hintpres):
            for hint in hintpres:
                if hint in treeexts and len(treeexts[hint]) == selfbeds:
                        exteqoccs.append(hint)
        if not exteqoccs:
            for treefrag, dfsgwise in treeexts.items():
                if len(dfsgwise) == selfbeds:
                    exteqoccs.append(treefrag)
        return exteqoccs

    def specificity_early_termination(self):
        nodelabrange, edgelabrange = self.specific_label_range()
        for ltype, lrange in [('node', nodelabrange), ('edge', edgelabrange)]:
            for ltag, labls in enumerate(lrange):
                for labl in labls[1:]:
                    matches = self.match_label(labl, ltag, self.embeds, ltype)
                    if len(matches) == len(self.embeds):
                        print "Early Termination"
                        return True
        return False

    def specific_label_range(self):
        nodelabs = [[label] for label in self.embeds[0].code._tags]
        edgelabs = [[frag.edgelab] for frag in self.code if not frag.isnode]
        for index, labels in enumerate([nodelabs, edgelabs]):
            labelslist = self.run.domain[index]
            ltype = 'edge' if index else 'node'
            for atag, taglabels in enumerate(labels):
                main_label = taglabels[0]
                labelneeded = labelslist[main_label.name]
                for speclab in labelslist.specific_to(main_label):
                    beds = self.match_label(speclab, atag, self.embeds, ltype)
                    if len({abed.id for abed in beds}) >= self.run.minsup:
                        taglabels.append(speclab)
        return nodelabs, edgelabs

    def specific_label_range1(self):
        nodelabs = [[label] for label in self.embeds[0].code._tags]
        edgelabs = [[frag.edgelab] for frag in self.code if not frag.isnode]
        # beds = [(mbed.code._tags, list(list(zip(*mbed.edges)[2])), mbed.id) for mbed in self.embeds]
        # for index, typelabs in enumerate([nodelabs, edgelabs]):
        #     labelslist = self.run.domain[index]
        #     indexbeds = [(abed[index], abed[-1]) for abed in beds]
        #     for ltag, (baselab,) in enumerate(typelabs):
        #         if baselab.is_specific:
        #             continue
        #         count = defaultdict(set)
        #         for mbed in indexbeds:
        #             count[mbed[ltag].value].add(mbed.id)
        #         for speclab in labelslist.specific_to(baselab):
        #             for eqval in speclab.eqvals - {speclab.value}:
        #                 count[speclab.value].update(count[eqval])
        #         for labl, byid in count.items():
        #             if len(byid) >= self.run.minsup:
        #                 typelabs[ltag].append(labl)
        return nodelabs, edgelabs

    @classmethod
    def matching_subset(cls, bedlist, frag, edgex):
        if not frag.isnode:
            embeds = cls.match_label(frag.edgelab, edgex, bedlist, 'edge')
            if frag.intree:
                embeds = cls.match_label(frag.edonlab, frag.edontag, embeds, 'node')
            return embeds
        return cls.match_label(frag[1], frag[0], bedlist, 'node')

    @staticmethod
    def match_label(label, tagex, bedlist, ltype='node'):
        embeds = list()
        saved = None
        for embed in bedlist:
            if ltype == 'node':
                nodeid = embed.nodeids[tagex]
                saved = nodeid
                embed_label = embed._graphdata.node_labels[nodeid]
            else:
                embed_label = embed.edges[tagex][2]
            if label.covers(embed_label):
                embeds.append(embed)
        return embeds

    def extensions(self, treemode, oldfrags):
        kwargs = {'treeedges': treemode, 'backedges': not treemode}
        isatroot = len(self.run.base) == len(self.code) # <=> all are rootnodes
        extensions = defaultdict(list) # lambda: defaultdict(list))
        for dfsgen in self.embeds: # extend only along the latest new node
            lastnodes = dfsgen.nodeids[:] if isatroot else [dfsgen.nodeids[-1]]
            for child in dfsgen.children(lastnodes, **kwargs):
                extensions[child].append(dfsgen)
        extensions = self.add_new_altcodes(extensions, treemode)
        extensions = self.extend_inherited(extensions, treemode, oldfrags)
        extensions = self.abort_infrequent(extensions)
        if treemode:
            extensions = self.abort_at_pattern(extensions) #TODO : decide on indentation
        return extensions

    def make_altkids(self, child, treemode):
        genfrags = list()
        genlabls = [[child[4]], [child[3]]]
        label = self.run.domain.edge[child[3]]
        genlabls[1] += self.run.domain.edge.general_to(label)
        if treemode:
            label = self.run.domain.node[child[4]]
            genlabls[0] += self.run.domain.node.general_to(label)
        for nlab, elab in product(genlabls[0], genlabls[1]):
            genfrags.append(child.change_labels(child[2], elab, nlab))
        return genfrags[1:]

    def add_new_altcodes(self, extensions, treemode):
        for newfrag in list(extensions):
            generals = self.make_altkids(newfrag, treemode)
            for general in generals:
                extensions[general] += extensions[newfrag]
            for frag1, frag2 in combinations(generals, 2):
                if len(extensions[frag1]) == len(extensions[frag2]):
                    if frag1.covers_labels(frag2):
                        extensions.pop(frag1, None)
                    elif frag2.covers_labels(frag1):
                        extensions.pop(frag2, None)
        return extensions

    def extend_inherited(self, extensions, treemode, oldfrags):
        # enlist inherited frequent edges
        newfrags = list()
        for oldfrag in oldfrags:
            if treemode:
                oldfrag = (oldfrag[0], oldfrag[1]+1) + oldfrag[2:]
            oldfrag = DFSCodeFragment.remult(oldfrag, self.run.domain)
            newfrags.append(oldfrag)
        # check embeddings of inherited frequent edges
        for dfsgen in self.embeds:
            for newfrag in newfrags:
                if dfsgen.step(newfrag, extend=False):
                    extensions[newfrag.unmult()].append(dfsgen)
        return extensions

    def abort_infrequent(self, extensions):
        # ensure only frequent extensions
        for code in list(extensions.keys()):
            if DFSGen.graph_count(extensions[code]) < self.run.minsup:
                del extensions[code]
        return extensions

    def abort_at_pattern(self, extensions):
        newfrags = list(extensions)
        for newfrag in newfrags:
            newfrag = DFSCodeFragment.remult(newfrag, self.run.domain)
            newcode = self.code.copy_append(newfrag)
            newtsgf = RNAStructureGraph.from_code(newcode)
            for abortat in self.run.abortats:
                if len(newtsgf.embed_code(abortat)):
                    del extensions[newfrag]
                    break
        return extensions
    

class GraphEmbeds(object):

    # def hatch_from_egg(self, eggself, treeexts, backexts, hintexts):
    def __init__(self, eggself, treeexts, backexts, hintexts):
        self.run = eggself.run
        self.code = eggself.code
        self.base = DFSCode(eggself.code[:len(self.run.base)])
        self.locus = eggself.locus
        self.embeds = eggself.embeds
        self.treeexts = treeexts
        self.backexts = backexts
        self.hintexts = hintexts
        assert(isinstance(self.base[0][1], Label))

    @property
    def occurrence(self):
        support = len(set(embed.id for embed in self.embeds))
        return (len(self.embeds), support)

    def childtrees(self):
        # gSpan redundancy removal
        newcodes = list()
        for newfrag in self.treeexts:
            newfrag = DFSCodeFragment.remult(newfrag, self.run.domain)
            newcode = self.code.copy_append(newfrag)
            newtsgf = RNAStructureGraph.from_code(newcode)
            mincode = newtsgf.mindfscode(self.base).unmulted
            if mincode == newcode.unmulted: # will not prevent alttrees in other paths
                newcodes.append((newfrag, self.run.encode(newcode)))

        legacy = (self.treeexts.keys(), self.backexts.keys(), self.hintexts)
        childtrees = list()
        broodsize = len(newcodes)
        for whichth, (newfrag, newcode) in enumerate(sorted(newcodes)):
            locus = self.locus + [(whichth+1, broodsize)]
            newchild = GraphEmbedsEgg(self.run, newcode, locus, legacy)
            childtrees.append(newchild)
        return childtrees

    def maximal_subgraph_expansion(self):
        if len(self.hintexts): # external edge pruning
            return list()

        associative = set() # tail shrinking
        for loop in self.backexts:
            if len(self.backexts[loop]) == len(self.embeds):
                associative.add(loop)

        # DONOT use lethal associative edge set for DFSCode canonical trees
        # mincode = self.plus_loops(associative)
        # if not self.in_class(mincode):
        #     return list()
        
        remainder = list(set(self.backexts).difference(associative))
        loopsets = sorted(self.search_graphs(associative, remainder), key=len)
        maxgraphs = list()
        for index, set1 in enumerate(loopsets):
            mincode = self.plus_loops(set1)
            if self.in_class(mincode): # remove not in_class
                for set2 in loopsets[index+1:]:
                    if len(set1) < len(set2) and \
                      not len(set1.difference(set2)): #TODO : issubset?
                        break # remove non-maximal
                else:
                    if not self.graph_associative_external_edge(set1):
                        maxgraphs.append(self.run.encode(mincode))
        return maxgraphs

    def search_graphs(self, loops, candidates):
        # maximality is not ensured by search graphs, only frequency
        if not candidates:
            return [loops] # recursion terminal case

        loopsbeds = set(self.embeds)
        for loopfrag in loops:
            loopsbeds = loopsbeds.intersection(self.backexts[loopfrag])

        # bottom-up pruning
        embedsall = set(loopsbeds)
        for candidate in candidates:
            embedsall = embedsall.intersection(self.backexts[candidate])
        if len(set([embed.id for embed in embedsall])) >= self.run.minsup:
            return [loops.union(candidates)]

        # dynamic reordering
        frags = list()
        for candidate in candidates:
            candbeds = loopsbeds.intersection(self.backexts[candidate])
            support = len(set([candbed.id for candbed in candbeds]))
            if support >= self.run.minsup: # removing infrequent ones
                frags.append((support, len(candbeds), candidate))
        if not frags:
            return [loops]

        frags = [frag for _, _, frag in sorted(frags)] # sorting by support
        loopsets = list()
        for index, frag in enumerate(frags):
            loopsets += self.search_graphs(loops.union([frag]), frags[index+1:])
        return loopsets

    def plus_loops(self, loopset): #TODO: possible misses when edge is mult
        for loop1, loop2 in combinations(list(loopset), 2):
            frag1 = DFSCodeFragment.remult(loop1, self.run.domain)
            frag2 = DFSCodeFragment.remult(loop2, self.run.domain)
            if frag1.compares_with(frag2):
                if frag1.covers_labels(frag2) and loop1 in loopset:
                    loopset.remove(loop1)
                if frag2.covers_labels(frag1) and loop2 in loopset:
                    loopset.remove(loop2)

        withloops = list(self.code.unmulted) + list(loopset)
        withloops = DFSCode.remulted(withloops, self.run.domain)
        newgraph = RNAStructureGraph.from_code(withloops)
        return newgraph.mindfscode(self.base).unmulted

    def in_class(self, graphmincode):
        graphmincode = DFSCode(graphmincode)
        newtree = graphmincode[:len(self.base)]
        for frag in graphmincode[len(self.base):]:
            if frag.intree: # i.e. only front edge
                newtree.append(frag)
        newtree = DFSCode.remulted(newtree, self.run.domain)
        newtsgf = RNAStructureGraph.from_code(newtree)
        treecode = newtsgf.mindfscode(self.base).unmulted #TODO: why unmulted?
        return tuple(treecode) == self.code.unmulted

    def graph_associative_external_edge(self, loopset):
        embedsall = set(self.embeds)
        for loop in loopset:
            embedsall = embedsall.intersection(self.backexts[loop])
        for treefrag, treeexts in self.treeexts.items():
            if embedsall == set(treeexts):
                return True
        return False


class Run(object):
    POPDELAY = 0.1
    
    def __init__(self, base, graphs, minsup, abortats, domain):
        self.base = base
        self.graphs = graphs
        self.minsup = minsup
        self.abortats = abortats
        self.domain = domain

    def encode(self, newcode):
        return DFSCode(newcode)

    def embed(self, code):
        return GraphEmbedsEgg(self, code)


class SharedLogger(object):

    LOGDELAY = 10.0

    def __init__(self, todo, done, last, lock, file, level, procid):
        self._todo = todo   # lock0
        self._done = done   # lock0
        self._last = last   # lock1
        self._lock = lock
        self._file = file
        self._level = level
        self.procid = procid
        self._count = 0

    @classmethod
    def create(cls, inittodo, logfile, procid):
        todo = mp.RawValue(ctypes.c_ulonglong, inittodo)
        done = mp.RawValue(ctypes.c_ulonglong, 0)
        last = mp.RawValue(ctypes.c_double, time.time())
        lock = (mp.Lock(), mp.Lock())
        file = logfile
        level = mp.RawArray(ctypes.c_ulonglong, [inittodo]+[0]*99)
        return cls(todo, done, last, lock, file, level, procid)

    def done_one(self):
        with self._lock[0]:
            self._done.value += 1

    def todo_more(self, newtodo, level):
        if newtodo > 0:
            with self._lock[0]:
                self._todo.value += newtodo
            with self._lock[1]:
                self._level[level] += newtodo

    def status(self):
        todo, done = self._todo.value, self._done.value
        if todo > done:
            return todo, done, todo - done
        with self._lock[0]:
            todo, done = self._todo.value, self._done.value
        return todo, done, todo - done

    def level(self):
        levellist = list(self._level)
        if levellist[-1] == 0:
            while levellist[-2] == 0:
                levellist.pop()
        return levellist

    def time_to_log(self):
        if self.procid[1] > 1:
            self._count = (self._count % self.procid[1]) + 1
            if self._count != self.procid[0]:
                return False
        timenow = float(time.time())
        if timenow - self._last.value > self.LOGDELAY:
            self._last.value = timenow
            return True
        return False

    def is_on(self):
        if self.status()[2] == 0:
            return False
        return True

    def update_logfile(self):
        if self._file:
            toout = list(self.status()) + self.level()
            toout = " ".join([str(count).rjust(6) for count in toout])
            with open(self._file, 'a') as logfile:
                logfile.write(str(self.procid[0]).ljust(4)+":"+toout+"\n")
            timenow = float(time.time())
            self._last.value = timenow


class CSpinStarter(object):
    
    def __init__(self, run):
        baseegg = run.embed(run.base)
        baseegg.warm()
        self.run = run
        self.logs = SharedLogger.create(1, None, 0)
        self.labrange = baseegg.specific_label_range()
        self.altbases = list()
        self.altlist = mp.Queue()
        self.altlist.put((0, None, baseegg.embeds))

    def alternate_labels(self, nowlen, labels):
        frag = self.run.base[nowlen]
        nlabrange, elabrange = self.labrange
        if frag.isnode:
            return [[[nlab], []] for nlab in nlabrange[frag.nodetag]]
        elabs = elabrange[nowlen-1]
        if not frag.intree:
            return [[labels[0], labels[1] + [elab]] for elab in elabs]
        nlabs = [frag[4]] if frag[0] > frag[1] else nlabrange[frag[1]]
        altlabels = [[labels[0]+[nlab], labels[1]+[elab]] 
                     for nlab, elab in product(nlabs, elabs)]
        return altlabels

    def frequent_specifics(self, nowlen, altlabls, bedlist):
        altcodes = dict()
        base = self.run.base
        for index1, labels1 in enumerate(altlabls):
            code1 = base.relabel(*labels1)
            frag1 = code1[nowlen]
            beds1 = GraphEmbedsEgg.matching_subset(bedlist, frag1, nowlen-1)
            strlabels1 = [" ".join([i.name for i in k]) for k in labels1]
            if DFSGen.graph_count(beds1) >= self.run.minsup:
                for index2, (_, code2, beds2) in list(altcodes.items()):
                    if len(beds1) == len(beds2):
                        frag2 = code2[nowlen]
                        if frag1.compares_with(frag2):
                            if frag1.covers_labels(frag2):
                                break
                            elif frag2.covers_labels(frag1):
                                del altcodes[index2]
                else:
                    altcodes[index1] = (labels1, code1, beds1)
        return altcodes

    def add_new_altbase(self, altbase, occurrence): # TODO: becareful with mindfscode
        unmult = tuple(altbase.unmulted)
        codegf = RNAStructureGraph.from_code(DFSCode(unmult))
        minnow = tuple(codegf.mindfscode())
        self.altbases.append((minnow, unmult, occurrence))

    def pop_altlist(self):
        while self.logs.is_on():
            try:
                depth, labs, beds = self.altlist.get(True, 1.0)
            except Empty:
                time.sleep(self.run.POPDELAY)
            else:
                yield depth, labs, beds

    def base_trees(self):
        base = self.run.base
        alts = {unmulted: occce for _, unmulted, occce in self.altbases}
        print "Alternate Bases\t", len(alts)
        for key1, occ1 in list(alts.items()):
            if key1 in alts:
                graph1 = RNAStructureGraph.from_code(self.run.encode(key1))
                for key2, occ2 in list(alts.items()):
                    if key1 != key2 and occ1 == occ2 and\
                      len(graph1.embed_code(self.run.encode(key2))):
                        del alts[key2]
        inittrees = [self.run.embed(self.run.base)]
        alts.pop(base.unmulted, None)
        for index, altbase in enumerate(sorted(alts)):
            altbase = self.run.encode(altbase)
            inittrees.append(self.run.embed(altbase))
        return inittrees

    def update(self, results):
        alts = dict()
        for mincode, unmult, occurrence in results:
            if mincode in alts:
                unmult = min([alts[mincode][0], unmult])
            alts[mincode] = (unmult, occurrence)
        self.altbases = [(code,) + value for code, value in alts.items()]


def explore_bases(args):
    shared.logs.procid = args
    print "exploring bases with", args[0], os.getpid()
    for len0, labs0, beds0 in shared.pop_altlist():
        altlabels = shared.alternate_labels(len0, labs0)
        altcodes = shared.frequent_specifics(len0, altlabels, beds0)
        for _, (labsi, codei, bedsi) in sorted(altcodes.items()):
            if len0+1 == len(shared.run.base):
                occurrence = (len(bedsi), DFSGen.graph_count(bedsi))
                shared.add_new_altbase(codei, occurrence)
            else:
                labsi = labsi
                shared.altlist.put((len0+1, labsi, bedsi))
                shared.logs.todo_more(1, len0+1)
        shared.logs.done_one()
    return shared.altbases


class CSpinRun(object):

    DELETABLE = ["N", "p"]
    ATTACHABLE = ["N", "R", "Y", "A", "C", "G", "U", "p", "WWC"]

    def __init__(self, run, trees, logfile):
        self.run = run
        self.logs = SharedLogger.create(len(trees), logfile, 0)
        self.treelist = mp.Queue()
        self.graphlist = list()
        for tree in trees:
            self.treelist.put(tree)
  
    def trim(self):
        trimmed = set()
        for dfscode in self.graphlist:
            unmult = dfscode.unmulted
            graph = RNAStructureGraph.from_code(DFSCode(dfscode.unmulted))
            nodes = [node for node in graph.nodes if node[1] in self.DELETABLE]
            terminals = None
            while terminals is None or terminals:
                terminals = list()
                for node in nodes:
                    elist = [edge for edge in graph.edges if node in edge]
                    edons = set([edge[1-edge.index(node)] for edge in elist])
                    if len(edons) == 1 and list(edons)[0][1] in self.ATTACHABLE:
                        terminals.append(node)
                for node in terminals:
                    graph.delete_node(node)
                    nodes.remove(node)
            trimmed.add(tuple(graph.mindfscode(unmult[:len(self.run.base)])))
        self.graphlist = [self.run.encode(code) for code in trimmed]

    def maximal(self):
        graphs = {code.unmulted: (RNAStructureGraph.from_code(code), code)\
                                        for code in self.graphlist}
        for key1, (graph1, code1) in sorted(graphs.items()):
            for key2, (graph2, code2) in list(graphs.items()):
                if key1 != key2 and len(code1) <= len(code2)\
                  and len(graph2.embed_code(code1)):
                    del graphs[key1]
                    break
        self.graphlist = [code for _, (_, code) in graphs.items()]

    def pop_tree(self):
        while self.logs.is_on():
            try:
                initegg = self.treelist.get(True, 1.0)
            except Empty:
                time.sleep(self.run.POPDELAY)
            else:
                if self.logs.time_to_log():
                    self.logs.update_logfile()
                initegg.warm(self.run)
                yield initegg.hatch()
        self.logs.update_logfile()

    def add_new_trees(self, level, kidtrees, inittree):
        for kidtree in kidtrees:
            self.treelist.put(kidtree)
        self.logs.todo_more(len(kidtrees), level+1)
        self.logs.done_one()

    def update(self, results):
        self.graphlist = [DFSCode(code) for code in results]


def explore_trees(args):
    def explore():
        shared.logs.procid = args
        print "exploring trees with", args[0], os.getpid()
        for inittree in shared.pop_tree():
            with DomainChange(inittree.run.domain):
                kidtrees = inittree.childtrees()
                kidtrees.sort(key=lambda akid: akid.code.unmulted)
                shared.add_new_trees(len(inittree.locus), kidtrees, inittree)
                maxgraphs = inittree.maximal_subgraph_expansion()
                shared.graphlist += maxgraphs
    logfile = shared.logs._file.replace(".run", ".time")+str(args[0])
    cProfile.runctx('explore()', globals(), locals(), logfile)
    return shared.graphlist


def initmaker(toshare):
    global shared
    shared = toshare


class StopCodeFragment(DFSCodeFragment):

    def covers(self, fles):
        return self.compares_with(fles) and\
               self.covers_labels(fles) and fles.covers_labels(self)


class StopCode(DFSCode):

    def __init__(self, code):
        if isinstance(code, DFSCode):
            code = [StopCodeFragment(frag) for frag in code]
        super(StopCode, self).__init__(code)


class CSpinManager(object):


    def __init__(self, graphs, logfile, maxprocs=1):
        self.runs = 0
        self.graphs = GraphList(graphs)
        self.logfile = os.path.abspath(logfile)
        self.maxprocs = maxprocs
        self.domain = self.graphs.get_domain()

    def create_abort_codes(self, stoplen, domain):
        stopcodes = list()
        # 
        symmetric = 1
        graph0 = RNAStructureGraph("")
        sequence = [RNAStructureGraph.N.name]*stoplen
        sequence = list(zip(map(str, range(stoplen)), sequence))
        graph0.add_sequence(sequence)
        engrapher = RNAStructureEngrapher()
        
        p = graph0.PHOSPHATE
        phosnodes = [("%s.%s" % (nucid, p.name), p) for nucid, _ in sequence]
        if symmetric:
            graph0.delete_node("0.%s" % p.name)

        edgelabels = [domain.edge[label.name] for label in graph0.SEQUENCE, graph0.ISPARTOF]
        edgelabels = (edgelabels,)*(len(sequence)-symmetric)
        for edgeset in product(*edgelabels):
            graphi = graph0.copy("")
            for index, label in enumerate(edgeset, symmetric):
                nucid, phosid = "%d" % index, "%d.P" % index
                if label.name == RNAStructureGraph.ISPARTOF.name:
                    graphi.edges.remove((nucid, phosid, label))
                else:
                    graphi.edges.remove((phosid, nucid, label))
            engrapher.graph = graphi
            mincode = graphi.mindfscode()
            stopcodes.append(StopCode(mincode))
        return stopcodes

    def distribute_work(self, method, toshare):
        if self.maxprocs > 2:
            pool = Pool(self.maxprocs-1, initializer=initmaker, initargs=(toshare,))
            args = [(indx, self.maxprocs) for indx in range(1, self.maxprocs)]
            results = sum(pool.map(method, args), list())
            pool.close()
            pool.join()
        else:
            initmaker(toshare)
            results = method((0, 1))
        toshare.update(results)

    def add_general_labels(self, multlabs):
        nodelabels, edgelabels = self.domain
        for index in (0, 1):
            multlabs[index] = sorted(set(multlabs[index]))
            labelslist = RNAStructureGraph.LABELS[index]
            multlabs[index] = labelslist.create_string_labels(multlabs[index])
        nodelabels = nodelabels.create_superset(multlabs[0])
        edgelabels = edgelabels.create_superset(multlabs[1])
        return Domain(nodelabels, edgelabels)

    def update_multlabs(self, multlabs, base):
        for frag in base:
            if frag.nodelab.is_general:
                multlabs[0].append(frag.nodelab.name)
            if not frag.isnode:
                if frag.edonlab.is_general:
                    multlabs[0].append(frag.edonlab.name)
                if frag.edgelab.is_general:
                    multlabs[1].append(frag.edgelab.name)

    def run(self, base, minsup, multlabs, stoplen):
        self.update_multlabs(multlabs, base)
        domain = self.add_general_labels(multlabs)
        with DomainChange(domain) as dchange:
            runbase = DFSCode(base).copy()
            dchange.code(runbase)
            abortats = self.create_abort_codes(stoplen, domain)
            self.runs += 1
            graphs = list()
            for graph in self.graphs:
                changed = graph.copy()
                dchange.graph(changed)
                graphs.append(changed)
            graphs = GraphList(graphs)
        return Run(runbase, graphs, minsup, abortats, domain)

    def start(self, base, minsup, multlabs, stoplen):
        logfile = "%s.run%s" % (self.logfile, self.runs)
        thisrun = self.run(base, minsup, multlabs, stoplen)
        starter = CSpinStarter(thisrun)
        self.distribute_work(explore_bases, starter)
        spinrun = CSpinRun(thisrun, starter.base_trees(), logfile)
        self.distribute_work(explore_trees, spinrun)
        spinrun.trim()
        spinrun.maximal()
        supports = dict()
        for graphcode in spinrun.graphlist:
            embeds = self.graphs.embed_code(graphcode)
            support = DFSGen.graph_count(embeds)
            supports[graphcode.unmulted] = (support, len(embeds))
        return supports


