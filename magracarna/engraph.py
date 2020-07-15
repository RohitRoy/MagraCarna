import os
import re
import sys
import time
import traceback
import numpy as np

try:
    import graph_tool as gt
    import graph_tool.draw as draw
except:
    pass

from itertools import product, combinations, combinations_with_replacement
from collections import defaultdict, deque, namedtuple

from .engrid import Residue
from .aminoacids import aminoacids
from .heteroatoms import hetresdus, someothers
from .nucleotides import allnucs, purines, pyrimidines, RNAnucs,\
                         variants, phosphate, atoms as nucatoms,\
                         edges as nucedges

MAX_SEQCE_DIST = 7

class Label(object):

    def __init__(self, value, eqvals=None, name=None):
        if eqvals:
            eqvals_, eqvals = list(eqvals), list()
            for eqval in eqvals_:
                eqvals += eqval.eqvals if self.islabel(eqval) else [eqval]
        self.eqvals = frozenset(eqvals if eqvals else [value])
        self.value = value
        self.name = str(value if name is None else name)

    @classmethod
    def islabel(cls, label):
        return isinstance(label, cls)

    def __str__(self):
        return self.name

    def __repr__(self):
        return "Label(%s)" % self.name

    @property
    def is_general(self):
        return len(self.eqvals) > 1

    @property
    def is_specific(self):
        return self.eqvals == {self.value}

    @property
    def is_redundant(self):
        return not (self.is_general or self.is_specific)

    def __hash__(self):
        return id(self)
        return tuple([self.value]+sorted(self.eqvals))

    def __eq__(self, fles):
        return self.value == fles.value

    def __ne__(self, fles):
        return self.value != fles.value

    def covers(self, fles):
        return fles.eqvals.issubset(self.eqvals)

    def __lt__(self, fles):
        return self.name < fles.name

    def __gt__(self, fles):
        return self.name > fles.name

    def __ge__(self, fles):
        return self.name >= fles.name

    def __le__(self, fles):
        return self.name <= fles.name


class StrLabel(Label):
    def __init__(self, value, eqvals=None):
        for aval in [value] + (list(eqvals) if eqvals else []):
            aval = aval.value if Label.islabel(aval) else aval
            if not isinstance(aval, str):
                raise ValueError("StrLabel accepts only 'str' values")
        super(StrLabel, self).__init__(value, eqvals, value)


class LabelsList(object):

    def __init__(self, strlabels):
        """Creates LabelsList from iterable of StrLabels."""
        arguments = [args[::-1] for args in enumerate(strlabels)]
        self.names = tuple(self.enlist_label_names(strlabels))
        self.labels = tuple([self._create_label(*args) for args in arguments])
        for label in self.labels:
            self._check_sanity(label)
        self._general_to = dict()
        self._specific_to = dict()
        self._name_to_index = {slab.name: index for slab, index in arguments}

    def general_to(self, label):
        if 0 <= label.value < len(self):
            if label.value not in self._general_to:
                generals = list()
                for lebal in self.labels:
                    if lebal.covers(label):
                        generals.append(lebal.value)
                generals.remove(label.value)
                generals.remove(self["*"].value)
                self._general_to[label.value] = generals
            generals = self._general_to[label.value]
            return [self.labels[index] for index in generals]
        raise ValueError("Label %s is not in LabelsList" % (label))

    def specific_to(self, label):
        if 0 <= label.value < len(self):
            if label.value not in self._specific_to:
                specifics = list()
                for lebal in self.labels:
                    if label.covers(lebal):
                        specifics.append(lebal.value)
                specifics.remove(label.value)
                self._specific_to[label.value] = specifics
            specifics = self._specific_to[label.value]
            return [self.labels[index] for index in specifics]
        raise ValueError("Label %s is not in LabelsList" % (label))

    def _check_sanity(self, label):
        for eqval in label.eqvals:
            if not self.labels[eqval].is_specific:
                message = "Label '%s' refers non-specific labels"
                raise ValueError(message % label.name)

    def _create_label(self, strlabel, value=None):
        """Creates LabelsList compatible label from given StrLabel """
        if not isinstance(strlabel, StrLabel):
            raise ValueError("Labels are created only from StrLabel")
        label_eqvals = list()
        label_value = strlabel.value if value is None else value
        for eqval in strlabel.eqvals:
            if eqval in self.names:
                label_eqvals.append(self.names.index(eqval))
        if not label_eqvals:
            raise ValueError("LabelsList labels can't define '%s'" % strlabel)
        return Label(label_value, label_eqvals, strlabel.name)

    def create_label(self, strlabel, value=None):
        label = self._create_label(strlabel, value)
        self._check_sanity(label)
        return label

    def create_subset(self, labels):
        # subset is not comparable to the superset
        names = self.enlist_label_names(labels)
        relabels = list()
        for label in self.labels:
            if label.name in names:
                relabels.append(self.create_string_label(label.value))
        if len(relabels) != len(names):
            raise ValueError("Labels are not subset of LabelsList")
        return LabelsList(relabels)

    def create_superset(self, strlabels):
        # superset is comparable to subset
        strlabelslist = list()
        for index in range(len(self)):
            strlabelslist.append(self.create_string_label(index))
        strlabelslist += strlabels
        return type(self)(strlabelslist)

    def create_string_label(self, index):
        try:
            label = self.labels[index]
        except IndexError:
            raise IndexError("Label index %s out of range" % str(index))
        if label.is_specific:
            return StrLabel(label.name)
        eqvals = [self.names[index] for index in label.eqvals]
        return StrLabel(label.name, eqvals)

    def create_string_labels(self, names):
        # not an actual LabelsList object, just a list
        strlabels = list()
        for name in names:
            strlabel = self.create_string_label(self[name].value)
            strlabels.append(strlabel)
        return strlabels

    @staticmethod
    def enlist_label_names(labels):
        names = [str(label) for label in labels]
        if len(set(names)) < len(names):
            for name in set(names):
                if names.count(name) > 1:
                    repeat_label = name
                    break
            raise ValueError("Label '%s' is repeated" % repeat_label)
        return names

    def __getitem__(self, label):
        if isinstance(label, Label):
            return self[label.name]
        elif isinstance(label, str):
            try:
                return self[self._name_to_index[label]]
            except:
                raise ValueError("No Label named %s in LabelList" % label)
        elif isinstance(label, int):
            try:
                return self.labels[label]
            except:
                raise IndexError("Index '%d' is out of range" % label)
        raise ValueError("Retrieves labels from only int, str or Label keys")

    def index(self, label):
        try:
            label = label.name if isinstance(label, Label) else label
            return self.names.index(label)
        except ValueError:
            raise ValueError("No Label named '%s' in LabelList" % label)

    def __contains__(self, label):
        label = label.name if isinstance(label, Label) else label
        return label in self.names

    def __iter__(self):
        return iter(self.labels)

    def __len__(self):
        return len(self.labels)


Domain = namedtuple('Domain', ['node', 'edge'])


class DFSCodeFragment(tuple):

    class InvalidFragmentError(ValueError):
        MESG = 'invalid DFSCodeFragment "%s"'
        def __init__(self, fragment):
            self.fragment = fragment
            message = self.MESG % (fragment)
            cls = DFSCodeFragment.InvalidFragmentError
            super(cls, self).__init__(self, message)

    def __init__(self, fragment):
        self.isnode = len(self) == 2
        self.intree = self.isnode or self[0] < self[1]

    @property
    def nodetag(self):
        return self[0]

    @property
    def nodelab(self):
        return self[1] if self.isnode else self[2]

    @property
    def edontag(self):
        return None if self.isnode else self[1]

    @property
    def edonlab(self):
        return None if self.isnode else self[4]

    @property
    def edgelab(self):
        return None if self.isnode else self[3]

    def check_validity(self):
        if len(self) == 2:
            if not isinstance(self[0], int) or not Label.islabel(self[1]):
                raise self.InvalidFragmentError(str(self))
        elif len(self) == 5:
            for value, thetype in zip(self, [int, int, Label, Label, Label]):
                if not isinstance(value, thetype):
                    raise self.InvalidFragmentError(str(self))
        else:
            raise self.InvalidFragmentError(str(self))

    def unmult(self):
        if self.isnode:
            return (self[0], self[1].name)
        return (self[0], self[1], self[2].name, self[3].name, self[4].name)

    @classmethod
    def remult(cls, fragment, domain):
        if len(fragment) == 2:
            return cls((fragment[0], domain.node[fragment[1]]))
        nodelab = domain.node[fragment[2]]
        edgelab = domain.edge[fragment[3]]
        edonlab = domain.node[fragment[4]]
        return cls((fragment[0], fragment[1], nodelab, edgelab, edonlab))

    def compares_with(self, fles):
        if self.isnode == fles.isnode and self[0] == fles[0]\
          and (self.isnode or self[1] == fles[1]):
            return True
        return False

    def covers_labels(self, fles):
        if fles.isnode: # fles is DFSCodeFragment and NodeFragment
            return self[1].covers(fles[1]) # self is NodeFragment
        return self[2].covers(fles[2]) and self[3].covers(fles[3])\
          and self[4].covers(fles[4]) # self, fles are EdgeFragment

    def covers(self, fles):
        return self.compares_with(fles) and self.covers_labels(fles)

    def coverage_with(self, fles):
        """ Returns whether it covers, iscoveredby given DFSCodeFragment."""
        if not self.compares_with(fles):
            return False, False
        return self.covers_labels(fles), fles.covers_labels(self)

    def change_labels(self, *labels):
        if self.isnode:
            return DFSCodeFragment((self[0], labels[0]))
        return DFSCodeFragment(self[:2] + tuple(labels[:3]))

    def change_tags(self, *tags):
        if self.isnode:
            return DFSCodeFragment((tags[0], self[1]))
        return DFSCodeFragment((tags[0], tags[1]) + self[2:])

    def lexical_form(self):
        if self.isnode:
            return (True, -(-1), self[0], "", "", self[1].name)
        fragment = (self.intree, -self[0], self[1])
        fragment += (self[2].name, self[3].name, self[4].name)
        return fragment

    def __lt__(self, fles):
        return self.lexical_form() < fles.lexical_form()

    def __gt__(self, fles):
        return self.lexical_form() > fles.lexical_form()

    def __le__(self, fles):
        return self.lexical_form() <= fles.lexical_form()

    def __ge__(self, fles):
        return self.lexical_form() >= fles.lexical_form()


class DFSCode(object):

    def __init__(self, code):
        if isinstance(code, tuple):
            code = DFSCode([DFSCodeFragment(frag) for frag in code])
        tags = dict()
        for fragment in code:
            tags[fragment.nodetag] = fragment.nodelab
            if not fragment.isnode:
                tags[fragment.edontag] = fragment.edonlab
        self._tags = [label for _, label in sorted(tags.items())]
        self._code = code

    @property
    def unmulted(self):
        return tuple([fragment.unmult() for fragment in self._code])

    @classmethod
    def remulted(cls, dfscode, domain):
        remulted = list()
        for frag in dfscode:
            remulted.append(DFSCodeFragment.remult(frag, domain))
        return cls(remulted)

    def __getitem__(self, index):
        value = self._code[index]
        if isinstance(value, slice):
            return DFSCode(value)
        return value

    def to_string(self, header="", indent="", indenthead=False, outfile=None):
        #TODO: change to "to graph string"
        code_string = "".join([indent+line+"\n" for line in self.to_lines()])
        out_string = (indent if indenthead else "") + header + "\n"
        out_string = ("\n"+out_string if header else out_string) + code_string
        if outfile:
            with open(outfile, 'a') as outfile:
                outfile.write(out_string+"\n")
        else:
            return out_string

    def to_lines(self):
        code_lines = list()
        for nodetag, nodelab in enumerate(self._tags):
            code_lines.append([str(nodetag), nodelab.name])
        for frag in self:
            nodetag, edontag, edgelab = frag.nodetag, frag.edontag, frag.edgelab
            if not frag.intree:
                code_lines[nodetag] += [str(edontag), frag.edgelab.name]
            elif not frag.isnode:
                code_lines[edontag] += [str(nodetag), frag.edgelab.name]
        return ["\t".join(line) for line in code_lines]

    @classmethod
    def from_lines(cls, lines, domain=None, base=None):
        domain = RNAGraph.LABELS if not domain else domain
        make = lambda frag: DFSCodeFragment.remult(frag, domain)
        code = base.copy() if base else cls([])
        for line in lines:
            line = [field.strip() for field in line.strip().split("\t")]
            nodetag, nodelab = int(line[0]), line[1]
            if not line[2:]: # implies isnode
                code.append(make((nodetag, nodelab)))
                continue
            edges = zip(line[2::2], line[3::2])
            for edgex, (edontag, edgelab) in enumerate(edges):
                edontag, edonlab = int(edontag), code._tags[int(edontag)]
                frag = (nodetag, edontag, nodelab, edgelab, edonlab) if edgex\
                       else (edontag, nodetag, edonlab, edgelab, nodelab)
                code.append(make(frag))
        return code

    def relabel(self, nlabs, elabs):
        nlabs, elabs = nlabs[:], elabs[:]
        # nodes = sorted(Graph.from_code(self).nodes)
        nlabs += [self._tags[ntag] for ntag in range(len(nlabs), len(self._tags))]
        elabs += [self[etag][3] for etag in range(len(elabs)+1, len(self))]
        created = list()
        for edgex, frag in enumerate(self):
            if frag.isnode:
                created.append(frag.change_labels(nlabs[frag[0]]))
            else:
                labels = nlabs[frag[0]], elabs[edgex-1], nlabs[frag[1]]
                created.append(frag.change_labels(*labels))
        return DFSCode(created)

    def append(self, newfrag):
        newfrag.check_validity()
        self._code.append(newfrag)
        if newfrag.intree:
            if newfrag.isnode:
                nodetag, nodelab = newfrag.nodetag, newfrag.nodelab
            else:
                nodetag, nodelab = newfrag.edontag, newfrag.edonlab
            if nodetag == len(self._tags):
                self._tags.append(nodelab)

    def copy(self):
        return DFSCode(self._code[:])

    def copy_append(self, newfrag):
        copied = self.copy()
        copied.append(newfrag)
        return copied

    def __iadd__(self, newcode):
        self._code += newcode

    def __lt__(self, fles):
        return self._code < fles._code

    def __le__(self, fles):
        return self._code <= fles._code

    def __gt__(self, fles):
        return self._code > fles._code

    def __ge__(self, fles):
        return self._code >= fles._code

    def __iter__(self):
        return iter(self._code)

    def __len__(self):
        return len(self._code)

    def __str__(self):
        return str(self._code)

    def __eq__(self, fles):
        if len(self) != len(fles):
            return False
        for fragA, fragB in zip(self._code, fles._code):
            if not fragA.covers(fragB) or not fragB.covers(fragA):
                return False
        return True


GraphData = namedtuple("GraphData", ["name", "node_labels", "incidence"])


class DFSGenerator(object):

    def __init__(self, code, nodeids, labels, edges, queue, graphdata):
        """
        Args:
            code: DFSCode upto previous step in the DFS
            nodeids: IDs of nodes embedded as per 'code'
            labels: node labels assigned as per 'code' #TODO: replace with code._tags
            edges: edges embedded as per 'code'
            queue: DFS queue as list of node IDs upto the
              previous step. Must be subset of 'nodeids'.
            graphdata: a GraphData object that completely
              describes the Graph on which DFS is run.
              This object should remain unchanged through
              the steps of DFS.
        """
        self.code = code
        self.queue = queue
        self.edges = edges
        self.labels = labels
        self.nodeids = nodeids
        self._node_exts = dict()
        self._graphdata = graphdata

    @property
    def id(self):
        return self._graphdata.name
        # return tuple(self.nodeids)

    def nodelabel(self, nodetag):
        nodeid = self.nodeids[nodetag]
        nodelabel = self._graphdata.node_labels[nodeid]
        return nodelabel

    @classmethod
    def graph_count(cls, dfsgens):
        return len({dfsg._graphdata.name for dfsg in dfsgens})

    def children(self, parents=None, backedges=True, treeedges=True):
        if parents == 'rightmost':
            parents = self.queue
        elif parents == 'uncovered':
            parents = set(self.nodeids) - set(self.queue)
        elif parents is None:
            parents = self.nodeids
        children = defaultdict(lambda: list())
        for nodeid in parents:
            node_exts = self._filtered_exts(nodeid, backedges, treeedges)
            for edge, frag in node_exts:
                children[frag].append((nodeid, edge))
        return children

    def minstep(self, extend=True, backedges=True, treeedges=True):
        #TODO: write explanation for how centering works without order
        queue = self.queue[:]
        while len(queue):
            nodeid = queue.pop()
            node_exts = self._filtered_exts(nodeid, backedges, treeedges)
            if node_exts:
                next_frag = min(frag for _, frag in node_exts)
                return self.step(next_frag, extend, backedges, treeedges)
        return list()

    def step(self, next_frag, extend=True, backedges=True, treeedges=True):
        extensions = list()
        next_nodetag = next_frag.nodetag
        if not (0 <= next_nodetag < len(self.nodeids))\
          or next_frag.nodelab.name != self.code._tags[next_nodetag].name:
            return extensions
        # 
        nodeid = self.nodeids[next_nodetag]
        for edge, frag in self._filtered_exts(nodeid, backedges, treeedges):
            if next_frag.covers(frag):
                extensions.append((nodeid, edge))
        if extend:
            return [self._extend(*args+(next_frag,)) for args in extensions]
        return extensions

    def as_embedding(self):
        embedding = list()
        for edgex, fragment in enumerate(self.code, -1):
            if fragment.isnode:
                embedding.append(self.nodeids[fragment.nodetag])
            else:
                embedding.append(self.edges[edgex][2].name)
                if fragment.intree:
                    embedding.append(self.nodeids[fragment.edontag])
        return embedding

    def _extend(self, nodeid, edge, next_frag):
        edonid = edge[1-edge.index(nodeid)]
        args = nodeids, labels, edges, queue = \
          self.nodeids[:], self.labels[:], self.edges[:], self.queue[:]
        if next_frag.intree: # i.e., a new node, edon, is being added
            queue.append(edonid)
            labels.append(next_frag.edonlab)
            nodeids.append(edonid)
        edges.append(edge)
        args = (self.code.copy_append(next_frag),) + args + (self._graphdata,)
        return DFSGenerator(*args)

    def _filtered_exts(self, nodeid, backedges, treeedges):
        filtered = list()
        for edge, frag in self._node_extensions(nodeid):
            if frag.intree and treeedges or not frag.intree and backedges:
                filtered.append((edge, frag))
        return filtered

    def _node_extensions(self, nodeid):
        if nodeid not in self._node_exts:
            incidence = self._graphdata.incidence[nodeid]
            node_elist = list(incidence.difference(self.edges))
            edge_frags = [self._fragment(nodeid, edge) for edge in node_elist]
            self._node_exts[nodeid] = list(zip(node_elist, edge_frags))
        return self._node_exts[nodeid]

    def _fragment(self, startid, edge):
        endid = edge[1-edge.index(startid)]
        nodetag = self._tag_nodeid(startid)
        edontag = self._tag_nodeid(endid)
        nodelab = self.code._tags[nodetag]
        if edontag <= nodetag:
            edonlab = self.code._tags[edontag]
        else:
            edonlab = self._graphdata.node_labels[endid]
        return DFSCodeFragment((nodetag, edontag, nodelab, edge[2], edonlab))

    def _tag_nodeid(self, nodeid):
        try:
            return self.nodeids.index(nodeid)
        except ValueError:
            return len(self.nodeids)


class Graph(object):

    def __init__(self, name=None):
        self.name = name
        self.nodes = dict()
        self.edges = set()

    def __str__(self):
        return self.to_string()

    def to_string(self):
        rows = "\n\t".join(["\t".join(arow) for arow in self.to_rows()])
        header = "\n%s\n\t" % (self.name if self.name else "")
        return "%s%s\n\n" % (header, rows)

    def to_rows(self, sorter=None):
        rows = defaultdict(list)
        for nodeid, nodelab in sorted(self.nodes.items()):
            rows[nodeid] = [str(nodeid), nodelab.name]
        for edge in sorted(self.edges):
            (nodeid, edonid), edgelab = sorted(edge[:2]), edge[2]
            rows[edonid] += [str(nodeid), edgelab.name]
        rows = [rows[nodeid] for nodeid in sorted(rows, key=sorter)]
        return rows

    def copy(self, name=None):
        graph = Graph(name)
        graph.nodes = self.nodes.copy()
        graph.edges = self.edges.copy()
        return graph

    def add_node(self, nodeid, nodelabel):
        if nodeid not in self.nodes:
            self.nodes[nodeid] = nodelabel
            return True
        return False

    def add_edge(self, nodeid1, nodeid2, edgelabel):
        """
        Args:
            node1, 
            node2: 2-tuples in format (nodeid, label)
            label: edge label
        """
        edge = (nodeid1, nodeid2, edgelabel)
        if edge not in self.edges:
            self.edges.add(edge)
            return True
        return False

    def delete_node(self, nodeid):
        if self.nodes.pop(nodeid, None) is not None:
            for edge in self.incident_edges(nodeid):
                self.edges.remove(edge)

    def incident_edges(self, nodeid):
        return {edge for edge in self.edges if nodeid in edge[:2]}

    def mindfscode(self, base=None):
        if base is None:
            base = DFSCode([DFSCodeFragment((0, min(self.nodes.values())))])
        # 
        mincode = None
        dfslist = self.embed_code(base) # self.embed_center(base)
        while dfslist:
            mincode = min([dfsg.code for dfsg in dfslist])
            minlist = [dfsg for dfsg in dfslist if dfsg.code == mincode]
            dfslist = sum((dfsg.minstep() for dfsg in minlist), list())
        return mincode

    def embed_code(self, code):
        newstep = self.embed_node(code[0][1])
        for frag in code[1:]:
            newstep = sum((dfsg.step(frag) for dfsg in newstep), list())
        if not newstep:
            return list()
        return newstep

    def match_code(self, code):
        queue = [self.embed_node(code[0][1])]
        while queue:
            if not queue[-1]:
                queue.pop()
            elif len(queue) == len(code):
                return queue[-1][-1]
            else:
                dfsg = queue[-1].pop()
                queue.append(dfsg.step(code[len(queue)]))
        return None

    def embed_node(self, nlabel):
        dfscode = DFSCode([DFSCodeFragment((0, nlabel)),])
        incidence = {node: self.incident_edges(node) for node in self.nodes}
        graphdata = GraphData(self.name, self.nodes, incidence)
        candidates = list()
        for nodeid, nodelab in self.nodes.items():
            if nlabel.covers(nodelab):
                dfsg = DFSGenerator(code=dfscode.copy(), nodeids=[nodeid],
                                    labels=[nlabel], edges=[], queue=[nodeid],
                                    graphdata=graphdata)
                candidates.append(dfsg)
        return candidates

    def edges_with_label_in(self, labels, from_edges=None):
        filtered = list()
        if from_edges is None:
            from_edges = self.edges
        for edge in from_edges:
            for label in labels:
                if label.covers(edge[2]):
                    filtered.append(edge)
                    break
        return filtered

    def nodes_with_label_in(self, labels, from_nodeids=None):
        filtered = list()
        if from_nodeids is None:
            from_nodeids = self.nodes.keys()
        for nodeid in from_nodeids:
            nodelab = self.nodes[nodeid]
            for label in labels:
                if label.covers(nodelab):
                    filtered.append(nodeid)
                    break
        return filtered

    @classmethod
    def from_code(cls, dfscode, name=None):
        graph = cls(name)
        for frag in dfscode:
            if frag.isnode:
                graph.add_node(frag[0], frag[1])
            else:
                nodeid, edonid, nodelab, edgelab, edonlab = frag
                graph.add_node(edonid, edonlab)
                graph.add_edge(nodeid, edonid, edgelab)
        return graph


class RNALabels(object):

    # ====================================================
    # =================== EDGE LABELS ====================
    # ====================================================

    # Labels for edges (and their atoms) in nucleobase edge-wise interactions

    # Labels for general relationships and interactions
    IN_SEQUENCE_WITH = "SEQ" # relationship between entities in sequence
    HYDROGEN_BOND = "HB" # relationship between hydrogen-bonding atoms
    IS_PART_OF = "IPO" # relationship of atom/moiety with larger entity
    CONTIGUOUS = "CON" # relationship between contiguous base pairs

    NUCLEOBASE_EDGES = HOOGSTEEN, SUGAREDGE, WATSCRICK = "H", "S", "W"
    ANY_BASE_EDGE = "E"
    HBOND_VIA_C = "HBC"
    PROTONATED = "+"
    DEPROTONATE = {"+": WATSCRICK, "g": HOOGSTEEN, "z": SUGAREDGE}

    # Labels for magnesium ion interactions
    ION_LIGAND_CONTACT = LIGATION, WATERMED = "L", "LW1"
    BASE_PAIR_CONTACT = "LBP"
    EITHER_SHELL = "LW*"
    CIS, TRANS = "CIS", "TRS"

    EDGELABELS = [HYDROGEN_BOND, IS_PART_OF, IN_SEQUENCE_WITH, CONTIGUOUS,
                    HBOND_VIA_C, PROTONATED, HOOGSTEEN, SUGAREDGE, WATSCRICK,
                    LIGATION, WATERMED, BASE_PAIR_CONTACT, CIS, TRANS]
    EDGELABELS = map(StrLabel, EDGELABELS)
    EDGELABELS += [    StrLabel(ANY_BASE_EDGE, NUCLEOBASE_EDGES),
                    StrLabel(EITHER_SHELL, ION_LIGAND_CONTACT)]
    EDGELABELS.append(StrLabel("*", EDGELABELS))
    EDGELABELS = LabelsList(EDGELABELS)

    # ====================================================
    # =================== NODE LABELS ====================
    # ====================================================

    # Labels for different kinds of nucleotides and residues
    WATER = "HOH"
    WATER_NAMES = {"HOH", "O"}
    UNKNOWN_RESIDUE = "UNK"
    AMINOACID_RESIDUES = sorted(aminoacids)
    NUCLEOTIDE_RESIDUES = sorted(set(allnucs) - {"N", "R", "Y"})
    ANY_NUCLEOTIDE, PURINE, PYRIMIDINE = "5Nu", "R", "Y"
    
    RESIDUES = NUCLEOTIDE_RESIDUES + AMINOACID_RESIDUES
    RESIDUES += hetresdus + someothers

    # Labels for different nucleotide/hetero-residue moieities and atoms
    PHOSPHATE = "P"

    NUCLEOTIDE_ATOMS = sorted(set().union(*nucatoms.values()) - {"P"})
    NUCLEOBASE_ATOMS = set(NUCLEOTIDE_ATOMS)-set(nucatoms['ribose']+phosphate)
    NUCLEOBASE_ATOMS = sorted(NUCLEOBASE_ATOMS)
    NONBRIDGING_PHOSPHATE_O = sorted(set(phosphate) - {"P"})

    ANY_NUCLEOTIDE_ATOM = "nt" # any nucleotide atom
    ANY_NUCLEOSIDE_ATOM = "ns" # any nucleoside atom
    Xb = "Xb" # any nucleobase atom
    Nb, Ob, Or, Oph = "Nb", "Ob", "Or", "Oph" 
    
    Ox = "O*" # any oxygen atom not in NUCLEOTIDE_ATOMS
    Nx = "N*" # any nitrogen atom not in NUCLEOTIDE_ATOMS
    Xx = "X*" # any atom not in NUCLEOTIDE_ATOMS

    ATOMS = NUCLEOTIDE_ATOMS + [Ox, Nx]

    # Labels for base pair families
    BASEPAIR_FAMILIES = combinations_with_replacement(NUCLEOBASE_EDGES, 2)
    BASEPAIR_FAMILIES = product(map("".join, BASEPAIR_FAMILIES), ("C", "T"))
    BASEPAIR_FAMILIES = map("".join, BASEPAIR_FAMILIES)
    ANY_BASE_PAIR = "BP*" # matches any of the BASEPAIR_FAMILIES

    NODELABELS = RESIDUES + [PHOSPHATE] + ATOMS + BASEPAIR_FAMILIES
    NODELABELS = map(StrLabel, NODELABELS)
    NODELABELS += [    
        StrLabel(PURINE, purines),
        StrLabel(PYRIMIDINE, pyrimidines),
        StrLabel(ANY_NUCLEOTIDE, NUCLEOTIDE_RESIDUES),
        # 
        StrLabel(Nb, filter(lambda atom: "N" in atom, NUCLEOBASE_ATOMS)),
        StrLabel(Ob, filter(lambda atom: "O" in atom, NUCLEOBASE_ATOMS)),
        StrLabel(Or, filter(lambda atom: "O" in atom, nucatoms['ribose'])),
        StrLabel(Oph, NONBRIDGING_PHOSPHATE_O),
        # 
        StrLabel(ANY_NUCLEOTIDE_ATOM, NUCLEOTIDE_ATOMS),
        StrLabel(ANY_NUCLEOSIDE_ATOM, NUCLEOBASE_ATOMS + nucatoms['ribose']),
        StrLabel(Xb, NUCLEOBASE_ATOMS),
        StrLabel(Xx, [Ox, Nx]),
        # 
        StrLabel(ANY_BASE_PAIR, BASEPAIR_FAMILIES)
    ]
    NODELABELS.append(StrLabel("*", NODELABELS))

    RESIDUES += [ANY_NUCLEOTIDE, PURINE, PYRIMIDINE]
    ATOMS += [ANY_NUCLEOTIDE_ATOM, ANY_NUCLEOSIDE_ATOM]
    ATOMS += [Nb, Ob, Xb, Or, Oph, Xx]

    # Fuzzy labels for atoms in a particular nucleobase edge
    for _nuc, _edge in product(RNAnucs, NUCLEOBASE_EDGES):
        _name = "%s.%s" % (_nuc, _edge)
        NODELABELS.append(StrLabel(_name, nucedges[_nuc][_edge.lower()]))
        ATOMS.append(_name)
    for _edge in NUCLEOBASE_EDGES:
        _name = "%s.%s" % (PURINE, _edge)
        _atoms = nucedges["A"][_edge.lower()] + nucedges["G"][_edge.lower()]
        NODELABELS.append(StrLabel(_name, _atoms))
        ATOMS.append(_name)
    for _edge in NUCLEOBASE_EDGES:
        _name = "%s.%s" % (PYRIMIDINE, _edge)
        _atoms = nucedges["C"][_edge.lower()] + nucedges["U"][_edge.lower()]
        NODELABELS.append(StrLabel(_name, _atoms))
        ATOMS.append(_name)
    # 
    del _name, _atoms
    NODELABELS = LabelsList(NODELABELS)

    # ====================================================
    
    LABELS = Domain(NODELABELS, EDGELABELS)


class RNAGraph(Graph, RNALabels):

    def __init__(self, name=None):
        super(RNAGraph, self).__init__(name)
        self.pairs = set()

    def add_edge(self, nodeid1, nodeid2, edgelabel):
        """
        edgelabel: A Label or str. The final Label assigned
            to the edge is found in RNAGraph.LABELS.edge,
            using 'edgelabel' to query the LabelsList.
        """
        edgelabel = self.LABELS.edge[edgelabel]
        super(RNAGraph, self).add_edge(nodeid1, nodeid2, edgelabel)

    def add_node(self, nodeid, nodelabel):
        """
        nodelabel: A Label or str. The final Label assigned
            to the edge is found in RNAGraph.LABELS.edge,
            using 'nodelabel' to query the LabelsList.
        """
        nodelabel = self.LABELS.node[nodelabel]
        super(RNAGraph, self).add_node(nodeid, nodelabel)

    def draw_structure(self, output):
        try:
            graph = gt.Graph(directed=False)
            nlabs = graph.new_vertex_property("string")
            elabs = graph.new_edge_property("string")
            graph.vertex_properties["node_labels"] = nlabs
            graph.edge_properties["edge_labels"] = elabs
            
            graph.add_vertex(len(self.nodes))
            for node in self.nodes:
              vertex = graph.vertex(int(node[0]))
              nlabs[vertex] = node[1]
            for edge in self.edges:
              edgetag = (int(edge[0][0]), int(edge[1][0]))
              edgeobj = graph.add_edge(*edgetag)
              elabs[edgeobj] = edge[2]

            pos = draw.radial_tree_layout(graph, graph.vertex(0))
            draw.graph_draw(graph, pos, output=output, output_size=(800, 800),\
                            bg_color=(1,1,1,1), vertex_size=35,\
                            vertex_text=nlabs, edge_text=elabs,\
                            vertex_font_size=18, edge_font_size=18)
        except:
            raise NotImplementedError("graph_tool could not be imported")

    def name_residue_node(self, resnode):
        resid, resname = resnode
        if resname == "N":
            return resid, self.ANY_NUCLEOTIDE
        if resname in self.WATER_NAMES:
            return resid, self.WATER
        if resname in self.RESIDUES:
            return resid, resname
        return resid, self.UNKNOWN_RESIDUE

    @staticmethod
    def resolve_residue_id(nodeid):
        return Residue.from_string(nodeid).entuple()

    def name_residue_part_node(self, resnode, part):
        resid, resname = resnode
        partid = "%s.%s" % (resid, part)
        if resname in self.WATER_NAMES:
            return resid, self.WATER
        if part == self.PHOSPHATE or part in self.NUCLEOTIDE_ATOMS:
            return partid, part
        if "O" in part:
            return partid, self.Ox
        if "N" in part:
            return partid, self.Nx
        return partid, self.Xx
    
    def connect_part_to_residue(self, resnode, part):
        """Adds an atom/resnode as part of a phosphate or nucleotide.
        Args:
            resnode: 2-tuple of (resid, nuctype)
            part: str, name of atom or moiety
        """
        resid, reslab = self.name_residue_node(resnode)
        self.add_node(resid, reslab)
        if reslab in self.WATER_NAMES:
            return # since water has no subparts in this representation.
        # 
        partid, partlab = self.name_residue_part_node(resnode, part)
        self.add_node(partid, partlab)
        self.add_edge(resid, partid, self.IS_PART_OF)
        if part in self.NONBRIDGING_PHOSPHATE_O:
            pid, plab = self.name_residue_part_node(resnode, self.PHOSPHATE)
            self.add_node(pid, plab)
            self.add_edge(resid, pid, self.IS_PART_OF)
            self.add_edge(pid, partid, self.IS_PART_OF)
        return resid, partid

    def name_basepair_node(self, nucnode1, nucnode2, bptype):
        nucid1, nucid2 = nucnode1[0], nucnode2[0]
        edge1, edge2, orient = bptype[0], bptype[1], bptype[2]
        if (edge1, nucid1) > (edge2, nucid2):
            nucid1, nucid2, edge1, edge2 = nucid2, nucid1, edge2, edge1
        nodeid = "%s.%s=%s.%s%s" % (nucid1, edge1, nucid2, edge2, orient)
        nodelab = "".join((edge1, edge2, orient))
        if nodelab not in self.BASEPAIR_FAMILIES:
            raise ValueError("Invalid base pair type: %s" % nodelab)
        return nodeid, nodelab

    def decode_bpfind_edgecodes(self, bptype):
        edges1, edges2 = list(), list()
        edgecode1, edgecode2, _ = bptype
        if edgecode1 in self.DEPROTONATE:
            edges1.append(self.PROTONATED)
            edgecode1 = self.DEPROTONATE[edgecode1]
        if edgecode2 in self.DEPROTONATE:
            edges2.append(self.PROTONATED)
            edgecode2 = self.DEPROTONATE[edgecode2]
        if edgecode1.islower() or edgecode2.islower():
            edges1.append(self.HBOND_VIA_C)
            edges2.append(self.HBOND_VIA_C)
        edges1.insert(0, edgecode1.upper())
        edges2.insert(0, edgecode2.upper())
        return edges1, edges2

    def add_basepair(self, nucnode1, nucnode2, bptype):
        """
        Args:
            nucnode1, nucnode2: 2-tuples of (resid, nuctype)
            bptype: 3-tuple of (edge1, edge2, orientation)
        """
        edges1, edges2 = self.decode_bpfind_edgecodes(bptype)
        bptype = (edges1[0], edges2[0], bptype[2])
        # 
        nucid1, nuclab1 = self.name_residue_node(nucnode1)
        nucid2, nuclab2 = self.name_residue_node(nucnode2)
        bpid, bplab = self.name_basepair_node(nucnode1, nucnode2, bptype)
        self.add_node(nucid1, nuclab1)
        self.add_node(nucid2, nuclab2)
        self.add_node(bpid, bplab)
        # 
        for edgelabel in edges1:
            self.add_edge(nucid1, bpid, edgelabel)
        for edgelabel in edges2:
            self.add_edge(nucid2, bpid, edgelabel)
        if bpid.startswith(nucid1):
            self.pairs.add((nucid1, nucid2, bpid))
        else:
            self.pairs.add((nucid2, nucid1, bpid))
        return nucid1, nucid2, bpid

    def add_hbond(self, resnode1, atom1, resnode2, atom2):
        """
        Args:
            resnode1, resnode2: 2-tuples of (resid, restype)
            atom1, atom2: string names of the atoms from
              residues 1 and 2 respectively that are in
              hydrogen bonding interaction.
        """
        atom1id, atom1lab = self.name_residue_part_node(resnode1, atom1)
        atom2id, atom2lab = self.name_residue_part_node(resnode2, atom2)
        self.connect_part_to_residue(resnode1, atom1)
        self.connect_part_to_residue(resnode2, atom2)
        self.add_edge(atom1id, atom2id, self.HYDROGEN_BOND)
        return atom1id, atom2id

    def add_sequence(self, nucnode_sequence):
        """
        Args:
            sequence: list of nodes in 2-tuples, (reside, nuctype)
        """
        nodeids = list()
        for nucnode in nucnode_sequence:
            nodeid, nodelab = self.name_residue_node(nucnode)
            self.connect_part_to_residue(nucnode, self.PHOSPHATE)
            partid, _ = self.name_residue_part_node(nucnode, self.PHOSPHATE)
            nodeids += [partid, nodeid]
        for nodeid1, nodeid2 in zip(nodeids, nodeids[1:]):
            self.add_edge(nodeid1, nodeid2, self.IN_SEQUENCE_WITH)
        return nodeids

    def add_contiguous_basepairs(self, bpid1, bpid2):
        self.add_edge(bpid1, bpid2, self.CONTIGUOUS)


class RNAStructureGraph(RNAGraph):

    def __init__(self, name=None):
        super(RNAStructureGraph, self).__init__(name)
        self.segts = list()
        self.stems = list()
        self._residue_locus = dict()

    def record_nucleotide_sequence(self, sequence):
        nodeids = self.add_sequence(sequence)
        if len(sequence) > 1:
            nucids = nodeids[1::2]
            # 
            merged = nucids[:]
            added_to = list()
            for seqex, seqce in enumerate(self.segts):
                post_merger = self.merge_sequences(merged, seqce)
                if post_merger:
                    added_to.append(seqex)
                    self.segts[seqex] = merged = post_merger[:]
            # the last added to will have the complete 'merged'
            if added_to:
                for seqex in added_to[::-1]:
                    del self.segts[seqex] # have partial mergers only
                if self.resolve_residue_id(merged[0])\
                  > self.resolve_residue_id(merged[-1]):
                    merged.reverse()
            self.segts.append(merged)
    
    @staticmethod
    def merge_sequences(seq1, seq2):
        common = sorted([(seq1.index(nodeid), seq2.index(nodeid))
                         for nodeid in set(seq1).intersection(seq2)])
        if common:
            len2, len1 = len(seq2), len(seq1)
            if common[0][1] > common[-1][1]: # seq2 indices are decreasing
                seq2.reverse()
                common = [(indx1, len2-indx2-1) for indx1, indx2 in common]
            add_beg, add_end = list(), list()
            beg_com, end_com = common[0], common[-1]
            if beg_com[0] == 0 and beg_com[1] > 0:
                add_beg = seq2[:beg_com[1]]
            if end_com[0] == len1-1 and end_com[1] < len2-1:
                add_end = seq2[end_com[1]+1:]
            return add_beg + seq1 + add_end
        return list()

    def record_contiguous_basepairs(self):
        for segts, bp_across in self._locate_basepairs_by_segtpairs().items():
            stems = [list()]
            segtA, segtB = segts
            for pairA, pairB in combinations(bp_across, 2):
                (idxA1, idxA2, bpidA), (idxB1, idxB2, bpidB) = pairA, pairB
                if idxB1 - idxA1 == 1 and abs(idxB2 - idxA2) == 1:
                    self.add_contiguous_basepairs(bpidA, bpidB)
                    if stems[-1]:
                        (idxY1, idxY2, bpidY), pairZ = stems[-1][-1]
                        if pairZ != pairA or idxB2 == idxY2:
                            stems.append([(pairA, pairB)])
                            continue
                    stems[-1].append((pairA, pairB))
            self.stems += [stem for stem in stems if len(stem) > 1]

    def residue_locus(self, resid):
        try:
            return self._residue_locus[resid]
        except KeyError:
            for segex, segt in enumerate(self.segts):
                try:
                    locus = (segex, segt.index(resid))
                    self._residue_locus[resid] = locus
                    return locus
                except ValueError:
                    continue
            return None

    def _locate_basepairs_by_segtpairs(self):
        bplocs = defaultdict(list)
        for nucid1, nucid2, bpid in self.pairs:
            loc1 = self.residue_locus(nucid1)
            loc2 = self.residue_locus(nucid2)
            if loc1 and loc2:
                if loc1 > loc2:
                    loc1, loc2 = loc2, loc1
                (segex1, index1), (segex2, index2) = loc1, loc2
                bplocs[(segex1, segex2)].append((index1, index2, bpid))
        bplocs = {segts: sorted(bpairs) for segts, bpairs in bplocs.items()}
        return bplocs

    def sort_residue_sequences(self):
        if self.segts:
            sortby = [self.resolve_residue_id(segt[0]) for segt in self.segts]
            self.segts = zip(*sorted(zip(sortby, self.segts)))[1]

    def derive_structure(self):
        self.sort_residue_sequences()
        self.record_contiguous_basepairs()


class RNAStructureEngrapher(object):

    _TYPE = "StructureGraph"
    _NODE_FORMAT = "%% %dd\t%% %ds\t%% %ds\t%%s"
    _EDGE_FORMAT = "%% %dd\t%% %ds"

    HBOND_ELEMENTS = ("O", "N")
    NUCLEOTIDE_RESNAMES = RNALabels.NUCLEOTIDE_RESIDUES + ["N", "R", "Y"]

    def __init__(self, hbfinder=None):
        self._hbonds = hbfinder is not None
        self.hbfinder = hbfinder
        self.structure = None
        self.set_structure(None)

    def set_structure(self, structure=None):
        """
        Args:
            structure: An engrid.Structure instance or None
        """
        self.graph = None
        if self.structure != structure:
            self.structure = structure
            if self.hbfinder:
                self.hbfinder.set_structure(structure)
        self.residues_added = set()

    def engraph_sequences(self, residue_set):
        residue_set = residue_set.union(self.residues_added)
        residues = sorted(residue_set, key=lambda res_: res_.index)
        for begex, endex, _ in self.structure.get_connected_segments():
            while residues and residues[0].index < endex:
                segment = [residues.pop(0)]
                for index, residue in enumerate(residues[:]):
                    if residue.index-1 == segment[-1].index\
                      and residue.index < endex:
                        segment.append(residues.pop(0))
                    else:
                        break
                if len(segment) > 1:
                    sequence = [(str(res_), res_.name) for res_ in segment]
                    self.graph.record_nucleotide_sequence(sequence)
                    self.residues_added.update(segment)

    def engraph_basepairs(self, residue_set, expand):
        """
        Args:
            residue_set: set of residues whose base pairs
              are to be included
            expand: A boolean value indicating whether the
              base pairing information should include
              residues not in the residue_set
        """
        #TODO: needs to be faster
        if self.structure.pairs:
            residue_set = residue_set.union(self.residues_added)
            for nuc1 in residue_set.intersection(self.structure.pairs):
                node1 = (str(nuc1), nuc1.name)
                for nuc2, bptype in self.structure.pairs[nuc1].items():
                    if expand or nuc2 in residue_set:
                        node2 = (str(nuc2), nuc2.name)
                        self.graph.add_basepair(node1, node2, bptype)
                        self.residues_added.update({nuc1, nuc2})

    def engraph_hbondings(self, residue_set, expand):
        """
        Args:
            residue_set: set of residues whose base pairs
              are to be included
            expand: A boolean value indicating whether the
              base pairing information should include
              residues not in the residue_set
        """
        add_hbond = self.graph.add_hbond
        if self._hbonds:
            residue_set = residue_set.union(self.residues_added)
            hbonds = self.hbfinder.find_all_hbonds(residue_set, expand)
            for residue1, atom1, residue2, atom2 in hbonds:
                if atom1.atype in self.HBOND_ELEMENTS\
                  and atom2.atype in self.HBOND_ELEMENTS:
                    resnode1 = (str(residue1), residue1.name)
                    resnode2 = (str(residue2), residue2.name)
                    add_hbond(resnode1, atom1.name, resnode2, atom2.name)
                    self.residues_added.update({residue1, residue2})

    def engraph_structure(self, residue_set=None, expand=False):
        if residue_set is None:
            residue_set = set(structure.residues)
        if not self.graph:
            self.graph = RNAStructureGraph(structure.structid)
        self.engraph_basepairs(residue_set, expand)
        self.engraph_hbondings(residue_set, expand)
        self.engraph_sequences(residue_set)
        if self.residues_added - residue_set:
            self.engraph_basepairs(residue_set, False)
            self.engraph_hbondings(residue_set, False)
        self.graph.derive_structure()

    @classmethod
    def extract_sequences(cls, infile, graphnames=None):
        IN_SEQUENCE_WITH = RNALabels.EDGELABELS[RNALabels.IN_SEQUENCE_WITH]
        PHOSPHATE = RNALabels.NODELABELS[RNALabels.PHOSPHATE]
        # 
        graphs = RNAGraphList.load(infile, graphnames)
        sequences = list()
        for graph in graphs:
            nodes = graph.nodes
            fragments = list()
            for edge in graph.edges_with_label_in([IN_SEQUENCE_WITH]):
                nodeid, edonid, _ = edge
                if PHOSPHATE in (nodes[nodeid], nodes[edonid]):
                    if PHOSPHATE == nodes[nodeid]: # if edonid is nucleotide,
                        nodeid, edonid = edonid, nodeid # swap so nodeid is.
                    edonid = edonid.split(".")[0]
                    if nodeid != edonid:
                        fragments.append([nodeid, edonid])
            nodeid_sequences = cls._merge_fragments(fragments)
            base_sequences = cls._base_sequences(nodeid_sequences)
            sequences.append(list(zip(base_sequences, nodeid_sequences)))
        return sequences

    @classmethod
    def _base_sequences(cls, nodeid_sequences):
        base_sequences = list()
        for nodeid_sequence in nodeid_sequences:
            base_sequence = list()
            for nodeid in nodeid_sequence:
                base = nodeid.split("]")[0].strip("[")
                if base in variants:
                    base_sequence.append(base)
                else:
                    for stdbase in variants:
                        if base in variants[stdbase]:
                            base_sequence.append(stdbase.lower())
                            break
                    else:
                        base_sequence.append("N")
            base_sequences.append("".join(base_sequence))
        return base_sequences

    @classmethod
    def _merge_fragments(cls, fragments):
        segments = list()
        for fragment in fragments:
            assert(fragment[0].split(".")[0] != fragment[1].split(".")[0])
            added_to = list()
            for segex, segment in enumerate(segments):
                if fragment[1] == segment[0]:
                    segment.insert(0, fragment[0])
                    added_to.append(segex)
                elif fragment[0] == segment[-1]:
                    segment.append(fragment[1])
                    added_to.append(segex)
            if not added_to:
                segments.append(fragment)
            elif len(added_to) == 2:
                segB = segments.pop(added_to[1])
                segA = segments.pop(added_to[0])
                if segA[0] == segB[-2]:
                    segments.append(segB[:-2]+segA)
                elif segA[-2] == segB[0]:
                    segments.append(segA[:-2]+segB)
        for segment in segments:
            for idA, idB in zip(segment, segment[1:]):
                if idA == idB:
                    raise ValueError()
        return segments

    def encode_graph(self):
        nodes = {key_: label.name for key_, label in self.graph.nodes.items()}
        edges = defaultdict(set)
        for nodeid, edonid, label in self.graph.edges:
            nodeid, edonid = sorted((nodeid, edonid), key=self._id_sorter)
            edges[edonid].add((nodeid, label.name))
        nformat, eformat = self._formatters(nodes, edges)
        # actual output
        tags = dict()
        outstring = [self._TYPE+":\t"+ self.graph.name.replace(" ", "\t")]
        for index, nodeid in enumerate(sorted(nodes, key=self._id_sorter), 1):
            tags[nodeid] = index
            elist = [(tags[edonid], lab_) for edonid, lab_ in edges[nodeid]]
            elist = "\t".join([eformat % args for args in sorted(elist)])
            outstring.append(nformat % (index, nodeid, nodes[nodeid], elist))
        return "\n".join(outstring)+"\n"

    def _formatters(self, nodes, edges):
        tag_length = len(str(len(nodes)))+1
        id_length = max(map(len, nodes))+1
        nlab_len = max(map(len, nodes.values()))+1
        elab_len = max(map(len, list(zip(*set.union(*edges.values())))[1]))+1
        nformat = self._NODE_FORMAT % (tag_length, -id_length, -nlab_len)
        eformat = self._EDGE_FORMAT % (tag_length, -elab_len)
        return nformat, eformat

    def _id_sorter(self, nodeid):
        restuple = Residue.from_string(nodeid.split("=")[0]).entuple()
        chain, resloc, resname = restuple[1], restuple[3:5], restuple[5]
        isnotnuc = resname not in self.NUCLEOTIDE_RESNAMES
        isbpnode = "=" in nodeid
        return (isbpnode, isnotnuc, chain, resloc, nodeid)


class RNAGraphList(object):

    FILE_FORMAT = """
    ! The header marks the start of a site-graph, and follows the format below:
    ! SiteGraph:	<structure ID>	<magnesium ID><"!" if unreliable site>
    ! 
    ! The graph is represented as a tab-separated adjacency list, padded with
    ! empty spaces. The format is as described below:
    ! 
    !     Column 1 = node serial number
    !     Column 2 = node ID (derived from RasMol primitive expressions)
    !     Column 3 = node label (see legend at the end of file)
    ! 
    !     Columns 4 onwards describes edges incident on the node.
    !     Columns 2n and 2n+1 together describe a particular edge, where
    !         Column 2n gives the serial number of the adjacent node,
    !         Column 2n+1 denotes the edge-label (labels legend at the end)
    !\n!\n\n
    """.replace("\n    ", "\n")

    _TYPE = RNAStructureEngrapher._TYPE

    def __init__(self, filepath, engrapher=None, encoded=False):
        self.graphs = list()
        self.encoded = encoded
        self.filepath = filepath
        self.engrapher = engrapher
        self.node_labels = set()
        self.edge_labels = set()

    @classmethod
    def dump_in(cls, filepath, structure, engrapher, **kwargs):
        graphslist = cls(filepath, engrapher)
        graphslist.add(structure, **kwargs)
        graphslist.dump()

    def add(self, structure, **kwargs):
        self.engrapher.set_structure(structure)
        self.engrapher.engraph_structure(**kwargs)
        self.update()

    def update(self):
        graph = self.engrapher.graph
        name = lambda label: label.name
        self.node_labels.update(map(name, graph.nodes.values()))
        self.edge_labels.update(map(name, zip(*list(graph.edges))[2]))
        if self.encoded:
            self.graphs.append(self.engrapher.encode_graph())
        else:
            self.graphs.append(graph)

    def encode(self):
        assert(self.encoded == False)
        for index in range(len(self.graphs)):
            self.engrapher.graph = self.graphs[index]
            self.graphs[index] = self.engrapher.encode_graph()
        self.encoded = True

    def dump(self):
        if not self.encoded:
            self.encode()
        self.graphs.sort(key=self.sorter)
        legends = LabelsLegend.generate_legends(self.node_labels, self.edge_labels)
        with open(self.filepath, 'w') as outfile:
            outfile.write(self.FILE_FORMAT)
            outfile.write("\n\n".join(self.graphs))
            outfile.write(legends)

    @staticmethod
    def sorter(encoded_graph):
        name = " ".join(encoded_graph.split("\n")[0].split("\t")[1:])
        return name

    @classmethod
    def load(cls, infile, graphnames=None):
        graphs = list()
        untag = None
        with open(infile, 'r') as thefile:
            ingraph = False
            for line in thefile:
                line = line.strip()
                if line.startswith("!"):
                    continue
                if line:
                    line = [field.strip() for field in line.split("\t")]
                    if line[0].startswith(cls._TYPE):
                        gfname = " ".join(line[1:])
                        ingraph = False
                        if not graphnames or gfname in graphnames:
                            graphs.append(RNAGraph(gfname))
                            untag = dict()
                            ingraph = True
                    elif ingraph:
                        graph = graphs[-1]
                        nodetag, nodeid = int(line[0]), line[1]
                        untag[int(nodetag)] = nodeid
                        graph.add_node(nodeid, graph.LABELS.node[line[2]])
                        for edontag, edgelab in zip(line[3::2], line[4::2]):
                            edonid = untag[int(edontag)]
                            edgelab = graph.LABELS.edge[edgelab]
                            graph.add_edge(nodeid, edonid, edgelab)
        if graphnames:
            graphs.sort(key=lambda graph: graphnames.index(graph.name))
        return graphs


class LabelsLegend(object):

    SPECIAL_NODES = { 
        "MG": "Magnesium ion", RNALabels.WATER: "Water oxygen atom",
        RNALabels.PHOSPHATE: "Phosphate moiety in a nucleotide",
        RNALabels.ANY_BASE_PAIR : "Any base pair node (fuzzy)",
        RNALabels.UNKNOWN_RESIDUE: "Default for any unrecognised residue"
    }

    FUZZY_NUCLEOTIDE_ATOMS = {
        RNALabels.ANY_NUCLEOTIDE_ATOM : "Any nucleotide atom",
        RNALabels.ANY_NUCLEOSIDE_ATOM : "Any nucleoside atom",
        RNALabels.Xb : "Any nucleobase atom",
        RNALabels.Nb : "Any nucleobase nitrogen atom",
        RNALabels.Ob : "Any nucleobase oxygen atom",
        RNALabels.Oph : "Any non-bridging phosphate oxygen atom"
    }

    NON_NUCLEOTIDE_ATOMS = {RNALabels.Xx : "Any non-nucleotide atom (fuzzy)",
                            RNALabels.Nx : "Any non-nucleotide nitrogen atom",
                            RNALabels.Ox : "Any non-nucleotide oxygen atom"}

    FUZZY_NUCLEOTIDES = { RNALabels.PURINE : "Purine Nucleotides",
                          RNALabels.PYRIMIDINE : "Pyrimidine Nucleotides",
                          RNALabels.ANY_NUCLEOTIDE : "Any Nucleotide"}

    MODIFIED_NUCLEOTIDES = set(RNALabels.NUCLEOTIDE_RESIDUES)
    MODIFIED_NUCLEOTIDES -= (set(RNAnucs) | set(FUZZY_NUCLEOTIDES))

    OTHER_RESIDUES = set(RNALabels.RESIDUES) - set(FUZZY_NUCLEOTIDES)
    OTHER_RESIDUES -= set(RNALabels.NUCLEOTIDE_RESIDUES)

    NODE_LABELS = [ ("Special Node Labels", SPECIAL_NODES),
                    ("Typical Nucleotide Atoms", RNALabels.NUCLEOTIDE_ATOMS),
                    ("Non-Nucleotide, Non-Water Atoms", NON_NUCLEOTIDE_ATOMS),
                    ("Fuzzy Nucleotide Atom Labels", FUZZY_NUCLEOTIDE_ATOMS),
                    ("Typical RNA Nucleotides", RNAnucs),
                    ("Fuzzy Nucleotide Labels", FUZZY_NUCLEOTIDES),
                    ("Modified Nucleotides", MODIFIED_NUCLEOTIDES),
                    ("Other Residues", OTHER_RESIDUES),
                    ("Base Pair Node Labels", RNALabels.BASEPAIR_FAMILIES)]

    SPECIAL_EDGES = {
      RNALabels.IS_PART_OF : "Smaller node IS PART OF the larger",
      RNALabels.CONTIGUOUS : "Base pairs are CONTIGUOUS",
      RNALabels.HYDROGEN_BOND : "The atoms interact via a HYDROGEN BOND",
      RNALabels.IN_SEQUENCE_WITH : "Phosphate is IN SEQUENCE WITH nucleoside",
    }

    MAGNESIUM_EDGES = {
      RNALabels.LIGATION : "Direct LIGATION (ligand in ion inner-shell)",
      RNALabels.WATERMED : "WATER MEDIATED (ligand in ion outer-shell)",
      RNALabels.EITHER_SHELL : "Ligand in EITHER inner or outer shell (fuzzy)",
      RNALabels.BASE_PAIR_CONTACT : "Ion interacts with a BASE PAIR",
    }

    SHELLGEOMETRY_EDGES = {
      RNALabels.CIS : "Ligands are in CIS orientation w.r.t each other",
      RNALabels.TRANS : "Ligands are in TRANS orientation w.r.t each other"
    }

    BASEPAIRING_EDGES = {
      RNALabels.HOOGSTEEN : "Base interacts via HOOGSTEEEN edge in pair",
      RNALabels.WATSCRICK : "Base interacts via WATSON-CRICK edge in pair",
      RNALabels.SUGAREDGE : "Base interacts via SUGAR edge in pair",
      RNALabels.ANY_BASE_EDGE : "Base interacts via H, W or S edges (fuzzy)",
      RNALabels.PROTONATED : "Base edge involved in base pair is PROTONATED",
      RNALabels.HBOND_VIA_C : "Base pair has HYDROGEN BOND with CARBON DONOR"
    }

    EDGE_LABELS = [ ("Special Edge Labels", SPECIAL_EDGES),
                    ("Edges Involving Magnesium", MAGNESIUM_EDGES),
                    ("Edges for Inner-Shell Geometry", SHELLGEOMETRY_EDGES),
                    ("Edges Describing Base Pairs", BASEPAIRING_EDGES)]

    @classmethod
    def generate_legends(cls, node_labels, edge_labels):
        if not node_labels or not edge_labels:
            return "\n\n"
        # 
        legends = ["", "%sNODE LABELS%s" % (" "*34, " "*34), ""]
        for header, check_labels in cls.NODE_LABELS:
            legends += cls.legends_for(check_labels, header, node_labels)
        # 
        legends += ["", "%sEDGE LABELS%s" % (" "*34, " "*34), ""]
        for header, check_labels in cls.EDGE_LABELS:
            legends += cls.legends_for(check_labels, header, edge_labels)
        # 
        return "\n\n!%s\n" % "\n!".join(legends)


    @classmethod
    def legends_for(cls, check_labels, header, main_labels):
        append = list()
        padded = lambda label: "% 5s" % label
        if main_labels & set(check_labels):
            labels = sorted(map(padded, main_labels & set(check_labels)))
            append.append("    %s" % header)
            if isinstance(check_labels, dict):
                for label in sorted(labels):
                    description = check_labels[label.strip()]
                    append.append("        %s : %s" % (label, description))
            else:
                while labels:
                    append.append("    %s" % ";".join(labels[:12]))
                    labels = labels[12:]
            main_labels -= set(check_labels)
            append.append("")
        return append