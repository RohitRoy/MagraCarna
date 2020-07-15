import os

from time import time, sleep
from itertools import combinations, product, count
from collections import defaultdict

from .motifs import MotifFilesParser
from .engraph import RNAGraph, DFSCode, DFSCodeFragment, RNAStructureEngrapher


class StopCodeFragment(DFSCodeFragment):

    def covers(self, fles):
        return self.compares_with(fles) and\
               self.covers_labels(fles) and fles.covers_labels(self)


class StopCode(DFSCode):

    def __init__(self, code):
        if isinstance(code, DFSCode):
            code = [StopCodeFragment(frag) for frag in code]
        super(StopCode, self).__init__(code)


class InfrequentError(ValueError):
    """
    """


class BadLabelError(ValueError):
    """
    """


class MotifGraphModifier(object):

    HEADERS = ("S", "E", "L", "T")
    ENDERS = ("L",)

    ORIENTS = [RNAGraph.LABELS.edge[name] for name in ("C", "T")]
    NLABEL = RNAGraph.LABELS.node["N"]
    NTLABEL = RNAGraph.LABELS.node["nt"]
    PHLABEL = RNAGraph.LABELS.node["P"]
    OPHLABEL = RNAGraph.LABELS.node["Oph"]
    LWXLABEL = RNAGraph.LABELS.edge["LW*"]
    LIGATIONS = [RNAGraph.LABELS.edge[name] for name in ("L", "LW1")]

    # class DisconnectedError(ValueError):
    #     """
    #     """

    @staticmethod
    def trim_nodes(code, nodetags):
        rnagf = RNAGraph.from_code(code)
        for nodetag in nodetags:
            rnagf.delete_node(nodetag)
        tagmap = {old_: new_ for new_, old_ in enumerate(sorted(rnagf.nodes))}
        edgeset = set()
        for nodetag, edontag, edgelab in rnagf.edges:
            edgeset.add((tagmap[nodetag], tagmap[edontag], edgelab))
        nodedict = {new_: rnagf.nodes[old_] for old_, new_ in tagmap.items()}
        rnagf.edges = edgeset
        rnagf.nodes = nodedict
        code = DFSCode.from_lines(rnagf.to_string().strip().split("\n"))
        return code

    @classmethod
    def trim_leaves(cls, code):
        rnagf = RNAGraph.from_code(code)
        totrim = list()
        ispartlabel = RNAGraph.LABELS.edge["I"]
        for nodetag in range(len(code._tags)):
            edges = list(rnagf.incident_edges(nodetag))
            if (cls.NTLABEL.covers(rnagf.nodes[nodetag]) and len(edges) == 1\
              or cls.OPHLABEL.covers(rnagf.nodes[nodetag]) and len(edges) == 2)\
              and set(zip(*edges)[2]) == {ispartlabel}:
                totrim.append(nodetag)
        return cls.trim_nodes(code, totrim)

    @classmethod
    def _trim_sequence_ends(cls, code, nodelab, edonlab):
        rnagf = RNAGraph.from_code(code)
        nodes = rnagf.nodes_with_label_in([nodelab])
        nodes = [node for node in nodes if rnagf.nodes[node].covers(nodelab)]
        totrim = set()
        for endnode in nodes:
            connects, hosts = set(), set()
            for nodeid, edonid, label in rnagf.incident_edges(endnode):
                edonid = edonid if nodeid == endnode else nodeid
                if edonlab.covers(rnagf.nodes[edonid]):
                    connects.add(edonid)
                else:
                    hosts.add(edonid)
            if not hosts and len(connects) == 1:
                totrim.add(endnode)
        return cls.trim_nodes(code, totrim)

    @classmethod
    def trim_sequence_ends(cls, code):
    	previous = float("inf")
    	while previous > len(code):
    		previous = len(code)
    		code = cls._trim_sequence_ends(code, cls.PHLABEL, cls.NLABEL)
    		code = cls._trim_sequence_ends(code, cls.NLABEL, cls.PHLABEL)
    	return code

    @classmethod
    def trim_redundant(cls, code):
        code = cls.trim_leaves(code)
        code = cls.trim_sequence_ends(code)
        return cls.trim_islands(code)

    @classmethod
    def trim_islands(cls, code):
        rnagf = RNAGraph.from_code(code)
        totrim = list()
        for nodetag in range(len(code._tags)):
            if not rnagf.incident_edges(nodetag):
                totrim.append(nodetag)
        return cls.trim_nodes(code, totrim)

    @classmethod
    def merge_into_motifsfile(cls, modifieds, motifsfile):
        out_string = list()
        for source in modifieds:
            basename = "_".join(os.path.basename(source).split("_")[:-1])
            for _, header, lines in MotifFilesParser._read_motiffile(source):
                # if header.endswith(cls.HEADERS):
                #     out_string.append("Graph:\t%s_%s" % (basename, header))
                #     out_string.append("\t%s\n" % "\n\t".join(lines))
                if header.endswith(cls.ENDERS):
                    round_no = int(header[:-1])
                    code = cls.trim_leaves(DFSCode.from_lines(lines))
                    name = "Graph:\t%s_%dT" % (basename, round_no)
                    out_string.append(code.to_string(name, "\t").lstrip("\n"))
            out_string.append("")
        with open(motifsfile, 'w') as outfile:
            outfile.write("\n".join(out_string))

    @classmethod
    def retrim_trimmed(cls, rnagf):
        if rnagf.name.endswith("T"):
            codelines = list(map("\t".join, rnagf.to_rows()))
            code = cls.trim_leaves(DFSCode.from_lines(codelines))
            rnagf = RNAGraph.from_code(code, rnagf.name)
        codelines = list(map("\t".join, rnagf.to_rows()))
        code = cls.trim_islands(DFSCode.from_lines(codelines))
        rnagf = RNAGraph.from_code(code, rnagf.name)
        if len(rnagf.mindfscode()._tags) != len(rnagf.nodes):
            raise cls.DisconnectedError
        # code = DFSCode.from_lines(rnagf.to_string().strip().split("\n")[1:])
        # mincode = rnagf.mindfscode().to_string(header=rnagf.name)
        return rnagf.to_string().rstrip("\n")

    @classmethod
    def generalise_ligation(cls, rnagf):
        rnagf.edges -= set(rnagf.edges_with_label_in(cls.ORIENTS))
        ligation_edges = rnagf.edges_with_label_in(cls.LIGATIONS)
        rnagf.edges -= set(ligation_edges)
        for nodeid, edonid, edgelab in ligation_edges:
            rnagf.edges.add((nodeid, edonid, cls.LWXLABEL))
        for nodeid, nodelab in list(rnagf.nodes.items()):
            if cls.OPHLABEL.covers(nodelab):
                rnagf.nodes[nodeid] = cls.OPHLABEL
        return rnagf

    @classmethod
    def ligation_nonspecific(cls, source, destination):
        out_string = list()
        for _, header, codelines in MotifFilesParser._read_motiffile(source):
            header = "Graph:\t%s" % header
            rnagf = RNAGraph.from_code(DFSCode.from_lines(codelines), header)
            rnagf = cls.generalise_ligation(rnagf)
            out_string.append(cls.retrim_trimmed(rnagf))
        with open(destination, 'w') as outfile:
            outfile.write("\n".join(out_string))

    # @classmethod
    # def disconnecting_groups(cls, rnagf, edges):
    #     print(rnagf.to_string())
    #     nodecount = len(rnagf.nodes)
    #     minbuilds = set() # defaultdict(list)
    #     rnagf = rnagf.copy()
    #     rnagf.edges -= set(edges)
    #     for length in range(1, len(edges)+1):
    #         for edgeset in combinations(edges, length):
    #             rnagf_copy = rnagf.copy()
    #             rnagf_copy.edges.update(edgeset)
    #             if len(rnagf_copy.mindfscode()._tags) == nodecount:
    #                 minbuilds.add(edgeset)
    #         if minbuilds:
    #             break
    #     return minbuilds

    # @classmethod
    # def ligation_variable(cls, source, destination):
    #     lwxlabel = RNAGraph.LABELS.edge["LW*"]
    #     lbplabel = RNAGraph.LABELS.edge["LBP"]
    #     out_string = list()
    #     for _, header, codelines in MotifFilesParser._read_motiffile(source):
    #         header = "Graph:\t%s" % header
    #         name_format = "%s_%%d%s" % (header[:-1], header[-1])
    #         rnagf = RNAGraph.from_code(DFSCode.from_lines(codelines), header)
    #         lwx_edges = set(rnagf.edges_with_label_in([lwxlabel]))
    #         lbp_edges = set(rnagf.edges_with_label_in([lbplabel]))
    #         groups = cls.disconnecting_groups(rnagf, lwx_edges|lbp_edges)
    #         groups = {tuple(sorted(set(group)-lbp_edges)) for group in groups}
    #         rnagf.edges -= lwx_edges
    #         copy_no = 1
    #         for edgeset in sorted(groups): # product(*groups):
    #             rnagf_copy = rnagf.copy(name_format % copy_no)
    #             rnagf_copy.edges.update(edgeset)
    #             out_string.append(cls.retrim_trimmed(rnagf_copy))
    #             copy_no += 1
    #     with open(destination, 'w') as outfile:
    #         outfile.write("\n".join(out_string))

    @classmethod
    def merge_and_modify(cls, sources, merged_file, general_file):
        cls.merge_into_motifsfile(sources, merged_file)
        cls.ligation_nonspecific(merged_file, general_file)


class MotifGraphExtender(MotifGraphModifier):

    @staticmethod
    def relabel_frag(newlabels, oldfrag):
        newfrags = list()
        newfragdict = dict(enumerate(oldfrag))
        for index, label in newlabels.items():
            newfragdict[index] = label
        newfrag = tuple(zip(*sorted(newfragdict.items()))[1])
        return DFSCodeFragment(newfrag)

    @classmethod
    def generalise_frag(cls, fraglabels, oldfrag, mode="max"):
        genlabels = dict()
        for (index, labeltype), label in fraglabels.items():
            labelslist = RNAGraph.LABELS[labeltype]
            maxgeneral = list([label])
            for general in labelslist.general_to(label):
                if mode == "all" or not maxgeneral[0].covers(general):
                    if general.covers(maxgeneral[0]):
                        if mode == "max":
                            maxgeneral.pop(0)
                        maxgeneral.append(general)
            genlabels[index] = maxgeneral
        generalised = list()
        indices, labels = zip(*sorted(genlabels.items()))
        for labelcase in product(*labels):
            genlabelcase = dict(zip(indices, labelcase))
            generalised.append(cls.relabel_frag(genlabelcase, oldfrag))
        return generalised[0] if mode == "max" else generalised

    @classmethod
    def generalise_child(cls, frag, mode="max"):
        NODE, EDGE = 0, 1
        fraglabels = dict()
        if frag.isnode:
            fraglabels[(1, NODE)] = frag.nodelab
        else:
            fraglabels[(3, EDGE)] = frag.edgelab
            if frag.intree:
                fraglabels[(4, NODE)] = frag.edonlab
        return cls.generalise_frag(fraglabels, frag, mode)

    @staticmethod
    def relabel_code_by_child(oldcode, child):
        nodelabels = oldcode._tags[:]
        if child.intree:
            nodelabels[child.nodetag] = child.nodelab
        else:
            nodelabels[child.edontag] = child.edonlab
        newcode = oldcode.relabel(nodelabels, [])
        newcode.append(child)
        return newcode

    @staticmethod
    def relabel_code_by_label(oldcode, index, labeltype, label):
        nodelabs = list()
        edgelabs = list()
        NODE, EDGE = 0, 1
        if labeltype == NODE:
            for nodetag, nodelab in enumerate(oldcode._tags):
                if nodetag == index:
                    nodelabs.append(label)
                    break
                nodelabs.append(nodelab)
        if labeltype == EDGE:
            for fragex, frag in enumerate(oldcode):
                if fragex == index:
                    edgelabs.append(label)
                edgelabs.append(frag.edgelab)
        newcode = oldcode.relabel(nodelabs, edgelabs)
        return newcode

    @staticmethod
    def retag_child_by_code(code, child):
        if child.intree:
            child = list(child)
            child[1] = len(code._tags)
            child = DFSCodeFragment(tuple(child))
        code = code.copy_append(child)
        return code


class RNAGraphListAnalyser(MotifGraphExtender):
    
    DELETE_LABELS = ("N_K", "N_M", "N_S", "N_W", "N_B", "N_D", "N_H", "N_V")
    BAD_LABELS = ("ns", "nt") + DELETE_LABELS
    BAD_LABELS = tuple([RNAGraph.LABELS.node[label] for label in BAD_LABELS])

    def __init__(self, graphs):
        self.graphs = graphs
        self.graphsnum = len(graphs)
        self.graphnames = {graph.name for graph in graphs}

    @classmethod
    def specty(cls, extension):
        if isinstance(extension, DFSCode):
            specty = 0.0
            for frag in extension:
                specty += cls.specty(frag)
        elif isinstance(extension, DFSCodeFragment):
            if extension.isnode:
                specty = (1. / len(extension.nodelab.eqvals))
            elif extension.intree:
                specty = (1. / len(extension.edonlab.eqvals))
                specty += (1. / len(extension.edgelab.eqvals))
            else:
                specty = (1. / len(extension.edgelab.eqvals))
        else:
            specty = (1. / len(extension.eqvals))
        return specty

    def embed(self, code, names, covered=0, frequent=True, first=False):
        for graph in self.graphs:
            graph_name = graph.name
            if graph_name not in names:
                continue
            graph_embeds = graph.embed_code(code)
            if frequent and not graph_embeds:
                raise InfrequentError()
            while graph_embeds:
                embed = graph_embeds.pop(0)
                embed.queue = embed.nodeids[covered:]
                yield (graph_name, embed)
                del embed
                if first:
                    break

    def specific_fraglabels(self, code, names):
        fragwise = dict()
        NODE, EDGE = 0, 1
        for name, embed in self.embed(code, names=names):
            edgex = -1
            for frag in code:
                if not frag.isnode:
                    edgex += 1
                if not frag in fragwise:
                    fragwise[frag] = defaultdict(lambda: defaultdict(set))
                labels = fragwise[frag]
                if frag.isnode:
                    labels[(1, NODE)][embed.nodelabel(frag.nodetag)].add(name)
                else:
                    labels[(3, EDGE)][embed.edges[edgex][2]].add(name)
                    labels[(2, NODE)][embed.nodelabel(frag.nodetag)].add(name)
                    labels[(4, NODE)][embed.nodelabel(frag.edontag)].add(name)
        return fragwise

    def specify_label_max_coverage(self, selected, labeltype):
        eqvals = set(sum([list(label.eqvals) for label in selected], []))
        by_eqval_len = defaultdict(list)
        for label in RNAGraph.LABELS[labeltype]:
            if eqvals.issubset(label.eqvals)\
              and label.name not in self.DELETE_LABELS:
                by_eqval_len[len(label.eqvals)].append(label)
        return by_eqval_len[min(by_eqval_len)][0]

    def specify_child(self, fraglabels, oldfrag, names):
        bestlabels = dict()
        for (index, labeltype), labels in sorted(fraglabels.items()):
            coverage = set()
            selected = list()
            for label in sorted(labels, key=lambda alab: -len(labels[alab])):
                coverage |= (labels[label] & names)
                selected.append(label)
                if len(coverage) >= len(names):
                    break
            bestlabel = self.specify_label_max_coverage(selected, labeltype)
            if bestlabel != oldfrag[index]:
                bestlabels[index] = bestlabel # forced hierarchy
        return self.relabel_frag(bestlabels, oldfrag)

    def specify(self, oldcode, names):
        names = self.graphnames if not names else names
        newcode = None
        fragwise = self.specific_fraglabels(oldcode, names=names)
        for frag in oldcode:
            newfrag = self.specify_child(fragwise.pop(frag), frag, names)
            if not newcode:
                newcode = DFSCode([newfrag])
            else:
                newcode = self.relabel_code_by_child(newcode, newfrag)
        embeds = self.embed(newcode, names=names, first=True, frequent=True)
        try:
            for _ in embeds:
                pass
            return newcode
        except InfrequentError:
            return oldcode


class ConsensusFinder(RNAGraphListAnalyser):

    EDGE_LIMIT = -1
    NODE_LIMIT = -1
    ITERATIONS = {'tree': 0, 'loop': 0}
    HEADER = "Graph:\t%d%s"

    def __init__(self, graphs, base):
        self.base = DFSCode(base)
        self.began = None
        super(ConsensusFinder, self).__init__(graphs)

    def embed(self, code, covered=0, **kwargs):
        superembed = super(ConsensusFinder, self).embed
        kwargs.update({"names": self.graphnames, "frequent": True})
        return superembed(code, covered=covered, **kwargs)

    def specify(self, code):
        superspecify = super(ConsensusFinder, self).specify
        return superspecify(code, self.graphnames)

    def common_extensions_old(self, code, covered, mode):
        if mode == 'tree':
            # tocover = None
            arguments = ('rightmost', False, True)
        else:
            covered = 0 ; # tocover = None
            arguments = (None, True, False)
        common = None
        print("\n\n\t%d : \t%.3f" % (covered, time() - self.began))
        common = set()
        for _, embed in self.embed(code, covered): #, tocover):
            children = set()
            for frag in embed.children(*arguments).keys():
                 children.add(self.generalise_child(frag))
            common = children if common is None else common & children
        print("\t\t%.3f\n" % (time() - self.began))
        return common

    def common_extensions(self, code, covered, mode):
        if mode == 'tree':
            arguments = ('rightmost', False, True)
        else:
            arguments, covered = (None, True, False), 0
        print("\n\n\t%d : \t%.3f" % (covered, time() - self.began))
        children = defaultdict(set)
        for _, embed in self.embed(code, covered):
            children[embed.id].update(embed.children(*arguments).keys())
        common = set.intersection(*children.values())
        for name, frags in children.items():
            children[name] = {self.generalise_child(frag) for frag in frags}
        common.update(set.intersection(*children.values()))
        common = sorted(common, key=lambda frag: (-self.specty(frag), frag))
        print("\t%d\t%.3f\n" % (len(common), time() - self.began))
        return common

    def extend_old(self, code, covered, mode):
        common = self.common_extensions(code, covered, mode)
        for iteration in count():
            print(len(common))
            if not common or iteration == self.ITERATIONS[mode]:
                break
            legacy = list()
            for index, frag in enumerate(common, 1):
                try:
                    newcode = self.retag_child_by_code(code, frag)
                    newcode = self.specify(newcode)
                    newfrag = newcode[-1]
                    if newfrag.intree and newfrag.edonlab in self.BAD_LABELS:
                        raise BadLabelError
                    code = newcode
                    legacy.append(frag)
                    printargs = (index, time() - self.began, newfrag)
                    print("\t%d\t%.3f\t%s" % printargs)
                except (InfrequentError, BadLabelError):
                    print("\t%d\t%.3f" % (index, time() - self.began))
                    pass
            common = legacy
        return code

    def extend(self, code, covered, mode):
        common = self.common_extensions(code, covered, mode)
        for index, frag in enumerate(common, 1):
            for iteration in count(1):
                try:
                    newcode = self.retag_child_by_code(code, frag)
                    newcode = self.specify(newcode)
                    newfrag = newcode[-1]
                    if newfrag.intree and newfrag.edonlab in self.BAD_LABELS:
                        raise BadLabelError
                    code = newcode
                    args = (index, time() - self.began, iteration, newfrag)
                    print("\t%d\t%.3f\t%d\t%s" % args)
                except (InfrequentError, BadLabelError):
                    print("\t%d\t%.3f\n" % (index, time() - self.began))
                    break
                if iteration == self.ITERATIONS[mode]:
                    print("")
                    break
        return code

    def start(self, outfile=None, rounds=1, covered=0):
        self.began = time()
        code = self.base
        # 
        def output(code, iteration, mode):
            header = self.HEADER % (iteration, mode)
            if outfile:
                code.to_string(header, "\t", outfile=outfile)
        # 
        output(code, 0, "I")
        if not covered:
            code = self.specify(code) ; output(code, 0, "S")
            code = self.extend(code, 0, 'loop') ; output(code, 0, "L")
        if not rounds:
            output(code, 0, "O")
            return
        for iteration in count(1):
            print("\n\nI:%d" % iteration)
            try:
                newcode = self.extend(code, covered, 'tree')
                output(newcode, iteration, "E")
                newcode = self.extend(newcode, covered, 'loop')
                output(newcode, iteration, "L")
            except InfrequentError:
                print("infreq")
                break
            if newcode == code:
                print(len(newcode), len(code))
                break
            covered = len(code._tags)
            code = newcode
            if iteration == rounds:
                break
        output(code, iteration, "O")
        return code


class Discriminator(RNAGraphListAnalyser):

    MAX_SPLIT = 10

    def __init__(self, graphs, base):
        super(Discriminator, self).__init__(graphs)
        self.base = DFSCode(base)
        self.began = None
        self.queue = list()
        self.minimum = self.specty(base)

    def fracty(self, names, nameset, mode):
        fracty = float(len(names)) / len(nameset)
        if mode == "cut":
            return min(fracty, 1-fracty)
        return fracty

    def best_specified(self, code, nameset, mode):
        best = list()
        fragwise = self.specific_fraglabels(code, nameset)
        repeating = set()
        NODE, EDGE = 0, 1
        for fragex, frag in enumerate(code):
            for labeltype, labels in fragwise[frag].items():
                if labeltype == (2, NODE) and not frag.isnode\
                  or labeltype == (4, NODE) and not frag.intree:
                    continue
                labeltag = fragex if labeltype[1] else frag[1]\
                                  if labeltype[0] == 4 else frag[0]
                for label, names in labels.items():
                    entry = (labeltag, labeltype[1], label)
                    if entry in repeating or label == frag[labeltype[0]]:
                        continue
                    specty = self.specty(label)
                    fracty = self.fracty(names, nameset, mode) * specty
                    best.append((fracty, entry, names))
                    repeating.add(entry)
                    break
        best = sorted(best)[::-1]
        print("\tSpecified")
        return best

    def select_expanse(self, children):
        best = list()
        newchildren = defaultdict(set, children)
        for child, names in children.items():
            if (child.edgelab.name == "I" and self.NTLABEL.covers(child.edonlab))\
              or self.graphsnum == len(names):
                newchildren.pop(child)
                continue
            for general_child in self.generalise_child(child, "all"):
                if not child.covers(general_child)\
                  and general_child.edonlab not in self.BAD_LABELS:
                    newchildren[general_child].update(names)
        return newchildren

    def best_expanded(self, code, nameset, mode, covered=0):
        best = list()
        children = defaultdict(set)
        arguments = ('rightmost', False, True)
        for name, embed in self.embed(code, nameset, covered, False):
            for child in embed.children(*arguments).keys():
                children[child].add(name)
        children = self.select_expanse(children)
        for child, names in children.items():
            fracty = self.fracty(names, nameset, mode) * self.specty(child)
            best.append((fracty, child, names))
        best = sorted(best)[::-1]
        print("\tExpanded")
        return best

    def best_looped(self, code, nameset, mode):
        best = list()
        children = defaultdict(set)
        arguments = (None, True, False)
        for name, embed in self.embed(code, nameset, 0, False):
            for child in embed.children(*arguments).keys():
                children[self.generalise_child(child)].add(name)
        for child, names in children.items():
            fracty = self.fracty(names, nameset, mode) * self.specty(child)
            best.append((fracty, child, names))
        best = sorted(best)[::-1]
        print("\tLooped")
        return best

    def best(self, code, nameset, mode, covered=0):
        best = self.best_specified(code, nameset, mode)
        best += self.best_expanded(code, nameset, mode, covered)
        best += self.best_looped(code, nameset, mode)
        cutsorter = lambda entry: (entry[0], not isinstance(entry, tuple))
        if mode == "shed":
            sorter = lambda entry: (len(entry[2]), cutsorter(entry))
            return sorted(best, key=sorter)[::-1]
        else:
            return sorted(best, key=cutsorter)[::-1]

    def extend_code(self, code, extension):
        if isinstance(extension, DFSCodeFragment):
            code = self.relabel_code_by_child(code, extension)
            code = self.retag_child_by_code(DFSCode(code[:-1]), extension)
        else:
            code = self.relabel_code_by_label(code, *extension)
        return code

    def extend(self, code, extension, nameset, cutoff):
        code = self.extend_code(code, extension)
        total = len(nameset)
        common = set()
        covered = 0
        while cutoff > 0.5:
            best = self.best(code, nameset, "shed", covered)
            for fracty, case, names in best:
                print("\t\t%.3f\t%s\t%d" % (fracty, case, len(names)))
            print("\n")
            for fracty, case, names in best:
                names &= nameset
                if len(names) < (cutoff * total):
                    continue
                newcode = self.extend_code(code, case)
                embeds = self.embed(newcode, names, frequent=False)
                newnameset = {name for name, _ in embeds}
                if len(newnameset) >= (cutoff * total):
                    code, nameset = newcode, newnameset
                    printargs = (fracty, self.specty(code), len(nameset), case)
                    print("\t\t%.3f\t%.3f\t%d\t%s" % printargs)
                    if self.specty(code) >= 2*self.minimum:
                        return (code, nameset)
                else:
                    print("\t\t%.3f\t%s" % (fracty, case))
            cutoff *= cutoff
            covered = len(code._tags)
        return None

    def extend_base(self):
        code = self.base
        covered = 0
        for iteration in range(3):
            grown = self.best(code, self.graphnames, "shed", covered)
            for fracty, extension, names in grown:
                if len(names) == self.graphsnum:
                    try:
                        newcode = self.extend_code(code, extension)
                        for _, _ in self.embed(newcode, names, frequent=True):
                            pass
                        code = newcode
                    except InfrequentError:
                        continue
            covered = len(code._tags)
        print(code.to_string("\n\nSetup\t%d" % len(self.graphnames), "\t"))
        self.base = code
        self.minimum = self.specty(code)
        return code

    def start(self):
        nameset = self.graphnames
        code = self.extend_base()
        while nameset:
            cuts = self.best(code, nameset, "cut")[:self.MAX_SPLIT]
            for count, (fracty, extension, names) in enumerate(cuts):
                print(code.to_string("\nExtend %d" % count, "\t"))
                print("\t%.3f\t%s\t%d\n" % (fracty, extension, len(names)))
                extended = self.extend(code, extension, names, cutoff=0.99)
                if extended:
                    newcode, newnameset = extended
                    print(newcode.to_string("First", "\t"))
                    nameset -= newnameset
                    break
            else:
                break
