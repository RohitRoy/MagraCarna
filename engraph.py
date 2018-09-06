import os
import re
import sys
import time
import traceback
import numpy as np
import itertools as it

try:
    import graph_tool as gt
    import graph_tool.draw as draw
except:
    pass

from collections import defaultdict, deque, namedtuple

from .engrid import Residue
from .aminoacids import aminoacids
from .heteroatoms import hetresdus, someothers
from .nucleotides import allnucs, purines, pyrimidines, nucleotides,\
                         phosphate, atoms as nucatoms, edges as nucedges


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

    def __eq__(self, fles):
        print("label equality")
        traceback.print_stack()
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

    def create_string_label(self, index):
        try:
            label = self.labels[index]
        except IndexError:
            raise IndexError("Label index %s out of range" % str(index))
        if label.is_specific:
            return StrLabel(label.name)
        eqvals = [self.names[index] for index in label.eqvals]
        return StrLabel(label.name, eqvals)

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


LabelsRange = namedtuple('LabelsRange', ['node', 'edge'])


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
    def remult(cls, fragment, labelsrange):
        if len(fragment) == 2:
            return cls((fragment[0], labelsrange.node[fragment[1]]))
        nodelab = labelsrange.node[fragment[2]]
        edgelab = labelsrange.edge[fragment[3]]
        edonlab = labelsrange.node[fragment[4]]
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
        tags = dict()
        for fragment in code:
            tags[fragment.nodetag] = fragment.nodelab
            if not fragment.isnode:
                tags[fragment.edontag] = fragment.edonlab
        self._tags = [label for _, label in sorted(tags.items())]
        self._code = code

    @property
    def unmulted(self):
        return tuple(fragment.unmult() for fragment in self._code)

    def __getitem__(self, index):
        value = self._code[index]
        if isinstance(value, slice):
            return DFSCode(value)
        return value

    def to_string(self, header="", indent="", indenthead=False, outfile=None):
        #TODO: change to "to graph string"
        code_string = (indent if indenthead else "") + header + "\n"
        code_string = "\n"+code_string if header else code_string
        code_string = "".join([indent+line+"\n" for line in self.to_lines()])
        if outfile:
            with open(outfile, 'a') as outfile:
                outfile.write(outstr+"\n")
        else:
            return outstr

    def to_lines(self):
        code_lines = list()
        node_format = ""
        for nodetag in range(len(self._tags)):
            code_lines.append([str(nodetag), self._tags[nodetag].name])
        for frag in self:
            nodetag, edontag, edgelab = frag.nodetag, frag.edontag, frag.edgelab
            if not frag.intree:
                code_lines[nodetag] += [node_format%edontag, frag.edgelab.name]
            elif not frag.isnode:
                code_lines[edontag] += [str(nodetag), frag.edgelab.name]
        return ["\t".join(line) for line in code_lines]

    @classmethod
    def from_lines(cls, lines, domain, base=None):
        code = base.copy() if base else cls([])
        make = lambda frag: DFSCodeFragment.remult(frag, domain)
        for line in lines:
            line = [field.strip() for field in line.strip().split("\t")]
            nodetag, nodelab = int(line[0]), line[1]
            if not line[2:]: # implies isnode
                code.append(make((0, nodelab)))
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
        nodes = sorted(Graph.from_code(self).nodes)
        nlabs += [nodes[ntag][1] for ntag in range(len(nlabs), len(nodes))]
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


class DFSGenerator(object):

    GraphData = namedtuple("GraphData", ["name", "node_labels", "incidence"])

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

    @classmethod
    def graph_count(cls, dfsgens):
        return len({dfsg._graphdata.name for dfsg in dfsgens})

    def children(self, parents=None, backedges=True, treeedges=True):
        if parents == 'rightmost':
            parents = self.queue
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
        while len(self.queue):
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
          or next_frag.nodelab.name != self.labels[next_nodetag].name:
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
        nodelab = self.labels[nodetag]
        if edontag <= nodetag:
            edonlab = self.labels[edontag]
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
        nodestr = lambda node: str(node[1]).rjust(3)+" "+str(node[0]).ljust(3)
        to_str = lambda nodes: [nodestr(node) for node in nodes]
        outstr = to_str(sorted(zip(*list(zip(*self.nodes.items()))[::-1])))
        for edge in sorted(self.edges, key=lambda edge: (edge[2], edge[:2])):
            outstr.append("%s %s %s" % ((edge[2],) + edge[:2]))
        return "\n".join([self.name] + outstr)+"\n"

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
            min_node_label = min(node[1] for node in self.nodes)
            base = DFSCode([DFSCodeFragment((0, min_node_label))])
        
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
        graphdata = DFSGenerator.GraphData(self.name, self.nodes, incidence)
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


class RNAGraph(Graph):

    # Labels for general relationships and interactions
    ISPARTOF = StrLabel("I") # relationship of atom/moiety with larger entity
    SEQUENCE = StrLabel("Q") # relationship between entities in sequence
    PROXIMUM = StrLabel("X") # relationship between nearby atoms
    # STACKING = StrLabel("S")
    # HBONDING = StrLabel("HB")

    # Labels for different kinds of nucleotides and residues
    _NUCNAMES = sorted(set(allnucs) - {"N", "R", "Y"})
    N = StrLabel("N", _NUCNAMES)
    R = StrLabel("R", purines)
    Y = StrLabel("Y", pyrimidines)
    _RESNAMES = _NUCNAMES + aminoacids + hetresdus + someothers
    RESLABELS = [N, R, Y] + [StrLabel(restype) for restype in _RESNAMES]
    _RESMAP = dict(zip(_RESNAMES, RESLABELS[3:]))

    # Labels for different nucleotide/hetero-residue moieities and atoms
    Oph = StrLabel("OP", phosphate[:3])
    _NUCATOMNAMES = sorted(set(sum(nucatoms.values(), [])))
    _RIBOATOMNAMES = nucatoms['ribose'][:]
    _PHOSATOMNAMES = phosphate + [Oph.name]
    _BASEATOMNAMES = set(_NUCATOMNAMES) - set(_RIBOATOMNAMES+_PHOSATOMNAMES)
    _BASEATOMNAMES = sorted(_BASEATOMNAMES)
    n = StrLabel("n", _NUCATOMNAMES)
    r = StrLabel("r", _RIBOATOMNAMES)
    b = StrLabel("b", _BASEATOMNAMES)
    Ob = StrLabel("Ob", [atom for atom in _BASEATOMNAMES if "O" in atom])
    Nb = StrLabel("Nb", [atom for atom in _BASEATOMNAMES if "N" in atom])
    Or = StrLabel("Or", [atom for atom in _RIBOATOMNAMES if "O" in atom])
    hO = StrLabel("hO")
    hN = StrLabel("hN")
    PHOSPHATE = StrLabel("p")
    ATOMLABELS = [StrLabel(atom) for atom in _NUCATOMNAMES]
    ATOMLABELS += [n, b, r, Oph, Ob, Nb, Or, hO, hN]
    _NUCATOMMAP = dict(zip(_NUCATOMNAMES, ATOMLABELS))

    # Labels for edges (and their atoms) in nucleobase edge-wise interactions
    HBONDONC = StrLabel("HBC")
    PROTONATE = StrLabel("+")
    HOOGSTEEN = StrLabel("H")
    SUGAREDGE = StrLabel("S")
    WATSCRICK = StrLabel("W")
    _NUCEDGECODES = {"W": (WATSCRICK,), "H": (HOOGSTEEN,), "S": (SUGAREDGE,),
                     "w": (WATSCRICK, HBONDONC), "+": (WATSCRICK, PROTONATE),
                     "h": (HOOGSTEEN, HBONDONC), "g": (HOOGSTEEN, PROTONATE),
                     "s": (SUGAREDGE, HBONDONC), "z": (SUGAREDGE, PROTONATE)}
    NUCEDGELABELS = [HOOGSTEEN, SUGAREDGE, WATSCRICK, PROTONATE, HBONDONC]
    _NUCEDGENAMES = [_label.name for _label in NUCEDGELABELS]
    ANUCEDGE = StrLabel("e", _NUCEDGENAMES)

    _NUCEDGEATOMSMAP = dict()
    for _lab, _name in zip(NUCEDGELABELS[:3], ["h", "w", "s"]):
        _atoms = list()
        for _nuc in nucleotides:
            _atoms += nucedges[_nuc][_name]
        _NUCEDGEATOMSMAP[_lab] = sorted(set().union(_atoms))

    # Labels for kinds of ligand interactions
    LIGATION = StrLabel("L")
    WATERMED = StrLabel("M")
    LIGLABELS = [LIGATION, WATERMED]
    TOLIGAND = StrLabel("i", LIGLABELS) #TODO: change to 'l'
    
    # Labels for relationships w.r.t. relative geometric orientation
    INCISGEO = StrLabel("C")
    TRANSGEO = StrLabel("T")
    # WEAKLYCIS = 
    # WEAKTRANS = 
    GEOLABELS = [INCISGEO, TRANSGEO] #, WEAKLYCIS, WEAKTRANS]

    # Labels for base pair families and interactions with base pairs
    _BPFAMNAMES = list()
    for _bpfam in it.product(*[NUCEDGELABELS[:3]]*2 + [GEOLABELS]):
        _BPFAMNAMES.append("".join([str(part) for part in _bpfam]))
    BPFAMLABELS = [StrLabel(_bpfam) for _bpfam in _BPFAMNAMES]
    _BPFAMMAP = dict(zip(_BPFAMNAMES, BPFAMLABELS))
    ABASEPAIR = StrLabel("eeo", _BPFAMNAMES)
    TOABPAIR = StrLabel("ip") #TODO: change to something else

    NODELABELS = RESLABELS+[PHOSPHATE]+ATOMLABELS+BPFAMLABELS+[ABASEPAIR]
    EDGELABELS = [ISPARTOF, SEQUENCE, PROXIMUM]+NUCEDGELABELS+[ANUCEDGE]+\
                 LIGLABELS+[TOLIGAND]+GEOLABELS+[TOABPAIR]
    NODELABELS = LabelsList(NODELABELS)
    EDGELABELS = LabelsList(EDGELABELS)
    LABELS = LabelsRange(NODELABELS, EDGELABELS)

    def __init__(self, name=None):
        super(RNAGraph, self).__init__(name)
        self.pairs = set()

    def add_edge(self, nodeid1, nodeid2, edgelabel):
        edgelabel = self.EDGELABELS[edgelabel]
        super(RNAGraph, self).add_edge(nodeid1, nodeid2, edgelabel)

    def add_node(self, nodeid, nodelabel):
        nodelabel = self.NODELABELS[nodelabel]
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
        nodelabel = self._RESMAP[resnode[1]]
        return resnode[0], nodelabel

    @staticmethod
    def resolve_residue_id(nodeid):
        return Residue.from_string(nodeid).entuple()

    def name_residue_part_node(self, resnode, part):
        partname = str(part)
        nodeid = resnode[0]+"."+partname
        if partname == self.PHOSPHATE.name:
            return nodeid, self.PHOSPHATE
        if partname in self._NUCATOMNAMES:
            if partname in self._PHOSATOMNAMES[:3]:
                return nodeid, self.Oph
            return nodeid, self._NUCATOMMAP[partname]
        if "O" in partname:
            return nodeid, self.hO
        if "N" in partname:
            return nodeid, self.hN
        raise ValueError("Cannot decipher atom name %s" % partname)
    
    def connect_part_to_residue(self, resnode, part):
        """Adds an atom/resnode as part of a phosphate or nucleotide.
        Args:
            resnode: 2-tuple of (resid, nuctype)
            part: str, name of atom or moiety
        """
        resid, reslab = self.name_residue_node(resnode)
        partid, partlab = self.name_residue_part_node(resnode, part)
        self.add_node(resid, reslab)
        self.add_node(partid, partlab)
        if part in self._PHOSATOMNAMES:
            phos_part = self.PHOSPHATE.name
            phosid, phoslab = self.name_residue_part_node(resnode, phos_part)
            self.add_node(phosid, phoslab)
            self.add_edge(resid, phosid, self.ISPARTOF)
            self.add_edge(phosid, partid, self.ISPARTOF)
        else:
            self.add_edge(resid, partid, self.ISPARTOF)

    def name_basepair_node(self, nucnode1, nucnode2, bptype):
        nucid1, nucid2 = nucnode1[0], nucnode2[0]
        edge1 = self._NUCEDGECODES[bptype[0]][0]
        edge2 = self._NUCEDGECODES[bptype[1]][0]
        if (edge1.name, nucid1) > (edge2.name, nucid2):
            nucid1, nucid2, edge1, edge2 = nucid2, nucid1, edge2, edge1
        nodeid = "%s.%s=%s.%s" % (nucid1, edge1.name, nucid2, edge2.name)
        bpfamily = "".join((edge1.name, edge2.name, bptype[2]))
        nodelabel = self._BPFAMMAP[bpfamily]
        return nodeid, nodelabel

    def add_basepair(self, nucnode1, nucnode2, bptype):
        """
        Args:
            nucnode1, nucnode2: 2-tuples of (resid, nuctype)
            bptype: 3-tuple of (edge1, edge2, orientation)
        """
        bpid, bplab = self.name_basepair_node(nucnode1, nucnode2, bptype)
        nucid1, nuclab1 = self.name_residue_node(nucnode1)
        nucid2, nuclab2 = self.name_residue_node(nucnode2)
        self.add_node(nucid1, nuclab1)
        self.add_node(nucid2, nuclab2)
        self.add_node(bpid, bplab)
        for edgelabel in self._NUCEDGECODES[bptype[0]]:
            self.add_edge(nucid1, bpid, edgelabel)
        for edgelabel in self._NUCEDGECODES[bptype[1]]:
            self.add_edge(nucid2, bpid, edgelabel)
        if not bpid.startswith(nucid1):
            self.pairs.add((nucid1, nucid2, bpid))
        else:
            self.pairs.add((nucid2, nucid1, bpid))


class RNAStructureGraph(RNAGraph):

    def __init__(self, name=None, plinked=True):
        super(RNAStructureGraph, self).__init__(name)
        self.segts = list()
        self.stems = list()
        self.plinked = plinked
        self._residue_locus = dict()

    def add_sequence(self, sequence):
        """
        Args:
            sequence: list of nodes in 2-tuples, (reside, nuctype)
        """
        nodeids = list()
        phoslab_name = self.PHOSPHATE.name
        for resnode in sequence:
            nodeid, nodelab = self.name_residue_node(resnode)
            if self.plinked: # phosphate and residue nodes are added
                self.connect_part_to_residue(resnode, phoslab_name)
                partid = self.name_residue_part_node(resnode, phoslab_name)[0]
                nodeids.append(partid)
            else: # only residues node is added
                self.add_node(nodeid, nodelab)
            nodeids.append(nodeid)
        for nodeid1, nodeid2 in zip(nodeids, nodeids[1:]):
            self.add_edge(nodeid1, nodeid2, self.SEQUENCE)
        if len(sequence) > 1:
            nucids = nodeids[1::2] if self.plinked else nodeids[:]
            self.record_residue_sequence(nucids)

    def record_residue_sequence(self, newseq):
        merged = newseq[:]
        added_to = list()
        for seqex, seqce in enumerate(self.segts):
            post_merger = self.merge_sequences(merged, seqce)
            if post_merger:
                added_to.append(seqex)
                self.segts[seqex] = merged = post_merger[:]
        # the last added to will have the complete 'merged'
        if added_to:
            for seqex in added_to[:-1][::-1]:
                del self.segts[seqex] # have partial mergers only
            if self.resolve_residue_id(merged[0])\
              > self.resolve_residue_id(merged[-1]):
                merged.reverse()
    
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
        return []

    def record_basepair_sequences(self):
        for bp_across_segts in self._locate_basepairs_by_segtpairs().values():
            stems = [list()]
            for pairA, pairB in zip(bp_across_segts, bp_across_segts[1:]):
                (idxA1, idxA2, bpidA), (idxB1, idxB2, bpidB) = pairA, pairB
                if idxB1 - idxA1 == 1 and abs(idxB2 - idxA2) == 1:
                    self.add_edge((bpidA, bpidB, self.SEQUENCE))
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
                except IndexError:
                    continue
            return None

    def _locate_basepairs_by_segtpairs(self):
        bplocs = dict()
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
        self.record_basepair_sequences()


class RNAStructureEngrapher(object):

    _TYPE = "StructureGraph"
    _NODE_FORMAT = "%% %dd\t%% %ds\t%% %ds\t%%s"
    _EDGE_FORMAT = "%% %dd\t%% %ds"
    NUCLEOTIDE_RESNAMES = allnucs

    def __init__(self, plinked=True):
        self.reset()
        self.plinked = plinked

    def reset(self):
        self.graph = None
        self.residues_added = set()

    def engraph_sequences(self, structure, residue_set):
        residue_set = residue_set.union(self.residues_added)
        residues = sorted(residue_set, key=lambda res_: res_.index)
        for begex, endex, _ in structure.get_connected_segments():
            if not residues:
                break
            segment = list()
            for index, residue in enumerate(residues):
                if residue.index >= endex:
                    break
                if begex <= residue.index:
                    segment.append(residue)
            del residues[:index]
            if len(segment) > 1:
                sequence = [(str(res_), res_.name) for res_ in segment]
                self.graph.add_sequence(sequence)
                self.residues_added.update(segment)

    def engraph_basepairs(self, structure, residue_set, expand):
        #TODO: needs to be faster
        if structure.pairs:
            residue_set = residue_set.union(self.residues_added)
            for nuc1 in residue_set.intersection(structure.pairs):
                node1 = (str(nuc1), nuc1.name)
                for nuc2, bptype in structure.pairs[nuc1].items():
                    if expand or nuc2 in residue_set:
                        node2 = (str(nuc2), nuc2.name)
                        self.graph.add_basepair(node1, node2, bptype)
                        self.residues_added.update({nuc1, nuc2})

    def engraph_structure(self, structure, residue_set=None, expand=False):
        if residue_set is None:
            residue_set = set(structure.residues)
        if not self.graph:
            self.graph = RNAStructureGraph(structure.structid, self.plinked)
        self.engraph_basepairs(structure, residue_set, expand)
        self.engraph_sequences(structure, residue_set)
        self.graph.derive_structure()

    def encode_graph(self, outfile=None):
        outstring = self._encode_graph()
        if outfile:
            outfile.write(outstring)
        else:
            return outstring

    def _encode_graph(self):
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

    @classmethod
    def decode_file(cls, infile):
        graphs = list()
        untag = None
        with open(infile, 'r') as thefile:
            for lineno, line in enumerate(thefile):
                line = line.strip()
                if line:
                    line = [field.strip() for field in line.split("\t")]
                    if line[0].startswith(cls._TYPE):
                        graphs.append(RNAGraph(" ".join(line[1:])))
                        untag = dict()
                    elif lineno:
                        graph = graphs[-1]
                        nodetag, nodeid = int(line[0]), line[1]
                        untag[int(nodetag)] = nodeid
                        graph.add_node(nodeid, graph.LABELS.node[line[2]])
                        for edontag, edgelab in zip(line[3::2], line[4::2]):
                            edonid = untag[int(edontag)]
                            edgelab = graph.LABELS.edge[edgelab]
                            graph.add_edge(nodeid, edonid, edgelab)
        return graphs


