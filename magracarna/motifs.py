from __future__ import print_function

import os
import sys
import math

from io import StringIO
from itertools import combinations, chain, product, count
from collections import defaultdict

from .engraph import RNALabels, RNAStructureGraph, DFSCode, DFSCodeFragment


class MotifGraphCode(DFSCode):

    def __init__(self, code):
        super(MotifGraphCode, self).__init__(code)
        self.symmetry = self._count_symmetry()
        self._complexity = self._count_complexity()

    def append(self, newfrag):
        super(MotifGraphCode, self).append(newfrag)
        self.symmetry = self._count_symmetry()

    def relabel(self, nlabs, elabs):
        relabelled = super(MotifGraphCode, self).relabel(nlabs, elabs)
        return MotifGraphCode(relabelled._code)

    def copy(self):
        return MotifGraphCode(super(MotifGraphCode, self).copy()._code)

    @classmethod
    def from_lines(cls, lines, domain=None, base=None):
        dfscode = super(MotifGraphCode, cls).from_lines(lines, domain, base)
        return MotifGraphCode(dfscode._code)

    def __iadd__(self, newcode):
        super(MotifGraphCode, self).__iadd__(newcode)
        self.symmetry = self._count_symmetry()

    def _count_symmetry(self):
        if len(self):
            graph = RNAStructureGraph.from_code(self)
            return len(graph.embed_code(self))
        return None

    def _count_complexity(self):
        specificity = 0
        for nodetag, nodelab in enumerate(self._tags):
            specificity += (1.0 / len(nodelab.eqvals))
        for frag in self:
            if not frag.isnode:
                specificity += (1.0 / len(frag.edgelab.eqvals))
        return specificity

    def _iter_tagstr(self):
        """
        Returns:
            2-tuple (isnode, tagstr)
              isnode : Boolean for whether the s
        """
        for frag in self:
            if frag.isnode:
                nodeargs = (frag.nodetag, frag.nodelab)
                yield (True, "%d %s" % nodeargs, nodeargs)
            else:
                edgeargs = (frag.edgelab, frag.nodetag, frag.edontag)
                edgetuple = (frag.nodetag, frag.edontag, frag.edgelab)
                yield (False, "%s (%d->%d)" % edgeargs, edgetuple)
                if frag.intree:
                    nodeargs = (frag.edontag, frag.edonlab)
                    yield (True, "%d %s" % nodeargs, nodeargs)

    def __abs__(self):
        return self._complexity


class MotifAssignments(object):
    
    def __init__(self, assigner, limited=set()):
        self.assigner = assigner
        self.by_case = defaultdict(set)
        self.by_graph = defaultdict(dict)
        self.by_formula = defaultdict(dict)
        self.by_type = {"Graph": self.by_graph, "Formula": self.by_formula}
        self._motifs = dict()
        self._loaded = set()
        self._header = dict()
        self.limit_to(limited)

    def limit_to(self, limited):
        if limited:
            limited = set(limited)
            self._motifs["Graph"] = set(self.assigner.graphs) & limited
            self._motifs["Formula"] = set(self.assigner.formulae) & limited
            for mtype in self.by_type:
                for motif in self.by_type[mtype]:
                    if motif not in self._motifs[mtype]:
                        del self.by_type[mtype][motif]
            for case in self.by_case:
                self.by_case[case] &= limited
        else:
            self._motifs["Graph"] = set(self.assigner.graphs)
            self._motifs["Formula"] = set(self.assigner.formulae)

    def motif_iterator(self):
        limited = self._motifs["Graph"]|self._motifs["Formula"]
        return self.assigner.motif_iterator(limited)

    def add_from(self, rnagfs):
        assigner = self.assigner
        for rnagf in rnagfs:
            preassigns, prerejects = set(), set()
            if rnagf.name in self.by_case:
                preassigns = set(self.by_case[rnagf.name])
                prerejects = self._loaded - preassigns
            assigns = assigner.assign_motifs(rnagf, preassigns, prerejects)
            self.by_case[rnagf.name].update(assigns)
            for motif in assigns: # will not contain preassigns
                if motif in self._motifs["Graph"]:
                    self.by_graph[motif][rnagf.name] = assigns[motif]
                elif motif in self._motifs["Formula"]:
                    self.by_formula[motif][rnagf.name] = assigns[motif]

    @classmethod
    def load(cls, assigner, filepath, detail=False, limited=set()):
        assignments = cls(assigner, limited)
        for isheader, mtype, motif, line in cls._assigns_reader(filepath):
            if isheader:
                if motif in assignments._motifs[mtype]:
                    assignments._loaded.add(motif)
                    assignments._header[motif] = line.split("\t")[3:]
            elif motif in assignments._motifs[mtype]:
                line = [value.strip() for value in line.split("\t")]
                case, record = " ".join(line[:2]), line[2:]
                assignments.by_case[case].add(motif)
                if case not in assignments.by_type[mtype][motif]:
                    details = assignments.by_type[mtype][motif][case] = list()
                    if not detail:
                        details.append(list())
                if detail:
                    if mtype == "Formula":
                        details += assignments._load_truth_details(record)
                    else:
                        details.append(record[1:])
        return assignments

    def dump(self, filepath, detail=False, append=False):
        outstr = list()
        for mtype, motif in self.motif_iterator():
            if motif in self.by_type[mtype] and self.by_type[mtype][motif]:
                outstr += self._dump_graph(motif, detail) if mtype == "Graph"\
                          else self._dump_formula(motif, detail)
                outstr.append("\n")
        with open(filepath, 'a' if append else 'w') as outfile:
            outfile.write("\n".join(outstr))

    @classmethod
    def merge_files(cls, assigner, assignsfiles, outputfile):
        lines = defaultdict(list)
        for assignsfile in assignsfiles:
            for isheader, _, motif, line in cls._assigns_reader(assignsfile):
                if (isheader and motif not in lines) or not isheader:
                    lines[motif].append(line)
        out_string = list()
        for _, motif in assigner.motif_iterator(limited=set(lines)):
            out_string += [""]+lines[motif]
        with open(outputfile, 'w') as assignsfile:
            assignsfile.write("\n".join(out_string))

    @classmethod
    def quick_formula(cls, motif, assigner, assignsfiles, outputfile=None):
        motifassigns = list()
        limits = set(list(zip(*sum(assigner.formulae[motif], [])))[1])
        for assignsfile in assignsfiles:
            assigns = cls.load(assigner, assignsfile, False, limits)
            for case in sorted(assigns.by_case):
                assignset = set(assigns.by_case[case])
                if assigner._truth_values(assignset, motif)[0]:
                    motifassigns.append(case)
        if outputfile:
            cls._write_quick_formula(motif, motifassigns, outputfile)
        else:
            return motifassigns
    
    @classmethod
    def _write_quick_formula(cls, motif, cases, outputfile):
        out_string = ["Formula\t%s" % motif]
        for case in cases:
            out_string.append(case.replace(" ", "\t"))
        with open(outputfile, 'a') as quickassigns:
            quickassigns.write("\n".join(out_string)+"\n\n")

    @classmethod
    def limit_file(cls, assignsfile, limited=set(), outputfile=None):
        out_string = list()
        keep_motif = False
        if_headers = list()
        cases = limited.get("cases", set())
        motifs = limited.get("motifs", set())
        for isheader, _, motif, line in cls._assigns_reader(assignsfile):
            if isheader:
                keep_motif = motif in motifs if motifs else True
                if if_headers and if_headers[-1]:
                    out_string.pop()
                    if_headers.pop()
                    keep_motif = False
                elif keep_motif:
                    if_headers.append(True)
                    out_string.append("\n"+line)
            elif keep_motif: # checking case
                if not cases or line.split("\t")[1] in cases:
                    if_headers.append(False)
                    out_string.append(line)
        if outputfile:
            with open(outputfile, 'w') as outfile:
                outfile.write("\n".join(out_string))
        else:
            return "\n".join(out_string)

    @staticmethod
    def _assigns_reader(assignsfile):
        mtype, motif = None, None
        is_stream = isinstance(assignsfile, StringIO)
        with (assignsfile if is_stream else open(assignsfile, 'r')) as stream:
            for line in stream:
                line = line.strip()
                if line.startswith(("Formula", "Graph")):
                    mtype, motif = line.split("\t")[:2]
                    yield (True, mtype, motif, line)
                elif line:
                    yield (False, mtype, motif, line)

    @staticmethod
    def _load_truth_details(record):
        truth_details = list()
        for index, value in enumerate(record, -len(record)):
            if not value and index != -1:
                clause = list()
                truth_details.append(clause)
            elif value == "T":
                clause.append(True)
            elif value == "F":
                clause.append(False)
            else:
                clause = list() # one that will be ignored
        while not truth_details[-1]:
            truth_details.pop()
        return truth_details

    def _dump_graph(self, graph, detail=False):
        cases = sorted(self.by_graph[graph])
        count = self.assigner.graphs[graph].symmetry
        records = [["Graph", graph, str(count)]]
        records += [case.split(" ") for case in cases]
        if detail:
            detailed = [records[0] + self._graph_header_details(graph)]
            for case, record in zip(cases, records[1:]):
                for embed in self.by_graph[graph][case]:
                    detailed.append(record[:] + [""] + embed)
            records = detailed
        records = ["\t".join(record) for record in records]
        return records

    def _graph_header_details(self, graph):
        if graph not in self._header:
            header = list()
            for _, tagstr, _ in self.assigner.graphs[graph]._iter_tagstr():
                header.append(tagstr)
            self._header[graph] = header
        return self._header[graph]

    def _dump_formula(self, formula, detail=False):
        _HAS = self.assigner._HAS
        cases = sorted(self.by_formula[formula])
        records = [["Formula", formula]] + [case.split(" ") for case in cases]
        if detail:
            for clause in self.assigner.formulae[formula]:
                records[0].append("")
                records[0] += [_HAS[has_]+term for has_, term in clause]
            for case, record in zip(cases, records[1:]):
                record += self._formula_case_details(formula, case)
        records = ["\t".join(record) for record in records]
        return records

    def _formula_case_details(self, formula, case):
        details = list()
        truths = defaultdict(list, enumerate(self.by_formula[formula][case]))
        for clausex, clause in enumerate(self.assigner.formulae[formula]):
            details.append("")
            details += ["T" if truth else "F" for truth in truths[clausex]]
            details += ["?" for _ in clause[len(truths[clausex]):]]
        return details


class MotifAssigner(object):

    _ROOT = ""
    MODES = ["all", "none", "first"]
    _HAS = {True: "", False: "no "}

    class StopBranchIteration(Exception):
        """
        """

    def __init__(self, definedin, graphs, formulae, mode):
        if mode not in self.MODES:
            raise ValueError("Invalid 'mode' argument '%s'" % mode)
        self.mode = mode
        self.graphs = graphs
        self.formulae = formulae
        self.definedin = definedin # motif to motif-file
        self.graph_tree = self._resolve_graph_tree()
        self.formula_list = self._resolve_formula_order()

    @classmethod
    def extract(cls, *motiffiles, **kwargs):
        mode = kwargs.pop("mode", "all")
        parser = MotifFilesParser()
        parser.extract(*motiffiles)
        return cls(parser.definedin, parser.graphs, parser.formulae, mode)

    def print_motifs(self):
        for graph, code in self._by_complexity(self.graphs):
            print(code.to_string(header="Graph:\t%s" % graph, indent="\t"))
        for formula in self.formula_list:
            print("Formula:\t%s" % formula)
            for clause in self.formulae[formula]:
                line = [self._HAS[has_]+part for has_, part in clause]
                print("\t"+(" and ".join(line)))
            print("")

    def print_graph_tree(self):
        print("\n")
        for depth, pattern in self.graph_iterator():
            print("%s%s" % ("\t"*depth, pattern))
        print("\n")

    def motif_iterator(self, limited=set()):
        graph_iterator = self.graph_iterator()
        for _, graph in graph_iterator:
            if not limited or graph in limited:
                try:
                    yield ("Graph", graph)
                except self.StopBranchIteration:
                   graph_iterator.throw(self.StopBranchIteration)
        for formula in self.formula_list:
            if not limited or formula in limited:
                yield ("Formula", formula)

    def assign_graphs(self, rnagf, preassigns=set(), prerejects=set()):
        """

        Arguments 'preassigns' and 'prerejects' are taken
        as fact, and are not checked again.

        Embeddings are represented as a list of nodeids
        and in-graph edge labels, following the DFSCode of
        the graph. Embeddings begin with a nodeid of the
        first node in the code. Thereafter, for each tree
        edge, an edge label and a nodeid are added. For
        each back edge, only the apt edge label is added.

        Args:
            rnagf: RNAStructureGraph object.
            prerejects: 'set' of motifs that are known to
              be present in 'rnagf'. See above.
            prerejects: 'set' of motifs that are known to
              be absent in 'rnagf'. See above.

        Returns:
            A 'dict' mapping motif-graphs assigned to 
            'rnagf' (excluding 'preassigns'), to values
            depending on the 'mode' attribute as follows:
                "all": 'list' of all corresponding
                  embeddings, as described above
                "first": 'list' with only the first
                  embedding, if any, discovered.
                "none": boolean 'True'.
        """
        preassigns = preassigns.intersection(self.graphs)
        prerejects = prerejects.intersection(self.graphs)
        assigns = dict()
        iterator = self.graph_iterator()
        for _, graph in iterator:
            if graph in prerejects:
                iterator.throw(self.StopBranchIteration)
            elif graph not in preassigns:
                embeddings = self._embeddings(graph, rnagf)
                if embeddings:
                    assigns[graph] = embeddings
        return assigns

    def assign_motifs(self, rnagf, preassigns=set(), prerejects=set()):
        """
        Arguments 'preassigns' and 'prerejects' are taken
        as fact, and are not checked again. If their union
        doesn't account for all motif-graphs in 'graphs'
        attribute, then 'assign_graphs' method is run with
        these arguments, before assigning motif-formulae.

        Args:
            rnagf: RNAStructureGraph object.
            prerejects: 'set' of motifs that are known to
              be present in 'rnagf'. See above.
            prerejects: 'set' of motifs that are known to
              be absent in 'rnagf'. See above.

        Returns:
            A 'dict' mapping motifs assigned to 'rnagf'
            (excluding 'preassigns') to motif-graph
            embeddings, or term-wise truth-values in
            motif-formulae. 
            See 'assign_graphs' for more about embeddings.
            Values mapped to motif-formule depend on
            'mode' attribute as follows:
                "all": a 'list' of clause-wise 'list's of
                  boolean values, indicating truth value
                  for every term in the formula.
                "first": a 'list' of 'list' objects, as in
                  the previous 'mode', but only containing
                  truth values for those terms which were
                  used in deciding the assignmnent of the
                  motif-formula.
                "none": boolean 'True'.
        """
        assigns = dict()
        assignset = preassigns.copy()
        considered = set(self.graphs).union(self.formulae)
        considered &= (preassigns|prerejects)
        if considered.intersection(self.graphs) != set(self.graphs):
            assigns.update(self.assign_graphs(rnagf, preassigns, prerejects))
            assignset.update(assigns)
        for formula in self.formula_list:
            if formula not in considered:
                assigned, truths = self._truth_values(assignset, formula)
                if assigned:
                    assigns[formula] = truths
                    assignset.add(formula)
        return assigns

    def graph_iterator(self):
        queue = [(0, self._ROOT)]
        considered = set()
        while queue:
            try:
                entry = depth, node = queue.pop()
                if node in considered:
                    continue
                kids = set(self.graph_tree[node]) - considered
                if node != self._ROOT:
                    considered.add(node)
                    yield entry
            except self.StopBranchIteration:
                considered.update(kids)
                yield 
            else:
                kids = {kid_: self.graphs[kid_] for kid_ in kids}
                kids = self._by_complexity(kids)[::-1]
                queue += [(depth+1, kid_) for kid_, _ in kids]

    def type_of(self, motif):
        if motif in self.graphs:
            return "Graph"
        elif motif in self.formulae:
            return "Formula"
        raise ValueError("No motif named '%s'" % motif)

    def generate_motiffile(self, motiffile):
        parser = MotifFilesParser()
        order = list()
        graphs = dict()
        formulae = dict()
        definedin = dict()
        for motif in chain(sorted(self.graphs), self.formula_list):
            if self.definedin[motif][1] == motiffile:
                order.append(motif)
                parser.definedin[motif] = (1, motiffile)
        for graph in set(order) & set(self.graphs):
            parser.graphs[graph] = self.graphs[graph]
        for formula in set(order) & set(self.formulae):
            parser.formulae[formula] = self.formulae[formula]
        parser.files.append(motiffile)
        parser.generate(motiffile, order)

    def _embeddings(self, graph, rnagf):
        dfscode = self.graphs[graph]
        if self.mode in ("none", "first"):
            embed = rnagf.match_code(dfscode)
            if self.mode == "none":
                return bool(embed)
            return [embed.as_embedding()] if embed else []
        embeds = [embed.as_embedding() for embed in rnagf.embed_code(dfscode)]
        return sorted(embeds, key=self._embedding_sorter)

    @staticmethod
    def _embedding_sorter(embedding):
        sort_by = list()
        for element in embedding:
            try:
                restuple = Residue.from_string(nodeid.split("=")[0]).entuple()
                chain, resloc = restuple[1], restuple[3:5]
                sort_by.append((chain, resloc, element))
            except:
                sort_by.append(element)
        return sort_by

    def _truth_values(self, assigns, formula):
        """
        Returns:
            A 2-tuple of (assigned, value), where:
              assigned: boolean, indicates if the formula
                is assigned, based on 'assigns' argument.
              value: depends on the 'mode' attribute, as:
                "all": See 'assign_motifs' under "all"
                "first": See 'assign_motifs' under "first"
                "none": boolean, equal to 'assigned'.
        """
        truth_values = list()
        for clause in self.formulae[formula]:
            truth_values.append(list())
            for hasterm, term in clause:
                truth_values[-1].append(hasterm == (term in assigns))
                if not truth_values[-1][-1] and self.mode != "all":
                    assigned = False
                    break
            else:
                if self.mode != "all":
                    assigned = True
                    break
        if self.mode == "all":
            assigned = any([all(values) for values in truth_values])
        return assigned, assigned if self.mode == "none" else truth_values

    def _resolve_graph_tree(self):
        children = defaultdict(list, {self._ROOT: list()})
        parentage = self._resolve_graph_parentage()
        for child in parentage:
            for parent in parentage[child]:
                children[parent].append(child)
        return children

    def _resolve_graph_parentage(self):
        ancestry = self._resolve_graph_ancestry()
        for graph in sorted(ancestry.keys()):
            for mom1, mom2 in combinations(sorted(ancestry[graph]), 2):
                if mom1 not in ancestry[graph] or mom2 not in ancestry[graph]:
                    continue
                granny1 = mom1 in ancestry[mom2] and len(ancestry[mom2]) > 1
                granny2 = mom2 in ancestry[mom1] and len(ancestry[mom1]) > 1
                if granny1 or granny2:
                    if granny1 and granny2:
                        ancestry[mom1].remove(mom2) # granny1 only
                    ancestry[graph].remove(mom1 if granny1 else mom2)
        return ancestry

    def _resolve_graph_ancestry(self):
        ancestors = defaultdict(list)
        graphs = self._by_complexity(self.graphs)
        for name1, code1 in graphs:
            graph1 = RNAStructureGraph.from_code(code1, name1)
            complexity1 = abs(code1)
            ancestors[name1].append(self._ROOT)
            for name2, code2 in graphs:
                if abs(code2) > complexity1:
                    break
                if name2 != name1 and graph1.embed_code(code2):
                    ancestors[name1].append(name2)
        return ancestors

    def _resolve_formula_order(self):
        if len(self.formulae) == 1:
            return list(self.formulae)
        dependency = self._inter_formula_dependencies()
        formula_order = sorted(set(self.formulae) - set(dependency))
        while dependency:
            for formula, dependencies in self._by_complexity(dependency):
                dependencies.difference_update(formula_order)
                if not dependencies:
                    formula_order.append(formula)
                    dependency.pop(formula)
        return formula_order

    def _inter_formula_dependencies(self):
        dependency = defaultdict(set)
        for formula, definition in self.formulae.items():
            for clause in definition:
                dependency[formula].update(refers for _, refers in clause)
            dependency[formula].difference_update(self.graphs)
            if not dependency[formula]:
                dependency.pop(formula)
        for formula, dependencies in dependency.items():
            if len(dependencies - set(self.formulae)):
                raise ValueError("Formula Dependency Parsing Error")
        return dependency

    @staticmethod
    def _by_complexity(mapping):
        if not mapping:
            return mapping.items()
        complexity = len
        if isinstance(list(mapping.values())[0], MotifGraphCode):
            complexity = abs
        sorter = lambda item: (complexity(item[1]), item[0])
        return sorted(mapping.items(), key=sorter)


class MotifFilesParser(object):

    class ParsingError(ValueError):
        """Base class for parsing errors"""

    class FileFormatError(ParsingError):
        MESG = "line %d in file '%s'"
        def __init__(self, filename, lineno):
            self.filename = filename
            self.lineno = lineno
            message = self.MESG % (lineno, filename)
            cls = MotifFilesParser.FileFormatError
            super(cls, self).__init__(message)

    class RepeatPatternError(ParsingError):
        MESG = "'%s' defined in files '%s' and '%s'"
        def __init__(self, defname, file1, file2):
            self.filenames = [file1, file2]
            message = self.MESG % (defname, file1, file2)
            cls = MotifFilesParser.RepeatPatternError
            super(cls, self).__init__(message)

    class BadReferenceError(ParsingError):
        MESG = "'%s' not defined before defining '%s'"
        def __init__(self, undefined, referrer):
            self.undefined = undefined
            self.referrer = referrer
            message = self.MESG % (undefined, referrer)
            cls = MotifFilesParser.BadReferenceError
            super(cls, self).__init__(message)

    class BadDefinitionError(ParsingError):
        MESG = "Can't parse '%s' in definition for '%s'"
        def __init__(self, unparsed, found_in):
            self.unparsed = unparsed
            self.found_in = found_in
            message = self.MESG % (unparsed, found_in)
            cls = MotifFilesParser.BadReferenceError
            super(cls, self).__init__(message)

    def __init__(self):
        self.files = list()
        self.graphs = dict()
        self.formulae = dict()
        self.definedin = dict()
        self.domain = RNALabels.LABELS

    def copy(self):
        copy = MotifFilesParser()
        copy.files = self.files[:]
        copy.domain = self.domain
        copy.graphs, copy.formulae, copy.definedin =\
          self.graphs.copy(), self.formulae.copy(), self.definedin.copy()
        return copy

    def extract(self, *motiffiles):
        copy = self.copy()
        for file in motiffiles:
            try:
                self._extract_from(file)
            except self.ParsingError as e:
                self.files, self.graphs, self.formulae, self.definedin\
                  = copy.files, copy.graphs, copy.formulae, copy.definedin
                raise e

    def generate(self, motiffile, order):
        out_string = list()
        for motif in order:
            if motif in self.graphs:
                header = "Graph:\t%s" % motif
                out_string.append(self.graphs[motif].to_string(header, "\t"))
            else:
                out_string.append("Formula:\t%s" % motif)
                lines = [" and ".join(line) for line in self.formulae[motif]]
                out_string += ["\t"+line for line in lines]
            out_string.append("")
        with open(motiffile, 'w') as outfile:
            outfile.write("\n".join(out_string))

    def _extract_from(self, file):
        defbytype = {"Graph": list(), "Formula": list()}
        for mtype, motif, deflist in self._read_motiffile(file):
            if motif in self.definedin:
                _, definedin = self.definedin[motif]
                raise self.RepeatPatternError(motif, file, definedin)
            self.definedin[motif] = (len(self.files)+1, file)
            defbytype[mtype].append((motif, deflist))
        self._parse_graphs(defbytype["Graph"])
        self._parse_formulae(defbytype["Formula"])
        self.files.append(file)

    @classmethod
    def _read_motiffile(cls, file):
        with open(file, 'r') as infile:
            mtype, motif, definition = None, None, None
            for lineno, line in chain(enumerate(infile, 1), [(-1, "")]):
                try:
                    line = line.split("!")[0].strip() # ignore comments
                    if line.startswith(("Graph:", "Formula:")):
                        mtype, motif = map(str.strip, line.split(":")[:2])
                        definition = list()
                    elif line:
                        definition.append(line)
                    else:
                        if mtype and motif and definition:
                            yield (mtype, motif, definition)
                        motif, mtype, definition = None, None, None
                except:
                    raise cls.FileFormatError(file, lineno)

    def _parse_graphs(self, graphdefs):
        for motif, refers, deflist in self._parse_graph_references(graphdefs):
            try:
                base = self.graphs[refers] if refers else None
                code = MotifGraphCode.from_lines(deflist, self.domain, base)
                self.graphs[motif] = code
            except KeyError as e:
                raise self.BadReferenceError(refers, motif)
            except DFSCodeFragment.InvalidFragmentError as e:
                raise self.BadDefinitionError(e.args[0], motif)
            except BaseException as e:
                raise ValueError("Error while parsing %s:\n%s" % (motif, e.message))

    def _parse_formulae(self, formuladefs):
        for motif, deflist in formuladefs:
            formula = self._parse_formula_clauses(deflist)
            for term in set(list(zip(*sum(formula, list())))[1]):
                if term not in chain(self.graphs, self.formulae):
                    raise self.BadReferenceError(term, motif)
            self.formulae[motif] = formula

    @staticmethod
    def _parse_graph_references(graphdefs):
        refs = list()
        for motif, deflist in graphdefs:
            if deflist[0].startswith("+"):
                refs.append((motif, deflist[0].split("\t")[1], deflist[1:]))
            else:
                refs.append((motif, None, deflist))
        return refs

    @staticmethod
    def _parse_formula_clauses(deflist):
        formula = list()
        for defline in deflist:
            defline = [term.strip() for term in defline.split(" and ")]
            clause = list()
            for term in defline:
                if term.startswith("not "):
                   clause.append((False, term[4:].strip()))
                else:
                    clause.append((True, term))
            formula.append(clause)
        return formula


class MotifGraphAnalyser(object):

    IS_PART_OF = RNALabels.EDGELABELS[RNALabels.IS_PART_OF]
    IN_SEQUENCE_WITH = RNALabels.EDGELABELS[RNALabels.IN_SEQUENCE_WITH]
    SEQUENCE_EDGES = [IS_PART_OF, IN_SEQUENCE_WITH]
    # 
    SEQUENCE_NODES = [RNALabels.ANY_NUCLEOTIDE, RNALabels.PHOSPHATE]
    SEQUENCE_NODES = map(RNALabels.NODELABELS.__getitem__, SEQUENCE_NODES)
    ANY_NUCLEOTIDE, PHOSPHATE = SEQUENCE_NODES

    def __init__(self, motifgraph, assigner):
        self.code = assigner.graphs[motifgraph]
        self.graph = motifgraph
        self.index = None
        self.sequence = dict()
        self.assigner = assigner
        self._residues = dict()
        self._segments = list()
        self._direction = dict()
        self._general_nodes = dict()
        self._general_edges = dict()
        self._sorted_items = list()
        self.__setup()

    @staticmethod
    def _read_limit_cases(limitfile):
        limited = set()
        if limitfile:
            with open(limitfile, 'r') as infile:
                for index, line in enumerate(infile):
                    if not index and line.startswith(("Formula", "Graph")):
                        continue
                    case = " ".join(line.strip("\n").split("\t")[:2])
                    limited.add(case)
        return limited

    def _set_index(self, assigns):
        header = assigns._graph_header_details(self.graph)
        self.index = {str_: indx for indx, str_ in enumerate(header)}

    @classmethod
    def _parse_nodename(cls, embed, hideresidue=True):
        if "=" in embed:
            parts = list()
            for part in embed.split("="):
                parts.append(cls._parse_nodename(part, False).split("."))
            return " ".join([":".join(pair) for pair in zip(*parts)])
        elif "." in embed:
            residue, part = embed.split(".")
            residue = cls._parse_nodename(residue)
            if hideresidue and residue in ["A", "C", "G", "U"]:
                return part
            return "%s.%s" % (residue, part)
        else:
            return embed.split("]")[0].strip("[")

    @classmethod
    def _parse_residue(cls, embed):
        if "=" in embed:
            pair = list()
            for part in embed.split("="):
                pair.append(cls._parse_nodename(part.split(".")[0]))
            return ":".join(pair)
        elif "." in embed:
            return cls._parse_nodename(embed.split(".")[0])
        else:
            return cls._parse_nodename(embed)

    def __setup(self):
        redundant = self.__remove_redundant()
        for isnode, tagstr, tagtuple in self.code._iter_tagstr():
            if isnode:
                nodetag, nodelab = tagtuple
                if nodetag not in redundant and\
                  nodelab.name in RNALabels.RESIDUES:
                    self._residues[tagstr] = defaultdict(list)
                if nodelab.is_general:
                    self._general_nodes[tagstr] = defaultdict(list)
            else:
                nodetag, edontag, edgelab = tagtuple
                if edgelab.is_general:
                    self._general_edges[tagstr] = defaultdict(list)
        self.__set_segments()
        for segment in self.__directionless():
            segment = segment if segment[0] < segment[-1] else segment[::-1]
            segstr = " - ".join(map(str, segment))
            self._direction[segstr] = defaultdict(list)
        self.__sort_items()
       
    def __remove_redundant(self):
        graph = RNAStructureGraph.from_code(self.code)
        def check(key):
            return (key in RNALabels.RESIDUES, 
                    key == RNALabels.PHOSPHATE, 
                    key in RNALabels.ATOMS)
        redundant = set()
        for nodetag, edontag, _ in graph.edges_with_label_in([self.IS_PART_OF]):
            nodelab = graph.nodes[nodetag].name
            edonlab = graph.nodes[edontag].name
            nodeis, edonis = check(nodelab), check(edonlab)
            redundant.add(edontag if nodeis > edonis else nodetag)
        for nodetag, nodelab in enumerate(self.code._tags):
            if nodelab.name in RNALabels.RESIDUES or\
              nodelab.name in RNALabels.BASEPAIR_FAMILIES:
                redundant.add(nodetag)
        return redundant

    def __set_segments(self):
        graph = RNAStructureGraph.from_code(self.code)
        sequence_nodes = set(graph.nodes_with_label_in(self.SEQUENCE_NODES))
        seqedgelabels = self.SEQUENCE_EDGES
        fragments = set()
        for nodetag, edontag, _ in graph.edges_with_label_in(seqedgelabels):
            if nodetag in sequence_nodes and edontag in sequence_nodes:
                if graph.nodes[nodetag] == self.PHOSPHATE:
                    nodetag, edontag = edontag, nodetag # node is nucleotide
                fragments.add((nodetag, edontag))
        fragments = sorted(fragments)
        self._merge_fragments(fragments)
        sequence_nodes -= set(sum(fragments, tuple()))
        for nodetag in sequence_nodes:
            self._segments.append([nodetag])

    def _merge_fragments(self, fragments):
        for fragment in fragments:
            added_to = list()
            for segex, segment in enumerate(self._segments):
                if segment[0] == fragment[0]:
                    segment.insert(0, fragment[1])
                    added_to.append(segex)
                elif segment[0] == fragment[1]:
                    segment.insert(0, fragment[0])
                    added_to.append(segex)
                if segment[-1] == fragment[0]:
                    segment.append(fragment[1])
                    added_to.append(segex)
                elif segment[-1] == fragment[1]:
                    segment.append(fragment[0])
                    added_to.append(segex)
            if len(added_to) == 0:
                self._segments.append(list(fragment))
            elif len(added_to) == 2:
                segB = self._segments.pop(added_to[1])
                segA = self._segments.pop(added_to[0])
                if segA[0] == segB[-2]:
                    self._segments.append(segB[:-2]+segA)
                elif segA[1] == segB[-1]:
                    self._segments.append(segB+segA[2:])
                elif segA[-2] == segB[0]:
                    self._segments.append(segA[:-2]+segB)
                elif segA[-1] == segB[1]:
                    self._segments.append(segA+segB[2:])

    def __directionless(self):
        directionless = self._segments[:]
        directioned = set()
        graph = RNAStructureGraph.from_code(self.code)
        sequence_nodes = set(graph.nodes_with_label_in(self.SEQUENCE_NODES))
        for nodetag, edontag, _ in graph.edges_with_label_in([self.IS_PART_OF]):
            if nodetag in sequence_nodes and edontag in sequence_nodes:
                directioned.add(nodetag)
                directioned.add(edontag)
        hasdirection = list()
        for segex, segment in enumerate(directionless):
            if set(segment) & directioned or len(segment) == 1:
                hasdirection.append(segex)
        for segex in reversed(hasdirection):
            directionless.pop(segex)
        return directionless

    def __sort_items(self):
        keyvalues = list()
        for segstr in self._direction:
            tags = map(int, [atag.strip() for atag in segstr.split("-")])
            key_ = ("sequence", "goes") + tuple(tags)
            keyvalues.append((key_, segstr))
        for edgestr in self._general_edges:
            tags = edgestr.split(" ")[1].strip("()").split("->")
            key_ = ("edge", "is") + tuple(sorted(map(int, tags)))
            keyvalues.append((key_, edgestr))
        for nodestr in self._general_nodes:
            key_ = ("node", "is", int(nodestr.split(" ")[0]))
            keyvalues.append((key_, nodestr))
        for nodestr in self._residues:
            key_ = ("node", "of", int(nodestr.split(" ")[0]))
            keyvalues.append((key_, nodestr))
        self._sorted_items = list(zip(*sorted(keyvalues)))[1]

    def _sort_itemset(self, itemset):
        sortedset = list()
        for typestr in self._sorted_items:
            for item in itemset:
                if item.startswith(typestr):
                    sortedset.append(item)
        return tuple(sortedset)


class MotifDivideFile(object):

    HEADER = "Division %d\t-\t%d\t%d"

    @classmethod
    def write_divisions(cls, outfile, divisions, mode='w'):
        out_string = ["="*50]
        for division in sorted(divisions):
            details, cases = divisions[division]
            uniques = len(set(zip(*cases)[0]))
            out_string += ["", cls.HEADER % (division, len(cases), uniques)]
            if details:
                out_string += ["\t%s" % line for line in details]
            out_string += ["\t\t%s" % case for case in map("\t".join, cases)]
        with open(outfile, mode) as outfile:
            outfile.write("\n".join(out_string)+"\n")

    @classmethod
    def read_divisions(cls, dividefile, only=None):
        with open(dividefile, 'r') as infile:
            line = ""
            try:
                while not (line.startswith("=") and line.endswith("=")):
                    line = next(infile).strip("\n")
                divisions = cls._read_divisions(infile, only)
            except StopIteration:
                divisions = dict()
        if only:
            details, cases = divisions.pop(only)
            return tuple(details), cases
        return divisions

    @classmethod
    def _read_divisions(cls, fileobj, only):
        divisions = dict()
        indivision = None
        for line in fileobj:
            line = line.strip("\n")
            stripped = line.strip()
            if line.startswith("Division "):
                if indivision and only: # previous indivision
                    break
                indivision = int(line.split()[1])
                if only and indivision != only:
                    indivision = None
                    continue
                divisions[indivision] = details, cases = list(), list()
            elif not line.strip():
                continue
            elif indivision and line.startswith("\t\t"):
                cases.append(stripped.split("\t"))
            elif indivision and (not stripped or line.startswith("\t")):
                details.append(stripped)
        for divno, (details, cases) in list(divisions.items()):
            divisions[divno] = (tuple(details), cases)
        return divisions

    @classmethod
    def example(cls, dividefile, resolutionfile):
        resolutions = cls.read_resolutions(resolutionfile)
        divisions = cls.read_divisions(dividefile)
        return cls._example(divisions, resolutions)

    @classmethod
    def _example(cls, divisions, resolutions):
        examples = dict()
        for division in divisions:
            cases = zip(*divisions[division][1])[0]
            best = min(cases, key=lambda case: resolutions[case.split()[0]])
            resolution = resolutions[best.split()[0]]
            examples[division] = (best, resolution)
        return examples
    
    @classmethod
    def read_resolutions(cls, resolutionfile):
        resolutions = list()
        with open(resolutionfile, 'r') as infile:
            for line in infile:
                pdbid, resolution = line.split()
                resolutions.append((pdbid, float(resolution)))
        resolutions = dict(resolutions)
        return resolutions


class MotifRedivideFile(MotifDivideFile):

    WARN = True
    CONVERT = dict()
    CONVERT["R"] = {"R", "G", "A", "g", "a"}
    CONVERT["Y"] = {"Y", "U", "C", "u", "c"}
    for N in ("A", "C", "G", "U"):
        CONVERT[N] = {N}
        CONVERT[N.lower()] = {N, N.lower()}
    CONVERT["N"] = set().union(*CONVERT.values())

    @classmethod
    def ensure_redividefile(cls, filepath):
        for details, caselines in cls.read_divisions(filepath).values():
            for line in details:
                if line.startswith("Organisms:\t"):
                    break
            else:
                raise ValueError("Not in redivide file format")

    @staticmethod
    def get_cutoff(minseqs, length, cutoff):
        if not length and cutoff >= 1:
            return 15
        if cutoff == 1:
            return length
        cutoff = sum(map(len, minseqs)) / float(cutoff)
        for length in count(length):
            minseqs = filter(lambda seqce: len(seqce) >= length, minseqs)
            if length ** 4 > (cutoff - length*len(minseqs)):
                return length
    
    @classmethod
    def common_sequences(cls, redividefile, length=0, cutoff=1):
        cls.ensure_redividefile(redividefile)
        divisions = cls.read_divisions(redividefile)
        bydivno = defaultdict(list)
        minseqs = set()
        for divno, (details, _) in divisions.items():
            for line in details:
                if line.startswith("Organisms"):
                    break
                molecule, sequence = line.split("\t")[:2]
                bydivno[divno].append((molecule, sequence))
                minseqs.add((molecule, sequence))
        length = cls.get_cutoff(minseqs, length, cutoff)
        if cls.WARN:
            print("Using cut-off length %d" % length, file=sys.stderr)
        mergers = defaultdict(set)
        for molecule, minseq in sorted(minseqs, key=len, reverse=True):
            if len(minseq) < length:
                break
            for endex in range(length, len(minseq)):
                stretch = minseq[endex-length:endex]
                divnos = set()
                for divno, div_sequences in bydivno.items():
                    for div_molecule, div_minseq in div_sequences:
                        if molecule == div_molecule and stretch in div_minseq:
                            divnos.add(divno)
                if len(divnos) >= 2:
                    for merger in combinations(sorted(divnos), 2):
                        mergers[merger].add(stretch)
        for merger in set(mergers):
            mergers[merger] = cls._min_common(mergers.pop(merger))
        return mergers

    @classmethod
    def _min_common(cls, sequences):
        minsequences = list()
        for sequence in sequences:
            for minsequence in minsequences[:]:
                if cls.matches(minsequence, sequence):
                    break
                if cls.matches(sequence, minsequence)\
                  and sequence != minsequence:
                    minsequences.remove(minsequence)
            else:
                minsequences.append(sequence)
        minsequences.sort(key=lambda sequence: (len(sequence), sequence))
        return minsequences

    @classmethod
    def matches(cls, shorter, longer):
        if len(shorter) <= len(longer):
            longer = map(lambda char: cls.CONVERT[char], longer)
            shorter = map(lambda char: cls.CONVERT[char], shorter)
            for index in range(1+len(longer)-len(shorter)):
                for shortchar, longchar in zip(shorter, longer[index:]):
                    if not longchar.issubset(shortchar):
                        break
                else:
                    return True
        return False


class MotifClusterFile(MotifRedivideFile):

    @classmethod
    def ensure_clusterfile(cls, filepath):
        for details, caselines in cls.read_divisions(filepath).values():
            for line1, line2 in zip(details, details[1:]):
                if line1.startswith("Organisms:\t")\
                  and line2.startswith("Merges:\t"):
                    break
            else:
                raise ValueError("Not in cluster file format")
    
    @classmethod
    def from_redividefile(cls, redividefile, clusterfile):
        cls.ensure_redividefile(redividefile)
        clusters = dict()
        divisions = cls.read_divisions(redividefile)
        for divno, (details, caselines) in divisions.items():
            for lineno, line in enumerate(details):
                if line.startswith("Organisms:\t"):
                    break
            details = list(details[:lineno+1]) + ["Merges:\t%d" % divno]
            clusters[divno] = (tuple(details), caselines)
        cls.write_divisions(clusterfile, clusters, 'w')

    @classmethod
    def cluster_divisions(cls, clusterfile, clusters):
        cls.ensure_clusterfile(clusterfile)
        clusters = cls.simplify_clusters(clusters)
        divisions = cls.read_divisions(clusterfile)
        for newdivno, cluster in enumerate(clusters, 1):
            newcaselines = list()
            newdetails = list()
            organisms = set()
            merges = set()
            for divno in cluster:
                details, caselines = divisions.pop(divno)
                newcaselines += caselines
                newdetails += list(details[:-2])
                organisms.update(details[-2].split("\t")[1].split(" ; "))
                merges.update(details[-1].split("\t")[1].split(" ; "))
            merges = map(str, sorted(map(int, merges)))
            newcaselines = sorted(newcaselines)
            newdetails.append("Organisms:\t%s" % " ; ".join(sorted(organisms)))
            newdetails.append("Merges:\t%s" % " ; ".join(merges))
            divisions[-newdivno] = (tuple(newdetails), newcaselines)
        divisions = cls.sort_divisions(divisions)
        cls.write_divisions(clusterfile, divisions, 'w')

    @staticmethod
    def simplify_clusters(clusters):
        clusters = list(map(set, clusters))
        previous = float("inf")
        while len(clusters) < previous:
            previous = len(clusters)
            todelete = set()
            for index1, cluster1 in enumerate(clusters, 1):
                if index1-1 in todelete:
                    continue
                for index2, cluster2 in enumerate(clusters[index1:], index1):
                    if cluster1 & cluster2:
                        cluster1.update(cluster2)
                        todelete.add(index2)
            for index in sorted(todelete, reverse=True):
                clusters.pop(index)
        return [sorted(cluster) for cluster in clusters]

    @classmethod
    def sort_divisions(cls, divisions):
        by_count = list()
        for division in divisions.values():
            cases = zip(*division[1])[0]
            uniques, count = len(set(cases)), len(cases)
            by_count.append((uniques, count, division))
        by_count = sorted(by_count, reverse=True)
        divisions = dict(enumerate(zip(*by_count)[2], 1))
        return divisions


class MotifDivider(MotifGraphAnalyser):

    def __init__(self, motifgraph, assigner, restrict=None, limited=None):
        super(MotifDivider, self).__init__(motifgraph, assigner)
        self.cases = list()
        self.embeds = list()
        self.itemsets = defaultdict(list)
        self.restrict = restrict
        self.limited = self._read_limit_cases(limited)
        self.results = None

    def add_embeddings(self, assignsfiles, keepembeds=False):
        # if not (self._general_edges or self._general_nodes):
        #     return
        load = MotifAssignments.load
        limit = {self.graph, self.restrict} - {None}
        assigner = self.assigner
        for assignsfile in assignsfiles:
            assigns = load(assigner, assignsfile, True, limit)
            for case, embeddings in assigns.by_graph[self.graph].items():
                if not self.index:
                    header = assigns._graph_header_details(self.graph)
                    self.index = {str_: indx for indx, str_ in enumerate(header)}
                if self.restrict and self.restrict not in assigns.by_case[case]:
                    continue
                if self.limited and case not in self.limited:
                    continue
                for embedding in embeddings:
                    self._add_embedding(case, embedding, keepembeds)
            del(assigns)

    def add_sequences(self, sitesfiles, engraphertype):
        if not self.index:
            return
        indices = self._sequence_indices()
        for pdbid, sitesfile in self._sitesfiles(sitesfiles):
            casexes, uniques = list(), set()
            for casex, case in enumerate(self.cases):
                if case.startswith(pdbid):
                    casexes.append(casex)
                    uniques.add(case)
            uniques = sorted(uniques)
            seqss = engraphertype.extract_sequences(sitesfile, uniques)
            seqss = dict(zip(uniques, seqss))
            for casex in casexes:
                embed = [self.embeds[casex][index] for index in indices]
                caseseqs = seqss[self.cases[casex]]
                self._add_case_sequences(casex, embed, caseseqs)

    def _sitesfiles(self, sitesfiles):
        sites = dict()
        pdbids = {case.split()[0] for case in self.cases}
        for sitesfile in sitesfiles:
            pdbid = os.path.basename(sitesfile).split("_")[0]
            if pdbid in pdbids:
                if pdbid in sites:
                    raise ValueError("Multiple files for PDB ID %s" % pdbid)
                sites[pdbid] = sitesfile
        if len(pdbids) != len(sites):
            difference = ", ".join(pdbids - set(sites))
            raise ValueError("No sites for PDB IDs %s" % difference)
        return [(pdbid, sites[pdbid]) for pdbid in sorted(sites)]

    def _sequence_indices(self):
        indices = list()
        for segment in self._segments:
            nodetag = segment[0]
            nodelab = self.code._tags[nodetag]
            indices.append(self.index["%d %s" % (nodetag, nodelab)])
        return indices

    def _add_case_sequences(self, casex, embed, caseseqs):
        self.sequence[casex] = list()
        for nodeid in embed:
            nodeid = nodeid.split(".")[0]
            chainid = nodeid.split(":")[1]
            for baseseq, idseq in caseseqs:
                if nodeid in idseq:
                    self.sequence[casex].append("%s:%s" % (chainid, baseseq))
                    break
            else:
                base = nodeid.split("]")[0].strip("[")
                self.sequence[casex].append("%s:%s" % (chainid, base))

    def summarise(self, outfile=None):
        if not self.index:
            return
        out_string = list()
        self._summarise(out_string, "Direction", self._direction)
        self._summarise(out_string, "Specific Edges", self._general_edges)
        self._summarise(out_string, "Specific Nodes", self._general_nodes)
        self._summarise(out_string, "Residue", self._residues)
        out_string = "\n".join(out_string)+"\n"
        if outfile:
            with open(outfile, 'a') as outfile:
                outfile.write(out_string+"\n")
        else:
            print(out_string)

    def divide(self, outfilepath, support=-1, show=True, sequence=False):
        if not len(self.itemsets):
            print("No instances")
            return
        self.results = self._divide(support)
        total, uniques = len(self.cases), len(set(self.cases))
        out_string = "Divisions at %s support for %d cases (%d unique)"
        out_string = [out_string % (support, total, uniques)]
        out_string = "\n".join(out_string)
        with open(outfilepath, 'a') as outfile:
            outfile.write(out_string+"\n")
        divisions = self._summarise_divisions(out_string, show, sequence)
        MotifDivideFile.write_divisions(outfilepath, divisions, "a")

    def _add_embedding(self, case, embedding, keepembeds):
        itemset = list()
        casex = len(self.cases)
        self.cases.append(case)
        if keepembeds:
            self.embeds.append(embedding)
        self._add_embed_segments(embedding, itemset, casex)
        self._add_embed_generals(embedding, itemset, casex)
        self._add_embed_residues(embedding, itemset, casex)
        itemset = self._sort_itemset(itemset)
        self.itemsets[itemset].append(casex)

    def _summarise(self, out_string, header, attribute):
        if not attribute:
            return
        sorter = lambda item: item[1]
        out_string.append(header)
        if header != "Direction":
            for tagstr, _ in sorted(self.index.items(), key=sorter):
                if tagstr not in attribute:
                    continue
                self._details(out_string, tagstr, attribute[tagstr])
        else:
            for segstr in self._sorted_items:
                if " - " in segstr:
                    self._details(out_string, segstr, attribute[segstr])
        out_string.append("")

    def _summarise_divisions(self, out_string, show, sequence):
        divisions = dict()
        for divno, (itemset, _) in enumerate(self.results, 1):
            itemset = self._sort_itemset(itemset)
            cases = list()
            if show:
                for casex in self.itemsets[itemset]:
                    cases.append([self.cases[casex]])
                    if sequence:
                        cases[-1].append("=".join(self.sequence[casex]))
                assert cases, "Divno %d" % divno
            divisions[divno] = (tuple(itemset), cases)
        return divisions

    def _divide(self, support):
        results = list()
        divisor = 1 if support <= 0 else float(len(self.cases))
        support = abs(support)
        for itemset in self.itemsets:
            netsupport = len(set(self.itemsets[itemset])) / divisor
            if netsupport >= support:
                casexes = self.itemsets[itemset]
                uniques = len(set(self.cases[casex] for casex in casexes))
                results.append((itemset, (netsupport, uniques)))
        results.sort(key=lambda entry: (entry[1], entry[0]), reverse=True)
        return results

    def _details(self, out_string, tagstr, details):
        freqs = dict()
        for spec, casexes in details.items():
            cases = [self.cases[index] for index in casexes]
            freqs[spec] = (len(set(cases)), len(cases))
        sorter = lambda item: item[1]
        out_string.append("\t%s" % tagstr)
        for case, freq in sorted(freqs.items(), key=sorter, reverse=True):
            out_string.append("\t\t%s\t%d\t%d" % (case, freq[0], freq[1]))
        out_string.append("")

    def _add_embed_segments(self, embedding, itemset, casex):
        for segstr in self._direction:
            tags = map(int, [atag.strip() for atag in segstr.split("-")])
            nodestr = "%s %s" % (tags[0], self.code._tags[int(tags[0])])
            edonstr = "%s %s" % (tags[1], self.code._tags[int(tags[1])])
            residue1 = embedding[self.index[nodestr]].split(".")
            residue2 = embedding[self.index[edonstr]].split(".")
            if len(residue1) == 2:
                mode = ">" if residue1[0] == residue2[0] else "<"
            else:
                mode = "<" if residue1[0] == residue2[0] else ">"
            self._direction[segstr][mode].append(casex)
            itemset.append("%s goes %s" % (segstr, mode))

    def _add_embed_generals(self, embedding, itemset, casex):
        for edgestr in self._general_edges:
            edgelab = embedding[self.index[edgestr]]
            self._general_edges[edgestr][edgelab].append(casex)
            itemset.append("%s is %s" % (edgestr, edgelab))
        for nodestr in self._general_nodes:
            embedpart = embedding[self.index[nodestr]]
            nodename = self._parse_nodename(embedpart).split(" ")[-1]
            self._general_nodes[nodestr][nodename].append(casex)
            itemset.append("%s is %s" % (nodestr, nodename))

    def _add_embed_residues(self, embedding, itemset, casex):
        for nodestr in self._residues:
            embedpart = embedding[self.index[nodestr]]
            residue = self._parse_residue(embedpart)
            self._residues[nodestr][residue].append(casex)
            itemset.append("%s of %s" % (nodestr, residue))
            residueid = embedpart.split(".")[0]
            if residueid in embedding:
                for index, debmepart in enumerate(embedding):
                    if debmepart == residueid:
                        for tagstr, header_index in self.index.items():
                            if header_index == index:
                                break
                        else:
                            raise ValueError("Residue not in header")
                        nodetag = tagstr.split(" ")[0]
                        self._residues[nodestr][nodetag].append(casex)
                        itemset.append("%s of %s" % (nodestr, nodetag))


class MotifRedivider(object):

    def __init__(self, entities):
        self.chainnames = dict()
        self._chain_details_folder = entities

    def columnwise_minimum_sequences(self, subdivisions):
        """
        Returns:
            A list, where each entry corresponds to a
            sequence column in the original divide file.
            Each entry in the list is a list of minimum
            sequences associated with that column.
        """
        subdivisions = sorted(subdivisions)
        sequences = [list(zip(*division)[1]) for division in subdivisions]
        sequences = list(zip(*sequences)) # row-wise to column-wise sequences
        minimum_sequences = list()
        for column in range(len(sequences)):
            column_minseqs = MotifRedivideFile._min_common(sequences[column])
            minimum_sequences.append(column_minseqs)
        return minimum_sequences

    def regroup_by_minimum_sequences(self, subdivisions):
        """
        Returns:
            A mapping from an n-tuple corresponding to
            a minimum sequence for each of the n sequence
            columns in the original divide file, to a set
            of subdivisions. Subdivisions are defined by
            n-tuple of (molecule_name, sequence) 2-tuples.
        """
        simplified = defaultdict(set)
        minseqss = self.columnwise_minimum_sequences(subdivisions)
        for subdivision in subdivisions:
            minrows = list()
            for row_minseqs in product(*minseqss):
                for (_, seq), minseq in zip(subdivision, row_minseqs):
                    if not MotifRedivideFile.matches(minseq, seq):
                        break
                else:
                    lengths = [len(minseq) for minseq in row_minseqs]
                    swapped = [minseq.swapcase() for minseq in row_minseqs]
                    minrows.append((tuple(lengths), tuple(swapped)))
            _, swapped = min(minrows)
            unswapped = tuple([minseq.swapcase() for minseq in swapped])
            simplified[unswapped].add(subdivision)
        return dict(simplified)

    def simplify_groupset(self, subdivisions):
        """
        Returns:
            A mapping from an n-tuple of 2-tuples, each of
            the form (molecule_name, minimum_sequence),
            where n is the number of sequence columns in
            the original divide file - to a set of
            subdivisions defined by similar n-tuples with
            the original sequences.
        """
        groups = dict()
        molecules = [molecule for molecule, _ in list(subdivisions)[0]]
        simplified = self.regroup_by_minimum_sequences(subdivisions)
        for minrow in simplified:
            group = list()
            for molecule, minseq in zip(molecules, minrow):
                group.append((molecule, minseq))
            groups[tuple(group)] = simplified[minrow]
        return groups

    def regroup_subdivisions(self, subdivisions):
        groups = defaultdict(set)
        for subdivision in subdivisions:
            group = [molecule for molecule, _ in subdivision]
            for (_, seq1), (_, seq2) in combinations(subdivision, 2):
                group.append(seq1 == seq2)
            groups[tuple(group)].add(subdivision)
        mingroups = dict()
        for group in groups:
            mingroups.update(self.simplify_groupset(groups[group]))
        # 
        groups = defaultdict(lambda: defaultdict(set))
        for mingroup, mingroup_subdivisions in mingroups.items():
            for subdivision in mingroup_subdivisions:
                for division, caselines in subdivisions[subdivision].items():
                    groups[mingroup][division].update(caselines)
            assert groups[mingroup], "Empty %s" % mingroup
        return groups

    def subdivide(self, division_details):
        """
        Returns:
        """
        subdivisions = defaultdict(lambda: defaultdict(set))
        for division, caselines in division_details:
            for case, chain_details in self._get_chain_details(caselines):
                subdivision = list()
                for _, sequence, molecule, organism in chain_details:
                    subdivision.append((molecule, sequence))
                subdivision = tuple(subdivision)
                subdivisions[subdivision][division].add((case, chain_details))
        subdivisions = self.regroup_subdivisions(subdivisions)
        return subdivisions

    def redivide(self, dividefile, redividefile):
        divisions = MotifDivideFile.read_divisions(dividefile)
        redivisions = defaultdict(dict)
        for divno in sorted(divisions):
            subdivisions = self.subdivide([divisions[divno]])
            assert(subdivisions)
            for subdivision, divisionwise in subdivisions.items():
                for itemset, cases_detailed in divisionwise.items():
                    redivisions[subdivision][itemset] = cases_detailed
        redivisions = self.regroup_subdivisions(redivisions)
        self.write_redivisions(redivisions, redividefile)

    def write_redivisions(self, redivisions, redividefile):
        out_string = ["="*50+"\n"]
        redivisions = self._finalise_redivisions(redivisions)
        sorted_divisions = self._sort_redivisions(redivisions)
        for divno, (key, counts) in enumerate(sorted_divisions, 1):
            out_string.append("Division %d\t-\t%d\t%d" % ((divno,)+counts))
            out_string += ["\t%s" % line.strip() for line in key]
            # 
            all_organisms = set()
            for division, cases_detailed in redivisions[key].items():
                organisms = zip(*sum(zip(*cases_detailed)[1], tuple()))[3]
                all_organisms.update(organisms)
            all_organisms = " ; ".join(sorted(all_organisms))
            out_string.append("\tOrganisms:\t%s" % all_organisms)
            # 
            subdivisions_items = list(redivisions[key].items())
            subdivisions_items.sort(key=lambda entry: -len(entry[1]))
            for division, cases_detailed in subdivisions_items:
                cases_detailed = sorted(cases_detailed)
                out_string.append("\t")
                out_string += ["\t%s" % line.strip() for line in division]
                out_string += list(map(self._to_caseline, cases_detailed))
            out_string.append("\n")
        with open(redividefile, 'w') as outfile:
            outfile.write("\n".join(out_string))

    def _finalise_redivisions(self, subdivisions):
        redivisions = subdivisions.copy()
        for subdivision in subdivisions:
            details = list()
            for index, (molecule, sequence) in enumerate(subdivision):
                line = [molecule, sequence]
                for subindex in range(index):
                    isequal = sequence == subdivision[subindex][1]
                    line.append("=" if isequal else "!")
                details.append("\t%s" % "\t".join(line))
            redivisions[tuple(details)] = redivisions.pop(subdivision)
        return redivisions

    def _sort_redivisions(self, redivisions):
        bycount = defaultdict(list)
        for subdivision in redivisions:
            cases = list()
            for itemset, cases_detailed in redivisions[subdivision].items():
                cases += list(list(zip(*cases_detailed))[0])
            bycount[subdivision] = (len(set(cases)), len(cases))
        bycount = sorted(bycount.items(), key=lambda entry: entry[::-1])[::-1]
        return bycount

    def _get_chain_details(self, caselines):
        chain_details = list()
        for caseline in caselines:
            case = caseline[0]
            pdbid = case.split(" ")[0]
            chainids, sequences = list(), list()
            organisms, molecules = list(), list()
            for chain_sequence in caseline[1].split("="):
                chainid, sequence = chain_sequence.split(":")
                molecule, organism = self._get_chainname(pdbid, chainid)
                chainids.append(chainid)
                sequences.append(sequence)
                molecules.append(molecule)
                organisms.append(organism)
            chain_data = tuple(zip(chainids, sequences, molecules, organisms))
            chain_details.append((case, chain_data))
        return sorted(chain_details)

    @staticmethod
    def read_chains_file(chainsfile):
        entities = dict()
        with open(chainsfile, 'r') as infile:
            for line in infile:
                line = line.strip("\n").split("\t")
                chainid, _, molecule, organism = line[:4]
                entities[chainid] = (molecule, organism)
        return entities

    def _get_all_chainnames(self, pdbid):
        prefix = "%s_" % pdbid
        for filename in os.listdir(self._chain_details_folder):
            if filename.startswith(prefix):
                break
        else:
            mesg = "%%s entry not in %s" % self._chain_details_folder
            raise ValueError(mesg % pdbid)
        filepath = os.path.join(self._chain_details_folder, filename)
        self.chainnames[pdbid] = self.read_chains_file(filepath)

    def _get_chainname(self, pdbid, chainid):
        if pdbid not in self.chainnames:        
            self._get_all_chainnames(pdbid)
        return self.chainnames[pdbid][chainid]

    def _to_caseline(self, case_chain_details):
        case, chain_details = case_chain_details
        chainids, sequences = list(zip(*chain_details))[:2]
        caseline = [":".join(entry) for entry in zip(chainids, sequences)]
        caseline = "\t\t%s\t%s" % (case, "=".join(caseline))
        return caseline


class MotifCombiner(object):

    class AddToAugMap(object):

        def __init__(self, mapping, backmap, code, common, extras):
            self.to = mapping
            self.fro = backmap
            self.code = code
            self.common = common
            self.extras = extras

        def _update_dictionary(self, totagstr_pairs, dictionary, remapping):
            for dtagstr, totagstr in dictionary.items():
                if "->" not in dtagstr:
                    totag, nodelab = totagstr.split(" ")
                    totag = int(totag)
                    if totag in remapping:
                        totagstr = "%d %s" % (remapping[totag], nodelab)
                        dictionary[dtagstr] = totagstr
                    continue
                nodetags = totagstr.split(" ")[1].strip("()")
                nodetag, edontag = map(int, nodetags.split("->"))
                if nodetag in remapping:
                    newstr = "(.%d.-" % remapping[nodetag]
                    totagstr = totagstr.replace("(%d-" % nodetag, newstr)
                    totagstr = totagstr.replace("(.", "(").replace(".-", "-")
                if edontag in remapping:
                    newstr = ">.%d.)" % remapping[edontag]
                    totagstr = totagstr.replace(">%d)" % edontag, newstr)
                    totagstr = totagstr.replace(">.", ">").replace(".)", ")")
                if totagstr != dictionary[dtagstr]:
                    dictionary[dtagstr] = totagstr
                for right, wrong in totagstr_pairs:
                    if totagstr in (right, wrong):
                        if totagstr == wrong:
                            dictionary[dtagstr] = right
                        break
                else:
                    raise ValueError("AddTagStr '%s' not found" % dtagstr)

        def _update_tagstr(self, remapping):
            totagstr_pairs = list()
            tagstr_iter = self.code._iter_tagstr()
            for isnode, tagstr, args in tagstr_iter:
                if not isnode:
                    nodetag, edontag, edgelab = args
                    totagstr = "%s (%%d->%%d)" % edgelab.name
                    totagstr_right = totagstr % (nodetag, edontag)
                    totagstr_wrong = totagstr % (edontag, nodetag)
                    totagstr_pairs.append((totagstr_right, totagstr_wrong))
            self._update_dictionary(totagstr_pairs, self.extras, remapping)
            self._update_dictionary(totagstr_pairs, self.common, remapping)

        def update_code(self, code, remap):
            self.code = code
            remapping = dict()
            frodict = self.fro.copy()
            for nodetag, nodeid in enumerate(remap.nodeids):
                if nodetag != nodeid:
                    remapping[nodeid] = nodetag
                    if nodeid in frodict:
                        fromtag = frodict[nodeid]
                        self.fro[nodetag] = fromtag
                        self.to[fromtag] = nodetag
            self._update_tagstr(remapping)

        def __eq__(self, fles):
            return self.code == fles.code

    
    def __init__(self, name, base, augend, addend, assigner, restrict=None):
        """
            base must be a part of augend
            if base occurs in augend or addend, it can occur only once
        """
        self.name = name
        self.base = base
        self.augend = MotifGraphAnalyser(augend, assigner)
        self.addend = MotifGraphAnalyser(addend, assigner)
        self.combined = dict()
        self.restrict = None
        self.assigner = None
        self.assigns = None
        self.index = None # note this is many to one, not bijective like Analyser
        self.__extend()

    def __extend(self):
        mappable, basemap, addedges, augedges = self.__get_mappable()
        auggf = RNAStructureGraph.from_code(self.augend.code)
        firstid = max(auggf.nodes)+1
        mappings = self.__generate_mappings(mappable, basemap, firstid)
        combined = list()
        for mapping in mappings:
            combinedgf  = auggf.copy("")
            self.__extend_nodes(combinedgf, mapping, firstid)
            self.__extend_edges(combinedgf, mapping, addedges, augedges)
            try:
                code = MotifGraphCode(combinedgf.mindfscode(self.augend.code))
                mapping.update_code(code, combinedgf.match_code(code))
                if mapping not in combined:
                    combined.append(mapping)
            except ValueError:
                continue
        formatter = "%%s_%%0%dd" % len(str(len(combined)))
        for mapnumber, mapping in enumerate(combined, 1):
            mapping_name = formatter % (self.name, mapnumber)
            self.combined[mapping_name] = mapping

    def generate_assigner(self, motiffile=None, mode=None):
        if motiffile and mode == "load":
            self._limit_combinations(motiffile)
        self._generate_assigner(motiffile)
        if motiffile and mode == "dump":
            self.assigner.generate_motiffile(motiffile)

    def _generate_assigner(self, motiffile):
        graphs = {motif: dgmap.code for motif, dgmap in self.combined.items()}
        formulae = dict()
        assigner = self.augend.assigner
        group_no = max(zip(*assigner.definedin.values())[0])+1
        definedin = {motif: (group_no, motiffile) for motif in self.combined}
        for motif in {self.augend.graph, self.addend.graph, self.restrict}:
            if motif:
                definedin[motif] = assigner.definedin[motif]
                if assigner.type_of(motif) == "Graph":
                    graphs[motif] = assigner.graphs[motif]
                else:
                    formulae[motif] = assigner.formulae[motif]
        self.assigner = MotifAssigner(definedin, graphs, formulae, "all")

    def _limit_combinations(self, loadfile):
        loaded = MotifAssigner.extract(loadfile)
        augcode = self.augend.code
        combined = dict()
        for loadmotif, code in loaded.graphs.items():
            code = RNAStructureGraph.from_code(code).mindfscode(augcode)
            for motif, mapping in self.combined.items():
                if mapping.code == code:
                    combined[loadmotif] = mapping
                    break
            else:
                raise ValueError("Couldn't find %s" % loadmotif)
        self.combined = combined

    def generate_assigns(self, assignsfiles, outputfile):
        if not self.assigner:
            self.generate_assigner()
        self.assigns = MotifAssignments(self.assigner)
        self.augend._set_index(self.assigns)
        self.addend._set_index(self.assigns)
        self._set_index()
        self.add_embeddings(assignsfiles)
        self.assigns.dump(outputfile, detail=True)

    def add_embeddings(self, assignsfiles):
        load = MotifAssignments.load
        limit = {self.augend.graph, self.addend.graph, self.restrict} - {None}
        loaded = None
        assigner = self.augend.assigner
        filterkey = lambda case: self.restrict in loaded.by_case[case]
        for started, assignsfile in enumerate(assignsfiles):
            loaded = load(assigner, assignsfile, True, limit)
            cases = set(loaded.by_graph[self.augend.graph])
            cases &= set(loaded.by_graph[self.addend.graph])
            if self.restrict:
                cases = filter(filterkey, cases)
            cases = sorted(cases)
            for casecount, case in enumerate(cases, 1):
                augbeds = loaded.by_graph[self.augend.graph][case]
                addbeds = loaded.by_graph[self.addend.graph][case]
                for augbed, addbed in product(augbeds, addbeds):
                    self._add_embedding(case, augbed, addbed)

    def _add_embedding(self, case, augbed, addbed):
        for graph, dgmap in self.combined.items():
            if self._check_embedding(dgmap, augbed, addbed):
                embed = dict()
                for dtagstr, totagstr in dgmap.extras.items():
                    toindex = self.index[graph][totagstr]
                    embed[toindex] = addbed[self.addend.index[dtagstr]]
                for totagstr, augex in self.augend.index.items():
                    embed[self.index[graph][totagstr]] = augbed[augex]
                embed = list(list(zip(*sorted(embed.items())))[1])
                self.assigns.by_case[case].add(graph)
                if case not in self.assigns.by_graph[graph]:
                    self.assigns.by_graph[graph][case] = list()
                self.assigns.by_graph[graph][case].append(embed)

    def _check_embedding(self, dgmap, augbed, addbed):
        common_nodes = dict()
        for dtagstr, gtagstr in dgmap.common.items():
            if addbed[self.addend.index[dtagstr]]\
              != augbed[self.augend.index[gtagstr]]:
                return False
            if "->" not in dtagstr:
                dnodetag = int(dtagstr.split(" ")[0])
                gnodetag = int(gtagstr.split(" ")[0])
                common_nodes[dnodetag] = gnodetag
        for dtagstr, totagstr in dgmap.extras.items():
            for gtagstr, index in self.augend.index.items():
                if "->" in dtagstr and "->" in gtagstr\
                  and dtagstr.split(" ")[0] == gtagstr.split(" ")[0]\
                  and addbed[self.addend.index[dtagstr]] == augbed[index]:
                    dnodetags = dtagstr.split(" ")[1].strip("()").split("->")
                    try:
                        gnodetag = common_nodes[int(dnodetags[0])]
                        gedontag = common_nodes[int(dnodetags[1])]
                    except KeyError:
                        continue
                    if gtagstr.endswith("(%d->%d)" % (gnodetag, gedontag))\
                      or gtagstr.endswith("(%d->%d)" % (gedontag, gnodetag)):
                        return False
                elif "->" not in dtagstr and "->" not in gtagstr:
                    if addbed[self.addend.index[dtagstr]] == augbed[index]:
                        return False
        return True

    def _set_index(self):
        self.index = dict()
        for graph in self.combined:
            head_iter = enumerate(self.assigns._graph_header_details(graph))
            self.index[graph] = {tagstr: index for index, tagstr in head_iter}
            for totagstr, augex in self.augend.index.items():
                if totagstr in self.index[graph]:
                    continue
                edgelab, nodetags = totagstr.split(" ")
                nodetags = tuple(map(int, nodetags.strip("()").split("->")))
                totagstr_right = "%s (%d->%d)" % ((edgelab,)+nodetags[::-1])
                if totagstr_right not in self.index[graph]:
                    raise ValueError("AugTagStr not found")
                self.index[graph][totagstr] = self.index[graph][totagstr_right]

    def __remove_base(self, code):
        graph = RNAStructureGraph.from_code(code)
        embed = graph.match_code(self.base)
        if embed:
            for edge in embed.edges:
                graph.edges.remove(edge)
            for nodeid in embed.nodeids:
                if not len(graph.incident_edges(nodeid)):
                    graph.delete_node(nodeid)
        return graph, embed

    def __get_mappable(self):
        auggf, augbed = self.__remove_base(self.augend.code)
        addgf, addbed = self.__remove_base(self.addend.code)
        basemap = self.AddToAugMap(dict(), dict(), self.base, None, None)
        addnodes = set(addgf.nodes)
        if addbed:
            addnodes -= set(addbed.nodeids)
            for addid, augid in zip(addbed.nodeids, augbed.nodeids):
                basemap.to[addid] = augid
                basemap.fro[augid] = addid
        augnodes = set(auggf.nodes) - set(augbed.nodeids)
        mappable = {addid: list([-1]) for addid in addnodes}
        for addid, augid in product(addnodes, augnodes):
            if addgf.nodes[addid].covers(auggf.nodes[augid]):
                mappable[addid].append(augid)
        return mappable, basemap, addgf.edges, auggf.edges

    def __generate_mappings(self, mappable, basemap, firstid):
        mappings = list()
        mapkeys, mapvalues = zip(*sorted(mappable.items()))
        for mapping in product(*mapvalues):
            if (len(mapping) - mapping.count(-1)) != len(set(mapping)-{-1}):
                continue
            mapping = dict(zip(mapkeys, mapping))
            mapping = self.AddToAugMap(mapping, dict(), None, dict(), dict())
            counter = count(firstid)
            for addid in mapkeys:
                if mapping.to[addid] == -1:
                    mapping.to[addid] = counter.next()
                else:
                    dnodelab = self.addend.code._tags[addid]
                    gnodelab = self.augend.code._tags[mapping.to[addid]]
                    dtagstr = "%d %s" % (addid, dnodelab)
                    gtagstr = "%d %s" % (mapping.to[addid], gnodelab)
                    mapping.common[dtagstr] = gtagstr
                mapping.fro[mapping.to[addid]] = addid
            mapping.to.update(basemap.to)
            mapping.fro.update(basemap.fro)
            mappings.append(mapping)
        return mappings

    def __extend_nodes(self, combinedgf, mapping, firstid):
        for toid in count(firstid):
            if toid in mapping.fro:
                addid = mapping.fro[toid]
                nodelab = self.addend.code._tags[addid]
                dtagstr = "%d %s" % (addid, nodelab.name)
                totagstr = "%d %s" % (toid, nodelab.name)
                mapping.extras[dtagstr] = totagstr
                combinedgf.add_node(toid, nodelab)
            else:
                break

    def __extend_edges(self, combinedgf, mapping, addedges, augedges):
        for dnodetag, dedontag, edgelab in addedges:
            dtagstr = "%s (%d->%d)" % (edgelab.name, dnodetag, dedontag)
            tonodetag, toedontag = mapping.to[dnodetag], mapping.to[dedontag]
            if (tonodetag, toedontag, edgelab) in augedges:
                gtagstr = "%s (%d->%d)" % (edgelab, tonodetag, toedontag)
                mapping.common[dtagstr] = gtagstr
            elif (toedontag, tonodetag, edgelab) in augedges:
                gtagstr = "%s (%d->%d)" % (edgelab, toedontag, tonodetag)
                mapping.common[dtagstr] = gtagstr
            else:
                combinedgf.add_edge(tonodetag, toedontag, edgelab)
                totagstr = "%s (%d->%d)" % (edgelab, tonodetag, toedontag)
                mapping.extras[dtagstr] = totagstr


class CaseReader(object):

    @classmethod
    def read_limited(cls, limited=None, restrict=None, division=None):
        if (not limited) or (division and restrict):
            raise ValueError("Values set for both limited and division")
        if division:
            _, cases = MotifDivideFile.read_divisions(limited, int(division))
            cases = set(list(zip(*cases))[0])
        else:
            cases = set()
            reader = MotifAssignments._assigns_reader(limited)
            for isheader, _, motif, line in reader:
                if not restrict or (motif == restrict):
                    cases.add(" ".join(line.split("\t")[:2]))
        return cases

