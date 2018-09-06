import os

from io import StringIO
from itertools import combinations, chain
from collections import defaultdict

from .engraph import RNAStructureGraph, DFSCode, DFSCodeFragment


class MotifGraphCode(DFSCode):

    def __init__(self, code):
        super(MotifGraphCode, self).__init__(code)
        self.symmetry = self._count_symmetry()

    def append(self, newfrag):
        super(MotifGraphCode, self).append(newfrag)
        self.symmetry = self._count_symmetry()

    def relabel(self, nlabs, elabs):
        relabelled = super(MotifGraphCode, self).relabel(nlabs, elabs)
        return MotifGraphCode(relabelled._code)

    def copy(self):
        return MotifGraphCode(super(MotifGraphCode, self).copy()._code)

    @classmethod
    def from_lines(cls, lines, domain, base=None):
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


class MotifAssignments(object):

    _HAS = {True: "", False: "no "}
    
    def __init__(self, assigner, limited=set()):
        self.assigner = assigner
        self.by_case = defaultdict(set)
        self.by_graph = defaultdict(dict)
        self.by_formula = defaultdict(dict)
        self.by_type = {"Graph": self.by_graph, "Formula": self.by_formula}
        self._motifs = {"Graph": set(self.assigner.graphs),
                        "Formula": set(self.assigner.formulae)}
        self.limit_to(limited)

    def type_of(self, motif):
        if motif in self.assigner.graphs:
            return "Graph"
        elif motif in self.assigner.formulae:
            return "Formula"
        raise ValueError("No motif named '%s'" % motif)

    def limit_to(self, limited):
        if limited:
            limited = set(limited)
            self._motifs["Graph"] = set(self.assigner.graphs) & limited
            self._motifs["Formula"] = set(self.assigner.formulae) & limited

    def motif_iterator(self):
        limited = self._motifs["Graph"]|self._motifs["Formula"]
        return self.assigner.motif_iterator(limited)

    def add_from(self, rnagfs):
        assigner = self.assigner
        considered = (self._motifs["Graph"] | self._motifs["Formula"])
        for rnagf in rnagfs:
            preassigns, prerejects = set(), set()
            if rnagf.name in self.by_case:
                preassigns = set(self.by_case[rnagf.name])
                prerejects = considered - preassigns
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
            if not isheader:
                line = [value.strip() for value in line.split("\t")]
                case, record = " ".join(line[:2]), line[2:]
                assignments.by_case[case].add(motif)
                if case not in assignments.by_type[mtype][motif]:
                    assignments.by_type[mtype][motif][case] = list()
                    if not detail:
                        assignments.by_type[mtype][motif][case].append(list())
                if detail:
                    if mtype == "Formula":
                        record = assignments._load_truth_details(record)
                        assignments.by_type[mtype][motif][case] += record
                    else:
                        assignments.by_type[mtype][motif][case].append(record)
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
    def merge_files(cls, assigner, assignfiles, outputfile):
        lines = defaultdict(list)
        for assignfile in assignfiles:
            for isheader, _, motif, line in cls._assigns_reader(assignfile):
                if (isheader and motif not in lines) or not isheader:
                    lines[motif].append(line)
        outstring = list()
        for _, motif in assigner.motif_iterator(limited=set(lines)):
            outstring += [""]+lines[motif]
        with open(outputfile, 'w') as assignfile:
            assignfile.write("\n".join(outstring))

    @classmethod
    def limit_file(cls, assignfile, limited=set(), outputfile=None):
        outstring = list()
        keep_motif = False
        if_headers = list()
        cases = limited.get("cases", set())
        motifs = limited.get("motifs", set())
        for isheader, _, motif, line in cls._assigns_reader(assignfile):
            if isheader:
                keep_motif = motif in motifs if motifs else True
                if if_headers and if_headers[-1]:
                    outstring.pop()
                    if_headers.pop()
                    keep_motif = False
                elif keep_motif:
                    if_headers.append(True)
                    outstring.append("\n"+line)
            elif keep_motif: # checking case
                if not cases or line.split("\t")[1] in cases:
                    if_headers.append(False)
                    outstring.append(line)
        if outputfile:
            with open(outputfile, 'w') as outfile:
                outfile.write("\n".join(outstring))
        else:
            return "\n".join(outstring)

    @staticmethod
    def _assigns_reader(assignfile):
        mtype, motif = None, None
        is_stream = isinstance(assignfile, StringIO)
        with (assignfile if is_stream else open(assignfile, 'r')) as stream:
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
        header = list()
        for frag in self.assigner.graphs[graph]:
            if frag.isnode:
                header.append("%d %s" % (frag.nodetag, frag.nodelab.name))
            else:
                edgeargs = (frag.edgelab.name, frag.nodetag, frag.edontag)
                header.append("%s (%d->%d)" % edgeargs)
                if frag.intree:
                    header.append("%d %s" % (frag.edontag, frag.edonlab.name))
        return header

    def _dump_formula(self, formula, detail=False):
        cases = sorted(self.by_formula[formula])
        records = [["Formula", formula]] + [case.split(" ") for case in cases]
        if detail:
            for clause in self.assigner.formulae[formula]:
                records[0].append("")
                records[0] += [self._HAS[has_]+term for has_, term in clause]
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
            print(code.to_string(header="Graph:\t%s" % name, indent="\t"))
        for formula in self.formula_list:
            print("Formula:\t%s" % formula)
            for clause in self.formulae[formula]:
                line = [self._HAS[has_] for has_, part in clause]
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
                granny1 = mom1 in ancestry[mom2] and mom1 in ancestry[graph]
                granny2 = mom2 in ancestry[mom1] and mom2 in ancestry[graph]
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
            code1len = len(code1)
            ancestors[name1].append(self._ROOT)
            for name2, code2 in graphs:
                if len(code2) > code1len:
                    break
                if name2 != name1 and graph1.embed_code(code2):
                    ancestors[name1].append(name2)
        return ancestors

    def _resolve_formula_order(self):
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
        return dependency

    @staticmethod
    def _by_complexity(mapping):
        sorter = lambda item: (len(item[1]), item[0])
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
        self.domain = RNAStructureGraph.LABELS

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

    def _extract_from(self, file):
        defbytype = {"Graph": list(), "Formula": list()}
        for mtype, motif, deflist in self._extract_definitions_from(file):
            if motif in self.definedin:
                definedin = self.definedin[motif]
                raise self.RepeatPatternError(motif, file, definedin)
            self.definedin[motif] = file
            defbytype[mtype].append((motif, deflist))
        self._parse_graphs(defbytype["Graph"])
        self._parse_formulae(defbytype["Formula"])
        self.files.append(file)

    @classmethod
    def _extract_definitions_from(cls, file):
        definitions = list()
        with open(file, 'r') as infile:
            for lineno, line in enumerate(infile, 1):
                try:
                    line = line.split("!")[0].strip() # ignore comments
                    if line.startswith("Graph:")\
                      or line.startswith("Formula:"):
                        mtype, motif = line.split(":")[:2]
                        defentry = (mtype, motif.strip(), list())
                        definitions.append(defentry)
                    elif line:
                        definitions[-1][2].append(line)
                except:
                    raise cls.FileFormatError(file, lineno)
        return definitions

    def _parse_graphs(self, graphdefs):
        for motif, refers, deflist in self._parse_graph_references(graphdefs):
            try:
                base = self.graphs[refers] if refers else None
                code = MotifGraphCode.from_lines(deflist, self.domain, base)
                self.graphs[motif] = code
            except KeyError as e:
                raise BadReferenceError(refers, motif)
            except DFSCodeFragment.InvalidFragmentError as e:
                raise self.BadDefinitionError(e.args[0], motif)
            except AssertionError as e:
                print(motif)
                raise e

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
