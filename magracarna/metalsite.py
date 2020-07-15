from __future__ import print_function

import os
import sys

import re
import numpy as np
import subprocess

from decimal import Decimal
from warnings import warn
from itertools import combinations
from collections import defaultdict

from .engrid import HydrogenBondFinder
from .nucleotides import phosphate, backbone, PO_DISTANCE
from .engraph import RNAStructureGraph, RNAStructureEngrapher, RNAGraphList
from .motifs import MotifAssignments, MotifAssigner, MotifGraphCode,\
                    MotifDivider, CaseReader
from .vectorise import MotifVectorList


PHOSPHATE_GROUP = phosphate[:3]+["O5'", "O3'"]
DONORS = HydrogenBondFinder.DONORS


class RNAMetalSiteGraph(RNAStructureGraph):
    
    INNER_SHELL_NO = 1
    OUTER_SHELL_NO = 2
    SHELL_NAME = {INNER_SHELL_NO: RNAStructureGraph.LIGATION,
                  OUTER_SHELL_NO: RNAStructureGraph.WATERMED}
    ORIENTATION = {"cis": RNAStructureGraph.CIS,
                   "trans": RNAStructureGraph.TRANS}

    def __init__(self, metnode):
        super(RNAMetalSiteGraph, self).__init__(metnode[0])
        metid, metlab = self.name_residue_node(metnode)
        self.metid = metid
        self.metlab = metlab
        self.add_node(metid, metlab)
        self.nucleobases = set()

    def add_ligand_interaction(self, resnode, atom, shell_no):
        if shell_no not in self.SHELL_NAME:
            raise ValueError("Invalid shell")
        nucid, nuclab = self.name_residue_node(resnode)
        atomid, atomlab = self.name_residue_part_node(resnode, atom)
        self.connect_part_to_residue(resnode, atom)
        edgelabel = self.SHELL_NAME[shell_no]
        self.add_edge(self.metid, atomid, edgelabel)
        if nuclab in self.NUCLEOTIDE_RESIDUES\
          or nuclab in (self.ANY_NUCLEOTIDE, self.PURINE, self.PYRIMIDINE):
            if atom in self.NUCLEOBASE_ATOMS:
                self.nucleobases.add(nucid)

    def add_inter_ligand_orientation(self, lig1, lig2, atom1, atom2, orient):
        atomid1, atomlab1 = self.name_residue_part_node(lig1, atom1)
        atomid2, atomlab2 = self.name_residue_part_node(lig2, atom2)
        self.add_node(atomid1, atomlab1)
        self.add_node(atomid2, atomlab2)
        self.add_edge(atomid1, atomid2, self.ORIENTATION[orient])

    def record_metal_basepair_interactions(self):
        for nucid1, nucid2, bpid in self.pairs:
            if nucid1 in self.nucleobases or nucid2 in self.nucleobases:
                self.add_edge(bpid, self.metid, self.BASE_PAIR_CONTACT)

    def derive_structure(self):
        super(RNAMetalSiteGraph, self).derive_structure()
        self.record_metal_basepair_interactions()


class MgSiteValidator(object):
    
    _VERBOSE = False
    MIN_COORDINATION = 4
    MAX_COORDINATION = 6
    _Q_E_CUTOFF = Decimal("0.5")
    _Q_S_CUTOFF = Decimal("0.6")
    _Q_V_CUTOFF = Decimal("0.5")

    # Bond valence values are from two different sources
    # For O, N, F, S, Cl and Br, it is from 
    #   "Data Mining..." by Zheng, et al.
    # For P, I, Se, Te and As, it is from
    #   "Bond Valence..." by Brese, et al.
    #

    VALENCE_PARAMS = {"O": Decimal("1.67"), "N": Decimal("1.78"),
                      "F": Decimal("1.64"), "S": Decimal("2.20"),
                      "CL": Decimal("2.11"), "BR": Decimal("2.26"),
                      "P": Decimal("2.290"), "I": Decimal("2.460"),
                      "SE": Decimal("2.32"), "TE": Decimal("2.53"),
                      "AS": Decimal("2.38")}

    @classmethod
    def check(cls, site):
        coordination = len(site.inner_shell)
        valid_coordination = cls.MIN_COORDINATION <= coordination
        valid_coordination &= coordination <= cls.MAX_COORDINATION
        valencies, net_valence = cls.valencies(site)
        Q_e_pass, Q_e = cls.check_environment(site, valencies, net_valence)
        Q_s_pass, Q_s = cls.check_symmetry(site, valencies, net_valence)
        Q_v_pass, Q_v = cls.check_valence(site, valencies, net_valence)
        passes = valid_coordination & Q_e_pass & Q_v_pass & Q_s_pass
        if cls._VERBOSE:
            details = (site.graph.name, coordination, Q_v, Q_s, Q_e)
            details = details + ("" if passes else site.INVALIDCHAR,)
            print("%s\t%d\t%s\t%s\t%s\t%s" % details)
        return passes

    @classmethod
    def valencies(cls, site):
        valencies = dict()
        net_valence = 0.0
        for entry in site.environment:
            _, atom, distance = entry
            valencies[entry] = cls.calculate_valence(atom, distance)
        for entry in site.inner_shell:
            net_valence += valencies[entry]
        return valencies, net_valence

    @classmethod
    def check_valence(cls, site, valencies, net_valence):
        Q_v_value = 1 - abs(net_valence - 2) / 2.0
        Q_v_value = Decimal("%.3f" % Q_v_value)
        if Q_v_value == Decimal("-0.00"):
            Q_v_value = Decimal("0.00")
        passes = Q_v_value >= cls._Q_V_CUTOFF
        return passes, Q_v_value

    @classmethod
    def check_symmetry(cls, site, valencies, net_valence):
        try:
            in_sh_vector = np.zeros((3,))
            for entry in site.inner_shell:
                valence = valencies[entry]
                _, atom, distance = entry
                vector = site._atom.vector_to(atom)
                in_sh_vector += (atom.occy * valence * vector) / distance
            Q_s_value = 1 - (np.linalg.norm(in_sh_vector) / net_valence)
            Q_s_value = Decimal("%.3f" % Q_s_value)
            if Q_s_value == Decimal("-0.00"):
                Q_s_value = Decimal("0.00")
        except (TypeError, ZeroDivisionError):
            Q_s_value = Decimal("0")
        passes = Q_s_value >= cls._Q_S_CUTOFF
        return passes, Q_s_value

    @classmethod
    def check_environment(cls, site, valencies, net_valence):
        env_occy, env_bfac = cls.environment_values(site, valencies)
        try:
            if env_occy == site._atom.occy == 1.0:
                min_bfac, max_bfac  = sorted((site._atom.bfac, env_bfac))
                Q_e_value = float(min_bfac) / float(max_bfac)
            else:
                env_ratio = site._atom.bfac / site._atom.occy
                env_ratio /= float(env_bfac) / float(env_occy)
                env_ratio = min(env_ratio, 1/env_ratio)
                Q_e_value = min(env_occy, site._atom.occy) * env_ratio
            Q_e_value = Decimal("%.3f" % Q_e_value)
            if Q_e_value == Decimal("-0.00"):
                Q_e_value = Decimal("0.00")
        except (TypeError, ZeroDivisionError):
            Q_e_value = Decimal("0")
        passes = Q_e_value >= cls._Q_E_CUTOFF
        return passes, Q_e_value

    @classmethod
    def environment_values(cls, site, valencies):
        env_valence = 0.0
        env_bfactor = 0.0
        env_occupancy = 0.0
        for entry in site.environment:
            valence = valencies[entry]
            _, atom, _, = entry
            env_valence += valence
            env_bfactor += atom.bfac * valence
            env_occupancy += atom.occy * valence
        return env_occupancy/env_valence, env_bfactor/env_valence

    @classmethod
    def calculate_valence(cls, atom, distance):
        valence = (float(cls.VALENCE_PARAMS[atom.atype]) - distance) / 0.37
        valence = np.power(np.e, valence)
        return valence


class RNAMgSiteEngrapher(RNAStructureEngrapher):
        
    _TYPE = "SiteGraph"
    INVALIDCHAR = "!"
    _RELAX_BY = Decimal("0.5")
    _RELAX_ANGLE = Decimal("50.0")
    HBOND_ELEMENTS = ("O", "N")
    WATER_NAMES = RNAMetalSiteGraph.WATER_NAMES
    ENVIRONMENT = MgSiteValidator.VALENCE_PARAMS.keys()

    MAX_COORDINATION = MgSiteValidator.MAX_COORDINATION

    ENV_BOUND = 4.0
    CIS_LIMIT = 135.0
    LIGDISTS = {"O": Decimal("2.58"), "N": Decimal("2.70")}

    class NoRNAInteraction(ValueError):
        MESG = "%s in structure %s does not interact with RNA"
        def __init__(self, siteid):
            self.structid, self.residueid = siteid.split(" ")
            message = self.MESG % (self.residueid, self.structid)
            super(RNAMgSiteEngrapher.NoRNAInteraction, self).__init__(message)
    
    class FailsValidation(ValueError):
        MESG = "%s in structure %s does not pass validation criteria"
        def __init__(self, siteid):
            self.structid, self.residueid = siteid.split(" ")
            message = self.MESG % (self.residueid, self.structid)
            super(RNAMgSiteEngrapher.FailsValidation, self).__init__(message)

    def __init__(self, hbfinder, depth=3, hbonds=False, ligdists=None,
                 relaxed=False, invalids=True, modbases=True, verbose=False):
        super(RNAMgSiteEngrapher, self).__init__(hbfinder)
        self.depth = depth
        self._hbonds = hbonds
        self._ligdists = ligdists if ligdists else self.LIGDISTS.copy()
        self._relaxed = relaxed
        self._invalids = invalids
        self._verbose = verbose
        if not modbases:
            self.NUCLEOTIDE_RESNAMES = ("A", "C", "G", "U")
        self.set_structure()

    @classmethod
    def limit_file(cls, filepath, limited=dict(), outputfile=None):
        cases = limited.get("cases", set())
        invalids = limited.get("invalids", True)
        to_keep = False
        outstring = list()
        with open(filepath, 'r') as infile:
            for line in infile:
                line = line.strip()
                if line:
                    isheader = line.startswith(cls._TYPE)
                    if isheader:
                        case = line.split("\t")[2]
                        isvalid = not case.endswith(cls.INVALIDCHAR)
                        incases = case.strip("!") in cases if cases else True
                        to_keep = (invalids or isvalid) and incases
                    if to_keep:
                        outstring.append(("\n"+line if isheader else line))
        if outputfile:
            with open(outputfile, 'w') as outfile:
                outfile.write("\n".join(outstring))
        else:
            return "\n".join(outstring)

    def set_structure(self, structure=None):
        super(RNAMgSiteEngrapher, self).set_structure(structure)
        self._node = None
        self._atom = None
        self.residue = None
        self.environment = set()
        self.inner_shell = set()
        self.outer_shell = set()

    def set_residue(self, residue, structure):
        self.set_structure(structure)
        self.residue = residue
        self._atom = residue.atoms["MG"]
        self._node = (str(residue), "MG")
        self.graph = RNAMetalSiteGraph(self._node)
        self.graph.name = "%s %s" % (structure.structid, self.graph.name)

    def detect_mgsite(self, residue, engraph):
        self.set_residue(residue, self.structure)
        self.detect_environment()
        self.detect_inner_shell(engraph=engraph)
        if self._relaxed and len(self.inner_shell) < self.MAX_COORDINATION:
            self._relaxed_inner_shell()
        self.detect_outer_shell(engraph=engraph)
        if not engraph:
            for entry in self.inner_shell | self.outer_shell:
                residue = entry[0]
                if residue.name in self.NUCLEOTIDE_RESNAMES:
                    self.residues_added.add(residue)
        if not self.residues_added:
            raise self.NoRNAInteraction(self.graph.name)
        is_validated = MgSiteValidator.check(self)
        if not self._invalids and not is_validated:
            raise self.FailsValidation(self.graph.name)
        if not is_validated:
            self.graph.name = self.graph.name+self.INVALIDCHAR

    def engraph_mgsite_context(self, residue):
        self.detect_mgsite(residue, engraph=True)
        self.engraph_sequence_context()
        self.engraph_innershell_orientations()

    def _relaxed_inner_shell(self, engraph=True):
        strict = self._ligdists
        relaxed = {atom: dist+self._RELAX_BY for atom, dist in strict.items()}
        environment = self.environment
        old_inner_shell = sorted(self.inner_shell, key=lambda entry: entry[2])
        self.set_residue(self.residue, self.structure)
        self.environment = environment
        self._ligdists = relaxed
        self.detect_inner_shell(engraph=False)
        self._ligdists = strict
        self.inner_shell = self._validate_relaxed_inner_shell(old_inner_shell)
        if engraph:
            self.engraph_shell(self.inner_shell, self.graph.INNER_SHELL_NO)

    def _validate_relaxed_inner_shell(self, old_inner_shell):
        """ Validate new inner-shell ligands and add to old ligands.

        Args:
            old_inner_shell: inner-shell entries from the stricter
                inner-shell detection step.
        """
        phosphates = set() # set of all phosphate groups in inner-shell
        inner_shell = old_inner_shell[:]
        angle_between = self._atom.angle_between
        new_inner_shell = sorted(self.inner_shell, key=lambda entry: entry[2])
        for entry in new_inner_shell:
            residue, atom, distance = entry
            if self._include_in_relaxed_inner_shell(entry, phosphates):
                for _, oldatom, _ in inner_shell:
                    if oldatom is atom:
                        break
                    angle = Decimal("%.2f" % angle_between(atom, oldatom))
                    if angle < self._RELAX_ANGLE:
                        break
                else:
                    inner_shell.append(entry)
        return set(inner_shell[:self.MAX_COORDINATION])

    def _include_in_relaxed_inner_shell(self, entry, phosphates):
        residue, atom, distance = entry
        include = distance > self._ligdists[atom.atype]
        if include and atom.atype == "N" and residue.name in DONORS:
            include = atom.name not in DONORS[residue.name]
        elif atom.name in PHOSPHATE_GROUP\
          and residue.name in self.NUCLEOTIDE_RESNAMES:
            rindex = residue.index
            begex, endex, _ = self.structure.chains[residue.chain]
            if atom.name == "O3'":
                get_nearby_atoms = self.structure.grid.get_nearby_atoms
                for entry_ in get_nearby_atoms(PO_DISTANCE, atom):
                    rindex_, atom_ = entry_[0].index, entry_[1]
                    if atom_.name == "P" and begex <= rindex_ < endex:
                        rindex = rindex_
                        break
            if rindex in phosphates:
                include = False
            else:
                phosphates.add(rindex)
        return include

    def detect_environment(self):
        get_nearby_atoms = self.structure.grid.get_nearby_atoms
        envmt = get_nearby_atoms(self.ENV_BOUND, self._atom)
        for residue, atom, distance in envmt:
            if atom.atype in self.ENVIRONMENT:
                self.environment.add((residue, atom, distance))

    def detect_inner_shell(self, engraph=True):
        for entry in self.environment:
            _, atom, distance = entry
            if atom.atype in self._ligdists\
              and distance <= self._ligdists[atom.atype]:
                self.inner_shell.add(entry)
        if engraph:
            self.engraph_shell(self.inner_shell, self.graph.INNER_SHELL_NO)

    def detect_outer_shell(self, engraph=True):
        add_hbond = self.graph.add_hbond
        find_hbonds = self.hbfinder.find_hbonds
        inner_shell = set(zip(*list(zip(*self.inner_shell))[:2]))
        for residue, atom in inner_shell:
            if residue.name in self.WATER_NAMES\
              and atom.atype in self.HBOND_ELEMENTS:
                for entry in find_hbonds(atom, residue):
                    bond_residue, bond_atom = entry
                    if entry not in inner_shell\
                      and bond_atom.atype in self.HBOND_ELEMENTS:
                        self.outer_shell.add(entry)
                        if self._hbonds and engraph:
                            res1 = (str(residue), residue.name)
                            res2 = (str(bond_residue), bond_residue.name)
                            add_hbond(res1, atom.name, res2, bond_atom.name)
        if engraph:
            self.engraph_shell(self.outer_shell, self.graph.OUTER_SHELL_NO)

    def engraph_shell(self, shell, shell_no):
        add_ligand_interaction = self.graph.add_ligand_interaction
        for entry in shell:
            residue, atom = entry[:2]
            if self._hbonds or residue.name not in self.WATER_NAMES:
                resnode = (str(residue), residue.name)
                add_ligand_interaction(resnode, atom.name, shell_no)
                if residue.name in self.NUCLEOTIDE_RESNAMES:
                    self.residues_added.add(residue)

    def engraph_innershell_orientations(self):
        angle_between = self._atom.angle_between
        for entry1, entry2 in combinations(self.inner_shell, 2):
            (res1, atom1), (res2, atom2) = entry1[:2], entry2[:2]
            if self._hbonds or not {res1.name, res2.name} & self.WATER_NAMES:
                res1, res2 = (str(res1), res1.name), (str(res2), res2.name)
                angle = Decimal("%.2f" % angle_between(atom1, atom2))
                orient = "cis" if angle <= self.CIS_LIMIT else "trans"
                args = (res1, res2, atom1.name, atom2.name, orient)
                self.graph.add_inter_ligand_orientation(*args)

    def engraph_sequence_context(self):
        chains = self.structure.chains
        rindices = set()
        for residue in self.residues_added:
            if residue.chain in chains:
                begex, endex, _ = chains[residue.chain]
                rindex = residue.index
                if begex <= rindex < endex:
                    minex = max(rindex-self.depth, begex)
                    maxex = min(rindex+self.depth, endex-1)
                    rindices.update(range(minex, maxex+1))
        residue_set = {self.structure.residues[rindex] for rindex in rindices}
        self.engraph_structure(residue_set, True)

    def _id_sorter(self, nodeid):
        sorter = super(RNAMgSiteEngrapher, self)._id_sorter(nodeid)
        return (nodeid != self._node[0],) + sorter


class RNAMgSiteList(RNAGraphList):

    _TYPE = RNAMgSiteEngrapher._TYPE

    def add(self, structure):
        """ Add Mg-RNA sites in given structure to itself.

        Args:
            structure: engrid.Structure object of the RNA
        """
        for residue in structure.residues[structure.chainend:]:
            if residue.name == "MG":
                try:
                    self.engrapher.set_structure(structure)
                    self.engrapher.engraph_mgsite_context(residue)
                except RNAMgSiteEngrapher.NoRNAInteraction:
                    pass
                except RNAMgSiteEngrapher.FailsValidation:
                    pass
                else:
                    self.update()

    @staticmethod
    def sorter(encoded_graph):
        header = encoded_graph.split("\n")[0].split("\t")
        sortkey = [header[1]]
        resid, chain = header[2].split(":")
        resno, insco = resid.replace("[MG]", ""), ""
        if "^" in resid:
            resno, insco = resid.split("^")
        sortkey += [int(resno), insco, chain]
        return tuple(sortkey)

    @classmethod
    def limited(cls, sitesfiles, limited, dumpfile=None):
        limited = sorted(limited)
        rnagfs = list()
        for (count, sitesfile) in enumerate(sitesfiles, 1):
            rnagfs += cls.load(sitesfile, limited)
            print("% 3d / %d" % (count, len(sitesfiles)))
        if not dumpfile:
            return rnagfs
        engrapher = RNAMgSiteEngrapher(None, None)
        siteslist = cls(dumpfile, engrapher, encoded=True)
        for rnagf in rnagfs:
            engrapher._node = (rnagf.name.split()[1], "MG")
            engrapher.graph = rnagf
            siteslist.update()
        siteslist.dump()


class RNAMgSiteContextParser(object):

    SITE = ["\t0\tMG", "\t1\t*\t0\tLW*"]
    SITE = MotifGraphCode.from_lines(SITE, RNAStructureGraph.LABELS, None)
    SITE_ASSIGNER = MotifAssigner({"SITE": ""}, {"SITE": SITE}, {}, "all")
    _CHAINS_HEADER = "Chains"

    @classmethod
    def site_selectors(cls, rnagfs):
        sites = dict()
        assignments = MotifAssignments(cls.SITE_ASSIGNER)
        assignments.add_from(rnagfs)
        for rnagf in rnagfs:
            residues = set()
            assigned_to = assignments.by_case[rnagf.name]
            for graph in assigned_to: # only graphs in that file
                for embed in assignments.by_graph[graph][rnagf.name]:
                    residues.update(cls.graph_embed_residues(embed))
            case, _ = cls.parse_site(rnagf.name)
            residues = [resdu.replace(":", ":\"")+"\"" for resdu in residues]
            sites[case] = " or ".join(residues)
        return sites

    @classmethod
    def context_selectors(cls, rnagfs):
        contexts = dict()
        for rnagf in rnagfs:
            case, _ = cls.parse_site(rnagf.name)
            residues = cls.graph_embed_residues(rnagf.nodes)
            residues = [resdu.replace(":", ":\"")+"\"" for resdu in residues]
            contexts[case] = " or ".join(residues)
        return contexts

    @classmethod
    def create_chainsfile(cls, rnagfs, chainsfile):
        sites = cls.site_selectors(rnagfs)
        sites = {case: value.split(" or ") for case, value in sites.items()}
        chain_wise = defaultdict(set)
        for case, identifiers in sites.items():
            for chain in cls.parse_chains(identifiers):
                chain_wise[chain].add(case)
        outstring = list()
        for chain, cases in chain_wise.items():
            cases = sorted(cases, key=lambda case: int(case.split(":")[0]))
            outstring.append("%s:\t%s" % (cls._CHAINS_HEADER, chain))
            outstring += ["\t"+case for case in cases] + [""]
        with open(chainsfile, 'w') as outfile:
            outfile.write("\n".join(outstring))

    @classmethod
    def read_chainsfile(cls, chainsfile):
        header = cls._CHAINS_HEADER+":"
        chain_wise = defaultdict(list)
        with open(chainsfile, 'r') as infile:
            for line in infile:
                line = [part.strip() for part in line.split("\t")]
                if len(line) > 1:
                    if line[0] == header:
                        chain = line[1]
                    else:
                        chain_wise[chain].append(line[1])
        return sorted(chain_wise.items())

    @staticmethod
    def graph_embed_residues(identifiers):
        residues = set()
        for identifier in identifiers:
            if identifier.startswith("["): # ensure only node identifiers
                for nodeid_part in identifier.split("="):
                    residues.add(nodeid_part.split(".")[0])
        return residues

    @staticmethod
    def parse_site(case):
        case = case.split(" ")[-1]
        case = case.replace(re.findall(r'\[.*\]', case)[0], "")
        case, isvalid = case.strip("!"), not case.endswith("!")
        return (case, isvalid)

    @staticmethod
    def parse_chains(identifiers):
        chains = set()
        for identifier in identifiers:
            if identifier.startswith("["):
                chains.add(identifier.split(".")[0].split(":")[1])
        return chains

    @staticmethod
    def parse_category(category_file):
        category = os.path.basename(category_file).split(".")[0].split("_")
        try:
            order_by, category = int(category[0]), " ".join(category[1:])
        except ValueError:
            order_by, category = 0, " ".join(category)
        return order_by, category


class RNAMgMotifCounts(RNAMgSiteContextParser):

    _SORTERS = ("category", "frequency", "dictionary")
    _MODIFIERS = (str.rjust, str.ljust, str.rjust, str.rjust)
    _ROW_LEGEND = ("category", "motif", "cindex", "mindex", "all", "valids")
    _HEADER_ROW = "Sites Present"
    
    def __init__(self, table):
        self.table = table

    @classmethod
    def create_from_assigns(cls, assigns):
        valids = sum(int(cls.parse_site(case)[1]) for case in assigns.by_case)
        table = [["", cls._HEADER_ROW, 0, 0, len(assigns.by_case), valids]]
        for mindex, (mtype, motif) in enumerate(assigns.motif_iterator(), 1):
            cindex, category = assigns.assigner.definedin[motif]
            _, category = cls.parse_category(category)
            allcases = len(assigns.by_type[mtype][motif])
            table.append([category, motif, cindex, mindex, allcases, 0])
            for case in assigns.by_type[mtype][motif]:
                _, isvalid = cls.parse_site(case)
                table[-1][-1] += int(isvalid)
        return cls(table)

    def dump(self, countsfile, sort_by="category"):
        table = sorted(self.table, key=self._row_sorter(sort_by))
        table = [map(str, entry[:2]+entry[4:]) for entry in table if entry[4]]
        maxlens = [max(map(len, column)) for column in zip(*table)]
        outable = list()
        for entry in table:
            outable.append(list())
            entry[1] = entry[1].replace("_", " ")
            for field, maxlen, modify in zip(entry, maxlens, self._MODIFIERS):
                outable[-1].append(modify(field, maxlen))
        with open(countsfile, 'w') as outfile:
            outstr = "\n".join(["\t".join(outrow) for outrow in outable])
            outfile.write(outstr)

    @classmethod
    def _row_sorter(cls, sort_by):
        if sort_by == "category":
            return lambda entry: (entry[2], entry[3], -entry[5], -entry[4])
        if sort_by == "frequency":
            return lambda entry: (-entry[5], -entry[4], entry[2], entry[3])
        if sort_by == "dictionary":
            return lambda entry: (entry[2], entry[1], -entry[5], entry[4])

    @classmethod
    def load(cls, assigner, countsfile):
        counts = cls.create_from_assigns(MotifAssignments(assigner))
        motifs = {entry[1]: entry for entry in counts.table}
        with open(countsfile, 'r') as infile:
            for index, line in enumerate(infile):
                line = [field.strip() for field in line.split("\t")]
                if not index and line[1] != "Sites Present":
                    message = '"Sites Present" row is absent'
                    message = "Cannot load '%s': %s" % (countsfile, message)
                    raise ValueError(message)
                if line:
                    line[1] = line[1].replace(" ", "_") if index else line[1]
                    if line[1] in motifs and line[0] == motifs[line[1]][0]:
                        motifs[line[1]][4] = int(line[2])
                        motifs[line[1]][5] = int(line[3])
        return counts

    @classmethod
    def merge(cls, assigner, countsfiles, mergedfile):
        merged = cls.create_from_assigns(MotifAssignments(assigner))
        motifs = {entry[1]: entry for entry in merged.table}
        for countsfile in countsfiles:
            counts = cls.load(assigner, countsfile)
            for entry in counts.table:
                motifs[entry[1]][4] += entry[4]
                motifs[entry[1]][5] += entry[5]
        merged.dump(mergedfile)


class RNAMgSiteContextVisualiser(object):

    def __init__(self, structure, rnagfs):
        self.sites = RNAMgSiteContextParser.site_selectors(rnagfs)
        self.rnagfs = rnagfs
        self.contexts = RNAMgSiteContextParser.context_selectors(rnagfs)
        self.structure = structure

    def _display(self, script):
        echo = subprocess.Popen(["echo", script], stdout=subprocess.PIPE)
        command = ["jmol", self.structure, "-s", "-"]
        with open(os.devnull, 'w') as NULL:
            subprocess.Popen(command, stdin=echo.stdout, stdout=NULL)

    def interactive(self):
        while True:
            try:
                case = raw_input()
                if case in self.sites:
                    self.display(case)
                else:
                    print("'%s' does not match any site identifier" % case)
            except KeyboardInterrupt:
                break

    def display(self, site, context=True, waters=False, structure=False):
        site = site.split("]")[-1].strip("!")
        mg_selector = site.replace(":", ":\"")+"\""
        site_selector = self.sites[site]
        context_selector = self.contexts[site]
        script = ["color background white",
                  "select not (metal or HOH)",
                  "trace only", "trace "+("0.01" if structure else "off"),
                  "select (metal or HOH)", "trace only", "trace off",
                  #
                  "select ("+context_selector+") and not (metal or HOH)",
                  "trace off" if context else "",
                  "wireframe "+("0.01" if context else "off"),
                  "spacefill "+("0.05" if context else "off"),
                  "select ("+context_selector+") and (metal or HOH)",
                  "spacefill 0.1",
                  #
                  "select (" +site_selector+ ")",
                  "trace off", "wireframe 0.15", "spacefill 0.25",
                  #
                  "zoomto ("+site_selector+") 0", "center "+mg_selector,
                  #
                  # "select ("+context_selector+") and nucleic and within(10.0, MG and clickable)", "trace 0.1", "color trace red",
                  # "select (" +site_selector+ ") and (*.P or *.O5' or *.O3' or *.OP1 or *.OP2 or not nucleic)",
                  # "wireframe 0.15", "spacefill 0.25",
                  # "select HOH and not within(3.1, MG and clickable) and protein; spacefill 0 ; wireframe 0 ; select none",
                  # 
                  # # bone
                  # "select nucleic and clickable ; color cartoon purple ; cartoon 0.1", "spacefill 0 ; wireframe 0",
                  # "select ("+context_selector+") and not (metal or HOH)", "spacefill 0", "wireframe 0",
                  # "select ("+site_selector+") and nucleic", "spacefill 0.25", "wireframe 0.15", 
                  # "select protein and clickable; spacefill 0 ; wireframe 0", color halo chain ; halo 0.5"
                  
                  # # site
                  # "select nucleic and not ("+site_selector+")", "wireframe 0 ; spacefill 0",
                  # # pair
                  # "select nucleic and clickable ; wireframe 0.15 ; spacefill 0.25"
                  # 
                  "select HOH and within(3.1, MG and clickable)",
                  "select selected "+("or HOH and within(3.8, selected)" if waters else ""),
                  "select HOH and not selected",
                  "spacefill 0" if not (waters and context) else "spacefill 0.25",
                  "select none",
                  ]
        self._display(";".join(script))


class RNAMgMotifDivider(MotifDivider):

    def __init__(self, motifgraph, assigner,
                 restrict=None, limited=None, invalids=True):
        superclass = super(RNAMgMotifDivider, self)
        superclass.__init__(motifgraph, assigner, restrict, limited)
        self.invalids = invalids

    def _add_embedding(self, case, embedding, keepembeds):
        if self.invalids or not case.endswith(RNAMgSiteEngrapher.INVALIDCHAR):
            args = (case, embedding, keepembeds)
            super(RNAMgMotifDivider, self)._add_embedding(*args)

    def add_sequences(self, sitesfiles):
        args = (sitesfiles, RNAMgSiteEngrapher)
        super(RNAMgMotifDivider, self).add_sequences(*args)


class RNAMgMotifVectorList(MotifVectorList):

    LIGATION = RNAStructureGraph.LIGATION
    WATERMED = RNAStructureGraph.WATERMED
    LIGLABELS = (LIGATION, WATERMED)
    
    def __init__(self, motifgraph, assigner, invalids=True, **kwargs):
        superclass = super(RNAMgMotifVectorList, self)
        superclass.__init__(motifgraph, assigner, **kwargs)
        self.invalids = invalids

    def _load_embeds(self, assignsfile):
        embeds = super(RNAMgMotifVectorList, self)._load_embeds(assignsfile)
        if not self.invalids:
            embeds, old_embeds = list(), embeds
            for case, embed in old_embeds:
                if not case.endswith(RNAMgSiteEngrapher.INVALIDCHAR):
                    embeds.append((case, embed))
        return embeds

    @classmethod
    def _embed_sorter(cls, embed):
        super_key = super(RNAMgMotifVectorList, cls)._embed_sorter(embed)
        embed_key = [part for part in embed if part in cls.LIGLABELS]
        embed_key = [int(part == cls.WATERMED) for part in embed_key]
        return (embed_key, super_key)

