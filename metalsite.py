
from decimal import Decimal
from warnings import warn
from itertools import combinations

import numpy as np

from .engraph import RNAStructureGraph, RNAStructureEngrapher


class RNAMetalContextGraph(RNAStructureGraph):
    
    INNER_SHELL_NO = 1
    OUTER_SHELL_NO = 2
    SHELL_NAME = {INNER_SHELL_NO: RNAStructureGraph.LIGATION,
                  OUTER_SHELL_NO: RNAStructureGraph.WATERMED}
    ORIENTATION = {"cis": RNAStructureGraph.INCISGEO,
                   "trans": RNAStructureGraph.TRANSGEO}

    def __init__(self, metnode, plinked=True):
        super(RNAMetalContextGraph, self).__init__(metnode[0], plinked)
        metid, metlab = self.name_residue_node(metnode)
        self.metid = metid
        self.metlab = metlab
        self.add_node(metid, metlab)

    def add_ligand_interaction(self, resnode, atom, shell_no):
        if shell_no not in self.SHELL_NAME:
            raise ValueError("Invalid shell number")
        nucid, nuclab = self.name_residue_node(resnode)
        atmid, atmlab = self.name_residue_part_node(resnode, atom)
        self.connect_part_to_residue(resnode, atom)
        edgelabel = self.SHELL_NAME[shell_no]
        self.add_edge(self.metid, atmid, edgelabel)
        if nuclab.name in self._NUCNAMES:
            for nucedge, edgeatoms in self._NUCEDGEATOMSMAP.items():
                if atom in edgeatoms:
                    self.add_edge(nucid, self.metid, nucedge)

    def add_inter_ligand_orientation(self, lig1, lig2, atom1, atom2, orient):
        atmid1, atmlab1 = self.name_residue_part_node(lig1, atom1)
        atmid2, atmlab2 = self.name_residue_part_node(lig2, atom2)
        self.add_node(atmid1, atmlab1)
        self.add_node(atmid2, atmlab2)
        edgelabel = self.ORIENTATION[orient]
        self.add_edge(atmid1, atmid2, edgelabel)

    def record_metal_basepair_interactions(self):
        for nucid1, nucid2, bpid in self.pairs:
            self.add_edge(bpid, self.metid, self.TOABPAIR)

    def derive_structure(self):
        super(RNAMetalContextGraph, self).derive_structure()
        self.record_metal_basepair_interactions()


class RNAMgSiteEngrapher(RNAStructureEngrapher):

    _TYPE = "ContextGraph"
    INVALIDCHAR = "!"
    HBOND_ELEMENTS = ("O", "N")
    WATER_RESIDUE_NAMES = ("HOH", "O")

    ENV_BOUND = 4.0
    CIS_LIMIT = 135.0
    LIGDISTS = {"O": Decimal("2.58"), "N": Decimal("2.70")}
    VALCEPRMS = {"O": 1.693, "N": 1.850, "P": 2.290, "CL": 2.08, "F": 1.581,
                 "BR": 2.28, "I": 2.460, "S": 2.180, "SE": 2.32, "TE": 2.53,
                 "AS": 2.38}

    class NoRNAInteraction(ValueError):
        MESG = "%s in structure %s does not interact with RNA"
        def __init__(self, residueid, structid):
            self.structid = structid
            self.residueid = residueid
            message = self.MESG % (residueid, structid)
            super(RNAMgSiteEngrapher.NoRNAInteraction, self).__init__(message)
    
    class FailsValidation(ValueError):
        MESG = "%s in structure %s does not pass validation criteria"
        def __init__(self, residueid, structid):
            self.structid = structid
            self.residueid = residueid
            message = self.MESG % (residueid, structid)
            super(RNAMgSiteEngrapher.FailsValidation, self).__init__(message)

    def __init__(self, hbfinder, depth=7, plinked=True, 
                 ligdists=None, modbases=True, invalids=True):
        self.depth = depth
        self.plinked = plinked
        self.hbfinder = hbfinder
        self._ligdists = ligdists if ligdists else self.LIGDISTS.copy()
        self._invalids = invalids
        if not modbases:
            self.NUCLEOTIDE_RESNAMES = ("A", "C", "G", "U")
        self.reset()

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

    def reset(self):
        super(RNAMgSiteEngrapher, self).reset()
        self._node = None
        self._atom = None
        self.residue = None
        self._valence = 0.0
        self.structure = None
        self.environment = set()
        self.inner_shell = set()
        self.outer_shell = set()

    @property
    def is_valid(self):
        valid_coordination = 4 <= len(self.inner_shell) <= 6
        return (valid_coordination and self._validate_environment()
                and self._validate_symmetry() and self._validate_valence())

    def set_residue(self, residue, structure):
        self.reset()
        self.structure = structure
        self.residue = residue
        self._atom = residue.atoms["MG"]
        self._node = (str(residue), "MG")
        self.graph = RNAMetalContextGraph(self._node, self.plinked)

    def engraph_mgsite_context(self, residue, structure):
        self.set_residue(residue, structure)
        self.detect_environment()
        self.detect_inner_shell(engraph=True)
        self.detect_outer_shell(engraph=True)
        structid = self.structure.structid
        residueid = self.graph.name
        if not self.residues_added:
            raise self.NoRNAInteraction(residueid, structid)
        is_validated = self.is_valid
        if not self._invalids and not is_validated:
            raise self.FailsValidation(residueid, structid)
        self.engraph_sequence_context()
        self.engraph_innershell_orientations()
        invalidchar = "" if is_validated else self.INVALIDCHAR
        self.graph.name = "%s %s%s" % (structid, residueid, invalidchar)

    def calculate_valence(self, atom, distance):
        valence = (self.VALCEPRMS[atom.atype] - distance) / 0.37
        valence = np.power(np.e, valence)
        return valence

    def detect_environment(self):
        get_nearby_atoms = self.structure.grid.get_nearby_atoms
        envmt = get_nearby_atoms(self.ENV_BOUND, self._atom)
        for residue, atom, distance in envmt:
            if atom.atype in self.VALCEPRMS:
                valence = self.calculate_valence(atom, distance)
                entry = (residue, atom, distance, valence)
                self.environment.add(entry)

    def detect_inner_shell(self, engraph=True):
        for entry in self.environment:
            residue, atom, distance, valence = entry
            if atom.atype in self._ligdists\
              and distance <= self._ligdists[atom.atype]:
                self._valence += valence
                self.inner_shell.add(entry)
        if engraph:
            self.engraph_shell(self.inner_shell, self.graph.INNER_SHELL_NO)

    def detect_outer_shell(self, engraph=True):
        find_hbonds = self.hbfinder.find_hbonds
        inner_shell = set(zip(*list(zip(*self.inner_shell))[:2]))
        for residue, atom in inner_shell:
            if residue.name in self.WATER_RESIDUE_NAMES:
                hbonds = find_hbonds(atom, residue, self.structure)
                for entry in hbonds:
                    _, bond_atom = entry
                    if entry not in inner_shell\
                      and bond_atom.atype in self.HBOND_ELEMENTS:
                        self.outer_shell.add(entry)
        if engraph:
            self.engraph_shell(self.outer_shell, self.graph.OUTER_SHELL_NO)

    def engraph_shell(self, shell, shell_no):
        add_ligand_interaction = self.graph.add_ligand_interaction
        for entry in shell:
            residue, atom = entry[:2]
            if residue.name not in self.WATER_RESIDUE_NAMES:
                resnode = (str(residue), residue.name)
                add_ligand_interaction(resnode, atom.name, shell_no)
                if residue.name in self.NUCLEOTIDE_RESNAMES:
                    self.residues_added.add(residue)

    def engraph_innershell_orientations(self):
        for entry1, entry2 in combinations(self.inner_shell, 2):
            (res1, atom1), (res2, atom2) = entry1[:2], entry2[:2]
            if res1.name not in self.WATER_RESIDUE_NAMES\
              and res2.name not in self.WATER_RESIDUE_NAMES:
                res1, res2 = (str(res1), res1.name), (str(res2), res2.name)
                angle = self._atom.angle_between(atom1, atom2)
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
        self.engraph_structure(self.structure, residue_set, True)

    def _id_sorter(self, nodeid):
        sorter = super(RNAMgSiteEngrapher, self)._id_sorter(nodeid)
        return (nodeid != self._node[0],) + sorter

    def _validate_valence(self):
        Q_v_value = 1 - abs(self._valence - 2) / 2.0
        return Q_v_value >= 0.5

    def _validate_symmetry(self):
        try:
            in_sh_vector = np.zeros((3,))
            for _, atom, distance, valence in self.inner_shell:
                vector = self._atom.vector_to(atom)
                in_sh_vector += (atom.occy * valence * vector) / distance
            Q_s_value = 1 - (np.linalg.norm(in_sh_vector) / self._valence)
            return Q_s_value >= 0.6
        except (TypeError, ZeroDivisionError):
            return False

    def _validate_environment(self):
        env_occy, env_bfac = self._env_occupancy_and_bfactor()
        try:

            if env_occy == self._atom.occy == 1.0:
                min_bfac, max_bfac  = sorted((self._atom.bfac, env_bfac))
                Q_e_value = min_bfac / max_bfac
            else:
                env_ratio = self._atom.bfac / self._atom.occy
                env_ratio /= env_bfac / env_occy
                env_ratio = min(env_ratio, 1/env_ratio)
                Q_e_value = min(env_occy, self._atom.occy) * env_ratio
            return Q_e_value >= 0.5
        except (TypeError, ZeroDivisionError):
            return False

    def _env_occupancy_and_bfactor(self):
        env_valence = 0.0
        env_bfactor = 0.0
        env_occupancy = 0.0
        for _, atom, _, valence in self.environment:
            env_valence += valence
            env_bfactor += atom.bfac * valence
            env_occupancy += atom.occy * valence
        return env_occupancy/env_valence, env_bfactor/env_valence


class RNAMgSiteList(object):

    def __init__(self, engrapher, encode=True):
        self.graphs = list()
        self.gcodes = list()
        self.encode = encode
        self._engrapher = engrapher

    def construct_site(self, residue, structure):
        self._engrapher.engraph_mgsite_context(residue, structure)
        self.graphs.append(self._engrapher.graph)
        if self.encode:
            self.gcodes.append(self._engrapher.encode_graph())
        self._engrapher.reset()

    def add_sites_from(self, structure):
        """ Add Mg-RNA sites in given structure to itself.

        Args:
            structure: engrid.Structure object of the RNA

        Returns:
            3-tuple (total sites, rna sites, shown sites)
            of integer counts, as described below.
              total sites: Total number of magnesium ions
                in the structure.
              rna sites: Total number of all sites where
                the ion interacts with RNA (regardless of
                validation criteria).
              shown sites: If engrapher excludes sites not
                passing validation criteria, these sites
                are not 'shown' in the output. Hence, the
                value may differ from the 'rna sites'.
        """
        rna_count = 0
        nonrna_count = 0
        invalid_count = 0
        for residue in structure.residues[structure.chainend:]:
            if residue.name == "MG":
                try:
                    self.construct_site(residue, structure)
                    rna_count += 1
                except RNAMgSiteEngrapher.NoRNAInteraction:
                    nonrna_count += 1
                except RNAMgSiteEngrapher.FailsValidation:
                    rna_count += 1
                    invalid_count += 1
        print(rna_count, nonrna_count, invalid_count)
        return rna_count+nonrna_count, rna_count, rna_count-invalid_count

    def dump_graphs(self, filepath, fileformat="tsv"):
        with open(filepath, 'w') as outfile:
            outfile.write("\n\n".join(self.gcodes))

