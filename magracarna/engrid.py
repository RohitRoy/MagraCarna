#!/usr/bin/env python

"""Module for representing coordinate files as a grid.

Quick Reference:
    Example:
        $ python engrid.py [DataFolder] [PDBList] [OutFolder] [LogFile]

Requirements:
    - Arguements:
        * Data Folder: contains *.out and *.nhx.co Files
        * PDB List   : has PDB IDs of required Files
        * Out Folder : to contain output files
        * Log File   : to contain runtime logs
    - In the directory where the Script is run:
        * residuify.py

Dependencies:
    * Python 2.7.3
    * Libraries:
        * numpy

Purpose of the Script:
    Merges *.cor and *.pdb information for each residue,
    Arranges *.cor residues by chains and residue numbers,
    Grids all *.pdb residues into cubic cells.

Actions Performed:
    * For each pdb ID * in [PDBList]
        * Reads *.ctx.cor and *.ctx.pdb in [DataFolder]
    * For each residue in *.ctx.cor
        * Retrieves residue identification from *.ctx.pdb
            * i.e. residue number and insertion code
        * Categorises residue based on chain ID
        * Further sorts residue based on residue number
    * Appends remaining *.ctx.pdb residues to residues-list
    * Builds empty grid with max, min coordinate info
        * For each residue, corresponding cell is found
        * Residue index in residues-list is added to cell
    * Builds structure with residues-list, chains, grid

Suggested Improvements:
    * Faster implementation for searching residues by
      using chains' information.
"""


import os
import sys

import re
import numpy as np
import warnings

from decimal import Decimal
from itertools import chain, product
from subprocess import call
from collections import namedtuple, defaultdict

from .aminoacids import donors as aa_donors, accors as aa_accors,\
                       hnames as aa_hnames
from .heteroatoms import donors as het_donors, accors as het_accors,\
                        hnames as het_hnames
from .nucleotides import donors as nuc_donors, accors as nuc_accors,\
                        hnames as nuc_hnames
from .nucleotides import allnucs, PO_DISTANCE


def showwarning_(message, category, filename, lineno, file=None, line=None):
    file = sys.stdout if file is None else file
    filename = os.path.basename(filename)
    warning = warnings.formatwarning(message, category, filename, str(lineno), line)
    sys.stderr.write(warning)

warnings.showwarning = showwarning_


class Atom(object):
    """Contains all information regarding atom.

    Attributes:
        name (str): atom name of atom.
        resnm (str): residue name.
        atmno (int): atom number.
        coos (numpy.array): (3,) shape array of atom
            coordinates.

    """

    def __init__(self, name, resname, atmno, atype, coos, occy, bfac):
        self.name = name
        self.resnm = resname
        self.atmno = atmno
        self.atype = atype
        self.coos = coos
        self.occy = occy
        self.bfac = bfac

    @classmethod
    def from_pdb_record(cls, line, iscor=False):
        if iscor:
            name = line[12:17].strip().replace("*", "'")
            atmno = int(line[4:11].strip())
            atype = None
        else:
            name = line[12:16].strip()
            atmno = int(line[6:11].strip())
            atype = line[76:78].strip()
        resname = line[17:20].strip()
        coos = np.array([line[30:38], line[38:46], line[46:54]], np.float64)
        occy = float(line[54:60]) if line[54:60].strip() else None
        bfac = float(line[60:66]) if line[60:66].strip() else None
        return cls(name, resname, atmno, atype, coos, occy, bfac)

    @classmethod
    def from_mmcif_row(cls, table, row_no):
        name = table["_atom_site.label_atom_id"][row_no].strip('"')
        resname = table["_atom_site.label_comp_id"][row_no]
        atmno = table["_atom_site.id"][row_no]
        atype = table["_atom_site.type_symbol"][row_no]
        coos = [float(table["_atom_site.Cartn_%s" % axis][row_no])
                                         for axis in ("x", "y", "z")]
        # limiting to the third decimal due to BPFind compatibility issues
        coos = np.array(["%.3f" % coo_ for coo_ in coos], dtype=np.float64)
        try:
            occy = float(table["_atom_site.occupancy"][row_no])
        except KeyError:
            occy = None
        try:
            bfac = float(table["_atom_site.B_iso_or_equiv"][row_no])
        except KeyError:
            bfac = None
        return cls(name, resname, atmno, atype, coos, occy, bfac)

    def __str__(self):
        return str((self.name, self.resnm, self.atype, self.atmno))

    def to_pdb_record(self, residue):
        args = [str(self.atmno), "% 4s" % ("% -3s" % self.name)]
        args += [residue.name, residue.chain, residue.resno, residue.insco]
        for value in self.coos:
            intval, decval = str(round(value, 3)).split(".")
            args.append("%s.%s" % (intval.rjust(4), decval.ljust(3, "0")))
        args.append("" if self.occy is None else "%.2f" % self.occy)
        args.append("" if self.bfac is None else "%.2f" % self.bfac)
        args = tuple(args + [self.atype])
        #            6  11 16    20  22  26  27---384654  60  66   78
        return "ATOM  % 5s %s % -3s% 2s% 4s% 1s   %s%s%s% 6s% 6s% 12s  " % args

    def vector_to(self, atom2):
        return atom2.coos - self.coos

    def distance_to(self, atom2):
        return np.linalg.norm(self.vector_to(atom2))

    def angle_between(self, atom1, atom2):
        """Calculates angle formed at the atom by two atoms.

        Args:
            atom1, atom2,: Atom objects (or None) in sequence
                to calculate angle between them

        Returns:
            Angle formed at this atom, in degrees, between
            the argument atom objects, ranging from 0 to
            180. If any argument is None, returns
            numpy.nan instead.
        """
        if not (atom1 or atom2):
            return np.nan
        bond1 = self.vector_to(atom1)
        bond2 = self.vector_to(atom2)
        norms = np.linalg.norm(bond1) * np.linalg.norm(bond2)
        cosine = (bond1.dot(bond2) / norms).clip(-1, 1)
        angle = np.degrees(np.arccos(cosine))
        return angle

    def dihedral_towards(self, atom1, atom2, atom3):
        bond1 = self.vector_to(atom1)
        bond2 = atom1.vector_to(atom2)
        bond3 = atom2.vector_to(atom3)
        for bond in (bond1, bond2, bond3):
            bond /= np.linalg.norm(bond)
        normal1 = np.cross(bond1, bond2)
        normal2 = np.cross(bond2, bond3)
        normal3 = np.cross(normal1, normal2)
        sin_part = bond2.dot(normal3)
        cos_part = normal2.dot(normal1)
        dihedral = np.degrees(np.arctan2(sin_part, cos_part))
        return dihedral


class Residue(object):
    """Contains all information regarding residue.

    Attributes:
        name: string residue name
        chain: string chain ID
        resno: integer author assigned residue number in
          structure file. If absent in MMCIF format, uses
          label_seq_id if it isn't '.' else _atom_site.id.
        serno: serial number in BPFind output's .cor file.
        insco: insertion code. If absent or '?', assigned
          an empty string. Assigned None for .cor files.
        index: used to represent unique ID in a structure,
          usually as an index to a reference list.
        atoms: maps atom names to Atom objects.
        structid: string ID of corresponding Structure.
    """

    REGEX = {"name": r"^\[(\w+)\]",
             "chain": r"^.*:(\w+)",
             "resno": r"^.*[\]](-?\d+)",
             "insco": r"^.*\^(\w+)"}

    def __init__(self, name=None, chain=None, structid=None, 
                       resno=None, insco=None, serno=None):
        """Initialises Residue object with given args.

        Args:
            All are optional, None if not assigned.
            name: string residue name.
            resno: integer residue number.
            serno: integer serial number in BPFind output.
            insco: string insertion code.
            chain: string chain ID.
            structid: string ID of corresponding Structure.
        """
        self.name = name
        self.resno = resno
        self.serno = serno
        self.insco = insco
        self.chain = chain
        self.index = None
        self.atoms = dict()
        self.structid = structid

    @classmethod
    def from_pdb_record(cls, line, iscor=False, asdict=False):
        """Creates Residue object from a PDB file record.

        Args:
            line: line of the file corresponding to the
              PDB record from which to create Residue.
            iscor: (default False) If True, treats 'line' 
              as a record in BPFind's .cor output file. 
              If False, treats 'line' as a record in a
              standard .pdb file format structure file.
            asdict: (default False) If True, returns a
              residue properties dict. If False, returns a
              Residue object.

        Returns:
            Residue object or a Residue properties dict,
            based on the value of 'asdict', where the name
            chain, resno and insco are assigned. 
            See Residue for details.
        """
        name = line[17:20].strip()
        props = {"name": name}
        if iscor:
            chain = None
            # New BPFind format does not keep chain IDS in .cor files
            isnew = int(SecondaryStructureManager.NEW_BPFIND)
            serno = props["serno"] = int(line[22-2*(isnew):26])
            insco = ""
        else:
            chain = props["chain"] = line[20:22].strip()
            resno = props["resno"] = int(line[22:26])
            insco = props["insco"] = line[26].strip()
        if asdict:
            return props
        return cls(**props)

    @classmethod
    def from_mmcif_row(cls, table, row_no, asdict=False):
        """Creates Residue object from an MMCIF table row.

        Args:
            table: a dict representing an MMCIF _atom_site
              category (table), mapping field names to a
              list of corresponding values in the column.
            row_no: row number in 'table' from which to 
              create a Residue object.
            asdict: (default False) If True, returns a
              residue properties dict. If False, returns a
              Residue object.

        Returns:
            Residue object or a Residue properties dict,
            based on the value of 'asdict', where the name
            chain, resno and insco are assigned. 
            See Residue for details.
        """
        try:
            resno = table["_atom_site.auth_seq_id"][row_no]
        except KeyError:
            resno = table["_atom_site.label_seq_id"][row_no]
            resno = table["_atom_site.id"][row_no] if resno == "." else resno
        name = table["_atom_site.label_comp_id"][row_no]
        chain = table["_atom_site.auth_asym_id"][row_no]
        resno = int(resno)
        try:
            insco = table["_atom_site.pdbx_PDB_ins_code"][row_no]
            insco = "" if insco == "?" else insco
        except KeyError:
            insco = ""
        props = {"name": name, "resno": resno, "insco": insco, "chain": chain}
        if asdict:
            return props
        return cls(**props)

    def __str__(self):
        """Returns string representation of residue.

        Returns:
            String representation to specify residue. 
            If 'resno' is assigned, uses the format used
            by RasMol to specify residues.
            Otherwise, returns string of the output of the
            entuple method.
        """
        if self.resno is not None:
            if self.insco:
                properties = (self.name, self.resno, self.insco, self.chain)
                return "[%s]%d^%s:%s" % properties
            else:
                return "[%s]%d:%s" % (self.name, self.resno, self.chain)
        return str(self.entuple())

    @classmethod
    def from_string(cls, string, asdict=False):
        """Creates Residue object from RasMol.

        The format used by RasMol can only specify residue
        properties mentioned in PDB file format, that are:
        'name', 'chain', 'resno' and 'insco'.

        Args:
            string: string in the format used by RasMol to
              specify residues.
            asdict: (default False) If True, returns a
              residue properties dict. If False, returns a
              Residue object.

        Returns:
            Residue object or a Residue properties dict,
            based on the value of 'asdict', where the name
            chain, resno and insco are assigned. 
            See Residue for details.
        """
        props = dict()
        for prop, regex in cls.REGEX.items():
            try:
                props[prop] = re.findall(regex, string)[0]
            except IndexError:
                continue
        if "resno" in props:
            props["resno"] = int(props["resno"])
            if "insco" not in props:
                props["insco"] = ""
        if asdict:
            return props
        return cls(**props)

    def entuple(self):
        """Returns unique residue identifier as tuple."""
        noserno = self.serno is None # with serno > without serno
        return (noserno, self.chain, self.serno,
                self.resno, self.insco, self.name)

    def matches(self, name=None, chain=None,
                      serno=None, resno=None, insco=None):
        if (name is not None and name != self.name)\
          or (chain is not None and chain != self.chain)\
          or (serno is not None and serno != self.serno):
            return False
        if resno is not None:
            insco = insco if insco is not None else ""
            if resno != self.resno or insco != self.insco:
                return False
        return True

    def __eq__(self, fles):
        """Checks if this Residue equals a given object.

        Checks for following attribute sets (in order):
            {index}, {serno}, {resno, insco},
            {atom.coos} of the first common atom name
            under their atoms attribute.
        If first set of attributes are None valued in
        either object, checks next set.

        Returns:
            If the attributes checked are not None:
                True if they are equal, False if not.
            False, if other object lacks attribute being
            checked, or no common atom name (in last set).
        """
        try:
            if None not in (self.index, fles.index):
                return self.index == fles.index
            elif None not in (self.serno, fles.serno):
                return self.serno == fles.serno
            elif None not in (self.resno, fles.resno, self.insco, fles.insco):
                return (self.resno == fles.resno and self.insco == fles.insco)
            else:
                for name in set(self.atoms).intersection(fles.atoms):
                    self_atom = self.atoms[name]
                    fles_atom = fles.atoms[name]
                    return (self_atom.coos == fles_atom.coos).all()
                return False
        except AttributeError:
            return False

    def __hash__(self):
        return id(self)


class Grid(object):
    """Construct to place atoms in a 3-D grid.

    Each cell in the Grid is a list of (residue, atom)
    corresponding to each atom, where rindex uniquely
    identifies the Residue to which the atom belongs.

    Attributes:
        size: size of each cell in Angstroms.
        cell: maps cell position tuples to cells.
    """

    def __init__(self, residues, size):
        """Creates Grid of given cell size for residues.

        Args:
            residues: list of Residue objects
            size: integer size of a cell, in Angstroms.

        Returns:
            Grid object initialised with all atoms in the
            given resides, at appropriate cell positions.
            An atom is identified by (residue, atom),
            rindex being the Residue index in given list.
        """
        limits = maxmm, minmm = self._limiting_coordinates(residues)
        limits = [(max_-min_)/size for max_, min_ in zip(*limits)]
        self._limit = np.array(np.ceil(limits), dtype=int)
        # limits hold for even maximum coordinate, since max_ +0.5 shifted
        self._min = minmm
        self.size = size
        self.cell = {pos_: list() for pos_ in self}
        for residue in residues:
            for atom in residue.atoms.values():
                cell = self.cell[self.get_cell_position(atom.coos)]
                cell.append((residue, atom))

    def get_cell_position(self, coos):
        """Retrieves cell position for given coordinates.

        Args:
            coos: (3, ) shape np.array with coordinates.

        Returns:
            Position of the cell as a 3-tuple.
        """        
        cell_position = np.floor((coos - self._min) / self.size)
        cell_position = tuple(int(index) for index in cell_position)
        return cell_position

    def get_cell(self, coos):
        cell_position = self.get_cell_position(coos)
        if self._cell_exists(cell_position):
            return self.cell[cell_position]
        return None

    def get_neighbours(self, positions, delta=1):
        """Returns concatenated info from adjacent cells.

        Args:
            positions: sequence of 3-tuple cell positions.

        Returns:
            Sorted, concatenated cells containing atoms in
            cells including and adjacent to 'positions'.
        """
        positions = np.array(sorted(set(positions)), dtype=int)
        neighbours = set()
        for addend in product(*(range(-delta, delta+1),)*3):
            adjacents = [tuple(pos_) for pos_ in positions + addend]
            neighbours.update(filter(self._cell_exists, adjacents))
        neighbour_atoms = sum([self.cell[pos_] for pos_ in neighbours], [])
        return neighbour_atoms

    def get_nearby_atoms(self, cutoff, entity):
        atomlist, positions = self._get_atoms_and_positions(entity)
        delta = int(np.ceil(float(cutoff) / self.size))
        nearby = list()
        for residue2, atom2 in self.get_neighbours(positions, delta):
            for atom1 in atomlist:
                distance = atom1.distance_to(atom2)
                if distance <= cutoff:
                    nearby.append((residue2, atom2, distance))
        return nearby

    @staticmethod
    def _limiting_coordinates(residues):
        minmm, maxmm = np.array([np.inf]*3), np.array([-np.inf]*3)
        for residue in residues:
            for atom in residue.atoms.values():
                minmm = np.min([minmm, atom.coos], 0)
                maxmm = np.max([maxmm, atom.coos], 0)
        maxmm = tuple(np.array(np.ceil(maxmm+[0.1]*3), dtype=int))
        minmm = tuple(np.array(np.floor(minmm-[0.1]*3), dtype=int))
        return maxmm, minmm

    def _cell_exists(self, position):
        """Returns True iff given cell position is valid.
        """
        return all(0 <= pos_ for pos_ in position)\
          and all(position < self._limit)

    def _get_atoms_and_positions(self, entity):
        try:
            atom_coos = entity.coos
            positions = {self.get_cell_position(atom_coos)}
            return [entity], positions
        except AttributeError:
            atomlist = entity.atoms.values()
            positions = set()
            for atom in atomlist:
                positions.update(self.get_cell_position(atom.coos))
            return atomlist, positions

    def __iter__(self):
        """Returns iterator over all cell positions.
        """
        return product(*[range(limit) for limit in self._limit])


class Structure(object):
    """Constructs structure from coordinates.

    The 'chains' attribute accounts for the nucleotide 
    residues in the structure that have the same chain ID.
    The beginning and ending residue indices are recorded,
    as per their order in the 'residues' attribute.
    Optionally, when constructed with 'breaks' argument as
    True, 'chains' only accounts for those nucleotides
    that are linked by P-O3' linkage with a neighbouring
    nucleotide of the same chain. These nucleotides may
    also be connected in fragments, and a breaks in the 
    chain is recorded as (startex, jump). There is a break
    entry for each fragment except for the first.
        startex: Index at which a fragment starts.
        jump: Estimated number of missing residues between
          this fragment and the previous.
    When constructed with 'breaks' argument as True,
    'chains' excludes nucleotides which are not in a 
    phospho-diester linkage, or have unusual atom names in
    the phosphate moiety.

    The index attribute of each residue in 'residues' is
    matches its index in the 'residues' list.

    Attributes:
        structid: string ID for structure (usually PDB ID)
        residues: sequence of all residues.
        grid: Grid object which grids all atoms.
        chains: Maps chain ID to (begex, endex, breaks)
          begex, endex: beginning and ending indices of
            the chain in 'residues', respectively.
          breaks: sequence of (startex, jump) break tuples
            as described above. It is empty if Structure
            is constructed with 'breaks' argument as False.
        chainend: Index of the first residue in 'residues'
          which is not accounted for by 'chains'.
        pairs: None, or a mapping from residues to dict of
          pairing-residue mapped with base pair type.
    """

    _WARN_ONE_RESIDUE_SEGT = "Removed one residue segment %%s at %%d / %d"

    def __init__(self, structid, residues, grid, pairs=None, breaks=False):
        """Initialises Structure with provided 

        During construction, the residues that can be
        accounted for by the 'chains' attribute are placed
        at the start of the list in 'residues' attribute. 
        So, the order of the residues in the 'residues'
        argument may not be the retained in the 'residues'
        attribute. The index attribute of each Residue is
        appropriately reassigned to match its new index.

        Args:
            structid: string ID for structure.
            residues: sequence of all residues.
            grid: Grid object which grids all atoms in the
              provided 'residues'.
            pairs: (optional) a dict where (chain, serno)
              of a residue is mapped to a dict object, 
              that in turn maps (chain, serno) of its 
              pairing residues to a 5-tuple representing
              the type of the base pair they form.
            breaks: If True, the 'chains' attribute will
              only account for phosphate-linked segments
              of the nucleotide chain, and the breaks in
              the chains are aptly recorded (for details,
              see Structure). If False, every nucleotide
              is accounted for by the 'chains' attribute.
        """
        self.grid = grid
        self.structid = structid
        self.residues = residues
        self.chains = None
        self.chainend = 0
        self._chain_nucleotides()
        if breaks:
            self._chain_nucleotides_with_breaks()
        if pairs:
            pairs = self._match_residues_to_base_pairs(pairs)
        self.pairs = pairs

    def find(self, **props):
        if "chain" in props:
            if props["chain"] in self.chains:
                props_ = props.copy()
                begex, endex, _ = self.chains[props_.pop("chain")]
                for residue in self.residues[begex:endex]:
                    if residue.matches(**props_):
                        return residue
            for residue in self.residues[self.chainend:]:
                if residue.matches(**props):
                    return residue
        for residue in self.residues:
            if residue.matches(**props):
                return residue
        return None

    @staticmethod
    def linked(nuc_at5, nuc_at3):
        if "O3'" in nuc_at5.atoms and "P" in nuc_at3.atoms:
            patom = nuc_at3.atoms["P"]
            o3atom = nuc_at5.atoms["O3'"]
            return patom.distance_to(o3atom) <= PO_DISTANCE
        return False

    def get_connected_segments(self):
        segments = list()
        for chain, (begex, endex, breaks) in self.chains.items():
            if breaks:
                breakpoints = (begex,) + tuple(zip(*breaks))[0] + (endex,)
            else:
                breakpoints = (begex, endex)
            for seg_begex, seg_endex in zip(breakpoints, breakpoints[1:]):
                segments.append((seg_begex, seg_endex, chain))
        return sorted(segments)

    def _chain_nucleotides(self):
        self.chains = dict()
        self.chainend = 0
        hetlist, nuclist = list(), list()
        chains = defaultdict(set)
        for residue in self.residues:
            (nuclist if residue.name in allnucs else hetlist).append(residue)
        # sorter = lambda residue: (residue.chain, residue.resno, residue.insco)
        # nuclist.sort(key=sorter)
        for res_index, residue in enumerate(nuclist):
            chains[residue.chain].add(res_index)
            if residue.chain not in self.chains:
                self.chains[residue.chain] = [res_index, res_index+1, []]
            else:
                self.chains[residue.chain][1] = res_index+1
        self.chainend = len(nuclist)
        self.residues = nuclist + hetlist
        for res_index, residue in enumerate(self.residues):
            residue.index = res_index

    def _chain_nucleotides_with_breaks(self):
        chained = list()
        unchained = list()
        het_residues = self.residues[self.chainend:]
        for chain, breaks, chained_, unchained_ in self._find_chain_breaks():
            begex = len(chained)
            breaks = [(begex+start, jump) for start, jump in breaks]
            chained += chained_
            unchained += unchained_
            self.chains[chain] = (begex, len(chained), breaks)
        self.residues = chained + unchained + het_residues
        for index, residue in enumerate(self.residues[:self.chainend]):
            residue.index = index
        self.chainend = len(chained)

    def _match_residues_to_base_pairs(self, pairs):
        matched = defaultdict(dict)
        for (chain1, serno1), pairs_with in pairs.items():
            residue1 = self.find(chain=chain1, serno=serno1)
            if residue1:
                for (chain2, serno2), bptype in pairs_with.items():
                    residue2 = self.find(chain=chain2, serno=serno2)
                    if residue2:
                        matched[residue1][residue2] = bptype
        return matched

    def _find_chain_breaks(self):
        sorter = lambda res_: res_.index
        chain_breaks = list()
        for chain, (begex, endex, _) in list(self.chains.items()):
            chain_residues = self.residues[begex:endex]
            segments = self._find_segments(chain_residues)
            self._remove_single_residue_segments(segments)
            segments.sort(key=lambda segment: segment[0].index)
            breaks = self._calculate_chain_breaks(segments)
            chained = sum(segments, list())
            unchained = sorted(set(chain_residues) - set(chained), key=sorter)
            chain_breaks.append((chain, breaks, chained, unchained))
        return chain_breaks

    def _find_segments(self, chain_residues):
        segments = [[chain_residues[0], ]]
        for residue in chain_residues[1:]:
            for segment in reversed(segments):
                if self.linked(segment[-1], residue):
                    segment.append(residue)
                    break
                elif self.linked(residue, segment[0]):
                    segment.insert(0, residue)
                    break
            else:
                segments.append([residue])
        for segment in segments:
            if segment[0].index > segment[-1].index:
                segment.reverse()
        return segments

    @classmethod
    def _remove_single_residue_segments(cls, segments):
        segments_number = len(segments)
        warn_disconnect = cls._WARN_ONE_RESIDUE_SEGT % segments_number
        for index in range(segments_number-1, -1, -1):
            segment = segments[index]
            if len(segment) == 1:
                residue = segments.pop(index)[0]
                # warnings.warn(warn_disconnect % (residue, index+1))

    @classmethod
    def _calculate_chain_breaks(cls, segments):
        breaks = list()
        break_start = 0
        duplicate_residue = "Duplicate residue numbering at %s and %s"
        for seg1, seg2 in zip(segments, segments[1:]):
            jump = cls._calculate_jump_value(seg1[-1], seg2[0])
            if jump == -1:
                raise ValueError(duplicate_residue % (seg1[-1], seg2[0]))
            break_start += len(seg1)
            breaks.append((break_start, jump))
        return breaks

    @staticmethod
    def _calculate_jump_value(residue1, residue2):
        if residue1.resno == residue2.resno:
            if residue1.insco:
                return ord(residue2.insco) - ord(residue1.insco) - 1
            return ord(residue2.insco) - (ord('A')-1) - 1
        return residue2.resno - residue1.resno - 1


class HydrogenBondFinder(object):

    MAX_DA = Decimal('3.80')
    MIN_DA = Decimal('2.10')
    MAX_HA = Decimal('2.70')
    MIN_HA = Decimal('1.00')
    MAX_DHA = Decimal('180')
    MIN_DHA = Decimal('110')

    DONORS = chain(nuc_donors.items(), aa_donors.items(), het_donors.items())
    ACCORS = chain(nuc_accors.items(), aa_accors.items(), het_accors.items())
    HNAMES = chain(nuc_hnames.items(), aa_hnames.items(), het_hnames.items())
    DONORS, ACCORS, HNAMES = dict(DONORS), dict(ACCORS), dict(HNAMES)

    class InvalidHBondParameters(ValueError):
        """
        """

    def __init__(self, ignore_H=True,
                 max_da=Decimal('3.80'), min_da=Decimal('2.10'),
                 max_ha=Decimal('2.70'), min_ha=Decimal('1.00'),
                 max_dha=Decimal('180'), min_dha=Decimal('110')):
        self._recall = None
        self._structure = None
        self._max_da, self._min_da = max_da, min_da
        self._max_ha, self._min_ha = max_ha, min_ha
        self._max_dha, self._min_dha = max_dha, min_dha
        self._set_verify_hbond(ignore_H)
        self._validate_parameters()

    def set_structure(self, structure=None):
        if not structure:
            self._structure = None
            self._recall = None
            self._get_nearby_atoms = None
        elif self._structure != structure:
            self._structure = structure
            self._recall = dict()
            self._near = structure.grid.get_nearby_atoms

    def find_hbonds(self, atom, residue, exclude=dict()):
        if atom in self._recall:
            return self._exclude(self._recall[atom], exclude)
        hbonds = list()
        is_donor = self._is_donor(residue, atom)
        is_accor = self._is_accor(residue, atom)
        if is_donor or is_accor:
            for residue2, atom2, da_dist in self._near(self._max_da, atom):
                if residue2 != residue:
                    if is_donor and self._is_accor(residue2, atom2)\
                      and self.verify_hbond(atom, residue, atom2, da_dist):
                        hbonds.append((residue2, atom2))
                    if is_accor and self._is_donor(residue2, atom2)\
                      and self.verify_hbond(atom2, residue2, atom, da_dist):
                        hbonds.append((residue2, atom2))
        self._recall[atom] = hbonds
        return self._exclude(hbonds, exclude)

    @staticmethod
    def _exclude(hbonds, exclude):
        if not exclude:
            return hbonds[:]
        excluded = list()
        for residue, atom in hbonds:
            if residue not in exclude or atom not in exclude[residue]:
                excluded.append((residue, atom))
        return excluded

    def find_all_hbonds(self, residues, expand=False):
        all_hbonds = list()
        found_hbonds = defaultdict(lambda: defaultdict(set))
        for residue1 in residues:
            for name1, atom1 in residue1.atoms.items():
                exclude = found_hbonds[(residue1, atom1)]
                hbonds = self.find_hbonds(atom1, residue1, exclude)
                for residue2, atom2 in hbonds:
                    if expand or residue2 in residues:
                        found_hbonds[(residue2, atom2)][residue1].add(atom1)
                        all_hbonds.append((residue1, atom1, residue2, atom2))
        return all_hbonds
    
    def _set_verify_hbond(self, ignore_H=True):
        if ignore_H:
            self.verify_hbond = self._verify_hbond_ignoring_H
        else:
            self.verify_hbond = self._verify_hbond_that_has_H

    def _validate_parameters(self):
        if not (self.MIN_DHA <= self._min_dha <= self._max_dha <= self.MAX_DHA
          and self.MIN_DA <= self._min_da <= self._max_da <= self.MAX_DA
          and self.MIN_HA <= self._min_ha <= self._max_ha <= self.MAX_HA):
            raise self.InvalidHBondParameters()

    def _verify_hbond_ignoring_H(self, donor, donor_residue, accor, da_dist):
        return self._min_da <= da_dist <= self._max_da

    def _verify_hbond_that_has_H(self, donor, donor_residue, accor, da_dist):
        hatom_names = set(self.HNAMES[donor_residue.name][donor.name])
        for hatom_name in (hatom_names & set(donor_residue.atoms)):
            hatom = donor_residue.atoms[hatom_name]
            ha_distance = accor.distance_to(hatom)
            dha_angle = hatom.angle_between(donor, accor)
            if self._min_dha <= dha_angle <= self._max_dha\
              and self._min_ha <= ha_distance <= self._max_ha\
              and self._min_da <= da_distance <= self._max_da:
                return True
        return False

    @classmethod
    def _is_donor(cls, residue, atom):
        return (residue.name in cls.DONORS 
                and atom.name in cls.DONORS[residue.name])

    @classmethod
    def _is_accor(cls, residue, atom):
        return (residue.name in cls.ACCORS 
                and atom.name in cls.ACCORS[residue.name])


class StructureParser(object):

    class FileFormatError(ValueError):
        """
        """


class PDBParser(StructureParser):

    _COORDINATE_SECTION_RECORDS = ("ATOM", "HETATM", "MODEL", "ENDMDL")
    _COORDINATE_RECORDS = ("ATOM", "HETATM")
    
    def __init__(self, lines, iscor=False):
        self.lines = self.retrieve_coordinate_section(lines)
        self.iscor = iscor

    def ensure_single_model(self):
        if self.iscor:
            return
        newlines = list()
        found_first_model = False
        ignore_current_model = False
        for record in self.lines:
            if record.startswith("MODEL"):
                if not found_first_model:
                    found_first_model = True
                else:
                    ignore_current_model = True
            if not ignore_current_model:
                newlines.append(record)
            if record.startswith("ENDMDL"):
                ignore_current_model = False
        self.lines = newlines

    def ensure_no_alternate_locations(self):
        if self.iscor:
            return
        chosen_locs = dict()
        delete_lines = list()
        model_number = 0
        for lineno, line in enumerate(self.lines):
            if line.startswith("MODEL"):
                model_number = int(line.split()[1])
            elif line.startswith(self._COORDINATE_RECORDS):
                altloc = line[16].strip()
                if altloc:
                    occy = float(line[54:60])
                    atomkey = (model_number, line[12:16], line[17:27])
                    if atomkey in chosen_locs:
                        occy0, altloc0, lineno0 = chosen_locs[atomkey]
                        if (occy0 > occy)\
                          or (occy0 == occy and altloc0 < altloc):
                            delete_lines.append(lineno)
                            continue # in case the previous was better
                        delete_lines.append(lineno0)
                    chosen_locs[atomkey] = (occy, altloc, lineno)
        for _, _, lineno in chosen_locs.values():
            newline = self.lines[lineno][:16]+" "+self.lines[lineno][17:]
            self.lines[lineno] = newline
        for lineno in sorted(delete_lines, reverse=True):
            self.lines.pop(lineno)

    @classmethod
    def retrieve_coordinate_section(cls, lines):
        newlines = list()
        for record in lines:
            if record.startswith(cls._COORDINATE_SECTION_RECORDS):
                newlines.append(record)
        return newlines

    def convert_back_to_file_format(self):
        return "".join(self.lines)

    def extract_atoms(self):
        for line in self.lines:
            if line.startswith(self._COORDINATE_RECORDS):
                atom = Atom.from_pdb_record(line, self.iscor)
                resprops = Residue.from_pdb_record(line, self.iscor, True)
                yield (atom, resprops)


class MMCIFParser(StructureParser):

    def __init__(self, lines):
        self.tables = list()
        self._columns = list()
        self._retrieve_coordinate_tables(lines)

    def _retrieve_coordinate_tables(self, lines):
        for block in self.partition_into_blocks(lines)[1:]:
            prefix, table, columns = self.convert_block_to_table(block)
            table["rows"] = len(table[columns[0]])
            self._columns.append(columns)
            self.tables.append([prefix, table])

    @staticmethod
    def convert_block_to_table(block):
        prefix = [block[0]]
        columns = list()
        tablerows = list()
        looping_over_atoms = False
        for lineno, line in enumerate(block):
            if line.startswith("_atom_site."):
                if block[lineno-1].startswith("loop_"):
                    prefix.append(block[lineno-1])
                    looping_over_atoms = True
                columns.append(line.strip())
            elif looping_over_atoms:
                if line.startswith(("_", "loop_", "stop_")):
                    looping_over_atoms = False
                else:
                    tablerows.append(line.split())
        table = dict(zip(columns, zip(*tablerows)))
        return prefix, table, columns

    @classmethod
    def partition_into_blocks(cls, lines):
        blocks = [list()]
        try:
            for line in lines:
                if line.startswith("data_"):
                    blocks.append(list())
                if not line.startswith("#"):
                    blocks[-1].append(line)
        except:
            raise cls.FileFormatError("Uncommented lines outside data block")
        return blocks

    def ensure_no_alternate_locations(self):
        for (_, table), columns in zip(self.tables, self._columns):
            try:
                alternate_locs = table["_atom_site.label_alt_id"]
            except KeyError:
                continue
            try:
                model_numbers = table["_atom_site.pdbx_PDB_model_num"]
            except:
                model_numbers = len(alternate_locs)*[""]

            atom_names = table["_atom_site.label_atom_id"]
            occupancies = table["_atom_site.occupancy"]
            delete_rows = list()
            chosen_locs = dict()
            for row_no, altloc in enumerate(alternate_locs):
                if altloc != ".":
                    occy = float(occupancies[row_no])
                    atom_name = atom_names[row_no]
                    model_number = model_numbers[row_no]
                    restuple = Residue.from_mmcif_row(table, row_no).entuple()
                    atomkey = (model_number, atom_name, restuple)
                    if atomkey in chosen_locs:
                        occy0, altloc0, row_no0 = chosen_locs[atomkey]
                        if (occy0 > occy)\
                          or (occy0 == occy and altloc0 < altloc):
                            delete_rows.append(row_no)
                            continue # in case the previous was better
                        delete_rows.append(row_no0)
                    chosen_locs[atomkey] = (occy, altloc, row_no)
            for _, _, row_no in chosen_locs.values():
                alternate_locs[row_no] = "."
            for row_no in sorted(delete_rows, reverse=True):
                for column in columns:
                    table[column].pop(row_no)
            table["rows"] -= len(delete_rows)

    def ensure_single_model(self):
        for (_, table), columns in zip(self.tables, self._columns):
            try:
                model_numbers = table["_atom_site.pdbx_PDB_model_num"]
            except:
                continue
            model_numbers = [int(model_no) for model_no in model_numbers]
            first_model = min(model_numbers)
            retain_rows = list()
            for row_no, model_number in enumerate(model_numbers):
                if first_model == model_number:
                    retain_rows.append(row_no)
            for column in columns:
                col_values = table[column]
                table[column] = [col_values[row_no] for row_no in retain_rows]
            table["rows"] = len(retain_rows)

    def ensure_single_data_block(self):
        self.tables = self.tables[:1]

    def convert_back_to_file_format(self):
        blocks = list()
        for (prefix, table), columns in zip(self.tables, self._columns):
            block = prefix[:]
            tablecols = list()
            for column in columns:
                block.append("%s \n" % column)
                col_values = table[column]
                string_len = max(len(value) for value in col_values) + 1
                col_values = [value.ljust(string_len) for value in col_values]
                tablecols.append(col_values)
            block += ["".join(row)+"\n" for row in zip(*tablecols)]
            blocks.append(block)
        return "".join(["".join(block)+"#\n" for block in blocks])

    def extract_atoms(self):
        for _, table in self.tables:
            for row_no in range(table["rows"]):
                resprops = Residue.from_mmcif_row(table, row_no, True)
                atom = Atom.from_mmcif_row(table, row_no)
                yield (atom, resprops)


class StructureFile(object):

    def __init__(self, filepath, scratch=None, is_cor=False):
        self.filepath = filepath
        self.structid = os.path.basename(filepath).split(".")[0]
        self.write_to = None
        self._parser = None
        self.format = None
        self.is_cor = is_cor
        if scratch:
            self.write_to = os.path.join(scratch, os.path.basename(filepath))
            if not os.path.exists(scratch) or not os.path.isdir(scratch):
                raise ValueError("Folder %s does not exist" % scratch)

    def __enter__(self):
        self.read()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self.write_to:
            self.write()

    def read(self):
        with open(self.filepath, 'r') as original_file:
            lines = original_file.readlines()
        for line in lines:
            if not line.startswith("#"):
                if line.startswith("data_"):
                    self._parser = MMCIFParser(lines)
                    self.format = "CIF"
                else:
                    self._parser = PDBParser(lines, self.is_cor)
                    self.format = "PDB"
                break

    def modify(self, multimodel=False, altlocs=False, multiblock=False):
        if not multimodel:
            self._parser.ensure_single_model()
        if not altlocs:
            self._parser.ensure_no_alternate_locations()
        if not multiblock and self.format == "CIF":
            self._parser.ensure_single_data_block()

    def write(self, write_to=None):
        write_to = self.write_to if write_to is None else write_to
        if write_to is None:
            raise ValueError("Destination file path not provided")
        with open(write_to, 'w') as scratch_file:
            scratch_file.write(self._parser.convert_back_to_file_format())

    def extract_residues(self, only_elements=None, 
                         only_residues=None, only_ids=None):
        """Converts file into list of Residue objects.

        Args:
            only_elements: an optional set of elements.
                If not None, only those Atom objects are 
                read which have elements in this set.
            only_residues: optional set of residue names.
                If not None, Atom objects from only these
                residue types are read.
            only_ids: optional set of residue IDs.
                If not None, Atom objects with only these
                residue IDs are read. Should not be used
                for cor files.

        Yields:
            Residue objects read from the file.
        """
        if self.is_cor and only_ids:
            raise ValueError("Residue IDs not accepted for COR files.")
        keep_residue = False
        current_residue = Residue() # empty residue
        for atom, resprops in self._parser.extract_atoms():
            if (only_elements and atom.atype not in only_elements)\
              or (only_residues and atom.resnm not in only_residues):
                continue
            new_residue = None
            if not current_residue.matches(**resprops):
                if keep_residue:
                    yield current_residue
                new_residue = Residue(structid=self.structid, **resprops)
                keep_residue =(not only_ids) or (str(new_residue) in only_ids)
                current_residue = new_residue
            if keep_residue:
                current_residue.atoms[atom.name] = atom
        if keep_residue:
            yield current_residue

    def extract_structure(self, size=7.0, breaks=True, multimodel=False,
                             altlocs=False, multiblock=False, **limits):
        self.modify(multimodel, altlocs, multiblock) # is it necessary?
        residues = list(self.extract_residues(**limits))
        grid = Grid(residues, size)
        return Structure(self.structid, residues, grid, None, breaks)


class SecondaryStructureManager(object):

    # New BPFind format does not keep chain IDs in .cor files
    NEW_BPFIND = True

    class BPFindParams(object):
        """Specify parameters for a BPFind run.
        """
        def __init__(self, scratch=".", cleanup=True, verbose=False,
                     hbond_params=HydrogenBondFinder()):
            self.scratch = os.path.abspath(scratch)
            self.da_dist = hbond_params._max_da
            self.cleanup = cleanup
            self.verbose = verbose

    BPFIND_COMMAND = ["bpfind.linux", "-HT"]
    _OUTPUT_EXTS = ("nup", "dat", "hlx", "fasta", "out", "cor")

    WARN = True

    class InvalidFileNameError(ValueError):
        """
        """

    def __init__(self, structurepath, run_bpfind=None, run_nuparm=None):
        if run_bpfind:
            self._scratch = run_bpfind.scratch
            self._da_dist = run_bpfind.da_dist
            self.clean_up = run_bpfind.cleanup
            self._verbose = run_bpfind.verbose
        else:
            self._scratch = self._da_dist = None
            self.clean_up = False
            self._verbose = False
        self.validate_filename(structurepath)
        with StructureFile(structurepath, self._scratch) as structfile:
            structfile.modify(False, False, False)
            self.structfile = structfile
        self.corfile = None
        self.pairs = None

    @classmethod
    def validate_filename(cls, filepath):
        filename = os.path.basename(filepath)
        message = None
        if len(filename) > 50:
            message = "'%s' is longer than 50 characters" % filename
            raise cls.InvalidFileNameError(message)
        elif "/" in filename:
            message = "'%s' contains '/'"
        elif filename.startswith(" ") or filename.startswith("-"):
            message = "'%s' starts with invalid character '%s'"
            message = message % (filename, filename[0])
        elif "." not in filename:
            message = "'%s' has no file extension" % filename
        else:
            extension = ".".join(filename.split(".")[1:])
            if len(extension) > 4:
                message = "'%s' has file extension longer than 4 characters"
            elif extension in cls._OUTPUT_EXTS:
                message = "'%s' has an invalid file extension\n"
                message += "List of invalid file extensions: %s"
                message = message % (filename, ", ".join(cls._OUTPUT_EXTS))
        if message:
            raise cls.InvalidFileNameError(message)

    def run_bpfind(self):
        original_folder = os.path.abspath(os.curdir)
        os.chdir(self._scratch)
        structfilename = os.path.basename(self.structfile.write_to)
        command = self.BPFIND_COMMAND[:]
        if self.structfile.format == "CIF":
            command.append("-cif")
        command += ["-HD ", str(self._da_dist), structfilename]
        with open(os.devnull, 'w') as NULL:
            error_out = sys.stderr if self._verbose else NULL
            exit_code = call(command, stdout=NULL, stderr=error_out)
        os.chdir(original_folder)

    def read_bpfind_output_files(self):
        input_folder = os.path.dirname(self.structfile.filepath)
        if self._scratch:
            input_folder = self._scratch
        input_files = os.path.join(input_folder, self.structfile.structid)
        corfilepath = "%s.cor" % input_files
        outfilepath = "%s.out" % input_files
        self.corfile = StructureFile(corfilepath, scratch=None, is_cor=True)
        self.corfile.read()
        self.pairs = self._extract_base_pairs(outfilepath)

    def __enter__(self):
        if self._scratch:
            self.run_bpfind()
        self.read_bpfind_output_files()
        return self

    def clean_up_bpfind_output(self):
        input_folder = os.path.dirname(self.structfile.filepath)
        if self._scratch:
            input_folder = self._scratch
        input_files = os.path.join(input_folder, self.structfile.structid)
        for extension in self._OUTPUT_EXTS:
            expected_file = "%s.%s" % (input_files, extension)
            if os.path.exists(expected_file):
                os.remove(expected_file)

    def __exit__(self, exc_type, exc_value, traceback):
        if self.clean_up:
            self.clean_up_bpfind_output()

    def extract_paired_structure(self, size=7.0, breaks=True, **limits):
        if None in [self.corfile, self.pairs]:
            self.read_bpfind_output_files()
        structresidues = list(self.structfile.extract_residues(**limits))
        limits.pop("residueids", None)
        corresidues = self.corfile.extract_residues(**limits) # why not limits?
        structid = self.structfile.structid
        grid = self._grid_merged_residues(corresidues, structresidues, size)
        return Structure(structid, structresidues, grid, self.pairs, breaks)

    def _extract_base_pairs(self, outfilepath):
        pairs = defaultdict(dict)
        for line in self._split_into_pairs(outfilepath):
            line = line.strip()
            if line:
                line = line.split()
                serno1, resno1, resname1, chain1 = line[:4]
                serno2, resno2, resname2, chain2 = line[4:8]
                bptype = tuple(char for char in line[8] if char != ":")
                pairs[(chain1, int(serno1))][(chain2, int(serno2))] = bptype
        return dict(pairs)

    @staticmethod
    def _split_into_pairs(outfilepath):
        rewritten = list()
        with open(outfilepath, 'r') as bpfind_outfile:
            for line in bpfind_outfile:
                if line.startswith(" ") and line.strip():
                    pairs_line = line.split()
                    main_base = "\t"+"\t".join(pairs_line[:4])
                    pairs_line = pairs_line[4:]
                    while pairs_line:
                        pairs_with = "\t"+"\t".join(pairs_line[:7])
                        pairs_line = pairs_line[7:]
                        rewritten.append("%s%s\n" % (main_base, pairs_with))
        return rewritten

    def _grid_merged_residues(self, corresidues, structresidues, size):
        grid = Grid(structresidues, size)
        cor_resnames = set()
        for corres in corresidues:
            cor_resnames.add(corres.name)
            for _, atom in corres.atoms.items():
                cell = grid.get_cell(atom.coos)
                if cell and self._merge_residue_into_cell(cell, corres, atom):
                    break
            else:
                if self.WARN:
                    warnings.warn("Can't find COR Residue %s in PDB" % corres)
        for residue in structresidues:
            if residue.name in cor_resnames and residue.serno is None:
                if self.WARN:
                    warnings.warn("PDB Residue not found in COR %s" % residue)
        return grid
        
    @staticmethod
    def _merge_residue_into_cell(cell, corres, coratom):
        for residue, atom in cell:
            if (atom.coos == coratom.coos).all():
                corres.resno = residue.resno
                corres.insco = residue.insco
                corres.atoms = residue.atoms
                residue.serno = corres.serno
                return True
        return False

