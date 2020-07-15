from __future__ import print_function

import os
import sys
import argparse

import numpy as np

from ...magracarna.engrid import Structure, HydrogenBondFinder as HBFinder,\
							  SecondaryStructureManager as SSManager
from ...magracarna.heteroatoms import cations
from ...magracarna.nucleotides import rings


SSManager.WARN = False
SSManager.NEW_BPFIND = False

required = {("U", "A", "W", "H", "T"): [0, 0]}


def get_basepairs(structure):
	basepairs = list()
	required_residues = {bptype[0] for bptype in required}
	for residue in structure.pairs:
		if residue.name in required_residues:
			itspairs = structure.pairs[residue].items()
			for pairresidue, bptype in itspairs:
				bptype = (residue.name, pairresidue.name)+bptype
				if bptype in required:
					basepairs.append((residue, pairresidue, bptype))
	return basepairs


def print_triplets(basepairs, structure):
	for residue, pairresidue, bptype in basepairs:
		print("%s\t%s\t%s" % (str(residue), str(pairresidue), " ".join(bptype)))
		residue_pairs = structure.pairs[residue].items()
		if len(residue_pairs) > 1:
			print("\tPairs Of:\t%s" % str(residue))
			for tripresidue, triptype in residue_pairs:
				print("\t\t%s\t%s" % (tripresidue, " ".join(triptype)))
		pairresidue_pairs = structure.pairs[pairresidue].items()
		if len(pairresidue_pairs) > 1:
			print("\tPairs Of:\t%s" % str(pairresidue))
			for tripresidue, triptype in pairresidue_pairs:
				print("\t\t%s\t%s" % (tripresidue, " ".join(triptype)))


def uawht_analyse(basepairs, structure):
	sorter = lambda entry: entry[2]
	structid = structure.structid
	get_nearby_atoms = structure.grid.get_nearby_atoms
	for u_residue, a_residue, bptype in basepairs:
		if bptype[:2] != ("U", "A"):
			continue
		site = "%s %s %s" % (structid, u_residue, a_residue)
		try:
			a_p = a_residue.atoms["P"]
			a_n9 = a_residue.atoms["N9"]
			a_n7 = a_residue.atoms["N7"]
			a_c1_ = a_residue.atoms["C1'"]
			# 
			u_o4 = u_residue.atoms["O4"]
			a_oph = set(map(a_residue.atoms.get, ("OP1", "OP2", "OP3")))
			a_oph = min(a_oph-{None}, key=lambda atom: u_o4.distance_to(atom))
		except KeyError as e:
			print("\t%s\t%s" % (site, e.message))
			continue
		# 
		in_plane = a_c1_.angle_between(a_oph, a_n9)
		distance = u_o4.distance_to(a_oph)
		bond1 = (a_n7.coos - a_n9.coos)
		bond2 = (u_o4.coos - a_n9.coos)
		normal = np.cross(bond1 / np.linalg.norm(bond1), bond2 / np.linalg.norm(bond2))
		off_plane = np.dot(normal, a_oph.coos-a_n9.coos)
		# 
		near_o4 = sorted(get_nearby_atoms(3.8, u_o4), key=sorter)
		near_o4 = {(residue, atom.name) for residue, atom, _ in near_o4}
		near_oph = sorted(get_nearby_atoms(3.8, a_oph), key=sorter)
		near_oph = {(residue, atom.name) for residue, atom, _ in near_oph}
		nearest = ("", "", np.inf)
		near_O, near_water = False, False
		near_N, near_amino = False, False
		near_MG, near_cation = False, False
		near_lysarg = False
		for residue, atom in sorted(near_o4 & near_oph):
			if residue not in (u_residue, a_residue) and residue.atoms[atom].atype != "H":
				near_O |= (residue.atoms[atom].atype == "O")
				near_N |= (residue.atoms[atom].atype == "N")
				near_MG |= (residue.atoms[atom].atype == "MG")
				near_water |= (residue.name == "HOH")
				near_cation |= (residue.name in cations)
				near_lysarg |= (residue.name in ("LYS", "ARG"))
				if residue.name == "HOH":
					for hetresidue, _, _ in get_nearby_atoms(2.58, residue.atoms[atom]):
						near_MG |= (hetresidue.name == "MG")
						near_cation |= (hetresidue.name in cations)
				if HBFinder._is_donor(residue, residue.atoms[atom])\
				  and residue.name in rings and atom not in rings[residue.name]:
					near_amino = True
		# 
		near = (near_MG, near_cation, near_N, near_amino, near_O, near_water, near_lysarg)
		unusual = (not sum(near, False))
		# 
		site = "%s.%s" % (site, a_oph.name)
		near = "\t".join(map(str, map(int, near)))
		print("% 9.3f\t% 9.3f\t% 9.3f\t%s\t%s" % (in_plane, off_plane, distance, near, site))


def detect_bptypes(files):
	for file in files:
		try:
			with SSManager(file, None) as manager:
				structure = manager.extract_paired_structure(breaks=False)
		except AttributeError:
			print("Error: %s" % file)
			continue
		basepairs = get_basepairs(structure)
		# print_triplets(basepairs, structure)
		uawht_analyse(basepairs, structure)
		del structure


def get_filepaths(folder, extension, ids=None):
	filepaths = list()
	if not isinstance(extension, str):
		for each_extension in extension:
			filepaths += get_filepaths(folder, each_extension, ids)
	else:
		len_ext = len(extension)
		for filename in sorted(os.listdir(folder)):
			if filename.endswith(extension)\
			  and (not ids or filename[:-len_ext] in ids):
				filepaths.append(os.path.join(folder, filename))
	return filepaths

def get_idslist(ids):
    if ids and os.path.exists(ids):
        with open(ids, 'r') as idsfile:
            ids = [line.strip() for line in idsfile.readlines()]
    return ids


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("folder", metavar="FOLDER")
    parser.add_argument("--ids", metavar="IDSLIST")
    # parser.add_argument("outfile", metavar="OUTPUT")
    return parser.parse_args()


if __name__ == "__main__":
	args = parse_args()
	idslist = get_idslist(args.ids)
	files = get_filepaths(args.folder, (".pdb", ".cif"), idslist)
	detect_bptypes(files)

# III {('G', 'G', 'W', 'H', 'C'): [123, 0], ('G', 'C', 'W', 'W', 'T'): [116, 12]}
# II  {('G', 'G', 'W', 'H', 'C'): [219, 68], ('G', 'C', 'W', 'W', 'T'): [86, 9]}