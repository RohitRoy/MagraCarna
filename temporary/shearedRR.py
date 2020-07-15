import os
import sys
from itertools import permutations, product

from ...utilities.args import ArgParser
from ...magracarna.nucleotides import edges, phosphate
from ...magracarna.engrid import HydrogenBondFinder, SecondaryStructureManager

hbfinder = HydrogenBondFinder()

FOLDER = sys.argv[1]
IDS = ArgParser.read_idslist(sys.argv[2])

REQUIRED = product(["A", "G"], ["A", "G"])
REQUIRED = product(REQUIRED, [("S", "H", "T"), ("s", "h", "T")])
REQUIRED = {purines+bptype for purines, bptype in REQUIRED}
REQUIRED = {bptype: [0, 0, 0, 0, 0] for bptype in REQUIRED}

SecondaryStructureManager.WARN = False

def get_basepairs(structure):
	basepairs = list()
	required_residues = {bptype[0] for bptype in REQUIRED}
	for residue in structure.pairs:
		if residue.name in required_residues:
			itspairs = structure.pairs[residue].items()
			for pairresidue, bptype in itspairs:
				bptype = (residue.name, pairresidue.name)+bptype
				if bptype in REQUIRED:
					yield (residue, pairresidue, bptype)

def filtered_basepairs(structure, hbfinder, conditions=(True, True)):
	basepairs = list()
	for residue, pairresidue, bptype in get_basepairs(structure):
		donors = set()
		bondresidues = set()
		for atomname in edges[residue.name]['w']:
			atom = residue.atoms[atomname]
			for bondresidue, bondatom in hbfinder.find_hbonds(atom, residue):
				condition1 = bondatom.name in phosphate
				condition2 = structure.linked(bondresidue, pairresidue)
				if (condition1, condition2) == conditions:
					donors.add(atomname)
					bondresidues.add(bondresidue)
		if donors:
			yield (residue, pairresidue, bptype, donors)


def yield_near(residue, atoms, checker, checkfor, within):
	for atomname in sorted(atoms):
		atom = residue.atoms[atomname]
		for nearresidue, _, _ in checker(within, atom):
			if nearresidue.name == checkfor:
				yield nearresidue


def check_near(residue, atoms, checker, checkfor, within):
	for _ in yield_near(residue, atoms, checker, checkfor, within):
		return True
	return False

def check_hoogsteen_mg(sresidue, checker):
	hoogsteen = set(edges[sresidue.name]["H"])
	nearmg = check_near(sresidue, hoogsteen, checker, "MG", 7.0)
	mgbound = check_near(sresidue, hoogsteen & {"O6"}, checker, "MG", 3.08)
	mgbound |= check_near(sresidue, hoogsteen & {"N7"}, checker, "MG", 3.2)
	nearhoh, nearmghoh = False, False
	for hohresidue in yield_near(sresidue, hoogsteen, checker, "HOH", 3.8):
		nearhoh = True
		if check_near(hohresidue, ["O"], checker, "MG", 3.08) and not mgbound:
			nearmghoh = True
	return list(map(int, [True, nearmg, mgbound, nearmghoh, nearhoh]))


try:
	print len(ArgParser.get_sources(FOLDER, IDS, ".pdb", ".cif"))
	for filepath in ArgParser.get_sources(FOLDER, IDS, ".pdb", ".cif"):
		with SecondaryStructureManager(filepath, None) as manager:
			structure = manager.extract_paired_structure()
		hbfinder.set_structure(structure)
		get_nearby_atoms = structure.grid.get_nearby_atoms
		print(structure.structid)
		for sresidue, hresidue, bptype, watoms in filtered_basepairs(structure, hbfinder):
			results = check_hoogsteen_mg(sresidue, get_nearby_atoms)
			for index in xrange(5):
			    REQUIRED[bptype][index] += results[index]
			# if check_near(sresidue, edges[sresidue.name]["H"], get_nearby_atoms, "MG", 7.0):
			# 	REQUIRED[bptype][1] += 1
except KeyboardInterrupt:
	pass

print("% 10s\t% 8s\t %8s\t %8s\t% 8s\t% 8s" % ("BPTYPE", "TOTAL", "NEAR", "L", "LW1", "HOH"))
for bptype in sorted(REQUIRED):
	if REQUIRED[bptype][0] != 0:
		counts = "\t".join(["% 8d" % count for count in REQUIRED[bptype]])
		printargs = (" ".join(bptype),)+tuple(REQUIRED[bptype])
		print("%s\t%s" % (" ".join(bptype), counts))
