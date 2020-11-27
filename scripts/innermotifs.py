from __future__ import print_function

import os
import sys
import argparse

from math import floor
from decimal import Decimal
from numpy.linalg import norm
from collections import defaultdict
from numpy import array, flatnonzero

from ..bin.homology import LigandsFile, MoleculeDataset, ClusterFile
from ..magracarna.engrid import Residue, HydrogenBondFinder
from ..magracarna.metalsite import RNAMgSiteEngrapher
from ..magracarna.nucleotides import allnucs as NUCLEOTIDES, purines as PURINES

DONORS = HydrogenBondFinder.DONORS
OPH = ("OP1", "OP2", "OP3")
UOB = ("O4", "O2")

def read_ligands(filepath):
	ligands = defaultdict(dict)
	for _, siteid, ligandslist in LigandsFile.read(filepath):
		if not ligandslist:
			continue
		siteid = siteid.strip()
		pdbid, residueid = siteid.split(" ")
		chainid = residueid.split(":")[1].strip("!*")
		ligandslist = [ligand.strip().split(".") for ligand in ligandslist]
		ligandslist = [(int(locus), atom) for locus, atom in ligandslist]
		ligands[(pdbid, chainid)][siteid] = ligandslist
	return dict(ligands)


def read_clusters(filepath):
	clusterfile = ClusterFile(filepath)
	clusternames = clusterfile.read_cluster_names()
	clusterof = dict()
	sitesof = defaultdict(list)
	for cluster in clusternames:
		for _, siteid in clusterfile.read_cluster_sites(cluster):
			clusterof[siteid.strip()] = cluster
			sitesof[cluster].append(siteid)
	return clusterof, dict(sitesof)


def by_chain(pdbid, chainid, structure, aligned, ligands, extra):
	subset = set()
	if (pdbid, chainid) in ligands:
		subset.update(ligands[(pdbid, chainid)])
	get_nearby_atoms = structure.grid.get_nearby_atoms
	for pdbid_, chainid_ in extra:
		if pdbid_ != pdbid:
			continue
		for siteid, ligandslist in ligands[(pdbid, chainid_)].items():
			locus, atom = ligandslist[0]
			if not aligned[locus-1] or atom not in aligned[locus-1].atoms:
				continue
			residueid = siteid.strip("!* ").split(" ")[1]
			nearby = get_nearby_atoms(7.0, aligned[locus-1].atoms[atom])
			for residue, _, _ in nearby:
				if str(residue) == residueid:
					subset.add(siteid)
	return subset


def get_shells(detector, residue, structure):
	detector.set_residue(residue, structure)
	detector.detect_environment()
	detector.detect_inner_shell(engraph=False)
	if len(detector.inner_shell) < detector.MAX_COORDINATION:
		detector._relaxed_inner_shell()
	detector.detect_outer_shell(engraph=False)
	#
	inner = sorted(detector.inner_shell, key=lambda entry: entry[2])
	inner = [(entry[0], entry[1].name) for entry in detector.inner_shell]
	outer = [(entry[0], entry[1].name) for entry in detector.outer_shell]
	inner = [ligand for ligand in inner if ligand[0].name != "HOH"]
	outer = [ligand for ligand in outer if ligand[0].name != "HOH"]
	return inner, outer


def assign_motifs(inner, outer, structure, aligned):
	CN = len(inner)
	if CN == 0 or CN > 3:
		return None
	#
	nonrna = [_ for residue, _ in inner if residue.name not in NUCLEOTIDES]
	if len(nonrna) > 0:
		return None
	#
	inner = sorted(inner, key=lambda entry: entry[1])
	nonOph = [(residue, atom) for residue, atom in inner if atom[:2] != "OP"]
	Oph = CN - len(nonOph)
	#
	if CN == 1 and inner[0][1] == "N7" and inner[0][0].name == "G":
		Pout = {residue for residue, atom in outer if atom[:2] == "OP"}
		Hout = [_ for residue, _ in outer if residue.name not in NUCLEOTIDES]
		if not len(Hout) and len(Pout) == 1 and Pout.pop() == inner[0][0]:
			Glocus = aligned.index(inner[0][0])+1
			return ("G-N7 MacroChelate I", (Glocus,))
	#
	if CN == 2 and Oph == 1:
		Presidue = [residue for residue, atom in inner if atom[:2] == "OP"][0]
		if nonOph[0][0].name == "G" and nonOph[0][1] in ("O6",):
			Glocus = aligned.index(nonOph[0][0])+1
			Plocus = aligned.index(Presidue)+1
			return ("G-Phosphate", (Glocus, Plocus))
		elif nonOph[0][0].name == "U" and nonOph[0][1] in ("O4", "O2"):
			Ulocus = aligned.index(nonOph[0][0])+1
			Plocus = aligned.index(Presidue)+1
			return ("U-Phosphate", (Ulocus, Plocus))
	#
	if CN == 2 and Oph == 0 and (inner[0][1], inner[1][1]) == ("N7", "N7")\
	  and inner[0][0].name in PURINES and inner[1][0].name in PURINES:
		R1locus = aligned.index(inner[0][0])+1
		R2locus = aligned.index(inner[1][0])+1
		R1locus, R2locus = sorted([R1locus, R2locus])
		return ("Purine-N7 Seat", (R1locus, R2locus))
	#
	if CN == 3 and Oph == 2:
		if nonOph[0][1] == "N7" and nonOph[0][0].name in PURINES:
			Oph = [residue for residue, atom in inner if atom[:2] == "OP"]
			Ophlocus = [aligned.index(residue)+1 for residue in Oph]
			first, second = sorted(zip(Ophlocus, Oph))
			Flocus, first = first
			Slocus, second = second
			if not structure.linked(first, second):
				return None
			if structure.linked(second, nonOph[0][0]):
				Rlocus = aligned.index(nonOph[0][0])+1
				return ("10-Member Ring Purine-N7", (Flocus, Slocus, Rlocus))
	#
	return None


def well_conserved(cluster_sites, allchains):
	cluster_size = len(cluster_sites)
	if cluster_size < max(int(floor(allchains * 0.50)), 2):
		return False
	reliable_size = 0
	for siteid in cluster_sites:
		residue = siteid.split()[1]
		if not residue.startswith(("[NA]", "[MG]")):
			reliable_size += 1
		elif residue.startswith("[MG]") and not siteid.endswith(("!", "*")):
			reliable_size += 1
	if reliable_size < max(int(floor(allchains * 0.10)), 2):
		return False
	return True


def motif_counts(dataset, structuredir, ligands, clusterof):
	allchains = len(dataset.chains)
	extra = set(ligands.keys()) - set(dataset.chains)
	hbfinder = HydrogenBondFinder(max_da=Decimal("3.50"))
	detector = RNAMgSiteEngrapher(hbfinder, relaxed=True)
	#
	motifs = defaultdict(lambda: defaultdict(lambda: 0))
	dataiter = dataset.iterate(True, structuredir, align=True)
	for pdbid, chain, count, structure, aligned in dataiter:
		message = "% 4d / % 4d    \t%s:%s" % (count, allchains, pdbid, chain)
		print(message, file=sys.stderr)
		subset = by_chain(pdbid, chain, structure, aligned, ligands, extra)
		for siteid in sorted(subset):
			if siteid.endswith(("!", "*")):
				continue
			residueid = siteid.split(" ")[1].strip("!* ")
			if not residueid.startswith("[MG]"):
				continue
			residue = structure.find(**Residue.from_string(residueid, True))
			inner, outer = get_shells(detector, residue, structure)
			try:
				motif = assign_motifs(inner, outer, structure, aligned)
			except ValueError:
				continue
			if motif:
				motif, loci = motif
				motifs[motif][(loci, clusterof[siteid])] += 1
	return motifs


def conserved_motif_clusters(motifcounts, sitesof, allchains):
	conserved = defaultdict(set)
	for motif in motifcounts:
		for (loci, cluster), count in motifcounts[motif].items():
			if cluster != 0 and count >= 2 \
			  and well_conserved(sitesof[cluster], allchains):
				conserved[motif].add((loci, cluster))
	for motif in sorted(conserved):
		print("\t%s" % motif)
		for loci, cluster in sorted(conserved[motif]):
			print("\t\t%s\t%s" % (cluster, loci))
	return conserved


def distance_to(metal, residue, atom):
	try:
		return round(norm(metal.coos - residue.atoms[atom].coos), 2)
	except KeyError:
		return float("inf")

def get_vector(motif, residues, metal, structure):
	siteid = structure.find()
	if motif == "G-N7 MacroChelate I":
		toN7 = distance_to(metal, residues[0], "N7")
		toOph = min([distance_to(metal, residues[0], atom) for atom in OPH])
		return (toN7, toOph)
	elif motif in ("G-Phosphate", "U-Phosphate"):
		Bresidue, Presidue = residues
		toOph = min([distance_to(metal, Presidue, atom) for atom in OPH])
		if motif == "G-Phosphate":
			toO6 = distance_to(metal, Bresidue, "O6")
			return (toO6, toOph)
		if motif == "U-Phosphate":
			toOb = min([distance_to(metal, Bresidue, atom) for atom in UOB])
			return (toOb, toOph)
	elif motif == "Purine-N7 Seat":
		toN71 = distance_to(metal, residues[0], "N7")
		toN72 = distance_to(metal, residues[1], "N7")
		return (toN71, toN72)
	elif motif == "10-Member Ring Purine-N7":
		P1residue, P2residue, N7residue = residues
		toP1 = min([distance_to(metal, P1residue, atom) for atom in OPH])
		toP2 = min([distance_to(metal, P2residue, atom) for atom in OPH])
		toN7 = distance_to(metal, N7residue, "N7")
		return (toP1, toP2, toN7)


def getdistances(dataset, structuredir, conserved, sitesof):
	distances = defaultdict(lambda: defaultdict(list))
	allchains = len(dataset.chains)
	dataiter = dataset.iterate(True, structuredir, align=True)
	for pdbid, chain, count, structure, aligned in dataiter:
		message = "% 4d / % 4d    \t%s:%s" % (count, allchains, pdbid, chain)
		print(message, file=sys.stderr)
		for motif in conserved:
			for loci, cluster in sorted(conserved[motif]):
				residues = [aligned[locus-1] for locus in loci]
				if None in residues:
					continue
				for siteid in sitesof[cluster]:
					if not siteid.startswith(pdbid):
						continue
					metal = siteid.split(" ")[1].strip("!*")
					chain_ = metal.split(":")[1]
					if (pdbid, chain_) in dataset.chains and chain_ != chain:
						continue
					metal = structure.find(**Residue.from_string(metal, True))
					for atom in metal.atoms.values():
						if atom.atype not in ("O", "N"):
							break
					vector = get_vector(motif, residues, atom, structure)
					if float("inf") not in vector:
						if any([distance <= 7.0 for distance in vector]):
							distances[motif][(loci, cluster)].append(vector)
	return distances

def print_distances(distances):
	for motif in sorted(distances):
		print("\t%s" % motif)
		for (loci, cluster) in sorted(distances[motif]):
			header = [cluster] + list(map(str, loci))
			header = ['% 5s' % part for part in header]
			print("\t\t%s" % "\t".join(header))
			for line in distances[motif][(loci, cluster)]:
				line = list(map(lambda part: "% 5s" % ("%.2f" % part), line))
				print("\t\t\t%s" % "\t".join(line))


def getmotifs(structuredir, msafile, ligandsfile, clusterfile):
	ligands = read_ligands(ligandsfile)
	dataset = MoleculeDataset("*", "*", msafile=msafile)
	allchains = len(dataset.chains)
	clusterof, sitesof = read_clusters(clusterfile)
	motifcounts = motif_counts(dataset, structuredir, ligands, clusterof)
	conserved = conserved_motif_clusters(motifcounts, sitesof, allchains)
	distances = getdistances(dataset, structuredir, conserved, sitesof)
	print_distances(distances)


def read_conserved_motifs(conservedmotifs):
	byMotif, motif = defaultdict(list), None
	with open(conservedmotifs, 'r') as infile:
		for line in infile:
			if not line.strip():
				continue
			elif "(" in line:
				cluster, loci = line.strip().split("\t")
				loci = [each for each in loci.strip("()").split(",") if each]
				loci = tuple(map(int, loci))
				byMotif[motif].append((int(cluster), loci))
			elif line.startswith("\t"):
				motif = line.strip()
	return dict(byMotif)


def bin_motif_distances(conservedmotifs, motif, column):
	BINS = (1.5, 1.9, 2.3, 2.4, 2.6, 3.2, 3.8, 4.6, 5.0, 6.0, 7.0, 8.0)
	values = defaultdict(list)
	molnow, ismotif, locusnow = None, None, None
	with open(conservedmotifs, 'r') as infile:
		for line in infile:
			if not line.strip():
				continue
			if not line.startswith("\t"):
				molnow = line.strip()
			if line.startswith("\t") and not line.startswith("\t\t"):
				ismotif = (line.strip() == motif)
			elif line.startswith("\t\t") and not line.startswith("\t\t\t"):
				if ismotif and not "(" in line:
					locusnow = (molnow,) + tuple(map(int, line.split())[1:])
			elif line.startswith("\t\t\t") and ismotif:
				values[locusnow].append(float(line.split()[column - 1]))
	for akey in sorted(values):
		string = [" ".join(["% 6s" % value for value in map(str, akey)]+[":"])]
		row_ = array(values[akey])
		for min_, max_ in zip(BINS, BINS[1:]):
			string.append(len(flatnonzero((row_ >= min_) & (row_ < max_))))
		print("\t".join(map(str, string)))


def variation(motif, loci, aligned, inner, outer):
	classes = set()
	CN = len(inner)
	Ophres = [ligand[0] for ligand in inner if ligand[1][:2] == "OP"]
	Pout = {ligand[0] for ligand in outer if ligand[1][:2] == "OP"}
	if motif == "G-N7 MacroChelate I":
		Hout = [_ for residue, _ in outer if residue.name not in NUCLEOTIDES]
		if CN > 1 or len(Pout) > 1 or len(Hout):
			classes.add("E")
		if (aligned[loci[0]-1], "N7") not in inner:
			classes.add("O" if (aligned[loci[0]-1], "N7") in outer else "M")
			if CN > 0:
				classes.add("E")
		if aligned[loci[0]-1] not in Pout:
			classes.add("I" if aligned[loci[0]-1] in Ophres else "M")
			if len(Pout) > 0:
				classes.add("E")
		return classes
	#
	innermissed = 0
	if motif == "Purine-N7 Seat":
		for locus in loci:
			if (aligned[locus-1], "N7") not in inner:
				classes.add("O" if (aligned[locus-1], "N7") in outer else "M")
				innermissed += 1
	elif motif in ("G-Phosphate", "U-Phosphate"):
		if aligned[loci[1]-1] not in Ophres:
			classes.add("O" if aligned[loci[1]-1] in Pout else "M")
			innermissed += 1
		if motif[0] == "G" and (aligned[loci[0]-1], "O6") not in inner:
			classes.add("O" if (aligned[loci[0]-1], "O6") in outer else "M")
			innermissed += 1
		if motif[0] == "U":
			if (aligned[loci[0]-1], "O4") not in inner\
			  and (aligned[loci[0]-1], "O2") not in inner:
				in_outer = (aligned[loci[0]-1], "O4") in outer
				in_outer |= (aligned[loci[0]-1], "O2") in outer
				classes.add("O" if in_outer else "M")
				innermissed += 1
	elif motif == "10-Member Ring Adjacent":
		for locus in loci[:2]:
			if aligned[locus-1] not in Ophres:
				classes.add("O" if aligned[locus-1] in Pout else "M")
				innermissed += 1
		if (aligned[locus[2]-1], "N7") not in inner:
			classes.add("O" if (aligned[locus[2]-1], "N7") in outer else "M")
			innermissed += 1
	#
	print("\t%s\t%d\t%d" % (innermissed, len(inner), len(outer)))
	try:
		print("\t%s" % str(zip(map(str, zip(*inner)[0]), zip(*inner)[1])))
	except:
		pass
	try:
		print("\t%s" % str(zip(map(str, zip(*outer)[0]), zip(*outer)[1])))
	except:
		pass
	if CN > (len(loci)-innermissed):
		classes.add("E")
	return classes


def variability(structuredir, msafile, clusterfile, conservedmotifs):
	byMotif = read_conserved_motifs(conservedmotifs)
	_, sitesof = read_clusters(clusterfile)
	allclusters = sum([list(zip(*each)[0]) for each in byMotif.values()], [])
	allpdbids = set()
	for cluster in allclusters:
		allpdbids.update((siteid.split()[0] for siteid in sitesof[cluster]))
	#
	dataset = MoleculeDataset("*", "*", msafile=msafile)
	dataiter = dataset.iterate(True, structuredir, align=True)
	allchains = len(dataset.chains)
	hbfinder = HydrogenBondFinder(max_da=Decimal("3.50"))
	detector = RNAMgSiteEngrapher(hbfinder, relaxed=True)
	#
	classes = defaultdict(list)
	for pdbid, chain, count, structure, aligned in dataiter:
		message = "% 4d / % 4d    \t%s:%s" % (count, allchains, pdbid, chain)
		print(message, file=sys.stderr)
		for motif in byMotif:
			for cluster, loci in sorted(byMotif[motif]):
				for siteid in sitesof[cluster][:]:
					if not siteid.startswith(pdbid):
						continue
					metal = siteid.split(" ")[1].strip("!*")
					chain_ = metal.split(":")[1]
					if (pdbid, chain_) in dataset.chains and chain_ != chain:
						continue
					if (pdbid, chain_) not in dataset.chains:
						cation = Residue.from_string(metal, True)
						cation = structure.find(**cation)
						nearby = structure.grid.get_nearby_atoms(7.0, cation)
						nearby = set(zip(*nearby)[0])
						expected = {aligned[locus-1] for locus in loci}
						expected.discard(None)
						if not len(expected.intersection(nearby)):
							continue
					sitesof[cluster].remove(siteid)
					if not metal.startswith("[MG]"):
						classes[(motif, cluster, loci)].append((siteid, "N"))
						continue
					mgion = structure.find(**Residue.from_string(metal, True))
					inner, outer = get_shells(detector, mgion, structure)
					print(siteid)
					assigned = variation(motif, loci, aligned, inner, outer)
					assigned = "".join(sorted(assigned)) if assigned else "X"
					classes[(motif, cluster, loci)].append((siteid, assigned))
	#
	for cluster in allclusters:
		if len(sitesof[cluster]):
			print("Missed in %d:\t" % cluster)
			for siteid in sitesof[cluster]:
				print("\t%s" % siteid)
	for motif, cluster, loci in sorted(classes):
		print("%s\t%s\t%s" % (motif, cluster, loci))
		for siteid, classname in sorted(classes[motif, cluster, loci]):
			print("\t% 5s\t%s" % (classname, siteid))


def parse_args():
	parser = argparse.ArgumentParser()
	commands = parser.add_subparsers(dest="command")
	#
	getmotifs = commands.add_parser('getmotifs')
	getmotifs.add_argument("structuredir", metavar="STRUCTURE-FOLDER")
	getmotifs.add_argument("msafile", metavar="MSA-OUTPUT-FILE")
	getmotifs.add_argument("ligandsfile", metavar="LIGANDS-FILE")
	getmotifs.add_argument("clusterfile", metavar="CLUSTER-FILE")
	#
	variability = commands.add_parser('variability')
	variability.add_argument("structuredir", metavar="STRUCTURE-FOLDER")
	variability.add_argument("msafile", metavar="MSA-OUTPUT-FILE")
	variability.add_argument("clusterfile", metavar="CLUSTER-FILE")
	variability.add_argument("conservedmotifs", metavar="CONSERVED-MOTIFS")
	#
	bin = commands.add_parser('bin')
	bin.add_argument("conservedmotifs", metavar="CONSERVED-MOTIFS")
	bin.add_argument("motif", metavar="MOTIF-NAME")
	bin.add_argument("column", metavar="COLUMN-NUMBER", type=int)
	return parser.parse_args()


COMMANDS = {"getmotifs": getmotifs, "variability": variability,
			"bin": bin_motif_distances}


if __name__ == "__main__":
	args = vars(parse_args())
	COMMANDS[args.pop("command")](**args)
