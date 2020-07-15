import os
import sys

from collections import defaultdict

import numpy as np

from ..utilities.args import ArgParser
from ..magracarna.engrid import StructureFile
from ..magracarna.motifs import MotifDivideFile, MotifRedivider
from ..magracarna.heteroatoms import cations as CATIONS


RESIDUES = {"MG", "HOH"}
def count(idsfile, structuredir):
	for pdbid in ArgParser.read_idslist(idsfile):
		count = {residue: 0 for residue in RESIDUES}
		path = ArgParser.get_sources(structuredir, [pdbid], ".pdb", ".cif")[0]
		with StructureFile(path) as structure:
			for residue in structure.extract_residues(only_residues=RESIDUES):
				count[residue.name] += 1
		print("%s\t% 5d\t% 5d" % (pdbid, count["MG"], count["HOH"]))

def check_rRNA(chainsfile):
	molecules = zip(*MotifRedivider.read_chains_file(chainsfile).values())[0]
	molecules = map(str.split, molecules)
	molecules = [molecule for molecule in molecules if "rRNA" in molecule]
	return "other" if not len(molecules) else "rRNA"

def read_counts(countsfile):
	counts = dict()
	with open(countsfile, 'r') as infile:
		for line in infile:
			pdbid, mgcount, hohcount = line.split("\t")
			counts[pdbid] = (int(mgcount), int(hohcount))
	return counts

def read_sites(sitesfile):
	sites = dict()
	with open(sitesfile, 'r') as infile:
		for line in infile:
			pdbid, allsites, onlyvalids = line.split("\t")
			sites[pdbid] = (int(allsites), int(onlyvalids))
	return sites

def sites(idsfile, sitesdir):
	idslist = ArgParser.read_idslist(idsfile)
	sitesfiles = ArgParser.get_sources(sitesdir, idslist, "_sites.tsv")
	for pdbid, sitesfile in zip(idslist, sitesfiles):
		count = 0
		valids = 0
		with open(sitesfile, 'r') as infile:
			for line in infile:
				if line.startswith("ContextGraph:"):
					count += 1
					line = line.strip()
					if not line.endswith("!"):
						valids += 1
		print("%s\t%s\t%s" % (pdbid, count, valids))

def table(idsfile, countsfile, sitesfile, resolutionsfile, chainsdir):
	idslist = ArgParser.read_idslist(idsfile)
	counts = read_counts(countsfile)
	sites = read_sites(sitesfile)
	resolutions = MotifDivideFile.read_resolutions(resolutionsfile)
	resolutions = [resolutions.get(pdbid, None) for pdbid in idslist]
	chainsfiles = ArgParser.get_sources(chainsdir, idslist, "_entities.txt")
	for pdbid, rvalue, chainsfile in zip(idslist, resolutions, chainsfiles):
		hasrRNA = check_rRNA(chainsfile)
		line = (pdbid, rvalue, hasrRNA) + sites[pdbid] + counts[pdbid]
		print("\t".join(map(str, line)))

def read_table(idsfile, tablefile):
	data = list()
	idslist = ArgParser.read_idslist(idsfile)
	with open(tablefile, 'r') as infile:
		for line in map(str.split, infile):
			if line[0] in idslist:
				pdbid, rvalue, rRNA, sites, valids, mgcount, hohcount = line
				hohtomg = round(float(hohcount) / float(mgcount), 2)
				sites, valids, rvalue = int(sites), float(valids), float(rvalue)
				data.append((sites, valids, rvalue, hohtomg, rRNA == "rRNA"))
	return np.array(data)

def summary(idsfile, tablefile):
	data = read_table(idsfile, tablefile)
	all_ = float(data[:, 0].sum())
	rRNA = float(np.array(data[data[:, 4] == 1][:, 0], dtype=int).sum())
	other = all_ - rRNA
	# 
	valids = float(data[:, 1].sum())
	valid_rRNA = float(np.array(data[data[:, 4] == 1][:, 1], dtype=int).sum())
	valid_other = valids - valid_rRNA
	# 
	frac_valid = round(valids / all_, 2) * 100
	frac_valid_rRNA = round(valid_rRNA / rRNA, 2) * 100
	frac_valid_other = round(valid_other / other, 2) * 100
	print("% 2.1f\t% 2.1f\t% 2.1f" % (frac_valid, frac_valid_rRNA, frac_valid_other))
	#  
	bad_rRNA = data[(data[:, 3] <= 3) & (data[:, 4] == 1)][:, 0].sum()
	bad_other = data[(data[:, 3] <= 3) & (data[:, 4] == 0)][:, 0].sum()
	frac_rRNA = round(rRNA / all_, 3) * 100
	frac_bad_rRNA = round(bad_rRNA / rRNA, 3) * 100
	frac_bad_other = round(bad_other / other, 3) * 100
	print("% 2.1f\t% 2.1f\t% 2.1f" % (frac_rRNA, frac_bad_rRNA, frac_bad_other))
	print("")
	# 

def wbyr(idsfile, tablefile):
	output = list()
	data = read_table(idsfile, tablefile)
	print("   \t  \t   \t%s" % "\t".join(map(lambda x: "%4d" % x, range(8))))
	resgroups = map(lambda x: x/10, map(float, range(15, 36, 2)))
	for minres, maxres in zip([0]+resgroups, resgroups):
		row = list()
		subdata = data[(data[:, 2] > minres) & (data[:, 2] <= maxres)]
		for minW, maxW in enumerate(range(7) + [float("inf")], -1):
			count = subdata[(subdata[:, 3] > minW) & (subdata[:, 3] <= maxW)]
			row.append(int(count[:, 0].sum()))
		print("%.1f\tto\t%.1f\t%s" % (minres, maxres, "\t".join(map(str, row))))
		output.append(row)
	print("")
	output = np.array(output)
	netsites = float(output.sum())
	for upto in range(0, 4):
		frac = sum([output[:, j].sum() for j in range(upto+1)]) / netsites
		print("%d\t% 2.1f" % (upto, round(frac * 100, 1)))
	print("")

def print_ligands(pdbid, isvalid, ligands):
	if isvalid is None or ligands is None:
		return
	rna, nonrna, water = 0, 0, 0
	for ligand in ligands:
		if ligand in ("Oh", "Nh"):
			nonrna += 1
		elif ligand == "HOH":
			water += 1
		else:
			rna += 1
	print("%s\t%d\t%d\t%d\t%d" % (pdbid, isvalid, rna, nonrna, water))

def ligands(idsfile, sitesdir):
	idslist = ArgParser.read_idslist(idsfile)
	sitesfiles = ArgParser.get_sources(sitesdir, idslist, "_sites.tsv")
	for pdbid, sitesfile in zip(idslist, sitesfiles):
		isvalid, ligands = None, None
		with open(sitesfile, 'r') as infile:
			for line in infile:
				line = line.strip().split()
				if not line:
					continue
				elif line[0] == "ContextGraph:":
					print_ligands(pdbid, isvalid, ligands)
					isvalid = int(not line[2].endswith("!"))
					ligands = list()
				elif "L" in line:
					ligands.append(line[2])
			print_ligands(pdbid, isvalid, ligands)

def rnaligands(idsfile, ligandsfile):
	idslist = ArgParser.read_idslist(idsfile)
	with open(ligandsfile, 'r') as infile:
		data = [line.split() for line in infile]
		data = [map(int, line[1:]) for line in data if line[0] in idslist]
	data = np.array(data)
	valid = data[data[:, 0] == 1]
	ligands = "\t".join(["%1.2f" % valid[:, i].mean().round(2) for i in range(1, 4)])
	average_coo = (valid[:, 1] + valid[:, 2] +  valid[:, 3]).mean().round(2)
	print("%s\t%1.2f" % (ligands, average_coo))
	#
	invalid = data[data[:, 0] == 0]
	ligands = "\t".join(["%1.2f" % invalid[:, i].mean().round(2) for i in range(1, 4)])
	average_coo = (invalid[:, 1] + invalid[:, 2] +  invalid[:, 3]).mean().round(2)
	print("%s\t%1.2f" % (ligands, average_coo))
	# 
	ligands = "\t".join(["%1.2f" % data[:, i].mean().round(2) for i in range(1, 4)])
	average_coo = (data[:, 1] + data[:, 2] +  data[:, 3]).mean().round(2)
	print("%s\t%1.2f" % (ligands, average_coo))

def read_coordination(idsfile, ligandsfile):
	idslist = ArgParser.read_idslist(idsfile)
	with open(ligandsfile, 'r') as infile:
		data = [line for line in map(str.split, infile) if line[0] in idslist]
		data = [(line[0], sum(map(int, line[2:]))) for line in data]
	data = [case + (data.count(case),) for case in set(data)]
	coordination = defaultdict(lambda: defaultdict(lambda: 0))
	for pdbid, coono, count in data:
		if coono > 6:
			coordination[pdbid][str("> 6")] += count
		else:
			coordination[pdbid][str(coono)] += count
	return coordination

def wbyc(idsfile, tablefile, ligandsfile):
	idslist = ArgParser.read_idslist(idsfile)
	ratio = dict()
	with open(tablefile, 'r') as infile:
		for line in map(str.split, infile):
			if line[0] in idslist:
				pdbid, _, _, _, _, mgcount, hohcount = line
				ratio[pdbid] = round(float(hohcount) / float(mgcount), 3)
	# 
	coordination = read_coordination(idsfile, ligandsfile)
	# 
	row = ["% 5s" % coono for coono in map(str, range(7))+["> 6"]]
	print("ratio\t%s" % "\t".join(row))
	#
	for minW, maxW in enumerate(range(7) + [float("inf")], -1):
		structures = [pdbid for pdbid in ratio if minW < ratio[pdbid] <= maxW]
		row = defaultdict(lambda: 0)
		for pdbid in structures:
			for coono in coordination[pdbid]:
				row[coono] += coordination[pdbid][coono]
		row = ["% 5d" % row[coono] for coono in map(str, range(7))+["> 6"]]
		print("%s\t%s" % (maxW, "\t".join(row)))

def main():
	if sys.argv[1] == "count":
		count(sys.argv[2], sys.argv[3])
	if sys.argv[1] == "sites":
		sites(sys.argv[2], sys.argv[3])
	elif sys.argv[1] == "table":
		table(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
	elif sys.argv[1] == "summary":
		summary(sys.argv[2], sys.argv[3])
	elif sys.argv[1] == "wbyr":
		wbyr(sys.argv[2], sys.argv[3])
	elif sys.argv[1] == "ligands":
		ligands(sys.argv[2], sys.argv[3])
	elif sys.argv[1] == "rnaligands":
		rnaligands(sys.argv[2], sys.argv[3])
	elif sys.argv[1] == "wbyc":
		wbyc(sys.argv[2], sys.argv[3], sys.argv[4])

if __name__=="__main__":
	main()
