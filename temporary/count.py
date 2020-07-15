import os
import sys

from numpy import argsort
from collections import defaultdict
from ...magracarna.engrid import StructureFile
from ...magracarna.motifs import MotifDivideFile
from ...magracarna.engraph import RNAGraph
from ...utilities.args import ArgParser

PAIR_INDICES = (4, 6, 1, 2, 5)
RNANUCS = ("A", "C", "G", "U")
ATOMTYPES = [RNAGraph.Oph, RNAGraph.Or, RNAGraph.b]

def count_divides(dividefile):
	def as_pairs(details):
		base_pairs = list()
		for divdetails in details:
			pair = list()
			for index in PAIR_INDICES:
				pair.append(divdetails[index].split()[-1])
			pair[-1] = pair[-1][-1]
			base_pairs.append(tuple(pair))
		return base_pairs
	pairings = MotifDivideFile.read_divisions(dividefile).values()
	details, caselines = zip(*pairings)
	base_pairs = as_pairs(details)
	cases = [set(zip(*divcaselines)[0]) for divcaselines in caselines]
	withoutatom = defaultdict(set)
	for pair, caseset in zip(base_pairs, cases):
		withoutatom[pair].update(caseset)
	for pair in sorted(withoutatom, key=lambda value: -len(withoutatom[value])):
		print("%s:%s %s:%s%s\t% 5d" % (pair + (len(withoutatom[pair]),)))

def count_residues(structuredir):
	print("    \t%s" % "\t".join(RNANUCS))
	for filepath in ArgParser.get_sources(structuredir, None, ".cif", ".pdb"):
		count = {nuc: 0 for nuc in RNANUCS}
		with StructureFile(filepath) as manager:
			for residue in manager.extract_residues(only_residues=RNANUCS):
				count[residue.name] += 1
			row = [manager.structid] + [count[nuc] for nuc in RNANUCS]
		print("\t".join(map(str, row)))

def count_hbonds_divides(dividefile):
	bytype = defaultdict(set)
	hbondings = MotifDivideFile.read_divisions(dividefile).values()
	for detail, divcaselines in hbondings:
		divcases = set(zip(*divcaselines)[0])
		detail = [line.split()[-1] for line in detail]
		if detail.pop() not in RNANUCS or detail[2] not in RNANUCS\
		  or detail[1] not in RNAGraph.b.eqvals\
		  or (detail[3] not in RNAGraph.b.eqvals and detail[3] != "O2'"):
			continue
		for atomtype in ATOMTYPES:
			if detail[-1] in atomtype.eqvals:
				detail[-1] = atomtype.name
				break
		bondtype = "-".join(detail)
		bytype[bondtype].update(divcases)
	for bondtype in bytype:
		bytype[bondtype] = len(bytype[bondtype])
	for bondtype, count in sorted(bytype.items(), key=lambda item: -item[1]):
		print("% 20s\t% 5s" % (bondtype, count))

def count_multiplet_divides(*dividefiles):
	counts = list()
	for dividefile in dividefiles:
		bytype = defaultdict(set)
		multiplets = MotifDivideFile.read_divisions(dividefile).values()
		for detail, divcaselines in multiplets:
			divcases = set(zip(*divcaselines)[0])
			detail = [line for line in detail if not line.startswith("1 b is ")]
			detail = "-".join([line.split()[-1] for line in detail][1:])
			bytype[detail].update(divcases)
		bytype["*"] = set.union(*(bytype.values() + [set()]))
		bytype = {key: len(cases) for key, cases in bytype.items()}
		counts.append(defaultdict(lambda: 0, bytype))
	multiplets = set.union(*map(set, counts))
	multiplets = [[-bytype[key] for bytype in counts]+[key] for key in multiplets]
	for multiplet in ["*"]: # sorted(multiplets):
		row = multiplet[-1].split("-")
		row += ["% 3d" % bytype[multiplet[-1]] for bytype in counts[::-1]]
		print("\t".join(row))

def count_multiplet_format(countsfile):
	with open(countsfile, 'r') as infile:
		header = None
		order = None
		cases = list()
		for line in infile:
			if line.startswith("!"):
				continue
			elif line.startswith("Multiplet:"):
				header = line.split("\t")[1].strip()
				print("\n%s" % header)
			elif line.startswith("Format:"):
				order = map(int, map(str.strip, line.split("\t")[1:]))
				order = argsort(order) + 1
			elif line.strip() and header and format:
				line = map(str.strip, line.split("\t"))
				for index, value in enumerate(line):
					if value in RNANUCS:
						line[index] = "*"+value
						break
				key = [line[index] for index in order]
				key = [value[-1] if ":" in value else value for value in key]
				for index, value in enumerate(key):
					if value.strip("*") not in RNANUCS:
						first = index
						break
				else:
					first = 100
				pairs = zip(key[first::3], key[first+1::3], key[first+2::3])
				key = " ".join(["".join(key[:first])] + map("".join, pairs))
				values = "\t".join(line[len(order)+1:])
				print("\t%s\t%s" % (key, values))

def main():
	if sys.argv[1] == "residue":
		if sys.argv[2] == "structure":
			count_residues(sys.argv[3])
	elif sys.argv[1] == "pair":
		if sys.argv[2] == "divides":
			count_divides(sys.argv[3])
	elif sys.argv[1] == "hbond":
		if sys.argv[2] == "divides":
			count_hbonds_divides(sys.argv[3])
	elif sys.argv[1] == "multiplet":
		if sys.argv[2] == "divides":
			count_multiplet_divides(*sys.argv[3:])
		elif sys.argv[2] == "format":
			count_multiplet_format(sys.argv[3])

main()