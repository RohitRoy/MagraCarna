import os
import sys
import time

from collections import defaultdict
from ...magracarna.nucleotides import nucleotides, allnucs
from ...magracarna.aminoacids import aminoacids
from ...utilities.args import ArgParser

def ligation_type(line):
	line = line.split()
	if "L" in line:
		return "  L"
	elif "LW1" in line:
		return "LW1"
	else:
		return ""

unknowns = set()

def residue_type(residue):
	# if residue in nucleotides:
	# 	return "std-RNA"
	if residue in allnucs:
		return "any-RNA"
	elif residue in aminoacids:
		return "protein"
	else:
		unknowns.add(residue)
		return "unknown"

def nonrna(sitesdir, idsfile):
	counts = defaultdict(lambda: 0)
	idslist = ArgParser.read_idslist(sys.argv[2])
	sitess = ArgParser.get_sources(sitesdir, idslist, "_sites.tsv")
	print("% 3d" % len(sitess))
	lastlog = time.time()
	for index, sitesfile in enumerate(sitess, 1):
		if time.time() - lastlog > 15:
			print("% 3d" % index)
			lastlog = time.time()
		with open(sitesfile, 'r') as infile:
			ligands = None
			siteid = None
			for line in infile:
				if line.startswith("ContextGraph:") and "!" not in line:
					ligands = set()
					siteid = " ".join(line.split()[1:3])
				elif line.strip() and ligands is not None:
					ligation = ligation_type(line)
					if ligation:
						residue = line.split()[1].split("]")[0].strip("[")
						if residue != "HOH":
							ligands.add((ligation, residue_type(residue)))
				elif ligands:
					counts[tuple(sorted(ligands))] += 1
					ligands = None
			if ligands:
				counts[tuple(sorted(ligands))] += 1
	summation = 0
	for key in sorted(counts, key=len):
		print("% 5d\t%s" % (counts[key], "\t".join(map(" ".join, key))))
		summation += counts[key] if key != "*" else 0
	print("% 5d\n"% summation)
	print("%s\n" % "\t".join(sorted(unknowns)))

	summary = defaultdict(lambda: 0)
	for key in counts:
		for subkey in key:
			if "RNA" not in subkey[1]:
				break
		else:
			summary["only any-RNA"] += counts[key]
			for subkey in key:
				if "  L" == subkey[0]:
					summary["only L-any-RNA"] += counts[key]
					break
			else:
				summary["only LW1-any-RNA"] += counts[key]
		for subkey in key:
			if "unknown" in subkey[1]:
				summary["unknowns"] += counts[key]
				break
		for subkey in key:
			if "protein" in subkey[1]:
				summary["protein + RNA"] += counts[key]
				break
	for key in summary:
		print("% 5d\t%s" % (summary[key], key))

def main():
	nonrna(sys.argv[1], sys.argv[2])

main()