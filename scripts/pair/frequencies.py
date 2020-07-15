from __future__ import print_function

import os
import sys
import time
import argparse

from collections import defaultdict

from ...magracarna.engrid import Structure, HydrogenBondFinder as HBFinder,\
								 SecondaryStructureManager as SSManager
from ...utilities.args import ArgParser


SSManager.WARN = False
SSManager.NEW_BPFIND = False

def get_basepairs(structure):
	basepairs = list()
	for residue in structure.pairs:
		itspairs = structure.pairs[residue].items()
		for pairresidue, bptype in itspairs:
			bptype = (residue.name, pairresidue.name) + bptype
			basepairs.append((residue, pairresidue, bptype))
	return basepairs


def update_frequencies(basepairs, frequencies):
	if not frequencies:
		frequencies = defaultdict(lambda: 0)
	for base1, base2, bptype in basepairs:
		if bptype[0] != bptype[1] or str(base1) <= str(base2):
			frequencies[bptype] += 1
	return frequencies

def print_frequencies(frequencies):
	frequencies = sorted(frequencies.items(), key=lambda entry: -entry[1])
	for bptype, count in frequencies:
		bptype = " ".join([part.rjust(3) for part in bptype])
		print("%s\t% 6d" % (bptype, count))


def detect_bptypes(files):
	totalfiles = len(files)
	timelogged = time.time()
	frequencies = None
	print(len(files))
	for fileno, file in enumerate(files, 1):
		if time.time() - timelogged > 50:
			print("%d / %d" % (fileno, totalfiles), file=sys.stderr)
			timelogged = time.time()
		print("entered")
		try:
			with SSManager(file, None) as manager:
				structure = manager.extract_paired_structure(breaks=False)
			print("constructed")
		except AttributeError:
			print("Error: %s" % file, file=sys.stderr)
			continue
		basepairs = get_basepairs(structure)
		print("got pairs")
		frequencies = update_frequencies(basepairs, frequencies)
		print("updated")
	print_frequencies(frequencies)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("folder", metavar="FOLDER")
    parser.add_argument("--ids", metavar="IDSLIST")
    # parser.add_argument("outfile", metavar="OUTPUT")
    return parser.parse_args()


if __name__ == "__main__":
	args = parse_args()
	idslist = ArgParser.read_idslist(args.ids)
	files = ArgParser.get_sources(args.folder, idslist, ".pdb", ".cif")
	detect_bptypes(files)

# III {('G', 'G', 'W', 'H', 'C'): [123, 0], ('G', 'C', 'W', 'W', 'T'): [116, 12]}
# II  {('G', 'G', 'W', 'H', 'C'): [219, 68], ('G', 'C', 'W', 'W', 'T'): [86, 9]}