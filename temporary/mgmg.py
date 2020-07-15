import os
import sys
import time

from itertools import combinations
from collections import defaultdict

import numpy as np

from ...utilities.args import ArgParser
from ...magracarna.engrid import StructureFile
from ...magracarna.motifs import MotifDivideFile, MotifRedivider
from ...magracarna.heteroatoms import cations as CATIONS


RESIDUES = {"MG"}
def mindistance(idsfile, structuredir):
	minpair = ""
	lasttime = 0
	mindistance = float("inf")
	for count, pdbid in enumerate(ArgParser.read_idslist(idsfile), 1):
		if time.time() - lasttime > 60:
			print("% 5d %1.3f %s" % (count, mindistance, minpair))
			lasttime = time.time()
		path = ArgParser.get_sources(structuredir, [pdbid], ".pdb", ".cif")[0]
		with StructureFile(path) as structure:
			mgresidues = structure.extract_residues(only_residues=RESIDUES)
			for mg1, mg2 in combinations(mgresidues, 2):
				distance = mg1.atoms["MG"].distance_to(mg2.atoms["MG"])
				if distance < mindistance and distance != 0:
					mindistance = distance
					minpair = "%s %s %s" % (pdbid, mg1, mg2)
	print("% 5d %1.3f %s" % (count, mindistance, minpair))

if __name__=="__main__":
	mindistance(sys.argv[1], sys.argv[2])
