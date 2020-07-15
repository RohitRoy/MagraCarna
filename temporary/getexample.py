import os
import sys

import time

from itertools import count
from collections import defaultdict

from ..utilities.args import ArgParser
from ..magracarna.metalsite import RNAMgSiteList
from ..magracarna.aminoacids import aminoacids
from ..magracarna.nucleotides import backbone

hetatoms = {"Nh", "Oh"}

def getexample(structuredir, idsfile):
	graphex = 0 # one-indexed
	lastlog = time.time()
	idslist = ArgParser.read_idslist(idsfile)
	filepaths = ArgParser.get_sources(structuredir, idslist, "_sites.tsv")
	for filepath in filepaths:
		pdbid = filepath.split("/")[-1].split("_")[0]
		for graph in RNAMgSiteList.load(filepath):
			graphex += 1
			if time.time() - lastlog > 60:
				print "% 10d\t%s" % (graphex, graph.name)
				lastlog = time.time()
			try:
				assert(not graph.name.endswith("!"))
				inner, outer, bpair = list(), list(), list()
				for edge in graph.edges:
					node0lab = graph.nodes[edge[0]].name
					node1lab = graph.nodes[edge[1]].name
					nonmg = node0lab if node0lab != "MG" else node1lab
					if edge[2].name == "L":
						inner.append(nonmg)
					elif edge[2].name == "LW1":
						outer.append(nonmg)
					elif edge[2].name == "LBP":
						bpair.append(nonmg)
				assert(len(inner) == 6)
				assert(0 < len(outer) <= 7)
				assert(0 < len(bpair) <= 3)
				assert(0 < sum((atom != "HOH" for atom in inner)))
				assert(0 < sum((atom != "HOH" for atom in outer)))
				assert(0 < sum((atom in hetatoms for atom in outer), 0))
				assert(sum((atom in hetatoms for atom in inner), 0) == 0)
			except AssertionError:
				continue
			else:
				print "Candidate:\t%s" % graph.name

if __name__ == "__main__":
	getexample(sys.argv[1], sys.argv[2])
