
from __future__ import print_function

import os
import sys
import math

import numpy as np

from itertools import combinations, product
from collections import defaultdict
from scipy.spatial.distance import pdist, squareform

from ..bin.homology import MgRNAProfile, ClusterFile, MoleculeDataset
from ..magracarna.engrid import StructureFile


def getfile(code, filetype):
	return os.path.join("Output/Homology/Molecules/"+code+"/"+code+filetype)

def radius():
	codes = ["05DR", "16TT"]
	savefile = "savefile.txt"
	for code in codes:
		sitexwise = dict()
		clusterfile = ClusterFile(getfile(code, ".clusters.txt"))
		for cluster in clusterfile.read_cluster_names():
			for sitex in zip(*clusterfile.read_cluster_sites(cluster))[0]:
				sitexwise[sitex] = cluster
		# 
		with StructureFile(getfile(code, ".pdb"), None) as manager:
			structure = manager.extract_structure(breaks=False, multimodel=True)
		begex, endex, _ = structure.chains["A"]
		# 
		coos = defaultdict(list)
		for residue in structure.residues[endex:]:
			sitex = residue.resno + int(residue.chain)*10000
			coos[sitexwise[sitex]].append(residue.atoms.values()[0].coos)
		del structure
		# 
		radii = list()
		for cluster in set(coos) - {0}:
			clustercoos = np.array(coos[cluster])
			# radii.append(max(pdist(clustercoos)) / 2.0)
			# 
			centroid = clustercoos.sum(axis=0) / clustercoos.shape[0]
			# 
			# distances = squareform(pdist(clustercoos))
			# indices = np.unravel_index(np.argmax(distances), distances.shape)
			# centroid = (clustercoos[indices[0]] + clustercoos[indices[1]]) / 2
			# distance = np.round(distances[indices] / 2.0, 3)
			# # 
			# farones = np.linalg.norm(clustercoos-centroid, axis=1)
			# farones = clustercoos[np.round(farones, 3) >= distance]
			# centroid = farones.sum(axis=0) / farones.shape[0]
			# if farones.shape[0] < 2:
			# 	print("%s\t%d" % (code, cluster))
			# 	continue
			radii.append(max(np.linalg.norm(clustercoos-centroid, axis=1)))
		radii = "\t".join(["% 7s" % ("%.3f" % radius) for radius in radii])
		with open(savefile, 'a') as outfile:
			outfile.write("%s\t%s\n" % (code, radii))

if __name__ == "__main__":
	radius()
