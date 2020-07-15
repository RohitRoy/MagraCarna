from __future__ import print_function

import os
import sys
import math

import numpy as np

from itertools import combinations, product
from collections import defaultdict

from ..bin.homology import MgRNAProfile, ClusterFile, MoleculeDataset
from ..magracarna.engrid import StructureFile
from ..magracarna.nucleotides import backbone, phosphate
from ..magracarna.heteroatoms import cations

backbone = set(backbone) | set(phosphate)
backboneO = {atom for atom in backbone if "O" in atom}

RHO = 90

class ClusterAnalyser(object):

	RELIABILITY = ("N", "U", "R", "B", "M", "H")

	def __init__(self, profilefile, clusterfile, chainscount):
		with StructureFile(profilefile, None) as manager:
			self.structure = manager.extract_structure(breaks=True, multimodel=True)
		self.clusterfile = ClusterFile(clusterfile)
		self.chainscount = chainscount
		self.molecule = os.path.basename(profilefile).split(".")[0]
		self.clusters = dict()
		self.category = dict()
		self.centroids = None
		self._metalcoos = dict()
		self._chainresidue = dict()
		# 
		self._min10 = max(int(math.floor(chainscount * 0.1)), 2)
		self._min50 = max(int(math.floor(chainscount * 0.5)), 2)
		self._min80 = max(int(math.floor(chainscount * 0.8)), 2)
		# 
		for cluster in self.clusterfile.read_cluster_names() - {0}:
			sitexes = zip(*self.clusterfile.read_cluster_sites(cluster))[0]
			self.clusters[cluster] = sitexes
		self.__get_centroids()

	def chainresidue(self, locus):
		if locus not in self._chainresidue:
			residue = self.structure.find(resno=locus, chain="A")
			self._chainresidue[locus] = residue
		return self._chainresidue[locus]

	def metalcoos(self, sitex):
		if sitex not in self._metalcoos:
			resno, chain = sitex % 10000, str(sitex / 10000)
			residue = self.structure.find(resno=resno, chain=chain)
			metalatom = list(residue.atoms.values())[0]
			self._metalcoos[sitex] = metalatom.coos
		return self._metalcoos[sitex]

	def pop(self, cluster):
		for dictionary in [self.clusters, self.centroids, self.category]:
			if dictionary and cluster in dictionary:
				dictionary.pop(cluster)

	def restrict_to(self, reliability):
		acceptable = self.RELIABILITY[self.RELIABILITY.index(reliability):]
		for cluster in self.clusterfile.read_cluster_names() - {0}:
			category = self.categorise(cluster)
			if category not in acceptable:
				self.pop(cluster)

	def categorise(self, cluster):
		sitexes = set()
		reliables = 0
		for sitex, siteid in self.clusterfile.read_cluster_sites(cluster):
			sitexes.add(sitex)
			residue = siteid.split()[1]
			if not residue.startswith(("[NA]", "[MG]")):
				isreliable = True
			elif residue.startswith("[NA]") or siteid.endswith(("!", "*")):
				isreliable = False
			elif residue.startswith("[MG]") and siteid.endswith(" "):
				isreliable = True
			reliables += int(isreliable)
		self.category[cluster] = self._categorise(len(sitexes), reliables)
		return self.category[cluster]

	def _categorise(self, sites, reliables):
		if sites < self._min10:
			return "N" # not conserved
		if reliables < 2:
			return "U" # unreliably conserved
		if reliables >= self._min50 and sites >= self._min80:
			return "H"	# highly conserved
		if reliables >= self._min10 and sites >= self._min50:
			return "M" # moderately well conserved
		if reliables >= self._min10 or sites >= self._min50:
			return "B" # borderline well conserved
		return "R" # rarely conserved

	def __get_centroids(self):
		sitexwise = defaultdict(lambda: 0)
		self.centroids = dict()
		for cluster, sitexes in self.clusters.items():
			sitexwise.update({sitex: cluster for sitex in sitexes})
			self.centroids[cluster] = np.zeros(3)
		# 
		_, endex, _ = self.structure.chains["A"]
		for residue in self.structure.residues[endex:]:
			if residue.chain != "A":
				sitex = int(residue.chain)*10000 + residue.resno
				if sitexwise[sitex]:
					centroid = self.centroids[sitexwise[sitex]]
					centroid += residue.atoms.values()[0].coos
		for cluster in self.centroids:
			self.centroids[cluster] /= len(self.clusters[cluster])


def bends(profilefile, clusterfile, chainscount, reliability, verbose=True):
	analyser = ClusterAnalyser(profilefile, clusterfile, chainscount)
	molecule = os.path.basename(profilefile).split(".")[0]
	structure = analyser.structure
	analyser.restrict_to(reliability)
	# 
	bends = set()
	begex, endex, breaks = analyser.structure.chains["A"]
	breakexes = [begex] + (list(zip(*breaks)[0]) if breaks else []) + [endex]
	for begex, endex in zip(breakexes, breakexes[1:]):
		window = list()
		for residue in analyser.structure.residues[begex:endex]:
			window.append(residue)
			window = window[-3:]
			for index in range(0, len(window)-2):
				try:
					endP = window[-1].atoms["P"].coos
					startP = window[index].atoms["P"].coos
					postend = endP - window[-2].atoms["P"].coos
					prestart = window[index+1].atoms["P"].coos - startP
				except KeyError:
					continue
				postend /= np.linalg.norm(postend)
				prestart /= np.linalg.norm(prestart)
				angle = np.clip(np.dot(prestart, postend), -1.0, +1.0)
				angle = 180.0 - np.degrees(np.arccos(angle))
				if angle <= RHO:
					bends.add(tuple(window[index:])+(angle,))
	# 
	bends = sorted(bends, key=lambda window: window[0].resno)
	segments = list()
	if bends:
		segments.append(list(bends.pop(0)))
		selector = list()
		for window in map(list, bends):
			if set(window[:-1]) & set(segments[-1][:-1]):
				rhoangle = min(segments[-1].pop(-1), window.pop(-1))
				segments[-1] += window + [rhoangle]
			else:
				segments.append(window)
	if verbose:
		for segment in segments:
			rhoangle = "%.1f" % segment[-1]
			begno, endno = segment[0].resno, segment[-2].resno
			print("% 5d-% -5d\t% 5s" % (begno, endno, rhoangle))
	return analyser, segments


def phosphateO(residue, residue5end=None):
	atomcoos = list()
	if residue5end:
		try:
			atomcoos.append(residue5end.atoms["O3'"].coos)
		except KeyError:
			pass
	for atomname in ["OP1", "OP2", "O5'"]:
		try:
			atomcoos.append(residue.atoms[atomname].coos)
		except KeyError:
			continue
	return atomcoos


def kinktype(cutoff, centroids, res1, res2, res3):
	atoms5P = phosphateO(res1)
	atomsCP = phosphateO(res2, res1)
	atoms3P = phosphateO(res3, res2)
	# 
	cluster5P, cluster3P, clusterCP = set(), set(), set()
	for cluster, centroid in centroids.items():
		for coos in atoms5P:
			if np.linalg.norm(centroid - coos) <= cutoff:
				cluster5P.add(cluster)
		for coos in atoms3P:
			if np.linalg.norm(centroid - coos) <= cutoff:
				cluster3P.add(cluster)
		if cluster in (cluster5P & cluster3P):
			break
		for coos in atomsCP:
			if np.linalg.norm(centroid - coos) <= cutoff:
				clusterCP.add(cluster)
	# 
	kinktype = dict()
	kinktype["any"] = bool(cluster5P) or bool(clusterCP) or bool(cluster3P)
	kinktype["one"] = bool(cluster5P) or bool(cluster3P)
	kinktype["two"] = bool(cluster5P) and bool(cluster3P)
	kinktype["both"] = bool((cluster5P & cluster3P))
	return kinktype


def longkinktype(cutoff, centroids, segment):
	atoms5P = phosphateO(segment[0])
	atoms3P = phosphateO(segment[-1], segment[-2])
	for cluster, centroid in centroids.items():
		for coos in atoms5P:
			if np.linalg.norm(centroid - coos) <= cutoff:
				for coos in atoms3P:
					if np.linalg.norm(centroid - coos) <= cutoff:
						return True
	return False


def cutoff(profilefile, clusterfile, chainscount, cutoff, reliability):
	FORMAT = "% 6s\t% 12s\t% 1d\t% 1d\t% 1d\t% 1d\t% 3.1f\t% 1d\t% 5s"
	# molecule code, segment name,
	# whether ion interacts at any, one-side, two-sides, both-sides of kink,
	# (sharpest) rhoangle, kink-length, kink P-P distance
	# 
	analyser, segments = bends(profilefile, clusterfile, chainscount, reliability, False)
	for segment in segments:
		rhoangle = segment.pop(-1)
		contactat = {"any": False, "one": False, "two": False, "both": False}
		for res123 in zip(segment[0::3], segment[1::3], segment[2::3]):
			for item in kinktype(cutoff, analyser.centroids, *res123).items():
				contactat[item[0]] |= item[1]
			if contactat["both"]:
				break
		if len(segment) > 3 and not contactat["both"]:
			straddles = longkinktype(cutoff, analyser.centroids, segment)
			contactat["two"] |= straddles
			contactat["both"] |= straddles
		# 
		vector = segment[0].atoms["P"].coos - segment[-1].atoms["P"].coos
		segname = "%s-%s:A" % (segment[0].resno, segment[-1].resno)
		args = [analyser.molecule, segname]
		args += [contactat[atype] for atype in ("any", "one", "two", "both")]
		args += [rhoangle, len(set(segment)), "%.1f" % np.linalg.norm(vector)]
		print(FORMAT % tuple(args))


def closest(profilefile, clusterfile, chainscount, reliability):
	FORMAT = "% 6s\t% 12s\t% 3.1f\t% 1d\t% 5s\t% 5s"
	# molecule code, segment name, (sharpest) rhoangle, kink-length,
	# kink P-P distance, distance to closest ion cluster centroid
	# 
	analyser, segments = bends(profilefile, clusterfile, chainscount, reliability, False)
	for segment in segments:
		rhoangle = segment.pop(-1)
		residues = sorted(set(segment), key=lambda residue: residue.resno)
		residue5end = None
		dist = float("inf")
		for residue in residues:
			for atomcoos in phosphateO(residue, residue5end):
				for cluster, centroid in analyser.centroids.items():
					dist = min(np.linalg.norm(centroid - atomcoos), dist)
			residue5end = residue
		# 
		vector = segment[0].atoms["P"].coos - segment[-1].atoms["P"].coos
		segname = "%s-%s:A" % (segment[0].resno, segment[-1].resno)
		args = [analyser.molecule, segname, rhoangle, len(residues)]
		args += ["%.1f" % value for value in (np.linalg.norm(vector), dist)]
		print(FORMAT % tuple(args))


def control(profilefile, clusterfile, chainscount, cutoff, reliability):
	FORMAT = "% 6s\t% 5d / % 5d\t% 5d / % 5d\t% 5d / % 5d"
	# molecule, number of 3-length segments with:
	# no kinking residues and how many contact a cluster,
	# some kinking residue and how many contact a cluster,
	# only kinking residues and how many contact a cluster.
	analyser, segments = bends(profilefile, clusterfile, chainscount, reliability, False)
	kinkresidues = set().union(*[segment[:-1] for segment in segments])
	# 
	totals = {"no": 0, "all": 0, "part": 0}
	counts = {"no": 0, "all": 0, "part": 0}
	begex, endex, breaks = analyser.structure.chains["A"]
	breakexes = [begex] + (list(zip(*breaks)[0]) if breaks else []) + [endex]
	for begex, endex in zip(breakexes, breakexes[1:]):
		window, contacts = list(), list()
		for residue in analyser.structure.residues[begex:endex]:
			window.append(residue)
			window = window[-3:]
			contacts = contacts[-2:]
			residue5end = None if len(window) < 2 else window[-1]
			for coos in phosphateO(residue, window[-1]):
				for cluster, centroid in analyser.centroids.items():
					if np.linalg.norm(centroid - coos) <= cutoff:
						contacts.append(True)
						break
				if len(contacts) == len(window):
					break
			if len(contacts) != len(window):
				contacts.append(False)
			if len(window) == 3:
				kinks = len(set(window) & kinkresidues)
				thetype = "all" if kinks == 3 else "part" if kinks else "no"
				totals[thetype] += 1
				counts[thetype] += int(all(contacts))
	# 
	args = [analyser.molecule]
	ORDER = ("no", "part", "all")
	args += sum([[counts[atype], totals[atype]] for atype in ORDER], [])
	print(FORMAT % tuple(args))


def types(profilefile, clusterfile, chainscount):
	analyser = ClusterAnalyser(profilefile, clusterfile, chainscount)
	category = defaultdict(list)
	sitexes = defaultdict(list)
	for cluster in sorted(set(analyser.clusters) - {0}):
		code = analyser.categorise(cluster)[0]
		sitexes[code] += list(analyser.clusters[cluster])
		category[code].append(cluster)
	category["N"].append(0)
	sitexes["N"] += list(zip(*analyser.clusterfile.read_cluster_sites(0))[0])
	# 
	for code in ClusterAnalyser.RELIABILITY:
		print("% 2s\t% 4d\t% 7d" % (code, len(category[code]), len(sitexes[code])))


def check(kinkstable, msafile, structuredir):
	molecule = os.path.basename(msafile).split(".")[0]
	# 
	kinks = list()
	with open(kinkstable, "r") as infile:
		for line in infile:
			if line.strip().startswith(molecule):
				angle = line.split()[6]
				begex, endex = line.split()[1].split(":")[0].split("-")
				kinks.append(tuple(map(int, (begex, endex))) + (angle,))
	counts = {kink: 0 for kink in kinks}
	allkinks = len(kinks)
	# 
	dataset = MoleculeDataset("*", "*", None, msafile=msafile)
	total = len(dataset.chains)
	iterator = dataset.iterate(count=True, structuredir=structuredir, align=True)
	for pdbid, chainid, chainex, structure, alignment in iterator:
		found = 0
		for kink in kinks:
			begex, endex, _ = kink
			residues = [alignment[index] for index in range(begex-1, endex)]
			while len(residues) >= 3:
				bend, residues = residues[:3], residues[1:]
				try:
					atom5P, atomCP, atom3P = [each.atoms["P"] for each in bend]
				except (KeyError, AttributeError):
					continue
				if atomCP.angle_between(atom5P, atom3P) <= RHO:
					counts[kink] += 1
					found += 1
					break
		print("% 3d / % 3d\t%s\t%s\t% 3d / % 3d" % (chainex, total, pdbid, chainid, found, allkinks))
	# 
	missed = 0
	badmissed = 0
	print(molecule, file=sys.stderr)
	for kink in kinks:
		if counts[kink] < total:
			frac = counts[kink] / float(total)
			isbad = "Almost" if frac > 0.95 else "Missed"
			missed += 1
			if frac <= 0.95:
				badmissed += 1
			# 
			strfrac = "%.1f" % (100 * frac)
			args = (isbad, counts[kink], total, strfrac, kink[0], kink[1], kink[2])
			print("\t\t%s %d / %d % 5s %d-%d\t% 5s" % args, file=sys.stderr)
	print("\tMissed %d\t%d\t%d" % (missed, badmissed, allkinks), file=sys.stderr)


if __name__ == "__main__":
	if sys.argv[1] == "bends":
		bends(sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5])
	elif sys.argv[1] == "closest":
		closest(sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5])
	elif sys.argv[1] == "cutoff":
		cutoff(sys.argv[2], sys.argv[3], int(sys.argv[4]), float(sys.argv[5]), sys.argv[6])
	elif sys.argv[1] == "control":
		control(sys.argv[2], sys.argv[3], int(sys.argv[4]), float(sys.argv[5]), sys.argv[6])
	elif sys.argv[1] == "types":
		types(sys.argv[2], sys.argv[3], int(sys.argv[4]))
	elif sys.argv[1] == "check":
		check(sys.argv[2], sys.argv[3], sys.argv[4])

