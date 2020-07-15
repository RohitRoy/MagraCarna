from __future__ import print_function

import os
import sys
import scipy
import argparse
import numpy as np

from decimal import Decimal
from itertools import count, combinations
from collections import defaultdict

from scipy.sparse import csgraph
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

from ..utilities.args import ArgParser
from ..magracarna.engraph import RNAGraph
from ..magracarna.motifs import MotifDivideFile, MotifRedivideFile,\
								MotifClusterFile, MotifGraphCode
from ..magracarna.metalsite import RNAMgSiteEngrapher, RNAMgSiteList,\
									RNAMgSiteContextVisualiser,\
									RNAMgMotifVectorList
from ..magracarna.vectorise import RNAVectorAnalyser


def merge(outfile, *infiles):
	divisionslist = list()
	allkeys = set()
	for file in infiles:
		divisions = MotifDivideFile.read_divisions(file).values()
		divisions = {tuple(key): map(tuple, value) for key, value in divisions}
		divisionslist.append(divisions)
		allkeys.update(divisions.keys())

	newdivisions = dict()
	by_counts = dict()
	for key in allkeys:
		newcaselist = set()
		for divisions in divisionslist:
			if key in divisions:
				newcaselist.update(divisions[key])
		newdivisions[key] = sorted(newcaselist)
		count_cases = len(set(zip(*newdivisions[key])[0]))
		by_counts[key] = (-count_cases, -len(newcaselist))

	by_counts = zip(*sorted(by_counts.items(), key=lambda entry: entry[1]))[0]
	for divno, key in enumerate(by_counts, 1):
		newdivisions[divno] = (key, newdivisions.pop(key))

	MotifDivideFile.write_divisions(outfile, newdivisions, 'w')

def limit(infile, outfile, ids=None):
	divisions = MotifDivideFile.read_divisions(infile)
	if ids:
		olddivisions = divisions.copy()
		idslist = ArgParser.read_idslist(ids)
		inlist = lambda entry: entry[0].split(" ")[0] in idslist
		for division in olddivisions:
			details, caselines = olddivisions[division]
			caselines = filter(inlist, caselines)
			if caselines:
				divisions[division] = (details, caselines)
			else:
				del divisions[division]
	MotifDivideFile.write_divisions(outfile, divisions)

def example(dividefile, resolutionsfile):
	examples = MotifDivideFile.example(dividefile, resolutionsfile)
	order = sorted(examples, key=lambda key: (key <= 0, abs(key)))
	for division in order:
		best, resolution = examples[division]
		print("% 3d\t%- 20s\t%.2f" % (division, best, resolution))

def disambiguate_cutoff(cutoff):
	kwargs = dict()
	try:
		cutoff = float(cutoff)
		if cutoff:
			if cutoff >= 1:
				kwargs["length"] = int(cutoff)
			else:
				kwargs["cutoff"] = float(cutoff)
	except ValueError:
		pass
	return kwargs

def common(redividefile, cutoff=1, show=False):
	kwargs = disambiguate_cutoff(cutoff)
	length = 0 if 0 < cutoff <= 1 else int()
	mergers = MotifRedivideFile.common_sequences(redividefile, **kwargs)
	for candidate, common_sequences in sorted(mergers.items()):
		print("%d\t%d" % candidate)
		if show:
			for sequence in common_sequences:
				print("\t%s" % sequence)

def copy(redividefile, clusterfile):
	MotifClusterFile.from_redividefile(redividefile, clusterfile)

def group(clusterfile, divisions):
	MotifClusterFile.cluster_divisions(clusterfile, [divisions])


def plot(dividefile, vectors, limit=4, show=False, save=None, label=None):
	divisions = MotifDivideFile.read_divisions(dividefile)
	cases = _cases_list(divisions.values())
	data = _obtain_data(cases, vectors)
	_, data, assigned = _data_assigned(divisions, cases, data)
	# 
	limit = int(limit) if int(limit) == limit else limit
	pca = PCA(n_components=limit)
	reduced = pca.fit(data).transform(data)
	components = pca.explained_variance_.shape[0]
	print(pca.explained_variance_ratio_)
	
	count = {divno: sum(1 for divno_ in assigned if divno_ == divno) for divno in set(assigned)}
	mask = np.array([count[divno] >= 2 for divno in assigned], dtype=bool)
	cases = [case for keep, case in zip(mask, cases) if keep]
	reduced = reduced[mask]
	assigned = assigned[mask]

	if isinstance(limit, float):
		print(components)
	if show:
		RNAVectorAnalyser.plot_reduced(cases, reduced, assigned, components, None, label)
	if save:
		RNAVectorAnalyser.plot_reduced(cases, reduced, assigned, components, save, label)


def neighbours(dividefile, vectors, show=False, n_components=4, k_neighbours=3):
	divisions = MotifDivideFile.read_divisions(dividefile)
	for division in sorted(divisions):
		if len(divisions[division][1]) <= k_neighbours:
			divisions.pop(division)
	cases = _cases_list(divisions.values())
	data = _obtain_data(cases, vectors)
	_, data, assigned = _data_assigned(divisions, cases, data)
	# 
	pca = PCA(n_components=n_components)
	reduced = pca.fit(data).transform(data)
	components = pca.explained_variance_.shape[0]
	if isinstance(n_components, float):
		print(components)
	if show:
		RNAVectorAnalyser.plot_reduced([], reduced, assigned, components, None)
	# 
	nn = NearestNeighbors(n_neighbors=k_neighbours).fit(reduced)
	adjacency = nn.kneighbors_graph(reduced).toarray()
	_, labels = csgraph.connected_components(adjacency)
	casewise = defaultdict(set)
	clusterwise = defaultdict(set)
	for case, label, division in zip(cases, labels, assigned):
		casewise[label].add(case)
		clusterwise[label].add(division)
	clusters = list()
	for clusterno, clusterdivs in clusterwise.items():
		clusters.append(clusterdivs)
		clusterdivs = " ; ".join(map(str, sorted(clusterdivs)))
	clusters = MotifClusterFile.simplify_clusters(clusters)
	for clusterno, clusterdivs in enumerate(clusters, 1):
		clusterdivs = " ; ".join(map(str, sorted(clusterdivs)))
		print("Cluster %d:\t%s" % (clusterno, clusterdivs))
	# for component, (_, cases) in enumerate(sorted(casewise.items()), 1):
	# 	print("\nComponent\t%d" % component)
	# 	for case in sorted(cases):
	# 		print("\t".join(case.split()))

def _data_assigned(divisions, cases, data, limit=None):
	if limit:
		newcases, newdata = list(), list()
		for case, vector in zip(cases, data):
			if case in limit:
				newcases.append(case)
				newdata.append(vector)
		cases, data = newcases, newdata
	assigned = dict()
	for divno, (_, divcases) in divisions.items():
		divcases = sum(divcases, [])
		for case in divcases:
			if case in cases:
				assigned[case] = divno
	for case in cases:
		if case not in assigned:
			assigned[case] = 0
	assigned = np.array([assigned[case] for case in cases])
	return cases, np.array(data, dtype=np.float64), np.array(assigned)

def _cases_list(division_values):
	caselines = sum(zip(*list(division_values))[1], [])
	return sorted(set(zip(*caselines)[0]))

def _obtain_data(cases, vectors):
	data = list()
	datadict = _data_dict(vectors)
	for case in cases:
		try:
			data.append(datadict[case])
		except KeyError:
			print("Missing %s" % case, sys.stderr)
	return data

def _data_dict(datafile):
	data = RNAMgMotifVectorList.read_dataframe(datafile)
	cases, data = RNAVectorAnalyser.prepare_data(data, angle=False)
	data = dict(zip(cases, map(list, data.to_records(index=False))))
	return data


# 

class MotifClusterMerger(object):

	ANSWERS = {"y": True, "n": False, "": False}

	def __init__(self, clusterfile, resolutionsfile,
		         structuredir, sitesfolder, redividefile=None):
		if redividefile:
			copy(redividefile, clusterfile)
		self.file = clusterfile
		self.resolutions = MotifClusterFile.read_resolutions(resolutionsfile)
		self.structuredir = structuredir
		self.sitesfolder = sitesfolder

	def structure(self, pdbid):
		exts = (".cif", ".pdb")
		paths = ArgParser.get_sources(self.structuredir, [pdbid], *exts)
		assert len(paths) == 1, "Multiple structures: %s" % ";".join(paths)
		return paths[0]

	def sites(self, pdbid):
		exts = ("_sites.tsv",)
		paths = ArgParser.get_sources(self.sitesfolder, [pdbid], *exts)
		assert len(paths) == 1, "Multiple sitess: %s" % ";".join(paths)
		return paths[0]

	def visualise(self, site):
		pdbid = site.split()[0]
		sites = self.sites(pdbid)
		structure = self.structure(pdbid)
		rnagfs = RNAMgSiteList.limited([sites], [site])
		visualiser = RNAMgSiteContextVisualiser(structure, rnagfs)
		site = site.split(" ")[1].split("]")[1]
		visualiser.display(site.strip("!"))

	def merge_by_sequences(self):
		kwargs = disambiguate_cutoff(raw_input("Cutoff "))
		divisions = MotifClusterFile.read_divisions(self.file)
		examples = MotifDivideFile._example(divisions, self.resolutions)
		mergers = MotifClusterFile.common_sequences(self.file, **kwargs)
		tomerge = list()
		for (div1, div2), common_sequences in sorted(mergers.items()):
			site1, site2 = examples[div1][0], examples[div2][0]
			print("Candidate:\n\t%d\t%s\n\t%d\t%s" % (div1, site1, div2, site2))
			if self.confirm("View"):
				self.visualise(site1)
				self.visualise(site2)
			if self.confirm("Merge"):
				tomerge.append((div1, div2))
		self.merge(tomerge)

	def merge_by_viewing(self):
		divisions = MotifClusterFile.read_divisions(self.file)
		examples = MotifClusterFile._example(divisions, self.resolutions)
		print("Current clusters %s" % " , ".join(map(str, sorted(divisions))))
		commands = "Commands:\n\tView 2 4...\n\tMerge 3 5...\n\tQuit\n\tUndo"
		tomerge = list()
		while True:
			answer = raw_input()
			try:
				if answer.startswith("View "):
					for divno in map(int, answer.split()[1:]):
						site = examples[divno][0]
						print("\t%s" % site)
						self.visualise(site)
				elif answer.startswith("Merge "):
					tomerge.append(map(int, answer.split()[1:]))
				elif answer.startswith("Undo") and len(tomerge):
					cluster = tomerge.pop(-1)
					print("%s merge undone" % "\t".join(map(str, cluster)))
				elif answer.startswith("Quit"):
					break
				else:
					print(commands)
			except ValueError:
				print(commands)
		self.merge(tomerge)

	def merge(self, clusters):
		MotifClusterFile.cluster_divisions(self.file, clusters)

	def interactive(self):
		while True:
			if self.confirm("Merge by seqeuence"):
				self.merge_by_sequences()
			if self.confirm("Merge by vectors PCA"):
				pass
			if self.confirm("Merge by viewing"):
				self.merge_by_viewing()
			if not self.confirm("\nContinue"):
				break


	def confirm(self, question):
		answer = None
		while answer not in self.ANSWERS:
			answer = raw_input("%s? [Y/n] " % question).lower()
		return self.ANSWERS[answer]


def interactive(**kwargs):
	clusterer = MotifClusterMerger(**kwargs)
	clusterer.interactive()


def parse_args():
	parser = argparse.ArgumentParser()
	commands = parser.add_subparsers(dest="command")
	# 
	merge = commands.add_parser('merge')
	merge.add_argument("outfile", metavar="MERGED-DIVIDE-FILE")
	merge.add_argument("infiles", nargs='+', metavar="INPUT-DIVIDE-FILE")
	# 
	limit = commands.add_parser('limit')
	limit.add_argument("infile", metavar="INPUT-DIVIDE-FILE")
	limit.add_argument("outfile", metavar="LIMITED-DIVIDE-FILE")
	limit.add_argument("--ids", default=None)
	# 
	example = commands.add_parser('example')
	example.add_argument("dividefile", metavar="DIVIDE-FILE")
	example.add_argument("resolutionsfile", metavar="RESOLUTIONS-FILE")
	# 
	common = commands.add_parser('common')
	common.add_argument("redividefile", metavar="DIVIDE-FILE")
	common.add_argument("--cutoff", "-c", type=float, default=1)
	common.add_argument("--length", "-l", type=int, default=0)
	common.add_argument("--show", "-s", action="store_true")
	# 
	copy = commands.add_parser('copy')
	copy.add_argument("redividefile", metavar="REDIVIDE-FILE")
	copy.add_argument("clusterfile", metavar="CLUSTER-FILE")
	#
	group = commands.add_parser('group')
	group.add_argument("clusterfile", metavar="CLUSTER-FILE")
	group.add_argument("divisions", nargs='+', type=int, metavar="DIVISION")
	#
	interactive = commands.add_parser('interactive')
	interactive.add_argument("--redividefile", metavar="REDIVIDE-FILE")
	interactive.add_argument("clusterfile", metavar="CLUSTER-FILE")
	interactive.add_argument("resolutionsfile", metavar="RESOLUTIONS-FILE")
	interactive.add_argument("structuredir", metavar="STRUCTURE-FOLDER")
	interactive.add_argument("sitesfolder", metavar="CONTEXT-FOLDER")
	# 
	LABEL = ("case", "division", "centre", "outlier")
	plot = commands.add_parser('plot')
	plot.add_argument("dividefile", metavar="DIVIDE-FILE")
	plot.add_argument("vectors", metavar="VECTORS")
	plot.add_argument("--show", "-s", action="store_true")
	plot.add_argument("--save", "-f", default=None, metavar="IMAGE-FILE")
	plot.add_argument("--limit", "-n", metavar="N", type=float, default=4)
	plot.add_argument("--label", "-l", choices=LABEL, default=None)
	# 
	neighbours = commands.add_parser('neighbours')
	neighbours.add_argument("dividefile", metavar="DIVIDE-FILE")
	neighbours.add_argument("vectors", metavar="VECTORS")
	neighbours.add_argument("--show", "-s", action="store_true")
	neighbours.add_argument("--k-neighbours", "-k", metavar="K", type=int)
	neighbours.add_argument("--n-components", "-n", metavar="N", type=float, default=4)
	# 
	return parser.parse_args()

COMMAND = {"merge": merge, "limit": limit, "example": example, "plot": plot,
			"common": common, "copy": copy, "group": group, 
			"interactive": interactive, "neighbours": neighbours}

if __name__ == "__main__":
	args = vars(parse_args())
	COMMAND[args.pop("command")](**args)
