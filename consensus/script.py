
import os
import sys
import time
import argparse

import numpy as np

from itertools import count
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA

from ..utilities.args import ArgParser
from ..magracarna.motifs import CaseReader, MotifDivideFile
from ..magracarna.engraph import DFSCode, RNAGraph
from ..magracarna.metalsite import RNAMgSiteList, RNAMgMotifVectorList
from ..magracarna.consensus import ConsensusFinder
from ..magracarna.vectorise import Vectoriser, RNAVectorAnalyser
# from ..magracarna.vectorise import Vectoriser


def contexts(infolder, dividefile, division, outfile):
	cases = CaseReader.read_limited(limited=dividefile, division=division)
	pdbids = tuple(sorted({case.split(" ")[0] for case in cases}))
	contexts = ArgParser.get_sources(infolder, pdbids, "_context.tsv")
	RNAMgSiteList.limited(contexts, cases, dumpfile=outfile)

def merge(sources, merged_file, general_file):
	ConsensusFinder.merge_and_modify(sources, merged_file, general_file)


class ConsensusVectoriser(Vectoriser):
	"""
	"""

	def __init__(self, motifgraph, basegraph):
		self.base = basegraph
		self._code2base = None
		super(ConsensusVectoriser, self).__init__(motifgraph)

	def _setup(self):
		self._code2base = self._map_node_tagstrs()
		super(ConsensusVectoriser, self)._setup()

	def _set_tags(self):
		for isnode, tagstr, tagtuple in self.code._iter_tagstr():
			self.tags.append(tagstr)
			if tagstr in self._code2base:
				nodetag = tagtuple[0]
				label = self.base._tags[nodetag]
				for typelabel, tagtype in self.NODELABEL_TYPE.items():
					if typelabel.covers(label):
						self.types.append(tagtype)
						break
			if len(self.tags) != len(self.types):
				empty_type = self.NODE if isnode else self.EDGE
				self.types.append(empty_type)

	def _set_nucleosides(self):
		codegf = RNAGraph.from_code(self.base)
		for index, tagstr, tagtype in zip(count(), self.tags, self.types[:]):
			if tagtype == self.NUCLEOTIDE and tagstr in self._code2base:
				nucid = int(tagstr.split(" ")[0])
				if self._has_phosphate_node(nucid, codegf):
					self.types[index] = self.NUCLEOSIDE
		for tagex, tagstr in enumerate(self.tags):
			if tagstr in self._code2base:
				self.tags[tagex] = self._code2base[tagstr]

	def _map_node_tagstrs(self):
		""" Mapping of node tagstrs from code to base.
		"""
		self._check_base()
		tagstrs = [entry[1] for entry in self.base._iter_tagstr() if entry[0]]
		tagstr_map = dict()
		for isnode, tagstr, _ in self.code._iter_tagstr():
			if isnode:
				tagstr_map[tagstr] = tagstrs.pop(0)
			if not tagstrs:
				break
		return tagstr_map

	def _check_base(self):
		embeds = RNAGraph.from_code(self.code).embed_code(self.base)
		for embed in embeds:
			consensus_tags, base_tags = list(zip(*enumerate(embed.nodeids)))
			if consensus_tags == base_tags:
				return embed
		raise ValueError("Base not placed correctly in consensus graph")


class ConsensusVectorList(RNAMgMotifVectorList):

	def __init__(self, motifgraph, basegraph, assigner, **kwargs):
		self.base = basegraph
		superinit = super(ConsensusVectorList, self).__init__
		superinit(motifgraph, assigner, **kwargs)

	def _set_vectoriser(self):
		base = self.assigner.graphs[self.base]
		code = self.assigner.graphs[self.graph]
		self._vectoriser = ConsensusVectoriser(code, base)

	def dump_vectors(self, outfilepath):
		header = self._vectoriser.header
		args = (outfilepath, self.base, header, self.cases, self.vectors)
		self.write_vectors(*args)


def vectorise(name, base, structuredir, assignsdir, outfile, motifs):
	Vectoriser.EDGELABEL_TYPE = {}
	assigner = ArgParser.get_assigner(motifs)
	vectorlist = ConsensusVectorList(name, base, assigner, mode="min")
	structures = ArgParser.get_sources(structuredir, None, ".pdb", ".cif")
	assigns = ArgParser.get_sources(assignsdir, None, "_assigns.tsv")
	sources = ArgParser.map_sources(structures, assigns)
	formatting = "%% 5d / %d\t%%s" % len(sources)
	lastlog = time.time()
	for index, (structurefile, assignsfile) in enumerate(sources, 1):
		if (time.time() - lastlog) > 60:
			print(formatting % (index, ArgParser.basename(structurefile, ".")))
			lastlog = time.time()
		vectorlist.add_vectors_from(assignsfile, structurefile)
	vectorlist.dump_vectors(outfile)


def finalise(outfile, default, clusters):
	read = RNAMgMotifVectorList.read_vectors
	cases = list()
	vectors = list()
	graph, header, casesD, vectorsD = read(default)
	undone, done = set(casesD), set()
	for filepath in ArgParser.get_filepaths(clusters):
		_, headerX, casesX, vectorsX = read(filepath)
		assert(header == headerX)
		assert(len(done & set(casesX)) == 0)
		done |= set(casesX)
		undone -= set(casesX)
		cases += casesX
		vectors += vectorsX
	print(len(undone), len(done))
	for case in sorted(undone):
		cases.append(case)
		vectors.append(vectorsD[casesD.index(case)])
	RNAMgMotifVectorList.write_vectors(outfile, graph, header, cases, vectors)



def plot(dividefile, vectors, n_components=4, save=None):
	divisions = MotifDivideFile.read_divisions(dividefile)
	cases = _cases_list(divisions.values())
	data = _obtain_data(cases, vectors)
	test = _cases_list(divisions[divno] for divno in divisions if divno != 0)
	train = _cases_list(divisions[divno] for divno in divisions if divno > 0)
	# extra = _cases_list(divisions[divno] for divno in divisions if divno == 0)
	_, traindata, trainclasses = _data_assigned(divisions, cases, data, train)
	_, testdata, testclasses = _data_assigned(divisions, cases, data, test)
	# 
	_plot_PCA(n_components, traindata, trainclasses, testdata, testclasses)
	# _plot_LDA(n_components, traindata, trainclasses, testdata, testclasses, test)

def _plot_PCA(n_components, traindata, trainclasses, testdata, testclasses):
	transformer = PCA(n_components=n_components)
	transformer.fit(traindata)
	n_components = transformer.explained_variance_.shape[0]
	print(sum(transformer.explained_variance_ratio_[:n_components]))
	testout = transformer.transform(testdata)
	RNAVectorAnalyser.plot_reduced([], testout, testclasses, n_components, None)

# def _plot_LDA(n_components, traindata, trainclasses, testdata, testclasses, testcases):
# 	transformer = LDA(n_components=n_components)
# 	transformer.fit(traindata, trainclasses)
# 	testout = transformer.transform(testdata)
# 	RNAVectorAnalyser.plot_reduced(testcases, testout, testclasses, n_components, None)

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
	# return sorted(set(sum(, [])))

def _obtain_data(cases, vectors):
	datalist = map(_data_dict, vectors)
	data = list()
	for case in cases:
		for datadict in datalist:
			if case in datadict:
				data.append(datadict[case])
				break
		else:
			print("Missing %s" % case, sys.stderr)
	return data

def _data_dict(datafile):
	data = RNAMgMotifVectorList.read_dataframe(datafile)
	cases, data = RNAVectorAnalyser.prepare_data(data, angle=False)
	data = dict(zip(cases, map(list, data.to_records(index=False))))
	return data


def parse_args():
	parser = argparse.ArgumentParser()
	commands = parser.add_subparsers(dest="command")
	# 
	contexts = commands.add_parser('contexts')
	contexts.add_argument("infolder", metavar="IN-FOLDER")
	contexts.add_argument("dividefile", metavar="DIVIDE-FILE")
	contexts.add_argument("division", type=int, metavar="DIVISION")
	contexts.add_argument("outfile", metavar="OUT-FILE")
	# 
	merge = commands.add_parser('merge')
	merge.add_argument("merged_file", metavar="MERGED-FILE")
	merge.add_argument("general_file", metavar="GENERAL-FILE")
	merge.add_argument("sources", nargs="+", metavar="IN-FILE")
	# 
	vectorise = commands.add_parser('vectorise')
	vectorise.add_argument("name", metavar="CONSENSUS-GRAPH")
	vectorise.add_argument("base", metavar="BASE-GRAPH")
	vectorise.add_argument("structuredir", metavar="STRUCTURE-DIR")
	vectorise.add_argument("assignsdir", metavar="ASSIGNS-FOLDER")
	vectorise.add_argument("outfile", metavar="OUT-FILE")
	vectorise.add_argument("motifs", nargs='+', metavar="MOTIFS")
	# 
	finalise = commands.add_parser('finalise')
	finalise.add_argument("outfile", metavar="FINAL-VECTORS")
	finalise.add_argument("default", metavar="DEFAULT-VECTORS")
	finalise.add_argument("clusters", nargs='+', metavar="CLUSTER-VECTORS")
	# 
	plot = commands.add_parser('plot')
	plot.add_argument("dividefile", metavar="DIVIDE-FILE")
	plot.add_argument("vectors", nargs="+", metavar="VECTORS")
	plot.add_argument("--n-components", "-n", metavar="N", type=float)
	# 
	return parser.parse_args()

COMMAND = {"contexts": contexts, "merge": merge, "vectorise": vectorise,
			"finalise": finalise, "plot": plot}

if __name__ == "__main__":
	args = vars(parse_args())
	command = args.pop("command")
	COMMAND[command](**args)
