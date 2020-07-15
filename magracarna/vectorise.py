
import numpy as np

from itertools import combinations, count, product
from collections import defaultdict
from decimal import Decimal
from pandas import DataFrame
from scipy import stats

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from sklearn.preprocessing import scale
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA


from .engrid import Residue, SecondaryStructureManager as SSManager
from .engraph import RNAGraph, RNALabels
from .motifs import MotifAssignments, CaseReader, MotifDivideFile

SSManager.WARN = False


class RNAProperties(object):

	class NotInChain(ValueError):
		"""
		"""

	class MissingAtom(ValueError):
		""" An atom is missing, thus property could not 
		    be measured.
		"""

	METHODS = dict()

	ANGLE = "A"
	TORSION = "T"
	DISTANCE = "D"
	PROPERTY_TYPE = dict()
	for __property in ("rho", "sigma", "cpc", "pcp"):
		PROPERTY_TYPE[__property] = ANGLE
	for __property in ("eta4", "theta4", "eta1", "theta1"):
		PROPERTY_TYPE[__property] = TORSION

	@staticmethod
	def residue(identifier, structure):
		residueid = identifier.split(".")[0]
		resprops = Residue.from_string(residueid, True)
		residue = structure.find(**resprops)
		if residue:
			return residue
		raise cls.NotInChain("%s not in nucleotide chains" % residueid)

	@classmethod
	def segment(cls, rindex, structure):
		for begex, endex, chain in structure.get_connected_segments():
			if begex <= rindex < endex:
				return (begex, endex, chain)
		residueid = str(structure.residues[rindex])
		raise cls.NotInChain("%s not in nucleotide chains" % residueid)

	@classmethod
	def previous_residue(cls, residue, structure):
		rindex = residue.index
		begex, endex, chain = cls.segment(rindex, structure)
		if begex != rindex:
			return structure.residues[rindex-1]
		residueid = str(structure.residues[rindex])
		raise cls.NotInChain("No residue before %s" % residueid)

	@classmethod
	def following_residue(cls, residue, structure):
		rindex = residue.index
		begex, endex, chain = cls.segment(rindex, structure)
		if rindex != endex-1:
			return structure.residues[rindex+1]
		residueid = str(structure.residues[rindex])
		raise cls.NotInChain("No residue after %s" % residueid)
		

	@classmethod
	def residue_arguments(cls, resid, structure, *args):
		this_residue = cls.residue(resid, structure)
		arguments = list()
		residue = this_residue
		for argument in args:
			if argument in ("this", "prev", "next"):
				if argument == "this":
					residue = this_residue
				elif argument == "prev":
					residue = cls.previous_residue(this_residue, structure)
				elif argument == "next":
					residue = cls.following_residue(this_residue, structure)
				arguments.append(residue)
			else:
				try:
					arguments.append(residue.atoms[argument])
				except KeyError:
					errargs = (str(residue), argument)
					raise cls.MissingAtom("%s has no %s atom " % errargs)
		return arguments

	@classmethod
	def measure_eta4(cls, identifiers, structure):
		args = ["prev", "C4'", "this", "P", "C4'", "next", "P"] 
		args = cls.residue_arguments(identifiers[0], structure, *args)
		prev_residue, prev_catom = args[:2]
		this_residue, this_patom, this_catom = args[2:5]
		next_residue, next_patom = args[5:]
		return prev_catom.dihedral_towards(this_patom, this_catom, next_patom)

	@classmethod
	def measure_theta4(cls, identifiers, structure):
		args = ["this", "P", "C4'", "next", "P", "C4'"] 
		args = cls.residue_arguments(identifiers[0], structure, *args)
		this_residue, this_patom, this_catom = args[:3]
		next_residue, next_patom, next_catom = args[3:]
		return this_patom.dihedral_towards(this_catom, next_patom, next_catom)

	@classmethod
	def measure_eta1(cls, identifiers, structure):
		args = ["prev", "C1'", "this", "P", "C1'", "next", "P"]
		args = cls.residue_arguments(identifiers[0], structure, *args)
		prev_residue, prev_catom = args[:2]
		this_residue, this_patom, this_catom = args[2:5]
		next_residue, next_patom = args[5:]
		return prev_catom.dihedral_towards(this_patom, this_catom, next_patom)

	@classmethod
	def measure_theta1(cls, identifiers, structure):
		args = ["this", "P", "C1'", "next", "P", "C1'"]
		args = cls.residue_arguments(identifiers[0], structure, *args)
		this_residue, this_patom, this_catom = args[:3]
		next_residue, next_patom, next_catom = args[3:]
		return this_patom.dihedral_towards(this_catom, next_patom, next_catom)

	@classmethod
	def measure_rho(cls, identifiers, structure):
		args = ["prev", "P", "this", "P", "next", "P"]
		args = cls.residue_arguments(identifiers[0], structure, *args)
		prev_residue, prev_patom = args[:2]
		this_residue, this_patom = args[2:4]
		next_residue, next_patom = args[4:]
		return this_patom.angle_between(next_patom, prev_patom)

	@classmethod
	def measure_sigma(cls, identifiers, structure):
		args = ["prev", "C1'", "this", "C1'", "next", "C1'"]
		args = cls.residue_arguments(identifiers[0], structure, *args)
		prev_residue, prev_catom = args[:2]
		this_residue, this_catom = args[2:4]
		next_residue, next_catom = args[4:]
		return this_catom.angle_between(next_catom, prev_catom)

	@classmethod
	def measure_cpc(cls, identifiers, structure):
		args = ["prev", "C1'", "this", "P", "C1'"]
		args = cls.residue_arguments(identifiers[0], structure, *args)
		prev_residue, prev_catom = args[:2]
		this_residue, this_patom, this_catom = args[2:5]
		return this_patom.angle_between(prev_catom, this_catom)

	@classmethod
	def measure_pcp(cls, identifiers, structure):
		args = ["this", "P", "C1'", "next", "P"]
		args = cls.residue_arguments(identifiers[0], structure, *args)
		this_residue, this_patom, this_catom = args[:3]
		next_residue, next_patom = args[3:]
		return this_catom.angle_between(this_patom, next_patom)

	@classmethod
	def measure(cls, field, identifiers, structure):
		if field not in cls.METHODS:
			method_name = "measure_%s" % field
			try:
				cls.METHODS[field] = getattr(RNAProperties, method_name)
			except AttributeError:
				raise ValueError("Unknown property %s" % field)
		try:
			return cls.METHODS[field](identifiers, structure)
		except (cls.NotInChain, cls.MissingAtom):
			return np.nan


class Vectoriser(object):

	NODE, EDGE, PHOSPHATE, NUCLEOSIDE, NUCLEOTIDE = range(5)
	EMPTY_TYPES = (NODE, EDGE)
	NODELABEL_TYPE = {RNALabels.NODELABELS[RNALabels.ANY_NUCLEOTIDE]: NUCLEOTIDE,
					  RNALabels.NODELABELS[RNALabels.PHOSPHATE]: PHOSPHATE}
	EDGELABEL_TYPE = {}

	HEADER = {PHOSPHATE: ("pcp", "rho", "eta4", "theta4", "eta1", "theta1"),
			  NUCLEOSIDE: ("cpc", "sigma")}
	HEADER[NUCLEOTIDE] = HEADER[NUCLEOSIDE] + HEADER[PHOSPHATE]

	PHLABEL = RNALabels.NODELABELS[RNALabels.PHOSPHATE]
	ISPARTOF = RNALabels.LABELS.edge[RNALabels.IS_PART_OF]

	def __init__(self, motifgraph):
		self.code = motifgraph
		self.tags = list()
		self.types = list()
		self.header = None
		self.indices = None
		self._setup()

	def measure(self, field, identifiers, structure):
		value = RNAProperties.measure(field, identifiers, structure)
		return Decimal("%.2f" % round(value, 2))

	def vectorise(self, embed, structure):
		vector = list()
		for index, part in enumerate(embed):
			tagtype = self.types[index]
			if tagtype not in self.HEADER:
				continue
			identifiers = [embed[index_] for index_ in self.indices[index]]
			for field in self.HEADER[tagtype]:
				vector.append(self.measure(field, identifiers, structure))
		return vector

	def _setup(self):
		self._set_tags()
		self._set_indices()
		self._set_nucleosides()
		self._set_header()

	def _set_tags(self):
		self.tags = list()
		for isnode, tagstr, tagtuple in self.code._iter_tagstr():
			self.tags.append(tagstr)
			label = tagtuple[1] if isnode else tagtuple[2]
			TOTYPE = self.NODELABEL_TYPE if isnode else self.EDGELABEL_TYPE
			for typelabel, tagtype in TOTYPE.items():
				if tagtype != self.NUCLEOSIDE and typelabel.covers(label):
					self.types.append(tagtype)
					break
			else:
				self.types.append(self.NODE if isnode else self.EDGE)

	def _set_header(self):
		self.header = list()
		for tagstr, tagtype in zip(self.tags, self.types):
			if tagtype in self.HEADER:
				for field in self.HEADER[tagtype]:
					self.header.append("%s %s" % (tagstr, field))

	def _set_indices(self):
		self.indices = list()
		for index, tagstr in enumerate(self.tags):
			if "->" not in tagstr:
				self.indices.append((index,))
				continue
			nodetags = map(int, tagstr.split(" ")[1].strip("()").split("->"))
			nodestr = "%d %s" % (nodetags[0], self.code._tags[nodetags[0]])
			edonstr = "%d %s" % (nodetags[1], self.code._tags[nodetags[1]])
			nodex, edonex = map(self.tags.index, (nodestr, edonstr))
			self.indices.append((nodex, edonex, index))

	def _set_nucleosides(self):
		codegf = RNAGraph.from_code(self.code)
		for index, tagstr, tagtype in zip(count(), self.tags, self.types[:]):
			if tagtype == self.NUCLEOTIDE:
				nucid = int(tagstr.split(" ")[0])
				if self._has_phosphate_node(nucid, codegf):
					self.types[index] = self.NUCLEOSIDE

	def _has_phosphate_node(self, nucid, codegf):
		edges = codegf.incident_edges(nucid)
		edges = codegf.edges_with_label_in([self.ISPARTOF], edges)
		for nodeid, edonid, _ in edges:
			partid = edonid if nucid == nodeid else nodeid
			if self.PHLABEL.covers(codegf.nodes[partid]):
				return True
		return False


class MotifVectorList(object):

	def __init__(self, motifgraph, assigner,
		         restrict=None, division=None, limited=None, mode="all"):
		self.graph = motifgraph
		self.cases = list()
		self._mode = mode
		self.vectors = list()
		self.assigner = assigner
		self._limited = None
		self._restrict = None
		self._set_vectoriser()
		if limited:
			limited = CaseReader.read_limited(limited, restrict, division)
			self._limited = limited
		elif restrict in assigner.definedin:
			self._restrict = restrict

	def _set_vectoriser(self):
		self._vectoriser = Vectoriser(self.assigner.graphs[self.graph])

	def _check_case(self, case, assigns):
		if self._restrict and self._restrict not in assigns.by_case[case]:
			return False
		if self._limited and case not in self._limited:
			return False
		return True

	def _load_assigns(self, assignsfile):
		limited_motifs = {self.graph, self._restrict} - {None}
		arguments = (self.assigner, assignsfile, True, limited_motifs)
		assigns = MotifAssignments.load(*arguments)
		if self.graph in assigns.by_graph:
			return assigns
		return None

	def _get_assigns_order(self, assigns):
		motifgraph = self._vectoriser.code
		header = assigns._graph_header_details(self.graph)
		tagstrs = [tagstr for _, tagstr, _ in motifgraph._iter_tagstr()]
		indices = [header.index(tagstr) for tagstr in tagstrs]
		return indices

	def _load_embeds(self, assignsfile):
		embeds = list()
		assigns = self._load_assigns(assignsfile)
		if assigns:
			order = self._get_assigns_order(assigns)
			casewise = assigns.by_graph[self.graph]
			for case, caseembeds in sorted(casewise.items()):
				if self._check_case(case, assigns):
					for bedex, embed in enumerate(caseembeds):
						caseembeds[bedex] = [embed[index] for index in order]
					embeds += self._case_embeds(case, caseembeds)
		return embeds

	def _case_embeds(self, case, caseembeds):
		if self._mode == "min":
			caseembeds = [sorted(caseembeds, key=self._embed_sorter)[0]]
		return [(case, embed) for embed in caseembeds]

	@classmethod
	def _embed_sorter(cls, embed):
		embed_key = [part for part in embed if part.startswith("[")]
		embed_key = map(Residue.from_string, embed_key)
		return [(part.resno, part.insco, part.chain) for part in embed_key]

	def add_vectors_from(self, assignsfile, structurepath, *runargs):
		embeds = self._load_embeds(assignsfile)
		if embeds:
			with SSManager(structurepath, *runargs) as manager:
				structure = manager.extract_paired_structure()
			for case, embed in embeds:
				vector = self._vectoriser.vectorise(embed, structure)
				self.cases.append(case)
				self.vectors.append(vector)

	def dump_vectors(self, outfilepath):
		header = self._vectoriser.header
		args = (outfilepath, self.graph, header, self.cases, self.vectors)
		self.write_vectors(*args)

	@classmethod
	def write_vectors(cls, outfilepath, graph, header, cases, vectors):
		out_string = ["Graph\t%s\t%s" % (graph, "\t".join(header))]
		for case, vector in zip(cases, vectors):
			entry = "\t".join(map(str, vector))
			entry = "%s\t%s" % (case.replace(" ", "\t"), entry)
			out_string.append(entry)
		with open(outfilepath, 'w') as outfile:
			outfile.write("\n".join(out_string)+"\n")

	@classmethod
	def read_vectors(cls, infilepath):
		vectors = list()
		header = None
		cases = list()
		graph = None
		with open(infilepath, 'r') as infile:
			for line in infile:
				line = line.strip("\n")
				if line.startswith("Graph"):
					header = [field.strip() for field in line.split("\t")]
					graph, header = header[1], header[2:]
				elif line:
					entry = [field.strip() for field in line.split("\t")]
					cases.append(" ".join(entry[:2]))
					vectors.append(list(map(Decimal, entry[2:])))
				else:
					break
		return graph, header, cases, vectors

	@classmethod
	def read_dataframe(cls, infilepath):
		graph, header, cases, vectors = cls.read_vectors(infilepath)
		header.insert(0, "SITE ID")
		header.insert(0, "PDB ID")
		data = list()
		for case, vector in zip(cases, vectors):
			data.append(case.split(" ")+vector)
		data = dict(zip(header, list(zip(*data))))
		return DataFrame(data)

	@classmethod
	def rewrite_limited(cls, infilepath, outfilepath,
						limited, restrict=None, division=None):
		limited_cases = CaseReader.read_limited(limited, restrict, division)
		graph, header, cases, vectors = cls.read_vectors(infilepath)
		filtered = list()
		for case, vector in zip(cases, vectors):
			if case in limited_cases:
				filtered.append((case, vector))
		cases, vectors = map(list, zip(*filtered)) if filtered else ([], [])
		cls.write_vectors(outfilepath, graph, header, cases, vectors)


class RNAVectorAnalyser(object):

	COLORS = ['silver', 'red', 'yellow', 'blue', 'orange', 'green', 'purple',
			  'maroon', 'gold', 'navy', 'brown', 'olive', 'indigo', 'fuchsia',
			  'lime', 'aqua', 'teal', 'gray', 'black']
	
	@classmethod
	def prepare_data(cls, data, distance=True, angle=True, torsion=True):
		PROPERTY_TYPE = RNAProperties.PROPERTY_TYPE
		header = ["PDB ID", "SITE ID"]
		cases = zip(data["PDB ID"], data["SITE ID"])
		cases = [" ".join(case) for case in cases]
		transformed = dict()
		indices = set()
		scale = False
		for field in data.drop(["PDB ID", "SITE ID"], axis=1):
			fieldtype = PROPERTY_TYPE[field.split(" ")[-1]]
			column = np.array(data[field], dtype=float)
			indices.update(np.nonzero(np.isnan(column))[0])
			if fieldtype == RNAProperties.DISTANCE and distance:
				scale = True
				transformed[field] = column
			elif fieldtype == RNAProperties.ANGLE and angle:
				continue
				transformed["%s_cos" % field] = np.cos(np.radians(column))
			elif fieldtype == RNAProperties.TORSION and torsion:
				transformed["%s_sin" % field] = np.sin(np.radians(column))
				transformed["%s_cos" % field] = np.cos(np.radians(column))
		transformed = DataFrame(transformed)
		indices = sorted(indices, reverse=True)
		transformed = transformed.drop(indices, axis=0)
		for index in indices:
			cases.pop(index)
		if scale:
			transformed = scale(transformed)
		return cases, transformed

	@classmethod
	def calculate_mean(cls, field, column):
		field = field.split(" ")[-1]
		fieldtype = RNAProperties.PROPERTY_TYPE[field]
		if fieldtype == RNAProperties.DISTANCE:
			return np.mean(column)
		elif fieldtype == RNAProperties.ANGLE:
			return stats.circmean(column, low=0., high=180.)
		elif fieldtype == RNAProperties.TORSION:
			return stats.circmean(column, low=-180., high=180.)
		raise ValueError("Bad field %s" % field)

	@classmethod
	def calculate_std(cls, field, column):
		field = field.split(" ")[-1]
		fieldtype = RNAProperties.PROPERTY_TYPE[field]
		if fieldtype == RNAProperties.DISTANCE:
			return np.std(column)
		elif fieldtype == RNAProperties.ANGLE:
			return stats.circstd(column, low=0., high=180.)
		elif fieldtype == RNAProperties.TORSION:
			return stats.circstd(column, low=-180., high=180.)
		raise ValueError("Bad field %s" % field)

	@classmethod
	def calculate(cls, vectorsfile, fields, measures):
		methods = list()
		for value in measures:
			method_name = "calculate_%s" % value
			methods.append(getattr(cls, method_name))
		data = MotifVectorList.read_dataframe(vectorsfile)
		output = defaultdict(dict)
		for field in fields:
			for value, method in zip(measures, methods):
				column = np.array(data[field], dtype=float)
				output[value][field] = method(field, column)
		return DataFrame(dict(output))

	@classmethod
	def assign_divisions(cls, cases, dividefile):
		if not dividefile:
			return np.zeros(len(cases), dtype=int)
		caseset = set(cases)
		assigned = dict()
		divisions = MotifDivideFile.read_divisions(dividefile)
		max_div = max(divisions)
		max_len = len(cls.COLORS)
		for division, (_, divcases) in divisions.items():
			divcases = list(zip(*divcases))[0]
			for case in set(divcases) & caseset:
				assigned[case] = division if division < max_len-2 else\
								 -2 if division != max_div else -1
		aslist = [assigned[case] if case in assigned else 0 for case in cases]
		return np.array(aslist)

	@classmethod
	def perform_pca(cls, vectorsfile, n_components, 
		            dividefile=None, imagefile=None):
		data = MotifVectorList.read_dataframe(vectorsfile)
		cases, data = cls.prepare_data(data)
		pca = PCA(n_components=n_components)
		reduced = pca.fit(data).transform(data)
		print("Explained: %s" % pca.explained_variance_ratio_)
		assigned = cls.assign_divisions(cases, dividefile)
		cls.plot_reduced(cases, reduced, assigned, n_components, imagefile, mode)
		return pca, data, assigned

	@classmethod
	def perform_lda(cls, vectorsfile, n_components, 
		            dividefile=None, imagefile=None):
		data = MotifVectorList.read_dataframe(vectorsfile)
		cases, data = cls.prepare_data(data)
		lda = LDA(n_components=n_components)
		assigned = cls.assign_divisions(cases, dividefile)
		reduced = lda.fit(data, assigned).transform(data)
		cls.plot_reduced(cases, reduced, assigned, n_components, imagefile, mode)
		return lda, data, assigned

	@classmethod
	def plot_reduced(cls, cases, reduced, assigned, n_components, imagefile, label=None):
		cases = np.array(cases) if cases else None
		combos = Plotter.generate_combos(n_components, mode='min')
		divisions = list(set(assigned))
		fig, axarray, grid = Plotter.generate_subplots(len(combos))
		colours = Plotter.assign_colours(divisions)
		handles = dict()
		axlabel = "Component %d"
		annotated = 0
		divisions.sort(key=lambda divno: -sum(map(int, assigned==divno)))
		for divno in sorted(divisions, ):
			for (i, j), ax in zip(combos, axarray):
				mask = assigned == divno
				x_coos = reduced[mask, i]
				y_coos = reduced[mask, j]
				marker = "." if divno == 0 else 'x' if divno > 0 else "+"
				style = {"color": colours[divno], "marker": marker, "s": 30}
				handle = ax.scatter(x_coos, y_coos, **style)
				if len(handles) != len(divisions):
					handles[divno] = handle
				ax.set(xlabel= axlabel % (i+1), ylabel= axlabel % (j+1))
				ax.set_aspect(np.diff(ax.get_xlim())/np.diff(ax.get_ylim()))
				if label == "case" and cases is not None:
					for case, x_coo, y_coo in zip(cases[mask], x_coos, y_coos):
						ax.annotate(case, (x_coo, y_coo), fontsize=8)
				elif label == "division":
					for x_coo, y_coo in zip(x_coos, y_coos):
						ax.annotate(str(divno), (x_coo, y_coo), fontsize=8)
				elif label == "centre":
					x_coo, y_coo = np.mean(x_coos), np.mean(y_coos)
					ax.annotate(str(divno), (x_coo, y_coo), fontsize=5)
				elif label == "outlier":
					allmask = np.ones(len(reduced), dtype=bool)
					for case, x_coo, y_coo in zip(cases[mask], x_coos, y_coos):
						l1dist = np.square(reduced[allmask, i] - x_coo)
						l1dist += np.square(reduced[allmask, j] - y_coo)
						if assigned[np.argsort(l1dist)[1]] != divno:
							if annotated == 0:
								print(case, cases[np.argsort(l1dist)[1]])
								print(divno, assigned[np.argsort(l1dist)[1]])
							ax.annotate(case, (x_coo, y_coo), fontsize=8)
							annotated += 1
		print(annotated, len(assigned)*len(combos))
		divisions.sort(key=lambda divno: (divno < 0, abs(divno)))
		handles = [handles[divno] for divno in divisions]
		fig.legend(handles, divisions, bbox_to_anchor=(0.98, 0.9), loc="upper right")
		grid.tight_layout(fig, rect=[0, 0, 0.85, 1])
		if imagefile:
			plt.savefig(imagefile, dpi=800)
			plt.close()
		else:
			plt.show()


class Plotter(object):

	ARC = 255
	WHEEL = ARC*3

	@staticmethod
	def best_grid(number):
		assert(number > 0)
		for columns in count(1):
			for rows in range(1, columns+1):
				size = columns * rows
				if number <= size:
					return (rows, columns)
			assert(number > columns)

	@classmethod
	def generate_combos(cls, n_components, mode="max"):
		if mode == "max":
			return list(combinations(range(n_components), 2))
		if mode == "min":
			combos = list(range(n_components))
			if n_components % 2:
				combos.append(combos[0])
			return list(zip(combos[::2], combos[1::2]))
		raise ValueError("Could not recognise mode : %s" % str(mode))

	@classmethod
	def generate_subplots(cls, number):
		rows, columns = cls.best_grid(number)
		fig = plt.figure()
		grid = gridspec.GridSpec(rows, columns)
		grid.update(wspace=0.4, hspace=0.4)
		axarray = list()
		for figno, cell in enumerate(product(range(rows), range(columns))):
			if figno == number:
				break
			axarray.append(plt.subplot(grid[cell[0], cell[1]]))
		return fig, axarray, grid

	@classmethod
	def generate_colours(cls, number):
		step = cls.WHEEL / number
		code = 0
		colours = list()
		while len(colours) != number:
			arcvalue = cls.ARC - (code % cls.ARC)
			remainder = (cls.ARC - arcvalue)
			inarc = code / cls.ARC
			colour = [0, 0, 0]
			colour[inarc] = arcvalue
			colour[(inarc+1)%3] = remainder
			colours.append(cls.hexcode(colour))
			code += step
		return colours

	@classmethod
	def assign_colours(cls, divisions):
		colour_wheel = cls.generate_colours(max(divisions))
		colour_dict = dict()
		for divno in divisions:
			if divno == 0:
				colour_dict[0] = "#c0c0c0"
			elif divno < 0:
				colour_dict[divno] = cls.opposite_colour(colour_wheel[abs(divno)-1])
			else:
				colour_dict[divno] = colour_wheel[divno-1]
		return colour_dict

	@classmethod
	def darken_colour(cls, hexcode):
		return cls.hexcode([value / 2 for value in cls.rgbtuple(hexcode)])

	@classmethod
	def opposite_colour(cls, hexcode):
		rvalue, gvalue, bvalue = cls.rgbtuple(hexcode)
		rvalue_ = (bvalue + gvalue) / 2
		gvalue_ = (rvalue + bvalue) / 2
		bvalue_ = (rvalue + gvalue) / 2
		return cls.hexcode([rvalue_, gvalue_, bvalue_])


	@staticmethod
	def hexcode(rgbtuple):
		code = [("%x" % value).zfill(2) for value in rgbtuple]
		return "#%s%s%s" % tuple(code)

	@staticmethod
	def rgbtuple(hexcode):
		rgbtuple = list()
		hexcode = hexcode.strip("#").lower()
		for start in range(0, 6, 2):
			rgbtuple.append(int(hexcode[start:(start+2)], 16))
		return tuple(rgbtuple)
