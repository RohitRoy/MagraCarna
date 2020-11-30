from __future__ import print_function

import re
import os
import sys
import time
import argparse
import numpy as np
from numpy.linalg import norm
from subprocess import call, Popen, PIPE

from scipy.cluster.hierarchy import linkage
# from sklearn.cluster import AgglomerativeClustering as AC
from scipy.spatial.distance import pdist, squareform
from scipy.sparse.csgraph import connected_components
from scipy.sparse import csr_matrix
from collections import defaultdict
from itertools import count, combinations, product

from ..utilities.args import ArgParser
from ..magracarna.motifs import MotifRedivider
from ..magracarna.engrid import SecondaryStructureManager as SSManager,\
								HydrogenBondFinder as HBFinder, Structure,\
								StructureFile, Residue, Atom
from ..magracarna.metalsite import RNAMgSiteEngrapher
from ..magracarna.vectorise import Plotter
from ..magracarna.nucleotides import variants, Borgs, Bopps, Bfeet, backbone, phosphate
from ..magracarna.heteroatoms import cations as CATIONS

MSA_COMMAND = ["clustalo", "--seqtype=RNA", "--infmt=fa", "--force"]
BACKBONE = [atom for atom in set(phosphate) | set(backbone) if "O" in atom]
WARN = True

class MoleculeDataset(object):

	def __init__(self, molecule, organism, ids=None, chainsdir=None, msafile=None):
		self.molecule = molecule
		self.organism = organism
		self.chainsdir = chainsdir
		self.msafile = msafile
		#
		if self.chainsdir:
			self.chains = self.chains_from_folder(ids)
		elif self.msafile:
			self.chains = self.chains_from_msa(ids)
		else:
			raise ValueError("Either 'chainsdir' or 'msafile' must be input")

	def iterate(self, count=False, structuredir=None, sitesdir=None,
				fastadir=None, align=False, motif=None, start=None):
		structure = Structure(None, [], None)
		started = True if not start else False
		for index, (pdbid, chainid) in enumerate(self.chains, 1):
			if not started:
				if pdbid == start:
					started = True
				else:
					continue
			toyield = [pdbid, chainid]
			if count:
				toyield.append(index)
			if structuredir:
				if structure.structid != pdbid:
					structure, _ = self.get_structure(structuredir, pdbid)
				toyield.append(structure)
			if sitesdir:
				path = self.get_file(sitesdir, pdbid, "_sites.tsv")
				toyield.append(path)
			if fastadir:
				fastafile = self.get_file(fastadir, pdbid, ".fasta")
				fastachain = self.read_fasta(fastafile, pdbid, chainid)
				toyield.append(fastachain)
			if structuredir and align:
				toyield.append(self.align(pdbid, chainid, structure))
			if structuredir and motif:
				toyield.append(self.embed(motif, structure, chainid))
			yield tuple(toyield)

	@classmethod
	def get_file(cls, folder, pdbid, *extensions):
		folderlist = os.listdir(folder)
		for extension in extensions:
			filename = "%s%s" % (pdbid, extension)
			if filename in folderlist:
				return os.path.join(folder, filename)
		raise ValueError("No File Found")

	@classmethod
	def get_structure(cls, structuredir, pdbid):
		path = cls.get_file(structuredir, pdbid, ".pdb", ".cif")
		with StructureFile(path, None) as manager:
			structure = manager.extract_structure(breaks=True)
		return structure, path

	def chains_from_folder(self, ids):
		chains = list()
		pdbids = ArgParser.read_idslist(ids)
		chainsfiles = ArgParser.get_sources(self.chainsdir, pdbids, "_entities.txt")
		for filepath in chainsfiles:
			pdbid = ArgParser.basename(filepath, "_.")
			entities = MotifRedivider.read_chains_file(filepath)
			for chainid, (molecule, organism) in entities.items():
				if molecule == self.molecule and self.organism in (organism, "*"):
					chains.append((pdbid, chainid))
		return chains

	def chains_from_msa(self, ids):
		chains = list()
		with open(self.msafile, 'r') as infile:
			for line in infile:
				if line.startswith(">"):
					pdbid, chainid = line.strip().strip(">").split(":")
					if not (ids and pdbid not in ids):
						chains.append((pdbid, chainid))
		return chains

	@staticmethod
	def generalise(motif):
		motif = motif.replace("N", ".")
		motif = motif.replace("Y", "[CU]").replace("R", "[AG]")
		return motif

	@staticmethod
	def simplify(base):
		for rnabase in ("A", "C", "G", "U"):
			if base in variants[rnabase]:
				return rnabase
		if base == "PSU":
			return "U"
		return "N"

	def representatives(self, structuredir, resolutionsfile):
		pdbids = tuple(set(list(zip(*self.chains))[0]))
		resolutions = dict()
		with open(resolutionsfile, 'r') as infile:
			for line in infile:
				if line.startswith(pdbids):
					pdbid, resolution = line.split()
					resolutions[pdbid] = float(resolution)
		chains = defaultdict(list)
		for pdbid, chain in self.chains:
			chains[pdbid].append(chain)
		order = sorted(resolutions.items(), key=lambda item:(item[1],item[0]))
		inorder = list()
		for pdbid, resolution in order:
			inorder.append([resolution, pdbid] + chains[pdbid])
		return inorder

	def align(self, pdbid, chainid, structure):
		alignment = list()
		begex, endex, _ = structure.chains[chainid]
		iterator = iter(structure.residues[begex:endex])
		chainmsa = "".join(self.read_fasta(self.msafile, pdbid, chainid))
		for locus, basemsa in enumerate(chainmsa, 1):
			if basemsa in ("-", "P"):
				alignment.append(None)
			else:
				try:
					residue = next(iterator)
					assert(self.simplify(residue.name) == basemsa)
					alignment.append(residue)
				except StopIteration:
					if WARN:
						message = "Locus %d: MSA base not in structure"
						print(message % locus, file=sys.stderr)
					alignment.append(None)
				except AssertionError:
					if WARN:
						message = "Locus %d: Matching MSA %s with %s"
						message = message % (locus, basemsa, str(residue))
						print(message, file=sys.stderr)
					alignment.append(residue)
		return alignment

	@staticmethod
	def read_fasta(fastafile, pdbid, chainid):
		chain = list()
		with open(fastafile, 'r') as inputfile:
			for line in inputfile:
				if line.startswith(">%s" % pdbid)\
				  and line.split(":")[1].strip() == chainid:
					for line in inputfile:
						if line.startswith(">"):
							break
						chain.append(line.strip())
					break
		return chain

	@classmethod
	def embed(cls, motif, structure, chainid):
		length = len(motif)
		begex, endex, _ = structure.chains[chainid]
		chainsequence = map(cls.simplify, structure.residues[begex:endex])
		for match in re.finditer(cls.generalise(motif), chainsequence):
			embedding = list()
			start = match.start() + begex
			for rindex in range(start, start+length):
				embedding.append(structure.residues[rindex])
			embeddings.append(embedding)
		return embeddings


def sequences(structuredir, chainsdir, molecule, organism, output, ids=None, msa=None):
	fasta = list()
	dataset = MoleculeDataset(molecule, organism, ids, chainsdir=chainsdir)
	for pdbid, chainid, fastachain in dataset.iterate(fastadir=structuredir):
		fastachain.insert(0, "\n>%s:%s" % (pdbid, chainid))
		fasta.append(fastachain)
	with open(output, 'w') as outputfile:
		outputfile.write("\n".join(sum(fasta, [])))
	if msa:
		command = MSA_COMMAND[:]
		command += ["--infile=%s" % output, "--outfile=%s" % msa]
		with open(os.devnull, 'w') as NULL:
			exit_code = call(command, stdout=NULL, stderr=sys.stderr)


def check(msafile, structuredir):
	dataset = MoleculeDataset("*", "*", None, msafile=msafile)
	allchains = len(dataset.chains)
	iterator = dataset.iterate(count=True, structuredir=structuredir, align=True)
	for pdbid, chainid, count, _, _ in iterator:
		print("% 3d / % 3d\t%s\t%s" % (count, allchains, pdbid, chainid))


class EngraphedMetal(object):

	@classmethod
	def get_metal_id(cls, sitesfile, pdbid, metalresidue):
		metalid = str(metalresidue)
		if metalresidue.name != "MG":
			return "%s " % metalid
		with open(sitesfile, 'r') as infile:
			for line in infile:
				if line.startswith("ContextGraph"):
					graphid = line.split()[2]
					if graphid.rstrip("!") != metalid:
						continue
					return graphid if graphid.endswith("!") else graphid+" "
		return "%s*" % metalid

	@classmethod
	def get_ligands(cls, sitesfile, pdbid, metalid):
		ligands = list()
		if metalid.endswith("*") or not metalid.startswith("[MG]"):
			return ligands
		siteid = "%s %s" % (pdbid, metalid.strip())
		rnagf = RNAMgSiteList.load(sitesfile, [siteid])[0]
		metalid = metalid.strip("! ")
		ligands = list()
		for nodeid, edonid, edgelabel in rnagf.incident_edges(metalid):
			if edgelabel.name in ("L", "LW1"):
				edonid = edonid if edonid != metalid else nodeid
				ligands.append((edgelabel.name, edonid))
		ligands = filter(lambda item: not item[1].startswith("[HOH]"), ligands)
		return sorted(ligands)

	@classmethod
	def print_ligands(cls, sitesfile, pdbid, structure, residue, atom=None):
		"""
		For a residue, locates nearby metal ions and prints ligands if Mg(2+).
		"""
		metal = MetalLocator.nearest_mg(structure, residue, atom)
		if not metal:
			metal = MetalLocator.nearest_metal(structure, residue, atom)
		siteid = EngraphedMetal.get_metal_id(sitesfile, pdbid, metal)
		sitestring = "%s % 15s % 15s" % (pdbid, siteid, str(residue))
		#
		ligands = list()
		for label, atom in cls.get_ligands(sitesfile, pdbid, siteid):
			ligands.append("%s-%s" % (label, atom))
		print("%s\t%s" % (sitestring, "\t".join(ligands)))

	@classmethod
	def ligands(cls, structuredir, sitesdir, chainsdir, sequence,
				molecule, organism, ids=None, atom=None, position=None):
		"""
		"""
		position = len(sequence) / 2 if not position else position-1
		dataset = MoleculeDataset(molecule, organism, ids, chainsdir=chainsdir)
		iterator = dataset.iterate(count=True, structuredir=structuredir, sitesdir=sitesdir, motif=sequence)
		#
		logged = time.time()
		print("Logging", file=sys.stderr)
		for pdbid, chainid, count, structure, sitesfile, embeddings in iterator:
			if time.time() - logged > 120:
				print("\t%d / %d" % (count, len(molecule_chains)), file=sys.stderr)
				logged = time.time()
			try:
				assert(len(embeddings) == 1)
				residue = embeddings[0][position]
				cls.print_ligands(sitesfile, pdbid, structure, residue, atom)
			except BaseException as e:
				print("%s %s\tMissed" % (pdbid, chainid), file=sys.stderr)


class NucleotideFrame(object):

	def __init__(self, residue):
		self.residue = residue
		self.p_origin, self.p_frame = self.create_p_frame()
		self.b_origin, self.b_frame = self.create_b_frame()

	def create_p_frame(self):
		try:
			p_coos = self.residue.atoms["P"].coos
			o3r_coos = self.residue.atoms["O3'"].coos
			c1r_coos = self.residue.atoms["C1'"].coos
			borg_coos = self.residue.atoms[Borgs[self.residue.name]].coos
			#
			rB_bond = (borg_coos - c1r_coos) # glycosidic bond
			#
			p3_axis = (p_coos - o3r_coos)
			p3_axis /= norm(p3_axis)
			rB_axis = rB_bond - (np.dot(rB_bond, p3_axis) * p3_axis)
			rB_axis /= norm(rB_axis)
			pX_axis = np.cross(p3_axis, rB_axis)
			p_frame = [p3_axis, rB_axis, pX_axis]
			return p_coos, p_frame
		except:
			return np.array([np.nan]*3), np.array([[np.nan]*3]*3)

	def create_b_frame(self):
		try:
			borg_coos = self.residue.atoms[Borgs[self.residue.name]].coos
			bopp_coos = self.residue.atoms[Bopps[self.residue.name]].coos
			foot_coos = self.residue.atoms[Bfeet[self.residue.name]].coos
			#
			b_coos = (borg_coos + bopp_coos) / 2
			bW_edge = foot_coos - bopp_coos
			#
			HS_axis = bopp_coos - borg_coos
			HS_axis /= norm(HS_axis)
			WR_axis = bW_edge - (np.dot(bW_edge, HS_axis) * HS_axis)
			WR_axis /= norm(WR_axis)
			bX_axis = np.cross(HS_axis, WR_axis)
			b_frame = [HS_axis, WR_axis, bX_axis]
			return b_coos, b_frame
		except:
			return np.array([np.nan]*3), np.array([[np.nan]*3]*3)

	def transform(self, absolutecoos, mode):
		if mode not in ("P", "B"):
			raise ValueError("Mode must be 'P' or 'B'")
		if "P" == mode:
			return np.dot(self.p_frame, absolutecoos - self.p_origin)
		if "B" == mode:
			return np.dot(self.b_frame, absolutecoos - self.b_origin)

	def reverse(self, relativecoos, mode):
		if mode not in ("P", "B"):
			raise ValueError("Mode must be 'P' or 'B'")
		if "P" == mode:
			return np.dot(relativecoos, self.p_frame) + self.p_origin
		if "B" == mode:
			return np.dot(relativecoos, self.b_frame) + self.b_origin

	def transform_self(self, mode):
		atoms = [atom for atom in self.residue.atoms.values()]
		print_list = list()
		for atom in sorted(atoms, key=lambda atom: atom.atmno):
			coos = self.transform(atom.coos, mode)
			print_list.append(["% 8.3f" % coo_ for coo_ in coos]+[atom.name])
		for line in print_list:
			print("\t\t" + "\t".join(line))


class MetalLocator(object):

	def __init__(self, reference, ligandsfile, max_da):
		self.reference = reference
		self.ligands = LigandsFile(ligandsfile)
		self.ligands.new()
		self.metals = None
		hbfinder = HBFinder(max_da=max_da)
		self.mgion = RNAMgSiteEngrapher(hbfinder, invalids=True, relaxed=True)
		self.coos = None
		#
		self._refframes = dict()

	def refframes(self, refresidue):
		if refresidue not in self._refframes:
			self._refframes[refresidue] = NucleotideFrame(refresidue)
		return self._refframes[refresidue]

	def locate(self, structuredir, msafile):
		starttime = time.time()
		chainwiseiterator = self.iterate_chains(structuredir, msafile)
		for pdbid, chainid, structure, alignment in chainwiseiterator:
			self.mgion.set_structure(structure)
			locusof = dict()
			for locus, residue in enumerate(alignment, 1):
				if residue:
					locusof[str(residue)] = locus
			locuswise = enumerate(zip(alignment, self.reference.alignment), 1)
			for locus, (residue, refresidue) in locuswise:
				if not residue or not refresidue:
					continue
				metals = self.interacting_metals(residue, structure, locusof)
				metals = self.locus_can_frame(metals, locus)
				if metals:
					frame = NucleotideFrame(residue)
					refframe = self.refframes(refresidue)
				for metal in metals:
					metalid = self.metals[str(metal)][1]
					for _, metalatom in metal.atoms.items():
						if metalatom.atype not in ("N", "O"):
							break
					for mode in ("B", "P"):
						coos = frame.transform(metalatom.coos, mode)
						profilecoos = refframe.reverse(coos, mode)
						if not np.isnan(profilecoos).any():
							self.coos[metalid].append(profilecoos)
			#
			order = sorted(self.metals, key=self.metal_resno)
			self.ligands.update(self.metals, order)
			self.reference.update(self.coos, order)
		self.reference.finish()
		print(time.time() - starttime)

	def locus_can_frame(self, metals, locus):
		canframe = list()
		for metal in metals:
			metal_residueid = str(metal)
			if metal_residueid in self.metals\
			  and len(self.metals[metal_residueid][2:])\
			  and locus not in dict(self.metals[metal_residueid][2:]):
				continue
			canframe.append(metal)
		return canframe

	@staticmethod
	def metal_resno(metalid):
		return (Residue.from_string(metalid).resno, metalid)

	def iterate_chains(self, structuredir, msafile, start=None):
		dataset = MoleculeDataset("*", "*", None, msafile=msafile)
		allchains = len(dataset.chains)
		chainiterator = dataset.iterate(start=start, count=True,
										structuredir=structuredir, align=True)
		pdbid_old = None
		for pdbid, chainid, count, structure, alignment in chainiterator:
			if pdbid_old != pdbid:
				self.metals = dict()
				self.coos = defaultdict(list)
			print("% 3d / % 3d\t%s\t%s" % (count, allchains, pdbid, chainid))
			yield (pdbid, chainid, structure, alignment)

	def interacting_metals(self, residue, structure, locusof):
		metals = self.nearby_metals(structure, residue)
		if not metals:
			return list()
		metals = list(zip(*metals)[0])
		self.detect_ligands(metals, structure, locusof)
		return metals

	def detect_ligands(self, metals, structure, locusof):
		for metal in metals:
			if self.metals.has_key(str(metal)):
				continue
			ligands = [structure.structid]
			if metal.name == "MG":
				try:
					self.mgion.detect_mgsite(metal, structure)
				except self.mgion.NoRNAInteraction:
					ligands.append("%s*" % metal)
				else:
					metalid = self.mgion.graph.name.split()[1]
					metalid = metalid + ("" if metalid.endswith("!") else " ")
					ligands.append(metalid)
					shells = self.mgion.inner_shell | self.mgion.outer_shell
					for entry in shells:
						residue, atom = entry[0], entry[1].name
						try:
							ligands.append((locusof[str(residue)], atom))
						except KeyError:
							continue
			else:
				ligands.append("%s " % str(metal))
			self.metals[str(metal)] = ligands

	@classmethod
	def nearby_metals(cls, structure, residue, atom=None):
		"""
		Returns list of metal residues at increasing distances.
		"""
		nearby = list()
		entity = residue.atoms[atom] if atom else residue
		nearbyatoms = structure.grid.get_nearby_atoms(7.0, entity)
		for cation, cationatom, distance in nearbyatoms:
			if cation.name in CATIONS and cationatom.atype not in ("O", "N"):
				nearby.append((cation, distance))
		metals = dict()
		for metal, distance in nearby:
			if metal not in metals or metals[metal] > distance:
				metals[metal] = distance
		return sorted(metals.items(), key=lambda item: item[1])

	@classmethod
	def nearest_metal(cls, structure, residue, atom):
		metals = cls.nearby_metals(structure, residue, atom)
		return None if not metals else metals[0][0]

	@classmethod
	def nearest_mg(cls, structure, residue, atom=None):
		metals = cls.nearby_metals(structure, residue, atom)
		mgresidues = filter(lambda item: item[0].name == "MG", metals)
		return None if not mgresidues else mgresidues[0][0]


class LigandsFile(object):

	def __init__(self, filepath):
		self.path = filepath
		self.ligands = dict()

	def new(self):
		with open(self.path, 'w') as outfile:
			outfile.write("")
		self.ligands = dict()

	def update(self, metals, order):
		lines = list()
		for metalid in order:
			ligands = [metals[metalid][0], "% 15s" % metals[metalid][1]]
			ligands += self.sorted_ligands(metals[metalid][2:])
			lines.append("%s\n" % "\t".join(ligands))
		with open(self.path, 'a') as outfile:
			outfile.write("".join(lines))

	def sorted_ligands(self, ligands):
		sorted_ligands = list()
		for locus, atom in sorted(ligands):
			sorted_ligands.append("% 4s.% -3s" % (locus, atom))
		return sorted_ligands

	@classmethod
	def read(self, filepath):
		with open(filepath, 'r') as infile:
			for resno, line in enumerate(infile, 1):
				if line.strip():
					line = line.strip("\n\t").split("\t")
					siteid = "%s %s" % (line[0], line[1].lstrip())
					yield (resno, siteid, tuple(line[2:]))


class MgRNAProfile(object):

	"""

	locations: a list of 5-tuples, in the format as follows
		(pdbid, chainid, locus, locuscoos, metalid), where
		locus is the one-indexed residue number in the MSA
			for the reference nucleotide used as a frame
		locuscoos are the coordinates relative to the frame
			calculated according to NucleotideFrame
		metalid is the siteid for the metal ion, ending in
			a '!' if unreliable magnesium ion,
			a '*' if non-RNA-interacting magnesium ion,
			a ' ' if a reliable magnesium ion, or other ion.
	"""

	def __init__(self, filepath, alignment, name):
		self.path = filepath
		self.name = name
		self.minresno = 0 # max resno for metal id written as of now
		self.model = 0
		self.atomno = count(1)
		self.alignment = alignment
		self.new()

	@classmethod
	def create(cls, pdbid, chainid, msafile, structuredir, filepath):
		dataset = MoleculeDataset("*", "*", msafile=msafile)
		structure, structurefile = dataset.get_structure(structuredir, pdbid)
		alignment = dataset.align(pdbid, chainid, structure)
		return cls(filepath, alignment, "%s:%s" % (pdbid, chainid))

	def new(self):
		lines = ["MODEL% 9d%s" % (self.model, " "*66)]
		for locus, residue in enumerate(self.alignment, 1):
			if residue:
				residue.chain, residue.resno, residue.insco = "A", locus, ""
				atmno = lambda atom: residue.atoms[atom].atmno
				for atom in sorted(residue.atoms, key=atmno):
					residue.atoms[atom].atmno = self.atomno.next()
					lines.append(residue.atoms[atom].to_pdb_record(residue))
		with open(self.path, 'w') as outfile:
			outfile.write("%s\n" % "\n".join(lines))

	def finish(self):
		with open(self.path, 'a') as outfile:
			outfile.write("ENDMDL%s\n" % (" "*74,))

	def update(self, metalcoos, order):
		if metalcoos:
			byresno = dict()
			for metalid, cooslist in metalcoos.items():
				resno = self.minresno + order.index(metalid.strip("*! "))
				byresno[resno] = (metalid, cooslist, resno)
			lines = list()
			for resno in sorted(byresno):
				model_expected = (resno / 100000) + 1
				if model_expected != self.model:
					self.model = model_expected
					self.atomno = count(1)
					lines.append("ENDMDL%s" % (" "*74,))
					lines.append("MODEL% 9d%s" % (self.model, " "*66))
				lines.append(self.metal_pdb_record(*(byresno[resno])))
			with open(self.path, 'a') as outfile:
				outfile.write("%s\n" % "\n".join(lines))
			self.minresno += len(order)

	def metal_pdb_record(self, metalid, cooslist, resno):
		atmno = self.atomno.next() % 100000
		coos, rmsd = self.agreement(cooslist)
		residue = Residue.from_string(metalid)
		residue.resno = resno % 10000
		residue.chain = str((resno / 10000) % 10)
		residue.insco = ""
		#
		atomname = "OS" if residue.name == "OHX" else\
		           "IR" if residue.name == "IRI" else residue.name
		atom = Atom(atomname, residue.name, atmno, atomname, coos, None, rmsd)
		return atom.to_pdb_record(residue)

	def agreement(self, cooslist):
		full = len(cooslist)
		half = full / 2.0
		cooslist = np.array(cooslist)
		meancoos = cooslist.mean(axis=0)
		mean_rmsd = np.sqrt(np.square(cooslist - meancoos).sum(axis=1).mean())
		return meancoos, mean_rmsd


def get_all_leaves(linkarray, node, verbose=False):
	size = len(linkarray) + 1
	queue = {node}
	nonleaves = set()
	#
	def check_nonleaves(nonleaves):
		nonleaves.clear()
		nonleaves.update(filter(lambda node: node >= size, queue))
		return bool(len(nonleaves))
	#
	while check_nonleaves(nonleaves):
		for nonleaf in map(int, nonleaves):
			queue.update(linkarray[nonleaf-size][:2])
		queue -= nonleaves
	if verbose:
		print(node, queue)
	return sorted(map(int, queue))


class MetalCluster(object):

	"""
	Minimum Mg-Mg distance in PDB-XRAY is 1.135 A, between
	[MG]1677:AA and [MG]1788:AA in the structure 4v9a.
	"""

	def __init__(self, ligandsfile, profilefile, chainscount):
		self.coos = list()
		self.rmsd = list()
		self.sites = list()
		self.sitex = dict()
		self.length = None
		self.ligands = list()
		self.clusters = dict()
		self.min10 = int(np.ceil(max(2, chainscount * 0.1)))
		self.min01 = int(np.ceil(self.min10 * 0.1))
		#
		self.read_ligands(ligandsfile)
		self.read_profile(profilefile)

	def read_ligands(self, ligandsfile):
		for resno, siteid, ligands in LigandsFile.read(ligandsfile):
			self.sites.append(siteid)
			self.sitex[siteid] = resno-1
			self.ligands.append(ligands)
		self.length = len(self.sites)

	def read_profile(self, profilefile):
		atom = dict()
		nanvalued = 0
		with StructureFile(profilefile, None) as manager:
			structure = manager.extract_structure(breaks=False, multimodel=True)
		begex, endex, _ = structure.chains["A"]
		for residue in structure.residues[endex:]:
			model, chain, resno = residue.model, residue.chain, residue.resno
			sitex = (model - 1)*100000 + int(chain) * 10000 + resno
			atom[sitex] = list(residue.atoms.values())[0]
		for sitex in range(self.length):
			if sitex in atom:
				self.coos.append(atom[sitex].coos)
				self.rmsd.append(atom[sitex].bfac)
			else:
				self.coos.append(np.array([np.nan, np.nan, np.nan]))
				self.rmsd.append(np.nan)
				nanvalued += 1
		self.coos = np.array(self.coos)
		self.rmsd = np.array(self.rmsd)
		print("% 10d\t% 10d" % (self.length, nanvalued))

	def cluster_sites(self):
		cutoff = 3.08
		self.exclude_unclusterable(cutoff)
		self.print_check(0)
		self.agglomerative_clusters(cutoff)
		self.print_check(1)
		self.knowledge_based_unmerge(2.08)
		self.print_check(2)
		self.remove_small_clusters()
		self.print_check(3)

	def print_check(self, message):
		singletons, clusters = self.singletons()
		print("%d\t%d\t%s" % (len(sum(clusters.values(), [])), len(clusters), message))
		return len(clusters)

	def singletons(self):
		singletons = set()
		clusters = self.clusters.copy()
		for cluster in set(self.clusters) - {0}:
			if len(clusters[cluster]) < 2:
				singletons.add(clusters.pop(cluster)[0])
		return singletons, clusters

	def exclude_unclusterable(self, cutoff):
		components = dict()
		logtime = time.time()
		for sitex in xrange(self.length):
			for key, component in components.items():
				if sitex in component:
					known = component
					break
			else:
				known = set()
			#
			dists = norm(self.coos[sitex:]-self.coos[sitex], axis=1)
			incidence = set(sitex + np.flatnonzero(dists <= cutoff))
			# <= automatically excludes all nan-valued sites,
			# since nan < x or nan > x is False for all x
			if not len(incidence - known):
				continue
			#
			addedto = list()
			for key in sorted(components):
				if len(components[key] & incidence):
					addedto.append(key)
					components[key].update(incidence)
			for key in addedto[1:]:
				components[addedto[0]].update(components.pop(key))
			if not addedto:
				components[sitex+1] = incidence
		coverage = set().union(*components.values())
		#
		badrmsd = set(np.flatnonzero(self.rmsd > 1.0))
		for key in set(components):
			components[key] -= badrmsd
			if len(components[key]) < self.min10:
				components.pop(key)
			else:
				sitexes = np.array(sorted(components[key]))
				exclude = self._cluster_breaking(sitexes, cutoff)
				components[key] -= set(exclude)
		#
		coverage = set().union(*components.values())
		self.clusters = {key: sorted(components[key]) for key in components}
		self.clusters[0] = sorted(set(range(self.length)) - coverage)
		self.renumber_clusters()

	def _cluster_breaking(self, sitexes, cutoff):
		"""
		"""
		size = len(sitexes)
		# identifying clusters of average link > cutoff
		linkarray = linkage(self.coos[sitexes], method='average')
		indices = np.flatnonzero(linkarray[:, 2] > cutoff)
		# clusters with average link > cutoff are parents of other nodes,
		# so the direct children of each cluster is noted in "nodes"
		nodes = np.array(linkarray[indices][:, (0, 1)].flatten(), dtype=int)
		# all those clusters that do not have any direct children with
		# average link <= cutoff are removed, so only clusters whose
		# children have average linkage <= cutoff remain in "nodes"
		nodes = np.array(list(set(nodes - size) - set(indices))) + size
		# "exclude" clusters smaller than "min01", starting with leaves
		exclude = set(nodes[nodes < size])
		for node in set(nodes) - exclude:
			if linkarray[node - size][3] < self.min01:
				exclude.update(get_all_leaves(linkarray, node))
		exclude = sorted(set(sitexes[sorted(exclude)]))
		return exclude

	def remove_small_clusters(self):
		for cluster in set(self.clusters) - {0}:
			if len(self.clusters[cluster]) < self.min10:
				self.clusters[0] += self.clusters.pop(cluster)
		self.renumber_clusters()

	def knowledge_based_unmerge(self, cutoff):
		LARGE_FLOAT = 1000000.
		new_cluster = count(len(self.clusters)+1)
		for cluster in sorted(set(self.clusters) - {0}):
			merged = self._merged_clusters(self.clusters[cluster], cutoff)
			if not merged:
				continue
			sitexes = np.array(self.clusters.pop(cluster))
			matrix = squareform(pdist(self.coos[sitexes]))
			size = len(sitexes)
			for pair in merged:
				matrix[pair] = matrix[pair[::-1]] = LARGE_FLOAT
			linkarray = linkage(squareform(matrix), method='average')
			parents = np.flatnonzero(linkarray[:, 2] > cutoff)
			kids = np.array(linkarray[parents][:, (0, 1)].flatten(), dtype=int)
			kids = np.array(list(set(kids - size) - set(parents))) + size
			for node in kids:
				new_sitexes = list(sitexes[get_all_leaves(linkarray, node)])
				self.clusters[next(new_cluster)] = new_sitexes
		self.renumber_clusters()

	def _merged_clusters(self, sitexes, cutoff):
		merged = list()
		sitexes = np.array(sitexes)
		siteids = [self.sites[sitex][:-1] for sitex in sitexes]
		chains = np.array([re.sub(" .+:", ":", each) for each in siteids])
		if len(set(chains)) != len(sitexes):
			for chain in set(chains):
				repeating = np.flatnonzero(chains == chain)
				if len(repeating) > 1:
					for index1, index2 in combinations(repeating, 2):
						coos1, coos2 = self.coos[sitexes[[index1, index2]]]
						if norm(coos1 - coos2) > cutoff:
							merged.append((index1, index2))
		return merged

	def agglomerative_clusters(self, cutoff):
		disappeared = set()
		new_cluster = count(1)
		clusters = self.clusters.copy()
		self.clusters = {0: clusters[0]}
		for component in sorted(set(clusters) - {0}):
			sitexes = np.array(clusters[component])
			size = len(sitexes)
			linkarray = linkage(self.coos[sitexes], method='average')
			# all nodes with average link > cutoff are called "parents"
			parents = np.flatnonzero(linkarray[:, 2] > cutoff)
			if not len(parents): # no parents => all nodes have link < cutoff
				self.clusters[next(new_cluster)] = list(sitexes)
				continue
			# the direct children of each parent are noted in "kids"
			kids = np.array(linkarray[parents][:, (0, 1)].flatten(), dtype=int)
			# if a kid is also a parent => kid also has average link > cutoff
			# kids-who-are-parents are removed => all kids have link <= cutoff
			kids = np.array(list(set(kids - size) - set(parents))) + size
			for node in kids:
				new_sitexes = list(sitexes[get_all_leaves(linkarray, node)])
				self.clusters[next(new_cluster)] = new_sitexes
		self.renumber_clusters()

	def renumber_clusters(self):
		singletons, clusters = self.singletons()
		order = lambda cluster: (-len(clusters[cluster]), clusters[cluster][0])
		order = sorted(set(clusters) - {0}, key=order)
		clusters = {0: self.clusters.pop(0)}
		for renumbered, cluster in enumerate(order, 1):
			clusters[renumbered] = self.clusters.pop(cluster)
		for cluster, sitex in enumerate(sorted(singletons), max(clusters)+1):
			clusters[cluster] = [sitex]
		self.clusters = clusters

	def shared_ligands(self):
		ligands = dict()
		for cluster in set(self.clusters) - {0}:
			sitexes = self.clusters[cluster]
			mgrna = [sitex for sitex in sitexes if self.sites[sitex][-1]!="*"]
			mgrna = [sitex for sitex in mgrna if self.sites[sitex][:4]=="[MG]"]
			valid = [sitex for sitex in sitexes if self.sites[sitex][-1]==" "]
			valid = [sitex for sitex in valid if self.sites[sitex][:4]!="[NA]"]
			ligands[cluster] = self.frequent_ligands(valid)
			ligands[cluster] |= self.frequent_ligands(mgrna)
		return ligands

	def frequent_ligands(self, sitexes):
		ligands = set()
		counts = defaultdict(lambda: 0)
		half = len(sitexes) / 2.0
		for sitex in sitexes:
			for ligand in self.ligands[sitex]:
				counts[ligand] += 1
		ligands = {ligand for ligand, count in counts.items() if count >= half}
		return ligands

	def write_clusters(self, clusterfile):
		clusterfile = ClusterFile(clusterfile)
		singletons, clusters = self.singletons()
		self.clusters = clusters
		self.clusters[0] = sorted(list(clusters.pop(0)) + list(singletons))
		clusterfile.update(self.clusters, self.sites, self.shared_ligands())


class ClusterFile(object):

	LIGAND_HEADERS = ("Reliable", "Frequent", "Uncommon")
	Representative = ("structurefile", "chainid")

	def __init__(self, filepath):
		self.path = filepath

	def update(self, clusters, sites, ligands):
		lines = list()
		for cluster in clusters:
			sitexes = clusters[cluster]
			lines.append("Cluster:\t%d\t%d\n" % (cluster, len(sitexes)))
			if cluster != 0:
				atoms = list(ligands[cluster])
				atoms.sort(key=lambda atom: int(atom.split(".")[0]))
				lines.append("\tLigands:\t%s\n" % " ; ".join(atoms))
			for sitex in sorted(clusters[cluster]):
				lines.append("\t\t% 5d\t% 20s\n" % (sitex, sites[sitex]))
			lines.append("\n")
		with open(self.path, 'w') as outfile:
			outfile.write("".join(lines))

	def read_cluster_names(self):
		clusters = set()
		with open(self.path, 'r') as infile:
			for line in infile:
				if line.startswith("Cluster:"):
					clusters.add(int(line.split()[1]))
		return clusters

	def read_cluster_sites(self, cluster):
		sites = list()
		include = False
		with open(self.path, 'r') as infile:
			for line in infile:
				if line.startswith("Cluster:"):
					include = (int(line.split()[1]) == cluster)
				elif line.startswith("\t\t") and include:
					line = line.strip("\t\n").split("\t")
					resno, siteid = int(line[0]), line[1].lstrip()
					sites.append((resno, siteid))
		return sites

	def read_all_clusters_sites(self):
		current = None
		clusters = defaultdict(list)
		with open(self.path, 'r') as infile:
			for line in infile:
				if line.startswith("Cluster:"):
					current = int(line.split()[1])
				elif line.startswith("\t\t"):
					line = line.strip("\t\n").split("\t")
					resno, siteid = int(line[0]), line[1].lstrip()
					clusters[current].append((resno, siteid))
		return clusters

	def read_cluster_ligands(self):
		bycluster = dict()
		with open(self.path, 'r') as infile:
			for line in infile:
				if line.startswith("Cluster:"):
					cluster = int(line.split()[1])
					if cluster == 0:
						continue
					nextline = next(infile).strip().split("\t")
					if nextline[0] == "Ligands:" and len(nextline) > 1:
						ligands = list()
						for ligand in nextline[1].split(" ; "):
							frame, atom = ligand.split(".")
							ligands.append((int(frame), atom))
						bycluster[cluster] = ligands
		return bycluster


class ClusterViewer(object):
	CLUSTER_LABEL = " \n ".join(['set echo c%d {%.3f %.3f %.3f}',
						'echo "%d"', 'color echo %s', 'font echo 12 bold'])
	SITE_LABEL = " \n ".join(['select [MG]%d:"#".%s', 'label %s'])

	def __init__(self, profilefile, clusterfile):
		self.structure = profilefile
		self.clusterfile = ClusterFile(clusterfile)
		self.clusters = self.clusterfile.read_cluster_names()
		#
		colours = Plotter.generate_colours(len(self.clusters) - 1)
		colours = map(lambda c: "[%s]" % c.upper().replace("#", "x"), colours)
		self.colour = dict(zip(sorted(set(self.clusters) - {0}), colours))
		self.colour[0] = "[x000000]"

	def interactive(self):
		commands = "Commands:\n\tView All\n\tView 1\n\tQuit"
		print("Clusters: %d - %d" % (min(self.clusters), max(self.clusters)))
		while True:
			answer = raw_input().strip()
			if answer.startswith("View All"):
				self.display_all()
			elif answer.startswith("View "):
				clusters = set(map(int, answer.split()[1:]))
				clusters &= self.clusters
				if len(clusters):
					self.display_clusters(clusters)
			elif answer.startswith("Quit"):
				break
			else:
				print(commands)

	@staticmethod
	def get_selector(sitex):
		resno = (sitex % 10000)
		chain = (sitex / 10000) % 10
		model = (sitex / 100000) + 1
		return "%s:%s/%s" % (resno, chain, model)

	def display_clusters(self, clusters):
		script = list()
		selectors = list()
		for cluster in clusters:
			selection = list()
			for sitex, siteid in self.clusterfile.read_cluster_sites(cluster):
				selector = self.get_selector(sitex)
				script.append('select %s; label "%s"' % (selector, siteid))
				selection.append(selector)
			script.append("select %s" % " or ".join(selection))
			script.append("color %s" % self.colour[cluster])
			selectors += selection
		selectors = " or ".join(selectors)
		script += ['select %s' % selectors, 'color label black', 'font label 9',
					'hide not selected', 'zoom 0', # 'display all',
					'select :A and within(15.0, true, %s)' % selectors,
					'display selected or %s' % selectors,
					'cartoon 0.1', 'wireframe 0.1', 'spacefill 0.15',
					'select :A and not selected', 'cartoon 0', #'wireframe 0.01', 'spacefill 0.05',
					'select not :A', 'spacefill 0.1', 'wireframe 0']
		self._display(script)

	def display_all(self):
		script = list()
		for cluster in self.clusters:
			selection = list()
			for sitex, siteid in self.clusterfile.read_cluster_sites(cluster):
				selection.append(self.get_selector(sitex))
			script.append("select %s" % " or ".join(selection))
			script.append("color %s" % self.colour[cluster])
			if cluster:
				script.append('label "%s"' % cluster)
		script += ['select :"1"', 'font label 9', 'zoom 0',
					'select :A', 'cartoon 0.1', 'wireframe 0.1', 'spacefill 0.15',
					'select not :A', 'spacefill 0.1', 'wireframe 0']
		self._display(script)

	def _display(self, script):
		script += ["color background white", "select none"]
		prescript = ["load MODELS {0 -1 1} %s" % self.structure,
					 "frame all", "display *"]
		script = "; \n".join(prescript + script)
		command = ["jmol", "-I"]
		with open(os.devnull, 'w') as NULL:
			process = Popen(command, stdin=PIPE, stdout=PIPE, stderr=PIPE)
			process.stdin.write(script)


def example(structuredir, resolutionsfile, msafile):
	dataset = MoleculeDataset("*", "*", msafile=msafile)
	ordered = dataset.representatives(structuredir, resolutionsfile)
	for entry in ordered:
		toprint = ["%.2f" % entry[0], entry[1]]
		toprint += ["% 2s" % chainid for chainid in entry[1:]]
		print("\t".join(toprint))

def locate(structuredir, msafile, pdbid, chainid, outfile, ligandsfile, max_da=3.8):
	profile = MgRNAProfile.create(pdbid, chainid, msafile, structuredir, outfile)
	MetalLocator(profile, ligandsfile, max_da).locate(structuredir, msafile)

def cluster(ligandsfile, profilefile, clusterfile, chainscount):
	clusterer = MetalCluster(ligandsfile, profilefile, chainscount)
	clusterer.cluster_sites()
	clusterer.write_clusters(clusterfile)

def view(profilefile, clusterfile):
	viewer = ClusterViewer(profilefile, clusterfile)
	viewer.interactive()

#

def parse_args():
	parser = argparse.ArgumentParser()
	commands = parser.add_subparsers(dest="command")
	#
	ligands = commands.add_parser('ligands')
	ligands.add_argument("structuredir", metavar="STRUCTURE-FOLDER")
	ligands.add_argument("sitesdir", metavar="CONTEXT-FOLDER")
	ligands.add_argument("chainsdir", metavar="CHAINS-FOLDER")
	ligands.add_argument("sequence", metavar="SEQUENCE")
	ligands.add_argument("molecule", metavar="MOLECULE")
	ligands.add_argument("organism", metavar="ORGANISM")
	ligands.add_argument("--ids", metavar="IDS-FILE", default=None)
	ligands.add_argument("--atom", metavar="ATOM-NAME", default=None)
	ligands.add_argument("--position", metavar="BASE-POSITION", default=None, type=int)
	#
	sequences = commands.add_parser('sequences')
	sequences.add_argument("structuredir", metavar="BPFIND-OUTPUT-FOLDER")
	sequences.add_argument("chainsdir", metavar="CHAINS-FOLDER")
	sequences.add_argument("molecule", metavar="MOLECULE")
	sequences.add_argument("organism", metavar="ORGANISM")
	sequences.add_argument("output", metavar="OUTPUT-FILE")
	sequences.add_argument("--ids", default=None)
	sequences.add_argument("--msa", metavar="MSA-OUTPUT-FILE")
	#
	check = commands.add_parser('check')
	check.add_argument("msafile", metavar="MSA-OUTPUT")
	check.add_argument("structuredir", metavar="STRUCTURE-FOLDER")
	#
	locate = commands.add_parser('locate')
	locate.add_argument("structuredir", metavar="STRUCTURE-FOLDER")
	locate.add_argument("msafile", metavar="MSA-OUTPUT")
	locate.add_argument("pdbid", metavar="REPRESENTATIVE-STRUCTURE")
	locate.add_argument("chainid", metavar="REPRESENTATIVE-CHAIN")
	locate.add_argument("outfile", metavar="OUTPUT-STRUCTURE-FILE")
	locate.add_argument("ligandsfile", metavar="LIGANDS-FILE")
	locate.add_argument("--max-da", default="3.8", type=float)
	#
	view = commands.add_parser('view')
	view.add_argument("profilefile", metavar="PROFILE-FILE")
	view.add_argument("clusterfile", metavar="CLUSTER-FILE")
	#
	example = commands.add_parser('example')
	example.add_argument("structuredir", metavar="STRUCTURE-FOLDER")
	example.add_argument("resolutionsfile", metavar="RESOLUTIONS-FILE")
	example.add_argument("msafile", metavar="MSA-OUTPUT")
	#
	cluster = commands.add_parser('cluster')
	cluster.add_argument("ligandsfile", metavar="LIGANDS-FILE")
	cluster.add_argument("profilefile", metavar="PROFILE-FILE")
	cluster.add_argument("clusterfile", metavar="CLUSTER-FILE")
	cluster.add_argument("chainscount", metavar="CHAINS-COUNT", type=int)
	#
	return parser.parse_args()


COMMAND = {"ligands": EngraphedMetal.ligands, "sequences": sequences, "check": check,
			"locate": locate, "cluster": cluster, "view": view, "example": example}

if __name__ == "__main__":
	args = vars(parse_args())
	COMMAND[args.pop("command")](**args)
