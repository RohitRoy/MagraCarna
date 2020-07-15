import os
import sys

from collections import defaultdict

from ..magracarna.motifs import MotifAssignments
from ..magracarna.engraph import RNAGraph
from ..utilities.args import ArgParser
from ..magracarna.aminoacids import aminoacids as AMINOACIDS
from ..magracarna.nucleotides import backbone as BACKBONE

HETERO = ("1Oh_Graph", "1Nh_Graph", "Hetero-Outer")
NUCLEOTIDES = ("1Oph_Graph", "1Or_Graph", "1Ob_Graph", "1Nb_Graph", "1Pout_Graph", "1Bout_Graph", "1Rout_Graph")
RNANUCS = ("A", "C", "G", "U")

def add_to_list(resset, siteid, line):
	mode = line[3]
	resname = line[4].split("]")[0].strip("[")
	atomname = line[4].split(".")[1]
	resset.add((siteid, resname, mode, atomname))

def stats_aa_contact(heteros, onlysites=None):
	sites = defaultdict(set)
	for siteid, resname, _, _ in heteros:
		if resname in AMINOACIDS:
			if resname not in ["ARG", "LYS", "GLY", "HIS", "SER", "ASP", "ASN", "ALA", "GLU", "GLN", "THR"]:
				resname = "Other"
			sites[resname].add(siteid)
	if onlysites:
		for resname in sites:
			sites[resname] &= onlysites
	return {resname: len(sites[resname]) for resname in sites}

def stats_aaa_contact(heteros, onlysites=None):
	sites = defaultdict(set)
	for siteid, resname, mode, atomname in heteros:
		if not resname in AMINOACIDS:
			continue
		if atomname in ("N", "O"):
			sites["aa.%s %s" % (atomname, mode)].add(siteid)
		else:
			sites["%s.%s %s" % (resname, atomname, mode)].add(siteid)
	if onlysites:
		for key in sites:
			sites[key] &= onlysites
	# 		print(key)
	# 		for site in sorted(sites[key]):
	# 			print("\t%s" % site)
	# print('\n')
	return {key: len(sites[key]) for key in sites}

def stats_Naa_contact(heteros, nucleotides, onlysites=None):
	bysiteid = defaultdict(set)
	for siteid, resname, mode, atomname in nucleotides:
		if atomname in RNAGraph._RIBOATOMNAMES+RNAGraph._PHOSATOMNAMES:
			if atomname in ["O3'", "O5'"]:
				atomname = "O3'/O5'"
			if atomname in RNAGraph._PHOSATOMNAMES:
				atomname = "OP1/2/3"
			resname = "N"
		# bysiteid[siteid].add("%s.%s %s" % (resname, atomname, mode))
		bysiteid[siteid].add("%s.%s" % (resname, atomname))
	sites = defaultdict(set)
	for siteid, resname, mode, atomname in heteros:
		if (resname not in AMINOACIDS) or (onlysites and siteid not in onlysites):
			continue
		if resname == "ARG" and atomname in ("NH1", "NH2"):
			atomname = "NH1/2"
		# key = "%s.%s %s" % (resname, atomname, mode)
		key = "%s.%s" % (resname, atomname)
		for subkey in bysiteid[siteid]:
			sites["%s %s" % (subkey, key)].add(siteid)
	return {key: len(sites[key]) for key in sites}

def stats_Na_aaa_contact(heteros, nucleotides, onlysites=None):
	bysiteid = defaultdict(set)
	for siteid, resname, mode, atomname in nucleotides:
		if resname not in RNANUCS:
			continue
		if atomname in RNAGraph._RIBOATOMNAMES+RNAGraph._PHOSATOMNAMES:
			if atomname in ["O3'", "O5'"]:
				atomname = "O3'/O5'"
			if atomname in RNAGraph._PHOSATOMNAMES:
				atomname = "OP1/2/3"
			resname = "N"
		bysiteid[siteid].add("%s.%s" % (resname, atomname))
	sites = defaultdict(lambda: defaultdict(lambda : set()))
	keysortby = dict()
	subkeys = set()
	for siteid, resname, mode, atomname in heteros:
		if (resname not in AMINOACIDS) or (onlysites and siteid not in onlysites):
			continue
		if resname == "ARG" and atomname in ("NH1", "NH2"):
			atomname = "NH1/2"
		key = "%s.%s" % (resname, atomname)
		for subkey in bysiteid[siteid]:
			sites[key][subkey].add(siteid)
			subkeys.add(subkey)
		keysortby[key] = len(set.union(*sites[key].values()))
	subkeys = sorted(subkeys)
	theformat = "%% %ds" % max(map(len, subkeys))
	print(theformat % 5)
	header = [" "*9] + [theformat % subkey for subkey in subkeys]
	print("\t".join(header))
	for key in sorted(sites, key=lambda key: -keysortby[key]):
		row = ["% 9s" % key]
		for subkey in subkeys:
			row.append(theformat % len(sites[key][subkey]))
		print("\t".join(row))
	# return {key: len(sites[key]) for key in sites}

def combo_stats(*data):
	# resnames = set()
	data = [defaultdict(lambda: 0, datadict) for datadict in data]
	resnames = set.union(*map(set, data))
	sortby = lambda name : tuple([datadict[name] for datadict in data])
	resnames = sorted(resnames, reverse=True, key=sortby)
	for resname in resnames:
		row = [datadict[resname] if resname in datadict else 0 for datadict in data]
		if sum(row) != 0:
			row = ["% 21s" % value for value in [resname]+row]
			print("\t".join(row))

def aminoacids(assignsdir, idsfiles):
	statslist = list()
	for idsfile in idsfiles:
		idslist = ArgParser.read_idslist(idsfile)
		heteros = set()
		nucleotides = set()
		for pdbid in idslist:
			filepath = ArgParser.get_sources(assignsdir, [pdbid], "_assigns.tsv")[0]
			assignsiter = MotifAssignments._assigns_reader(filepath)
			for header, _, motif, line in assignsiter:
				if not header:
					line = line.split()
					siteid = " ".join(line[:2])
					if not siteid.endswith("!"):
						if motif in HETERO:
							add_to_list(heteros, siteid, line)
						if motif in NUCLEOTIDES:
							add_to_list(nucleotides, siteid, line)
		nucset = {siteid for siteid, resname, _, _ in nucleotides if resname in RNANUCS}
		# nomodcount = stats_aa_contact(heteros, onlysites=nucset)
		# nomodcount = stats_aaa_contact(heteros, onlysites=nucset)
		nomodcount = stats_Naa_contact(heteros, nucleotides, onlysites=nucset)
		# stats_Na_aaa_contact(heteros, nucleotides, onlysites=nucset)
		statslist.append(nomodcount)
	combo_stats(*statslist)

if __name__ == "__main__":
	aminoacids(sys.argv[1], sys.argv[2:])
