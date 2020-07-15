from __future__ import print_function

import os
import sys
from magracarna.metalsite import RNAMgSiteEngrapher, RNAMgSiteList

# folder = os.path.abspath(sys.argv[1])

# for filename in os.listdir(folder):
# 	if filename.endswith("_sites.tsv"):
# 		filepaths.append(os.path.join(folder, filename))
rnagfs = list()
for sitesfile in sys.argv[1:]:
	rnagfs.append(RNAMgSiteList.decode_file(sitesfile)[0])

common_edges = set(rnagfs[0].edges)
for rnagf in rnagfs[1:]:
	common_edges &= rnagf.edges
for rnagf, filepath in zip(rnagfs, sys.argv[1:]):
	nodeids = rnagf.nodes.keys()
	rnagf.edges -= common_edges
	hbedges = set()
	for edge in rnagf.edges:
		if str(edge[2]) == "HB":
			res1 = edge[0].split(".")[0]
			res2 = edge[1].split(".")[0]
			if res1 == res2:
				hbedges.add(edge)
	print(len(hbedges))
	for hbedge in hbedges.copy():
		for nodeid in hbedge[:2]:
			resnodeid = nodeid.split(".")[0]
			alledges = rnagf.incident_edges(nodeid)
			alledges = list(set(alledges) - hbedges)
			if len(alledges) == 1:
				lastedge = alledges[0]
				if lastedge[:2] == (resnodeid, nodeid)\
				  or lastedge[:2] == (nodeid, resnodeid):
					if str(lastedge[2]) == "I":
						hbedges.add(alledges[0])
	print(len(hbedges))
	rnagf.edges -= hbedges
	# 
	# hohedges = list()
	# for edge in rnagf.edges:
	# 	if edge[0].startswith("[HOH]") and edge[1].startswith("[HOH]")\
	# 	  and str(edge[2]) == "I":
	# 		hohedges.append(edge)
	# 		print(edge)
	# print(len(hohedges))
	# # 
	for nodeid in nodeids:
		if not rnagf.incident_edges(nodeid):
			del rnagf.nodes[nodeid]
	print(len(rnagf.nodes), len(rnagf.edges))
	engrapher = RNAMgSiteEngrapher(None, True)
	engrapher.graph = rnagf
	engrapher._node = ("[MG]3001:AA", "MG")
	print(filepath)
	print(engrapher.encode_graph())
	print("\n\n")

