from __future__ import print_function

import os
import sys

from magracarna.engraph import StructureGraphFile
from magracarna.metalsite import RNAMetalContextGraph, RNAMgSiteEngrapher, RNAMgSiteList
from magracarna.motifs import MotifAssigner, MotifAssignments
ISPARTOF = RNAMetalContextGraph.EDGELABELS[RNAMetalContextGraph.ISPARTOF]


def rewrite_phosphate_edges(graph):
    for nodeid, nlabel in graph.nodes.items():
        if nlabel.name == graph.Oph.name:
            resnodeid = nodeid.split(".")[0]
            graph.add_edge(resnodeid, nodeid, ISPARTOF)
    return graph

def rewrite_reliability(graph, reference):
    name = graph.name.strip("!")
    graph.name = "%s%s" % (name, reference[name])

def rewrite_basepair_node(graph):
    mapping = dict()
    for nodeid in graph.nodes.copy():
        if "=" in nodeid:
            orient = graph.nodes[nodeid].name[-1]
            mapping[nodeid] = "%s%s" % (nodeid, orient)
            graph.nodes[mapping[nodeid]] = graph.nodes.pop(nodeid)
    for edge in graph.edges.copy():
        if edge[0] in mapping or edge[1] in mapping:
            graph.edges.remove(edge)
            nodeid = mapping[edge[0]] if edge[0] in mapping else edge[0]
            edonid = mapping[edge[1]] if edge[1] in mapping else edge[1]
            graph.edges.add((nodeid, edonid, edge[2]))

def rewrite_graphs(engrapher, infile, outfile, reference=None):
    graphs = RNAMgSiteList.load(infile)
    mgsitelist = RNAMgSiteList(outfile, engrapher, encoded=True)
    mgsitelist.graphs = graphs
    for graph in mgsitelist.graphs:
        rewrite_basepair_node(graph)
        metnode = (graph.name.split(" ")[1].strip("!"), "MG")
        mgsitelist.engrapher.graph = graph
        mgsitelist.engrapher._node = metnode
        mgsitelist.graphs.append(mgsitelist.engrapher.encode_graph())
    mgsitelist.dump()

def _read_reliablity(infile):
    reference = dict()
    with open(infile, 'r') as reffile:
        for line in reffile:
            line = line.strip("\n").split("\t")
            reference[line[0]] = line[5]
    return reference

def _read_sites(reffile):
    mapping = dict()
    with open(reffile, 'r') as sitesfile:
        for line in sitesfile:
            line = [part.strip() for part in line.split("\t")]
            if len(line) >= 3 and "=" in line[1] and line[1] not in mapping:
                mapping[line[1]] = line[1].strip("CT")
    return {value: key for key, value in mapping.items()}

def rewrite_assigns(assigner, infile, outfile, sitesfile):
    assigns = MotifAssignments.load(assigner, infile, detail=True)
    mapping = _read_sites(sitesfile)
    for graph, casewise in assigns.by_graph.items():
        for case, embeds in casewise.items():
            for embed in embeds:
                for index, partid in enumerate(embed[:]):
                    if "=" in partid:
                        embed[index] = mapping[partid]
    assigns.dump(outfile, detail=True)

def rewrite(mode, infolder, outfolder, *references):
    if mode == "sites":
        reference = _read_reliablity(references[0]) if references else None
        engrapher = RNAMgSiteEngrapher(None)
        for infile, outfile in iterate(infolder, outfolder, "_sites.tsv"):
            rewrite_graphs(engrapher, infile, outfile, reference)
    if mode == "assigns":
        sitesfolder, motiffiles = references[0], _motiffiles(references[1:])
        assigner = MotifAssigner.extract(*motiffiles, mode="all")
        for iteritem in iterate(infolder, outfolder, "_assigns.tsv"):
            pdbid, infile, outfile = iteritem
            sitesfile = os.path.join(sitesfolder, pdbid+"_sites.tsv")
            rewrite_assigns(assigner, infile, outfile, sitesfile)

def _motiffiles(motifsargs):
    motifsfiles = list()
    for motifsarg in motifsargs:
        if os.path.isdir(motifsarg):
            folderfiles = list()
            for filename in os.listdir(motifsarg):
                if filename.split("_")[0].isdigit():
                    order_by = int(filename.split("_")[0])
                    filepath = os.path.join(motifsarg, filename)
                    folderfiles.append((order_by, filepath))
            motifsfiles += [filepath for _, filepath in sorted(folderfiles)]
        else:
            motifsfiles.append(motifsarg)
    return motifsfiles


def iterate(infolder, outfolder, extension):
    for file in sorted(os.listdir(infolder)):
        if file.endswith(extension):
            pdbid = file.split("_")[0]
            infile = os.path.join(infolder, file)
            outfile = os.path.join(outfolder, file)
            print(pdbid)
            yield (pdbid, infile, outfile)
            

if __name__ == "__main__":
    mode = "assigns"
    infolder = sys.argv[1]
    outfolder = sys.argv[2]
    rewrite(mode, infolder, outfolder, *sys.argv[3:])
