#!/usr/bin/python

from __future__ import print_function

import os
import sys
import argparse

from decimal import Decimal

from ..utilities.args import ArgParser
from ..magracarna.engrid import HydrogenBondFinder as HBFinder
from ..magracarna.engrid import SecondaryStructureManager as SSManager
from ..magracarna.engraph import DFSCode, RNAStructureGraph
from ..magracarna.metalsite import RNAMgSiteEngrapher, RNAMgSiteList,\
                                 RNAMgMotifCounts, RNAMgSiteContextVisualiser,\
                                 RNAMgSiteContextParser, RNAMgMotifDivider,\
                                 RNAMgMotifVectorList
from ..magracarna.motifs import MotifFilesParser, MotifAssigner, CaseReader,\
                              MotifAssignments, MotifCombiner, MotifRedivider

_TYPES = ["sites", "assigns", "counts", "chains"]

def _filepath(filepath, is_structure, to_type):
    assert(to_type in _TYPES)
    delim = "." if is_structure else "_"
    folder = os.path.dirname(filepath)
    filename = delim.join(os.path.basename(filepath).split(delim)[:-1])
    extension = "_"+to_type+(".txt" if to_type == "chains" else ".tsv")
    return os.path.join(folder, filename+extension)

def _get_motifsfiles(args):
    motifsargs = args.pop("motifs")
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
    args["motifsfiles"] = motifsfiles

# 

def bpfind(structurefile, rundir, max_da=HBFinder.MAX_DA, verbose=False):
    hbfinder = HBFinder(max_da=max_da)
    bpfind_run = SSManager.BPFindParams(rundir, False, verbose, hbfinder)
    with SSManager(structurefile, bpfind_run) as manager:
        manager.extract_paired_structure()

#

def _parameters(ignore_H=True, hbonds=False, relaxed=False, invalids=True,
                modbases=True, depth=3, verbose=False,
                ligdists=RNAMgSiteEngrapher.LIGDISTS,
                max_da=HBFinder.MAX_DA, min_da=HBFinder.MIN_DA,
                max_ha=HBFinder.MAX_HA, min_ha=HBFinder.MIN_HA,
                max_dha=HBFinder.MAX_DHA, min_dha=HBFinder.MIN_DHA):
    # hydrogen bonding parameters
    hbfinder = HBFinder(max_da=max_da, max_ha=max_ha, max_dha=max_dha,
                        min_da=min_da, min_ha=min_ha, min_dha=min_dha,
                        ignore_H=ignore_H)
    # sites-graph parameters
    engrapher = RNAMgSiteEngrapher(hbfinder, depth=depth, hbonds=hbonds,
                                   ligdists=ligdists, relaxed=relaxed,
                                   invalids=invalids, modbases=modbases,
                                   verbose=verbose)
    return hbfinder, engrapher


def sites(structurefile, rundir=None, **kwargs):
    hbfinder, engrapher = _parameters(**kwargs)
    sitesfile = _filepath(structurefile, True, "sites")
    chainsfile = _filepath(structurefile, True, "chains")
    mgsitelist = RNAMgSiteList(sitesfile, engrapher, encoded=False)
    if rundir:
        bpfind_run = SSManager.BPFindParams(rundir, False, False, hbfinder)
    else:
        bpfind_run = None
    with SSManager(structurefile, bpfind_run) as manager:
        mgsitelist.add(manager.extract_paired_structure())
        if rundir:
            structurefile = manager.structfile.write_to
            sitesfile = _filepath(structurefile, True, "sites")
            mgsitelist.filepath = sitesfile
    RNAMgSiteContextParser.create_chainsfile(mgsitelist.graphs, chainsfile)
    mgsitelist.dump()

#

def assigns(sitesfile, motifsfiles, mode="first", reassign=False,
            detail=False, sort_by=RNAMgMotifCounts._SORTERS):
    rnagfs = RNAMgSiteList.load(sitesfile)
    assigner = MotifAssigner.extract(*motifsfiles, mode=mode)
    assignsfile = _filepath(sitesfile, False, "assigns")
    if reassign and os.path.exists(assignsfile):
        assignments = MotifAssignments.load(assigner, assignsfile, detail)
    else:
        assignments = MotifAssignments(assigner)
    assignments.add_from(rnagfs)
    assignments.dump(assignsfile, detail)
    countsfile = _filepath(sitesfile, False, "counts")
    counts = RNAMgMotifCounts.create_from_assigns(assignments)
    counts.dump(countsfile, sort_by)

#

def counts(assignsfile, countsfile, motifsfiles):
    assigner = MotifAssigner.extract(*motifsfiles, mode="first")
    assigns = MotifAssignments.load(assigner, assignsfile, detail=False)
    RNAMgMotifCounts.create_from_assigns(assigns).dump(countsfile)

# 

def merge_sites(sources, destination, only_valid=False):
    header = RNAMgSiteEngrapher._TYPE
    with open(destination, 'w') as outfile:
        for index, source in enumerate(sources):
            ignore = False
            if index:
                outfile.write("\n\n")
            with open(source, 'r') as infile:
                for line in infile:
                    if only_valid and line.startswith(header):
                        ignore = line.split("\t")[2].endswith("!")
                    if not ignore:
                        outfile.write(line)

def merge_assigns(sources, destination, motifsfiles):
    assigner = MotifAssigner.extract(*motifsfiles, mode="all")
    MotifAssignments.merge_files(assigner, sources, destination)

def merge_counts(sources, destination, motifsfiles):
    assigner = MotifAssigner.extract(*motifsfiles, mode="none")
    RNAMgMotifCounts.merge(assigner, sources, destination)

MERGE_TYPE = {"sites": merge_sites,
              "assigns": merge_assigns,
              "counts": merge_counts}

def merge(filetype, folder, dataset, ids=None, **kwargs):
    extension = "_"+filetype+".tsv"
    sources = ArgParser.get_sources(folder, ids, extension)
    destination = os.path.join(folder, dataset) + extension
    MERGE_TYPE[filetype](sources, destination, **kwargs)

#

def visualise(structure, sitesfile, site=None, assignsfile=None, context=False, waters=False):
    names = [site] if site else None
    rnagfs = RNAMgSiteEngrapher.decode_file(sitesfile, names)
    visualiser = RNAMgSiteContextVisualiser(structure, rnagfs)
    if site:
        visualiser.display(site, context, waters)
    else:
        print("\t".join(sorted(visualiser.site)))
        visualiser.interactive()

#

def divide(name, assigns, outfile, motifsfiles,
           ids=None, show=False, sitesfolder=None,
           only_valid=False, limited=None, restrict=None):
    kwargs = {"limited": limited, "restrict": restrict}
    kwargs["invalids"] = not only_valid
    assigner = MotifAssigner.extract(*motifsfiles, mode="all")
    divider = RNAMgMotifDivider(name, assigner, **kwargs)
    assigns = [assigns]
    if os.path.isdir(assigns[0]):
        assigns = ArgParser.get_sources(assigns[0], ids, "_assigns.tsv")
    divider.add_embeddings(assigns, bool(sitesfolder))
    divider.summarise(outfile)
    if sitesfolder:
        sites = ArgParser.get_sources(sitesfolder, ids, "_sites.tsv")
        if len(assigns) > 1:
            assert(len(assigns) == len(sites))
        divider.add_sequences(sites)
        show = True
    divider.divide(outfile, -1, show, bool(sitesfolder))

def redivide(dividefile, outfile, chainsfolder):
    redivider = MotifRedivider(chainsfolder)
    redivider.redivide(dividefile, outfile)

# 

def formula(report, formula, folder, motifsfiles, 
            ids=None, only_valid=False, **kwargs):
    name = kwargs.pop("name", "formula")
    parser = MotifFilesParser()
    parser.extract(*motifsfiles)
    parser._parse_formulae([(name, formula.split(" or "))])
    sources = ArgParser.get_sources(folder, ids, "_assigns.tsv")
    initargs = (parser.definedin, parser.graphs, parser.formulae, "all")
    assigner = MotifAssigner(*initargs)
    cases = MotifAssignments.quick_formula(name, assigner, sources, None)
    if only_valid:
        cases = filter(lambda case: not case.endswith("!"), cases)
    if report == "assigns":
        MotifAssignments._write_quick_formula(name, cases, kwargs["outfile"])
    elif report == "counts":
        print(len(set(cases)))
        

# 

def combine_assigns(combiner, folder, outfile, motiffile=None, limit=None, ids=None):
    assignsfiles = ArgParser.get_sources(folder, ids, "_assigns.tsv")
    combiner.restrict = limit
    combiner.generate_assigner(motiffile=motiffile, mode="load")
    combiner.generate_assigns(assignsfiles, outfile)

def combine_motifs(combiner, outfile):
    combiner.generate_assigner(motiffile=outfile, mode="dump")

COMBINE_TYPE = {"motifs": combine_motifs, "assigns": combine_assigns}

def combine(name, augend, addend, motifsfiles, filetype, **kwargs):
    base = DFSCode.from_lines(["0\tMG"], RNAStructureGraph.LABELS)
    assigner = MotifAssigner.extract(*motifsfiles, mode="all")
    combiner = MotifCombiner(name, base, augend, addend, assigner)
    COMBINE_TYPE[filetype](combiner, **kwargs)

# 

def vectorise(name, assignsfolder, structurefolder, outfile,
              motifsfiles, ids=None, only_valid=False, mode="all",
              restrict=None, limited=None, division=None):
    kwargs = {"limited": limited, "division": division, "mode": mode}
    kwargs["invalids"] = not only_valid
    assigner = MotifAssigner.extract(*motifsfiles, mode="all")
    vectorlist = RNAMgMotifVectorList(name, assigner, **kwargs)
    structures = ArgParser.get_sources(structurefolder, ids, ".pdb", ".cif")
    assignsfiles = ArgParser.get_sources(assignsfolder, ids, "_assigns.tsv")
    sources = ArgParser.map_sources(structures, assignsfiles)
    formatting = "%% 5d / %d\t%%s" % len(sources)
    for index, (structurefile, assignsfile) in enumerate(sources, 1):
        print(formatting % (index, ArgParser.basename(structurefile, "._")))
        vectorlist.add_vectors_from(assignsfile, structurefile)
    vectorlist.dump_vectors(outfile)

# 

def test(site, contextfile):
    for baseseqs, idseqs in RNAMgSiteEngrapher.extract_sequences(sitesfile, [site]):
        for baseseq, idseq in zip(baseseqs, idseqs):
            print("%s\n\t%s" % (baseseq, " ".join(idseq)))

def parse_args():
    parser = argparse.ArgumentParser()
    commands = parser.add_subparsers(dest="command")
    # 
    bpfind = commands.add_parser('bpfind')
    bpfind.add_argument("structure", metavar="STRUCTURE")
    bpfind.add_argument("--rundir", "-d", default=os.path.realpath("."))
    bpfind.add_argument("--verbose", "-v", action="store_false")
    bpfind.add_argument("--max-da", type=Decimal, default=HBFinder.MAX_DA)
    #
    sites = commands.add_parser('sites')
    sites.add_argument("structure", metavar="STRUCTURE")
    sites.add_argument("--rundir", "-d", default=None)
    sites.add_argument("--has-H", "-H", action="store_true")
    sites.add_argument("--hbonds", "-B", action="store_true")
    sites.add_argument("--relaxed", "-R", action="store_true")
    sites.add_argument("--verbose", "-v", action="store_true")
    sites.add_argument("--only-valid", "-V", action="store_true")
    sites.add_argument("--unmodified", "-U", action="store_true")
    sites.add_argument("--depth", "-D", default=3, type=int)
    LIGDISTS = RNAMgSiteEngrapher.LIGDISTS
    sites.add_argument("--mgo", "-O", type=Decimal, default=LIGDISTS["O"])
    sites.add_argument("--mgn", "-N", type=Decimal, default=LIGDISTS["N"])
    sites.add_argument("--max-da", type=Decimal, default=HBFinder.MAX_DA)
    sites.add_argument("--min-da", type=Decimal, default=HBFinder.MIN_DA)
    sites.add_argument("--max-ha", type=Decimal, default=HBFinder.MAX_HA)
    sites.add_argument("--min-ha", type=Decimal, default=HBFinder.MIN_HA)
    sites.add_argument("--max-dha", type=Decimal, default=HBFinder.MAX_DHA)
    sites.add_argument("--min-dha", type=Decimal, default=HBFinder.MIN_DHA)
    #
    assigns = commands.add_parser('assigns')
    assigns.add_argument("sites", metavar="CONTEXT")
    assigns.add_argument("motifs", nargs="+", metavar="MOTIFS")
    MODES = MotifAssigner.MODES
    assigns.add_argument("--mode", "-m", choices=MODES, default="first")
    assigns.add_argument("--detail", "-d", action="store_true")
    assigns.add_argument("--reassign", "-r", action="store_true")
    SORTERS = RNAMgMotifCounts._SORTERS
    assigns.add_argument("--sort-by", choices=SORTERS, default="dictionary")
    #
    merge = commands.add_parser('merge')
    merge.add_argument("folder")
    merge.add_argument("--ids", default=None)
    filetypes = merge.add_subparsers(dest="filetype")
    merge_sites = filetypes.add_parser("sites")
    merge_sites.add_argument("--only-valid", "-V", action="store_true")
    merge_assigns = filetypes.add_parser("assigns")
    merge_assigns.add_argument("motifs", nargs="+", metavar="MOTIFS")
    merge_counts = filetypes.add_parser("counts")
    merge_counts.add_argument("motifs", nargs="+", metavar="MOTIFS")
    #
    visualise = commands.add_parser('visualise')
    visualise.add_argument("structure", metavar="STRUCTURE")
    visualise.add_argument("sitesfile", metavar="CONTEXT-FILE")
    visualise.add_argument("--assignsfile", metavar="ASSIGNS", default=None)
    visualise.add_argument("--context", "-C", action="store_true")
    visualise.add_argument("--waters", "-W", action="store_true")
    visualise.add_argument("--site", metavar="SITE", default=None)
    #
    divide = commands.add_parser('divide')
    divide.add_argument("name", metavar="MOTIF-GRAPH-NAME")
    divide.add_argument("--limited", default=None, metavar="SITES-FILE")
    divide.add_argument("--restrict", default=None, metavar="MOTIF")
    divide.add_argument("assigns", metavar="ASSIGNS")
    divide.add_argument("--sitesfolder", metavar="CONTEXT-FOLDER")
    divide.add_argument("--ids", default=None, metavar="IDSLIST")
    divide.add_argument("outfile", metavar="OUTPUT")
    divide.add_argument("--show", action="store_true")
    divide.add_argument("--only-valid", "-V", action="store_true")
    divide.add_argument("motifs", nargs="+", metavar="MOTIFS")
    #
    redivide = commands.add_parser('redivide')
    redivide.add_argument("dividefile", metavar="DIVIDE-FILE")
    redivide.add_argument("outfile", metavar="OUTPUT")
    redivide.add_argument("chainsfolder", metavar="CHAIN-DETAILS-FOLDER")
    # 
    formula = commands.add_parser('formula')
    formula.add_argument("formula", metavar="MOTIF-FORMULA")
    formula.add_argument("folder", metavar="FOLDER")
    formula.add_argument("--ids", default=None, metavar="IDSLIST")
    formula.add_argument("--only-valid", "-V", action="store_true")
    formula_report = formula.add_subparsers(dest="report")
    formula_assigns = formula_report.add_parser("assigns")
    formula_assigns.add_argument("name", metavar="MOTIF-NAME")
    formula_assigns.add_argument("outfile", metavar="OUTPUT")
    formula_assigns.add_argument("motifs", nargs="+", metavar="MOTIFS")
    formula_counts = formula_report.add_parser("counts")
    formula_counts.add_argument("motifs", nargs="+", metavar="MOTIFS")
    # 
    combine = commands.add_parser('combine')
    combine.add_argument("name", metavar="MOTIF-NAME")
    combine.add_argument("augend", metavar="MOTIF-GRAPH-1")
    combine.add_argument("addend", metavar="MOTIF-GRAPH-2")
    combine_types = combine.add_subparsers(dest="filetype")
    combine_motifs = combine_types.add_parser("motifs")
    combine_motifs.add_argument("outfile", metavar="OUTPUT")
    combine_motifs.add_argument("motifs", nargs="+", metavar="MOTIFS")
    combine_assigns = combine_types.add_parser("assigns")
    combine_assigns.add_argument("--motiffile", metavar="MOTIF-FILE")
    combine_assigns.add_argument("--limit", default=None, metavar="LIMIT_TO")
    combine_assigns.add_argument("folder", metavar="FOLDER")
    combine_assigns.add_argument("--ids", default=None, metavar="IDSLIST")
    combine_assigns.add_argument("outfile", metavar="OUTPUT")
    combine_assigns.add_argument("motifs", nargs="+", metavar="MOTIFS")
    #
    counts = commands.add_parser('counts')
    counts.add_argument("assignsfile", metavar="ASSIGNS")
    counts.add_argument("countsfile", metavar="OUTPUT")
    counts.add_argument("motifs", nargs='+', metavar="MOTIFS")
    # 
    vectorise = commands.add_parser('vectorise')
    vectorise.add_argument("name", metavar="MOTIF-GRAPH-NAME")
    vectorise.add_argument("assignsfolder", metavar="ASSIGNS-FOLDER")
    vectorise.add_argument("structurefolder", metavar="STRUCTURE-FOLDER")
    vectorise.add_argument("--ids", default=None, metavar="IDSLIST")
    vectorise.add_argument("--only-valid", "-V", action="store_true")
    vectorise.add_argument("--limited", default=None, metavar="SITES-FILE")
    vectorise.add_argument("--division", default=None, metavar="DIVISION")
    vectorise.add_argument("--restrict", default=None, metavar="MOTIF")
    vectorise.add_argument("--mode", default="all", metavar="MODE")
    vectorise.add_argument("outfile", metavar="OUTPUT")
    vectorise.add_argument("motifs", nargs='+', metavar="MOTIFS")
    # 
    test = commands.add_parser('test')
    test.add_argument("site")
    test.add_argument("sitesfile")
    return parser.parse_args()

def _get_args(args, foldertype, filetype, extension=""):
    if os.path.isdir(args[foldertype]):
        folder = args.pop(foldertype)
        for filename in os.listdir(folder):
            if filename.endswith(extension):
                filepath = os.path.join(folder, filename)
                args[filetype] = filepath
                yield args
    else:
        args[filetype] = args.pop(foldertype)
        yield args


COMMAND = {"divide": divide, "formula": formula, "vectorise": vectorise,
           "combine": combine, "counts": counts}

def main():
    args = vars(parse_args())
    command = args.pop("command")
    if command == "bpfind":
        for args in _get_args(args, "structure", "structurefile", ""):
            bpfind(**args)
    elif command == "sites":
        args["ignore_H"] = not args.pop("has_H")
        args["invalids"] = not args.pop("only_valid")
        args["modbases"] = not args.pop("unmodified")
        args["ligdists"] = {"O": args.pop("mgo"), "N": args.pop("mgn")}
        for args in _get_args(args, "structure", "structurefile", ""):
            sites(**args)
    elif command == "assigns":
        _get_motifsfiles(args)
        for args in _get_args(args, "sites", "sitesfile", "_sites.tsv"):
            assigns(**args)
    elif command == "merge":
        if "motifs" in args:
            _get_motifsfiles(args)
        if "ids" in args and args["ids"]:
            args["dataset"] = args["ids"]
            args["ids"] = ArgParser.read_idslist(args.pop("ids"))
        else:
            args["dataset"] = args["folder"]
        args["dataset"] = ArgParser.basename(args["dataset"], "._")
        merge(**args)
    elif command in COMMAND:
        _get_motifsfiles(args)
        if "ids" in args and args["ids"]:
            args["ids"] = ArgParser.read_idslist(args.pop("ids"))
        COMMAND[command](**args)
    elif command == "visualise":
        visualise(**args)
    elif command == "redivide":
        redivide(**args)
    elif command == "test":
        test(**args)


if __name__=="__main__":
    main()
