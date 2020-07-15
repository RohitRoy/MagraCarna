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
from ..magracarna.consensus import ConsensusFinder, Discriminator

_TYPES = ["context", "assigns", "counts", "chains"]

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
    # context-graph parameters
    engrapher = RNAMgSiteEngrapher(hbfinder, plinked=True, depth=depth,
                                   hbonds=hbonds, ligdists=ligdists,
                                   relaxed=relaxed, invalids=invalids,
                                   modbases=modbases, verbose=verbose)
    return hbfinder, engrapher


def context(structurefile, rundir=None, **kwargs):
    hbfinder, engrapher = _parameters(**kwargs)
    mgsitelist = RNAMgSiteList(engrapher)
    if rundir:
        bpfind_run = SSManager.BPFindParams(rundir, False, False, hbfinder)
    else:
        bpfind_run = None
    with SSManager(structurefile, bpfind_run) as manager:
        structure = manager.extract_paired_structure()
        mgsitelist.add_sites_from(structure)
        if rundir:
            structurefile = manager.structfile.write_to
    contextfile = _filepath(structurefile, True, "context")
    mgsitelist.dump_graphs(contextfile)
    chainsfile = _filepath(structurefile, True, "chains")
    RNAMgSiteContextParser.create_chainsfile(mgsitelist.graphs, chainsfile)

#

def assigns(contextfile, motifsfiles, mode="first", reassign=False,
            detail=False, sort_by=RNAMgMotifCounts._SORTERS):
    rnagfs = RNAMgSiteList.load(contextfile)
    assigner = MotifAssigner.extract(*motifsfiles, mode=mode)
    assignsfile = _filepath(contextfile, False, "assigns")
    if reassign and os.path.exists(assignsfile):
        assignments = MotifAssignments.load(assigner, assignsfile, detail)
    else:
        assignments = MotifAssignments(assigner)
    assignments.add_from(rnagfs)
    assignments.dump(assignsfile, detail)
    countsfile = _filepath(contextfile, False, "counts")
    counts = RNAMgMotifCounts.create_from_assigns(assignments)
    counts.dump(countsfile, sort_by)

#

def counts(assignsfile, countsfile, motifsfiles):
    assigner = MotifAssigner.extract(*motifsfiles, mode="first")
    assigns = MotifAssignments.load(assigner, assignsfile, detail=False)
    RNAMgMotifCounts.create_from_assigns(assigns).dump(countsfile)

# 

def merge_context(sources, destination, only_valid=False):
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

MERGE_TYPE = {"context": merge_context,
              "assigns": merge_assigns,
              "counts": merge_counts}

def merge(filetype, folder, dataset, ids=None, **kwargs):
    extension = "_"+filetype+".tsv"
    sources = ArgParser.get_sources(folder, ids, extension)
    destination = os.path.join(folder, dataset) + extension
    MERGE_TYPE[filetype](sources, destination, **kwargs)

#

def visualise(structure, contextfile, site=None, assignsfile=None, context=False, waters=False):
    names = [site] if site else None
    rnagfs = RNAMgSiteEngrapher.decode_file(contextfile, names)
    visualiser = RNAMgSiteContextVisualiser(structure, rnagfs)
    if site:
        visualiser.display(site, context, waters)
    else:
        print("\t".join(sorted(visualiser.site)))
        visualiser.interactive()

#

def divide(name, assigns, outfile, motifsfiles,
           ids=None, show=False, contextfolder=None,
           only_valid=False, limited=None, restrict=None):
    kwargs = {"limited": limited, "restrict": restrict}
    kwargs["invalids"] = not only_valid
    assigner = MotifAssigner.extract(*motifsfiles, mode="all")
    divider = RNAMgMotifDivider(name, assigner, **kwargs)
    assigns = [assigns]
    if os.path.isdir(assigns[0]):
        assigns = ArgParser.get_sources(assigns[0], ids, "_assigns.tsv")
    divider.add_embeddings(assigns, bool(contextfolder))
    divider.summarise(outfile)
    if contextfolder:
        contexts = ArgParser.get_sources(contextfolder, ids, "_context.tsv")
        if len(assigns) > 1:
            assert(len(assigns) == len(contexts))
        divider.add_sequences(contexts)
        show = True
    divider.divide(outfile, -1, show, bool(contextfolder))

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

def consensus(base, contexts, limited, outfile, motifsfiles, ids=None,
                    division=None, restrict=None, rounds=3, covered=0):
    assigner = MotifAssigner.extract(*motifsfiles)
    basecode = assigner.graphs[base]
    cases = None
    if limited:
        cases = CaseReader.read_limited(limited, restrict, division)
        caseids = {case.split(" ")[0] for case in cases}
        ids = sorted(caseids & set(ids)) if ids else sorted(caseids)
    if os.path.isdir(contexts):
        contexts = ArgParser.get_sources(contexts, ids, "_context.tsv")
    else:
        assert(os.path.isfile(contexts))
        contexts = [contexts]
    print(len(cases))
    rnagfs = RNAMgSiteList.limited(contexts, cases, dumpfile=None)
    print(len(rnagfs))
    ConsensusFinder.MAX_DEPTH = rounds
    finder = ConsensusFinder(rnagfs, basecode)
    finder.start(outfile, rounds, covered)

def discriminate(base, contexts, limited, outfile, motifsfiles, ids=None):
	assigner = MotifAssigner.extract(*motifsfiles)
	basecode = assigner.graphs[base]
	contexts = ArgParser.get_sources(contexts, ids, "_context.tsv")
	rnagfs = RNAMgSiteList.limited(contexts, CaseReader.read_limited(limited))
	print(len(rnagfs))
	discer = Discriminator(rnagfs, basecode)
	discer.start()

# 

def test(site, contextfile):
    for baseseqs, idseqs in RNAMgSiteEngrapher.extract_sequences(contextfile, [site]):
        for baseseq, idseq in zip(baseseqs, idseqs):
            print("%s\n\t%s" % (baseseq, " ".join(idseq)))

def parse_args():
    parser = argparse.ArgumentParser()
    commands = parser.add_subparsers(dest="command")
    # 
    context = commands.add_parser('bpfind')
    context.add_argument("structure", metavar="STRUCTURE")
    context.add_argument("--rundir", "-d", default=os.path.realpath("."))
    context.add_argument("--verbose", "-v", action="store_false")
    context.add_argument("--max-da", type=Decimal, default=HBFinder.MAX_DA)
    #
    context = commands.add_parser('context')
    context.add_argument("structure", metavar="STRUCTURE")
    context.add_argument("--rundir", "-d", default=None)
    context.add_argument("--has-H", "-H", action="store_true")
    context.add_argument("--hbonds", "-B", action="store_true")
    context.add_argument("--relaxed", "-R", action="store_true")
    context.add_argument("--verbose", "-v", action="store_true")
    context.add_argument("--only-valid", "-V", action="store_true")
    context.add_argument("--unmodified", "-U", action="store_true")
    context.add_argument("--depth", "-D", default=3, type=int)
    LIGDISTS = RNAMgSiteEngrapher.LIGDISTS
    context.add_argument("--mgo", "-O", type=Decimal, default=LIGDISTS["O"])
    context.add_argument("--mgn", "-N", type=Decimal, default=LIGDISTS["N"])
    context.add_argument("--max-da", type=Decimal, default=HBFinder.MAX_DA)
    context.add_argument("--min-da", type=Decimal, default=HBFinder.MIN_DA)
    context.add_argument("--max-ha", type=Decimal, default=HBFinder.MAX_HA)
    context.add_argument("--min-ha", type=Decimal, default=HBFinder.MIN_HA)
    context.add_argument("--max-dha", type=Decimal, default=HBFinder.MAX_DHA)
    context.add_argument("--min-dha", type=Decimal, default=HBFinder.MIN_DHA)
    #
    assigns = commands.add_parser('assigns')
    assigns.add_argument("context", metavar="CONTEXT")
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
    merge_context = filetypes.add_parser("context")
    merge_context.add_argument("--only-valid", "-V", action="store_true")
    merge_assigns = filetypes.add_parser("assigns")
    merge_assigns.add_argument("motifs", nargs="+", metavar="MOTIFS")
    merge_counts = filetypes.add_parser("counts")
    merge_counts.add_argument("motifs", nargs="+", metavar="MOTIFS")
    #
    visualise = commands.add_parser('visualise')
    visualise.add_argument("structure", metavar="STRUCTURE")
    visualise.add_argument("contextfile", metavar="CONTEXT-FILE")
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
    divide.add_argument("--contextfolder", metavar="CONTEXT-FOLDER")
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
    consensus = commands.add_parser('consensus')
    consensus.add_argument("base", metavar="MOTIF-GRAPH-NAME")
    consensus.add_argument("contexts", metavar="CONTEXT")
    consensus.add_argument("--ids", default=None, metavar="IDSLIST")
    consensus.add_argument("--limited", required=True, metavar="SITES-FILE")
    consensus.add_argument("--division", default=None, metavar="DIVISION")
    consensus.add_argument("--restrict", default=None, metavar="MOTIF")
    consensus.add_argument("--rounds", default=3, type=int, metavar="ROUNDS")
    consensus.add_argument("--covered", default=0, type=int, metavar="COVERED")
    consensus.add_argument("outfile", metavar="OUTPUT")
    consensus.add_argument("motifs", nargs='+', metavar="MOTIFS")
    # 
    discriminate = commands.add_parser('discriminate')
    discriminate.add_argument("base", metavar="MOTIF-GRAPH-NAME")
    discriminate.add_argument("contexts", metavar="CONTEXT")
    discriminate.add_argument("--ids", default=None, metavar="IDSLIST")
    discriminate.add_argument("--limited", required=True, metavar="SITES-FILE")
    discriminate.add_argument("outfile", metavar="OUTPUT")
    discriminate.add_argument("motifs", nargs='+', metavar="MOTIFS")
    # 
    test = commands.add_parser('test')
    test.add_argument("site")
    test.add_argument("contextfile")
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
           "combine": combine, "counts": counts, "consensus": consensus,
           "discriminate": discriminate}

def main():
    args = vars(parse_args())
    command = args.pop("command")
    if command == "bpfind":
        for args in _get_args(args, "structure", "structurefile", ""):
            bpfind(**args)
    elif command == "context":
        args["ignore_H"] = not args.pop("has_H")
        args["invalids"] = not args.pop("only_valid")
        args["modbases"] = not args.pop("unmodified")
        args["ligdists"] = {"O": args.pop("mgo"), "N": args.pop("mgn")}
        for args in _get_args(args, "structure", "structurefile", ""):
            context(**args)
    elif command == "assigns":
        _get_motifsfiles(args)
        for args in _get_args(args, "context", "contextfile", "_context.tsv"):
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
