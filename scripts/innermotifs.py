from __future__ import print_function

import os
import sys
import argparse

from math import ceil
from decimal import Decimal
from numpy.linalg import norm
from collections import defaultdict
from numpy import array, flatnonzero
from itertools import count as itercount, product

from ..bin.homology import LigandsFile, MoleculeDataset, ClusterFile
from ..magracarna.engrid import Residue, HydrogenBondFinder
from ..magracarna.metalsite import RNAMgSiteEngrapher
from ..magracarna.nucleotides import allnucs as NUCLEOTIDES, purines as PURINES

DONORS = HydrogenBondFinder.DONORS
OPH = ("OP1", "OP2", "OP3")
UOB = ("O4", "O2")


class InnerMotifsAnalyzer(object):

    def __init__(self, structure, alignment, max_da=Decimal("3.50")):
        hbfinder = HydrogenBondFinder(max_da=max_da)
        self.detector = RNAMgSiteEngrapher(hbfinder, relaxed=True)
        self.detector.set_structure(structure)
        self.structure = structure
        self.alignment = alignment

    MOTIFS = {"Purine-N7 Seat": (("I", "N7", "R", 2),),
      "G-N7 MacroChelate I": (("I", "N7", "G", 1), ("O", "Oph", "G", 1)),
      "10-Member Ring Purine-N7": (("I", "Oph", "*", 2), ("I", "N7", "R", 1)),
      "G-Phosphate": (("I", "Ob", "G", 1), ("I", "Oph", "*", 1)),
      "U-Phosphate": (("I", "Ob", "U", 1), ("I", "Oph", "*", 1))}

    DICTIONARY = {"*": NUCLEOTIDES, "R": PURINES,
                "Oph": ("OP1", "OP2", "OP3"), "Ob": ("O2", "O4", "O6")}

    BINS = {"O": (0.0, 1.9, 2.3, 2.6, 3.2, 3.8, 4.6, 5.0, float("inf")),
            "N": (0.0, 1.9, 2.4, 2.6, 3.2, 3.8, 4.6, 5.0, float("inf"))}
    for element, breaks in BINS.items():
        BINS[element] = list(zip(breaks, breaks[1:], itercount(1)))

    @classmethod
    def check_definition(cls, motif, inner, outer):
        definition = cls.MOTIFS[motif]
        expected_CN = sum([part[3] for part in definition if part[0] == "I"])
        if len(inner) != expected_CN:
            return None
        if "O" in zip(*definition)[0]:
            for residue, atom in outer:
                if residue.name not in NUCLEOTIDES:
                    return None
        #
        matches = defaultdict(list)
        for index, (shell, atomtype, resname, count) in enumerate(definition):
            for residue, atom in (inner if shell == "I" else outer):
                if atom in cls.DICTIONARY.get(atomtype, (atomtype,)):
                    if residue.name in cls.DICTIONARY.get(resname, (resname,)):
                        matches[index].append(residue)
                    else:
                        return None
            if len(matches[index]) != count:
                return None
            matches[index].sort(key=lambda each: (each.resno, each.insco))
        return tuple(sum(zip(*sorted(matches.items()))[1], list()))

    def check_motif(self, motif, inner, outer):
        match = self.check_definition(motif, inner, outer)
        if motif == "G-N7 MacroChelate I":
            if match[0] != match[1]:
                return None
        elif motif == "10-Member Ring Purine-N7":
            for isat5end, isat3end in zip(match, match[1:]):
                if not self.structure.linked(isat5end, isat3end):
                    return None
        return match

    def get_shells(self, residue):
        self.detector.detect_shells(residue)
        shells = list()
        for detected in (self.detector.inner_shell, self.detector.outer_shell):
            shells.append(list())
            for entry in detected:
                if entry[0].name != "HOH":
                    shells[-1].append((entry[0], entry[1].name))
        return tuple(shells)

    def assign_motifs(self, residue):
        inner_shell, outer_shell = self.get_shells(residue)
        for motif in sorted(self.MOTIFS):
            try:
                match = self.check_motif(motif, inner_shell, outer_shell)
                loci = [self.alignment.index(each) + 1 for each in match]
                return motif, tuple(loci)
            except TypeError: # for when match == None, and can't iterate
                continue
        return None

    def _ligand_contact(self, motif, residues, metal):
        variant = list()
        inner, outer = self.get_shells(metal)
        ligands = [("I", residue, atom) for residue, atom in inner]
        ligands += [("O", residue, atom) for residue, atom in outer]
        for shell, atomtype, _, count in self.MOTIFS[motif]:
            for _ in range(count):
                residue, residues = residues[0], residues[1:]
                contact = list()
                res_ligands = [each for each in ligands if each[1] == residue]
                for ligand in res_ligands:
                    shell_, residue_, atomname_ = ligand
                    if atomname_ in self.DICTIONARY.get(atomtype, (atomtype,)):
                        contacttype_ = "1" if shell_ == shell else "2"
                        contact.append((contacttype_, ligand))
                variant.append(min(contact)[0] if contact else "4")
                if contact:
                    ligands.remove(min(contact)[1])
        return variant, ligands

    def _additional_ligands(self, motif, extra_ligands):
        if any([each for each in extra_ligands if each[0] == "I"]):
            return True
        definition = self.MOTIFS[motif]
        if "O" in zip(*definition)[0]:
            outer = [each for each in extra_ligands if each[0] == "O"]
            if not outer:
                return False
            atoms = set(zip(*outer)[2])
            types = [each[1] for each in definition if each[0] == "O"]
            types = [list(self.DICTIONARY.get(each, (each))) for each in types]
            if atoms & set(sum(types, [])):
                return True
            resnames = [each.name for each in zip(*outer)[1]]
            if set(resnames) - set(self.DICTIONARY["*"]):
                return True
        return False

    def interaction_variants(self, motif, loci, siteids):
        """
            1 : interaction as expected
            2 : interaction in different shell
            3 : additional interactions
            4 : missing interactions
            5 : different metal ion
        """
        variations = list()
        residues = [self.alignment[locus - 1] for locus in loci]
        if all(residues):
            metalids = [siteid.split()[1].strip("!* ") for siteid in siteids]
            for metal in self.structure.residues[self.structure.chainend:]:
                if str(metal) not in metalids:
                    continue
                if metal.name != "MG":
                    variations.append("5")
                else:
                    types, extra = self._ligand_contact(motif, residues, metal)
                    if self._additional_ligands(motif, extra):
                        types.append("3")
                    variations.append(max(types))
        return variations

    def metal_to_ligands(self, motif, residues, metal):
        distances = list()
        definition = self.MOTIFS[motif]
        for shell, atomtype, resname, count in definition:
            for _ in range(count):
                toligand = list()
                residue, residues = residues[0], residues[1:]
                for atomname, atom in residue.atoms.items():
                    if atomname in self.DICTIONARY.get(atomtype, (atomtype,)):
                        toligand.append(round(norm(metal.coos - atom.coos), 2))
                toligand = min(toligand) if toligand else float("inf")
                distances.append(toligand)
        return distances

    def distances_to_ligands(self, motif, loci, siteids):
        distances = list()
        residues = [self.alignment[locus - 1] for locus in loci]
        if all(residues):
            metalids = [siteid.split()[1].strip("!* ") for siteid in siteids]
            for residue in self.structure.residues[self.structure.chainend:]:
                if str(residue) in metalids:
                    for metal in residue.atoms.values():
                        if metal.atype not in ("O", "N"):
                            break
                    vector = self.metal_to_ligands(motif, residues, metal)
                    if float("inf") not in vector:
                        distances.append(vector)
        return distances


class MoleculeSitesDataset(MoleculeDataset):

    def __init__(self, msafile, ligands, structuredir, ids=None):
        kwargs = {"ids": ids, "msafile": msafile}
        super(MoleculeSitesDataset, self).__init__("*", "*", **kwargs)
        self.ligands = ligands
        self.structures = structuredir

    def iterate(self, valids=False, onlymg=False):
        message = "%% 4d / % 4d    \t%%s:%%s" % len(self.chains)
        args = {"count": True, "structuredir": self.structures, "align": True}
        dataiter = super(MoleculeSitesDataset, self).iterate(**args)
        for entry in dataiter:
            pdbid, chainid, count, structure, aligned = entry
            sites = set(self.ligands.get((pdbid, chainid), []))
            sites.update(self.off_chain_sites(structure, aligned))
            sites = self.only_mg(sites) if onlymg else sites
            sites = self.mgrna_valids(sites) if valids else sites
            #
            print(message % (count, pdbid, chainid), file=sys.stderr)
            analyzer = InnerMotifsAnalyzer(structure, aligned, Decimal("3.50"))
            yield chainid, structure, analyzer, sites

    def off_chain_sites(self, structure, aligned):
        offchain = set()
        pdbid = structure.structid
        extra = set(self.ligands.keys()) - set(self.chains)
        extra = [chain_ for pdbid_, chain_ in extra if pdbid_ == pdbid]
        within = structure.grid.get_nearby_atoms
        for chainid in extra:
            for siteid, ligandslist in self.ligands[(pdbid, chainid)].items():
                try:
                    locus, atom = ligandslist[0]
                    ligand = aligned[locus-1].atoms[atom]
                    residueid = siteid.strip("!* ").split(" ")[1]
                    for residue, _, _ in within(7.0, ligand):
                        if str(residue) == residueid:
                            offchain.add(siteid)
                            break
                except:
                    continue
        return offchain

    @staticmethod
    def mgrna_valids(siteids):
        valids = set()
        for siteid in siteids:
            if "[MG]" not in siteid and "[NA]" not in siteid:
                valids.add(siteid)
            elif "[MG]" in siteid:
                if not siteid.endswith(("!", "*")):
                    valids.add(siteid)
        return valids

    @staticmethod
    def only_mg(siteids):
        return {siteid for siteid in siteids if "[MG]" in siteid}


class Atlas(object):

    def __init__(self, folder, structuredir):
        self.folder = folder
        self.structures = structuredir

    def filepaths(self):
        filepaths = list()
        for subfolder in sorted(os.listdir(self.folder)):
            if subfolder not in ("05TT", "16EC", "16TT", "23EC", "23HM", "23TT", "D2GrS", "FMNrS3"):
                continue
            clusterfile, ligandsfile, msafile = None, None, None
            subfolder = os.path.join(self.folder, subfolder)
            for filename in os.listdir(subfolder):
                if filename.endswith("clusters.txt"):
                    clusterfile = os.path.join(subfolder, filename)
                if filename.endswith(".ligands.txt"):
                    ligandsfile = os.path.join(subfolder, filename)
                if filename.endswith(".msa.fasta"):
                    msafile = os.path.join(subfolder, filename)
            if clusterfile and ligandsfile and msafile:
                molecule = os.path.basename(subfolder)
                yield molecule, clusterfile, ligandsfile, msafile

    def profilewise(self):
        for molecule, clusterfile, ligandsfile, msafile in self.filepaths():
            ligands = self.read_ligands(ligandsfile)
            dataset = MoleculeSitesDataset(msafile, ligands, self.structures)
            chains = set(ligands.keys())
            extra = chains - set(dataset.chains)
            for chain in chains:
                if chain not in extra:
                    ligands[chain] = list(ligands[chain].keys())
            clusterof, sitesof = self.read_clusters(clusterfile)
            yield molecule, dataset, clusterof, sitesof

    @staticmethod
    def read_ligands(ligandsfile):
        ligands = defaultdict(dict)
        for _, siteid, ligandslist in LigandsFile.read(ligandsfile):
            siteid = siteid.strip()
            pdbid, residueid = siteid.split(" ")
            chainid = residueid.split(":")[1].strip("!*")
            ligandslist = [ligand.strip().split(".") for ligand in ligandslist]
            ligandslist = [(int(locus), atom) for locus, atom in ligandslist]
            ligands[(pdbid, chainid)][siteid] = ligandslist
        return dict(ligands)

    @staticmethod
    def read_clusters(clusterfile):
        clusterfile = ClusterFile(clusterfile)
        clusternames = clusterfile.read_cluster_names()
        clusterof = dict()
        sitesof = defaultdict(list)
        for cluster in clusternames:
            for _, siteid in clusterfile.read_cluster_sites(cluster):
                clusterof[siteid.strip()] = cluster
                sitesof[cluster].append(siteid.strip())
        return clusterof, dict(sitesof)


class ConservedMotifsFile(object):

   def __init__(self, filepath):
      self.path = filepath

   def read(self, molecule):
       byMotif = defaultdict(list)
       motif, ismolecule = None, False
       with open(self.path, 'r') as infile:
           for line in infile:
               if not line.strip():
                   continue
               elif not line.startswith("\t"):
                   ismolecule = (molecule == line.strip())
               elif "(" in line and ismolecule:
                   cluster, loci = line.strip().split("\t")
                   loci = [each for each in loci.strip("()").split(",") if each]
                   loci = tuple(map(int, loci))
                   byMotif[motif].append((int(cluster), loci))
               elif line.startswith("\t") and ismolecule:
                   motif = line.strip()
       return dict(byMotif)

   def new(self):
      with open(self.path, "w") as motifsfile:
         motifsfile.write("")

   def write(self, molecule, conservedmotifs):
      with open(self.path, "a") as cmfile:
          print(molecule, file=cmfile)
          for motif in sorted(conservedmotifs):
              print("\t%s" % motif, file=cmfile)
              for cluster, loci in sorted(conservedmotifs[motif]):
                  print("\t\t%s\t%s" % (cluster, loci), file=cmfile)


def count_reliable(dataset, structuredir, clusterof):
    motifs = defaultdict(lambda: defaultdict(lambda: 0))
    dataiter = dataset.iterate(valids=True, onlymg=True)
    for chain, structure, analyzer, inchain in dataiter:
        reliable = dict()
        for siteid in inchain:
            reliable[siteid.strip("!* ").split(" ")[1]] = clusterof[siteid]
        for residue in structure.residues[structure.chainend:]:
            if str(residue) in reliable:
                try:
                    motif, loci = analyzer.assign_motifs(residue)
                    cluster = reliable[str(residue)]
                    if cluster:
                        motifs[motif][(cluster, loci)] += 1
                except TypeError: # for when match == None, as it can't iterate
                    continue
    return motifs

def conserved_clusters(motifcounts, sitesof, dataset):
    chains = len(dataset.chains)
    conserved = defaultdict(set)
    for motif in motifcounts:
        for (cluster, loci), count in motifcounts[motif].items():
            if cluster != 0 and count >= 2:
                if len(sitesof[cluster]) >= max(int(ceil(chains * 0.50)), 2):
                    conserved[motif].add((cluster, loci))
    return conserved

def conserved(profilesdir, structuredir, conservedfile, allmotifsfile):
    atlas = Atlas(profilesdir, structuredir)
    allmotifsfile = ConservedMotifsFile(allmotifsfile)
    conservedfile = ConservedMotifsFile(conservedfile)
    allmotifsfile.new()
    conservedfile.new()
    for entry in atlas.profilewise():
        molecule, dataset, clusterof, sitesof = entry
        instances = count_reliable(dataset, structuredir, clusterof)
        conserved = conserved_clusters(instances, sitesof, dataset)
        allmotifsfile.write(molecule, instances)
        conservedfile.write(molecule, conserved)


def print_binned_distances(bins, filepath):
    with open(filepath, "w") as outfile:
        for (motif, locusex), bybin in bins.items():
            print("\n% -20s\t%d" % (motif, locusex), file=outfile)
            contexts = sorted(bybin["max"])
            header = ["% -6s % 4d" % context for context in contexts]
            print("\t".join([" "*3] + header), file=outfile)
            for name, bycontext in sorted(bybin.items()):
                row_ = ["% 11.2f" % bycontext[context] for context in contexts]
                print("\t".join(["% 3s" % name] + row_), file=outfile)

def print_binned_combinations(combos, filepath):
    with open(filepath, "w") as outfile:
        for motif, bycombo in combos.items():
            print("\n% -20s" % motif, file=outfile)
            contexts = sorted(bycombo[(1,) * len(bycombo.keys()[0])])
            header = ["% -6s % 4d" % context for context in contexts]
            print("\t".join([" "*3] + header), file=outfile)
            for combo, bycontext in sorted(bycombo.items()):
                combo = "".join(map(str, combo))
                row_ = [bycontext[context] for context in contexts]
                if sum(row_) > 0:
                    row_ = [combo] + ["% 11.2f" % part for part in row_]
                    print("\t".join(row_), file=outfile)

def bin_distances(distances):
    BINS = InnerMotifsAnalyzer.BINS
    bins = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 0)))
    combos = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 0)))
    for molecule, bymolecule in distances.items():
        for motif, bymotif in sorted(bymolecule.items()):
            elements = dict()
            for _, atomtype, _, recurs in InnerMotifsAnalyzer.MOTIFS[motif]:
                for index in range(recurs):
                    elements[len(elements)] = "O" if "O" in atomtype else "N"
            for (cluster, loci), vectors in bymotif.items():
                context = (molecule, cluster)
                binsets = defaultdict(set)
                for locusex, darray in enumerate(map(array, zip(*vectors))):
                    locus = (motif, locusex)
                    for min_, max_, name in BINS[elements[locusex]]:
                        inbin = flatnonzero((min_ <= darray) * (darray < max_))
                        binsets[(locusex, name)] = set(inbin)
                        bins[locus][name][context] = len(inbin)
                    bins[locus]["max"][context] = max(darray)
                for name in product(*[range(1, 9)] * len(loci)): # name = combo
                    sets = [binsets[each] for each in enumerate(name)]
                    combos[motif][name][context] = len(set.intersection(*sets))
    return bins, combos

def print_binned_interactions(variants, filepath):
    with open(filepath, "w") as outfile:
        for motif, byvariant in variants.items():
            print("\n% -20s" % motif, file=outfile)
            contexts = set()
            for varianttype, bycontext in byvariant.items():
                contexts.update(bycontext.keys())
            contexts = sorted(contexts)
            header = ["% -6s % 4d" % context for context in contexts]
            print("\t".join([" "*3] + header), file=outfile)
            for varianttype, bycontext in sorted(byvariant.items()):
                row_ = ["% 11.2f" % bycontext[context] for context in contexts]
                print("\t".join(["% 3s" % varianttype] + row_), file=outfile)

def variations(profilesdir, structuredir, motifsfile, distance_tables=None,
                distance_combos=None, interaction_tables=None):
    distances = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    variants = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 0)))
    atlas = Atlas(profilesdir, structuredir)
    conserved_motif_sites = ConservedMotifsFile(motifsfile)
    for entry in atlas.profilewise():
        molecule, dataset, clusterof, sitesof = entry
        conserved = conserved_motif_sites.read(molecule)
        if not conserved:
            continue
        for chain, structure, analyzer, inchain in dataset.iterate():
            for motif in conserved:
                for cluster, loci in sorted(conserved[motif]):
                    siteids = sorted(set(sitesof[cluster]) & set(inchain))
                    toadd = analyzer.distances_to_ligands(motif, loci, siteids)
                    distances[molecule][motif][(cluster, loci)] += toadd
                    tobin = analyzer.interaction_variants(motif, loci, siteids)
                    for varianttype in tobin:
                        variants[motif][varianttype][(molecule, cluster)] += 1
    bins, combos = bin_distances(distances)
    if distance_tables:
        print_binned_distances(bins, distance_tables)
    if distance_combos:
        print_binned_combinations(combos, distance_combos)
    if interaction_tables:
        print_binned_interactions(variants, interaction_tables)


def parse_args():
    parser = argparse.ArgumentParser()
    commands = parser.add_subparsers(dest="command")
    #
    conserved = commands.add_parser('conserved')
    conserved.add_argument("profilesdir", metavar="PROFILES-FOLDER")
    conserved.add_argument("structuredir", metavar="STRUCTURE-FOLDER")
    conserved.add_argument("conservedfile", metavar="CONSERVED-MOTIFS-FILE")
    conserved.add_argument("allmotifsfile", metavar="ALL-MOTIFS-FILE")
    #
    variations = commands.add_parser('variations')
    variations.add_argument("profilesdir", metavar="PROFILES-FOLDER")
    variations.add_argument("structuredir", metavar="STRUCTURE-FOLDER")
    variations.add_argument("motifsfile", metavar="CONSERVED-MOTIFS-SITES")
    variations.add_argument("--distance-tables", "-d", default=None)
    variations.add_argument("--distance-combos", "-c", default=None)
    variations.add_argument("--interaction-tables", "-i", default=None)
    #
    return parser.parse_args()


COMMANDS = {"conserved": conserved, "variations": variations}


if __name__ == "__main__":
    args = vars(parse_args())
    COMMANDS[args.pop("command")](**args)
