from __future__ import print_function

import os
import sys
import argparse
from itertools import combinations
from collections import defaultdict

from numpy.linalg import norm, det
from numpy.random import randint
from numpy import mean, array, dot, cross, sqrt, ceil, flatnonzero
from numpy import all as npall, round as npround
import matplotlib.pyplot as plt

from ..magracarna.engrid import StructureFile
from ..bin.homology import ClusterFile
from ..utilities.args import ArgParser


def diff(clusterfile1, clusterfile2):
    clusters1 = ClusterFile(clusterfile1).read_all_clusters_sites()
    clusters2 = ClusterFile(clusterfile2).read_all_clusters_sites()
    differences = dict()
    uncommon1 = set(clusters1)
    uncommon2 = set(clusters2)
    for cluster1 in sorted(clusters1):
        cluster1set = set(clusters1[cluster1])
        for cluster2 in sorted(clusters2):
            if 0 in (cluster1, cluster2) and cluster1 != cluster2:
                continue
            common =  set(clusters2[cluster2]) & cluster1set
            if len(common):
                uncommon1 -= {cluster1}
                uncommon2 -= {cluster2}
                diffin1 = sorted(cluster1set - common)
                diffin2 = sorted(set(clusters2[cluster2]) - common)
                cluster1size = float(len(cluster1set))
                cluster2size = float(len(clusters2[cluster2]))
                if (len(diffin1) or len(diffin2)) \
                  and (len(common) / cluster1size > 0.1 \
                  or len(common) / cluster2size > 0.1):
                    differences[(cluster1, cluster2)] = (diffin1, diffin2)
    for cluster1, cluster2 in sorted(differences):
        diffin1, diffin2 = differences[(cluster1, cluster2)]
        print("Clusters\t:\t% 5d\t% 5d" % (cluster1, cluster2))
        for site in diffin1:
            print("\t< %d %s" % site)
        for site in diffin2:
            print("\t> %d %s" % site)
    for character, clusters, uncommon in (("<", clusters1, uncommon1), (">", clusters2, uncommon2)):
        for cluster in sorted(uncommon):
            print("Cluster \t:\t% 5d" % cluster)
            for site in clusters[cluster]:
                print("\t%s %d %s" % (character, site[0], site[1]))


class Circumsphere(object):

    def __init__(self, coos):
        self.coos = array(coos)

    def _circumscribes_triangle(self, points):
        pointA, pointB, pointC = points = [self.coos[point] for point in points]
        vecAC = pointA - pointC
        vecBC = pointB - pointC
        crossAC_BC = cross(vecAC, vecBC)
        centre = (norm(vecAC)**2)*vecBC - (norm(vecBC)**2)*vecAC
        centre = cross(centre, crossAC_BC) / (2 * norm(crossAC_BC)**2) + pointC
        centre = npround(centre, 3)
        # radius = norm(vecBC) * norm(vecAC) * norm(pointA - pointB)
        # radius = npround(radius / (2 * norm(crossAC_BC)), 3)
        radius = max([npround(norm(point - centre), 3) for point in points])
        return centre, radius

    def _circumscribes_tetrahedron(self, points):
        points = [self.coos[point] for point in points]
        main = array([[dot(each, each)] + list(each) + [1] for each in points])
        det_a = det(main[:, (1, 2, 3, 4)])
        det_c = det(main[:, (0, 1, 2, 3)])
        det_Dx = det(main[:, (0, 2, 3, 4)])
        det_Dy = -det(main[:, (0, 1, 3, 4)])
        det_Dz = det(main[:, (0, 1, 2, 4)])
        centre = npround(array([det_Dx, det_Dy, det_Dz]) / (2 * det_a), 3)
        # radius = sqrt(det_Dx**2 + det_Dy**2 + det_Dz**2 - (4 *det_a * det_c))
        # radius = npround(radius / (2 * abs(det_a)), 3)
        radius = max([npround(norm(point - centre), 3) for point in points])
        return centre, radius

    def _trivial_circumsphere(self, boundary):
        boundpoints = len(boundary)
        if boundpoints == 0:
            return [0, 0, 0], 0.0
        if boundpoints == 1:
            return self.coos[list(boundary)[0]], 0.0
        if boundpoints == 2:
            centre = npround(mean(self.coos[list(boundary)], axis=0), 3)
            radius = max([norm(self.coos[each] - centre) for each in boundary])
            return centre, npround(radius, 3)
        if boundpoints == 3:
            return self._circumscribes_triangle(list(boundary))
        if boundpoints == 4:
            return self._circumscribes_tetrahedron(list(boundary))

    def welzl(self, points, boundary=set()):
        if len(boundary) == 4 or len(points) == 0:
            if len(boundary) > 4:
                print(boundary, points)
                exit()
            return self._trivial_circumsphere(boundary)
        without = points[randint(0, len(points))]
        points_without = sorted(set(points) - {without})
        centre, radius = self.welzl(points_without, boundary)
        distance_to_without = npround(norm(self.coos[without] - centre), 3)
        if distance_to_without <= radius:
            return centre, radius
        elif distance_to_without > radius:
            return self.welzl(points_without, boundary | {without})

    @classmethod
    def circumsphere(cls, coos):
        return cls(coos).welzl(range(len(coos)))

    def verify(self, centre, radius):
        distances = npround(norm(self.coos - centre, axis=1), 3)
        outside = len(flatnonzero(distances > radius))
        if outside:
            print("Outside range: % 5d out of % 5d" % (outside, len(self.coos)))


class DatasetClusters(object):

    def __init__(self, folder, suffix, resolutions=None):
        self.folder = folder
        self.suffix = suffix
        if resolutions:
            self.resolution = self._read_resolutions(resolutions)
        else:
            self.resolution = dict()
        self._min = dict()

    def _read_resolutions(self, resolutionsfile):
        resolution = dict()
        with open(resolutionsfile, 'r') as infile:
            for line in infile:
                if line:
                    pdbid, value = line.split()
                    resolution[pdbid] = float(value)
        return resolution

    def min(self, chains, percentage):
        if chains not in self._min:
            self._min[chains] = dict()
        if percentage not in self._min[chains]:
            value = max(ceil((percentage / 100.) * chains), 2)
            self._min[chains][percentage] = value
        return self._min[chains][percentage]

    def filepaths(self):
        filepaths = list()
        for subfolder in sorted(os.listdir(self.folder)):
            if subfolder == "kTHM":
                continue
            print(subfolder, file=sys.stderr)
            clusterfile, profile, chains = None, None, 0
            subfolder = os.path.join(self.folder, subfolder)
            for filename in os.listdir(subfolder):
                if filename.endswith("clusters.txt"+self.suffix):
                    clusterfile = os.path.join(subfolder, filename)
                if filename.endswith(".pdb"):
                    profile = os.path.join(subfolder, filename)
                if filename.endswith(".msa.fasta"):
                    with open(os.path.join(subfolder, filename)) as msafile:
                        for line in msafile:
                            chains += int(line.startswith(">"))
            if clusterfile and profile and chains:
                yield profile, clusterfile, chains

    def cluster_atoms(self, profile, clusterof, siteidof):
        cluster_atoms = defaultdict(list)
        with StructureFile(profile) as manager:
            structure = manager.extract_structure(breaks=False, multimodel=True)
            _, endex, _ = structure.chains["A"]
            for each in structure.residues[endex:]:
                model, chain, resno = each.model, int(each.chain), each.resno
                sitex = (model - 1) * 100000 + int(chain) * 10000 + resno
                pdbid = siteidof[sitex].split(" ")[0]
                resolution = self.resolution.get(pdbid, None)
                atom = each.atoms.values()[0]
                entry = (atom.coos, atom.bfac, resolution, siteidof[sitex])
                cluster_atoms[clusterof[sitex]].append(entry)
        return cluster_atoms

    @staticmethod
    def mgrna_valid(siteid):
        residueid = siteid.split(" ")[1]
        if residueid.startswith("[MG]"):
            if residueid.endswith(("*", "!")):
                return False
            return True
        elif residueid.startswith("[NA]"):
            return False
        return True

    def cluster_type(self, cluster_sites, chains):
        sites = len(cluster_sites)
        if sites >= self.min(chains, 80):
            return "H"
        if sites >= self.min(chains, 50):
            return "M"
        if sites >= self.min(chains, 10):
            return "R"
        return "N" # sites < min10

    def profilewise(self, category=False, ):
        for profile, clusterfile, chains in self.filepaths():
            siteidof = dict()
            clusterof = dict()
            cluster_type = dict()
            clusters = ClusterFile(clusterfile).read_all_clusters_sites()
            for cluster in clusters:
                cluster_sites = clusters[cluster]
                for sitex, siteid in cluster_sites:
                    siteidof[sitex] = siteid
                    clusterof[sitex] = cluster
                cluster_type[cluster] = self.cluster_type(cluster_sites, chains)
            cluster_atoms = self.cluster_atoms(profile, clusterof, siteidof)
            yield cluster_type, cluster_atoms


def radius(folder, suffix):
    circumradii = list()
    ccircumradii = list()
    maxtocentroid = list()
    dataset = DatasetClusters(folder, suffix)
    for cluster_type, cluster_atoms in dataset.profilewise():
        for cluster in cluster_atoms:
            if cluster == 0:
                continue
            coos = zip(*cluster_atoms[cluster])[0]
            centre, radius = Circumsphere.circumsphere(coos)
            circumradii.append(radius)
            if cluster_type[cluster] in ("M", "H"):
                ccircumradii.append(radius)
            centroid = mean(coos, axis=0)
            tocentroid = npround(norm(coos - centroid, axis=1), 3)
            maxtocentroid.append(max(tocentroid))
    bincircradii = defaultdict(lambda: 0)
    bincentradii = defaultdict(lambda: 0)
    cbincircradii = defaultdict(lambda: 0)
    for circumradius, farthestpoint in zip(circumradii, maxtocentroid):
        bincircradii[npround(circumradius, 1)] += 1
        bincentradii[npround(farthestpoint, 1)] += 1
    for ccircumradius in ccircumradii:
        cbincircradii[npround(ccircumradius, 1)] += 1
    for binvalue in sorted(set(bincircradii) | set(bincentradii)):
        output = [binvalue]
        circradius = bincircradii[binvalue] ; output.append(circradius)
        centradius = bincentradii[binvalue] ; output.append(centradius)
        cradius = cbincircradii[binvalue] ; output.append(cradius)
        print("% 7.3f\t% 7d\t% 7d\t% 7d" % tuple(output))


def rmsd_by_resolution(folder, suffix, resolutionsfile):
    dataset = DatasetClusters(folder, suffix, resolutionsfile)
    binned = defaultdict(lambda: defaultdict(lambda: 0))
    rescolumns = set()
    for _, cluster_atoms in dataset.profilewise():
        for cluster in cluster_atoms:
            for _, rmsd, resolution, _ in cluster_atoms[cluster]:
                binned[int(rmsd * 10) / 10.][int(resolution * 10) / 10.] += 1
                rescolumns.add(int(resolution * 10) / 10.)
    rescolumns = range(int(min(rescolumns)*10), int(max(rescolumns)*10)+1)
    rescolumns = [value / 10. for value in rescolumns]
    #
    string = [" "*7]
    for value in rescolumns:
        string.append("% 7s" % value)
    print("\t".join(string))
    for rmsd in sorted(binned):
        string = ["% 7s" % rmsd]
        for value in rescolumns:
            string.append("% 7s" % binned[rmsd][value])
        print("\t".join(string))

def rmsd(folder, suffix):
    rmsds = list()
    dataset = DatasetClusters(folder, suffix)
    for _, cluster_atoms in dataset.profilewise():
        for cluster in cluster_atoms:
            rmsds += list(zip(*cluster_atoms[cluster])[1])
    for rmsdvalue in rmsds:
        print(rmsdvalue)


def conserved(folder, suffix):
    conserved = 0
    bytype = defaultdict(lambda: 0)
    dataset = DatasetClusters(folder, suffix)
    for cluster_type, cluster_atoms in dataset.profilewise():
        for cluster in sorted(set(cluster_atoms) - {0}):
            bytype[cluster_type[cluster]] += 1
            if cluster_type[cluster] in ("M", "H"):
                conserved += len(cluster_atoms[cluster])
    print(conserved, bytype)


def within(folder, suffix):
    PERCENTS = (80, 90, 95, 98, 99)
    percentwise = {pct_: defaultdict(lambda: 0.0) for pct_ in PERCENTS}
    dataset = DatasetClusters(folder, suffix)
    for cluster_type, cluster_atoms in dataset.profilewise():
        for cluster in cluster_atoms:
            if cluster == 0:
                continue
            coos = zip(*cluster_atoms[cluster])[0]
            tocentroid = npround(norm(coos - mean(coos, axis=0), axis=1), 3)
            for pct_ in PERCENTS:
                radius = sorted(tocentroid)[int(ceil(len(coos) * pct_/100.)-1)]
                percentwise[pct_][ceil(radius * 10) / 10.] += 1
    for binvalue in sorted(set.union(*map(set, percentwise.values()))):
        string = ["% 7d" % percentwise[pct_][binvalue] for pct_ in PERCENTS]
        print("% 5.1f\t%s" % (binvalue, "\t".join(string)))


def fromcentroid(folder, suffix):
    distancewise = defaultdict(lambda: 0)
    cdistancewise = defaultdict(lambda: 0)
    dataset = DatasetClusters(folder, suffix)
    for cluster_type, cluster_atoms in dataset.profilewise():
        for cluster in cluster_atoms:
            if cluster == 0:
                continue
            coos = zip(*cluster_atoms[cluster])[0]
            tocentroid = npround(norm(coos - mean(coos, axis=0), axis=1), 1)
            for distance in tocentroid:
                distancewise[distance] += 1
                if cluster_type[cluster] in ("M", "H"):
                    cdistancewise[distance] += 1
    for howfar, count in sorted(distancewise.items()):
        print("% 5.1f\t% 7d\t% 7d" % (howfar, count, cdistancewise[howfar]))


def badclusters(folder, suffix):
    dataset = DatasetClusters(folder, suffix)
    for cluster_type, cluster_atoms in dataset.profilewise():
        for cluster in cluster_atoms:
            if cluster == 0:
                continue
            coos, _, _, siteids = zip(*cluster_atoms[cluster])
            pdbids = [each.split(" ")[0] for each in siteids]
            chains = [each.split(":")[1][:-1] for each in siteids]
            chains = array([":".join(each) for each in zip(pdbids, chains)])
            if len(set(chains)) == len(chains):
                continue
            print("\tCluster : % 5d" % cluster)
            for chainid in sorted(set(chains)):
                badsites = flatnonzero(chains == chainid)
                badids = [siteids[each].strip("*! ")[5:] for each in badsites]
                if len(badsites) > 1:
                    distances = list()
                    for site1, site2 in combinations(badsites, 2):
                        distance = npround(norm(coos[site1] - coos[site2]), 3)
                        distances.append("% 7.3f" % distance)
                    print("\t\t% 7s\t%s" % (chainid, "\t".join(distances)))
                    for badid in badids:
                        print("\t\t\t%s" % badid)


def analyze_badclusters(badclustersfile, structuresdir):
    clusters = defaultdict(dict)
    with open(badclustersfile, 'r') as infile:
        for line in infile:
            if not line.strip():
                continue
            if not line.startswith("\t"):
                molecule = line.strip()
            if line.startswith("\tCluster"):
                cluster = line.split()[2]
            if line.startswith("\t\t") and not line.startswith("\t\t\t"):
                chainid = line.split()[0]
                thelist = clusters[chainid][(molecule, cluster)] = list()
            if line.startswith("\t\t\t"):
                thelist.append(line.strip())
    bypdbid = defaultdict(dict)
    for chainid in clusters.keys():
        pdbid = chainid.split(":")[0]
        for (molecule, cluster), metallist in clusters[chainid].items():
            bypdbid[pdbid][(molecule, cluster)] = tuple(metallist)
    total = len(bypdbid)
    bycluster = defaultdict(list)
    for count, pdbid in enumerate(sorted(bypdbid), 1):
        print("% 5d / % 5d\t%s" % (count, total, pdbid), file=sys.stderr)
        path = ArgParser.get_sources(structuresdir, [pdbid], ".pdb", ".cif")[0]
        residues = set(sum(bypdbid[pdbid].values(), tuple()))
        byresidue = dict()
        with StructureFile(path) as structure:
            for residue in structure.extract_residues():
                residueid = str(residue).split("/")[0]
                if residueid in residues:
                    for _, atom in residue.atoms.items():
                        if atom.atype not in ("N", "O"):
                            break
                    byresidue[residueid] = atom.coos
        for clusterkey in bypdbid[pdbid]:
            distances = list()
            cooslist = [byresidue[each] for each in bypdbid[pdbid][clusterkey]]
            for coos1, coos2 in combinations(cooslist, 2):
                distances.append(npround(norm(coos1 - coos2), 3))
            bycluster[clusterkey].append(max(distances))
    distances = [max(bycluster[each]) for each in bycluster]
    print(array(distances))


def parse_args():
    parser = argparse.ArgumentParser()
    #
    commands = parser.add_subparsers(dest="command")
    diff = commands.add_parser('diff')
    diff.add_argument("clusterfile1", metavar="CLUSTERFILE-1")
    diff.add_argument("clusterfile2", metavar="CLUSTERFILE-2")
    #
    radius = commands.add_parser('radius')
    radius.add_argument("folder")
    radius.add_argument("suffix")
    #
    resolution = commands.add_parser('resolution')
    resolution.add_argument("folder")
    resolution.add_argument("suffix")
    resolution.add_argument("resolutionsfile")
    #
    rmsd = commands.add_parser('rmsd')
    rmsd.add_argument("folder")
    rmsd.add_argument("suffix")
    #
    conserved = commands.add_parser('conserved')
    conserved.add_argument("folder")
    conserved.add_argument("suffix")
    #
    within = commands.add_parser('within')
    within.add_argument("folder")
    within.add_argument("suffix")
    #
    fromcentroid = commands.add_parser('fromcentroid')
    fromcentroid.add_argument("folder")
    fromcentroid.add_argument("suffix")
    #
    badclusters = commands.add_parser("badclusters")
    badclusters.add_argument("folder")
    badclusters.add_argument("suffix")
    #
    analyze = commands.add_parser("analyze")
    analyze.add_argument("badclustersfile")
    analyze.add_argument("structuresdir")
    return parser.parse_args()

COMMAND = {"diff": diff, "radius": radius, "resolution": rmsd_by_resolution,
            "rmsd": rmsd, "conserved": conserved, "within": within,
            "fromcentroid": fromcentroid, "badclusters": badclusters,
            "analyze": analyze_badclusters}

if __name__ == "__main__":
    args = vars(parse_args())
    COMMAND[args.pop("command")](**args)
