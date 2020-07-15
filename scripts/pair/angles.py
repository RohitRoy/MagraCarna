from __future__ import print_function

import os
import sys

from magracarna.engrid import Structure, SecondaryStructureManager as SSManager
from magracarna.heteroatoms import cations as metals

folder = "/home/rohit/B/Data/BPFind/III/"

files = list()

for file in sorted(os.listdir(folder)):
    if file.endswith(".pdb"):
        files.append(os.path.join(folder, file))

required = {("G", "C", "W", "W", "T"), ("G", "G", "W", "H", "C")}

aminoacids = ["ARG", "LYS", "HIS", "CYS", "SER", "PHE", "TYR", "GLN", "ASN"]

cations = metals + aminoacids

TEMPLATE = "\t\t%3s\t%3s\t%.3f\t%.3f\t%.1f"

for file in files:
    print(file.split("/")[-1].split(".")[0].lower())
    try:
        with SSManager(file, None) as manager:
            if None in [manager.corfile, manager.pairs]:
                manager.read_bpfind_output_files()
            structresidues = manager.structfile.extract_residues()
            corresidues = manager.corfile.extract_residues() # why not limits?
            structid = manager.structfile.structid
            grid = manager._grid_merged_residues(corresidues, structresidues, 7.0)
            structure = Structure(structid, structresidues, grid, manager.pairs, False)
    except AttributeError:
        print("Error")
        continue
    get_nearby_atoms = structure.grid.get_nearby_atoms
    for residue in structure.pairs:
        if residue.name == "G":
            itspairs = structure.pairs[residue].items()
            for index0, (pairresidue, bptype) in enumerate(itspairs):
                bptype = (residue.name, pairresidue.name)+bptype
                if bptype in required:
                	print("\t%s\t%s" % (str(residue), " ".join(bptype)))
                	n3atm, n7atm = residue.atoms["N7"], residue.atoms["N3"]
                	for hetres, hetatm, n3_d in get_nearby_atoms(10.0, n3atm):
                		hetname = hetres.name
                		if hetres.name in cations:
                			toprint = [hetatm.name, hetres.name, n3_d]
                			toprint.append(n7atm.distance_to(hetatm))
                			toprint.append(n3atm.angle_between(n7atm, hetatm))
                			print(TEMPLATE % tuple(toprint))
                	break
    del structure


# III {('G', 'G', 'W', 'H', 'C'): [123, 0], ('G', 'C', 'W', 'W', 'T'): [116, 12]}
# II  {('G', 'G', 'W', 'H', 'C'): [219, 68], ('G', 'C', 'W', 'W', 'T'): [86, 9]}