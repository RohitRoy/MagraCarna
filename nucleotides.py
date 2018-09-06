from itertools import product
from collections import deque, defaultdict

PO_DISTANCE = 2.00
DA_DISTANCE = (2.2, 3.8)
HA_DISTANCE = (1.2, 2.7)
DHA_ANGLE = (110, 180)

atoms = dict()
rings = dict()
heavy = dict()
edges = dict()
bonds = dict()
angles = dict()
hnames = dict()
dihedrals = dict()

purines = ['A', 'G']
pyrimidines = ['C', 'U']
nucleotides = sorted(purines+pyrimidines)
phosphate = ["OP1", "OP2", "OP3", "P"]

backbone = ["P", "O5'", "C5'", "C4'", "C3'", "O3'"]
ribose = deque(["C1'", "C2'", "C3'", "C4'", "O4'"])
imidazole = deque(['C4', 'C5', 'N7', 'C8', 'N9'])
pyrimidine = deque(['N1', 'C2', 'N3', 'C4', 'C5', 'C6'])


rings['R'] = [pyrimidine, imidazole]
rings['Y'] = [pyrimidine]
rings['A'] = rings['R']
rings['G'] = rings['R']
rings['C'] = rings['Y']
rings['U'] = rings['Y']


def heavy_atoms(alist):
    alist = alist[:]
    for atom in alist[:]:
        if 'C' in atom:
            alist.remove(atom)
    return alist


def measures(atoms, n, cyclic=True):
    ring = list()
    if cyclic:
        for i in range(n):
            ring.append(list(atoms))
            atoms.rotate(1)
        ring.reverse()
    else:
        for i in range(n):
            ring.append(atoms[i:])
    return list(zip(*tuple(ring)))

heavy['ribose'] = heavy_atoms(list(ribose))+phosphate[:3]
bonds['ribose'] = measures(ribose, 2)
bonds['ribose'] += [("C2'", "O2'"), ("C3'", "O3'"),\
                    ("C4'", "C5'"), ("C5'", "O5'")]
angles['ribose'] = measures(ribose, 3)
angles['ribose'] += [("C1'", "C2'", "O2'"), ("C3'", "C2'", "O2'"),\
                     ("C2'", "C3'", "O3'"), ("C4'", "C3'", "O3'"),\
                     ("C3'", "C4'", "C5'"), ("O4'", "C4'", "C5'"),\
                     ("C4'", "C5'", "O5'")]

dihedrals['ribose'] = measures(ribose, 4)

heavy['Y'] = heavy_atoms(list(pyrimidine))+heavy_atoms(backbone)+["O2'"]
heavy['Y'].remove("P")
heavy['R'] = heavy_atoms(list(imidazole))+heavy['Y']
heavy['Y'] += ["O2"]+heavy['ribose']
heavy['R'] += heavy['ribose']

heavy['C'] = heavy['Y']+["N4"]
heavy['U'] = heavy['Y']+["O4"]
heavy['A'] = heavy['R']+["N6"]
heavy['G'] = heavy['R']+["O6", "N2"]


bonds['Y'] = measures(pyrimidine, 2)
bonds['R'] = bonds['Y']+measures(imidazole, 2)
bonds['R'].remove(('C4', 'C5'))
bonds['Y'].append(('C2', 'O2'))
bonds['R'] += bonds['ribose']+[("C1'", "N9")]
bonds['Y'] += bonds['ribose']+[("C1'", "N1")]

angles['Y'] = measures(pyrimidine, 3)
angles['R'] = angles['Y']+measures(imidazole, 3)
angles['R'] += angles['ribose']+[('C6', 'C5', 'N7'), ('N2', 'C4', 'N9')]
angles['Y'] += angles['ribose']+[('N1', 'C2', 'O2'), ('N3', 'C2', 'O2')]

dihedrals['R'] = dihedrals['ribose'] + [("O4'", "C1'", "N9", "C4")]
dihedrals['Y'] = dihedrals['ribose'] + [("O4'", "C1'", "N1", "C2")]

bonds['C'] = bonds['Y']+[("C4", "N4")]
bonds['U'] = bonds['Y']+[("C4", "O4")]
bonds['A'] = bonds['R']+[("C6", "N6")]
bonds['G'] = bonds['R']+[("C6", "O6"), ("C2", "N2")]

angles['C'] = angles['Y']+[("N3", "C4", "N4"), ("C5", "C4", "N4")]
angles['U'] = angles['Y']+[("N3", "C4", "O4"), ("C5", "C4", "O4")]
angles['A'] = angles['R']+[("N1", "C6", "N6"), ("C5", "C6", "N6")]
angles['G'] = angles['R']+[("N1", "C6", "O6"), ("C5", "C6", "O6")]\
                         +[("N1", "C2", "N2"), ("N3", "C2", "N2")]

dihedrals['C'] = dihedrals['Y']
dihedrals['U'] = dihedrals['Y']
dihedrals['A'] = dihedrals['R']
dihedrals['G'] = dihedrals['R']

variants = dict()

variants['A'] = ["A", "ADE", "DA5", "DA3", "1MA", "MIA", "+A", "AMP", "AMO",\
                 "12A", "AET", "PSD", "AVC", "APC", "GOM", "MAD", "A23",\
                 "ATP", "2MA", "A2M", "T6A", "RIA", "6MZ", "6IA", "DA", "ADP",\
                 "5AA", "PR5", "2AD", "3DA", "ANZ", "AVC", "TSB", "QSI", "VAA",\
                 "RA", "RA5", "RA3", "SAM",\
                 "0A"]
variants['C'] = ["C", "CYT", "DC5", "DC3", "5MC", "+C", "OMC", "S4C", "CB2",\
                 "5IC", "CCC", "1SC", "DC", "CBV", "DCZ", "CSL", "CBR", "C38",\
                 "BLS", "RC", "RC5", "RC3",\
                 "0C"]
variants['G'] = ["G", "DG", "GUA", "DG5", "DG3", "2MG", "OMG", "G7M", "GNP",\
                 "7MG", "1MG", "5CG", "+G", "GDP", "M2G", "I", "GMP", "GTP",\
                 "PGP", "YG", "YYG", "2PR", "XUG", "RG", "RG5", "RG3", "6OG",\
                 "0G"]
variants['U'] = ["T", "DT", "THY", "DT5", "DT3", "U", "URA", "H2U", "5MU",\
                 "2MU", "4SU", "FMU", "CMO", "OMU", "70U", "+U", "DHU", "UR3",\
                 "RT", "5BU", "S4U", "MTU", "MNU", "UMS", "IU", "UD5", "PYO",\
                 "SUR", "SSU", "RU", "RU5", "RU3", "FMN",\
                 "0U", "DU"]

spcnucs = ['PSU', 'FHU', 'QUO', 'N6G']
varnucs = sorted(set(sum(variants.values(), list())))
modnucs = spcnucs + varnucs
allnucs = nucleotides + spcnucs + varnucs + ["N"]

hnames['ribose'] = {"C5'": ["H5'", "H5''"], "O5'": ["H5T"], "C4'": ["H4'"],\
                    "C3'": ["H3'"], "O3'": ["H3T"], "C2'": ["H2''"],\
                    "O2'": ["H2'"], "C1'": ["H1'"]}
hnames['A'] = {"C2": ["H2"], "N6": ["H61", "H62"], "C8": ["H8"]}
hnames['C'] = {"N4": ["H41", "H42"], "C5": ["H5"], "C6": ["H6"]}
hnames['G'] = {"N1": ["H1"], "N2": ["H21", "H22"], "C8": ["H8"]}
hnames['U'] = {"N3": ["H3"], "C5": ["H5"], "C6": ["H6"]}

for N in nucleotides:
    for ringclass in [purines, pyrimidines]:
        if N in ringclass:
            ringclass += variants[N]
    for base in variants[N]:
        rings[base] = rings[N]
        bonds[base] = bonds[N]
        hnames[base] = hnames[N]
        angles[base] = angles[N]
        dihedrals[base] = dihedrals[N]

purines = sorted(set(purines))
pyrimidines = sorted(set(pyrimidines))

rings['PSU'] = rings['U']
bonds['PSU'] = bonds['U'][:]
bonds['PSU'][bonds['PSU'].index(("C1'", "N1"))] = ("C1'", "C4")
angles['PSU'] = angles['U']
dihedrals['PSU'] = dihedrals['ribose'] + [("O4'", "C1'", "C4", "C2")]

rings['FHU'] = rings['PSU']
bonds['FHU'] = bonds['PSU']
angles['FHU'] = angles['PSU']
dihedrals['FHU'] = dihedrals['PSU']

rings['N6G'] = rings['G']
bonds['N6G'] = bonds['G'][:]
bonds['N6G'][bonds['N6G'].index(("C6", "O6"))] = ("C6", "N6")
angles['N6G'] = angles['G']
angles['N6G'][angles['N6G'].index(("N1", "C6", "O6"))] = ("N1", "C6", "N6")
angles['N6G'][angles['N6G'].index(("C5", "C6", "O6"))] = ("C5", "C6", "N6")
dihedrals['N6G'] = dihedrals['G']

rings['QUO'] = rings['G']
bonds['QUO'] = bonds['G'][:]
bonds['QUO'][bonds['QUO'].index(("C5", "N7"))] = ("C5", "C7")
bonds['QUO'][bonds['QUO'].index(("N7", "C8"))] = ("C7", "C8")
angles['QUO'] = angles['G'][:]
angles['QUO'][angles['QUO'].index(('C6', 'C5', 'N7'))] = ('C6', 'C5', 'C7')
angles['QUO'][angles['QUO'].index(('C4', 'C5', 'N7'))] = ('C4', 'C5', 'C7')
angles['QUO'][angles['QUO'].index(('C5', 'N7', 'C8'))] = ('C5', 'C7', 'C8')
angles['QUO'][angles['QUO'].index(('N7', 'C8', 'N9'))] = ('C7', 'C8', 'N9')
dihedrals['QUO'] = dihedrals['G']

distance = dict()
distance[('C', 'C')] = 1.388
for bond in [('C', 'N'), ('N', 'C')]: distance[bond] = 1.334

def bondtype(bond):
    atm1, atm2 = bond
    atm1 = 'C' if 'C' in atm1 else 'N'
    atm2 = 'C' if 'C' in atm2 else 'N'
    return (atm1, atm2)


homa_alpha = dict()
homa_alpha[('C', 'C')] = 257.7
for bond in [('C', 'N'), ('N', 'C')]: homa_alpha[bond] = 93.52

deprotonated = {'+': ['W', 'w'], 'z': ['S', 's']}

for N in nucleotides+spcnucs:
    edges[N] = dict()

edges['A']['h'] = ['C5', 'C6', 'N6', 'N7']
edges['A']['w'] = ['N1', 'C2', 'C6', 'N6']
edges['A']['s'] = ['C2', 'N3', 'C4', "C2'", "O2'"]
edges['G']['h'] = ['C5', 'C6', 'O6', 'N7']
edges['G']['w'] = ['N1', 'C2', 'N2', 'C6', 'O6']
edges['G']['s'] = ['C2', 'N2', 'N3', 'C4', "C2'", "O2'"]
edges['C']['h'] = ['C4', 'N4', 'C5', 'C6']
edges['C']['w'] = ['C2', 'O2', 'N3', 'C4', 'N4']
edges['C']['s'] = ['C2', 'O2', "C1'", "O2'"]
edges['U']['h'] = ['C4', 'O4', 'C5', 'C6']
edges['U']['w'] = ['C2', 'O2', 'N3', 'C4', 'O4']
edges['U']['s'] = ['C2', 'O2', "C1'", "O2'"]

edges['PSU']['h'] = ['N1', 'C2', 'O2', 'C6']
edges['PSU']['s'] = ['C2', 'O2', 'N3', 'C4', 'O4']
edges['PSU']['w'] = ['C4', 'O4', "C1'", "O2'"]

edges['FHU'] = {E: edges['PSU'][E][:] for E in edges['PSU'].keys()}
edges['FHU']['h'].append('O6')

edges['N6G'] = {E: edges['G'][E][:] for E in edges['G'].keys()}
edges['N6G']['w'][edges['N6G']['w'].index('O6')] = 'N6'
edges['N6G']['h'][edges['N6G']['h'].index('O6')] = 'N6'

edges['QUO'] = {E: edges['G'][E][:] for E in edges['G'].keys()}
edges['QUO']['h'][edges['QUO']['h'].index('N7')] = 'C7'

for N in nucleotides+spcnucs:
    for E in ['h', 'w', 's']:
        edges[N][E.upper()] = heavy_atoms(edges[N][E])
    edges[N]['+'] = edges[N]['w']
    edges[N]['z'] = edges[N]['s']


def base_pair_bonds(bptype):
    n1 = bptype[0]
    n2 = bptype[1]
    for N in nucleotides+spcnucs:
        if bptype[0] in variants[N]: n1 = N
        if bptype[1] in variants[N]: n2 = N
    return [("C1'", "C1'")] + list(product(edges[n1][bptype[2]],\
                                           edges[n2][bptype[3]]))


atoms['ribose'] = sorted(set(list(sum(bonds['ribose'], tuple()))), key=lambda x: x[1::-1])

for N in nucleotides+spcnucs:
    atoms[N] = sorted(set(list(sum(bonds[N], tuple()))+backbone+phosphate[:3]), key=lambda x: ("P" in x, "'" in x, x[1::-1]))

for N in nucleotides:
    for base in variants[N]:
        atoms[base] = atoms[N]

donors = dict()
accors = dict()

# for hetatm in hetatoms[1:]:
#     donors[hetatm] = [hetatm]

donors["ribose"] = ["O2'", "C1'", "C2'", "C3'", "C4'", "C5'"]
accors["ribose"] = ["O2'", "O3'", "O4'", "O5'"]

donors['N'] = donors["ribose"]
accors['N'] = accors["ribose"]+phosphate[:3]

donors['R'] = donors['N'] + ["C8"]
donors['A'] = donors['R'] + ["C2", "N6"]
donors['G'] = donors['R'] + ["N1", "N2"]

accors['R'] = accors['N'] + ["N3", "N7"]
accors['A'] = accors['R'] + ["N1", "N6"]
accors['G'] = accors['R'] + ["N2", "O6"]


donors['Y'] = donors['N'] + ["C5", "C6"]
donors['C'] = donors['Y'] + ["N4"]
donors['U'] = donors['Y'] + ["N3"]

accors['Y'] = accors['N'] + ["O2"]
accors['C'] = accors['Y'] + ["N3", "N4"]
accors['U'] = accors['Y'] + ["O4"]


hnames = defaultdict(dict)
hnames['ribose'] = {"C5'": ["H5'", "H5''"], "O5'": ["H5T"], "C4'": ["H4'"],\
                    "C3'": ["H3'"], "O3'": ["H3T"], "C2'": ["H2''"],\
                    "O2'": ["H2'"], "C1'": ["H1'"]}
hnames['N'] = hnames["ribose"]
hnames['R'] = {"C8": ["H8"]}
hnames['A'] = {"C2": ["H2"], "N6": ["H61", "H62"]}
hnames['G'] = {"N2": ["H21", "H22"], "N1": ["H1"]}
hnames['Y'] = {"C5": ["H5"], "C6": ["H6"]}
hnames['C'] = {"N4": ["H41", "H42"]}
hnames['U'] = {"N3": ["H3"]}

hnames['R'].update(hnames['N'])
hnames['A'].update(hnames['R'])
hnames['G'].update(hnames['R'])
hnames['Y'].update(hnames['N'])
hnames['C'].update(hnames['Y'])
hnames['U'].update(hnames['Y'])

for N, varlist in zip(['A', 'C', 'G', 'U'], [['ADE', '0A'], ['CYT', '0C'], ['GUA', '0G'], ['URA', '0U']]):
    for variant in varlist:
        donors[variant] = donors[N][:]
        accors[variant] = accors[N][:]
        hnames[variant] = hnames[N].copy()

for N, variant in zip(['A', 'C', 'G',], ['DA', 'DC', 'DG']):
    donors[variant] = donors[N][:]
    accors[variant] = accors[N][:]
    hnames[variant] = hnames[N].copy()
    donors[variant].remove("O2'")
    accors[variant].remove("O2'")
    del hnames[variant]["O2'"]
