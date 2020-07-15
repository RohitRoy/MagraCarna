from collections import defaultdict

donors = defaultdict(list)
accors = defaultdict(list)
hnames = defaultdict(dict)

aminoacids = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY",\
              "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",\
              "THR", "TRP", "TYR", "VAL"]

donors['aminoacid'] = ["N"]
accors['aminoacid'] = ["O"]

# REF : http://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/charge/#hydrogen

donors['ARG'] = ["NE", "NH1", "NH2"]
donors['ASN'] = ["ND2"]
# donors['ASP'] = []
donors['GLN'] = ["NE2"]
# donors['GLU'] = []
donors['HIS'] = ["ND1", "NE2"]
donors['LYS'] = ["NZ"]
donors['SER'] = ["OG"]
donors['THR'] = ["OG1"]
donors['TRP'] = ["NE1"]
donors['TYR'] = ["OH"]

# accors['ARG'] = []
accors['ASN'] = ["OD1"]
accors['ASP'] = ["OD1", "OD2"]
accors['GLN'] = ["OE1"]
accors['GLU'] = ["OE1", "OE2"]
accors['HIS'] = ["ND1", "NE2"]
# accors['LYS'] = []
accors['SER'] = ["OG"]
accors['THR'] = ["OG1"]
# accors['TRP'] = []
accors['TYR'] = ["OH"]

for aacid in aminoacids:
    donors[aacid] += donors['aminoacid']
    accors[aacid] += accors['aminoacid']


hnames['aminoacid'] = {"N": ["H1", "H2", "H3", "H"], "CA": ["HA"]}

# hnames['ALA'] = {""}
hnames['ARG'] = {"NE": ["HE"], "NH1": ["HH11", "HH12"], "NH2": ["HH21", "HH22"]}
# hnames['ASP']
hnames['ASN'] = {"ND2": ["HD21", "HD22"]}
# hnames['CYS']
# hnames['GLU']
hnames['GLN'] = {"NE2": ["HE21", "HE22"]}
# hnames['GLY']
hnames['HIS'] = {"ND1": ["HD1"], "NE2": ["HE2"]}
# hnames['ILE']
# hnames['LEU']
hnames['LYS'] = {"NZ": ["HZ1", "HZ2", "HZ3"]}
# hnames['MET']
# hnames['PHE']
# hnames['PRO']
hnames['SER'] = {"OG": ["HG"]}
hnames['THR'] = {"OG1": ["HG1"]}
hnames['TRP'] = {"NE1": ["HE1"]}
hnames['TYR'] = {"OH": ["HH"]}
# hnames['VAL']

for aacid in aminoacids:
    hnames[aacid].update(hnames['aminoacid'])

donors = dict(donors)
accors = dict(accors)
