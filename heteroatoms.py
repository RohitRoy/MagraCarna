cations = ["CA", "SR", "MG", "K", "NA", "BA", "CD", "PT", "OHX", "ZN", "MN", "IRI"]
anions = ["SO4", "CL"]

neutral = ["HOH", "O", "NH3", "SF4"]

organic = ["BTN", "004", "04X"]

hetresdus = cations + anions + organic + neutral

donors = {"HOH": "O", "O": "O"}
accors = {"HOH": "O", "O": "O"}
hnames = {"HOH": {"O": ["H1", "H2"]}, "O": {"O": []}}

someothers = ["0EC",  "0G6", "0TD",  "0U1", "10C", "13T", "16D", "1CC", "1F2", "1F3", "1PE", 
                "1RN", "1TU",
              "20V", "218", "23G", "29G", "29H", "2AU", "2BA", "2BP", "2HP", "2IA", "2OP", "2PE",
                "2QB", "2QC", "2QY", "2QZ", "2R1", "2R3", "2SG", "2TB", "2ZE", "2ZY", "2ZZ", 
              "31H", "31M", "34G", "365", "38E", "3AD", "3AT", "3AW", "3AY", "3CO", "3DR", "3H3",
                "3HE", "3J2", "3J6", "3K5", "3K8", "3KD", "3KF", "3L2", "3TD", "3TS", "3V6",
              "42B", "4BW", "4D4", "4L4", "4M2", "4OC", "50L", "50N", "574", "5AZ", "5CF", "5CM",
                "5FU", "5GP", "5GS", "5IU", "5OH", "5PC",
              "6AP", "6FC", "6FU", "6GO", "6GS", "6GU", "6HA", "6HC", "6HG", "6HS", "6HT", "6MN",
                "773", "7AT", "7DG", "84T", "8AN", "8XA", "8XC", "8XG", "8XU", "9DG",
              "A2F", "A2L", "A2P", "A3P", "A44", "A5A", "A5L", "A5M", "A5O", "A6A", "A6C", "A6G",
                "A6U", "A7E", "A9Z", "AB6", "AB9", "ACA", "ACE", "ACP", "ACT", "ACY", "ADN", 
                "ADS", "AF2", "AF3", "AG2", "AG9", "AGS", "AKN", "ALF", "AM2", "AMZ", "ANM",
                "ANP", "AP7", "APN", "ARF",  "AS", "AT7", "ATL", "AU3", 
              "B12", "B1Z", "B3P",  "BB9", "BCM", "BDG", "BDR", "BEF", "BFT", "BGC", "BGM", "BGR",
                "BMA", "BME", "BO4",  "BR", "BRU", "BU1", 
              "C2E", "C2L", "C43", "C5L", "C5P", "C6G", "C7P", "CAC", "CAI", "CAR", "CDP", "CFL",
                "CFZ", "CG1",  "CH", "CH1", "CIR", "CIT", "CLM", "CLY", "CM0", "CMY", "CNC",
                "CNY",  "CO", "CPN",  "CS", "CSG", "CSX", "CTC", "CTP", "CTY",  "CU", "CYY",
              "D2T", "D2X", "DAI", "DAL", "DAR", "DBB", "DBU",  "DCT", "DCY",  "DGP", "DGT",
                "DHA", "DIO", "DJF", "DOC", "DOL", "DPO", "DPP", "DPR", "DST",  "DTP", "DUT",
                "DX4",
              "EDE", "EDO", "EEM", "EFZ", "EM1", "EMK", "ENX", "EOH", "EPE", "ERY", 
              "F", "F3N", "F3O",  "FB", "FE2", "FFO", "FGA", "FLC", "FME", "FMT", "FOU",
                "FOZ", "FPD", "FUA", "FUC", "FUL", "FYA",
              "G0B", "G19", "G2L", "G2P", "G46", "G48", "G6P", "GAL", "GAO", "GAU", "GCP", "GDO",
                "GE1", "GE2", "GE3", "GET", "GF2", "GFL", "GH3", "GIR", "GLA", "GLP", "GND",
                "GNG", "GOL", "GPN", "GRB",  "GS", "GSU", "GUN",
              "H4B", "HFA",  "HG", "HMT", "HPA", "HRG", "HYG", "HYP", "I2A",  "IC", "IDG", "IEL",
                "IG", "ILA", "ILX", "IMD", "IOD", "IPA", "IPH", "IR3", "ISH", "ISI",  "IUM",
              "JOS", "JS4", "JS5", "JS6", "KAN", "KBE", "KIR", "KSG", 
              "L8H", "L94", "LC2", "LCA", "LCC", "LCG", "LHA", "LIV", "LKC", "LLL", "LMA", "LMS",
                "LU", "LYA",
              "M2M", "M3O", "M5M", "M5Z", "MA6", "MAN", "MEA", "MEQ", "MES", "MGR", "MGT", "MH6",
                "MHT", "MHU", "MHV", "MHW", "MIX", "MLI", "MLZ", "MMC", "MMT", "MPD", "MRC",
                "MRD", "MSE", "MSP", "MTT", "MUB", "MUL", "MVA", "MYN",
              "N30", "N33", "N5C", "N5M", "N79", "NAG", "NCO", "NCU", "NDG", "NEB", "NEG", "NF2",
                "NH2", "NH4", "NHE",  "NI", "NME", "NMY", "NMZ", "NO3", "NPM", "NTT", "NVA",
                "NVP", "O2C", "OCS", "OLZ", "ON0", "ONE",  "OS",
              "P12", "P13", "P14", "P1P", "P24", "P5P", "P6G", "PA1", "PA2", "PAE", "PAR",  "PB", 
                "PCY", "PDI", "PDU", "PEG", "PG4", "PG6", "PGE", "PGN", "PHA", "PLR", "PMZ",
                "PO2", "PO4", "POP", "PPU", "PPV", "PQ0", "PQ1", "PRF", "PRI", "PRL", "PST", 
                "PT4", "PUP", "PUT", "PYI", "PYY", "QUA",
              "R14",  "RB", "RBF", "RHD", "RIB", "RIO", "ROS", "RPC", "RPO", "RS3", "RSP", 
                "RSQ", "RTP", "RUS",
              "S9L", "SAH", "SAR",  "SC", "SCM", "SCY", "SDG", "SE4", "SEP", "SFG", "SIN", "SIS",
                "SJP", "SLD", "SLZ", "SPD", "SPM", "SPR", "SPS", "SRA", "SRY", "SS0", "SSA",
                "STD", "SUC", "SVN",
              "T1C", "T2T", "T39", "T3P", "T8B", "TAC", "TAF", "TAO", "TAR",  "TB", "TCP", "TEL",
                "TEP", "THF",  "TL", "TLN", "TM2", "TOA", "TOB", "TOC", "TOY", "TPN", "TPO",
                "TPP", "TPS", "TRS", "TRX", "TS6", "TS9", "TSE", "TSP",  "TT", "TYE", "TYK",
              "U2L", "U33", "U34", "U36", "U37", "U3H", "U5P", "U8U", "UAL", "UAM", "UAR", "UBD", 
                "UDP", "UFT", "UNK", "UNX", "UPV", "URE", "URI", "URU", "US5", "UTP", "UZR",
              "VIB", "VIF", "VIR", "VO4", "WCP", "WIN", "WO2", "XAN", "XAR", "XCR", "XGR", "XTR",
                "XXX", "Y5P", "YMP", "ZBA", "ZIT", "ZLD", "ZUK"]