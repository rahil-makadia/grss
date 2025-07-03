"""ADES (Astrometric Data Exchange Standard) data handling for the GRSS orbit determination code"""

__all__ = [ 'ades_keep_columns',
            'ades_add_columns',
            'ades_grss_columns',
            'ades_column_types',
            'ades_catalog_map',
            'pack_letters',
            'unpack_letters',
            'prog_codes',
]

# from table 16 of ADES_Description.pdf from the ADES-Master git repo
# comments correspond to future GRSS capabilities.
ades_keep_columns = {
    'permID':'str', 'provID': 'str', 'mode': 'str', 'stn': 'str', 'prog': 'str',
    'obsTime': 'str', 'rmsTime': 'float',
    'ra': 'float', 'dec': 'float',
    'rmsRA': 'float', 'rmsDec': 'float', 'rmsCorr': 'float', 'astCat': 'str',
    'trx': 'str', 'rcv': 'str', 'delay': 'float', 'doppler': 'float',
    'rmsDelay': 'float', 'rmsDoppler': 'float', 'com': 'Int64', 'frq': 'float',
    'raStar': 'float', 'decStar': 'float', 'deltaRA': 'float', 'deltaDec': 'float',
    'sys': 'str', 'ctr': 'Int64', 'pos1': 'float', 'pos2': 'float', 'pos3': 'float',
    'vel1': 'float', 'vel2': 'float', 'vel3': 'float', 'posCov11': 'float', 'posCov12': 'float',
    'posCov13': 'float', 'posCov22': 'float', 'posCov23': 'float', 'posCov33': 'float'
}
ades_add_columns = {
    'resRA': 'float', 'resDec': 'float', 'selAst': 'str',
    'sigRA': 'float', 'sigDec': 'float', 'sigCorr': 'float', 'sigTime': 'float',
    'biasRA': 'float', 'biasDec': 'float', 'biasTime': 'float',
    'resDelay': 'float', 'resDoppler': 'float',
    'sigDelay': 'float', 'sigDoppler': 'float',
}
ades_grss_columns = {
    'obsTimeMJD': 'float', 'obsTimeMJDTDB': 'float', 'cosDec': 'float',
    'resChi': 'float',
}
ades_column_types = ades_keep_columns | ades_add_columns | ades_grss_columns

# from ADES-Master/Python/bin/packUtil.py
ades_catalog_map = {
    ' ': 'UNK',
    'a': 'USNOA1',
    'b': 'USNOSA1',
    'c': 'USNOA2',
    'd': 'USNOSA2',
    'e': 'UCAC1',
    'f': 'Tyc1',
    'g': 'Tyc2',
    'h': 'GSC1.0',
    'i': 'GSC1.1',
    'j': 'GSC1.2',
    'k': 'GSC2.2',
    'l': 'ACT',
    'm': 'GSCACT',
    'n': 'SDSS8',
    'o': 'USNOB1',
    'p': 'PPM',
    'q': 'UCAC4',
    'r': 'UCAC2',
    's': 'USNOB2',  # USNOB2 missing on ADES web page
    't': 'PPMXL',
    'u': 'UCAC3',
    'v': 'NOMAD',
    'w': 'CMC14',
    'x': 'Hip2',
    'y': 'Hip1',
    'z': 'GSC',
    'A': 'AC',
    'B': 'SAO1984',
    'C': 'SAO',
    'D': 'AGK3',
    'E': 'FK4',
    'F': 'ACRS',
    'G': 'LickGas',
    'H': 'Ida93',
    'I': 'Perth70',
    'J': 'COSMOS',
    'K': 'Yale',
    'L': '2MASS',
    'M': 'GSC2.3',
    'N': 'SDSS7',
    'O': 'SSTRC1',
    'P': 'MPOSC3',
    'Q': 'CMC15',
    'R': 'SSTRC4',
    'S': 'URAT1',
    'T': 'URAT2',  # URAT2 missing on ADES web page
    'U': 'Gaia1',
    'V': 'Gaia2',
    'W': 'Gaia3',
    'X': 'Gaia3E',  
    'Y': 'UCAC5',  # UCAC5 mission on ADES web page
    'Z': 'ATLAS2', 
    '0': 'IHW', 
    '1': 'PS1_DR1', 
    '2': 'PS1_DR2', 
    '3': 'Gaia_Int', 
    '4': 'GZ', 
    '5': 'UBSC', 
    '6': 'Gaia_2016', 
}
# flip the dictionary to get the catalog codes
ades_catalog_map = {v: k for k, v in ades_catalog_map.items()}

# from ADES-Master/Python/bin/packUtil.py
pack_letters = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'
unpack_letters = {pack_letters[i]: i for i in range(len(pack_letters))}
prog_codes = R"""0123456789!"#$%&'()*+,-./[\]^_`{|}~:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz£"""

special_codes = {
    'gaia': {'258'},
    'occultation': {'244', '275'},
    'spacecraft': {'S/C', 'S_C', '245', '249', '250', '273', '274',
                   'C49', 'C50', 'C51', 'C52', 'C53', 'C54', 'C55', 'C56', 'C57', 'C58', 'C59', },
    'roving': {'247', '270'},
}
