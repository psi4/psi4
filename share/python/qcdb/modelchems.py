

class QCEssential(object):
    """Class to link literature and external representation of some
    aspect of quantum chemistry (basis set, method, etc.) with a
    shorthand and indexed representation of same.

    """
    def __init__(self, name, fullname=None, latex=None, citation=None, pdfdatabase=None, comment=None):
        """

        """
        self.name = name.lower()
        self.fullname = fullname
        if fullname is not None and latex is None:
            self.latex = fullname
        else:
            self.latex = latex
        self.citation = citation
        self.dsdbid = pdfdatabase
        self.comment = comment

    def __str__(self):
        text = ''
        text += """  ==> %s BasisSet <==\n\n""" % (self.name)
        text += """  Formal name:          %s\n""" % (self.fullname)
        text += """  LaTeX representation: %s\n""" % (self.latex)
        text += """  PDF database id:      %s\n""" % (self.dsdbid)
        text += """  Literature citation:  %s\n""" % (self.citation)
        text += """  Comment:              %s\n""" % (self.comment)
        text += """\n"""
        return text


class BasisSet(QCEssential):
    """Specialization of :pyclass:`QCEssential` for basis sets.

    """
    def __init__(self, name, fullname=None, latex=None, citation=None, pdfdatabase=None, comment=None, zeta=None, build=None):
        QCEssential.__init__(self, name, fullname, latex, citation, pdfdatabase, comment)
        self.name = name.lower()
        self.zeta = zeta
        self.build = [[self.name]] if build is None else build

    def __str__(self):
        text = ''
        text += """  ==> %s BasisSet Treatment <==\n\n""" % (self.name)
        text += """  Formal name:          %s\n""" % (self.fullname)
        text += """  LaTeX representation: %s\n""" % (self.latex)
        text += """  Zeta:                 %s\n""" % (self.zeta)
        text += """  CBS build:            %s\n""" % (self.build)
        text += """  PDF database id:      %s\n""" % (self.dsdbid)
        text += """  Literature citation:  %s\n""" % (self.citation)
        text += """  Comment:              %s\n""" % (self.comment)
        text += """\n"""
        return text


class Method(QCEssential):
    """Specialization of :pyclass:`QCEssential` for quantum chemical methods.

    """
    def __init__(self, name, fullname=None, latex=None, citation=None, pdfdatabase=None, comment=None):
        QCEssential.__init__(self, name, fullname, latex, citation, pdfdatabase, comment)
        self.name = name.upper()

    def __str__(self):
        text = ''
        text += """  ==> %s Method <==\n\n""" % (self.name)
        text += """  Formal name:          %s\n""" % (self.fullname)
        text += """  LaTeX representation: %s\n""" % (self.latex)
        text += """  PDF database id:      %s\n""" % (self.dsdbid)
        text += """  Literature citation:  %s\n""" % (self.citation)
        text += """  Comment:              %s\n""" % (self.comment)
        text += """\n"""
        return text


class Error(QCEssential):
    """Specialization of :pyclass:`QCEssential` for measures of error.

    """
    def __init__(self, name, fullname=None, latex=None, citation=None, pdfdatabase=None, comment=None):
        QCEssential.__init__(self, name, fullname, latex, citation, pdfdatabase, comment)
        self.name = name.lower()

    def __str__(self):
        text = ''
        text += """  ==> %s Error Measure <==\n\n""" % (self.name)
        text += """  Formal name:          %s\n""" % (self.fullname)
        text += """  LaTeX representation: %s\n""" % (self.latex)
        text += """  PDF database id:      %s\n""" % (self.dsdbid)
        text += """  Literature citation:  %s\n""" % (self.citation)
        text += """  Comment:              %s\n""" % (self.comment)
        text += """\n"""
        return text


_tlist = [
    BasisSet('dz',         fullname='cc-pVDZ'),
    BasisSet('jadz',       fullname='jun-cc-pVDZ'),
    BasisSet('hadz',       fullname='heavy-aug-cc-pVDZ'),
    BasisSet('adz',        fullname='aug-cc-pVDZ'),
    BasisSet('addz',       fullname='aug-cc-pV(D+d)Z'),
    BasisSet('tz',         fullname='cc-pVTZ'),
    BasisSet('matz',       fullname='may-cc-pVTZ'),
    BasisSet('jatz',       fullname='jun-cc-pVTZ'),
    BasisSet('hatz',       fullname='heavy-aug-cc-pVTZ'),
    BasisSet('atz',        fullname='aug-cc-pVTZ'),
    BasisSet('qz',         fullname='cc-pVQZ'),
    BasisSet('aaqz',       fullname='apr-cc-pVQZ'),
    BasisSet('maqz',       fullname='may-cc-pVQZ'),
    BasisSet('jaqz',       fullname='jun-cc-pVQZ'),
    BasisSet('haqz',       fullname='heavy-aug-cc-pVQZ'),
    BasisSet('aqz',        fullname='aug-cc-pVQZ'),
    BasisSet('a5z',        fullname='aug-cc-pV5Z'),
    BasisSet('dtz',        fullname='cc-pVDTZ', build=[None, ['tz', 'dtz']]),
    BasisSet('jadtz',      fullname='jun-cc-pVDTZ', build=[None, ['jatz', 'jadtz']]),
    BasisSet('hadtz',      fullname='heavy-aug-cc-pVDTZ', build=[None, ['hatz', 'hadtz']]),
    BasisSet('adtz',       fullname='aug-cc-pVDTZ', build=[['adtz'], ['atz', 'adtz']]),
    BasisSet('tqz',        fullname='cc-pVTQZ', build=[None, ['qz', 'tqz']]),
    BasisSet('matqz',      fullname='may-cc-pVTQZ', build=[None, ['maqz', 'matqz']]),
    BasisSet('jatqz',      fullname='jun-cc-pVTQZ', build=[None, ['jaqz', 'jatqz']]),
    BasisSet('hatqz',      fullname='heavy-aug-cc-pVTQZ', build=[None, ['haqz', 'hatqz']]),
    BasisSet('atqz',       fullname='aug-cc-pVTQZ', build=[['atqz'], ['aqz', 'atqz']]),
    BasisSet('aq5z',       fullname='aug-cc-pVQ5Z', build=[['aq5z'], ['a5z', 'aq5z']]),
    BasisSet('a6z',        fullname='aug-cc-pV6Z'),
    BasisSet('a56z',       fullname='aug-cc-pV56Z', build=[['a56z'], ['a6z', 'a56z']]),
    BasisSet('atzdz',      fullname='[aTZ; D:DZ]', latex="""[aTZ; $\delta$:DZ]""",
        build=[None, None, ['atz', 'atz', 'dz']]),
    BasisSet('adtzdz',     fullname='[aDTZ; D:DZ]', latex="""[aDTZ; $\delta$:DZ]""",
        build=[None, None, ['atz', 'adtz', 'dz']]),
    BasisSet('atqzdz',     fullname='[aTQZ; D:DZ]', latex="""[aTQZ; $\delta$:DZ]""",
        build=[None, None, ['aqz', 'atqz', 'dz']]),
    BasisSet('atzjadz',    fullname='[aTZ; D:jaDZ]', latex="""[aTZ; $\delta$:jaDZ]""",
        build=[None, None, ['atz', 'atz', 'jadz']]),
    BasisSet('adtzjadz',   fullname='[aDTZ; D:jaDZ]', latex="""[aDTZ; $\delta$:jaDZ]""",
        build=[None, None, ['atz', 'adtz', 'jadz']]),
    BasisSet('atqzjadz',   fullname='[aTQZ; D:jaDZ]', latex="""[aTQZ; $\delta$:jaDZ]""",
        build=[None, None, ['aqz', 'atqz', 'jadz']]),
    BasisSet('atzhadz',    fullname='[aTZ; D:haDZ]', latex="""[aTZ; $\delta$:haDZ]""",
        build=[None, None, ['atz', 'atz', 'hadz']]),
    BasisSet('adtzhadz',   fullname='[aDTZ; D:haDZ]', latex="""[aDTZ; $\delta$:haDZ]""",
        build=[None, None, ['atz', 'adtz', 'hadz']]),
    BasisSet('atqzhadz',   fullname='[aTQZ; D:haDZ]', latex="""[aTQZ; $\delta$:haDZ]""",
        build=[None, None, ['aqz', 'atqz', 'hadz']]),
    BasisSet('atzadz',     fullname='[aTZ; D:aDZ]', latex="""[aTZ; $\delta$:aDZ]""",
        build=[None, None, ['atz', 'atz', 'adz']]),
    BasisSet('adtzadz',    fullname='[aDTZ; D:aDZ]', latex="""[aDTZ; $\delta$:aDZ]""",
        build=[None, None, ['atz', 'adtz', 'adz']]),
    BasisSet('atqzadz',    fullname='[aTQZ; D:aDZ]', latex="""[aTQZ; $\delta$:aDZ]""",
        build=[None, None, ['aqz', 'atqz', 'adz']]),
    BasisSet('aq5zadz',    fullname='[aQ5Z; D:aDZ]', latex="""[aQ5Z; $\delta$:aDZ]""",
        build=[None, None, ['a5z', 'aq5z', 'adz']]),
    BasisSet('atzdtz',     fullname='[aTZ; D:DTZ]', latex="""[aTZ; $\delta$:DTZ]""",
        build=[None, None, ['atz', 'atz', 'dtz']]),
    BasisSet('atqzdtz',    fullname='[aTQZ; D:DTZ]', latex="""[aTQZ; $\delta$:DTZ]""",
        build=[None, None, ['aqz', 'atqz', 'dtz']]),
    BasisSet('atzjadtz',   fullname='[aTZ; D:jaDTZ]', latex="""[aTZ; $\delta$:jaDTZ]""",
        build=[None, None, ['atz', 'atz', 'jadtz']]),
    BasisSet('atqzjadtz',  fullname='[aTQZ; D:jaDTZ]', latex="""[aTQZ; $\delta$:jaDTZ]""",
        build=[None, None, ['aqz', 'atqz', 'jadtz']]),
    BasisSet('atzhadtz',   fullname='[aTZ; D:haDTZ]', latex="""[aTZ; $\delta$:haDTZ]""",
        build=[None, None, ['atz', 'atz', 'hadtz']]),
    BasisSet('atqzhadtz',  fullname='[aTQZ; D:haDTZ]', latex="""[aTQZ; $\delta$:haDTZ]""",
        build=[None, None, ['aqz', 'atqz', 'hadtz']]),
    BasisSet('atzadtz',    fullname='[aTZ; D:aDTZ]', latex="""[aTZ; $\delta$:aDTZ]""",
        build=[None, None, ['atz', 'atz', 'adtz']]),
    BasisSet('atqzadtz',   fullname='[aTQZ; D:aDTZ]', latex="""[aTQZ; $\delta$:aDTZ]""",
        build=[None, None, ['aqz', 'atqz', 'adtz']]),
    BasisSet('aq5zadtz',   fullname='[aQ5Z; D:aDTZ]', latex="""[aQ5Z; $\delta$:aDTZ]""",
        build=[None, None, ['a5z', 'aq5z', 'adtz']]),
    BasisSet('atqztz',     fullname='[aTQZ; D:TZ]', latex="""[aTQZ; $\delta$:TZ]""",
        build=[None, None, ['aqz', 'atqz', 'tz']]),
    BasisSet('atqzjatz',   fullname='[aTQZ; D:jaTZ]', latex="""[aTQZ; $\delta$:jaTZ]""",
        build=[None, None, ['aqz', 'atqz', 'jatz']]),
    BasisSet('atqzmatz',   fullname='[aTQZ; D:maTZ]', latex="""[aTQZ; $\delta$:maTZ]""",
        build=[None, None, ['aqz', 'atqz', 'matz']]),
    BasisSet('atqzhatz',   fullname='[aTQZ; D:haTZ]', latex="""[aTQZ; $\delta$:haTZ]""",
        build=[None, None, ['aqz', 'atqz', 'hatz']]),
    BasisSet('atqzatz',    fullname='[aTQZ; D:aTZ]', latex="""[aTQZ; $\delta$:aTZ]""",
        build=[None, None, ['aqz', 'atqz', 'atz']]),
    BasisSet('aq5zatz',    fullname='[aQ5Z; D:aTZ]', latex="""[aQ5Z; $\delta$:aTZ]""",
        build=[None, None, ['a5z', 'aq5z', 'atz']]),
    BasisSet('dzf12',      fullname='cc-pVDZ-F12'),
    BasisSet('tzf12',      fullname='cc-pVTZ-F12'),
    BasisSet('qzf12',      fullname='cc-pVQZ-F12'),
    BasisSet('dtzf12',     fullname='cc-pVDTZ-F12', build=[['dtzf12'], ['tzf12', 'dtzf12']]),
    BasisSet('tqzf12',     fullname='cc-pVTQZ-F12', build=[['tqzf12'], ['qzf12', 'tqzf12']]),
    BasisSet('hill1_adtz', build=[['hillcc_adtz'], ['atz', 'hillcc_adtz']]),  # TODO should have None or non-xtpl first element?
    BasisSet('hill1_atqz', build=[['hillcc_atqz'], ['aqz', 'hillcc_atqz']]),
    BasisSet('hill1_aq5z', build=[['hillcc_aq5z'], ['a5z', 'hillcc_aq5z']]),
    BasisSet('hill1_dtzf12', build=[['hillcc_dtzf12'], ['tzf12', 'hillcc_dtzf12']]),
    BasisSet('hill1_tqzf12', build=[['hillcc_tqzf12'], ['qzf12', 'hillcc_tqzf12']]),
    BasisSet('hill2_dtzf12', build=[None, None, ['tzf12', 'hillcc_dtzf12', 'hillt_dtzf12']]),
    BasisSet('hill2_tqzf12', build=[None, None, ['qzf12', 'hillcc_tqzf12', 'hillt_tqzf12']]),
    BasisSet('hill2_adtz', build=[None, None, ['atz', 'hillcc_adtz', 'hillt_adtz']]),
    BasisSet('hill2_atqz', build=[None, None, ['aqz', 'hillcc_atqz', 'hillt_atqz']]),
    BasisSet('hill2_aq5z', build=[None, None, ['a5z', 'hillcc_aq5z', 'hillt_aq5z']]),
    BasisSet('dadz',       fullname='double-aug-cc-pVDZ'),
    BasisSet('datz',       fullname='double-aug-cc-pVTZ'),
    BasisSet('631pgs',     fullname='6-31+G(d)'),
    BasisSet('6311pg_3df_2p_', fullname='6-311+G(3df,2p)'),
    BasisSet('6311ppg_3df_2p_', fullname='6-311++G(3df,2p)'),
]
bases = {}
for item in _tlist:
    bases[item.name] = item

# Key name must be [A-Z], [0-9], and _, being either all upper or all lowercase according to Essential
# fullname can be anything on the keyboard, no ascii codes
# latex can contain escape codes for LaTeX
_tlist = [
    Method('SAPT0',           fullname='SAPT0'),
    Method('SAPT0S',          fullname='sSAPT0', latex="""\\textit{s}SAPT0"""),
    Method('SAPTSCS',         fullname='SCS-SAPT0'),
    Method('SAPTDFT',         fullname='DFT-SAPT'),
    Method('SAPT2',           fullname='SAPT2'),
    Method('SAPT2P',          fullname='SAPT2+'),
    Method('SAPT3',           fullname='SAPT2+(3)'),
    Method('SAPT3F',          fullname='SAPT2+3'),
    Method('SAPT2PC',         fullname='SAPT2+(CCD)'),
    Method('SAPT3C',          fullname='SAPT2+(3)(CCD)'),
    Method('SAPT3FC',         fullname='SAPT2+3(CCD)'),
    Method('SAPT2PM',         fullname='SAPT2+dMP2', latex="""SAPT2+$\delta$MP2"""),
    Method('SAPT3M',          fullname='SAPT2+(3)dMP2', latex="""SAPT2+(3)$\delta$MP2"""),
    Method('SAPT3FM',         fullname='SAPT2+3dMP2', latex="""SAPT2+3$\delta$MP2"""),
    Method('SAPT2PCM',        fullname='SAPT2+(CCD)dMP2', latex="""SAPT2+(CCD)$\delta$MP2"""),
    Method('SAPT3CM',         fullname='SAPT2+(3)(CCD)dMP2', latex="""SAPT2+(3)(CCD)$\delta$MP2"""),
    Method('SAPT3FCM',        fullname='SAPT2+3(CCD)dMP2', latex="""SAPT2+3(CCD)$\delta$MP2"""),
    Method('SAPT2LCM',        fullname='MP2(CCD)', comment="""Identical to SAPT2+(CCD)dMP2"""),
    Method('HF',              fullname='HF'),
    Method('MP2',             fullname='MP2'),
    Method('SCSMP2',          fullname='SCS-MP2'),
    Method('SCSNMP2',         fullname='SCS(N)-MP2'),
    Method('SCSMIMP2',        fullname='SCS(MI)-MP2'),
    Method('DWMP2',           fullname='DW-MP2'),
    Method('MP2C',            fullname='MP2C'),
    Method('MP3',             fullname='MP3'),
    Method('MP25',            fullname='MP2.5'),
    Method('CCSD',            fullname='CCSD'),
    Method('SCSCCSD',         fullname='SCS-CCSD'),
    Method('SCSMICCSD',       fullname='SCS(MI)-CCSD'),
    Method('CCSDT',           fullname='CCSD(T)'),
    Method('HFCABS',          fullname='HF-CABS'),
    Method('MP2F12',          fullname='MP2-F12'),
    Method('SCSMP2F12',       fullname='SCS-MP2-F12'),
    Method('SCSNMP2F12',      fullname='SCS(N)-MP2-F12'),
    Method('SCSMIMP2F12',     fullname='SCS(MI)-MP2-F12'),
    Method('DWMP2F12',        fullname='DW-MP2-F12'),
    Method('MP2CF12',         fullname='MP2C-F12'),
    Method('CCSDAF12',        fullname='CCSD-F12a'),
    Method('CCSDBF12',        fullname='CCSD-F12b'),
    Method('CCSDCF12',        fullname='CCSD-F12c'),
    Method('SCSCCSDAF12',     fullname='SCS-CCSD-F12a'),
    Method('SCSCCSDBF12',     fullname='SCS-CCSD-F12b'),
    Method('SCSCCSDCF12',     fullname='SCS-CCSD-F12c'),
    Method('SCMICCSDAF12',    fullname='SCS(MI)-CCSD-F12a'),
    Method('SCMICCSDBF12',    fullname='SCS(MI)-CCSD-F12b'),
    Method('SCMICCSDCF12',    fullname='SCS(MI)-CCSD-F12c'),
    Method('CCSDTAF12',       fullname='CCSD(T**)-F12a'),
    Method('CCSDTBF12',       fullname='CCSD(T**)-F12b'),
    Method('CCSDTCF12',       fullname='CCSD(T**)-F12c'),
    Method('DWCCSDTF12',      fullname='DW-CCSD(T**)-F12'),
#        build=lambda: ['DW-CCSD(T**)-F12 TOTAL ENERGY'],
#        ['HF-CABS TOTAL ENERGY', 'DW-CCSD(T**)-F12 CORRELATION ENERGY'],
#        ['HF-CABS TOTAL ENERGY', 'MP2-F12 CORRELATION ENERGY', 'DW-CCSD(T**)-F12 CC CORRECTION ENERGY'],
#        ['HF-CABS TOTAL ENERGY', 'MP2-F12 CORRELATION ENERGY', 'DW-CCSD-F12 CC CORRECTION ENERGY', 'DW-(T**)-F12 CORRECTION ENERGY'])
    Method('B97',             fullname='B97'),
    Method('B97D2',           fullname='B97-D2'),
    Method('B97D3',           fullname='B97-D3'),
    Method('B97D3BJ',         fullname='B97-D3bj'),
    Method('B3LYP',           fullname='B3LYP'),
    Method('B3LYPD2',         fullname='B3LYP-D2'),
    Method('B3LYPD3',         fullname='B3LYP-D3'),
    Method('B3LYPD3',         fullname='B3LYP-D3'),
    Method('B3LYPXDM',        fullname='B3LYP-XDM'),
    Method('B2PLYP',          fullname='B2PLYP'),
    Method('B2PLYPD2',        fullname='B2PLYP-D2'),
    Method('B2PLYPD3',        fullname='B2PLYP-D3'),
    Method('B2PLYPD3BJ',      fullname='B2PLYP-D3bj'),
    Method('M052X',           fullname='M05-2X'),
    Method('M052XD3',         fullname='M05-2X-D3'),
    Method('M062X',           fullname='M06-2X'),
    Method('M062XD3',         fullname='M06-2X-D3'),
    Method('M08HX',           fullname='M08-HX'),
    Method('M08SO',           fullname='M08-SO'),
    Method('M11',             fullname='M11'),
    Method('M11L',            fullname='M11L'),
    Method('XYG3',            fullname='XYG3'),
    Method('DLDFD',           fullname='dlDF+D'),
    Method('DSDPBEP86',       fullname='DSD-PBEP86'),
    Method('VV10',            fullname='VV10'),
    Method('LCVV10',          fullname='LC-VV10'),
    Method('WB97XD',          fullname='wB97X-D', latex="""$\omega$B97X-D"""),
    Method('WB97X2',          fullname='wB97X-2', latex="""$\omega$B97X-2"""),
    Method('PBE',             fullname='PBE'),
    Method('PBED2',           fullname='PBE-D2'),
    Method('PBED3',           fullname='PBE-D3'),
    Method('PBED3BJ',         fullname='PBE-D3bj'),
    Method('PBE0',            fullname='PBE0'),
    Method('PBE0D2',          fullname='PBE0-D2'),
    Method('PBE0D3',          fullname='PBE0-D3'),
    Method('PBE0D3BJ',        fullname='PBE0-D3BJ'),
    Method('PBE02',           fullname='PBE0-2'),
    Method('CCSDTNSAF12',     fullname='CCSD(T)-F12a'),
    Method('CCSDTNSBF12',     fullname='CCSD(T)-F12b'),
    Method('CCSDTNSCF12',     fullname='CCSD(T)-F12c'),
    Method('B970',            fullname='B970'),
    Method('B970D2',          fullname='B970-D2'),
    Method('BP86',            fullname='BP86'),
    Method('BP86D2',          fullname='BP86-D2'),
    Method('BP86D3',          fullname='BP86-D3'),
    Method('BP86D3BJ',        fullname='BP86-D3bj'),
    Method('CCSDTQ',          fullname='CCSDT(Q)'),
    Method('CCSDFULLT',       fullname='CCSDT'),
    Method('CCSDTSAF12',      fullname='CCSD(T*)-F12a'),
    Method('CCSDTSBF12',      fullname='CCSD(T*)-F12b'),
    Method('CCSDTSCF12',      fullname='CCSD(T*)-F12c'),
    Method('DWCCSDTNSF12',    fullname='DW-CCSD(T)-F12'),
    Method('DWCCSDTSF12',     fullname='DW-CCSD(T*)-F12'),
    Method('DELTQ',           fullname='d(TQ)', latex="""$\delta$(TQ)"""),  # TODO kill this once non-IE impl in reap-DB
    Method('DEL2T',           fullname='d(T)', latex="""$\delta$(T)"""),  # TODO kill this once non-IE impl in reap-DB
    ]

methods = {}
for item in _tlist:
    methods[item.name] = item

_tlist = [
    Error('maxe',            fullname='maxE'),
    Error('mine',            fullname='minE'),
    Error('me',              fullname='ME'),
    Error('mae',             fullname='MAE'),
    Error('rmse',            fullname='rmsE'),
    Error('stde',            fullname='stdE'),
    Error('maxpe',           fullname='maxPE'),
    Error('minpe',           fullname='minPE'),
    Error('mpe',             fullname='MPE'),
    Error('mape',            fullname='MAPE', latex="""MA\%E"""),
    Error('rmspe',           fullname='rmsPE'),
    Error('stdpe',           fullname='stdPE'),
    Error('maxpbe',          fullname='maxPBE'),
    Error('minpbe',          fullname='minPBE'),
    Error('mpbe',            fullname='MPBE'),
    Error('mapbe',           fullname='MAPBE', latex="""MA\%BE"""),
    Error('rmspbe',          fullname='rmsPBE'),
    Error('stdpbe',          fullname='stdPBE'),
]
errors = {}
for item in _tlist:
    errors[item.name] = item
