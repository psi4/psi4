#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2017 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

from __future__ import absolute_import
from __future__ import print_function
import struct


def getrec(reclabelarray, verbose=False):
    """Reads binary files JOBARC and JAINDX and returns contents
    of each record in *reclabelarray*.

    """
    knownlabels = {
        "AU_LENGT": 'DOUBLE',
        "CHARGE_E": 'DOUBLE',
        "AMU     ": 'DOUBLE',
        "NUC_MAGN": 'DOUBLE',
        "MASS_ELE": 'DOUBLE',
        "MASS_PRO": 'DOUBLE',
        "HBAR    ": 'DOUBLE',
        "AU_MASSP": 'DOUBLE',
        "SP_LIGHT": 'DOUBLE',
        "AU_EV   ": 'DOUBLE',
        "AVOGADRO": 'DOUBLE',
        "AU_ENERG": 'DOUBLE',
        "AU_CM-1 ": 'DOUBLE',
        "CM-1_KCA": 'DOUBLE',
        "CM-1_KJ ": 'DOUBLE',
        "AU_DIPOL": 'DOUBLE',
        "AU_VELOC": 'DOUBLE',
        "AU_TIME ": 'DOUBLE',
        "EL_GFACT": 'DOUBLE',
        "EA_IRREP": 'INTEGER',
        "UHFRHF  ": 'INTEGER',
        "IFLAGS  ": 'INTEGER',
        "IFLAGS2 ": 'INTEGER',
        "OCCUPYA ": 'INTEGER',
        "NUMDROPA": 'INTEGER',
        "JODAFLAG": 'INTEGER',
        "TITLE   ": 'CHARACTER',
        "NCNSTRNT": 'INTEGER',
        "ICNSTRNT": 'INTEGER',
        "VCNSTRNT": 'DOUBLE',
        "NMPROTON": 'INTEGER',
        "NREALATM": 'INTEGER',
        "COORDINT": 'DOUBLE',
        "VARNAINT": 'DOUBLE',
        "COORD000": 'DOUBLE',
        "ROTCONST": 'DOUBLE',
        "ORIENT2 ": 'DOUBLE',  # input orientation into interial frame
        "LINEAR  ": 'INTEGER',
        "NATOMS  ": 'INTEGER',
        "COORD   ": 'DOUBLE',
        "ORIENTMT": 'DOUBLE',  # input orientation from ZMAT (mostly useful for Cartesians) to Cfour standard orientation
        "ATOMMASS": 'DOUBLE',
        "ORIENT3 ": 'DOUBLE',
        "FULLPTGP": 'CHARACTER',
        "FULLORDR": 'INTEGER',
        "FULLNIRR": 'INTEGER',
        "FULLNORB": 'INTEGER',
        "FULLSYOP": 'DOUBLE',
        "FULLPERM": 'INTEGER',
        "FULLMEMB": 'INTEGER',
        "FULLPOPV": 'INTEGER',
        "FULLCLSS": 'INTEGER',
        "FULLSTGP": 'CHARACTER',
        "ZMAT2MOL": 'INTEGER',
        "COMPPTGP": 'CHARACTER',
        "COMPORDR": 'INTEGER',
        "COMPNIRR": 'INTEGER',
        "COMPNORB": 'INTEGER',
        "COMPSYOP": 'DOUBLE',
        "COMPPERM": 'INTEGER',
        "COMPMEMB": 'INTEGER',
        "COMPPOPV": 'INTEGER',
        "COMPCLSS": 'INTEGER',
        "COMPSTGP": 'CHARACTER',
        "BMATRIX ": 'DOUBLE',
        "NUCREP  ": 'DOUBLE',
        "TIEDCORD": 'INTEGER',
        "MPVMZMAT": 'INTEGER',
        "ATOMCHRG": 'INTEGER',
        "NTOTSHEL": 'INTEGER',
        "NTOTPRIM": 'INTEGER',
        "BASISEXP": 'DOUBLE',
        "BASISCNT": 'DOUBLE',
        "SHELLSIZ": 'INTEGER',
        "SHELLPRM": 'INTEGER',
        "SHELLANG": 'INTEGER',
        "SHELLLOC": 'INTEGER',
        "SHOFFSET": 'INTEGER',
        "SHELLORB": 'INTEGER',
        "PROFFSET": 'INTEGER',
        "PRIMORBT": 'INTEGER',
        "FULSHLNM": 'INTEGER',
        "FULSHLTP": 'INTEGER',
        "FULSHLSZ": 'INTEGER',
        "FULSHLAT": 'INTEGER',
        "JODAOUT ": 'INTEGER',
        "NUMIIII ": 'INTEGER',
        "NUMIJIJ ": 'INTEGER',
        "NUMIIJJ ": 'INTEGER',
        "NUMIJKL ": 'INTEGER',
        "NBASTOT ": 'INTEGER',
        "NAOBASFN": 'INTEGER',
        "NUMBASIR": 'INTEGER',
        "FAOBASIR": 'DOUBLE',
        "AO2SO   ": 'DOUBLE',
        "FULLSOAO": 'DOUBLE',
        "FULLAOSO": 'DOUBLE',
        "AO2SOINV": 'DOUBLE',
        "CART3CMP": 'DOUBLE',
        "CART2CMP": 'DOUBLE',
        "CMP3CART": 'DOUBLE',
        "CMP2CART": 'DOUBLE',
        "ANGMOMBF": 'INTEGER',
        "NBASATOM": 'INTEGER',
        "NAOBFORB": 'INTEGER',
        "MAP2ZMAT": 'INTEGER',
        "CENTERBF": 'INTEGER',
        "CNTERBF0": 'INTEGER',
        "ANMOMBF0": 'INTEGER',
        "CMP2ZMAT": 'DOUBLE',
        "ZMAT2CMP": 'DOUBLE',
        "OVERLAP ": 'DOUBLE',
        "ONEHAMIL": 'DOUBLE',
        "AOOVRLAP": 'DOUBLE',
        "SHALFMAT": 'DOUBLE',
        "SCFEVCA0": 'DOUBLE',
        "RPPBMAT ": 'DOUBLE',
        "OCCUPYA0": 'INTEGER',
        "SYMPOPOA": 'INTEGER',
        "SYMPOPVA": 'INTEGER',
        "SCFEVLA0": 'DOUBLE',
        "SCFDENSA": 'DOUBLE',
        "FOCKA   ": 'DOUBLE',
        "SMHALF  ": 'DOUBLE',
        "EVECOAOA": 'DOUBLE',
        "ONEHMOA ": 'DOUBLE',
        "NOCCORB ": 'INTEGER',
        "NVRTORB ": 'INTEGER',
        "SCFENEG ": 'DOUBLE',
        "TOTENERG": 'DOUBLE',
        "IRREPALP": 'INTEGER',
        "OMEGA_A ": 'DOUBLE',
        "EVECAOXA": 'DOUBLE',
        "EVALORDR": 'DOUBLE',
        "EVECAO_A": 'DOUBLE',
        "EVCSYMAF": 'CHARACTER',
        "EVCSYMAC": 'CHARACTER',
        "TESTVECT": 'DOUBLE',
        "MODROPA ": 'INTEGER',
        "VRHARMON": 'DOUBLE',
        "NEWRECRD": 'INTEGER',
        "VRCORIOL": 'DOUBLE',
        "VRQUADRA": 'DOUBLE',
        "VRANHARM": 'DOUBLE',
        "REFINERT": 'DOUBLE',
        "DIDQ    ": 'DOUBLE',
        "REFCOORD": 'DOUBLE',
        "REFDIPOL": 'DOUBLE',
        "REFGRADI": 'DOUBLE',
        "REFDIPDR": 'DOUBLE',
        "REFNORMC": 'DOUBLE',
        "REFD2EZ ": 'DOUBLE',
        "REFFREQS": 'DOUBLE',
        "REFORIEN": 'DOUBLE',
        "NUSECORD": 'INTEGER',
        "NZMATANH": 'INTEGER',
        "ISELECTQ": 'INTEGER',
        "NEXTGEOM": 'DOUBLE',
        "NEXTGEO1": 'DOUBLE',
        "FCMDISPL": 'DOUBLE',
        "GRDDISPL": 'DOUBLE',
        "DPMDISPL": 'DOUBLE',
        "DIPDISPL": 'DOUBLE',
        "NMRDISPL": 'DOUBLE',
        "SRTDISPL": 'DOUBLE',
        "CHIDISPL": 'DOUBLE',
        "POLDISPL": 'DOUBLE',
        "EFGDISPL": 'DOUBLE',
        "THEDISPL": 'DOUBLE',
        "JFCDISPL": 'DOUBLE',
        "JSDDISPL": 'DOUBLE',
        "JSODISPL": 'DOUBLE',
        "JDSODISP": 'DOUBLE',
        "CUBCOUNT": 'INTEGER',
        "FCMMAPER": 'DOUBLE',
        "QPLSMINS": 'INTEGER',
        "CUBCOORD": 'INTEGER',
        "PASS1   ": 'INTEGER',
        "REFFORDR": 'INTEGER',
        "REFFSYOP": 'DOUBLE',
        "REFFPERM": 'INTEGER',
        "REFNUMIC": 'INTEGER',
        "REFAMAT ": 'DOUBLE',
        "REFTTEN ": 'DOUBLE',
        "REFLINER": 'INTEGER',
        "DIPOLMOM": 'DOUBLE',
        "POLARTEN": 'DOUBLE',
        "CHITENSO": 'DOUBLE',
        "EFGTENSO": 'DOUBLE',
        "IRREPPOP": 'INTEGER',
        "REORDERA": 'INTEGER',
        "IRREPBET": 'INTEGER',
        "SCFEVLB0": 'DOUBLE',
        "SCFEVCB0": 'DOUBLE',
        "IRREPCOU": 'INTEGER',
        "IDROPA  ": 'INTEGER',
        "OCCSCF  ": 'INTEGER',
        "VRTSCF  ": 'INTEGER',
        "SCFEVECA": 'DOUBLE',
        "NCOMPA  ": 'INTEGER',
        "NBASCOMP": 'INTEGER',
        "SCFEVALA": 'DOUBLE',
        "SCFEVALB": 'DOUBLE',
        "SVAVA0  ": 'INTEGER',
        "SVAVA0X ": 'INTEGER',
        "SVAVA0I ": 'INTEGER',
        "SVBVB0  ": 'INTEGER',
        "SVBVB0X ": 'INTEGER',
        "SVBVB0I ": 'INTEGER',
        "SOAOA0  ": 'INTEGER',
        "SOAOA0X ": 'INTEGER',
        "SOAOA0I ": 'INTEGER',
        "SOBOB0  ": 'INTEGER',
        "SOBOB0X ": 'INTEGER',
        "SOBOB0I ": 'INTEGER',
        "SVAVA1  ": 'INTEGER',
        "SVAVA1X ": 'INTEGER',
        "SVAVA1I ": 'INTEGER',
        "SVBVB1  ": 'INTEGER',
        "SVBVB1X ": 'INTEGER',
        "SVBVB1I ": 'INTEGER',
        "SOAOA1  ": 'INTEGER',
        "SOAOA1X ": 'INTEGER',
        "SOAOA1I ": 'INTEGER',
        "SOBOB1  ": 'INTEGER',
        "SOBOB1X ": 'INTEGER',
        "SOBOB1I ": 'INTEGER',
        "SVAOA2  ": 'INTEGER',
        "SVAOA2X ": 'INTEGER',
        "SVAOA2I ": 'INTEGER',
        "SVBOB2  ": 'INTEGER',
        "SVBOB2X ": 'INTEGER',
        "SVBOB2I ": 'INTEGER',
        "SOBVA2  ": 'INTEGER',
        "SOBVA2X ": 'INTEGER',
        "SOBVA2I ": 'INTEGER',
        "SVBOA2  ": 'INTEGER',
        "SVBOA2X ": 'INTEGER',
        "SVBOA2I ": 'INTEGER',
        "SVAVB2  ": 'INTEGER',
        "SVAVB2X ": 'INTEGER',
        "SVAVB2I ": 'INTEGER',
        "SOAOB2  ": 'INTEGER',
        "SOAOB2X ": 'INTEGER',
        "SOAOB2I ": 'INTEGER',
        "SOAVA2  ": 'INTEGER',
        "SOAVA2X ": 'INTEGER',
        "SOAVA2I ": 'INTEGER',
        "SOBVB2  ": 'INTEGER',
        "SOBVB2X ": 'INTEGER',
        "SOBVB2I ": 'INTEGER',
        "SOAVB2  ": 'INTEGER',
        "SOAVB2X ": 'INTEGER',
        "SOAVB2I ": 'INTEGER',
        "SVAVA2  ": 'INTEGER',
        "SVAVA2X ": 'INTEGER',
        "SVAVA2I ": 'INTEGER',
        "SVBVB2  ": 'INTEGER',
        "SVBVB2X ": 'INTEGER',
        "SVBVB2I ": 'INTEGER',
        "SOAOA2  ": 'INTEGER',
        "SOAOA2X ": 'INTEGER',
        "SOAOA2I ": 'INTEGER',
        "SOBOB2  ": 'INTEGER',
        "SOBOB2X ": 'INTEGER',
        "SOBOB2I ": 'INTEGER',
        "SYMPOPOB": 'INTEGER',
        "SYMPOPVB": 'INTEGER',
        "T2NORM  ": 'DOUBLE',
        "MOIOVEC ": 'INTEGER',
        "MOIOWRD ": 'INTEGER',
        "MOIOSIZ ": 'INTEGER',
        "MOIODIS ": 'INTEGER',
        "MOIOFIL ": 'INTEGER',
        "ISYMTYP ": 'INTEGER',
        "TOTRECMO": 'INTEGER',
        "TOTWRDMO": 'INTEGER',
        "RELDENSA": 'DOUBLE',
        "IINTERMA": 'DOUBLE',
        "OCCNUM_A": 'DOUBLE',
        "SCRATCH ": 'DOUBLE',
        "SETUP2  ": 'INTEGER',
        "MOLHES2 ": 'INTEGER',
        "GRAD2   ": 'INTEGER',
        "COORDMAS": 'INTEGER',
        "NUCMULT ": 'INTEGER',
        "SYMCOORD": 'DOUBLE',
        "SYMCOOR2": 'DOUBLE',
        "SYMCOOR3": 'DOUBLE',
        "SYMMLENG": 'INTEGER',
        "SKIP    ": 'INTEGER',
        "NSYMPERT": 'INTEGER',
        "NPERTB  ": 'INTEGER',
        "TRANSINV": 'INTEGER',
        "IBADNUMB": 'INTEGER',
        "IBADINDX": 'INTEGER',
        "IBADIRRP": 'INTEGER',
        "IBADPERT": 'INTEGER',
        "IBADSPIN": 'INTEGER',
        "TREATPER": 'INTEGER',
        "MAXAODSZ": 'INTEGER',
        "PERTINFO": 'INTEGER',
        "GRADIENT": 'DOUBLE',
        "HESSIANM": 'DOUBLE',
        "GRDZORDR": 'DOUBLE',
        "D2EZORDR": 'DOUBLE',
        "REALCORD": 'DOUBLE',
        "DUMSTRIP": 'INTEGER',
        "BMATRIXC": 'DOUBLE',
        "REALATOM": 'INTEGER',
        "NORMCORD": 'DOUBLE',
        "DIPDERIV": 'DOUBLE',
        "I4CDCALC": 'DOUBLE',
        "FREQUENC": 'DOUBLE',
        "RATMMASS": 'DOUBLE',
        "RATMPOSN": 'INTEGER',
        "DEGENERT": 'INTEGER',
        "REFSHILD": 'DOUBLE',
        "CORIZETA": 'DOUBLE',
        "NMPOINTX": 'INTEGER',
        "REFD3EDX": 'DOUBLE',
        "BPPTOB  ": 'DOUBLE',
        "BPTOB   ": 'DOUBLE',
        "BSRTOB  ": 'DOUBLE',
        "BARTOB  ": 'DOUBLE',
        "VRTOTAL ": 'DOUBLE',
        "D2DIPOLE": 'DOUBLE',
        "D3DIPOLE": 'DOUBLE',
        "D1DIPOLE": 'DOUBLE',
        "REFNORM2": 'DOUBLE',
        "NUSECOR2": 'INTEGER',
        "FCMDISP2": 'DOUBLE',
        "RGTDISPL": 'DOUBLE',
        "CUBCOOR1": 'INTEGER',
        "CUBCOOR2": 'INTEGER',
        "REFFPEM2": 'INTEGER',
        "RGTTENSO": 'DOUBLE',
        "REFFPER2": 'INTEGER',
        "REFD4EDX": 'DOUBLE',
        "ZPE_ANHA": 'DOUBLE',
        "OPENSLOT": 'INTEGER',

        "BOLTZMAN": 'DOUBLE',
        "MRCCOCC ": 'INTEGER',
        "ABELPTGP": 'CHARACTER',
        "ABELORDR": 'INTEGER',
        "ABELNIRR": 'INTEGER',
        "ABELNORB": 'INTEGER',
        "ABELSYOP": 'DOUBLE',
        "ABELPERM": 'INTEGER',
        "ABELMEMB": 'INTEGER',
        "ABELPOPV": 'INTEGER',
        "ABELCLSS": 'INTEGER',
        "ABELSTGP": 'CHARACTER',
        "REALCHRG": 'INTEGER',      # atom/mol? charge taking into acct edp
        "NSOSCF  ": 'INTEGER',      # whether is spin orbital calc?
        "SCFVCFLA": 'DOUBLE',       # scf vector expanded from sph to cart basis for symm anal - determin orb sym
        "EFG_SYM1": 'INTEGER',       # symmetry property of components of electric field gradient  integrals
        "EFG_SYM2": 'INTEGER',       # symm prop of comp of EFG

        "DCTDISPL": 'DOUBLE',
        "DANGERUS": 'INTEGER',   #?
        "FULLCHAR": 'CHARACTER', #?
        "FULLDEGN": 'CHARACTER', #?
        "FULLLABL": 'CHARACTER', #?
        "FULLNIRX": 'CHARACTER', #?
        "COMPCHAR": 'CHARACTER', #?
        "COMPDEGN": 'CHARACTER', #?
        "COMPLABL": 'CHARACTER', #?
        "COMPNIRX": 'CHARACTER', #?
        "ROTVECX ": 'CHARACTER', #?
        "ROTVECY ": 'CHARACTER', #?
        "ROTVECZ ": 'CHARACTER', #?
        "COMPNSYQ": 'CHARACTER', #?
        "COMPSYQT": 'CHARACTER', #?
        "COMPSYMQ": 'CHARACTER', #?
        "TRAVECX ": 'CHARACTER', #?
        "TRAVECY ": 'CHARACTER', #?
        "TRAVECZ ": 'CHARACTER', #?
        "NVIBSYM ": 'CHARACTER', #?
        "NUMVIBRT": 'CHARACTER', #?
        "SBGRPSYM": 'CHARACTER', #?
        "ORDERREF": 'CHARACTER', #?
        "OPERSREF": 'CHARACTER', #?
        "NVIBSYMF": 'CHARACTER', #?
        "FULLNSYQ": 'CHARACTER', #?
        "FULLSYQT": 'CHARACTER', #?
        "FULLSYMQ": 'CHARACTER', #?
        "INVPSMAT": 'CHARACTER', #?
        "FDCOORDS": 'CHARACTER', #?
        "FDCALCTP": 'CHARACTER', #?
        "NUMPOINT": 'CHARACTER', #?
        "NPTIRREP": 'CHARACTER', #?
        "GRDPOINT": 'CHARACTER', #?
        "DIPPOINT": 'CHARACTER', #?
        "ENGPOINT": 'CHARACTER', #?
        "PASS1FIN": 'CHARACTER', #?
        "REFENERG": 'CHARACTER', #?
        "NEXTCALC": 'CHARACTER', #?
        "PRINSPIN": 'CHARACTER', #?
        "PRINFROM": 'CHARACTER', #?
        "PRININTO": 'CHARACTER', #?
        "NEXTGEOF": 'CHARACTER', #?
        "ZPE_HARM": 'DOUBLE', #?
        "NDROPPED": 'INTEGER',
        "REFCPTGP": 'INTEGER', #?
        "REFFPTGP": 'INTEGER', #?
        }

    with open('JAINDX', mode='rb') as file:  # b is important -> binary
        fileContent = file.read()
        fileLength = len(fileContent)

    if fileLength == 16012:
        srcints = 4
        srcrecs = 4
    elif fileLength == 16020:
        srcints = 4
        srcrecs = 8
    elif fileLength == 24016:
        srcints = 8
        srcrecs = 4
    elif fileLength == 24024:
        srcints = 8
        srcrecs = 8

    # fixed number of slots for options
    nopt = 1000

    type2len = {
        'DOUBLE': 8,
        'INTEGER': srcints,
        'CHARACTER': 1,
        }

    intlen2format = {
        4: 'i',
        8: 'l',
        }

    type2format = {
        'DOUBLE': 'd',
        'INTEGER': intlen2format[type2len['INTEGER']],
        'CHARACTER': 'c',
        }

    if verbose:
        print('\n<<<  JAINDX  >>>\n')

    posf = srcrecs
    istr = intlen2format[srcrecs]
    jastart = struct.unpack(istr, fileContent[:posf])
    if verbose:
        print('%10s%10d%10d' % ('start', 0, posf))

    poss = posf
    posf = poss + 8 * nopt
    istr = '8s' * nopt
    jaindx = struct.unpack(istr, fileContent[poss:posf])
    if verbose:
        print('%10s%10d%10d' % ('jaindx', poss, posf))

    poss = posf
    posf = poss + srcints * nopt
    istr = intlen2format[srcints] * nopt
    jaindx2 = struct.unpack(istr, fileContent[poss:posf])
    if verbose:
        print('%10s%10d%10d' % ('jaindx2', poss, posf))

    poss = posf
    posf = poss + srcints * nopt
    istr = intlen2format[srcints] * nopt
    jaindx3 = struct.unpack(istr, fileContent[poss:posf])
    if verbose:
        print('%10s%10d%10d' % ('jaindx3', poss, posf))

    poss = posf
    posf = poss + srcints
    istr = intlen2format[srcints]
    jamid = struct.unpack(istr, fileContent[poss:posf])
    if verbose:
        print('%10s%10d%10d' % ('mid', poss, posf))

    poss = posf
    posf = poss + srcrecs
    istr = intlen2format[srcrecs]
    jaend = struct.unpack(istr, fileContent[poss:posf])
    if verbose:
        print('%10s%10d%10d' % ('end', poss, posf))

    nrecs = jaindx.index('OPENSLOT')  # number of active records

    if verbose:
        print('\n')
        print('%20s%10d' % ('File Length:', fileLength))
        print('%20s%10d' % ('srcints Int Length:', srcints))
        print('%20s%10d' % ('srcrecs Int Length:', srcrecs))
        print('%20s%10d' % ('First Rec:', jastart[0]))
        print('%20s%10d' % ('Second Rec:', jamid[0]))
        print('%20s%10d' % ('Last Rec:', jaend[0]))
        print('%20s%10d' % ('Full Records:', nrecs))
        print('\n')

        print('\n<<<  JOBARC  >>>\n')

    with open('JOBARC', mode='rb') as file:  # b is important -> binary
        fileContent = file.read()

    returnRecords = {}
    poss = 0
    for item in range(nrecs):
        posf = poss + type2len[knownlabels[jaindx[item]]] * jaindx3[item]
        istr = type2format[knownlabels[jaindx[item]]] * jaindx3[item]
        if knownlabels[jaindx[item]] == 'CHARACTER':
            bound = type2len[knownlabels[jaindx[item]]] * jaindx3[item] * 8
            posf = poss + bound
            istr = str(bound) + 's'
        jobarc = struct.unpack(istr, fileContent[poss:posf])

        if verbose:
            #print item, istr, poss, posf, '\t', jaindx[item], jaindx2[item], jaindx3[item], jobarc
            if jaindx3[item] < 120:
                print(jaindx[item], jaindx2[item], jaindx3[item], jobarc)

        poss = posf
        if jaindx[item] in reclabelarray:
            returnRecords[jaindx[item]] = jobarc

    return returnRecords

#if __name__ == "__main__":
#    want = ['NATOMS  ', 'AU_LENGT', 'COORD   ', 'HBAR    ', 'ATOMCHRG']
##    got = get_jajo_record(want)
#    got = getrec(want)
#    for item in got.keys():
#        print item, got[item]
