import re
import logging
from decimal import Decimal

import numpy as np
import qcelemental as qcel
from qcelemental.models import Molecule
from qcelemental.molparse import regex

from ..util import PreservingDict, load_hessian

logger = logging.getLogger(__name__)


def harvest_output(outtext):
    """Function to separate portions of a CFOUR output file *outtest*,
    divided by xjoda.

    """
    pass_psivar = []
    pass_coord = []
    pass_grad = []

    # for outpass in re.split(r'--invoking executable xjoda', outtext, re.MULTILINE):
    for outpass in re.split(r"JODA beginning optimization cycle", outtext, re.MULTILINE):
        psivar, c4coord, c4grad, version, module, error = harvest_outfile_pass(outpass)
        pass_psivar.append(psivar)
        pass_coord.append(c4coord)
        pass_grad.append(c4grad)

        # print('\n\nXXXXXXXXXXXXXXXXXXXXXXXXXX\n\n')
        # print(outpass)
        # print(psivar, c4coord, c4grad, version, error)
        # print('\n\nxxxxxxxxxxxxxxxxxxxxxxxxxx\n\n')

    retindx = -1 if pass_coord[-1] else -2

    #    print '    <<<  C4 PSIVAR  >>>'
    #    for item in pass_psivar[retindx]:
    #        print('       %30s %16.8f' % (item, pass_psivar[retindx][item]))
    #    print '    <<<  C4 COORD   >>>'
    #    for item in pass_coord[retindx]:
    #        print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))
    #    print '    <<<   C4 GRAD   >>>'
    #    for item in pass_grad[retindx]:
    #        print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))

    return pass_psivar[retindx], pass_coord[retindx], pass_grad[retindx], version, module, error


def harvest_outfile_pass(outtext):
    """Function to read CFOUR output file *outtext* and parse important
    quantum chemical information from it in

    """
    psivar = PreservingDict()
    psivar_coord = None
    psivar_grad = None
    version = ""
    module = None
    error = ""

    #    TODO: BCC
    #          CI
    #          other ROHF tests
    #          vcc/ecc

    NUMBER = r"(?x:" + regex.NUMBER + ")"
    DECIMAL = r"(?x:" + regex.DECIMAL + ")"

    # Process version
    mobj = re.search(r"^\s*" + r"Version" + r"\s+" + r"(?P<version>[\w.]+)" + r"\s*$", outtext, re.MULTILINE)
    if mobj:
        print("matched version")
        version = mobj.group("version")

    # Process NRE
    mobj = re.search(
        r"^\s+" + r"(?:Nuclear repulsion energy :)" + r"\s+" + NUMBER + r"\s+a\.u\.\s*$",
        outtext,
        re.MULTILINE | re.IGNORECASE,
    )
    if mobj:
        print("matched nre")
        psivar["NUCLEAR REPULSION ENERGY"] = mobj.group(1)

    # Process calcinfo
    mobj = re.search(
        r"^\s*" + r"There are" + r"\s+" + r"(?P<nbf>\d+)" + r"\s+" + r"functions in the AO basis." + r"\s*$",
        outtext,
        re.MULTILINE,
    )
    if mobj:
        print("matched nbf", mobj.groups())
        psivar["N BASIS FUNCTIONS"] = mobj.group("nbf")
        psivar["N MOLECULAR ORBITALS"] = mobj.group("nbf")  # TODO BAD

    mobj = re.search(
        # fmt: off
        r"^\s*" + "Alpha population by irrep:" + r"(?P<aocc>[\d\s]+)" + r"\s*" +
        r"^\s*" + "Beta population by irrep:" + r"(?P<bocc>[\d\s]+)" + r"\s*",
        # fmt: on
        outtext,
        re.MULTILINE,
    )
    if mobj:
        print("matched occupied", mobj.groups())
        psivar["N ALPHA ELECTRONS"] = sum([int(d) for d in mobj.group("aocc").split()])
        psivar["N BETA ELECTRONS"] = sum([int(d) for d in mobj.group("bocc").split()])

    # Process SCF
    mobj = re.search(r"^\s+" + r"(?:E\(SCF\))" + r"\s+=\s+" + NUMBER + r"\s+a\.u\.\s*$", outtext, re.MULTILINE)
    if mobj:
        print("matched scf1")
        psivar["SCF TOTAL ENERGY"] = mobj.group(1)

    mobj = re.search(r"^\s+" + r"(?:E\(SCF\)=)" + r"\s+" + NUMBER + r"\s+" + NUMBER + r"\s*$", outtext, re.MULTILINE)
    if mobj:
        print("matched scf2")
        psivar["SCF TOTAL ENERGY"] = mobj.group(1)

    if "SCF TOTAL ENERGY" not in psivar:
        # can be too greedy and match across scf cycles
        mobj = re.search(
            # fmt: off
            r'^\s+' + r'(?:SCF has converged.)' + r'\s*$' +
            r'(?:.*?)' +
            r'^\s+' + r'(?:\d+)' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$',
            # fmt: on
            outtext,
            re.MULTILINE | re.DOTALL,
        )
        if mobj:
            print("matched scf3")
            psivar["SCF TOTAL ENERGY"] = mobj.group(1)

    mobj = re.search(r"^\s+" + r"(?:E\(ROHF\)=)" + r"\s+" + NUMBER + r"\s+" + NUMBER + r"\s*$", outtext, re.MULTILINE)
    if mobj:
        print("matched scf4")
        psivar["SCF TOTAL ENERGY"] = mobj.group(1)

    # Process MP2
    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?:E2\(AA\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:E2\(AB\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:E2\(TOT\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:Total MP2 energy)' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*$',
        # fmt: on
        outtext,
        re.MULTILINE,
    )
    if mobj:
        logger.debug("matched mp2r")
        psivar["MP2 SAME-SPIN CORRELATION ENERGY"] = 2 * Decimal(mobj.group(1))
        psivar["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = mobj.group(2)
        psivar["MP2 CORRELATION ENERGY"] = 2 * Decimal(mobj.group(1)) + Decimal(mobj.group(2))
        psivar["MP2 TOTAL ENERGY"] = mobj.group(4)

    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?:E2\(AA\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:E2\(BB\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:E2\(AB\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:E2\(TOT\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:Total MP2 energy)' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*$',
        # fmt: on
        outtext,
        re.MULTILINE,
    )
    if mobj:
        print("matched mp2u")
        psivar["MP2 SAME-SPIN CORRELATION ENERGY"] = Decimal(mobj.group(1)) + Decimal(mobj.group(2))
        psivar["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = mobj.group(3)
        psivar["MP2 CORRELATION ENERGY"] = Decimal(mobj.group(1)) + Decimal(mobj.group(2)) + Decimal(mobj.group(3))
        psivar["MP2 TOTAL ENERGY"] = mobj.group(5)

    mobj = re.search(
        # particularly, want to avoid capture when following line present:
        #  "MP2 energies are correct only for semicanonical orbitals."
        # fmt: off
        r'Singles contribution will be calculated.' + r'\s*' +
        r'^\s+' + r'-*' + r'\s*' +
        r'^\s+' + r'(?:E\(SCF\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:E2\(AA\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:E2\(BB\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:E2\(AB\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:E2\(SINGLE\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:E2\(TOT\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:Total MP2 energy)' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*$',
        # fmt: on
        outtext,
        re.MULTILINE,
    )
    if mobj:
        print("matched mp2ro")
        psivar["MP2 SAME-SPIN CORRELATION ENERGY"] = Decimal(mobj.group(2)) + Decimal(mobj.group(3))
        psivar["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = mobj.group(4)
        psivar["MP2 SINGLES ENERGY"] = mobj.group(5)
        psivar["MP2 CORRELATION ENERGY"] = (
            Decimal(mobj.group(5)) + Decimal(mobj.group(2)) + Decimal(mobj.group(3)) + Decimal(mobj.group(4))
        )
        psivar["MP2 TOTAL ENERGY"] = mobj.group(7)

    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?:S-MBPT\(2\))' + r'\s+' + r'(?P<sgl>' + NUMBER + r')' + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:D-MBPT\(2\))' + r'\s+' + r'(?P<dbl>' + NUMBER + r')' + r'\s+' +
                                                r'(?P<mp2tot>' + NUMBER + r')' + r'\s*$',
        # fmt: on
        outtext,
        re.MULTILINE,
    )
    if mobj:
        print("matched mp2ro2")
        # psivar['MP2 SAME-SPIN CORRELATION ENERGY'] = Decimal(mobj.group(1)) + Decimal(mobj.group(2))
        # psivar['MP2 OPPOSITE-SPIN CORRELATION ENERGY'] = mobj.group(3)
        psivar["MP2 SINGLES ENERGY"] = mobj.group("sgl")
        psivar["MP2 CORRELATION ENERGY"] = Decimal(mobj.group("sgl")) + Decimal(mobj.group("dbl"))
        psivar["MP2 TOTAL ENERGY"] = mobj.group("mp2tot")

    # Process MP3
    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?:D-MBPT\(2\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:D-MBPT\(3\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        print("matched mp3r")
        dmp2 = Decimal(mobj.group(1))
        dmp3 = Decimal(mobj.group(3))
        psivar["MP2 CORRELATION ENERGY"] = dmp2
        psivar["MP2 TOTAL ENERGY"] = mobj.group(2)
        psivar["MP3 CORRELATION ENERGY"] = dmp2 + dmp3
        psivar["MP3 TOTAL ENERGY"] = mobj.group(4)
        psivar["MP2.5 CORRELATION ENERGY"] = dmp2 + Decimal("0.500000000000") * dmp3
        psivar["MP2.5 TOTAL ENERGY"] = psivar["MP2.5 CORRELATION ENERGY"] + psivar["SCF TOTAL ENERGY"]
        psivar["MP3 SINGLES ENERGY"] = Decimal("0.0")
        psivar["MP3 DOUBLES ENERGY"] = dmp2 + dmp3
        module = "vcc"

    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?:S-MBPT\(2\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:D-MBPT\(2\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:S-MBPT\(3\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:D-MBPT\(3\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        print("matched mp3ro")
        dmp2 = Decimal(mobj.group(1)) + Decimal(mobj.group(3))
        dmp3 = Decimal(mobj.group(5)) + Decimal(mobj.group(7))
        psivar["MP3 CORRELATION ENERGY"] = dmp2 + dmp3
        psivar["MP3 TOTAL ENERGY"] = mobj.group(8)
        psivar["MP2.5 CORRELATION ENERGY"] = dmp2 + Decimal("0.500000000000") * dmp3
        psivar["MP2.5 TOTAL ENERGY"] = psivar["MP2.5 CORRELATION ENERGY"] + psivar["SCF TOTAL ENERGY"]
        psivar["MP3 SINGLES ENERGY"] = Decimal(mobj.group(1)) + Decimal(mobj.group(5))
        psivar["MP3 DOUBLES ENERGY"] = Decimal(mobj.group(3)) + Decimal(mobj.group(7))
        module = "vcc"

    mobj = re.search(
        # fmt: off
        r"^\s*" + r"(?:MP2 correlation energy:)\s+" + r"(?P<mp2corl>" + NUMBER + ")" + r"\s*" +
        r"^\s*" + r"(?:MP3 correction:)\s+" +         r"(?P<mp3corr>" + NUMBER + ")" + r"\s*" +
        r"^\s*" + r"(?:MP3 correlation energy:)\s+" + r"(?P<mp3corl>" + NUMBER + ")" + r"\s*" +
        r"(?:.*?)" +
        r"^\s*" + r"(?:Non-iterative calculation of MP3)" + r".*" +
        r"(?:.*?)" +
        r"^\s*" + r"(?:Total MP3 energy:)" + r"\s+" + r"(?P<mp3tot>" + NUMBER + ")" + r"\s*$",
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        logger.debug("matched mp3 ncc")
        # psivar["MP2 CORRELATION ENERGY"] = mobj.group("mp2corl")
        psivar["MP3 CORRELATION ENERGY"] = mobj.group("mp3corl")
        psivar["MP3 CORRECTION ENERGY"] = mobj.group("mp3corr")
        psivar["MP3 TOTAL ENERGY"] = mobj.group("mp3tot")
        # looks like ncc is rhf-only
        # psivar["MP2 DOUBLES ENERGY"] = mobj.group("mp2corl")
        psivar["MP3 DOUBLES ENERGY"] = mobj.group("mp3corl")
        module = "ncc"

    # Process MP4
    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?:D-MBPT\(2\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:D-MBPT\(3\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:D-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:Q-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:S-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        print("matched mp4r")
        dmp2 = Decimal(mobj.group(1))
        dmp3 = Decimal(mobj.group(3))
        dmp4sdq = Decimal(mobj.group(5)) + Decimal(mobj.group(7)) + Decimal(mobj.group(9))
        psivar["MP2 CORRELATION ENERGY"] = dmp2
        psivar["MP2 TOTAL ENERGY"] = mobj.group(2)
        psivar["MP3 CORRELATION ENERGY"] = dmp2 + dmp3
        psivar["MP3 TOTAL ENERGY"] = mobj.group(4)
        psivar["MP2.5 CORRELATION ENERGY"] = dmp2 + Decimal("0.500000000000") * dmp3
        psivar["MP2.5 TOTAL ENERGY"] = psivar["MP2.5 CORRELATION ENERGY"] + psivar["SCF TOTAL ENERGY"]
        psivar["MP4(SDQ) CORRELATION ENERGY"] = dmp2 + dmp3 + dmp4sdq
        psivar["MP4(SDQ) TOTAL ENERGY"] = mobj.group(10)

    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?:S-MBPT\(2\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:D-MBPT\(2\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:S-MBPT\(3\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:D-MBPT\(3\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:L-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:NL-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        print("matched mp4ro")
        module = "vcc"
        dmp2 = Decimal(mobj.group(1)) + Decimal(mobj.group(3))
        dmp3 = Decimal(mobj.group(5)) + Decimal(mobj.group(7))
        dmp4sdq = Decimal(mobj.group(9)) + Decimal(mobj.group(11))
        psivar["MP2 CORRELATION ENERGY"] = dmp2
        psivar["MP2 TOTAL ENERGY"] = mobj.group(4)
        psivar["MP3 CORRELATION ENERGY"] = dmp2 + dmp3
        psivar["MP3 TOTAL ENERGY"] = mobj.group(8)
        psivar["MP2.5 CORRELATION ENERGY"] = dmp2 + Decimal("0.500000000000") * dmp3
        psivar["MP2.5 TOTAL ENERGY"] = psivar["MP2.5 CORRELATION ENERGY"] + psivar["SCF TOTAL ENERGY"]
        psivar["MP4(SDQ) CORRELATION ENERGY"] = dmp2 + dmp3 + dmp4sdq
        psivar["MP4(SDQ) TOTAL ENERGY"] = mobj.group(12)

    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?:D-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:Q-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:S-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:T-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        print("matched mp4tr")
        dmp4sdq = Decimal(mobj.group(1)) + Decimal(mobj.group(3)) + Decimal(mobj.group(5))
        dmp4t = Decimal(mobj.group(7))
        psivar["MP4(SDQ) CORRELATION ENERGY"] = psivar["MP3 CORRELATION ENERGY"] + dmp4sdq
        psivar["MP4(SDQ) TOTAL ENERGY"] = mobj.group(6)
        psivar["MP4(T) CORRECTION ENERGY"] = dmp4t
        psivar["MP4 CORRELATION ENERGY"] = psivar["MP3 CORRELATION ENERGY"] + dmp4sdq + dmp4t
        psivar["MP4 TOTAL ENERGY"] = mobj.group(8)

    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?:L-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:NL-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:WT12-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:T-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        print("matched mp4tro")
        module = "vcc"
        dmp4sdq = Decimal(mobj.group(1)) + Decimal(mobj.group(3))
        dmp4t = Decimal(mobj.group(5)) + Decimal(mobj.group(7))  # WT12 with T, not SDQ
        psivar["MP4(SDQ) CORRELATION ENERGY"] = psivar["MP3 CORRELATION ENERGY"] + dmp4sdq
        psivar["MP4(SDQ) TOTAL ENERGY"] = mobj.group(4)
        psivar["MP4(T) CORRECTION ENERGY"] = dmp4t
        psivar["MP4 CORRELATION ENERGY"] = psivar["MP3 CORRELATION ENERGY"] + dmp4sdq + dmp4t
        psivar["MP4 TOTAL ENERGY"] = mobj.group(8)

    mobj = re.search(
        # fmt: off
        r"^\s*" + r"(?:MP2 correlation energy:)\s+" + r"(?P<mp2corl>" + NUMBER + ")" + r"\s*" +
        r"^\s*" + r"(?:MP3 correction:)\s+" +         r"(?P<mp3corr>" + NUMBER + ")" + r"\s*" +
        r"^\s*" + r"(?:MP3 correlation energy:)\s+" + r"(?P<mp3corl>" + NUMBER + ")" + r"\s*" +
        r"^\s*" + r"(?:SDQ-MP4 correction:)\s+" + r"(?P<mp4sdqcorr>" + NUMBER + ")" + r"\s*" +
        r"^\s*" + r"(?:SDQ-MP4 correlation energy:)\s+" + r"(?P<mp4sdqcorl>" + NUMBER + ")" + r"\s*" +
        r"(" +
        r"^\s*" + r"(?:T-MP4 correction:)\s+" + r"(?P<mp4tcorr>" + NUMBER + ")" + r"\s*" +
        r"^\s*" + r"(?:Total MP4 correction:)\s+" + r"(?P<mp4corr>" + NUMBER + ")" + r"\s*" +
        r"^\s*" + r"(?:MP4 correlation energy:)\s+" + r"(?P<mp4sdtqcorl>" + NUMBER + ")" + r"\s*" +
        r")?" +
        r"(?:.*?)" +
        r"^\s*" + r"(?:Non-iterative calculation of (MP4|SDQ-MP4))" + r".*" +
        r"(?:.*?)" +
        r"^\s*" + r"(?:Total (?P<mp4flavor>(MP4|SDQ-MP4)) energy:)" + r"\s+" + r"(?P<mp4flavortot>" + NUMBER + ")" + r"\s*$",
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        logger.debug("matched mp4 ncc")
        # psivar["MP2 CORRELATION ENERGY"] = mobj.group("mp2corl")
        module = "ncc"
        mtd = {"MP4": "MP4", "SDQ-MP4": "MP4(SDQ)"}[mobj.group("mp4flavor")]
        psivar["MP3 CORRELATION ENERGY"] = mobj.group("mp3corl")
        psivar["MP3 CORRECTION ENERGY"] = mobj.group("mp3corr")
        psivar["MP4(SDQ) CORRELATION ENERGY"] = mobj.group("mp4sdqcorl")
        # looks like ncc is rhf-only
        # psivar["MP2 DOUBLES ENERGY"] = mobj.group("mp2corl")
        psivar["MP3 DOUBLES ENERGY"] = mobj.group("mp3corl")
        if mtd == "MP4(SDQ)":
            psivar["MP4(SDQ) TOTAL ENERGY"] = mobj.group("mp4flavortot")
        elif mtd == "MP4":
            psivar["MP4(T) CORRECTION ENERGY"] = mobj.group("mp4tcorr")
            psivar["MP4 CORRECTION ENERGY"] = mobj.group("mp4corr")
            psivar["MP4 TOTAL ENERGY"] = mobj.group("mp4flavortot")
            psivar["MP4 CORRELATION ENERGY"] = mobj.group("mp4sdtqcorl")

    # Process CI Iterations
    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?P<fullCI>(?P<iterCI>Q?CI(?:\w+))(?:\(T\))?)' + r'\s+(?:energy will be calculated.)\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'(?:\d+)' + r'\s+' + r"(?P<corl>" + NUMBER + r")" + r'\s+' + r"(?P<tot>" + NUMBER + r")" + r'\s+DIIS\s*' +
        r'^\s*(?:-+)\s*' +
        # CI iterations for CISD, CC iterations for QCISD
        r'^\s*(?:A miracle (?P<ccprog>has come|come) to pass. The (CI|CC) iterations have converged.)\s*$',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        logger.debug("matched ci with full %s iterating %s" % (mobj.group("fullCI"), mobj.group("iterCI")))
        module = {"has come": "vcc", "come": "ecc"}[mobj.group("ccprog")]

        mtd = mobj.group("iterCI").upper()
        psivar[f"{mtd} CORRELATION ENERGY"] = mobj.group("corl")
        psivar[f"{mtd} TOTAL ENERGY"] = mobj.group("tot")

        mobj2 = re.search(
            # fmt: off
            r"^\s+" + r"E\(QCISD\)\s+=\s+" + r"(?P<qcisd>" + NUMBER + r")" + r"\s*" +
            r"^\s+" + r"E\(QCISD\(T\)\)\s+=\s+" + r"(?P<qcisdt>" + NUMBER + r")" + r"\s*$",
            # fmt: on
            outtext,
            re.MULTILINE | re.DOTALL,
        )
        if mobj2 and mobj.group("fullCI") == "QCISD(T)":
            logger.debug("matched qcisd(t)")
            psivar["QCISD(T) TOTAL ENERGY"] = mobj2.group("qcisdt")
            psivar["QCISD(T) CORRECTION ENERGY"] = Decimal(mobj2.group("qcisdt")) - Decimal(mobj2.group("qcisd"))
            psivar["QCISD(T) CORRELATION ENERGY"] = psivar["QCISD(T) TOTAL ENERGY"] - psivar["SCF TOTAL ENERGY"]

    # Process CC Iterations
    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?P<fullCC>(?P<iterCC>L?CC(?:\w+(?:-(?:1|1b|2|3))?))(?:\(T\))?)' + r'\s+(?:energy will be calculated.)\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'(?:\d+)' + r'\s+' + r"(?P<corl>" + NUMBER + r")" + r'\s+' + r"(?P<tot>" + NUMBER + r")" + r'\s+DIIS\s*' +
        r'^\s*(?:-+)\s*' +
        r'^\s*(?:A miracle (?P<ccprog>has come|come) to pass. The CC iterations have converged.)\s*$',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        logger.debug("matched cc with full %s iterating %s" % (mobj.group("fullCC"), mobj.group("iterCC")))
        module = {"has come": "vcc", "come": "ecc"}[mobj.group("ccprog")]

        mobj4 = re.search(r"CALCLEVEL\s+ICLLVL\s+CCSDT-1b", outtext)
        mtd = mobj.group("iterCC").upper()
        if mtd == "CCSDT-1":
            if mobj4 and module == "vcc":
                mtd = "CCSDT-1B"
            else:
                mtd = "CCSDT-1A"
        elif mtd == "CCSDT-1b":
            mtd = "CCSDT-1B"
        psivar[f"{mtd} CORRELATION ENERGY"] = mobj.group("corl")
        psivar[f"{mtd} TOTAL ENERGY"] = mobj.group("tot")

        mobj3 = re.search(r"SCF reference function:  RHF", outtext)
        if mobj3 and mtd not in ["CCSDT-1A", "CCSDT-1B", "CCSDT-2", "CCSDT-3", "CCSDT"]:
            psivar[f"{mtd} DOUBLES ENERGY"] = mobj.group("corl")

    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?:\d+)' + r'\s+' + r'(?P<corl>' + NUMBER + r')\s+' +
                  NUMBER + r'\s+' + NUMBER + r'\s+' +
                  NUMBER + r'\s+' + NUMBER +  r'\s+' +
                  r'(' + NUMBER +  r')?' + r'(' + r'\s+' + NUMBER + r')?' + r'\s*' +
        r'^\s*' +
        r'^\s*' + r'(?:\w+(?:-(1a|1b|2|3))? iterations converged .*?)' +
        r'^\s*' +
        r'^\s*' + r'(?:Total (?P<iterCC>\w+(?:-(1a|1b|2|3))?) energy:)' + r'\s+' + r'(?P<tot>' + NUMBER + r')\s*$',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        print("matched ncc cc iter")
        # looks like ncc is rhf-only
        mtd = mobj.group("iterCC").upper()
        psivar[f"{mtd} CORRELATION ENERGY"] = mobj.group("corl")
        if mtd not in ["CCSDT-1A", "CCSDT-1B", "CCSDT-2", "CCSDT-3", "CCSDT"]:
            psivar[f"{mtd} DOUBLES ENERGY"] = mobj.group("corl")
        psivar[f"{mtd} TOTAL ENERGY"] = mobj.group("tot")
        module = "ncc"

    mobj = re.search(
        # fmt: off
        r'^\s+' + r"Beginning iterative solution of (?P<iterCC>\w+(?:-\d)?) equations" + r"\s*" +
        r'(?P<iterations>.*)' +
        r'^\s*' + r"It\." + r"\s+" + "Correlation Energy" + r".*" +
        r'^\s*(?:-+)\s*' +
        r'^\s*' +
        r'^\s*' + r'(?:\w+(?:-\d)? iterations converged .*?)' +
        r'^\s*' +
        r'^\s*' + r'(?:Total \1 energy:)' + r'\s+' + r'(?P<tot>' + NUMBER + r')\s*$',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        mobj2 = re.findall(
            # fmt: off
            r"(\d+)" + r"\s+" + r"(?P<corl>" + DECIMAL + r")\s+" +
                    DECIMAL + r"\s+" + DECIMAL + r"\s+" + DECIMAL + r"\s+" + DECIMAL + r"\s+" + r"(" + DECIMAL + r")?" + r"(" + r"\s+" + DECIMAL + r")?",
            # fmt: on
            mobj.group("iterations"),
        )
        if mobj2:
            print("matched ncc cc iter mod5", mobj.groupdict(), mobj2[-1])
            mtd = mobj.group("iterCC").upper()
            psivar[f"{mtd} CORRELATION ENERGY"] = mobj2[-1][2]
            if mtd not in ["CCSDT-1A", "CCSDT-1B", "CCSDT-2", "CCSDT-3", "CCSDT"]:
                psivar[f"{mtd} DOUBLES ENERGY"] = mobj2[-1][2]
            psivar[f"{mtd} TOTAL ENERGY"] = mobj.group("tot")
            module = "ncc"

    # Process CC(T)
    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?:E\(SCF\))' + r'\s+=\s+' + NUMBER + r'\s+a\.u\.\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'(?:E\(CCSD\))' + r'\s+=\s+' + NUMBER + r'\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'(?:E\(CCSD\(T\)\))' + r'\s+=\s+' + NUMBER + r'\s*$',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        print("matched ccsd(t) vcc")
        psivar["SCF TOTAL ENERGY"] = mobj.group(1)
        psivar["CCSD TOTAL ENERGY"] = mobj.group(2)
        psivar["(T) CORRECTION ENERGY"] = Decimal(mobj.group(3)) - Decimal(mobj.group(2))
        psivar["CCSD(T) CORRELATION ENERGY"] = Decimal(mobj.group(3)) - Decimal(mobj.group(1))
        psivar["CCSD(T) TOTAL ENERGY"] = mobj.group(3)
        module = "vcc"

    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?:E\(CCSD\))' + r'\s+=\s+' + NUMBER + r'\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'(?:E\(CCSD\(T\)\))' + r'\s+=\s+' + NUMBER + r'\s*$',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        print("matched ccsd(t) vcc v2")
        psivar["CCSD TOTAL ENERGY"] = mobj.group(1)
        psivar["(T) CORRECTION ENERGY"] = Decimal(mobj.group(2)) - Decimal(mobj.group(1))
        psivar["CCSD(T) TOTAL ENERGY"] = mobj.group(2)
        module = "vcc"

    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?:E\(SCF\))' + r'\s+=\s*' + NUMBER + r'\s+a\.u\.\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'(?:CCSD energy)' + r'\s+' + NUMBER + r'\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'(?:Total perturbative triples energy:)' + r'\s+' + NUMBER + r'\s*' +
        r'^\s*(?:-+)\s*' +
        r'^\s+' + r'(?:CCSD\(T\) energy)' + r'\s+' + NUMBER + r'\s*$',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        print("matched ccsd(t) ecc")
        psivar["SCF TOTAL ENERGY"] = mobj.group(1)
        psivar["CCSD TOTAL ENERGY"] = mobj.group(2)
        psivar["(T) CORRECTION ENERGY"] = mobj.group(3)
        psivar["CCSD(T) CORRELATION ENERGY"] = Decimal(mobj.group(4)) - Decimal(mobj.group(1))
        psivar["CCSD(T) TOTAL ENERGY"] = mobj.group(4)
        module = "ecc"

    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?:HF-SCF energy)' + r'\s+' + r"(?P<hf>" + NUMBER + r")" + r'\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'(?:CCSD energy)' + r'\s+' + r"(?P<ccsd>" + NUMBER + r")" + r'\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'(?:E4T  to CCSD\(T\))' + r'\s+' + r"(?P<e4t>" + NUMBER + r")" + r'\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'(?:E4T \+ E5ST)' + r'\s+' + r"(?P<e4te5st>" + NUMBER + r")" + r'\s*' +
        r'(?:.*?)' +
        r'^\s*(?:-+)\s*' +
        r'^\s+' + r'(?:CCSD\(T\) energy)' + r'\s+' + r"(?P<ccsd_t_>" + NUMBER + r")" + r'\s*$',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        print("matched ccsd(t) ecc v2")
        psivar["SCF TOTAL ENERGY"] = mobj.group("hf")
        psivar["CCSD TOTAL ENERGY"] = mobj.group("ccsd")
        psivar["T(CCSD) CORRECTION ENERGY"] = mobj.group("e4t")
        psivar["(T) CORRECTION ENERGY"] = mobj.group("e4te5st")
        psivar["CCSD(T) CORRELATION ENERGY"] = Decimal(mobj.group("ccsd_t_")) - Decimal(mobj.group("hf"))
        psivar["CCSD(T) TOTAL ENERGY"] = mobj.group("ccsd_t_")
        psivar["CCSD+T(CCSD) CORRELATION ENERGY"] = (
            psivar["CCSD CORRELATION ENERGY"] + psivar["T(CCSD) CORRECTION ENERGY"]
        )
        psivar["CCSD+T(CCSD) TOTAL ENERGY"] = psivar["CCSD TOTAL ENERGY"] + psivar["T(CCSD) CORRECTION ENERGY"]
        module = "ecc"

    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?:CCSD energy)' + r'\s+' + NUMBER + r'\s*' +
        r'^\s*(?:-+)\s*' +
        r'^\s+' + r'(?:CCSD\(T\) energy)' + r'\s+' + NUMBER + r'\s*$',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        print("matched ccsd(t) lamb")
        psivar["CCSD TOTAL ENERGY"] = mobj.group(1)
        psivar["(T) CORRECTION ENERGY"] = Decimal(mobj.group(2)) - Decimal(mobj.group(1))
        psivar["CCSD(T) CORRELATION ENERGY"] = Decimal(mobj.group(2)) - psivar["SCF TOTAL ENERGY"]
        psivar["CCSD(T) TOTAL ENERGY"] = mobj.group(2)

    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?:CCSD\(T\) contribution:)\s+' + r'(?P<tcorr>' + NUMBER + ')' + r'\s*'
        r'^\s*' + r'(?:CCSD\[T\] contribution:)\s+' + r'(?P<bkttcorr>' + NUMBER + ')' + r'\s*'
        r'^\s*' + r'(?:Total CCSD\(T\) energy:)\s+' + r'(?P<ttot>' + NUMBER + ')' + r'\s*$',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        print("matched ccsd(t) ncc")
        psivar["(T) CORRECTION ENERGY"] = mobj.group("tcorr")
        psivar["T(CCSD) CORRECTION ENERGY"] = mobj.group("bkttcorr")
        psivar["CCSD(T) TOTAL ENERGY"] = mobj.group("ttot")
        module = "ncc"

    mobj = re.search(
        # fmt: off
        r'^\s*' + r'(?:CCSD\[T\] correlation energy:)\s+' + r'(?P<bkttcorr>' + NUMBER + ')' + r'\s*'
        r'^\s*' + r'(?:CCSD\(T\) correlation energy:)\s+' + r'(?P<tcorr>' + NUMBER + ')' + r'\s*$',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        logger.debug("matched ccsd(t) ncc v2")
        psivar["(T) CORRECTION ENERGY"] = mobj.group("tcorr")
        psivar["T(CCSD) CORRECTION ENERGY"] = mobj.group("bkttcorr")
        psivar["CCSD+T(CCSD) TOTAL ENERGY"] = psivar["T(CCSD) CORRECTION ENERGY"] + psivar["CCSD TOTAL ENERGY"]
        psivar["CCSD+T(CCSD) CORRELATION ENERGY"] = (
            psivar["T(CCSD) CORRECTION ENERGY"] + psivar["CCSD CORRELATION ENERGY"]
        )
        psivar["CCSD(T) TOTAL ENERGY"] = psivar["(T) CORRECTION ENERGY"] + psivar["CCSD TOTAL ENERGY"]
        psivar["CCSD(T) CORRELATION ENERGY"] = psivar["(T) CORRECTION ENERGY"] + psivar["CCSD CORRELATION ENERGY"]
        module = "ncc"

    mobj = re.search(
        # fmt: off
        r"^\s*" + r"(?:Lambda-CCSD iterations converged .*)" +
        r'(?:.*?)' +
        r"^\s*" + r"(?:CCSD\[T\] correlation energy:)\s+" + r"(?P<bkttcorr>" + NUMBER + r")" + r"\s*" +
        r"^\s*" + r"(?:CCSD\(T\) correlation energy:)\s+" + r"(?P<tcorr>" + NUMBER + r")" + r"\s*" +
        r"^\s*" + r"(?:CCSD\[T\]_L correlation energy:)\s+" + r"(?P<abkttcorr>" + NUMBER + r")" + r"\s*" +
        r"^\s*" + r"(?:CCSD\(T\)_L correlation energy:)\s+" + r"(?P<atcorr>" + NUMBER + r")" + r"\s*" +
        r'(?:.*?)' +
        r"^\s*" + r"(?:Non-iterative calculation of CCSD\(T\)_L .*)" +
        r'(?:.*?)' +
        r"^\s*" + r"(?:Total CCSD\(T\)_L energy:)\s+" + r"(?P<accsdttot>" + NUMBER + r")" + r"\s*$",
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        logger.debug("matched a-ccsd(t) ncc")
        psivar["(T) CORRECTION ENERGY"] = mobj.group("tcorr")
        psivar["T(CCSD) CORRECTION ENERGY"] = mobj.group("bkttcorr")
        psivar["A-(T) CORRECTION ENERGY"] = mobj.group("atcorr")
        psivar["CCSD+T(CCSD) TOTAL ENERGY"] = psivar["T(CCSD) CORRECTION ENERGY"] + psivar["CCSD TOTAL ENERGY"]
        psivar["CCSD+T(CCSD) CORRELATION ENERGY"] = (
            psivar["T(CCSD) CORRECTION ENERGY"] + psivar["CCSD CORRELATION ENERGY"]
        )
        psivar["CCSD(T) TOTAL ENERGY"] = psivar["(T) CORRECTION ENERGY"] + psivar["CCSD TOTAL ENERGY"]
        psivar["CCSD(T) CORRELATION ENERGY"] = psivar["(T) CORRECTION ENERGY"] + psivar["CCSD CORRELATION ENERGY"]
        psivar["A-CCSD(T) TOTAL ENERGY"] = mobj.group("accsdttot")
        psivar["A-CCSD(T) CORRELATION ENERGY"] = psivar["A-(T) CORRECTION ENERGY"] + psivar["CCSD CORRELATION ENERGY"]
        module = "ncc"

    mobj = re.search(
        # fmt: off
        r"^\s*" + r"(?:A miracle come to pass. The CC iterations have converged.)" + r"\s*" +
        r'(?:.*?)' +
        r"^\s*" + r"(?:HF-SCF energy          )\s+" + r"(?P<hf>" + NUMBER + r")" + r"\s*" +
        r"^\s*" + r"(?:MP2 correlation energy )\s+" + r"(?P<mp2corl>" + NUMBER + r")" + r"\s*" +
        r"^\s*" + r"(?:MP2 energy             )\s+" + r"(?P<mp2>" + NUMBER + r")" + r"\s*" +
        r"^\s*" + r"(?:CCSD correlation energy)\s+" + r"(?P<ccsdcorl>" + NUMBER + r")" + r"\s*" +
        r"^\s*" + r"(?:CCSD energy            )\s+" + r"(?P<ccsd>" + NUMBER + r")" + r"\s*" +
        r"^\s*" + r"(?:Delta ET to CCSD\[T\]_L)\s+" + r"(?P<abkttcorr>" + NUMBER + r")" + r"\s*" +
        r"^\s*" + r"(?:CCSD\[T\]_L energy     )\s+" + r"(?P<accsdbktt>" + NUMBER + r")" + r"\s*" +
        r"^\s*" + r"(?:Delta ET to CCSD\(T\)_L)\s+" + r"(?P<atcorr>" + NUMBER + r")" + r"\s*" +
        r"^\s*" + r"(?:CCSD\(T\)_L energy     )\s+" + r"(?P<accsdt>" + NUMBER + r")" + r"\s*" +
        r'(?:.*?)' +
        r"^\s*" + r"(?:Non-iterative perturbative treatment of triple)" + r"\s*" +
        r"^\s*" + r"(?:excitations using the CCSD\(T\)_L method:)" + r"\s*" +
        r'(?:.*?)' +
        r"^\s*" + r"(?:HF-SCF energy)\s+" + r"(?P<hf_again>" + NUMBER + r")" + r"\s*" +
        r"^\s*" + r"(?:MP2 correlation energy)\s+" + r"(?P<mp2corl_again>" + NUMBER + r")" + r"\s*" +
        r"^\s*" + r"(?:MP2 energy)\s+" + r"(?P<mp2_again>" + NUMBER + r")" + r"\s*" +
        r"^\s*" + r"(?:CCSD correlation energy)\s+" + r"(?P<lccsdcorl>" + NUMBER + r")" + r"\s*" +
        r"^\s*" + r"(?:CCSD energy)\s+" + r"(?P<lccsd>" + NUMBER + r")" + r"\s*" +
        r"^\s*" + r"(?:Delta ET to CCSD\(T\)_L)\s+" + r"(?P<atcorr_again>" + NUMBER + r")" + r"\s*" +
        r"^\s*" + r"(?:CCSD\(T\)_L energy)\s+" + r"(?P<laccsdt>" + NUMBER + r")" + r"\s*$",
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        logger.debug("matched a-ccsd(t) ecc")
        psivar["HF TOTAL ENERGY"] = mobj.group("hf")
        psivar["MP2 CORRELATION ENERGY"] = mobj.group("mp2corl")
        psivar["CCSD CORRELATION ENERGY"] = mobj.group("ccsdcorl")
        psivar["A-(T) CORRECTION ENERGY"] = mobj.group("atcorr")
        psivar["A-CCSD(T) TOTAL ENERGY"] = mobj.group("accsdt")
        psivar["A-CCSD(T) CORRELATION ENERGY"] = psivar["A-(T) CORRECTION ENERGY"] + psivar["CCSD CORRELATION ENERGY"]

        mobj3 = re.search(r"SCF reference function:  (R|U)HF", outtext)
        if mobj3:
            psivar["CCSD SINGLES ENERGY"] = Decimal("0.0")
        module = "ecc"

    mobj = re.search(
        # fmt: off
           r"^\s*" + r"(?:E\(CCSD\))" + r"\s+=\s+" + r"(?P<ccsdtot>" + NUMBER + ")" + r"\s*" +
           r"^\s*" + r"(?:E\(CCSD \+ T\(CCSD\)\))" + r"\s*=\s*" + r"(?P<ccsdtccsdtot>" + NUMBER + ")" + r"\s*$",
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        logger.debug("matched ccsd+t(ccsd) vcc")
        psivar["CCSD TOTAL ENERGY"] = mobj.group("ccsdtot")
        psivar["CCSD+T(CCSD) TOTAL ENERGY"] = mobj.group("ccsdtccsdtot")
        psivar["CCSD+T(CCSD) CORRELATION ENERGY"] = psivar["CCSD+T(CCSD) TOTAL ENERGY"] - psivar["SCF TOTAL ENERGY"]
        mobj3 = re.search(r"Reference function is (R|U)HF Hartree-Fock", outtext)
        if mobj3:
            psivar["CCSD SINGLES ENERGY"] = Decimal("0.0")
        module = "vcc"

    mobj = re.search(
        # fmt: off
        r'^\s*' + r'(?:CCSDT\[Q\] correlation energy:)\s+' + r'(?P<bkttcorr>' + NUMBER + ')' + r'\s*'
        r'^\s*' + r'(?:CCSDT\(Q\) correlation energy:)\s+' + r'(?P<tcorr>' + NUMBER + ')' + r'\s*$',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        logger.debug("matched ccsdt(q) ncc")
        psivar["(Q) CORRECTION ENERGY"] = mobj.group("tcorr")
        psivar["[Q] CORRECTION ENERGY"] = mobj.group("bkttcorr")
        psivar["CCSDT(Q) TOTAL ENERGY"] = psivar["(Q) CORRECTION ENERGY"] + psivar["CCSDT TOTAL ENERGY"]
        psivar["CCSDT(Q) CORRELATION ENERGY"] = psivar["(Q) CORRECTION ENERGY"] + psivar["CCSDT CORRELATION ENERGY"]
        module = "ncc"

    # Process DBOC
    mobj = re.search(
        # fmt: off
        r'^\s*' + r'(?:The total diagonal Born-Oppenheimer correction \(DBOC\) is:)\s+' +
        r'(?P<dboc>' + NUMBER + ')' + r'\s*a.u.\s*',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        print("matched dboc ecc")
        psivar["CCSD DBOC ENERGY"] = mobj.group("dboc")
        module = "ecc"

    # Process SCS-CC
    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?P<fullCC>(?P<iterCC>CC(?:\w+))(?:\(T\))?)' + r'\s+(?:energy will be calculated.)\s*' +
        r'(?:.*?)' +
        r'^\s*' + r'(?:@CCENRG-I, Correlation energies.)' + r'\s+(?:ECCAA)\s+' + NUMBER + r'\s*' +
        r'^\s+(?:ECCBB)\s+' + NUMBER + r'\s*' +
        r'^\s+(?:ECCAB)\s+' + NUMBER + r'\s*' +
        r'^\s+(?:Total)\s+' + NUMBER + r'\s*',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:  # PRINT=2 to get SCS-CC components
        print("matched scscc")
        if float(mobj.group(4)) == 0.0:
            ss = 2 * Decimal(mobj.group(3))
        else:
            ss = Decimal(mobj.group(3)) + Decimal(mobj.group(4))

        if not (
            re.search(r"executable xvcc finished", outtext)
            and re.search(r"The reference state is a ROHF wave function.", outtext)
        ):
            psivar["%s SAME-SPIN CORRELATION ENERGY" % (mobj.group("iterCC"))] = ss
        psivar["%s OPPOSITE-SPIN CORRELATION ENERGY" % (mobj.group("iterCC"))] = mobj.group(5)
        psivar["%s CORRELATION ENERGY" % (mobj.group("iterCC"))] = mobj.group(6)

    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?P<fullCC>(?P<iterCC>CC(?:\w+))(?:\(T\))?)' + r'\s+(?:energy will be calculated.)\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'Amplitude equations converged in' + r'\s*\d+\s*' + r'iterations.\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'The AA contribution to the correlation energy is:\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'The BB contribution to the correlation energy is:\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'The AB contribution to the correlation energy is:\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'The total correlation energy is\s+' + NUMBER + r'\s+a.u.\s*' +
        r'(?:.*?)' +
        # r'^\s+' + r'The CC iterations have converged.' + r'\s*$',
        r'^\s+' + r'(?:A miracle come to pass. )?' + r'The CC iterations have converged.' + r'\s*$',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:  # PRINT=2 to get SCS components
        logger.debug("matched scscc2")
        iterCC = mobj.group("iterCC")
        mobj3 = re.search(r"The reference state is a ROHF wave function.", outtext)
        mobj4 = re.search(r"executable xvcc finished", outtext)
        if mobj4:  # vcc
            if not (iterCC == "CCD" and mobj3):
                # uncertain if ROHF CCD correct
                psivar[f"{iterCC} OPPOSITE-SPIN CORRELATION ENERGY"] = mobj.group(5)
            if not mobj3:
                psivar[f"{iterCC} SAME-SPIN CORRELATION ENERGY"] = Decimal(mobj.group(3)) + Decimal(mobj.group(4))
        else:  # ecc
            psivar[f"{iterCC} SAME-SPIN CORRELATION ENERGY"] = Decimal(mobj.group(3)) + Decimal(mobj.group(4))
            if not mobj3:
                psivar[f"{iterCC} OPPOSITE-SPIN CORRELATION ENERGY"] = mobj.group(5)
        psivar[f"{iterCC} CORRELATION ENERGY"] = mobj.group(6)

    mobj = re.search(
        # fmt: off
        #r'^\s+' + r'(?P<fullCC>(?P<iterCC>L?CC(?:\w+))(?:\(T\))?)' + r'\s+(?:energy will be calculated.)\s*' +
        # better one for LCC and one for CC, right?
        r'^\s+' + r'(?P<fullCC>(?P<iterCC>LCC(?:\w+))(?:\(T\))?)' + r'\s+(?:energy will be calculated.)\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'Amplitude equations converged in' + r'\s*\d+\s*' + r'iterations.\s*' +
        r'(?:.*?)' +
         r'^\s+' + r'The AA contribution to the correlation energy is:\s+' + r"(?P<AA>" + NUMBER + r")" + r'\s+a.u.\s*' +
        r'(^\s+' + r'The BB contribution to the correlation energy is:\s+' + r"(?P<BB>" + NUMBER + r")" + r'\s+a.u.\s*' + r")?" +
         r'^\s+' + r'The AB contribution to the correlation energy is:\s+' + r"(?P<AB>" + NUMBER + r")" + r'\s+a.u.\s*' +
        r'^\s+' + r'The total correlation energy is\s+' + r"(?P<corl>" + NUMBER + r")" + r'\s+a.u.\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'(?:A miracle come to pass. )?' + r'The CC iterations have converged.' + r'\s*$',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:  # PRINT=2 to get SCS components
        logger.debug("matched scslccd")
        mobj3 = re.search(r"The reference state is a ROHF wave function.", outtext)
        mobj4 = re.search(r"executable xvcc finished", outtext)
        iterCC = mobj.group("iterCC")
        if mobj4:  # vcc
            if mobj.group("BB"):
                aabb = Decimal(mobj.group("AA")) + Decimal(mobj.group("BB"))
            else:
                aabb = Decimal("2") * Decimal(mobj.group("AA"))
            psivar[f"{iterCC} OPPOSITE-SPIN CORRELATION ENERGY"] = mobj.group("AB")
            if not mobj3:
                psivar[f"{iterCC} SAME-SPIN CORRELATION ENERGY"] = aabb
        psivar["%s CORRELATION ENERGY" % (mobj.group("iterCC"))] = mobj.group("corl")

    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?P<fullCC>(?P<iterCC>CC(?:\w+))(?:\(T\))?)' + r'\s+(?:energy will be calculated.)\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'Amplitude equations converged in' + r'\s*\d+\s*' + r'iterations.\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'The AA contribution to the correlation energy is:\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'The AB contribution to the correlation energy is:\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'The total correlation energy is\s+' + NUMBER + r'\s+a.u.\s*' +
        r'(?:.*?)' +
        # r'^\s+' + r'The CC iterations have converged.' + r'\s*$',
        r'^\s+' + r'(?:A miracle come to pass. )?' + r'The CC iterations have converged.' + r'\s*$',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:  # PRINT=2 to get SCS components
        print("matched scscc rhf", mobj.groups())
        psivar["%s SAME-SPIN CORRELATION ENERGY" % (mobj.group("iterCC"))] = 2 * Decimal(mobj.group(3))
        psivar["%s OPPOSITE-SPIN CORRELATION ENERGY" % (mobj.group("iterCC"))] = mobj.group(4)
        psivar["%s CORRELATION ENERGY" % (mobj.group("iterCC"))] = mobj.group(5)

    # Process gradient
    mobj = re.search(
        # fmt: off
        r'\s+' + r'Molecular gradient' + r'\s*' +
        r'\s+' + r'------------------' + r'\s*' +
        r'\s+' + r'\n' +
        r'(?:(?:\s+[A-Z]+\s*#\d+\s+[xyz]\s+[-+]?\d+\.\d+\s*\n)+)' +  # optional, it seems
        r'\n\n' +  # optional, it seems
        r'((?:\s+[A-Z]+\s*#\d+\s+\d?\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s*\n)+)' +
        r'\n\n' +
        r'\s+' + 'Molecular gradient norm',
        # fmt: on
        outtext,
        re.MULTILINE,
    )
    if mobj:
        print("matched molgrad")
        atoms = []
        psivar_grad = []
        for line in mobj.group(1).splitlines():
            lline = line.split()
            atoms.append(lline[0])
            # psivar_gradient.append([Decimal(lline[-3]), Decimal(lline[-2]), Decimal(lline[-1])])
            psivar_grad.append([float(lline[-3]), float(lline[-2]), float(lline[-1])])

    # Process geometry
    mobj = re.search(
        # fmt: off
        # r'\s+(?:-+)\s*' +
        # r'^\s+' + r'Z-matrix   Atomic            Coordinates (in bohr)' + r'\s*' +
        r'^\s+' + r'Symbol    Number           X              Y              Z' + r'\s*' +
        r'^\s+(?:-+)\s*' +
        r'((?:\s+[A-Z]+\s+([0-9]+|\*\*\*)\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s*\n)+)' +  # allows ghosts
        r'^\s+(?:-+)\s*',
        # fmt: on
        outtext,
        re.MULTILINE,
    )
    if mobj:
        print("matched geom")
        molxyz = "%d bohr\n\n" % len(mobj.group(1).splitlines())
        for line in mobj.group(1).splitlines():
            lline = line.split()
            if lline[1] == "***":
                tag = "@Xe"  # potentially dangerous bypass
            else:
                tag = lline[0]
            molxyz += "%s %16s %16s %16s\n" % (tag, lline[-3], lline[-2], lline[-1])
        # Rather a dinky Molecule as no ghost, charge, or multiplicity
        psivar_coord = Molecule(
            validate=False,
            **qcel.molparse.to_schema(
                qcel.molparse.from_string(molxyz, dtype="xyz+", fix_com=True, fix_orientation=True)["qm"], dtype=2
            ),
        )

    # Process atom geometry
    mobj = re.search(r"^\s+" + r"@GETXYZ-I,     1 atoms read from ZMAT." + r"\s*$", outtext, re.MULTILINE)
    mobj2 = re.search(
        r"^([A-Z]+)#1" + r"\s+" + NUMBER + r"\s+" + NUMBER + r"\s+" + NUMBER + r"\s*$", outtext, re.MULTILINE
    )
    if mobj and mobj2:
        print("matched atom2")  # unsavory for when atom never printed except for basis file
        # Dinky Molecule
        molxyz = "1 bohr\n\n%s 0.0 0.0 0.0\n" % (mobj2.group(1))
        psivar_coord = Molecule(
            validate=False,
            **qcel.molparse.to_schema(
                qcel.molparse.from_string(molxyz, dtype="xyz+", fix_com=True, fix_orientation=True)["qm"], dtype=2
            ),
        )

    mobj = re.search(
        # fmt: off
        r'^\s+' + r'@GETXYZ-I,     1 atoms read from ZMAT.' + r'\s*' +
        r'^\s+' + r'[0-9]+\s+([A-Z]+)\s+[0-9]+\s+' + NUMBER + r'\s*',
        # fmt: on
        outtext,
        re.MULTILINE,
    )
    if mobj:
        print("matched atom")
        # Dinky Molecule
        molxyz = "1 bohr\n\n%s 0.0 0.0 0.0\n" % (mobj.group(1))
        psivar_coord = Molecule(
            validate=False,
            **qcel.molparse.to_schema(
                qcel.molparse.from_string(molxyz, dtype="xyz+", fix_com=True, fix_orientation=True)["qm"], dtype=2
            ),
        )

    # Process error codes
    mobj = re.search(
        # fmt: off
        r"^\s*" + r"--executable " + r"(?P<c4exe>\w+)" + r" finished with status" + r"\s+" + r"(?P<errcode>[1-9][0-9]*)",
        # fmt: on
        outtext,
        re.MULTILINE,
    )
    if mobj:
        print("matched error")
        # psivar['CFOUR ERROR CODE'] = mobj.group(2)
        c4exe = mobj.group("c4exe")
        errcode = int(mobj.group("errcode"))
        if errcode != 0:
            error += f"--executable {c4exe} finished with status {errcode}"
            if c4exe in ["xvcc", "xecc", "xncc"]:
                module = c4exe[1:]

    # Process CURRENT energies (TODO: needs better way)
    if "SCF TOTAL ENERGY" in psivar:
        psivar["CURRENT REFERENCE ENERGY"] = psivar["SCF TOTAL ENERGY"]
        psivar["CURRENT ENERGY"] = psivar["SCF TOTAL ENERGY"]
        psivar["HF TOTAL ENERGY"] = psivar["SCF TOTAL ENERGY"]

    if "MP2 TOTAL ENERGY" in psivar and "MP2 CORRELATION ENERGY" in psivar:
        psivar["CURRENT CORRELATION ENERGY"] = psivar["MP2 CORRELATION ENERGY"]
        psivar["CURRENT ENERGY"] = psivar["MP2 TOTAL ENERGY"]

    if "MP3 TOTAL ENERGY" in psivar and "MP3 CORRELATION ENERGY" in psivar:
        psivar["CURRENT CORRELATION ENERGY"] = psivar["MP3 CORRELATION ENERGY"]
        psivar["CURRENT ENERGY"] = psivar["MP3 TOTAL ENERGY"]

    if "MP4(SDQ) TOTAL ENERGY" in psivar and "MP4(SDQ) CORRELATION ENERGY" in psivar:
        psivar["CURRENT CORRELATION ENERGY"] = psivar["MP4(SDQ) CORRELATION ENERGY"]
        psivar["CURRENT ENERGY"] = psivar["MP4(SDQ) TOTAL ENERGY"]

    if "MP4 TOTAL ENERGY" in psivar and "MP4 CORRELATION ENERGY" in psivar:
        psivar["CURRENT CORRELATION ENERGY"] = psivar["MP4 CORRELATION ENERGY"]
        psivar["CURRENT ENERGY"] = psivar["MP4 TOTAL ENERGY"]

    if "CISD TOTAL ENERGY" in psivar and "CISD CORRELATION ENERGY" in psivar:
        psivar["CURRENT CORRELATION ENERGY"] = psivar["CISD CORRELATION ENERGY"]
        psivar["CURRENT ENERGY"] = psivar["CISD TOTAL ENERGY"]

    if "QCISD TOTAL ENERGY" in psivar and "QCISD CORRELATION ENERGY" in psivar:
        psivar["CURRENT CORRELATION ENERGY"] = psivar["QCISD CORRELATION ENERGY"]
        psivar["CURRENT ENERGY"] = psivar["QCISD TOTAL ENERGY"]

    if "QCISD(T) TOTAL ENERGY" in psivar and "QCISD(T) CORRELATION ENERGY" in psivar:
        psivar["CURRENT CORRELATION ENERGY"] = psivar["QCISD(T) CORRELATION ENERGY"]
        psivar["CURRENT ENERGY"] = psivar["QCISD(T) TOTAL ENERGY"]

    if "LCCD TOTAL ENERGY" in psivar and "LCCD CORRELATION ENERGY" in psivar:
        psivar["CURRENT CORRELATION ENERGY"] = psivar["LCCD CORRELATION ENERGY"]
        psivar["CURRENT ENERGY"] = psivar["LCCD TOTAL ENERGY"]

    if "LCCSD TOTAL ENERGY" in psivar and "LCCSD CORRELATION ENERGY" in psivar:
        psivar["CURRENT CORRELATION ENERGY"] = psivar["LCCSD CORRELATION ENERGY"]
        psivar["CURRENT ENERGY"] = psivar["LCCSD TOTAL ENERGY"]

    #    if ('%s TOTAL ENERGY' % (mobj.group('fullCC')) in psivar) and \
    #       ('%s CORRELATION ENERGY' % (mobj.group('fullCC')) in psivar):
    #        psivar['CURRENT CORRELATION ENERGY'] = psivar['%s CORRELATION ENERGY' % (mobj.group('fullCC')]
    #        psivar['CURRENT ENERGY'] = psivar['%s TOTAL ENERGY' % (mobj.group('fullCC')]

    if "CC2 TOTAL ENERGY" in psivar and "CC2 CORRELATION ENERGY" in psivar:
        psivar["CURRENT CORRELATION ENERGY"] = psivar["CC2 CORRELATION ENERGY"]
        psivar["CURRENT ENERGY"] = psivar["CC2 TOTAL ENERGY"]

    if "CCD TOTAL ENERGY" in psivar and "CCD CORRELATION ENERGY" in psivar:
        psivar["CURRENT CORRELATION ENERGY"] = psivar["CCD CORRELATION ENERGY"]
        psivar["CURRENT ENERGY"] = psivar["CCD TOTAL ENERGY"]

    if "CCSD TOTAL ENERGY" in psivar and "CCSD CORRELATION ENERGY" in psivar:
        psivar["CURRENT CORRELATION ENERGY"] = psivar["CCSD CORRELATION ENERGY"]
        psivar["CURRENT ENERGY"] = psivar["CCSD TOTAL ENERGY"]

    if "CCSD+T(CCSD) TOTAL ENERGY" in psivar and "CCSD+T(CCSD) CORRELATION ENERGY" in psivar:
        psivar["CURRENT CORRELATION ENERGY"] = psivar["CCSD+T(CCSD) CORRELATION ENERGY"]
        psivar["CURRENT ENERGY"] = psivar["CCSD+T(CCSD) TOTAL ENERGY"]

    if "CCSD(T) TOTAL ENERGY" in psivar and "CCSD(T) CORRELATION ENERGY" in psivar:
        psivar["CURRENT CORRELATION ENERGY"] = psivar["CCSD(T) CORRELATION ENERGY"]
        psivar["CURRENT ENERGY"] = psivar["CCSD(T) TOTAL ENERGY"]

    if "A-CCSD(T) TOTAL ENERGY" in psivar and "A-CCSD(T) CORRELATION ENERGY" in psivar:
        psivar["CURRENT CORRELATION ENERGY"] = psivar["A-CCSD(T) CORRELATION ENERGY"]
        psivar["CURRENT ENERGY"] = psivar["A-CCSD(T) TOTAL ENERGY"]

    if "CC3 TOTAL ENERGY" in psivar and "CC3 CORRELATION ENERGY" in psivar:
        psivar["CURRENT CORRELATION ENERGY"] = psivar["CC3 CORRELATION ENERGY"]
        psivar["CURRENT ENERGY"] = psivar["CC3 TOTAL ENERGY"]

    if "CCSDT TOTAL ENERGY" in psivar and "CCSDT CORRELATION ENERGY" in psivar:
        psivar["CURRENT CORRELATION ENERGY"] = psivar["CCSDT CORRELATION ENERGY"]
        psivar["CURRENT ENERGY"] = psivar["CCSDT TOTAL ENERGY"]

    if "CCSDT-1A TOTAL ENERGY" in psivar and "CCSDT-1A CORRELATION ENERGY" in psivar:
        psivar["CURRENT CORRELATION ENERGY"] = psivar["CCSDT-1A CORRELATION ENERGY"]
        psivar["CURRENT ENERGY"] = psivar["CCSDT-1A TOTAL ENERGY"]

    if "CCSDT-1B TOTAL ENERGY" in psivar and "CCSDT-1B CORRELATION ENERGY" in psivar:
        psivar["CURRENT CORRELATION ENERGY"] = psivar["CCSDT-1B CORRELATION ENERGY"]
        psivar["CURRENT ENERGY"] = psivar["CCSDT-1B TOTAL ENERGY"]

    if "CCSDT-2 TOTAL ENERGY" in psivar and "CCSDT-2 CORRELATION ENERGY" in psivar:
        psivar["CURRENT CORRELATION ENERGY"] = psivar["CCSDT-2 CORRELATION ENERGY"]
        psivar["CURRENT ENERGY"] = psivar["CCSDT-2 TOTAL ENERGY"]

    if "CCSDT-3 TOTAL ENERGY" in psivar and "CCSDT-3 CORRELATION ENERGY" in psivar:
        psivar["CURRENT CORRELATION ENERGY"] = psivar["CCSDT-3 CORRELATION ENERGY"]
        psivar["CURRENT ENERGY"] = psivar["CCSDT-3 TOTAL ENERGY"]

    if "CCSDT(Q) TOTAL ENERGY" in psivar and "CCSDT(Q) CORRELATION ENERGY" in psivar:
        psivar["CURRENT CORRELATION ENERGY"] = psivar["CCSDT(Q) CORRELATION ENERGY"]
        psivar["CURRENT ENERGY"] = psivar["CCSDT(Q) TOTAL ENERGY"]

    if "CCSDTQ TOTAL ENERGY" in psivar and "CCSDTQ CORRELATION ENERGY" in psivar:
        psivar["CURRENT CORRELATION ENERGY"] = psivar["CCSDTQ CORRELATION ENERGY"]
        psivar["CURRENT ENERGY"] = psivar["CCSDTQ TOTAL ENERGY"]

    psivar[f"N ATOMS"] = len(psivar_coord.symbols)

    return psivar, psivar_coord, psivar_grad, version, module, error


def harvest(in_mol: Molecule, method: str, c4out, **largs):
    """Parses all the pieces of output from Cfour: the stdout in
    *c4out* and the contents of various scratch files like GRD stored
    in their namesake keys in *largs*. Since all Cfour output uses
    its own orientation and atom ordering for the given molecule,
    a qcdb.Molecule *in_mol*, if supplied, is used to transform the
    Cfour output back into consistency with *in_mol*.

    """
    # Collect results from output file and subsidiary files
    qcvars, out_mol, outGrad, version, module, error = harvest_output(c4out)

    if largs.get("GRD"):
        grdMol, grdGrad = harvest_GRD(largs["GRD"])
    else:
        grdMol, grdGrad = None, None

    if largs.get("FCMFINAL"):
        fcmHess = load_hessian(largs["FCMFINAL"], dtype="fcmfinal")
        if np.count_nonzero(fcmHess) == 0:
            fcmHess = None
    else:
        fcmHess = None

    if largs.get("DIPOL"):
        dipolDip = harvest_DIPOL(largs["DIPOL"])
    else:
        dipolDip = None

    # Sometimes the hierarchical setting of CURRENT breaks down
    if method == "ccsd+t(ccsd)":
        qcvars["CURRENT CORRELATION ENERGY"] = qcvars["CCSD+T(CCSD) CORRELATION ENERGY"]
        qcvars["CURRENT ENERGY"] = qcvars["CCSD+T(CCSD) TOTAL ENERGY"]

    if fcmHess is not None and method in ["hf", "scf"]:
        # MP2 available in HF Hessian so need to counteract
        qcvars.pop("CURRENT CORRELATION ENERGY")
        qcvars["CURRENT ENERGY"] = qcvars["HF TOTAL ENERGY"]

    # Reconcile the coordinate information: several cases
    #   Case                            p4Mol   GRD      Check consistency           Apply orientation?     ReturnMol (1-19-2014)
    #   sp with mol thru cfour {}       None    None              outMol             N.C.                   outMol
    #   opt with mol thru cfour {}      None    grdMol            outMol && grdMol   N.C.                   grdMol
    #   sp with mol thru molecule {}    p4Mol   None     p4Mol && outMol             p4Mol <-- outMol       p4Mol (same as input arg)
    #   opt with mol thru molecule {}   p4Mol   grdMol   p4Mol && outMol && grdMol   p4Mol <-- grdMol       p4Mol (same as input arg)
    # Jul 2021: above describes longtime orientation strategy. Now, mol through cfour {} no longer allowed, and fix_* signal whether input (T) or cfour native (F) frames for returned data.

    if out_mol:
        if grdMol:
            if abs(out_mol.nuclear_repulsion_energy() - grdMol.nuclear_repulsion_energy()) > 1.0e-3:
                raise ValueError(
                    f"""CFOUR outfile (NRE: {out_mol.nuclear_repulsion_energy()} inconsistent with CFOUR GRD (NRE: {grdMol.nuclear_repulsion_energy()})."""
                )
        if in_mol:
            if abs(out_mol.nuclear_repulsion_energy() - in_mol.nuclear_repulsion_energy()) > 1.0e-3:
                raise ValueError(
                    f"""CFOUR outfile (NRE: {out_mol.nuclear_repulsion_energy()}) inconsistent with AtomicInput.molecule (NRE: {in_mol.nuclear_repulsion_energy()})."""
                )
    else:
        raise ValueError("""No coordinate information extracted from CFOUR output.""")

    # Set up array reorientation object(s)
    if in_mol and out_mol and grdMol:
        # Jul 2021: apparently GRD and FCMFINAL can have different atom orderings :-)

        _, data = grdMol.align(out_mol, atoms_map=False, mols_align=True, verbose=0)
        g2o_mill = data["mill"]

        oriCoord = g2o_mill.align_coordinates(grdMol.geometry)
        oriGrad = g2o_mill.align_gradient(np.array(grdGrad))

        if dipolDip is None:
            oriDip = None
        else:
            oriDip = g2o_mill.align_vector(dipolDip)

        if fcmHess is None:
            oriHess = None
        else:
            oriHess = fcmHess

        # Frame considerations
        if in_mol.fix_com and in_mol.fix_orientation:
            # Impose input frame if important as signalled by fix_*=T
            return_mol = in_mol
            _, data = out_mol.align(in_mol, atoms_map=False, mols_align=True, verbose=0)
            mill = data["mill"]

        else:
            return_mol, _ = in_mol.align(out_mol, atoms_map=False, mols_align=True, verbose=0)
            mill = qcel.molutil.compute_scramble(
                len(in_mol.symbols), do_resort=False, do_shift=False, do_rotate=False, do_mirror=False
            )  # identity AlignmentMill

        #        _, data = out_mol.align(in_mol, atoms_map=False, mols_align=True, verbose=0)
        #        o2i_mill = data["mill"]
        #
        #        oriCoord = o2i_mill.align_coordinates(oriCoord)
        #        oriGrad = o2i_mill.align_gradient(oriGrad)
        #        if oriDip is not None:
        #            oriDip = o2i_mill.align_vector(oriDip)
        #        if oriHess is not None:
        #            oriHess = o2i_mill.align_hessian(oriHess)

        oriCoord = mill.align_coordinates(oriCoord)
        oriGrad = mill.align_gradient(oriGrad)
        if oriDip is not None:
            oriDip = mill.align_vector(oriDip)
        if oriHess is not None:
            oriHess = mill.align_hessian(oriHess)

    elif in_mol and out_mol:
        # TODO watch out - haven't seen atom_map=False yet

        if in_mol.fix_com and in_mol.fix_orientation:
            # Impose input frame if important as signalled by fix_*=T
            return_mol = in_mol
            _, data = out_mol.align(in_mol, atoms_map=True, mols_align=True, generic_ghosts=True, verbose=0)
            mill = data["mill"]

        else:
            return_mol, _ = in_mol.align(out_mol, atoms_map=False, mols_align=True, generic_ghosts=True, verbose=0)
            mill = qcel.molutil.compute_scramble(
                len(in_mol.symbols), do_resort=False, do_shift=False, do_rotate=False, do_mirror=False
            )  # identity AlignmentMill

        oriCoord = mill.align_coordinates(out_mol.geometry)  # (np_out=True))
        oriGrad = None
        oriHess = None  # I don't think we ever get FCMFINAL w/o GRAD
        if dipolDip is None:
            oriDip = None
        else:
            oriDip = mill.align_vector(np.array(dipolDip))
        # p4c4 = OrientMols(in_mol, out_mol)
        # oriCoord = p4c4.transform_coordinates2(out_mol)
        # oriGrad = None
        # oriDip = None if dipolDip is None else p4c4.transform_vector(dipolDip)

    elif out_mol:
        oriGrad = None
        oriHess = None
        oriDip = None if dipolDip is None else dipolDip

    # not sure of purpose but it interferes now that return_mol overwrites atres.mol
    # return_mol = None if in_mol else grdMol

    if oriDip is not None:
        qcvars["CURRENT DIPOLE"] = oriDip
        oriDip *= qcel.constants.dipmom_au2debye
        # outPsivar["CURRENT DIPOLE X"] = oriDip[0]
        # outPsivar["CURRENT DIPOLE Y"] = oriDip[1]
        # outPsivar["CURRENT DIPOLE Z"] = oriDip[2]
        # outPsivar['CURRENT DIPOLE X'] = str(oriDip[0] * psi_dipmom_au2debye)
        # outPsivar['CURRENT DIPOLE Y'] = str(oriDip[1] * psi_dipmom_au2debye)
        # outPsivar['CURRENT DIPOLE Z'] = str(oriDip[2] * psi_dipmom_au2debye)

    if oriGrad is not None:
        return_grad = oriGrad
    elif grdGrad is not None:
        return_grad = grdGrad
    else:
        return_grad = None

    if oriHess is not None:
        return_hess = oriHess
    else:
        return_hess = None

    # if oriCoord is not None:
    #     retCoord = oriCoord
    # else:
    #     retCoord = None

    return qcvars, return_hess, return_grad, return_mol, version, module, error


def harvest_GRD(grd):
    """Parses the contents *grd* of the Cfour GRD file into the gradient
    array and coordinate information. The coordinate info is converted
    into a rather dinky Molecule (no charge, multiplicity, or fragment),
    but this is these coordinates that govern the reading of molecule
    orientation by Cfour. Return qcel.models.Molecule and gradient array.

    """
    grd = grd.splitlines()
    Nat = int(grd[0].split()[0])
    molxyz = f"{Nat} bohr\n\n"

    grad = []
    for at in range(Nat):
        mline = grd[at + 1].split()
        el = "GH" if int(float(mline[0])) == 0 else qcel.periodictable.to_E(int(float(mline[0])))
        molxyz += "%s %16s %16s %16s\n" % (el, mline[-3], mline[-2], mline[-1])
        lline = grd[at + 1 + Nat].split()
        grad.append([float(lline[-3]), float(lline[-2]), float(lline[-1])])
    grad = np.array(grad).reshape((-1, 3))

    mol = Molecule(
        validate=False,
        **qcel.molparse.to_schema(
            qcel.molparse.from_string(molxyz, dtype="xyz+", fix_com=True, fix_orientation=True)["qm"], dtype=2
        ),
    )

    return mol, grad


def harvest_DIPOL(dipol):
    """Parses the contents *dipol* of the Cfour DIPOL file into a dipol vector."""
    dipol = dipol.splitlines()
    lline = dipol[0].split()
    dip = np.array([float(lline[0]), float(lline[1]), float(lline[2])])

    # return None if empty else dip
    return dip
