"""Compute quantum chemistry using Iowa State's GAMESS executable."""

import logging
import pprint
import re
from decimal import Decimal
from typing import Dict, Optional, Tuple

import numpy as np
import qcelemental as qcel
from qcelemental.models import Molecule
from qcelemental.molparse import regex

from ..util import PreservingDict, load_hessian

pp = pprint.PrettyPrinter(width=120, compact=True, indent=1)
logger = logging.getLogger(__name__)


def harvest(
    in_mol: Molecule, method: str, gamessout: str, **outfiles
) -> Tuple[PreservingDict, Optional[np.ndarray], Optional[np.ndarray], Molecule]:
    """Parses all the pieces of output from gamess: the stdout in
    *gamessout* Scratch files are not yet considered at this moment.
    """
    qcvars, calc_mol, calc_grad, module = harvest_output(gamessout)

    datasections = {}
    if outfiles.get("gamess.dat"):
        datasections = harvest_datfile(outfiles["gamess.dat"])

    calc_hess = None
    if "$HESS" in datasections:
        calc_hess = load_hessian(datasections["$HESS"], dtype="gamess")
        if np.count_nonzero(calc_hess) == 0:
            calc_hess = None

    # Sometimes the hierarchical setting of CURRENT breaks down
    if method == "ccsd+t(ccsd)":
        qcvars["CURRENT CORRELATION ENERGY"] = qcvars["CCSD+T(CCSD) CORRELATION ENERGY"]
        qcvars["CURRENT ENERGY"] = qcvars["CCSD+T(CCSD) TOTAL ENERGY"]

    if calc_mol:
        qcvars["NUCLEAR REPULSION ENERGY"] = str(round(calc_mol.nuclear_repulsion_energy(), 8))
        if in_mol:
            if abs(calc_mol.nuclear_repulsion_energy() - in_mol.nuclear_repulsion_energy()) > 1.0e-3:
                raise ValueError(
                    f"""GAMESS outfile (NRE: {calc_mol.nuclear_repulsion_energy()}) inconsistent with AtomicInput.molecule (NRE: {in_mol.nuclear_repulsion_energy()})."""
                )

        # Frame considerations
        # * `in_mol` built with deliberation and with all fields accessible.
        # * `calc_mol` has the internally consistent geometry frame but otherwise dinky (geom & symbols & maybe chgmult).
        if in_mol.fix_com and in_mol.fix_orientation:
            # Impose input frame if important as signalled by fix_*=T
            return_mol = in_mol
            _, data = calc_mol.align(in_mol, atoms_map=False, mols_align=True, run_mirror=True, verbose=0)
            mill = data["mill"]

        else:
            return_mol, _ = in_mol.align(calc_mol, atoms_map=False, mols_align=True, run_mirror=True, verbose=0)
            mill = qcel.molutil.compute_scramble(
                len(in_mol.symbols), do_resort=False, do_shift=False, do_rotate=False, do_mirror=False
            )  # identity AlignmentMill

        return_grad = None
        if calc_grad is not None:
            return_grad = mill.align_gradient(calc_grad)

        return_hess = None
        if calc_hess is not None:
            return_hess = mill.align_hessian(np.array(calc_hess))

    else:
        raise ValueError("""No coordinate information extracted from gamess output.""")

    return qcvars, return_hess, return_grad, return_mol, module


def harvest_datfile(datfile: str) -> Dict[str, str]:
    sections = datfile.split(r"$END")
    goodies = {}

    for i, sec in enumerate(sections):
        lsec = sec.split("\n")
        for iln, ln in enumerate(lsec):
            if ln.strip():
                key = ln.strip()
                break
        goodies[key] = "\n".join(lsec[iln + 1 :])

    return goodies


def harvest_output(outtext):
    """Function to separate portions of a gamess output file *outtext*,
    divided by "Step".
    """
    pass_qcvar = []
    pass_coord = []
    pass_grad = []

    for outpass in re.split(
        # fmt: off
        r"^\s+" + r"--------" + r"NSERCH:" + r"([1-9][0-9][0-9][0-9]*)" + r"\s*" +
        r"^\s+" + r"--------",
        # fmt: on
        outtext,
        re.MULTILINE,
    ):
        qcvar, gamesscoord, gamessgrad, module = harvest_outfile_pass(outpass)
        pass_qcvar.append(qcvar)
        pass_coord.append(gamesscoord)
        pass_grad.append(gamessgrad)

    retindx = -1 if pass_coord[-1] else -2
    return pass_qcvar[retindx], pass_coord[retindx], pass_grad[retindx], module


def harvest_outfile_pass(outtext):
    """Function to read gamess output file *outtext* and parse important
    quantum chemical information from it in
    """
    qcvar = PreservingDict()
    qcvar_coord = None
    qcvar_grad = None
    module = None

    NUMBER = r"(?x:" + regex.NUMBER + ")"

    # If calculation fail to converge
    mobj = re.search(r"^\s+" + r"(?:GAMESS TERMINATED ABNORMALLY)" + r"\s*$", outtext, re.MULTILINE)
    if mobj:
        logger.debug("GAMESS TERMINATED ABNORMALLY")

    # If calculation converged
    else:
        mobj = re.search(
            r"^\s+" + r"(?:            TOTAL ENERGY)" + r"\s+=\s*" + NUMBER + r"s*$", outtext, re.MULTILINE
        )
        if mobj:
            logger.debug("matched gamess_RHF energy")
            qcvar["HF TOTAL ENERGY"] = mobj.group(1)
            qcvar["SCF TOTAL ENERGY"] = mobj.group(1)

        # Process NRE
        mobj = re.search(
            r"^\s+" + r"(?:   NUCLEAR REPULSION ENERGY)" + r"\s+=\s*" + NUMBER + r"\s*$", outtext, re.MULTILINE
        )
        if mobj:
            logger.debug("matched NRE")
            qcvar["NUCLEAR REPULSION ENERGY"] = mobj.group(1)

        # Process calcinfo
        mobj = re.search(
            # fmt: off
            r"^\s+" + r"NUMBER OF OCCUPIED ORBITALS \(ALPHA\)" + r"\s+=\s+" + r"(?P<nao>\d+)" + r"\s*" +
            r"^\s+" + r"NUMBER OF OCCUPIED ORBITALS \(BETA \)" + r"\s+=\s+" + r"(?P<nbo>\d+)" + r"\s*" +
            r"^\s+" + r"TOTAL NUMBER OF ATOMS" + r"\s+=\s+" + r"(?P<nat>\d+)" + r"\s*",
            # fmt: on
            outtext,
            re.MULTILINE,
        )
        if mobj:
            logger.debug("matched calcinfo")
            print("matched calcinfo", mobj.groups())
            qcvar["N ALPHA ELECTRONS"] = mobj.group("nao")
            qcvar["N BETA ELECTRONS"] = mobj.group("nbo")

        mobj = re.search(
            # fmt: off
            r"^\s+" + r"TOTAL NUMBER OF MOS IN VARIATION SPACE=" + r"\s+" + r"(?P<nmo>\d+)" + r"\s*$",
            # fmt: on
            outtext,
            re.MULTILINE,
        )
        if mobj:
            logger.debug("matched calcinfo 2")
            print("matched calcinfo 2", mobj.groups())
            qcvar["N MOLECULAR ORBITALS"] = mobj.group("nmo")
            qcvar["N BASIS FUNCTIONS"] = mobj.group("nmo")  # TODO BAD
        else:
            mobj2 = re.search(
                # fmt: off
                r"^\s+" + r"NUMBER OF CORE -A-  ORBITALS" + r"\s+=\s+" + r"(?P<core_nao>\d+)" + r"\s*" +
                r"^\s+" + r"NUMBER OF CORE -B-  ORBITALS" + r"\s+=\s+" + r"(?P<core_nbo>\d+)" + r"\s*" +
                r"^\s+" + r"NUMBER OF OCC. -A-  ORBITALS" + r"\s+=\s+" + r"(?P<occ_nao>\d+)" + r"\s*" +
                r"^\s+" + r"NUMBER OF OCC. -B-  ORBITALS" + r"\s+=\s+" + r"(?P<occ_nbo>\d+)" + r"\s*" +
                r"^\s+" + r"NUMBER OF MOLECULAR ORBITALS" + r"\s+=\s+" + r"(?P<nmo>\d+)" + r"\s*" +
                r"^\s+" + r"NUMBER OF   BASIS  FUNCTIONS" + r"\s+=\s+" + r"(?P<nbf>\d+)" + r"\s*",
                # fmt: on
                outtext,
                re.MULTILINE,
            )
            if mobj2:
                logger.debug("matched calcinfo 3")
                print("matched calcinfo 3", mobj2.groups())
                qcvar["N ALPHA ELECTRONS"] = mobj2.group("occ_nao")
                qcvar["N BETA ELECTRONS"] = mobj2.group("occ_nbo")
                qcvar["N MOLECULAR ORBITALS"] = mobj2.group("nmo")
                qcvar["N BASIS FUNCTIONS"] = mobj2.group("nbf")

        # Process MP2
        mobj = re.search(
            # fmt: off
            r"^\s+" + r"MP2 CONTROL INFORMATION" + r"\s*" +
            r"^\s+" + r"-+" + r"\s*" +
            r"^\s+" + r"NACORE" + r"\s+=\s+" + r"\d+" + r"\s+" + r"NBCORE" + r"\s+=\s+" + r"\d+" + r"\s*" +
            r"^\s+" + r"LMOMP2" + r"\s+=\s+" + r"\w+" + r"\s+" + r"AOINTS" + r"\s+=\s+" + r"\w+" + r"\s*" +
            r"^\s+" + r"METHOD" + r"\s+=\s+" + r"\d+" + r"\s+" + r"NWORD" + r"\s+=\s+" + r"\d+" + r"\s*" +
            r"^\s+" + r"MP2PRP" + r"\s+=\s+" + r"\w+" + r"\s+" + r"OSPT" + r"\s+=\s+" + r"\w+"+ r"\s*" +
            r"^\s+" + r"CUTOFF" + r"\s+=\s+" + NUMBER + r"\s+" + r"CPHFBS" + r"\s+=\s+" + r"\w+"+ r"\s*" +
            r"^\s+" + r"CODE" + r"\s+=\s+" + r"(?P<code>\w+)"  + r"\s*$",
            # fmt: on
            outtext,
            re.MULTILINE,
        )
        if mobj:
            logger.debug("matched mp2 module")
            module = mobj.group("code").lower()

        mobj = re.search(
            # fmt: off
            r'^\s+' + r'RESULTS OF MOLLER-PLESSET 2ND ORDER CORRECTION ARE\n'
            r'^\s+' + r'E\(0\)' + r'=\s+' + NUMBER + r'\s*' +
            r'^\s+' + r'E\(1\)' + r'=\s+' + NUMBER + r'\s*' +
            r'^\s+' + r'E\(2\)' + r'=\s+' + NUMBER + r'\s*' +
            r'^\s+' + r'E\(MP2\)' + r'=\s+' + NUMBER + r'\s*$',
            # fmt: on
            outtext,
            re.MULTILINE,
        )
        if mobj:
            logger.debug("matched mp2 a")
            print("matched mp2 a", mobj.groups())
            qcvar["MP2 CORRELATION ENERGY"] = mobj.group(3)
            qcvar["MP2 TOTAL ENERGY"] = mobj.group(4)
            mobj3 = re.search(r"\s+[RU]HF\s*SCF\s*CALCULATION", outtext, re.MULTILINE)
            if mobj3:
                qcvar[f"MP2 DOUBLES ENERGY"] = mobj.group(3)

        mobj = re.search(
            # fmt: off
            r'^\s+' + r'RESULTS OF 2ND-ORDER ZAPT CORRECTION' + r'\s*' +
            r'^\s+' + r'E\(HF\)   ' + r'=\s+' + NUMBER + r'\s*' +
            r'^\s+' + r'E\(ZAPT\) ' + r'=\s+' + NUMBER + r'\s*' +
            r'^\s+' + r'-----------------------------------' + r'\s*' +
            r'^\s+' + r'E\(MP2\)  ' + r'=\s+' + NUMBER + r'\s*$',
            # fmt: on
            outtext,
            re.MULTILINE,
        )
        if mobj:
            logger.debug("matched mp2 b")
            print("matched mp2 b")
            qcvar["MP2 CORRELATION ENERGY"] = mobj.group(2)
            qcvar["MP2 TOTAL ENERGY"] = mobj.group(3)

        mobj = re.search(
            # fmt: off
            r'^\s+' + r'RESULTS OF MOLLER-PLESSET 2ND ORDER CORRECTION ARE\n'
            r'^\s+' + r'E\(0\)=' + r'\s+' + NUMBER + r'\s*' +
            r'^\s+' + r'E\(1\)=' + r'\s+' + NUMBER + r'\s*' +
            r'^\s+' + r'E\(2\)=' + r'\s+' + NUMBER + r'\s*' +
            r'^\s+' + r'E\(MP2\)=' + r'\s+' + NUMBER + r'\s*' +
            r'^\s+' + r'SPIN-COMPONENT-SCALED MP2 RESULTS ARE' + r'\s*' +
            r'^\s+' + r'E\(2S\)=' + r'\s+' + NUMBER + r'\s*' +
            r'^\s+' + r'E\(2T\)=' + r'\s+' + NUMBER + r'\s*' +
            r'^\s+' + r'E\(2ST\)=' + r'\s+' + NUMBER + r'\s*',
            # fmt: on
            outtext,
            re.MULTILINE,
        )
        if mobj:
            logger.debug("matched mp2 rhf c")
            print("matched mp2 c", mobj.groups())
            qcvar["HF TOTAL ENERGY"] = mobj.group(1)
            qcvar["MP2 SINGLES ENERGY"] = mobj.group(2)
            qcvar["MP2 DOUBLES ENERGY"] = mobj.group(3)
            qcvar["MP2 TOTAL ENERGY"] = mobj.group(4)
            qcvar["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = mobj.group(5)
            qcvar["MP2 SAME-SPIN CORRELATION ENERGY"] = mobj.group(6)

        mobj = re.search(
            # fmt: off
            r'^\s+' + r'SINGLE EXCITATION CONTRIBUTION' + r'\s*' +
            r'^\s+' + r'ALPHA' + r'\s*' + NUMBER + r'\s*' +
            r'^\s+' + r'BETA' + r'\s*' + NUMBER + r'\s*' +
            r'^\s+' + r'DOUBLE EXCITATION CONTRIBUTION' + r'\s*' +
            r'^\s+' + NUMBER + r'\s*' +
            r'^\s+' +
            r'^\s+' + r'RESULTS OF MOLLER-PLESSET 2ND ORDER CORRECTION ARE' + r'\s*$',
            # fmt: on
            outtext,
            re.MULTILINE,
        )
        if mobj:
            logger.debug("matched mp2 rohf d")
            print("matched mp2 rohf d", mobj.groups())
            qcvar["MP2 SINGLES ENERGY"] = Decimal(mobj.group(1)) + Decimal(mobj.group(2))
            qcvar["MP2 DOUBLES ENERGY"] = mobj.group(3)

        mobj = re.search(r"^\s+" + "UHF-MP2 CALCULATION", outtext, re.MULTILINE)
        if mobj:
            logger.debug("matched mp2 uhf e")
            print("matched mp2 uhf e")
            qcvar["MP2 SINGLES ENERGY"] = "0.0"

        mobj = re.search(
            # fmt: off
            r"^\s+" + r"E\(SCF\)" + r"\s*=\s+" + NUMBER + r"\s*" +
            r"^\s+" + r"ZAPT E\(2\)" + r"\s*=\s+" + NUMBER + r"\s*" +
            r"^\s+" + r"E\(MP2\)" + r"\s*=\s+" + NUMBER + r"\s*",
            # fmt: on
            outtext,
            re.MULTILINE,
        )
        if mobj:
            logger.debug("matched mp2 rohf f")
            print("matched mp2 rohf f", mobj.groups())
            qcvar["HF TOTAL ENERGY"] = mobj.group(1)
            qcvar["MP2 CORRELATION ENERGY"] = mobj.group(2)
            qcvar["MP2 TOTAL ENERGY"] = mobj.group(3)

        mobj = re.search(
            # fmt: off
            r'^\s+' + r'RESULTS OF MOLLER-PLESSET 2ND ORDER CORRECTION ARE\n'
            r'^\s+' + r'E\(SCF\)=' + r'\s+' + NUMBER + r'\s*' +
            r'^\s+' + r'E\(1\)=' + r'\s+' + NUMBER + r'\s*' +
            r'^\s+' + r'E\(2\)=' + r'\s+' + NUMBER + r'\s*' +
            r'^\s+' + r'E\(MP2\)=' + r'\s+' + NUMBER + r'\s*',
            # fmt: on
            outtext,
            re.MULTILINE,
        )
        if mobj:
            logger.debug("matched mp2 g llel")
            qcvar["HF TOTAL ENERGY"] = mobj.group(1)
            qcvar["MP2 SINGLES ENERGY"] = mobj.group(2)
            qcvar["MP2 DOUBLES ENERGY"] = mobj.group(3)
            qcvar["MP2 CORRELATION ENERGY"] = mobj.group(3)
            qcvar["MP2 TOTAL ENERGY"] = mobj.group(4)

        mobj = re.search(
            # fmt: off
            r'^\s+' + r'RESULTS OF MOLLER-PLESSET 2ND ORDER CORRECTION ARE\n'
            r'^\s+' + r'E\(SCF\)=' + r'\s+' + NUMBER + r'\s*' +
            r'^\s+' + r'E\(2\)=' + r'\s+' + NUMBER + r'\s*' +
            r'^\s+' + r'E\(MP2\)=' + r'\s+' + NUMBER + r'\s*' +
            r'^\s+' + r'SPIN-COMPONENT-SCALED MP2 RESULTS ARE' + r'\s*' +
            r'^\s+' + r'E\(2S\)=' + r'\s+' + NUMBER + r'\s*' +
            r'^\s+' + r'E\(2T\)=' + r'\s+' + NUMBER + r'\s*' +
            r'^\s+' + r'E\(2ST\)=' + r'\s+' + NUMBER + r'\s*',
            # fmt: on
            outtext,
            re.MULTILINE,
        )
        if mobj:
            logger.debug("matched mp2 rhf h llel")
            qcvar["HF TOTAL ENERGY"] = mobj.group(1)
            qcvar["MP2 SINGLES ENERGY"] = "0.0"
            qcvar["MP2 DOUBLES ENERGY"] = mobj.group(2)
            qcvar["MP2 CORRELATION ENERGY"] = mobj.group(2)
            qcvar["MP2 TOTAL ENERGY"] = mobj.group(3)
            qcvar["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = mobj.group(4)
            qcvar["MP2 SAME-SPIN CORRELATION ENERGY"] = mobj.group(5)

        # Process CCSD
        mobj = re.search(
            # fmt: off
            r'^\s+' + r'SUMMARY OF RESULTS' + r'\s+' + r'\n' +
            r'^\s+' + r'REFERENCE ENERGY:' + r'\s+' + NUMBER + r'\s*' +
            r'^\s+' + r'MBPT\(2\) ENERGY:' + r'\s+' + NUMBER + r'\s*' + r'CORR.E=\s+' + r'\s+' + NUMBER + r'\s*' +
            r'^\s+' + r"(?P<mtd>CCSD|LCCD|CCD)" + r"\s+" r'ENERGY:' + r'\s+' + NUMBER + r'\s*' + r'CORR.E=\s+' + r'\s+' + NUMBER + r'\s*$',
            # fmt: on
            outtext,
            re.MULTILINE,
        )
        if mobj:
            logger.debug("matched rhf ccsd/lccd/ccd")
            mtd = mobj.group("mtd")
            qcvar["HF TOTAL ENERGY"] = mobj.group(1)
            qcvar["SCF TOTAL ENERGY"] = mobj.group(1)
            qcvar["MP2 CORRELATION ENERGY"] = mobj.group(3)
            qcvar["MP2 DOUBLES ENERGY"] = mobj.group(3)
            qcvar["MP2 TOTAL ENERGY"] = mobj.group(2)
            qcvar[f"{mtd} DOUBLES ENERGY"] = mobj.group(6)
            qcvar[f"{mtd} CORRELATION ENERGY"] = mobj.group(6)
            qcvar[f"{mtd} TOTAL ENERGY"] = mobj.group(5)

        mobj = re.search(
            # fmt: off
            r'^\s+' + r'SUMMARY OF CCSD RESULTS' + r'\s+' + r'\n' +
            r'^\s+' + r'REFERENCE ENERGY:' + r'\s+' + NUMBER + r'\s*' +
            r'^\s+' + r'CCSD ENERGY:'   + r'\s+' + NUMBER + r'\s*' + r'CORR. E=\s+' + r'\s+' + NUMBER + r'\s*$',
            # fmt: on
            outtext,
            re.MULTILINE,
        )
        if mobj:
            logger.debug("matched rohf ccsd")
            qcvar["SCF TOTAL ENERGY"] = mobj.group(1)
            qcvar["CCSD CORRELATION ENERGY"] = mobj.group(3)
            qcvar["CCSD TOTAL ENERGY"] = mobj.group(2)

        # Process CR-CC(2,3)
        mobj = re.search(
            # fmt: off
            r'^\s+' + r'CCSD                       ENERGY:'       + r'\s+' + NUMBER + r'\s*' + r'CORR.E=\s+' + r'\s+' + NUMBER + r'\s*' +
            r'^\s+' + r'CR-CC\(2,3\),A OR CCSD\(2\)_T  ENERGY:'   + r'\s+' + NUMBER + r'\s*' + r'CORR.E=\s+' + r'\s+' + NUMBER + r'\s*' +
            r'^\s+' + r'CR-CC\(2,3\) OR CR-CCSD\(T\)_L ENERGY:'   + r'\s+' + NUMBER + r'\s*' + r'CORR.E=\s+' + r'\s+' + NUMBER + r'\s*',
            # fmt: on
            outtext,
            re.MULTILINE,
        )
        if mobj:
            logger.debug("matched cc-cr(2,3)")
            qcvar["CCSD TOTAL ENERGY"] = mobj.group(1)
            qcvar["CCSD CORRELATION ENERGY"] = mobj.group(2)
            qcvar["CR-CC(2,3),A TOTAL ENERGY"] = mobj.group(3)
            qcvar["CR-CC(2,3),A CORRELATION ENERGY"] = mobj.group(4)
            qcvar["CR-CC(2,3) TOTAL ENERGY"] = mobj.group(5)
            qcvar["CR-CC(2,3) CORRELATION ENERGY"] = mobj.group(6)

        # Process CCSD(T)
        mobj = re.search(
            # fmt: off
            r'^\s+' + r'SUMMARY OF RESULTS' + r'\s+' + r'\n' +
            r'^\s+' + r'REFERENCE ENERGY:' + r'\s+' + NUMBER + r'\s*' +
            r'^\s+' + r'MBPT\(2\) ENERGY:' + r'\s+' + NUMBER + r'\s*' + r'CORR.E=\s+' + r'\s+' + NUMBER + r'\s*' +
            r'^\s+' + r'CCSD    ENERGY:'   + r'\s+' + NUMBER + r'\s*' + r'CORR.E=\s+' + r'\s+' + NUMBER + r'\s*' +
            r'^\s+' + r'CCSD\[T\] ENERGY:' + r'\s+' + NUMBER + r'\s*' + r'CORR.E=\s+' + r'\s+' + NUMBER + r'\s*' +
            r'^\s+' + r'CCSD\(T\) ENERGY:' + r'\s+' + NUMBER + r'\s*' + r'CORR.E=\s+' + r'\s+' + NUMBER + r'\s*$',
            # fmt: on
            outtext,
            re.MULTILINE,
        )
        if mobj:
            logger.debug("matched ccsd(t)")
            qcvar["HF TOTAL ENERGY"] = mobj.group(1)
            qcvar["SCF TOTAL ENERGY"] = mobj.group(1)
            qcvar["CCSD CORRELATION ENERGY"] = mobj.group(5)
            qcvar["CCSD TOTAL ENERGY"] = mobj.group(4)
            qcvar["(T) CORRECTION ENERGY"] = Decimal(mobj.group(8)) - Decimal(mobj.group(4))
            # qcvar['CCSD(T) CORRELATION ENERGY'] = Decimal(mobj.group(5)) - qcvar['SCF TOTAL ENERGY']
            qcvar["CCSD+T(CCSD) CORRELATION ENERGY"] = mobj.group(7)
            qcvar["CCSD+T(CCSD) TOTAL ENERGY"] = mobj.group(6)
            qcvar["CCSD(T) CORRELATION ENERGY"] = mobj.group(9)
            qcvar["CCSD(T) TOTAL ENERGY"] = mobj.group(8)

        # Process CISD
        mobj = re.search(
            # fmt: off
            r"^\s+" + r"PROPERTY VALUES FOR THE " + r"(RHF|ROHF)" + r"\s+" + r"SELF-CONSISTENT FIELD WAVEFUNCTION" + r"\s*" +
            r"(?:.*?)" +
            r"^\s+" + r"ENERGY COMPONENTS" + r"\s*" +
            r"(?:.*?)" +
            r"^\s*" + r"(TOTAL ENERGY =)\s+" + r"(?P<hf>" + NUMBER + r")" + r"\s*" +
            r"(?:.*?)" +
            r"^\s+" + r"(?P<module>(FSOCI|GUGA))" + r"\s+" + r"CI PROPERTIES...FOR THE WAVEFUNCTION OF STATE    1" + r"\s*" +
            r"(?:.*?)" +
            r"^\s+" + r"ENERGY COMPONENTS" + r"\s*" +
            r"(?:.*?)" +
            r"^\s*" + r"(TOTAL ENERGY =)\s+" + r"(?P<ci>" + NUMBER + r")" + r"\s*$",
            # fmt: on
            outtext,
            re.MULTILINE | re.DOTALL,
        )
        if mobj:
            logger.debug("matched cisd fsoci/guga", mobj.groupdict())
            print("matched cisd fsoci/guga", mobj.groupdict())
            module = mobj.group("module").lower()
            qcvar["CISD CORRELATION ENERGY"] = Decimal(mobj.group("ci")) - Decimal(mobj.group("hf"))
            qcvar["CISD TOTAL ENERGY"] = mobj.group("ci")
            qcvar["CI TOTAL ENERGY"] = mobj.group("ci")

        # Process FCI
        mobj = re.search(
            # fmt: off
            r'^\s+' + r'ALDET CI PROPERTIES...FOR THE WAVEFUNCTION OF STATE' + r'\s+' + r'(?P<state>\d+)' + r'\s*' +
            r'^\s+' + r'USING THE EXPECTATION VALUE DENSITY' + r'\s*' +
            r'(?:.*?)' +
            r'^\s+' + r'TOTAL ENERGY =' + r'\s+' + NUMBER + r'\s*$',
            # fmt: on
            outtext,
            re.MULTILINE | re.DOTALL,
        )
        if mobj:
            logger.debug("matched fci")
            qcvar["FCI TOTAL ENERGY"] = mobj.group(2)
            qcvar["CI TOTAL ENERGY"] = mobj.group(2)

        # Process EOM-CCSD
        mobj = re.search(
            # fmt: off
            r'^\s+' + r'---- SUMMARY OF EOM-CCSD CALCULATIONS ----' +
            r'\s+' + r'\n' +
            r'^\s+' + r'EXCITATION      EXCITATION      TOTAL' + r'\s*' +
            r'^\s+' + r'SYMMETRY     ENERGY \(H\)      ENERGY \(EV\)     ENERGY \(H\)          ITERATIONS' + r'\s*' +
            r'^\s+' + r'A' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + r'CONVERGED\s+' +
            r'^\s+' + r'A' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + r'CONVERGED\s+' +
            r'^\s+' + r'A' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + r'CONVERGED\s+' +
            r'^\s+' + r'A' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + r'CONVERGED\s+',
            # fmt: on
            outtext,
            re.MULTILINE,
        )
        if mobj:
            logger.debug("matched eom-ccsd")
            qcvar["EOM-CCSD ROOT 1 EXCITATION ENERGY"] = mobj.group(1)
            qcvar["EOM-CCSD ROOT 2 EXCITATION ENERGY"] = mobj.group(4)
            qcvar["EOM-CCSD ROOT 3 EXCITATION ENERGY"] = mobj.group(7)
            qcvar["EOM-CCSD ROOT 4 EXCITATION ENERGY"] = mobj.group(10)
            qcvar["EOM-CCSD ROOT 1 TOTAL ENERGY"] = mobj.group(3)
            qcvar["EOM-CCSD ROOT 2 TOTAL ENERGY"] = mobj.group(6)
            qcvar["EOM-CCSD ROOT 3 TOTAL ENERGY"] = mobj.group(9)
            qcvar["EOM-CCSD ROOT 4 TOTAL ENERGY"] = mobj.group(12)

        # Process DFT
        # mobj = re.search(
        #     r'^\s+' + r'DFT EXCHANGE + CORRELATION ENERGY' + r'=\s+' + NUMBER + r'\s*$'
        #     ,outtext, re.MULTILINE)
        mobj = re.search(
            r"^\s+" + r"DFT EXCHANGE \+ CORRELATION ENERGY" + r"\s+=\s*" + NUMBER + r"\s*$", outtext, re.MULTILINE
        )
        if mobj:
            logger.debug("matched dft xc")
            qcvar["DFT XC ENERGY"] = mobj.group(1)

        # Process Geometry
        mobj = re.search(
            # fmt: off
            r'^\s+' + r'ATOM      ATOMIC                      COORDINATES \(BOHR\)' + r'\s*' +
            r'^\s+' + r'CHARGE         X                   Y                   Z'+ r'\s*' +
            r'((?:\s+([A-Za-z]\w*)\s+\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s*\n)+)' +
            r"\s*$",
            # fmt: on
            outtext,
            re.MULTILINE | re.IGNORECASE,
        )
        if mobj:
            logger.debug("matched geom")
            molxyz = "%d bohr\n\n" % len(mobj.group(1).splitlines())
            for line in mobj.group(1).splitlines():
                lline = line.split()
                chg_on_center = int(float(lline[-4]))
                if chg_on_center > 0:
                    tag = f"{lline[-5]}"
                else:
                    tag = f"@{lline[-5].strip()}"
                molxyz += "%s %16s %16s %16s\n" % (tag, lline[-3], lline[-2], lline[-1])
            qcvar_coord = Molecule(
                validate=False,
                **qcel.molparse.to_schema(
                    qcel.molparse.from_string(molxyz, dtype="xyz+", fix_com=True, fix_orientation=True)["qm"], dtype=2
                ),
            )

        # Process Gradient
        mobj = re.search(
            # fmt: off
            r'^\s+' + r'GRADIENT OF THE ENERGY' + r'\s*'+
            r'^\s+' + r'----------------------' + r'\s*'+
            r'\s+' + r'\n'+
            r'^\s+' + r'UNITS ARE HARTREE/BOHR    E\'X               E\'Y               E\'Z' + r'\s*' +
            r'((?:\s+([1-9][0-9]*)+\s+([A-Za-z]\w*)+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s*\n)+)' +
            r'\s*$',
            # fmt: off
                outtext, re.MULTILINE
        )
        if mobj:
            logger.debug("matched gradient - after")
            atoms = []
            qcvar_grad = []
            for line in mobj.group(1).splitlines():
                lline = line.split()
                if lline == []:
                    pass
                else:
                    logger.debug("printing gradient")
                    atoms.append(lline[1])
                    qcvar_grad.append([float(lline[-3]), float(lline[-2]), float(lline[-1])])
            qcvar_grad = np.array(qcvar_grad).reshape((-1, 3))

        # Process SCF blocks
        mobj = re.search(
            # fmt: off
            r'^\s+' + r'PROPERTIES FOR THE .* DFT FUNCTIONAL .* DENSITY MATRIX' +
            r'(.*)' +
            r'\.\.\.\.\.\. END OF PROPERTY EVALUATION \.\.\.\.\.\.',
            # fmt: on
            outtext,
            re.MULTILINE | re.DOTALL,
        )
        if mobj:
            logger.debug("matched dft props")
            prop_block = mobj.group(1)
            mobj_energy = re.search(
                # fmt: off
                r'\s+ENERGY COMPONENTS\s*\n' +
                r'\s+-----------------\s*\n' +
                r'\s*\n' +
                r'\s+WAVEFUNCTION NORMALIZATION =.*\n' +
                r'\s*\n' +
                r'\s+ONE ELECTRON ENERGY =\s+' + NUMBER + r'\n' +
                r'\s+TWO ELECTRON ENERGY =\s+' + NUMBER + r'\n' +
                r'\s+NUCLEAR REPULSION ENERGY =\s+' + NUMBER + r'\n' +
                r'\s+------------------\s*\n' +
                r'\s+TOTAL ENERGY =\s+' + NUMBER+r'\n',
                # fmt: on
                prop_block
            )
            if mobj_energy:
                qcvar["ONE-ELECTRON ENERGY"] = mobj_energy.group(1)
                qcvar["TWO-ELECTRON ENERGY"] = mobj_energy.group(2)
                qcvar["SCF TOTAL ENERGY"] = mobj_energy.group(4)
                qcvar["DFT TOTAL ENERGY"] = mobj_energy.group(4)

            mobj_dipole = re.search(
                # fmt: off
                r'\s+ELECTROSTATIC MOMENTS\s*\n'
                r'\s+---------------------\s*\n'
                r'\s*\n'
                r'\s*POINT\s+1\s+X\s+Y\s+Z\s+\(BOHR\)\s+CHARGE\s*\n'
                r'.*\n'
                r'\s*DX\s+DY\s+DZ\s+/D/\s+\(DEBYE\)\s*\n'
                r'\s*' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*\n',
                # fmt: on
                prop_block
            )
            if mobj_dipole:
                d2au = Decimal(qcel.constants.conversion_factor("debye", "e * bohr"))
                qcvar["SCF DIPOLE"] = d2au * np.array(
                    [Decimal(mobj_dipole.group(1)), Decimal(mobj_dipole.group(2)), Decimal(mobj_dipole.group(3))]
                )

    # Process CURRENT Energies
    if "HF TOTAL ENERGY" in qcvar:
        qcvar["CURRENT REFERENCE ENERGY"] = qcvar["HF TOTAL ENERGY"]
        qcvar["CURRENT ENERGY"] = qcvar["HF TOTAL ENERGY"]

    if "MP2 TOTAL ENERGY" in qcvar and "MP2 CORRELATION ENERGY" in qcvar:
        qcvar["CURRENT CORRELATION ENERGY"] = qcvar["MP2 CORRELATION ENERGY"]
        qcvar["CURRENT ENERGY"] = qcvar["MP2 TOTAL ENERGY"]

    if "CISD TOTAL ENERGY" in qcvar and "CISD CORRELATION ENERGY" in qcvar:
        qcvar["CURRENT CORRELATION ENERGY"] = qcvar["CISD CORRELATION ENERGY"]
        qcvar["CURRENT ENERGY"] = qcvar["CISD TOTAL ENERGY"]

    if "LCCD TOTAL ENERGY" in qcvar and "LCCD CORRELATION ENERGY" in qcvar:
        qcvar["CURRENT CORRELATION ENERGY"] = qcvar["LCCD CORRELATION ENERGY"]
        qcvar["CURRENT ENERGY"] = qcvar["LCCD TOTAL ENERGY"]

    if "CCD TOTAL ENERGY" in qcvar and "CCD CORRELATION ENERGY" in qcvar:
        qcvar["CURRENT CORRELATION ENERGY"] = qcvar["CCD CORRELATION ENERGY"]
        qcvar["CURRENT ENERGY"] = qcvar["CCD TOTAL ENERGY"]

    if "CCSD TOTAL ENERGY" in qcvar and "CCSD CORRELATION ENERGY" in qcvar:
        qcvar["CURRENT CORRELATION ENERGY"] = qcvar["CCSD CORRELATION ENERGY"]
        qcvar["CURRENT ENERGY"] = qcvar["CCSD TOTAL ENERGY"]

    if "CR-CC(2,3) TOTAL ENERGY" in qcvar and "CR-CC(2,3) CORRELATION ENERGY" in qcvar:
        qcvar["CURRENT CORRELATION ENERGY"] = qcvar["CR-CC(2,3) CORRELATION ENERGY"]
        qcvar["CURRENT ENERGY"] = qcvar["CR-CC(2,3) TOTAL ENERGY"]

    if "CCSD+T(CCSD) TOTAL ENERGY" in qcvar and "CCSD+T(CCSD) CORRELATION ENERGY" in qcvar:
        qcvar["CURRENT CORRELATION ENERGY"] = qcvar["CCSD+T(CCSD) CORRELATION ENERGY"]
        qcvar["CURRENT ENERGY"] = qcvar["CCSD+T(CCSD) TOTAL ENERGY"]

    if "CCSD(T) TOTAL ENERGY" in qcvar and "CCSD(T) CORRELATION ENERGY" in qcvar:
        qcvar["CURRENT CORRELATION ENERGY"] = qcvar["CCSD(T) CORRELATION ENERGY"]
        qcvar["CURRENT ENERGY"] = qcvar["CCSD(T) TOTAL ENERGY"]

    if "DFT TOTAL ENERGY" in qcvar:
        qcvar["CURRENT REFERENCE ENERGY"] = qcvar["DFT TOTAL ENERGY"]
        qcvar["CURRENT ENERGY"] = qcvar["DFT TOTAL ENERGY"]

    if "FCI TOTAL ENERGY" in qcvar:  # and 'FCI CORRELATION ENERGY' in qcvar:
        qcvar["CURRENT ENERGY"] = qcvar["FCI TOTAL ENERGY"]

    qcvar[f"N ATOMS"] = len(qcvar_coord.symbols)

    return qcvar, qcvar_coord, qcvar_grad, module
