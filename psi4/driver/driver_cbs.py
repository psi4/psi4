#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2025 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#
"""Plan, run, and assemble QC tasks to obtain composite method, basis, & options treatments.

========
CBS Flow
========
Bullet points are major actions
Lines of dashes denote function calls
stage: scf, corl, delta1, delta2, ...
e/d/dd=dg/g/h := energy, dipole, dipole derivative = dipole gradient, gradient, Hessian

cbs_text_parser()
-----------------
* called from task_planner() only if "/" in method

    _parse_cbs_gufunc_string()
    --------------------------
    * break user string into paired method and basis stages

* transform user string into cbs kwargs inc'l basic cbs_metadata
  cbs kwargs may signal simple method/basis single point -or- a modelchem requiring CompositeComputer


----------------------------
CompositeComputer.__init__()
----------------------------

    _process_cbs_kwargs()
    ---------------------
    * if input is cbs_metadata dict, skip to _validate_cbs_inputs()
    * otherwise, transform user kwargs into trial cbs_metadata format (aka dict spec)

        _validate_cbs_inputs()
        ----------------------

            _get_default_xtpl()
            -------------------
            * supply default xtpl fn for stage and basis conditions

            _expand_bracketed_basis()
            -------------------------
            * parse and validate user bases

        * check and supply defaults for cbs_metadata format (various calls to above two fns)

* BaseComputer.__init__()

    _build_cbs_compute()
    --------------------

        _expand_scheme_orders()
        -----------------------
        * form f_fields dict of entries for each zeta in a scheme (single NEED; entries related by nonlinear fn
         (that is, constructing the CBS energy from the component energies is nonlinear))

        _contract_bracketed_basis()
        ---------------------------
        * form basis abbr. string from basis seq

    * form d_fields list of stages or stage halves from NEEDs (GRAND_NEED; items related linearly to form final val)
    * form list of entries (entry:= mtd-bas-opt specification) mentioned in GRAND_NEED (MODELCHEM; redundant, naive)
    * form subset of MODELCHEM with minimal list of jobs (job:= entry on which to call QC) to satisfy CBS (JOBS; minimal, enlightened)
    * form superset of JOBS with maximal list of entries resulting from JOBS (TROVE)
    * return GRAND_NEED/cbsrec, JOBS/compute_list, TROVE/trove

* form task_list of AtomicComputers 1:1 from JOBS/compute_list

-------------------------------
CompositeComputer.build_tasks()
-------------------------------
* pass

---------------------------
CompositeComputer.compute()
---------------------------
* compute() for each job in task list

-----------------------------------
CompositeComputer.get_psi_results()
-----------------------------------

    Computer.get_results()
    ----------------------

        Computer._prepare_results()
        ---------------------------
        * get_results() for each job in task list
        * arrange atomicresult data into e/d/g/h fields in compute_list and copy them into cbs tables

            _assemble_cbs_components()
            --------------------------
            * fill in results from TROVE/trove into GRAND_NEED/cbsrec

                _contract_scheme_orders()
                -------------------------
                * prepare arguments for xtpl fns based on desired E/D/G/H quantity

            * form extrapolated values for all available E/D/G/H quantities
            * return structure of extrapolated values and filled-in GRAND_NEED/cbsrec

            _summary_table()
            ----------------
            * build string table of cbs results

    * form cbs qcvars, inc'l number, E, DG, G, H as available
    * form model, including detailed dict at atomicresult.extras["cbs_record"]

* convert result to psi4.core.Matrix (non-energy)

    _cbs_schema_to_wfn()
    --------------------
    * build wfn from cbs mol and basis (always def2-svp) and module (if present)
    * push qcvars to P::e and wfn

* return e/g/h and wfn

"""

import copy
import logging
import re
import sys
from typing import TYPE_CHECKING, Any, Callable, Dict, List, Optional, Tuple, Union

import numpy as np

try:
    from pydantic.v1 import Field, validator
except ImportError:
    from pydantic import Field, validator

from qcelemental.models import AtomicResult, DriverEnum

from psi4 import core

from . import driver_util, p4util, qcdb
from .constants import pp
from .driver_cbs_helper import (  # lgtm[py/unused-import]
    composite_procedures,
    register_composite_function,
    register_xtpl_function,
    xtpl_procedures,
)
from .driver_util import UpgradeHelper
from .p4util.exceptions import ValidationError
from .procrouting.interface_cfour import cfour_psivar_list
from .task_base import AtomicComputer, BaseComputer, EnergyGradientHessianWfnReturn

if TYPE_CHECKING:
    import qcportal

logger = logging.getLogger(__name__)

zeta_values = 'dtq5678'
_zeta_val2sym = {k + 2: v for k, v in enumerate(zeta_values)}
_zeta_sym2val = {v: k for k, v in _zeta_val2sym.items()}
_addlremark = {'energy': '', 'gradient': ', GRADIENT', 'hessian': ', HESSIAN'}
_f_fields = ['f_wfn', 'f_basis', 'f_zeta', 'f_options', 'f_energy', 'f_gradient', 'f_hessian', 'f_dipole', 'f_dipder']
_lmh_labels = {
    1: ['HI'],
    2: ['LO', 'HI'],
    3: ['LO', 'MD', 'HI'],
    4: ['LO', 'MD', 'M2', 'HI'],
    5: ['LO', 'MD', 'M2', 'M3', 'HI']
}
CBSMetadata = List[Dict[str, Any]]


# remove in 1.8
# these get input files to the point where they raise an UpgradeHelper
def xtpl_highest_1():
    pass


def scf_xtpl_helgaker_2():
    pass


def scf_xtpl_truhlar_2():
    pass


def scf_xtpl_karton_2():
    pass


def scf_xtpl_helgaker_3():
    pass


def corl_xtpl_helgaker_2():
    pass


def _expand_bracketed_basis(basisstring: str, molecule: Union["qcdb.Molecule", core.Molecule] = None) -> Tuple[List[str], List[int]]:
    """Function to transform and validate basis series specification for cbs().

    Parameters
    ----------
    basisstring
        A string containing the basis sets to be expanded.
        A basis set with no paired square brackets is passed through
        with zeta level 0 (e.g., ``'6-31+G(d,p)'`` is returned as
        ``(["6-31+G(d,p)"], [0])``). A basis set with square brackets is checked
        for sensible sequence and returned as separate basis sets
        (e.g., ``'cc-pV[Q5]Z'` is returned as ``(["cc-pVQZ", "cc-pV5Z"], [4, 5])``).
        Allows out-of-order zeta specification (e.g., ``[qtd]``) and numeral for
        number (e.g., ``[23]``). Does not allow skipped zetas (e.g., ``[dq]``), zetas
        outside the [2,8] range, non-Dunning, non-Ahlrichs, or non-Jensen sets,
        or non-findable .gbs sets.
    molecule
        This function checks that the basis is valid by trying to build
        the qcdb.BasisSet object for *molecule* or for H2 if None.

    Returns
    -------
    tuple
        Tuple in the ``([basis set names], [basis set zetas])`` format.

    """
    BSET = []
    ZSET = []
    legit_compound_basis = re.compile(
        r'^(?P<pre>.*cc-.*|def2-|.*pcs+eg-|.*)\[(?P<zeta>[dtq2345678,s1]*)\](?P<post>.*z.*|)$', re.IGNORECASE)
    pc_basis = re.compile(r'.*pcs+eg-$', re.IGNORECASE)
    def2_basis = re.compile(r'def2-', re.IGNORECASE)
    zapa_basis = re.compile(r'.*zapa.*', re.IGNORECASE)

    if legit_compound_basis.match(basisstring):
        basisname = legit_compound_basis.match(basisstring)
        # handle def2-svp* basis sets as double-zeta
        if def2_basis.match(basisname.group('pre')):
            bn_gz = basisname.group('zeta').replace("s", "d")
        # handle pc-n basis set polarisation -> zeta conversion
        elif pc_basis.match(basisname.group('pre')):
            bn_gz = basisname.group('zeta').replace("4", "5").replace("3", "4").replace("2", "3").replace("1", "2")
        else:
            bn_gz = basisname.group('zeta')
        # filter out commas and be forgiving of e.g., t5q or 3q
        zetas = [z for z in zeta_values if (z in bn_gz or str(zeta_values.index(z) + 2) in bn_gz)]
        for b in zetas:
            if ZSET and (int(ZSET[len(ZSET) - 1]) - zeta_values.index(b)) != 1:
                raise ValidationError("""Basis set '%s' has skipped zeta level '%s'.""" %
                                      (basisstring, _zeta_val2sym[_zeta_sym2val[b] - 1]))
            # reassemble def2-svp* properly instead of def2-dzvp*
            if def2_basis.match(basisname.group('pre')) and b == "d":
                BSET.append(basisname.group('pre') + "s" + basisname.group('post')[1:])
            # reassemble pc-n basis sets properly
            elif pc_basis.match(basisname.group('pre')):
                BSET.append(basisname.group('pre') + "{0:d}".format(_zeta_sym2val[b] - 1))
            # assemble nZaPa basis sets
            elif zapa_basis.match(basisname.group('post')):
                bzapa = b.replace("d", "2").replace("t", "3").replace("q", "4")
                BSET.append(basisname.group('pre') + bzapa + basisname.group('post'))
            else:
                BSET.append(basisname.group('pre') + b + basisname.group('post'))
            ZSET.append(zeta_values.index(b) + 2)
    elif re.match(r'.*\[.*\].*$', basisstring, flags=re.IGNORECASE):
        raise ValidationError(
            """Basis series '%s' invalid. Specify a basis series matching"""
            """ '*cc-*[dtq2345678,]*z*'. or 'def2-[sdtq]zvp*' or '*pcs[s]eg-[1234]' or '[1234567]ZaPa' """ %
            (basisstring))
    else:
        BSET.append(basisstring)
        ZSET.append(0)

    if molecule is None:
        molecule = """\nH\nH 1 1.00\n"""
    elif isinstance(molecule, core.Molecule):
        molecule = qcdb.Molecule(molecule.to_dict())

    for basis in BSET:
        try:
            qcdb.BasisSet.pyconstruct(molecule, "BASIS", basis)
        except qcdb.BasisSetNotFound:
            e = sys.exc_info()[1]
            raise ValidationError(f"""Basis set '{basis}' not available for molecule.""")

    return (BSET, ZSET)


def _contract_bracketed_basis(basisarray: List[str]) -> str:
    """Function to re-form a bracketed basis set string from a sequential series
    of basis sets. Essentially the inverse of _expand_bracketed_basis(). Used to
    print a nicely formatted basis set string in the results table.

    Parameters
    ----------
    basisarray
        Basis set names, differing by zeta level, e.g. ``["cc-pvqz", "cc-pv5z"]``.

    Returns
    -------
    str
        A nicely formatted basis set string, e.g. ``"cc-pv[q5]z"`` for the above example.

    """

    if len(basisarray) == 1:
        return basisarray[0]

    else:
        zetaindx = [i for i in range(len(basisarray[0])) if basisarray[0][i] != basisarray[1][i]][0]
        ZSET = [bas[zetaindx] for bas in basisarray]
        pre = basisarray[1][:zetaindx]
        post = basisarray[1][zetaindx + 1:]

        return "".join([pre, "[", *ZSET, "]", post])


def return_energy_components():
    """Define some quantum chemical knowledge, namely what methods are subsumed in others."""

    # yapf: disable
    VARH = {}
    VARH['scf'] = {
                            'scf': 'SCF TOTAL ENERGY'}
    VARH['hf'] = {
                             'hf': 'HF TOTAL ENERGY'}
    VARH['mp2'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY'}
    VARH['dlpno-mp2'] = {
                             'hf': 'HF TOTAL ENERGY',
                      'dlpno-mp2': 'MP2 TOTAL ENERGY'}
    VARH['mp2d'] = {
                             'hf': 'HF TOTAL ENERGY',
                           'mp2d': 'MP2D TOTAL ENERGY'}
    VARH['mp2.5'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                          'mp2.5': 'MP2.5 TOTAL ENERGY',
                            'mp3': 'MP3 TOTAL ENERGY'}
    VARH['mp3'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                          'mp2.5': 'MP2.5 TOTAL ENERGY',
                            'mp3': 'MP3 TOTAL ENERGY'}
    VARH['mp4(sdq)'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                          'mp2.5': 'MP2.5 TOTAL ENERGY',
                            'mp3': 'MP3 TOTAL ENERGY',
                       'mp4(sdq)': 'MP4(SDQ) TOTAL ENERGY'}
    VARH['mp4'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                          'mp2.5': 'MP2.5 TOTAL ENERGY',
                            'mp3': 'MP3 TOTAL ENERGY',
                       'mp4(sdq)': 'MP4(SDQ) TOTAL ENERGY',
                            'mp4': 'MP4 TOTAL ENERGY'}
    VARH['omp2'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                           'omp2': 'OMP2 TOTAL ENERGY'}
    VARH['omp2.5'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                          'mp2.5': 'MP2.5 TOTAL ENERGY',
                         'omp2.5': 'OMP2.5 TOTAL ENERGY'}
    VARH['omp3'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                            'mp3': 'MP3 TOTAL ENERGY',
                           'omp3': 'OMP3 TOTAL ENERGY'}
    VARH['olccd'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                          'olccd': 'OLCCD TOTAL ENERGY'}
    VARH['oremp2'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                         'oremp2': 'OREMP2 TOTAL ENERGY'}
    VARH['lccd'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                           'lccd': 'LCCD TOTAL ENERGY'}
    VARH['remp2'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                          'remp2': 'REMP2 TOTAL ENERGY'}
    VARH['lccsd'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                          'lccsd': 'LCCSD TOTAL ENERGY'}
    VARH['cepa(0)'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                        'cepa(0)': 'CEPA(0) TOTAL ENERGY'}
    VARH['cepa(1)'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                        'cepa(1)': 'CEPA(1) TOTAL ENERGY'}
    VARH['cepa(3)'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                        'cepa(3)': 'CEPA(3) TOTAL ENERGY'}
    VARH['acpf'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                           'acpf': 'ACPF TOTAL ENERGY'}
    VARH['aqcc'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                           'aqcc': 'AQCC TOTAL ENERGY'}
    VARH['qcisd'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                          'mp2.5': 'MP2.5 TOTAL ENERGY',
                            'mp3': 'MP3 TOTAL ENERGY',
                       'mp4(sdq)': 'MP4(SDQ) TOTAL ENERGY',
                          'qcisd': 'QCISD TOTAL ENERGY'}
    VARH['cc2'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                            'cc2': 'CC2 TOTAL ENERGY'}
    VARH['ccsd'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                           'ccsd': 'CCSD TOTAL ENERGY'}
    VARH['bccd'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                           'bccd': 'CCSD TOTAL ENERGY'}
    VARH['cc3'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                            'cc3': 'CC3 TOTAL ENERGY'}
    VARH['fno-ccsd'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                       'fno-ccsd': 'CCSD TOTAL ENERGY'}
    VARH['fno-ccsd(t)'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                       'fno-ccsd': 'CCSD TOTAL ENERGY',
                    'fno-ccsd(t)': 'CCSD(T) TOTAL ENERGY'}
    VARH['dlpno-ccsd'] = {
                             'hf': 'HF TOTAL ENERGY',
                      'dlpno-mp2': 'MP2 TOTAL ENERGY',
                     'dlpno-ccsd': 'CCSD TOTAL ENERGY'}
    VARH['dlpno-ccsd(t0)'] = {
                             'hf': 'HF TOTAL ENERGY',
                      'dlpno-mp2': 'MP2 TOTAL ENERGY',
                     'dlpno-ccsd': 'CCSD TOTAL ENERGY',
                 'dlpno-ccsd(t0)': 'CCSD(T) TOTAL ENERGY'}
    VARH['dlpno-ccsd(t)'] = {
                             'hf': 'HF TOTAL ENERGY',
                      'dlpno-mp2': 'MP2 TOTAL ENERGY',
                     'dlpno-ccsd': 'CCSD TOTAL ENERGY',
                  'dlpno-ccsd(t)': 'CCSD(T) TOTAL ENERGY'}
    VARH['qcisd(t)'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                          'mp2.5': 'MP2.5 TOTAL ENERGY',
                            'mp3': 'MP3 TOTAL ENERGY',
                       'mp4(sdq)': 'MP4(SDQ) TOTAL ENERGY',
                          'qcisd': 'QCISD TOTAL ENERGY',
                       'qcisd(t)': 'QCISD(T) TOTAL ENERGY'}
    VARH['ccsd(t)'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                           'ccsd': 'CCSD TOTAL ENERGY',
                        'ccsd(t)': 'CCSD(T) TOTAL ENERGY'}
    VARH['ccsd(at)'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                           'ccsd': 'CCSD TOTAL ENERGY',
                       'ccsd(at)': 'CCSD(AT) TOTAL ENERGY'}
    VARH["a-ccsd(t)"] = VARH["ccsd(at)"]
    VARH['bccd(t)'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                           'ccsd': 'CCSD TOTAL ENERGY',
                        'bccd(t)': 'CCSD(T) TOTAL ENERGY'}
    VARH['cisd'] = {
                             'hf': 'HF TOTAL ENERGY',
                           'cisd': 'CISD TOTAL ENERGY'}
    VARH['cisdt'] = {
                             'hf': 'HF TOTAL ENERGY',
                          'cisdt': 'CISDT TOTAL ENERGY'}
    VARH['cisdtq'] = {
                             'hf': 'HF TOTAL ENERGY',
                         'cisdtq': 'CISDTQ TOTAL ENERGY'}
    VARH['fci'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'fci': 'FCI TOTAL ENERGY'}
    VARH['ccsdt'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                          'ccsdt': 'CCSDT TOTAL ENERGY'}
    VARH['ccsdt(q)'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                          'ccsdt': 'CCSDT TOTAL ENERGY',
                       'ccsdt(q)': 'CCSDT(Q) TOTAL ENERGY'}

    for cilevel in range(2, 99):
        VARH[f'ci{cilevel}'] = {
                             'hf': 'HF TOTAL ENERGY',
                   f'ci{cilevel}': 'CI TOTAL ENERGY'}

    for mplevel in range(5, 99):
        VARH[f'mp{mplevel}'] = {
                             'hf': 'HF TOTAL ENERGY',
                   f'mp{mplevel}': f'MP{mplevel} TOTAL ENERGY'}
        for mplevel2 in range(2, mplevel):
            VARH[f'mp{mplevel}'][f'mp{mplevel2}'] = f'MP{mplevel2} TOTAL ENERGY'

    # Integrate CFOUR methods
    VARH.update(cfour_psivar_list())
    return VARH
    # yapf: enable


VARH = return_energy_components()


def _get_default_xtpl(nbasis: int, xtpl_type: str) -> Callable:
    """ A helper function to determine default extrapolation type.

    Parameters
    ----------
    nbasis
        Number of basis sets
    xtpl_type
        {'scf', 'corl'}
        Extrapolation type: 'scf' for the total energy, 'corl' for just the
        correlation component.

    Returns
    -------
    Callable
        Extrapolation function to be used.
    """

    if nbasis == 1 and xtpl_type in ["scf", "corl"]:
        return "xtpl_highest_1"
    elif xtpl_type == "scf":
        if nbasis == 2:
            return "scf_xtpl_helgaker_2"
        elif nbasis == 3:
            return "scf_xtpl_helgaker_3"
        else:
            raise ValidationError(f"Wrong number of basis sets supplied to scf_xtpl: {nbasis}")
    elif xtpl_type == "corl":
        if nbasis == 2:
            return "corl_xtpl_helgaker_2"
        else:
            raise ValidationError(f"Wrong number of basis sets supplied to corl_xtpl: {nbasis}")
    else:
        raise ValidationError(f"Stage treatment must be 'corl' or 'scf', not '{xtpl_type}'")


def _validate_cbs_inputs(cbs_metadata: CBSMetadata, molecule: Union["qcdb.Molecule", core.Molecule]) -> CBSMetadata:
    """ A helper function which validates the ``cbs_metadata`` format,
    expands basis sets, and provides sensible defaults for optional arguments.

    Parameters
    ----------
    cbs_metadata
        List of dicts containing CBS stage keywords.
    molecule
        Molecule to be passed to _expand_bracketed_basis()

    Returns
    -------
    list
        Validated list of dictionaries, with each item consisting of an extrapolation
        stage. All validation takes place here.
    """

    # TODO: split options into mixable (qc_module=ccenergy/"") or non-mixable (freeze_core=true/false)

    metadata = []
    for iitem, item in enumerate(cbs_metadata):
        # 1a) all items must have wfn
        if "wfn" not in item:
            raise ValidationError(f"Stage {iitem} doesn't have defined level of theory!")
    # 1b) all items must have basis set
        if "basis" not in item:
            raise ValidationError(f"Stage {iitem} doesn't have defined basis sets!")
    # 2a) process required stage parameters and assign defaults
        stage = {}
        stage["wfn"] = item["wfn"].lower()
        stage["basis"] = _expand_bracketed_basis(item["basis"].lower(), molecule)
        # 2b) if first item is not HF, generate it
        if len(metadata) == 0 and stage["wfn"] not in ["hf", "c4-hf", "scf", "c4-scf"]:
            scf = {}
            if stage["wfn"].startswith("c4"):
                scf["wfn"] = "c4-hf"
            else:
                scf["wfn"] = "hf"
            scf["basis"] = ([stage["basis"][0][-1]], [stage["basis"][1][-1]])
            scf["treatment"] = "scf"
            scf["stage"] = "scf"
            scf["scheme"] = _get_default_xtpl(len(scf["basis"][1]), scf["treatment"])
            scf["alpha"] = None
            scf["options"] = False
            scf["options_lo"] = False
            metadata.append(scf)
    # 2c) keep processing current stage
        stage["treatment"] = item.get("treatment", "scf" if len(metadata) == 0 else "corl")
        stage["stage"] = item.get("stage", False)
        if not stage["stage"]:
            if len(metadata) == 0:
                stage["stage"] = "scf"
            elif len(metadata) == 1:
                stage["stage"] = "corl"
            else:
                stage["stage"] = f"delta{len(metadata) - 1}"
        stage["scheme"] = item.get("scheme", _get_default_xtpl(len(stage["basis"][1]), stage["treatment"]))
        if len(metadata) > 0:
            stage["wfn_lo"] = item.get("wfn_lo", metadata[-1].get("wfn")).lower()
            stage["basis_lo"] = _expand_bracketed_basis(item.get("basis_lo", item["basis"]).lower(), molecule)
            if len(stage["basis"][0]) != len(stage["basis_lo"][0]):
                raise ValidationError("""Number of basis sets inconsistent
                                            between high ({}) and low ({}) levels.""".format(
                    len(stage["basis"][0]), len(stage["basis_lo"][0])))
        stage["alpha"] = item.get("alpha", None)
        stage["options"] = item.get("options", False)
        stage["options_lo"] = item.get("options_lo", False)
        metadata.append(stage)
    return metadata


def _process_cbs_kwargs(kwargs: Dict) -> CBSMetadata:
    """ A helper function which translates supplied kwargs into the
    ``cbs_metadata`` format and passes it for validation.

    Parameters
    ----------
    kwargs
        kwargs containing the CBS function specification.

    Returns
    -------
    cbs_metadata
        List of dictionaries, with each item consisting of an extrapolation
        stage. All validation takes place here.
    """

    molecule = kwargs.get('molecule', core.get_active_molecule())

    if "cbs_metadata" in kwargs:
        # if we passed in a dict, validate it right away
        cbs_metadata = kwargs["cbs_metadata"]
    else:
        # if we passed in options, check for consecutive correlations first
        if "delta_wfn" in kwargs and "corl_wfn" not in kwargs:
            raise ValidationError("Delta function supplied without corl_wfn defined.")
        if "delta2_wfn" in kwargs and "delta_wfn" not in kwargs:
            raise ValidationError("Second delta function supplied without delta_wfn defined.")
        cbs_metadata = []
        possible_stages = ["scf", "corl"]
        while len(possible_stages) > 0:
            sn = possible_stages.pop(0)
            if f"{sn}_wfn" in kwargs and f"{sn}_basis" in kwargs:
                # either both *_wfn and *_basis have to be specified
                stage = {"wfn": kwargs[f"{sn}_wfn"], "basis": kwargs[f"{sn}_basis"]}
            elif sn == "scf" and f"{sn}_basis" in kwargs:
                # or we're at a scf stage which can be implied with a provided scf_basis
                stage = {"wfn": "hf", "basis": kwargs[f"{sn}_basis"]}
            else:
                # otherwise go to the next possible stage
                continue
            # if we made it here, stage exists - parse other keywords
            if f"{sn}_scheme" in kwargs:
                stage["scheme"] = kwargs[f"{sn}_scheme"]
            if f"{sn}_wfn_lesser" in kwargs:
                stage["wfn_lo"] = kwargs[f"{sn}_wfn_lesser"]
            if f"cbs_{sn}_alpha" in kwargs:
                stage["alpha"] = kwargs[f"cbs_{sn}_alpha"]
            elif f"{sn}_alpha" in kwargs:
                stage["alpha"] = kwargs[f"{sn}_alpha"]
            cbs_metadata.append(stage)
            if sn == "corl":
                possible_stages.append("delta")
            elif sn == "delta":
                possible_stages.append("delta2")

    return _validate_cbs_inputs(cbs_metadata, molecule)


###################################
##  Start of Complete Basis Set  ##
###################################


def cbs(func, label, **kwargs):
    r"""Function to define a multistage energy method from combinations of
    basis set extrapolations and delta corrections and condense the
    components into a minimum number of calculations.

    :aliases: complete_basis_set()

    :returns: (*float*) -- Total electronic energy in Hartrees

    :PSI variables:

    .. hlist::
       :columns: 1

       * :psivar:`CBS TOTAL ENERGY`
       * :psivar:`CBS REFERENCE ENERGY`
       * :psivar:`CBS CORRELATION ENERGY`
       * :psivar:`CURRENT ENERGY`
       * :psivar:`CURRENT REFERENCE ENERGY`
       * :psivar:`CURRENT CORRELATION ENERGY`

    .. caution:: Some features are not yet implemented. Buy a developer a coffee.

       - No way to tell function to boost fitting basis size for all calculations.

       - Need to add more extrapolation schemes

    As represented in the equation below, a CBS energy method is defined in several
    sequential stages (scf, corl, delta1, delta2, ... ) covering treatment
    of the reference total energy, the correlation energy, a delta correction to the
    correlation energy, and a second delta correction, etc.. Each is activated by its
    stage_wfn keyword, or as a field in the ```cbs_metadata``` list, and is only
    allowed if all preceding stages are active.

    .. include:: /cbs_eqn.rst

    * Energy Methods
        The presence of a stage_wfn keyword is the indicator to incorporate
        (and check for stage_basis and stage_scheme keywords) and compute
        that stage in defining the CBS energy.

        The cbs() function requires, at a minimum, ``name='scf'`` and ``scf_basis``
        keywords to be specified for reference-step only jobs and ``name`` and
        ``corl_basis`` keywords for correlated jobs.

        The following energy methods have been set up for cbs().

        .. hlist::
           :columns: 5

           * scf
           * hf
           * mp2
           * mp2.5
           * mp3
           * mp4(sdq)
           * mp4
           * mp\ *n*
           * omp2
           * omp2.5
           * omp3
           * olccd
           * lccd
           * lccsd
           * cepa(0)
           * cepa(1)
           * cepa(3)
           * acpf
           * aqcc
           * qcisd
           * cc2
           * ccsd
           * fno-ccsd
           * bccd
           * cc3
           * qcisd(t)
           * ccsd(t)
           * fno-ccsd(t)
           * bccd(t)
           * cisd
           * cisdt
           * cisdtq
           * ci\ *n*
           * fci
           * mrccsd
           * mrccsd(t)
           * mrccsdt
           * mrccsdt(q)

    :type name: str
    :param name: ``'scf'`` || ``'ccsd'`` || etc.

        First argument, usually unlabeled. Indicates the computational method
        for the correlation energy, unless only reference step to be performed,
        in which case should be ``'scf'``. Overruled if stage_wfn keywords supplied.

    :type scf_wfn: str
    :param scf_wfn: |dl| ``'scf'`` |dr| || ``'c4-scf'`` || etc.

        Indicates the energy method for which the reference energy is to be
        obtained. Generally unnecessary, as 'scf' is *the* scf in |PSIfour| but
        can be used to direct lone scf components to run in |PSIfour| or Cfour
        in a mixed-program composite method.

    :type corl_wfn: str
    :param corl_wfn: ``'mp2'`` || ``'ccsd(t)'`` || etc.

        Indicates the energy method for which the correlation energy is to be
        obtained. Can also be specified with ``name`` or as the unlabeled
        first argument to the function.

    :type delta_wfn: str
    :param delta_wfn: ``'ccsd'`` || ``'ccsd(t)'`` || etc.

        Indicates the (superior) energy method for which a delta correction
        to the correlation energy is to be obtained.

    :type delta_wfn_lesser: str
    :param delta_wfn_lesser: |dl| ``corl_wfn`` |dr| || ``'mp2'`` || etc.

        Indicates the inferior energy method for which a delta correction
        to the correlation energy is to be obtained.

    :type delta2_wfn: str
    :param delta2_wfn: ``'ccsd'`` || ``'ccsd(t)'`` || etc.

        Indicates the (superior) energy method for which a second delta correction
        to the correlation energy is to be obtained.

    :type delta2_wfn_lesser: str
    :param delta2_wfn_lesser: |dl| ``delta_wfn`` |dr| || ``'ccsd(t)'`` || etc.

        Indicates the inferior energy method for which a second delta correction
        to the correlation energy is to be obtained.

    * Basis Sets
        Currently, the basis set set through ``set`` commands have no influence
        on a cbs calculation.

    :type scf_basis: :ref:`basis string <apdx:basisElement>`
    :param scf_basis: |dl| ``corl_basis`` |dr| || ``'cc-pV[TQ]Z'`` || ``'jun-cc-pv[tq5]z'`` || ``'6-31G*'`` || etc.

        Indicates the sequence of basis sets employed for the reference energy.
        If any correlation method is specified, ``scf_basis`` can default
        to ``corl_basis``.

    :type corl_basis: :ref:`basis string <apdx:basisElement>`
    :param corl_basis: ``'cc-pV[TQ]Z'`` || ``'jun-cc-pv[tq5]z'`` || ``'6-31G*'`` || etc.

        Indicates the sequence of basis sets employed for the correlation energy.

    :type delta_basis: :ref:`basis string <apdx:basisElement>`
    :param delta_basis: ``'cc-pV[TQ]Z'`` || ``'jun-cc-pv[tq5]z'`` || ``'6-31G*'`` || etc.

        Indicates the sequence of basis sets employed for the delta correction
        to the correlation energy.

    :type delta2_basis: :ref:`basis string <apdx:basisElement>`
    :param delta2_basis: ``'cc-pV[TQ]Z'`` || ``'jun-cc-pv[tq5]z'`` || ``'6-31G*'`` || etc.

        Indicates the sequence of basis sets employed for the second delta correction
        to the correlation energy.

    * Schemes
        Transformations of the energy through basis set extrapolation for each
        stage of the CBS definition. A complaint is generated if number of basis
        sets in stage_basis does not exactly satisfy requirements of stage_scheme.
        An exception is the default, ``'xtpl_highest_1'``, which uses the best basis
        set available. See :ref:`sec:cbs_xtpl` for all available schemes.

    :type scf_scheme: str
    :param scf_scheme: |dl| ``'xtpl_highest_1'`` |dr| || ``'scf_xtpl_helgaker_3'`` || etc.

        Indicates the basis set extrapolation scheme to be applied to the reference energy.
        Defaults to :py:func:`~psi4.driver.driver_cbs_helper.scf_xtpl_helgaker_3` if three valid basis sets
        present in ``psi4.driver.driver_cbs.scf_basis``, :py:func:`~psi4.driver.driver_cbs_helper.scf_xtpl_helgaker_2` if two valid basis
        sets present in ``scf_basis``, and :py:func:`~psi4.driver.driver_cbs_helper.xtpl_highest_1` otherwise.

        .. hlist::
           :columns: 1

           * :py:func:`~psi4.driver.driver_cbs_helper.xtpl_highest_1`
           * :py:func:`~psi4.driver.driver_cbs_helper.scf_xtpl_helgaker_3`
           * :py:func:`~psi4.driver.driver_cbs_helper.scf_xtpl_helgaker_2`
           * :py:func:`~psi4.driver.driver_cbs_helper.scf_xtpl_truhlar_2`
           * :py:func:`~psi4.driver.driver_cbs_helper.scf_xtpl_karton_2`

    :type corl_scheme: str
    :param corl_scheme: |dl| ``'xtpl_highest_1'`` |dr| || ``'corl_xtpl_helgaker_2'`` || etc.

        Indicates the basis set extrapolation scheme to be applied to the correlation energy.
        Defaults to :py:func:`~psi4.driver.driver_cbs_helper.corl_xtpl_helgaker_2` if two valid basis sets
        present in ``corl_basis`` and :py:func:`~psi4.driver.driver_cbs_helper.xtpl_highest_1` otherwise.

        .. hlist::
           :columns: 1

           * :py:func:`~psi4.driver.driver_cbs_helper.xtpl_highest_1`
           * :py:func:`~psi4.driver.driver_cbs_helper.corl_xtpl_helgaker_2`

    :type delta_scheme: str
    :param delta_scheme: |dl| ``'xtpl_highest_1'`` |dr| || ``'corl_xtpl_helgaker_2'`` || etc.

        Indicates the basis set extrapolation scheme to be applied to the delta correction
        to the correlation energy.
        Defaults to :py:func:`~psi4.driver.driver_cbs_helper.corl_xtpl_helgaker_2` if two valid basis sets
        present in ``delta_basis`` and :py:func:`~psi4.driver.driver_cbs_helper.xtpl_highest_1` otherwise.

        .. hlist::
           :columns: 1

           * :py:func:`~psi4.driver.driver_cbs_helper.xtpl_highest_1`
           * :py:func:`~psi4.driver.driver_cbs_helper.corl_xtpl_helgaker_2`

    :type delta2_scheme: str
    :param delta2_scheme: |dl| ``'xtpl_highest_1'`` |dr| || ``'corl_xtpl_helgaker_2'`` || etc.

        Indicates the basis set extrapolation scheme to be applied to the second delta correction
        to the correlation energy.
        Defaults to :py:func:`~psi4.driver.driver_cbs_helper.corl_xtpl_helgaker_2` if two valid basis sets
        present in ``delta2_basis`` and :py:func:`~psi4.driver.driver_cbs_helper.xtpl_highest_1` otherwise.

        .. hlist::
           :columns: 1

           * :py:func:`~psi4.driver.driver_cbs_helper.xtpl_highest_1`
           * :py:func:`~psi4.driver.driver_cbs_helper.corl_xtpl_helgaker_2`

    :type scf_alpha: float
    :param scf_alpha: |dl| ``1.63`` |dr|

        Overrides the default \alpha parameter used in the listed SCF extrapolation procedures.
        Has no effect on others, including :py:func:`~psi4.driver.driver_cbs_helper.xtpl_highest_1` and :py:func:`~psi4.driver.driver_cbs_helper.scf_xtpl_helgaker_3`.

        .. hlist::
           :columns: 1

           * :py:func:`~psi4.driver.driver_cbs_helper.scf_xtpl_helgaker_2`
           * :py:func:`~psi4.driver.driver_cbs_helper.scf_xtpl_truhlar_2`
           * :py:func:`~psi4.driver.driver_cbs_helper.scf_xtpl_karton_2`

    :type corl_alpha: float
    :param corl_alpha: |dl| ``3.00`` |dr| 

        Overrides the default \alpha parameter used in the listed :py:func:`~psi4.driver.driver_cbs_helper.corl_xtpl_helgaker_2` correlation
        extrapolation to the corl stage. The supplied \alpha does not impact delta or any further stages.

        .. hlist::
           :columns: 1

           * :py:func:`~psi4.driver.driver_cbs_helper.corl_xtpl_helgaker_2`

    :type delta_alpha: float
    :param delta_alpha: |dl| ``3.00`` |dr| 

        Overrides the default \alpha parameter used in the listed
        :py:func:`~psi4.driver.driver_cbs_helper.corl_xtpl_helgaker_2` correlation extrapolation for the delta correction. Useful when
        delta correction is performed using smaller basis sets for which a different \alpha might
        be more appropriate.

        .. hlist::
           :columns: 1

           * :py:func:`~psi4.driver.driver_cbs_helper.corl_xtpl_helgaker_2`

    * Combined interface
    
    :type cbs_metadata: List[Dict]
    :param cbs_metadata: |dl| autogenerated from above keywords |dr| || ``[{"wfn": "hf", "basis": "cc-pv[TQ5]z"}]`` || etc.

        This is the interface to which all of the above calls are internally translated. The first item in
        the array is always defining the SCF contribution to the total energy. The required items in the
        dictionary are:

        * ```wfn```: typically ```HF```, which is subsumed in correlated methods anyway.
        * ```basis```: basis set, can be in a bracketed form (eg. ```cc-pv[tq]z```)

        |  Other supported arguments for the first dictionary are:

        * ```scheme```: scf extrapolation scheme function, by default it is worked out from the number of basis sets (1 - 3) supplied as ```basis```.
        * ```alpha```: alpha for the above scheme, if the default is to be overriden
        * ```options```: if special options are required for a step, they should be entered as a dict here. If some options should be used for both parts of the stage, they should be entered in both ```options``` and ```options_lo```. This is helpful for calculating all electron corrections in otherwise frozen core calculations, or relativistic (DKH) Hamiltionian corrections for otherwise nonrelativistic.
        * ```options_lo```: special options for lower method in a given stage. This is useful to calculate a direct stage in an otherwise density-fitted calculation, or similar.
        * ```treatment```: treat extrapolation stage as ```scf``` or ```corl```, by default only the first stage is ```scf``` and every later one is ```corl```.
        * ```stage```: tag for the stage used in tables.

        |  The next items in the ```cbs_metadata``` array extrapolate correlation. All of the above parameters are available, with only the ```wfn``` and ```basis``` keywords required. Other supported parameters are:

        * ```wfn_lo```: the lower method from which the delta correction is to be calculated. By default, it is set to ```wfn``` from the previous field in the ```cbs_metadata``` array.
        * ```basis_lo```: basis set to be used for the delta correction. By default, it is the same as the ```basis``` specified above.


    * Others

    :type molecule: :ref:`molecule <op_py_molecule>`
    :param molecule: ``h2o`` || etc.

        The target molecule, if not the last molecule defined.

    :examples:


    >>> # [1] replicates with cbs() the simple model chemistry scf/cc-pVDZ: set basis cc-pVDZ energy('scf')
    >>> energy(cbs, scf_wfn='scf', scf_basis='cc-pVDZ')

    >>> # [2] replicates with cbs() the simple model chemistry mp2/jun-cc-pVDZ: set basis jun-cc-pVDZ energy('mp2')
    >>> energy(cbs, corl_wfn='mp2', corl_basis='jun-cc-pVDZ')

    >>> # [3] DTQ-zeta extrapolated scf reference energy
    >>> energy('cbs', scf_wfn='scf', scf_basis='cc-pV[DTQ]Z', scf_scheme='scf_xtpl_helgaker_3')

    >>> # [4] DT-zeta extrapolated mp2 correlation energy atop a T-zeta reference
    >>> energy('cbs', corl_wfn='mp2', corl_basis='cc-pv[dt]z', corl_scheme='corl_xtpl_helgaker_2')

    >>> # [5] a DT-zeta extrapolated coupled-cluster correction atop a TQ-zeta extrapolated mp2 correlation energy atop a Q-zeta reference (both equivalent)
    >>> energy('cbs', corl_wfn='mp2', corl_basis='aug-cc-pv[tq]z', delta_wfn='ccsd(t)', delta_basis='aug-cc-pv[dt]z')
    >>> energy('cbs', corl_wfn='mp2', corl_basis='aug-cc-pv[tq]z', corl_scheme='corl_xtpl_helgaker_2', delta_wfn='ccsd(t)', delta_basis='aug-cc-pv[dt]z', delta_scheme='corl_xtpl_helgaker_2')

    >>> # [6] a D-zeta ccsd(t) correction atop a DT-zeta extrapolated ccsd cluster correction atop a TQ-zeta extrapolated mp2 correlation energy atop a Q-zeta reference
    >>> energy('cbs', corl_wfn='mp2', corl_basis='aug-cc-pv[tq]z', corl_scheme=corl_xtpl_helgaker_2, delta_wfn='ccsd', delta_basis='aug-cc-pv[dt]z', delta_scheme='corl_xtpl_helgaker_2', delta2_wfn='ccsd(t)', delta2_wfn_lesser='ccsd', delta2_basis='aug-cc-pvdz')

    >>> # [7] a Q5-zeta MP2 calculation, corrected by CCSD(T) at the TQ-zeta extrapolated level, and all-electron CCSD(T) correlation at T-zeta level
    >>> energy(cbs, cbs_metadata=[{"wfn": "hf", "basis": "cc-pv5z"}, {"wfn": "mp2", "basis": "cc-pv[q5]z"}, {"wfn": "ccsd(t)", "basis": "cc-pv[tq]z"}, {"wfn": "ccsd(t)", "basis": "cc-pvtz", "options": {"freeze_core": "False"}}])

    >>> # [8] cbs() coupled with database()
    >>> TODO database('mp2', 'BASIC', subset=['h2o','nh3'], symm='on', func=cbs, corl_basis='cc-pV[tq]z', corl_scheme='corl_xtpl_helgaker_2', delta_wfn='ccsd(t)', delta_basis='sto-3g')

    >>> # [9] cbs() coupled with optimize()
    >>> TODO optimize('mp2', corl_basis='cc-pV[DT]Z', corl_scheme='corl_xtpl_helgaker_2', func=cbs)

    """


##  Aliases  ##
complete_basis_set = cbs

#   LAB: below is a piece of pre-class cbs() that didn't make the transition. it has details, so preserving for future revival
#
#    #psioh = core.IOManager.shared_object()
#    #psioh.set_specific_retention(psif.PSIF_SCF_MOS, True)
#    # projection across point groups not allowed and cbs() usually a mix of symm-enabled and symm-tol calls
#    #   needs to be communicated to optimize() so reset by that optstash
#    core.set_local_option('SCF', 'GUESS_PERSIST', True)
#
#    # Run necessary computations
#    for mc in JOBS:
#        kwargs['name'] = mc['f_wfn']
#
#        # Build string of molecule and commands that are dependent on the database
#        commands = '\n'
#        commands += """\ncore.set_global_option('BASIS', '%s')\n""" % (mc['f_basis'])
#        commands += """core.set_global_option('WRITER_FILE_LABEL', '%s')\n""" % \
#            (user_writer_file_label + ('' if user_writer_file_label == '' else '-') + mc['f_wfn'].lower() + '-' + mc['f_basis'].lower())
#        exec(commands)
#
#    psioh.set_specific_retention(psif.PSIF_SCF_MOS, False)


def _expand_scheme_orders(scheme: str, basisname: List[str], basiszeta: List[int], wfnname: str, options: Dict) -> Dict[str, Dict[str, Any]]:
    """Check that the length of *basiszeta* array matches the implied degree of
    extrapolation in *scheme* name. Return a dictionary of same length as
    basiszeta, with *basisname* and *basiszeta* distributed therein.

    """
    Nxtpl = len(basiszeta)

    try:
        scheme.split()
    except AttributeError:
        raise UpgradeHelper(scheme, repr(scheme.__name__), 1.6, ' Replace extrapolation function with function name.')

    if scheme not in xtpl_procedures:
        raise ValidationError(f"Extrapolation function ({scheme}) not among registered extrapolation schemes: {list(xtpl_procedures.keys())}. Use 'register_xtpl_function' function.")

    if int(scheme.split("_")[-1]) != Nxtpl:
        raise ValidationError(f"""Call to '{scheme}' not valid with '{len(basiszeta)}' basis sets.""")

    NEED = {}
    for idx in range(Nxtpl):
        NEED[_lmh_labels[Nxtpl][idx]] = dict(
            zip(_f_fields, [wfnname, basisname[idx], basiszeta[idx], options, 0.0, None, None, None, None]))
    return NEED


def _contract_scheme_orders(needdict, datakey: str = 'f_energy') -> Dict[str, Any]:
    """Prepared named arguments for extrapolation functions by
    extracting zetas and values (which one determined by *datakey*) out
    of *needdict* and returning a dictionary whose keys are constructed
    from _lmh_labels.

    """
    largs = {}
    largs['functionname'] = needdict['HI']['f_wfn']
    Nxtpl = len(needdict)
    zlabels = _lmh_labels[Nxtpl]  # e.g., ['LO', 'HI']

    for zeta in range(Nxtpl):
        zlab = zlabels[zeta]  # e.g., LO
        largs['z' + zlab] = needdict[zlab]['f_zeta']
        largs['value' + zlab] = needdict[zlab][datakey]

    return largs


def _parse_cbs_gufunc_string(method_name: str):
    """ A helper function that parses a ``"method/basis"`` input string
    into separate method and basis components. Also handles delta corrections.

    Parameters
    ----------
    method_name
        A ``"method/basis"`` style string defining the calculation.

    Returns
    -------
    tuple
        Tuple in the ``(method_list, basis_list)`` format, where ``method_list``
        is the list of the component methods, and ``basis_list`` is the list of
        basis sets forming the extrapolation for each specified method.
        E.g. ``"mp2/cc-pv[tq]z+D:ccsd(t)/cc-pvtz"`` would return:
        ``(["mp2", "ccsd(t)"], ["cc-pv[tq]z", "cc-pvtz"])``.
    """

    method_name_list = re.split(r"""\+(?=\s*[Dd]:)""", method_name)
    if len(method_name_list) > 2:
        raise ValidationError(
            "CBS gufunc: Text parsing is only valid for a single delta, please use the CBS wrapper directly")

    method_list = []
    basis_list = []
    for num, method_str in enumerate(method_name_list):
        if (method_str.count("[") > 1) or (method_str.count("]") > 1):
            raise ValidationError(f"""CBS gufunc: Too many brackets given! {method_str}""")

        if method_str.count('/') != 1:
            raise ValidationError(f"""CBS gufunc: All methods must specify a basis with '/'. {method_str}""")

        if num > 0:
            method_str = method_str.strip()
            if method_str[:2].lower() != 'd:':
                raise ValidationError("""CBS gufunc: Delta method must start with 'D:'.""")
            else:
                method_str = method_str[2:]
        method, basis = method_str.split('/')
        method_list.append(method)
        basis_list.append(basis)
    return method_list, basis_list


def cbs_text_parser(total_method_name: str, **kwargs) -> Dict:
    """
    A text based parser of the CBS method string. Provided to handle "method/basis"
    specification of the requested calculations. Also handles "simple" (i.e.
    one-method and one-basis) calls.

    Parameters
    ----------
    total_method_name
        String in a ``"method/basis"`` syntax. Simple calls (e.g. ``"blyp/sto-3g"``) are
        bounced out of CBS. More complex calls (e.g. ``"mp2/cc-pv[tq]z"`` or
        ``"mp2/cc-pv[tq]z+D:ccsd(t)/cc-pvtz"``) are expanded by `_parse_cbs_gufunc_string()`
        and pushed through :py:func:`~psi4.driver.cbs`.

    Returns
    -------
    dict of updated CBS keyword arguments
    """

    ptype = kwargs.pop('ptype', None)

    # Sanitize total_method_name
    total_method_name = total_method_name.lower()
    total_method_name = total_method_name.replace(' ', '')

    # Split into components
    method_list, basis_list = _parse_cbs_gufunc_string(total_method_name)

    # Single energy call?
    single_call = len(method_list) == 1
    single_call &= '[' not in basis_list[0]
    single_call &= ']' not in basis_list[0]

    if single_call:
        method_name = method_list[0]
        basis = basis_list[0]

        return {'method': method_name, 'basis': basis}

    # Drop out for unsupported calls
    if ptype is None:
        raise ValidationError("A CBS call was detected, but no ptype was passed in. Please alert a dev.")
    elif ptype not in ["energy", "gradient", "hessian"]:
        raise ValidationError(f"{ptype.title()}: Cannot extrapolate or delta correct {ptype} yet.")

    # Catch kwarg issues for CBS methods only
    user_dertype = kwargs.pop('dertype', None)
    cbs_verbose = kwargs.pop('cbs_verbose', False)

    # If we are not a single call, let CBS wrapper handle it!
    cbs_kwargs = {}
    cbs_kwargs['ptype'] = ptype
    cbs_kwargs['verbose'] = cbs_verbose

    if user_dertype is not None:
        cbs_kwargs['dertype'] = user_dertype

    # Find method and basis
    metadata = []
    if method_list[0] in ['scf', 'hf', 'c4-scf', 'c4-hf']:
        stage = {}
        stage['wfn'] = method_list[0]
        stage['basis'] = basis_list[0]
        if 'scf_scheme' in kwargs:
            stage['scheme'] = kwargs.pop('scf_scheme')
        stage['stage'] = "scf"
        stage['treatment'] = "scf"
    else:
        # _validate_cbs_inputs will produce scf stage automatically
        stage = {}
        stage['wfn'] = method_list[0]
        stage['basis'] = basis_list[0]
        if 'corl_scheme' in kwargs:
            stage['scheme'] = kwargs.pop('corl_scheme')
        stage['stage'] = "corl"
        stage['treatment'] = "corl"
    metadata.append(stage)

    # "method/basis" syntax only allows for one delta correction
    # via "method/basis+D:delta/basis". Maximum length of method_list is 2.
    if len(method_list) == 2:
        stage = {}
        stage['wfn'] = method_list[1]
        stage['basis'] = basis_list[1]
        if 'delta_scheme' in kwargs:
            stage['scheme'] = kwargs.pop('delta_scheme')
        stage['stage'] = "delta1"
        stage['treatment'] = "corl"
        metadata.append(stage)

    cbs_kwargs["cbs_metadata"] = metadata

    return cbs_kwargs


def _build_cbs_compute(metameta: Dict[str, Any], metadata: CBSMetadata):
    label = metameta['label']
    ptype = metameta['ptype']
    verbose = metameta['verbose']

    # Build string of title banner
    instructions = "\n" + p4util.banner(f" CBS Setup{':' + label if label else ''} ", strNotOutfile=True) + "\n"

    # Call schemes for each portion of total energy to 'place orders' for calculations needed
    d_fields = [
        'd_stage', 'd_scheme', 'd_basis', 'd_wfn', 'd_alpha', 'd_need', 'd_coef', 'd_energy', 'd_gradient', 'd_hessian', 'd_dipole', 'd_dipder'
    ]
    GRAND_NEED = []

    NEED = _expand_scheme_orders(metadata[0]["scheme"], metadata[0]["basis"][0], metadata[0]["basis"][1],
                                 metadata[0]["wfn"], metadata[0]["options"])
    GRAND_NEED.append(
        dict(
            zip(d_fields, [
                'scf', metadata[0]["scheme"],
                _contract_bracketed_basis(metadata[0]["basis"][0]), metadata[0]["wfn"], metadata[0]["alpha"], NEED, +1,
                0.0, None, None, None, None
            ])))
    if len(metadata) > 1:
        for delta in metadata[1:]:
            NEED = _expand_scheme_orders(delta["scheme"], delta["basis"][0], delta["basis"][1], delta["wfn"],
                                         delta["options"])
            GRAND_NEED.append(
                dict(
                    zip(d_fields, [
                        delta["stage"], delta["scheme"],
                        _contract_bracketed_basis(delta["basis"][0]), delta["wfn"], delta["alpha"], NEED, +1, 0.0,
                        None, None, None, None
                    ])))
            NEED = _expand_scheme_orders(delta["scheme"], delta["basis_lo"][0], delta["basis_lo"][1], delta["wfn_lo"],
                                         delta["options_lo"])
            GRAND_NEED.append(
                dict(
                    zip(d_fields, [
                        delta["stage"], delta["scheme"],
                        _contract_bracketed_basis(delta["basis_lo"][0]), delta["wfn_lo"], delta["alpha"], NEED, -1,
                        0.0, None, None, None, None
                    ])))

    # MODELCHEM is unordered, possibly redundant list of single result *entries* needed to satisfy full CBS
    # JOBS is subset of MODELCHEM with minimal list of single result *jobs* needed to satisfy full CBS
    # TROVE is superset of JOBS with maximal list of single result *entries* resulting from JOBS
    # "entry" here is a mtd-bas-opt spec that can support E/G/H data
    # "job" here is an entry on which to sic Psi4 that, through VARH, may fill in multiple entries

    MODELCHEM = []
    for stage in GRAND_NEED:
        for lvl in stage['d_need'].values():
            MODELCHEM.append(lvl)

    # Apply chemical reasoning to choose the minimum computations to run
    JOBS = MODELCHEM[:]
    listfmt = """   {:>12} / {:24} for  {}{}\n"""

    # TODO: In the "naive" and "enlightened" loops below, I had to remove condition `and (job['f_options'] is not False))`
    #   to get them working, and I feel like they were added to fix the same thing. someday, seek to understand.

    #     Remove duplicate modelchem portion listings
    for mc in MODELCHEM:
        dups = -1
        for indx_job, job in enumerate(JOBS):
            if ((job['f_wfn'] == mc['f_wfn']) and (job['f_basis'] == mc['f_basis'])
                    and (job['f_options'] == mc['f_options'])):
                dups += 1
                if dups >= 1:
                    del JOBS[indx_job]

    instructions += """    Naive listing of computations required.\n"""
    for mc in JOBS:
        instructions += listfmt.format(mc['f_wfn'], mc['f_basis'] + " + options" * bool(mc['f_options']),
                                       VARH[mc['f_wfn']][mc['f_wfn']], _addlremark[ptype])

    #     Remove chemically subsumed modelchem portion listings
    if ptype == 'energy':
        for mc in MODELCHEM:
            for wfn in VARH[mc['f_wfn']]:
                for indx_job, job in enumerate(JOBS):
                    if ((VARH[mc['f_wfn']][wfn] == VARH[job['f_wfn']][job['f_wfn']])
                            and (mc['f_basis'] == job['f_basis'])
                            and not (mc['f_wfn'] == job['f_wfn'])
                            and (mc['f_options'] == job['f_options'])):
                        del JOBS[indx_job]

    instructions += """\n    Enlightened listing of computations required.\n"""
    for mc in JOBS:
        instructions += listfmt.format(mc['f_wfn'], mc['f_basis'] + " + options" * bool(mc['f_options']),
                                       VARH[mc['f_wfn']][mc['f_wfn']], _addlremark[ptype])

    #     Expand listings to all that will be obtained
    TROVE = []
    for job in JOBS:
        for wfn in VARH[job['f_wfn']]:
            TROVE.append(dict(zip(_f_fields, [wfn, job['f_basis'], job['f_zeta'], job['f_options'], 0.0, None, None, None, None])))

    instructions += """\n    Full listing of computations to be obtained (required and bonus).\n"""
    for mc in TROVE:
        instructions += listfmt.format(mc['f_wfn'], mc['f_basis'] + " + options" * bool(mc['f_options']),
                                       VARH[mc['f_wfn']][mc['f_wfn']], _addlremark[ptype])
    if verbose:
        core.print_out(instructions)
        logger.info(instructions)

    return GRAND_NEED, JOBS, TROVE


def _assemble_cbs_components(metameta, TROVE, GRAND_NEED):
    """Absorb job E/G/H results from `TROVE` into `GRAND_NEED`. Process
    those into stage E/G/H in `GRAND_NEED`, returning the latter.
    Accumulate into final E/G/H quantities, returning them in dict.

    """
    label = metameta['label']
    nat = metameta['molecule'].natom()
    ptype = metameta['ptype']
    verbose = metameta['verbose']

    # Build string of title banner
    instructions = "\n" + p4util.banner(f" CBS Results{':' + label if label else ''} ", strNotOutfile=True) + "\n"
    core.print_out(instructions)
    logger.info(instructions)

    # Insert obtained energies into the array that stores the cbs stages
    for stage in GRAND_NEED:
        for lvl in stage['d_need'].values():
            for job in TROVE:
                # Don't ask
                if (((lvl['f_wfn'] == job['f_wfn']) or
                     ((lvl['f_wfn'][3:] == job['f_wfn']) and lvl['f_wfn'].startswith('c4-')) or
                     ((lvl['f_wfn'] == job['f_wfn'][3:]) and job['f_wfn'].startswith('c4-')) or
                     (('c4-' + lvl['f_wfn']) == job['f_wfn']) or (lvl['f_wfn'] == ('c4-' + job['f_wfn'])))
                        and (lvl['f_basis'] == job['f_basis']) and (lvl['f_options'] == job['f_options'])):
                    lvl['f_energy'] = job['f_energy']
                    lvl['f_gradient'] = job['f_gradient']
                    lvl['f_hessian'] = job['f_hessian']
                    lvl['f_dipole'] = job['f_dipole']
                    lvl['f_dipder'] = job['f_dipder']

    # Make xtpl() call
    finalenergy = 0.0
    finalgradient = None
    finalhessian = None
    finaldipole = None
    finaldipder = None

    for stage in GRAND_NEED:
        hiloargs = {'alpha': stage['d_alpha'], 'verbose': verbose}

        grad_available = all([lmh['f_gradient'] is not None for lmh in stage['d_need'].values()])
        hess_available = all([lmh['f_hessian'] is not None for lmh in stage['d_need'].values()])
        dipole_available = all([lmh['f_dipole'] is not None for lmh in stage['d_need'].values()])
        dipder_available = all([lmh['f_dipder'] is not None for lmh in stage['d_need'].values()])

        hiloargs.update(_contract_scheme_orders(stage['d_need'], 'f_energy'))
        stage['d_energy'] = xtpl_procedures[stage['d_scheme']](**hiloargs)
        finalenergy += stage['d_energy'] * stage['d_coef']

        if ptype == 'gradient' or grad_available:
            if finalgradient is None:
                finalgradient = np.zeros((nat, 3))
            hiloargs.update(_contract_scheme_orders(stage['d_need'], 'f_gradient'))
            stage['d_gradient'] = xtpl_procedures[stage['d_scheme']](**hiloargs)
            finalgradient += stage['d_gradient'] * stage['d_coef']

        if ptype == 'hessian' or hess_available:
            if finalhessian is None:
                finalhessian = np.zeros((3 * nat, 3 * nat))
            hiloargs.update(_contract_scheme_orders(stage['d_need'], 'f_hessian'))
            stage['d_hessian'] = xtpl_procedures[stage['d_scheme']](**hiloargs)
            finalhessian += stage['d_hessian'] * stage['d_coef']

        if dipole_available:
            if finaldipole is None:
                finaldipole = np.zeros((3))
            hiloargs.update(_contract_scheme_orders(stage['d_need'], 'f_dipole'))
            stage['d_dipole'] = xtpl_procedures[stage['d_scheme']](**hiloargs)
            finaldipole += stage['d_dipole'] * stage['d_coef']

        if dipder_available:
            if finaldipder is None:
                finaldipder = np.zeros((3 * nat, 3))
            hiloargs.update(_contract_scheme_orders(stage['d_need'], 'f_dipder'))
            stage['d_dipder'] = xtpl_procedures[stage['d_scheme']](**hiloargs)
            finaldipder += stage['d_dipder'] * stage['d_coef']

    cbs_results = {
        'ret_ptype': {
            'energy': finalenergy,
            'gradient': finalgradient,
            'hessian': finalhessian,
        }[ptype],
        'energy': finalenergy,
        'gradient': finalgradient,
        'hessian': finalhessian,
        'dipole': finaldipole,
        'dipole gradient': finaldipder,
    }

    return cbs_results, GRAND_NEED


def _summary_table(metadata, TROVE, GRAND_NEED) -> str:
    """Build string of results table"""

    delimit = '  ' + '-' * 105 + '\n'
    blckfmt = """\n   ==> {} <==\n\n"""
    headfmt = """     {:>6} {:>20} {:1} {:26} {:>3} {:>16}   {}\n"""
    linefmt = """     {:>6} {:>20} {:1} {:27} {:2} {:16.8f}   {}\n"""

    tables = ''
    tables += blckfmt.format('Components')
    tables += delimit

    required = []
    finalenergy = 0.0
    for stage in GRAND_NEED:
        finalenergy += stage['d_energy'] * stage['d_coef']
        for lvl in stage['d_need'].values():
            required.append((lvl['f_wfn'], lvl['f_basis'], lvl['f_options']))

    tables += headfmt.format('', 'Method', '/', 'Basis', 'Rqd', 'Energy [Eh]', 'Variable')
    tables += delimit
    for job in TROVE:
        star = ''
        for mc in required:
            if (job['f_wfn'], job['f_basis'], job['f_options']) == mc:
                star = '*'
        tables += linefmt.format('', job['f_wfn'], '/', job['f_basis'] + " + options" * bool(job['f_options']), star,
                                 job['f_energy'], VARH[job['f_wfn']][job['f_wfn']])
    tables += delimit

    tables += blckfmt.format('Stages')
    tables += delimit
    tables += headfmt.format('Stage', 'Method', '/', 'Basis', 'Wt', 'Energy [Eh]', 'Scheme')
    tables += delimit
    for stage in GRAND_NEED:
        tables += linefmt.format(stage['d_stage'], stage['d_wfn'], '/', stage['d_basis'], stage['d_coef'],
                                 stage['d_energy'], stage['d_scheme'])
    tables += delimit

    tables += blckfmt.format('CBS')
    tables += delimit
    tables += headfmt.format('Stage', 'Method', '/', 'Basis', '', 'Energy [Eh]', 'Scheme')
    tables += delimit
    tables += linefmt.format(GRAND_NEED[0]['d_stage'], GRAND_NEED[0]['d_wfn'], '/', GRAND_NEED[0]['d_basis'], '',
                             GRAND_NEED[0]['d_energy'], GRAND_NEED[0]['d_scheme'])

    if len(metadata) > 1:
        dc = 1
        for delta in metadata[1:]:
            mtdstr = GRAND_NEED[dc]['d_wfn']
            if dc != 1:
                mtdstr += ' - ' + GRAND_NEED[dc + 1]['d_wfn']
            tables += linefmt.format(GRAND_NEED[dc]['d_stage'], mtdstr, '/', GRAND_NEED[dc]['d_basis'], '',
                                     GRAND_NEED[dc]['d_energy'] - GRAND_NEED[dc + 1]['d_energy'],
                                     GRAND_NEED[dc]['d_scheme'])
            dc += 2

    tables += linefmt.format('total', 'CBS', '', '', '', finalenergy, '')
    tables += delimit

    return tables


class CompositeComputer(BaseComputer):

    molecule: Any
    basis: str = "(auto)"
    method: str = "(auto)"
    driver: DriverEnum
    keywords: Dict[str, Any] = {}
    metadata: Any
    metameta: Dict[str, Any] = {}

    verbose: int = 1

    # List of model chemistries with extrapolation scheme applied. Can reconstruct CBS. Keys are d_fields. Formerly GRAND_NEED.
    cbsrec: List[Dict[str, Any]] = []

    # Maximal list of model chemistries extractable from running `compute_list`. Keys are _f_fields. Formerly JOBS_EXT.
    trove: List[Dict[str, Any]] = []

    # Minimal (enlightened) list of jobs to run to satisfy full CBS. Keys are _f_fields. Formerly JOBS.
    compute_list: List[Dict[str, Any]] = []

    # One-to-One list of AtomicComputer-s corresponding to `compute_list`.
    task_list: List[AtomicComputer] = []

    # One-to-One list of QCSchema corresponding to `task_list`.
    results_list: List[Any] = []

    @validator('molecule')
    def set_molecule(cls, mol):
        mol.update_geometry()
        mol.fix_com(True)
        mol.fix_orientation(True)
        return mol

    def __init__(self, **data):
        data = p4util.kwargs_lower(data)
        data["metadata"] = _process_cbs_kwargs(data)
        BaseComputer.__init__(self, **data)

        self.metameta = {
            'kwargs': data,
            'ptype': self.driver,
            'verbose': self.verbose,
            'label': None,
            'molecule': self.molecule,
        }
        # logger.debug("METAMETA\n" + pp.pformat(self.metameta))

        if data['metadata']:
            if data['metadata'][0]["wfn"] not in VARH.keys():
                raise ValidationError(
                    """Requested SCF method '%s' is not recognized. Add it to VARH in driver_cbs.py to proceed.""" %
                    (metadata[0]["wfn"]))

            if len(self.metadata) > 1:
                for delta in self.metadata[1:]:
                    if delta["wfn"] not in VARH.keys():
                        raise ValidationError(
                            f"""Requested higher {delta["treatment"]} method '{delta["wfn"]}' is not recognized. Add it to VARH in driver_cbs.py to proceed."""
                        )
                    if delta["wfn_lo"] not in VARH.keys():
                        raise ValidationError(
                            f"""Requested lesser {delta["treament"]} method '{delta["wfn_lo"]}' is not recognized. Add it to VARH in driver_cbs.py to proceed."""
                        )

            self.cbsrec, self.compute_list, self.trove = _build_cbs_compute(self.metameta, self.metadata)

            for job in self.compute_list:
                keywords = copy.deepcopy(self.metameta['kwargs']['keywords'])
                if job["f_options"] is not False:
                    stage_keywords = dict(job["f_options"].items())
                    keywords = {**keywords, **stage_keywords}
                task = AtomicComputer(
                    **{
                        "molecule": self.molecule,
                        "driver": self.driver,
                        "method": job["f_wfn"],
                        "basis": job["f_basis"],
                        "keywords": keywords or {},
                    })
                self.task_list.append(task)

                # logger.debug("TASK\n" + pp.pformat(task.dict()))

    def build_tasks(self, obj, **kwargs):
        # permanently a dummy function
        pass

    def plan(self):
        # uncalled function
        return [t.plan() for t in self.task_list]

    def compute(self, client: Optional["qcportal.FractalClient"] = None):
        label = self.metameta['label']
        instructions = "\n" + p4util.banner(f" CBS Computations{':' + label if label else ''} ",
                                            strNotOutfile=True) + "\n"
        logger.debug(instructions)
        core.print_out(instructions)

        with p4util.hold_options_state():
            for t in reversed(self.task_list):
                t.compute(client=client)

    def _prepare_results(self, client: Optional["qcportal.FractalClient"] = None):
        results_list = [x.get_results(client=client) for x in self.task_list]

        modules = [getattr(v.provenance, "module", None) for v in results_list]
        if self.driver != "energy" and len(set(modules)) == 2 and modules.count("scf") == len(modules) / 2:
            # signature of "MP2 GRAD" - "HF GRAD" implementation detail
            # * avoid having post-scf single-method gradients/Hessians show up as "(mixed)" module just because an outright HF call in the jobs list
            modules = set(modules) - {"scf"}
        modules = list(set(modules))
        modules = modules[0] if len(modules) == 1 else "(mixed)"

        # load results_list numbers into compute_list (task_list is AtomicComputer-s)
        for itask, mc in enumerate(self.compute_list):
            task = results_list[itask]
            response = task.return_result

            if self.metameta['ptype'] == 'energy':
                mc['f_energy'] = response

            elif self.metameta['ptype'] == 'gradient':
                mc['f_gradient'] = response
                mc['f_energy'] = task.extras['qcvars']['CURRENT ENERGY']

            elif self.metameta['ptype'] == 'hessian':
                mc['f_hessian'] = response
                mc['f_energy'] = task.extras['qcvars']['CURRENT ENERGY']
                if 'CURRENT GRADIENT' in task.extras['qcvars']:
                    mc['f_gradient'] = task.extras['qcvars']['CURRENT GRADIENT']

            if 'CURRENT DIPOLE' in task.extras['qcvars']:
                mc['f_dipole'] = task.extras['qcvars']['CURRENT DIPOLE']

            if 'CURRENT DIPOLE GRADIENT' in task.extras['qcvars']:
                mc['f_dipder'] = task.extras['qcvars']['CURRENT DIPOLE GRADIENT']

            # Fill in energies for subsumed methods
            if self.metameta['ptype'] == 'energy':
                for wfn in VARH[mc['f_wfn']]:
                    for job in self.trove:
                        if ((wfn == job['f_wfn']) and (mc['f_basis'] == job['f_basis'])
                                and (mc['f_options'] == job['f_options'])):
                            job['f_energy'] = task.extras['qcvars'][VARH[wfn][wfn]]

            # Copy data from 'run' to 'obtained' table
            for mce in self.trove:
                if ((mc['f_wfn'] == mce['f_wfn']) and (mc['f_basis'] == mce['f_basis'])
                        and (mc['f_options'] == mce['f_options'])):
                    mce['f_energy'] = mc['f_energy']
                    mce['f_gradient'] = mc['f_gradient']
                    mce['f_hessian'] = mc['f_hessian']
                    mce['f_dipole'] = mc['f_dipole']
                    mce['f_dipder'] = mc['f_dipder']

            # logger.debug("MC\n" + pp.pformat(mc))

        cbs_results, self.cbsrec = _assemble_cbs_components(self.metameta, self.trove, self.cbsrec)

        instructions = _summary_table(self.metadata, self.trove, self.cbsrec)
        core.print_out(instructions)
        logger.info(instructions)

        # logger.debug('CBS_RESULTS\n' + pp.pformat(cbs_results))
        # logger.debug('GRAND_NEED\n' + pp.pformat(self.cbsrec))

        cbs_results["module"] = modules
        return cbs_results

    def get_results(self, client: Optional["qcportal.FractalClient"] = None) -> AtomicResult:
        """Return results as Composite-flavored QCSchema."""

        assembled_results = self._prepare_results(client=client)
        E0 = assembled_results["energy"]

        # load QCVariables & properties
        qcvars = {
            'CBS NUMBER': len(self.compute_list),
            'NUCLEAR REPULSION ENERGY': self.molecule.nuclear_repulsion_energy(),
        }

        properties = {
            "calcinfo_natom": self.molecule.natom(),
            "nuclear_repulsion_energy": self.molecule.nuclear_repulsion_energy(),
            "return_energy": E0,
        }

        for qcv in ['CBS', 'CURRENT']:
            qcvars[qcv + ' REFERENCE ENERGY'] = self.cbsrec[0]['d_energy']
            qcvars[qcv + ' CORRELATION ENERGY'] = E0 - self.cbsrec[0]['d_energy']
            qcvars[qcv + ('' if qcv == 'CURRENT' else ' TOTAL') + ' ENERGY'] = E0

        for idelta in range(int(len(self.cbsrec) / 2)):
            if idelta == 0:
                continue
            dc = idelta * 2 + 1
            qcvars[f"CBS {self.cbsrec[dc]['d_stage'].upper()} TOTAL ENERGY"] = self.cbsrec[dc]["d_energy"] - self.cbsrec[dc + 1]["d_energy"]

        G0 = assembled_results["gradient"]
        if G0 is not None:
            qcvars["CURRENT GRADIENT"] = G0
            qcvars["CBS TOTAL GRADIENT"] = G0
            properties["return_gradient"] = G0

        H0 = assembled_results["hessian"]
        if H0 is not None:
            qcvars["CURRENT HESSIAN"] = H0
            qcvars["CBS TOTAL HESSIAN"] = H0
            properties["return_hessian"] = H0

        D0 = assembled_results["dipole"]
        if D0 is not None:
            qcvars["CURRENT DIPOLE"] = D0
            qcvars["CBS DIPOLE"] = D0

        DD0 = assembled_results["dipole gradient"]
        if DD0 is not None:
            qcvars["CURRENT DIPOLE GRADIENT"] = DD0
            qcvars["CBS DIPOLE GRADIENT"] = DD0

        cbs_model = AtomicResult(
            **{
                'driver': self.driver,
                #'keywords': self.keywords,
                'model': {
                    'method': self.method,
                    'basis': self.basis,
                },
                'molecule': self.molecule.to_schema(dtype=2),
                'properties': properties,
                'provenance': p4util.provenance_stamp(__name__, module=assembled_results["module"]),
                'extras': {
                    'qcvars': qcvars,
                    'cbs_record': copy.deepcopy(self.cbsrec),
                },
                'return_result': assembled_results['ret_ptype'],
                'success': True,
            })

        logger.debug('CBS QCSchema:\n' + pp.pformat(cbs_model.dict()))

        return cbs_model

    def get_psi_results(
        self,
        client: Optional["qcportal.FractalClient"] = None,
        *,
        return_wfn: bool = False) -> EnergyGradientHessianWfnReturn:
        """Called by driver to assemble results into Composite-flavored QCSchema,
        then reshape and return them in the customary Psi4 driver interface: ``(e/g/h, wfn)``.

        Parameters
        ----------
        return_wfn
            Whether to additionally return the dummy :py:class:`~psi4.core.Wavefunction`
            calculation result as the second element of a tuple. Contents are:

            - molecule
            - dummy basis, def2-svp
            - e/g/h member data
            - QCVariables
            - module if simple

        Returns
        -------
        ret
            Energy, gradient, or Hessian according to self.driver.
        wfn
            Wavefunction described above when *return_wfn* specified.

        """
        cbs_model = self.get_results(client=client)

        if cbs_model.driver == 'energy':
            ret_ptype = cbs_model.return_result
        else:
            ret_ptype = core.Matrix.from_array(cbs_model.return_result)
        wfn = _cbs_schema_to_wfn(cbs_model)

        if return_wfn:
            return (ret_ptype, wfn)
        else:
            return ret_ptype


def _cbs_schema_to_wfn(cbs_model):
    """Helper function to produce Wavefunction from a Composite-flavored AtomicResult."""

    mol = core.Molecule.from_schema(cbs_model.molecule.dict())
    basis = core.BasisSet.build(mol, "ORBITAL", 'def2-svp', quiet=True)
    wfn = core.Wavefunction(mol, basis)
    if hasattr(cbs_model.provenance, "module"):
        wfn.set_module(cbs_model.provenance.module)

    # wfn.set_energy(cbs_model['extras'['qcvars'].get('CBS TOTAL ENERGY'))  # catches Wfn.energy_
    for qcv, val in cbs_model.extras['qcvars'].items():
        for obj in [core, wfn]:
            obj.set_variable(qcv, val)

    return wfn
