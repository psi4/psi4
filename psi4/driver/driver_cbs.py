#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2022 The Psi4 Developers.
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

import math
import re
import sys
import copy
import pprint
from typing import Any, Callable, Dict, List, Union

import numpy as np
import pydantic

from psi4 import core
from psi4.driver import qcdb
from psi4.driver import p4util
from psi4.driver import driver_util
from psi4.driver import psifiles as psif
from psi4.driver.driver_cbs_helper import xtpl_procedures, register_xtpl_scheme
from psi4.driver.p4util.exceptions import ValidationError
from psi4.driver.procrouting.interface_cfour import cfour_psivar_list
from psi4.driver.task_base import BaseTask, SingleResult, unnp, plump_qcvar

zeta_values = 'dtq5678'
_zeta_val2sym = {k + 2: v for k, v in enumerate(zeta_values)}
_zeta_sym2val = {v: k for k, v in _zeta_val2sym.items()}
_addlremark = {'energy': '', 'gradient': ', GRADIENT', 'hessian': ', HESSIAN'}
_f_fields = ['f_wfn', 'f_basis', 'f_zeta', 'f_options', 'f_energy', 'f_gradient', 'f_hessian']
_lmh_labels = {
    1: ['HI'],
    2: ['LO', 'HI'],
    3: ['LO', 'MD', 'HI'],
    4: ['LO', 'MD', 'M2', 'HI'],
    5: ['LO', 'MD', 'M2', 'M3', 'HI']
}


def _expand_bracketed_basis(basisstring: str, molecule=None):
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
    molecule : qcdb.molecule or psi4.core.Molecule
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


def _contract_bracketed_basis(basisarray: List):
    """Function to reform a bracketed basis set string from a sequential series
    of basis sets. Essentially the inverse of _expand_bracketed_basis(). Used to
    print a nicely formatted basis set string in the results table.

    Parameters
    ----------
    basisarray
        Basis set names, differing by zeta level, e.g. ``["cc-pvqz", "cc-pv5z"]``.

    Returns
    -------
    string
        A nicely formatted basis set string, e.g. ``"cc-pv[q5]z"`` for the above example.

    """

    if len(basisarray) == 1:
        return basisarray[0]

    else:
        zetaindx = [i for i in range(len(basisarray[0])) if basisarray[0][i] != basisarray[1][i]][0]
        ZSET = [bas[zetaindx] for bas in basisarray]

        pre = basisarray[1][:zetaindx]
        post = basisarray[1][zetaindx + 1:]
        basisstring = pre + '[' + ''.join(ZSET) + ']' + post
        return basisstring


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
    VARH['lccd'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                           'lccd': 'LCCD TOTAL ENERGY'}
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
    VARH['mrccsd'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                         'mrccsd': 'CCSD TOTAL ENERGY'}
    VARH['mrccsd(t)'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                         'mrccsd': 'CCSD TOTAL ENERGY',
                      'mrccsd(t)': 'CCSD(T) TOTAL ENERGY'}
    VARH['mrccsdt'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                        'mrccsdt': 'CCSDT TOTAL ENERGY'}
    VARH['mrccsdt(q)'] = {
                             'hf': 'HF TOTAL ENERGY',
                            'mp2': 'MP2 TOTAL ENERGY',
                        'mrccsdt': 'CCSDT TOTAL ENERGY',
                     'mrccsdt(q)': 'CCSDT(Q) TOTAL ENERGY'}

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
        return 'xtpl_highest_1'
    elif xtpl_type == "scf":
        if nbasis == 2:
            return 'scf_xtpl_helgaker_2'
        elif nbasis == 3:
            return 'scf_xtpl_helgaker_3'
        else:
            raise ValidationError(f"Wrong number of basis sets supplied to scf_xtpl: {nbasis}")
    elif xtpl_type == "corl":
        if nbasis == 2:
            return 'corl_xtpl_helgaker_2'
        else:
            raise ValidationError(f"Wrong number of basis sets supplied to corl_xtpl: {nbasis}")
    else:
        raise ValidationError(f"Stage treatment must be 'corl' or 'scf', not '{xtpl_type}'")


def _validate_cbs_inputs(cbs_metadata, molecule):
    """ A helper function which validates the ``cbs_metadata`` format,
    expands basis sets, and provides sensible defaults for optional arguments.

    Parameters
    ----------
    cbs_metadata : list
        List of dicts containing CBS stage keywords.
    molecule : qcdb.molecule or psi4.core.Molecule
        Molecule to be passed to _expand_bracketed_basis()

    Returns
    -------
    list
        Validated list of dictionaries, with each item consisting of an extrapolation
        stage. All validation takes place here.
    """

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
    return (metadata)


def _process_cbs_kwargs(kwargs):
    """ A helper function which translates supplied kwargs into the
    ``cbs_metadata`` format and passes it for validation.

    Parameters
    ----------
    kwargs : dict
        kwargs containing the CBS function specification.

    Returns
    -------
    list
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

    :type scf_scheme: string
    :param scf_scheme: |dl| ``'xtpl_highest_1'`` |dr| || ``'scf_xtpl_helgaker_3'`` || etc.

        Indicates the basis set extrapolation scheme to be applied to the reference energy.
        Defaults to :py:func:`~psi4.driver.driver_cbs.scf_xtpl_helgaker_3` if three valid basis sets
        present in ``psi4.driver.driver_cbs.scf_basis``, :py:func:`~psi4.driver.driver_cbs.scf_xtpl_helgaker_2` if two valid basis
        sets present in ``scf_basis``, and :py:func:`~psi4.driver.driver_cbs.xtpl_highest_1` otherwise.

        .. hlist::
           :columns: 1

           * :py:func:`~psi4.driver.driver_cbs.xtpl_highest_1`
           * :py:func:`~psi4.driver.driver_cbs.scf_xtpl_helgaker_3`
           * :py:func:`~psi4.driver.driver_cbs.scf_xtpl_helgaker_2`
           * :py:func:`~psi4.driver.driver_cbs.scf_xtpl_truhlar_2`
           * :py:func:`~psi4.driver.driver_cbs.scf_xtpl_karton_2`

    :type corl_scheme: string
    :param corl_scheme: |dl| ``'xtpl_highest_1'`` |dr| || ``'corl_xtpl_helgaker_2'`` || etc.

        Indicates the basis set extrapolation scheme to be applied to the correlation energy.
        Defaults to :py:func:`~psi4.driver.driver_cbs.corl_xtpl_helgaker_2` if two valid basis sets
        present in ``corl_basis`` and :py:func:`~psi4.driver.driver_cbs.xtpl_highest_1` otherwise.

        .. hlist::
           :columns: 1

           * :py:func:`~psi4.driver.driver_cbs.xtpl_highest_1`
           * :py:func:`~psi4.driver.driver_cbs.corl_xtpl_helgaker_2`

    :type delta_scheme: string
    :param delta_scheme: |dl| ``'xtpl_highest_1'`` |dr| || ``'corl_xtpl_helgaker_2'`` || etc.

        Indicates the basis set extrapolation scheme to be applied to the delta correction
        to the correlation energy.
        Defaults to :py:func:`~psi4.driver.driver_cbs.corl_xtpl_helgaker_2` if two valid basis sets
        present in ``delta_basis`` and :py:func:`~psi4.driver.driver_cbs.xtpl_highest_1` otherwise.

        .. hlist::
           :columns: 1

           * :py:func:`~psi4.driver.driver_cbs.xtpl_highest_1`
           * :py:func:`~psi4.driver.driver_cbs.corl_xtpl_helgaker_2`

    :type delta2_scheme: string
    :param delta2_scheme: |dl| ``'xtpl_highest_1'`` |dr| || ``'corl_xtpl_helgaker_2'`` || etc.

        Indicates the basis set extrapolation scheme to be applied to the second delta correction
        to the correlation energy.
        Defaults to :py:func:`~psi4.driver.driver_cbs.corl_xtpl_helgaker_2` if two valid basis sets
        present in ``delta2_basis`` and :py:func:`~psi4.driver.driver_cbs.xtpl_highest_1` otherwise.

        .. hlist::
           :columns: 1

           * :py:func:`~psi4.driver.driver_cbs.xtpl_highest_1`
           * :py:func:`~psi4.driver.driver_cbs.corl_xtpl_helgaker_2`

    :type scf_alpha: float
    :param scf_alpha: |dl| ``1.63`` |dr|

        Overrides the default \alpha parameter used in the listed SCF extrapolation procedures.
        Has no effect on others, including :py:func:`~psi4.driver.driver_cbs.xtpl_highest_1` and :py:func:`~psi4.driver.driver_cbs.scf_xtpl_helgaker_3`.

        .. hlist::
           :columns: 1

           * :py:func:`~psi4.driver.driver_cbs.scf_xtpl_helgaker_2`
           * :py:func:`~psi4.driver.driver_cbs.scf_xtpl_truhlar_2`
           * :py:func:`~psi4.driver.driver_cbs.scf_xtpl_karton_2`

    :type corl_alpha: float
    :param corl_alpha: |dl| ``3.00`` |dr| 

        Overrides the default \alpha parameter used in the listed :py:func:`~psi4.driver.driver_cbs.corl_xtpl_helgaker_2` correlation
        extrapolation to the corl stage. The supplied \alpha does not impact delta or any further stages.

        .. hlist::
           :columns: 1

           * :py:func:`~psi4.driver.driver_cbs.corl_xtpl_helgaker_2`

    :type delta_alpha: float
    :param delta_alpha: |dl| ``3.00`` |dr| 

        Overrides the default \alpha parameter used in the listed
        :py:func:`~psi4.driver.driver_cbs.corl_xtpl_helgaker_2` correlation extrapolation for the delta correction. Useful when
        delta correction is performed using smaller basis sets for which a different \alpha might
        be more appropriate.

        .. hlist::
           :columns: 1

           * :py:func:`~psi4.driver.driver_cbs.corl_xtpl_helgaker_2`

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
    pass


##  Aliases  ##
complete_basis_set = cbs




######### OUTSOURCE
######### PARSE / PREPARE
######### PREPARE / COMPUTE

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


def _expand_scheme_orders(scheme, basisname, basiszeta, wfnname, options, natom):
    """Check that the length of *basiszeta* array matches the implied degree of
    extrapolation in *scheme* name. Return a dictionary of same length as
    basiszeta, with *basisname* and *basiszeta* distributed therein.

    """
    Nxtpl = len(basiszeta)

    try:
        scheme.split()
    except AttributeError:
        raise UpgradeHelper(scheme, repr(scheme.__name__), 1.4, ' Replace extrapolation function with function name.')

    if scheme not in xtpl_procedures:
        raise ValidationError('Extrapolation function ({}) not among registered extrapolation schemes: {}'.format(
            scheme, list(xtpl_procedures.keys())))

    if int(scheme.split('_')[-1]) != Nxtpl:
        raise ValidationError(f"""Call to '{scheme}' not valid with '{len(basiszeta)}' basis sets.""")

    NEED = {}
    for idx in range(Nxtpl):
        NEED[_lmh_labels[Nxtpl][idx]] = dict(
            zip(_f_fields, [
                wfnname, basisname[idx], basiszeta[idx], options, 0.0,
                np.zeros((natom, 3)),
                np.zeros((3 * natom, 3 * natom))
            ]))
    return NEED


def _contract_scheme_orders(needdict, datakey='f_energy'):
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


def _cbs_wrapper_methods(**kwargs):
    """ A helper function for the driver to enumerate methods used in the
    stages of a cbs calculation.

    Parameters
    ----------
    kwargs : dict
        kwargs containing cbs specification either in the ``cbs_metadata``
        format, or in separate keywords (``scf_wfn``, ``corl_wfn`` etc.).

    Returns
    -------
    list
        List containing method name for each active stage.
    """

    cbs_methods = []
    if "cbs_metadata" in kwargs:
        for item in kwargs["cbs_metadata"]:
            cbs_methods.append(item.get("wfn"))
    else:
        cbs_method_kwargs = ['scf_wfn', 'corl_wfn', 'delta_wfn']
        cbs_method_kwargs += [f'delta{x}_wfn' for x in range(2, 6)]
        for method in cbs_method_kwargs:
            if method in kwargs:
                cbs_methods.append(kwargs[method])
    return cbs_methods


def _parse_cbs_gufunc_string(method_name):
    """ A helper function that parses a ``"method/basis"`` input string
    into separate method and basis components. Also handles delta corrections.

    Parameters
    ----------
    method_name : str
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


def _cbs_text_parser(total_method_name, **kwargs):
    """
    A text based parser of the CBS method string. Provided to handle "method/basis"
    specification of the requested calculations. Also handles "simple" (i.e.
    one-method and one-basis) calls.

    Parameters
    ----------
    total_method_name : str
        String in a ``"method/basis"`` syntax. Simple calls (e.g. ``"blyp/sto-3g"``) are
        bounced out of CBS. More complex calls (e.g. ``"mp2/cc-pv[tq]z"`` or
        ``"mp2/cc-pv[tq]z+D:ccsd(t)/cc-pvtz"``) are expanded by `_parse_cbs_gufunc_string()`
        and pushed through :py:func:`~psi4.cbs`.

    Returns
    -------
    dict of updated CBS keyword arguments
    """

    ptype = kwargs.pop('ptype', None)

    # Sanitize total_method_name
    label = total_method_name
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
        raise ValidationError("A CBS call was detected, but no ptyped was passed in. Please alert a dev.")
    elif ptype not in ["energy", "gradient", "hessian"]:
        raise ValidationError("%s: Cannot extrapolate or delta correct %s yet." % (ptype.title(), ptype))

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


def cbs_gufunc(func, total_method_name, **kwargs):
    """
    A text based wrapper of the CBS function. Provided to handle "method/basis"
    specification of the requested calculations. Also handles "simple" (i.e.
    one-method and one-basis) calls.

    Parameters
    ----------
    func : function
        Function to be called (energy, gradient, frequency or cbs).
    total_method_name : str
        String in a ``"method/basis"`` syntax. Simple calls (e.g. ``"blyp/sto-3g"``) are
        bounced out of CBS. More complex calls (e.g. ``"mp2/cc-pv[tq]z"`` or
        ``"mp2/cc-pv[tq]z+D:ccsd(t)/cc-pvtz"``) are expanded by `_parse_cbs_gufunc_string()`
        and pushed through :py:func:`~psi4.cbs`.

    Returns
    -------
    tuple or float
        Float, or if ``return_wfn`` is specified, a tuple of ``(value, wavefunction)``.

    """

    # Catch kwarg issues for all methods
    kwargs = p4util.kwargs_lower(kwargs)
    return_wfn = kwargs.pop('return_wfn', False)
    core.clean_variables()

    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop('molecule', core.get_active_molecule())
    molecule.update_geometry()

    cbs_kwargs = _cbs_text_parser(total_method_name, **kwargs)
    cbs_kwargs['molecule'] = molecule
    cbs_kwargs['return_wfn'] = True

    if 'cbs_metadata' not in cbs_kwargs:
        # Single call
        method_name = cbs_kwargs['method']
        basis = cbs_kwargs['basis']

        # Save some global variables so we can reset them later
        optstash = p4util.OptionsState(['BASIS'])
        core.set_global_option('BASIS', basis)
        ptype_value, wfn = func(method_name, return_wfn=True, molecule=molecule, **kwargs)
        if core.get_option("SCF", "DF_INTS_IO") != "SAVE":
            core.clean()

        optstash.restore()

    else:
        ptype_value, wfn = cbs(func, total_method_name, **cbs_kwargs)

    if return_wfn:
        return (ptype_value, wfn)
    else:
        return ptype_value


def _build_cbs_compute(metameta, metadata):
    label = metameta['label']
    natom = metameta['molecule'].natom()
    ptype = metameta['ptype']

    # Build string of title banner
    cbsbanners = ''
    cbsbanners += """core.print_out('\\n')\n"""
    cbsbanners += """p4util.banner(' CBS Setup: %s ' % label)\n"""
    cbsbanners += """core.print_out('\\n')\n\n"""
    exec(cbsbanners)

    # Call schemes for each portion of total energy to 'place orders' for calculations needed
    d_fields = [
        'd_stage', 'd_scheme', 'd_basis', 'd_wfn', 'd_need', 'd_coef', 'd_energy', 'd_gradient', 'd_hessian', 'd_alpha'
    ]
    GRAND_NEED = []

    NEED = _expand_scheme_orders(metadata[0]["scheme"], metadata[0]["basis"][0], metadata[0]["basis"][1],
                                 metadata[0]["wfn"], metadata[0]["options"], natom)
    GRAND_NEED.append(
        dict(
            zip(d_fields, [
                'scf', metadata[0]["scheme"],
                _contract_bracketed_basis(metadata[0]["basis"][0]), metadata[0]["wfn"], NEED, +1, 0.0, None, None,
                metadata[0]["alpha"]
            ])))
    if len(metadata) > 1:
        for delta in metadata[1:]:
            NEED = _expand_scheme_orders(delta["scheme"], delta["basis"][0], delta["basis"][1], delta["wfn"],
                                         delta["options"], natom)
            GRAND_NEED.append(
                dict(
                    zip(d_fields, [
                        delta["stage"], delta["scheme"],
                        _contract_bracketed_basis(delta["basis"][0]), delta["wfn"], NEED, +1, 0.0, None, None,
                        delta["alpha"]
                    ])))
            NEED = _expand_scheme_orders(delta["scheme"], delta["basis_lo"][0], delta["basis_lo"][1], delta["wfn_lo"],
                                         False, natom)
            GRAND_NEED.append(
                dict(
                    zip(d_fields, [
                        delta["stage"], delta["scheme"],
                        _contract_bracketed_basis(delta["basis_lo"][0]), delta["wfn_lo"], NEED, -1, 0.0, None, None,
                        delta["alpha"]
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

    #     Remove duplicate modelchem portion listings
    for mc in MODELCHEM:
        dups = -1
        for indx_job, job in enumerate(JOBS):
            if ((job['f_wfn'] == mc['f_wfn']) and (job['f_basis'] == mc['f_basis'])
                    and (job['f_options'] == mc['f_options'])):
                dups += 1
                if dups >= 1:
                    del JOBS[indx_job]

    instructions = ''
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
                            and (mc['f_basis'] == job['f_basis']) and not (mc['f_wfn'] == job['f_wfn'])):
                        del JOBS[indx_job]

    instructions += """\n    Enlightened listing of computations required.\n"""
    for mc in JOBS:
        instructions += listfmt.format(mc['f_wfn'], mc['f_basis'] + " + options" * bool(mc['f_options']),
                                       VARH[mc['f_wfn']][mc['f_wfn']], _addlremark[ptype])
    core.print_out(instructions)

    #     Expand listings to all that will be obtained
    TROVE = []
    for job in JOBS:
        for wfn in VARH[job['f_wfn']]:
            TROVE.append(
                dict(
                    zip(_f_fields, [
                        wfn, job['f_basis'], job['f_zeta'], job['f_options'], 0.0,
                        np.zeros((natom, 3)),
                        np.zeros((3 * natom, 3 * natom))
                    ])))

    instructions = """\n    Full listing of computations to be obtained (required and bonus).\n"""
    for mc in TROVE:
        instructions += listfmt.format(mc['f_wfn'], mc['f_basis'] + " + options" * bool(mc['f_options']),
                                       VARH[mc['f_wfn']][mc['f_wfn']], _addlremark[ptype])
    core.print_out(instructions)

    return GRAND_NEED, JOBS, TROVE


def _assemble_cbs_components(metameta, TROVE, GRAND_NEED):
    """Absorb job E/G/H results from `TROVE` into `GRAND_NEED`. Process
    those into stage E/G/H in `GRAND_NEED`, returning the latter.
    Accumulate into final E/G/H quantities, returning them in dict.

    """
    label = metameta['label']
    natom = metameta['molecule'].natom()
    ptype = metameta['ptype']

    # Build string of title banner
    cbsbanners = ''
    cbsbanners += """core.print_out('\\n')\n"""
    cbsbanners += """p4util.banner(' CBS Results: %s ' % label)\n"""
    cbsbanners += """core.print_out('\\n')\n\n"""
    exec(cbsbanners)

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

    # Make xtpl() call
    finalenergy = 0.0
    finalgradient = np.zeros((natom, 3))
    finalhessian = np.zeros((3 * natom, 3 * natom))

    for stage in GRAND_NEED:
        hiloargs = {'alpha': stage['d_alpha']}

        grad_available = all([np.count_nonzero(lmh['f_gradient']) for lmh in stage['d_need'].values()])
        hess_available = all([np.count_nonzero(lmh['f_hessian']) for lmh in stage['d_need'].values()])

        hiloargs.update(_contract_scheme_orders(stage['d_need'], 'f_energy'))
        stage['d_energy'] = xtpl_procedures[stage['d_scheme']](**hiloargs)
        finalenergy += stage['d_energy'] * stage['d_coef']

        if ptype == 'gradient' or grad_available:
            hiloargs.update(_contract_scheme_orders(stage['d_need'], 'f_gradient'))
            stage['d_gradient'] = xtpl_procedures[stage['d_scheme']](**hiloargs)
            finalgradient += stage['d_gradient'] * stage['d_coef']

        if ptype == 'hessian' or hess_available:
            hiloargs.update(_contract_scheme_orders(stage['d_need'], 'f_hessian'))
            stage['d_hessian'] = xtpl_procedures[stage['d_scheme']](**hiloargs)
            finalhessian += stage['d_hessian'] * stage['d_coef']

    cbs_results = {
        'ret_ptype': {
            'energy': finalenergy,
            'gradient': finalgradient,
            'hessian': finalhessian,
        }[ptype],
        'energy': finalenergy,
        'gradient': finalgradient,
        'hessian': finalhessian,
    }

    return cbs_results, GRAND_NEED


def _summary_table(metadata, TROVE, GRAND_NEED):
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


class CBSComputer(BaseTask):

    molecule: Any
    driver: str = "energy"
    metadata: Any
    metameta: Dict[str, Any] = {}

    return_wfn: bool = False
    verbose: int = 0

    # List of model chemistries with extrapolation scheme applied. Can reconstruct CBS. Keys are d_fields.
    grand_need: List[Dict[str, Any]] = []

    # Maximal list of model chemistries extractable from running `compute_list`. Keys are _f_fields. Formerly JOBS_EXT.
    trove: List[Dict[str, Any]] = []

    # Minimal (enlightened) list of jobs to run to satisfy full CBS. Keys are _f_fields. Formerly JOBS.
    compute_list: List[Dict[str, Any]] = []

    # One-to-One list of SingleResult-s corresponding to `compute_list`.
    task_list: List[Any] = []

    # One-to-One list of QCSchema corresponding to `task_list`.
    results_list: List[Any] = []

    @pydantic.validator('driver')
    def set_driver(cls, driver):
        if driver not in ['energy', 'gradient', 'hessian']:
            raise ValidationError(f"""Wrapper complete_basis_set is unhappy to be calling
                                     function '{driver}' instead of 'energy', 'gradient' or 'hessian'.""")
        return driver

    @pydantic.validator('molecule')
    def set_molecule(cls, mol):
        mol.update_geometry()
        mol.fix_com(True)
        mol.fix_orientation(True)
        return mol

    def __init__(self, **data):
        data = p4util.kwargs_lower(data)
        data["metadata"] = _process_cbs_kwargs(data)
        BaseTask.__init__(self, **data)

        self.metameta = {
            'kwargs': data,
            'ptype': self.driver,
            'verbose': self.verbose,
            'label': None,
            'molecule': self.molecule,
        }
        print('METAMETA')
        import pprint
        pprint.pprint(self.metameta)

        if data['metadata']:
            print('DATAMETADATA', data['metadata'])
            if data['metadata'][0]["wfn"] not in VARH.keys():
                raise ValidationError(
                    """Requested SCF method '%s' is not recognized. Add it to VARH in driver_cbs.py to proceed.""" %
                    (metadata[0]["wfn"]))

            if len(self.metadata) > 1:
                for delta in self.metadata[1:]:
                    if delta["wfn"] not in VARH.keys():
                        raise ValidationError(
                            """Requested higher %s method '%s' is not recognized. Add it to VARH in driver_cbs.py to proceed."""
                            % (delta["treatment"], delta["wfn"]))
                    if delta["wfn_lo"] not in VARH.keys():
                        raise ValidationError(
                            """Requested lesser %s method '%s' is not recognized. Add it to VARH in driver_cbs.py to proceed."""
                            % (delta["treatment"], delta["wfn_lo"]))

            self.grand_need, self.compute_list, self.trove = _build_cbs_compute(self.metameta, self.metadata)

            for job in self.compute_list:
                keywords = copy.deepcopy(self.metameta['kwargs']['keywords'])
                if job["f_options"] is not False:
                    stage_keywords = dict(job["f_options"].items())
                    print('STAGE KEY', stage_keywords)
                    keywords = {**keywords, **stage_keywords}
                task = SingleResult(
                    **{
                        "molecule": self.molecule,
                        "driver": self.driver,
                        "method": job["f_wfn"],
                        "basis": job["f_basis"],
                        "keywords": keywords or {},
                    })
                print(task)
                self.task_list.append(task)

    def build_tasks(self, obj, **kwargs):
        # permanently a dummy function
        pass

    def plan(self):
        return [x.plan() for x in self.task_list]

    def compute(self):
        all_options = p4util.prepare_options_for_modules(changedOnly=True, commandsInsteadDict=False)
        # gof = core.get_output_file()
        # core.close_outfile()

        for x in self.task_list:
            x.compute()

        # core.set_output_file(gof, True)
        p4util.reset_pe_options(all_options)

    def _prepare_results(self):
        results_list = [x.get_results() for x in self.task_list]
        print(results_list)
        energies = [x['properties']["return_energy"] for x in results_list]
        print('ENE', energies)
        for x in self.task_list:
            print('\nTASK')
            pprint.pprint(x)
        for x in self.results_list:
            print('\nRESULT')
            pprint.pprint(x)

        # load results_list numbers into compute_list (task_list is SingleResult-s)
        for itask, mc in enumerate(self.compute_list):
            task = results_list[itask]
            response = task['return_result']
            print('ITASK:', itask, '\nRETURN', response, '\nTASK:')
            pprint.pprint(task)

            if self.metameta['ptype'] == 'energy':
                mc['f_energy'] = response

            elif self.metameta['ptype'] == 'gradient':
                mc['f_gradient'] = np.array(response).reshape((-1, 3))
                mc['f_energy'] = task['psi4:qcvars']['CURRENT ENERGY']

            elif self.metameta['ptype'] == 'hessian':
                mc['f_hessian'] = plump_qcvar(response, 'hessian')
                mc['f_energy'] = task['psi4:qcvars']['CURRENT ENERGY']
                if 'CURRENT GRADIENT' in task['psi4:qcvars']:
                    mc['f_gradient'] = plump_qcvar(task['psi4:qcvars']['CURRENT GRADIENT'], 'gradient')

            # Fill in energies for subsumed methods
            if self.metameta['ptype'] == 'energy':
                for wfn in VARH[mc['f_wfn']]:
                    for job in self.trove:
                        if ((wfn == job['f_wfn']) and (mc['f_basis'] == job['f_basis'])
                                and (mc['f_options'] == job['f_options'])):
                            job['f_energy'] = task['psi4:qcvars'][VARH[wfn][wfn]]

            # Copy data from 'run' to 'obtained' table
            for mce in self.trove:
                if ((mc['f_wfn'] == mce['f_wfn']) and (mc['f_basis'] == mce['f_basis'])
                        and (mc['f_options'] == mce['f_options'])):
                    mce['f_energy'] = mc['f_energy']
                    mce['f_gradient'] = mc['f_gradient']
                    mce['f_hessian'] = mc['f_hessian']

            print('MC:\n', mc)

        cbs_results, self.grand_need = _assemble_cbs_components(self.metameta, self.trove, self.grand_need)

        core.print_out(_summary_table(self.metadata, self.trove, self.grand_need))

        print('\nCBS_RESULTS', cbs_results)
        print('\nGRAND_NEED', self.grand_need)

        return cbs_results

    def get_results(self):
        """Return CBS results as CBS-flavored QCSchema."""

        assembled_results = self._prepare_results()

        # load QCVariables
        qcvars = {
            'CBS NUMBER': len(self.compute_list),
            'NUCLEAR REPULSION ENERGY': self.molecule.nuclear_repulsion_energy(),
        }

        for qcv in ['CBS', 'CURRENT']:
            qcvars[qcv + ' REFERENCE ENERGY'] = self.grand_need[0]['d_energy']
            qcvars[qcv + ' CORRELATION ENERGY'] = assembled_results['energy'] - self.grand_need[0]['d_energy']
            qcvars[qcv + ('' if qcv == 'CURRENT' else ' TOTAL') + ' ENERGY'] = assembled_results['energy']

        for idelta in range(int(len(cbs.grand_need) / 2)):
            if idelta == 0:
                continue
            dc = idelta * 2 + 1
            qcvars[f"CBS {cbs.grand_need[dc]['d_stage'].upper()} TOTAL ENERGY"] = self.grand_need[dc]["d_energy"] - self.grand_need[dc + 1]["d_energy"]

        if np.count_nonzero(assembled_results['gradient']):
            for qcv in ['CURRENT GRADIENT', 'CBS TOTAL GRADIENT']:
                qcvars[qcv] = assembled_results['gradient']

        if np.count_nonzero(assembled_results['hessian']):
            for qcv in ['CURRENT HESSIAN', 'CBS TOTAL HESSIAN']:
                qcvars[qcv] = assembled_results['hessian']

        cbsrec = {
            #'cbs_jobs': copy.deepcopy(self.compute_list),
            'cbs_record': copy.deepcopy(self.grand_need),
            'driver': self.driver,
            #'keywords': self.keywords,
            #'model': {
            #    'method': self.method,
            #    'basis': self.basis
            #},
            'molecule': self.molecule.to_schema(dtype=1)['molecule'],
            'properties': {
                #'calcinfo_nalpha': wfn.nalpha(),
                #'calcinfo_nbeta': wfn.nbeta(),
                'calcinfo_natom': self.molecule.natom(),
                'nuclear_repulsion_energy': self.molecule.nuclear_repulsion_energy(),
                'return_energy': assembled_results['energy'],
            },
            'provenance': p4util.provenance_stamp(__name__),
            'psi4:qcvars': qcvars,
            #'raw_output': None,
            'return_result': assembled_results['ret_ptype'],
            'schema_name': 'qc_schema_output',
            'schema_version': 1,
            'success': True,
        }

        cbsrec = unnp(cbsrec, flat=True)
        print('\nCBS QCSchema:')
        pprint.pprint(cbsrec)
        return cbsrec

    def get_psi_results(self, return_wfn=False):
        """Return CBS results in usual E/wfn interface."""

        cbs_schema = self.get_results()

        ret_ptype = plump_qcvar(cbs_schema['return_result'], shape_clue=cbs_schema['driver'], ret='psi4')
        wfn = self._cbs_schema_to_wfn(cbs_schema)

        if return_wfn:
            return (ret_ptype, wfn)
        else:
            return ret_ptype

    @staticmethod
    def _cbs_schema_to_wfn(cbs_schema):
        """Helper function to keep Wavefunction dependent on CBS-flavored QCSchemus."""

        # new skeleton wavefunction w/mol, highest-SCF basis (just to choose one), & not energy
        mol = core.Molecule.from_schema(cbs_schema)
        basis = core.BasisSet.build(mol, "ORBITAL", 'def2-svp')
        wfn = core.Wavefunction(mol, basis)

        for qcv, val in cbs_schema['psi4:qcvars'].items():
            for obj in [core, wfn]:
                obj.set_variable(qcv, plump_qcvar(val, qcv))

        flat_grad = cbs_schema['psi4:qcvars'].get('CBS TOTAL GRADIENT')
        if flat_grad is not None:
            finalgradient = plump_qcvar(flat_grad, 'gradient', ret='psi4')
            wfn.set_gradient(finalgradient)

            if finalgradient.rows(0) < 20:
                core.print_out('CURRENT GRADIENT')
                finalgradient.print_out()

        flat_hess = cbs_schema['psi4:qcvars'].get('CBS TOTAL HESSIAN')
        if flat_hess is not None:
            finalhessian = plump_qcvar(flat_hess, 'hessian', ret='psi4')
            wfn.set_hessian(finalhessian)

            if finalhessian.rows(0) < 20:
                core.print_out('CURRENT HESSIAN')
                finalhessian.print_out()

        return wfn
