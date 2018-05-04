#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2018 The Psi4 Developers.
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
"""Module with utility function for dumping the Hamiltonian to file in FCIDUMP format."""
from __future__ import division

import mmap
import re
from datetime import datetime

import numpy as np
from deepdiff import DeepDiff

from psi4.driver import constants
from psi4.driver.p4util.util import success
from psi4.driver.procrouting.proc_util import check_iwl_file_from_scf_type

from .exceptions import *


def fcidump(wfn, fname='INTDUMP', oe_ints=None):
    """Save integrals to file in FCIDUMP format as defined in Comp. Phys. Commun. 54 75 (1989)
    Additional one-electron integrals, including orbital energies, can also be saved.
    This latter format can be used with the HANDE QMC code but is not standard.

    :returns: None

    :raises: ValidationError when SCF wavefunction is not RHF

    :type wfn: :py:class:`~psi4.core.Wavefunction`
    :param wfn: set of molecule, basis, orbitals from which to generate cube files
    :param fname: name of the integrals file, defaults to INTDUMP
    :param oe_ints: list of additional one-electron integrals to save to file.
    So far only EIGENVALUES is a valid option.

    :examples:

    >>> # [1] Save one- and two-electron integrals to standard FCIDUMP format
    >>> E, wfn = energy('scf', return_wfn=True)
    >>> fcidump(wfn)

    >>> # [2] Save orbital energies, one- and two-electron integrals.
    >>> E, wfn = energy('scf', return_wfn=True)
    >>> fcidump(wfn, oe_ints=['EIGENVALUES'])

    """
    # Get some options
    reference = core.get_option('SCF', 'REFERENCE')
    ints_tolerance = core.get_global_option('INTS_TOLERANCE')
    # Some sanity checks
    if reference not in ['RHF', 'UHF']:
        raise ValidationError('FCIDUMP not implemented for {} references\n'.format(reference))
    if oe_ints is None:
        oe_ints = []

    molecule = wfn.molecule()
    docc = wfn.doccpi()
    frzcpi = wfn.frzcpi()
    frzvpi = wfn.frzvpi()
    active_docc = docc - frzcpi
    active_socc = wfn.soccpi()
    active_mopi = wfn.nmopi() - frzcpi - frzvpi

    nbf = active_mopi.sum() if wfn.same_a_b_orbs() else 2 * active_mopi.sum()
    nirrep = wfn.nirrep()
    nelectron = 2 * active_docc.sum() + active_socc.sum()

    core.print_out('Writing integrals in FCIDUMP format to ' + fname + '\n')
    # Generate FCIDUMP header
    header = '&FCI\n'
    header += 'NORB={:d},\n'.format(nbf)
    header += 'NELEC={:d},\n'.format(nelectron)
    header += 'MS2={:d},\n'.format(wfn.nalpha() - wfn.nbeta())
    header += 'UHF=.{}.,\n'.format(not wfn.same_a_b_orbs()).upper()
    orbsym = ''
    for h in range(active_mopi.n()):
        for n in range(frzcpi[h], frzcpi[h] + active_mopi[h]):
            orbsym += '{:d},'.format(h + 1)
            if not wfn.same_a_b_orbs():
                orbsym += '{:d},'.format(h + 1)
    header += 'ORBSYM={}\n'.format(orbsym)
    header += '&END\n'
    with open(fname, 'w') as intdump:
        intdump.write(header)

    # Get an IntegralTransform object
    check_iwl_file_from_scf_type(core.get_option('SCF', 'SCF_TYPE'), wfn)
    spaces = [core.MOSpace.all()]
    trans_type = core.IntegralTransform.TransformationType.Restricted
    if not wfn.same_a_b_orbs():
        trans_type = core.IntegralTransform.TransformationType.Unrestricted
    ints = core.IntegralTransform(wfn, spaces, trans_type)
    ints.transform_tei(core.MOSpace.all(), core.MOSpace.all(), core.MOSpace.all(), core.MOSpace.all())
    core.print_out('Integral transformation complete!\n')

    DPD_info = {'instance_id': ints.get_dpd_id(), 'alpha_MO': ints.DPD_ID('[A>=A]+'), 'beta_MO': 0}
    if not wfn.same_a_b_orbs():
        DPD_info['beta_MO'] = ints.DPD_ID("[a>=a]+")
    # Write TEI to fname in FCIDUMP format
    core.fcidump_tei_helper(nirrep, wfn.same_a_b_orbs(), DPD_info, ints_tolerance, fname)

    # Read-in OEI and write them to fname in FCIDUMP format
    # Indexing functions to translate from zero-based (C and Python) to
    # one-based (Fortran)
    mo_idx = lambda x: x + 1
    alpha_mo_idx = lambda x: 2 * x + 1
    beta_mo_idx = lambda x: 2 * (x + 1)

    with open(fname, 'a') as intdump:
        core.print_out('Writing frozen core operator in FCIDUMP format to ' + fname + '\n')
        if reference == 'RHF':
            PSIF_MO_FZC = 'MO-basis Frozen-Core Operator'
            moH = core.Matrix(PSIF_MO_FZC, wfn.nmopi(), wfn.nmopi())
            moH.load(core.IO.shared_object(), constants.PSIF_OEI)
            mo_slice = core.Slice(frzcpi, active_mopi)
            MO_FZC = moH.get_block(mo_slice, mo_slice)
            offset = 0
            for h, block in enumerate(MO_FZC.nph):
                il = np.tril_indices(block.shape[0])
                for index, x in np.ndenumerate(block[il]):
                    row = mo_idx(il[0][index] + offset)
                    col = mo_idx(il[1][index] + offset)
                    if (abs(x) > ints_tolerance):
                        intdump.write('{:29.20E} {:4d} {:4d} {:4d} {:4d}\n'.format(x, row, col, 0, 0))
                offset += block.shape[0]
            # Additional one-electron integrals as requested in oe_ints
            # Orbital energies
            core.print_out('Writing orbital energies in FCIDUMP format to ' + fname + '\n')
            if 'EIGENVALUES' in oe_ints:
                eigs_dump = write_eigenvalues(wfn.epsilon_a().get_block(mo_slice).to_array(), mo_idx)
                intdump.write(eigs_dump)
        else:
            PSIF_MO_A_FZC = 'MO-basis Alpha Frozen-Core Oper'
            moH_A = core.Matrix(PSIF_MO_A_FZC, wfn.nmopi(), wfn.nmopi())
            moH_A.load(core.IO.shared_object(), constants.PSIF_OEI)
            mo_slice = core.Slice(frzcpi, active_mopi)
            MO_FZC_A = moH_A.get_block(mo_slice, mo_slice)
            offset = 0
            for h, block in enumerate(MO_FZC_A.nph):
                il = np.tril_indices(block.shape[0])
                for index, x in np.ndenumerate(block[il]):
                    row = alpha_mo_idx(il[0][index] + offset)
                    col = alpha_mo_idx(il[1][index] + offset)
                    if (abs(x) > ints_tolerance):
                        intdump.write('{:29.20E} {:4d} {:4d} {:4d} {:4d}\n'.format(x, row, col, 0, 0))
                offset += block.shape[0]
            PSIF_MO_B_FZC = 'MO-basis Beta Frozen-Core Oper'
            moH_B = core.Matrix(PSIF_MO_B_FZC, wfn.nmopi(), wfn.nmopi())
            moH_B.load(core.IO.shared_object(), constants.PSIF_OEI)
            mo_slice = core.Slice(frzcpi, active_mopi)
            MO_FZC_B = moH_B.get_block(mo_slice, mo_slice)
            offset = 0
            for h, block in enumerate(MO_FZC_B.nph):
                il = np.tril_indices(block.shape[0])
                for index, x in np.ndenumerate(block[il]):
                    row = beta_mo_idx(il[0][index] + offset)
                    col = beta_mo_idx(il[1][index] + offset)
                    if (abs(x) > ints_tolerance):
                        intdump.write('{:29.20E} {:4d} {:4d} {:4d} {:4d}\n'.format(x, row, col, 0, 0))
                offset += block.shape[0]
            # Additional one-electron integrals as requested in oe_ints
            # Orbital energies
            core.print_out('Writing orbital energies in FCIDUMP format to ' + fname + '\n')
            if 'EIGENVALUES' in oe_ints:
                alpha_eigs_dump = write_eigenvalues(wfn.epsilon_a().get_block(mo_slice).to_array(), alpha_mo_idx)
                beta_eigs_dump = write_eigenvalues(wfn.epsilon_b().get_block(mo_slice).to_array(), beta_mo_idx)
                intdump.write(alpha_eigs_dump + beta_eigs_dump)
        # Dipole integrals
        #core.print_out('Writing dipole moment OEI in FCIDUMP format to ' + fname + '\n')
        # Traceless quadrupole integrals
        #core.print_out('Writing traceless quadrupole moment OEI in FCIDUMP format to ' + fname + '\n')
        # Frozen core + nuclear repulsion energy
        core.print_out('Writing frozen core + nuclear repulsion energy in FCIDUMP format to ' + fname + '\n')
        e_fzc = ints.get_frozen_core_energy()
        e_nuc = molecule.nuclear_repulsion_energy(wfn.get_dipole_field_strength())
        intdump.write('{: 29.20E} {:4d} {:4d} {:4d} {:4d}\n'.format(e_fzc + e_nuc, 0, 0, 0, 0))
    core.print_out('Done generating {} with integrals in FCIDUMP format.\n'.format(fname))


def write_eigenvalues(eigs, mo_idx):
    """Prepare multi-line string with one-particle eigenvalues to be written to the FCIDUMP file.
    """
    eigs_dump = ''
    iorb = 0
    for h, block in enumerate(eigs):
        for idx, x in np.ndenumerate(block):
            eigs_dump += '{: 29.20E} {:4d} {:4d} {:4d} {:4d}\n'.format(x, mo_idx(iorb), 0, 0, 0)
            iorb += 1
    return eigs_dump


def readin_fcidump(fname):
    """Function to read in a FCIDUMP file.

    :returns: a tuple with the header and the data in the FCIDUMP file

    :param fname: FCIDUMP file name
    """
    with open(fname, 'r+') as f:
        data = mmap.mmap(f.fileno(), 0)
        # Get the string in between &FCI and &END.
        header_string = re.match(b'(?:&FCI((?:.*?\r?\n?)*)&END)', data).group(1)
    header = _parse_fcidump_header(header_string.decode('utf-8'))
    integrals = _parse_fcidump_integrals(fname)
    return header, integrals


def compare_fcidump_headers(expected, computed, label):
    """Function to compare the headers of two FCIDUMP files. Prints :py:func:`util.success`
    when value *computed* matches value *expected*.
    Performs a system exit on failure. Used in input files in the test suite.

    :returns: a tuple of energies computed from the MO integrals (SCF 1-el, SCF 2-el, SCF total, MP2 correlation)

    :param expected: reference FCIDUMP file
    :param computed: computed FCIDUMP file
    :param label: string labelling the test
    """
    # Grab expected header and integrals
    ref_header, ref_intdump = readin_fcidump(expected)
    header, intdump = readin_fcidump(computed)

    # Compare headers
    header_diff = DeepDiff(ref_header, header, ignore_order=True)
    if header_diff:
        message = ("\tComputed FCIDUMP file header does not match expected header.\n")
        raise TestComparisonError(header_diff)
    success(label)

    nbf = header['norb']
    # Skip the last line, i.e. Enuc + Efzc
    ints = intdump[:-1, 0]
    # Get indices
    idxs = intdump[:, 1:].astype(np.int) - 1
    # Slices
    sl = slice(ints.shape[0] - nbf, ints.shape[0])
    # Extract orbital energies
    epsilon = np.zeros(nbf)
    epsilon[idxs[sl, 0]] = ints[sl]
    # Count how many 2-index integrals we have
    sl = slice(sl.start - nbf * nbf, sl.stop - nbf)
    two_index = np.all(idxs[sl, 2:] == -1, axis=1).sum()
    sl = slice(sl.stop - two_index, sl.stop)
    # Extract Hcore
    Hcore = np.zeros((nbf, nbf))
    Hcore[(idxs[sl, 0], idxs[sl, 1])] = ints[sl]
    Hcore[(idxs[sl, 1], idxs[sl, 0])] = ints[sl]
    # Extract ERIs
    sl = slice(0, sl.start)
    eri = np.zeros((nbf, nbf, nbf, nbf))
    eri[(idxs[sl, 0], idxs[sl, 1], idxs[sl, 2], idxs[sl, 3])] = ints[sl]
    eri[(idxs[sl, 0], idxs[sl, 1], idxs[sl, 3], idxs[sl, 2])] = ints[sl]
    eri[(idxs[sl, 1], idxs[sl, 0], idxs[sl, 2], idxs[sl, 3])] = ints[sl]
    eri[(idxs[sl, 1], idxs[sl, 0], idxs[sl, 3], idxs[sl, 2])] = ints[sl]
    eri[(idxs[sl, 2], idxs[sl, 3], idxs[sl, 0], idxs[sl, 1])] = ints[sl]
    eri[(idxs[sl, 3], idxs[sl, 2], idxs[sl, 0], idxs[sl, 1])] = ints[sl]
    eri[(idxs[sl, 2], idxs[sl, 3], idxs[sl, 1], idxs[sl, 0])] = ints[sl]
    eri[(idxs[sl, 3], idxs[sl, 2], idxs[sl, 1], idxs[sl, 0])] = ints[sl]
    # Compute SCF energy
    scf_1el_e, scf_2el_e = _scf_energy_from_fcidump(Hcore, eri, np.where(epsilon < 0)[0], header['uhf'])
    # Compute MP2 energy
    mp2_e = _mp2_energy_from_fcidump(eri, epsilon, header['uhf'])
    return scf_1el_e, scf_2el_e, (scf_1el_e + scf_2el_e), mp2_e


def _scf_energy_from_fcidump(Hcore, ERI, occ_sl, unrestricted):
    scf_1el_e = np.einsum('ii->', Hcore[np.ix_(occ_sl, occ_sl)])
    if not unrestricted:
        scf_1el_e *= 2
    coulomb = np.einsum('iijj->', ERI[np.ix_(occ_sl, occ_sl, occ_sl, occ_sl)])
    exchange = np.einsum('ijij->', ERI[np.ix_(occ_sl, occ_sl, occ_sl, occ_sl)])
    if unrestricted:
        scf_2el_e = 0.5 * (coulomb - exchange)
    else:
        scf_2el_e = 2.0 * coulomb - exchange
    return scf_1el_e, scf_2el_e


def _mp2_energy_from_fcidump(ERI, epsilon, unrestricted):
    # Occupied and virtual slices
    occ_sl = np.where(epsilon < 0)[0]
    vir_sl = np.where(epsilon > 0)[0]
    eocc = epsilon[occ_sl]
    evir = epsilon[vir_sl]
    denom = 1 / (eocc.reshape(-1, 1, 1, 1) - evir.reshape(-1, 1, 1) + eocc.reshape(-1, 1) - evir)
    MO = ERI[np.ix_(occ_sl, vir_sl, occ_sl, vir_sl)]
    if unrestricted:
        mp2_e = 0.5 * np.einsum("abrs,abrs,abrs->", MO, MO - MO.swapaxes(1, 3), denom)
    else:
        mp2_e = np.einsum('iajb,iajb,iajb->', MO, MO, denom) + np.einsum('iajb,iajb,iajb->', MO - MO.swapaxes(1, 3),
                                                                         MO, denom)
    return mp2_e


def _parse_fcidump_header(header_string):
    """Parse FCIDUMP preamble, using regexes.

    :returns: a dictionary with the contents of the header

    :param fname: the multiline header string
    """
    # Set up a dictionary of regexes for the stuff we want.
    # The value in the (key, value) pair is a tuple.
    # The first element is the regex, the second is a function to process the result of the match.
    # yapf: disable
    parser = {
        # Get NORB, transform to integer
          'norb' : ('(?:NORB=(\d+))', int)
        # Get NELEC, transform to integer
        , 'nelec' : ('(?:NELEC=(\d+))', int)
        # Get MS2, transform to integer
        , 'ms2' : ('(?:MS2=(\d+))', int)
        # Get ISYM, if present, transform to integer
        , 'isym' : ('(?:ISYM=(\d+))?', (lambda x: int(x) if x is not None else 1))
        # Get ORBSYM, transform to list of integers
        , 'orbsym' : ('(?:ORBSYM=(\d.+))', (lambda x: [int(y.strip()) for y in x.rstrip(',').split(',')]))
        # Get UHF, transform to True or False
        , 'uhf' : ('(?:UHF=((?:.TRUE.|.FALSE.)))', (lambda x: True if x.lower() == '.true.' else False))
    }
    # yapf: enable
    # Generate preamble dictionary.
    # The keys are shared with the parser dictionary, the values result
    # from searching, case-insensitively, with the corresponding regex (v[0])
    # and applying the processing function (v[1]) on the result.
    preamble = {k: v[1](re.search(v[0], header_string, re.I).group(1)) for k, v in parser.items()}
    return preamble


def _parse_fcidump_integrals(fname):
    """Parse FCIDUMP integrals, using a regex.

    :returns: an array with the integrals and their indices as floats, i.e. the body of the FCIDUMP

    :param fname: FCIDUMP file name
    """
    regex = r'(-?[0-9]+\.?[0-9]*(?:[Ee][+-]?[0-9]+)?)\s*(\d+)\s*(\d+)\s*(\d+)\s*(\d+)'
    return np.fromregex(fname, regex, dtype=float)
