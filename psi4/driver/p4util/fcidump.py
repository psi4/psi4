#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2023 The Psi4 Developers.
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

__all__ = [
    "compare_fcidumps",
    "energies_from_fcidump",
    "fcidump",
    "fcidump_from_file",
]

from typing import Any, Dict, List, Optional

import numpy as np

from psi4.driver import psifiles as psif
from psi4.driver.p4util.testing import compare_integers, compare_values, compare_recursive
from psi4.driver.procrouting.proc_util import check_iwl_file_from_scf_type

from psi4 import core
from .exceptions import ValidationError, TestComparisonError


def fcidump(wfn: core.Wavefunction, fname: str = 'INTDUMP', oe_ints: Optional[List] = None):
    """Save integrals to file in FCIDUMP format as defined in Comp. Phys. Commun. 54 75 (1989),
    https://doi.org/10.1016/0010-4655(89)90033-7 .
    Additional one-electron integrals, including orbital energies, can also be saved.
    This latter format can be used with the HANDE QMC code but is not standard.

    Parameters
    ----------
    wfn
        Set of molecule, basis, orbitals from which to generate FCIDUMP file.
    fname
        Name of the integrals file, defaults to INTDUMP.
    oe_ints
        List of additional one-electron integrals to save to file. So far only
        EIGENVALUES is a valid option.

    Raises
    ------
    ValidationError
        When SCF wavefunction is not RHF.

    Examples
    --------

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
    irrep_map = _irrep_map(wfn)

    wfn_irrep = 0
    for h, n_socc in enumerate(active_socc):
        if n_socc % 2 == 1:
            wfn_irrep ^= h

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
            orbsym += '{:d},'.format(irrep_map[h])
            if not wfn.same_a_b_orbs():
                orbsym += '{:d},'.format(irrep_map[h])
    header += 'ORBSYM={}\n'.format(orbsym)
    header += 'ISYM={:d},\n'.format(irrep_map[wfn_irrep])
    header += '&END\n'
    with open(fname, 'w') as intdump:
        intdump.write(header)

    # Get an IntegralTransform object
    check_iwl_file_from_scf_type(core.get_global_option('SCF_TYPE'), wfn)
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
            moH.load(core.IO.shared_object(), psif.PSIF_OEI)
            mo_slice = core.Slice(frzcpi, frzcpi+active_mopi)
            MO_FZC = moH.get_block(mo_slice, mo_slice)
            offset = 0
            for h, block in enumerate(MO_FZC.nph):
                il = np.tril_indices(block.shape[0])
                for index, x in np.ndenumerate(block[il]):
                    row = mo_idx(il[0][index] + offset)
                    col = mo_idx(il[1][index] + offset)
                    if (abs(x) > ints_tolerance):
                        intdump.write('{:28.20E}{:4d}{:4d}{:4d}{:4d}\n'.format(x, row, col, 0, 0))
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
            moH_A.load(core.IO.shared_object(), psif.PSIF_OEI)
            mo_slice = core.Slice(frzcpi, active_mopi)
            MO_FZC_A = moH_A.get_block(mo_slice, mo_slice)
            offset = 0
            for h, block in enumerate(MO_FZC_A.nph):
                il = np.tril_indices(block.shape[0])
                for index, x in np.ndenumerate(block[il]):
                    row = alpha_mo_idx(il[0][index] + offset)
                    col = alpha_mo_idx(il[1][index] + offset)
                    if (abs(x) > ints_tolerance):
                        intdump.write('{:28.20E}{:4d}{:4d}{:4d}{:4d}\n'.format(x, row, col, 0, 0))
                offset += block.shape[0]
            PSIF_MO_B_FZC = 'MO-basis Beta Frozen-Core Oper'
            moH_B = core.Matrix(PSIF_MO_B_FZC, wfn.nmopi(), wfn.nmopi())
            moH_B.load(core.IO.shared_object(), psif.PSIF_OEI)
            mo_slice = core.Slice(frzcpi, active_mopi)
            MO_FZC_B = moH_B.get_block(mo_slice, mo_slice)
            offset = 0
            for h, block in enumerate(MO_FZC_B.nph):
                il = np.tril_indices(block.shape[0])
                for index, x in np.ndenumerate(block[il]):
                    row = beta_mo_idx(il[0][index] + offset)
                    col = beta_mo_idx(il[1][index] + offset)
                    if (abs(x) > ints_tolerance):
                        intdump.write('{:28.20E}{:4d}{:4d}{:4d}{:4d}\n'.format(x, row, col, 0, 0))
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
        intdump.write('{:28.20E}{:4d}{:4d}{:4d}{:4d}\n'.format(e_fzc + e_nuc, 0, 0, 0, 0))
    core.print_out('Done generating {} with integrals in FCIDUMP format.\n'.format(fname))


def write_eigenvalues(eigs, mo_idx):
    """Prepare multi-line string with one-particle eigenvalues to be written to the FCIDUMP file.
    """
    eigs_dump = ''
    iorb = 0
    for h, block in enumerate(eigs):
        for idx, x in np.ndenumerate(block):
            eigs_dump += '{:28.20E}{:4d}{:4d}{:4d}{:4d}\n'.format(x, mo_idx(iorb), 0, 0, 0)
            iorb += 1
    return eigs_dump


def _irrep_map(wfn):
    """Returns an array of irrep indices that maps from Psi4's ordering convention to the standard FCIDUMP convention.
    """
    symm = wfn.molecule().point_group().symbol()
    psi2dump = {'c1' : [1],               # A
                'ci' : [1,2],             # Ag Au
                'c2' : [1,2],             # A  B
                'cs' : [1,2],             # A' A"
                'd2' : [1,4,3,2],         # A  B1  B2  B3
                'c2v' : [1,4,2,3],        # A1 A2  B1  B2
                'c2h' : [1,4,2,3],        # Ag Bg  Au  Bu
                'd2h' : [1,4,6,7,8,5,3,2] # Ag B1g B2g B3g Au B1u B2u B3u
                }

    irrep_map = psi2dump[symm]
    return np.array(irrep_map, dtype='int')


def fcidump_from_file(fname: str) -> Dict[str, Any]:
    """Function to read in a FCIDUMP file.

    :returns: a dictionary with FCIDUMP header and integrals

      - 'norb' : number of basis functions
      - 'nelec' : number of electrons
      - 'ms2' : spin polarization of the system
      - 'isym' : symmetry of state (if present in FCIDUMP)
      - 'orbsym' : list of symmetry labels of each orbital
      - 'uhf' : whether restricted or unrestricted
      - 'enuc' : nuclear repulsion plus frozen core energy
      - 'epsilon' : orbital energies
      - 'hcore' : core Hamiltonian
      - 'eri' : electron-repulsion integrals

    :param fname: FCIDUMP file name

    """
    intdump = {}
    with open(fname, 'r') as handle:
        firstline = handle.readline().strip()
        assert '&FCI' == firstline, firstline

        skiplines = 1
        read = True
        while True:
            skiplines += 1
            line = handle.readline()
            if 'END' in line:
                break

            key, value = line.split('=')
            value = value.strip().rstrip(',')
            if key == 'UHF':
                value = 'TRUE' in value
            elif key == 'ORBSYM':
                value = [int(x) for x in value.split(',')]
            else:
                value = int(value.replace(',', ''))

            intdump[key.lower()] = value

    # Read the data and index, skip header
    raw_ints = np.genfromtxt(fname, skip_header=skiplines)

    # Read last line, i.e. Enuc + Efzc
    intdump['enuc'] = raw_ints[-1, 0]

    # Read in integrals and indices
    ints = raw_ints[:-1, 0]

    # Get dimensions and indices
    nbf = intdump['norb']
    idxs = raw_ints[:, 1:].astype(int) - 1

    # Slices
    sl = slice(ints.shape[0] - nbf, ints.shape[0])

    # Extract orbital energies
    epsilon = np.zeros(nbf)
    epsilon[idxs[sl, 0]] = ints[sl]
    intdump['epsilon'] = epsilon

    # Count how many 2-index intdump we have
    sl = slice(sl.start - nbf * nbf, sl.stop - nbf)
    two_index = np.all(idxs[sl, 2:] == -1, axis=1).sum()
    sl = slice(sl.stop - two_index, sl.stop)

    # Extract Hcore
    Hcore = np.zeros((nbf, nbf))
    Hcore[(idxs[sl, 0], idxs[sl, 1])] = ints[sl]
    Hcore[(idxs[sl, 1], idxs[sl, 0])] = ints[sl]
    intdump['hcore'] = Hcore

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
    intdump['eri'] = eri

    return intdump


def compare_fcidumps(expected: str, computed: str, label: str):
    """Comparison function for FCIDUMP files.
    Compares the first six below, then computes energies from MO integrals and compares the last four.

      - 'norb' : number of basis functions
      - 'nelec' : number of electrons
      - 'ms2' : spin polarization of the system
      - 'isym' : symmetry of state (if present in FCIDUMP)
      - 'orbsym' : list of symmetry labels of each orbital
      - 'uhf' : whether restricted or unrestricted
      - 'ONE-ELECTRON ENERGY' : SCF one-electron energy
      - 'TWO-ELECTRON ENERGY' : SCF two-electron energy
      - 'SCF TOTAL ENERGY' : SCF total energy
      - 'MP2 CORRELATION ENERGY' : MP2 correlation energy

    :param expected: Reference FCIDUMP file against which `computed` is compared.
    :param computed: Input FCIDUMP file to compare against `expected`.
    :param label: string labeling the test
    """

    # Grab expected header and integrals
    ref_intdump = fcidump_from_file(expected)
    intdump = fcidump_from_file(computed)

    # Compare headers
    compare_recursive(
        ref_intdump,
        intdump,
        'FCIDUMP header',
        forgive=['enuc', 'hcore', 'eri', 'epsilon'])

    ref_energies = energies_from_fcidump(ref_intdump)
    energies = energies_from_fcidump(intdump)

    pass_1el = compare_values(ref_energies['ONE-ELECTRON ENERGY'], energies['ONE-ELECTRON ENERGY'], 7,
                              label + '. 1-electron energy')
    pass_2el = compare_values(ref_energies['TWO-ELECTRON ENERGY'], energies['TWO-ELECTRON ENERGY'], 7,
                              label + '. 2-electron energy')
    pass_scf = compare_values(ref_energies['SCF TOTAL ENERGY'], energies['SCF TOTAL ENERGY'], 10,
                              label + '. SCF total energy')
    pass_mp2 = compare_values(ref_energies['MP2 CORRELATION ENERGY'], energies['MP2 CORRELATION ENERGY'], 10,
                              label + '. MP2 correlation energy')

    compare_integers(True, (pass_1el and pass_2el and pass_scf and pass_mp2), label)


def energies_from_fcidump(intdump) -> Dict[str, float]:
    """From integrals dictionary generated from :py:func:`fcidump_from_file`,
    compute energies.

    :returns: a dictionary with energies

      - 'NUCLEAR REPULSION ENERGY'
      - 'ONE-ELECTRON ENERGY'
      - 'TWO-ELECTRON ENERGY'
      - 'SCF TOTAL ENERGY'
      - 'MP2 CORRELATION ENERGY'

    """
    energies = {}
    energies['NUCLEAR REPULSION ENERGY'] = intdump['enuc']
    epsilon = intdump['epsilon']
    Hcore = intdump['hcore']
    eri = intdump['eri']

    # Compute SCF energy
    energies['ONE-ELECTRON ENERGY'], energies['TWO-ELECTRON ENERGY'] = _scf_energy(Hcore, eri,
                                                                                   np.where(epsilon < 0)[0],
                                                                                   intdump['uhf'])
    # yapf: disable
    energies['SCF TOTAL ENERGY'] = energies['ONE-ELECTRON ENERGY'] + energies['TWO-ELECTRON ENERGY'] + energies['NUCLEAR REPULSION ENERGY']
    # yapf: enable

    # Compute MP2 energy
    energies['MP2 CORRELATION ENERGY'] = _mp2_energy(eri, epsilon, intdump['uhf'])

    return energies


def _scf_energy(Hcore, ERI, occ_sl, unrestricted):
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


def _mp2_energy(ERI, epsilon, unrestricted):
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
