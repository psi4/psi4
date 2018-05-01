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

import re
import sys
import itertools

import numpy as np

from ..exceptions import *


try:
    long(1)
except NameError:
    long = int

try:
    from itertools import ifilterfalse as filterfalse  # py2
except ImportError:
    from itertools import filterfalse  # py3


def _unique_everseen(iterable, key=None):
    "List unique elements, preserving order. Remember all elements ever seen."
    # unique_everseen('AAAABBBCCDAABBB') --> A B C D
    # unique_everseen('ABBCcAD', str.lower) --> A B C D
    # straight from the docs, https://docs.python.org/3/library/itertools.html#itertools-recipes
    seen = set()
    seen_add = seen.add
    if key is None:
        for element in filterfalse(seen.__contains__, iterable):
            seen_add(element)
            yield element
    else:
        for element in iterable:
            k = key(element)
            if k not in seen:
                seen_add(k)
                yield element


def _apply_default(llist, default):
    return [default if (c is None) else c for c in llist]


def _high_spin_sum(mult_list):
    mm = 1
    for m in mult_list:
        mm += m - 1
    return mm


def _mult_ok(m):
    return isinstance(m, (int, np.int64, long)) and m >= 1


def _sufficient_electrons_for_mult(z, c, m):
    """Require sufficient electrons in total: total mult ({}) - 1 > raw electrons ({}) - total chg ({})"""
    return m - 1 <= z - c


def _parity_ok(z, c, m):
    """Check total electrons (neutral protons `z` and charge `c`) is (un)paired-compatible with multiplicity `m`"""
    return (m % 2) != ((z - c) % 2)


def _alpha_beta_allocator(z, c, m):
    nbeta = (z - c - m + 1) // 2
    nalpha = nbeta + m - 1
    return nalpha, nbeta


def validate_and_fill_chgmult(zeff,
                              fragment_separators,
                              molecular_charge,
                              fragment_charges,
                              molecular_multiplicity,
                              fragment_multiplicities,
                              zero_ghost_fragments=False,
                              verbose=1):
    """
    Applies physical constraints and sensible defaults to reconciling and
    completing the molecular and fragment charge and multiplicity
    specification.

    Parameters
    ----------
    zeff : ndarray of float
        (nat,) electron counts for neutral atoms, generally Z nuclear charge.
        0 indicates ghosts such that a full fragment of 0s will be constained
        to `0 1` charge & multiplicity.
    fragment_separators : ndarray of int
        (nfr - 1, ) indices splitting `zeff` into nfr fragments.
    molecular_charge : float or None
        Total charge for molecular system.
    fragment_charges : list of float or None
        (nfr,) known fragment charges with `None` as placeholder for
        unknown. Expected pre-defaulted so even if nothing known if
        `fragment_separators` breaks `zeff` into `nfr=2` fragments, input
        value should be `fragment_charges=[None, None]`.
    molecular_multiplicity : int or None
        Total multiplicity for molecular system.
    fragment_multiplicity : list of int or None
        (nfr,) known fragment charges with `None` as placeholder for
        unknown. Expected pre-defaulted so even if nothing known if
        `fragment_separators` breaks `zeff` into `nfr=2` fragments, input
        value should be `fragment_multiplicities=[None, None]`.
    zero_ghost_fragments : bool, optional
        Fragments composed entirely of ghost atoms (Zeff=0) are required to have
        chgmult `0 1`. When `False`, violations of this will cause a
        ValidationError. When `True`, treat ghost fragments indicated by `zeff` to
        contain superior information over chgmult arguments that might still
        correspond to full-real molecule. Clears information from
        `molecular_charge` and `molecular_multiplicity` and sets ghost fragments
        to `0 1`, leaving other positions free to readjust. Unused (prefer to set
        up such manipulations outside function call) but works.
    verbose : int, optional
        Amount of printing.

    Notes
    -----
    Returns combination of total & fragment charge & multiplicity among
    values of S1-7 that fulfill rules R1-9. A few derived implications in I1-3.

    * Constraints
    R1 * require all chg & mult exist
    R2 * require total charge to be the sum of frag chg
    R3 * require mult is positive int
    R4 * require sufficient tot electrons for mult: mult - 1 <= neut_el - chg
    R5 * require total parity consistent among tot electrons and mult: (mult % 2) != ((neut_el - chg) % 2)
    R6 * require chg match input argument values
    R7 * require mult match input argument values
    R8 * require that tot = sum(frag) mult follow high spin addition unless tot & frag mult fully specified
    R9 * require that ghost fragments (zeff all 0) be neutral singlet

    * Allowed values
    S1 * suggest input argument values for tot chg, frag chg, tot mult or frag mult
    S2 * suggest sum frag chg for tot chg, allowing for indiv frag chg defaulting to 0
    S3 * suggest distributing unallocated chg onto frag chg
    S4 * suggest 0 default for frag chg
    S5 * suggest range of high-spin sum frag mult for tot mult, allowing for indiv frag mult defaulting to 1 or 2
    S6 * suggest range of unallocated mult = tot - high_spin_sum(frag - 1), allowing for all indiv but self defaulting to 1 or 2.
    S7 * suggest 1 or 2 default for frag mult

    * Implications
    I1 * won't form an ion just to be closed shell (would require choosing +1 vs. -1)
    I2 * unallocated chg or mult lands on the first unspecified fragment able to
         bear it (enforced by returning first match encountered; subsequent
         matches distribute charge to later frags)
    I3 * missing chg or mult from tot - frags will always be allocated as a block, not distributed

    Examples
    --------
    >>> validate_and_fill_chgmult(*sys('He'), 0, [0], 1, [1])
    0, [0], 1, [1]
    >>> validate_and_fill_chgmult(*sys('He'), None, [None], None, [None])
    0, [0], 1, [1]
    >>> validate_and_fill_chgmult(*sys('He/He'), None, [None, None], None, [None, None])
    0, [0, 0], 1, [1, 1])
    >>> validate_and_fill_chgmult(*sys('He/He'), 2, [None, None], None, [None, None])
    2, [2, 0], 1, [1, 1])
    >>> validate_and_fill_chgmult(*sys('He/He'), None, [2, None], None, [None, None])
    2, [2, 0], 1, [1, 1])
    >>> validate_and_fill_chgmult(*sys('He/He'), 0, [2, None], None, [None, None])
    0, [2, -2], 1, [1, 1])
    >>> validate_and_fill_chgmult(*sys('Ne/He/He'), -2, [None, 2, None], None, [None, None, None])
    -2, [-4, 2, 0], 1, [1, 1, 1]
    >>> validate_and_fill_chgmult(*sys('Ne/He/He'), 2, [None, -2, None], None, [None, None, None])
    2, [4, -2, 0], 1, [1, 1, 1]
    #  9 - residual +4 distributes to first fragment able to wholly accept it (He+4 is no-go)
    >>> validate_and_fill_chgmult(*sys('He/He/Ne'), 2, [None, -2, None], None, [None, None, None])
    2, [0, -2, 4], 1, [1, 1, 1]
    # 10 - residual +4 unsuited for only open fragment, He, so irreconcilable
    >>> validate_and_fill_chgmult(*sys('He/He/Ne'), 2, [None, -2, 0], None, [None, None, None])
    ValidationError
    # 11 - non-positive multiplicity
    >>> validate_and_fill_chgmult(*sys('He/He/Ne'), 2, [2, -2, None], None, [None, None, None])
    2, [2, -2, 2], 1, [1, 1, 1])
    >>> validate_and_fill_chgmult(*sys('He/He'), None, [-2, 2], None, [None, None])
    0, [-2, 2], 1, [1, 1]
    >>> validate_and_fill_chgmult(*sys('He/He'), None, [None, -2], None, [None, None])
    -2, [0, -2], 1, [1, 1]
    >>> validate_and_fill_chgmult(*sys('Ne/Ne'), 0, [None, 4], None, [None, None])
    0, [-4, 4], 1, [1, 1]
    >>> validate_and_fill_chgmult(*sys('He/He/He'), 4, [2, None, None], None, [None, None, None])
    4, [2, 2, 0], 1, [1, 1, 1]
    >>> validate_and_fill_chgmult(*sys('He/He'), 0, [-2, 2], None, [None, None])
    0, [-2, 2], 1, [1, 1]
    >>> validate_and_fill_chgmult(*sys('He/He'), 0, [-2, -2], None, [None, None])
    ValidationError
    >>> validate_and_fill_chgmult(*sys('He'), None, [None], 0, [None])
    ValidationError
    >>> validate_and_fill_chgmult(*sys('He'), None, [None], None, [1])
    0, [0], 1, [1]
    # 20 - doublet non consistent with closed-shell, neutral default charge
    >>> validate_and_fill_chgmult(*sys('He'), None, [None], None, [2])
    ValidationError
    >>> validate_and_fill_chgmult(*sys('He'), None, [None], None, [3])
    0, [0], 3, [3]
    # 22 - insufficient electrons for pentuplet
    >>> validate_and_fill_chgmult(*sys('He'), None, [None], None, [5])
    ValidationError
    >>> validate_and_fill_chgmult(*sys('He'), None, [-1], None, [2])
    -1, [-1], 2, [2]
    # 24 - doublet not consistent with even charge
    >>> validate_and_fill_chgmult(*sys('He'), None, [-2], None, [2])
    ValidationError
    >>> validate_and_fill_chgmult(*sys('He/He'), None, [None, None], None, [1, 1])
    0, [0, 0], 1, [1, 1]
    >>> validate_and_fill_chgmult(*sys('He/He'), None, [None, None], None, [3, 1])
    0, [0, 0], 3, [3, 1]
    >>> validate_and_fill_chgmult(*sys('He/He'), None, [None, None], None, [1, 3])
    0, [0, 0], 3, [1, 3]
    >>> validate_and_fill_chgmult(*sys('He/He'), None, [None, None], None, [3, 3])
    0, [0, 0], 5, [3, 3]
    >>> validate_and_fill_chgmult(*sys('He/He'), None, [None, None], 3, [3, 3])
    0, [0, 0], 3, [3, 3]
    # 30 - bad parity btwn mult and total # electrons
    >>> validate_and_fill_chgmult(*sys('He/He'), None, [None, None], 2, [3, 3])
    ValidationError
    >>> validate_and_fill_chgmult(*sys('H'), None, [None], None, [None])
    0, [0], 2, [2]
    >>> validate_and_fill_chgmult(*sys('H'), 1, [None], None, [None])
    1, [1], 1, [1]
    >>> validate_and_fill_chgmult(*sys('H'), None, [-1], None, [None])
    -1, [-1], 1, [1]
    >>> validate_and_fill_chgmult(*sys('funnyH'), None, [None], None, [None])
    0, [0], 1, [1]
    # 35 - insufficient electrons
    >>> validate_and_fill_chgmult(*sys('funnierH'), None, [None], None, [None])
    ValidationError
    >>> validate_and_fill_chgmult(*sys('H/H'), None, [None, None], None, [None, None])
    0, [0, 0], 3, [2, 2]
    >>> validate_and_fill_chgmult(*sys('H/He'), None, [None, None], None, [None, None])
    0, [0, 0], 2, [2, 1]
    >>> validate_and_fill_chgmult(*sys('H/He'), None, [1, 1], None, [None, None])
    2, [1, 1], 2, [1, 2]
    >>> validate_and_fill_chgmult(*sys('H/He'), -2, [-1, None], None, [None, None])
    -2, [-1, -1], 2, [1, 2]
    >>> validate_and_fill_chgmult(*sys('H/He/Na/Ne'), None, [1, None, 1, None], None, [None, None, None, None])
    2, [1, 0, 1, 0], 1, [1, 1, 1, 1]
    >>> validate_and_fill_chgmult(*sys('H/He/Na/Ne'), None, [-1, None, 1, None], None, [None, None, None, None])
    0, [-1, 0, 1, 0], 1, [1, 1, 1, 1]
    >>> validate_and_fill_chgmult(*sys('H/He/Na/Ne'), 2, [None, None, 1, None], None, [None, None, None, None])
    2, [1, 0, 1, 0], 1, [1, 1, 1, 1]
    >>> validate_and_fill_chgmult(*sys('H/He/Na/Ne'), 3, [None, None, 1, None], None, [None, None, None, None])
    3, [0, 2, 1, 0], 2, [2, 1, 1, 1]
    >>> validate_and_fill_chgmult(*sys('H/He'), None, [1, None], None, [2, None])
    ValidationError
    >>> validate_and_fill_chgmult(*sys('H/He'), None, [None, 0], None, [None, 2])
    ValidationError
    >>> validate_and_fill_chgmult(*sys('H/He'), None, [None, -1], None, [None, 3])
    ValidationError
    >>> validate_and_fill_chgmult(*sys('H/He/Na/Ne'), None, [None, 1, 0, 1], None, [None, None, None, None])
    2, [0, 1, 0, 1], 5, [2, 2, 2, 2]
    >>> validate_and_fill_chgmult(*sys('H/He/Na/Ne'), None, [None, 1, 0, None], None, [None, None, None, None])
    1, [0, 1, 0, 0], 4, [2, 2, 2, 1]
    >>> validate_and_fill_chgmult(*sys('H/He/Na/Ne'), None, [None, 1, 0, None], None, [None, None, 4, None])
    1, [0, 1, 0, 0], 6, [2, 2, 4, 1]
    >>> validate_and_fill_chgmult(*sys('He/He/He'), 0, [None, None, 1], None, [1, None, 2])
    0, [0, -1, 1], 3, [1, 2, 2]
    >>> validate_and_fill_chgmult(*sys('N/N/N'), None, [1, 1, 1], 3, [None, 3, None])
    3, [1, 1, 1], 3, [1, 3, 1]
    >>> validate_and_fill_chgmult(*sys('N/N/N'), None, [1, 1, 1], 3, [None, None, None])
    3, [1, 1, 1], 3, [3, 1, 1]
    >>> validate_and_fill_chgmult(*sys('N/N/N'), None, [None, None, None], 3, [None, None, 2])
    ValidationError
    >>> validate_and_fill_chgmult(*sys('N/N/N'), 1, [None, -1, None], 3, [None, None, 2])
    1, [2, -1, 0], 3, [2, 1, 2]
    # 55 - both (1, (1, 0.0, 0.0), 4, (1, 3, 2)) and (1, (0.0, 0.0, 1), 4, (2, 3, 1)) plausible
    >>> validate_and_fill_chgmult(*sys('N/Ne/N'), 1, [None, None, None], 4, [None, 3, None])
    1, [1, 0, 0], 4, [1, 3, 2]
    >>> validate_and_fill_chgmult(*sys('N/Ne/N'), None, [None, None, 1], 4, [None, 3, None])
    1, [0, 0, 1], 4, [2, 3, 1]
    >>> validate_and_fill_chgmult(*sys('He/He'), None, [-1, 1], None, [None, None])
    0, [-1, 1], 3, [2, 2]
    >>> validate_and_fill_chgmult(*sys('Gh'), 1, [None], None, [None])
    ValidationError
    >>> validate_and_fill_chgmult(*sys('Gh'), -1, [None], None, [None])
    ValidationError
    >>> validate_and_fill_chgmult(*sys('Gh'), None, [None], 3, [None])
    ValidationError
    >>> validate_and_fill_chgmult(*sys('He/Gh'), None, [2, None], None, [None, None])
    2, [2, 0], 1, [1, 1]
    >>> validate_and_fill_chgmult(*sys('Gh/He'), None, [2, None], None, [None, None])
    ValidationError
    >>> validate_and_fill_chgmult(*sys('Gh/He/Ne'), 2, [None, -2, None], None, [None, None, None])
    2, [0, -2, 4], 1, [1, 1, 1]
    >>> validate_and_fill_chgmult(*sys('Gh/He/Gh'), 1, [None, None, None], None, [None, None, None])
    1, [0, 1, 0], 2, [1, 2, 1]
    >>> sys = {
        'He': (np.array([2]), np.array([])),
        'He/He': (np.array([2, 2]), np.array([1])),
        'Ne/He/He': (np.array([10, 2, 2]), np.array([1, 2])),
        'He/He/Ne': (np.array([2, 2, 10]), np.array([1, 2])),
        'Ne/Ne': (np.array([10, 10]), np.array([1])),
        'He/He/He': (np.array([2, 2, 2]), np.array([1, 2])),
        'H': (np.array([1]), np.array([])),
        'funnyH': (np.array([0]), np.array([])),  # has no electrons
        'funnierH': (np.array([-1]), np.array([])),  # has positron
        'H/H': (np.array([1, 1]), np.array([1])),
        'H/He': (np.array([1, 2]), np.array([1])),
        'H/He/Na/Ne': (np.array([1, 2, 11, 10]), np.array([1, 2, 3])),
        'N/N/N': (np.array([7, 7, 7]), np.array([1, 2])),
        'N/Ne/N': (np.array([7, 10, 7]), np.array([1, 2])),
        'He/Gh': (np.array([2, 0]), np.array([1])),
        'Gh/He': (np.array([0, 2]), np.array([1])),
        'Gh': (np.array([0, 0]), np.array([])),
        'Gh/He/Ne': (np.array([0, 0, 2, 10]), np.array([2, 3])),
        'Gh/He/Gh': (np.array([0, 2, 0]), np.array([1, 2]))}

    """
    text = []

    felez = np.split(zeff, fragment_separators)
    nfr = len(felez)
    text.append('felez: {}'.format(felez))

    cgmp_exact_c = []  # exact_* are candidates for the final value
    cgmp_exact_fc = [[] for f in range(nfr)]
    cgmp_exact_m = []
    cgmp_exact_fm = [[] for f in range(nfr)]

    cgmp_range = []  # tests that the final value must pass to be valid
    cgmp_rules = []  # key to what rules in cgmp_range are T/F

    real_fragments = np.array([not all(f == 0 for f in felez[ifr]) for ifr in range(nfr)])
    all_fc_known = all(f is not None for f in fragment_charges)
    all_fm_known = all(f is not None for f in fragment_multiplicities)

    if zero_ghost_fragments and not all(real_fragments):
        print('possibly adjusting charges')
        molecular_charge = None
        fragment_charges = [(fr if real_fragments[ifr] else 0.0) for ifr, fr in enumerate(fragment_charges)]
        molecular_multiplicity = None
        fragment_multiplicities = [(fr if real_fragments[ifr] else 1) for ifr, fr in enumerate(fragment_multiplicities)]

    # <<< assert broad physical requirements

    #   * (R1) require all chg & mult exist
    cgmp_range.append(lambda c, fc, m, fm: c is not None and
                                           all(f is not None for f in fc) and
                                           m is not None and
                                           all(f is not None for f in fm))
    cgmp_rules.append('1')

    #   * (R2) require total charge to be the sum of fragment charges
    cgmp_range.append(lambda c, fc, m, fm: c == sum(fc))
    cgmp_rules.append('2')

    #   * (R3) require mult is positive int
    cgmp_range.append(lambda c, fc, m, fm: _mult_ok(m) and all(_mult_ok(f) for f in fm))
    cgmp_rules.append('3')

    # <<< assert electron count requirements

    zel = np.sum(zeff)  # note: number electrons in neutral species, not number total electrons
    fzel = [np.sum(f) for f in felez]
    text.append('zel: {}'.format(zel))
    text.append('fzel: {}'.format(fzel))

    #   * (R4) require sufficient electrons for mult: mult - 1 <= neutral_electrons - chg
    cgmp_range.append(lambda c, fc, m, fm: _sufficient_electrons_for_mult(zel, c, m))
    cgmp_rules.append('4')
    for ifr in range(nfr):
        cgmp_range.append(lambda c, fc, m, fm, ifr=ifr: _sufficient_electrons_for_mult(fzel[ifr], fc[ifr], fm[ifr]))
        cgmp_rules.append('4-' + str(ifr))

    #   * (R5) require total parity consistent among neutral_electrons, chg, and mult
    cgmp_range.append(lambda c, fc, m, fm: _parity_ok(zel, c, m))
    cgmp_rules.append('5')
    for ifr in range(nfr):
        cgmp_range.append(lambda c, fc, m, fm, ifr=ifr: _parity_ok(fzel[ifr], fc[ifr], fm[ifr]))
        cgmp_rules.append('5-' + str(ifr))

    # <<< (R6, R7, S1) assert & suggest input values

    if molecular_charge is not None:
        cgmp_exact_c.append(molecular_charge)
        cgmp_range.append(lambda c, fc, m, fm: c == molecular_charge)
        cgmp_rules.append('6')
    for ifr, chg in enumerate(fragment_charges):
        if chg is not None:
            cgmp_exact_fc[ifr].append(chg)
            cgmp_range.append(lambda c, fc, m, fm, ifr=ifr, chg=chg: fc[ifr] == chg)
            cgmp_rules.append('6-' + str(ifr))
    if molecular_multiplicity is not None:
        cgmp_exact_m.append(molecular_multiplicity)
        cgmp_range.append(lambda c, fc, m, fm: m == molecular_multiplicity)
        cgmp_rules.append('7')
    for ifr, mult in enumerate(fragment_multiplicities):
        if mult is not None:
            cgmp_exact_fm[ifr].append(mult)
            cgmp_range.append(lambda c, fc, m, fm, ifr=ifr, mult=mult: fm[ifr] == mult)
            cgmp_rules.append('7-' + str(ifr))

    # <<< assert high-spin-rule and suggest "missing quantity" and default values

    #   * (S2) suggest net frag charge for total charge, allowing for indiv frag defaulting to 0
    cgmp_exact_c.append(sum(filter(None, fragment_charges)))

    missing_frag_chg = 0. if molecular_charge is None else molecular_charge
    missing_frag_chg -= sum(filter(None, fragment_charges))

    #   * (S3) suggest distributing unallocated charge onto fragment
    #   * (S4) suggest 0 default charge for fragment
    for ifr in range(nfr):
        if fragment_charges[ifr] is None:  # unneeded, but shortens the exact lists
            cgmp_exact_fc[ifr].append(missing_frag_chg)
            cgmp_exact_fc[ifr].append(0.)

    #   * (R8) require that frag mult follow high spin addition unless fully specified
    if molecular_multiplicity is None or any(f is None for f in fragment_multiplicities):
        cgmp_range.append(lambda c, fc, m, fm: m == _high_spin_sum(fm))
        cgmp_rules.append('8')

    #   * (S5) suggest range of net frag mult for total mult, allowing for indiv frag defaulting to 1 or 2.
    #          many in range may be unphysical, but those will be caught by physical rules.
    if molecular_multiplicity is None:  # unneeded, but shortens the exact lists
        frag_mult_hi = _high_spin_sum(_apply_default(fragment_multiplicities, 2))
        frag_mult_lo = _high_spin_sum(_apply_default(fragment_multiplicities, 1))
        for m in range(frag_mult_lo, frag_mult_hi + 1):
            cgmp_exact_m.append(m)

    #   * (S6) suggest range of missing mult = tot - high_spin_sum(frag - 1),
    #          allowing for all indiv but self defaulting to 1 or 2.  Many in range
    #          may be unphysical, but those will be caught by physical rules.
    #   * (S7) suggest 1 or 2 default multiplicity for fragment
    if molecular_multiplicity is not None and any(f is None for f in fragment_multiplicities):
        frag_mult_less_one_none = fragment_multiplicities[:]
        frag_mult_less_one_none.remove(None)  # "missing" slot to solve for
        frag_mult_hi = _high_spin_sum(_apply_default(frag_mult_less_one_none, 2))
        frag_mult_lo = _high_spin_sum(_apply_default(frag_mult_less_one_none, 1))
        missing_mult_hi = molecular_multiplicity - frag_mult_lo + 1
        missing_mult_lo = molecular_multiplicity - frag_mult_hi + 1
    else:
        missing_mult_hi = 0
        missing_mult_lo = 0

    for ifr in range(nfr):
        if fragment_multiplicities[ifr] is None:  # unneeded, but shortens the exact lists
            for m in reversed(range(max(missing_mult_lo, 1), missing_mult_hi + 1)):
                cgmp_exact_fm[ifr].append(m)
            cgmp_exact_fm[ifr].append(1)
            cgmp_exact_fm[ifr].append(2)

    #   * (R9) require that ghost fragments be neutral singlet
    for ifr in range(nfr):
        if all(f == 0 for f in felez[ifr]):
            cgmp_range.append(lambda c, fc, m, fm, ifr=ifr: fc[ifr] == 0 and fm[ifr] == 1)
            cgmp_rules.append('9-' + str(ifr))

    # <<< reconcile and report

    def reconcile(exact_c, exact_fc, exact_m, exact_fm):
        """Returns a member from all combinations of `exact` that passes all tests in cgmp_range, else raises error."""

        # remove duplicates
        uniq_c = _unique_everseen(exact_c)
        uniq_fc = [_unique_everseen(f) for f in exact_fc]
        uniq_m = _unique_everseen(exact_m)
        uniq_fm = [_unique_everseen(f) for f in exact_fm]
        text.append('c: {}'.format(list(exact_c)))
        for f in exact_fc:
            text.append('fc:'.format(list(f)))
        text.append('m: {}'.format(list(exact_m)))
        for f in exact_fm:
            text.append('fm:'.format(list(f)))

        header = True
        for candidate in itertools.product(*[uniq_c, itertools.product(*uniq_fc),
                                             uniq_m, itertools.product(*uniq_fm)]):  # yapf: disable
            cc, cfc, cm, cfm = candidate
            if header:
                text.append(
                    """Assess candidate {}: {}""".format(candidate, ' '.join(('{:3}'.format(r) for r in cgmp_rules))))
                header = False
            assessment = [fn(cc, cfc, cm, cfm) for fn in cgmp_range]
            sass = ['{:3}'.format('T' if b else '') for b in assessment]
            text.append("""Assess candidate {:}: {} --> {}""".format(candidate, ' '.join(sass), all(assessment)))
            if all(assessment):
                return candidate
        else:
            err = """Inconsistent or unspecified chg/mult: sys chg: {}, frag chg: {}, sys mult: {}, frag mult: {}""".format(
                molecular_charge, fragment_charges, molecular_multiplicity, fragment_multiplicities)
            if verbose > -1:
                print('\n\n' + '\n'.join(text))
            raise ValidationError(err)

    def stringify(start, final):
        fcgmp = '{:^4}'
        return fcgmp.format(final) if final == start else fcgmp.format('(' + str(int(final)) + ')')

    # TODO could winnow down the exact_* lists a bit by ruling out
    #      independent values. do this if many-frag molecular systems take too
    #      long in the itertools.product

    c_final, fc_final, m_final, fm_final = reconcile(cgmp_exact_c, cgmp_exact_fc, cgmp_exact_m, cgmp_exact_fm)

    c_text = stringify(molecular_charge, c_final)
    fc_text = ', '.join((stringify(fs, ff) for fs, ff in zip(fragment_charges, fc_final)))
    m_text = stringify(molecular_multiplicity, m_final)
    fm_text = ', '.join((stringify(fs, ff) for fs, ff in zip(fragment_multiplicities, fm_final)))

    brief = []
    brief.append('    {:26} {}'.format('      charge = ' + c_text, 'fragments = ' + fc_text))
    brief.append('    {:26} {}'.format('multiplicity = ' + m_text, 'fragments = ' + fm_text))

    been_defaulted = []
    if c_text.count('(') + fc_text.count('(') > 1:
        been_defaulted.append('charge')
    if '(' in m_text or '(' in fm_text:
        been_defaulted.append('multiplicity')

    if been_defaulted:
        brief.append('    Note: Default values have been applied for {}. Specify intentions in molecule input block'.
                     format(' and '.join(been_defaulted)))

    if m_final != _high_spin_sum(fm_final):
        brief.append(
            '    Warning: Total multiplicity is not high-spin sum of fragments; may be clobbered by psi4.core.Molecule.update_geometry().'
        )

    if verbose >= 2:
        print('\n'.join(text))
    if verbose >= 1:
        # TODO add back when printing worked out
        #print('\n'.join(brief))
        pass

    return {'molecular_charge': float(c_final),
            'fragment_charges': [float(f) for f in fc_final],
            'molecular_multiplicity': m_final,
            'fragment_multiplicities': list(fm_final)}
