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

import warnings

import qcelemental as qcel


def B787(cgeom,
         rgeom,
         cuniq,
         runiq,
         do_plot=False,
         verbose=1,
         atoms_map=False,
         run_resorting=False,
         mols_align=False,
         run_to_completion=False,
         uno_cutoff=1.e-3,
         run_mirror=False):
    """Use Kabsch algorithm to find best alignment of geometry `cgeom` onto
    `rgeom` while sampling atom mappings restricted by `runiq` and `cuniq`.

    Parameters
    ----------
    rgeom : ndarray of float
        (nat, 3) array of reference/target/unchanged geometry. Assumed [a0]
        for RMSD purposes.
    cgeom : ndarray of float
        (nat, 3) array of concern/changeable geometry. Assumed [a0] for RMSD
        purposes. Must have same nat, units, and atom content as rgeom.
    runiq : ndarray of str
        (nat,) array indicating which rows (atoms) in `rgeom` are shuffleable
        without changing the molecule. Generally hashes of element symbol and
        mass are used, but could be as simple as ['C', 'H', 'H', 'D', 'H'] for
        monodeuterated methane.
    cuniq : ndarray of str
        (nat,) array indicating which rows (atoms) in `cgeom` are shuffleable.
        See `runiq` for more details. Strings and count in `cuniq` must match
        `runiq`. That is, `sorted(cuniq) == sorted(runiq)`.
    do_plot : bool, optional
        Pops up a mpl plot showing before, after, and ref geometries.
    verbose : int, optional
        Quantity of printing. 0 to silence.
    atoms_map : bool, optional
        Whether atom1 of rgeom already corresponds to atom1 of cgeom and so on.
        If `True`, no resorting will be run, parameters `runiq` and `cuniq`
        may be passed as `None`, and much time will be saved.
    run_resorting : bool, optional
        Run the resorting machinery even if unnecessary because `atoms_map=True`.
    mols_align : bool or float, optional
        Whether ref_mol and concern_mol have identical geometries by eye
        (barring orientation or atom mapping) and expected final RMSD = 0.
        If `True`, procedure is truncated when RMSD condition met, saving time.
        If float, convcrit at which search for minimium truncates.
    run_to_completion : bool, optional
        Run reorderings to completion (past RMSD = 0) even if unnecessary because
        `mols_align=True`. Used to test worst-case timings.
    uno_cutoff : float, optional
        TODO
    run_mirror : bool, optional
        Run alternate geometries potentially allowing best match to `rgeom`
        from mirror image of `cgeom`. Only run if system confirmed to
        be nonsuperimposable upon mirror reflection.

    Returns
    -------
    float, tuple
        First item is RMSD [A] between `rgeom` and the optimally aligned
        geometry computed.
        Second item is a AlignmentMill namedtuple with fields
        (shift, rotation, atommap, mirror) that prescribe the transformation
        from `cgeom` and the optimally aligned geometry.

    """
    warnings.warn(
        "Using `qcdb.align.B787` instead of `qcelemental.molutil.B787` is deprecated, and as soon as 1.5 it will stop working\n",
        category=FutureWarning,
        stacklevel=2)

    return qcel.molutil.B787(cgeom,
                             rgeom,
                             cuniq,
                             runiq,
                             do_plot=do_plot,
                             verbose=verbose,
                             atoms_map=atoms_map,
                             run_resorting=run_resorting,
                             mols_align=mols_align,
                             run_to_completion=run_to_completion,
                             uno_cutoff=uno_cutoff,
                             run_mirror=run_mirror)


def compute_scramble(nat, do_resort=True, do_shift=True, do_rotate=True, deflection=1.0, do_mirror=False):
    """Generate a random or directed translation, rotation, and atom shuffling.

    Parameters
    ----------
    nat : int
        Number of atoms for which to prepare an atom mapping.
    do_resort : bool or array-like, optional
        Whether to randomly shuffle atoms (`True`) or leave 1st atom 1st, etc. (`False`)
        or shuffle according to specified (nat, ) indices (e.g., [2, 1, 0])
    do_shift : bool or array-like, optional
        Whether to generate a random atom shift on interval [-3, 3) in each
        dimension (`True`) or leave at current origin (`False`) or shift along
        specified (3, ) vector (e.g., np.array([0., 1., -1.])).
    do_rotate : bool or array-like, optional
        Whether to generate a random 3D rotation according to algorithm of Arvo (`True`)
        or leave at current orientation (`False`) or rotate with specified (3, 3) matrix.
    deflection : float, optional
        If `do_rotate`, how random a rotation: 0.0 is no change, 0.1 is small
        perturbation, 1.0 is completely random.
    do_mirror : bool, optional
        Whether to set mirror reflection instruction. Changes identity of
        molecule so off by default.

    Returns
    -------
    tuple
        AlignmentMill namedtuple with fields (shift, rotation, atommap, mirror)
        as requested: identity, random, or specified.

    """
    warnings.warn(
        "Using `qcdb.align.compute_scramble` instead of `qcelemental.molutil.compute_scramble` is deprecated, and as soon as 1.5 it will stop working\n",
        category=FutureWarning,
        stacklevel=2)

    return qcel.molutil.compute_scramble(nat, do_resort=do_resort, do_shift=do_shift, do_rotate=do_rotate, deflection=deflection, do_mirror=do_mirror)
