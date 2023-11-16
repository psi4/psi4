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

__all__ = [
    "get_qm_atoms_opts",
    "modify_Fock_induced",
    "modify_Fock_permanent",
]

from typing import Any, Dict, List, Tuple

import numpy as np

from psi4 import core


def get_qm_atoms_opts(mol: core.Molecule) -> Tuple[List[float], List[float], Dict[str, Any]]:
    """Provides list of coordinates of quantum mechanical atoms from
    psi4.core.Molecule `mol` to pylibefp.core.efp() `efpobj`. Also
    converts from `read_options("EFP"` to pylibefp opts dictionary.

    """
    efpobj = mol.EFP

    ptc = []
    coords = []
    for iat in range(mol.natom()):
        ptc.append(mol.charge(iat))
        coords.append(mol.x(iat))
        coords.append(mol.y(iat))
        coords.append(mol.z(iat))

    # set options
    # * 'chtr', 'qm_exch', 'qm_disp', 'qm_chtr' may be enabled in a future libefp release
    opts = {}
    for opt in ['elst', 'exch', 'ind', 'disp',
                'elst_damping', 'ind_damping', 'disp_damping']:
        psiopt = 'EFP_' + opt.upper()
        if core.has_option_changed('EFP', psiopt):
            opts[opt] = core.get_option('EFP', psiopt)
    for opt in ['elst', 'ind']:
        psiopt = 'EFP_QM_' + opt.upper()
        if core.has_option_changed('EFP', psiopt):
            opts['qm_' + opt] = core.get_option('EFP', psiopt)

    return ptc, coords, opts


def modify_Fock_permanent(mol: core.Molecule, mints: core.MintsHelper, verbose: int = 1) -> np.ndarray:
    """Computes array of the EFP contribution to the potential felt by
    QM atoms due to permanent EFP moments. Used for SCF procedure.

    Parameters
    ----------
    mol
        Source of quantum mechanical atoms. As its `EFP` member data, contains
        a :py:class:`pylibefp.core.efp` object that is the source and computer
        of EFP fragments.
    mints
        Integral computer.
    verbose
        Whether to print out multipole coordinates and values. 0: no printing.
        1: print charges and dipoles. 2: additionally print quadrupoles and octupoles.

    Returns
    -------
    ~numpy.ndarray
        (nbf, nbf) EFP charge through octupole contribution to the potential

    """
    # get composition counts from libefp
    efpobj = mol.EFP
    nfr = efpobj.get_frag_count()
    natoms = efpobj.get_frag_atom_count()

    # get multipoles count, pos'n, values from libefp
    #   charge + dipoles + quadrupoles + octupoles = 20
    nmp = efpobj.get_multipole_count()
    xyz_mp = np.asarray(efpobj.get_multipole_coordinates(verbose=verbose)).reshape(nmp, 3)
    val_mp = np.asarray(efpobj.get_multipole_values(verbose=verbose)).reshape(nmp, 20)

    #                    0  X  Y  Z  XX   YY   ZZ   XY   XZ   YZ
    prefacs = np.array([ 1, 1, 1, 1, 1/3, 1/3, 1/3, 2/3, 2/3, 2/3,
        1/15, 1/15, 1/15, 3/15, 3/15, 3/15, 3/15, 3/15, 3/15, 6/15])
    #   XXX   YYY   ZZZ   XXY   XXZ   XYY   YYZ   XZZ   YZZ   XYZ

    # EFP permanent moment contribution to the Fock Matrix
    nbf = mints.basisset().nbf()
    V2 = np.zeros((nbf, nbf))

    # Cartesian basis one-electron EFP perturbation
    efp_ints = np.zeros((20, nbf, nbf))

    for imp in range(nmp):
        origin = xyz_mp[imp]

        # get EFP multipole integrals from Psi4
        p4_efp_ints = mints.ao_efp_multipole_potential(origin=origin)
        for pole in range(20):
            efp_ints[pole] = np.asarray(p4_efp_ints[pole])

        # add frag atom Z into multipole charge (when pos'n of atom matches mp)
        for ifr in range(nfr):
            atoms = efpobj.get_frag_atoms(ifr)
            for iat in range(natoms[ifr]):
                xyz_atom = [atoms[iat]['x'], atoms[iat]['y'], atoms[iat]['z']]
                if np.allclose(xyz_atom, origin, atol=1e-10):
                    val_mp[imp, 0] += atoms[iat]['Z']

        # scale multipole integrals by multipole magnitudes. result goes into V
        for pole in range(20):
            efp_ints[pole] *= -prefacs[pole] * val_mp[imp, pole]
            V2 += efp_ints[pole]

    return V2


def modify_Fock_induced(efpobj: "pylibefp.core.efp", mints: core.MintsHelper, verbose: int = 1) -> np.ndarray:
    """Returns shared matrix containing the EFP contribution to the potential
    felt by QM atoms due to EFP induced dipoles. Used in SCF procedure.

    Parameters
    ----------
    efpobj
        Source of EFP induced dipole information.
    mints
        Integral computer.
    verbose
        Whether to print out induced dipole coordinates and values.
        0: no printing. 1: print induced dipole info.

    Returns
    -------
    ndarray
        (nbf, nbf) EFP contribution to potential.

    """
    # get induced dipoles count, pos'n, values from libefp
    #   dipoles = 3
    nid = efpobj.get_induced_dipole_count()
    xyz_id = np.asarray(efpobj.get_induced_dipole_coordinates(verbose=verbose)).reshape(nid, 3)
    val_id = np.asarray(efpobj.get_induced_dipole_values(verbose=verbose)).reshape(nid, 3)
    val_idt = np.asarray(efpobj.get_induced_dipole_conj_values(verbose=verbose)).reshape(nid, 3)

    # take average of induced dipole and conjugate
    val_id = (val_id + val_idt) * 0.5

    # EFP induced dipole contribution to the Fock Matrix
    coords = core.Matrix.from_array(xyz_id)
    V_ind = mints.induction_operator(coords, core.Matrix.from_array(val_id)).np
    return V_ind
