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

import os
import sys
import hashlib
import warnings
import itertools
import collections

import qcelemental as qcel

from .exceptions import *
from .psiutil import search_file
from .molecule import Molecule
from .libmintsgshell import ShellInfo
from .libmintsbasissetparser import Gaussian94BasisSetParser
from .basislist import corresponding_basis, corresponding_zeta


basishorde = {}

class BasisSet(object):
    """Basis set container class
    Reads the basis set from a checkpoint file object. Also reads the molecule
    from the checkpoint file storing the information in an internal Molecule class
    which can be accessed using molecule().

    """

    # <<< Globals >>>

    # Has static information been initialized?
    initialized_shared = False
    # Global arrays of x, y, z exponents (Need libmint for max ang mom)
    LIBINT_MAX_AM = 6  # TODO
    exp_ao = [[] for l in range(LIBINT_MAX_AM)]

    def __init__(self, *args):

        # <<< Basic BasisSet Information >>>

        # The name of this basis set (e.g. "BASIS", "RI BASIS")
        self.name = None
        # Array of gaussian shells
        self.shells = None
        # Molecule object.
        self.molecule = None
        # Shell information
        self.atom_basis_shell = None

        # <<< Scalars >>>

        # Number of atomic orbitals (Cartesian)
        self.PYnao = None
        # Number of basis functions (either cartesian or spherical)
        self.PYnbf = None
        # The number of unique primitives
        self.n_uprimitive = None
        # The number of shells
        self.n_shells = None
        # The number of primitives
        self.PYnprimitive = None
        # The maximum angular momentum
        self.PYmax_am = None
        # The maximum number of primitives in a shell
        self.PYmax_nprimitive = None
        # Whether the basis set is uses spherical basis functions or not
        self.puream = None

        # <<< Arrays >>>

        # The number of primitives (and exponents) in each shell
        self.n_prim_per_shell = None
        # The first (Cartesian) atomic orbital in each shell
        self.shell_first_ao = None
        # The first (Cartesian / spherical) basis function in each shell
        self.shell_first_basis_function = None
        # Shell number to atomic center.
        self.shell_center = None
        # Which shell does a given (Cartesian / spherical) function belong to?
        self.function_to_shell = None
        # Which shell does a given Cartesian function belong to?
        self.ao_to_shell = None
        # Which center is a given function on?
        self.function_center = None
        # How many shells are there on each center?
        self.center_to_nshell = None
        # What's the first shell on each center?
        self.center_to_shell = None

        # The flattened lists of unique exponents
        self.uexponents = None
        # The flattened lists of unique contraction coefficients (normalized)
        self.ucoefficients = None
        # The flattened lists of unique contraction coefficients (as provided by the user)
        self.uoriginal_coefficients = None
        # The flattened lists of ERD normalized contraction coefficients
        self.uerd_coefficients = None
        # The flattened list of Cartesian coordinates for each atom
        self.xyz = None
        # label/basis to number of core electrons mapping for ECPs
        self.ecp_coreinfo = None

        # Divert to constructor functions
        if len(args) == 0:
            self.constructor_zero_ao_basis()
        elif len(args) == 2 and \
            isinstance(args[0], BasisSet) and \
            isinstance(args[1], int):
            self.constructor_basisset_center(*args)
        elif len(args) == 3 and \
            isinstance(args[0], str) and \
            isinstance(args[1], Molecule) and \
            isinstance(args[2], collections.OrderedDict):
            self.constructor_role_mol_shellmap(*args)
        elif len(args) == 4 and \
            isinstance(args[0], str) and \
            isinstance(args[1], Molecule) and \
            isinstance(args[2], collections.OrderedDict) and \
            isinstance(args[3], bool):
            self.constructor_role_mol_shellmap(*args)
        else:
            raise ValidationError('BasisSet::constructor: Inappropriate configuration of constructor arguments')

    def __eq__(self, other):
        """Naive equality test. Haven't considered exp/coeff distribution among shells or AM"""

        if isinstance(other, self.__class__):
            if ((self.name == other.name) and
                (self.puream == other.puream) and
                (self.PYnao == other.PYnao) and
                (self.PYnbf == other.PYnbf) and
                (self.n_prim_per_shell == other.n_prim_per_shell) and
                (self.ucoefficients == other.ucoefficients) and
                (self.uexponents == other.uexponents)):
                return True
            else:
                return False
        return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def allclose(self, other, atol: float=1.e-8, verbose: int=1):
        """Equality test. Sorts the coefficients so handles different shell orderings. Print any failed exp/coeff differences if verbose > 1."""
        sc, se = (list(t) for t in zip(*sorted(zip(self.uoriginal_coefficients, self.uexponents))))
        oc, oe = (list(t) for t in zip(*sorted(zip(other.uoriginal_coefficients, other.uexponents))))

        if isinstance(other, self.__class__):
            if ((self.name == other.name) and
                (self.puream == other.puream) and
                (self.PYnao == other.PYnao) and
                (self.PYnbf == other.PYnbf) and
                (self.n_prim_per_shell == other.n_prim_per_shell) and
                (all(abs(isc - ioc) < atol for isc, ioc in zip(sc, oc))) and
                (all(abs(ise - ioe) < atol for ise, ioe in zip(se, oe)))):
                return True
            else:
                if verbose > 1:
                    print("")
                    for idx in range(len(sc)):
                        if not((abs(sc[idx] - oc[idx]) < atol) and (abs(se[idx] - oe[idx]) < atol)):
                            print(f"{sc[idx]:20.12f} {oc[idx]:20.12f}\t\t{sc[idx] - oc[idx]:8.1E}\t\t\t{se[idx]:20.12f} {oe[idx]:20.12f}\t\t{se[idx] - oe[idx]:8.1E}")
                return False
        return False


    # <<< Methods for Construction >>>
    @classmethod
    def initialize_singletons(cls):
        """Initialize singleton values that are shared by all basis set objects."""
        # Populate the exp_ao arrays
        for l in range(cls.LIBINT_MAX_AM):
            for i in range(l + 1):
                x = l - i
                for j in range(i + 1):
                    y = i - j
                    z = j
                    cls.exp_ao[l].append([x, y, z])
        cls.initialized_shared = True

    def constructor_zero_ao_basis(self):
        """Constructs a zero AO basis set"""

        if not self.initialized_shared:
            self.initialize_singletons()

        # Add a dummy atom at the origin, to hold this basis function
        self.molecule = Molecule()
        self.molecule.add_atom(0, 0.0, 0.0, 0.0, 'X')
        # Fill with data representing a single S function, at the origin, with 0 exponent
        self.n_uprimitive = 1
        self.n_shells = 1
        self.PYnprimitive = 1
        self.PYnao = 1
        self.PYnbf = 1
        self.uerd_coefficients = [1.0]
        self.n_prim_per_shell = [1]
        self.uexponents = [0.0]
        self.ucoefficients = [1.0]
        self.uoriginal_coefficients = [1.0]
        self.shell_first_ao = [0]
        self.shell_first_basis_function = [0]
        self.ao_to_shell = [0]
        self.function_to_shell = [0]
        self.function_center = [0]
        self.shell_center = [0]
        self.center_to_nshell = [0]
        self.center_to_shell = [0]
        self.puream = False
        self.PYmax_am = 0
        self.PYmax_nprimitive = 1
        self.xyz = [0.0, 0.0, 0.0]
        self.name = '(Empty Basis Set)'
        self.shells = []
        self.shells.append(ShellInfo(0, self.uoriginal_coefficients,
            self.uexponents, 'Cartesian', 0, self.xyz, 0))

    def constructor_role_mol_shellmap(self, role, mol, shell_map, is_ecp = False):
        """The most commonly used constructor. Extracts basis set name for *role*
        from each atom of *mol*, looks up basis and role entries in the
        *shell_map* dictionary, retrieves the ShellInfo objects and returns
        the BasisSet.

        """
        self.molecule = mol
        self.name = role
        self.xyz = self.molecule.geometry()  # not used in libmints but this seems to be the intent
        self.atom_basis_shell = shell_map
        natom = self.molecule.natom()

        # Singletons
        if not self.initialized_shared:
            self.initialize_singletons()

        # These will tell us where the primitives for [basis][symbol] start and end in the compact array
        primitive_start = {}
        primitive_end = {}

        # First, loop over the unique primitives, and store them
        uexps = []
        ucoefs = []
        uoriginal_coefs = []
        uerd_coefs = []
        rpowers = []
        self.n_uprimitive = 0

        for symbolfirst, symbolsecond in shell_map.items():
            label = symbolfirst
            basis_map = symbolsecond
            primitive_start[label] = {}
            primitive_end[label] = {}
            for basisfirst, basissecond in basis_map.items():
                basis = basisfirst
                shells = basis_map[basis]  # symbol --> label
                primitive_start[label][basis] = self.n_uprimitive  # symbol --> label
                for i in range(len(shells)):
                    shell = shells[i]
                    for prim in range(shell.nprimitive()):
                        rpowers.append(shell.rpower(prim))
                        uexps.append(shell.exp(prim))
                        ucoefs.append(shell.coef(prim))
                        uoriginal_coefs.append(shell.original_coef(prim))
                        uerd_coefs.append(shell.erd_coef(prim))
                        self.n_uprimitive += 1
                primitive_end[label][basis] = self.n_uprimitive  # symbol --> label

        # Count basis functions, shells and primitives
        self.n_shells = 0
        self.PYnprimitive = 0
        self.PYnao = 0
        self.PYnbf = 0
        for n in range(natom):
            atom = self.molecule.atom_entry(n)
            basis = atom.basisset(role)
            label = atom.label()  # symbol --> label
            shells = shell_map[label][basis]  # symbol --> label
            for i in range(len(shells)):
                shell = shells[i]
                nprim = shell.nprimitive()
                self.PYnprimitive += nprim
                self.n_shells += 1
                self.PYnao += shell.ncartesian()
                self.PYnbf += shell.nfunction()

        # Allocate arrays
        self.n_prim_per_shell = [0] * self.n_shells
        # The unique primitives
        self.uexponents = [0.0] * self.n_uprimitive
        self.ucoefficients = [0.0] * self.n_uprimitive
        self.uoriginal_coefficients = [0.0] * self.n_uprimitive
        self.uerd_coefficients = [0.0] * self.n_uprimitive
        for i in range(self.n_uprimitive):
            self.uexponents[i] = uexps[i]
            self.ucoefficients[i] = ucoefs[i]
            self.uoriginal_coefficients[i] = uoriginal_coefs[i]
            self.uerd_coefficients[i] = uerd_coefs[i]

        self.shell_first_ao = [0] * self.n_shells
        self.shell_first_basis_function = [0] * self.n_shells
        self.shells = [None] * self.n_shells
        self.ao_to_shell = [0] * self.PYnao
        self.function_to_shell = [0] * self.PYnbf
        self.function_center = [0] * self.PYnbf
        self.shell_center = [0] * self.n_shells
        self.center_to_nshell = [0] * natom
        self.center_to_shell = [0] * natom

        # Now loop over all atoms, and point to the appropriate unique data
        shell_count = 0
        ao_count = 0
        bf_count = 0
        xyz_ptr = [0.0, 0.0, 0.0]  # libmints seems to be always passing ShellInfo zeros, so following suit
        self.puream = False
        self.PYmax_am = 0
        self.PYmax_nprimitive = 0
        for n in range(natom):
            atom = self.molecule.atom_entry(n)
            basis = atom.basisset(role)
            label = atom.label()  # symbol --> label
            shells = shell_map[label][basis]  # symbol --> label
            ustart = primitive_start[label][basis]  # symbol --> label
            uend = primitive_end[label][basis]  # symbol --> label
            nshells = len(shells)
            self.center_to_nshell[n] = nshells
            self.center_to_shell[n] = shell_count
            atom_nprim = 0
            for i in range(nshells):
                thisshell = shells[i]
                self.shell_first_ao[shell_count] = ao_count
                self.shell_first_basis_function[shell_count] = bf_count
                shell_nprim = thisshell.nprimitive()
                am = thisshell.am()
                self.PYmax_nprimitive = max(shell_nprim, self.PYmax_nprimitive)
                self.PYmax_am = max(am, self.PYmax_am)
                self.shell_center[shell_count] = n
                self.puream = thisshell.is_pure()
                tst = ustart + atom_nprim
                tsp = ustart + atom_nprim + shell_nprim
                self.shells[shell_count] = ShellInfo(am,
                    self.uoriginal_coefficients[tst:tsp],
                    self.uexponents[tst:tsp],
                    'Pure' if self.puream else 'Cartesian',
                    n, xyz_ptr, bf_count, pt='Normalized' if is_ecp else 'Unnormalized',
                    rpowers=rpowers[tst:tsp])
                for thisbf in range(thisshell.nfunction()):
                    self.function_to_shell[bf_count] = shell_count
                    self.function_center[bf_count] = n
                    bf_count += 1
                for thisao in range(thisshell.ncartesian()):
                    self.ao_to_shell[ao_count] = shell_count
                    ao_count += 1
                atom_nprim += shell_nprim
                shell_count += 1

            if atom_nprim != uend - ustart:
                raise ValidationError("Problem with nprimitive in basis set construction!")

    def constructor_basisset_center(self, bs, center):
        """
        * Creates a new basis set object for an atom, from an existing basis set
        * bs: the basis set to copy data from
        * center: the atom in bs to copy over

        """
        # Singletons; these should've been initialized by this point, but just in case
        if not self.initialized_shared:
            self.initialize_singletons()

        # First, find the shells we need, and grab the data
        uexps = []
        ucoefs = []
        uoriginal_coefs = []
        uerd_coefs = []
        self.name = bs.name
        self.n_shells = 0
        self.n_uprimitive = 0
        self.PYnao = 0
        self.PYnbf = 0
        for shelln in range(bs.nshell()):
            shell = bs.shell(shelln)
            if shell.ncenter() == center:
                nprim = shell.nprimitive()
                for prim in range(nprim):
                    uexps.append(shell.exp(prim))
                    ucoefs.append(shell.coef(prim))
                    uoriginal_coefs.append(shell.original_coef(prim))
                    uerd_coefs.append(shell.erd_coef(prim))
                    self.n_uprimitive += 1
                self.n_shells += 1
                self.PYnao += shell.ncartesian()
                self.PYnbf += shell.nfunction()
        self.PYnprimitive = self.n_uprimitive

        # Create a "molecule", i.e., an atom, with 1 fragment
        mol = bs.molecule
        self.molecule = Molecule.from_arrays(elem=[mol.symbol(center)],
                                             geom=mol.xyz(center),
                                             units='Bohr',
                                             fix_com=True,
                                             verbose=0)
        # Allocate arrays
        self.n_prim_per_shell = [0] * self.n_shells
        # The unique primitives
        self.uexponents = [0.0] * self.n_uprimitive
        self.ucoefficients = [0.0] * self.n_uprimitive
        self.uoriginal_coefficients = [0.0] * self.n_uprimitive
        self.uerd_coefficients = [0.0] * self.n_uprimitive
        for i in range(self.n_uprimitive):
            self.uexponents[i] = uexps[i]
            self.ucoefficients[i] = ucoefs[i]
            self.uoriginal_coefficients[i] = uoriginal_coefs[i]
            self.uerd_coefficients[i] = uerd_coefs[i]

        self.shell_first_ao = [0] * self.n_shells
        self.shell_first_basis_function = [0] * self.n_shells
        self.shells = [None] * self.n_shells
        self.ao_to_shell = [0] * self.PYnao
        self.function_to_shell = [0] * self.PYnbf
        self.function_center = [0] * self.PYnbf
        self.shell_center = [0] * self.n_shells
        self.center_to_nshell = [0]
        self.center_to_shell = [0]
        self.xyz = mol.xyz(center)

        # Now loop over shell for this atom, and point to the appropriate unique data
        shell_count = 0
        ao_count = 0
        bf_count = 0
        self.puream = False
        self.PYmax_am = 0
        self.PYmax_nprimitive = 0
        prim_count = 0
        for shelln in range(bs.nshell()):
            shell = bs.shell(shelln)
            if shell.ncenter() == center:
                self.center_to_nshell[0] = self.n_shells
                #self.center_to_shell[0] = shell_count  # diff from libmints
                self.shell_first_ao[shell_count] = ao_count
                self.shell_first_basis_function[shell_count] = bf_count
                shell_nprim = shell.nprimitive()
                am = shell.am()
                self.PYmax_nprimitive = shell_nprim if shell_nprim > self.PYmax_nprimitive else self.PYmax_nprimitive
                self.PYmax_am = max(self.PYmax_am, am)
                self.shell_center[shell_count] = center
                self.puream = shell.is_pure()
                tst = prim_count
                tsp = prim_count + shell_nprim
                self.shells[shell_count] = ShellInfo(am,
                    self.uoriginal_coefficients[tst:tsp],
                    self.uexponents[tst:tsp],
                    'Pure' if self.puream else 'Cartesian',
                    center, self.xyz, bf_count, pt='Unnormalized', rpowers=None)
                self.shells[shell_count].pyprint()
                for thisbf in range(shell.nfunction()):
                    self.function_to_shell[bf_count] = shell_count
                    self.function_center[bf_count] = center
                    bf_count += 1
                for thisao in range(shell.ncartesian()):
                    self.ao_to_shell[ao_count] = shell_count
                    ao_count += 1
                shell_count += 1
                prim_count += shell_nprim

    # <<< Methods for Construction by Another Name >>>

    @staticmethod
    def zero_ao_basis_set():
        """Returns an empty basis set object.
        Returns a BasisSet object that actually has a single s-function
        at the origin with an exponent of 0.0 and contraction of 1.0.
        *  @return A new empty BasisSet object.

        """
        # In the new implementation, we simply call the default constructor
        return BasisSet()

    def atomic_basis_set(self, center):
        """Return a BasisSet object containing all shells at center i (0-index)
        * Used for Atomic HF computations for SAD Guesses
        * @param center Atomic center to provide a basis object for.
        * @returns A new basis set object for the atomic center.

        """
        return BasisSet(self, center)

    @staticmethod
    def build(molecule, shells):
        """Builder factory method
        * @param molecule the molecule to build the BasisSet around
        * @param shells array of *atom-numbered* ShellInfo to build the BasisSet from
        * @return BasisSet corresponding to this molecule and set of shells

        """
        raise FeatureNotImplemented('BasisSet::build')

    @staticmethod
    def pyconstruct_combined(mol, keys, targets, fitroles, others):
        """Builds a BasisSet object for *mol* (either a qcdb.Molecule or
        a string that can be instantiated into one) from two basis set
        specifications passed in as python functions or as a string that
        names a basis to be applied to all atoms. The union of the two
        basis sets can be used as an auxiliary basis set for F12 methods.

        Parameters
        ----------
        mol : :py:class:`qcdb.Molecule` or dict or str
            If not a :py:class:`qcdb.Molecule`, something that can be converted into one.
            If the latter, the basisset dict is returned rather than the BasisSet object.
            If you've got a psi4.core.Molecule, pass `qcdb.Molecule(psimol.to_dict())`.
        keys : {'BASIS', 'CABS_BASIS'}
            Labels to append to the two basis sets.
        targets : str or function
            Defines the basis sets to be combined.
        fitroles : {'ORBITAL', 'F12'}
        others :
            Like `target` only supplies the definitions
            for the orbital basis.

        Returns
        -------
        bas : :py:class:`qcdb.BasisSet`
            Built BasisSet object for `mol`.
        dbas : dict, optional
            Only returned if `mol` is a qcdb.Molecule. Suitable for feeding to libmints.
        """
        # make sure the lengths are all the same
        if len(keys) != len(targets) or len(keys) != len(fitroles):
            raise ValidationError("""Lengths of keys, targets, and fitroles must be equal""")

        # Create (if necessary) and update qcdb.Molecule
        if isinstance(mol, str):
            mol = Molecule(mol)
            returnBasisSet = False
        elif isinstance(mol, Molecule):
            returnBasisSet = True
        else:
            raise ValidationError("""Argument mol must be psi4string or qcdb.Molecule""")
        mol.update_geometry()

        # load in the basis sets
        sets = [BasisSet.pyconstruct(mol, keys[at], targets[at], fitroles[at], others[at]) for at in range(len(keys))]
        name = " + ".join(targets)
        keywords = " + ".join(keys)
        blends = " + ".join([bas.name.upper() for bas in sets])

        # merges the two maps of ShellInfo into one combined map
        combined_atom_basis_shell = collections.OrderedDict()
        for basis in sets:
            atom_basis_shell = basis.atom_basis_shell

            for label, basis_map in atom_basis_shell.items():
                if label not in combined_atom_basis_shell:
                    combined_atom_basis_shell[label] = collections.OrderedDict.fromkeys([name], [])
                for shells in basis_map.values():
                    combined_atom_basis_shell[label][name].extend(shells)

        # sort the shells by angular momentum
        for label, basis_map in combined_atom_basis_shell.items():
            combined_atom_basis_shell[label][name] = sorted(combined_atom_basis_shell[label][name],
                    key=lambda shell: shell.l)

        # Molecule and parser prepped, call the constructor
        mol.set_basis_all_atoms(name, "CABS")

        # Construct the grand BasisSet for mol
        basisset = BasisSet("CABS", mol, combined_atom_basis_shell)

        # Construct all the one-atom BasisSet-s for mol's CoordEntry-s
        for at in range(mol.natom()):
            oneatombasis = BasisSet(basisset, at)
            oneatombasishash = hashlib.sha1(oneatombasis.print_detail(numbersonly=True).encode('utf-8')).hexdigest()
            mol.set_shell_by_number(at, oneatombasishash, role="CABS")
        mol.update_geometry()  # re-evaluate symmetry taking basissets into account

        text = """   => Creating Basis Set <=\n\n"""
        text += """    Role: %s\n""" % (fitroles)
        text += """    Keyword: %s\n""" % (keys)
        text += """    Name: %s\n""" % (name)

        if returnBasisSet:
            print(text)
            return basisset
        else:
            bsdict = {}
            bsdict['message'] = text
            bsdict['name'] = basisset.name
            bsdict['key'] = keywords
            bsdict['blend'] = blends
            bsdict['puream'] = int(basisset.has_puream())
            bsdict['shell_map'] = basisset.export_for_libmints("CABS")
            return bsdict

    @staticmethod
    def pyconstruct(mol, key, target, fitrole='ORBITAL', other=None, return_atomlist=False, return_dict=False, verbose=1):
        """Builds a BasisSet object for *mol* (either a qcdb.Molecule or
        a string that can be instantiated into one) from basis set
        specifications passed in as python functions or as a string that
        names a basis to be applied to all atoms. Always required is the
        keyword *key* and string/function *target* of the basis to be
        constructed. For orbital basis sets, *key* will likely be 'BASIS'
        and, together with *target*, these arguments suffice.
        ``pyconstruct(smol, "BASIS", basisspec_psi4_yo_631pg_d_p_)``
        ``pyconstruct(mol, "BASIS", "6-31+G**")``
        When building an auxiliary basis, *key* is again the keyword,
        *target* is the string or function for the fitting basis (this
        may also be an empty string). In case the fitting basis isn't
        fully specified, also provide a *fitrole* and the string/function
        of the orbital basis as *other*, so that orbital hints can be
        used to look up a suitable default basis in BasisFamily.
        ``pyconstruct(smol, "DF_BASIS_MP2", basisspec_psi4_yo_ccpvdzri, 'RIFIT', basisspec_psi4_yo_631pg_d_p_)``
        ``pyconstruct(mol, "DF_BASIS_MP2", "", "RIFIT", "6-31+G(d,p)")``

        Parameters
        ----------
        mol : :py:class:`qcdb.Molecule` or dict or str
            If not a :py:class:`qcdb.Molecule`, something that can be converted into one.
            If the latter, the basisset dict is returned rather than the BasisSet object.
            If you've got a psi4.core.Molecule, pass `qcdb.Molecule(psimol.to_dict())`.
        key : {'BASIS', 'DF_BASIS_SCF', 'DF_BASIS_MP2', 'DF_BASIS_CC'}
            Label (effectively Psi4 keyword) to append the basis on the molecule.
        target : str or function
            Defines the basis set to be constructed. Can be a string (naming a
            basis file) or a function (multiple files, shells).
            Required for `fitrole='ORBITAL'`. For auxiliary bases to be built
            entirely from orbital default, can be empty string.
        fitrole : {'ORBITAL', 'JKFIT', 'RIFIT'}
        other : 
            Only used when building fitting bases. Like `target` only supplies
            the definitions for the orbital basis.
        return_atomlist : bool, optional
            Build one-atom basis sets (for SAD) rather than one whole-Mol basis set
        return_dict : bool, optional
            Additionally return the dictionary representation of built BasisSet

        Returns
        -------
        bas : :py:class:`qcdb.BasisSet`
            Built BasisSet object for `mol`.
        dbas : dict, optional
            Only returned if `return_dict=True`. Suitable for feeding to libmints.
            
        """
        orbonly = True if (fitrole == 'ORBITAL' and other is None) else False
        if orbonly:
            orb = target
            aux = None
        else:
            orb = other
            aux = target

        if verbose >= 2:
            print('BasisSet::pyconstructP', 'key =', key, 'aux =', aux, 'fitrole =', fitrole, 'orb =', orb, 'orbonly =', orbonly) #, mol)

        # Create (if necessary) and update qcdb.Molecule
        if not isinstance(mol, Molecule):
            mol = Molecule(mol)
        mol.update_geometry()

        # Apply requested basis set(s) to the molecule
        #   - basstrings only a temp object so using fitrole as dict key instead of psi4 keyword
        #   - error checking not needed since C-side already checked for NULL ptr
        mol.clear_basis_all_atoms()
        # TODO now need to clear shells, too
        basstrings = collections.defaultdict(dict)
        if orb is None or orb == '':
            raise ValidationError("""Orbital basis argument must not be empty.""")
        elif callable(orb):
            basstrings['BASIS'] = orb(mol, 'BASIS')
            callby = orb.__name__.replace('basisspec_psi4_yo__', '')
        elif orb in basishorde:
            basstrings['BASIS'] = basishorde[orb](mol, 'BASIS')
            callby = orb
        elif isinstance(orb, str):
            mol.set_basis_all_atoms(orb, role='BASIS')
            callby = orb
        else:
            raise ValidationError("""Orbital basis argument must be function that applies basis sets to Molecule or a string of the basis to be applied to all atoms.""")

        if not orbonly:
            if aux is None or aux == '':
                callby = '({} Aux)'.format(callby)
            elif callable(aux):
                basstrings[fitrole] = aux(mol, fitrole)
                callby = aux.__name__.replace('basisspec_psi4_yo__', '')
            elif isinstance(aux, str):
                mol.set_basis_all_atoms(aux, role=fitrole)
                callby = aux
            else:
                raise ValidationError("""Auxiliary basis argument must be function that applies basis sets to Molecule or a string of the basis to be applied to all atoms.""")

        # Not like we're ever using a non-G94 format
        parser = Gaussian94BasisSetParser()

        # Molecule and parser prepped, call the constructor
        bs, msg, ecp = BasisSet.construct(parser, mol,
                                          'BASIS' if fitrole == 'ORBITAL' else fitrole,
                                          None if fitrole == 'ORBITAL' else fitrole,
                                          basstrings['BASIS'] if fitrole == 'ORBITAL' else basstrings[fitrole],
                                          return_atomlist=return_atomlist,
                                          verbose=verbose)

        text = """   => Loading Basis Set <=\n\n"""
        text += """    Name: %s\n""" % (callby.upper())
        text += """    Role: %s\n""" % (fitrole)
        text += """    Keyword: %s\n""" % (key)
        text += msg

        if return_atomlist:
            if return_dict:
                atom_basis_list = []
                for atbs in bs:
                    bsdict = {}
                    bsdict['message'] = text
                    bsdict['key'] = key
                    bsdict['name'] = callby.upper()
                    bsdict['blend'] = atbs.name.upper()
                    bsdict['puream'] = int(atbs.has_puream())
                    bsdict['shell_map'] = atbs.export_for_libmints('BASIS' if fitrole == 'ORBITAL' else fitrole)
                    if ecp:
                        bsdict['ecp_shell_map'] = ecp.export_for_libmints('BASIS')
                    bsdict['molecule'] = atbs.molecule.to_dict(force_c1=True)
                    atom_basis_list.append(bsdict)
                return bs, atom_basis_list
            else:  
                return bs

        if return_dict:
            bsdict = {}
            bsdict['message'] = text
            bsdict['key'] = key
            bsdict['name'] = callby.upper()
            bsdict['blend'] = bs.name.upper()
            bsdict['puream'] = int(bs.has_puream())
            bsdict['shell_map'] = bs.export_for_libmints('BASIS' if fitrole == 'ORBITAL' else fitrole)
            if ecp:
                bsdict['ecp_shell_map'] = ecp.export_for_libmints('BASIS')
            return bs, bsdict
        else:
            return bs

    @classmethod
    def construct(cls, parser, mol, key, deffit=None, basstrings=None, return_atomlist=False, verbose=1):
        """Returns a new BasisSet object configured from the *mol*
        Molecule object for *key* (generally a Psi4 keyword: BASIS,
        DF_BASIS_SCF, etc.). Fails utterly if a basis has not been set for
        *key* for every atom in *mol*, unless *deffit* is set (JFIT,
        JKFIT, or RIFIT), whereupon empty atoms are assigned to *key*
        from the :py:class:`~BasisFamily`. This function is significantly
        re-worked from its libmints analog.

        Parameters
        ----------
        mol : qcdb.Molecule
            A molecule object for which every atom has had a basisset set for `key`

        basstrings : dict, optional
            Additional source for basis data. Keys are regularized basis names and values are gbs strings.
        return_atomlist
            Return list of one-atom BasisSet-s, rather than single whole-mol BasisSet.
    
        """
        # Update geometry in molecule, if there is a problem an exception is thrown.
        mol.update_geometry()

        # Paths to search for gbs files: here + PSIPATH + library
        try:
            from psi4 import core
        except ImportError:
            libraryPath = ''
        else:
            psidatadir = core.get_datadir()
            libraryPath = os.pathsep + os.path.join(os.path.abspath(psidatadir), 'basis')
        #nolongerenvvar psidatadir = os.environ.get('PSIDATADIR', None)
        #nolongerpredicatble psidatadir = __file__ + '/../../..' if psidatadir is None else psidatadir

        basisPath = os.path.abspath('.') + os.pathsep
        basisPath += os.pathsep.join([os.path.abspath(x) for x in os.environ.get('PSIPATH', '').split(os.pathsep)])
        basisPath += libraryPath

        # Validate deffit for key
        univdef_zeta = 4
        univdef = {'JFIT': ('def2-universal-jfit', 'def2-universal-jfit', None),
                   'JKFIT': ('def2-universal-jkfit', 'def2-universal-jkfit', None),
                   'RIFIT': ('def2-qzvpp-ri', 'def2-qzvpp-ri', None),
                   'DECON': (None, None, BasisSet.decontract),
                   'F12': ('def2-qzvpp-f12', 'def2-qzvpp-f12', None)}

        if deffit is not None:
            if deffit not in univdef.keys():
                raise ValidationError("""BasisSet::construct: deffit argument invalid: %s""" % (deffit))

        # Map of ShellInfo
        atom_basis_shell = collections.OrderedDict()
        ecp_atom_basis_shell = collections.OrderedDict()
        ecp_atom_basis_ncore = collections.OrderedDict()
        names = {}
        summary = []
        bastitles = []

        for at in range(mol.natom()):
            symbol = mol.atom_entry(at).symbol()  # O, He
            label = mol.atom_entry(at).label()  # O3, C_Drot, He
            basdict = mol.atom_entry(at).basissets()  # {'BASIS': 'sto-3g', 'DF_BASIS_MP2': 'cc-pvtz-ri'}

            if label not in atom_basis_shell:
                atom_basis_shell[label] = collections.OrderedDict()
            if label not in ecp_atom_basis_shell:
                ecp_atom_basis_shell[label] = collections.OrderedDict()

            # Establish search parameters for what/where basis entries suitable for atom
            seek = {}
            try:
                requested_basname = basdict[key]
            except KeyError:
                if key == 'BASIS' or deffit is None:
                    raise BasisSetNotDefined("""BasisSet::construct: No basis set specified for %s and %s.""" %
                        (symbol, key))
                else:
                    # No auxiliary / decon basis set for atom, so try darnedest to find one.
                    #   This involves querying the BasisFamily for default and
                    #   default-default and finally the universal default (defined
                    #   in this function). Since user hasn't indicated any specifics,
                    #   look for symbol only, not label.
                    tmp = []
                    tmp.append(corresponding_basis(basdict['BASIS'], deffit))
                    #NYI#tmp.append(corresponding_basis(basdict['BASIS'], deffit + '-DEFAULT'))
                    orbital_zeta = corresponding_zeta(basdict['BASIS'])
                    if orbital_zeta is None or orbital_zeta <= univdef_zeta:
                        tmp.append(univdef[deffit])
                    seek['basis'] = [item for item in tmp if item != (None, None, None)]
                    seek['entry'] = [symbol]
                    seek['path'] = basisPath
                    seek['strings'] = ''
            else:
                # User (I hope ... dratted has_changed) has set basis for atom,
                #   so look only for basis (don't try defaults), look for label (N88)
                #   or symbol (N) (in that order; don't want to restrict use of atom
                #   labels to basis set spec), look everywhere (don't just look
                #   in library)
                if requested_basname.lower().endswith('-decon'):
                    bas_recipe = requested_basname, requested_basname[:-6], BasisSet.decontract
                else:
                    bas_recipe = requested_basname, requested_basname, None
                seek['basis'] = [bas_recipe]
                seek['entry'] = [symbol] if symbol == label else [label, symbol]
                seek['path'] = basisPath
                seek['strings'] = '' if basstrings is None else list(basstrings.keys())

            if verbose >= 2:
                print("""  Shell Entries: %s""" % (seek['entry']))
                print("""  Basis Sets: %s""" % (seek['basis']))
                print("""  File Path: %s""" % (', '.join(map(str, seek['path'].split(os.pathsep)))))
                print("""  Input Blocks: %s\n""" % (', '.join(seek['strings'])))

            # Search through paths, bases, entries
            for bas in seek['basis']:
                (bastitle, basgbs, postfunc) = bas

                filename = cls.make_filename(basgbs)

                # If name is prefixed with "bse:", then load from basis set exchange library
                if bastitle.lower().startswith('bse:'):
                    if not qcel.util.which_import('basis_set_exchange', return_bool=True):
                        raise ModuleNotFoundError("Python module 'basis_set_exchange' not found. Solve by installing it: `conda install -c conda-forge basis_set_exchange` or `pip install basis_set_exchange`")

                    import basis_set_exchange as bse

                    index = 'bse %s' % (bastitle)

                    if index not in names:
                        basparts = bastitle.split(":")

                        # If len=3, then bse:name:version
                        # else it's just bse:name
                        if len(basparts) == 3:
                            _, basname, basver = basparts
                        else:
                            _, basname = basparts
                            basver = None

                        names[index] = bse.get_basis(basname, version=basver, fmt='psi4')

                # -- seek bas string in input file strings
                elif filename[:-4] in seek['strings']:
                    index = 'inputblock %s' % (filename[:-4])
                    # Store contents
                    if index not in names:
                        names[index] = basstrings[filename[:-4]].split('\n')
                else:
                    # -- Else seek bas.gbs file in path
                    fullfilename = search_file(_basis_file_warner_and_aliaser(filename), seek['path'])
                    if fullfilename is None:
                        # -- Else skip to next bas
                        continue
                    # Store contents so not reloading files
                    index = 'file %s' % (fullfilename)
                    if index not in names:
                        names[index] = parser.load_file(fullfilename)

                lines = names[index]

                for entry in seek['entry']:

                    # Seek entry in lines, else skip to next entry
                    shells, msg, ecp_shells, ecp_msg, ecp_ncore = parser.parse(entry, lines)
                    if shells is None:
                        continue

                    # Found!
                    # -- Post-process
                    if postfunc:
                        shells = postfunc(shells)
                        fmsg = 'func {}'.format(postfunc.__name__)
                    else:
                        fmsg = ''
                    # -- Assign to Molecule
                    atom_basis_shell[label][bastitle] = shells
                    ecp_atom_basis_shell[label][bastitle] = ecp_shells
                    if key == 'BASIS':
                        ecp_atom_basis_ncore[label] = ecp_ncore
                    mol.set_basis_by_number(at, bastitle, role=key)
                    bastitles.append(bastitle.upper())
                    if ecp_msg:
                        summary.append("""entry %-10s %s (ECP: %s) %s %s""" % (entry, msg, ecp_msg, index, fmsg))
                    else:
                        summary.append("""entry %-10s %s %s %s""" % (entry, msg, index, fmsg))
                    break

                # Break from outer loop if inner loop breaks
                else:
                    continue
                break

            else:
                text2  = """  Shell Entries: %s\n""" % (seek['entry'])
                text2 += """  Basis Sets: %s\n""" % (seek['basis'])
                text2 += """  File Path: %s\n""" % (', '.join(map(str, seek['path'].split(os.pathsep))))
                text2 += """  Input Blocks: %s\n""" % (', '.join(seek['strings']))

                # basis/NOTES [53]
                deprecated_basis_sets = ["APR-CC-PCV(5+D)Z", "APR-CC-PCV(6+D)Z", "APR-CC-PCV(Q+D)Z", "AUG-CC-PCV(5+D)Z", "AUG-CC-PCV(6+D)Z", "AUG-CC-PCV(D+D)Z", "AUG-CC-PCV(Q+D)Z", "AUG-CC-PCV(T+D)Z", "CC-PCV(5+D)Z", "CC-PCV(6+D)Z", "CC-PCV(D+D)Z", "CC-PCV(Q+D)Z", "CC-PCV(T+D)Z", "FEB-CC-PCV(6+D)Z", "HEAVY-AUG-CC-PCV(5+D)Z", "HEAVY-AUG-CC-PCV(6+D)Z", "HEAVY-AUG-CC-PCV(D+D)Z", "HEAVY-AUG-CC-PCV(Q+D)Z", "HEAVY-AUG-CC-PCV(T+D)Z", "JUN-CC-PCV(5+D)Z", "JUN-CC-PCV(6+D)Z", "JUN-CC-PCV(D+D)Z", "JUN-CC-PCV(Q+D)Z", "JUN-CC-PCV(T+D)Z", "MAR-CC-PCV(5+D)Z", "MAR-CC-PCV(6+D)Z", "MAY-CC-PCV(5+D)Z", "MAY-CC-PCV(6+D)Z", "MAY-CC-PCV(Q+D)Z", "MAY-CC-PCV(T+D)Z", "APR-CC-PWCV(5+D)Z", "APR-CC-PWCV(Q+D)Z", "AUG-CC-PWCV(5+D)Z", "AUG-CC-PWCV(D+D)Z", "AUG-CC-PWCV(Q+D)Z", "AUG-CC-PWCV(T+D)Z", "CC-PWCV(5+D)Z", "CC-PWCV(D+D)Z", "CC-PWCV(Q+D)Z", "CC-PWCV(T+D)Z", "HEAVY-AUG-CC-PWCV(5+D)Z", "HEAVY-AUG-CC-PWCV(D+D)Z", "HEAVY-AUG-CC-PWCV(Q+D)Z", "HEAVY-AUG-CC-PWCV(T+D)Z", "JUN-CC-PWCV(5+D)Z", "JUN-CC-PWCV(D+D)Z", "JUN-CC-PWCV(Q+D)Z", "JUN-CC-PWCV(T+D)Z", "MAR-CC-PWCV(5+D)Z", "MAY-CC-PWCV(5+D)Z", "MAY-CC-PWCV(Q+D)Z", "MAY-CC-PWCV(T+D)Z"]
                if  seek['basis'][0][0] in deprecated_basis_sets:
                    raise BasisSetNotFoundDeprecated(
                        f'BasisSet::construct: Unable to find a basis set for atom {at + 1} for key {key} among:\n{text2}',
                        f"This basis was ill-advised. Access if you must via {seek['basis'][0][0].lower()}-deprecated. See note [54] for details in https://github.com/psi4/psi4/blob/master/psi4/share/psi4/basis/NOTES.\n"
                    )
                else:
                    # Ne'er found :-(
                    raise BasisSetNotFound(f'BasisSet::construct: Unable to find a basis set for atom {at + 1} for key {key} among:\n{text2}')

        # Construct the grand BasisSet for mol
        basisset = BasisSet(key, mol, atom_basis_shell)

        # If an ECP was detected, and we're building BASIS, process it now
        ecpbasisset = None
        ncore = 0
        for atom in ecp_atom_basis_ncore:
            ncore += ecp_atom_basis_ncore[atom]
        if ncore and key == 'BASIS':
            ecpbasisset = BasisSet(key, mol, ecp_atom_basis_shell, True)
            ecpbasisset.ecp_coreinfo = ecp_atom_basis_ncore

        # Construct all the one-atom BasisSet-s for mol's CoordEntry-s
        atom_basis_list = []
        for at in range(mol.natom()):
            oneatombasis = BasisSet(basisset, at)
            oneatombasishash = hashlib.sha1(oneatombasis.print_detail(numbersonly=True).encode('utf-8')).hexdigest()
            if return_atomlist:
                oneatombasis.molecule.set_shell_by_number(0, oneatombasishash, role=key)
                atom_basis_list.append(oneatombasis)
            mol.set_shell_by_number(at, oneatombasishash, role=key)

        mol.update_geometry()  # re-evaluate symmetry taking basissets into account

        bastitles = list(set(bastitles))
        bastitles.sort()
        basisset.name = ' + '.join(bastitles)

        # Summary printing
        tmp = collections.defaultdict(list)
        for at, v in enumerate(summary):
            tmp[v].append(at + 1)
        tmp2 = collections.OrderedDict()
        maxsats = 0
        for item in sorted(tmp.values()):
            for msg, ats in tmp.items():
                if item == ats:
                    G = (list(x) for _, x in itertools.groupby(ats, lambda x, c=itertools.count(): next(c) - x))
                    sats = ", ".join("-".join(map(str, (g[0], g[-1])[:len(g)])) for g in G)
                    maxsats = max(maxsats, len(sats))
                    tmp2[sats] = msg
        text = ''
        for ats, msg in tmp2.items():
            text += """    atoms %s %s\n""" % (ats.ljust(maxsats), msg)
        text += '\n'

        if return_atomlist:
            return atom_basis_list, text, ecpbasisset
        else:
            return basisset, text, ecpbasisset

    # <<< Simple Methods for Basic BasisSet Information >>>

    def name(self):
        """Returns the name of this basis set"""
        return self.name

    def set_name(self, name):
        """Sets the name of this basis set"""
        self.name = name

# JET added but I think should fail
#+    def atom_shell_map(self):
#+        return self.atom_shell_map

    def nprimitive(self):
        """Number of primitives.
        *  @return The total number of primitives in all contractions.

        """
        return self.PYnprimitive

    def max_nprimitive(self):
        """Maximum number of primitives in a shell.
        *  Examines each shell and find the shell with the maximum number of primitives returns that
        *  number of primitives.
        *  @return Maximum number of primitives.

        """
        return self.PYmax_nprimitive

    def nshell(self):
        """Number of shells.
        *  @return Number of shells.

        """
        return self.n_shells

    def nao(self):
        """Number of atomic orbitals (Cartesian).
        * @return The number of atomic orbitals (Cartesian orbitals, always).

        """
        return self.PYnao

    def nbf(self):
        """Number of basis functions (Spherical).
        *  @return The number of basis functions (Spherical, if has_puream() == true).

        """
        return self.PYnbf

    def max_am(self):
        """Maximum angular momentum used in the basis set.
        *  @return Maximum angular momentum.

        """
        return self.PYmax_am

    def has_puream(self):
        """Spherical harmonics?
        *  @return true if using spherical harmonics

        """
        return self.puream

    def max_function_per_shell(self):
        """Compute the maximum number of basis functions contained in a shell.
        *  @return The max number of basis functions in a shell.

        """
        return 2 * self.PYmax_am + 1 if self.puream else (self.PYmax_am + 1) * (self.PYmax_am + 2) / 2

    def molecule(self):
        """Molecule this basis is for.
        *  @return Shared pointer to the molecule for this basis set.

        """
        return self.molecule

    def shell_to_ao_function(self, i):
        """Given a shell what is its first AO function
        *  @param i Shell number
        *  @return The function number for the first function for the i'th shell.

        """
        return self.shell_first_ao[i]

    def shell_to_center(self, i):
        """Given a shell what is its atomic center
        *  @param i Shell number
        *  @return The atomic center for the i'th shell.

        """
        return self.shell_center[i]

    def shell_to_basis_function(self, i):
        """Given a shell what is its first basis function (spherical) function
        *  @param i Shell number
        *  @return The function number for the first function for the i'th shell.

        """
        return self.shell_first_basis_function[i]

    def function_to_shell(self, i):
        """Given a function number what shell does it correspond to."""
        return self.function_to_shell[i]

    def function_to_center(self, i):
        """Given a function what is its atomic center
        *  @param i Function number
        *  @return The atomic center for the i'th function.

        """
        return self.function_center[i]

    def ao_to_shell(self, i):
        """Given a Cartesian function (AO) number what shell does it correspond to."""
        return self.ao_to_shell[i]

    def shell(self, si, center=None):
        """Return the si'th Gaussian shell on center
        *  @param i Shell number
        *  @return A shared pointer to the ShellInfo object for the i'th shell.

        """
        if center is not None:
            si += self.center_to_shell[center]
        if si < 0 or si > self.nshell():
            text = """BasisSet::shell(si = %d), requested a shell out-of-bound.\n   Max shell size: %d\n   Name: %s\n""" % \
                (si, self.nshell(), self.name)
            raise ValidationError("BasisSet::shell: requested shell is out-of-bounds:\n%s" % (text))
        return self.shells[si]

    def nshell_on_center(self, i):
        """Return the number of shells on a given center."""
        return self.center_to_nshell[i]

    def shell_on_center(self, center, shell):
        """Return the overall shell number"""
        return self.center_to_shell[center] + shell

    # <<< Methods for Printing >>>

    def print_by_level(self, out=None, level=2):
        """Print basis set information according to the level of detail in print_level
        @param out The file stream to use for printing. Defaults to outfile.
        @param print_level: defaults to 2
        *  < 1: Nothing
        *    1: Brief summary
        *    2: Summary and contraction details
        *  > 2: Full details

        """
        if level < 1:
            return
        elif level == 1:
            text = self.pyprint(out=None)
        elif level == 2:
            text = self.print_summary(out=None)
        elif level > 2:
            text = self.print_detail(out=None)

        if out is None:
            print(text)
        else:
            with open(out, mode='w') as handle:
                handle.write(text)

    def pyprint(self, out=None):
        """Print the basis set.
        *  @param out The file stream to use for printing. Defaults to outfile.

        """
        text = ''
        text += """  Basis Set: %s\n""" % (self.name)
        text += """    Number of shells: %d\n""" % (self.nshell())
        text += """    Number of basis function: %d\n""" % (self.nbf())
        text += """    Number of Cartesian functions: %d\n""" % (self.nao())
        text += """    Spherical Harmonics?: %s\n""" % ('true' if self.has_puream() else 'false')
        text += """    Max angular momentum: %d\n\n""" % (self.max_am())
        #text += """    Source:\n%s\n""" % (self.source())  # TODO

        if out is None:
            return text
        else:
            with open(outfile, mode='w') as handle:
                handle.write(text)

    def print_summary(self, out=None):
        """Prints a short string summarizing the basis set
        *  @param out The file stream to use for printing. Defaults to outfile.

        """
        text = ''
        text += """  -AO BASIS SET INFORMATION:\n"""
        text += """    Name                   = %s\n""" % (self.name)
        text += """    Total number of shells = %d\n""" % (self.nshell())
        text += """    Number of primitives   = %d\n""" % (self.nprimitive())
        text += """    Number of AO           = %d\n""" % (self.nao())
        text += """    Number of SO           = %d\n""" % (self.nbf())
        text += """    Maximum AM             = %d\n""" % (self.max_am())
        text += """    Spherical Harmonics    = %s\n""" % ('TRUE' if self.puream else 'FALSE')
        text += """\n"""
        text += """  -Contraction Scheme:\n"""
        text += """    Atom   Type   All Primitives // Shells:\n"""
        text += """   ------ ------ --------------------------\n"""

        for A in range(self.molecule.natom()):
            nprims = [0] * (self.PYmax_am + 1)
            nunique = [0] * (self.PYmax_am + 1)
            nshells = [0] * (self.PYmax_am + 1)
            amtypes = [None] * (self.PYmax_am + 1)

            text += """    %4d    """ % (A + 1)
            text += """%2s     """ % (self.molecule.symbol(A))

            first_shell = self.center_to_shell[A]
            n_shell = self.center_to_nshell[A]

            for Q in range(n_shell):
                shell = self.shells[Q + first_shell]
                nshells[shell.am()] += 1
                nunique[shell.am()] += shell.nprimitive()
                nprims[shell.am()] += shell.nprimitive()
                amtypes[shell.am()] = shell.amchar()

            # All Primitives
            for l in range(self.PYmax_am + 1):
                if nprims[l] == 0:
                    continue
                text += """%d%c """ % (nprims[l], amtypes[l])

            # Shells
            text += """// """
            for l in range(self.PYmax_am + 1):
                if nshells[l] == 0:
                    continue
                text += """%d%c """ % (nshells[l], amtypes[l])
            text += """\n"""
        text += """\n"""

        if out is None:
            return text
        else:
            with open(out, mode='w') as handle:
                handle.write(text)

    def print_detail(self, out=None, numbersonly=False):
        """Prints a detailed PSI3-style summary of the basis (per-atom)
        *  @param out The file stream to use for printing. Defaults to outfile.

        """
        text = ''
        if not numbersonly:
            text += self.print_summary(out=None)
            text += """  ==> AO Basis Functions <==\n"""
            text += '\n'
            text += """    [ %s ]\n""" % (self.name)
        text += """    spherical\n""" if self.has_puream() else """    cartesian\n"""
        text += """    ****\n"""

        for uA in range(self.molecule.nunique()):
            A = self.molecule.unique(uA)
            if not numbersonly:
                text += """   %2s %3d\n""" % (self.molecule.symbol(A), A + 1)
            first_shell = self.center_to_shell[A]
            n_shell = self.center_to_nshell[A]

            for Q in range(n_shell):
                text += self.shells[Q + first_shell].pyprint(outfile=None)
            text += """    ****\n"""
        text += """\n"""

        if out is None:
            return text
        else:
            with open(out, mode='w') as handle:
                handle.write(text)


    def export_for_libmints(self, role):
        """From complete BasisSet object, returns array where
        triplets of elements are each unique atom label, the hash
        of the string shells entry in gbs format and the
        shells entry in gbs format for that label. This packaging is
        intended for return to libmints BasisSet::construct_from_pydict for
        instantiation of a libmints BasisSet clone of *self*.

        """
        info = []
        for A in range(self.molecule.natom()):
            label = self.molecule.label(A)
            first_shell = self.center_to_shell[A]
            n_shell = self.center_to_nshell[A]

            atominfo = [label]
            atominfo.append(self.molecule.atoms[A].shell(key=role))
            if self.ecp_coreinfo:
                # If core information is present, this is an ECP so we add the
                # number of electrons this atom's ECP basis accounts for here.
                try:
                   atominfo.append(self.ecp_coreinfo[label])
                except KeyError:
                    raise ValidationError("Problem with number of cores in ECP constuction!")
            for Q in range(n_shell):
                atominfo.append(self.shells[Q + first_shell].aslist())
            info.append(atominfo)
        return info

    def print_detail_gamess(self, out=None, numbersonly=False):
        """Prints a detailed PSI3-style summary of the basis (per-atom)
        *  @param out The file stream to use for printing. Defaults to outfile.

        """
        text = ''
        if not numbersonly:
            text += self.print_summary(out=None)
            text += """  ==> AO Basis Functions <==\n"""
            text += '\n'
            text += """    [ %s ]\n""" % (self.name)
        text += """    spherical\n""" if self.has_puream() else """    cartesian\n"""
        text += """    ****\n"""

        for uA in range(self.molecule.nunique()):
            A = self.molecule.unique(uA)
            if not numbersonly:
                text += """%s\n""" % (qcel.periodictable.to_element(self.molecule.Z(A)))
            first_shell = self.center_to_shell[A]
            n_shell = self.center_to_nshell[A]

            for Q in range(n_shell):
                text += self.shells[Q + first_shell].pyprint_gamess(outfile=None)
            #text += """    ****\n"""
        text += """\n"""

        if out is None:
            return text
        else:
            with open(out, mode='w') as handle:
                handle.write(text)

    def print_detail_cfour(self, out=None):
        """Returns a string in CFOUR-style of the basis (per-atom)
        *  Format from https://web.archive.org/web/20221130013041/http://slater.chemie.uni-mainz.de/cfour/index.php?n=Main.OldFormatOfAnEntryInTheGENBASFile

        """
        text = ''

        for uA in range(self.molecule.nunique()):
            A = self.molecule.unique(uA)
            text += """%s:P4_%d\n""" % (self.molecule.symbol(A), A + 1)
            text += """Psi4 basis %s for element %s atom %d\n\n""" % \
                (self.name.upper(), self.molecule.symbol(A), A + 1)

            first_shell = self.center_to_shell[A]
            n_shell = self.center_to_nshell[A]

            max_am_center = 0
            for Q in range(n_shell):
                max_am_center = self.shells[Q + first_shell].am() if \
                self.shells[Q + first_shell].am() > max_am_center else max_am_center

            shell_per_am = [[] for i in range(max_am_center + 1)]
            for Q in range(n_shell):
                shell_per_am[self.shells[Q + first_shell].am()].append(Q)

            # Write number of shells in the basis set
            text += """%3d\n""" % (max_am_center + 1)

            # Write angular momentum for each shell
            for am in range(max_am_center + 1):
                text += """%5d""" % (am)
            text += '\n'

            # Write number of contracted basis functions for each shell
            for am in range(max_am_center + 1):
                text += """%5d""" % (len(shell_per_am[am]))
            text += '\n'

            exp_per_am = [[] for i in range(max_am_center + 1)]
            coef_per_am = [[] for i in range(max_am_center + 1)]
            for am in range(max_am_center + 1):
                # Collect unique exponents among all functions
                for Q in range(len(shell_per_am[am])):
                    for K in range(self.shells[shell_per_am[am][Q] + first_shell].nprimitive()):
                        if self.shells[shell_per_am[am][Q] + first_shell].exp(K) not in exp_per_am[am]:
                            exp_per_am[am].append(self.shells[shell_per_am[am][Q] + first_shell].exp(K))

                # Collect coefficients for each exp among all functions, zero otherwise
                for Q in range(len(shell_per_am[am])):
                    K = 0
                    for ep in range(len(exp_per_am[am])):
                        if abs(exp_per_am[am][ep] - self.shells[shell_per_am[am][Q] + first_shell].exp(K)) < 1.0e-8:
                            coef_per_am[am].append(self.shells[shell_per_am[am][Q] + first_shell].original_coef(K))
                            if (K + 1) != self.shells[shell_per_am[am][Q] + first_shell].nprimitive():
                                K += 1
                        else:
                            coef_per_am[am].append(0.0)

            # Write number of exponents for each shell
            for am in range(max_am_center + 1):
                text += """%5d""" % (len(exp_per_am[am]))
            text += '\n\n'

            for am in range(max_am_center + 1):
                # Write exponents for each shell
                for ep in range(len(exp_per_am[am])):
                    text += """%14.7f""" % (exp_per_am[am][ep])
                    if ((ep + 1) % 5 == 0) or ((ep + 1) == len(exp_per_am[am])):
                        text += '\n'
                text += '\n'

                # Write contraction coefficients for each shell
                for ep in range(len(exp_per_am[am])):
                    for bf in range(len(shell_per_am[am])):
                        text += """%10.7f """ % (coef_per_am[am][bf * len(exp_per_am[am]) + ep])
                    text += '\n'
                text += '\n'

        if out is None:
            return text
        else:
            with open(out, mode='w') as handle:
                handle.write(text)

    # <<< Misc. Methods >>>

    def refresh(self):
        """Refresh internal basis set data. Useful if someone has pushed
        to shells. Pushing to shells happens in the BasisSetParsers, so
        the parsers will call refresh(). This function is now defunct.

        """
        raise FeatureNotImplemented('BasisSet::refresh')

    @staticmethod
    def make_filename(name):
        """Converts basis set name to a compatible filename.
        * @param basisname Basis name
        * @return Compatible file name.

        """
        # Modify the name of the basis set to generate a filename: STO-3G -> sto-3g
        basisname = name

        # First make it lower case
        basisname = basisname.lower()

        # Replace all '(' with '_'
        basisname = basisname.replace('(', '_')

        # Replace all ')' with '_'
        basisname = basisname.replace(')', '_')

        # Replace all ',' with '_'
        basisname = basisname.replace(',', '_')

        # Replace all '*' with 's'
        basisname = basisname.replace('*', 's')

        # Replace all '+' with 'p'
        basisname = basisname.replace('+', 'p')

        # Add file extension
        basisname += '.gbs'

        return basisname

    @staticmethod
    def decontract(shells):
        """Procedure applied to list to ShellInfo-s *shells* that returns
        another list of shells, one for every AM and exponent pair in the input
        list. Decontracts the shells.

        """
        # vector of uncontracted shells to return
        shell_list = []

        # map of AM to a vector of exponents for duplicate basis functions check
        exp_map = collections.defaultdict(list)

        for shell in shells:
            am = shell.am()
            pure = shell.is_pure()
            nc = shell.ncenter()
            center = shell.center
            start = shell.start

            for prim in range(shell.nprimitive()):
                exp = shell.exp(prim)
                unique = True
                for _exp in exp_map[am]:
                    if abs(exp - _exp) < 1.0e-6:
                        unique = False
                if unique:
                    us = ShellInfo(am, [1.0], [exp],
                        'Pure' if pure else 'Cartesian',
                        nc, center, start, 'Unnormalized')
                    shell_list.append(us)
                    exp_map[am].append(exp)

        return shell_list

    # <<< Methods not Implemented >>>

    def zero_so_basis_set(cls, factory):
        """ **NYI** Returns an empty SO basis set object.
        *  Returns an SOBasis object that actually has a single s-function
        *  at the origin with an exponent of 0.0 and contraction of 1.0.
        *  @return A new empty SOBasis object.

        """
        raise FeatureNotImplemented('BasisSet::zero_so_basis_set')  # FINAL

    @staticmethod
    def test_basis_set(max_am):
        """Returns a shell-labeled test basis set object
        * @param max_am maximum angular momentum to build
        * @return pair containing shell labels and four-center
        * test basis for use in benchmarking
        * See libmints/benchmark.cc for details
        The libmints version seems not to have been updated along with the classes.

        """
        raise FeatureNotImplemented('BasisSet::test_basis_set')

    def get_ao_sorted_shell(self, i):
        """Returns the value of the sorted shell list. Defunct"""
        raise FeatureNotImplemented('BasisSet::get_ao_sorted_shell')

    def get_ao_sorted_list(self):
        """Returns the vector of sorted shell list. Defunct"""
        raise FeatureNotImplemented('BasisSet::get_ao_sorted_list')

    def compute_phi(self, phi_ao, x, y, z):
        """Returns the values of the basis functions at a point"""

        phi_ao = [0.0] * self.nao()
        ao = 0
        for ns in range(self.nshell()):
            shell = self.shells[ns]
            am = shell.am()
            nprim = shell.nprimitive()
            a = shell.exps()
            c = shell.coefs()

            xyz = shell.center()
            dx = x - xyz[0]
            dy = y - xyz[1]
            dz = z - xyz[2]
            rr = dx * dx + dy * dy + dz * dz

            cexpr = 0
            for np in range(nprim):
                cexpr += c[np] * math.exp(-a[np] * rr)

            for l in range(INT_NCART(am)):
                components = exp_ao[am][l]
                phi_ao[ao + l] += pow(dx, components[0]) * \
                                pow(dy, components[1]) * \
                                pow(dz, components[2]) * \
                                cexpr

            ao += INT_NCART(am)

    def concatenate(self, b):
        """Concatenates two basis sets together into a new basis without
        reordering anything. Unless you know what you're doing, you should
        use the '+' operator instead of this method. Appears defunct.

        """
        raise FeatureNotImplemented('BasisSet::concatenate')

    def add(self, b):
        """Adds this plus another basis set and returns the result.
        Equivalent to the '+' operator. Appears defunct.

        """
        raise FeatureNotImplemented('BasisSet::add')

    @staticmethod
    def shell_sorter_ncenter(d1, d2):
        return d1.ncenter() < d2.ncenter()

    @staticmethod
    def shell_sorter_am(d1, d2):
        return d1.am() < d2.am()


def _basis_file_warner_and_aliaser(filename):
    aliased_in_1p4 = {
        "def2-qzvp-jkfit": "def2-universal-jkfit",
        "def2-qzvpp-jkfit": "def2-universal-jkfit",
        "def2-sv_p_-jkfit": "def2-universal-jkfit",
        "def2-svp-jkfit": "def2-universal-jkfit",
        "def2-tzvp-jkfit": "def2-universal-jkfit",
        "def2-tzvpp-jkfit": "def2-universal-jkfit",

        "def2-qzvp-jfit": "def2-universal-jfit",
        "def2-qzvpp-jfit": "def2-universal-jfit",
        "def2-sv_p_-jfit": "def2-universal-jfit",
        "def2-svp-jfit": "def2-universal-jfit",
        "def2-tzvp-jfit": "def2-universal-jfit",
        "def2-tzvpp-jfit": "def2-universal-jfit",
    }
    for k, v in aliased_in_1p4.items():
        if filename.endswith(k + ".gbs"):
            warnings.warn(
                f"Using basis set `{k}` instead of its generic name `{v}` is deprecated, and as soon as 1.5 it will stop working\n",
                category=FutureWarning,
                stacklevel=2)
            return filename.replace(k, v)
    else:
        return filename
