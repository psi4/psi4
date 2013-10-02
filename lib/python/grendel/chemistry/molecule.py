from __future__ import print_function
from collections import Iterable
from copy import deepcopy, copy
from functools import total_ordering
from numbers import Number
import os
import re
import itertools
import math
from math import sin, cos, pi

import numpy as np

from grendel import type_checking_enabled, sanity_checking_enabled, caching_enabled


from grendel.gmath.geometry import angle_between_vectors, rotate_point_about_axis
from grendel.gmath.vector import magnitude, Vector, LightVector, cross
from grendel.gmath.matrix import Matrix

from grendel.symmetry.identity_operation import IdentityOperation
from grendel.symmetry.improper_rotation import ImproperRotation
from grendel.symmetry.inversion import Inversion
from grendel.symmetry.point_group import PointGroup
from grendel.symmetry.reflection import Reflection
from grendel.symmetry.rotation import Rotation
from grendel.symmetry.symmetry_operation import SymmetryOperation

from grendel.util import cirpy
from grendel.util.aliasing import function_alias
from grendel.util.decorators import with_flexible_arguments, cached_method, CachedProperty, typecheck, typechecked, IterableOf
from grendel.util.exceptions import raises_error, min_arg_count
from grendel.util.iteration import grouper
from grendel.util.misc import distinct, full_path
from grendel.util.parsing import re_float_or_int_or_scientific
from grendel.util.overloading import overloaded, OverloadedFunctionCallError
from grendel.util.units import Angstroms, Degrees, Radians, EnergyUnit, AngularUnit, DistanceUnit, ValueWithUnits, Kilograms, Meters, Joules, MassUnit
from grendel.util.units.physical_constants import ReducedPlanckConstant
from grendel.util.units.unit import isunit
from grendel.util.units.value_with_units import  hasunits, strip_units
from grendel.util.strings import classname, indented
from grendel.util.web_getter import MoleculeNotFoundError, input_identifiers

__all__ = [
    "Molecule",
    "InvalidZMatrixException",
    "InvalidXYZFormatError"
]

class InvalidZMatrixException(Exception):
    """ Raised if a line of a z-matrix has the wrong number of entries.
    """
    def __init__(self, line):
        super(InvalidZMatrixException, self).__init__("Invalid z-matrix entry:  " + repr(line))

class InvalidXYZFormatError(Exception):
    """ Raised if an invalid string in the XYZ format is passed into the Molecule class constructor
    """

#TODO Units for atom positions
@total_ordering
class Molecule(object):
    """
    Encapsulates all of the functionality and attributes of a Molecule itself.

    **Signatures :**
        * ``Molecule(xyz_string)``
        * ``Molecule(atoms)``
        * ``Molecule(atom_names, cart_mat)``


    :Parameters:

    xyz_string : str
        a string in the format of a standard .xyz file
    atoms : list of Atom
        a list of `Atom` objects
    atom_names : list of str
        a list of atomic symbols corresponding to the rows of the `cart_mat` parameter
    cart_mat : Matrix
        an Nx3 `Matrix` of positions


    :Other Parameters:

    description : str
        optional keyword argument that works with all forms.  See the `description` attribute


    :Attributes:

    atoms : list of Atom
    internal_representations : list of InternalRepresentation
    normal_representation : NormalRepresentation
    cartesian_representation : `CartesianRepresentation`
    description : str


    :Examples:


    *Constructor*

    >>> Molecule(\"""
    ...     5
    ...
    ...     C        0.000000        0.000000        0.000000
    ...     H        0.000000        0.000000        1.089000
    ...     H        1.026719        0.000000       -0.363000
    ...     H       -0.513360       -0.889165       -0.363000
    ...     H       -0.513360        0.889165       -0.363000
    ... \""")
    Molecule([
        Atom('C', [  0.00000000,  0.00000000,  0.00000000 ] ),
        Atom('H', [  0.00000000,  0.00000000,  1.08900000 ] ),
        Atom('H', [  1.02671900,  0.00000000, -0.36300000 ] ),
        Atom('H', [ -0.51336000, -0.88916500, -0.36300000 ] ),
        Atom('H', [ -0.51336000,  0.88916500, -0.36300000 ] )
    ])
    >>> Molecule([
    ...     Atom('H', [0,0,0]),
    ...     Atom('H', [0,0,0.75])
    ... ])
    Molecule([
        Atom('H', [  0.00000000,  0.00000000,  0.00000000 ] ),
        Atom('H', [  0.00000000,  0.00000000,  0.75000000 ] )
    ])
    >>> Molecule(['H','C'], Matrix([[0,0,0],[1,0.5,0]]))
    Molecule([
        Atom('H', [  0.00000000,  0.00000000,  0.00000000 ] ),
        Atom('C', [  1.00000000,  0.50000000,  0.00000000 ] )
    ])
    >>> mol = Molecule(\"""
    ...     5
    ...
    ...     C        0.000000        0.000000        0.000000
    ...     H        0.000000        0.000000        1.089000
    ...     H        1.026719        0.000000       -0.363000
    ...     H       -0.513360       -0.889165       -0.363000
    ...     H       -0.513360        0.889165       -0.363000
    ... \""")
    >>> mol
    Molecule([
        Atom('C', [  0.00000000,  0.00000000,  0.00000000 ] ),
        Atom('H', [  0.00000000,  0.00000000,  1.08900000 ] ),
        Atom('H', [  1.02671900,  0.00000000, -0.36300000 ] ),
        Atom('H', [ -0.51336000, -0.88916500, -0.36300000 ] ),
        Atom('H', [ -0.51336000,  0.88916500, -0.36300000 ] )
    ])

    *Iteration*

    >>> from __future__ import print_function
    >>> for atom in mol:
    ...     print(repr(atom))
    ...
    Atom('C', [  0.00000000,  0.00000000,  0.00000000 ] )
    Atom('H', [  0.00000000,  0.00000000,  1.08900000 ] )
    Atom('H', [  1.02671900,  0.00000000, -0.36300000 ] )
    Atom('H', [ -0.51336000, -0.88916500, -0.36300000 ] )
    Atom('H', [ -0.51336000,  0.88916500, -0.36300000 ] )

    """

    ####################
    # Class Attributes #
    ####################

    linear_cutoff = 5 * Degrees
    default_multiplicity = 1
    default_charge = 0
    global_result_getters = []

    ##############
    # Attributes #
    ##############

    internal_representations = None
    """ A list of all the InternalRepresentation objects associated with self.
    Usually there is just one (or 0), and it may be more convenient to access it using the `internal_representation()` method.
    """

    normal_representation = None
    """ The NormalRepresentation associated with the molecule, or None if one has not been computed yet. """

    atoms = None
    """ The atoms that `self` is composed of. """

    description = None
    """ A string describing the molecule (optional) """

    computations = None
    """ The list of Computation objects that Molecule is a part of. """

    cartesian_units = None
    """ The units to use for the `CartesianRepresentation` objects automatically created by the Molecule instance.
    Defaults to Angstroms.
    """

    result_getters = None
    """ `list` of `ResultGetter` objects associated with `self`.
    """

    displacement = None
    """ If the Molecule was generated as a displacement, the Displacement object is stored here.
    """

    ######################
    # Private Attributes #
    ######################

    _cirmol = None             # Used if the molecule is initialized from the NIH CIR service.  Can be used to retrieve may properties.
    _point_group = None
    _cartesian_representation = None
    # Variables for cached methods and properties
    _center_of_mass = None
    _pmi = None
    _principal_axes = None
    _A_e = _B_e = _C_e = None
    _is_linear = None
    _properties = None


    ##################
    # Initialization #
    ##################

    @overloaded
    @with_flexible_arguments(
        optional=[
            ('units', 'cartesian_units'),
            ('description', 'comment'),
            ('charge',),
            ('multiplicity',),
        ],
        what_to_call_it='Molecule constructor'
    )
    def __init__(self, *args, **kwargs):
        """
        """
        raise OverloadedFunctionCallError

    # TODO Raise a ValueError for a incorrectly formatted string
    @__init__.overload_with(
        xyz_string=basestring)
    def __init__(self, xyz_string, **kwargs):
        orig_string = xyz_string
        # First see if the header is there
        m = re.match(
            r'\s*'                       # zero or more leading blank lines
            r'[ \t]*(\d+)[ \t]*\n'       # A number on a line by itself
            r'(.*)\n'                    # A description line (. matches anything _but_ the newline by default)
            r'([ \t]*.*\d(.*\n?)+)\Z',   # The rest, starting with a line containing at least one number
                                         #   if there are problems with this, they will show up later.
            xyz_string,
            re.MULTILINE
        )
        if m:
            # the headed version is being used; get the number of atoms and the description
            num_atoms = int(m.group(1))
            if 'description' in kwargs:
                kwargs['description'] += m.group(2)
            else:
                kwargs['description'] = m.group(2)
            xyz_string = m.group(3)
        else:
            # we don't care how many atoms there are
            num_atoms = None

        # Now parse the main part
        m = re.match(
            r'\s*'                                                     # zero or more leading blank lines
            r'((?:[ \t]*\w{{1,3}}(?:[ \t]+{number}){{3}}[ \t]*\n?)+)'  # an entry in an xyz file.  (The doubled {{ and }}
                                                                       #   are escaped curly brackets in the format
                                                                       #   string syntax, so that the format method on the line
                                                                       #   below this doesn't try to replace the {} expressions)
            r'\s*\Z'                                                   # zero or more trailing blank lines
            r''.format(number=re_float_or_int_or_scientific),          # insert the canned float/int/scientific regex
            xyz_string,
            re.MULTILINE
        )
        if m:
            atoms = []
            natoms=0
            for line in m.group(1).splitlines():
                symb, x, y, z = line.strip().split()
                atom = Atom(symb, float(x), float(y), float(z))
                atom.parent_molecule = self
                atoms.append(atom)
                natoms+=1
            if not num_atoms is None and natoms != num_atoms:
                raise InvalidXYZFormatError("expecting {} atoms but found {} atoms"
                                            " in XYZ string:\n{}".format(
                    num_atoms,
                    natoms,
                    indented(orig_string)
                ))
            Molecule.__init__(self, atoms, copy_atoms=False, **kwargs)
        else:
            raise InvalidXYZFormatError("invalid XYZ string:\n{}".format(indented(orig_string)))

    @__init__.overload_with(
        atom_names=IterableOf(basestring),
        cart_mat=Matrix,
        copy_atoms=bool)
    def __init__(self, atom_names, cart_mat, copy_atoms=False, **kwargs):
        if sanity_checking_enabled:
            if len(atom_names) != cart_mat.shape[0]:
                raise ValueError("dimension mismatch in Molecule constructor ({} != {})".format(
                    len(atom_names),
                    cart_mat.shape[0]
                ))
            if cart_mat.shape[1] != 3:
                raise ValueError("non-cartesian matrix passed to Molecule constructor")
        atoms = []
        for sym, vect in zip(atom_names, cart_mat.rows):
            atoms.append(Atom(sym, deepcopy(vect) if copy_atoms else vect))
        Molecule.__init__(self, atoms, copy_atoms=False, **kwargs)

    @__init__.overload_with(
        molecule='Molecule',
        charge=(int, None),
        multiplicity=(int, None),
        description=(basestring, None))
    def __init__(self,
            molecule,
            charge=None,
            multiplicity=None,
            description=None):
        charge = charge or molecule.charge
        multiplicity = multiplicity or molecule.multiplicity
        if description is None and molecule.description is not None:
            description = '(Copy of) ' + molecule.description
        Molecule.__init__(self,
            atoms=molecule.atoms,
            copy_atoms=True,
            description=description,
            charge=charge,
            multiplicity=multiplicity,
            units=molecule.cartesian_units
        )

    # TODO @usability as is, this won't print a very helpful error message if an optional argument has the wrong type and a different initializer is called (need the common initializer arguments improvement to overloaded)
    # Principal initializer (all other __init__ methods call this one)
    @__init__.overload_with(
        atoms=IterableOf('Atom'),
        copy_atoms=bool,
        description=(basestring, None),
        units=isunit,
        charge=int,
        multiplicity=int)
    def __init__(
            self,
            atoms,
            copy_atoms=False,
            description=None,
            units=DistanceUnit.default,
            charge=0,
            multiplicity=1):
        #----------------------------------------#
        self.atoms = []
        self.internal_representations = []
        self.computations = []
        self.result_getters = []
        self.result_getters.extend(Molecule.global_result_getters)
        self.description = description
        self.cartesian_units = units
        self.charge = charge
        self.multiplicity = multiplicity
        self._properties = []
        #----------------------------------------#
        if False in map(isinstance, atoms, [Atom]*len(atoms)):
            raise TypeError("Expecting a list of Atom objects.")
        if copy_atoms:
            self.atoms = [deepcopy(atom) for atom in atoms]
        else:
            self.atoms = [atom for atom in atoms]
        for atom in self.atoms:
            atom.parent_molecule = self
            atom._cartesian_units = self.cartesian_units
        # Add the current default cartesian representation
        self.update_cartesian_representation()



    ##############
    # Properties #
    ##############

    @property
    def natoms(self):
        return len(self.atoms)

    @property
    def internal_representation(self):
        """ The first internal representations (of type InternalRepresentation), or None if it does not exist yet.
        """
        if len(self.internal_representations) == 0:
            return None
        else:
            return self.internal_representations[0]

    @internal_representation.setter
    def internal_representation(self, other_rep):
        self.internal_representations[0] = other_rep
        other_rep.molecule = self

    @property
    def mass(self):
        ret_val = 0.0
        for atom in self:
            ret_val += atom.mass
        return ret_val

    @property
    def xyz(self):
        """ The molecule's position as a natoms x 3 `Matrix`.
        The ordering of rows is (as expected) the same as the ordering of the `Molecule.atoms` list attribute.  Aliased
        to `xyz_mat` and `position` (the latter to be "consistant-ish" with the naming in `Atom`.
        """
        return Matrix([atom.pos for atom in self])
    xyz_mat = position = xyz

    @xyz.setter
    def xyz(self, new_xyz):
        if type_checking_enabled:
            if not isinstance(new_xyz, Matrix):
                raise TypeError
        if sanity_checking_enabled:
            if not new_xyz.shape == (self.natoms, 3):
                raise ValueError("Invalid dimensions for matrix to set as molecules xyz. Expected ({}, 3), got ({})".format(
                    str(self.natoms),
                    ','.join(str(s) for s in new_xyz.shape)
                ))
        for i, row in enumerate(new_xyz.row_iter):
            self.atoms[i].position = row
            self.atoms[i].position.units = self.cartesian_units
        self.update_cartesian_representation()


    @CachedProperty
    def A_e(self):
        """ The (non-zero-point-corrected) rotational constant A_e.
        .. note::
           This property is cached, and the cached value gets flushed in
           `update_cartesian_representation()`.  If you change an atom's position
           (or mass) and forget to call `update_cartesian_representation()`, you may
           get some funny results for this property. You can detect whether caching
           is causing your problems by setting the environment variable
           `GRENDEL_NO_CACHE` to 1 and rerunning your tests.  If tests that were failing
           subsequently succeed, you probably forgot to call `update_cartesian_representation()`
           somewhere, or you were assuming that it was automatically called somewhere when
           in fact it was not getting called.
        """
        I = self.pmi()[0]
        units = element_data.MASS_UNIT * self.cartesian_units**2
        conv_factor = units.to(Kilograms * Meters**2)
        I *= conv_factor
        ret_val = ReducedPlanckConstant**2 / (2 * I)
        ret_val *= Joules.to(EnergyUnit.default)
        return ValueWithUnits(ret_val, EnergyUnit.default)

    @CachedProperty
    def B_e(self):
        """ The (non-zero-point-corrected) rotational constant B_e.

        .. note::
           This property is cached.  See discussion of caching in `A_e`

        """
        I = self.pmi()[1]
        units = MassUnit.default * self.cartesian_units**2
        conv_factor = units.to(Kilograms * Meters**2)
        I *= conv_factor
        ret_val = ReducedPlanckConstant**2 / (2 * I)
        ret_val *= Joules.to(EnergyUnit.default)
        return ValueWithUnits(ret_val, EnergyUnit.default)

    @CachedProperty
    def C_e(self):
        """ The (non-zero-point-corrected) rotational constant C_e.

        .. note::
           This property is cached.  See discussion of caching in `A_e`

        """
        I = self.pmi()[2]
        units = MassUnit.default * self.cartesian_units**2
        conv_factor = units.to(Kilograms * Meters**2)
        I *= conv_factor
        ret_val = ReducedPlanckConstant**2 / (2 * I)
        ret_val *= Joules.to(EnergyUnit.default)
        return ValueWithUnits(ret_val, EnergyUnit.default)

    @property
    def a(self):
        """ The first principal axis of rotation, a.  The positive phase is always chosen.

        .. note::
           This property depends directly on `Molecule.principal_axes()`, which is cached.  See discussion of
           caching in `A_e`

        """
        rv = self.principal_axes()[:,0]
        # Always choose the positive phase.
        #if sum(rv) < 0:
        #    return -rv
        #else:
        #    return rv
        return rv

    @property
    def b(self):
        """ The second principal axis of rotation, b.  The positive phase is always chosen.


        .. note::
           This property depends directly on `Molecule.principal_axes()`, which is cached.  See discussion of
           caching in `A_e`

        """
        rv = self.principal_axes()[:,1]
        # Always choose the positive phase.
        #if sum(rv) < 0:
        #    return -rv
        #else:
        #    return rv
        return rv

    @property
    def c(self):
        """ The third principal axis of rotation, c.  The positive phase is always chosen.

        .. note::
           This property depends directly on `Molecule.principal_axes()`, which is cached.  See discussion of
           caching in `A_e`

        """
        rv = self.principal_axes()[:,2]
        # Always choose the positive phase.
        #if sum(rv) < 0:
        #    return -rv
        #else:
        #    return rv
        return rv

    @property
    def point_group(self):
        """
        """
        if not self._point_group is None:
            return self._point_group
        self._compute_point_group()
        return self._point_group

    @CachedProperty
    def nelectrons(self):
        return self.Z - self.charge

    @CachedProperty
    def Z(self):
        """
        @return: The total nuclear charge of the molecule.  Molecule.charge is not taken into account.
        """
        return sum(atom.element.atomic_number for atom in self)

    @property
    def cartesian_representation(self):
        """ The current default cartesian representation associated with the molecule.
        If the molecule has multiple cartesian representations, the first in the list `self.cartesian_representations`
        is returned.  The setter for this property pushes the cartesian representation onto the front of the list.

        .. note::
            The setter does not check for uniqueness (unless it's exactly the same object as another `CartesianRepresentation`).
            Thus, you should check to make sure the `CartesianRepresentation` you are assigning to the molecule is not
            already part of the molecule's list of representations (to the degree of accuracy required for your particular
            application) before assigning.  Repeatedly failing to do this check could cause runaway memory usage.

        """
        # The molecule must have at least one cartesian_representation, the one created when the molecule was created.
        return self._cartesian_representation

    @cartesian_representation.setter
    @typechecked(new_rep='CartesianRepresentation')
    def cartesian_representation(self, new_rep):
        self._cartesian_representation = new_rep
        # At this point, we've released the old representation.  However, if any properties depend on a representation,
        #  those properties will still have a reference to the representation they depend on.  The molecule can
        #  always be reverted back to the original (for use with the property in question) if needed simply by setting the
        #  cartesian_reperesentation property as we have done here!
        self._clear_cartesian_dependent_cached_variables()
        # Now adapt atoms to the new cartesian representation...
        for i, coords in enumerate(grouper(3, self._cartesian_representation.value)):
            coords = [strip_units(c, convert_to=self.cartesian_units, assume_units=new_rep.units)
                            for c in coords]
            self.atoms[i].position = Vector(coords, units=self.cartesian_units)

    @property
    def ninternals(self):
        if self.is_linear():
            return 3 * self.natoms - 5
        else:
            return 3 * self.natoms - 6

    @CachedProperty
    def inverse_mass_matrix(self):
        return np.diag(sum(([1.0 / atom.mass]*3 for atom in self), [])).view(Matrix)

    @CachedProperty
    def inverse_sqrt_mass_matrix(self):
        return np.diag(sum(([1.0 / math.sqrt(atom.mass)]*3 for atom in self), [])).view(Matrix)

    ###################
    # Special Methods #
    ###################

    def __add__(self, other):
        if isinstance(other, Molecule):
            new_atoms = list(self.atoms) + list(other.atoms)
            if self.description is None or other.description is None:
                new_desc = None
            else:
                new_desc = self.description + " + " + other.description
            if self.cartesian_units != other.cartesian_units:
                raise NotImplementedError()
            if self.multiplicity != 1 or other.multiplicity != 1:
                raise NotImplementedError()
            new_charge = self.charge + other.charge
            return Molecule(
                atoms=new_atoms,
                description=new_desc,
                units=self.cartesian_units,
                charge=new_charge
            )
        else:
            return NotImplemented

    def __contains__(self, item):
        for atom in self:
            if atom == item:
                return True
        return False

    def __getitem__(self, item):
        return self.atoms[item]

    def __iter__(self):
        return itertools.chain(self.atoms)

    def __eq__(self, other):
        my_atoms=sorted(self.atoms)
        other_atoms=sorted(other.atoms)
        return my_atoms == other_atoms

    def __lt__(self, other):
        my_atoms = sorted(self.atoms)
        other_atoms = sorted(other.atoms)
        return my_atoms < other_atoms

    def __copy__(self):
        copied = self.copy_with_atoms(self.atoms)
        # copied Molecules have their own cartesian_representations, internal_representations,
        # and normal_representation, so don't copy these by default.
        copied.description = self.description
        # computations should not be copied over, since they do not refer to the new instance
        copied.multiplicity = self.multiplicity
        copied.charge = self.charge
        # The methods of getting results for the copied molecule should be the same as for the
        #   original, even though, say, the Computation instances themselves are not.
        copied.result_getters = copy(self.result_getters)
        return copied

    def __deepcopy__(self, memo):
        copied = self.copy_with_atoms(self.atoms, deep_copy=True, deep_copy_memo=memo)
        # copy over cached variables, since the atom's positions will not have changed (yet, at least)
        if caching_enabled:
            copied._pmi = self._pmi
            copied._principal_axes = self._principal_axes
            copied._center_of_mass = self._center_of_mass
            copied._A_e, copied._B_e, copied._C_e = self._A_e, self._B_e, self._C_e
        # No need to deepcopy, ResultGetter instances are not molecule-specific
        copied.result_getters = copy(self.result_getters)
        return copied

    #------------------------#
    # Output Representations #
    #------------------------#

    def __repr__(self):
        ret_val = "Molecule([\n    " + ',\n    '.join(repr(atom) for atom in self.atoms)
        other_stuff_added = False
        # Put any other information in here as kwargs...(such as representations)
        # Then...
        if len(self.atoms) == 0 and not other_stuff_added:
            ret_val = ret_val[:-5] + '])'
        else:
            ret_val += "\n])"
        return ret_val

    def __short_str__(self):
        if self.description is not None:
            return self.description
        else:
            # TODO full symbol here for non-principal isotopes
            return ''.join(a.symbol for a in self.atoms)


    #################
    # Class Methods #
    #################

    #TODO make this take a string argument
    #TODO better error checking, better type checking
    #TODO dummy atoms?
    #TODO set units based on global units or kwarg
    #TODO allow zero_based indexes (one_based still as default)
    #TODO helpful error messages
    #TODO check the handedness of the coordinate system
    @classmethod
    @with_flexible_arguments(
        optional=[
            ('one_based',),
            ('dist_units', 'distance_units', 'length_units'),
            ('ang_units', 'angular_units', 'angle_units'),
            ('create_representation',),
        ]
    )
    def from_z_matrix(cls, *args, **kwargs):
        '''
        TODO Document this more

        **Signatures :**
            * Molecule.from_z_matrix(atom1, atom2, ...)
            * Molecule.from_z_matrix(atoms, create_representation = False)


        :Examples:


        >>> Molecule.from_z_matrix(\"""
        ... O
        ... H1 O 1.0
        ... H2 O 1.0 H1 90.0
        ... """)
        Molecule([
            Atom('O', [  0.00000000,  0.00000000,  0.00000000 ] ),
            Atom('H', [  0.00000000,  0.00000000,  1.00000000 ] ),
            Atom('H', [  0.00000000, -1.00000000, -0.00000000 ] )
        ])
        >>> Molecule.from_z_matrix(
        ...    [
        ...       ['H'],
        ...       ['O', 1, 0.9],
        ...       ['O', 2, 1.4, 1, 105.0],
        ...       ['H', 3, 0.9, 2, 105.0, 1, 120.0]
        ...   ]
        ... )
        Molecule([
            Atom('H', [  0.00000000,  0.00000000,  0.00000000 ] ),
            Atom('O', [  0.00000000,  0.00000000,  0.90000000 ] ),
            Atom('O', [  0.00000000,  1.35229616,  1.26234666 ] ),
            Atom('H', [ -0.75286467,  1.46479616,  1.74249095 ] )
        ])
        >>> Molecule.from_z_matrix(
        ...       ['H'],
        ...       ['O', 1, 0.9]
        ... )
        Molecule([
            Atom('H', [  0.00000000,  0.00000000,  0.00000000 ] ),
            Atom('O', [  0.00000000,  0.00000000,  0.90000000 ] )
        ])

        '''
        # Pop off optional keyword arguments
        create_representation = kwargs.pop("create_representation", False)
        one_based = kwargs.pop("one_based", True)
        offset = 1 if one_based else 0
        dist_units = kwargs.pop('dist_units', DistanceUnit.default)
        ang_units = kwargs.pop('ang_units', Degrees)
        ang_conv = ang_units.to(Radians) * Radians
        #--------------------------------------------------------------------------------#
        # Make the *args into a uniform iterable
        itlist = []
        if len(args) == 1:
            # Remember strings are iterable, so we have to do that first
            if isinstance(args[0], basestring):
                for line in args[0].splitlines():
                    # ignore blank lines
                    if re.match(r'^\s*$', line): continue
                    itlist.append(line)
            elif isinstance(args[0], Iterable):
                for item in args[0]:
                    itlist.append(item)
            else:
                raise TypeError
        elif len(args) == 0:
            raise min_arg_count(1, 0)
        else:
            itlist = list(args)
        #----------------------------------------#
        # now parse the uniform iterable into two arrays, zatoms and zlines (one is
        #   the original lines, the other is the lines split into lists)
        zatoms = []
        zlines = []
        for zatom in itlist:
            if isinstance(zatom, basestring):
                zlines.append(zatom)
                zatoms.append(zatom.split())
            elif isinstance(zatom, Iterable):
                zatoms.append(zatom)
                zlines.append(repr(zatom))
            else:
                raise TypeError
        #--------------------------------------------------------------------------------#
        # Helper function to uniformly get the atom from a string
        def _labeled_atom(atomstr, cartvect):
            m = re.match(r'^([A-Z][a-z]?)(\d*)$', str(atomstr))
            if m:
                ret_val = Atom(m.group(1), cartvect, units=dist_units)
                if len(m.group(2)) > 0:
                    ret_val.zmat_label = m.group(2)
            else:
                raise ValueError("Invalid atom label in z-matrix: " + str(atomstr))
            return ret_val
        #--------------------------------------------------------------------------------#
        # Parse the first line
        if len(zatoms[0]) != 1:
            raise InvalidZMatrixException(zlines[0])
        atoms = [_labeled_atom(zatoms[0][0], Vector(0.0, 0.0, 0.0, units=dist_units))]
        # If there's only one atom, we're done
        if len(zatoms) == 1:
            if create_representation:
                raise ValueError("Monoatomics are only representable as cartesians.  Cannot create an internal representation.")
            return Molecule(atoms)
        #--------------------------------------------------------------------------------#
        # Define another helper that handles all of the possible formats uniformly
        def _get_zatom(zatom):
            if not raises_error(int, zatom, errors=(ValueError, TypeError)):
                idx = int(zatom)
                if idx - offset > len(atoms) - 1:
                    raise IndexError("z-matrix index '{}' is out of range.".format(idx))
                elif idx - offset < 0:
                    raise IndexError("z-matrix index '{}' is out of range."
                                     " to use zero-based indices, give"
                                     " the keyword option one_based=False"
                                     " to Molecule.from_z_matrix".format(idx))
                return atoms[idx - offset]
            else:
                # Convert it to a string, and try and use it
                m = re.match(r'^([A-Z][a-z]?)(\d*)$', zatom)
                if m:
                    el = m.group(1)
                    z_id = m.group(2) if len(m.group(2)) > 0 else None
                    # I know, this is a mess.  Sorry
                    #  All it does is prevent z_id from storing '' as its id string, and checks for this (which needs
                    #  to use the 'is' syntax to avoid 'some_var == None')
                    possiblities = [atom for atom in atoms
                                    if ((atom.zmat_label is None and z_id is None)
                                           or (z_id is not None and atom.zmat_label == z_id))
                                        and atom.symbol == el
                                   ]
                    if len(possiblities) == 0:
                        raise ValueError("Unknown atom identifier in z-matrix: " + zatom)
                    elif len(possiblities) > 1:
                        raise ValueError("Ambiguous atom identifier in z-matrix: " + zatom)
                    return possiblities[0]
                else:
                    raise ValueError("Invalid atom identifier in z-matrix: " + zatom)
        #--------------------------------------------------------------------------------#
        # Parse the second atom
        if len(zatoms[1]) != 3:
            raise InvalidZMatrixException(zlines[1])
        # If create_representation is given, initialize the list of coordinates
        coords_to_make = None
        if create_representation:
            coords_to_make = []
        # and utilize the helper functions to get the specifics of the atom
        dist_atom = _get_zatom(zatoms[1][1])
        dist = float(zatoms[1][2])
        # place the atom...
        pos_vect = dist_atom.position + Vector([0.0, 0.0, dist], units=dist_units)
        atoms.append(_labeled_atom(zatoms[1][0], pos_vect))
        # create a BondLength coordinate, if necessary
        if create_representation:
            coords_to_make.append((BondLength, dist_atom, atoms[-1]))
        # if there are only two atoms, we're done
        if len(zatoms) == 2:
            ret_val = Molecule(atoms)
            if create_representation:
                # Now that the atoms are not orphaned, make the coordinates
                cargs = coords_to_make[0]
                coords = [cargs[0](cargs[1:], units=dist_units)]
                # The __ indicates that we ignore the created object, since the InternalRepresentation
                #   constructor itself sets the internal_representation of the Molecule
                __ = InternalRepresentation(
                         ret_val,
                         coords,
                         units={DistanceUnit:dist_units, AngularUnit:ang_units})
            return ret_val
        #--------------------------------------------------------------------------------#
        # Parse the third atom
        if len(zatoms[2]) != 5:
            raise InvalidZMatrixException(zlines[2])
        # get the data about the atom using the helper functions...
        dist_atom, dist = _get_zatom(zatoms[2][1]), float(zatoms[2][2])
        ang_atom, ang = _get_zatom(zatoms[2][3]), float(zatoms[2][4]) * ang_conv
        if not distinct(dist_atom, ang_atom):
            raise InvalidZMatrixException("atom indices in z-matrix line '{}' must" \
                                          " all be distinct".format(
                zlines[2].strip()
            ))
        # place the atom...
        n = (dist_atom.position - ang_atom.position).normalize() * dist
        pos_vect = dist_atom.position + Vector(n.x, math.sin(math.pi - ang) * n.z, math.cos(math.pi - ang) * n.z, units=dist_units)
        atoms.append(_labeled_atom(zatoms[2][0], pos_vect))
        # create BondLength and BondAngle coordinates, if necessary
        if create_representation:
            coords_to_make.append((BondLength, dist_atom, atoms[-1]))
            coords_to_make.append((BondAngle, ang_atom, dist_atom, atoms[-1]))
        # if there are only three atoms, we're done
        if len(zatoms) == 3:
            ret_val = Molecule(atoms)
            if create_representation:
                # Now that the atoms are not orphaned, make the coordinates
                coords = []
                for cargs in coords_to_make:
                    if len(cargs) == 3:
                        coords.append(cargs[0](cargs[1:], units=dist_units))
                    else:
                        coords.append(cargs[0](cargs[1:], units=ang_units))
                # The __ indicates that we ignore the created object, since the InternalRepresentation
                #   constructor itself sets the internal_representation of the Molecule
                __ = InternalRepresentation(
                         ret_val,
                         coords,
                         units={DistanceUnit:dist_units, AngularUnit:ang_units})
            return ret_val
        #--------------------------------------------------------------------------------#
        # Finally, parse any remaining lines, all of which should have 7 parts
        #TODO check for linear part in torsion angle raise an appropriate error
        for idx, zatom in enumerate(zatoms[3:]):
            # Allow the "alternate z-matrix format" from Gaussian
            if len(zatom) not in (7, 8):
                raise InvalidZMatrixException(zlines[idx])
            if len(zatom) == 8 and int(zatom[7]) != 0:
                alternate_format = True
            else:
                alternate_format = False
            #----------------------------------------#
            # get the data about the atom using the helper functions...
            # common to both normal and alternate format
            dist_atom, dist = _get_zatom(zatom[1]), float(zatom[2])
            dvect = dist_atom.position
            ang_atom,  ang  = _get_zatom(zatom[3]), float(zatom[4]) * ang_conv
            avect = ang_atom.position
            #----------------------------------------#
            if not alternate_format:
                tors_atom, tors = _get_zatom(zatom[5]), float(zatom[6]) * ang_conv
                tvect = tors_atom.position
                # Check and make sure they are all different
                if not distinct(tors_atom, dist_atom, ang_atom):
                    raise InvalidZMatrixException("atom indices in z-matrix line '{}' must"\
                                                  " all be distinct".format(
                        zlines[3].strip()
                    ))
                # place the atom: first displace from origin by bond length, then apply rotations
                #   for both the angle and the torsion, then translate it into position
                axis = cross(tvect-avect, dvect-avect)
                pos = (avect - dvect).normalize() * dist
                pos = rotate_point_about_axis(pos, axis, ang)
                tors_axis = dvect - avect
                pos = dvect + rotate_point_about_axis(pos, tors_axis, tors) * dist_units
                atoms.append(_labeled_atom(zatom[0], pos))
                # create BondLength, BondAngle, and Torsion coordinates, if necessary
                if create_representation:
                    coords_to_make.append((BondLength, dist_atom, atoms[-1]))
                    coords_to_make.append((BondAngle, ang_atom, dist_atom, atoms[-1]))
                    coords_to_make.append((Torsion, tors_atom, ang_atom, dist_atom, atoms[-1]))
            #----------------------------------------#
            else:
                raise NotImplementedError("Gaussian alternate Z-matrix format not yet implemented")
        #--------------------------------------------------------------------------------#
        ret_val = Molecule(atoms)
        if create_representation:
            # Now that the atoms are not orphaned, make the coordinates
            coords = []
            for cargs in coords_to_make:
                if len(cargs) == 3:
                    coords.append(cargs[0](cargs[1:], units=dist_units))
                else:
                    coords.append(cargs[0](cargs[1:], units=ang_units))
            # The __ indicates that we ignore the created object, since the InternalRepresentation
            #   constructor itself sets the internal_representation of the Molecule
            __ = InternalRepresentation(
                     ret_val,
                     coords,
                     units={DistanceUnit:dist_units, AngularUnit:ang_units})
        return ret_val

    @classmethod
    def from_zmatrix_with_labels(cls, string, *args, **kwargs):
        """

        @param string: String containing a z-matrix with labels defined at the end using e.g. a234=8.5
        @param args: Passed on to Molecule.from_z_matrix
        @param kwargs: Passed on to Molecule.from_z_matrix
        @return: Molecule object corresponding to the z-matrix passed in
        """

        string = string.strip()
        lines = string.splitlines()
        vars = dict()
        zmat_lines = []
        for line in reversed(lines):
            splitline = re.split(r'\s*,?\s*', line.strip())
            if '=' in line:
                varname, val = line.split('=')
                vars[varname.strip()] = val.strip()
            elif len(splitline) == 2:
                vars[splitline[0]] = splitline[1]
            elif "ariables:" in line:
                continue
            elif len(line.strip()) > 0:
                zmat_lines.insert(0, line)


        new_zmat = ''
        for line in zmat_lines:
            new_vals = []
            for i, entry in enumerate(re.split(r'\s*,?\s*', line.strip())):
                if i == 7:
                    # skip the eighth entry; I don't know what it does
                    continue
                elif i == 0 or i % 2 == 1:
                    new_vals.append(entry.title())
                else:
                    entry_no_minus = entry[1:] if entry[0] == '-' else entry
                    if entry in vars:
                        new_vals.append(vars[entry])
                    elif entry_no_minus in vars:
                        new_vals.append(entry.replace(entry_no_minus, vars[entry_no_minus]))
                    else:
                        new_vals.append(entry)
            new_zmat += ' '.join(new_vals) + "\n"

        # remove double negatives
        new_zmat = new_zmat.replace('--', '')

        return Molecule.from_z_matrix(new_zmat, *args, **kwargs)

    @classmethod
    def linear_alkane(cls, length,
            rCC=1.5299919220184575 * Angstroms,  # C50H102 average using NIH structure
            rCH=1.0900147940154403 * Angstroms,  # C50H102 average for interior CH bonds using NIH structure
            aCCC=109.4721093826136 * Degrees,    # C50H102 average using NIH structure
            aCCH=109.47073572803828 * Degrees,   # C50H102 average using NIH structure
            #tCCCH=60.00042875324615 * Degrees,   # C50H102 average using NIH structure
    ):
        if length <= 2:
            raise NotImplementedError
        zmat = [['H']]
        zmat.append(['C', 0, rCH])
        zmat.append(['H', 1, rCH, 0, aCCH])
        zmat.append(['H', 1, rCH, 2, aCCH, 0, -120.0])
        zmat.append(['C', 1, rCC, 0, aCCC, 2, -120.0])
        zmat.append(['H', 4, rCH, 1, aCCH, 0, 60.0])
        zmat.append(['H', 4, rCH, 1, aCCH, 0, -60.0])
        for i in xrange(1,length):
            zmat.append(['C', 3*i+1, rCC, 3*(i-1)+1, aCCC, max(3*(i-2)+1, 0), 180.0])
            zmat.append(['H', 3*i+4, rCH, 3*i+1, aCCH, max(3*(i-2)+1, 0), 60.0])
            zmat.append(['H', 3*i+4, rCH, 3*i+1, aCCH, max(3*(i-2)+1, 0),-60.0])
        zmat.append(['H', 3*length+1, rCH, 3*(length-1)+1, aCCC, 3*(length-2)+1, 180.0])
        return Molecule.from_z_matrix(zmat, one_based=False)

    @classmethod
    def acene(cls, length,
        rx = 2.42260241888759,
        rcc_cross = 1.41, # approximately
        rcc_side = 1.398821025,
        rch = 1.0800004490739805
    ):
        # sorry for the hard to read code.  I was in a hurry...
        v = 0.5*rcc_cross
        w = v + rcc_side*sin(pi / 6)
        y_coords = [v, w, w+rch]
        x = 0
        atoms = []
        def a(s, ynum):
            atoms.extend([Atom(s, [x, y_coords[ynum], 0]), Atom(s, [x, -y_coords[ynum], 0])])

        y_coords = [v, v+rch*sin(pi / 6), w+rch]
        a("H", 1)
        x += rch*cos(pi/6)
        y_coords = [v, w, w+rch]
        a("C", 0)
        x += rcc_side*cos(pi/6)
        a("C", 1); a("H", 2)
        for _ in xrange(1, length):
            x += rcc_side*cos(pi/6)
            a("C", 0)
            x += rcc_side*cos(pi/6)
            a("C", 1); a("H", 2)
        x += rcc_side*cos(pi/6)
        a("C", 0)
        x += rch*cos(pi/6)
        y_coords = [v, v+rch*sin(pi / 6), w+rch]
        a("H", 1)

        return Molecule(atoms)

    # TODO a "RequiresWeb" decorator to prevent tests from failing when not connected to internet (or perhaps a decorator for the tests themselves?)
    @classmethod
    def from_identifier(cls, *args, **kwargs):
        """ Create a Molecule object using only the name, SMILES, InChIKey, etc.

        The easiest way to call `Molecule.from_identifier()` is to give it a single argument
        that is one of the following properties:
        * smiles
        * stdinchikey
        * stdinchi
        * ncicadd_identifier      # (for FICTS, FICuS, uuuuu)
        * hashisy
        * cas_number
        * chemspider_id           # input must be chemspider_id=1234567
        * opsin_name
        * cir_name
        (These are resolved in roughly this order.)  Visit the NIH CIR documentation for (some) explaination of
        what these identifiers are:
        http://cactus.nci.nih.gov/chemical/structure/documentation

        `Molecule.from_identifier()` can be called with any (single) keyword argument in the
        `grendel.util.web_getter.input_identifiers` list.  In this form, there must be only


        This can also be called using an argument that is the value of a property and a second argument that is
        a list of fields to search for that property in.

        Note:  `Molecule.get()` is a *very* useful alias for this.

        :Examples:

        TODO

        """
        cir_mol = None
        idents = None
        val = None
        # Remove special kwargs from kwargs array here:
        if 'proxy' in kwargs:
            raise NotImplementedError
        if len(args) == 1 and len(kwargs) == 0:
            cir_mol = cirpy.CIRMolecule(args[0])
            idents = '__any__'
            val = args[0]
        elif len(kwargs) > 1:
            raise TypeError("Cannot resolve multiple types of keyword identifiers.")
        elif len(kwargs) == 1:
            if not kwargs.keys()[0] in input_identifiers:
                raise TypeError("Invalid identifier \"" + kwargs.keys()[0] + "\" for contructing molecule.")
            cir_mol = cirpy.CIRMolecule(kwargs[kwargs.keys()[0]], kwargs.keys())
            idents = kwargs.keys()
            val = kwargs[kwargs.keys()[0]]
        elif len(args) == 2 and len(kwargs) == 0 and isinstance(args[0], basestring):
            if isinstance(args[1], basestring):
                cir_mol = cirpy.CIRMolecule(args[0], [args[1]])
                idents = args[1]
                val = args[0]
            elif isinstance(args[1], (list, tuple)):
                cir_mol = cirpy.CIRMolecule(args[0], list(args[1]))
                idents = list(args[1])
                val = args[0]
            else:
                raise TypeError("Invalid arguments for Molecule.from_identifier().  See documentation.")
        else:
            raise TypeError("Invalid arguments for Molecule.from_identifier().  See documentation.")

        # There must be a sdf structure available
        sdf_text = cir_mol.sdf
        if sdf_text is None:
            raise MoleculeNotFoundError(idents, val, cir_mol.sdf_url)
            # TODO:  Utilize the rest of the information available in the SDF file?
        ret_val = Molecule(cir_mol.xyz)
        ret_val._cirmol = cir_mol
        return ret_val
    # Alias
    get = from_identifier


    #############
    #  Methods  #
    #############

    #-------------------------------------#
    # Inquiry methods (which return bool) #
    #-------------------------------------#

    def is_centered(self, tol=1e-8*Angstroms, cartesian_representation=None):
        """ True if the center of mass is at the origin.
        This actually computes the vector from the origin to the center of mass and then determines if the magnitude of
        that vector is less than `tol`.

        :Parameters:

        tol : float or ValueWithUnits
            The maximum 'off-centeredness' that will be tolerated and still return True.  If a `float` is given,
            the units are assumed to be DistanceUnit.default
        cartesian_representation : Representation or None
            Determine if the molecule is centered when represented in `cartesian_representation`.  If
            `None`, just use the molecule's current cartesian representation.

        """
        if cartesian_representation is None:
            tol = strip_units(tol, self.cartesian_units)
            return magnitude(self.center_of_mass()) < tol
        else:
            com = sum(LightVector([p.value for p in pos]) * a.mass
                for pos, a in zip(grouper(3, cartesian_representation), self))
            com /= sum(atom.mass for atom in self)
            return magnitude(com) < tol

    def is_linear(self, tol=None):
        """ True if the molecule is linear to within `tol`.  All diatomics should return True.

        If `tol` is a ValueWithUnits and `tol.units` is an `AngularUnit`, then this method returns True only
        if *all* angles in the molecule are within `tol` of 180 Degrees.  If `tol.units` is a unit-compatible with
        a moment of inertia (i.e. MassUnit * DistanceUnit**2 ), then this method returns True if the smallest principal
        moment of inertia is less than `tol` *and* the difference between the two largest principal moments of
        inertia is less than `tol`.  If no units are given (i.e. `tol` is a `float` or other float-compatible unit),
        `tol` is assumed to have units of AngularUnit.default and the method proceeds as if `tol.units` was an
        `AngularUnit` subclass.

        :Parameters:

        tol : float or ValueWithUnits
            The linearity tolerance.

        :Examples:

        >>> from grendel.chemistry import SampleMolecules, init_sample_molecules
        >>> init_sample_molecules()
        >>> SampleMolecules['water'].is_linear()
        False
        >>> SampleMolecules['CO2'].is_linear()
        True
        >>> SampleMolecules['Benzene'].is_linear()
        False

        """
        if self.natoms == 2:
            return True
        tol = tol or Molecule.linear_cutoff
        tol = (tol * AngularUnit.default) if not isinstance(tol, ValueWithUnits) else tol
        thresh = strip_units(tol, Radians)
        def _tmp(a, b, c):
            rv = angle_between_vectors(b.pos - a.pos, b.pos - c.pos)
            return min(abs(rv-math.pi), abs(rv))
        return all(_tmp(a, b, c) < thresh for a in self for b in self for c in self if a != b != c != a)

    def has_same_elements(self, other):
        """ Returns True if the elements of self correspond directly to the elements of other
        (i.e. `self.atoms[0].element == other.atoms[0].element and self.atoms[1].element == other.atoms[1].element and ...`)
        """
        if self.natoms != other.natoms:
            return False
        for atom, other_atom in zip(self, other):
            if atom.element != other_atom.element:
                return False
        return True

    @typecheck(other='Molecule', tol=(None, Number))
    def has_same_geometry(self, other, tol=None):
        """
        Returns True if all of the elements are the same (including isotope and nuclear spin)
        and the `reoriented()` versions of `self` and `other` have no atoms whose pairwise
        position difference has a magnitude greater than `tol` (which defaults to 1e-8 Angstroms)

        :Examples:

        >>> mol = Molecule('''
        ...     O 1.5 0.0 0.0
        ...     H 0.2 0.0 0.0
        ... ''')
        >>> m1 = Molecule('''
        ...     O 0.0 1.5 0.0
        ...     H 0.0 0.2 0.0
        ... ''')
        >>> m2 = Molecule('''
        ...     O 0.00000 0.00000 0.200000
        ...     H 0.00000 0.00000 1.500001
        ... ''')
        >>> m3 = Molecule('''
        ...     H 0.0  0.0 0.0
        ...     H 0.0 -1.3 0.0
        ... ''')
        >>> mol.has_same_geometry(m1)
        True
        >>> m1.has_same_geometry(m2)
        False
        >>> m1.has_same_geometry(m2, 1e-5)
        True
        >>> m2.has_same_geometry(mol, 1e-5)
        True
        >>> m1.has_same_geometry(m3)
        False

        """
        if tol is None:
            tol = 1e-8*Angstroms
        failed = False
        # If tol has units, turn it into our units.  If not, assume it's already in our units.
        if hasunits(tol):
            tol = strip_units(tol, self.cartesian_units)
        # Check for same element composition
        if not self.has_same_elements(other):
            return False
        # Create reoriented versions of self and other:
        me = self.reoriented()
        other = other.reoriented()
        # Compare these reoriented versions atom-wise
        for atom, other_atom in zip(me, other):
            apos = strip_units(atom.pos, DistanceUnit.default)
            opos = strip_units(other_atom.pos, DistanceUnit.default)
            if magnitude(apos - opos) > tol:
                failed = True
                break
        # Since reorienting is only specified up to a phase factor (i.e. up to a sign), we have to
        #   check the negative version as well:
        if failed:
            failed = False
            for atom, other_atom in zip(me, other):
                apos = strip_units(atom.pos, DistanceUnit.default)
                opos = strip_units(other_atom.pos, DistanceUnit.default)
                if magnitude(apos + opos) > tol:
                    failed = True
                    break
            if failed:
                return False
        # No failures, they must be the same
        return True

    @typecheck(op=SymmetryOperation)
    def has_symmetry(self, op):
        """ True if `op` is a valid symmetry operation on `self`
        Each coordinate of the transformed atoms must differ from the original by less than
        """
        new_atoms = []
        for atom in self:
            new_atom = Atom(atom.symbol, op.matrix * atom.position)
            new_atom.same_coordinate_tolerance = PointGroup.symmetry_tolerance
            new_atoms.append(new_atom)
            atom.same_coordinate_tolerance = PointGroup.symmetry_tolerance
        retval = (Molecule(new_atoms) == self)
        # Change things back
        for atom in self:
            atom.same_coordinate_tolerance = Atom.same_coordinate_tolerance
        return retval

    @typecheck(details=('ComputationDetails', None))
    def has_energy(self, details=None):
        return self.has_property(Energy, details)

    def has_property(self, property, details=None):
        return any(rg.has_property_for_molecule(self, property, details) for rg in self.result_getters)

    def can_get_energy(self, details=None):
        return self.can_get_property(Energy, details)

    def can_get_property(self, property, details=None):
        return any(rg.can_get_property_for_molecule(self, property, details) for rg in self.result_getters)

    #------------------------------#
    # Computed and derived results #
    #------------------------------#

    #TODO Make this get_computations_for_properties and find minimal number of computations that get everything needed
    def get_computation_for_property(self, property, details=None):
        for rg in self.result_getters:
            if isinstance(rg, ComputationResultGetter):
                if rg.can_get_property_for_molecule(self, property, details):
                    return rg.get_computation_for_property(self, property, details)
        # If we've gotten here, we've failed.
        raise RuntimeError("Don't know how to generate a computation to get {prop}"
                           " for molecule:\n{0}\nwith details:\n{1}".format(
            self, details,
            prop=MolecularProperty.property_type_of(prop).__name__
        ))

    def get_energy(self, details=None, run_computation=False):
        return self.get_property(Energy, details, run_computation)

    def get_property(self, property, details=None, run_computation=False):
        # Check if we have the property already...
        for p in self._properties:
            if MolecularProperty.is_same_property(property, p):
                if ComputationDetails.is_compatible_details(details, p.details):
                    return p
        # Check if any of the ResultGetters have the property already:
        for rg in self.result_getters:
            if rg.has_property_for_molecule(self, property, details):
                rv = rg.get_property_for_molecule(self, property, details)
                self._properties.append(rv)
                return rv
        # Don't have it already...darn...see if we can get it
        if run_computation:
            comp = self.get_computation_for_property(property, details)
            if comp.started and not comp.completed:
                raise RuntimeError("requested computation has been started, but is not completed."
                                   "  Make sure you wait for all computations to finish before"
                                   " querying them for properties.")
            else:
                comp.run()
                ret_val = comp.get_property(property, details)
                if ret_val is None:
                    raise RuntimeError("Did not successfully get property {} for Molecule:\n{}\nwith"
                                       " details:\n{}\n(Computation.get_property()"
                                       " returned 'None')".format(
                        MolecularProperty.property_type_of(property).__name__, self, details
                    ))
                else:
                    self._properties.append(ret_val)
                return ret_val
        else:
            raise RuntimeError("Don't already have {prop} for molecule:\n{0}\nwith"
                               " details:\n{1},\n and automatic running of computations"
                               " to obtain this property is disabled".format(
                self, details,
                prop=MolecularProperty.property_type_of(property).__name__
            ))
        # If we've gotten here, we've failed.
        raise RuntimeError("Don't know how to get {prop} for molecule:\n{0}\nwith"
                           " details:\n{1}".format(
            self, details,
            prop=MolecularProperty.property_type_of(property).__name__
        ))

    def get_optimized_geometry(self, details=None, property=None):
        raise NotImplementedError

    #----------------------------------------#
    # Geometric properties and modifications #
    #----------------------------------------#

    @cached_method
    def center_of_mass(self):
        """ Returns a Vector giving the center of mass of the molecule in the current Cartesian representation.
        The units of the returned value are self.cartesian_units.
        .. note::
            This result of this method is cached, and the cached value gets flushed in
            `update_cartesian_representation()`.  If you change an atom's position
            (or mass) and forget to call `update_cartesian_representation()`, you may
            get some funny results for this method or any methods that depend on it,
            including `recenter()`, `reorient()`, and `principal_moments_of_inertia()`.  You
            can detect whether caching is causing problems by setting the environment variable
            `GRENDEL_NO_CACHE` to 1 and rerunning your tests.  If tests that were failing
            subsequently succeed, you probably forgot to call `update_cartesian_representation()`
            somewhere, or you were assuming that it was automatically called somewhere when
            in fact it was not getting called.
        """
        center = Vector([0.0, 0.0, 0.0])
        for atom in self:
            center += atom.position * atom.mass
        center /= self.mass
        return center

    def recenter(self):
        """ Recenters the molecule about the center of mass
        This modifies the molecule in place.

        :Examples:


        >>> m = Molecule([Atom("O", [1.1, 1.3, 1.7])])
        >>> m
        Molecule([
            Atom('O', [  1.10000000,  1.30000000,  1.70000000 ] )
        ])
        >>> m.recenter()
        >>> m
        Molecule([
            Atom('O', [  0.00000000,  0.00000000,  0.00000000 ] )
        ])
        >>> m = Molecule(['H','H'], Matrix([[0,0,0],[1,0.0,0]]))
        >>> m
        Molecule([
            Atom('H', [  0.00000000,  0.00000000,  0.00000000 ] ),
            Atom('H', [  1.00000000,  0.00000000,  0.00000000 ] )
        ])
        >>> m.recenter()
        >>> m
        Molecule([
            Atom('H', [ -0.50000000,  0.00000000,  0.00000000 ] ),
            Atom('H', [  0.50000000,  0.00000000,  0.00000000 ] )
        ])

        """
        # This takes a little longer to make this check first, but it
        # lessens the severity of round-off error
        if self.is_centered(): return
        center = self.center_of_mass()
        for atom in self:
            atom.position -= center
        self.update_cartesian_representation()
        return

    def recentered(self):
        """ Same as `recenter`, but makes returns a copy.  `self` is not modified.

        :Examples:


        >>> foo = Molecule([Atom("O", [1.3, 1.7, 2.0])])
        >>> bar = foo.recentered()
        >>> foo.recenter()
        >>> foo
        Molecule([
            Atom('O', [  0.00000000,  0.00000000,  0.00000000 ] )
        ])
        >>> foo == bar
        True
        >>> foo is bar
        False

        """
        ret_val = deepcopy(self)
        ret_val.recenter()
        return ret_val

    # TODO Make this work consistently for symmetric tops
    # TODO Standard Nuclear Orientation (using nuclear charge) from CPL 209:506 (1993)
    def reorient(self, representation = "II"):
        """ Reorient the molecule to align the x, y, and z axes with the principal axes of rotation.

        :Parameters:

        representation : str, optional
            Must be one of "I", "II", or "III".  "I" means {x,y,z} = {b,c,a}.  "II" (the default) means
            {x,y,z} = {c,a,b}.  "III" means {x,y,z} = {a,b,c}

        :Examples:


        >>> from grendel.chemistry import SampleMolecules, init_sample_molecules
        >>> init_sample_molecules()
        >>> from grendel.gmath import chopped
        >>> h2o = SampleMolecules['quantum water']
        >>> h2o.principal_axes()
        Matrix([[-0.        ,  0.        ,  1.        ],
                [-0.53654222, -0.84387348,  0.        ],
                [-0.84387348,  0.53654222,  0.        ]])
        >>> h2o.reorient()
        >>> # Use chopped to get rid of very small numerical artifacts...
        >>> # using abs to get rid of phase factor
        >>> abs(chopped(h2o.principal_axes()))
        Matrix([[ 0.,  0.,  1.],
                [ 1.,  0.,  0.],
                [ 0.,  1.,  0.]])
        >>> h2o.reorient("I")
        >>> from grendel.gmath import chopped
        >>> abs(chopped(h2o.principal_axes()))
        Matrix([[ 0.,  1.,  0.],
                [ 0.,  0.,  1.],
                [ 1.,  0.,  0.]])
        >>> h2o.reorient("III")
        >>> abs(chopped(h2o.principal_axes()))
        Matrix([[ 1.,  0.,  0.],
                [ 0.,  1.,  0.],
                [ 0.,  0.,  1.]])

        """
        nx, ny, nz = Vector(1,0,0), Vector(0,1,0), Vector(0,0,1)
        i = self.principal_axes()
        a, b, c = self.a, self.b, self.c
        if representation == "II":
            x = c
            y = a
            z = b
        elif representation == "I":
            x = b
            y = c
            z = a
        elif representation == "III":
            x = a
            y = b
            z = c
        else: # pragma: no cover
            raise ValueError("Argument to reorient must be one of \"I\", \"II\", or \"III\"")
            # See http://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions#Rotation_matrix_.E2.86.94_Euler_angles
        ca = lambda n, w: math.cos(angle_between_vectors(n, w))
        rot_mat = Matrix([
            [ca(nx, x), ca(ny, x), ca(nz, x)],
            [ca(nx, y), ca(ny, y), ca(nz, y)],
            [ca(nx, z), ca(ny, z), ca(nz, z)]
        ])
        for atom in self:
            atom.position = rot_mat * atom.position
        self.update_cartesian_representation()
        return

    def reoriented(self, representation = "II"):
        """ Same as `reorient`, but makes returns a copy.  `self` is not modified.

        :Examples:


        >>> from grendel.chemistry import SampleMolecules, init_sample_molecules
        >>> init_sample_molecules()
        >>> from grendel.gmath import chopped
        >>> h2o = Molecule([
        ...    Atom('O', [  0.00000000,  0.00000000,  0.00000000 ] ),
        ...    Atom('H', [  0.00000000,  0.00000000,  1.00000000 ] ),
        ...    Atom('H', [  0.00000000, -1.00000000, -0.00000000 ] )
        ... ])
        >>> new = h2o.reoriented()
        >>> # Use chopped to get rid of very small numerical artifacts...
        >>> chopped(new.principal_axes())
        Matrix([[ 0.,  0.,  1.],
                [-1.,  0.,  0.],
                [ 0.,  1.,  0.]])
        >>> new is h2o
        False
        >>> # original molecule is unchanged...
        >>> h2o
        Molecule([
            Atom('O', [  0.00000000,  0.00000000,  0.00000000 ] ),
            Atom('H', [  0.00000000,  0.00000000,  1.00000000 ] ),
            Atom('H', [  0.00000000, -1.00000000, -0.00000000 ] )
        ])

        """
        ret_val = deepcopy(self)
        ret_val.reorient(representation)
        return ret_val

    def rotate(self, axis, angle):
        """ Rotates the molecule about axis (a `Vector`) by angle

        """
        angle = strip_units(angle, Radians, assume_units=AngularUnit.default)
        for atom in self:
            atom.position = rotate_point_about_axis(atom.position, axis, angle)
        self.update_cartesian_representation()
        return

    def rotated(self, axis, angle):
        ret_val = deepcopy(self)
        ret_val.rotate(axis, angle)
        return ret_val

    def translate(self, translation):
        translation = strip_units(translation, self.cartesian_units, assume_units=self.cartesian_units)
        for atom in self:
            atom.position += translation
        self.update_cartesian_representation()
        return self

    def inertial_system(self):
        """ Returns a tuple of the principal moments of inertia vector and the principal axes matrix.

        mol.inertial_system()

        Computes the principal moments of inertia and the principal axes.

        .. note::
           This funtion `mol.recenter()` if the molecule is not centered, so any `CartesianRepresentations`
           that are not both `frozen` and referenced elsewhere (e.g. in a `RepresentationDependentProperty`)
           will be lost forever.

        .. note::
           This method is cached.  See discussion of the consequences of caching in `Molecule.center_of_mass()`

        :Returns:

        A tuple with types (`Vector`, `Matrix`) containing the principal moments of inertia and the principal axes,
        respectively.  These should be aligned (i.e. ret_val[0][1] corresponds to the vector ret_val[1][1])

        :Examples:


        >>> hnco = Molecule.from_z_matrix(\"""
        ... N
        ... C N 1.2145
        ... O C 1.1634 N 172.22
        ... H N 1.0030 C 123.34 O 180.0
        ... \"""
        ... )
        >>> i = hnco.inertial_system()
        >>> i[0]
        Vector([  0.60187342,  45.48378728,  46.0856607 ])
        >>> i[1]
        Matrix([[ 0.        ,  0.        ,  1.        ],
                [ 0.09870438, -0.9951168 ,  0.        ],
                [ 0.9951168 ,  0.09870438,  0.        ]])
        >>> j = hnco.inertial_system()
        >>> i == j
        True
        >>> hnco.A_e
        28.0085955528 Wavenumber
        >>> hnco.B_e
        0.370629408869 Wavenumber
        >>> hnco.C_e
        0.365789031483 Wavenumber
        >>> hnco.a
        Vector([ 0.        ,  0.09870438,  0.9951168 ])
        >>> hnco.b
        Vector([ 0.        , -0.9951168 ,  0.09870438])
        >>> hnco.c
        Vector([ 1.,  0.,  0.])

        """
        if caching_enabled and self._pmi is not None and self._principal_axes is not None:
            return self._pmi, self._principal_axes
        self.recenter()
        imat = Matrix(shape=(3,3))
        for i in xrange(self.natoms):
            a = self.atoms[i]
            v = a.position
            imat[1,2] -= a.mass * v.y * v.z
            imat[0,1] -= a.mass * v.x * v.y
            imat[0,2] -= a.mass * v.x * v.z
            imat[0,0] += a.mass * (v.y**2 + v.z**2)
            imat[1,1] += a.mass * (v.x**2 + v.z**2)
            imat[2,2] += a.mass * (v.x**2 + v.y**2)
        imat[2,1] = imat[1,2]
        imat[1,0] = imat[0,1]
        imat[2,0] = imat[0,2]
        evals, evecs = np.linalg.eigh(imat)
        evals = evals.view(Vector)
        evecs = evecs.view(Matrix)
        esys = [(evals[i], evecs[:,i]) for i in range(3)]
        esys.sort(key=lambda x: x[0])
        self._pmi = Vector([x[0] for x in esys])
        self._principal_axes = Matrix([x[1] for x in esys]).T
        return self._pmi, self._principal_axes

    def pmi(self):
        """ The principal moments of inertia, as a `Vector`

        .. note::
           This funtion `mol.recenter()` if the molecule is not centered, so any `CartesianRepresentations`
           that are not both `frozen` and referenced elsewhere (e.g. in a `RepresentationDependentProperty`)
           will be lost forever.

        .. note::
           This method is cached.  See discussion of the consequences of caching in `Molecule.center_of_mass()`

        :Examples:


        >>> hnco = Molecule.from_z_matrix(\"""
        ... N
        ... C N 1.2145
        ... O C 1.1634 N 172.22
        ... H N 1.0030 C 123.34 O 180.0
        ... \"""
        ... )
        >>> hnco.principal_moments_of_inertia()
        Vector([  0.60187342,  45.48378728,  46.0856607 ])
        >>> hnco.pmi()
        Vector([  0.60187342,  45.48378728,  46.0856607 ])


        :See Also:

        principal_axes, inertial_system, A_e, B_e, C_e

        """
        return self.inertial_system()[0]
    principal_moments_of_inertia = function_alias('principal_moments_of_inertia', pmi)

    @cached_method
    def principal_axes(self):
        """ The principal axes as column vectors in a `Matrix`.
        The `Vector` object `mol.principal_axes()[:,i]` corresponds to the `i`th moment of inertia,
        `mol.principal_moments_of_inertia()[i]`.

        mol.principal_axes()

        .. note::
           This funtion `mol.recenter()` if the molecule is not centered, so any `CartesianRepresentations`
           that are not both `frozen` and referenced elsewhere (e.g. in a `RepresentationDependentProperty`)
           will be lost forever.

        .. note::
           This method is cached.  See discussion of the consequences of caching in `Molecule.center_of_mass()`

        :Examples:


        >>> hnco = Molecule.from_z_matrix(\"""
        ... N
        ... C N 1.2145
        ... O C 1.1634 N 172.22
        ... H N 1.0030 C 123.34 O 180.0
        ... \"""
        ... )
        >>> hnco.principal_axes()
        Matrix([[ 0.        ,  0.        ,  1.        ],
                [ 0.09870438, -0.9951168 ,  0.        ],
                [ 0.9951168 ,  0.09870438,  0.        ]])


        :See Also:

        principal_moments_of_inertia, inertial_system, A_e, B_e, C_e
        """
        return self.inertial_system()[1]

    def largest_difference_with(self, other):
        me = self.reoriented()
        other = other.reoriented()
        largest_pos = float('-inf') * DistanceUnit.default
        for atom, other_atom in zip(me, other):
            apos = strip_units(atom.pos, DistanceUnit.default)
            opos = strip_units(other_atom.pos, DistanceUnit.default)
            mag = magnitude(apos - opos)
            if mag > largest_pos:
                largest_pos = mag
            # Try the opposite sign:
        largest_neg = float('-inf') * DistanceUnit.default
        for atom, other_atom in zip(me, other):
            apos = strip_units(atom.pos, DistanceUnit.default)
            opos = strip_units(other_atom.pos, DistanceUnit.default)
            mag = magnitude(apos + opos)
            if magnitude(apos + opos) > largest_neg:
                largest_neg = mag
        return (min(largest_pos, largest_neg) * DistanceUnit.default).in_units(self.cartesian_units)

    def update_cartesian_representation(self):
        self._clear_cartesian_dependent_cached_variables()
        self.cartesian_representation = CartesianRepresentation(self, units=self.cartesian_units)
        for atom in self.atoms:
            self.cartesian_representation.add_atom(atom)

    @typechecked(cart_vect=Vector)
    def displace(self, cart_vect):
        if sanity_checking_enabled:
            if len(cart_vect) != self.natoms * 3:
                raise ValueError("dimension mismatch.  Expected cartesian displacement vector "\
                                 "with {0} elements, got one with {1} elements".format(self.natoms*3, len(cart_vect)))
        for atom, disp in zip(self.atoms, grouper(3, cart_vect)):
            disp = LightVector(disp)
            atom.displace(disp)
        self.update_cartesian_representation()


    #----------------------------#
    # Other unclassified methods #
    #----------------------------#

    def copy_with_atoms(self,
            new_atoms,
            deep_copy=False,
            deep_copy_memo=None,
            new_charge=None,
            new_multiplicity=None,
            new_description=None
    ):
        if deep_copy:
            new_atoms = deepcopy(new_atoms, deep_copy_memo)
        else:
            new_atoms = [copy(atom) for atom in new_atoms]
        #----------------------------------------#
        copied = Molecule(
            new_atoms,
            units=self.cartesian_units
        )
        #----------------------------------------#
        for atom in copied:
            atom.parent_molecule = copied
        #----------------------------------------#
        # copy attributes...
        if new_description is None:
            copied.description = self.description
        else:
            copied.description = new_description
        if new_charge is None:
            copied.charge = self.charge
        else:
            copied.charge = new_charge
        if new_multiplicity is None:
            copied.multiplicity = self.multiplicity
        else:
            copied.multiplicity = new_multiplicity
        #----------------------------------------#
        return copied

    def fragment(self,
            atom_numbers,
            charge=None,
            multiplicity=None,
            description=None
    ):
        if description is None:
            description = ("fragment of " + self.description) if self.description is not None else None
        fatoms = [self.atoms[n] for n in atom_numbers]
        return self.copy_with_atoms(fatoms,
            new_charge=charge,
            new_multiplicity=multiplicity,
            new_description=description
        )

    def copy_without_ghost_atoms(self):
        return self.copy_with_atoms([a for a in self if not a.is_ghost()])

    def index(self, atom):
        """ Returns the index of `atom` in the atoms array of the molecule

        .. note::
           This method returns a cached property `Atom.index`.  If you reorder the atoms in a molecule,
           be sure and flush this cache by setting `atom._index` to `None` for all of the atoms in the
           reordered molecule.

        :Raises:

        IndexError : if `atom` is not found in `self`

        :Examples:

        >>> from grendel.chemistry import SampleMolecules, init_sample_molecules
        >>> init_sample_molecules()
        >>> h2o = SampleMolecules['water']
        >>> atom1 = h2o[0]
        >>> atom2 = h2o[1]
        >>> atom3 = h2o[2]
        >>> h2o
        Molecule([
            Atom('O', [  0.00000000,  0.00000000,  0.11815400 ] ),
            Atom('H', [  0.00000000,  0.75873400, -0.47261400 ] ),
            Atom('H', [  0.00000000, -0.75873400, -0.47261400 ] )
        ])
        >>> h2o.index(atom2)
        1
        >>> h2o.index(atom1)
        0
        >>> # It must be exactly the same instance to avoid raising an index error
        >>> h2o.index(Atom('O', [  0.00000000,  0.00000000,  0.11815400 ] )) # doctest: +ELLIPSIS
        Traceback (most recent call last):
            ...
        IndexError: ...

        """
        if sanity_checking_enabled:
            if atom.parent_molecule is not self:
                raise IndexError(repr(atom) + " is not in \n" + repr(self))
        return atom.index

    # TODO figure out how to test this
    def write_xyz(self, filename, overwrite=False, format_str = "%-3s %12.8f %12.8f %12.8f"):
        """ Writes the Molecule to the standard xyz format.

        See `Molecule.xyz_string() for more details.

        :See Also:

        xyz_string

        """
        filename = full_path(filename)
        if not overwrite and os.access(filename, os.F_OK):
            raise IOError("%s already exists. Use 'overwrite=True' to overwrite it." % filename)
        file = open(filename, "w+")
        file.write(self.xyz_string(format_str))
        file.close()
        return

    def xyz_string(self, format_str="%-3s %12.8f %12.8f %12.8f", header=True):
        """ The molecule, represented as a string in standard xyz format.

        The first line will always be the number of atoms.
        The second line (which is a comment in the xyz specification) is filled by the
        description first line of the `description` attribute.

        :Parameters:

        format_str : str, optional
            The format to apply to the lines of the xyz output.  The default is `"%-3s %12.8f %12.8f %12.8f"`, which
            should be fine for most purposes.
        header : bool, optional
            Whether or not to include the standard two-line header which is part of the standard xyz format (defaults to True)


        :Examples:


        >>> from __future__ import print_function
        >>> from grendel.chemistry import SampleMolecules, init_sample_molecules
        >>> init_sample_molecules()
        >>> print(SampleMolecules['water'].xyz_string())
        3
        CCSD(T)/aug-cc-pVTZ Water
        O     0.00000000   0.00000000   0.11815400
        H     0.00000000   0.75873400  -0.47261400
        H     0.00000000  -0.75873400  -0.47261400
        >>> print(Molecule("O 1.0 0.5 1.5").xyz_string("Atom: %2s, x: %3.1f, y: %4.2f, z: %5.3f"))
        1
        <BLANKLINE>
        Atom:  O, x: 1.0, y: 0.50, z: 1.500
        >>> print(Molecule("O 1.0 0.5 1.5").xyz_string("Atom: %2s, x: %3.1f, y: %4.2f, z: %5.3f", False))
        Atom:  O, x: 1.0, y: 0.50, z: 1.500

        """
        ret_val = ""
        if header:
            ret_val = str(self.natoms) + "\n"
            if not self.description is None:
                ret_val += self.description + "\n"
            else:
                ret_val += "\n"
        for atom in self:
            ret_val += (format_str % (atom.symbol, atom.x, atom.y, atom.z)) + "\n"
        return ret_val[:-1]

    def convert_units(self, new_units):
        for atom in self:
            atom.convert_units(new_units)
        self.cartesian_units = new_units
        self.update_cartesian_representation()

    def displacement_description(self,
            eq_name="eq",
            divider="_",
            delta_name='',
            one_based=True,
            include_zeros=False):
        ret_val = ''
        if self.displacement is None:
            return eq_name
        incs = self.displacement.increments()
        if all(d == 0 for d in incs):
            return eq_name
        #----------------------------------------#
        disp_descs = []
        for disp, coord in zip(incs, self.displacement.representation):
            if disp != 0 or include_zeros:
                disp_descs.append(coord.generate_name(one_based) + "{:+d}".format(disp) + delta_name)
        return divider.join(disp_descs)


    def use_result_getter(self, rg):
        #TODO set a flag that generates a warning if a different ResultGetter must be used.
        self.result_getters.insert(0, rg)


    def get_stub(self):
        raise NotImplementedError

    ###################
    # Private Methods #
    ###################

    def _clear_cartesian_dependent_cached_variables(self):
        self._pmi = None
        self._principal_axes = None
        self._center_of_mass = None
        self._A_e = self._B_e = self._C_e = None
        self._a = self._b = self._c = None
        self._is_linear = None
        if hasattr(self, '_reoriented_matrix'):
            self._reoriented_matrix = None
        if hasattr(self, '_pseudo_hash'):
            self._pseudo_hash = None
        if hasattr(self, '_pseudo_hash_lower'):
            self._pseudo_hash = None
        if hasattr(self, '_pseudo_hash_upper'):
            self._pseudo_hash = None

    def _compute_point_group(self):
        """

        >>> from grendel.chemistry import SampleMolecules, init_sample_molecules
        >>> init_sample_molecules()
        >>> pg = SampleMolecules['benzene'].point_group   # doctest: +SKIP
        >>> len(pg) # doctest: +SKIP
        24
        >>> len(pg.classes) # doctest: +SKIP
        12
        >>> pg.classes # doctest: +SKIP
        [{E}, {2 C_6}, {2 C_6^2}, {C_6^3}, {3 C_2(z)}, {3 C_2'}, {i}, {2 S_6(x)}, {2 S_3}, {sigma_h(yz)}, {3 sigma_v(xy)}, {3 sigma_v(xz)}]


        """
        if self.is_linear(2.0 * Degrees):
            raise NotImplementedError("Symmetry of linear molecules not yet impmlemented.")

        self.recenter()
        self.reorient()
        self.inertial_system()

        ops = []

        # First, make a set of trial axes
        trial_axes = [
            [self._principal_axes[:,0], 2],
            [self._principal_axes[:,1], 2],
            [self._principal_axes[:,2], 2]
        ]

        def append_if_unique(vect, max_order):
            if magnitude(vect) < SymmetryOperation.zero_vector_cutoff:
                return
            same_axes = [axis for axis in trial_axes if SymmetryOperation.is_same_axis(vect, axis[0])]
            if len(same_axes) == 0:
                trial_axes.append([vect.normalized(), max_order])
            elif max_order > same_axes[0][1]:
                same_axes[0][1] = max_order
            return

        natoms =self.natoms
        for a1 in self:
            for a2 in self:
                if a1 == a2: continue
                # Append the bond midpoints...
                b1 = a1.pos - a2.pos
                mid1 = (a1.pos + a2.pos)/2.0
                append_if_unique(mid1, 2)
                for a3 in self:
                    if a1 == a3 or a2 == a3: continue
                    b2 = a3.pos - a2.pos
                    ang = angle_between_vectors(b1, b2)
                    mid2 = (a2.pos + a3.pos)/2.0
                    could_be_ring = False
                    max_ring = 2
                    for x in xrange(3, natoms):
                        ring_angle = ((x-2)*math.pi)/x
                        if abs(ring_angle - ang) < PointGroup.ring_angle_tolerance:
                            could_be_ring = True
                            max_ring = x
                            break
                        elif ang < ring_angle:
                            break
                    if could_be_ring:
                        n = cross(b1, b2)
                        if magnitude(n) > SymmetryOperation.zero_vector_cutoff:
                            n.normalize()
                            append_if_unique(n, max_ring)
                            # See http://mathforum.org/library/drmath/view/62814.html for simple derivation of the following
                            #v1, v2 = cross(n, b1), cross(n, b2)
                            #a = magnitude(cross(mid2 - mid1, v2)) / magnitude(cross(v1, v2))
                            #vect = mid1 + (a * v1)
                            #append_if_unique(vect, max_ring)
                    if a2.pos.magnitude() < SymmetryOperation.same_axis_tolerance:
                        # Fixes methane...
                        append_if_unique(cross(b1, b2), 2)

        # Now test them...
        ops.append(IdentityOperation())
        if Inversion.exists_for_molecule(self):
            ops.append(Inversion())
        for pair in trial_axes:
            axis = pair[0]
            max_ring = pair[1]
            rots = Rotation.about_axis(self, axis, max_ring)
            for op in rots:
                if not any(op == existing for existing in ops):
                    ops.append(op)
            refl = Reflection.with_normal(self, axis)
            if refl is not None and not any(refl == existing for existing in ops):
                ops.append(refl)
            improper_rots = ImproperRotation.about_axis(self, axis, max_ring)
            for op in improper_rots:
                if not any(op == existing for existing in ops):
                    ops.append(op)
        self._point_group = PointGroup(self, ops)


class MoleculeStub(Molecule):
    """ The exact same thing as the parent Molecule class (for now, anyway), but used for holding
    specifications of molecules that can later be replaced by a full molecule class (so things like,
    for instance, the `calculations` attribute won't be filled in, but later if the geometry is
    found to match, the class using the `MoleculeStub` in some list can replace the stub with the full
    `Molecule` instance upon matching.
    """

    ####################
    # Class Attributes #
    ####################

    eq_precision = 7
    hash_precision = 5

    ##############
    # Attributes #
    ##############

    multiplicity_specified = None
    charge_specified = None
    units_specified = None

    ######################
    # Private Attributes #
    ######################

    _reoriented_matrix = None
    _hash = None

    ##################
    # Initialization #
    ##################

    @with_flexible_arguments(
        optional=[
            ('units', 'cartesian_units'),
            ('charge',),
            ('multiplicity',),
        ],
        what_to_call_it='MoleculeStub constructor'
    )
    def __init__(self, *args, **kwargs):
        super(MoleculeStub, self).__init__(*args, **kwargs)
        # Keep track of what attributes were specified.  If an attribute is
        # specified, it must be matched in "is_valid_stub_for()"
        self.multiplicity_specified = 'multiplicity' in kwargs
        self.charge_specified = 'charge' in kwargs
        self.units_specified = 'units' in kwargs
        if self.multiplicity_specified or self.charge_specified:
            raise NotImplementedError

    #################
    # Class Methods #
    #################

    @classmethod
    def stub_for(cls, molecule):
        return MoleculeStub(
            atom_names=[a.symbol for a in molecule.atoms],
            cart_mat=molecule.position,
            copy_atoms=True
        )

#####################
# Dependent Imports #
#####################

from grendel.chemistry.molecular_properties import Energy, MolecularProperty
from grendel.chemistry import element_data
from grendel.chemistry.atom import Atom
from grendel.chemistry.molecule_dict import MoleculeDict

from grendel.interface.result_getter import ComputationResultGetter
from grendel.interface.computation_details import ComputationDetails

from grendel.representations.cartesian_representation import CartesianRepresentation
from grendel.representations.internal_representation import InternalRepresentation

from grendel.coordinates.bond_angle import BondAngle
from grendel.coordinates.bond_length import BondLength
from grendel.coordinates.torsion import Torsion

