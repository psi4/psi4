from copy import deepcopy
from functools import total_ordering
from numbers import Real
import numpy as np

from grendel import type_checking_enabled, sanity_checking_enabled
from grendel.chemistry.element import Element, Isotope

from grendel.gmath.vector import Vector, LightVector, magnitude

from grendel.util.decorators import with_flexible_arguments, CachedProperty, IterableOf, typechecked
from grendel.util.metaprogramming import ReadOnlyAttribute
from grendel.util.overloading import overloaded, OverloadedFunctionCallError
from grendel.util.units.unit import DistanceUnit, isunit, Angstroms
from grendel.util.units.value_with_units import hasunits, strip_units

__all__ = [
    "Atom"
]

#TODO Decide if we want to make the parent_molecule a weakref or just control its pickling/copying behavour
#TODO Figure out what the most reasonable way to handle copying should be with respect to the parent_molecule attribute (and implement it)
@total_ordering
class Atom(object):
    """
    Encapsulates an atom.  Molecules are made up of Atoms.

    :Attributes:

    position : `Vector`
        The atom's Cartesian coordinates [x y z] as a :py:class:`~grendel.math.vector.Vector` object.
    parent_molecule : `Molecule`
        Reference back to the parent molecule of the atom.  (Not set in initializer, but set when the atom gets
        included in a `Molecule`).  Note that for purposes of pickling and copying this may become a weakref in the
        future.  Do not implement behavior that depends on it being a strong reference.
    zmat_label : `str`
        Used in the construction of a z-matrix to differentiate atoms of the same Element.  Currently not used anywhere
        else, and not set when any other Atom construction method is used.



    """

    ####################
    # Class Attributes #
    ####################

    same_coordinate_tolerance = 1e-6  # Affects the __lt__ and __eq__ functions

    ################
    #  Attributes  #
    ################

    position = None
    parent_molecule = None
    zmat_label = None

    base_atom = None
    """
    If the `Atom` instance is a displaced copy of another `Atom`, the original is held here.
    """

    cartesian_units = ReadOnlyAttribute('cartesian_units')


    ########################
    #  Private Attributes  #
    ########################

    _element = None
    _isotope = None
    _orphaned_cart_coords = None


    ##################
    # Initialization #
    ##################

    @overloaded
    @with_flexible_arguments(
        optional=[
            ('units', 'cartesian_units')
        ],
        what_to_call_it='Atom constructor'
    )
    def __init__(self, *args, **kwargs):
        """
        Atom(*args, **kwargs)

        **Signatures**
            * ``Atom(symbol, position)``
            * ``Atom(symbol, x, y, z)``


        :Parameters:

        symbol : str
            The atomic symbol for the atom object
        position : Vector
            The atom's cartesian position
        x : float
            The atom's x posistion
        y : float
            The atom's y posistion
        z : float
            The atom's z posistion



        :Examples:

        >>> Atom('H', [0,0,0])
        Atom('H', [  0.00000000,  0.00000000,  0.00000000 ] )
        >>> Atom(position=(1.5, 0., 0.), symbol='H')
        Atom('H', [  1.50000000,  0.00000000,  0.00000000 ] )
        >>> Atom([1.5, 1.5, 0.], symbol='H')
        Atom('H', [  1.50000000,  1.50000000,  0.00000000 ] )

        """
        raise OverloadedFunctionCallError

    @__init__.overload_with(
        symbol=basestring,
        position=IterableOf(Real))
    def __init__(self, symbol, position, **kwargs):
        self._init_common(**kwargs)
        self.symbol = symbol
        if hasunits(position):
            if 'units' not in kwargs:
                self._cartesian_units = position.units
            elif self._cartesian_units != kwargs['units']:
                position = position.in_units(self._cartesian_units)
        self.position = Vector(position, copy=True, units=self.cartesian_units)

    @__init__.overload_with(
        symbol=basestring,
        x=Real, y=Real, z=Real)
    def __init__(self, symbol, x, y, z, **kwargs):
        self._init_common(**kwargs)
        self.symbol = symbol
        self.position = Vector([x, y, z], units=self.cartesian_units)

    @__init__.overload_with(
        element=Element,
        isotope=Isotope,
        position=Vector
    )
    def __init__(self, element, isotope, position, **kwargs):
        self._init_common(**kwargs)
        self._element = element
        self.isotope = isotope
        self.position = position

    def _init_common(self, **kwargs):
        self.position = None
        self.parent_molecule = None
        self._cartesian_units = kwargs.pop('units', DistanceUnit.default)
        if type_checking_enabled:
            if not isunit(self.cartesian_units) or self.cartesian_units.genre is not DistanceUnit:
                raise ValueError("Invalid units for Atom constructor: {}".format(self.cartesian_units))

    ###################
    # Special Methods #
    ###################

    def __copy__(self):
        """ Copy the atom.  Calls `deepcopy` (i.e. there's only one way to copy an Atom instance).
        """
        return deepcopy(self)

    def __deepcopy__(self, memo):
        """ Override deepcopy to prevent copying of the element object
        """
        # -> *Don't* copy the parent_molecule attribute.  Rationale: this will most often be called by the
        #       Molecule class's __deepcopy__, and thus the atom will be assigned to a new molecule
        dup = Atom(self.symbol, deepcopy(self.position, memo), units=self.cartesian_units)
        dup._isotope = self._isotope
        dup.zmat_label = self.zmat_label
        return dup

    def __eq__(self, other):
        if (self.symbol, self.isotope) != (other.symbol, other.isotope):
            return False
        elif abs(self.x - other.x) > self.same_coordinate_tolerance:
            return False
        elif abs(self.y - other.y) > self.same_coordinate_tolerance:
            return False
        elif abs(self.z - other.z) > self.same_coordinate_tolerance:
            return False
        return True

    def __lt__(self, other):
        """
        a < b as Atoms if a.symbol < b.symbol or (if a.symbol == b.symbol), a.isotope < b.isotope, or
        (if a.isotope == b.isotope) if a.x < b.x, or (if a.x == a.x), if a.y < b.y, or (if a.y == a.y),
        if a.z < b.z.
        This establishes a well-ordering of atoms in a molecule to make comparison of molecule objects easier.
        """
        if (self.symbol, self.isotope) < (other.symbol, other.isotope):
            return True
        elif (self.symbol, self.isotope) > (other.symbol, other.isotope):
            return False
        elif abs(self.x - other.x) > self.same_coordinate_tolerance:
            return self.x < other.x
        elif abs(self.y - other.y) > self.same_coordinate_tolerance:
            return self.y < other.y
        elif abs(self.z - other.z) > self.same_coordinate_tolerance:
            return self.z < other.z
        else:
            return False

    #------------------------#
    # Output Representations #
    #------------------------#

    def __repr__(self):
        return "Atom('" + self.symbol + "', [" + ','.join(["%12.8f" % num for num in self.position]) + " ] )"

    def __str__(self):
        symb_iso = self.symbol
        if not self.isotope.is_principal():
            symb_iso = self.isotope.symbol
        if self.zmat_label:
            ret_val = "Atom '{}'".format(self.zmat_label)
        elif not self.is_orphaned():

            ret_val = "Atom #{0}: '{1}'".format(self.index + 1, self.symbol)
        else:
            ret_val = "orphaned Atom '{0}'".format(self.symbol)
        ret_val += ' with position {0}'.format(self.pos)

        return ret_val

    ################
    #  Properties  #
    ################

    @property
    def symbol(self):
        """ The atomic symbol from the periodic table (as a :py:class:`~__builtin__.str` object)
        , with correct capitalization.  Updating this property transparently updates the
        :py:class:`~grendel.util.element.Element` object associated with ``self``.
        """
        return self.element.symbol

    @symbol.setter
    def symbol(self, symb):
        old = self._element
        try:
            self._element = Elements[symb]
        except KeyError:
            raise StandardError("Unknown element symbol " + repr(symb) + ".  Please add element information to grendel/util/ElementData.py")

        if old is None or not old == self._element:
            self.isotope = self._element.principal_isotope

    @CachedProperty
    def index(self):
        """ The index of the atom in its parent molecule.
        If the atom is orphaned but has a `base_atom`, that `Atom`'s `index` will be returned (and,
        of course, if `base_atom` is orphaned, proceed recursively until either the atom has no
        `base_atom` or a non-orphaned `Atom` is found).

        .. note::
           This is a cached property.  Be sure and reset it (by setting atom._index to None) if you
           reorder the atoms in the parent molecule (which is a disasterous thing to do for many other
           parts of the program as well.  If you really need to reorder atoms in a `Molecule`, you
           should create a new `Molecule` instance and copy the atoms to the new instance.)

        """
        if self.is_orphaned():
            if self.base_atom is not None:
                return self.base_atom.index
            else:
                # None will not be cached...
                return None
        else:
            return self.parent_molecule.atoms.index(self)

    @property
    def element(self):
        """ The :py:class:`~grendel.util.element.Element` object corresponding to `self`.
        Updated transparently using the :py:obj:`~grendel.atom.Atom.symbol` property
        """
        return self._element

    @property
    def mass(self):
        """ The atom's mass as a `float` (in AMU).  Automatically retrieved from the `element` and `isotope` attributes.
        """
        if self.isotope is None: return None
        return self.isotope.mass

    @property
    def isotope(self):
        """ The isotope of the element that `self` is composed of.  Defaults to `element.principal_isotope`
        """
        return self._isotope

    @isotope.setter
    def isotope(self, value):
        if self._element is None: raise StandardError("self.element is not set.  Programmer needs more coffee...")
        if value not in self._element.isotopes:
            raise StandardError("Unknown isotope " + str(self._isotope) + " of element " + str(self._element))
        self._isotope = value

    @CachedProperty
    def full_symbol(self):
        if self.isotope.is_principal():
            return self.symbol
        else:
            return self.symbol + str(int(round(self.isotope.mass)))

    @property
    def x(self):
        return self.position[0]

    @x.setter
    def x(self, val):
        self.position[0] = val

    @property
    def y(self):
        return self.position[1]

    @y.setter
    def y(self, val):
        self.position[1] = val

    @property
    def z(self):
        return self.position[2]

    @z.setter
    def z(self, val):
        self.position[2] = val

    @property
    def xyz(self):
        """ Alias for `self.position`.  `pos` is also an alias for `self.position`
        """
        return self.position

    @xyz.setter
    def xyz(self, val):
        self.position = val
    pos = xyz

    @property
    def cart_coords(self):
        return [coord for coord in self.iter_cart_coords()]

    @property
    def cart_indices(self):
        return [idx for coord, idx in self.iter_cart_coords(True)]

    @property
    def parent(self):
        """
        Alias for the `parent_molecule` attribute.
        """
        return self.parent_molecule

    @parent.setter
    def parent(self, new_val):
        self.parent_molecule = new_val

    ###########
    # Methods #
    ###########

    #-------------------------------------#
    # Inquiry methods (which return bool) #
    #-------------------------------------#

    def is_orphaned(self):
        return self.parent_molecule is None

    def is_ghost(self):
        return self._element == Elements['X']

    def is_bonded_to(self, other_atom, n_vdw_radii=1.2, default_vdw_radius=2.0*Angstroms):
        #TODO deal with units
        dist = magnitude(self.xyz - other_atom.xyz)
        vdw_me = self.element.vdw_radius
        if vdw_me is None:
            vdw_me = default_vdw_radius
        vdw_me = strip_units(vdw_me, convert_to=self.cartesian_units, assume_units=DistanceUnit.default)
        vdw_other = self.element.vdw_radius
        if vdw_other is None:
            vdw_other = default_vdw_radius
        vdw_other = strip_units(vdw_other, convert_to=self.cartesian_units, assume_units=DistanceUnit.default)
        return dist < (vdw_me * n_vdw_radii/2 + vdw_other * n_vdw_radii/2)

    #---------------#
    # Other methods #
    #---------------#

    @typechecked(vect=LightVector)
    def displace(self, vect):
        if sanity_checking_enabled:
            if len(vect) != 3:
                raise ValueError("displacement vector must have three dimensions")
        self.position += vect

    def displaced(self, disp_vect):
        """
        Creates an orphaned copy of `atom` and displaces it by disp_vect
        """
        ret_val = deepcopy(self)
        ret_val.displace(disp_vect)
        ret_val.base_atom = self
        return ret_val

    def convert_units(self, new_units):
        self.position = self.position * self._cartesian_units.to(new_units)
        self._cartesian_units = new_units

    def in_units(self, new_units):
        ret_val = deepcopy(self)
        ret_val.convert_units(new_units)
        return ret_val

    def iter_cart_coords(self, with_index=False):
        if not self.is_orphaned():
            cart_rep = self.parent_molecule.cartesian_representation
            idx = self.index
            for x in [0,1,2]:
                if with_index:
                    yield cart_rep[3*idx + x], 3*idx + x
                else:
                    yield cart_rep[3*idx + x]
        else:
            raise ValueError("cannot iterate over cartesian coordinates of an orphaned atom.")

#####################
# Dependent Imports #
#####################

from grendel.chemistry.element_data import Elements

