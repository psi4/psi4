from abc import ABCMeta, abstractmethod
from copy import copy
from functools import partial
from itertools import chain, product, combinations_with_replacement, permutations
from grendel import sanity_checking_enabled, type_checking_enabled
from grendel.coordinates.internal_coordinate import InternalCoordinate, CoordinateDependencyError
from grendel.differentiation.derivative_collection import DerivativeCollection
from grendel.differentiation.finite_difference import Differentiable, FloatShell, FiniteDifferenceDerivative, FiniteDifferenceFunction
from grendel.gmath import DeclareIndexRange
from grendel.gmath.matrix import Matrix
from grendel.gmath.tensor import Tensor
from grendel.gmath.vector import LightVector, Vector
from grendel.util.compatibility import abstractclassmethod
from grendel.util.containers import LRUDict
from grendel.util.decorators import CachedProperty, typechecked
from grendel.util.freezing import SetOnceAttribute
from grendel.util.iteration import flattened, grouper, first
from grendel.util.overloading import pop_multikwarg
from grendel.util.sentinal_values import All
from grendel.util.strings import classname, short_str, andjoin
from grendel.util.units import strip_units, Angstroms, Radians, AngularUnit, ValueWithUnits
from grendel.util.units.unit import AngularUnit

__all__ = [
    "SimpleInternalCoordinate"
]


class SimpleInternalCoordinate(InternalCoordinate):
    """
    Superclass for all types of simple (non-symmetry, non-normal) internal coordinates.

    :Attributes:

    atoms : list of `Atom`
        `list` of `Atom` objects that the coordinate refers to


    """
    __metaclass__ = ABCMeta

    ####################
    # Class Attributes #
    ####################

    # Deltas for computation of B vectors by finite difference
    xyz_delta = 0.001

    analytic_b_orders = [1]

    created_coords = {}
    created_coord_maps = {}

    ######################
    # Private Attributes #
    ######################

    _atoms = None
    _cart_to_idx = None

    ##################
    # Initialization #
    ##################

    def _init(self, **kwargs):
        self._bvector = None
        super(SimpleInternalCoordinate, self)._init(**kwargs)

    ###################
    # Special Methods #
    ###################

    def __copy_kwargs__(self):
        ret_val = super(SimpleInternalCoordinate, self).__copy_kwargs__()
        ret_val.update(
            atoms=copy(self.atoms)
        )
        return ret_val

    #------------------------#
    # Output Representations #
    #------------------------#

    def __format__(self, formspec):
        return format(str(self), formspec)

    def __str__(self):
        """
        """
        tyname = self.__class__.__name__
        if self.is_orphaned():
            return "(orphaned) {} for atoms {}".format(tyname, andjoin([str(atom) for atom in self.atoms]))
        else:
            mol_short_name = short_str(self.molecule)
            return "{} for atom indices {} of {}".format(tyname, andjoin(self.atom_indices), mol_short_name)

    def __short_str__(self):
        tyname = self.__class__.__name__
        if not self.is_orphaned():
            return "{} for atom indices {}".format(tyname, andjoin(self.atom_indices))
        else:
            return "(orphaned) {} with atoms {}".format(tyname, andjoin(self.atom_indices))

    __repr__ = __str__ # for now...

    ##############
    # Properties #
    ##############

    @property
    def value(self):
        if self.units.genre is AngularUnit:
            return self.value_for_positions(*[a.pos for a in self.atoms]) * Radians.to(self.units)
        else:
            return self.value_for_positions(*[a.pos for a in self.atoms]) * self.molecule.cartesian_units.to(self.units)

    @property
    def xyz(self):
        return Vector([atom.pos for atom in self.atoms])


    ######################
    # Private Properties #
    ######################

    @property
    def _l_xyz(self):
        return LightVector([atom.pos for atom in self.atoms])

    ###################
    # Special Methods #
    ###################

    def __add__(self, other):
        if isinstance(other, InternalCoordinate):
            return SumOfInternalCoordinates(
                [self, other],
                coefficients = [1.0, 1.0]
            )
        else:
            return NotImplemented
    __radd__ = __add__

    def __sub__(self, other):
        if isinstance(other, InternalCoordinate):
            return SumOfInternalCoordinates(
                [self, other],
                coefficients = [1.0, -1.0]
            )
        else:
            return NotImplemented

    def __rsub__(self, other):
        if isinstance(other, InternalCoordinate):
            return SumOfInternalCoordinates(
                [other, self],
                coefficients = [1.0, -1.0]
            )
        else:
            return NotImplemented

    # TODO __mul__ for giving coefficients in function coordinates

    ##################
    # Static Methods #
    ##################

    @staticmethod
    def b_tensor_element_reindexed(current_coord, *coords):
        coords = tuple(coords) if not isinstance(coords[0], tuple) else coords[0]
        internal = coords[0]
        carts = coords[1:]
        #========================================#
        if internal.is_orphaned():
            #----------------------------------------#
            # avoid the gratuitous creation of coordinates
            key = (type(internal), tuple(internal.atoms))
            internal = SimpleInternalCoordinate.created_coords.get(key, None)
            if internal is None:
                SimpleInternalCoordinate.created_coords[key] = internal = coords[0]
                # TODO UNITS FOR CREATED COORDINATES!!!
            #----------------------------------------#
            # Check to make sure their all integers...
            if sanity_checking_enabled:
                if any(isinstance(cart, CartesianCoordinate) for cart in carts):
                    cart = first(carts, lambda c: isinstance(c, CartesianCoordinate))
                    raise ValueError("ambiguous B tensor element retrieval for orphaned internal"
                                     " coordinate '{}' using cartesian coordinate '{}'".format(
                        short_str(internal), short_str(cart)
                    ))
                if not all(isinstance(cart, int) for cart in carts):
                    raise ValueError(
                        "b_tensor elements may not be retrieved for "
                        " orphaned coordinates by anything but integer indices"
                    )
            #----------------------------------------#
            # get the offset indices
            # what we have now is integers in current_coord's internal indexing scheme
            new_idxs = current_coord.molecule_indices_for(carts)
            # This method uses the molecular indexing scheme.  This is why we needed to
            #   do the above index translation.
            return internal.b_tensor_element(*new_idxs)
        #========================================#
        else: # internal is not orphaned
            #----------------------------------------#
            # Tidy things up and convert internal indices
            #   to molecular indices
            if sanity_checking_enabled:
                if any(isinstance(cart, CartesianCoordinate) for cart in carts) \
                        and not all(isinstance(cart, CartesianCoordinate) for cart in carts):
                    raise ValueError("mixing of CartesianCoordinates and integers in b tensor"
                                     " element retrieval is confusing and thus no longer allowed."
                                     " coordinates were ('{}')".format(
                        "', '".join(str(c) for c in carts)
                    ))
            if all(isinstance(c, int) for c in carts):
                # These are integers in current_coord's internal scheme.  Translate them to
                #   the molecular indexing scheme.
                carts = current_coord.molecule_indices_for(carts)
            elif type_checking_enabled:
                if not all(isinstance(cart, CartesianCoordinate) for cart in carts):
                    raise TypeError("arguments beyond the first two to b_tensor_element_reindexed"
                                    " must be either all ints or all CartesianCoordinates."
                                    " coordinates were ('{}')".format(
                        "', '".join(str(c) for c in carts)
                    ))
            #----------------------------------------#
            # This method uses the molecular indexing scheme.  This is why we needed to
            #   do the above index translation.
            return internal.b_tensor_element(*carts)

    ##########################
    # Abstract Class Methods #
    ##########################

    @abstractclassmethod
    def value_for_xyz(cls, xyz):
        """ Value for the xyz matrix constructed from the coordinates of the atoms for the coordinate.
        For instance, for a
        """
        return NotImplemented

    @abstractclassmethod
    def b_vector_for_positions(cls, *args):
        """ Returns the b vector corresponding to a coordinate with atoms in the positions
        given as arguments (assumed to be `Vector` or at least `LightVector` instances).
        Note that the returned value uses the indexing scheme corresponding to the order of
        the argument's passed in, not the indexing scheme corresponding to a `Molecule` instance
        somewhere.
        """
        return NotImplemented

    #################
    # Class Methods #
    #################

    @classmethod
    def value_for_positions(cls, *posvects):
        """ Get the value of the coordinate position vectors xyz.

        .. note::
           For angular coordinates, this *always* returns a value in Radians.  `value`, `value_with_units`, and
           `value_for_molecule`, however, return values in `self.units`

        """
        return cls.value_for_xyz(posvects)

    ###########
    # Methods #
    ###########

    def atom_num(self, i):
        """
        Returns the number of the atom (as an index in `self.molecule`) that the i'th atom composing the
        coordinate `self` refers to.
        """
        return self.molecule.index(self.atoms[i])

    def get_coord(self, ctype, *atoms):
        if not self.is_orphaned():
            self.parent_representation._gen_icoord_map()
            key = (ctype, atoms)
            ret_val = self.parent_representation._icoord_map.get(key, None)
            if ret_val is None:
                return ctype(*atoms)
            else:
                return ret_val
        else:
            return ctype(*atoms)

    #--------------------------------#
    # Methods abstract in Coordinate #
    #--------------------------------#

    def copy_for_representation(self, rep, **kwargs):
        copykw = self.__copy_kwargs__()
        copykw.update(kwargs)
        copykw.update(parent=rep)
        copykw.update(atoms=[rep.molecule[i] for i in self.atom_indices])
        return self.__class__(**copykw)

    def value_for_molecule_matrix(self, mat):
        i = self.atom_indices
        new_mat = Matrix([mat[i[j]] for j in range(len(self.atoms))])
        return self.__class__.value_for_xyz(new_mat)

    #----------------------------------------#
    # Methods abstract in InternalCoordinate #
    #----------------------------------------#

    def copy_with_atoms(self, atoms, **kwargs):
        copykw = self.__copy_kwargs__()
        copykw.update(kwargs)
        if 'units' not in kwargs:
            copykw.update(units=self.units)
        if 'base_analog' not in kwargs:
            copykw.update(base_analog=self)
        copykw.update(atoms=atoms)
        return self.__class__(**copykw)


    ##########################
    # Private Static Methods #
    ##########################

    @staticmethod
    def _check_b_tensor_cache(clstype, order, *idx_n_vects):
        try:
            cache_val = SimpleInternalCoordinate.b_tensor_cache[
                   (clstype, order) + tuple(tuple(v) for v in idx_n_vects)]
            return cache_val
        except KeyError:
            return None

    @staticmethod
    def _set_b_tensor_cache_entry(value, clstype, order, *idx_n_vects):
        SimpleInternalCoordinate.b_tensor_cache[
            (clstype, order) + tuple(tuple(v) for v in idx_n_vects)] = value
        return None

    ###################
    # Private Methods #
    ###################

    def _get_state(self):
        # for all simple internal coordinates, the B vector depends only on the positions of the atoms...
        ret_val = []
        for atom in self.atoms:
            ret_val.extend(atom.pos)
        return tuple(ret_val)





#####################
# Dependent Imports #
#####################

from grendel.chemistry.molecule import Molecule

from grendel.representations.cartesian_representation import X, Y, Z

from grendel.coordinates.cartesian_coordinate import CartesianCoordinate

from grendel.differentiation.derivative_tensor import  DerivativeTensor

from grendel.coordinates.function_coordinate import SumOfInternalCoordinates

# needed for dynamic typechecking; do not delete
if type_checking_enabled:
    from grendel.representations.internal_representation import InternalRepresentation
    from grendel.representations.cartesian_representation import CartesianRepresentation

