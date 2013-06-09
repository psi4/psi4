from abc import abstractmethod, ABCMeta, abstractproperty
from copy import deepcopy
from functools import partial
from itertools import product, permutations, combinations_with_replacement as symmetric_product
from grendel import sanity_checking_enabled
import math
import numpy as np
from grendel.chemistry.molecule import Molecule
from grendel.coordinates.cartesian_coordinate import CartesianCoordinate
from grendel.coordinates.coordinate import Coordinate, OrphanedCoordinateError, CoordinateDependencyError
from grendel.differentiation.finite_difference import FiniteDifferenceFunction, FloatShell, FiniteDifferenceDerivative
from grendel.gmath.matrix import Matrix
from grendel.gmath.tensor import Tensor
from grendel.gmath.vector import LightVector, Vector
from grendel.util.containers import LRUDict
from grendel.util.decorators import typechecked
from grendel.util.exceptions import ProgrammerNeedsMoreCoffeeError
from grendel.util.iteration import grouper, flattened
from grendel.util.overloading import pop_multikwarg
from grendel.util.sentinal_values import All
from grendel.util.strings import andjoin
from grendel.util.units import DistanceUnit, AngularUnit, Radians, Angstroms, isunit, ValueWithUnits, strip_units
from grendel.util.units.value_with_units import hasunits

__all__ = [
    "InternalCoordinate"
]


class InternalCoordinate(Coordinate, FiniteDifferenceFunction):
    """ General superclass for all types internal coordinates.
    Internal coordinates may only refer to one representation.
    """
    __metaclass__ = ABCMeta

    ####################
    # Class Attributes #
    ####################

    max_b_tensor_cache_entries = int(1e7)
    b_tensor_cache = LRUDict(max_b_tensor_cache_entries)
    b_tensor_finite_difference_rigor = 10
    b_tensor_finite_difference_delta = 0.001 * Angstroms
    analytic_b_orders = [1]

    ############################
    # Private Class Attributes #
    ############################

    _b_finite_difference_managers = {}

    ######################
    # Private Attributes #
    ######################

    _atoms = abstractproperty()
    _bvector = None
    _btensor = None
    _state = None
    _created_cart_coords = None

    ##############
    # Properties #
    ##############

    @property
    def atoms(self):
        return self._atoms

    @property
    def b_vector(self):
        if self._bvector is None or self._state != self._get_state():
            self._recompute_b_vector()
        return self._bvector
    # Aliases
    bvector = b_vector

    @property
    def cart_coords(self):
        if self.is_orphaned():
            raise OrphanedCoordinateError("'cart_coords' property is only valid for non-orphaned"
                                          " coordinates.  Coordinate '{}' is orphaned.".format(self))
        else:
            return [coord for coord in self.iter_cart_coords()]

    @property
    def cart_indices(self):
        return [idx for coord, idx in self.iter_cart_coords(True)]

    #-------------------------------------------------------------------#
    # Properties needed to conform to FiniteDifferenceFunction protocol #
    #-------------------------------------------------------------------#

    @property
    def variables(self):
        return [f for f in flattened(
            [[self.molecule.cartesian_representation[3*i + d] for d in [0,1,2]] for i in self.atom_indices]
        )]

    ####################
    # Abstract Methods #
    ####################

    @abstractmethod
    def copy_with_atoms(self, atoms):
        return NotImplemented

    @abstractmethod
    def b_vector_for_positions(self, *pos):
        return NotImplemented

    ###########
    # Methods #
    ###########

    def internal_indices_for_coordinates(self, *cart_coords):
        mol_idxs = []
        for coord in cart_coords:
            if isinstance(coord, CartesianCoordinate):
                if not coord.is_orphaned():
                    mol_idxs.append(coord.index)
                else:
                    raise ValueError("can't get internal index on coordinate '{}' corresponding"
                                     " to orphaned CartesianCoordinate '{}'".format(
                        self,
                        coord
                    ))
            else:
                raise TypeError("all arguments to 'internal_indices_for_coordinates' must be"
                                " CartesianCoordinate instances (got at least one that"
                                " was a '{}')".format(type(coord).__name__))
        try:
            ret_val = self.internal_indices_for(mol_idxs)
        except CoordinateDependencyError:
            for coord in cart_coords:
                try:
                    if isinstance(coord, CartesianCoordinate):
                        self.internal_indices_for([coord.index])
                    else:
                        self.internal_indices_for(coord)
                except CoordinateDependencyError:
                    if isinstance(coord, CartesianCoordinate):
                        raise CoordinateDependencyError("coordinate '{}' does not depend on"
                                                        " CartesianCoordinate '{}'".format(
                            self, coord
                        ))
                    else:
                        # reraise the original error
                        raise
            # if we've gotten here, something has gone terribly wrong.  Just reraise the error we got.
            raise
        return tuple(ret_val)

    def generate_name(self, one_based=True):
        off = 1 if one_based else 0
        if not self.is_orphaned() or all(not a.is_orphaned() for a in self.atoms):
            return self.coordinate_symbol + "".join(str(a.index+off) for a in self.atoms)
        else:
            return self.coordinate_symbol + "(o)"

    #--------------------------#
    # B tensor related methods #
    #--------------------------#

    def analytic_b_tensor_for_order(self, order):
        """
        Compute the analytic B tensor of a given order and return a tensor indexed by the
        coordinate's own indexing scheme.
        """
        if order == 1:
            # We can't use the `b_vector` attribute since that vector is indexed in the parent
            #   molecule's indexing scheme
            my_b = [self.__class__.b_vector_for_positions(*[a.pos for a in self.atoms])]
            if NotImplemented in my_b:
                return NotImplemented
            # Return now, to avoid double-caching
            return Vector(my_b)
        else:
            return NotImplemented

    def get_b_tensor(self, max_order):
        """
        Get the coordinate's b tensor as a `Coordinate`-based `DerivativeCollection`.  If it has not
        already done so, the coordinate will fill in the b tensor up to `max_order`.  The returned
        value uses the `Coordinate`'s internal indexing scheme.
        """
        if self._btensor is None:
            self._init_btensor()
        for order in xrange(1, max_order+1):
            if self._btensor.has_order(order):
                continue
            value = self.fill_b_tensor_analytically(order=order)
            if value is NotImplemented:
                self._btensor[...] = self.finite_difference_b_tensor(
                    order=order,
                    use_parent_indices=False)
        return self._btensor

    def b_tensor_element(self, *cart_coords_or_indices):
        """
        Retrieve an element of the Coordinate's B tensor using CartesianCoordinate instances
        and/or integer indices *in the molecule's indexing scheme*.
        """
        if self._btensor is None:
            self._init_btensor()
        if all(isinstance(c, CartesianCoordinate) for c in cart_coords_or_indices):
            idxs = cart_coords_or_indices
        elif all(isinstance(c, int) for c in cart_coords_or_indices):
            try:
                idxs = self.internal_indices_for(cart_coords_or_indices)
            except CoordinateDependencyError:
                #TODO raise an error here instead?
                raise
                #return 0.0
        else:
            raise TypeError(
                "arguments to b_tensor_element must be all integers (in the molecule's"
                " indexing scheme) or all CartesianCoordinate instances."
                " Mixing of CartesianCoordinates and integers in b_tensor_element"
                " arguments is confusing and thus no longer allowed."
                " indices were ('{}')".format(
                    "', '".join(str(c) for c in cart_coords_or_indices)
                )
            )
        if not self.is_orphaned():
            try:
                if self._btensor.has_order(len(cart_coords_or_indices)):
                    return self._btensor[idxs]
                else:
                    self.fill_b_tensor_analytically(order=len(cart_coords_or_indices))
                    return self._btensor[idxs]
            except CoordinateDependencyError:
                #TODO raise an error here instead?
                raise
                #return 0.0
        else:
            if sanity_checking_enabled:
                if not all(isinstance(cart, int) for cart in cart_coords_or_indices):
                    raise ValueError(
                        "b_tensor elements may not be retrieved for orphaned coordinates by"
                        " anything but integer indices (got indices of types {})".format(
                            andjoin("'" + type(item).__name__ + "'" for item in cart_coords_or_indices)
                        )
                    )
            order = len(idxs)
            if not self._btensor.has_order(order):
                response = self.fill_b_tensor_analytically(order=order)
                if response is NotImplemented:
                    self._btensor[...] = self.finite_difference_b_tensor(order=order)
            try:
                return self._btensor.for_order(order)[idxs]
            except IndexError:
                raise

    @typechecked(
        B=('DerivativeCollection', None),
        order=int)
    def fill_b_tensor_analytically(self, B=None, order=None):
        """
        Fill the coordinate's own `_btensor` attribute for a given order.  (The `_btensor`
        attribute should be retrieved using the `get_b_tensor()` method, which calls this
        if the _btensor has not already been computed.)  If the b tensor cannot be computed
        analytically at the given order, this method returns `NotImplemented`.  Otherwise,
        it computes the B tensor to the given order and returns `None`.  If the optional `B`
        argument is given, this Coordinate's part of the `Representation`-based `TensorCollection`
        that `B` refers to is filled as well.
        """
        if order is None:
            raise TypeError
        if self._btensor is None:
            self._init_btensor()
        if not self._btensor.has_order(order):
            fill_val = self.analytic_b_tensor_for_order(order)
            if fill_val is not NotImplemented:
                self._btensor[...] = fill_val
            else:
                return NotImplemented
        if B is not None:
            # Fill in the parts of the collection corresponding to the atoms in self
            #   This unwrapping needs to be done because within the coordinate itself,
            #   the atoms are indexed by number within the coordinate, whereas in the
            #   representation, the atoms are indexed by number within the molecule.
            #   This is somewhat confusing, but necessary to allow "detached coordinates"
            #   to have B tensors associated with them (which is necessary, for instance,
            #   if a BondAngle(a, b, c) coordinate is given with one of BondLength(a,b) or
            #   BondLength(b, c) not being a coordinate in the parent representation)
            bfill = Tensor(shape=(3*B.representation.molecule.natoms,) * order)
            for grp in product(enumerate(self.iter_molecule_indices()), repeat=order):
                # 'Unzip' the index groups
                my_idxs, mol_idxs = zip(*grp)
                bfill[mol_idxs] = self._btensor[my_idxs]
            B.for_order(order)[self] = bfill
        # Return nothing, signaling success
        return

    @typechecked(order=int)
    def finite_difference_b_tensor(
            self,
            order,
            using_order=None,
            robustness=None,
            forward=False,
            use_parent_indices=True):
        """ Computes the B tensor for order `order` by finite difference
        of B tensors of order `order-1`, taking into account permutational symmetry.
        If `use_parent_indices` is `False`, a `Tensor` is returned that is indexed
        by the coordinate's internal indices rather than the parent molecule's atom indices.
        Otherwise, a `DerivativeTensor` is returned with its `representation` attribute set to the
        coordinate's parent molecule's `cartesian_representation` attribute.
        """
        if using_order is None:
            using_order = order - 1
        else:
            raise NotImplementedError
        robustness = robustness or self.b_tensor_finite_difference_rigor
        if use_parent_indices:
            ret_val = DerivativeTensor(
                representation=self.molecule.cartesian_representation,
                shape=(3*self.molecule.natoms,) * order
            )
        else:
            ret_val = Tensor(
                shape=(len(self.atoms) * 3,) * order
            )
        for coords_and_indices in symmetric_product(self.iter_cart_coords(with_index=True), order):
            coords, indices = zip(*coords_and_indices)
            val_to_spread = self.finite_difference_b_tensor_element(
                *coords,
                robustness=robustness,
                forward=forward
            )
            if use_parent_indices:
                for perm in permutations(indices):
                    ret_val[tuple(perm)] = val_to_spread
            else:
                for perm in permutations(coords):
                    perm_indices = self.internal_indices_for_coordinates(*perm)
                    ret_val[perm_indices] = val_to_spread
        return ret_val

    def finite_difference_b_tensor_element(self, *cart_coords, **kwargs):
        rigor = pop_multikwarg(kwargs, 'rigor', 'robustness') or self.b_tensor_finite_difference_rigor
        forward = pop_multikwarg(kwargs, 'forward', 'use_forward') or False
        key = (self, self.molecule.cartesian_units, rigor, forward)
        manager = InternalCoordinate._b_finite_difference_managers.get(key, None)
        if manager is None:
            manager = _BTensorFiniteDifferenceManager(
                robustness=rigor,
                units=self.molecule.cartesian_units,
                forward=forward
            )
            InternalCoordinate._b_finite_difference_managers[key] = manager
        return manager.get_b_tensor_entry(self, cart_coords)

    #-------------------#
    # Iteration methods #
    #-------------------#

    def iter_cart_coords(self, with_index=False):
        cart_rep = self.molecule.cartesian_representation
        for i in self.atom_indices:
            for x in [X, Y, Z]:
                if with_index:
                    yield cart_rep[3*i + x], 3*i + x
                else:
                    yield cart_rep[3*i + x]

    #----------------------------------------------#
    # Methods abstract in FiniteDifferenceFunction #
    #----------------------------------------------#

    def value_for_displacements(self, pairs):
        # get the position vector
        xyzvect = LightVector([a.position for a in self.atoms]).ravel()
        # make pairs into a list in case it's a general iterator
        pairs = list(pairs)
        int_idxs = self.internal_indices_for_coordinates(*[pair[0] for pair in pairs])
        for i, (__, ndeltas) in enumerate(pairs):
            xyzvect[int_idxs[i]] += self.b_tensor_finite_difference_delta * ndeltas
        return self.__class__.value_for_xyz(
            Matrix(
                list(grouper(3, xyzvect))
            )
        )

    def deltas_for_variables(self, vars):
        return (self.b_tensor_finite_difference_delta,) * len(set(vars))

    ############################
    # Private Abstract Methods #
    ############################

    #TODO benchmark caching of state versus simply recomputing every time (particularly for torsion angles)
    #TODO generalize and/or remove this (it is a bit redundant now that there is a B tensor cache)
    @abstractmethod
    def _get_state(self):
        """ Get any data that, if changed, would invalidate the b vector calculation.
        The B vector will be recomputed if self._get_state() != self._state
        """
        return NotImplemented

    # TODO B vector computing functions for all internal coordinate subclasses
    @abstractmethod
    def _recompute_b_vector(self):
        """Pure virtual method which must be implemented in subclasses.
        """
        # pragma: no cover
        return NotImplemented

    ###################
    # Private Methods #
    ###################

    def _init_btensor(self):
        if hasattr(self.__class__, '_btens_idx_range_set'):
            self._btensor = DerivativeCollection(
                coordinate=self,
                uncomputed=False,
                einsum_index='v',
                index_range_set=self.__class__._btens_idx_range_set
            )
        else:
            self._btensor = DerivativeCollection(
                coordinate=self,
                uncomputed=False
            )


    def _cache_state(self):
        """ Cache copies of any data that, if changed, would invalidate the b vector calculation
        """
        # do a deepcopy just to be on the safe side
        self._state = deepcopy(self._get_state())

class _BTensorFiniteDifferenceManager(object):

    ##############
    # Attributes #
    ##############

    cartesian_representation = None
    robustness = None
    delta = None
    forward = None
    displaced_molecules = None
    displaced_coordinates = None
    b_tensor_collection = None
    fdiff_instances = None

    ##################
    # Initialization #
    ##################

    @typechecked(
        robustness=int,
        delta=(float, ValueWithUnits),
        forward=bool,
        units=(None, isunit),
    )
    def __init__(self,
            robustness=InternalCoordinate.b_tensor_finite_difference_rigor,
            delta=InternalCoordinate.b_tensor_finite_difference_delta,
            units=None,
            forward=False):
        self.robustness = robustness
        self.delta = delta
        if self.cartesian_representation is not None:
            self.delta = strip_units(self.delta, units)
        self.forward = forward
        self.displaced_molecules = {}
        self.displaced_coordinates = {}
        self.fdiff_instances = {}

    ###########
    # Methods #
    ###########

    def get_b_tensor_entry(self, internal_coord, cart_coords):
        fdiff_key = (internal_coord, cart_coords)
        fdiff = self.fdiff_instances.get(fdiff_key, None)
        if fdiff is None:
            fdiff = FiniteDifferenceDerivative(
                #--------------#
                # function =
                internal_coord,
                #--------------#
                # *variables =
                cart_coords[0],
                #--------------#
                value_function=partial(
                    self.value_for_coords_and_displacements,
                    internal_coord=internal_coord,
                    cart_coords=cart_coords[1:]),
                forward=self.forward,
                delta=self.delta,
                robustness=self.robustness
            )
            self.fdiff_instances[fdiff_key] = fdiff
        ret_val = fdiff.value
        return ret_val

    def value_for_coords_and_displacements(self, pairs, internal_coord, cart_coords):
        cart_coords = tuple(cart_coords)
        disps = [0] * (len(internal_coord.atoms) * 3)
        for coord, ndeltas in pairs:
            idx = internal_coord.internal_indices_for_coordinates(coord)[0]
            disps[idx] = ndeltas
        coord_idxs = internal_coord.internal_indices_for_coordinates(*cart_coords)
        dispcoord = self.coordinate_with_cartesian_displacements(internal_coord, disps)
        if len(cart_coords) == 0:
            # For some reason, BondAngle thinks it's a subclass of PeriodicCoordinate (!?)
            if hasattr(dispcoord, 'possible_values_for_xyz'):
                baseval = internal_coord.value_for_xyz(
                    Matrix([a.pos for a in internal_coord.atoms])
                )
                pos_vals = dispcoord.possible_values_for_xyz(
                    Matrix([a.pos for a in dispcoord.atoms])
                )
                # Including the half-periods makes this work...
                # TODO fix and/or justify this
                val = min(*(pos_vals + tuple(p+math.pi for p in pos_vals)), key=lambda p: abs(p - baseval))
                # No need to convert angular units, since value_for_xyz always returns Radians
                return val
            else:
                val = dispcoord.value
                if hasunits(internal_coord) and internal_coord.units.genre is AngularUnit:
                    val *= dispcoord.units.to(Radians)
                return val
        else:
            B = dispcoord.get_b_tensor(max_order=len(cart_coords))
            return B[tuple(coord_idxs)]

    def coordinate_with_cartesian_displacements(self, base_coord, disps):
        disps = tuple(disps)
        if all(d == 0 for d in disps):
            return base_coord
        ret_val = self.displaced_coordinates.get((base_coord, disps), None)
        if ret_val is None:
            new_atoms = []
            for atom, disp in zip(base_coord.atoms, grouper(3, disps)):
                new_atoms.append(atom.displaced(LightVector(disp)*self.delta))
            ret_val = base_coord.copy_with_atoms(new_atoms)
            # Cache it for later to avoid gratuitous creation of displaced Molecule instances
            self.displaced_coordinates[(base_coord, disps)] = ret_val
        return ret_val


#####################
# Dependent Imports #
#####################

from grendel.differentiation.derivative_tensor import DerivativeTensor
from grendel.representations.cartesian_representation import X, Z, Y
from grendel.differentiation.derivative_collection import DerivativeCollection

