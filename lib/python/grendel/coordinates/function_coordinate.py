"""
Coordinates that are functions of other coordinates.
Warning:  This module is pretty sloppy right now.  Sorry about that.
"""
from collections import Iterable
from functools import partial
from inspect import ismethod
from grendel import sanity_checking_enabled, type_checking_enabled
from itertools import product
from grendel.coordinates.coordinate import Coordinate
from grendel.coordinates.internal_coordinate import InternalCoordinate
from grendel.gmath import Vector
from grendel.gmath.tensor import Tensor
from grendel.util.iteration import grouper
from grendel.util.sentinal_values import All
from grendel.util.units import strip_units, ValueWithUnits, isunit, Radians


class FunctionCoordinate(Coordinate):
    """
    A coordinate that is a function of one or more other `Coordinate` instances
    """

    ##############
    # Attributes #
    ##############

    subcoordinates = None
    function = None
    atoms = None
    default_delta = None
    index = None

    ######################
    # Private Attributes #
    ######################

    _idxmap = None

    # To make InternalCoordinate happy...
    _atoms = atoms

    ##################
    # Initialization #
    ##################

    def __init__(self,
            function,
            subcoordinates,
            index = None,
            default_delta = None,
            **kwargs):
        """
        ..note:
        If non-orphaned coordinates are passed in, orphaned copies are made.
        """
        self.function = function
        self.index = index
        #----------------------------------------#
        if 'units' in kwargs and kwargs['units'] is not None:
            units = kwargs['units']
        else:
            units = subcoordinates[0].units
            kwargs.update(units=units)
        #----------------------------------------#
        self.subcoordinates = []
        for sc in subcoordinates:
            self.subcoordinates.append(sc.orphaned_copy(units=units))
        #----------------------------------------#
        # Determine the atom ordering/numbering before we start
        # At the same time, make a list of coordinate_atom_index => function_coordinate_atom_index
        #   mappings for each subcoordinate.
        self.atoms = []
        self._idxmap = [{} for _ in xrange(len(self.subcoordinates))]
        for isc, sc in enumerate(self.subcoordinates):
            for iatom, atom in enumerate(sc.atoms):
                if atom not in self.atoms:
                    self.atoms.append(atom)
                    self._idxmap[isc][iatom] = len(self.atoms) - 1
                else:
                    self._idxmap[isc][iatom] = self.atoms.index(atom)
        #----------------------------------------#
        # default_delta is not required, but any
        #   method that requires a displacement will fail
        #   unless the delta is manually specified for that
        #   method.
        self.default_delta = default_delta
        #--------------------------------------------------------------------------------#
        self._init(**kwargs)

    ###################
    # Special Methods #
    ###################

    def __copy_kwargs__(self):
        rv = super(FunctionCoordinate, self).__copy_kwargs__()
        rv.update(
            function=self.function,
            subcoordinates=self.subcoordinates,
            default_delta=self.default_delta
        )
        return rv

    ###########
    # Methods #
    ###########

    #--------------------------------#
    # Methods abstract in Coordinate #
    #--------------------------------#

    def value_for_molecule_matrix(self, mat):
        vals = [sc.value_for_molecule_matrix(mat) for sc in self.subcoordinates]
        return self.function(*vals)

    def value_for_positions(self, *pos):
        vals = []
        for isc, sc in enumerate(self.subcoordinates):
            vals.append(sc.value_for_positions(*[pos[self._idxmap[isc][i]]
                                                 for i in xrange(len(sc.atoms))]))
        return self.function(*vals)

    def copy_for_representation(self, rep, **kwargs):
        copykw = self.__copy_kwargs__()
        copykw.update(kwargs)
        copykw.update(parent=rep)
        new_atoms = [rep.molecule[a.index] for a in self.atoms]
        #new_atoms = [rep.molecule[i[0]//3] for i in grouper(3, self.iter_molecule_indices)]
        new_sc = []
        for isc, sc in enumerate(self.subcoordinates):
            atmap = self._idxmap[isc]
            catoms = [new_atoms[a] for _, a in sorted(atmap.items())]
            new_sc.append(sc.copy_with_atoms(catoms))
        copykw.update(
            subcoordinates=new_sc
        )
        return self.__class__(**copykw)

class InternalFunctionCoordinate(FunctionCoordinate, InternalCoordinate):
    """
    """

    ##############
    # Attributes #
    ##############

    b_vector_same_function = None
    b_vector_function = None
    b_vector_takes_coords = None
    b_tensor_same_function = None
    b_tensor_function = None
    b_tensor_takes_coords = None
    value_takes_coords = None

    ##################
    # Initialization #
    ##################

    def __init__(self,
            b_vector_same_function=True,
            b_vector_function=None,
            b_vector_takes_coords=False,
            b_tensor_same_function=True,
            b_tensor_function=None,
            b_tensor_takes_coords=False,
            value_takes_coords=False,
            **kwargs):
        self.b_vector_same_function = b_vector_same_function
        self.b_vector_function = b_vector_function
        self.b_vector_takes_coords = b_vector_takes_coords
        self.b_tensor_same_function = b_tensor_same_function
        self.b_tensor_function = b_tensor_function
        self.b_tensor_takes_coords = b_tensor_takes_coords
        self.value_takes_coords = value_takes_coords
        if self.value_takes_coords:
            raise NotImplementedError
        super(InternalFunctionCoordinate, self).__init__(**kwargs)
        limiting_coordinates = [sc for sc in self.subcoordinates if sc.analytic_b_orders is not All]
        if len(limiting_coordinates) == 0:
            self.analytic_b_orders = All
        else:
            self.analytic_b_orders = set(limiting_coordinates[0].analytic_b_orders).intersection(*[
                sc.analytic_b_orders for sc in limiting_coordinates
            ])

    def _validate_initialization(self):
        # Check to make sure all of the subcoordinates are InternalCoordinate instances
        if not all(isinstance(c, InternalCoordinate) for c in self.subcoordinates):
            raise TypeError("InternalFunctionCoordinate instances may only be composed of"
                            " InternalCoordinate instances as subcoordinates.")
        super(InternalFunctionCoordinate, self)._validate_initialization()

    ###################
    # Special Methods #
    ###################

    def __copy_kwargs__(self):
        kw = InternalCoordinate.__copy_kwargs__(self)
        # Give precidence to FunctionCoordinate over InternalCoordinate
        kw.update(FunctionCoordinate.__copy_kwargs__(self))
        kw.update(
            b_vector_same_function = self.b_vector_same_function,
            b_vector_function = self.b_vector_function,
            b_vector_takes_coords = self.b_vector_takes_coords,
            b_tensor_same_function = self.b_tensor_same_function,
            b_tensor_function = self.b_tensor_function,
            b_tensor_takes_coords = self.b_tensor_takes_coords,
            value_takes_coords = self.value_takes_coords
        )
        return kw

    ###########
    # Methods #
    ###########

    #----------------------------------------#
    # Methods abstract in InternalCoordinate #
    #----------------------------------------#

    def copy_with_atoms(self, atoms, **kwargs):
        copykw = InternalCoordinate.__copy_kwargs__(self)
        # Give precidence to FunctionCoordinate over InternalCoordinate
        copykw.update(FunctionCoordinate.__copy_kwargs__(self))
        copykw.update(self.__copy_kwargs__())
        copykw['base_analog'] = self
        # allow manual override of any of the above
        copykw.update(kwargs)
        new_sc = []
        for isc, sc in enumerate(self.subcoordinates):
            atmap = self._idxmap[isc]
            catoms = [atoms[a] for _, a in sorted(atmap.items())]
            new_sc.append(sc.copy_with_atoms(catoms))
        copykw.update(subcoordinates=new_sc)
        return self.__class__(**copykw)

    def b_vector_for_positions(self, *pos):
        if not self.b_tensor_takes_coords:
            b_vect_vals = [sc.b_vector_for_positions(*pos) for sc in self.subcoordinates]
            reindexed = [Vector(shape=(len(self.atoms),)) for _ in self.subcoordinates]
            for isc, sc in enumerate(b_vect_vals):
                for idx in xrange(b_vect_vals[isc].shape[0]):
                    new_idx = self._idxmap[isc][idx//3] + idx % 3
                    reindexed[isc][new_idx] = b_vect_vals[isc][idx]
            if callable(self.b_vector_function):
                return self.b_vector_function(*reindexed)
            else:
                return self.function(*reindexed)
        elif callable(self.b_vector_function):
            if self.b_vector_takes_coords:
                return self.b_vector_function(*self.subcoordinates)
        elif self.b_vector_takes_coords:
            raise ValueError("b_tensor_function given to FunctionCoordinate is not callable")
        else:
            return NotImplemented

    def analytic_b_tensor_for_order(self, order):
        if not self.b_tensor_takes_coords:
            b_tens_vals = [sc.analytic_b_tensor_for_order(order) for sc in self.subcoordinates]
            if any(b is NotImplemented for b in b_tens_vals):
                return NotImplemented
            reindexed = [Tensor(shape=(len(self.atoms)*3,)*order) for _ in self.subcoordinates]
            for isc, sc in enumerate(b_tens_vals):
                for idxs in product(*map(xrange, b_tens_vals[isc].shape)):
                    new_idxs = tuple(self._idxmap[isc][idx//3]*3 + idx % 3 for idx in idxs)
                    reindexed[isc][new_idxs] = b_tens_vals[isc][idxs]
            if callable(self.b_tensor_function):
                return self.b_tensor_function(*reindexed)
            else:
                return self.function(*reindexed)
        elif callable(self.b_tensor_function):
            return self.b_tensor_function(*self.subcoordinates)
        elif self.b_tensor_takes_coords:
            raise ValueError("b_tensor_function given to FunctionCoordinate is not callable")
        else:
            if order == 1:
                return self.b_vector_for_positions([a.pos for a in self.atoms])
            return NotImplemented

    ###################
    # Private Methods #
    ###################

    def _get_state(self):
        # Effectively disable caching
        return float('nan')
        # We could use state caching like this:
        #     states = []
        #     for sc in self.subcoordinates:
        #         states.append(sc._get_state())
        #     return tuple(states)
        # but it's probably more efficient just
        # to use the existing caching in the subcoordinates,
        # unless the computation of the b vector for a given
        # function is expensive relative to the cost of
        # computing the subcoordinate's b vectors

    def _recompute_b_vector(self):
        if self.b_vector_same_function:
            self._bvector = self.function(*[sc.b_vector for sc in self.subcoordinates])
        elif callable(self.b_vector_function):
            if self.b_vector_takes_coords:
                self._bvector = self.b_vector_function(*self.subcoordinates)
            else:
                self._bvector = self.b_vector_function(*[sc.b_vector for sc in self.subcoordinates])
        else:
            raise ValueError("don't known how to get b vector for FunctionCoordinate '{}'")

class SumOfInternalCoordinates(InternalFunctionCoordinate):
    """
    """

    ##############
    # Attributes #
    ##############

    coefficients = []

    ##################
    # Initialization #
    ##################

    def __init__(self, subcoordinates, coefficients=None, **kwargs):
        if sanity_checking_enabled:
            if not all(c.units == subcoordinates[0].units for c in subcoordinates):
                raise ValueError("all subcoordinates of a sum coordinate must have the same units")
        if coefficients is None:
            self.coefficients = [1.0] * len(subcoordinates)
        else:
            self.coefficients = coefficients[:]
        kwargs.update(
            function=self.compute_sum,
            b_vector_same_function=True,
            b_tensor_same_function=True,
            subcoordinates=subcoordinates,
        )
        super(SumOfInternalCoordinates, self).__init__(**kwargs)

    def _validate_initialization(self):
        if len(self.subcoordinates) != len(self.coefficients):
            raise ValueError('dimension mismatch: {} != {}'.format(
                len(self.subcoordinates), len(self.coefficients)))
        super(SumOfInternalCoordinates, self)._validate_initialization()

    ###################
    # Special Methods #
    ###################

    def __copy_kwargs__(self):
        kw = super(SumOfInternalCoordinates, self).__copy_kwargs__()
        kw.update(
            coefficients=self.coefficients
        )
        return kw

    ###########
    # Methods #
    ###########

    def compute_sum(self, *items):
        if sanity_checking_enabled:
            if len(items) != len(self.coefficients):
                raise ValueError('dimension mismatch: {} != {}'.format(
                    len(items), len(self.coefficients)))
        return sum(c * v for c, v in zip(self.coefficients, items))

    def generate_name(self, one_based=True):
        off = 1 if one_based else 0
        rv = '('
        for isc, sc in enumerate(self.subcoordinates):
            if isc == 0 and self.coefficients[isc] < 0:
                rv += '-'
            if abs(self.coefficients[isc]) != 1:
                rv += '{:.1f}'.format(abs(self.coefficients[isc]))
            rv += sc.coordinate_symbol
            rv += ''.join(str(a.index + off) for a in sc.atoms)
            if isc != len(self.subcoordinates) - 1:
                if self.coefficients[isc+1] > 0:
                    rv += '+'
                elif self.coefficients[isc+1] < 0:
                    rv += '-'
        rv += ')'
        return rv

# TODO: SumOfPeriodicCoordinates which inherits from PeriodicCoordinate

######################
## Dependent Imports #
#####################

from grendel.coordinates.simple_internal_coordinate import SimpleInternalCoordinate

# needed for dynamic typechecking, do not delete
if type_checking_enabled:
    from grendel.coordinates.simple_internal_coordinate import SimpleInternalCoordinate
    from grendel.representations.internal_representation import InternalRepresentation
