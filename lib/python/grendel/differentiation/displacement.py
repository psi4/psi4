from __future__ import print_function, division
from collections import Iterable
from itertools import product
from numbers import Number
import numpy as np
from grendel import type_checking_enabled, sanity_checking_enabled
from grendel.coordinates.coordinate import Coordinate
from grendel.chemistry.molecule import Molecule
from grendel.differentiation.finite_difference import  FiniteDifferenceVariable
from grendel.gmath.vector import Vector
from grendel.representations.internal_representation import MaximumGeometryIterationsExceededError
from grendel.util.decorators import CachedProperty, typechecked, IterableOf, with_flexible_arguments
from grendel.util.units import AngularUnit, Radians
from grendel.util.metaprogramming import ReadOnlyAttribute
from grendel.util.units import IncompatibleUnitsError, compatible_units, hasunits
from grendel.util.units.value_with_units import strip_units


class DisplacementGenerator(object):
    """
    """

    #############
    # Constants #
    #############

    Forward = object()
    Backward = object()
    Both = object()

    ####################
    # Class Attributes #
    ####################

    verbose = False

    ##############
    # Attributes #
    ##############

    representation = None
    coordinates = None
    ranges = None
    npoints_for_coord = None
    displacements = None
    displaced_molecules = None
    direction = None

    ######################
    # Private Attributes #
    ######################

    _coord_idxs = None
    _reps = None


    ##################
    # Initialization #
    ##################

    @with_flexible_arguments(
        required=[
            ('representation', 'rep'),
            ('coordinates', 'coordinate', 'coord', 'coords'),
            ('ranges', 'range'),
            ('npoints', 'npts', 'npoints_for_coord', 'npts_for_coord'),
        ],
        optional=[
            ('direction',)
        ]

    )
    @typechecked(
        representation='Representation',
        coordinates=(Coordinate, IterableOf(Coordinate)),
        ranges=(Number, IterableOf((Number, IterableOf(Number)))),
        npoints=(Number, IterableOf(Number))
    )
    def __init__(self,
            representation,
            coordinates,
            ranges,
            npoints,
            direction=None):
        """
        TODO document this!
        """
        self.displacements = {}
        self.displaced_molecules = {}
        self._reps = {}
        self.representation = representation
        #----------------------------------------#
        self.coordinates = (coordinates,) if not isinstance(coordinates, Iterable) else tuple(coordinates)
        if sanity_checking_enabled:
            if any(c not in self.representation.coords for c in self.coordinates):
                raise ValueError("all coords must be in given representation")
        self._coord_idxs = [c.index for c in self.coordinates]
        #----------------------------------------#
        self.direction = direction
        if self.direction is None:
            self.direction = DisplacementGenerator.Both
        if sanity_checking_enabled and self.direction not in (
                DisplacementGenerator.Forward,
                DisplacementGenerator.Backward,
                DisplacementGenerator.Both):
            raise ValueError("Invalid value for DisplacementGenerator argument 'direction': {}".format(direction))
        #----------------------------------------#
        if isinstance(ranges, Number):
            self.ranges = (ranges,) * len(self.coordinates)
        elif all(isinstance(r, Number) for r in ranges):
            self.ranges = tuple(ranges)
            if sanity_checking_enabled and len(self.ranges) != len(self.coordinates):
                raise ValueError("dimension mismatch: {} != {}".format(len(self.ranges), len(self.coordinates)))
        else:
            raise NotImplementedError
        #----------------------------------------#
        if isinstance(npoints, Number):
            self.npoints_for_coord = (npoints,) * len(self.coordinates)
        elif all(isinstance(r, Number) for r in ranges):
            self.npoints_for_coord = tuple(npoints)
            if sanity_checking_enabled and len(self.npoints_for_coord) != len(self.coordinates):
                raise ValueError("dimension mismatch: {} != {}".format(
                    len(self.npoints_for_coord),
                    len(self.coordinates)
                ))
        else:
            raise NotImplementedError
        self.npoints_for_coord = map(int, self.npoints_for_coord)
        if sanity_checking_enabled:
            if not all(n > 1 for n in self.npoints_for_coord):
                raise ValueError("number of point for a coordinate must be greater than one (omit"
                                 " coordinate to only include equilibrium value).")
        #----------------------------------------#
        self.deltas = tuple(r / (n-1) for r, n in zip(self.ranges, self.npoints_for_coord))
        full_deltas = [0] * len(self.representation.coords)
        for delta, idx in zip(self.deltas, self._coord_idxs):
            full_deltas[idx] = delta
        #----------------------------------------#
        self._reps[(0,)*len(self.coordinates)] = self.representation
        for ndisps in product(*map(xrange, self.npoints_for_coord)):
            # Adjust the increments for the method we're using
            if self.direction is DisplacementGenerator.Forward:
                incrs = ndisps
            elif self.direction is DisplacementGenerator.Backward:
                incrs = tuple(-n for n in ndisps)
            else: # self.direction is DisplacementGenerator.Both
                # TODO start with equilibrium and move outwards
                incrs = tuple(n - npts//2 for n, npts in zip(ndisps, self.npoints_for_coord))
            # Find the representation we've already made that is closest to the displacement set
            min_rep_key = min(
                self._reps.keys(),
                key=lambda x: sum(abs(x[i]-ndisps[i]) for i in xrange(len(x)))
            )
            if self.verbose:
                print("Generating displacement <{}>".format(
                    ', '.join(c.name + "{:+d}".format(n) for c, n in zip(self.coordinates, incrs) if n != 0)
                ))
            rep = self._reps[min_rep_key]
            # See if we can generate the displaced molecule with this representation
            adj_incrs = [incrs[i] - min_rep_key[i] for i in xrange(len(ndisps))]
            adj_deltas = [inc * delta for inc, delta in zip(adj_incrs, self.deltas)]
            nominative_incrs = [0] * len(self.representation.coords)
            full_disps = [0] * len(self.representation)
            for incr, delt, idx in zip(incrs, adj_deltas, self._coord_idxs):
                nominative_incrs[idx] = incr
                full_disps[idx] = delt
            disp = Displacement(rep, full_disps, nominative_incrs)
            try:
                # Try to generate the displaced molecule...
                tmp = disp.displaced_molecule
                self.displaced_molecules[incrs] = tmp
                self.displacements[incrs] = disp
            except MaximumGeometryIterationsExceededError:
                # Create a new representation using the closest successful displacement:
                closest_disp_key = min(
                    self.displacements.keys(),
                    key=lambda x: sum(abs(x[i]-ndisps[i]) for i in xrange(len(x)))
                )
                closest_disp = self.displacements[closest_disp_key]
                rep = self.representation.copy_with_molecule(closest_disp.displaced_molecule)
                adj_incrs = [incrs[i] - closest_disp_key[i] for i in xrange(len(ndisps))]
                adj_deltas = [inc * delta for inc, delta in zip(adj_incrs, self.deltas)]
                nominative_incrs = [0] * len(self.representation.coords)
                full_disps = [0] * len(self.representation)
                for incr, delt, idx in zip(incrs, adj_deltas, self._coord_idxs):
                    nominative_incrs[idx] = incr
                    full_disps[idx] = delt
                disp = Displacement(rep, full_disps, nominative_incrs)
                self._reps[closest_disp_key] = rep
                # If this raises an error, let it be raised...we've done the best that we can
                self.displaced_molecules[incrs] = disp.displaced_molecule
                self.displacements[incrs] = disp

    ###################
    # Special Methods #
    ###################

    def __iter__(self):
        for _, disp in sorted(self.displacements.items()):
            yield disp


class Displacement(object):
    """ Encapsulates a displacement from a base molecule.

    :Attributes:

    base_molecule : `Molecule`
        The Molecule object that the displacement is relative to.
    representation : `Representation`
        The Representation object that the displacements are values in.
    disp_vect : `Vector`
        The displacement amounts in the given representation. (e.g. `disp_vect[i]` corresponds to the
        amount of displacement of the `i`th `Coordinate` in `representation`)

    """

    ####################
    # Class Attributes #
    ####################

    tolerance = 1e-12
    max_iterations = 35

    ##############
    # Attributes #
    ##############

    base_molecule = ReadOnlyAttribute('base_molecule',
        doc="""
        The Molecule object that the displacement is relative to.
        """)
    representation = None
    disp_vect = None

    ######################
    # Private Attributes #
    ######################

    _displaced_molecule = None
    _displaced_representation = None
    _increments = None
    _deltas = None

    ##################
    # Initialization #
    ##################

    @typechecked(
        representation='Representation',
        disps=IterableOf(Number),
        increments=(None, IterableOf(Number)),
        deltas=(None, IterableOf(float))
    )
    def __init__(self, representation, disps, increments=None, deltas=None):
        base = representation.molecule
        self._base_molecule = base
        self.representation = representation

        self.disp_vect = []
        for d in disps:
            self.disp_vect.append(d)
            # TODO Decide if I should handle units here or in displaced_by (probably there is better)
            #if hasunits(d):
            #    if isinstance(representation, InternalRepresentation):
            #        units = representation.units[d.units.genre]
            #    elif isinstance(representation, CartesianRepresentation):
            #        units = representation.units
            #    else:
            #        raise NotImplementedError
            #    self.disp_vect.append(strip_units(d, ) if hasunits(d) else d)
            #else:
            #    self.disp_vect.append(strip_units(d))
        self.disp_vect = Vector(self.disp_vect)
        self._increments = increments
        self._deltas = []
        # TODO Decide if I should handle units here or in displaced_by (probably there is better)
        #if deltas is not None:
        #    for d in deltas:
        #        self._deltas.append(strip_units(d, representation.units[d.units.genre]) if hasunits(d) else d)

    ##############
    # Properties #
    ##############

    @property
    def displaced_molecule(self):
        """ The displaced molecule object resulting from applying `self` to `Molecule`
        """
        if self._displaced_molecule is None:
            self._compute_displacement()
        return self._displaced_molecule

    @property
    def displaced_representation(self):
        """ The displaced molecule object resulting from applying `self` to `Molecule`
        """
        if self._displaced_representation is None:
            self._compute_displacement()
        return self._displaced_representation

    @CachedProperty
    def desired_values(self):
        return self.representation.values + self.disp_vect

    #################
    # Class Methods #
    #################

    @classmethod
    def get_default_deltas(cls, rep):
        if type_checking_enabled:
            if not isinstance(rep, Representation):
                raise TypeError
        del_list = []
        for coord in rep:
            #if coord.default_delta.units.genre is AngularUnits:
            #    del_list.append(coord.default_delta.in_units(Radians))
            #else:
            del_list.append(coord.default_delta)
        return del_list

    @classmethod
    def from_increments(cls, increments, rep, deltas=None):
        if deltas is None:
            deltas = cls.get_default_deltas(rep)
        if sanity_checking_enabled:
            if isinstance(increments, np.ndarray) and len(increments.shape) != 1:
                raise ValueError("'increments' should be 1-dimensional")
            if isinstance(deltas, np.ndarray) and len(deltas.shape) != 1:
                raise ValueError("'deltas' should be 1-dimensional")
        if type_checking_enabled:
            if not isinstance(increments, (np.ndarray, list, tuple)):
                raise TypeError
            if not isinstance(deltas, (np.ndarray, list, tuple)):
                raise TypeError
        if sanity_checking_enabled:
            if len(increments) != len(rep):
                raise ValueError("length of increments must match length of representation, but {0} != {1}".format(len(increments), len(rep)))
            if len(increments) != len(deltas):
                raise ValueError("length of increments must match deltas, but {0} != {1}".format(len(increments), len(deltas)))
        disps = []
        for inc, delta, coord in zip(increments, deltas, rep):
            if hasunits(delta):
                disps.append(inc * delta.in_units(coord.units))
            else:
                # Assume the default units for the genre
                disps.append(inc * (delta * coord.units.genre.default).in_units(coord.units))
        return Displacement(rep, disps, increments=increments, deltas=deltas)

    ###########
    # Methods #
    ###########

    def increments(self, deltas=None):
        if self._increments is not None:
            return self._increments
        if deltas is None:
            if self._deltas is not None:
                deltas = self._deltas
            else:
                deltas = Displacement.get_default_deltas(self.representation)
        increments = []
        for i in self.disp_vect/deltas:
            if abs(i - int(i)) > 1e-8:
                raise RuntimeError("the deltas given are not the original deltas used.")
            increments.append(int(i))
        self._increments = tuple(increments)
        return self._increments


    ###################
    # Private Methods #
    ###################

    def _compute_displacement(self, tol=None, maxiter=None):
        self._displaced_molecule, self._displaced_representation = self.representation.displaced_by(self, tol, maxiter)
        self._displaced_molecule.displacement = self
        return

FiniteDifferenceVariable.register(Coordinate)

#####################
# Dependent Imports #
#####################

from grendel.representations.representation import Representation
from grendel.representations.internal_representation import InternalRepresentation
from grendel.representations.cartesian_representation import CartesianRepresentation

