"""
"""
from __future__ import print_function, division
from itertools import product
import math
import numpy as np
from copy import copy, deepcopy
from grendel.differentiation.derivative_tensor import DerivativeTensor
from grendel.differentiation.forcefield import ForceField
from grendel.external.combinatorics import prod
from grendel.gmath.matrix import Matrix
from grendel.gmath.vector import Vector
from grendel.representations.representation import Representation
from grendel.util.decorators import typechecked
#noinspection PyUnresolvedReferences
from grendel.util.units import DistanceUnit, Bohr, Angstroms, Kilogram, Meter, Joules, Hartrees, \
                               AMU, Hertz, Wavenumbers, Attojoules, ReducedPlanckConstant, \
                               PlanckConstant, SpeedOfLight, AvogadrosNumber, Centimeters, Second


class NormalRepresentation(Representation):

    ##################
    # Initialization #
    ##################

    @typechecked(
        forcefield=ForceField
    )
    def __init__(self, forcefield):
        self.molecule = forcefield.base_molecule
        self.cartesian_representation = self.molecule.cartesian_representation
        if isinstance(forcefield.representation, InternalRepresentation):
            self.internal_representation = forcefield.representation
        _cart_ff = forcefield.in_representation(self.cartesian_representation)
        Fx = _cart_ff.for_order(2).view(Matrix)
        Fx_start_units = _cart_ff.property_units / self.cartesian_representation.units**2
        #----------------------------------------#
        #TODO Think this through...
        self._units = self.cartesian_representation.units / AMU**(1/2)
        #----------------------------------------#
        Fxm = Fx.transformed(self.molecule.inverse_sqrt_mass_matrix)
        # Symmetrize for numerical stability:
        Fxm = Fxm.symmetrized()
        evals, evecs = Fxm.eigensystem()
        Bmat = self.molecule.inverse_sqrt_mass_matrix * evecs
        to_si = (Fx_start_units / AMU).to(Joules/(Meter**2 * Kilogram))
        #TODO order by irrep in cannonical ordering
        self.frequencies = [(math.sqrt(e * to_si) if e > 0 else complex(0,math.sqrt(-e))) \
                            / (2.0 * math.pi) * Hertz.to(Wavenumbers) * Wavenumbers for e in evals]
        self.frequencies = self.frequencies[6:]
        #----------------------------------------#
        # Now make the Coordinate instances
        self.coords = []
        # TODO linear molecules
        # TODO Decide what to do with this stuff
        external = zip(self.frequencies[:6], list(Bmat.col_iter)[:6])
        for i, (freq, vec) in enumerate(zip(self.frequencies, list(Bmat.col_iter)[6:])):
            new_coord = NormalCoordinate(
                vec,
                freq,
                parent=self,
                index=i,
                units=self.units
            )
            self.coords.append(new_coord)
        #----------------------------------------#
        # Use the GF-matrix method to get the L matrix elements
        # (Disabled for now)
        # TODO re-enable this as a means of checking ourselves
        #F = hessian.view(Matrix)
        #G = rep.g_matrix
        #G_12 = G.sqrt_matrix()
        #GFG = F.transformed(G_12)
        #print(GFG.formatted_string(
        #    name="GFG",
        #    float_digits=3,
        #    float_type_string='f')
        #)
        ##TODO choose a phase unambiguously
        #evals, evecs = GFG.eigensystem()
        ##TODO order by irrep in cannonical ordering
        #old_evals = copy(evals)
        #evals, evecs = zip(*sorted((val, vec) for val, vec in zip(evals, evecs.col_iter)))
        #evecs = Matrix(evecs).T
        ##evecs = zip(*sorted((val, vec) for val, vec in zip(old_evals, evecs)))[1]
        #evecs = G_12 * Matrix(evecs)
        ###TODO figure out why this is the correct unit conversion
        #evals = Vector(evals) * (Hartrees**(1.0/2.0)/(AMU*Bohr)).to(Joules**(1.0/2.0)/(Kilograms*Meters))
        #evals = evals / (rep.units[DistanceUnit]**2).to(Bohr**2)
        #self._frequencies = Vector([math.sqrt(val) if val > 0 else -math.sqrt(-val) for val in evals])
        #self._frequencies *= Hertz.to(Wavenumbers) * Wavenumbers
        #----------------------------------------#
        # Now make the Coordinate instances

    ##############
    # Properties #
    ##############

    @property
    def b_matrix(self):
        return Matrix([c.b_vector for c in self.coords])


    ##################
    # Static Methods #
    ##################

    @staticmethod
    def convert_to_wavenumbers(tensor, energy_units):
        if not isinstance(tensor, DerivativeTensor):
            raise TypeError
        if not isinstance(tensor.representation, NormalRepresentation):
            raise ValueError
        order = len(tensor.shape)
        if order >= 2:
            cart_units = tensor.representation.cartesian_representation.units
            start_units = energy_units / (cart_units**order * AMU**(order/2))
            c = SpeedOfLight.in_units(Centimeters/Second).value
            h = PlanckConstant .value
            to_si = start_units.to(Joules/(Meter**order * Kilogram**(order/2)))
            ret_val = tensor * (to_si / ((4.0 * math.pi**2 * c / h)**(order/2) * h * c))
            for idxs in product(*map(xrange, tensor.shape)):
                ws = [tensor.representation.frequencies[i] for i in idxs]
                # TODO handle imaginary freqs
                ret_val[idxs] = ret_val[idxs] / math.sqrt(prod(ws))
            return ret_val
        else:
            raise NotImplementedError



    ###########
    # Methods #
    ###########


    #--------------------------------#
    # Methods abstract in Coordinate #
    #--------------------------------#

    #TODO decide what to do about this.  It doesn't make much sense to have it here, but the representation protocol requires it (for now, at least)
    def add_coordinate_copy(self, coordinate):
        raise NotImplementedError

    #TODO decide what to do about this.  It doesn't make much sense to have it here, but the representation protocol requires it (for now, at least)
    def copy_with_molecule(self, molecule):
        raise NotImplementedError

    def displaced_by(self, disp, tol=None, maxiter=None):
        raise NotImplementedError


#####################
# Dependent Imports #
#####################

from grendel.coordinates.normal_coordinate import NormalCoordinate

# Needed for dynamic typechecking, do not delete
from grendel.representations.internal_representation import InternalRepresentation

