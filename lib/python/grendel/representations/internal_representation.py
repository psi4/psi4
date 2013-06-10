from __future__ import print_function
from collections import defaultdict
from copy import deepcopy, copy
from itertools import chain, tee, product

import math
import sys

import numpy as np

from grendel import type_checking_enabled, sanity_checking_enabled
from grendel.differentiation.derivative_tensor import RepresentationDependentTensor
from grendel.external.combinatorics import prod

from grendel.gmath.matrix import Matrix
from grendel.gmath.vector import Vector, magnitude
from grendel.gmath.tensor import Vector, Tensor
from grendel.representations.cartesian_representation import CartesianRepresentation
from grendel.representations.representation import Representation, RepresentationError
from grendel.util.containers import CustomHashDict, LRUDict

from grendel.util.decorators import see_abstract_doc, IterableOf, typecheck
from grendel.util.iteration import stutter, grouper, I_lubar
from grendel.util.overloading import overloaded, OverloadedFunctionCallError
from grendel.util.strings import indented, shortstr
from grendel.util.units import AngularUnit, Radian, DistanceUnit
from grendel.util.units.errors import IncompatibleUnitsError
from grendel.util.units.value_with_units import hasunits

__all__ = [
    "InternalRepresentation"
]

class MaximumGeometryIterationsExceededError(RuntimeError): pass

class _BTensorCache(object):

    representation = None

    def __init__(self, representation):
        self.representation = representation
        self.value_dict = CustomHashDict(_BTensorCache.tensor_element_hash)

    @staticmethod
    def tensor_element_hash(tup):
        return sum(map(hash, tup))

#TODO make coordinates truely immutable by making copies
class InternalRepresentation(Representation):
    """
    An internal coordinate representations of a molecule.
    TODO: More thoroughly document what this means.

    """

    #############
    # Constants #
    #############

    STRE_NAMES = ["str", "stre", "bond"]
    BEND_NAMES = ["bend", "angle"]
    TORS_NAMES = ["tors"]

    ####################
    # Class Attributes #
    ####################

    zero_eigenvalue = 1e-8

    default_units = {
        AngularUnit : AngularUnit.default,
        DistanceUnit : DistanceUnit.default
    }

    ######################
    # Private Attributes #
    ######################

    _b_tensor_for_cart_rep = None
    _b_finite_difference_managers = None
    _icoord_map = None

    ##################
    # Initialization #
    ##################

    #TODO Named internal coordinates
    @overloaded
    def __init__(self, *args, **kwargs):
        """
        InternalRepresentation(*args, **kwargs)

        **Signatures**
            * ``InternalRepresentation(molecule, coords_string, one_based=True)``
            * ``InternalRepresentation(molecule, coords=None)``

        Parameters
        ----------
        molecule: `Molecule`
            The molecule to be represented
        coords_string : `str`
            A string of coordinates (newline separated) defined in the freeform notation allowed by grendel:
            TODO Describe notation here
        one_based : `bool`
            Whether the string given uses atom indices starting at 1 or 0 (defaults to True, meaning starts at 1)
        coords : `list` of `InternalCoordinate`
            The internal coordinates that together compose a representation.  If None, the coords attribute is initialized
            to an empty list.


        Examples
        --------

        >>> from grendel import SampleMolecules
        >>> mol = SampleMolecules["H2"]
        >>> i = InternalRepresentation(mol, '''
        ...     stre 1 2
        ... ''')
        >>> i
        InternalRepresentation(Molecule([
            Atom('H', [  0.00000000,  0.00000000,  0.00000000 ] ),
            Atom('H', [  0.00000000,  0.00000000,  0.74300000 ] )
        ]), [BondLength(0, 1)])
        >>> mol = grendel.chemistry.SampleMolecules["water"]
        >>> i = InternalRepresentation(mol,
        ...     '''
        ...     stre   1   2
        ...     stre   2   3
        ...     bend   1   2   3
        ...     '''
        ... )

        """
        raise OverloadedFunctionCallError

    @__init__.overload_with(
        molecule='Molecule',
        coords_string=basestring,
        one_based=bool,
        units=(dict, None))
    def __init__(self,
            molecule,
            coords_string,
            one_based=True,
            units=None):
        off = 1 if one_based else 0
        self.molecule = molecule
        atoms = self.molecule.atoms
        coords = []
        for line in coords_string.splitlines():
            # ignore blank lines
            if len(line.strip()) == 0: continue
            vals = line.split()
            name = vals[0]
            nums = vals[1:]
            #TODO More coordinate types
            if name.lower() in self.STRE_NAMES:
                if len(nums) != 2:
                    raise StandardError("Invalid coordinate specification string:  " + line)
                coord = BondLength([
                    atoms[int(nums[0])-off],
                    atoms[int(nums[1])-off]
                ], parent=self, index=len(coords))
                coords.append(coord)
            elif name.lower() in self.BEND_NAMES:
                if len(nums) != 3:
                    raise StandardError("Invalid coordinate specification string:  " + line)
                coord = BondAngle([
                    atoms[int(nums[0])-off],
                    atoms[int(nums[1])-off],
                    atoms[int(nums[2])-off]
                ], parent=self, index=len(coords))
                coords.append(coord)
            elif name.lower() in self.TORS_NAMES:
                if len(nums) != 4:
                    raise StandardError("Invalid coordinate specification string:  " + line)
                coord = Torsion([
                    atoms[int(nums[0])-off],
                    atoms[int(nums[1])-off],
                    atoms[int(nums[2])-off],
                    atoms[int(nums[3])-off]
                ], parent=self, index=len(coords))
                coords.append(coord)
            else:
                raise StandardError("Unknown or unimplemented internal coordinate string:  " + line)
        self.__init__(molecule=molecule, coords=coords, units=units)

    @__init__.overload_with(
        molecule='Molecule',
        coords=(IterableOf('InternalCoordinate'), None),
        units=(dict, None))
    def __init__(self, molecule, coords=None, units=None):
        self._b_tensor_for_cart_rep = {}
        self.molecule = molecule
        self._init_units(units)
        self.coords = []
        for i, coord in enumerate(coords or []):
            cunits=self.units[coord.default_delta.units.genre] if hasunits(coord.default_delta) else None
            new_coord = coord.copy_for_representation(self, index=i, units=cunits)
            self.coords.append(new_coord)
        self.molecule.internal_representations.append(self)

    def _init_units(self, units):
        self._units = units
        units_unset = False
        if self.units is None:
            units_unset = True
            self._units = copy(InternalRepresentation.default_units)
        else:
            if sanity_checking_enabled:
                error = ValueError("units passed into InternalRepresentation constructor must be a "
                                   " dictionary with AngularUnit and DistanceUnit as keys.")
                if not isinstance(self.units, dict):
                    raise error
                elif not AngularUnit in self.units:
                    raise error
                elif not DistanceUnit in self.units:
                    raise error
                elif len(self.units) != 2:
                    raise error
        if units_unset:
            self._units[DistanceUnit] = self.molecule.cartesian_units

    ##############
    # Properties #
    ##############

    @property
    def a_matrix(self):
        """ The El'yashevich--Wilson B matrix pseudoinverse, with the Sayvetz conditions used
        to generate the inverse mass external conditions.
        """
        u = self.molecule.inverse_mass_matrix
        B = self.b_matrix
        return u * B.T * (B * u * B.T).I

    @property
    def b_matrix(self):
        """ The El'yashevich--Wilson B matrix for the molecule in the representation.
        """
        rows = [coord.b_vector for coord in self]
        return Matrix(rows)

    @property
    def g_matrix(self):
        B = self.b_matrix
        # TODO G matrix units
        G = Matrix(shape=(len(self),)*2)
        for i, j in product(xrange(len(self)), repeat=2):
            G[i,j] = sum(
                B[i,k]*B[j,k] / self.molecule.atoms[k//3].mass
                    for k in xrange(3*self.molecule.natoms)
            )
        return G

    ###################
    # Special Methods #
    ###################

    def __repr__(self):
        # TODO find out why this doesn't always work
        try:
            return "InternalRepresentation(" + repr(self.molecule) + ", [" + ", ".join([c.short_repr() for c in self.coords]) + "])"
        except:
            return super(InternalRepresentation, self).__repr__()

    def __str__(self):
        return "InternalRepresentation of {} with coordinates\n{}".format(
            shortstr(self.molecule),
            indented('\n'.join('[{}] {}: {}'.format(
                i, shortstr(coord), str(coord.value_with_units)
            ) for i, coord in enumerate(self.coords)))
        )

    ###########
    # Methods #
    ###########

    def b_tensor(self, cartesian_representation=None):
        cartesian_representation = cartesian_representation or self.molecule.cartesian_representation
        if sanity_checking_enabled:
            if not cartesian_representation.molecule is self.molecule:
                raise ValueError("the molecule for the CartesianRepresentation that a B tensor is computed with"
                                 " repsect to must be the same exact instance as the molecule corresponding to the"
                                 " InternalRepresentation of which the B tensor is being computed.")
        if cartesian_representation in self._b_tensor_for_cart_rep:
            return self._b_tensor_for_cart_rep[cartesian_representation]
        #--------------------------------------------------------------------------------#
        # Looks like we don't already have it; so lets create it
        def b_tensor_calculator(tensor, indices, cartesian_representation=cartesian_representation, internal_representation=self):
            # we expect indices to always be integers
            if tensor.compute_function is None:
                # See if we already have it (we could have made the assignment on a previous call...
                return tensor[indices]
            else:
                #========================================#
                # build up some arrays of indices and coordinates
                internal_coord = internal_representation[indices[0]]
                internal_coord_idx = indices[0]
                cart_coords = [cartesian_representation[i] for i in indices[1:]]
                cart_coord_idxs = indices[1:]
                internal_carts = []
                for atom_idx in internal_coord.atom_indices:
                    internal_carts.extend([3*atom_idx, 3*atom_idx+1, 3*atom_idx+2])
                all_cart_coords = []
                for cart in cart_coord_idxs:
                    all_cart_coords.extend(CartesianRepresentation.same_atom_indices(cart))
                #========================================#
                # check the cache
                vects = [cartesian_representation[i].atom.position for i in cart_coord_idxs]
                if isinstance(internal_coord, SimpleInternalCoordinate):
                    cache_resp = SimpleInternalCoordinate._check_b_tensor_cache(type(internal_coord), len(vects), *vects)
                    if cache_resp is not None:
                        # TODO store all of the values for cache_resp in tensor
                        return cache_resp[tuple(idx % 3 for idx in cart_coord_idxs)]
                #========================================#
                # see if we know how to compute it analytically
                analytic_response = internal_coord.fill_b_tensor_analytically(
                    tensor.collection,
                    len(indices) - 1)
                if analytic_response is NotImplemented:
                    # compute by finite difference
                    return internal_coord.finite_difference_b_tensor_element(*cart_coords)
                else:
                    # success!
                    return tensor[indices]
        #--------------------------------------------------------------------------------#
        # Create a new b tensor
        b_tens = DerivativeCollection(
            representation=self,
            first_dimension_different=True,
            secondary_representation=cartesian_representation,
            uncomputed=True,
            compute_function=b_tensor_calculator,
            #permutational_symmetry=True
            permutational_symmetry=False  # For now
        )
        self._b_tensor_for_cart_rep[cartesian_representation] = b_tens
        return b_tens

    @see_abstract_doc
    def copy_with_molecule(self, molecule):
        ret_val = self.__class__(molecule, units=self.units)
        for coord in self:
            ret_val.add_coordinate_copy(coord)
        molecule.internal_representations.insert(0, ret_val)
        return ret_val

    @typecheck(coordinate='InternalCoordinate')
    @see_abstract_doc
    def add_coordinate_copy(self, coordinate):
        self.fail_if_frozen()
        self.coords.append(coordinate.copy_for_representation(self))

    @see_abstract_doc
    def displaced_by(self, disp, tol=None, maxiter=None):
        # TODO Units for tolerance
        if type_checking_enabled:
            if not isinstance(disp, Displacement):
                raise TypeError
        if sanity_checking_enabled:
            self.validate()
        if tol is None:
            tol = disp.tolerance
        if maxiter is None:
            maxiter = disp.max_iterations
        #--------------------------------------------------------------------------------#
        disp_mol = deepcopy(self.molecule)
        newrep = self.copy_with_molecule(disp_mol)
        #--------------------------------------------------------------------------------#
        ## Jeremy's method/optking method
        #if False:
        #  TODO Make this a fallback?
        #    if sanity_checking_enabled: progressive_vals = []
        #    total_disp_cart = Vector(shape=(3*self.molecule.natoms,), default_val=0.0)
        #    error = tol + 1
        #    niter = 0
        #    def_to_rad = AngularUnit.default.to(Radian)
        #    #--------------------------------------------------------------------------------#
        #    while error > tol and niter < maxiter:
        #        error = 0.0
        #        niter += 1
        #        total_disp_cart.fill(0.0)
        #        if sanity_checking_enabled: progressive_vals.append([])
        #        for i, coord in enumerate(self.coords):
        #            try:
        #                # Find the difference in the desired value and the current value
        #                current_value = newrep[i].value
        #                if sanity_checking_enabled: progressive_vals[-1].append(current_value)
        #                desired_value = disp.desired_values[i]
        #                diff = desired_value - current_value
        #                if coord.default_delta.units.genre is AngularUnit:
        #                    diff *= def_to_rad
        #                # accumulate the error
        #                error += abs(diff)
        #                # And find the normalization factor for the b-vector shift
        #                b_vect = coord.b_vector
        #                shift = b_vect.dot(b_vect)
        #                if shift < Vector.zero_cutoff:
        #                    raise ZeroDivisionError, "B vector magnitude of {0} is less than {1}, so displacement cannot be normalized.".format(shift, coord.b_vector)
        #                # ...then accumulate the normalized b vector into the total displacement vector
        #                total_disp_cart += b_vect * (diff / shift)
        #            except Exception as e:
        #                if sanity_checking_enabled:
        #                    print("Exception raised in InternalRepresentation.displaced_by().\n"
        #                        "The values of the coordinates by iteration were:\n{}".format(
        #                            '\n'.join("[{}]".format(', '.join(map(str, pvals))) for pvals in progressive_vals)
        #                        ), file=sys.stderr
        #                    )
        #                    sys.stderr.flush()
        #                    raise e
        #                else:
        #                    raise e
        #        if sanity_checking_enabled:
        #            if any(math.isnan(n) for n in progressive_vals[-1]):
        #                raise RuntimeError('Something went horribly wrong in internal representation displacement:\n'
        #                                   'One coordinate became not-a-number (nan).\n'
        #                                   'The values for the coordinates by iteration were:\n{}'.format(
        #                    '\n'.join("[{}]".format(', '.join(map(str, pvals))) for pvals in progressive_vals)
        #                ))
        #        disp_mol.displace(total_disp_cart)
        #        # Coordinate B vectors will automatically be recomputed.
        #        if niter >= maxiter:
        #            raise RuntimeError("Maximum geometry displacement iterations exceeded.  Most likely, this is "
        #                                "because your internal coordinates are nonsensical.  If you're sure that's"
        #                                " not the case, try increasing Displacement.max_iterations.")
        #    #--------------------------------------------------------------------------------#
        #    # Prepend the new representation
        #    disp_mol.internal_representations.insert(0, newrep)
        #    return disp_mol, newrep
        ##--------------------------------------------------------------------------------#
        # Dr. Allen's method (could be improved by including higher order A tensor contributions)
        # Debugging code...
        #def sanity_iters(t, p, s, c, n):
        #    print("Exception raised in InternalRepresentation.displaced_by().\n"
        #          "Target values were:\n[{}]\nThe values of the coordinates by iteration were:\n{}".format(
        #        ", ".join(str(des) for des in t),
        #        '\n'.join("[{}]\n   shifts: [{}]\n   cartesian shifts: [{}]\n   norm of shift: {}".format(
        #            ', '.join(map(str, pvals)),
        #            ', '.join(map(str, shiftvals)),
        #            ', '.join(map(str, cshiftvals)),
        #            str(norm)
        #        ) for pvals, shiftvals, cshiftvals, norm in zip(p, s, c, n))
        #    ), file=sys.stderr)
        #    sys.stderr.flush()
        #----------------------------------------#
        A = self.a_matrix
        desired = disp.desired_values
        shift = disp.disp_vect
        for i, c in enumerate(self.coords):
            if c.units.genre is AngularUnit:
                desired[i] *= c.units.to(Radian)
                shift[i] *= c.units.to(Radian)
        geom_vect = self.molecule.cartesian_representation.value
        niter = 0
        norm = math.sqrt(shift.dot(shift))
        # Debugging arrays...
        #progressive_vals = []; shifts = []; cart_shifts = []; norms = []
        while norm > tol and niter < maxiter:
            niter += 1
            try:
                shift_amt = A * shift
                geom_vect += shift_amt
                current = self.value_for_matrix(geom_vect.reshape((self.molecule.natoms, 3)))
                shift = desired - current
                norm = math.sqrt(shift.dot(shift))
                #if sanity_checking_enabled:
                #    progressive_vals.append(current)
                #    shifts.append(shift)
                #    cart_shifts.append(shift_amt)
                #    norms.append(norm)
            except Exception as e:
                #if sanity_checking_enabled:
                #    sanity_iters(desired, progressive_vals, shifts, cart_shifts, norms)
                raise e
        if niter >= maxiter:
            #if sanity_checking_enabled:
            #    sanity_iters(desired, progressive_vals, shifts, cart_shifts, norms)
            raise MaximumGeometryIterationsExceededError(
                "Maximum geometry displacement iterations exceeded.  Most likely, this is"
                " because your internal coordinates are nonsensical or you are trying to do too"
                " large of a displacement.  If you're sure that's not the case, try increasing"
                " Displacement.max_iterations.")
        disp_mol.xyz = geom_vect.reshape((self.molecule.natoms, 3))
        disp_mol.internal_representations.insert(0, newrep)
        #if sanity_checking_enabled:
        #    if any(math.isnan(n) for n in newrep.values):
        #        sanity_iters(desired, progressive_vals, shifts, cart_shifts, norms)
        #        raise RuntimeError('Something went horribly wrong in internal representation displacement: '
        #                           'One coordinate became not-a-number (nan).')
        return disp_mol, newrep

    def validate(self):
        bmat = self.b_matrix
        evals = np.linalg.eigvalsh(bmat * bmat.T)
        nvalid = len([e for e in evals if abs(e) > InternalRepresentation.zero_eigenvalue])
        if nvalid != self.molecule.ninternals:
            raise RepresentationError("Have {} nonredundant internal coordinates, but need {}.\nCoordnates are:\n{}".format(
                nvalid,
                self.molecule.ninternals,
                indented('\n'.join('[' + str(i) + '] ' + shortstr(c) for i, c in enumerate(self.coords)))
            ))

    def transform_tensor(self, tensor, to_representation):
        # TODO sanity checks similar to CartesianRepresentation.transform_tensor
        if to_representation is self:
            return copy(tensor)
        elif isinstance(to_representation, InternalRepresentation):
            # use an intermediary
            intermed = self.molecule.cartesian_representation
            return tensor.in_representation(intermed).in_representation(to_representation)
        elif isinstance(to_representation, CartesianRepresentation):
            self.freeze()
            to_representation.freeze()
            shape=(len(to_representation),)*len(tensor.shape)
            ret_val = RepresentationDependentTensor(
                shape=shape,
                representation=to_representation)
            if len(tensor.shape) == 1:
                # TODO Check and make sure the cartesian_representation is the one with which the b tensor was computed
                # TODO units, keeping in mind that B tensor always has angular units of radians.
                B = self.b_matrix
                ret_val[...] = B.T * tensor.view(Vector)
            else:
                raise NotImplementedError("use transform_forcefield instead")
            return ret_val
        elif isinstance(to_representation, NormalRepresentation):
            # use an intermediary
            intermed = self.molecule.cartesian_representation
            return tensor.in_representation(intermed).in_representation(to_representation)
        else:
            raise NotImplementedError

    def transform_forcefield(self, ff, to_representation):
        self.freeze()
        to_representation.freeze()
        if sanity_checking_enabled and ff.representation is not self:
            raise ValueError("int_rep.transform_forcefield must be given a force field whose"
                             " representation is int_rep")
        rv = ForceField(to_representation, ff.max_order, ff.molecular_property, ff.property_units)
        if isinstance(to_representation, CartesianRepresentation):
            B = self.b_tensor(cartesian_representation=to_representation)
            rv.for_order(1)[...] = self.transform_tensor(ff.for_order(1), to_representation)
            for o in xrange(1, ff.max_order+1):
                B.for_order(o).fill()
            if ff.max_order >= 2:
                # Do this by hand to speed things up...
                rv['i1,i2'] = ff['p1']*B['p1,i1,i2'] + ff['p1,p2']*B['p1,i1']*B['p2,i2']
                # TODO Permutational symmetry?
                #if ff.max_order >= 3:
                #    # Do this by hand also to speed things up...
                #    # The three middle terms could be done more efficiently by first summing
                #    #   over them and then doing the contraction, but this is not yet implemented
                #    #   for einstein summation
                #    rv['i1,i2,i3'] = ff['p1']*B['p1,i1,i2,i3'] + ff['p1,p2,p3']*B['p1,i1']*B['p2,i2']*B['p3,i3']
                #    for s in I_lubar(3, 2, ('i1','i2','i3')):
                #        rv['i1,i2,i3'] += ff['p1,p2']*B[('p1',)+s[0]]*B[('p2',) + s[1]]
                if ff.max_order >= 3:
                    def idxstrs(letter, num): return tuple(map(lambda x: letter + str(x), xrange(num)))
                    for order in xrange(3, ff.max_order + 1):
                        iidxs = idxstrs('i', order)
                        for k in xrange(1, order+1):
                            for cartsets in I_lubar(order, k, iidxs):
                                rv[iidxs] += ff[idxstrs('p', k)] * prod(
                                    B[('p'+str(m),) + cartsets[m]] for m in xrange(k)
                                )
            return rv
        elif isinstance(to_representation, InternalRepresentation):
            # use an intermediary
            intermed = self.molecule.cartesian_representation
            return ff.in_representation(intermed).in_representation(to_representation)
        elif isinstance(to_representation, NormalRepresentation):
            # use an intermediary
            intermed = self.molecule.cartesian_representation
            return ff.in_representation(intermed).in_representation(to_representation)
        else:
            raise NotImplementedError


    #-------------------------------------#
    # Inquiry methods (which return bool) #
    #-------------------------------------#

    def is_valid(self):
        try:
            self.validate()
            return True
        except RepresentationError:
            return False

    ###################
    # Private Methods #
    ###################

    def _gen_icoord_map(self):
        if isinstance(self._icoord_map, dict):
            return
        else:
            self._icoord_map = dict()
            for icoord in self:
                key = (icoord.__class__, tuple(icoord.atoms))
                self._icoord_map[key] = icoord





#####################
# Dependent Imports #
#####################

from grendel.chemistry.molecule import Molecule

from grendel.coordinates.bond_length import BondLength
from grendel.coordinates.bond_angle import BondAngle
from grendel.coordinates.torsion import Torsion
from grendel.coordinates.internal_coordinate import InternalCoordinate
from grendel.coordinates.simple_internal_coordinate import SimpleInternalCoordinate

from grendel.differentiation.displacement import Displacement
from grendel.differentiation.forcefield import ForceField
from grendel.differentiation.derivative_collection import DerivativeCollection
from grendel.representations.normal_representation import NormalRepresentation


