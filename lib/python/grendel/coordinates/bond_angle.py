from copy import copy
from functools import partial
from itertools import product, permutations, combinations_with_replacement as symmetric_product
from operator import mul
import numpy as np
from math import sin, cos, tan
from grendel import sanity_checking_enabled
from grendel.coordinates.bond_length import BondLength
from grendel.coordinates.simple_internal_coordinate import SimpleInternalCoordinate
from grendel.differentiation.derivative_collection import DerivativeCollection
from grendel.gmath import IndexRangeSet, DeclareIndexRange, IndexRange, ComputableTensor
from grendel.gmath.geometry import angle_between_vectors
from grendel.gmath.matrix import Matrix
from grendel.gmath.misc import delta
from grendel.gmath.tensor import Tensor
from grendel.gmath.vector import Vector, LightVector
from grendel.representations.cartesian_representation import X, Z, Y
from grendel.util.decorators import CachedProperty, SequenceOf, typechecked, IterableOf
from grendel.util.iteration import ordered_partitions_iter, brace_notation_iter, I_lubar
from grendel.util.sentinal_values import All
from grendel.util.units import Degrees, Radians, AngularUnit, isunit, Angstroms
from grendel.util.overloading import overloaded, OverloadedFunctionCallError

__all__ = [
    "BondAngle"
]

class BondAngle(SimpleInternalCoordinate):
    """

    """

    ####################
    # Class Attributes #
    ####################

    default_delta = 0.02*Radians
    analytic_b_orders = All
    coordinate_symbol = "phi"

    ############################
    # Private Class Attributes #
    ############################

    _btens_idx_range_set = IndexRangeSet()
    DeclareIndexRange('v', 9, index_range_set=_btens_idx_range_set).with_subranges(
        IndexRange('a', 0, 3, index_range_set=_btens_idx_range_set),
        IndexRange('b', 3, 6, index_range_set=_btens_idx_range_set),
        IndexRange('c', 6, 9, index_range_set=_btens_idx_range_set)
    )

    ##################
    # Initialization #
    ##################

    # TODO Write documentation for this analogous to BondLength
    # TODO move initialization methods up to simple_internal_coordinate, and typecheck 'self' in the overload_with to ensure correct number of atoms
    @overloaded
    def __init__(self, *args, **kwargs):
        """
        """
        raise OverloadedFunctionCallError

    @__init__.overload_with(atom1='Atom', atom2='Atom', atom3='Atom')
    def __init__(self, atom1, atom2, atom3, **kwargs):
        self.__init__([atom1, atom2, atom3], **kwargs)

    @__init__.overload_with(
        idx1=int, idx2=int, idx3=int,
        parent=('InternalRepresentation', None))
    def __init__(self, idx1, idx2, idx3, parent, one_based=True, **kwargs):
        off = 1 if one_based else 0
        self.__init__(
            (parent.molecule.atoms[i - off] for i in [idx1, idx2, idx3]),
            parent=parent,
            **kwargs
        )

    @__init__.overload_with(
        idx1=int, idx2=int, idx3=int,
        molecule='Molecule')
    def __init__(self, idx1, idx2, idx3, molecule, one_based=True, **kwargs):
        off = 1 if one_based else 0
        self.__init__(
            [molecule.atoms[i - off] for i in [idx1, idx2, idx3]],
            **kwargs
        )

    # Principal initializer (all other __init__ methods call this one)
    @__init__.overload_with(
        atoms=IterableOf('Atom'),
        parent=('InternalRepresentation', None),
        units=isunit)
    def __init__(self, atoms, parent=None, index=None, units=AngularUnit.default, **kwargs):
        self._atoms = list(atoms)
        self._index = index
        self._init(
            parent=parent,
            units=units,
            **kwargs)

    ##############
    # Properties #
    ##############

    @property
    def terminal_atoms(self):
        return self.atoms[0], self.atoms[2]

    #################
    # Class Methods #
    #################

    @classmethod
    def value_for_xyz(cls, xyz):
        v1 = LightVector.l_sub(xyz[0], xyz[1])
        v2 = LightVector.l_sub(xyz[2], xyz[1])
        return angle_between_vectors(v1, v2)

    @classmethod
    def b_vector_for_positions(cls, *args):
        if sanity_checking_enabled:
            if len(args) != 3:
                raise TypeError

        # Check the cache first
        cache_resp = SimpleInternalCoordinate._check_b_tensor_cache(cls, args)
        if cache_resp:
            return cache_resp

        # not in the cache, so compute it
        def e(j, k):
            ev = LightVector.l_sub(args[k-1], args[j-1]); ev.normalize()
            return ev
        def r(j, k):
            return LightVector.l_sub(args[k-1], args[j-1]).magnitude()
        phi = cls.value_for_xyz(list(args))
        cosph = cos(phi)
        sinph = sin(phi)
        e21, e23 = e(2,1), e(2,3)
        s1 = (e21 * cosph - e23)/(r(2,1) * sinph)
        s3 = (e23 * cosph - e21)/(r(2,3) * sinph)
        ret_val = [s1, -s1-s3, s3]
        SimpleInternalCoordinate.b_tensor_cache[(cls,) + tuple(args)] = ret_val
        return s1, -s1-s3, s3


    def analytic_b_tensor_for_order(self, order):
        #--------------------------------------------------------------------------------#
        # Now check the cache
        cache_key = (self.__class__, order) + tuple(a.pos for a in self.atoms)
        cache_resp = SimpleInternalCoordinate._check_b_tensor_cache(*cache_key)
        if cache_resp is not None:
            return cache_resp
        #--------------------------------------------------------------------------------#
        B = partial(SimpleInternalCoordinate.b_tensor_element_reindexed, self)
        #--------------------------------------------------------------------------------#
        # First order is already done elsewhere...
        if order == 1:
            # We can't use the `b_vector` attribute since that vector is indexed in the parent
            #   molecule's indexing scheme
            my_b = Vector(self.__class__.b_vector_for_positions(*[a.pos for a in self.atoms]))
            # Return now, to avoid double-caching
            return my_b
        #--------------------------------------------------------------------------------#
        # Second order terms
        if order == 2:
            # BondAngles are composed of 3 atoms (3 CartesianCoordinates each), and the
            #   order is 2, so the output Tensor will be a 9x9 Tensor
            my_b = np.ndarray(shape=(9,)*2)
            # some precomputed values
            phi = self.value * self.units.to(Radians)
            cotphi = 1.0 / tan(phi)
            cscphi = 1.0 / sin(phi)
            #========================================#
            # see comments in BondLength version
            # First, handle the terminal atom entries
            for (a_idx, a), (c_idx, c) in product(zip([0, 2], self.terminal_atoms), repeat=2):
                # Generate the temporary coordinates
                if a_idx == 0:
                    Rab = BondLength(self.atoms[0], self.atoms[1])
                    Rbc = BondLength(self.atoms[1], self.atoms[2])
                else: # a_idx == 2
                    Rab = BondLength(self.atoms[2], self.atoms[1])
                    Rbc = BondLength(self.atoms[1], self.atoms[0])
                #----------------------------------------#
                # Now iterate over the possible sets of cartesian coordinates
                for alpha, beta in product([X, Y, Z], repeat=2):
                    a_alpha, c_beta = a_idx*3 + alpha, c_idx*3 + beta
                    if a_idx == c_idx:
                        # From equation 19 in Allen, et al. Mol. Phys. 89 (1996), 1213-1221
                        a_beta = c_beta
                        other_terminal_index = 2 if a_idx == 0 else 0
                        my_b[a_alpha, a_beta] = (
                            -cotphi * B(self, a_alpha) * B(self, a_beta)
                            - cscphi * sum(
                                B(Rab, a_idx*3 + sigma, a_alpha, a_beta)
                                        * B(Rbc, other_terminal_index*3 + sigma)
                                    for sigma in [X, Y, Z]
                            )
                        )
                    else:
                        # From equation 20 in Allen, et al. Mol. Phys. 89 (1996), 1213-1221
                        my_b[a_alpha, c_beta] = (
                            -cotphi * B(self, a_alpha) * B(self, c_beta)
                            - cscphi * sum(
                                B(Rab, a_idx*3 + sigma, a_alpha)
                                        * B(Rbc, c_idx*3 + sigma, c_beta)
                                    for sigma in [X, Y, Z]
                            )
                        )
            # Now fill in the middle atom entries utilizing translational invariance
            # From equation 32 in Allen, et al. Mol. Phys. 89 (1996), 1213-1221
            for alpha, beta in product([X, Y, Z], repeat=2):
                b_beta = 3 + beta
                for a_idx, a in zip([0, 2], self.terminal_atoms):
                    a_alpha = 3*a_idx + alpha
                    my_b[b_beta, a_alpha] = my_b[a_alpha, b_beta] = \
                        -sum(my_b[a_alpha, t + beta]
                            for t in [0, 6] # 0 and 6 are the offsets for the terminal atoms
                        )
                # Now the b_alpha, b_beta entry:
                my_b[3 + alpha, 3 + beta] = sum(my_b[atom1*3 + alpha, atom2*3 + beta]
                                                for atom1, atom2 in product([0,2], repeat=2)
                                            )
        #--------------------------------------------------------------------------------#
        else:
            # behold, the general formula!
            # BondAngles are composed of 3 atoms (3 CartesianCoordinates each),
            #   so the output will be a 9x9x...x9 (`order`-dimensional) tensor
            my_b = Tensor(
                indices=','.join('v'*order),
                index_range_set=self.__class__._btens_idx_range_set
            )
            #Bphi = self._btensor
            Rab = self.get_coord(BondLength, self.atoms[0], self.atoms[1])
            Rbc = self.get_coord(BondLength, self.atoms[1], self.atoms[2])
            #Brab = Rab._btensor
            #Brbc = Rbc._btensor
            #for o in xrange(1, order):
            #    t = Bphi.for_order(o)
            #    rab = Brab.for_order(o)
            #    rbc = Rbc._btensor.for_order(o)
            #    if isinstance(t, ComputableTensor): t.fill()
            #    if isinstance(rab, ComputableTensor): rab.fill()
            #    if isinstance(rbc, ComputableTensor): rab.fill()
            #Brab.for_order(order).fill()
            #Brbc.for_order(order).fill()
            #remap_set = IndexRangeSet()
            #DeclareIndexRange('v', 6, index_range_set=remap_set).with_subranges(
            #    IndexRange('b', 3, 6, index_range_set=remap_set),
            #    IndexRange('c', 0, 3, index_range_set=remap_set)
            #)
            ##Brbc = DerivativeCollection(coordinate=Rbc, einsum_index='v', index_range_set=remap_set)
            #ba = Vector([Brab[(0,) + (Ellipsis,)*order]])
            #for o in xrange(1, order+1):
            #    Brbc.for_order(o).index_range_set = remap_set
            ## some precomputed values
            phi = self.value * self.units.to(Radians)
            cotphi = 1.0 / tan(phi)
            cscphi = 1.0 / sin(phi)
            #========================================#
            # First take care of the terminal atoms...
            def f_K(k):
                if k % 2 == 1:
                    if (k+1)/2 % 2 == 0:
                        return 1
                    else: # (k+1)/2 % 2 == 1
                        return -1
                elif k/2 % 2 == 0:
                    return cotphi
                else: # k/2 % 2 == 1
                    return -cotphi
            #----------------------------------------#
            a_idx, c_idx = 0, 2
            for num_a_coords in xrange(0, order+1):
                # outer loop over possible alphas
                for alphas in product([X, Y, Z], repeat=order):
                    a_alphas = tuple(3*a_idx + alpha for alpha in alphas[:num_a_coords])
                    c_alphas = tuple(3*c_idx + alpha for alpha in alphas[num_a_coords:])
                    # Now we have all of the specifics for the left-hand side, so compute the
                    #   right-hand side to go with it...
                    cum_sum = 0.0
                    for sigma in [X, Y, Z]:
                        cum_sum += B(Rab, 3*a_idx + sigma, *a_alphas) * B(Rbc, 3*c_idx+sigma, *c_alphas)
                    cum_sum *= -cscphi
                    for k in range(2, order+1):
                        inner_sum = 0.0
                        for s in I_lubar(order, k, a_alphas + c_alphas):
                            prod = 1.0
                            for i in range(k):
                                prod *= B(self, *s[i])
                            inner_sum += prod
                        cum_sum += f_K(k) * inner_sum
                    # Spread over permutations
                    for idxs in permutations(a_alphas + c_alphas):
                        my_b[idxs] = cum_sum
            #========================================#
            # now compute the terms involving the middle atom
            a1_idx, a2_idx, a3_idx = 0, 1, 2
            for num_a2s in xrange(1, order+1):
                # Now fill in the middle atom entries utilizing translational invariance
                # From equation 32 in Allen, et al. Mol. Phys. 89 (1996), 1213-1221
                for first_a2_position in xrange(0, order - num_a2s + 1):
                    for alphas in product([X, Y, Z], repeat=order):
                        a1_alphas = tuple(3*a1_idx + alpha for alpha in alphas[:first_a2_position])
                        middle_alphas = alphas[first_a2_position:first_a2_position+num_a2s]
                        a2_alphas = tuple(3*a2_idx + alpha for alpha in middle_alphas)
                        a3_alphas = tuple(3*a3_idx + alpha
                                              for alpha in alphas[first_a2_position+num_a2s:])
                        val_to_spread = 0.0
                        for midatoms in product([a1_idx, a3_idx], repeat=num_a2s):
                            idxs = a1_alphas
                            idxs += tuple(ai*3 + middle_alphas[i] for i, ai in enumerate(midatoms))
                            idxs += a3_alphas
                            val_to_spread += my_b[idxs]
                        if num_a2s % 2 == 1:
                            val_to_spread *= -1.0
                        for perm in permutations(a1_alphas + a2_alphas + a3_alphas):
                            my_b[perm] = val_to_spread
        #--------------------------------------------------------------------------------#
        # Cache the value we got
        SimpleInternalCoordinate._set_b_tensor_cache_entry(my_b, *cache_key)
        #--------------------------------------------------------------------------------#
        return my_b


    ###################
    # Private Methods #
    ###################

    def _recompute_b_vector(self):
        b = Matrix(shape=(self.molecule.natoms, 3))
        i = self.atom_indices
        b[i[0]], b[i[1]], b[i[2]] = self.b_vector_for_positions(*[a.pos for a in self.atoms])
        self._bvector = b.flatten('C').view(Vector)


#####################
# Dependent Imports #
#####################

from grendel.chemistry.atom import Atom
from grendel.chemistry.molecule import Molecule
from grendel.representations.internal_representation import InternalRepresentation

