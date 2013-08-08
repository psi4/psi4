from collections import Iterable
from copy import copy
import numpy as np
from functools import partial
from itertools import product, chain, permutations, combinations_with_replacement as symmetric_product
from types import NoneType
from grendel import sanity_checking_enabled
import operator
from grendel.coordinates.coordinate import Coordinate
from grendel.coordinates.simple_internal_coordinate import SimpleInternalCoordinate
from grendel.differentiation.derivative_collection import DerivativeCollection
from grendel.external.combinatorics import unlabeled_balls_in_unlabeled_boxes, labeled_balls_in_unlabeled_boxes
from grendel.gmath import IndexRangeSet, DeclareIndexRange, IndexRange
from grendel.gmath.misc import delta
from grendel.gmath.tensor import Tensor, ComputableTensor
from grendel.representations.cartesian_representation import Z, Y, X
from grendel.representations.internal_representation import InternalRepresentation
from grendel.chemistry.atom import Atom
from grendel.gmath.matrix import Matrix
from grendel.gmath.vector import Vector, LightVector
from grendel.util.decorators import typechecked, IterableOf
from grendel.util.iteration import brace_notation_iter, ordered_partitions, ordered_partitions_iter, I_lubar
from grendel.util.overloading import overloaded
from grendel.util.sentinal_values import All
from grendel.util.units import Angstroms
from grendel.util.units.unit import DistanceUnit, isunit

__all__ = [
    "BondLength"
]

class BondLength(SimpleInternalCoordinate):
    """

    """

    ####################
    # Class Attributes #
    ####################

    default_delta = 0.01*Angstroms
    analytic_b_orders = All
    coordinate_symbol = "R"

    ############################
    # Private Class Attributes #
    ############################

    _btens_idx_range_set = IndexRangeSet()
    DeclareIndexRange('v', 6, index_range_set=_btens_idx_range_set).with_subranges(
        IndexRange('a', 0, 3, index_range_set=_btens_idx_range_set),
        IndexRange('b', 3, 6, index_range_set=_btens_idx_range_set)
    )

    ##################
    # Initialization #
    ##################

    @overloaded
    def __init__(self, *args, **kwargs):
        """
        BondLength(*args, **kwargs)

        **Signatures**
            * ``BondLength(atom1, atom2)``
            * ``BondLength(idx1, idx2, parent, one_based=True)``
            * ``BondLength(atoms, parent, index)``
            * ``BondLength(atoms, index)``

        :Parameters:

        atom1 : `Atom`
            The first atom that the bond length coordinate is composed of.
        atom2 : `Atom`
            The second atom that the bond length coordinate is composed of.
        idx1 : `int`
            The index of the first atom that the bond length coordinate is composed of.
        idx2: `int`
            The index of the second atom that the bond length coordinate is composed of.
        one_based: `bool`
            Whether the indices refer to atom numbers beginning at 1 (default) or 0.
        atoms: `list` of `Atom`
            List of the atoms composing the coordinate
        parent: `InternalRepresentation`
            The representation containing the coordinate `self`.
        index: `int`
            Index of `self` in the parent representation

        """
        raise TypeError

    @__init__.overload_with(
        atom1=Atom, atom2=Atom,
        parent=(InternalRepresentation, None),
        units=isunit)
    def __init__(self, atom1, atom2, parent=None, units=DistanceUnit.default, **kwargs):
        self.__init__(
            atoms=[atom1, atom2],
            parent=parent,
            units=units,
            **kwargs)

    @__init__.overload_with(
        idx1=int, idx2=int,
        parent=InternalRepresentation,
        one_based=bool,
        units=isunit)
    def __init__(self, idx1, idx2, parent, one_based=True, units=DistanceUnit.default, **kwargs):
        off = 1 if one_based else 0
        self._parent = parent
        self.__init__(
            self.parent_representation.molecule.atoms[idx1 - off],
            self.parent_representation.molecule.atoms[idx2 - off],
            units=units,
            **kwargs
        )

    @__init__.overload_with(
        idx1=int, idx2=int,
        molecule='Molecule',
        one_based=bool,
        units=isunit)
    def __init__(self, idx1, idx2, molecule, one_based=True, units=DistanceUnit.default, **kwargs):
        off = 1 if one_based else 0
        self.__init__(
            (molecule.atoms[idx1 - off], molecule.atoms[idx2 - off]),
            units=units,
            **kwargs
        )

    @__init__.overload_with(
        atoms=IterableOf(Atom),
        parent=(InternalRepresentation, None),
        index=(int, None), units=isunit)
    def __init__(self, atoms, parent=None, index=None, units=DistanceUnit.default, **kwargs):
        self._atoms = copy(atoms)
        if index is not None:
            self._index = index
        self._init(
            parent=parent,
            units=units,
            **kwargs
        )

    ###########
    # Methods #
    ###########

    @classmethod
    def value_for_xyz(cls, xyz):
        """
        """
        return LightVector.l_sub(xyz[1], xyz[0]).magnitude()

    @classmethod
    def b_vector_for_positions(cls, *args):
        if sanity_checking_enabled:
            if len(args) != 2:
                raise TypeError
        #----------------------------------------#
        def e(j, k):
            ev = LightVector.l_sub(args[k-1], args[j-1]); ev.normalize()
            return ev
        e21 = e(2, 1)
        return e21, -e21

    def analytic_b_tensor_for_order(self, order):
        # First check the cache
        cache_key = (self.__class__, order) + tuple(a.pos for a in self.atoms)
        cache_resp = SimpleInternalCoordinate._check_b_tensor_cache(*cache_key)
        if cache_resp is not None:
            return cache_resp
        #--------------------------------------------------------------------------------#
        # Define a function that acts like the B tensor for getting lower-order terms
        B = partial(SimpleInternalCoordinate.b_tensor_element_reindexed, self)
        #--------------------------------------------------------------------------------#
        # First order is trivial...
        if order == 1:
            # We can't use the `b_vector` attribute since that vector is indexed in the parent
            #   molecule's indexing scheme
            my_b = Vector(self.__class__.b_vector_for_positions(*[a.pos for a in self.atoms]))
            # return early to avoid double-caching, since b_vector_for_positions caches
            #   the result in the B tensor cache already
            return my_b
        #--------------------------------------------------------------------------------#
        # Second order terms
        elif order == 2:
            # BondLengths are composed of 2 atoms (3 CartesianCoordinates each), and the
            #   order is 2, so the output Tensor will be a 6x6 Tensor
            my_b = np.ndarray(shape=(6,)*2)
            Rabinv = 1.0 / self.value
            # This uses the itertools.product function, which acts like
            #   a nested loop (with depth given by the 'repeat' parameter).
            #   Mostly, I used this here just to keep the code indentation
            #   from getting out of hand (particularly in the case of the
            #   Torsion coordinate).  Also, though, this is useful for the
            #   general coordinate formulas.
            # First do the "same atom" derivatives
            for a_idx, a in enumerate(self.atoms):
                # Now iterate over the possible sets of cartesian coordinates
                for alpha, beta in symmetric_product([X, Y, Z], 2):
                    a_alpha, a_beta = 3*a_idx + alpha, 3*a_idx + beta
                    if alpha != beta:
                        my_b[a_alpha, a_beta] = - Rabinv * (B(self, a_alpha) * B(self, a_beta))
                        my_b[a_beta, a_alpha] = - Rabinv * (B(self, a_alpha) * B(self, a_beta))
                    else:
                        my_b[a_alpha, a_beta] = - Rabinv * (B(self, a_alpha) * B(self, a_beta) - 1)
            # Now get the "different atom" derivatives from the "same atom ones"
            for alpha in [X, Y, Z]:
                for beta in [X, Y, Z]:
                    a_alpha = alpha
                    b_beta = 3 + beta
                    my_b[a_alpha, b_beta] = my_b[b_beta, a_alpha] = (-1) * my_b[alpha, beta]
        #--------------------------------------------------------------------------------#
        else:
            # Behold, the general formula!
            # BondLengths are composed of 2 atoms (3 CartesianCoordinates each),
            #   so the output Tensor will be a 6x6x...x6 (`order`-dimensional) Tensor
            my_b = Tensor(
                indices=','.join('v'*order),
                index_range_set=self.__class__._btens_idx_range_set
            )
            # New, fancy way with einstein summation (if I can get it to work)
            #B = self._btensor
            #for o in xrange(1, order):
            #    t = B.for_order(o)
            #    if isinstance(t, ComputableTensor):
            #        t.fill()
            # First do the "same atom" derivatives
            #def idxs(letter, num): return tuple(letter + '_' + str(x) for x in xrange(num))
            #aa = idxs('a', order)
            #bb = idxs('b', order)
            #for s in I_lubar(order, 2, aa):
            #    my_b[aa] += B[s[0]] * B[s[1]]
            #for s in I_lubar(order, 2, bb):
            #    my_b[bb] += B[s[0]] * B[s[1]]
            #Rabinv = 1.0 / self.value
            #my_b[aa] *= -Rabinv
            #my_b[bb] *= -Rabinv
            ## Now get the "different atom" derivatives from the "same atom ones"
            #for nacarts in xrange(1, order):
            #    aprefactor =  -1.0 if (order - nacarts) % 2 == 1 else 1.0
            #    for a_b_set in set(permutations("a"*nacarts + "b"*(order-nacarts))):
            #        ab = tuple(map(lambda x: x[0] + "_" + str(x[1]), zip(a_b_set, xrange(order))))
            #        my_b[ab] = aprefactor * my_b[aa]
            # old way...
            Rabinv = 1.0 / self.value
            for a_idx, a in enumerate(self.atoms):
                # Now iterate over the possible sets of cartesian coordinates
                for alphas in symmetric_product([X, Y, Z], order):
                    a_alphas = tuple(3*a_idx + alpha for alpha in alphas)
                    cumval = 0.0
                    for a_alps, a_bets in I_lubar(order, 2, a_alphas):
                        cumval += B(self, *a_alps) * B(self, *a_bets)
                    val_to_spread = -Rabinv * cumval
                    # and spread over permutations
                    for idxs in permutations(a_alphas):
                        my_b[idxs] = val_to_spread
            # Now get the "different atom" derivatives from the "same atom ones"
            #   (From eq. 31 in Allen, et al. Mol. Phys. 89 (1996), 1213-1221)
            for alphas in product([X, Y, Z], repeat=order):
                for nacarts in xrange(1, order):
                    acoords = tuple(acart for acart in alphas[:nacarts])
                    bcoords = tuple(3 + bcart for bcart in alphas[nacarts:])
                    aprefactor =  -1.0 if (order - nacarts) % 2 == 1 else 1.0
                    val_to_spread = my_b[alphas]
                    for perm in permutations(acoords + bcoords):
                        my_b[perm] = aprefactor * val_to_spread
        #--------------------------------------------------------------------------------#
        # Whew...that was tough.  Let's cache the value we got
        SimpleInternalCoordinate._set_b_tensor_cache_entry(my_b, *cache_key)
        #--------------------------------------------------------------------------------#
        return my_b

    ###################
    # Private Methods #
    ###################

    def _recompute_b_vector(self):
        b = Matrix(shape=(self.molecule.natoms, 3))
        i = self.atom_indices
        b[i[0]], b[i[1]] = self.__class__.b_vector_for_positions(*[a.pos for a in self.atoms])
        self._bvector = b.flatten('C').view(Vector)


#####################
# Dependent Imports #
#####################

from grendel.chemistry.molecule import Molecule

