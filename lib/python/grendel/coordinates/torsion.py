from functools import partial
from itertools import product, permutations, islice, combinations_with_replacement as symmetric_product
from operator import mul
from math import cos, sin, acos, atan
import math
import numpy as np
from grendel import sanity_checking_enabled
from grendel.coordinates.bond_angle import BondAngle
from grendel.coordinates.bond_length import BondLength
from grendel.coordinates.periodic_coordinate import PeriodicCoordinate

from grendel.gmath.geometry import angle_between_vectors, safe_acos, safe_asin
from grendel.coordinates.simple_internal_coordinate import SimpleInternalCoordinate
from grendel.gmath.matrix import Matrix
from grendel.gmath.tensor import Tensor
from grendel.gmath.vector import cross, Vector, l_cross, LightVector
from grendel.external.combinatorics import labeled_balls_in_labeled_boxes, labeled_balls_in_unlabeled_boxes, unlabeled_balls_in_unlabeled_boxes, unlabeled_balls_in_labeled_boxes
from grendel.representations.cartesian_representation import Z, Y, X
from grendel.util.decorators import IterableOf, SequenceOf
from grendel.util.iteration import ordered_partitions_iter, brace_notation_iter, all_partitions_iter, unique_permutations, I_ll, I_lubar, I_llbar
from grendel.util.overloading import overloaded
from grendel.util.sentinal_values import All
from grendel.util.strings import indented
from grendel.util.units import Degrees, Radians, AngularUnit
from grendel.util.units.unit import isunit

__all__ = [
    "Torsion"
]

class Torsion(PeriodicCoordinate):
    """

    """

    ####################
    # Class Attributes #
    ####################

    default_delta = 0.02*Radians
    analytic_b_orders = All
    preferred_range = [-math.pi, math.pi]
    coordinate_symbol = 'tau'

    ##################
    # Initialization #
    ##################

    @overloaded
    def __init__(self, *args, **kwargs):
        raise TypeError

    @__init__.overload_with(atom1='Atom', atom2='Atom', atom3='Atom', atom4='Atom')
    def __init__(self, atom1, atom2, atom3, atom4, **kwargs):
        self.__init__([atom1, atom2, atom3, atom4], **kwargs)

    @__init__.overload_with(
        idx1=int, idx2=int, idx3=int, idx4=int,
        parent=('InternalRepresentation', None))
    def __init__(self, idx1, idx2, idx3, idx4, parent, one_based=True, **kwargs):
        off = 1 if one_based else 0
        self.__init__(
            (parent.molecule.atoms[i - off] for i in [idx1, idx2, idx3, idx4]),
            parent=parent,
            **kwargs
        )

    @__init__.overload_with(
        idx1=int, idx2=int, idx3=int, idx4=int,
        molecule='Molecule',
        one_based=bool,
        units=isunit)
    def __init__(self, idx1, idx2, idx3, idx4, molecule, one_based=True, units=AngularUnit.default, **kwargs):
        off = 1 if one_based else 0
        self.__init__(
            (molecule.atoms[idx1-off], molecule.atoms[idx2-off], molecule.atoms[idx3-off], molecule.atoms[idx4-off]),
            units=units,
            **kwargs
        )

    # Principal initializer (all other __init__ methods call this one)
    @__init__.overload_with(
        atoms=IterableOf('Atom'),
        parent=('InternalRepresentation', None),
        index=(int, None),
        units=isunit)
    def __init__(self, atoms, parent=None, index=None, units=AngularUnit.default, **kwargs):
        self._atoms = list(atoms)
        self._index = index
        self._init(
            parent=parent,
            units=units,
            **kwargs
        )

    ##############
    # Properties #
    ##############

    @property
    def terminal_atoms(self):
        return self.atoms[0], self.atoms[3]

    #################
    # Class Methods #
    #################

    @classmethod
    def noncanonical_value_for_xyz(cls, xyz):
        def e(j, k):
            ev = LightVector.l_sub(xyz[k-1], xyz[j-1]); ev.normalize()
            return ev
        e12 = e(1,2)
        e21 = -e12
        e32 = e(3,2)
        e23 = -e32
        e34 = e(3,4)
        phi2 = angle_between_vectors(e21, e23)
        phi3 = angle_between_vectors(e32, e34)
        sinphi2, sinphi3 = sin(phi2), sin(phi3)
        sintau = e21.dot(cross(e32, e34)) / (sinphi2 * sinphi3)
        tau_0 = safe_asin(sintau)
        #----------------------------------------#
        # Deal with the sign
        sign_val = cross(e21, e23).dot(cross(e32, e34))
        #zero_cutoff = 1e-10
        zero_cutoff = 0
        if sign_val < -zero_cutoff:
            tau_e = math.pi - tau_0
            #tau_e = - tau_0
        elif sign_val > zero_cutoff:
            tau_e = tau_0
        else:
            tau_e = math.pi / 2.0
        #----------------------------------------#
        return tau_e

    @classmethod
    def b_vector_for_positions(cls, *args):
        if sanity_checking_enabled:
            if len(args) != 4:
                raise TypeError
        #--------------------------------------------------------------------------------#
        # Check the cache first
        cache_resp = SimpleInternalCoordinate._check_b_tensor_cache(cls, args)
        if cache_resp:
            return cache_resp
        #--------------------------------------------------------------------------------#
        # Not in the cache...so compute it
        def e(j, k):
            ev = LightVector.l_sub(args[k-1], args[j-1]); ev.normalize()
            return ev
        def r(j, k):
            return LightVector.l_sub(args[k-1], args[j-1]).magnitude()
        e12, e23, e43 = e(1,2), e(2,3), e(4,3)
        e32 = -e23
        r12, r23, r43 = r(1,2), r(2, 3), r(3, 4)
        r32 = r23
        sin2phi2 = sin(angle_between_vectors(e12, e32))**2
        cosphi2 = cos(angle_between_vectors(e12, e32))
        sin2phi3 = sin(angle_between_vectors(e23, e43))**2
        cosphi3 = cos(angle_between_vectors(e23, e43))
        ret_val = [0,0,0,0]
        ret_val[0] = cross(e12, e23) * (-1.0 / (r12 * sin2phi2))
        ret_val[1] = ((r23 - r12 * cosphi2)/(r23*r12*sin2phi2)) * cross(e12, e23)\
                        + (cosphi3/(r23*sin2phi3)) * cross(e43, e32)
        ret_val[2] = ((r32 - r43 * cosphi3)/(r32*r43*sin2phi3)) * cross(e43, e32)\
                        + (cosphi2/(r32*sin2phi2)) * cross(e12, e23)
        ret_val[3] = cross(e43, e32) * (-1.0 / (r43 * sin2phi3))
        SimpleInternalCoordinate.b_tensor_cache[(cls,) + tuple(args)] = ret_val
        return ret_val

    ###########
    # Methods #
    ###########

    def analytic_b_tensor_for_order(self, order):
        # First check the cache
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
        # TODO Explicit (matrix/tensor based) implementations of 2nd, 3rd, and 4th order to speed things up substantially
        else:
            # We only have general formulas for terminal atoms as of now...
            # So we can save a little bit of time by computing these terms and
            #   then doing only finite difference for everything else
            # Torsions are composed of 4 atoms (3 CartesianCoordinates each),
            #   so the output will be a 12x12x...x12 (`order`-dimensional) tensor
            my_b = np.ndarray(shape=(12,)*order)
            if sanity_checking_enabled:
                my_b = np.ones(shape=(12,)*order) * float('inf')
            # some precomputed values
            #========================================#
            # Helper function needed for derivative:
            def Bsin2phi(phi, *idx_alphas):
                # some precomputed values
                twophi = 2.0*phi.value*phi.units.to(Radians)
                sin2phi = sin(twophi)
                cos2phi = cos(twophi)
                def h_K(k):
                    if k % 2 == 1:
                        if ((k-1) / 2) % 2 == 0:
                            return cos2phi
                        else:
                            return -cos2phi
                    else: # k is even
                        if k/2 % 2 == 0:
                            return sin2phi
                        else:
                            return -sin2phi
                cumsum = 0.0
                n = len(idx_alphas)
                for k in xrange(1, n + 1):
                    inner_sum = 0.0
                    for s in I_lubar(n, k, idx_alphas):
                        inner_prod = 1.0
                        for i in range(k):
                            inner_prod *= B(phi, *s[i])
                        inner_sum += inner_prod
                    cumsum += 2.0**k * h_K(k) * inner_sum
                return cumsum
            #========================================#
            # The aaaa...bbbbb...cccc... terms
            phis = []
            rbcs = []
            for a_idx, a in zip([0, 3], self.terminal_atoms):
                # Determine whihc terminal atom we're handling
                if a_idx == 0:
                    b_idx, c_idx, d_idx = 1, 2, 3
                else: # a_idx == 3
                    b_idx, c_idx, d_idx = 2, 1, 0
                #----------------------------------------#
                phi_abc = self.get_coord(BondAngle,
                    self.atoms[a_idx],
                    self.atoms[b_idx],
                    self.atoms[c_idx],
                )
                # units of Radians
                ang_conv = phi_abc.units.to(Radians)
                # TODO figure out what 'units' should be here
                rbc = self.get_coord(BondLength, self.atoms[b_idx], self.atoms[c_idx])
                # Keep them for later
                phis.append(phi_abc)
                rbcs.append(rbc)
                #----------------------------------------#
                # some precomputed values
                sin2phi = sin(2.0 * phi_abc.value * ang_conv)
                #----------------------------------------#
                for num_a, num_b, num_c in unlabeled_balls_in_labeled_boxes(order-1, [order-1]*3):
                    num_a += 1
                    alph_iter = product(*map(lambda n: symmetric_product([X,Y,Z], n), (num_a, num_b, num_c)))
                    for alphas_a, alphas_b, alphas_c in alph_iter:
                        a_alphas = tuple(3*a_idx + alpha for alpha in alphas_a)
                        b_alphas = tuple(3*b_idx + alpha for alpha in alphas_b)
                        c_alphas = tuple(3*c_idx + alpha for alpha in alphas_c)
                        all_alphas = a_alphas + b_alphas + c_alphas
                        cum_sum = 0.0
                        for t1, t2 in I_ll(num_b + num_c, 2, b_alphas + c_alphas):
                            bphi = LightVector([
                                B(phi_abc, 3*a_idx + sigma, *(a_alphas[1:] + t1))
                                    for sigma in [X, Y, Z]
                            ])
                            brbc = LightVector([
                                B(rbc, 3*c_idx + sigma, *t2) for sigma in [X, Y, Z]]
                            )
                            cum_sum += cross(bphi, brbc)[alphas_a[0]]
                        cum_sum *= 2.0
                        cum_sum -= B(self, a_alphas[0]) * Bsin2phi(phi_abc, *all_alphas[1:])
                        for s1, s2 in I_llbar(num_a - 1 + num_b + num_c, 2, all_alphas[1:]):
                            cum_sum -= Bsin2phi(phi_abc, *s1) * B(self, *((a_alphas[0],) + s2))
                        cum_sum /= sin2phi
                        for perm in permutations(all_alphas):
                            my_b[perm] = cum_sum
            #========================================#
            # Note that the terminal-atom cross-derivatives are 0
            if sanity_checking_enabled:
                # Fill in the explicity zeros now, since we had infinity there to make sure
                #   uncomputed values weren't being used for something else.
                for a_idx, d_idx in permutations([0, 3]):
                    for remaining_idxs in product([0, 1, 2, 3], repeat=order-2):
                        for alphas in product([X, Y, Z], repeat=order):
                            idx_alphas = tuple(3*atom_idx + alpha
                                    for atom_idx, alpha in zip((a_idx, d_idx) + remaining_idxs, alphas))
                            for perm in permutations(idx_alphas):
                                my_b[perm] = 0.0
            #========================================#
            # Arbitrary order bbbbb.... terms
            a_idx, b_idx, c_idx, d_idx = 0, 1, 2, 3
            phi_abc = phis[0]
            phi_bcd = self.get_coord(BondAngle,
                self.atoms[b_idx],
                self.atoms[c_idx],
                self.atoms[d_idx],
            )
            ang_conv = phi_bcd.units.to(Radians)
            # TODO figure out what 'units' should be here
            r_ba = self.get_coord(BondLength,
                self.atoms[b_idx],
                self.atoms[a_idx]
            )
            r_cd = self.get_coord(BondLength,
                self.atoms[c_idx],
                self.atoms[d_idx]
            )
            #----------------------------------------#
            def Bcscphi(phi, *b_alphas):
                phi_val = phi.value * phi.units.to(Radians)
                sinphi = sin(phi_val)
                cscphi = 1.0 / sinphi
                if len(b_alphas) == 0:
                    return cscphi
                cotphi = cos(phi_val) / sinphi
                #------------------------------------#
                def dcsc_n(n):
                    def t(n_t, k):
                        if k == 0:
                            return 1
                        elif k <= n_t//2:
                            return (2*k + 1) * t(n_t-1, k) + (n_t - 2*k + 1) * t(n_t-1, k-1)
                        else:
                            return 0
                    #--------------------------------#
                    ret_val = 0.0
                    for kk in xrange(n//2 + 1):
                        ret_val += t(n, kk) * cotphi**(n - 2*kk) * cscphi**(2*kk + 1)
                    if n % 2 == 1:
                        return -ret_val
                    else:
                        return ret_val
                #------------------------------------#
                outer_sum = 0.0
                for k in xrange(1, len(b_alphas) + 1):
                    inner_sum = 0.0
                    for idx_sets in labeled_balls_in_unlabeled_boxes(len(b_alphas), [len(b_alphas)]*k):
                        if any(len(st) == 0 for st in idx_sets):
                            continue
                        b_idx_sets = tuple(tuple(b_alphas[i] for i in idxset) for idxset in idx_sets)
                        product = 1.0
                        for b_idxs in b_idx_sets:
                            product *= B(phi, *b_idxs)
                        inner_sum += product
                    outer_sum += dcsc_n(k) * inner_sum
                return outer_sum
            #----------------------------------------#
            for alphas in symmetric_product([X, Y, Z], order):
                # here we go...
                term_sum = first_term = second_term = third_term = 0.0
                iter_alphas = alphas[1:]
                b_alphas = tuple(3*b_idx + alpha for alpha in alphas)
                #----------------------------------------#
                for b_alphas1, b_alphas2, b_alphas3 in I_ll(order-1, 3, b_alphas[1:]):
                    #----------------------------------------#
                    # preconstruct the vectors we need for the cross products
                    b_a_phi_abc = LightVector([
                        B(phi_abc, 3*a_idx + sigma, *b_alphas2) for sigma in [X, Y, Z]
                    ])
                    b_c_phi_abc = LightVector([
                        B(phi_abc, 3*c_idx + sigma, *b_alphas2) for sigma in [X, Y, Z]
                    ])
                    b_a_r_ba = LightVector([
                        B(r_ba, 3*a_idx + sigma, *b_alphas3) for sigma in [X, Y, Z]
                    ])
                    #----------------------------------------#
                    # save a little bit of time by only computing this part once per iteration
                    Bcscphi_abc = Bcscphi(phi_abc, *b_alphas1)
                    #----------------------------------------#
                    # now add the contribution from this set of indices
                    first_term += Bcscphi_abc * cross(b_a_phi_abc, b_a_r_ba)[alphas[0]]
                    second_term += Bcscphi_abc * cross(b_c_phi_abc, b_a_r_ba)[alphas[0]]
                #----------------------------------------#
                term_sum -= first_term + second_term
                b_d_r_cd = LightVector([
                    B(r_cd, 3*d_idx + sigma) for sigma in [X, Y, Z]
                ])
                for b_alphas1, b_alphas2 in I_ll(order-1, 2, b_alphas[1:]):
                    #----------------------------------------#
                    b_b_phi_bcd = LightVector([
                        B(phi_bcd, 3*b_idx + sigma, *b_alphas2) for sigma in [X, Y, Z]
                    ])
                    #----------------------------------------#
                    third_term += Bcscphi(phi_bcd, *b_alphas1) * cross(b_b_phi_bcd, b_d_r_cd)[alphas[0]]
                term_sum += third_term
                #----------------------------------------#
                # and spread it across permutations
                for perm in permutations(tuple(3*b_idx + alpha for alpha in alphas)):
                    my_b[perm] = term_sum
            #========================================#
            # and fill in the bbbb...cccc... derivatives by translational invariance
            # Range from one c index to all
            for num_c in xrange(1, order+1):
                num_b = order - num_c
                alph_iter = product(*map(lambda n: symmetric_product([X,Y,Z], n), (num_b, num_c)))
                for alphas_b, alphas_c in alph_iter:
                    b_alphas = tuple(3*b_idx + alph for alph in alphas_b)
                    c_alphas = tuple(3*c_idx + alph for alph in alphas_c)
                    alphas_all = alphas_b + alphas_c
                    currsum = 0.0
                    for repl_atoms in product([a_idx, b_idx, d_idx], repeat=num_c):
                        repl_alphas = tuple(3*repl_atom + alph
                            for repl_atom, alph in zip(repl_atoms, alphas_c))
                        currsum += my_b[repl_alphas + b_alphas]
                        if sanity_checking_enabled and math.isinf(my_b[b_alphas + repl_alphas]):
                            raise IndexError("indices not filled in: {}, needed for {}".format(
                                b_alphas + repl_alphas,
                                b_alphas + c_alphas
                            ))
                    if num_c % 2 == 1:
                        currsum *= -1.0
                    #----------------------------------------#
                    # and spread this value over all permutations...
                    for perm in permutations(b_alphas + c_alphas):
                        my_b[perm] = currsum
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
        b[i[0]], b[i[1]], b[i[2]], b[i[3]] = self.__class__.b_vector_for_positions(*[a.pos for a in self.atoms])
        self._bvector = b.flatten('C').view(Vector)

#####################
# Dependent Imports #
#####################

from grendel.chemistry.atom import Atom
from grendel.chemistry.molecule import Molecule
from grendel.representations.internal_representation import InternalRepresentation

