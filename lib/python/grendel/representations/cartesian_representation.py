from copy import deepcopy, copy
from itertools import product

import numpy as np

from grendel import type_checking_enabled, sanity_checking_enabled
import math
from grendel.external.combinatorics import prod, m_way_ordered_combinations
from grendel.gmath.geometry import angle_between_vectors, rotate_point_about_axis
from grendel.gmath.matrix import Matrix
from grendel.gmath.tensor import Tensor, chop, chopped
from grendel.gmath.vector import LightVector, cross, Vector, norm
from grendel.representations.representation import Representation
from grendel.util.decorators import see_abstract_doc, typechecked, IterableOf
from grendel.util.freezing import Freezable
from grendel.util.iteration import grouper, I_lubar
from grendel.util.strings import shortstr, indent
from grendel.util.units import UnknownUnitError, IncompatibleUnitsError, isunit, DistanceUnit, Radians, Degrees, AngularUnit


__all__ = [
    "CartesianRepresentation"
]

X, Y, Z = 0, 1, 2

class CartesianRepresentation(Representation):
    """ A representations of a molecule in Cartesian coordinates.
    Since Atom objects already have Cartesian vectors
    associated with them, this class mostly consists of helper functions and convenience (accessor) functions that make
    CartesianRepresentation compatible with the parent Representation class.
    """

    ##################
    # Initialization #
    ##################

    @typechecked(
        molecule='Molecule',
        coords=(None, IterableOf('CartesianCoordinate')),
        units=(None, isunit)
    )
    def __init__(self, molecule, coords=None, units=None):
        """
        CartesianRepresentation(molecule, coords)

        Parameters
        ----------
        molecule : `Molecule`
            The molecule that the representation refers to.
        coords : list of `CartesianCoordinate`
            The coordinates represented, or `None` if the coordinates are to be extracted automatically
            from the current cartesian position of `molecule`

        """
        self.molecule = molecule
        if coords is None:
            self.coords = []
        else:
            self.coords = []
            for coord in coords:
                self.add_coordinate_copy(coord)
        if units is not None:
            if sanity_checking_enabled and not units.genre is DistanceUnit:
                raise IncompatibleUnitsError("{} is not a unit of distance.".format(units))
            self._units = units
        else:
            self._units = DistanceUnit.default

    ##############
    # Properties #
    ##############

    @property
    def atoms(self):
        """ Convenience property for referencing the atoms list in the molecule associated with self.
        """
        return self.molecule.atoms


    @property
    def l_xyz(self):
        """
        Returns the value of the representation as a light vector
        """
        return LightVector([c.value for c in self])

    ###################
    # Special Methods #
    ###################

    #------------------------#
    # Output Representations #
    #------------------------#

    def __str__(self):
        ret_val = "CartesianRepresentation of molecule {}, in {}:\n".format(
            shortstr(self.molecule), self.units)
        for i, (x, y, z) in enumerate(grouper(3, self)):
            ret_val += indent('{:>9} {:12.8f} {:12.8f} {:12.8f}\n'.format(
                '[{}-{}]'.format(3*i, 3*i+2),
                x.value, y.value, z.value
            ), 2)
        return ret_val

    ##################
    # Static Methods #
    ##################

    @staticmethod
    def same_atom_indices(index):
        atom_num = (int(index) / 3)
        return [3*atom_num, 3*atom_num + 1, 3*atom_num+2]

    ###########
    # Methods #
    ###########

    def atom(self, number):
        """ Convenience method for accessing atoms.  Returns self.molecule.atoms[number]
        """
        return self.atoms[number]

    def index(self, atom_or_coordinate):
        """ Returns the index of the parameter `atom` in the representation self.

        Parameters
        ----------
        atom_or_coordinate : `Atom` or `CartesianCoordinate`
            The atom or cartesian coordinate to obtain the index of.  If an atom is passed in, the index of the X coordinate of the atom is returned.

        Raises
        ------
        IndexError
            If the parameter `atom_or_coordinate` is not found.

        """
        if isinstance(atom_or_coordinate, Atom):
            if atom_or_coordinate not in self.molecule.atoms:
                raise IndexError(repr(atom_or_coordinate) + " not found in " + repr(self))
            return self.molecule.atoms.index(atom_or_coordinate) * 3
        elif isinstance(atom_or_coordinate, CartesianCoordinate):
            if atom_or_coordinate not in self.coords:
                raise IndexError(repr(atom_or_coordinate) + " not found in " + repr(self))
            return atom_or_coordinate.index
        else:
            raise TypeError


    @typechecked(coordinate='CartesianCoordinate')
    def add_coordinate_copy(self, coordinate, **kwargs):
        """ Add a CartesianCoordinate to the representation.
        """
        self.fail_if_frozen()
        new_coord = coordinate.copy_for_representation(
            rep=self,
            index=len(self.coords),
            **kwargs
        )
        self.coords.append(new_coord)

    def add_atom(self, atom):
        """ Add an atom to the representation by creating a CartesianCoordinate object cooresponding to the atom.
        """
        self.fail_if_frozen()
        if type_checking_enabled and not isinstance(atom, Atom):
            raise TypeError
        for direction in [0, 1, 2]:
            addval = CartesianCoordinate(
                atom=atom,
                direction=direction,
                parent=self)
            self.coords.append(addval)

    @typechecked(atom='Atom')
    def refresh_atom(self, atom):
        """ Renew the coordinates for an atom to coincide with the atom's current position
        """
        self.fail_if_frozen()
        spot = atom.index
        for direction in [0, 1, 2]:
            addval = CartesianCoordinate(
                atom=atom,
                direction=direction,
                parent=self)
            self.coords[spot*3 + direction] = addval

    @see_abstract_doc
    def copy_with_molecule(self, molecule):
        # TODO copy more stuff
        ret_val = self.__class__(
            molecule,
            units=self.units)
        for coord in self:
            ret_val.add_coordinate_copy(coord)
        molecule.cartesian_representation = ret_val
        return ret_val

    @see_abstract_doc
    def displaced_by(self, disp, tol=None, maxiter=None):
        # tol and maxiter are irrelevent here...
        dispvect = disp.desired_values - self.values
        disp_mol = deepcopy(self.molecule)
        newrep = self.copy_with_molecule(disp_mol)
        disp_mol.displace(dispvect)
        return disp_mol, newrep

    @typechecked(
        tensor='RepresentationDependentTensor',
        to_representation=Representation)
    def transform_tensor(self, tensor, to_representation):
        """
        """
        self.freeze()
        to_representation.freeze()
        shape=(len(to_representation),)*len(tensor.shape)
        ret_val = RepresentationDependentTensor(
            shape=shape,
            representation=to_representation)
        #--------------------------------------------------------------------------------#
        if isinstance(to_representation, CartesianRepresentation):
            if len(self) != len(to_representation):
                raise ValueError("incompatible representation sizes ({} != {})".format(
                    len(self),
                    len(to_representation)
                ))
            if tensor.representation is not self:
                raise ValueError("representation {} can only transform a tensor whose representation attribute is the same "
                                 " as itself (tensor.representation was {} ".format(self, tensor.representation
                ))
            if self.molecule.is_linear():
                raise NotImplementedError("linear molecule cartesian-to-cartesian transformation is not"
                                          " yet implemented.  Shouldn't be too hard...")
            if sanity_checking_enabled:
                #TODO use a recentered version of the 'from' representation if it is not centered
                pass
                # Make sure things are centered
                #if not self.molecule.is_centered(cartesian_representation=self):
                #    raise ValueError("CartesianRepresentation objects transformed from and to"
                #                     " must be centered at the center of mass ({} was not)".format(
                #        self
                #    ))
                #elif not self.molecule.is_centered(cartesian_representation=to_representation):
                #    raise ValueError("CartesianRepresentation objects transformed from and to"
                #                     " must be centered at the center of mass ({} was not)".format(
                #        to_representation
                #    ))
            #----------------------------------------#
            old = self.value
            new = to_representation.value
            unitconv = self.units.to(to_representation.units)
            oldmat = Matrix(old).reshape((len(self)/3, 3))
            newmat = Matrix(new).reshape((len(self)/3, 3))
            #----------------------------------------#
            # Check to see if the molecule is planar.  If so, append the cross product of any
            #   two atoms' positions not colinear with the origin
            if any(norm(oldmat[:dir].view(Vector)) < 1e-12 \
                    or norm(newmat[:dir].view(Vector)) < 1e-12 \
                        for dir in [X, Y, Z]):
                first_atom = None
                # TODO unit concious cutoff (e.g. this would fail if the units of the position matrix were meters)
                nonorigin_cutoff = 1e-2
                for i, v in enumerate(grouper(3, old)):
                    if norm(LightVector(v)) > nonorigin_cutoff:
                        first_atom = i
                        break
                second_atom = None
                for i in range(first_atom+1, len(self)//3):
                    #TODO make sure this works with Molecule.linear_cutoff
                    if norm(oldmat[i]) > nonorigin_cutoff and norm(newmat[i]) > nonorigin_cutoff \
                            and abs(angle_between_vectors(oldmat[first_atom], oldmat[i])) > 1e-5 \
                            and abs(angle_between_vectors(newmat[first_atom], newmat[i])) > 1e-5:
                        second_atom = i
                        break
                oldmat = Matrix(list(oldmat.rows) + [cross(oldmat[first_atom], oldmat[second_atom])])
                newmat = Matrix(list(newmat.rows) + [cross(newmat[first_atom], newmat[second_atom])])
            #----------------------------------------#
            rot_mat = newmat.T * np.linalg.pinv(oldmat.T)
            # Divide by unit conversion because the representation dependent tensor
            #   is [some units] *per* representation units
            rot_mat /= unitconv
            trans_mat = Matrix(shape=(len(self), len(self)))
            for idx in xrange(len(self)//3):
                off = idx * 3
                trans_mat[off:off+3,off:off+3] = rot_mat
            order = len(tensor.shape)
            # This is basically impossible to read what I'm doing here, but work it out
            #   looking at the numpy.einsum documentation and you'll see that this is correct.
            #   Basically, we want to contract the row indices of the transformation matrix
            #   with each axis of the tensor to be transformed.
            einsumargs = sum(([trans_mat, [2*i, 2*i+1]] for i in xrange(order)), [])
            einsumargs += [tensor, [i for i in xrange(2*order) if i % 2 != 0]]
            einsumargs += [[i for i in xrange(2*order) if i % 2 == 0]]
            np.einsum(*einsumargs, out=ret_val)
            return ret_val
        #--------------------------------------------------------------------------------#
        elif isinstance(to_representation, InternalRepresentation):
            if len(tensor.shape) == 1:
                A = to_representation.a_matrix
                ret_val[...] = A.T * tensor.view(Vector)
            else:
                raise NotImplementedError("use transform_forcefield instead")
            return ret_val
        #--------------------------------------------------------------------------------#
        elif isinstance(to_representation, NormalRepresentation):
            B = to_representation.b_matrix
            ret_val[...] = tensor.linearly_transformed(B)
            return ret_val
        else:
            raise NotImplementedError(
                "Transformation of arbitrary tensors from representation of type '{}' to "
                " representations of type '{}' is not implemented.".format(
                    self.__class__.__name__,
                    to_representation.__class__.__name__
                )
            )

    @typechecked(
        ff='ForceField',
        to_representation=Representation)
    def transform_forcefield(self, ff, to_representation):
        self.freeze()
        to_representation.freeze()
        if sanity_checking_enabled and ff.representation is not self:
            raise ValueError("cart_rep.transform_forcefield must be given a force field whose"
                             " representation is cart_rep")
        rv = ForceField(to_representation, ff.max_order, ff.molecular_property, ff.property_units)
        if isinstance(to_representation, CartesianRepresentation):
            for order in xrange(1, ff.max_order+1):
                rv[...] = self.transform_tensor(ff.for_order(order), to_representation)
            return rv
        elif isinstance(to_representation, NormalRepresentation):
            for order in xrange(1, ff.max_order+1):
                rv[...] = self.transform_tensor(ff.for_order(order), to_representation)
            return rv
        elif isinstance(to_representation, InternalRepresentation):
            # Units should take care of themselves here.  This is because
            #   the B tensor should have units (to-representation units)/(cart-representation units)**order
            # TODO write a test that checks the fact that units take care of themselves
            A = DerivativeCollection(
                representation=self,
                first_dimension_different=True,
                secondary_representation=to_representation,
                uncomputed=False
            )
            A.for_order(1)[...] = to_representation.a_matrix
            B = to_representation.b_tensor(cartesian_representation=self)
            rv.for_order(1)[...] = A.for_order(1).view(Matrix).T * ff.for_order(1).view(Vector)
            if ff.max_order >= 2:
                C = DerivativeCollection(representation=to_representation, uncomputed=False)
                C['q,p1,p2'] = B['q,i1,i2'] * A['i1,p1'] * A['i2,p2']
                rv['p1,p2'] = ff['i1,i2']*A['i1,p1']*A['i2,p2'] - rv['q']*C['q,p1,p2']
                #--------------------------------------------------------------------------------#
                if ff.max_order >= 3:
                    # Do third order by hand...
                    C['q,p1,p2,p3'] = B['q,i1,i2,i3'] * A['i1,p1'] * A['i2,p2'] * A['i3,p3']
                    rv['p1,p2,p3'] = ff['i1,i2,i3']*A['i1,p1']*A['i2,p2']*A['i3,p3']
                    rv['p1,p2,p3'] -= rv['q,p1']*C['q,p2,p3']
                    rv['p1,p2,p3'] -= rv['q,p2']*C['q,p1,p3']
                    rv['p1,p2,p3'] -= rv['q,p3']*C['q,p1,p2']
                    rv['p1,p2,p3'] -= rv['q']*C['q,p1,p2,p3']
                    if ff.max_order >= 4:
                        def idxstrs(letter, num): return tuple(map(lambda x: letter + str(x), xrange(num)))
                        for order in xrange(4, ff.max_order+1):
                            p_s = idxstrs('p', order)
                            i_s = idxstrs('i', order)
                            C[('q',)+p_s] = B[('q',)+i_s] * prod(A[i, p] for i, p in zip(i_s, p_s))
                            rv[p_s] = ff[i_s] * prod(A[i, p] for i, p in zip(i_s, p_s))
                            rv[p_s] -= rv['q'] * C[('q',)+p_s]
                            for k in xrange(2, order):
                                q_s = idxstrs('q', k)
                                for s in I_lubar(order, k, p_s):
                                    if any(len(seg) == 1 for seg in s):
                                        continue
                                    rv[p_s] -= rv[q_s] * prod(C[(q_s[m],) + s[m]] for m in xrange(k))
                                q_s = idxstrs('q', k-1)
                                for part1len in xrange(1, order - 2*(k-1) + 1):
                                    part2len = order - part1len
                                    for p1idxs, p2idxs in m_way_ordered_combinations(order, [part1len, part2len]):
                                        p1s = tuple(p_s[pidx] for pidx in p1idxs)
                                        p2s = tuple(p_s[pidx] for pidx in p2idxs)
                                        print(p1s, p2s, k)
                                        for s in I_lubar(part2len, k-1, p2s):
                                            rv[p_s] -= rv[q_s + p1s] * prod(
                                                C[(alpha,)+betas] for alpha, betas in zip(q_s, s)
                                            )
            return rv
        else:
            raise NotImplementedError


    def atom_coords(self, atom_or_index):
        if isinstance(atom_or_index, Atom):
            atom_or_index = atom_or_index.index()
        return self[atom_or_index*3], self[atom_or_index*3+1], self[atom_or_index*3+2]

    def iter_atom_coords(self, with_atom=False):
        for i, (a, b, c) in enumerate(grouper(3, self)):
            if with_atom:
                yield a, b, c, self.atoms[i]
            else:
                yield a, b, c

    def frozen_copy(self):
        new_rep = CartesianRepresentation(molecule=self.molecule, units=self.units)
        for coord in self:
            new_rep.add_coordinate_copy(
                coord,
                freeze_value=True,
                value=coord.value
            )
        new_rep.freeze()
        return new_rep

#####################
# Dependent Imports #
#####################

from grendel.chemistry.atom import Atom
from grendel.representations.internal_representation import InternalRepresentation
from grendel.differentiation.derivative_tensor import RepresentationDependentTensor
from grendel.differentiation.derivative_collection import DerivativeCollection
from grendel.coordinates.cartesian_coordinate import CartesianCoordinate
from grendel.differentiation.forcefield import ForceField
from grendel.representations.normal_representation import NormalRepresentation

if type_checking_enabled:
    from grendel.chemistry.molecule import Molecule

