from copy import copy
from grendel.chemistry.atom import Atom
from grendel.coordinates.bond_angle import BondAngle
from grendel.coordinates.bond_length import BondLength
from grendel.coordinates.simple_internal_coordinate import SimpleInternalCoordinate
from grendel.gmath.vector import LightVector, Vector
from grendel.gmath.matrix import Matrix
from grendel.representations.cartesian_representation import Y, Z, X
from grendel.representations.internal_representation import InternalRepresentation
from grendel.util.decorators import IterableOf
from grendel.util.overloading import overloaded
from grendel.util.units import isunit, DistanceUnit, Angstroms

def generate_fixed_fragment_coords(a, b, c, *others):
    """
    TODO Document this
    """
    ret_val = [BondLength(a, b), BondLength(a, c), BondAngle(b, a, c)]
    for d in others:
        ret_val.append(InternalCartesianX(a, b, c, d))
        ret_val.append(InternalCartesianY(a, b, c, d))
        ret_val.append(InternalCartesianZ(a, b, c, d))
    return ret_val

class InternalCartesian(SimpleInternalCoordinate):

    ####################
    # Class Attributes #
    ####################

    default_delta = 0.01*Angstroms

    ##################
    # Initialization #
    ##################

    @overloaded
    def __init__(self, *args, **kwargs):
        """
        """
        raise TypeError

    @__init__.overload_with(
        atom1=Atom, atom2=Atom, atom3=Atom, atom4=Atom,
        parent=(InternalRepresentation, None),
        units=isunit)
    def __init__(self, atom1, atom2, atom3, atom4, parent=None, units=DistanceUnit.default, **kwargs):
        self.__init__(
            atoms=[atom1, atom2, atom3, atom4],
            parent=parent,
            units=units,
            **kwargs)

    @__init__.overload_with(
        idx1=int, idx2=int, idx3=int, idx4=int,
        parent=InternalRepresentation,
        one_based=bool,
        units=isunit)
    def __init__(self, idx1, idx2, idx3, idx4, parent, one_based=True, units=DistanceUnit.default, **kwargs):
        off = 1 if one_based else 0
        self._parent = parent
        self.__init__(
            self.parent_representation.molecule.atoms[idx1 - off],
            self.parent_representation.molecule.atoms[idx2 - off],
            self.parent_representation.molecule.atoms[idx3 - off],
            self.parent_representation.molecule.atoms[idx4 - off],
            units=units,
            **kwargs
        )

    @__init__.overload_with(
        idx1=int, idx2=int, idx3=int, idx4=int,
        molecule='Molecule',
        one_based=bool,
        units=isunit)
    def __init__(self, idx1, idx2, idx3, idx4, molecule, one_based=True, units=DistanceUnit.default, **kwargs):
        off = 1 if one_based else 0
        self.__init__(
            (molecule.atoms[idx1 - off], molecule.atoms[idx2 - off], molecule.atoms[idx3 - off], molecule.atoms[idx4 - off]),
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



class InternalCartesianX(InternalCartesian):
    """
    """

    ####################
    # Class Attributes #
    ####################

    analytic_b_orders = [1]

    #################
    # Class Methods #
    #################

    @classmethod
    def value_for_xyz(cls, xyz):
        def e(j, k):
            ev = LightVector.l_sub(xyz[k-1], xyz[j-1]); ev.normalize()
            return ev
        return BondLength.value_for_positions(xyz[0], xyz[3]) * (e(1,2).dot(e(1,4)))

    @classmethod
    def b_vector_for_positions(cls, *args):
        def R(i, j):
            # This is ugly...
            # TODO figure out a better way to do this
            mol = Molecule([Atom('H', args[i-1]), Atom('H', args[j-1])])
            return BondLength(mol[0], mol[1])
        Rab, Rad = R(1,2), R(1,4)
        BRab, BRad = Rab.get_b_tensor(max_order=2), Rad.get_b_tensor(max_order=2)
        Bb = Rad.value * Vector([BRab[3+alpha,3:].dot(BRad[3:]) for alpha in [X, Y, Z]])
        Bd = Rad.value * Vector([BRad[3+alpha,3:].dot(BRab[3:]) for alpha in [X, Y, Z]])
        Bd += BRab[3:].dot(BRad[3:]) * BRad[3:]
        Ba = -Bb - Bd
        return Ba, Bb, Vector(0,0,0), Bd

    ###################
    # Private Methods #
    ###################

    def _recompute_b_vector(self):
        #TODO Fix units
        b = Matrix(shape=(self.molecule.natoms, 3))
        i = self.atom_indices
        b[i[0]], b[i[1]], b[i[2]], b[i[3]] = self.b_vector_for_positions(*[a.pos for a in self.atoms])
        self._bvector = b.flatten('C').view(Vector)

class InternalCartesianY(InternalCartesian):
    """
    """

    ####################
    # Class Attributes #
    ####################

    analytic_b_orders = []

    #################
    # Class Methods #
    #################

    @classmethod
    def value_for_xyz(cls, xyz):
        def e(j, k):
            ev = LightVector.l_sub(xyz[k-1], xyz[j-1]); ev.normalize()
            return ev
        eab, eac, ead = e(1,2), e(1,3), e(1,4)
        tmp = eab.l_cross(eac).normalized()
        tmp2 = tmp.l_cross(eab).normalized()
        return BondLength.value_for_positions(xyz[0], xyz[3]) * tmp2.dot(ead)

    @classmethod
    def b_vector_for_positions(cls, *args):
        def R(i, j):
            # This is ugly...
            # TODO figure out a better way to do this
            mol = Molecule([Atom('H', args[i-1]), Atom('H', args[j-1])])
            return BondLength(mol[0], mol[1])
        def phi(i, j, k):
            mol = Molecule([Atom('H', args[i-1]), Atom('H', args[j-1]), Atom('H', args[k-1])])
            return BondAngle(mol[0], mol[1], mol[2])
        Rab, Rad, phibac = R(1,2), R(1,4), phi(2,1,3)
        rab, rad = Rab.value, Rad.value
        BRab, BRad, Bphi = map(lambda x: x.get_b_tensor(max_order=2), (Rab, Rad, phibac))
        Bb = -Bphi[:3,:3].view(Matrix) * BRad[3:].view(Vector) * rab * rad
        Bb -= Bphi[:3].dot(BRad[3:]) * rad * BRab[3:]
        Bc = -Bphi[:3,6:].view(Matrix) * BRad[3:].view(Vector) * rab * rad
        Bd = -BRad[3:,3:].view(Matrix) * Bphi[:3].view(Vector) * rab * rad
        Bd -= Bphi[:3].dot(BRad[3:]) * rab * BRad[3:]
        Ba = -Bb - Bc - Bd
        return Ba, Bb, Bc, Bd

    ###################
    # Private Methods #
    ###################

    def _recompute_b_vector(self):
        #TODO Fix units
        b = Matrix(shape=(self.molecule.natoms, 3))
        i = self.atom_indices
        b[i[0]], b[i[1]], b[i[2]], b[i[3]] = self.b_vector_for_positions(*[a.pos for a in self.atoms])
        self._bvector = b.flatten('C').view(Vector)

class InternalCartesianZ(InternalCartesian):
    """
    """

    ####################
    # Class Attributes #
    ####################

    analytic_b_orders = [1]

    #################
    # Class Methods #
    #################

    @classmethod
    def value_for_xyz(cls, xyz):
        def e(j, k):
            ev = LightVector.l_sub(xyz[k-1], xyz[j-1]); ev.normalize()
            return ev
        eab, eac, ead = e(1,2), e(1,3), e(1,4)
        tmp = eab.l_cross(eac).normalized()
        return BondLength.value_for_positions(xyz[0], xyz[3]) * tmp.dot(ead)

    @classmethod
    def b_vector_for_positions(cls, *args):
        return NotImplemented

    ###################
    # Private Methods #
    ###################

    def _recompute_b_vector(self):
        self._bvector = self.finite_difference_b_tensor(1).view(Vector)


from grendel.chemistry.molecule import Molecule
