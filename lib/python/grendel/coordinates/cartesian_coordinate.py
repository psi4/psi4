"""
"""
from grendel import type_checking_enabled
from numbers import Real
from grendel.coordinates.coordinate import Coordinate
from grendel.util.decorators import typechecked
from grendel.util.freezing import SetOnceAttribute
from grendel.util.strings import classname
from grendel.util.metaprogramming import ReadOnlyAttribute
from grendel.util.units import DistanceUnit, isunit, Angstroms, ValueWithUnits, Unitized, strip_units
from grendel.util.units.value_with_units import hasunits

__all__ = [
    "CartesianCoordinate"
]

# TODO A cartesian coordinate should really be a single x, y, or z in a single atom (to make indexing work)
class CartesianCoordinate(Coordinate):
    """
    A cartesian coordinate compatible with the Coordinate class.
    CartesianCoordinate objects are immutable.  If you need a different cartesian coordinate, create a new
    CartesianCoordinate objects.

    Attributes
    ----------
    atom : `Atom`
        The atom to which the cartesian coordinate refers
    direction : `int`
        The direction (x, y, or z) that the coordinate describes. (0 for x, 1 for y, 2 for z, which are also class constants)

    """

    ####################
    # Class Attributes #
    ####################

    default_delta = 0.01 * Angstroms

    ##################
    # Initialization #
    ##################

    @typechecked(
        atom='Atom',
        direction=int,
        parent=('CartesianRepresentation', None),
        units=(None, isunit),
        value=(None, Real),
        freeze_value=(bool, None)
    )
    def __init__(self,
            atom,
            direction,
            parent=None,
            parent_internal_coordinate=None,
            index=None,
            value=None,
            freeze_value=False,
            units=DistanceUnit.default,
            **kwargs):
        """ Constructor

        Parameters
        ----------
        atom : `Atom`
            The atom represented by the coordinate
        index : `int`
            The index of the coordinate in the parent representation
        direction : `int`
            The direction (x, y, or z) that the coordinate describes. (0 for x, 1 for y, 2 for z, which are also class constants)
        parent : `CartesianRepresentation`
            The representation containing the coordinate `self`

        """
        self._atom = atom
        self._direction = direction
        if freeze_value:
            # TODO think through the implications of this when I'm a little more cogent
            if value is not None:
                self._value = strip_units(value, units)
            elif not self.is_orphaned():
                self._value = strip_units(atom.position[direction] if value is None else value, units)
            else:
                raise ValueError("don't know how to get value for orphaned CartesianCoordinate")
        elif value is not None:
            raise NotImplementedError("value given for CartesianCoordinate, but"
                                      " 'freeze_value' was not set to True; this sort"
                                      " of functionality is not yet implemented.")
        if parent_internal_coordinate is not None:
            self.parent_internal_coordinate = parent_internal_coordinate
            self._index = index
        super(CartesianCoordinate, self)._init(
            units=units,
            parent=parent,
            freeze_value=freeze_value,
            **kwargs
        )
        if not self.is_orphaned():
            self._index = self.molecule.index(atom) * 3 + self.direction

    ##############
    # Properties #
    ##############

    @property
    def direction(self):
        """
        The direction (x, y, or z) that the coordinate describes. (0 for X, 1 for Y, 2 for Z,
        which are also class constants)
        """
        return self._direction

    @property
    def atom(self):
        """The atom to which the cartesian coordinate refers. """
        return self._atom

    @property
    def index(self):
        """ The index in the parent representation.  This allows for the
        retrieval of the value for the coordinate on a different molecule,
        for instance.
        """
        return self._index

    #-----------------------------------#
    # Properties abstract in Coordinate #
    #-----------------------------------#

    @property
    def atoms(self):
        return [self.atom]

    ###################
    # Special Methods #
    ###################

    def __short_str__(self):
        if self.atom.zmat_label is None:
            return "'{xyz}' of {atom} (#{num})".format(
                xyz=['X', 'Y', 'Z'][self.direction],
                num = int(self.index/3),
                atom=self.atom.symbol
            )
        else:
            return "'{xyz}' of {atomlabel}".format(
                xyz=['X', 'Y', 'Z'][self.direction],
                atomlabel=self.atom.zmat_label
            )

    def __str__(self):
        return self.__short_str__() + " with value '{}'".format(self.value)

    __repr__ = __str__ # for now...

    ###########
    # Methods #
    ###########

    def generate_name(self, one_based=True):
        off = 1 if one_based else 0
        return self.atom.symbol + str(self.atom.index + off) + ['X', 'Y', 'Z'][self.direction]

    def iter_molecule_indices(self):
        yield self.index

    #--------------------------------#
    # Methods abstract in Coordinate #
    #--------------------------------#

    def value_for_molecule_matrix(self, mat):
        return mat[self.index/3, self.index %3]

    def value_for_positions(self, *pos):
        return pos[0][self.direction]

    def copy_for_representation(self, rep, **kwargs):
        copykw = self.__copy_kwargs__()
        copykw.update(
            atom=rep.molecule[self.index//3],
            direction=self.direction,
            freeze_value=self._frozen_value,
            value=self.value if self._frozen_value else None
        )
        copykw.update(kwargs)
        units=copykw.pop('units', None) or self.units
        return self.__class__(
            parent=rep,
            units=units,
            **copykw
        )



#####################
# Dependent Imports #
#####################

from grendel.representations.cartesian_representation import CartesianRepresentation
from grendel.chemistry.atom import Atom

