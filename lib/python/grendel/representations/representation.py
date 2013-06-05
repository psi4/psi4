from abc import ABCMeta, abstractmethod
from grendel.gmath.vector import Vector
from grendel.util.aliasing import function_alias
from grendel.util.freezing import Freezable, FreezableAttribute, FreezableListAttribute
from grendel.util.metaprogramming import ReadOnlyAttribute

__all__ = [
    "Representation"
]

class RepresentationError(ValueError):
    """ An error in the represetation of a molecule.
    """

class Representation(Freezable):
    """
    Superclass of all the representations types.

    Attributes
    ----------
        molecule : `Molecule`
            The `Molecule` object represented by `self`.
        coords : list of `Coordinate`
            The coordinates that make up the representation

    """

    __metaclass__ = ABCMeta

    ##############
    # Attributes #
    ##############

    molecule = FreezableAttribute(
        name="molecule",
        doc="""
        The `Molecule` object represented by `self`.
        """
    )

    coords = FreezableListAttribute(
        name="coords",
        doc="""
        The coordinates that make up the representation
        """
    )

    units = ReadOnlyAttribute(
        name='units',
        doc="""
        The units to use for the coordinates created, or if the created coordinates vary in units, a `dict` of
        `UnitGenre`, `Unit` pairs.  Must be passed into constructor; defaults to `genre.default`, where `genre` is
        the `UnitGenre` subclass of the applicable units.
        """
    )

    ###################
    # Special Methods #
    ###################

    def __getitem__(self, item):
        return self.coords[item]

    def __len__(self):
        return len(self.coords)

    def __iter__(self):
        for coord in self.coords:
            yield coord

    ##############
    # Properties #
    ##############

    @property
    def values(self):
        vals = [c.value for c in self.coords]
        return Vector(vals)
    value = values

    ####################
    # Abstract Methods #
    ####################

    @abstractmethod
    def add_coordinate_copy(self, coordinate):
        return NotImplemented

    @abstractmethod
    def copy_with_molecule(self, molecule):
        """ Make a copy of `self` that is the same in every way except for the `molecule` attribute.
        New Coordinate objects are created using the `Coordinate.copy_for_representation()` method for
        each element of `self.coords`.
        This is an abstract method that *must* be implemented by all Representation subclasses.
        """
        return NotImplemented

    @abstractmethod
    def displaced_by(self, disp, tol=None, maxiter=None):
        """ Apply the `Displacement` instance `disp` to the molecule and current representation,
        generating a new molecule and a new representation (which start as a `deepcopy` and the return value of
        `Representation.copy_with_molecule`, respectively) with the displacement applied.
        This is an abstract method that *must* be implemented by all Representation subclasses.
        """
        return NotImplemented

    ###########
    # Methods #
    ###########

    def values_for_molecule(self, mol):
        return Vector([c.value_for_molecule(mol) for c in self.coords])
    value_for_molecule = function_alias('value_for_molecule', values_for_molecule)

    def values_for_matrix(self, mat):
        return Vector([c.value_for_molecule_matrix(mat) for c in self.coords])
    value_for_matrix = function_alias('value_for_matrix', values_for_matrix)


