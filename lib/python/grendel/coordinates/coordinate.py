"""
`Coordinate` is the abstract base class of all coordinates.  This most important thing to understand
about coordinates is that they are completely immutable, like strings and tuples in Python.  If you
need to change something about a coordinate, you must make a copy of it.
"""
from abc import ABCMeta, abstractproperty, abstractmethod
from grendel import type_checking_enabled, sanity_checking_enabled
from grendel.differentiation.finite_difference import FiniteDifferenceVariable
from grendel.util.decorators import typechecked, IterableOf, CachedProperty
from grendel.util.exceptions import SuperCallMissingError
from grendel.util.freezing import SetOnceAttribute, Freezable, FreezableAttribute
from grendel.util.iteration import first
from grendel.util.metaprogramming import ReadOnlyAttribute
from grendel.util.strings import short_str
from grendel.util.units import Radians, Unitized, isunit, AngularUnit, UnknownUnitError, CompositeUnit


__all__ = [
    "Coordinate"
]

class OrphanedCoordinateError(RuntimeError):
    """
    Raised when an invalid operation is performed on an orphaned `Coordinate`.
    """
    pass


class CoordinateDependencyError(Exception):
    """
    Raised when retrieval of an internal ("`Coordinate`'s own") index fails because
    the `Coordinate` does not depend on the molecular index or `CartesianCoordinate` requested.
    """

# TODO a cached state with the atom's positions so the coordinate can make sure the atoms haven't moved
class Coordinate(Unitized):
    """ Abstract base class for a general coordinate of a `Molecule`.

    `Coordinate`s are the components in which `Molecule`s are represented.  The most important
    thing about `Coordinate`s is that they are immutible, like Python's `str` and `tuple`
    classes.[#f1]_  Let me say that again a bit louder:

    .. note::
        All coordinates are immutable!

    Parentage and the `molecule` attribute
    --------------------------------------

    `Coordinate` instances are typically (though not always) associated with a `Representation`,
    accessible through the `parent` read-only attribute.  `Coordinate` instances that do not
    have a `parent` associated with them are called "orphaned" coordinates, and can be identified
    using the `is_orphaned()` instance method.  Even orphaned coordinates must have a `Molecule`
    associated with them in one way or another.  This happens in one of three ways:
    1. If the `Coordinate` is not orphaned, then `parent`'s `molecule` attribute is used.  Note
       that in this case, the `Coordinate` protocol requires the `Coordinate`'s constituant
       atoms to have a parent molecule that is exactly the same instance as `parent`'s `molecule`
       attribute.  (If you aren't doing something weird, this shouldn't be an issue).
    2. If the `Coordinate` is orphaned, it still must be composed of one or more `Atom` instances,
       accessible through the `atoms` attribute.  These `Atom`s may themselves be orphaned, but
       if they are not, the first non-orphaned `Atom`'s `parent` is used.  Note that if one `Atom` is
       non-orphaned, the `Coordinate` protocol requires that all of them must not be orphaned (i. e.
       the `Coordinate` constructor enforces this when `sanity_checking_enabled` is `True`)
    3. If the `Coordinate` is orphaned and all of its `Atom`s are also orphaned (as is the
       case, for instance, with `Coordinate` instances created for the purpose of finite difference
       B tensor computation), the `Coordinate`'s `base_analog` attribute may be set to a
       `Coordinate` instance whose constituate `Atom` instances correspond to the derived
       `Coordinate`'s atoms' `base_atom` attributes.  In this case, the `base_analog`'s
       `molecule` attribute is used, though it should be utilized for index translation purposes only,
       since the positions of the atoms in the base molecule will not correspond to the positions
       of the `Coordinate`'s atoms.  (The `base_analog`, in turn, may be an orphaned coordinate
       of this third type, in which case *that* coordinate's `base_analog` will be used, and so
       on recursively.)
    4. If the `Coordinate` is orphaned and no `base_analog` is given, then all of the atoms must have
       a `base_atom` that either has a `parent_molecule` or also has a `base_atom` defined (and so
       on so that a `base_atom` eventually has a `parent_molecule` somewhere up the line).  In this
       case, the `Coordinate` protocol requires that all the atoms refer to the same `Molecule` obtained
       in this manner.  Coordinates created in the course of analytic computation of B tensors fall
       into this catagory.

    Indexing schemes
    ----------------

    Because all `Coordinate`s may be austensibly associated with a `Molecule` instance in one way
    or another, there are two different ways to index various tensor properties of a `Coordinate`:
    the `Molecule`'s indices (`3 * natoms` total indices) and the `Coordinate`'s indexing scheme
    (`3 * len(coord.atoms)` total indices, where `coord` is a `Coordinate` instance).  The internal
    indexing scheme makes a lot more sense in many instances, (see the documentation of
    `Coordinate.atoms` for a simple example), but for determining properties of the whole
    molecule, we need to use the `Molecule`'s indexing scheme.  All of this should be handled
    seemlessly behind the scenes, but it's a good distinction to be aware of.

    Getting a `Coordinate`'s value
    ------------------------------

    `Coordinate` has a number of methods that can be used to get the value of a `Coordinate`
    in a particular scenario.  They hierarchically call each other to determine the value of
    the `Coordinate`.  The basic scheme is as follows:

    .. graphviz::
       digraph value_getting {
           node [shape=box];
           value_with_units -> value -> value_for_molecule -> value_for_molecule_matrix
           value_with_units [label="value_with_units()"]
           value [label="@property\nvalue"]
           value_for_molecule [label="value_for_molecule()"]
           value_for_molecule_matrix [label="value_for_molecule_matrix()",style=dotted];
           value_for_molecule_matrix-> value_for_positions [style=dotted,label="(only InternalCoordinate)"]
           value_for_positions -> value_for_xyz [style=dotted,label="(only SimpleInternalCoordinate)"];
           value_for_positions [label="@classmethod\nvalue_for_positions()",style=dotted];
            value_for_xyz [label="@classmethod\nvalue_for_xyz()",style=dotted];
       }

    The reason for this relatively complex hierarchy is that there are instances in which I have needed
    each variety for some reason or another.  Nodes with dotted borders must be implemented in
    the relevant subclasses:  `SimpleInternalCoordinate` subclasses must implement `value_for_xyz`,
    `InternalCoordinate` subclasses must implement `value_for_positions`, and `Coordinate` subclasses
    must implement `value_for_molecule_matrix()`.

    Subclassing
    -----------

    TODO: write this part

    .. rubric:: Footnotes

    .. [#f1] Since Python doesn't really have private attributes, you *could in theory* change
             some of the attributes of a Coordinate instance after it is created; the designation
             of "immutible" here merely means that the class is not guarenteed to work if you do.
             While it's technically more of a "please don't" than a "you can't", I for one will not
             be held responsible for your code not working if you do something like that, and
             the way these private attributes are used and changed is not guarenteed to stay the
             same between different versions.  Bottom line: just don't do it.  If you don't know
             what I'm talking about, you're probably okay.  Just don't access or modify any
             attributes that start with a single underscore (you should never do this anyway
             from outside of a class).

    """
    __metaclass__ = ABCMeta

    ####################
    # Class Attributes #
    ####################

    coordinate_symbol = "chi"

    ##############
    # Attributes #
    ##############

    name = None

    ######################
    # Private Attributes #
    ######################

    _index = None
    _mol_to_internal = None
    _value = None
    _frozen_value = False

    ##################
    # Initialization #
    ##################

    #-----------------------------------------------------------------#
    # Initialization methods to be called internally or by subclasses #
    #-----------------------------------------------------------------#

    @typechecked(
        parent=('Representation', None),
        units=(isunit, None),
        base_analog=('Coordinate', None),
        freeze_value=(bool, None),
        name=(basestring, None))
    def _init(self,
            parent=None,
            units=None,
            base_analog=None,
            freeze_value=False,
            name=None):
        #--------------------------------------------------------------------------------#
        # Miscellanea
        self._parent = parent
        self._init_units(units)
        self._base_analog = base_analog
        #--------------------------------------------------------------------------------#
        # Naming
        self.name = name
        if name is None:
            self.name = self.generate_name()
        #--------------------------------------------------------------------------------#
        # Frozen values
        if freeze_value and type(self).__name__ != 'CartesianCoordinate':
            raise NotImplementedError("freeze_value not yet implemented for '{}'.".format(
                self.__class__.__name__
            ))
        self._frozen_value = freeze_value or False
        #--------------------------------------------------------------------------------#
        if sanity_checking_enabled:
            # Now that construction is complete, make sure we have all of the things we need.
            # And make sure that subclasses are behaving themselves and calling super:
            self._Coordinate_validate_called = False
            self._validate_initialization()
            if not self._Coordinate_validate_called:
                #noinspection PyExceptionInherit
                raise SuperCallMissingError(
                    "'{}' or one of its super classes failed to call"
                    " super()._validate_initialization()".format(
                    self.__class__.__name__))
        #--------------------------------------------------------------------------------#
        # Now call finalize.  Same drill as with _validate_coordinate(), requiring subclasses
        #   to behave themselves.
        self._Coordinate_finalize_called = False
        self._finalize_initialization()
        if not self._Coordinate_finalize_called:
            #noinspection PyExceptionInherit
            raise SuperCallMissingError(
                "'{}' or one of its super classes failed to call"
                " super()._finalize_initialization()".format(
                self.__class__.__name__))

    @typechecked(
        units=(isunit, dict, None)
    )
    def _init_units(self, units):
        if isinstance(units, dict):
            if self.__class__.default_delta.units.genre is CompositeUnit:
                raise NotImplementedError
            else:
                units = units[self.__class__.default_delta.units.genre]
        if type_checking_enabled:
            if not isunit(units):
                raise UnknownUnitError(units)
        if sanity_checking_enabled:
            # try to convert the default_delta's units to these
            # units to see if they are compatible:
            if self.__class__.default_delta is not None:
                self.__class__.default_delta.units.to(units)
        self._units = units

    def _validate_initialization(self):
        """ Make sure we have all of the things we need.  Subclasses may override this,
        but all subclass implementations should call `super()._validate_initialization()` at
        the end of this method (if they do not, a SuperCallMissingError is raised).

        .. note::
        This method should not make any modifications to the `Coordinate` itself, since it is
        only called when `sanity_checking_enabled` is `True`.  Use `_finilize_initialization()`
        for this purpose.
        """
        #--------------------------------------------------------------------------------#
        # Type 1 coordinate: non-orphaned
        if not self.is_orphaned():
            # Enforce the requirement that all atoms of a type 1 coordinate be non-orphaned
            if any(a.is_orphaned() for a in self.atoms):
                raise ValueError(
                    "non-orphaned Coordinates must be constructed entirely of non-orphaned"
                    " atoms, but atom '{}' on coordinate '{}' is orphaned.".format(
                        first(a for a in self.atoms if a.is_orphaned()),
                        self
                    ))
            # Enforce the requirement that all atoms of a type 1 coordinate should have
            #   the exact same parent molecule as the parent representation's molecule
            if any(a.parent is not self.molecule for a in self.atoms):
                raise ValueError(
                    "all atoms in a non-orphaned coordinate must have the same parent"
                    " molecule (i.e. the exact same instance) as the coordinate, but"
                    " atom '{}' does not have the same parent molecule as Coordinate"
                    " '{}' ({} != {}).".format(
                        first(a for a in self.atoms if a.parent is not self.molecule),
                        self,
                        id(first(a for a in self.atoms if a.parent is not self.molecule).parent),
                        id(self.molecule)
                    ))
        #--------------------------------------------------------------------------------#
        # Type 2 coordinate: orphaned with non-orphaned atoms
        elif any(not a.is_orphaned() for a in self.atoms):
            # Enforce the requirement that if one atom is non-orphaned, all of them must be.
            if any(a.is_orphaned() for a in self.atoms):
                raise ValueError(
                    "orphaned coordinates composed of non-orphaned atoms must have all"
                    " atoms non-orphaned, but atom '{}' on coordinate '{}' is orphaned.".format(
                        first(a for a in self.atoms if a.is_orphaned()),
                        self
                    ))
            # Enforce the requirement that all atoms be on the same molecule
            if any(a.parent is not self.molecule for a in self.atoms):
                raise ValueError(
                    "all atoms in a orphaned coordinate composed of non-orphaned atoms"
                    "must have the same parent molecule (i.e. the exact same instance),"
                    " but atom '{}' does not have the same parent molecule as '{}' for"
                    " coordinate '{}'.".format(
                        first(a for a in self.atoms if a.parent is not self.molecule),
                        self.atoms[0],
                        self
                    ))
        #--------------------------------------------------------------------------------#
        # Type 3 coordinate: orphaned coordinate with all orphaned atoms, base_analog not None
        elif self.base_analog is not None:
            # Enforce the requirement that the base_analog's atoms exactly match the
            #   coordinate's atoms' respective base_atom attributes
            if len(self.base_analog.atoms) != len(self.atoms):
                raise ValueError("orphaned coordinates with a base_analog must be composed"
                                 " of the same number of Atoms as their base_analog (got {}"
                                 " atoms for derived coordinate, {} atoms for base"
                                 " coordinate)".format(
                    len(self.atoms),
                    len(self.base_analog.atoms)))
            for my_atom, base_atom in zip(self.atoms, self.base_analog.atoms):
                if my_atom.base_atom is not base_atom:
                    raise ValueError("base_atom of displaced Atom must correspond to the"
                                     " analogous atom in the displaced coordinate's"
                                     " base_analog's atoms, but '{}' is not the same as"
                                     " '{}' for coordinate '{}'".format(
                        my_atom.base_atom, base_atom, short_str(self)
                    ))
        #--------------------------------------------------------------------------------#
        # Type 4 coordinate: orphaned coordinate with all orphaned atoms, base_analog is None
        else:
            # Enforce the requirement that all atom's eventual parent molecules should be
            #   the same.
            mol_found = None
            for atom in self.atoms:
                iter_atom = atom
                visited = []
                while iter_atom is not None and iter_atom.parent is None:
                    if iter_atom in visited:
                        raise ValueError("circular base_atom dependency for '{}'".format(atom))
                    iter_atom = iter_atom.base_atom
                if iter_atom is None:
                    raise ValueError("all orphaned atoms on an orphaned coordinate must have"
                                     " a `base_atom` that is not orphaned somewhere down the"
                                     " line.  '{}' in coordinate '{}' is not like that.".format(
                        atom, short_str(self)
                    ))
                if mol_found is None:
                    mol_found = iter_atom.parent
                elif mol_found is not iter_atom.parent:
                    raise ValueError("all orphaned atoms on an orphaned coordinate must have"
                                     " the same parent molecule of the first non-orphaned coordinate"
                                     " in the `base_atom` hierarchy.  The molecule obtained in this"
                                     " manner from '{}' in coordinate '{}' is not the same"
                                     " as the molecule obtained in this manner from '{}'.".format(
                        atom, self.atoms[0], short_str(self)
                    ))
        #--------------------------------------------------------------------------------#
        self._Coordinate_validate_called = True

    def _finalize_initialization(self):
        """ Do anything that needs to be done post-initialization to bring the coordinate to its
        final, immutible state.  This is called at the end of `Coordinate._init()` and allows
        subclasses to apply any special handling needed for the stuff set up by the hierarchy of
        `_init()` methods.  All subclass implementations should call
        `super()._finalize_initialization()` at the end of this method (if they do not,
        a SuperCallMissingError is raised).

        """
        #----------------------------------------#
        self._Coordinate_finalize_called = True

    ###################
    # Special Methods #
    ###################

    def __copy_kwargs__(self):
        return dict(
            units=self.units,
            name=self.name if self.name != self.generate_name() else None
        )

    #############################
    # Abstract Class Properties #
    #############################

    @abstractproperty
    def default_delta(self):
        """ Every coordinate must define a reasonable default amount for a finite displacement
        of that kind of coordinate.  This should be a class attribute.
        """
        return NotImplemented

    #######################
    # Abstract Properties #
    #######################

    @abstractproperty
    def atoms(self):
        """ List of the `Atom` instances that the `Coordinate` depends on.  In the case of
        a `CartesianCoordinate`, for instance, this is trivially one `Atom`, but in the case
        of other coordinates it can be much more substantial.  In the case of a `NormalCoordinate`,
        this is all the atoms in the parent `Molecule`.

        The presence of an atoms list is essential to establish an internal indexing scheme that
        can be converted to the indexing scheme of the parent molecule.  Besides making the
        B tensor code much more managable, the use of an internal indexing scheme make `Coordinate`
        instances much more self-contained, making "orphaned" coordinates much more viable.  Consider
        the following example:

            >>> from grendel.util.units import Degrees
            >>> from grendel.coordinates.bond_angle import BondAngle
            >>> ba = BondAngle(
            ...     Atom('O', [0.0, 0.0, 0.0]),
            ...     Atom('H', [0.0, 1.0, 0.0]),
            ...     Atom('H', [0.0, 0.0, 1.0]),
            ...     units=Degrees
            ... )
            >>> round(ba.value)
            90.

        If orphaned coordinates where not allowed to exist, we would have to create a `Molecule`
        instance containing those three atoms, construct an `InternalRepresentation` containing
        the `BondAngle` instance (which would need to be a valid, complete representation, a
        non-trivial task in the general case), pair that `InternalRepresentation` instance with the
        `Molecule` instance, get the right `Coordinate` instance from the `InternalRepresentation`
        corresponding to the `BondAngle` we want, and then get that coordinate's value.  If all of
        the coordinate's complicated methods were not implemented in terms of the Coordinate's own
        indexing scheme, it would have to have all of this structure in place just to get a value.

        .. note::
        Since all `Coordinate`s depend on `Atom`s and `Atom`s are not immutible, a `Coordinate`s
        `value` and other properties such as its `b_vector` *can* change even if the `Coordinate`
        itself cannot.  The one way around this is to specify `freeze_value` as `True` in the
        `Coordinate`'s constructor.  If that is done, the value of the `Coordinate` will not
        change over the life of the instance, even if the positions of the constituant atoms
        does change (so if you use this, make sure you know what you are doing).  This behavior
        is useful, for instance, in the use of `CartesianRepresentation`s as "snapshots" of
        a `Molecule` instance to be used for the parsing of a particular
        `RepresentationDependentProperty`.
        """
        return NotImplemented

    ################
    #  Properties  #
    ################

    @property
    def atom_indices(self):
        """ The indices of the atom(s) the `Coordinate` depends on in the `Molecule` the
        `Coordinate` depends on.

        This is guarenteed to work whether or not the `Coordinate` is orphaned.  For each
        of the various possible `Coordinate`--`Molecule` relationships detailed in the
        `Coordinate` class documentation, the reason the method works is as follows:
        1. If the `Coordinate` is not orphaned, then the `Coordinate`'s atoms must not be
           orphaned (the `Coordinate` constructor enforces this when `sanity_checking_enabled` is
           `True`).  Thus, each of the atoms must be non-orphaned and thus have an index.
        2. If the Coordinate is orphaned and at least one atom is non-orphaned, all atoms
           must be non-orphaned (enforced by the Coordinate protocol).  Thus, we can use
           the non-orphaned atoms' indices.
        3. If the coordinate is orphaned and all of its atoms are orphaned and `base_analog` is not
           None, the atoms must have the `base_atom` attribute set to the corresponding atom in the
           `Coordinate`'s `base_analog`.  Thus, the `Atom`'s `base_atom`'s index is used.
        4. If the coordinate is orphaned and all of its atoms are orphaned and `base_analog` is
           None, then each atom must have a `base_atom` somewhere down the line that has a parent
           and thus an index; this atom's index is used.

        Use this carefully.  Don't ever retrieve an atom's position by its index in the parent
        molecule, since that atom may not be the same instance (or even in the same position)
        as the `Coordinate`'s atom.
        """
        return [atom.index for atom in self.atoms]

    @property
    def molecule(self):
        """
        The molecule associated with the `Coordinate`.

        .. note::
        The molecule's atoms' positions may not correspond to the positions of the atoms the
        coordinate describes.  See the class documentation for `Coordinate` for details.

        """
        if not self.is_orphaned():
            # Type 1:  Non-orphaned coordinate
            return self.parent.molecule
        else:
            # The rest of this is a bit more expensive, so we cache it
            if hasattr(self, '_molecule'):
                return self._molecule
            for atom in self.atoms:
                if not atom.is_orphaned():
                    # Type 2: not all of the atoms are orphaned
                    self._molecule = atom.parent
                    return self._molecule
            if self.base_analog is not None:
                # Type 3: all Atoms are orphaned, base_analog is not None
                self._molecule = self.base_analog.molecule
                return self._molecule
            else:
                # Type 4: all Atoms orphaned, base_analog is None
                if self.base_analog is None and all(a.is_orphaned() for a in self.atoms):
                    atom = self.atoms[0]
                    while atom is not None and atom.parent is None:
                        atom = atom.base_atom
                    self._molecule = atom.parent
                    return self._molecule


    @CachedProperty
    def molecule_indices(self):
        """
        A list of the indices in the molecule's indexing scheme that the coordinate depends on.
        There should be `3 * len(coord.atoms)` of these (with the exception of `CartesianCoordinate`
        instances, which should have just one).
        """
        return list(self.iter_molecule_indices())

    @property
    def value(self):
        if self._value is not None:
            return self._value
        elif not self.is_orphaned():
            return self.value_for_molecule(self.molecule)
        else:
            return self.value_for_positions(*[a.pos for a in self.atoms])

    @property
    def value_with_units(self):
        return self.value * self.units

    @property
    def parent(self):
        """
        The :py:class:`~grendel.representations.internal_representation.InternalRepresentation`
        instance which the coordinate `self` is a part of.
        """
        return self._parent
    parent_representation = parent

    @property
    def base_analog(self):
        """ Analogous `Coordinate` object on the base molecule, if self is a `Coordinate`
        on a displaced `Molecule`, or the `Coordinate` from which an orphaned displaced `Coordinate`
        was created, if the coordinate was created this way.
        """
        return self._base_analog

    @property
    def units(self):
        """
        The units in which the coordinate's value is expressed.
        """
        return self._units

    @property
    def index(self):
        if self._index is not None:
            return self._index
        elif self.is_orphaned():
            raise OrphanedCoordinateError("'index' property is only valid for non-orphaned"
                                          " coordinates.  Coordinate '{}' is orphaned.".format(self))
        else:
            self._index = self.parent.coords.index(self)
            return self._index


    ####################
    # Abstract Methods #
    ####################

    @abstractmethod
    def value_for_molecule_matrix(self, mat):
        """ Get the value of the coordinate for the XYZ matrix `mat`,
        which is a (`self.molecule.natoms`) by 3 matrix.  In other words, the matrix argument uses
        the parent `Molecule`'s indexing scheme, not the `Coordinate`'s own indexing scheme.

        .. note::
        For angular coordinates, this *always* returns a value in Radians.  `value`,
        `value_with_units`, and `value_for_molecule`, however, return values in `self.units`

        """
        return NotImplemented

    @abstractmethod
    def value_for_positions(self, *pos):
        return NotImplemented

    @abstractmethod
    def copy_for_representation(self, rep, **kwargs):
        return NotImplemented

    ###########
    # Methods #
    ###########

    def value_for_molecule(self, mol):
        if self.units.genre is AngularUnit:
            return self.value_for_molecule_matrix(mol.xyz) * Radians.to(self.units)
        else:
            return self.value_for_molecule_matrix(mol.xyz)

    @typechecked(internal_indices=IterableOf(int))
    def molecule_indices_for(self, internal_indices):
        mol_idxs = self.molecule_indices
        return tuple(mol_idxs[internal] for internal in internal_indices)

    @typechecked(molecule_indices=IterableOf(int))
    def internal_indices_for(self, molecule_indices):
        if self._mol_to_internal is None:
            self._mol_to_internal = {}
        ret_val = []
        all_mol_indices = self.molecule_indices
        for mol_idx in molecule_indices:
            int_idx = self._mol_to_internal.get(mol_idx, None)
            if int_idx is None:
                try:
                    int_idx = all_mol_indices.index(mol_idx)
                except ValueError:
                    raise CoordinateDependencyError(
                        "Coordinate '{}' does not depend on molecular index {}".format(
                            self, mol_idx
                        ))
                self._mol_to_internal[mol_idx] = int_idx
            ret_val.append(int_idx)
        return tuple(ret_val)

    def generate_name(self, one_based=True):
        off = 1 if one_based else 0
        if not self.is_orphaned():
            return self.coordinate_symbol + str(self.index+off)
        else:
            return self.coordinate_symbol + '(o)'

    def orphaned_copy(self, **kwargs):
        copykw = self.__copy_kwargs__()
        copykw.update(kwargs)
        copykw.update(parent=None)
        return self.__class__(**copykw)

    #------------------------------#
    # Methods abstract in Unitized #
    #------------------------------#

    def in_units(self, new_units):
        return (self.value * self.units).in_units(new_units)

    #-----------#
    # Iteration #
    #-----------#

    def iter_molecule_indices(self):
        """
        Default implementation for coordinates dependent on atoms directly (basically anything
        but a `CartesianCoordinate` or some sort of symmetrized cartesian coordinate).
        """
        for atom in self.atoms:
            if not atom.is_orphaned():
                for direction in [X, Y, Z]:
                    yield 3*atom.index + direction
            else:
                visited = None
                if sanity_checking_enabled:
                    visited = []
                while atom.is_orphaned() and atom.base_atom is not None:
                    if sanity_checking_enabled:
                        if atom in visited:
                            raise ValueError("circular base_atom dependency")
                        visited.append(atom)
                    atom = atom.base_atom
                if sanity_checking_enabled:
                    if atom.is_orphaned():
                        raise ValueError("orphaned Atom instance {} does not have a"
                                         " base_atom".format(atom))
                for direction in [X, Y, Z]:
                    yield 3*atom.index + direction

    #-----------------#
    # Inquiry methods #
    #-----------------#

    def is_orphaned(self):
        """
        Whether or not the `Coordinate` instance has a parent `Representation`
        """
        return self.parent is None


#####################
# Dependent Imports #
#####################

from grendel.representations.cartesian_representation import X, Y, Z

# needed for dynamic typechecking, do not delete
if type_checking_enabled:
    from grendel.representations.representation import Representation
