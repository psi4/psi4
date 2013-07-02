"""
"""
from abc import ABCMeta, abstractmethod
from functools import partial
from numbers import Number
from grendel import type_checking_enabled, sanity_checking_enabled
from grendel.differentiation.finite_difference import Differentiable
from grendel.gmath.tensor import Tensor
from grendel.util.overloading import PartiallyConstructed
from grendel.util.parsing import MatrixRegexSequence
from grendel.util.strings import classname
from grendel.util.units import isunit, EnergyUnit, ValueWithUnits, IncompatibleUnitsError, DistanceUnit
from grendel.util.decorators import with_attributes, typechecked
from grendel.util.units.value_with_units import Unitized, hasunits

class ValueNotAvailableError(Exception):
    """ Raised when a MolecularProperty is requested but not available for whatever reason.
    """

class MolecularProperty(Differentiable):
    """ Abstract base class for all properties

    :Attributes:

    molecule : Molecule
    units : CompositeUnit or class with Unit metaclass


    """
    __metaclass__ = ABCMeta


    ####################
    # Class Attributes #
    ####################

    default_units = None

    ##############
    # Attributes #
    ##############

    molecule = None
    """ The molecule that `self` describes
    """

    units = None
    """ The units that the value is expressed in, or `None`
    if the property is unitless
    """

    details = None
    """ The ComputationDetails associated with the property.
    """

    getter = None
    """ A callable that takes the contents of a file as an argument and returns the value of the property (possibly
    as a `ValueWithUnits` instance).
    This replaces the deprecated from_groups/from_sequence motif.
    """



    ######################
    # Private Attributes #
    ######################

    _value = None


    ##################
    # Initialization #
    ##################

    def __init__(self, molecule, units=None, details=None):
        self._value = None
        self.details = details
        self.molecule = molecule
        self.getter = None
        if type_checking_enabled:
            if not isinstance(self.molecule, Molecule):
                raise TypeError("MolecularProperty initialization parameter"
                                " 'molecule' must be of type "
                                "Molecule (got {0})".format(
                                    type(molecule).__name__
                ))
            if units is not None and not isunit(units):
                raise TypeError("MolecularProperty initialization parameter"
                                " 'units' must be a unit or None")
        self.units = units or self.default_units


    ##############
    # Properties #
    ##############

    @property
    def value(self):
        """ The value of the molecular property.
        """
        return self._value

    @value.setter
    @typechecked(new_val=(Number, Tensor))
    def value(self, new_val):
        if hasunits(new_val):
            self._value = new_val.in_units(self.units).value
        else:
            self._value = new_val

    @property
    def value_with_units(self):
        if self.units is not None:
            return self.value * self.units
        else:
            return self.value

    ###################
    # Special Methods #
    ###################

    def __copy__(self):
        ret_val = self.__class__(self.molecule, self.units, self.details)
        ret_val.value = self.value
        return ret_val

    ##################
    # Static Methods #
    ##################

    @staticmethod
    def is_same_property(prop1, prop2):
        return MolecularProperty.property_type_of(prop1) == MolecularProperty.property_type_of(prop2)

    ##################
    # Static Methods #
    ##################

    @staticmethod
    def property_type_of(cls_or_inst):
        if isinstance(cls_or_inst, type):
            ret_val = cls_or_inst
        else:
            ret_val = type(cls_or_inst)
        if '_PartiallyConstructed__' in ret_val.__name__:
            return ret_val.__partial_parent__
        else:
            return ret_val

    #################
    # Class Methods #
    #################

    @classmethod
    def in_units(cls, units):
        return PartiallyConstructed(cls, units=units)


    ####################
    # Abstract Methods #
    ####################

    def from_sequence(self, seq, units):
        """ Retrieve the value of the property from a RegexSequence object that has already been used on a file.
        By default, this just calls self.from_groups() with `seq.groups(flatten)` where `flatten` comes from
        self.from_groups.flatten.
        """
        return self.from_groups(seq.groups(self.from_groups.flatten), units)

    @abstractmethod
    def from_groups(self, groups, units):
        """ Retrieve the value of the property from a tuple of groups matched in a `RegexSequence`.
        See `RegexSequence.groups()`
        """
        return NotImplemented


    ###########
    # Methods #
    ###########

    def has_value(self):
        return self._value is not None

    def clear_value(self):
        self._value = None

    def get_value(self, *args, **kwargs):
        """ The new paradigm for obtaining values.
        Call the callable instance attribute `getter` with the arguments and keyword arguments
        passed in.
        """
        # TODO More helpful error messages
        if self.getter is None:
            raise ValueNotAvailableError
        else:
            self.value = self.getter(*args, **kwargs)
            return self.value




    @property
    def property_type(self):
        ret_val = type(self)
        if '_PartiallyConstructed__' in ret_val.__name__:
            ret_val = ret_val.__partial_parent__
        return ret_val



class ScalarProperty(MolecularProperty):
    """ A property whose value is a scalar
    """

    ###################
    # Special Methods #
    ###################

    def __str__(self):
        if self._value is not None:
            return str(self.value_with_units)
        else:
            return '<{0} in {1} with no value yet>'.format(
                self.property_type.__name__,
                str(self.units or '[unitless]')
            )


    ###########
    # Methods #
    ###########

    @with_attributes(flatten=True)
    def from_groups(self, groups, units):
        """ By default, just get the last group matched.
        """
        if len(groups) == 0:
            raise ValueNotAvailableError(
                "Cannot find value of {0} for molecule:\n"
                "{1}\nand details:\n{2} in output "
                "file (no matching groups were found, "
                "check regular expressions)".format(
                    type(self).__name__,
                    self.molecule,
                    self.details
            ))
        if sanity_checking_enabled and len(groups) != 1:
            raise ValueError(
                "expecting exactly one value for ScalarProperty"
                 " '{0}', but got {1}".format(
                    classname(self.__class__), len(groups)
                ))
        if type_checking_enabled and not isinstance(groups[-1], basestring):
            raise TypeError("Got non-scalar last group: ('{1}') for ScalarProperty '{0}'".format(
                                classname(type(self)), "', '".join(groups)))
        ret_val = float(groups[-1]) * units
        self.value = ret_val.in_units(self.units)
        return ret_val


class TensorProperty(MolecularProperty):
    """ A property whose value is an n-dimensional tensor
    """


class VectorProperty(TensorProperty):
    """ A property whose value is a vector.
    """


class MatrixProperty(TensorProperty):
    """ A property whose value is a matrix.
    """

    ###########
    # Methods #
    ###########

    def from_sequence(self, seq, units):
        """ Default implementation for matrix properties; only works if `seq` is a MatrixRegexSequence
        """
        if isinstance(seq, MatrixRegexSequence):
            if units is not None:
                val = (seq.get_matrix() * units).in_units(self.units)
                val.units = None
                self._value = val
                return self.value_with_units
            elif self.units is None:
                return seq.get_matrix()
            else:
                raise IncompatibleUnitsError(self.units, units)

    @with_attributes(flatten=False)
    def from_groups(self, groups, units):
        """ Placeholder implementation to satisfy abstract class requirements.  Everything
        is already taken care of in `MatrixProperty.from_sequence`
        """
        return NotImplemented


class Energy(ScalarProperty):
    """ The energy of a molecule.
    """

    ####################
    # Class Attributes #
    ####################

    default_units = EnergyUnit.default


class OrbitalEnergies(VectorProperty):

    ####################
    # Class Attributes #
    ####################

    default_units = EnergyUnit.default

    ###########
    # Methods #
    ###########

    #@with_attributes(flatten=True)
    #def


#####################
# Dependent Imports #
#####################

from grendel.chemistry.molecule import Molecule

