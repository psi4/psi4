from abc import ABCMeta
from collections import Iterable
from copy import deepcopy, copy
from grendel.chemistry.molecular_properties import ScalarProperty, MolecularProperty
from grendel.differentiation.derivative_tensor import RepresentationDependentTensor
from grendel.interface.computation import ComputationUnavailableError
from grendel.interface.computation_details import ComputationDetails, DetailSpecificObject
from grendel.util.decorators import CachedProperty
from grendel.util.strings import classname
from grendel.util.units.unit import isunit

class OutputParsingError(RuntimeError): pass

class OutputParser(object):
    """ Abstract base class for classes that parse the output file for a computation.
    """
    __metaclass__ = ABCMeta


    #############################
    # Class/Instance Attributes #
    #############################

    property_parsers = None


    ##############
    # Attributes #
    ##############

    computation = None



    ##################
    # Initialization #
    ##################

    def __init__(self, computation, property_parsers=None):
        if property_parsers is not None:
            self.property_parsers = deepcopy(property_parsers)
        else:
            self.property_parsers = deepcopy(self.__class__.property_parsers)
        self.computation = computation

    ##############
    # Properties #
    ##############

    @property
    def needed_properties(self):
        if self.computation is None:
            return []
        else:
            return [prop for prop in self.computation.properties if not prop.has_value()]

    @property
    def needed_indices(self):
        if self.computation is None:
            return []
        else:
            return [i for i, prop in enumerate(self.computation.properties) if not prop.has_value()]

    @CachedProperty
    def valid_parsers(self):
        valid_parsers = {}
        for property in self.needed_properties:
            for parser in self.property_parsers:
                if parser.compatible_with_details(self.computation.details):
                    if MolecularProperty.is_same_property(property, parser.property):
                        valid_parsers[property.property_type] = parser
                        break
            if not any(MolecularProperty.is_same_property(property, p) for p in valid_parsers):
                raise ComputationUnavailableError("Could not find a valid parser for {0} in class {1} compatible with details: "
                                 "\n{2}".format(property.property_type.__name__, classname(self.__class__), self.computation.details))
        return valid_parsers

    @property
    def molecule(self):
        return self.computation.molecule

    ###########
    # Methods #
    ###########

    def validate(self):
        # Get the parsers.  A ComputationUnavailableError will be raised if things won't work out
        self.valid_parsers
        return True

    def reset(self):
        for p in self.property_parsers:
            p.reset()

    def parse_file(self, file):
        file_contents = None
        remaining_properties = []
        for property in self.needed_properties:
            p = self.valid_parsers[property.property_type]
            if p.getter is not None:
                # New style parser
                if file_contents is None:
                    # TODO @optimization this could slow things down a lot...
                    file_contents = open(file).read()
                p.get_value_for_molecule(property, file_contents)
            else:
                remaining_properties.append(property)
        # Old style, if any properties are left
        # DEPRECATED
        if len(remaining_properties) > 0:
            with open(file) as f:
                for line in f:
                    for prop in remaining_properties:
                        self.valid_parsers[prop.property_type].parse_line(line)
            # Now fill in the results
            for property in remaining_properties:
                p = self.valid_parsers[property.property_type]
                if p.custom_from is None:
                    property.from_sequence(p.sequence, p.units)
                else:
                    flatten = True
                    if hasattr(p.custom_from, 'flatten'):
                        flatten = p.custom_from.flatten
                    property.value = p.custom_from(p.sequence.groups(flatten), p.units)


class PropertyParser(DetailSpecificObject):
    """ Encapsulates a molecular property and the code to parse it from an output file.
    """

    ##############
    # Attributes #
    ##############

    getter = None
    """ New way of getting properties.  This is just a callable
    """

    sequence = None
    """ DEPRECATED
    The regular expression sequence to use for parsing.  Should have a `parse_line` method (typically,
    this is going to be an instance of RegexSequence"""

    property = None
    """ The class of the property to be parsed. """

    units = None
    """ The units the parsed property should be interpreted in."""

    custom_from = None
    """ DEPRECATED
    Advanced use only.
    A custom `MolecularProperty.from_sequence()` funtion to use instead of the one found in `self.property.from_sequence()`
    """


    ##################
    # Initialization #
    ##################

    def __init__(
            self,
            property,
            details,
            units,
            getter=None,
            sequence=None,
            unsupported_details=None,
            custom_from=None):
        self.getter = getter
        self.property = property
        self.sequence = sequence
        self.details = details
        self.units = units
        self.unsupported_details = unsupported_details or []
        self.custom_from = custom_from


    ###################
    # Special Methods #
    ###################

    def __deepcopy__(self, memo):
        return self.__class__(
            property=self.property,
            sequence=deepcopy(self.sequence, memo),
            details=deepcopy(self.details, memo),
            units=self.units,
            getter=self.getter,
            unsupported_details=deepcopy(self.unsupported_details, memo),
            custom_from=self.custom_from
        )

    def __copy__(self):
        # This is probably not what you want to do...
        return self.__class__(
            property=self.property,
            sequence=copy(self.sequence),
            details=copy(self.details),
            units=self.units,
            getter=self.getter,
            unsupported_details=copy(self.unsupported_details),
            custom_from=self.custom_from
        )

    ###########
    # Methods #
    ###########

    def get_value_for_molecule(self, property, file_contents):
        property.value = self.getter(file_contents)
        return property.value

    # deprecated
    def parse_line(self, line):
        self.sequence.parse_line(line)

    def reset(self):
        self.sequence.reset()


class RepresentationDependentPropertyParser(PropertyParser):

    ##############
    # Attributes #
    ##############

    representation_getter = None

    ##################
    # Initialization #
    ##################

    def __init__(self, representation_getter, **kwargs):
        self.representation_getter = representation_getter
        super(RepresentationDependentPropertyParser, self).__init__(**kwargs)

    def get_value_for_molecule(self, property, file_contents):
        representation = self.representation_getter(file_contents, property.molecule)
        gotten_val = self.getter(file_contents)
        tens = RepresentationDependentTensor(
            gotten_val,
            representation=representation,
            units=gotten_val.units
        )
        if isunit(property.units):
            tens = tens.in_units(property.units)
        # Now transform to the relevant representation of the molecule
        new_tens=tens.in_representation(property.representation)
        property.value = new_tens
        return property.value

    def get_value_and_representation(self, property, file_contents):
        """
        ..note ::
            This method ignores `property.units`
        """
        representation = self.representation_getter(file_contents, property.molecule)
        gotten_val = self.getter(file_contents)
        return gotten_val, representation



    # TODO Get rid of the need to call deepcopy by completely deprecating sequence-based PropertyParsers
    def __deepcopy__(self, memo):
        return self.__class__(
            property=self.property,
            sequence=deepcopy(self.sequence, memo),
            details=deepcopy(self.details, memo),
            units=self.units,
            getter=self.getter,
            unsupported_details=deepcopy(self.unsupported_details, memo),
            representation_getter=self.representation_getter,
            custom_from=self.custom_from
        )

    def __copy__(self):
        # This is probably not what you want to do...
        return self.__class__(
            property=self.property,
            sequence=copy(self.sequence),
            details=copy(self.details),
            units=self.units,
            getter=self.getter,
            unsupported_details=copy(self.unsupported_details),
            representation_getter=self.representation_getter,
            custom_from=self.custom_from
        )

class OptimizedGeometryParser(PropertyParser):

    def get_value_for_molecule(self, property, file_contents):
        raise NotImplementedError

    def __init__(self, geometry_getter, **kwargs):
        self.geometry_getter = geometry_getter
        super(OptimizedGeometryParser, self).__init__(**kwargs)

    def get_optimized_geometry(self, file_contents):
        return self.geometry_getter(file_contents)

