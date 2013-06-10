from abc import ABCMeta, abstractmethod
from collections import Iterable
from warnings import warn

import numpy
import grendel

from grendel.chemistry.molecular_properties import *
from grendel.interface.computation import ComputationUnavailableError
from grendel.interface.computation_details import DetailSpecificObject
from grendel.util.sentinal_values import Keyword, All
from grendel.util.strings import classname, camel_to_lower

class InputGenerator(object):
    """ Abstract base class for classes that generate input files.
    """
    __metaclass__ = ABCMeta

    ##############
    # Attributes #
    ##############

    computation = None
    """ The Computation object to generate an input for."""

    generated = None

    ##################
    # Initialization #
    ##################

    def __init__(self, computation):
        self.computation = computation
        self.generated = False


    ##############
    # Properties #
    ##############

    @property
    def needed_properties(self):
        if self.computation is None:
            return []
        if self.computation.started:
            return []
        else:
            return [prop for prop in self.computation.properties if isinstance(prop, type) or not prop.has_value()]

    ####################
    # Abstract Methods #
    ####################

    @abstractmethod
    def generate(self, filename):
        return NotImplemented

    @abstractmethod
    def validate(self):
        """ Determine if we will be able to generate an input.
        """
        return NotImplemented


class TemplateInputGenerator(InputGenerator):
    """ The simplest form of a generator.
    Takes a Mako-style template file, passes in the Computation object in as the argument "computation", and runs
    the Mako templete engine.
    """

    ####################
    # Class Attributes #
    ####################

    property_generators = []


    ###########
    # Methods #
    ###########

    def validate(self):
        prop_gen = None
        for gen in self.property_generators:
            if gen.property is All or all(
                    any(MolecularProperty.is_same_property(prop, p) for p in gen.property)
                        for prop in self.needed_properties):
                if gen.compatible_with_details(self.computation.details):
                    prop_gen = gen
                    break
        # and raise an error if we can't find a property generator that will work for us
        if prop_gen is None:
            raise ComputationUnavailableError("Could not find template file capable of generating an input to compute properties [{0}]"
                                              " with details given by {1}".format(
                ", ".join(classname(prop) for prop in self.needed_properties),
                self.computation.details))
        return True

    def generate(self, filename):
        # Find the first valid property generator that will work for what we need to do
        prop_gen = None
        for gen in self.property_generators:
            if gen.property is All or all(any(MolecularProperty.is_same_property(prop.property_type, p) for p in gen.property) for prop in self.needed_properties):
                if gen.compatible_with_details(self.computation.details):
                    prop_gen = gen
                    break
        # and raise an error if we can't find a property generator that will work for us
        if prop_gen is None:
            raise ComputationUnavailableError(
                "Could not find template file capable of generating"
                " an input to compute properties [{0}]"
                " with details given by {1}".format(
                    ", ".join(prop.property_type.__name__ for prop in self.needed_properties),
                    self.computation.details)
            )
        # generate a dict of the known keywords to include in the context...
        known_keywords = {}
        for kw in Keyword._known_roots:
            known_keywords[kw._name] = kw
        # individually pass the required_details variables in as variables at the template
        # scope to improve template readability (even though they could be accessed via
        # computation.details.<detail>
        dets = {}
        for det in prop_gen.details.required_details:
            # Pass both CamelCase and lower_case versions for now
            val = getattr(self.computation.details, det)
            if val is None:
                raise AttributeError("Missing required detail {0} needed for generation "
                                     "of input file {1}".format(det, filename))
            dets[det] = dets[camel_to_lower(det)] = val
        # combined kwargs to pass in...
        pass_in_kwargs = {}
        # give required_details precidence over keyword namespace.  If there is a clash,
        # the programmer/template creator really should rename the required_detail, the
        # keyword, or both.
        pass_in_kwargs.update(known_keywords)
        pass_in_kwargs.update(dets)
        # Only import mako once we need it
        try:
            from mako.template import Template
        except ImportError:
            raise ImportError("Need the Mako package to generate input files.  Try running 'easy_install mako' or 'pip install mako' and try again.")
            # finally, generate the input file
        template = Template(filename=prop_gen.filename)
        ofile = filename
        if isinstance(filename, basestring):
            ofile = open(filename, 'w+')
        ofile.write(template.render(
            computation=self.computation,
            **pass_in_kwargs
        ))
        # ofile will be garbage collected and thus closed if a new file was opened
        self.generated = True


class PropertyGenerator(DetailSpecificObject):
    """ Contains the filename of a template that generates the input file for the computation of a molecular property,
    optionally with a set of specific details for which the template is valid.  This class is mostly to mirror the
    `OutputParser` and `PropertyParser` from the output side of things.
    """

    ##############
    # Attributes #
    ##############

    filename = None
    """ Path to the Mako-style template file. """

    property = None
    """ The property or properties which can be calculated with this template given in `filename`, as a tuple. """

    details = None



    ##################
    # Initialization #
    ##################

    def __init__(self, filename, property, details=None, unsupported_details=None):
        self.filename, self.details = filename, details
        if details is None:
            # Initialize to empty...
            self.details = ComputationDetails()
        if property is All:
            self.property = All
        elif isinstance(property, Iterable):
            self.property = tuple(property)
        else:
            self.property = (property,)
        self.unsupported_details = unsupported_details or {}







class InputGeneratorFactory(object):
    """ A collection of convenience methods for generating InputGenerator classes
    """

    ####################
    # Class Attributes #
    ####################

    generated_number = 0

    #################
    # Class Methods #
    #################

    @classmethod
    def template_based_class(cls, filename, required_details = None, property=All, class_name=None):
        details = ComputationDetails()
        details.required_details = required_details or []
        prop_gens = [PropertyGenerator(filename, property, details)]
        if class_name is None:
            class_name = "_FactoryGeneratedInputGenerator__{0}".format(cls.generated_number)
            cls.generated_number += 1
        return type(class_name, (TemplateInputGenerator,), {'property_generators': prop_gens})


#####################
# Dependent Imports #
#####################

from .computation_details import ComputationDetails

