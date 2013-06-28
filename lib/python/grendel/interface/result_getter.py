"""
"""
from abc import ABCMeta, abstractmethod
from copy import copy

from grendel.chemistry.molecular_properties import MolecularProperty
from grendel.util.aliasing import function_alias


class ResultGetter(object):
    """ Abstract base class for all classes that encapsulate a method of retrieving a (set of) properties of
    a molecule from a computational chemistry package.  An instance of a ResultGetter object should know how to
    turn a request for a given molecular property into a Computation object (which, in turn, houses MolecularProperty
    objects for one or more Molecule objects).  **ResultGetter objects generate Computation objects.** *Any
    ResultGetter subclass should always check for an available computation on (exactly) the same molecule that
    already has the result available.*
    """
    __metaclass__ = ABCMeta

    ####################
    # Abstract Methods #
    ####################

    @abstractmethod
    def get_property_for_molecule(self, molecule, property, details=None):
        return NotImplemented

    @abstractmethod
    def can_get_property_for_molecule(self, molecule, property, details=None):
        return NotImplemented

    @abstractmethod
    def has_property_for_molecule(self, molecule, property, details=None, verbose=False):
        return NotImplemented


class ComputationResultGetter(ResultGetter):
    """
    """

    ####################
    # Class Attributes #
    ####################

    known_result_getters = {}
    directory_series = {}

    ##############
    # Attributes #
    ##############

    computations = None
    computation_kwargs = None
    use_directory_series = None

    ##################
    # Initialization #
    ##################

    def __new__(cls, use_directory_series=False, **kwargs):
        """
        """
        # First set ourselves up and parse our own kwargs...

        # Now remove the molecule property, if it's present, so the kwargs are usable for lots of stuff...
        kwargs.pop('molecule', None)
        # Check to see if an identical ComputationResultGetter exists
        try:
            key = frozenset((k, v) for k, v in kwargs.iteritems())
        except TypeError:
            key = frozenset((id(k), id(v)) for k, v in kwargs.iteritems())
        if key in ComputationResultGetter.known_result_getters:
            ret_val = ComputationResultGetter.known_result_getters[key]
        else:
            ret_val = object.__new__(cls)
        ret_val.use_directory_series = use_directory_series
        if use_directory_series:
            raise NotImplementedError
        # Then pass the remaining kwargs to computations when they are generated...
        ret_val.computation_kwargs = copy(kwargs)
        ComputationResultGetter.known_result_getters[key] = ret_val
        ret_val.computations = []
        return ret_val

    ###########
    # Methods #
    ###########

    def can_get_property_for_molecule(self, molecule, property, details=None):
        try:
            comp = self._get_comp(molecule, property, details, True)
            if comp.runner is not None:
                comp.runner.validate()
            if comp.input_generator is not None:
                comp.input_generator.validate()
            if comp.output_parser is not None:
                comp.output_parser.validate()
        except ComputationUnavailableError:
            return False
        return True
    can_get_property = function_alias('can_get_property', can_get_property_for_molecule)

    def has_property_for_molecule(self, molecule, property, details=None, verbose=False):
        comp = self._get_comp(molecule, property, details, False)
        if comp is None:
            return False
        return comp.has_property(property, details)
    has_property = function_alias('has_property', can_get_property_for_molecule)

    def get_computation_for_property(self, molecule, property, details=None):
        comp = self._get_comp(molecule, property, details, True)
        # Try to generate any errors that would cause can_get_property_for_molecule to fail...
        if comp.runner is not None:
            comp.runner.validate()
        if comp.input_generator is not None:
            comp.input_generator.validate()
        if comp.output_parser is not None:
            comp.output_parser.validate()
        # Return it, whether or not it has been run.
        return comp

    def get_property_for_molecule(self, molecule, property, details=None):
        comp = self.get_computation_for_property(molecule, property, details)
        # See if we already have the property...
        result = comp.get_property(property, details)
        if result is not None:
            return result
        else:
            # Now run the calculation
            comp.run()
            return comp.get_property(property, details)

    def add_computation(self, comp):
        self.computations.append(comp)

    ###################
    # Private Methods #
    ###################

    def _get_comp(self, molecule, property, details, generate):
        for comp in self.computations:
            # TODO offer both "is" and "==" (as well as "is_same_molecule") options
            if comp.molecule is molecule:
                for prop in comp.properties:
                    if MolecularProperty.is_same_property(prop, property):
                        if details is None or details.is_subset_of(prop.details):
                            return comp
        # well, we couldn't find one.  Should we make one?
        if generate:
            if 'property' in self.computation_kwargs:
                prop = self.computation_kwargs['property']
                if not MolecularProperty.is_same_property(prop, property):
                    raise ComputationUnavailableError("property mismatch.  Can't generate computation"
                                                      " of property '{}'".format(
                        MolecularProperty.property_type_of(property).__name__
                    ))
                # This will append the computation to the self.computations list automatically
                comp = Computation(molecule=molecule, result_getter=self, **self.computation_kwargs)
            else:
                # This will append the computation to the self.computations list automatically
                comp = Computation(molecule=molecule, property=property, result_getter=self, **self.computation_kwargs)
            return comp
        else:
            return None


class OptimizationResultGetter(ResultGetter):

    ####################
    # Class Attributes #
    ####################

    known_optimizers = {}

    ####################
    # Abstract Methods #
    ####################

    @abstractmethod
    def get_optimized_property_for_molecule(self, molecule, property, details=None):
        return NotImplemented

    @abstractmethod
    def can_optimize_property_for_molecule(self, molecule, property, details=None):
        return NotImplemented

    @abstractmethod
    def has_optimized_property_for_molecule(self, molecule, property, details=None, verbose=False):
        return NotImplemented


class ComputationOptimizationResultGetter(OptimizationResultGetter):


    ##############
    # Attributes #
    ##############

    computations = None

    ##################
    # Initialization #
    ##################

    def __new__(cls, **kwargs):
        """
        """
        # First set ourselves up and parse our own kwargs...

        # Now remove the molecule property, if it's present, so the kwargs are usable as a key
        kwargs.pop('molecule', None)
        # Check to see if an identical ComputationResultGetter exists
        try:
            key = frozenset((k, v) for k, v in kwargs.iteritems())
        except TypeError:
            key = frozenset((id(k), id(v)) for k, v in kwargs.iteritems())
        if key in ComputationResultGetter.known_result_getters:
            ret_val = ComputationResultGetter.known_optimizers[key]
        else:
            ret_val = object.__new__(cls)
        ret_val.computation_kwargs = copy(kwargs)
        ComputationResultGetter.known_result_getters[key] = ret_val
        ret_val.computations = []
        return ret_val

    ###########
    # Methods #
    ###########

    def get_property_for_molecule(self, molecule, property, details=None):
        return False

    def can_get_property_for_molecule(self, molecule, property, details=None):
        return False

    def has_property_for_molecule(self, molecule, property, details=None, verbose=False):
        return False

    def get_optimized_property_for_molecule(self, molecule, property, details=None):
        raise NotImplementedError

    def can_get_property_for_molecule(self, molecule, property, details=None):
        raise NotImplementedError

    def has_optimized_property_for_molecule(self, molecule, property, details=None, verbose=False):
        raise NotImplementedError


#####################
# Dependent Imports #
#####################

from grendel.interface.computation import Computation, ComputationUnavailableError

