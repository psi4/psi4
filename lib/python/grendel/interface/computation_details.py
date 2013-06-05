from collections import Iterable
from copy import deepcopy, copy
from exceptions import TypeError
import grendel
from grendel.util.overloading import listify_args

from grendel.util.sentinal_values import Keyword
from grendel.util.units import *

__all__ = [
    "ComputationDetails",
    "Methods",
    "Bases", "Basis"
]

Methods = Keyword('Methods')
Bases = Basis = Keyword('Basis')

class ComputationDetails(object):
    """ Encapsulates specifics about how the calculation was/is to be carried out,
    such as level of theory, basis set, program, etc.  Can be subclassed for specific
    programs which may require special detail setting methods.
    """

    ##############
    # Attributes #
    ##############

    required_details = []
    """ Options that must be specified before the computation can be run.  These can be given by the user, but more often
    they are automatically set when other options or properties are set.  For instance, if the user specifies Theory
    to be "DFT", then "Exchange" is added to required_options.  If Exchange is not specified, the calculation will not run.
    """

    ######################
    # Private Attributes #
    ######################

    _known_details = None


    ##################
    # Initialization #
    ##################

    def __init__(self, **kwargs):
        self._known_details = {}
        for key in kwargs:
            # Let __setattr__ handle it
            setattr(self, key, kwargs[key])


    ###################
    # Special Methods #
    ###################

    # Set up setattr and getattr to mimic an open struct
    def __setattr__(self, key, value):
        if key in type(self).__dict__:
            self.__dict__.update({key: value})
        else:
            name = Computation.standardize_attribute(key)
            if "_{0}_set".format(name) in self.__class__.__dict__:
                self._known_details[name] = getattr(self, "_{0}_set".format(name))(value)
            else:
                self._known_details[name] = value

    def __getattr__(self, item):
        # Remember, this is only called if item is not in self.__dict__
        name = Computation.standardize_attribute(item)
        if name not in self._known_details:
            return None
        return self._known_details[name]

    def __contains__(self, item):
        return Computation.standardize_attribute(item) in self._known_details

    def __getitem__(self, item):
        return getattr(self, item)

    def __setitem__(self, key, value):
        setattr(self, key, value)

    def __str__(self):
        return "<ComputationDetails:\n    {0}\n>".format("\n    ".join(
            (Computation.standardize_attribute(key) + " = " + str(self._known_details[key])) for key in self._known_details
        ))
    __repr__ = __str__

    def __deepcopy__(self, memo):
        # Since exactly one instance of every keyword may exist, there is no reason to copy the values of
        # known_details (they can also be strings, which are of course immutable).  There is also no reason to
        # copy the keys, since they are strings.  Thus, all we need to do is:
        ret_val = self.__class__(**self._known_details)
        # Note that deepcopy on a string is identical to copy since strings are immutable
        ret_val.required_details = deepcopy(self.required_details, memo)
        return ret_val

    def __copy__(self):
        # Since exactly one instance of every keyword may exist, there is no reason to copy the values of
        # known_details (they can also be strings, which are of course immutable).  There is also no reason to
        # copy the keys, since they are strings.  Thus, all we need to do is:
        ret_val = self.__class__(**self._known_details)
        ret_val.required_details = copy(self.required_details)
        return ret_val

    #----------------#
    # Detail Setting #
    #----------------#

    def _Method_set(self, value):
        # Make all methods upper case for now...
        v = value.upper()
        if v == 'DFT':
            self.required_details.append("Exchange")
            return Methods.DFT
        elif v == 'HF':
            return Methods.HF


    #################
    # Class Methods #
    #################

    @classmethod
    def is_compatible_details(cls, requested, details_computed_with):
        return requested is None or requested.is_subset_of(details_computed_with)

    @classmethod
    def from_details(cls, details):
        """ Effectively clones details.
        Used internally for making a details into a ComputationDetails subclass
        """
        ret_val = cls(details._known_details)
        ret_val.other_options = details.other_options
        ret_val.required_details = details.required_details
        return ret_val

    ###########
    # Methods #
    ###########

    # TODO call _compatible_values instead of this mess
    def is_superset_of(self, other):
        if grendel.type_checking_enabled and not isinstance(other, ComputationDetails):
            raise TypeError
        for key in self._known_details:
            value = self._known_details[key]
            if not hasattr(other, key):
                return False
            elif isinstance(value, Keyword):
                if not isinstance(getattr(other, key), Keyword):
                    return False
                elif not (value == getattr(other, key) or value in getattr(other, key).parents):
                    return False
            elif isinstance(value, Iterable):
                found = False
                for item in value:
                    if isinstance(item, Keyword):
                        if isinstance(getattr(other, key), Keyword) and (item == getattr(other, key) or item in getattr(other, key).parents):
                            found = True
                    elif not isinstance(getattr(other, key), Keyword) and item == getattr(other, key):
                            found = True
                if not found:
                    return False
            else:
                if value != getattr(other, key):
                    return False
        return True

    # NOTE:  THIS IS NOT THE INVERSE OF is_superset_of!!!
    # TODO Think this through more clearly
    def is_subset_of(self, other):
        if self._known_details is None or len(self._known_details) == 0:
            return True
        for key in self._known_details:
            value = self._known_details[key]
            if not key in other:
                return False
            elif not self._compatible_values(value, other._known_details[key]):
                return False
        return True




    #TODO optional 'num_per_line' kwarg
    #TODO optional 'separator' kwarg (i.e. for num_per_line > 1)
    #TODO optional 'indent' kwarg
    #TODO error checking on the format string (make sure it contains {key} and {value})
    def keywordify_if_available(self, format, *keys):
        """ Mostly for use by templates, this method returns a string of keyword-value pairs in using the given format
        for each argument in the argument list.
        `format` should use the str.format() protocol from the python standard library with the
        named slots `{key}` where the keyword name should go and `{value}` where the value of the
        detail should go (or, if the value of the detail is a `Keyword` object, the contents of the
        `value` attribute for that `Keyword`, or, if the `value` attribute for the keyword is not defined,
        the `_name` attribute for the `Keyword` object in question.

        Aliased as `keywordify`
        """
        ret_val = ""
        for key in listify_args(keys):
            if key in self:
                value = self[key]
                if isinstance(value, Keyword):
                    if value.value is not None:
                        ret_val += format.format(key, value.value) + "\n"
                    else:
                        ret_val += format.format(key, value._name) + "\n"
                else:
                    ret_val += format.format(key, value) + "\n"
        return ret_val
    keywordify = keywordify_if_available

    def available_pairs(self, *keys):
        for key in listify_args(keys):
            if key in self:
                value = self[key]
                if isinstance(value, Keyword):
                    if value.value is not None:
                        yield key, str(value.value)
                    else:
                        yield key, value._name
                else:
                    yield key, str(value)

    ###################
    # Private Methods #
    ###################

    def _compatible_values(self, parent, child):
        if isinstance(parent, Keyword):
            if not isinstance(child, Keyword):
                return False
            elif not parent == child or parent in child.parents:
                return False
        elif isinstance(parent, basestring):
            if isinstance(child, basestring):
                return parent == child
            else:
                return False
        elif isinstance(parent, Iterable):
            found = False
            for item in parent:
                if isinstance(item, Keyword):
                    if isinstance(child, Keyword) and (item == child or item in child.parents):
                        found = True
                elif not isinstance(child, Keyword) and item == child:
                    found = True
            if not found:
                return False
        else:
            if parent != child:
                return False
        return True



class DetailSpecificObject(object):
    """ Superclass for classes that are compatible with specific ComputationDetail objects
    or any ComputationDetail object more specifically specifice than the `details` variable.
    """

    ##############
    # Attributes #
    ##############

    details = None
    """ The ComputationDetails object describing the scope of validity for the object.
    An empty ComputationDetails object indicates that the subclass is valid for all properties (not likely to
    be the case)
    """

    unsupported_details = None
    """ Features not supported by the given template (i.e. features that must not be in a details object to match)."""


    ###########
    # Methods #
    ###########

    def compatible_with_details(self, details):
        if not self.details.is_superset_of(details):
            return False
        else:
            for key in (self.unsupported_details or []):
                val = self.unsupported_details[key]
                if key in details:
                    if isinstance(val, Keyword):
                        if isinstance(details[key], Keyword):
                            if val == details[key] or val in details.parents:
                                return False
                    elif val == details[key]:
                        return False
            # Otherwise, it's not a match; keep looking for unsupported details
            # No unsupported details found; we're good
            return True


#####################
# Dependent Imports #
#####################

from .computation import Computation