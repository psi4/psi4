from abc import ABCMeta, abstractmethod
from collections import Iterable
from functools import partial
from itertools import islice
import os
import re
from tempfile import TemporaryFile
import threading
from warnings import warn
import urllib

from grendel import show_warnings
from grendel.chemistry.molecular_properties import Energy, MolecularProperty
from grendel.util.iteration import grouper
from grendel.util.misc import full_path
from grendel.util.overloading import pop_kwarg
from grendel.util.decorators import with_flexible_arguments
from grendel.util.context_managers import working_directory, acquired_lock
from grendel.util.strings import classname, indent, make_safe_identifier

class ComputationFailureError(RuntimeError):
    """ Thrown when a computation fails for one reason or another
    """

class ComputationUnavailableError(RuntimeError):
    """ Raised when a computation can't be generated either because a compatible
    input template can't be found or a valid output parser can't be found
    """

# TODO __str__ and __repr__
# TODO default to not overwrite different input files (with global option to overwrite)
class Computation(object):
    """ Base class for all classes encapsulating runs of external programs
    """
    __metaclass__ = ABCMeta

    ####################
    # Class Attributes #
    ####################

    attribute_aliases = {
        "Theory" : ["Method", "LevelOfTheory"]
    }

    directory_template = None
    always_create_missing = False
    additional_safe_directory_characters = '+'

    ##############
    # Attributes #
    ##############

    molecule = None
    """ The Molecule object used to generate the input.  For single point computations, this (exact object) should
    also be the molecule that any properties refer to."""

    properties = []

    queued = False
    """ Whether or not the calculation is in a queue currently. """

    started = False
    """ True if and only if the calculation has been started. """

    runner = None
    """ The Runner object associated with the computation. """

    details = None
    """ ComputationDetails object that corresponds to this calculation.  Mostly maintained by the Calculation
    object itself, but can be modified to tweek data retrieval (if you really know what you're doing).
    """

    input_generator = None
    output_parser = None
    directory = None

    ######################
    # Private Attributes #
    ######################

    _completed = False


    ##################
    # Initialization #
    ##################

    # TODO Document this!!!!!!
    @with_flexible_arguments(
        required=[
            ('molecule',),
            ('details',),
            ('property', 'properties', 'calculate'),
            ('package', 'program', 'interface'),
        ],
        optional=[
            ('input_generator',),
            ('runner',),
            ('output_parser',),
            ('details_class', 'details_subclass'),
            ('computation_class', 'computation_subclass'),
            ('directory', 'path', 'folder', 'directory_template'),
            ('input_filename', 'input_file'),
            ('output_filename', 'output_file'),
            ('create_directories', 'create_missing', 'create_directory'),
            ('directory_creation_vars', 'directory_template_vars', 'directory_template_variables',
                'directory_creation_variables', 'directory_creation_context'),
            ('result_getter',),
        ],
        what_to_call_it="Computation constructor"
    )
    def __new__(cls, **kwargs):
        """
        """
        def get_kwarg(arg, cls, subclass=False):
            gotten = pop_kwarg(kwargs, arg)
            if not subclass:
                if gotten is not None and not isinstance(gotten, cls):
                    if not isinstance(cls, tuple):
                        raise TypeError("#{0} argument given to Computation() "\
                                        "must be an instance of the #{1} class.".format(arg, classname(cls)))
                    else:
                        raise TypeError("#{0} argument given to Computation() "\
                                        "must be one of (#{1}) .".format(arg, ', '.join(classname(cl) for cl in cls)))
            else:
                if gotten is not None and not issubclass(gotten, cls):
                    if not isinstance(cls, tuple):
                        raise TypeError("#{0} argument given to Computation() "\
                                        "must be an subclass of the #{1} class.".format(arg, classname(cls)))
                    else:
                        raise TypeError("#{0} argument given to Computation() "\
                                        "must be subclass of one of (#{1}) .".format(arg, ', '.join(classname(cl) for cl in cls)))
            return gotten
        #--------------------------------------------------------------------------------#
        # Handle molecule, make a ResultGetter
        molecule = get_kwarg('molecule', Molecule)
        res_getter = get_kwarg('result_getter', ComputationResultGetter)
        if res_getter is None:
            res_getter = ComputationResultGetter(**kwargs)
            molecule.result_getters.append(res_getter)
        #--------------------------------------------------------------------------------#
        # Now see what package we're using
        package = get_kwarg('package', (basestring, PackageInterface))
        comp_class = cls
        interface = None
        if isinstance(package, PackageInterface):
            interface = package
        elif isinstance(package, basestring):
            interface = PackageInterface.load_package(
                            package,
                            input_generator=get_kwarg('input_generator', InputGenerator, True),
                            runner=get_kwarg('runner', Runner, True),
                            output_parser=get_kwarg('output_parser', OutputParser, True),
                            details_class=get_kwarg('details_class', ComputationDetails, True),
                            computation_class=get_kwarg('computation_class', cls, True)
                        )
        if interface.computation_class is not None:
            comp_class = cls
        ret_val = object.__new__(comp_class)
        #--------------------------------------------------------------------------------#
        ret_val.molecule = molecule
        res_getter.add_computation(ret_val)
        ret_val.input_generator = interface.input_generator(ret_val)
        ret_val.output_parser = interface.output_parser(ret_val)
        #TODO accept option for directory to put input/output files in
        #TODO accept option for input/output filenames
        #TODO support for runners that don't follow the input/output pattern
        ret_val.runner = interface.runner(input_file='input.dat', output_file='output.dat')
        #--------------------------------------------------------------------------------#
        ret_val.details = get_kwarg('details', interface.details_class)
        ret_val.properties = pop_kwarg(kwargs, 'property')
        if isinstance(ret_val.properties, Iterable):
            ret_val.properties = list(ret_val.properties)
        else:
            ret_val.properties = [ret_val.properties]
        if not all(issubclass(prop, MolecularProperty) or isinstance(prop, MolecularProperty) for prop in ret_val.properties):
            raise TypeError("property argument given to Computation() must be a list of MolecularProperty subclasses or instances.".format(classname(cls)))
        fixed_properties = []
        for prop in ret_val.properties:
            if isinstance(prop, MolecularProperty):
                fixed_properties.append(prop)
                prop.details = ret_val.details
            elif issubclass(prop, MolecularProperty):
                fixed_properties.append(prop(molecule=ret_val.molecule, details=ret_val.details))
        ret_val.properties = fixed_properties
        #--------------------------------------------------------------------------------#
        # Get input_filename and output_filename, if specified.  Make sure that neither
        #   is a path.
        # first get input_filename
        input_filename = pop_kwarg(kwargs, 'input_filename')
        if input_filename:
            if os.path.split(input_filename)[0] != '':
                raise ValueError("File name '{0}' given for argument 'input_filename' of Computation constructor "
                                 "must be a single file, not a path.  Use the 'directory' argument to specify "
                                 "a directory.".format(input_filename))
            ret_val.runner.input_file = input_filename
        else:
            ret_val.runner.input_file = Runner.default_input_filename
        #----------------------------------------#
        # now get output_filename
        output_filename = pop_kwarg(kwargs, 'output_filename')
        if output_filename:
            if os.path.split(output_filename)[0] != '':
                raise ValueError("File name '{0}' given for argument 'output_filename' of Computation constructor "
                                 "must be a single file, not a path.  Use the 'directory' argument to specify "
                                 "a directory.".format(output_filename))
            ret_val.runner.output_file = output_filename
        else:
            ret_val.runner.output_file = Runner.default_output_filename
        #--------------------------------------------------------------------------------#
        # Find out the desired directory and make new ones if necessary and allowed.
        ret_val.directory = pop_kwarg(kwargs, 'directory')
        used_template = False
        if ret_val.directory is not None:
            ret_val.directory = full_path(ret_val.directory)
        elif cls.directory_template is not None:
            ret_val.directory = cls.directory_template
            used_template = True
        else:
            ret_val.directory = full_path(os.curdir)
        #----------------------------------------#
        # Do mako-template style replacement, if the user has mako
        try:
            from mako.template import Template
            from mako import exceptions
            # Start with some obvious standard things to add...
            variables = dict(
                mol = ret_val.molecule, molecule = ret_val.molecule,
                calc = ret_val, calculation = ret_val,
                rg = res_getter, result_getter = res_getter,
                details = ret_val.details
            )
            # Get all of the known details, as strings
            det_strs = dict()
            for k, v in ret_val.details._known_details.items():
                det_strs[str(k)] = str(v)
                det_strs[str(k).lower()] = str(v)
                det_strs[str(k).upper()] = str(v)
                if k in Computation.attribute_aliases:
                    for alias in Computation.attribute_aliases[k]:
                        alias = make_safe_identifier(alias)
                        det_strs[str(alias)] = str(v)
                        det_strs[str(alias).lower()] = str(v)
                        det_strs[str(alias).upper()] = str(v)
            variables.update(det_strs)
            # Get any additional keywords specified by the user (advanced feature)
            additional = kwargs.pop('directory_creation_vars', dict())
            variables.update(additional)
            # Now run mako on the template...
            try:
                ret_val.directory = Template(
                    text=ret_val.directory,
                    strict_undefined=True
                ).render(
                    **variables
                )
            except:
                raise ValueError(
                    "invalid MAKO template:\n{}\nMako came back with error:\n{}".format(
                        indent(ret_val.directory),
                        indent(exceptions.text_error_template().render())
                    )
                )
            # and make it into a "url safe string":
            ret_val.directory = urllib.quote_plus(
                ret_val.directory,
                safe='/' + cls.additional_safe_directory_characters
            )
        except ImportError:
            # Can't do mako parsing of path.  Oh well...
            # But if the user specified additional args, they probably thought they did have
            #   mako installed...
            additional = kwargs.pop('directory_creation_vars', None)
            if additional is not None:
                raise
            # and if they have a ${ in the directory name or a <% in the directory name,
            #   they probably also thought that they had mako.
            # (Note: this could be an issue in the unlikely scenario that the user is using
            #   e.g. zsh and specifying environment variables from within Python as ${ }.
            #   If someone is having that problem, they probably have enough intelligence
            #   also to figure out how to get around it.)
            if '${' in ret_val.directory or '<%' in ret_val.directory:
                raise
            # if they're using a class-level template, they also need to have mako installed
            if used_template:
                raise
        #----------------------------------------#
        ret_val.directory = full_path(ret_val.directory)
        #----------------------------------------#
        # Now create the directory if we're allowed to and need to
        create_dirs = kwargs.pop('create_directories', cls.always_create_missing)
        if create_dirs:
            if not os.path.isdir(ret_val.directory):
                if os.path.exists(ret_val.directory):
                    raise OSError(
                        "path '{}' exists and is not a directory".format(
                            ret_val.directory
                        )
                    )
                os.makedirs(ret_val.directory)
        else:
            if not os.path.isdir(ret_val.directory):
                raise OSError(
                    "'{0}' is not a valid directory, and will not be created.  Set the "
                    "'create_directories' argument of computation to True to do directory "
                    "creation automatically.".format(
                        ret_val.directory
                    ))
        #----------------------------------------#
        # now set the input_file and output_file values
        ret_val.runner.input_file = os.path.join(ret_val.directory, ret_val.runner.input_file)
        ret_val.runner.output_file = os.path.join(ret_val.directory, ret_val.runner.output_file)
        ret_val.runner.working_directory = ret_val.directory
        #--------------------------------------------------------------------------------#
        # Pass the remaining kwargs on to the initialization, which may have been customized for the particular package.
        ret_val.__init__(**kwargs)
        return ret_val


    ##############
    # Properties #
    ##############

    @property
    def completed(self):
        """ True if and only if the calculation has finished running.  (regardless of success or failure).
        Setting `completed` to True triggers the output parser.
        """
        return self._completed

    @completed.setter
    def completed(self, new_val):
        if not isinstance(new_val, bool):
            raise TypeError
        self._completed = new_val
        if new_val:
            self.output_parser.parse_file(self.runner.output_file)
        else:
            for p in self.properties:
                p.clear_value()

    @property
    def needed_properties(self):
        """ The types of the properties in the `self.properties` list.
        """
        return map(type, self.properties)


    ###################
    # Special Methods #
    ###################


    #################
    # Class Methods #
    #################

    @classmethod
    def standardize_attribute(cls, attr):
        """ Convert strings to standard form for attributes.
        Strings with spaces, underscore_joined_strings, Uppercase_Underscore_Strings, and
        lowerCamelCase strings are all converted to UpperCamelCase.

        Examples
        --------
        >>> Computation.standardize_attribute("basis_set")
        'BasisSet'
        >>> Computation.standardize_attribute("basis set")
        'BasisSet'
        >>> Computation.standardize_attribute("basis!Set")
        'BasisSet'
        >>> # Only _ and space serve as CamelCase separators:
        >>> Computation.standardize_attribute("basis!set")
        'Basisset'
        >>> Computation.standardize_attribute("Level___of!*theory")
        'LevelOftheory'
        >>> # "Accidental" caps are not fixed...
        >>> Computation.standardize_attribute("Level___of!*tHEory")
        'LevelOftHEory'
        >>> #All caps gets translated into title case
        >>> Computation.standardize_attribute("Level&$#of__THEORY")
        'LevelofTheory'

        """
        if len(attr) <= 1:
            ret_val = attr.title()
        else:
            parts = re.split(r'[\s_]+', attr)
            parts = [part for part in parts if len(part) > 0]
            if len(parts) == 0:
                raise NameError("'{0}' is not a valid identifier, and a valid one cannot be constructed from it".format(attr))
            # make ALLCAPS into title case, otherwise make only the first letter capital and ignore other capitalization
            ret_val = "".join(
                (part.title() if part.isupper() or len(part) == 1 else (part[0].upper() + part[1:])) for part in parts
            )
        # clear out any other junk
        ret_val = "".join(re.findall(r'[A-Za-z0-9]+', ret_val))
        if len(ret_val) == 0:
            raise NameError("'{0}' is not a valid identifier, and a valid one cannot be constructed from it".format(attr))
        for k in cls.attribute_aliases:
            v = cls.attribute_aliases[k]
            if isinstance(v, Iterable):
                if ret_val in v:
                    ret_val = k
                    break
            else:
                if ret_val == v:
                    ret_val = k
        return ret_val


    ####################
    # Abstract Methods #
    ####################

    # nothing here!

    ###########
    # Methods #
    ###########

    def needs_property(self, prop):
        if isinstance(prop, type):
            return any(MolecularProperty.is_same_property(prop, p) for p in self.needed_properties)
        else:
            prop = eval(prop)
            return any(MolecularProperty.is_same_property(prop, p) for p in self.needed_properties)

    def has_energy(self, details = None):
        return self.has_property(Energy, details) is not None

    def get_energy(self, details = None):
        return self.get_property(Energy, details)

    def has_property(self, prop_class, details=None):
        if not self.completed:
            return False
        for p in self.properties:
            if MolecularProperty.is_same_property(p, prop_class) \
                    and (details is None or details.is_subset_of(p.details)):
                return True
        return False

    def get_property(self, prop_class, details=None):
        ret_val = None
        if not self.completed:
            return None
        for p in self.properties:
            # TODO figure out which should be a superset of the other
            if MolecularProperty.is_same_property(p, prop_class) and (details is None or details.is_subset_of(p.details)):
                if not ret_val is None and not p.value == ret_val.value:
                    # TODO Raise a more helpful exception class
                    raise RuntimeError("Multiple conflicting values for property " + classname(prop_class))
                ret_val = p
        return ret_val

    def start(self):
        raise NotImplementedError
        # TODO make this work with the working_directory change thing...
        #self._input_check()
        #self.runner.start()
        #self.started = True
        #return self.runner

    def get_results(self):
        if not self.completed:
            raise RuntimeError("computation not yet completed.  Cannot get results.")
        else:
            self.output_parser.parse_file(self.runner.output_file)

    def run(self, force_rerun=False, lock=None):
        """
        """
        #TODO figure out if the lock is even necessary.  There should be at worst a 1-to-1 mapping between threads and Computation instances
        if lock is None:
            # Just create a dummy lock
            lock = threading.Lock()
        with acquired_lock(lock):
            already_run = self.already_run()
            completed = self.completed
            runner = self.runner
        if not force_rerun and already_run:
            with acquired_lock(lock):
                if show_warnings:
                    warn(
                        "Using already-run version of calculation in directory {0} since"
                        " generated input would be identical to the one already existing"
                        " and output appears to have completed successfully.  To force"
                        " rerun, set the 'force_rerun' argument of 'Computation.run()'"
                        " to True.".format(self.directory)
                    )
                self.started = True
                self.completed = True
                self.__output_check()
        elif not completed or force_rerun:
            with acquired_lock(lock):
                self.__input_check()
                self.started = True
            runner.run()
            with acquired_lock(lock):
                self.completed = True
                self.__output_check()
        else:
            with acquired_lock(lock):
                if show_warnings:
                    warn("calculation already ran and will not be rerun.  To force rerun, set the 'force_rerun' argument of `Computation.run()` to True")
        return runner

    def already_run(self):
        """ Checks to see if an identical input file with an identical name has already been run in the computation's directory
        to produce an output file of some sort (if that output file has errors, you probably don't want to rerun the *exact same* input file,
        so this doesn't do that by default.  If you need to do that, call `Computation.run()` with the `force_rerun` argument
        set to True).
        """
        if not os.path.isfile(self.runner.input_file):
            # Input file must exist, at the very least...
            return False
        # generate a temporary version of the input file.
        with TemporaryFile() as tmp:
            self.input_generator.generate(tmp)
            existing = open(self.runner.input_file, 'r').read()
            tmp.seek(0)
            new = tmp.read()
            if new == existing:
                if os.path.isfile(self.runner.output_file):
                    if hasattr(self.output_parser, "has_errors"):
                        return self.output_parser.has_errors(open(self.runner.output_file))
                    else:
                        return True
                else:
                    return False
            else:
                # otherwise, put the new input file in place of the old...
                #TODO option to backup the old input file (and global option to backup/systematically stash all files before overwriting, throughout the program)
                with open(self.runner.input_file, 'w+') as f:
                    f.write(new)
                return False

    def task_description(self):
        # TODO add some more helpful information to this, with an argument for verbosity level
        ret_val = ""
        basics = [
            "Path", self.directory,
            "Input File Name", self.runner.input_file.split(os.sep)[-1],
            "Output File Name", self.runner.output_file.split(os.sep)[-1],
        ]
        extras = [
            "Input Generator Class", self.input_generator.__class__.__name__,
            "Output Parser Class", self.output_parser.__class__.__name__,
            "Runner Class", self.runner.__class__.__name__,
        ]
        allstats = basics + extras
        lwidth = max(len(k) for k, v in grouper(2, allstats))
        rwidth = max(len(v) for k, v in grouper(2, allstats))
        props = []
        for name, value in grouper(2, allstats):
            props.append(('{:.<'+str(lwidth)+'}{:.>'+str(rwidth)+'}').format(
                name, value
            ))
        return '\n'.join(props)

    ###################
    # Private Methods #
    ###################

    def __input_check(self):
        with working_directory(self.directory):
            if not self.input_generator.generated:
                self.input_generator.generate(self.runner.input_file)

    def __output_check(self):
        if not os.path.isfile(self.runner.output_file):
            raise ComputationFailureError("Output file '{0}' not generated by run of {1}".format(
                self.runner.output_file,
                ' '.join(self.runner.command_sequence)
            ))
        elif hasattr(self.output_parser, 'failure_check'):
            # If the output parser knows how to check for errors, let it do its thing
            self.output_parser.failure_check(self.runner.output_file)

#####################
# Dependent Imports #
#####################

from grendel.chemistry.molecule import Molecule
from grendel.chemistry.molecular_properties import *
from grendel.chemistry.derivative_properties import *
from grendel.interface.computation_details import ComputationDetails
from grendel.interface.package_interface import PackageInterface
from grendel.interface.input_generator import InputGenerator
from grendel.interface.output_parser import OutputParser
from grendel.interface.runner import Runner
from grendel.interface.result_getter import ComputationResultGetter



