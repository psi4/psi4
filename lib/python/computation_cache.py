from __future__ import print_function, division
from abc import abstractmethod, abstractproperty, ABCMeta
from collections import Iterable, defaultdict
from inspect import getargspec
import inspect
from itertools import product
from os import makedirs, mkdir, remove, unlink
import os
from os.path import exists, isdir, split as path_split, sep as path_sep, join as path_join
import re
import shelve
import shutil
import tempfile
import traceback
from warnings import warn
import numpy as np
from numpy.lib.format import read_array_header_1_0, open_memmap, read_magic
import sys
from grendel.chemistry import Molecule
from grendel.util.containers import AliasedDict
from grendel.util.decorators import typecheck, LazyProperty
from functools import partial
from hashlib import md5
import psi4
import atexit
# Create a trivial molecule so psi4 doesn't freak out about not having a Molecule
#psi4.geometry("H 0 0 0\nH 0 0 1")

RegexType = type(re.compile(''))

#TODO test suite
#TODO move to grendel

xproduct = lambda *args: product(*[xrange(i) for i in args])
lambda_type = type(lambda x: x)

class CacheInconsistancyWarning(Warning):
    pass

class FileMissingWarning(CacheInconsistancyWarning):
    pass

class ShapeMismatchWarning(Warning):
    pass

class FileOverwriteWarning(Warning):
    pass

enabled_warnings = AliasedDict({
    (ShapeMismatchWarning, "ShapeMismatchWarning") : True,
    (FileOverwriteWarning, "FileOverwriteWarning") : True,
    (FileMissingWarning, "FileMissingWarning") : True
})
fail_warnings = AliasedDict({
    (ShapeMismatchWarning, "ShapeMismatchWarning") : False,
    (FileOverwriteWarning, "FileOverwriteWarning") : False,
    (FileMissingWarning, "FileMissingWarning") : False
})

#================================================================================#

#region | Warning handling |

#TODO move this to grendel, generalize

def enable_warning(w):
    enabled_warnings[w] = True

def disable_warning(w):
    enabled_warnings[w] = False

def enable_all_warnings():
    for w in enabled_warnings:
        enabled_warnings[w] = True

def disable_all_warnings():
    for w in enabled_warnings:
        enabled_warnings[w] = False

def raise_on_warning(w):
    fail_warnings[w] = True

def raise_warning(msg, warning_type):
    if enabled_warnings[warning_type]:
        warn(msg, warning_type)
    if fail_warnings[warning_type]:
        raise RuntimeError("Critical warning raised.")

#endregion

#================================================================================#

class UninitializedDependency(object):

    def __init__(self,
            init_func,
            owner,
            depends_on=None,
            requires_optional_arguments=None
    ):
        self.init_func = init_func
        self.owner = owner
        self.initialized = False
        if depends_on is None: depends_on = []
        self.depends_on = depends_on
        if requires_optional_arguments is None: requires_optional_arguments = []
        self.requires_optional_arguments = requires_optional_arguments
        for req in self.requires_optional_arguments:
            if not ComputationCache.optional_attributes[req] == ComputationCache.NoDefaultValue:
                raise ValueError("Can't require that a user-specified value be given for an optional"
                                 " attribute that has a default value.")


    def initialize(self, memo=None):
        if self.initialized:
            # Should be unreachable, unless someone is messing with things they shouldn't
            raise ValueError("Call to used UninitializedDependency initialize function")
        #----------------------------------------#
        # Check for circular dependencies
        if memo is None:
            memo = [self]
        else:
            if self in memo:
                raise ValueError("Circular Dependency")
            memo.append(self)
        #----------------------------------------#
        # Initialize dependencies
        for dep in self.depends_on:
            if not hasattr(self.owner, dep):
                raise AttributeError("Unknown dependent attribute: {0}".format(dep))
            attr = getattr(self.owner, dep)
            if isinstance(attr, UninitializedDependency):
                setattr(self.owner, dep, attr.initialize(memo))
        #----------------------------------------#
        # Check to make sure all of the necessary
        #   optional arguments (if any) have been
        #   specified by the user in the cached
        #   computation self.owner
        for req in self.requires_optional_arguments:
            if self.owner.optional_arguments[req] == ComputationCache.NoDefaultValue:
                # Try to get the name of the attribute we're looking for, for error purposes
                my_stack = inspect.stack()
                dep_name = None
                # First, see if we called from another uninitialized dependency
                if my_stack[1][1] == __file__ and my_stack[1][3] == "initialize":
                    caller_locals = inspect.getargvalues(my_stack[1][0]).locals
                    if 'dep' in caller_locals:
                        dep_name = caller_locals['dep']
                # Or check to see if we were called from get_lazy_attribute
                elif my_stack[1][1] == __file__ and my_stack[1][3] == "get_lazy_attribute":
                    argvals = inspect.getargvalues(my_stack[1][0])
                    argname = argvals.args[1]
                    caller_locals = argvals.locals
                    if argname in caller_locals:
                        dep_name = caller_locals[argname]
                # Then raise the (hopefully helpful) error
                raise ValueError("Initialization of attribute '{0}' requires optional CachedComputation"
                                 " argument '{1}' to have a user-defined value".format(
                    dep_name if dep_name is not None else "<unknown attribute>",
                    req
                ))
        #----------------------------------------#
        rv = self.init_func(self.owner)
        self.initialized = True
        return rv

    def __getstate__(self):
        raise TypeError("UninitializedDependency instances are not pickleable")

#================================================================================#

#region | DatumGetter subclasses |

class DatumGetter(object):
    __metaclass__ = ABCMeta

    @property
    def needs(self):
        if hasattr(self, "_needs"):
            return self._needs
        else:
            return getargspec(self.__call__).args[2:]
    @needs.setter
    def needs(self, value):
        self._needs = value

    @abstractmethod
    def __call__(self, *args, **kwargs):
        pass

    def __getstate__(self):
        raise TypeError("Instances of DatumGetter and its subclasses should never be pickled.")

class SimpleMethodCallGetter(DatumGetter):

    def __init__(self, obj, method_name, *args):
        self._needs = [obj]
        self.method_name = method_name
        self.method_args = args

    def __call__(self, **kwargs):
        if len(kwargs) > 1:
            raise TypeError("SimpleMethodCallGetter requires exactly one argument")
        obj = kwargs.values()[0]
        return getattr(obj, self.method_name)(*self.method_args)

class DimensionListGetter(SimpleMethodCallGetter):

    def __call__(self, **kwargs):
        if len(kwargs) > 1:
            raise TypeError("SimpleMethodCallGetter requires exactly one argument")
        obj = kwargs.values()[0]
        dimrv = getattr(obj, self.method_name)(*self.method_args)
        size = dimrv.n()
        rv = []
        for i in range(size):
            rv.append(dimrv[i])
        return rv

class MemmapArrayGetter(DatumGetter):
    @abstractmethod
    def get_shape(self, *args):
        pass

class TEIGetter(MemmapArrayGetter):
    """
    Abstract parent class of AOTEIGetter and MOTEIGetter
    """

    ####################
    # Class Attributes #
    ####################

    VALID_KERNELS = (
        "eri",
        "erf_eri",
        "f12",
        "f12_squared",
        "f12g12",
        "f12_double_commutator"
    )

    ##################
    # Initialization #
    ##################

    def __init__(self, kernel_name, physicist_notation=False):
        self.kernel_name = kernel_name
        if self.kernel_name not in TEIGetter.VALID_KERNELS:
            raise NameError("Don't know how to get two electron integrals for kernel '{0}'".format(
                self.kernel_name
            ))
        if self.kernel_name == "eri":
            # Make a special needs that doesn't include the correlation factor
            self.needs = [n for n in self.needs if n != "correlation_factor"]
        self.physicist_notation = physicist_notation

    ###################
    # Private Methods #
    ###################

    def _copy_in_matrix(self, out, psimat, physicist_notation=False):
        # TODO copy memory instead, if possible
        # TODO symmetry
        nbf1, nbf2, nbf3, nbf4 = out.value.shape
        ary = np.array(psimat)
        out.value[...] = ary.reshape(out.shape)
        if physicist_notation:
            out.value[...] = out.value.transpose([0,2,1,3])
        #if physicist_notation:
        #    for p, q, r, s in xproduct(*out.shape):
        #        out.value[p, r, q, s] = psimat[0, p*nbf2 + q, r*nbf4 + s]
        #else:
        #    for p, q, r, s in xproduct(*out.shape):
        #        out.value[p, q, r, s] = psimat[0, p*nbf2 + q, r*nbf4 + s]

    def _get_tei(self, out, func_name, mints, correlation_factor, physicist_notation, C1=None, C2=None, C3=None, C4=None):
        getter_args = []
        if C1 is not None:
            if C2 is None: C2 = C1
            if C4 is None:
                if C3 is not None: C4 = C3
                else: C4 = C1
            if C3 is None: C3 = C1
            getter_args.extend([C1, C2, C3, C4])
        mints_getter = getattr(mints, func_name)
        if "f12" in func_name:
            getter_args.insert(0, correlation_factor)
        elif "erf" in func_name:
            raise NotImplementedError()
        psimat = mints_getter(*getter_args)
        self._copy_in_matrix(out, psimat, physicist_notation)
        return out

class AOTEIGetter(TEIGetter):

    def get_shape(self, basis):
        nbf = basis.nbf()
        return (nbf,)*4

    def __call__(self, out, mints, correlation_factor=None):
        self._get_tei(out,
            "ao_" + self.kernel_name,
            mints, correlation_factor, self.physicist_notation
        )

class MOTEIGetter(TEIGetter):

    # TODO noncanonical transformation

    def get_shape(self, reference_wavefunction):
        nmo = reference_wavefunction.nmo()
        return (nmo,)*4

    def __call__(self, out, mints, mo_coefficients, correlation_factor=None):
        self._get_tei(out,
            "mo_" + self.kernel_name,
            mints, correlation_factor, self.physicist_notation,
            mo_coefficients
        )

class MatrixGetter(MemmapArrayGetter):

    def _copy_in_matrix(self, out, psimat):
        # TODO Symmetry
        #out.value[...] = np.array(psimat)
        for p, q in xproduct(*out.shape):
            out.value[p, q] = psimat[0, p, q]

class MOCoefficientsGetter(MatrixGetter):
    """
    """

    #TODO open shells

    ####################
    # Class Attributes #
    ####################

    ALLOWED_TYPES = [
        "canonical",
        "boys",
        "pipek_mezey"
    ]

    ##################
    # Initialization #
    ##################

    def __init__(self, orbital_type="canonical"):
        self.orbital_type = orbital_type
        if self.orbital_type not in MOCoefficientsGetter.ALLOWED_TYPES:
            raise ValueError("Unknown orbital type '{0}'".format(orbital_type))

    def get_shape(self, reference_wavefunction, basis):
        return basis.nbf(), reference_wavefunction.nmo()

    def _get_coeffs(self, out, basis, C):
        if self.orbital_type == "canonical":
            self._copy_in_matrix(out, C)
        elif self.orbital_type == "boys" or self.orbital_type == "pipek_mezey":
            localizer = psi4.Localizer.build(self.orbital_type.upper(), basis, C)
            localizer.localize()
            if not localizer.converged:
                raise RuntimeError("Localization didn't converge")
            self._copy_in_matrix(out, localizer.L)
        else:
            raise NotImplementedError()

    def __call__(self, out, basis, mo_coefficients):
        return self._get_coeffs(out, basis, mo_coefficients)

class OccupiedMOCoefficientsGetter(MOCoefficientsGetter):

    def get_shape(self, reference_wavefunction, basis):
        return basis.nbf(), sum(reference_wavefunction.doccpi())

    def __call__(self, out, basis, docc_space):
        self._get_coeffs(out, basis, docc_space.C())

class OEIGetter(MemmapArrayGetter):

    def __init__(self, kernel):
        self.kernel = kernel
        if not kernel in self.__class__.VALID_KERNELS:
            raise NameError("Don't know how to get one electron integrals for kernel '{0}'".format(
                self.kernel
            ))

#DEPRECATED: Chunks data incorrectly
class AOMultipoleGetter(MemmapArrayGetter):
    """
    NOTE: This getter chunks the data incorrectly
    """

    def __init__(self, max_order):
        self.max_order = max_order

    def get_shape(self, basis):
        o = self.max_order
        return basis.nbf(), basis.nbf(), int((o+1)*(o+2)*(o+3)/6 - 1)

    def __call__(self, out, factory, basis):
        intobj = factory.ao_multipoles(self.max_order)
        nsh = basis.nshell()
        for ish, jsh in xproduct(nsh, nsh):
            ioff, joff = map(basis.shell_to_basis_function, (ish, jsh))
            inbf, jnbf = map(lambda sh: basis.shell(sh).nfunction, (ish, jsh))
            intobj.compute_shell(ish, jsh)
            out.value[ioff:ioff+inbf, joff:joff+jnbf,:] = np.array(intobj.py_buffer).reshape(
                (inbf,jnbf,out.value.shape[2])
            )

class AOOEIGetter(OEIGetter, MatrixGetter):

    VALID_KERNELS = [
        "overlap",
        "kinetic",
        "potential"
    ]

    def get_shape(self, basis):
        nbf = basis.nbf()
        return nbf, nbf

    def __call__(self, out, mints):
        mints_getter = getattr(mints, "ao_" + self.kernel)
        psimat = mints_getter()
        self._copy_in_matrix(out, psimat)
        return out

class SOOEIGetter(OEIGetter, MatrixGetter):

    VALID_KERNELS = [
        "overlap",
        "kinetic",
        "potential"
    ]

    def get_shape(self, reference_wavefunction):
        nso = reference_wavefunction.nso()
        return nso, nso

    def __call__(self, out, mints):
        mints_getter = getattr(mints, "so_" + self.kernel)
        psimat = mints_getter()
        self._copy_in_matrix(out, psimat)
        return out

class AOTEIGetterWithComputer(MemmapArrayGetter):

    VALID_INTEGRAL_TYPES = [
        "eri", "f12", "f12g12",
        "f12_double_commutator",
        "f12_squared", "erf_eri",
        "erf_complement_eri"
    ]

    def __init__(self, integral_type):
        self.integral_type = integral_type
        if self.integral_type not in type(self).VALID_INTEGRAL_TYPES:
            raise ValueError("Don't know how to compute integrals of type '{0}'".format(
                integral_type
            ))
        # Make a special needs for the ERI doesn't include the correlation factor
        if self.integral_type == "eri":
            self.needs = [n for n in self.needs if n != "correlation_factor"]
        self.zero_basis_set = lambda: psi4.BasisSet.zero_ao_basis_set()

    def get_shape(self, basis):
        nbf = basis.nbf()
        return nbf, nbf, nbf, nbf

    def __call__(self, out, basis, correlation_factor=None):
        self._compute_ints(out,
            basis_sets=(basis, basis, basis, basis),
            correlation_factor=correlation_factor
        )

    @staticmethod
    def shell_slice(bs, ish):
        start = bs.shell_to_basis_function(ish)
        nbf = bs.shell(ish).nfunction
        return slice(start, start + nbf)

    def _compute_ints(self, out, basis_sets, correlation_factor=None, clip_zero_basis=True):
        # Only create the zero_basis_set instance once we get here and need it
        if isinstance(self.zero_basis_set, lambda_type): self.zero_basis_set = self.zero_basis_set()
        #----------------------------------------#
        #region | Make the TwoBodyInt object |
        factory = psi4.IntegralFactory(*basis_sets)
        computer = getattr(factory, self.integral_type)
        if "f12" in self.integral_type:
            if correlation_factor is None:
                raise NotImplementedError()
            computer = computer(correlation_factor)
        elif "erf" in self.integral_type:
            raise NotImplementedError()
        else:
            computer = computer()
        #endregion
        #----------------------------------------#
        #region | Compute the integrals and store them in out.value |
        computer.set_enable_pybuffer(True)
        pybuffer = computer.py_buffer_object
        for shell_nums in xproduct(*[bs.nshell() for bs in basis_sets]):
            computer.compute_shell(*shell_nums)
            slices = tuple(
                AOTEIGetterWithComputer.shell_slice(bs, ish)
                    for bs, ish in zip(basis_sets, shell_nums) if bs is not self.zero_basis_set
            )
            buff = np.array(pybuffer)
            if clip_zero_basis:
                squeeze_axes =tuple(i for i, bs in enumerate(basis_sets) if bs is self.zero_basis_set)
                out.value[slices] = buff.squeeze(axis=squeeze_axes)
            else:
                out.value[slices] = buff
        #endregion
        #----------------------------------------#
        return

class DFThreeCenterAOTEIGetter(AOTEIGetterWithComputer):

    def __init__(self, two_center_bra, integral_type='eri'):
        self.two_center_bra = two_center_bra
        super(DFThreeCenterAOTEIGetter, self).__init__(
            integral_type=integral_type
        )

    # noinspection PyMethodOverriding
    def get_shape(self, basis, df_basis):
        nbf = basis.nbf()
        dfnbf = df_basis.nbf()
        if self.two_center_bra:
            return nbf, nbf, dfnbf
        else:
            return dfnbf, nbf, nbf

    # noinspection PyMethodOverriding
    def __call__(self, out, basis, df_basis, correlation_factor=None):
        if isinstance(self.zero_basis_set, lambda_type): self.zero_basis_set = self.zero_basis_set()
        if self.two_center_bra:
            basis_sets = (basis, basis, df_basis, self.zero_basis_set)
        else:
            basis_sets = (df_basis, self.zero_basis_set, basis, basis)
        self._compute_ints(out, basis_sets, correlation_factor)

class DFTwoCenterAOTEIGetter(AOTEIGetterWithComputer):

    def get_shape(self, df_basis):
        dfnbf = df_basis.nbf()
        return dfnbf, dfnbf

    def __call__(self, out, df_basis, correlation_factor=None):
        if isinstance(self.zero_basis_set, lambda_type): self.zero_basis_set = self.zero_basis_set()
        self._compute_ints(out,
            basis_sets=(df_basis, self.zero_basis_set, df_basis, self.zero_basis_set),
            correlation_factor=correlation_factor
        )

class ArbitraryBasisDatumGetter(DatumGetter):

    def __init__(self, set_names):
        self.basis_set_names = tuple(bsname.lower() for bsname in set_names)

    @abstractproperty
    def name(self):
        return NotImplemented

class ArbitraryBasisAOTEIGetter(ArbitraryBasisDatumGetter, AOTEIGetterWithComputer):

    def __init__(self, integral_type, bs1, bs2=None, bs3=None, bs4=None):
        if bs2 is None: bs2 = bs1
        if bs3 is None: bs3 = bs1
        if bs4 is None: bs4 = bs2
        super(ArbitraryBasisAOTEIGetter, self).__init__((bs1, bs2, bs3, bs4))
        AOTEIGetterWithComputer.__init__(self, integral_type)

    @property
    def name(self):
        return (
            "ao__"
            + "__".join(self.basis_set_names[:2])
            + "___" + self.integral_type + "___"
            + "__".join(self.basis_set_names[2:])
        )

    def get_shape(self, basis_sets):
        return tuple(bs.nbf() for bs in basis_sets)

    def __call__(self, out, basis_sets, correlation_factor=None):
        self._compute_ints(out, basis_sets, correlation_factor)

class GeneralCorrelationFactorAOTEIGetter(ArbitraryBasisAOTEIGetter):

    def __init__(self, coefficients, exponents, *args, **kwargs):
        super(GeneralCorrelationFactorAOTEIGetter, self).__init__(*args, **kwargs)
        self.coefficients = coefficients
        self.exponents = exponents

    # noinspection PyMethodOverriding
    def __call__(self, out, basis_sets):
        cf = psi4.CorrelationFactor(len(self.coefficients))
        coefs = psi4.Vector(len(self.coefficients))
        expons = psi4.Vector(len(self.exponents))
        for (ic, c), (ie, e) in zip(*map(enumerate, (self.coefficients, self.exponents))):
            coefs.set(0, ic, c)
            expons.set(0, ie, e)
        cf.set_params(coefs, expons)
        self._compute_ints(out, basis_sets, cf)

class AOOEIGetterWithComputer(MemmapArrayGetter):

    VALID_INTEGRAL_TYPES = [
        "overlap", "kinetic",
        "potential"
    ]

    def __init__(self, integral_type):
        self.integral_type = integral_type
        if self.integral_type not in type(self).VALID_INTEGRAL_TYPES:
            raise ValueError("Don't know how to compute integrals of type '{0}'".format(
                integral_type
            ))

    # TODO this is very similar to the _compute_ints in AOTEIGetterWithComputer and could be combined in a common superclass
    def _compute_ints(self, out, basis_sets):
        #----------------------------------------#
        #region | Make the TwoBodyInt object |
        factory = psi4.IntegralFactory(*basis_sets)
        computer = getattr(factory, "ao_" + self.integral_type)()
        #endregion
        #----------------------------------------#
        #region | Compute the integrals and store them in out.value |
        computer.set_enable_pybuffer(True)
        pybuffer = computer.py_buffer_object
        for shell_nums in xproduct(*[bs.nshell() for bs in basis_sets]):
            computer.compute_shell(*shell_nums)
            slices = tuple(
                AOTEIGetterWithComputer.shell_slice(bs, ish) for bs, ish in zip(basis_sets, shell_nums)
            )
            out.value[slices] = np.array(pybuffer)
        #endregion
        #----------------------------------------#
        return

    def get_shape(self, basis):
        nbf = basis.nbf()
        return nbf, nbf

    def __call__(self, out, basis, correlation_factor=None):
        self._compute_ints(out,
            basis_sets=(basis, basis)
        )

class ArbitraryBasisAOOEIGetter(ArbitraryBasisDatumGetter, AOOEIGetterWithComputer):

    def __init__(self, integral_type, bs1, bs2=None):
        if bs2 is None: bs2 = bs1
        super(ArbitraryBasisAOOEIGetter, self).__init__((bs1, bs2))
        AOOEIGetterWithComputer.__init__(self, integral_type)

    @property
    def name(self):
        return "ao__" + ("___" + self.integral_type + "___").join(self.basis_set_names)

    def get_shape(self, basis_sets):
        return tuple(bs.nbf() for bs in basis_sets)

    def __call__(self, out, basis_sets):
        self._compute_ints(out, basis_sets)

class MOEigenvaluesGetter(MemmapArrayGetter):

    #TODO open shell
    #TODO symmetry

    def get_shape(self, reference_wavefunction):
        return reference_wavefunction.nmo(),

    def __call__(self, out, reference_wavefunction):
        vect = reference_wavefunction.epsilon_a()
        nbf = out.value.shape[0]
        for ibf in range(nbf):
            out.value[ibf] = vect[0, ibf]

#endregion

#================================================================================#

class ComputationCache(object):
    """
    A set of computations cached in a given directory
    """

    #========================================#
    #region | Class Attributes |

    # Always add new attributes to the end of these lists,
    #   so that older files attributes can be seen as a
    #   subset of newer files attributes (though right now
    #   this ordering property isn't used)
    # Note that "molecule" is understood to be needed
    required_attributes = [
        'basis'
    ]

    available_psi_options = [
        "scf_convergence"

    ]

    # Sentinel value for optional attributes to not include in the
    #   key unless they are set by the user.
    NoDefaultValue = "___NoDefaultValue___"

    optional_attributes = AliasedDict({
        'correlation_factor_exponent' : 1.0,
        ('df_basis', 'df_basis_name', 'density_fitting_basis', 'fitting_basis', 'dfbasis') : NoDefaultValue
    })

    case_sensative_attributes = [
    ]

    do_writeback_on_sync = False

    #========================================#

    #region | Initialization |

    def __init__(self,
            directory,
            make_directory=False,
            make_path=False,
            make_shelve=True,
            writeback=True
    ):
        #========================================#
        #region Check for the directory and make it if needed.
        self.directory = str(directory).rstrip(path_sep)
        pth, dir_name = path_split(self.directory)
        if not isdir(pth):
            if exists(pth):
                raise OSError("File '{0}' exists, but it is not a directory.".format(pth))
            if make_path:
                makedirs(pth)
            else:
                raise OSError("Path '{0}' does not exist.".format(pth))
        if not isdir(self.directory):
            if make_directory:
                mkdir(self.directory)
            else:
                raise OSError("Directory {0} does not exist.".format(self.directory))
        #endregion
        #----------------------------------------#
        #region Open the Shelf
        flag = 'c' if make_shelve else 'w'
        self.shelf = shelve.open(
            path_join(self.directory, "computation_shelf.db"),
            flag=flag,
            protocol=2,
            writeback=writeback
        )
        self.analogous_keys = shelve.open(
            path_join(self.directory, "analogous_keys.db"),
            flag=flag,
            protocol=2,
            writeback=writeback
        )
        self.writeback = writeback
        def close_shelf(shelf, other_shelf):
            try:
                shelf.close()
            except ValueError:
                # If it's already closed, it's okay
                pass
            try:
                other_shelf.close()
            except ValueError:
                # If it's already closed, it's okay
                pass
        atexit.register(close_shelf, self.shelf, self.analogous_keys)
        #endregion
        #----------------------------------------#

    #endregion

    #========================================#

    #region | Methods |

    def get_computation(self,
            molecule,
            needed_data=None,
            **kwargs
    ):
        """

        @param molecule:
        @type molecule: Molecule
        @param needed_data:
        @type needed_data: NoneType, Iterable
        @rtype: CachedComputation
        """
        molkey = self.get_molecule_key(molecule)
        key = [molkey]
        key_mandatory = [molkey]
        optional_part = []
        init_kwargs = dict(molecule=molecule)
        psi_options = dict()
        #----------------------------------------#
        #region | Gather required attributes |
        for arg in self.required_attributes:
            has_arg = arg in kwargs
            val = kwargs.pop(arg, None)
            if not has_arg:
                # Cheesy case insensitivity
                found = False
                for kw in kwargs:
                    if kw.lower() == arg.lower():
                        val = kwargs.pop(kw)
                        found = True
                        break
                if not found:
                    raise ValueError("Missing required argument {0}".format(arg))
            if isinstance(val, str) and arg not in ComputationCache.case_sensative_attributes:
                val=val.lower()
            key.append((arg,val))
            key_mandatory.append((arg,val))
            init_kwargs[arg] = val
        #endregion
        #----------------------------------------#
        #region | Gather optional attributes |
        # This isn't quite perfect, since data that do not depend on
        #   a given optional attribute will be recomputed when really
        #   that data could be reused in a new computation which has
        #   some extra optional argument.  I don't see an easy way
        #   around this right now, though; remember that iterating
        #   over existing Computations should be out of the question.
        optional_args = dict()
        for keyset, defaultval in self.optional_attributes.items():
            kw = [kw for kw in kwargs if kw in keyset]
            kw = None if len(kw) == 0 else kw[-1]
            val = kwargs.pop(kw, None)
            has_arg = kw is not None
            if not has_arg:
                # Cheesy case insensitivity
                lower_keyset = set(k.lower() for k in keyset)
                for kw in kwargs:
                    if kw.lower() in lower_keyset:
                        val = kwargs.pop(kw)
                        has_arg = True
                        break
            firstkey = self.optional_attributes.firstkey(keyset)
            if isinstance(val, str) and firstkey not in ComputationCache.case_sensative_attributes:
                val=val.lower()
            if has_arg:
                key.append((firstkey,val))
                optional_part.append((firstkey,val))
                optional_args[firstkey] = val
            elif not defaultval == ComputationCache.NoDefaultValue:
                key.append((firstkey,val))
                optional_part.append((firstkey,val))
                optional_args[firstkey] = defaultval
        if len(kwargs) > 0:
            raise ValueError("Unknown or duplicate attribute {0}".format(kwargs.keys()[0]))
        #endregion
        #----------------------------------------#
        #region | Gather psi options |
        # These get set before each UnitializedDependency is evaluated and restored to their
        #   previous values after the given uninitialized dependency is evaluated.  Options
        #   that carry around baggage should not be used with this (for instance, 'basis'
        #   needs for the global environment to have a molecule defined or it raises an
        #   exception)
        for arg in self.available_psi_options:
            has_arg = arg in kwargs
            val = kwargs.pop(arg, None)
            if not has_arg:
                # Cheesy case insensitivity
                for kw in kwargs:
                    if kw.lower() == arg.lower():
                        val = kwargs.pop(kw)
                        has_arg = True
                        break
            if isinstance(val, str) and arg not in ComputationCache.case_sensative_attributes:
                val=val.lower()
            if has_arg:
                key.append((arg,val))
                key_mandatory.append((arg,val))
                psi_options[arg] = val
        #endregion
        #----------------------------------------#
        # This is not ideal.  For instance, two computations
        #   with scf_convergence set to 8 and 8.0 would show
        #   up as unique.  But shelve apperently requires a
        #   string to hash with.
        key_tuple = tuple(key)
        key = str(key_tuple)
        key_mandatory = str(tuple(key_mandatory))
        if key not in self.shelf:
            # The directory name isn't all that important,
            #   since it won't (and shouldn't) be human
            #   readable anyway.  The complicated hashing
            #   is just to get a unique key.
            # The abs() is because I don't like directory
            #   names that start with "-"
            hash_key = abs(hash(key))
            hash_dir = "{0:016x}".format(hash_key)[:16]
            while exists(path_join(self.directory, hash_dir)):
                # Just choose the next available value.  This
                #   should pretty much never happen, since
                #   an existing directory should also have
                #   an existing key in the shelf.  But just in
                #   case...
                hash_key += 1
                hash_dir = "{0:016x}".format(hash_key)[:16]
            hash_dir = path_join(self.directory, hash_dir)
            # Now make the CachedComputation object.  No need
            #   to sync afterwards, since the assignment triggers
            #   an automatic sync.
            rv = CachedComputation(
                needed_data=needed_data,
                directory=hash_dir,
                owner=self,
                optional_arguments=optional_args,
                **init_kwargs
            )
            rv.shelf_key = key
            rv.minimal_shelf_key = key_mandatory
            # Add the key to the list of keys related to the
            #   mandatory minimal subset of keys
            self._map_analogous_keys(key_mandatory, key_tuple)
            # Make sure the object is constructed successfully
            #   before shelving it
            self.shelf[key] = rv
        else:
            rv = self._get_existing_computation(key, key_mandatory)
            # Add the key to the list of keys related to the
            #   mandatory minimal subset of keys if it isn't
            #   already there
            self._map_analogous_keys(key_mandatory, key_tuple)
            if needed_data is not None:
                rv.get_data(needed_data)
        return rv

    def sync_computation(self, comp, sync_writeback=None):
        self.shelf[comp.shelf_key] = comp
        if ((sync_writeback is None and self.do_writeback_on_sync) or sync_writeback) and self.writeback:
            self.shelf.sync()

    def clear_datum(self, datum_name, fail_if_missing=False):
        """
        Clear datum named `datum_name` from all computations in the cache

        @param datum_name: The name of the datum to clear from all known computations
        @type datum_name: str
        @param fail_if_missing: Whether or not to raise an exception if no computations have a datum named `datum_name`
        @return: The number of data deleted
        @rtype: int
        @see CachedComputation.clear_datum, clear_data_regex
        """
        num_found = 0
        for key in self.shelf:
            comp = self._get_existing_computation(key)
            cleared = comp.clear_datum(datum_name)
            if cleared:
                num_found += 1
        if fail_if_missing and num_found == 0:
            raise ValueError("No known computations have a datum named '{0}'".format(datum_name))
        return num_found

    def clear_data_regex(self, regex, fail_if_missing=False, flags=0):
        """
        Clear all data matching regex from known computations

        @param regex:  The regular expression to search for in the name of the datum
        @type regex: str or compiled regex
        @param fail_if_missing: Whether or not to raise an exception if no data match the regex
        @param flags: Flags to compile the regular expression with.  Ignored if regex is an instance of re.RegexObject.
        @type flags: int
        @raise ValueError: if no data names match the regex and fail_if_missing is True
        @return: Length 2 tuple of the number of data removed and the number of computations it was removed from.
        @rtype: tuple
        @see CachedComputation.clear_data_regex
        """
        num_found = 0
        num_comps = 0
        if hasattr(regex, "pattern") and hasattr(regex, "search"):
            re_string = regex.pattern
        else:
            re_string = str(regex)
            regex = re.compile(re_string, flags)
        for key in self.shelf:
            comp = self._get_existing_computation(key)
            found = comp.clear_data_regex(regex, flags=flags)
            if found > 0:
                num_found += found
                num_comps += 1
        if fail_if_missing and num_found == 0:
            raise ValueError("No known computations have a datum matching '{0}'".format(re_string))
        return num_found, num_comps

    #endregion

    #========================================#

    #region | Private Methods |

    def _get_existing_computation(self, key, key_mandatory=None):
        """
        Internal use only

        @type key: str
        @type key_mandatory: str
        @rtype : CachedComputation
        @raise KeyError : if key is not found in the computation
        """
        rv = self.shelf[key]
        rv.shelf_key = key
        if key_mandatory is not None:
            rv.minimal_shelf_key = key_mandatory
        rv.owner = self
        return rv

    def _map_analogous_keys(self, key_mandatory, key_tuple):
        if key_mandatory in self.analogous_keys:
            # Don't append, since mutation of the object will
            #   not get immediately registered.  Instead, assign
            #   the appended list
            known_keys = self.analogous_keys[key_mandatory]
            if key_tuple not in known_keys:
                self.analogous_keys[key_mandatory] = self.analogous_keys[key_mandatory] + [key_tuple]
        else:
            self.analogous_keys[key_mandatory] = [key_tuple]

    #endregion

    #========================================#

    #region | Class Methods |

    @classmethod
    def get_molecule_key(cls, molecule):
        """
        For now, this just returns the xyz string representation of the
        Molecule object.  A more sophisticated means could be devised if
        constant recomputation becomes a problem, but it would be a total
        mess (see my attempts in the MoleculeStub class in grendel)
        """
        if isinstance(molecule, str):
            return str
        return molecule.xyz_string(header=False)

    #endregion

    #========================================#

    # End of ComputationCache class
    pass

class CachedComputation(object):
    """
    A computation on a given molecule with various relevant integrals
    and such cached for later as tensors or numpy memory maps.
    """

    #========================================#

    #region | Class Attributes |

    PICKLE_VERSION = (2,1,1)

    parser = psi4.Gaussian94BasisSetParser()

    array_getters = dict()
    # Add the AO and MO TEI getters to the known getters
    for kernel in TEIGetter.VALID_KERNELS:
        array_getters["ao_" + kernel] = AOTEIGetter(kernel)
        array_getters["mo_" + kernel] = MOTEIGetter(kernel)
    for kernel in AOOEIGetter.VALID_KERNELS:
        array_getters["ao_" + kernel] = AOOEIGetter(kernel)
    # Add the MO Coefficient getters
    for orbital_type in MOCoefficientsGetter.ALLOWED_TYPES:
        if orbital_type == "canonical":
            array_getters["mo_coefficients"] = MOCoefficientsGetter("canonical")
            array_getters["occupied_mo_coefficients"] = OccupiedMOCoefficientsGetter("canonical")
        else:
            array_getters["mo_coefficients" + "_" + orbital_type] = MOCoefficientsGetter(orbital_type)
            array_getters["occupied_mo_coefficients" + "_" + orbital_type] = OccupiedMOCoefficientsGetter(orbital_type)
    array_getters["mo_eigenvalues"] = MOEigenvaluesGetter()

    other_getters = dict()
    # Simple getters that amount to nothing more than
    #   a method call on a single dependency
    # SimpleMethodCallGetters for molecule
    _molecule_methods = ["nuclear_repulsion_energy", "natom", "name"]
    for method in _molecule_methods:
        other_getters["molecule_" + method] = SimpleMethodCallGetter('psi_molecule', method)
    # SimpleMethodCallGetters for basis
    _basis_methods = ["nbf", "nao", "nprimitive", "nshell", "has_puream"]
    for method in _basis_methods:
        other_getters["basis_" + method] = SimpleMethodCallGetter('basis', method)
    # SimpleMethodCallGetters for wavefunction
    _reference_wavefunction_methods = [
        "nso",  "nmo", "nirrep",
        "energy", "nalpha", "nbeta"
    ]
    for method in _reference_wavefunction_methods:
        other_getters["reference_wavefunction_" + method] = SimpleMethodCallGetter('reference_wavefunction', method)
    _reference_wavefunction_dimension_methods = [
        "doccpi", "soccpi", "nsopi",
        "nalphapi", "nbetapi", "frzcpi", "frzvpi"
    ]
    for method in _reference_wavefunction_dimension_methods:
        other_getters["reference_wavefunction_" + method] = DimensionListGetter('reference_wavefunction', method)

    default_psi_options = dict(
        scf_type = "direct"
    )

    #endregion

    #========================================#

    #region | Initialization |

    @typecheck(
        molecule=(Molecule, str),
        cached_data=(dict, None)
    )
    def __init__(self,
            molecule,
            basis,
            directory,
            optional_arguments,
            needed_data=None,
            cached_data=None,
            owner=None,
            psi_options=None,
    ):
        #========================================#
        #region Set up arguments and attributes
        if isinstance(molecule, str):
            self.molecule = Molecule(molecule)
        else:
            self.molecule = molecule
        self.basis_name = str(basis)
        if needed_data is None:
            needed_data = []
        elif isinstance(needed_data, str):
            needed_data = [needed_data]
        elif isinstance(needed_data, Iterable):
            needed_data = list(needed_data)
        self.directory = directory
        if not isdir(self.directory):
            if exists(self.directory):
                raise IOError("File {0} exists, but it is not a directory.".format(self.directory))
            else:
                mkdir(self.directory)
        self.owner = owner
        self.optional_arguments = ComputationCache.optional_attributes.copy()
        self.optional_arguments.update(optional_arguments)
        if any(k not in ComputationCache.optional_attributes for k in self.optional_arguments):
            raise ValueError("Invalid optional_arguments key: {0}".format(
                tuple(next(k for k in self.optional_arguments if k not in ComputationCache.optional_attributes))
            ))
        self._custom_getters = dict()
        self._basis_registry = AliasedDict()
        self._optional_argument_dependencies = dict()
        self._analogously_loaded_data = set()
        self.parallel_ready = False
        self._dependency_aliases = defaultdict(lambda: [])
        # multiprocessing RLock object for locking
        #   psi to avoid PSIO errors
        self.psi_lock = None
        self._required_locks = dict()
        self._original_names = dict()
        #endregion
        #========================================#
        #region psi options
        # Psi4 options get set before the evaluation
        #   of any lazy attributes and set back to there
        #   former values after that evaluation.  This
        #   prevents a change of option in one computation
        #   from affecting the option value in another
        #   computation which assumes the default.
        self._psi_options = CachedComputation.default_psi_options.copy()
        if psi_options is not None:
            self._psi_options.update(psi_options)
        #endregion
        #========================================#
        #region Set up the dictionary of data bits to cache
        if cached_data is None:
            self.cached_data = dict()
        else:
            self.cached_data = cached_data
        #endregion
        #========================================#
        #region Set up the psi molecule and psi objects
        # These objects don't get pickled
        #----------------------------------------#
        # Lazily evaluated attributes
        self.register_dependency(
            "psi_molecule",
            UninitializedDependency(
                # Don't use symmetry for now
                # Note that psi4.geometry calls activate on the molecule
                lambda self: psi4.geometry("symmetry c1\n" + self.molecule.xyz_string(header=False).strip()),
                self, depends_on=[]
            )
        )
        # basis, df_basis, and add the unit basis (a.k.a. zero basis) to the registry
        self.register_dependency(
            "basis",
            self.uninitialized_basis(self.basis_name)
        )
        self._basis_registry[self.basis_name.lower()] = self.basis
        self.register_dependency(
            "df_basis",
            UninitializedDependency(
                partial(
                    CachedComputation.construct_basis,
                    basis_name=self.optional_arguments['df_basis']
                ),
                self, depends_on=['psi_molecule'],
                requires_optional_arguments=['df_basis']
            ),
            dependent_optional_arguments=['df_basis']
        )
        if self.optional_argument_given("df_basis"):
            self._basis_registry[self.optional_arguments['df_basis'].lower()] = self.df_basis
        self._basis_registry[('unit','1','zero','one')] = UninitializedDependency(
            CachedComputation.construct_unit_basis,
            self, depends_on=['psi_molecule']
        )
        # other stuff
        self.register_dependency(
            "mints",
            UninitializedDependency(
                lambda self: psi4.MintsHelper(self.basis),
                self, depends_on=['basis']
            )
        )
        self.register_dependency(
            "factory",
            UninitializedDependency(
                lambda self: self.mints.integral(),
                self, depends_on=['mints']
            )
        )
        def get_scf_wavefunction(self):
            # Check the basis to be sure it's the right one
            if psi4.options.basis.lower() != self.basis_name.lower():
                psi4.options.basis = self.basis_name
            # run SCF
            scratch_path = tempfile.mkdtemp()
            psi4.IOManager.shared_object().set_default_path(scratch_path)
            psi4.energy("scf")
            psi4.clean()
            shutil.rmtree(scratch_path)
            # grab the wavefunction of the most recently run computation (scf, in this case)
            return psi4.wavefunction()
        self.register_dependency(
            "reference_wavefunction",
            UninitializedDependency(
                get_scf_wavefunction,
                self, depends_on=['basis']
            ),
            required_locks=["psi_lock"]
        )
        self.register_dependency(
            "correlation_factor",
            UninitializedDependency(
                lambda self: psi4.FittedSlaterCorrelationFactor(
                    self.optional_arguments['correlation_factor_exponent']
                ),
                self, depends_on=['psi_molecule']
            ),
            dependent_optional_arguments=['correlation_factor_exponent']
        )
        # TODO handle open shells
        self.register_dependency(
            "docc_space",
            UninitializedDependency(
                lambda self: self.reference_wavefunction.alpha_orbital_space('i', 'SO', 'OCC'),
                self, depends_on=['reference_wavefunction']
            ),
            aliases=["doubly_occupied_space", "occupied_space", "alpha_occupied_space"]
        )
        self.register_dependency(
            "virt_space",
            UninitializedDependency(
                lambda self: self.reference_wavefunction.alpha_orbital_space('a', 'SO', 'VIR'),
                self, depends_on=['reference_wavefunction']
            ),
            aliases=["virtual_space", "alpha_virtual_space"]
        )
        self.register_dependency(
            "mo_coefficients",
            UninitializedDependency(
                lambda self: self.reference_wavefunction.Ca(),
                self, depends_on=['reference_wavefunction']
            ),
            aliases=["mo_coefficients_canonical", "canonical_mo_coefficients"]
        )
        #endregion
        #========================================#
        #region Get the data we need
        while len(needed_data) != 0:
            needed_name = needed_data.pop(0)
            self._fill_datum(needed_name)
        #endregion
        #========================================#

    #endregion

    #========================================#

    #region | Methods |

    def optional_argument_given(self, arg):
        return self.optional_arguments[arg] != ComputationCache.NoDefaultValue

    def construct_basis(self, basis_name=None):
        if basis_name is None:
            basis_name = self.basis_name
        old_value = psi4.options.basis
        psi4.options.basis = basis_name
        rv = psi4.BasisSet.construct(CachedComputation.parser, self.psi_molecule, 'BASIS')
        if old_value != '':
            psi4.options.basis = old_value
        self._basis_registry[basis_name.lower()] = rv
        return rv

    def construct_unit_basis(self):
        rv = psi4.BasisSet.zero_ao_basis_set()
        self._basis_registry['unit'] = rv
        return rv

    def uninitialized_basis(self, basis_name):
        return UninitializedDependency(
            partial(CachedComputation.construct_basis, basis_name=basis_name),
            self, depends_on=['psi_molecule']
        )

    def register_dependency(
            self,
            attribute_name,
            dependency_object_or_function,
            depends_on=None,
            aliases=None,
            dependent_optional_arguments=None,
            required_locks=tuple()
    ):
        if isinstance(dependency_object_or_function, UninitializedDependency):
            dependency_object = dependency_object_or_function
            if depends_on is not None:
                dependency_object.depends_on.extend(depends_on)
        else:
            if depends_on is None: depends_on = []
            dependency_object = UninitializedDependency(
                dependency_object_or_function, self, depends_on
            )
        setattr(self, attribute_name, dependency_object)
        if aliases is not None:
            # Register pointers to the original object under alternate names
            for alias in aliases:
                self.register_dependency_alias(alias, attribute_name)
        if dependent_optional_arguments is not None:
            self._optional_argument_dependencies[attribute_name] = dependent_optional_arguments
        else:
            self._optional_argument_dependencies[attribute_name] = []
        self._required_locks[attribute_name] = list(required_locks)

    def register_dependency_alias(self, alias, original_name):
        attr = getattr(self, original_name)
        self._dependency_aliases[attr].append(alias)
        self._original_names[alias] = original_name
        setattr(self, alias, attr)

    def register_basis(self, basis_name):
        if basis_name not in self._basis_registry:
            self._basis_registry[basis_name.lower()] = self.uninitialized_basis(basis_name)

    def get_datum(self,
            name=None,
            custom_getter=None,
            pre_computation_callback=None,
            post_computation_callback=None,
            analogous_load_callback=lambda akey, fname=None: print(
                "Loaded datum from analogous source file {0} with computation attributes {1}".format(
                    "(no file)" if fname is None else fname,
                    akey[1:]
            )),
            allow_analogous_load=True
    ):
        #----------------------------------------#
        if name is None:
            if not hasattr(custom_getter, "name"):
                raise TypeError()
            name = custom_getter.name
        #----------------------------------------#
        if not self.has_datum(name):
            if (name in self._custom_getters
                and name in self.cached_data
                and not self.cached_data[name].filled
            ):
                del self._custom_getters[name]
            self._fill_datum(
                name, custom_getter,
                pre_computation_callback=pre_computation_callback,
                post_computation_callback=post_computation_callback,
                allow_analogous_load=allow_analogous_load,
                analogous_load_callback=analogous_load_callback
            )
            # Sync the parent shelf since self has been modified
            if self.owner is not None and not self.parallel_ready:
                self.owner.sync_computation(self)
        return self.cached_data[name]

    def get_data(self, names):
        rv = []
        for name in names:
            rv.append(self.get_datum(name))
        return rv

    def clear_datum(self, name, fail_if_missing=False):
        """
        Clear a datum from the computation.  If the datum has a file
        associated with it, delete the file.

        @param name: The name of the datum to clear
        @type name: str
        @param fail_if_missing: Whether or not to raise an exception if the datum is not found
        @type fail_if_missing: bool
        @raise ValueError: if the datum named `name` is not found
        @return: True if the datum was found and deleted, False otherwise
        @rtype: bool
        """
        if name in self.cached_data:
            if name not in self._analogously_loaded_data:
                d = self.cached_data[name]
                if isinstance(d, CachedMemmapArray):
                    fname = d.filename
                    del self.cached_data[name]
                    if exists(fname):
                        os.remove(fname)
                    else:
                        raise_warning(
                            "File '{0}' associated with datum named '{1}' is missing and"
                            " cannot be deleted".format(
                                fname, name
                            ),
                            FileMissingWarning
                        )
                else:
                    del self.cached_data[name]
            else:
                # Just remove the reference from the dictionary
                del self.cached_data[name]
                self._analogously_loaded_data.remove(name)
            # Sync the parent shelf since self has been modified
            if self.owner is not None and not self.parallel_ready:
                self.owner.sync_computation(self)
            return True
        elif fail_if_missing:
            raise ValueError("Datum '{0}' does not exist.  Known data names are:\n{1}".format(
                name, "    " + "\n    ".join(self.cached_data.keys())
            ))
        return False

    def clear_data_regex(self, regex, fail_if_missing=False, flags=0):
        """
        Clear all data matching regex

        @param regex:  The regular expression to search for in the name of the datum
        @type regex: str or re.RegexObject
        @param fail_if_missing: Whether or not to raise an exception if no data match the regex
        @param flags: Flags to compile the regular expression with.  Ignored if regex is an instance of re.RegexObject.
        @type flags: int
        @raise ValueError: if no data names match the regex and fail_if_missing is True
        @return: The number of data removed
        @rtype: int
        """
        found_count = 0
        if hasattr(regex, "pattern") and hasattr(regex, "search"):
            re_string = regex.pattern
        else:
            re_string = str(regex)
            regex = re.compile(re_string, flags)
        for name in list(self.cached_data.keys()):
            if regex.search(name) is not None:
                found_count += 1
                self.clear_datum(name)
        if fail_if_missing and found_count == 0:
            raise ValueError("No data matching '{0}' found.  Known data names are:\n{1}".format(
                re_string, "    " + "\n    ".join(self.cached_data.keys())
            ))
        if found_count > 0:
            if self.owner is not None and not self.parallel_ready:
                self.owner.sync_computation(self)
        return found_count

    def clear_all_data(self):
        for datum_name in list(self.cached_data.keys()):
            self.clear_datum(datum_name)

    def get_lazy_attribute(self, needed_attr):
        # Make sure we know how to get the attribute
        if not hasattr(self, needed_attr):
            raise NameError("Don't know how to get attribute '{0}'".format(needed_attr))
        #----------------------------------------#
        if needed_attr in self._original_names:
            needed_attr = self._original_names[needed_attr]
        acquired = []
        try:
            for lock_name in self._required_locks[needed_attr]:
                lock = getattr(self, lock_name)
                if lock is not None:
                    lock.acquire(blocking=True)
                    acquired.append(lock)
                elif self.parallel_ready:
                    raise ValueError("Required lock named '{}' not set for"
                                     " parallel ready CachedComputation".format(
                        lock_name
                    ))
            if self.parallel_ready:
                psi4.clean()
            self._psi_options_use()
            attr = getattr(self, needed_attr)
            if isinstance(attr, UninitializedDependency):
                aliases = self._dependency_aliases[attr]
                attr = attr.initialize()
                setattr(self, needed_attr, attr)
                for alias in aliases:
                    setattr(self, alias, attr)
            if self.parallel_ready:
                psi4.clean()
            self._psi_options_restore()
        finally:
            for lock in acquired:
                lock.release()
        return attr

    def get_file_path(self, datum_name):
        return path_join(self.directory, datum_name + ".npy")

    def store_datum(self, name, value):
        """

        @param name: str
        @param value: object
        @raise: NameError
        """
        if name in self.cached_data:
            if not isinstance(value, np.ndarray) and not self.cached_data[name].filled:
                if self.cached_data[name] == value:
                    # Forgot to set the filled flag before
                    self.cached_data[name].filled = True
                else:
                    raise NameError("Already have cached datum named '{0}'".format(name))
        elif name in self.array_getters or name in self.other_getters:
            raise NameError("Naming conflict:  Data storage operation requested for"
                            " datum named '{0}', but a getter for this datum is "
                            " already known in CachedComputation.{1}".format(
                name,
                "array_getters" if name in self.array_getters else 'other_getters'
            ))
        elif name in self._custom_getters:
            raise NameError("Naming conflict:  Custom getter for datum named '{0}'"
                            " already exists".format(
                name
            ))
        #----------------------------------------#
        if isinstance(value, np.ndarray):
            out = CachedMemmapArray(
                filename=self.get_file_path(name),
                shape=value.shape,
                # Force overwrite because if the file exists at this point,
                #   it's not attached to self.cached_data, so it's useless
                #   anyway.
                force_overwrite=True
            )
            out.value[...] = value
            out.filled = True
            self.cached_data[name] = out
        else:
            # Store it as a regular Datum
            out = CachedDatum(value)
            out.filled = True
            self.cached_data[name] = out
        #----------------------------------------#
        # Sync the parent shelf since self has been modified
        if self.owner is not None and not self.parallel_ready:
            self.owner.sync_computation(self)

    def has_datum(self, name):
        return (name in self.cached_data
            and isinstance(self.cached_data[name], CachedDatum)
            and self.cached_data[name].filled
        )

    def save_computation(self):
        if self.owner is None:
            raise ValueError("Can't save computation because parent cache doesn't exist")
        self.owner.sync_computation(self)
        self.owner.shelf.sync()

    def begin_parallel_use(self):
        self.parallel_ready = True

    def end_parallel_use(self):
        self.parallel_ready = False

    #endregion

    #========================================#

    #region | Special Methods |

    def __reduce_ex__(self, protocol):
        unloader = CachedComputationUnloader(self)
        return unloader, tuple()

    def __del__(self):
        if self.owner is not None:
            self.owner.sync_computation(self)

    #endregion

    #========================================#

    #region | Private Methods |

    def _psi_options_use(self):
        if not isinstance(self.psi_molecule, UninitializedDependency):
            psi4.activate(self.psi_molecule)
        self._psi_saved_opts = dict()
        for opt, val in self._psi_options.items():
            self._psi_saved_opts[opt] = getattr(psi4.options, opt)
            setattr(psi4.options, opt, val)

    def _psi_options_restore(self):
        for opt, val in self._psi_saved_opts.items():
            setattr(psi4.options, opt, val)

    def _fill_datum(self,
            needed_name=None,
            getter=None,
            pre_computation_callback=None,
            post_computation_callback=None,
            analogous_load_callback=None,
            allow_analogous_load=True
    ):
        if needed_name is None:
            if not hasattr(getter, 'name'):
                raise TypeError()
            needed_name = getter.name
        #----------------------------------------#
        if "ao_multipole_" in needed_name and needed_name not in self.array_getters:
            max_order = int(needed_name[13:])
            self.array_getters[needed_name] = AOMultipoleGetter(max_order)
        #----------------------------------------#
        if getter is not None:
            if needed_name in self.array_getters or needed_name in self.other_getters:
                raise NameError("Naming conflict:  Custom getter requested for datum named"
                                " '{0}', but a getter for this datum is already known in"
                                " CachedComputation.{1}".format(
                    needed_name,
                    "array_getters" if needed_name in self.array_getters else 'other_getters'
                ))
            elif needed_name in self._custom_getters:
                raise NameError("Naming conflict:  Custom getter for datum named '{0}'"
                                " already exists".format(
                    needed_name
                ))
            self._custom_getters[needed_name] = getter
        #========================================#
        # First check if we've already filled the given datum
        if needed_name in self.cached_data and self.cached_data[needed_name].filled:
            return
        #----------------------------------------#
        #region | Handle MemmapArrayGetters |
        if ((getter is not None and isinstance(getter, MemmapArrayGetter))
            or needed_name in self.array_getters
        ):
            fname = self.get_file_path(needed_name)
            if getter is None:
                getter = self.array_getters[needed_name]
            shape, _, dtype = CachedMemmapArray.read_file_metadata(fname)
            load_successful = False
            # Try to load the memmap if we successfully retrieved a shape
            #   from the cached file.
            has_tb = False
            if shape is not None:
                try:
                    # The shape and dtype arguments present redundant information,
                    #   (which is retrieved when the file is loaded anyway)
                    #   but including them here does no harm, and they may be used
                    #   for error checking or something else in the future.
                    out = CachedMemmapArray(fname,
                        shape=shape,
                        dtype=dtype
                    )
                    load_successful = True
                except ValueError:
                    # Raised if the data is invalid.  This means the file is corrupted
                    #   and needs to be rewritten
                    load_successful = False
                    out = None
                    has_tb = True
                except IOError:
                    # Raised if the file cannot be opened correctly.  This should never
                    #   really happen unless the file changes between the call to
                    #   CachedMemmapArray.read_file_metadata() and here.
                    load_successful = False
                    out = None
                    has_tb = True
            # Only compute the shape using the getter if we need to.
            if not load_successful:
                loaded_key = False
                # For now, if we are parallel ready, don't borrow data
                if allow_analogous_load and not self.parallel_ready:
                    loaded_key = self._check_analogous_computation(needed_name, getter)
                if loaded_key is not False:
                    # Another computatation has a compatible datum.  We can use it.
                    #   The helper function has already put it in our cached_data
                    #   dictionary, so all we need to do is use it.
                    out = self.cached_data[needed_name]
                    # We should check that there isn't a file in our own directory,
                    #   and if so, delete it.
                    if exists(fname):
                        raise_warning(
                            "File '{0}' is corrupted (couldn't get shape), but an analogous"
                            " datum was found in '{1}'; the corrupted file will be deleted"
                            " and the analogous one will be used.".format(
                                fname,
                                out.filename
                            ),
                            FileOverwriteWarning
                        )
                        os.remove(fname)
                    self._analogously_loaded_data.add(needed_name)
                    if callable(analogous_load_callback):
                        analogous_load_callback(loaded_key, out.filename)
                else:
                    # We don't have the file, or the file is invalid,
                    #   so we need to first get the shape and then
                    #   create the file.
                    shape_needs = getargspec(getter.get_shape).args[1:]  # exclude "self"
                    shape_kwargs = self._fill_needed_kwargs(shape_needs, getter)
                    shape = getter.get_shape(**shape_kwargs)
                    # If the file exists but we've gotten here, this means that the file
                    #   somehow corrupted and needs to be overwritten.  Warn the user
                    #   of the situation.
                    if exists(fname):
                        if has_tb:
                            raise_warning(
                                "File '{0}' is corrupted and will be overwritten."
                                " Its data will be recomputed.  Traceback:\n    {1}".format(
                                    fname,
                                    "\n    ".join(traceback.format_exc())),
                                FileOverwriteWarning
                            )
                        else:
                            raise_warning(
                                "File '{0}' is corrupted (couldn't get shape) and will be overwritten."
                                " Its data will be recomputed.".format(
                                    fname),
                                FileOverwriteWarning
                            )
                    # Construct the CachedMemmapArray with force_overwrite=True since
                    #   we already checked for a valid file earlier and found that
                    #   we couldn't determine the shape.
                    out = CachedMemmapArray(
                        fname,
                        shape=shape,
                        force_overwrite=True
                    )
            #----------------------------------------#
            # Only fill it if we need to
            if not out.filled:
                if out.read_only:
                    # This shouldn't happen
                    raise_warning(
                        "File '{0}' exists, but was not marked as filled previously.  It's data"
                        " will be recomputed.".format(
                            fname
                        ),
                        FileOverwriteWarning
                    )
                    out = CachedMemmapArray(
                        fname,
                        shape=shape,
                        force_overwrite=True
                    )
                if pre_computation_callback is not None:
                    pre_computation_callback()
                getter_kwargs = self._fill_needed_kwargs(getter.needs, getter)
                # This is where the computation of the data actually happens
                try:
                    getter(out, **getter_kwargs)
                except Exception:
                    # Remove the file caching the failure
                    if exists(fname):
                        remove(fname)
                    raise
                else:
                    out.filled = True
                if post_computation_callback is not None:
                    post_computation_callback()
            self.cached_data[needed_name] = out
        #endregion
        #----------------------------------------#
        #region | Handle simpler getters |
        elif isinstance(getter, DatumGetter) or needed_name in CachedComputation.other_getters:
            if not isinstance(getter, DatumGetter):
                getter = CachedComputation.other_getters[needed_name]
            loaded_key = False
            # For now, if we are parallel ready, don't borrow data
            if allow_analogous_load and not self.parallel_ready:
                loaded_key = self._check_analogous_computation(needed_name, getter)
            if loaded_key is not False:
                self._analogously_loaded_data.add(needed_name)
                if callable(analogous_load_callback):
                    analogous_load_callback(loaded_key)
            else:
                if pre_computation_callback is not None:
                    pre_computation_callback()
                getter_kwargs = self._fill_needed_kwargs(getter.needs, getter)
                self.cached_data[needed_name] = CachedDatum(getter(**getter_kwargs))
                self.cached_data[needed_name].filled = True
                if post_computation_callback is not None:
                    post_computation_callback()
        #endregion
        #----------------------------------------#
        else:
            raise NotImplementedError("Don't know how to get '{}'".format(needed_name))

    def _check_analogous_computation(self, needed_datum_name, getter):
        if self.owner is not None and hasattr(self, "minimal_shelf_key"):
            analogs = self.owner.analogous_keys[self.minimal_shelf_key]
            needed_opt = []
            # remove the basis_sets attribute which is used by arbitrary basis getters
            #   and does not serve to make one computation distinct from another anyway
            needs_no_bs = [g for g in getter.needs if g != "basis_sets"]
            for needed_attr in needs_no_bs:
                needed_opt.extend(self._optional_argument_dependencies[needed_attr])
            for akey in analogs:
                if str(akey) == self.shelf_key:
                    # Don't check ourselves...
                    continue
                useable = True
                for opt in needed_opt:
                    # Skip the molecule key using akey[1:]; it's never optional
                    avalues = [v for k, v in akey[1:] if k == opt]
                    if len(avalues) == 0:
                        # A required optional argument was not given
                        #   for the analogou,s computation, so we
                        #   can't use it
                        useable = False
                        break
                    elif len(avalues) == 2:
                        # This should never happen
                        raise KeyError("Something went horribly wrong; multiply defined key")
                    else: # len(avalues) == 1
                        if avalues[0] != self.optional_arguments[opt]:
                            # It doesn't have the same value, so we
                            #   can't use it
                            useable = False
                            break
                if useable:
                    acomp = self.owner._get_existing_computation(str(akey), self.minimal_shelf_key)
                    if acomp.has_datum(needed_datum_name):
                        self.cached_data[needed_datum_name] = acomp.get_datum(needed_datum_name)
                        return akey
        return False



    def _fill_needed_kwargs(self, needed_kws, getter):
        rv = dict()
        for needed_attr in needed_kws:
            if isinstance(getter, ArbitraryBasisDatumGetter) and needed_attr == "basis_sets":
                sets = []
                for bsname in getter.basis_set_names:
                    if bsname.lower() not in self._basis_registry:
                        raise NameError("Unknown basis set '{0}'".format(bsname))
                    bs = self._basis_registry[bsname.lower()]
                    if isinstance(bs, UninitializedDependency):
                        bs = bs.initialize()
                        if bsname == self.basis_name:
                            self.basis = bs
                        elif self.optional_argument_given('df_basis') and bsname == self.optional_arguments['df_basis']:
                            self.df_basis = bs
                    sets.append(bs)
                rv[needed_attr] = tuple(sets)
            else:
                rv[needed_attr] = self.get_lazy_attribute(needed_attr)
        return rv

    #endregion

    #========================================#
    # end CachedComputation class
    pass

#================================================================================#

#region | Cached data classes and unloaders |

class CachedDatum(object):
    """
    Anything from a computation that might need to be cached.
    Assuming value is picklable, the pickling functionality of this
    class is trivial.  Subclasses may override for more complex
    unpickling behavior
    """
    def __init__(self, value):
        # Remember, this only gets called when constructed normally, not when unpickled
        self.value = value

class CachedMemmapArray(CachedDatum):
    """
    An instance of a numpy.ndarray subclass to be cached.  Only the filename is pickled
    """

    ####################
    # Class Attributes #
    ####################

    DEFAULT_MMAP_LOAD_MODE = 'r'
    files_needing_flush = set()

    ##################
    # Initialization #
    ##################

    # noinspection PyMissingConstructor
    def __init__(self,
            filename,
            dtype="float64",
            shape=None,
            force_overwrite=False,
            mmap_load_mode=DEFAULT_MMAP_LOAD_MODE
    ):
        self.filename = filename
        self.dtype = dtype
        # Note:  mmap_load_mode is the mode used when loading.  When creating,
        #   "w+" is always used.
        self.mmap_load_mode=mmap_load_mode
        self.shape = shape
        self.write_mode = False
        self.force_overwrite = force_overwrite
        self.read_only = not force_overwrite and exists(self.filename)
        self.filled = False

    ###################
    # Special Methods #
    ###################

    def __getstate__(self):
        state = self.__dict__.copy()
        # Just get rid of the value, leave everything else
        if 'value' in state:
            # remember, value only gets filled if it is accessed
            del state['value']
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        # Note that the 'filled' attribute will stick with the class
        #   from the pickled instance
        self.force_overwrite = False

    #################
    # Class Methods #
    #################

    @classmethod
    def read_file_metadata(cls, filename):
        """
        Try to get the shape, ordering, and dtype of the memmap cached in the
        file `filename`.  If the file doesn't exist, isn't readable, or the
        data is invalid, returns (None, None, None).

        See numpy.lib.format.read_array_header_1_0()
        """
        try:
            fp = open(filename)
            major, minor = read_magic(fp)
            shape, order, dtype = read_array_header_1_0(fp)
            fp.close()
            return shape, order, dtype
        except ValueError:
            return None, None, None
        except IOError:
            return None, None, None
        except TypeError:
            return None, None, None

    ##############
    # Properties #
    ##############

    @LazyProperty
    def value(self):
        return self._load_value(self.force_overwrite)

    ###########
    # Methods #
    ###########

    def flush_value(self):
        if self.filled:
            self.value.flush()

    ###################
    # Private Methods #
    ###################

    def _load_value(self, force_overwrite=False):
        load_successful = False
        if exists(self.filename) and not force_overwrite:
            try:
                value = np.load(
                    self.filename,
                    mmap_mode=self.mmap_load_mode
                )
                if self.shape is not None and self.shape != value.shape:
                    raise_warning("Loaded shape {0} does not match specified (or saved) shape {1}".format(
                        value.shape, self.shape
                    ), ShapeMismatchWarning)
                self.shape = value.shape
                self.filled = True
                self.write_mode = False
                load_successful = True
            except ValueError as e:
                raise_warning( "File '{0}' is corrupted and will be overwritten."
                      " Its data will be recomputed.  The following error was raised:\n{1}".format(
                    self.filename,
                    traceback.format_exc()),
                      FileOverwriteWarning
                )
                load_successful = False
                unlink(self.filename)
        if not load_successful:
            if self.shape is None:
                raise ValueError("Can't create cached array if I don't know the shape.")
            value = open_memmap(
                self.filename,
                shape=self.shape,
                mode='w+',
                dtype=self.dtype
            )
            self.write_mode = True
            # Make sure the value is flushed on exit
            if self.filename not in CachedMemmapArray.files_needing_flush:
                atexit.register(self.flush_value)
                CachedMemmapArray.files_needing_flush.add(self.filename)
            self.filled = False
        return value

class CachedComputationUnloader(object):
    """
    Callable, picklable class that CachedComputation.__reduce_ex__
    returns an instance of.  Pickling and unpickling proceeds by
    automatic mechanisms
    """
    def __init__(self, comp):
        """
        Argument `pickler_version` is the revision of the file pickling procedure
        used when the CachedComputation instance was written.
        """
        self.pickler_version = CachedComputation.PICKLE_VERSION
        self.molecule_xyz = comp.molecule.xyz_string(header=False)
        self.optional_arguments = comp.optional_arguments.get_simple_dict()
        self.init_kwargs = dict(
            basis=comp.basis_name,
            cached_data=comp.cached_data,
            psi_options=comp._psi_options,
            directory=comp.directory
        )
        self.basis_sets_to_register = tuple(comp._basis_registry.keys())
        if hasattr(comp, "minimal_shelf_key"):
            self.minimal_shelf_key = comp.minimal_shelf_key
        else:
            self.minimal_shelf_key = None
        if hasattr(comp, "shelf_key"):
            self.shelf_key = comp.shelf_key
        else:
            self.shelf_key = None

    def __call__(self):
        init_kwargs = dict(
            # Need to recreate the molecule from the xyz string
            molecule=self.molecule_xyz,
            optional_arguments=self.optional_arguments,
            **self.init_kwargs
        )
        #========================================#
        rv = CachedComputation(**init_kwargs)
        #========================================#
        if self.pickler_version >= (2, 0, 0):
            for bsname in self.basis_sets_to_register:
                if self.pickler_version >= (2, 0, 1):
                    # changed _basis_registry from a dict to an AliasedDict,
                    #   thus the keys will be frozenset objects
                    bsname = next(iter(bsname))
                if bsname not in rv._basis_registry:
                    rv.register_basis(bsname)
            #----------------------------------------#
            if self.pickler_version >= (2, 1, 0):
                if self.minimal_shelf_key is not None:
                    rv.minimal_shelf_key = self.minimal_shelf_key
                if self.shelf_key is not None:
                    rv.shelf_key = self.shelf_key
        #========================================#
        return rv

#endregion
