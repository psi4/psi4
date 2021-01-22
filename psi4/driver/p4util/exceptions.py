#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2021 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#
"""Module with non-generic exceptions classes."""

from psi4 import core
from psi4 import extras


class PsiException(Exception):
    """Error class for Psi."""
    extras._success_flag_ = False
    pass


class ValidationError(PsiException):
    """Error called for problems with the input file. Prints
    error message *msg* to standard output stream and output file.

    """

    def __init__(self, msg):
        PsiException.__init__(self, msg)
        # print("%s" % repr(msg))
        self.message = '\nPsiException: %s\n\n' % repr(msg)


class ParsingError(PsiException):
    """Error called for problems parsing a text file. Prints error message
    *msg* to standard output stream and output file.

    """

    def __init__(self, msg):
        PsiException.__init__(self, msg)
        self.message = '\nPsiException: %s\n\n' % msg


class PsiImportError(PsiException):
    """Error called for problems import python dependencies. Prints error message
    *msg* to standard output stream and output file.
    """

    def __init__(self, msg):
        PsiException.__init__(self, msg)
        self.message = '\nPsiException: %s\n\n' % msg


class TestComparisonError(PsiException):
    """Error called when a test case fails due to a failed
    compare_values() call. Prints error message *msg* to standard
    output stream and output file.

    """

    def __init__(self, msg):
        PsiException.__init__(self, msg)
        self.message = '\nPsiException: %s\n\n' % msg


class UpgradeHelper(PsiException):
    """Error called on previously valid syntax that now isn't and a
    simple syntax transition is possible.

    It is much preferred to leave the old syntax valid for a release
    cycle and have the old syntax raise a deprecation FutureWarning. For
    cases where the syntax just has to jump, this can be used to trap
    the old syntax at first error and suggest the new.

    """

    def __init__(self, old, new, version, elaboration):
        msg = "Using `{}` instead of `{}` is obsolete as of {}.{}".format(old, new, version, elaboration)
        PsiException.__init__(self, msg)
        core.print_out('\nPsiException: %s\n\n' % (msg))


class ConvergenceError(PsiException):
    """Error called for problems with converging an iterative method.

    Parameters
    ----------
    eqn_description : str
        Type of QC routine that has failed (e.g., SCF)
    iteration : int
        What iteration we failed on

    """

    def __init__(self, eqn_description, iteration, additional_info=None):
        msg = f"Could not converge {eqn_description:s} in {iteration:d} iterations."
        if additional_info is not None:
            msg += f"\n\n{additional_info}"
        PsiException.__init__(self, msg)
        self.iteration = iteration
        self.message = msg
        core.print_out(f'\nPsiException: {msg:s}\n\n')


class OptimizationConvergenceError(ConvergenceError):
    """Error called for problems with geometry optimizer."""

    def __init__(self, eqn_description, iteration, wfn):
        ConvergenceError.__init__(self, eqn_description, iteration)
        self.wfn = wfn


class SCFConvergenceError(ConvergenceError):
    """Error called for problems with SCF iterations.

    Parameters
    ----------
    wfn : psi4.core.Wavefunction
        Wavefunction at time of exception
    e_conv : float
        Change in energy for last iteration
    d_conv : float
        RMS change in density for last iteration

    """

    def __init__(self, eqn_description, iteration, wfn, e_conv, d_conv):
        ConvergenceError.__init__(self, eqn_description, iteration)
        self.e_conv = e_conv
        self.d_conv = d_conv
        self.wfn = wfn


class TDSCFConvergenceError(ConvergenceError):
    """Error called for problems with TDSCF iterations.

    Parameters
    ----------
    wfn : psi4.core.Wavefunction
        Wavefunction at time of exception
    what : str
        What we were trying to solve for (singlets/triplets, irrep) when we failed to converge
    stats : Dict
        Dictionary of convergence statistics of last iteration.
        Keys are:

          count : int, iteration number
          res_norm : np.ndarray (nroots, ), the norm of residual vector for each roots
          val : np.ndarray (nroots, ), the eigenvalue corresponding to each root
          delta_val : np.ndarray (nroots, ), the change in eigenvalue from the last iteration to this ones
          collapse : bool, if a subspace collapse was performed
          product_count : int, the running total of product evaluations that was performed
          done : bool, if all roots were converged
    """

    def __init__(self, iteration, wfn, what, stats):
        # prepare message, including excitation energies and residual norm
        conv_info = "==> Convergence statistics from last iteration <==\n\n"
        conv_info += "Excitation Energy".center(21) + f" {'D[value]':^15}" + "|R|".center(11) + "\n"
        conv_info += f"{'-':->20} {'-':->15} {'-':->15}\n"
        for e, diff, r_norm in zip(stats["val"], stats["delta_val"], stats["res_norm"]):
            conv_info += f"      {e:.6f}         {diff:-11.5e}    {r_norm:12.5e}\n"
        ConvergenceError.__init__(self,
                                  eqn_description=f"""TDSCF solver ({what})""",
                                  iteration=iteration,
                                  additional_info=conv_info)
        self.wfn = wfn
        self.what = what
        self.stats = stats


class CSXError(PsiException):
    """Error called when CSX generation fails.

    """

    def __init__(self, msg):
        PsiException.__init__(self, msg)
        self.message = '\nCSXException: %s\n\n' % msg


class MissingMethodError(ValidationError):
    """Error called when method not available.

    """

    def __init__(self, msg):
        ValidationError.__init__(self, msg)
        self.message = '\nMissingMethodError: %s\n\n' % msg


class ManagedMethodError(PsiException):
    def __init__(self, circs):
        if circs[5] == '':
            modulemsg = "not available"
        else:
            modulemsg = f"not directable to QC_MODULE '{circs[5]}'"

        if len(circs) == 7:
            msg = f"""{circs[0]}: Method '{circs[1]}' with {circs[2]} '{circs[3]}', FREEZE_CORE '{not circs[6]}', and REFERENCE '{circs[4]}' {modulemsg}"""
        else:
            msg = f"""{circs[0]}: Method '{circs[1]}' with {circs[2]} '{circs[3]}' and REFERENCE '{circs[4]}' {modulemsg}"""
        PsiException.__init__(self, msg)
        self.message = '\nPsiException: %s\n\n' % msg


class Dftd3Error(PsiException):
    """

    """

    def __init__(self, msg):
        PsiException.__init__(self, msg)
        self.message = '\nDftd3Error: %s\n\n' % msg


class PastureRequiredError(PsiException):
    """Error called when the specified value of *option* requires some
    module(s) from Psi4Pasture, but could not be imported.
    """
    msg_tmpl = """Psi4Pasture module(s) [{modlist}] are required to change the default value of {opt}

    """
    install_instructions = """
    Note: Psi4Pasture is currently in an experimental state with no reliable install
    procedure yet, but this is what it would look like.

    To Build Psi4Pasture and install the required modules within your current
    Psi4 installation

    >>> # clone the pasture repo
    >>> git clone https://github.com/psi4/psi4pasture.git

    >>> cmake -H. -Bobjdir -Dpsi4_DIR=$PSI4_INSTALL_PREFIX/share/cmake/psi4 {module_args}
    >>> # $PSI4_INSTALL_PREFIX is the $CMAKE_INSTALL_PREFIX for the psi4
    >>> # install you want to install pasture to

    >>> # build + install install location is detected automatically
    >>> cd objdir
    >>> make && make install

    See https://github.com/psi4/psi4pasture for more details

    Or to install using psi4's own build system add
         {module_args}
    to cmake command line when building psi4.
    """
    pasture_required_modules = {"RUN_CCTRANSORT": ["ccsort", "transqt2"]}

    def __init__(self, option):
        mods_str = ", ".join([m for m in PastureRequiredError.pasture_required_modules[option]])
        msg = PastureRequiredError.msg_tmpl.format(opt=option, modlist=mods_str)
        PsiException.__init__(self, msg)
        module_cmake_args = " ".join(
            ["-DENABLE_{}=ON".format(module) for module in PastureRequiredError.pasture_required_modules[option]])
        msg += PastureRequiredError.install_instructions.format(module_args=module_cmake_args)
        self.message = '\nPsiException: {}\n\n'.format(msg)
        core.print_out(self.message)
