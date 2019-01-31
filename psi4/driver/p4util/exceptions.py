#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2019 The Psi4 Developers.
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


class ConvergenceError(PsiException):
    """Error called for problems with converging and iterative method.

    Parameters
    ----------
    eqn_description : str
        Type of QC routine that has failed (e.g., SCF)
    iteration : int
        What iteration we failed on

    """

    def __init__(self, eqn_description, iteration):
        msg = "Could not converge %s in %d iterations." % (eqn_description, iteration)
        PsiException.__init__(self, msg)
        self.iteration = iteration
        self.message = msg
        core.print_out('\nPsiException: %s\n\n' % (msg))


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


class CSXError(PsiException):
    """Error called when CSX generation fails.

    """

    def __init__(self, msg):
        PsiException.__init__(self, msg)
        self.message = '\nCSXException: %s\n\n' % msg


class ManagedMethodError(PsiException):
    def __init__(self, circs):
        if circs[5] == '':
            msg = """{0}: Method '{1}' with {2} '{3}' and REFERENCE '{4}' not available{5}""".format(*circs)
        else:
            msg = """{0}: Method '{1}' with {2} '{3}' and REFERENCE '{4}' not directable to QC_MODULE '{5}'""".format(
                *circs)
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
