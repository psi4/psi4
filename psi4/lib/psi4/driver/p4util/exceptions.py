#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2024 The Psi4 Developers.
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

__all__ = [
    "ConvergenceError",
    "MissingMethodError",
    "ManagedMethodError",
    "OptimizationConvergenceError",
    "ParsingError",
    "PastureRequiredError",
    "PsiException",
    "SCFConvergenceError",
    "TDSCFConvergenceError",
    "TestComparisonError",
    "UpgradeHelper",
    "ValidationError",
    "docs_table_link",
    "sanitize_method",
]

from typing import Any, Dict, List, Optional

from psi4 import core, extras


class PsiException(Exception):
    """Error class for |PSIfour|. Flags success as False (triggering coffee)."""
    extras._success_flag_ = False


class ValidationError(PsiException):
    """Input specification has problems.

    Error message *msg* directed both to standard output stream and to outfile.

    Parameters
    ----------
    msg
        Human readable string describing the exception.

    Attributes
    ----------
    message
        Human readable string describing the exception.

    """
    message: str

    def __init__(self, msg: str):
        PsiException.__init__(self, msg)
        self.message = '\nPsiException: %s\n\n' % repr(msg)


class ParsingError(PsiException):
    """Error called for problems parsing a text file. Prints error message
    *msg* to standard output stream and output file.

    Only used by untested distributed CC response machinery.

    """

    def __init__(self, msg):
        PsiException.__init__(self, msg)
        self.message = '\nPsiException: %s\n\n' % msg


# PsiImportError ceased to be used by v1.1. Class removed by v1.7
# class PsiImportError(PsiException):


class TestComparisonError(PsiException):
    """Error called when a :func:`~psi4.compare_values` or other comparison
    function fails.

    Error message *msg* directed both to standard output stream and to outfile.

    Parameters
    ----------
    msg
        Human readable string describing the exception.

    Attributes
    ----------
    message
        Human readable string describing the exception.

    Example
    -------
    >>> psi4.compare_values(2, 3, 2, "asdf")
    asdf..................................................................................FAILED
    psi4.driver.p4util.exceptions.TestComparisonError:  asdf: computed value (3.0000) does not match (2.0000) to atol=0.01 by difference (1.0000).
    !----------------------------------------------------------------------------------!
    !                                                                                  !
    !         asdf: computed value (3.0000) does not match (2.0000) to atol=0.01 by    !
    !     difference (1.0000).                                                         !
    !                                                                                  !
    !----------------------------------------------------------------------------------!

    """
    message: str

    def __init__(self, msg: str):
        PsiException.__init__(self, msg)
        self.message = '\nPsiException: %s\n\n' % msg


class UpgradeHelper(PsiException):
    """Error called on previously valid syntax that now isn't and a
    simple syntax transition is possible.

    It is much preferred to leave the old syntax valid for a release
    cycle and have the old syntax raise a deprecation :class:`FutureWarning`.
    For cases where the syntax just has to jump, an UpgradeHelper can be used
    to trap the old syntax at first error and suggest the new.

    An UpgradeHelper can also be used after the :class:`FutureWarning`
    described above has expired. Then the body of the code can be deleted while
    the definition is preserved, and an UpgradeHelper called in place of the
    body to guide users with lagging syntax.

    Parameters
    ----------
    old
        Previously valid syntax.
    new
        Suggested replacement syntax.
    version
        First Major.minor version at which `old` syntax won't run. Generally
        the next release at time of commit.
    elaboration
        Any additional message to convey. Should start with a space.

    """
    def __init__(self, old: str, new: str, version: str, elaboration: str):
        msg = "Using `{}` instead of `{}` is obsolete as of {}.{}".format(old, new, version, elaboration)
        PsiException.__init__(self, msg)
        core.print_out('\nPsiException: %s\n\n' % (msg))


class ConvergenceError(PsiException):
    """Error called for problems with converging an iterative method.

    Parameters
    ----------
    eqn_description
        Type of QC routine that has failed (e.g., SCF, optimization).
    iteration
        Iteration number on which routine failed.
    additional_info
        Any additional message to convey.

    Attributes
    ----------
    message
        Human readable string describing the exception.
    iteration
        Iteration number on which routine failed.

    """
    message: str
    iteration: int

    def __init__(self, eqn_description: str, iteration: int, additional_info: Optional[str] = None):
        msg = f"Could not converge {eqn_description:s} in {iteration:d} iterations."
        if additional_info is not None:
            msg += f"\n\n{additional_info}"
        PsiException.__init__(self, msg)
        self.iteration = iteration
        self.message = msg
        core.print_out(f'\nPsiException: {msg:s}\n\n')


class OptimizationConvergenceError(ConvergenceError):
    """Error called for problems with geometry optimizer.

    Parameters
    ----------
    eqn_description
        Type of QC routine that has failed (e.g., geometry optimization).
    iteration
        Iteration number on which routine failed.
    wfn
        Wavefunction at time of exception.

    Attributes
    ----------
    message
        Human readable string describing the exception.
    iteration
        Iteration number on which routine failed.
    wfn
        Wavefunction at time of exception.

    """
    message: str
    iteration: int
    wfn: core.Wavefunction

    def __init__(self, eqn_description: str, iteration: int, wfn: core.Wavefunction):
        ConvergenceError.__init__(self, eqn_description, iteration)
        self.wfn = wfn


class SCFConvergenceError(ConvergenceError):
    """Error called for problems with SCF iterations.

    Parameters
    ----------
    eqn_description
        Type of QC routine that has failed (e.g., SCF preiterations).
    iteration
        Iteration number on which routine failed.
    wfn
        Wavefunction at time of exception.
    e_conv
        Change in energy for last iteration.
    d_conv
        RMS change in density for last iteration.

    Attributes
    ----------
    message
        Human readable string describing the exception.
    iteration
        Iteration number on which routine failed.
    wfn
        Wavefunction at time of exception.
    e_conv
        Change in energy for last iteration.
    d_conv
        RMS change in density for last iteration.

    """
    message: str
    iteration: int
    wfn: core.Wavefunction
    e_conv: float
    d_conv: float

    def __init__(
        self,
        eqn_description: str,
        iteration: int,
        wfn: core.Wavefunction,
        e_conv: float,
        d_conv: float,
    ):
        ConvergenceError.__init__(self, eqn_description, iteration)
        self.e_conv = e_conv
        self.d_conv = d_conv
        self.wfn = wfn


class TDSCFConvergenceError(ConvergenceError):
    """Error called for problems with TDSCF iterations.

    Parameters
    ----------
    wfn
        Wavefunction at time of exception
    what
        What we were trying to solve for (singlets/triplets, irrep) when we failed to converge
    stats
        Dictionary of convergence statistics of last iteration.
        Keys are:

          - count : int, iteration number
          - res_norm : np.ndarray (nroots, ), the norm of residual vector for each roots
          - val : np.ndarray (nroots, ), the eigenvalue corresponding to each root
          - delta_val : np.ndarray (nroots, ), the change in eigenvalue from the last iteration to this ones
          - collapse : bool, if a subspace collapse was performed
          - product_count : int, the running total of product evaluations that was performed
          - done : bool, if all roots were converged

    Attributes
    ----------
    message
        Human readable string describing the exception.
    iteration
        Iteration number on which routine failed.
    wfn
        Wavefunction at time of exception
    what
        What we were trying to solve for (singlets/triplets, irrep) when we
        failed to converge.
    stats
        Dictionary of convergence statistics of last iteration. See keys above.

    """
    message: str
    iteration: int
    wfn: core.Wavefunction
    what: str
    stats: Dict[str, Any]

    def __init__(self, iteration: int, wfn: core.Wavefunction, what: str, stats: Dict[str, Any]):
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


# CSXError ceased to be used by v1.4. Class removed by v1.7
# class CSXError(PsiException):


class MissingMethodError(ValidationError):
    """Error called when requested level or theory or derivative level are not
    available.

    Parameters
    ----------
    msg
        Human readable string describing the exception.

    Attributes
    ----------
    message
        Human readable string describing the exception.

    """
    message: str

    def __init__(self, msg: str):
        ValidationError.__init__(self, msg)
        self.message = '\nMissingMethodError: %s\n\n' % msg


class ManagedMethodError(PsiException):
    """Error called when a requested level of theory and derivative level are
    nominally available but not for the particular conditions (e.g., reference,
    algorithm, active orbitals, QC module, etc.) requested.

    Parameters
    ----------
    circs
        List providing calling function name, level of theory, algorithm,
        reference, QC module, and frozen-core/all-electron requested conditions.

    Attributes
    ----------
    message
        Human readable string describing the exception.
    stats
        Dictionary of conditions for which method request was made.
        Keys are:

          - driver : {"energy", "gradient", "hessian", "properties"}
          - derivative_int : {0, 1, 2, "prop"}, integer representation of driver
          - method : str, model/method/level-of-theory
          - link : str, link address to table in docs
          - method_type : (str, str), raw value and str repr of governing variable and algorithm
          - reference : (str, str), raw value and str repr of reference type
          - qc_module : (str, str), raw value and str repr of targeted implementation
          - fcae : (str, str), raw value and str repr of all-electron/frozen-core status

    """
    message: str
    stats: Dict[str, Any]

    def __init__(self, circs: List[str]):

        if circs[0].endswith("_property"):
            driver = "properties"
            derivative_int = "prop"  # differs from QCSchema Driver.derivative_int that uses `0` here
        elif circs[0].endswith("_hessian"):
            driver = "hessian"
            derivative_int = 2
        elif circs[0].endswith("_gradient"):
            driver = "gradient"
            derivative_int = 1
        else:
            driver = "energy"
            derivative_int = 0

        docslink = docs_table_link(circs[1], mode=("summary" if circs[5] == "" else "details"))
        all_electron = circs[6] if (len(circs) == 7) else "(n/a)"

        stats = {
            # use QCSchema AtomicInput keys where available
            "driver": driver,
            "derivative_int": derivative_int,
            "method": circs[1],
            "link": docslink,
            "method_type": (circs[3], f"{circs[2].upper()}={circs[3]}"),
            "reference": (circs[4], f"REFERENCE={circs[4]}"),
            "qc_module": (circs[5], f"QC_MODULE={circs[5] or '(auto)'}"),
            # Note on fcae: most select_* functions don't consider this circumstance. Index 0 of tuple documents this.
            #   Index 1 of tuple always returns the circumstance, regardless of whether considered, so user can look up availability in table.
            "fcae": (all_electron, f"FREEZE_CORE={core.get_global_option('FREEZE_CORE')}"),
        }

        # Note that this message is not generally seen, as driver_util catches and edits it from `stats`
        conditions2 = [stats[k][1] for k in ["method_type", "reference", "fcae", "qc_module"]]
        msg = f"Method={stats['method']} is not available for target derivative level (dertype={derivative_int}) under conditions {', '.join(conditions2)}. See {stats['link']}."

        PsiException.__init__(self, msg)
        self.message = msg
        self.stats = stats


# Dftd3Error ceased to be used by v1.4. Class removed by v1.7
# class Dftd3Error(PsiException):


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

    >>> cmake -S. -Bobjdir -Dpsi4_DIR=$PSI4_INSTALL_PREFIX/share/cmake/psi4 {module_args}
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


def sanitize_method(name: str) -> str:
    """Replace characters in method name so that suitable for function name or Sphinx table anchor
    (both as-is) or HTML table link (replace underscore in returned string by dash).

    """
    return name.lower(
        ).replace("(", "_pr"  # ccsd(t)
        ).replace(")", "_pr"
        ).replace(".", "p"    # mp2.5
        ).replace("+", "p"    # ccsd+t(ccsd)
        ).replace("-", ""     # ccsdt-1a
        )


def docs_table_link(name: str, mode: str) -> str:
    """Compose a link to *mode* documentation table.
    For modes {"summary", "details"}, link will be anchored at the table row for method *name*. (Method row not guaranteed to exist.)
    For other modes, link will be to table in general.

    """
    DOCS_BASE = "https://psicode.org/psi4manual/master/"

    if mode == "summary":
        anchor = "introduction.html#ss"
    elif mode == "details":
        anchor = "capabilities.html#dd"
    elif mode == "ccenergy":
        anchor = "cc.html#table-ccenergy-stdsuite"
    elif mode == "dfmp2":
        anchor = "dfmp2.html#table-dfmp2-stdsuite"
    elif mode == "fnocc":
        anchor = "fnocc.html#table-fnocc-stdsuite"
    elif mode == "occ_nonoo":
        anchor = "occ.html#table-occ-stdsuite-nonoo"
    elif mode == "occ_oo":
        anchor = "occ.html#table-occ-stdsuite-oo"
    elif mode == "scf":
        anchor = "dft.html#table-scf-stdsuite"
    else:
        raise KeyError("invalid table mode")

    # Sphinx reST anchors with underscore build into HTML docs with dash
    sanitized_method = sanitize_method(name).replace("_", "-")

    if mode in ["summary", "details"]:
        return f"{DOCS_BASE}{anchor}-{sanitized_method}"
    else:
        return f"{DOCS_BASE}{anchor}"
