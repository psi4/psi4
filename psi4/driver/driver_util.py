#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2023 The Psi4 Developers.
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

__all__ = [
    "highest_analytic_derivative_available",
    "highest_analytic_properties_available",
    "negotiate_convergence_criterion",
    "negotiate_derivative_type",
    "parse_cotton_irreps",
    "sort_derivative_type",
    "upgrade_interventions",
]

import math
from typing import Any, Dict, Optional, Tuple, Union

from psi4 import core

from . import p4util
from .p4util.exceptions import ManagedMethodError, MissingMethodError, UpgradeHelper, ValidationError, docs_table_link
from .procrouting import proc
from .procrouting.proc_table import procedures


def negotiate_convergence_criterion(dermode: Union[Tuple[str, str], Tuple[int, int]], method: str, return_optstash: bool = False):
    r"""
    This function will set local SCF and global energy convergence criterion
    to the defaults listed at:
    http://www.psicode.org/psi4manual/master/scf.html#convergence-and-
    algorithm-defaults. SCF will be converged more tightly (pscf_Ec and pscf_Dc)
    if a post-SCF method is selected. For a final SCF method, the looser
    (scf_Ec and scf_Dc) convergence criterion will be used.

    dermode -       Tuple of target and means derivative level or ("prop", "prop"). E.g., analytic gradient
                    is (1, 1); frequency by energies is (2, 0). Nearly always test on
                    procedures['energy'] since that's guaranteed to exist for a method.
    method -        Name of the method
    scf_Ec -        E convergence criterion for scf             target method
    pscf_Ec -       E convergence criterion for scf of post-scf target method
    scf_Dc -        D convergence criterion for scf             target method
    pscf_Dc -       D convergence criterion for scf of post-scf target method
    gen_Ec -        E convergence criterion for        post-scf target method

    """

    scf_Ec, pscf_Ec, scf_Dc, pscf_Dc, gen_Ec = {
        (0, 0): [6, 8, 6, 8, 6],
        (1, 0): [8, 10, 8, 10, 8],
        (2, 0): [10, 11, 10, 11, 10],
        (1, 1): [8, 10, 8, 10, 8],
        (2, 1): [8, 10, 8, 10, 8],
        (2, 2): [8, 10, 8, 10, 8],
        ("prop", "prop"): [6, 10, 6, 10, 8]
    }[dermode]

    # Set method-dependent scf convergence criteria, check against energy routines
    # Set post-scf convergence criteria (global will cover all correlated modules)
    cc = {}
    if procedures['energy'][method] in [proc.run_scf, proc.run_tdscf_energy]:
        if not core.has_option_changed('SCF', 'E_CONVERGENCE'):
            cc['SCF__E_CONVERGENCE'] = math.pow(10, -scf_Ec)
        if not core.has_option_changed('SCF', 'D_CONVERGENCE'):
            cc['SCF__D_CONVERGENCE'] = math.pow(10, -scf_Dc)
    else:
        if not core.has_option_changed('SCF', 'E_CONVERGENCE'):
            cc['SCF__E_CONVERGENCE'] = math.pow(10, -pscf_Ec)
        if not core.has_option_changed('SCF', 'D_CONVERGENCE'):
            cc['SCF__D_CONVERGENCE'] = math.pow(10, -pscf_Dc)
        if not core.has_global_option_changed('E_CONVERGENCE'):
            cc['E_CONVERGENCE'] = math.pow(10, -gen_Ec)

    if return_optstash:
        optstash = p4util.OptionsState(['SCF', 'E_CONVERGENCE'], ['SCF', 'D_CONVERGENCE'], ['E_CONVERGENCE'])
        p4util.set_options(cc)
        return optstash

    else:
        return cc


def upgrade_interventions(method):
    try:
        lowermethod = method.lower()
    except AttributeError as e:
        if method.__name__ == 'cbs':
            raise UpgradeHelper(method, repr(method.__name__), 1.6,
                                ' Replace cbs or complete_basis_set function with cbs string.')
        elif method.__name__ in ['sherrill_gold_standard', 'allen_focal_point']:
            raise UpgradeHelper(
                method.__name__, '"' + method.__name__ + '"', 1.6,
                f' Replace function `energy({method.__name__})` with string `energy("{method.__name__}")` or similar.')
        else:
            raise e

    if lowermethod.startswith("mrcc"):
        raise UpgradeHelper(method, method[2:], 1.7,
            f' Replace "mr"-prefixed methods in driver calls (e.g., `energy("{method}")`) by the plain method and specify the MRCC addon by option (e.g., `set qc_module mrcc; energy("{method[2:]}")`).')

    return lowermethod


def parse_cotton_irreps(irrep: Union[str, int], point_group: str) -> int:
    """Return validated Cotton ordering index of `irrep` within `point_group`.

    Parameters
    ----------
    irrep
        Irreducible representation. Either label (case-insensitive) or 1-based index (int or str).
    point_group
        Molecular point group label (case-insensitive).

    Returns
    -------
    int
        1-based index for **irrep** within **point_group** in Cotton ordering.

    Raises
    ------
    ValidationError
        If **irrep** out-of-bounds or invalid or if **point_group** doesn't exist.

    """
    cotton = {
        'c1': ['a'],
        'ci': ['ag', 'au'],
        'c2': ['a', 'b'],
        'cs': ['ap', 'app'],
        'd2': ['a', 'b1', 'b2', 'b3'],
        'c2v': ['a1', 'a2', 'b1', 'b2'],
        'c2h': ['ag', 'bg', 'au', 'bu'],
        'd2h': ['ag', 'b1g', 'b2g', 'b3g', 'au', 'b1u', 'b2u', 'b3u'],
    }

    boll = cotton[point_group.lower()]

    if str(irrep).isdigit():
        irrep = int(irrep)
        if irrep > 0 and irrep <= len(boll):
            return irrep
    else:
        if irrep.lower() in boll:
            return boll.index(irrep.lower()) + 1

    raise ValidationError(f"""Irrep '{irrep}' not valid for point group '{point_group}'.""")


def negotiate_derivative_type(
    target_dertype: str,
    method: str,
    user_dertype: Optional[int],
    verbose: int = 1,
    proc: Optional[Dict] = None,
) -> Union[Tuple[int, int], Tuple[str, str]]:
    r"""Find the best derivative level (0, 1, 2) and strategy (analytic, finite difference)
    for `method` to achieve `target_dertype` within constraints of `user_dertype`.

    Parameters
    ----------
    target_dertype: {"energy", "gradient", "hessian", "properties"}
        Type of calculation targeted by driver.
    method
        Quantum chemistry method targeted by driver. Should be correct case for procedures lookup.
    user_dertype
        User input on which derivative level should be employed to achieve `target_dertype`.
    verbose
        Control amount of output printing.
    proc
        For testing only! Procedures table to look up `method`. Default is psi4.driver.procedures .
    managed_keywords
        NYI Keywords that influence managed methods.

    Returns
    -------
    tuple : (int, int)
        "second" is highest accessible derivative level for `method` to achieve
        `target_dertype` "first" within constraints of `user_dertype`. When
        "first" == "second", analytic is the best strategy, otherwise finite
        difference of target "first" by means of "second".

    Raises
    ------
    ValidationError
        When input validation fails. When `user_dertype` exceeds `target_dertype`.
    MissingMethodError
        When `method` is unavailable at all. When `user_dertype` exceeds what available for `method`.

    """
    if target_dertype == "properties":
        if isinstance(method, list):
            # untested
            analytic = [highest_analytic_properties_available(mtd, proc) for mtd in method]
            if verbose > 1:
                print(method, analytic, '->', min(analytic, key=lambda t: t[0]))
            highest_analytic_dertype, proc_messages = min(analytic, key=lambda t: t[0])
        else:
            highest_analytic_dertype, proc_messages = highest_analytic_properties_available(method, proc)

        return (target_dertype[:4], highest_analytic_dertype)  # i.e., ("prop", "prop")

    else:
        if isinstance(method, list):
            analytic = [highest_analytic_derivative_available(mtd, proc) for mtd in method]
            if verbose > 1:
                print(method, analytic, '->', min(analytic, key=lambda t: t[0]))
            highest_analytic_dertype, proc_messages = min(analytic, key=lambda t: t[0])
        else:
            highest_analytic_dertype, proc_messages = highest_analytic_derivative_available(method, proc)

        return sort_derivative_type(target_dertype, highest_analytic_dertype, user_dertype, proc_messages, verbose)


def _alternative_methods_message(method_name: str, dertype: str, *, messages: Dict[int, Any], proc: Dict) -> str:
    """Format error message when *method_name* not available, whether at all in *proc* or simply not under current conditions.

    Parameters
    ----------
    method_name
        See caller.
    dertype
        Always "any"
    messages
        Dictionary of returned error statistics from ManagedMethodError.
        Entry for energy (key `0`) used if present. Alternate message composed if empty dict.
    proc
        See caller.

    Returns
    -------
    str
        Message saying not available and suggesting some alternatives in case of typo. If the
        method was probed under conditions and rejected (*messages* non-empty), the message
        includes the conditions and a link to docs table.

    """
    alt_method_name = p4util.text.find_approximate_string_matches(method_name, proc["energy"].keys(), 2)
    if method_name in alt_method_name:
        alt_method_name.remove(method_name)

    alternatives = ''
    query = "Or did you mean?" if messages else "Did you mean?"
    if len(alt_method_name) > 0:
        alternatives = f" {query} {' '.join(alt_method_name)}"

    assert dertype == "any"
    if messages:
        stats = messages[0]
        conditions2 = [stats[k][1] for k in ["method_type", "reference", "fcae", "qc_module"]]
        return f"Method={stats['method']} is not available for {dertype} derivative level under conditions {', '.join(conditions2)}. See {stats['link']}.{alternatives}"
    else:
        return f"Method={method_name} is not available for {dertype} derivative level.{alternatives}"


def highest_analytic_derivative_available(method: str,
                                          proc: Optional[Dict] = None,
                                          managed_keywords: Optional[Dict] = None) -> Tuple[int, Dict[int,str]]:
    """Find the highest dertype program can provide for method, as encoded in procedures and managed methods.

    Managed methods return finer grain "is available" info. For example, "is analytic ROHF DF HF gradient available?"
    from managed method, not just "is HF gradient available?" from procedures.

    Parameters
    ----------
    method
        Quantum chemistry method targeted by driver. Should be correct case for procedures lookup.
    proc
        For testing only! Procedures table to look up `method`. Default is psi4.driver.procedures .
    managed_keywords
        Keywords that influence managed methods.

    Returns
    -------
    int
        Highest available analytic derivative for `method`.
    dict
        Detailed error messages to be passed along. Keys are dertype.

    Raises
    ------
    MissingMethodError
        When `method` is unavailable at all. When `user_dertype` exceeds what available for `method`.

    """
    if managed_keywords is None:
        managed_keywords = {}

    if proc is None:
        proc = procedures

    dertype = "(auto)"
    proc_messages = {}

    if method in proc['hessian']:
        dertype = 2
        if proc['hessian'][method].__name__.startswith('select_'):
            try:
                proc["hessian"][method](method, probe=True, **managed_keywords)
            except ManagedMethodError as e:
                proc_messages[2] = e.stats
                dertype = 1
                if proc['gradient'][method].__name__.startswith('select_'):
                    try:
                        proc["gradient"][method](method, probe=True, **managed_keywords)
                    except ManagedMethodError as e:
                        proc_messages[1] = e.stats
                        dertype = 0
                        if proc['energy'][method].__name__.startswith('select_'):
                            try:
                                proc["energy"][method](method, probe=True, **managed_keywords)
                            except ManagedMethodError as e:
                                proc_messages[0] = e.stats
                                raise MissingMethodError(_alternative_methods_message(method, "any", messages=proc_messages, proc=proc))

    elif method in proc['gradient']:
        dertype = 1
        proc_messages[2] = {"method": method, "link": docs_table_link(method, mode="summary")}
        if proc['gradient'][method].__name__.startswith('select_'):
            try:
                proc["gradient"][method](method, probe=True, **managed_keywords)
            except ManagedMethodError as e:
                proc_messages[1] = e.stats
                dertype = 0
                if proc['energy'][method].__name__.startswith('select_'):
                    try:
                        proc["energy"][method](method, probe=True, **managed_keywords)
                    except ManagedMethodError as e:
                        proc_messages[0] = e.stats
                        raise MissingMethodError(_alternative_methods_message(method, "any", messages=proc_messages, proc=proc))

    elif method in proc['energy']:
        dertype = 0
        proc_messages[1] = {"method": method, "link": docs_table_link(method, mode="summary")}
        proc_messages[2] = proc_messages[1]
        if proc['energy'][method].__name__.startswith('select_'):
            try:
                proc["energy"][method](method, probe=True, **managed_keywords)
            except ManagedMethodError as e:
                proc_messages[0] = e.stats
                raise MissingMethodError(_alternative_methods_message(method, "any", messages=proc_messages, proc=proc))

    if dertype == '(auto)':
        raise MissingMethodError(_alternative_methods_message(method, "any", messages=proc_messages, proc=proc))

    if dertype == 2 and p4util.libint2_configuration()["eri"][2] is None:
        dertype = 1
        proc_messages[2] = {"method": method, "blame": "Libint2 build"}
        #core.print_out("  Warning: Analytical Hessians not available with this Libint2 library. Falling back to finite difference. Setting `points=5` may be needed for precision.\n")

    return dertype, proc_messages


def highest_analytic_properties_available(method: str,
                                          proc: Optional[Dict] = None,
                                          managed_keywords: Optional[Dict] = None) -> Tuple[str, Dict[int,str]]:
    """Find whether propgram can provide analytic properties for method, as encoded in procedures and managed methods.

    Managed methods return finer grain "is available" info. For example, "is analytic ROHF DF HF gradient available?"
    from managed method, not just "is HF gradient available?" from procedures.

    Parameters
    ----------
    method
        Quantum chemistry method targeted by driver. Should be correct case for procedures lookup.
    proc
        For testing only! Procedures table to look up `method`. Default is psi4.driver.procedures .
    managed_keywords
        Keywords that influence managed methods.

    Returns
    -------
    str
        "prop" if analytic properties available for `method`.
    dict
        Detailed error messages to be passed along. Keys are dertype.

    Raises
    ------
    MissingMethodError
        When `method` is unavailable at all. When `user_dertype` exceeds what available for `method`.

    """
    if managed_keywords is None:
        managed_keywords = {}

    if proc is None:
        proc = procedures

    dertype = "(auto)"
    proc_messages = {}

    if method in proc["properties"]:
        dertype = "prop"
        if proc["properties"][method].__name__.startswith('select_'):
            try:
                proc["properties"][method](method, probe=True, **managed_keywords)
            except ManagedMethodError as e:
                proc_messages[0] = e.message
                raise MissingMethodError(_alternative_methods_message(method, "any", messages=proc_messages, proc=proc))

    if dertype == '(auto)':
        raise MissingMethodError(_alternative_methods_message(method, "any", messages=proc_messages, proc=proc))

    return dertype, proc_messages


def sort_derivative_type(
    target_dertype,
    highest_analytic_dertype: int,
    user_dertype: Union[int, None],
    proc_messages: Dict[int, str],
    verbose: int = 1,
) -> Tuple[int, int]:
    r"""Find the best derivative level (0, 1, 2) and strategy (analytic, finite difference)
    to achieve `target_dertype` within constraints of `user_dertype` and `highest_analytic_dertype`.

    Parameters
    ----------
    target_dertype: {'energy', 'gradient', 'hessian'}
        Type of calculation targeted by driver.
    highest_analytic_dertype
        Highest derivative level program can provide analytically.
    user_dertype
        User input on which derivative level should be employed to achieve `target_dertype`.
    proc_messages
        Dertype-indexed detailed message to be appended to user error.
    verbose
        Control amount of output printing.

    Returns
    -------
    tuple : (int, int)
        "second" is highest accessible derivative level `highest_analytic_dertype` to achieve
        `target_dertype` "first" within constraints of `user_dertype`. When
        "first" == "second", analytic is the best strategy, otherwise finite
        difference of target "first" by means of "second".

    Raises
    ------
    ValidationError
        When input validation fails. When `user_dertype` exceeds `target_dertype`.
    MissingMethodError
        When `user_dertype` exceeds what available for `method`.

    """
    egh = ['energy', 'gradient', 'hessian']

    # Validate input dertypes
    if target_dertype not in egh:
        raise ValidationError(f"target_dertype ({target_dertype}) must be in {egh}.")

    if not (user_dertype is None or isinstance(user_dertype, int)):
        raise ValidationError(f"user_dertype ({user_dertype}) should only be None or int!")

    dertype = highest_analytic_dertype

    # Negotiations. In particular:
    # * don't return higher derivative than targeted by driver
    # * don't return higher derivative than spec'd by user. that is, user can downgrade derivative
    # * alert user to conflict between driver and user_dertype

    if user_dertype is not None and user_dertype > egh.index(target_dertype):
        raise ValidationError(f'User dertype ({user_dertype}) excessive for target calculation ({target_dertype})')

    if egh.index(target_dertype) < dertype:
        dertype = egh.index(target_dertype)

    if user_dertype is not None:
        if user_dertype <= dertype:
            dertype = user_dertype
        else:
            stats = proc_messages[user_dertype]
            msg = f"Method={stats['method']} is not available for requested derivative level (reqd={user_dertype} > avail={highest_analytic_dertype})"
            try:
                conditions2 = [stats[k][1] for k in ["method_type", "reference", "fcae", "qc_module"]]
                msg += f" under conditions {', '.join(conditions2)}. See {stats['link']}."
            except KeyError:
                if "blame" in stats:
                    msg += f" under {stats['blame']} conditions."
                else:
                    msg += f" under any conditions. See (possibly) {stats['link']}."
            raise MissingMethodError(msg)

    # hack section
    if core.get_global_option('PCM') and dertype != 0:
        core.print_out('\nPCM analytic gradients are not implemented yet, re-routing to finite differences.\n')
        dertype = 0

    if core.get_global_option("RELATIVISTIC") in ["X2C", "DKH"]:
        core.print_out("\nRelativistic analytic gradients are not implemented yet, re-routing to finite differences.\n")
        dertype = 0

    if verbose > 1:
        print(
            f'Derivative negotiations: target/driver={egh.index(target_dertype)}, best_available={highest_analytic_dertype}, user_dertype={user_dertype} -> FINAL={dertype}'
        )

    return (egh.index(target_dertype), dertype)
