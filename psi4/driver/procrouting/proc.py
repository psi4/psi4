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

"""Module with functions that encode the sequence of PSI module
calls for each of the *name* values of the energy(), optimize(),
response(), and frequency() function. *name* can be assumed lowercase by here.

"""
import re
import os
import sys
import shutil
import subprocess
import warnings
from typing import Dict, List, Union

import numpy as np
from qcelemental import constants
from qcelemental.util import which

from psi4 import extras
from psi4 import core
from psi4.driver import p4util
from psi4.driver import qcdb
from psi4.driver import psifiles as psif
from psi4.driver.p4util.exceptions import ManagedMethodError, PastureRequiredError, UpgradeHelper, ValidationError, docs_table_link
#from psi4.driver.molutil import *
from psi4.driver.qcdb.basislist import corresponding_basis
# never import driver, wrappers, or aliases into this file

from .proc_data import method_algorithm_type
from .roa import run_roa
from . import proc_util
from . import empirical_dispersion
from . import dft
from . import mcscf
from . import response
from . import solvent


# ADVICE on new additions:
# * two choices: basic `def run` or managed `def select`
# * consult http://psicode.org/psi4manual/master/proc_py.html  --or--  <psi4-repo>/doc/sphinxman/source/proc_py.rst


def select_scf_gradient(name, **kwargs):
    """Function selecting the algorithm for an SCF gradient call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    type_var, _, mtd_type = method_algorithm_type("scf")  # `"scf"` instead of `name` avoids adding every functional to governing dict in proc_data.py
    module = core.get_global_option('QC_MODULE')

    if mtd_type == 'CD':
        # manifestation of `"""No analytic derivatives for SCF_TYPE CD."""`.
        #   here, only hits upon `gradient("scf")` so above message also present in driver.py to catch e.g., mp2 gradient atop a cd reference.
        func = None
    else:
        func = run_scf_gradient

    if func is None:
        raise ManagedMethodError([__name__, name, type_var, mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_mp2(name, **kwargs):
    """Function selecting the algorithm for a MP2 energy call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    type_var, _, mtd_type = method_algorithm_type(name)
    module = core.get_global_option('QC_MODULE')

    func = None
    if reference == 'RHF':
        if mtd_type == 'CONV':
            if module == 'DETCI':
                func = run_detci
            elif module == 'FNOCC':
                func = run_fnocc
            elif module in ['', 'OCC']:
                func = run_occ
        elif mtd_type == 'DF':
            if module == 'OCC':
                func = run_dfocc
            elif module in ['', 'DFMP2']:
                func = run_dfmp2
        elif mtd_type == 'CD':
            if module in ['', 'OCC']:
                func = run_dfocc
    elif reference == 'UHF':
        if mtd_type == 'CONV':
            if module in ['', 'OCC']:
                func = run_occ
        elif mtd_type == 'DF':
            if module == 'OCC':
                func = run_dfocc
            elif module in ['', 'DFMP2']:
                func = run_dfmp2
        elif mtd_type == 'CD':
            if module in ['', 'OCC']:
                func = run_dfocc
    elif reference == 'ROHF':
        if mtd_type == 'CONV':
            if module == 'DETCI':
                raise UpgradeHelper("energy('mp2')", "energy('zapt2')", 1.7,
                    " Replace method MP with method ZAPT for ROHF reference. DETCI is orders-of-magnitude inefficient for perturbation theory.")
            elif module in ['', 'OCC']:
                func = run_occ
        elif mtd_type == 'DF':
            if module == 'OCC':
                func = run_dfocc
            elif module in ['', 'DFMP2']:
                func = run_dfmp2
        elif mtd_type == 'CD':
            if module in ['', 'OCC']:
                func = run_dfocc
    elif reference in ['RKS', 'UKS']:
        if mtd_type == 'DF':
            if module in ['', 'DFMP2']:
                func = run_dfmp2

    if module == 'DETCI':
        core.print_out("""\nDETCI is ill-advised for method MP2 as it is available inefficiently as a """
                       """byproduct of a CISD computation.\n""")

    if func is None:
        raise ManagedMethodError([__name__, name, type_var, mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_mp2_gradient(name, **kwargs):
    """Function selecting the algorithm for a MP2 gradient call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    type_var, _, mtd_type = method_algorithm_type(name)
    module = core.get_global_option('QC_MODULE')
    all_electron = (core.get_global_option('FREEZE_CORE') == "FALSE")

    func = None
    if reference == 'RHF':
        if mtd_type == 'CONV':
            if all_electron:
                if module in ['', 'OCC']:
                    func = run_occ_gradient
        elif mtd_type == 'DF':
                if module == 'OCC':
                    func = run_dfocc_gradient
                elif module in ['', 'DFMP2']:
                    func = run_dfmp2_gradient
    elif reference == 'UHF':
        if mtd_type == 'CONV':
            if all_electron:
                if module in ['', 'OCC']:
                    func = run_occ_gradient
        elif mtd_type == 'DF':
                if module in ['', 'OCC']:
                    func = run_dfocc_gradient

    if func is None:
        raise ManagedMethodError([__name__, name, type_var, mtd_type, reference, module, all_electron])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_mp2_property(name, **kwargs):
    """Function selecting the algorithm for a MP2 property call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    type_var, _, mtd_type = method_algorithm_type(name)
    module = core.get_global_option('QC_MODULE')

    func = None
    if reference == 'RHF':
        if mtd_type == 'DF':
            #if module == 'OCC':
            #    func = run_dfocc_property
            if module in ['', 'DFMP2']:
                func = run_dfmp2_property
    #elif reference == 'UHF':
    #    if mtd_type == 'DF':
    #        if module in ['', 'OCC']:
    #            func = run_dfocc_property

    if func is None:
        raise ManagedMethodError([__name__, name, type_var, mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_omp2(name, **kwargs):
    """Function selecting the algorithm for an OMP2 energy call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    type_var, _, mtd_type = method_algorithm_type(name)
    module = core.get_global_option('QC_MODULE')

    func = None
    if reference in ['RHF', 'UHF', 'ROHF', 'RKS', 'UKS']:
        if mtd_type == 'CONV':
            if module in ['', 'OCC']:
                func = run_occ
        elif mtd_type == 'DF':
            if module in ['', 'OCC']:
                func = run_dfocc
        elif mtd_type == 'CD':
            if module in ['', 'OCC']:
                func = run_dfocc

    if func is None:
        raise ManagedMethodError([__name__, name, type_var, mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_omp2_gradient(name, **kwargs):
    """Function selecting the algorithm for an OMP2 gradient call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    type_var, _, mtd_type = method_algorithm_type(name)
    module = core.get_global_option('QC_MODULE')

    func = None
    if reference in ['RHF', 'UHF', 'ROHF', 'RKS', 'UKS']:
        if mtd_type == 'CONV':
            if module in ['', 'OCC']:
                func = run_occ_gradient
        elif mtd_type == 'DF':
            if module in ['', 'OCC']:
                func = run_dfocc_gradient

    if func is None:
        raise ManagedMethodError([__name__, name, type_var, mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_omp2_property(name, **kwargs):
    """Function selecting the algorithm for an OMP2 property call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    type_var, _, mtd_type = method_algorithm_type(name)
    module = core.get_global_option('QC_MODULE')

    func = None
    if reference in ['RHF', 'UHF', 'ROHF', 'RKS', 'UKS']:
        if mtd_type == 'DF':
            if module in ['', 'OCC']:
                func = run_dfocc_property

    if func is None:
        raise ManagedMethodError([__name__, name, type_var, mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_omp2p5_property(name, **kwargs):
    """Function selecting the algorithm for an OMP2.5 property call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    type_var, _, mtd_type = method_algorithm_type(name)
    module = core.get_global_option('QC_MODULE')

    func = None
    if reference in ['RHF', 'UHF', 'ROHF', 'RKS', 'UKS']:
        if mtd_type == 'DF':
            if module in ['', 'OCC']:
                func = run_dfocc_property

    if func is None:
        raise ManagedMethodError([__name__, name, type_var, mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_omp3_property(name, **kwargs):
    """Function selecting the algorithm for an OMP3 property call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    type_var, _, mtd_type = method_algorithm_type(name)
    module = core.get_global_option('QC_MODULE')

    func = None
    if reference in ['RHF', 'UHF', 'ROHF', 'RKS', 'UKS']:
        if mtd_type == 'DF':
            if module in ['', 'OCC']:
                func = run_dfocc_property

    if func is None:
        raise ManagedMethodError([__name__, name, type_var, mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_olccd_property(name, **kwargs):
    """Function selecting the algorithm for an OLCCD property call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    type_var, _, mtd_type = method_algorithm_type(name)
    module = core.get_global_option('QC_MODULE')

    func = None
    if reference in ['RHF', 'UHF', 'ROHF', 'RKS', 'UKS']:
        if mtd_type == 'DF':
            if module in ['', 'OCC']:
                func = run_dfocc_property

    if func is None:
        raise ManagedMethodError([__name__, name, type_var, mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_mp3(name, **kwargs):
    """Function selecting the algorithm for a MP3 energy call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    type_var, _, mtd_type = method_algorithm_type(name)
    module = core.get_global_option('QC_MODULE')

    func = None
    if reference == 'RHF':
        if mtd_type == 'CONV':
            if module == 'DETCI':
                func = run_detci
            elif module == 'FNOCC':
                func = run_fnocc
            elif module in ['', 'OCC']:
                func = run_occ
        elif mtd_type == 'DF':
            if module in ['', 'OCC']:
                func = run_dfocc
        elif mtd_type == 'CD':
            if module in ['', 'OCC']:
                func = run_dfocc
    elif reference == 'UHF':
        if mtd_type == 'CONV':
            if module in ['', 'OCC']:
                func = run_occ
        elif mtd_type == 'DF':
            if module in ['', 'OCC']:
                func = run_dfocc
        elif mtd_type == 'CD':
            if module in ['', 'OCC']:
                func = run_dfocc
    elif reference == 'ROHF':
        if mtd_type == 'CONV':
            if module == 'DETCI':
                raise UpgradeHelper("energy('mp3')", "energy('zapt3')", 1.7,
                    " Replace method MP with method ZAPT for ROHF reference. DETCI is orders-of-magnitude inefficient for perturbation theory.")

    if module == 'DETCI':
        core.print_out("""\nDETCI is ill-advised for method MP3 as it is available inefficiently as a """
                       """byproduct of a CISD computation.\n""")

    if func is None:
        raise ManagedMethodError([__name__, name, type_var, mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_mp3_gradient(name, **kwargs):
    """Function selecting the algorithm for a MP3 gradient call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    type_var, _, mtd_type = method_algorithm_type(name)
    module = core.get_global_option('QC_MODULE')
    all_electron = (core.get_global_option('FREEZE_CORE') == "FALSE")

    func = None
    if reference == 'RHF':
        if mtd_type == 'CONV':
            if all_electron:
                if module in ['', 'OCC']:
                    func = run_occ_gradient
        elif mtd_type == 'DF':
                if module in ['', 'OCC']:
                    func = run_dfocc_gradient
    elif reference == 'UHF':
        if mtd_type == 'CONV':
            if all_electron:
                if module in ['', 'OCC']:
                    func = run_occ_gradient
        elif mtd_type == 'DF':
                if module in ['', 'OCC']:
                    func = run_dfocc_gradient

    if func is None:
        raise ManagedMethodError([__name__, name, type_var, mtd_type, reference, module, all_electron])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_omp3(name, **kwargs):
    """Function selecting the algorithm for an OMP3 energy call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    type_var, _, mtd_type = method_algorithm_type(name)
    module = core.get_global_option('QC_MODULE')

    func = None
    if reference in ['RHF', 'UHF', 'ROHF', 'RKS', 'UKS']:
        if mtd_type == 'CONV':
            if module in ['', 'OCC']:
                func = run_occ
        elif mtd_type == 'DF':
            if module in ['', 'OCC']:
                func = run_dfocc
        elif mtd_type == 'CD':
            if module in ['', 'OCC']:
                func = run_dfocc

    if func is None:
        raise ManagedMethodError([__name__, name, type_var, mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_omp3_gradient(name, **kwargs):
    """Function selecting the algorithm for an OMP3 gradient call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    type_var, _, mtd_type = method_algorithm_type(name)
    module = core.get_global_option('QC_MODULE')

    func = None
    if reference in ['RHF', 'UHF', 'ROHF', 'RKS', 'UKS']:
        if mtd_type == 'CONV':
            if module in ['', 'OCC']:
                func = run_occ_gradient
        elif mtd_type == 'DF':
            if module in ['', 'OCC']:
                func = run_dfocc_gradient

    if func is None:
        raise ManagedMethodError([__name__, name, type_var, mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_mp2p5(name, **kwargs):
    """Function selecting the algorithm for a MP2.5 energy call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    type_var, _, mtd_type = method_algorithm_type(name)
    module = core.get_global_option('QC_MODULE')

    func = None
    if reference in ['RHF', 'UHF']:
        if mtd_type == 'CONV':
            if module in ['', 'OCC']:
                func = run_occ
        elif mtd_type == 'DF':
            if module in ['', 'OCC']:
                func = run_dfocc
        elif mtd_type == 'CD':
            if module in ['', 'OCC']:
                func = run_dfocc

    if func is None:
        raise ManagedMethodError([__name__, name, type_var, mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_mp2p5_gradient(name, **kwargs):
    """Function selecting the algorithm for a MP2.5 gradient call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    type_var, _, mtd_type = method_algorithm_type(name)
    module = core.get_global_option('QC_MODULE')
    all_electron = (core.get_global_option('FREEZE_CORE') == "FALSE")

    func = None
    if reference in ['RHF', 'UHF']:
        if mtd_type == 'CONV':
            if all_electron:
                if module in ['', 'OCC']:
                    func = run_occ_gradient
        elif mtd_type == 'DF':
                if module in ['', 'OCC']:
                    func = run_dfocc_gradient

    if func is None:
        raise ManagedMethodError([__name__, name, type_var, mtd_type, reference, module, all_electron])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_omp2p5(name, **kwargs):
    """Function selecting the algorithm for an OMP2.5 energy call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    type_var, _, mtd_type = method_algorithm_type(name)
    module = core.get_global_option('QC_MODULE')

    func = None
    if reference in ['RHF', 'UHF', 'ROHF', 'RKS', 'UKS']:
        if mtd_type == 'CONV':
            if module in ['', 'OCC']:
                func = run_occ
        elif mtd_type == 'DF':
            if module in ['', 'OCC']:
                func = run_dfocc
        elif mtd_type == 'CD':
            if module in ['', 'OCC']:
                func = run_dfocc

    if func is None:
        raise ManagedMethodError([__name__, name, type_var, mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_omp2p5_gradient(name, **kwargs):
    """Function selecting the algorithm for an OMP2.5 gradient call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    type_var, _, mtd_type = method_algorithm_type(name)
    module = core.get_global_option('QC_MODULE')

    func = None
    if reference in ['RHF', 'UHF', 'ROHF', 'RKS', 'UKS']:
        if mtd_type == 'CONV':
            if module in ['', 'OCC']:
                func = run_occ_gradient
        elif mtd_type == 'DF':
            if module in ['', 'OCC']:
                func = run_dfocc_gradient

    if func is None:
        raise ManagedMethodError([__name__, name, type_var, mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_lccd(name, **kwargs):
    """Function selecting the algorithm for a LCCD energy call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    type_var, _, mtd_type = method_algorithm_type(name)
    module = core.get_global_option('QC_MODULE')

    func = None
    if reference == 'RHF':
        if mtd_type == 'CONV':
            if module == 'OCC':
                func = run_occ
            elif module in ['', 'FNOCC']:
                func = run_cepa
        elif mtd_type == 'DF':
            if module in ['', 'OCC']:
                func = run_dfocc
        elif mtd_type == 'CD':
            if module in ['', 'OCC']:
                func = run_dfocc
    elif reference == 'UHF':
        if mtd_type == 'CONV':
            if module in ['', 'OCC']:
                func = run_occ
        elif mtd_type == 'DF':
            if module in ['', 'OCC']:
                func = run_dfocc
        elif mtd_type == 'CD':
            if module in ['', 'OCC']:
                func = run_dfocc

    if func is None:
        raise ManagedMethodError([__name__, name, type_var, mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_lccd_gradient(name, **kwargs):
    """Function selecting the algorithm for a LCCD gradient call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    type_var, _, mtd_type = method_algorithm_type(name)
    module = core.get_global_option('QC_MODULE')
    all_electron = (core.get_global_option('FREEZE_CORE') == "FALSE")

    func = None
    if reference in ['RHF', 'UHF']:
        if mtd_type == 'CONV':
            if all_electron:
                if module in ['', 'OCC']:
                    func = run_occ_gradient
        elif mtd_type == 'DF':
                if module in ['', 'OCC']:
                    func = run_dfocc_gradient

    if func is None:
        raise ManagedMethodError([__name__, name, type_var, mtd_type, reference, module, all_electron])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_olccd(name, **kwargs):
    """Function selecting the algorithm for an OLCCD energy call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    type_var, _, mtd_type = method_algorithm_type(name)
    module = core.get_global_option('QC_MODULE')

    func = None
    if reference in ['RHF', 'UHF', 'ROHF', 'RKS', 'UKS']:
        if mtd_type == 'CONV':
            if module in ['', 'OCC']:
                func = run_occ
        elif mtd_type == 'DF':
            if module in ['', 'OCC']:
                func = run_dfocc
        elif mtd_type == 'CD':
            if module in ['', 'OCC']:
                func = run_dfocc

    if func is None:
        raise ManagedMethodError([__name__, name, type_var, mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_olccd_gradient(name, **kwargs):
    """Function selecting the algorithm for an OLCCD gradient call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    type_var, _, mtd_type = method_algorithm_type(name)
    module = core.get_global_option('QC_MODULE')

    func = None
    if reference in ['RHF', 'UHF', 'ROHF', 'RKS', 'UKS']:
        if mtd_type == 'CONV':
            if module in ['', 'OCC']:
                func = run_occ_gradient
        elif mtd_type == 'DF':
            if module in ['', 'OCC']:
                func = run_dfocc_gradient

    if func is None:
        raise ManagedMethodError([__name__, name, type_var, mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_fnoccsd(name, **kwargs):
    """Function selecting the algorithm for a FNO-CCSD energy call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    type_var, _, mtd_type = method_algorithm_type(name)
    module = core.get_global_option('QC_MODULE')

    func = None
    if reference == 'RHF':
        if mtd_type == 'CONV':
            if module in ['', 'FNOCC']:
                func = run_fnocc
        elif mtd_type == 'DF':
            if module in ['', 'FNOCC']:
                func = run_fnodfcc
        elif mtd_type == 'CD':
            if module in ['', 'FNOCC']:
                func = run_fnodfcc

    if func is None:
        raise ManagedMethodError([__name__, name, type_var, mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_ccsd(name, **kwargs):
    """Function selecting the algorithm for a CCSD energy call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    type_var, _, mtd_type = method_algorithm_type(name)
    module = core.get_global_option('QC_MODULE')

    # [Aug 2022] DF CCSD through CCENERGY for (RHF|ROHF) not enabled here since not advertised. It does run, though, see #2710

    func = None
    if reference == 'RHF':
        if mtd_type == 'CONV':
            if module == 'FNOCC':
                func = run_fnocc
            elif module == 'MRCC' and which("dmrcc", return_bool=True):
                func = run_mrcc
            elif module == 'CCT3' and extras.addons("cct3"):
                import cct3
                func = cct3.run_cct3
            elif module in ['', 'CCENERGY']:
                func = run_ccenergy
        elif mtd_type == 'DF':
            if module == 'OCC':
                func = run_dfocc
            elif module in ['', 'FNOCC']:
                func = run_fnodfcc
        elif mtd_type == 'CD':
            if module == 'OCC':
                func = run_dfocc
            elif module in ['', 'FNOCC']:
                func = run_fnodfcc
    elif reference == 'UHF':
        if mtd_type == 'CONV':
            if module == 'MRCC' and which("dmrcc", return_bool=True):
                func = run_mrcc
            elif module in ['', 'CCENERGY']:
                func = run_ccenergy
        elif mtd_type in ["DF", "CD"]:
            if module in ["", "OCC"]:
                func = run_dfocc
    elif reference == 'ROHF':
        if mtd_type == 'CONV':
            if module == 'CCT3' and extras.addons("cct3"):
                import cct3
                func = cct3.run_cct3
            elif module == 'MRCC' and which("dmrcc", return_bool=True):
                func = run_mrcc
            elif module in ['', 'CCENERGY']:
                func = run_ccenergy

    if func is None:
        raise ManagedMethodError([__name__, name, type_var, mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_ccsd_gradient(name, **kwargs):
    """Function selecting the algorithm for a CCSD gradient call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    type_var, _, mtd_type = method_algorithm_type(name)
    module = core.get_global_option('QC_MODULE')

    func = None
    if reference == 'RHF':
        if mtd_type == 'CONV':
            if module in ['', 'CCENERGY']:
                func = run_ccenergy_gradient
        elif mtd_type == 'DF':
            if module in ['', 'OCC']:
                func = run_dfocc_gradient
    elif reference == 'UHF':
        if mtd_type == 'CONV':
            if module in ['', 'CCENERGY']:
                func = run_ccenergy_gradient
        elif mtd_type == 'DF':
            if module in ['', 'OCC']:
                func = run_dfocc_gradient
    elif reference == 'ROHF':
        if mtd_type == 'CONV':
            if module in ['', 'CCENERGY']:
                func = run_ccenergy_gradient

    if func is None:
        raise ManagedMethodError([__name__, name, type_var, mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_fnoccsd_t_(name, **kwargs):
    """Function selecting the algorithm for a FNO-CCSD(T) energy call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    type_var, _, mtd_type = method_algorithm_type(name)
    module = core.get_global_option('QC_MODULE')

    func = None
    if reference == 'RHF':
        if mtd_type == 'CONV':
            if module in ['', 'FNOCC']:
                func = run_fnocc
        elif mtd_type == 'DF':
            if module in ['', 'FNOCC']:
                func = run_fnodfcc
        elif mtd_type == 'CD':
            if module in ['', 'FNOCC']:
                func = run_fnodfcc

    if func is None:
        raise ManagedMethodError([__name__, name, type_var, mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_ccsd_t_(name, **kwargs):
    """Function selecting the algorithm for a CCSD(T) energy call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    type_var, _, mtd_type = method_algorithm_type(name)
    module = core.get_global_option('QC_MODULE')

    func = None
    if reference == 'RHF':
        if mtd_type == 'CONV':
            if module == 'FNOCC':
                func = run_fnocc
            elif module == 'MRCC' and which("dmrcc", return_bool=True):
                func = run_mrcc
            elif module in ['', 'CCENERGY']:
                func = run_ccenergy
        elif mtd_type == 'DF':
            if module == 'OCC':
                func = run_dfocc
            elif module in ['', 'FNOCC']:
                func = run_fnodfcc
        elif mtd_type == 'CD':
            if module == 'OCC':
                func = run_dfocc
            elif module in ['', 'FNOCC']:
                func = run_fnodfcc
    elif reference == 'UHF':
        if mtd_type == 'CONV':
            if module == 'MRCC' and which("dmrcc", return_bool=True):
                func = run_mrcc
            elif module in ['', 'CCENERGY']:
                func = run_ccenergy
        elif mtd_type in ["DF", "CD"]:
            if module in ["OCC"]:  # SOON "",
                func = run_dfocc
    elif reference == 'ROHF':
        if mtd_type == 'CONV':
            if module == 'MRCC' and which("dmrcc", return_bool=True):
                func = run_mrcc
            elif module in ['', 'CCENERGY']:
                func = run_ccenergy

    if func is None:
        raise ManagedMethodError([__name__, name, type_var, mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_ccsd_t__gradient(name, **kwargs):
    """Function selecting the algorithm for a CCSD(T) gradient call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    type_var, _, mtd_type = method_algorithm_type(name)
    module = core.get_global_option('QC_MODULE')

    func = None
    if reference in ['RHF']:
        if mtd_type == 'CONV':
            if module in ['', 'CCENERGY']:
                func = run_ccenergy_gradient
        elif mtd_type == 'DF':
            if module in ['', 'OCC']:
                func = run_dfocc_gradient
    elif reference == 'UHF':
        if mtd_type == 'CONV':
            if module in ['', 'CCENERGY']:
                func = run_ccenergy_gradient
        elif mtd_type == 'DF':
            if module in ['OCC']:  # SOON "",
                func = run_dfocc_gradient

    if func is None:
        raise ManagedMethodError([__name__, name, type_var, mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_ccsd_at_(name, **kwargs):
    """Function selecting the algorithm for a a-CCSD(T) energy call
    and directing to specified or best-performance default modules.

    """
    if name.lower() == "a-ccsd(t)":
        pass
    elif name.lower() in ["ccsd(at)", "lambda-ccsd(t)", "ccsd(t)_l"]:
        core.print_out(f"""\nMethod "{name.lower()}" has been regularized to "a-ccsd(t)" for QCVariables.""")
        name = "a-ccsd(t)"

    reference = core.get_option('SCF', 'REFERENCE')
    type_var, _, mtd_type = method_algorithm_type(name)
    module = core.get_global_option('QC_MODULE')

    func = None
    if reference == 'RHF':
        if mtd_type == 'CONV':
            if module == 'MRCC' and which("dmrcc", return_bool=True):
                func = run_mrcc
            elif module in ['', 'CCENERGY']:
                func = run_ccenergy
        elif mtd_type == 'DF':
            if module in ['', 'OCC']:
                func = run_dfocc
        elif mtd_type == 'CD':
            if module in ['', 'OCC']:
                func = run_dfocc
    elif reference == "UHF":
        if mtd_type == 'CONV':
            if module in ['', 'MRCC'] and which("dmrcc", return_bool=True):
                func = run_mrcc
        elif mtd_type == "DF":
            if module in ["OCC"]:  # SOON "",
                func = run_dfocc
        elif mtd_type == "CD":
            if module in ["OCC"]:  # SOON "",
                func = run_dfocc
    elif reference == "ROHF":
        if mtd_type == 'CONV':
            if module in ['', 'MRCC'] and which("dmrcc", return_bool=True):
                func = run_mrcc

    if func is None:
        raise ManagedMethodError([__name__, name, type_var, mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_cisd(name, **kwargs):
    """Function selecting the algorithm for a CISD energy call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    type_var, _, mtd_type = method_algorithm_type(name)
    module = core.get_global_option('QC_MODULE')

    func = None
    if reference == 'RHF':
        if mtd_type == 'CONV':
            if module == 'DETCI':
                func = run_detci
            elif module in ['', 'FNOCC']:
                func = run_cepa
    elif reference == 'ROHF':
        if mtd_type == 'CONV':
            if module in ['', 'DETCI']:
                func = run_detci

    if func is None:
        raise ManagedMethodError([__name__, name, type_var, mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_mp4(name, **kwargs):
    """Function selecting the algorithm for a MP4 energy call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    type_var, _, mtd_type = method_algorithm_type(name)
    module = core.get_global_option('QC_MODULE')

    func = None
    if reference == 'RHF':
        if mtd_type == 'CONV':
            if module == 'DETCI':
                func = run_detci
            elif module in ['', 'FNOCC']:
                func = run_fnocc
    elif reference == 'ROHF':
        if mtd_type == 'CONV':
            if module == 'DETCI':
                raise UpgradeHelper("energy('mp4')", "energy('zapt4')", 1.7,
                    " Replace method MP with method ZAPT for ROHF reference. DETCI is orders-of-magnitude inefficient for perturbation theory.")

    if module == 'DETCI':
        core.print_out("""\nDETCI is ill-advised for method MP4 as it is available inefficiently as a """
                       """byproduct of a CISDT computation.\n""")

    if func is None:
        raise ManagedMethodError([__name__, name, type_var, mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_remp2(name, **kwargs):
    """Function selecting the algorithm for a REMP2 energy call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    type_var, _, mtd_type = method_algorithm_type(name)
    module = core.get_global_option('QC_MODULE')

    func = None
    if reference in ["RHF", "UHF"]:
        if mtd_type == 'CONV':
            if module in ['', 'OCC']:
                func = run_occ
        elif mtd_type == 'DF':
            if module in ['', 'OCC']:
                func = run_dfocc
        elif mtd_type == 'CD':
            if module in ['', 'OCC']:
                func = run_dfocc

    if func is None:
        raise ManagedMethodError([__name__, name, type_var, mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_ccd(name, **kwargs):
    """Function selecting the algorithm for a CCD energy call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option("SCF", "REFERENCE")
    type_var, _, mtd_type = method_algorithm_type(name)
    module = core.get_global_option("QC_MODULE")

    func = None
    if reference in ["RHF", "UHF"]:
        if mtd_type == "CONV":
            if module in [""]:
                core.print_out("""\nThis method is not available with conventional integrals. Add "set """
                               """cc_type df" or "set cc_type cd" to input to access this method.\n""")
        elif mtd_type == "DF":
            if module in ["", "OCC"]:
                func = run_dfocc
        elif mtd_type == "CD":
            if module in ["", "OCC"]:
                func = run_dfocc

    if func is None:
        raise ManagedMethodError([__name__, name, type_var, mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_ccd_gradient(name, **kwargs):
    """Function selecting the algorithm for a CCD gradient call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option("SCF", "REFERENCE")
    type_var, _, mtd_type = method_algorithm_type(name)
    module = core.get_global_option("QC_MODULE")

    func = None
    if reference in ["RHF", "UHF"]:
        if mtd_type == "CONV":
            if module in [""]:
                core.print_out("""\nThis method is not available with conventional integrals. Add "set """
                               """cc_type df" or "set cc_type cd" to input to access this method.\n""")
        elif mtd_type == "DF":
            if module in ["", "OCC"]:
                func = run_dfocc_gradient

    if func is None:
        raise ManagedMethodError([__name__, name, type_var, mtd_type, reference, module])

    if kwargs.pop("probe", False):
        return
    else:
        return func(name, **kwargs)


def select_cc2(name, **kwargs):
    """Function selecting the algorithm for a CC2 energy call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    type_var, _, mtd_type = method_algorithm_type(name)
    module = core.get_global_option('QC_MODULE')

    # [LAB Aug 2022] I'm leaving MRCC CC2 in as a route, but my c.2014 MRCC consistently yields:
    #   "Approximate CC methods are not implemented for excitation level 2!"
    # [LAB Aug 2022] DF CC2 enabled for test_gradient but only by deliberate `set qc_module ccenergy`
    #   since not advertised. See #2710.

    func = None
    if reference == 'RHF':
        if mtd_type == 'CONV':
            if module == 'MRCC' and which("dmrcc", return_bool=True):
                func = run_mrcc
            elif module in ['', 'CCENERGY']:
                func = run_ccenergy
        elif mtd_type == 'DF':
            if module in ['CCENERGY']:
                func = run_ccenergy
    elif reference == 'UHF':
        if mtd_type == 'CONV':
            if module == 'MRCC' and which("dmrcc", return_bool=True):
                func = run_mrcc
            elif module in ['', 'CCENERGY']:
                func = run_ccenergy
    elif reference == 'ROHF':
        if mtd_type == 'CONV':
            if module == 'MRCC' and which("dmrcc", return_bool=True):
                func = run_mrcc
            elif module in ['', 'CCENERGY']:
                func = run_ccenergy

    if func is None:
        raise ManagedMethodError([__name__, name, type_var, mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_cc2_gradient(name, **kwargs):
    """Function selecting the algorithm for a CC2 gradient call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    type_var, _, mtd_type = method_algorithm_type(name)
    module = core.get_global_option('QC_MODULE')

    # [LAB Aug 2022] Both UHF and ROHF gradients run in ccenergy but ROHF is slightly off (1.e-5)
    #   and UHF is more off (1.e-4). Moreover, manual only claims RHF are working, so restricting here.
    # [LAB Aug 2022] DF CC2 enabled for test_gradient but only by deliberate `set qc_module ccenergy`
    #   since not advertised. See #2710.

    func = None
    if reference == 'RHF':
        if mtd_type == 'CONV':
            if module in ['', 'CCENERGY']:
                func = run_ccenergy_gradient
        elif mtd_type == 'DF':
            if module in ['CCENERGY']:
                func = run_ccenergy_gradient

    if func is None:
        raise ManagedMethodError([__name__, name, type_var, mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_cc3(name, **kwargs):
    """Function selecting the algorithm for a CC3 energy call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    type_var, _, mtd_type = method_algorithm_type(name)
    module = core.get_global_option('QC_MODULE')

    func = None
    if reference == 'RHF':
        if mtd_type == 'CONV':
            if module == 'MRCC' and which("dmrcc", return_bool=True):
                func = run_mrcc
            elif module in ['', 'CCENERGY']:
                func = run_ccenergy
    elif reference == 'UHF':
        if mtd_type == 'CONV':
            if module == 'MRCC' and which("dmrcc", return_bool=True):
                func = run_mrcc
            elif module in ['', 'CCENERGY']:
                func = run_ccenergy
    elif reference == 'ROHF':
        if mtd_type == 'CONV':
            # ROHF MRCC CC3  CCn methods are not implemented for ROHF reference!
            if module in ['', 'CCENERGY']:
                func = run_ccenergy

    if func is None:
        raise ManagedMethodError([__name__, name, type_var, mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_mrcc(name, **kwargs):
    """Function selecting the algorithm for a CC* energy call
    and directing to specified MRCC module.

    This function is unusual among "select" functions in that it services multiple methods and a
    single module. This function could have been skipped and the methods associated directly with
    run_rmcc; however, routing through this function screens for conv only while
    providing uniform error messages with other select functions.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    type_var, _, mtd_type = method_algorithm_type(name)
    module = core.get_global_option('QC_MODULE')

    # todo fix table link anchor

    func = None
    if reference == 'RHF':
        if mtd_type == 'CONV':
            if module in ['', 'MRCC'] and which("dmrcc", return_bool=True):
                func = run_mrcc
    elif reference == 'UHF':
        if mtd_type == 'CONV':
            if module in ['', 'MRCC'] and which("dmrcc", return_bool=True):
                func = run_mrcc
    elif reference == 'ROHF':
        if mtd_type == 'CONV':
            if module in ['', 'MRCC'] and which("dmrcc", return_bool=True):
                func = run_mrcc

    if func is None:
        raise ManagedMethodError([__name__, name, type_var, mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def build_disp_functor(name, restricted, save_pairwise_disp=False, **kwargs):

    if core.has_option_changed("SCF", "DFT_DISPERSION_PARAMETERS"):
        modified_disp_params = core.get_option("SCF", "DFT_DISPERSION_PARAMETERS")
    else:
        modified_disp_params = None

    # Figure out functional
    superfunc, disp_type = dft.build_superfunctional(name, restricted)

    if disp_type:
        if isinstance(name, dict):
            # user dft_functional={} spec - type for lookup, dict val for param defs,
            #   name & citation discarded so only param matches to existing defs will print labels
            _disp_functor = empirical_dispersion.EmpiricalDispersion(name_hint='',
                                                                     level_hint=disp_type["type"],
                                                                     param_tweaks=disp_type["params"],
                                                                     save_pairwise_disp=save_pairwise_disp,
                                                                     engine=kwargs.get('engine', None))
        else:
            # dft/*functionals.py spec - name & type for lookup, option val for param tweaks
            _disp_functor = empirical_dispersion.EmpiricalDispersion(name_hint=superfunc.name(),
                                                                     level_hint=disp_type["type"],
                                                                     param_tweaks=modified_disp_params,
                                                                     save_pairwise_disp=save_pairwise_disp,
                                                                     engine=kwargs.get('engine', None))

        # [Aug 2018] there once was a breed of `disp_type` that quacked
        #   like a list rather than the more common dict handled above. if
        #   ever again sighted, make an issue so this code can accommodate.

        _disp_functor.print_out()
        return superfunc, _disp_functor

    else:
        return superfunc, None

def scf_wavefunction_factory(name, ref_wfn, reference, **kwargs):
    """Builds the correct (R/U/RO/CU HF/KS) wavefunction from the
    provided information, sets relevant auxiliary basis sets on it,
    and prepares any empirical dispersion.

    """
    # Figure out functional and dispersion
    superfunc, _disp_functor = build_disp_functor(name, restricted=(reference in ["RKS", "RHF"]), **kwargs)

    # Build the wavefunction
    core.prepare_options_for_module("SCF")
    if reference in ["RHF", "RKS"]:
        wfn = core.RHF(ref_wfn, superfunc)
    elif reference == "ROHF":
        wfn = core.ROHF(ref_wfn, superfunc)
    elif reference in ["UHF", "UKS"]:
        wfn = core.UHF(ref_wfn, superfunc)
    elif reference == "CUHF":
        wfn = core.CUHF(ref_wfn, superfunc)
    else:
        raise ValidationError("SCF: Unknown reference (%s) when building the Wavefunction." % reference)

    if _disp_functor and _disp_functor.engine != 'nl':
        wfn._disp_functor = _disp_functor

    # Set the DF basis sets
    df_needed = core.get_global_option("SCF_TYPE") in ["DF", "MEM_DF", "DISK_DF", "COSX", "LINK"]
    df_needed |= (core.get_global_option("SCF_TYPE") == "DIRECT" and core.get_option("SCF", "DF_SCF_GUESS"))
    if df_needed:
        aux_basis = core.BasisSet.build(wfn.molecule(), "DF_BASIS_SCF",
                                        core.get_option("SCF", "DF_BASIS_SCF"),
                                        "JKFIT", core.get_global_option('BASIS'),
                                        puream=wfn.basisset().has_puream())
        wfn.set_basisset("DF_BASIS_SCF", aux_basis)
    else:
        wfn.set_basisset("DF_BASIS_SCF", core.BasisSet.zero_ao_basis_set())

    # Set the relativistic basis sets
    if core.get_global_option("RELATIVISTIC") in ["X2C", "DKH"]:
        decon_basis = core.BasisSet.build(wfn.molecule(), "BASIS_RELATIVISTIC",
                                        core.get_option("SCF", "BASIS_RELATIVISTIC"),
                                        "DECON", core.get_global_option('BASIS'),
                                        puream=wfn.basisset().has_puream())
        wfn.set_basisset("BASIS_RELATIVISTIC", decon_basis)

    # Set the multitude of SAD basis sets
    if (core.get_option("SCF", "GUESS") in ["SAD", "SADNO", "HUCKEL"]):
        sad_basis_list = core.BasisSet.build(wfn.molecule(), "ORBITAL",
                                             core.get_global_option("BASIS"),
                                             puream=wfn.basisset().has_puream(),
                                             return_atomlist=True)
        wfn.set_sad_basissets(sad_basis_list)

        if ("DF" in core.get_option("SCF", "SAD_SCF_TYPE")):
            # We need to force this to spherical regardless of any user or other demands.
            optstash = p4util.OptionsState(['PUREAM'])
            core.set_global_option('PUREAM', True)
            sad_fitting_list = core.BasisSet.build(wfn.molecule(), "DF_BASIS_SAD",
                                                   core.get_option("SCF", "DF_BASIS_SAD"),
                                                   puream=True,
                                                   return_atomlist=True)
            wfn.set_sad_fitting_basissets(sad_fitting_list)
            optstash.restore()

    if hasattr(core, "EXTERN") and 'external_potentials' in kwargs:
        core.print_out("\n  Warning! Both an external potential EXTERN object and the external_potential" +
                       " keyword argument are specified. The external_potentials keyword argument will be ignored.\n")
        raise ValidationError("double extern")

    ep = kwargs.get("external_potentials", None)
    if ep is not None:
        _set_external_potentials_to_wavefunction(ep, wfn)

    return wfn


def _set_external_potentials_to_wavefunction(external_potential: Union[List, Dict[str, List]], wfn: "core.Wavefunction"):
    """Initialize :py:class:`psi4.core.ExternalPotential` object(s) from charges and locations and set on **wfn**.

    Parameters
    ----------
    external_potential
        List-like structure where each row corresponds to a charge. Lines can be composed of ``q, [x, y, z]`` or
        ``q, x, y, z``. Locations are in [a0].
        Or, dictionary where keys are FI-SAPT fragments A, B, or C and values are as above.

    """
    from psi4.driver.qmmm import QMMMbohr

    def validate_qxyz(qxyz):
        if len(qxyz) == 2:
            return qxyz[0], qxyz[1][0], qxyz[1][1], qxyz[1][2]
        elif len(qxyz) == 4:
            return qxyz[0], qxyz[1], qxyz[2], qxyz[3]
        else:
            raise ValidationError(f"Point charge '{qxyz}' not mapping into 'chg, [x, y, z]' or 'chg, x, y, z'")

    if isinstance(external_potential, dict):
        # For FSAPT, we can take a dictionary of external potentials, e.g.,
        # external_potentials={'A': potA, 'B': potB, 'C': potC} (any optional)
        # For the dimer SAPT calculation, we need to account for the external potential
        # in all of the subsystems A, B, C. So we add them all in total_external_potential
        # and set the external potential to the dimer wave function

        total_external_potential = core.ExternalPotential()

        for frag, frag_qxyz in external_potential.items():
            if frag.upper() in "ABC":
                chrgfield = QMMMbohr()
                for qxyz in frag_qxyz:
                    chrgfield.extern.addCharge(*validate_qxyz(qxyz))

                wfn.set_potential_variable(frag.upper(), chrgfield.extern)
                total_external_potential.appendCharges(chrgfield.extern.getCharges())

            else:
                core.print_out("\n  Warning! Unknown key for the external_potentials argument: %s" % frag)

        wfn.set_external_potential(total_external_potential)

    else:
        chrgfield = QMMMbohr()
        for qxyz in external_potential:
            chrgfield.extern.addCharge(*validate_qxyz(qxyz))
        wfn.set_potential_variable("C", chrgfield.extern)  # This is for the FSAPT procedure
        wfn.set_external_potential(chrgfield.extern)


def scf_helper(name, post_scf=True, **kwargs):
    """Function serving as helper to SCF, choosing whether to cast
    up or just run SCF with a standard guess. This preserves
    previous SCF options set by other procedures (e.g., SAPT
    output file types for SCF). Most run_* functions should call
    this function, common exceptions being when multireference
    SCF is needed or when restarting from converged SCF.

    """

    if post_scf:
        name = "scf"

    optstash = p4util.OptionsState(
        ['PUREAM'],
        ['BASIS'],
        ['QMEFP'],
        ['INTS_TOLERANCE'],
        ['DF_BASIS_SCF'],
        ['SCF', 'GUESS'],
        ['SCF', 'DF_INTS_IO'],
        ['SCF', 'ORBITALS_WRITE'],
        ['SCF_TYPE'],  # Hack: scope gets changed internally with the Andy trick
    )

    optstash2 = p4util.OptionsState(
        ['BASIS'],
        ['DF_BASIS_SCF'],
        ['SCF_TYPE'],
        ['SCF', 'DF_INTS_IO'],
    )

    # Make sure we grab the correctly scoped integral threshold for SCF
    core.set_global_option('INTS_TOLERANCE', core.get_option('SCF', 'INTS_TOLERANCE'))

    # Grab a few kwargs
    use_c1 = kwargs.get('use_c1', False)
    scf_molecule = kwargs.get('molecule', core.get_active_molecule())
    read_orbitals = core.get_option('SCF', 'GUESS') == "READ"
    do_timer = kwargs.pop("do_timer", True)
    ref_wfn = kwargs.pop('ref_wfn', None)
    if ref_wfn is not None:
        raise ValidationError("Cannot seed an SCF calculation with a reference wavefunction ('ref_wfn' kwarg).")

    # decide if we keep the checkpoint file
    _chkfile = kwargs.get('write_orbitals', True)
    write_checkpoint_file = False
    if isinstance(_chkfile, str):
        write_checkpoint_file = True
        filename = kwargs.get('write_orbitals')
        core.set_local_option("SCF", "ORBITALS_WRITE", filename)
    elif _chkfile is True:
        write_checkpoint_file = True

    # Continuum solvation needs to be run w/o symmetry
    if core.get_option("SCF", "PCM") or core.get_option("SCF", "DDX"):
        c1_molecule = scf_molecule.clone()
        c1_molecule.reset_point_group('c1')
        c1_molecule.update_geometry()

        scf_molecule = c1_molecule
        core.print_out("""  PCM or DDX continuum solvation does not make use of molecular symmetry: """
                       """further calculations in C1 point group.\n""")

    # PE needs to use exactly input orientation to correspond to potfile
    if core.get_option("SCF", "PE"):
        c1_molecule = scf_molecule.clone()
        if getattr(scf_molecule, "_initial_cartesian", None) is not None:
            c1_molecule._initial_cartesian = scf_molecule._initial_cartesian.clone()
            c1_molecule.set_geometry(c1_molecule._initial_cartesian)
            c1_molecule.reset_point_group("c1")
            c1_molecule.fix_orientation(True)
            c1_molecule.fix_com(True)
            c1_molecule.update_geometry()
        else:
            raise ValidationError("Set no_com/no_reorient/symmetry c1 by hand for PE on non-Cartesian molecules.")

        scf_molecule = c1_molecule
        core.print_out("""  PE does not make use of molecular symmetry: """
                       """further calculations in C1 point group.\n""")
        core.print_out("""  PE geometry must align with POTFILE keyword: """
                       """resetting coordinates with fixed origin and orientation.\n""")

    # SCF Banner data
    banner = kwargs.pop('banner', None)
    bannername = name

    # Did we pass in a DFT functional?
    dft_func = kwargs.pop('dft_functional', None)
    if dft_func is not None:
        if name.lower() != "scf":
            raise ValidationError("dft_functional was supplied to SCF, but method name was not SCF ('%s')" % name)
        name = dft_func
        bannername = name
        if isinstance(name, dict):
            bannername = name.get("name", "custom functional")


    # Setup the timer
    if do_timer:
        core.tstart()

    # Second-order SCF requires non-symmetric density matrix support
    if core.get_option('SCF', 'SOSCF'):
        proc_util.check_non_symmetric_jk_density("Second-order SCF")

    # sort out cast_up settings. no need to stash these since only read, never reset
    cast = False
    if core.has_option_changed('SCF', 'BASIS_GUESS'):
        cast = core.get_option('SCF', 'BASIS_GUESS')
        if p4util.yes.match(str(cast)):
            cast = True
        elif p4util.no.match(str(cast)):
            cast = False

    if cast:

        # A user can set "BASIS_GUESS" to True and we default to 3-21G
        if cast is True:
            guessbasis = corresponding_basis(core.get_global_option('BASIS'), 'GUESS')[0]
            if guessbasis is None:
                guessbasis = '3-21G'  # guess of last resort
        else:
            guessbasis = cast
        core.set_global_option('BASIS', guessbasis)

        castdf = 'DF' in core.get_global_option('SCF_TYPE')

        if core.has_option_changed('SCF', 'DF_BASIS_GUESS'):
            castdf = core.get_option('SCF', 'DF_BASIS_GUESS')
            if p4util.yes.match(str(castdf)):
                castdf = True
            elif p4util.no.match(str(castdf)):
                castdf = False

        if castdf:
            core.set_global_option('SCF_TYPE', 'DF')
            core.set_local_option('SCF', 'DF_INTS_IO', 'none')

            # Figure out the fitting basis set
            if castdf is True:
                core.set_global_option('DF_BASIS_SCF', '')
            elif isinstance(castdf, str):
                core.set_global_option('DF_BASIS_SCF', castdf)
            else:
                raise ValidationError("Unexpected castdf option (%s)." % castdf)


        # Switch to the guess namespace
        namespace = core.IO.get_default_namespace()
        guesspace = namespace + '.guess'
        if namespace == '':
            guesspace = 'guess'
        core.IO.set_default_namespace(guesspace)

        # Print some info about the guess
        core.print_out('\n')
        p4util.banner('Guess SCF, %s Basis' % (guessbasis))
        core.print_out('\n')

    if cast and read_orbitals:
        raise ValidationError("""Detected options to both cast and read orbitals""")

    if (core.get_option('SCF', 'STABILITY_ANALYSIS') == 'FOLLOW') and (core.get_option('SCF', 'REFERENCE') != 'UHF'):
        raise ValidationError(f"""Stability analysis root following is only available for unrestricted Hartree--Fock, not present {core.get_option('SCF', 'REFERENCE')}""")

    # If GUESS is auto guess what it should be
    if core.get_option('SCF', 'GUESS') == "AUTO":
        if (scf_molecule.natom() > 1):
            core.set_local_option('SCF', 'GUESS', 'SAD')
        else:
            core.set_local_option('SCF', 'GUESS', 'CORE')

    if core.get_global_option('BASIS') in ['', '(AUTO)']:
        if name in ['hf3c', 'hf-3c']:
            core.set_global_option('BASIS', 'minix')
        elif name in ['pbeh3c', 'pbeh-3c']:
            core.set_global_option('BASIS', 'def2-msvp')

    # the FIRST scf call
    if cast:
        # Cast is a special case
        base_wfn = core.Wavefunction.build(scf_molecule, core.get_global_option('BASIS'))
        core.print_out("\n         ---------------------------------------------------------\n")
        if banner:
            core.print_out("         " + banner.center(58))
        if cast:
            core.print_out("         " + "SCF Castup computation".center(58))
        ref_wfn = scf_wavefunction_factory(name, base_wfn, core.get_option('SCF', 'REFERENCE'), **kwargs)

        # Compute additive correction: dftd3, mp2d, dftd4, etc.
        if hasattr(ref_wfn, "_disp_functor"):
            disp_energy = ref_wfn._disp_functor.compute_energy(ref_wfn.molecule())
            ref_wfn.set_variable("-D Energy", disp_energy)
        ref_wfn.compute_energy()

    # cast clean-up
    if cast:

        # Move files to proper namespace
        core.IO.change_file_namespace(180, guesspace, namespace)
        core.IO.set_default_namespace(namespace)

        optstash2.restore()

        # Print the banner for the standard operation
        core.print_out('\n')
        p4util.banner(bannername.upper())
        core.print_out('\n')

    # the SECOND scf call
    base_wfn = core.Wavefunction.build(scf_molecule, core.get_global_option('BASIS'))
    if banner:
        core.print_out("\n         ---------------------------------------------------------\n")
        core.print_out("         " + banner.center(58))

    scf_wfn = scf_wavefunction_factory(name, base_wfn, core.get_option('SCF', 'REFERENCE'), **kwargs)

    # The wfn from_file routine adds the npy suffix if needed, but we add it here so that
    # we can use os.path.isfile to query whether the file exists before attempting to read
    read_filename = scf_wfn.get_scratch_filename(180) + '.npy'
    if ((core.get_option('SCF', 'GUESS') == 'READ') and os.path.isfile(read_filename)):
        old_wfn = core.Wavefunction.from_file(read_filename)

        Ca_occ = old_wfn.Ca_subset("SO", "OCC")
        Cb_occ = old_wfn.Cb_subset("SO", "OCC")

        if old_wfn.molecule().schoenflies_symbol() != scf_molecule.schoenflies_symbol():
            raise ValidationError("Cannot compute projection of different symmetries.")

        if old_wfn.basisset().name() == scf_wfn.basisset().name():
            core.print_out(f"  Reading orbitals from file {read_filename}, no projection.\n\n")
            scf_wfn.guess_Ca(Ca_occ)
            scf_wfn.guess_Cb(Cb_occ)
        else:
            core.print_out(f"  Reading orbitals from file {read_filename}, projecting to new basis.\n\n")
            core.print_out("  Computing basis projection from %s to %s\n\n" % (old_wfn.basisset().name(), scf_wfn.basisset().name()))

            pCa = scf_wfn.basis_projection(Ca_occ, old_wfn.nalphapi(), old_wfn.basisset(), scf_wfn.basisset())
            pCb = scf_wfn.basis_projection(Cb_occ, old_wfn.nbetapi(), old_wfn.basisset(), scf_wfn.basisset())
            scf_wfn.guess_Ca(pCa)
            scf_wfn.guess_Cb(pCb)

        # Strip off headers to only get R, RO, U, CU
        old_ref = old_wfn.name().replace("KS", "").replace("HF", "")
        new_ref = scf_wfn.name().replace("KS", "").replace("HF", "")
        if old_ref != new_ref:
            scf_wfn.reset_occ_ = True

    elif (core.get_option('SCF', 'GUESS') == 'READ') and not os.path.isfile(read_filename):
        core.print_out(f"\n !!!  Unable to find file {read_filename}, defaulting to SAD guess. !!!\n\n")
        core.set_local_option('SCF', 'GUESS', 'SAD')
        sad_basis_list = core.BasisSet.build(scf_wfn.molecule(), "ORBITAL",
                                             core.get_global_option("BASIS"),
                                             puream=scf_wfn.basisset().has_puream(),
                                             return_atomlist=True)
        scf_wfn.set_sad_basissets(sad_basis_list)

        if ("DF" in core.get_option("SCF", "SAD_SCF_TYPE")):
            sad_fitting_list = core.BasisSet.build(scf_wfn.molecule(), "DF_BASIS_SAD",
                                                   core.get_option("SCF", "DF_BASIS_SAD"),
                                                   puream=scf_wfn.basisset().has_puream(),
                                                   return_atomlist=True)
            scf_wfn.set_sad_fitting_basissets(sad_fitting_list)


    if cast:
        core.print_out("\n  Computing basis projection from %s to %s\n\n" % (ref_wfn.basisset().name(), base_wfn.basisset().name()))
        if ref_wfn.basisset().n_ecp_core() != base_wfn.basisset().n_ecp_core():
            raise ValidationError("Projecting from basis ({}) with ({}) ECP electrons to basis ({}) with ({}) ECP electrons will be a disaster. Select a compatible cast-up basis with `set guess_basis YOUR_BASIS_HERE`.".format(
                                  ref_wfn.basisset().name(), ref_wfn.basisset().n_ecp_core(), base_wfn.basisset().name(), base_wfn.basisset().n_ecp_core()))

        pCa = ref_wfn.basis_projection(ref_wfn.Ca(), ref_wfn.nalphapi(), ref_wfn.basisset(), scf_wfn.basisset())
        pCb = ref_wfn.basis_projection(ref_wfn.Cb(), ref_wfn.nbetapi(), ref_wfn.basisset(), scf_wfn.basisset())
        scf_wfn.guess_Ca(pCa)
        scf_wfn.guess_Cb(pCb)


    # Print basis set info
    if core.get_option("SCF", "PRINT_BASIS"):
        scf_wfn.basisset().print_detail_out()

    # Compute additive correction: dftd3, mp2d, dftd4, etc.
    if hasattr(scf_wfn, "_disp_functor"):
        disp_energy = scf_wfn._disp_functor.compute_energy(scf_wfn.molecule(), scf_wfn)
        scf_wfn.set_variable("-D Energy", disp_energy)

    # PCM preparation
    if core.get_option('SCF', 'PCM'):
        if core.get_option('SCF', 'PE'):
            raise ValidationError("""Error: 3-layer QM/MM/PCM not implemented.\n""")
        pcmsolver_parsed_fname = core.get_local_option('PCM', 'PCMSOLVER_PARSED_FNAME')
        pcm_print_level = core.get_option('SCF', "PRINT")
        scf_wfn.set_PCM(core.PCM(pcmsolver_parsed_fname, pcm_print_level, scf_wfn.basisset()))

    # DDPCM preparation
    if core.get_option('SCF', 'DDX'):
        if not solvent._have_ddx:
            raise ModuleNotFoundError('Python module ddx not found. Solve by installing it: `pip install pyddx`')
        ddx_options = solvent.ddx.get_ddx_options(scf_molecule)
        scf_wfn.ddx_state = solvent.ddx.DdxInterface(
            molecule=scf_molecule, options=ddx_options,
            basisset=scf_wfn.basisset()
        )

    # PE preparation
    if core.get_option('SCF', 'PE'):
        if not solvent._have_pe:
            raise ModuleNotFoundError('Python module cppe not found. Solve by installing it: `conda install -c psi4 pycppe`')
        # PE needs information about molecule and basis set
        pol_embed_options = solvent.pol_embed.get_pe_options()
        core.print_out(f""" Using potential file
                       {pol_embed_options["potfile"]}
                       for Polarizable Embedding calculation.\n""")
        scf_wfn.pe_state = solvent.pol_embed.CppeInterface(
            molecule=scf_molecule, options=pol_embed_options,
            basisset=scf_wfn.basisset()
        )

    e_scf = scf_wfn.compute_energy()
    for obj in [core, scf_wfn]:
        # set_variable("SCF TOTAL ENERGY")  # P::e SCF
        for pv in ["SCF TOTAL ENERGY", "CURRENT ENERGY", "CURRENT REFERENCE ENERGY"]:
            obj.set_variable(pv, e_scf)

    # We always would like to print a little property information
    if kwargs.get('scf_do_properties', True):
        oeprop = core.OEProp(scf_wfn)
        oeprop.set_title("SCF")

        # Figure our properties, if empty do dipole
        props = [x.upper() for x in core.get_option("SCF", "SCF_PROPERTIES")]
        if "DIPOLE" not in props:
            props.append("DIPOLE")

        proc_util.oeprop_validator(props)
        for x in props:
            oeprop.add(x)

        # Populate free-atom volumes
        # if we're doing MBIS
        if 'MBIS_VOLUME_RATIOS' in props:
            p4util.free_atom_volumes(scf_wfn)

        # Compute properties
        oeprop.compute()
        for obj in [core, scf_wfn]:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
            obj.set_variable("CURRENT DIPOLE", obj.variable("SCF DIPOLE"))  # P::e SCF

    # Write out MO's
    if core.get_option("SCF", "PRINT_MOS"):
        mowriter = core.MOWriter(scf_wfn)
        mowriter.write()

    # Write out a molden file
    if core.get_option("SCF", "MOLDEN_WRITE"):
        filename = core.get_writer_file_prefix(scf_molecule.name()) + ".molden"
        dovirt = bool(core.get_option("SCF", "MOLDEN_WITH_VIRTUAL"))

        occa = scf_wfn.occupation_a()
        occb = scf_wfn.occupation_a()

        mw = core.MoldenWriter(scf_wfn)
        mw.write(filename, scf_wfn.Ca(), scf_wfn.Cb(), scf_wfn.epsilon_a(),
                 scf_wfn.epsilon_b(), scf_wfn.occupation_a(),
                 scf_wfn.occupation_b(), dovirt)

    # Write checkpoint file (orbitals and basis); Can be disabled, e.g., for findif displacements
    if write_checkpoint_file and isinstance(_chkfile, str):
        filename = kwargs['write_orbitals']
        scf_wfn.to_file(filename)
        # core.set_local_option("SCF", "ORBITALS_WRITE", filename)
    elif write_checkpoint_file:
        filename = scf_wfn.get_scratch_filename(180)
        scf_wfn.to_file(filename)
        extras.register_numpy_file(filename) # retain with -m (messy) option

    if do_timer:
        core.tstop()

    optstash.restore()

    if (not use_c1) or (scf_molecule.schoenflies_symbol() == 'c1'):
        return scf_wfn
    else:
        # C1 copy quietly
        c1_optstash = p4util.OptionsState(['PRINT'])
        core.set_global_option("PRINT", 0)

        # If we force c1 copy the active molecule
        scf_molecule.update_geometry()
        core.print_out("""\n  A requested method does not make use of molecular symmetry: """
                           """further calculations in C1 point group.\n\n""")
        c1_molecule = scf_molecule.clone()
        c1_molecule.reset_point_group('c1')
        c1_molecule.fix_orientation(True)
        c1_molecule.fix_com(True)
        c1_molecule.update_geometry()
        c1_basis = core.BasisSet.build(c1_molecule, "ORBITAL", core.get_global_option('BASIS'), quiet=True)
        tmp = scf_wfn.c1_deep_copy(c1_basis)
        if not scf_wfn.has_variable("-D ENERGY"):
            tmp.del_variable("-D ENERGY")
        c1_jkbasis = core.BasisSet.build(c1_molecule, "DF_BASIS_SCF",
                                         core.get_global_option("DF_BASIS_SCF"),
                                         "JKFIT", core.get_global_option('BASIS'), quiet=True)
        tmp.set_basisset("DF_BASIS_SCF", c1_jkbasis)
        c1_optstash.restore()
        return tmp


def run_dct(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a density cumulant theory calculation.

    """

    if (core.get_global_option('FREEZE_CORE') == 'TRUE'):
        raise ValidationError('Frozen core is not available for DCT.')

    # Bypass the scf call if a reference wavefunction is given
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, **kwargs)

    if (core.get_global_option("DCT_TYPE") == "DF"):
        core.print_out("  Constructing Basis Sets for DCT...\n\n")
        aux_basis = core.BasisSet.build(ref_wfn.molecule(), "DF_BASIS_DCT",
                                        core.get_global_option("DF_BASIS_DCT"),
                                        "RIFIT", core.get_global_option("BASIS"))
        ref_wfn.set_basisset("DF_BASIS_DCT", aux_basis)

        scf_aux_basis = core.BasisSet.build(ref_wfn.molecule(), "DF_BASIS_SCF",
                                            core.get_option("SCF", "DF_BASIS_SCF"),
                                            "JKFIT", core.get_global_option('BASIS'),
                                            puream=ref_wfn.basisset().has_puream())
        ref_wfn.set_basisset("DF_BASIS_SCF", scf_aux_basis)
        dct_wfn = core.dct(ref_wfn)

    else:
        # Ensure IWL files have been written for non DF-DCT
        proc_util.check_iwl_file_from_scf_type(core.get_global_option('SCF_TYPE'), ref_wfn)
        dct_wfn = core.dct(ref_wfn)

    for k, v in dct_wfn.variables().items():
        core.set_variable(k, v)

    return dct_wfn


def run_dct_gradient(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    DCT gradient calculation.

    """
    optstash = p4util.OptionsState(
        ['GLOBALS', 'DERTYPE'])


    core.set_global_option('DERTYPE', 'FIRST')
    dct_wfn = run_dct_property(name, **kwargs)

    derivobj = core.Deriv(dct_wfn)
    derivobj.set_tpdm_presorted(True)
    if core.get_option('DCT', 'DCT_TYPE') == 'CONV':
        grad = derivobj.compute()
    else:
        grad = derivobj.compute_df('DF_BASIS_SCF', 'DF_BASIS_DCT')

    dct_wfn.set_gradient(grad)

    optstash.restore()
    return dct_wfn


def run_dct_property(name, **kwargs):
    """ Function encoding sequence of PSI module calls for
    DCT property calculation.

    """
    optstash = p4util.OptionsState(
        ['DCT', 'OPDM'])

    core.set_local_option('DCT', 'OPDM', 'true')
    dct_wfn = run_dct(name, **kwargs)

    # Run OEProp
    oe = core.OEProp(dct_wfn)
    oe.set_title("DCT")
    for prop in kwargs.get("properties", []):
        prop = prop.upper()
        if prop in core.OEProp.valid_methods or "MULTIPOLE(" in prop:
            oe.add(prop)
    oe.compute()
    dct_wfn.oeprop = oe

    for k, v in dct_wfn.variables().items():
        core.set_variable(k, v)

    optstash.restore()
    return dct_wfn


def run_dfocc(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a density-fitted or Cholesky-decomposed
    (non-)orbital-optimized MPN or CC computation.

    """
    dtl = docs_table_link("dummy", "occ_nonoo")

    optstash = p4util.OptionsState(
        ['SCF', 'DF_INTS_IO'],
        ['DFOCC', 'WFN_TYPE'],
        ['DFOCC', 'ORB_OPT'],
        ['DFOCC', 'DO_SCS'],
        ['DFOCC', 'DO_SOS'],
        ['DFOCC', 'READ_SCF_3INDEX'],
        ['DFOCC', 'CHOLESKY'],
        ['DFOCC', 'CC_LAMBDA'])

    def set_cholesky_from(corl_type):
        if corl_type == 'DF':
            core.set_local_option('DFOCC', 'CHOLESKY', 'FALSE')
            proc_util.check_disk_df(name.upper(), optstash)

        elif corl_type == 'CD':
            core.set_local_option('DFOCC', 'CHOLESKY', 'TRUE')
            # Alter default algorithm
            if not core.has_global_option_changed('SCF_TYPE'):
                optstash.add_option(['SCF_TYPE'])
                core.set_global_option('SCF_TYPE', 'CD')
                core.print_out("""    SCF Algorithm Type (re)set to CD.\n""")
            if core.get_global_option('SCF_TYPE') != 'CD':
                core.set_local_option('DFOCC', 'READ_SCF_3INDEX', 'FALSE')
        else:
            raise ValidationError(f"""Invalid type '{corl_type}' for DFOCC. See Capabilities Table at {dtl}""")

    if name in ["mp2.5", "mp3"] and not core.has_global_option_changed("MP_TYPE"):
        core.print_out(f"    Information: {name.upper()} default algorithm changed to DF in August 2020. Use `set mp_type conv` for previous behavior.\n")

    director = {
           "mp2":     {"wfn_type": "DF-OMP2",     "orb_opt": "FALSE", "nat_orbs": "FALSE",},
          "omp2":     {"wfn_type": "DF-OMP2",     "orb_opt": "TRUE",  "nat_orbs": "FALSE",},

           "mp2.5":   {"wfn_type": "DF-OMP2.5",   "orb_opt": "FALSE", "nat_orbs": "FALSE",},
          "omp2.5":   {"wfn_type": "DF-OMP2.5",   "orb_opt": "TRUE",  "nat_orbs": "FALSE",},

           "mp3":     {"wfn_type": "DF-OMP3",     "orb_opt": "FALSE", "nat_orbs": "FALSE",},
          "omp3":     {"wfn_type": "DF-OMP3",     "orb_opt": "TRUE",  "nat_orbs": "FALSE",},

         "remp2":     {"wfn_type": "DF-OREMP",    "orb_opt": "FALSE", "nat_orbs": "FALSE",},
        "oremp2":     {"wfn_type": "DF-OREMP",    "orb_opt": "TRUE",  "nat_orbs": "FALSE",},

          "lccd":     {"wfn_type": "DF-OLCCD",    "orb_opt": "FALSE", "nat_orbs": "FALSE",},
         "olccd":     {"wfn_type": "DF-OLCCD",    "orb_opt": "TRUE",  "nat_orbs": "FALSE",},

           "ccd":     {"wfn_type": "DF-CCD",      "orb_opt": "FALSE", "nat_orbs": "FALSE",},  # changes to DF-OCCD

           "ccsd":    {"wfn_type": "DF-CCSD",     "orb_opt": "FALSE", "nat_orbs": "FALSE",},

           "ccsd(t)": {"wfn_type": "DF-CCSD(T)",  "orb_opt": "FALSE", "nat_orbs": "FALSE",},

         "a-ccsd(t)": {"wfn_type": "DF-CCSD(AT)", "orb_opt": "FALSE", "nat_orbs": "FALSE",},
    }

    if name not in director:
        raise ValidationError(f"Invalid method {name} for DFOCC energy")

    # throw exception for CONV (approximately). run reference defaulting logic
    set_cholesky_from(method_algorithm_type(name).now)

    for k, v in director[name].items():
        core.set_local_option("DFOCC", k.upper(), v)

    core.set_local_option('DFOCC', 'DO_SCS', 'FALSE')
    core.set_local_option('DFOCC', 'DO_SOS', 'FALSE')
    core.set_local_option('SCF', 'DF_INTS_IO', 'SAVE')

    # Bypass the scf call if a reference wavefunction is given
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, use_c1=True, **kwargs)  # C1 certified
    else:
        if ref_wfn.molecule().schoenflies_symbol() != 'c1':
            raise ValidationError("""  DFOCC does not make use of molecular symmetry: """
                                  """reference wavefunction must be C1.\n""")

    if not core.get_local_option("DFOCC", "CHOLESKY"):
        core.print_out("  Constructing Basis Sets for DFOCC...\n\n")
        scf_aux_basis = core.BasisSet.build(ref_wfn.molecule(), "DF_BASIS_SCF",
                                           core.get_option("SCF", "DF_BASIS_SCF"),
                                           "JKFIT", core.get_global_option('BASIS'),
                                           puream=ref_wfn.basisset().has_puream())

        ref_wfn.set_basisset("DF_BASIS_SCF", scf_aux_basis)

        aux_basis = core.BasisSet.build(ref_wfn.molecule(), "DF_BASIS_CC",
                                            core.get_global_option("DF_BASIS_CC"),
                                            "RIFIT", core.get_global_option("BASIS"))
        ref_wfn.set_basisset("DF_BASIS_CC", aux_basis)

    if core.get_option('SCF', 'REFERENCE') == 'ROHF':
        ref_wfn.semicanonicalize()

    dfocc_wfn = core.dfocc(ref_wfn)

    # Shove variables into global space
    for k, v in dfocc_wfn.variables().items():
        core.set_variable(k, v)

    optstash.restore()
    return dfocc_wfn


def run_dfocc_gradient(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a density-fitted (non-)orbital-optimized MPN or CC computation.

    """
    dtl = docs_table_link("dummy", "occ_nonoo")

    optstash = p4util.OptionsState(
        ['SCF', 'DF_INTS_IO'],
        ['REFERENCE'],
        ['DFOCC', 'WFN_TYPE'],
        ['DFOCC', 'ORB_OPT'],
        ['DFOCC', 'CC_LAMBDA'],
        ['GLOBALS', 'DERTYPE'])

    if name in ["mp2.5", "mp3"] and not core.has_global_option_changed("MP_TYPE"):
        core.print_out(f"    Information: {name.upper()} default algorithm changed to DF in August 2020. Use `set mp_type conv` for previous behavior.\n")

    # CC_LAMBDA keyword was being set TRUE sporadically, but that's covered c-side

    director = {
           "mp2":     {"wfn_type": "DF-OMP2",     "orb_opt": "FALSE", "nat_orbs": "FALSE",},
          "omp2":     {"wfn_type": "DF-OMP2",     "orb_opt": "TRUE",  "nat_orbs": "FALSE",},

           "mp2.5":   {"wfn_type": "DF-OMP2.5",   "orb_opt": "FALSE", "nat_orbs": "FALSE",},
          "omp2.5":   {"wfn_type": "DF-OMP2.5",   "orb_opt": "TRUE",  "nat_orbs": "FALSE",},

           "mp3":     {"wfn_type": "DF-OMP3",     "orb_opt": "FALSE", "nat_orbs": "FALSE",},
          "omp3":     {"wfn_type": "DF-OMP3",     "orb_opt": "TRUE",  "nat_orbs": "FALSE",},

         "remp2":     {"wfn_type": "DF-OREMP",    "orb_opt": "FALSE", "nat_orbs": "FALSE",},
        "oremp2":     {"wfn_type": "DF-OREMP",    "orb_opt": "TRUE",  "nat_orbs": "FALSE",},

          "lccd":     {"wfn_type": "DF-OLCCD",    "orb_opt": "FALSE", "nat_orbs": "FALSE",},
         "olccd":     {"wfn_type": "DF-OLCCD",    "orb_opt": "TRUE",  "nat_orbs": "FALSE",},

           "ccd":     {"wfn_type": "DF-CCD",      "orb_opt": "FALSE", "nat_orbs": "FALSE",},  # changes to DF-OCCD

           "ccsd":    {"wfn_type": "DF-CCSD",     "orb_opt": "FALSE", "nat_orbs": "FALSE",},

           "ccsd(t)": {"wfn_type": "DF-CCSD(T)",  "orb_opt": "FALSE", "nat_orbs": "FALSE",},
    }

    if name not in director:
        raise ValidationError(f"Invalid method {name} for DFOCC gradient")

    # throw exception for CONV (approximately)
    if (corl_type := method_algorithm_type(name).now) not in ["DF", "CD"]:
        raise ValidationError(f"Invalid type {corl_type} for DFOCC gradient. See Capabilities Table at {dtl}")

    proc_util.check_disk_df(name.upper(), optstash)

    # throw exception for SCF_TYPE
    if core.get_global_option('SCF_TYPE') != 'DISK_DF':
        raise ValidationError('DFOCC gradients need DF-SCF reference.')

    for k, v in director[name].items():
        core.set_local_option("DFOCC", k.upper(), v)

    core.set_global_option('DERTYPE', 'FIRST')
    core.set_local_option('DFOCC', 'DO_SCS', 'FALSE')
    core.set_local_option('DFOCC', 'DO_SOS', 'FALSE')
    core.set_local_option('SCF', 'DF_INTS_IO', 'SAVE')

    # Bypass the scf call if a reference wavefunction is given
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, use_c1=True, **kwargs)  # C1 certified
    else:
        if ref_wfn.molecule().schoenflies_symbol() != 'c1':
            raise ValidationError("""  DFOCC does not make use of molecular symmetry: """
                                  """reference wavefunction must be C1.\n""")

    core.print_out("  Constructing Basis Sets for DFOCC...\n\n")
    scf_aux_basis = core.BasisSet.build(ref_wfn.molecule(), "DF_BASIS_SCF",
                                    core.get_option("SCF", "DF_BASIS_SCF"),
                                    "JKFIT", core.get_global_option('BASIS'),
                                    puream=ref_wfn.basisset().has_puream())
    ref_wfn.set_basisset("DF_BASIS_SCF", scf_aux_basis)

    aux_basis = core.BasisSet.build(ref_wfn.molecule(), "DF_BASIS_CC",
                                        core.get_global_option("DF_BASIS_CC"),
                                        "RIFIT", core.get_global_option("BASIS"))
    ref_wfn.set_basisset("DF_BASIS_CC", aux_basis)

    if core.get_option('SCF', 'REFERENCE') == 'ROHF':
        ref_wfn.semicanonicalize()
    dfocc_wfn = core.dfocc(ref_wfn)

    derivobj = core.Deriv(dfocc_wfn)
    derivobj.compute_df("DF_BASIS_SCF", "DF_BASIS_CC")

    dfocc_wfn.set_variable(f"{name.upper()} TOTAL GRADIENT", dfocc_wfn.gradient())

    # Shove variables into global space
    for k, v in dfocc_wfn.variables().items():
        core.set_variable(k, v)

    optstash.restore()
    return dfocc_wfn


def run_dfocc_property(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a density-fitted (non-)orbital-optimized MPN or CC computation.

    """
    optstash = p4util.OptionsState(
        ['SCF', 'DF_INTS_IO'],
        ['DFOCC', 'WFN_TYPE'],
        ['DFOCC', 'ORB_OPT'],
        ['DFOCC', 'OEPROP'])

    if name in ['mp2', 'omp2']:
        core.set_local_option('DFOCC', 'WFN_TYPE', 'DF-OMP2')
    elif name in ['omp3']:
        core.set_local_option('DFOCC', 'WFN_TYPE', 'DF-OMP3')
    elif name in ['omp2.5']:
        core.set_local_option('DFOCC', 'WFN_TYPE', 'DF-OMP2.5')
    elif name in ['olccd']:
        core.set_local_option('DFOCC', 'WFN_TYPE', 'DF-OLCCD')
    else:
        raise ValidationError('Unidentified method ' % (name))

    proc_util.check_disk_df(name.upper(), optstash)

    if name in ['mp2']:
        core.set_local_option('DFOCC', 'ORB_OPT', 'FALSE')
    elif name in ['omp2', 'omp3', 'omp2.5', 'olccd']:
        core.set_local_option('DFOCC', 'ORB_OPT', 'TRUE')

    core.set_local_option('DFOCC', 'OEPROP', 'TRUE')
    core.set_local_option('DFOCC', 'DO_SCS', 'FALSE')
    core.set_local_option('DFOCC', 'DO_SOS', 'FALSE')
    core.set_local_option('SCF', 'DF_INTS_IO', 'SAVE')

    # Bypass the scf call if a reference wavefunction is given
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, use_c1=True, **kwargs)  # C1 certified
    else:
        if ref_wfn.molecule().schoenflies_symbol() != 'c1':
            raise ValidationError("""  DFOCC does not make use of molecular symmetry: """
                                  """reference wavefunction must be C1.\n""")

    core.print_out("  Constructing Basis Sets for DFOCC...\n\n")
    scf_aux_basis = core.BasisSet.build(ref_wfn.molecule(), "DF_BASIS_SCF",
                                    core.get_option("SCF", "DF_BASIS_SCF"),
                                    "JKFIT", core.get_global_option('BASIS'),
                                    puream=ref_wfn.basisset().has_puream())
    ref_wfn.set_basisset("DF_BASIS_SCF", scf_aux_basis)

    aux_basis = core.BasisSet.build(ref_wfn.molecule(), "DF_BASIS_CC",
                                        core.get_global_option("DF_BASIS_CC"),
                                        "RIFIT", core.get_global_option("BASIS"))
    ref_wfn.set_basisset("DF_BASIS_CC", aux_basis)

    if core.get_option('SCF', 'REFERENCE') == 'ROHF':
        ref_wfn.semicanonicalize()
    dfocc_wfn = core.dfocc(ref_wfn)

    # Shove variables into global space
    # TODO: Make other methods in DFOCC update all variables, then add them to the list. Adding now, risks setting outdated information.
    if name in ['mp2', 'omp2']:
        for k, v in dfocc_wfn.variables().items():
            core.set_variable(k, v)

    optstash.restore()
    return dfocc_wfn


def run_qchf(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    an quadratically-convergent SCF computation.

    """
    dtl = docs_table_link("dummy", "occ_nonoo")

    optstash = p4util.OptionsState(
        ['SCF', 'DF_INTS_IO'],
        ['DF_BASIS_SCF'],
        ['SCF', 'FAIL_ON_MAXITER'],
        ['MAXITER'],
        ['DFOCC', 'ORB_OPT'],
        ['DFOCC', 'WFN_TYPE'],
        ['DFOCC', 'QCHF'],
        ['DFOCC', 'E_CONVERGENCE'])

    # throw exception for CONV
    if (corl_type := method_algorithm_type(name).now) not in ["DISK_DF", "DF", "CD"]:
        raise ValidationError(f"Invalid type {corl_type} for QCHF energy through `run_qchf`. See Capabilities Table at {dtl}")

    core.set_local_option('DFOCC', 'ORB_OPT', 'TRUE')
    core.set_local_option('DFOCC', 'WFN_TYPE', 'QCHF')
    core.set_local_option('DFOCC', 'QCHF', 'TRUE')
    core.set_local_option('DFOCC', 'E_CONVERGENCE', 8)

    core.set_local_option('SCF', 'DF_INTS_IO', 'SAVE')
    core.set_local_option('SCF', 'FAIL_ON_MAXITER', False)
    core.set_local_option('SCF', 'MAXITER', 1)

    # Bypass the scf call if a reference wavefunction is given
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, use_c1=True, **kwargs)  # C1 certified
    else:
        if ref_wfn.molecule().schoenflies_symbol() != 'c1':
            raise ValidationError("""  QCHF does not make use of molecular symmetry: """
                                  """reference wavefunction must be C1.\n""")

    if core.get_option('SCF', 'REFERENCE') == 'ROHF':
        ref_wfn.semicanonicalize()
    dfocc_wfn = core.dfocc(ref_wfn)

    # Shove variables into global space
    for k, v in dfocc_wfn.variables().items():
        core.set_variable(k, v)

    optstash.restore()
    return dfocc_wfn


def run_occ(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a conventional integral (O)MPN computation

    """
    # Stash these options so we can reload them at computation end.
    optstash = p4util.OptionsState(
        ['OCC', 'SPIN_SCALE_TYPE'],
        ['OCC', 'ORB_OPT'],
        ['OCC', 'WFN_TYPE'])

    director = {
                  "mp2":     {"wfn_type": "OMP2",   "orb_opt": "FALSE", "spin_scale_type": "NONE",  },
              "scs-mp2":     {"wfn_type": "OMP2",   "orb_opt": "FALSE", "spin_scale_type": "SCS",   },
           "scs(n)-mp2":     {"wfn_type": "OMP2",   "orb_opt": "FALSE", "spin_scale_type": "SCSN",  },
              "scs-mp2-vdw": {"wfn_type": "OMP2",   "orb_opt": "FALSE", "spin_scale_type": "SCSVDW",},
              "sos-mp2":     {"wfn_type": "OMP2",   "orb_opt": "FALSE", "spin_scale_type": "SOS",   },
           "sos-pi-mp2":     {"wfn_type": "OMP2",   "orb_opt": "FALSE", "spin_scale_type": "SOSPI", },
       "custom-scs-mp2":     {"wfn_type": "OMP2",   "orb_opt": "FALSE", "spin_scale_type": "CUSTOM",},

                 "omp2":     {"wfn_type": "OMP2",   "orb_opt": "TRUE",  "spin_scale_type": "NONE",  },
             "scs-omp2":     {"wfn_type": "OMP2",   "orb_opt": "TRUE",  "spin_scale_type": "SCS",   },
             "sos-omp2":     {"wfn_type": "OMP2",   "orb_opt": "TRUE",  "spin_scale_type": "SOS",   },
      "custom-scs-omp2":     {"wfn_type": "OMP2",   "orb_opt": "TRUE",  "spin_scale_type": "CUSTOM",},

                  "mp2.5":   {"wfn_type": "OMP2.5", "orb_opt": "FALSE", "spin_scale_type": "NONE",  },
       "custom-scs-mp2.5":   {"wfn_type": "OMP2.5", "orb_opt": "FALSE", "spin_scale_type": "CUSTOM",},

                 "omp2.5":   {"wfn_type": "OMP2.5", "orb_opt": "TRUE",  "spin_scale_type": "NONE",  },
      "custom-scs-omp2.5":   {"wfn_type": "OMP2.5", "orb_opt": "TRUE",  "spin_scale_type": "CUSTOM",},

                  "mp3":     {"wfn_type": "OMP3",   "orb_opt": "FALSE", "spin_scale_type": "NONE",  },
              "scs-mp3":     {"wfn_type": "OMP3",   "orb_opt": "FALSE", "spin_scale_type": "SCS",   },
       "custom-scs-mp3":     {"wfn_type": "OMP3",   "orb_opt": "FALSE", "spin_scale_type": "CUSTOM",},

                 "omp3":     {"wfn_type": "OMP3",   "orb_opt": "TRUE",  "spin_scale_type": "NONE",  },
             "scs-omp3":     {"wfn_type": "OMP3",   "orb_opt": "TRUE",  "spin_scale_type": "SCS",   },
             "sos-omp3":     {"wfn_type": "OMP3",   "orb_opt": "TRUE",  "spin_scale_type": "SOS",   },
      "custom-scs-omp3":     {"wfn_type": "OMP3",   "orb_opt": "TRUE",  "spin_scale_type": "CUSTOM",},

                "remp2":     {"wfn_type": "REMP",   "orb_opt": "FALSE", "spin_scale_type": "NONE",  },

               "oremp2":     {"wfn_type": "OREMP",  "orb_opt": "TRUE",  "spin_scale_type": "NONE",  },

                 "lccd":     {"wfn_type": "OCEPA",  "orb_opt": "FALSE", "spin_scale_type": "NONE",  },
      "custom-scs-lccd":     {"wfn_type": "OCEPA",  "orb_opt": "FALSE", "spin_scale_type": "CUSTOM",},

                "olccd":     {"wfn_type": "OCEPA",  "orb_opt": "TRUE",  "spin_scale_type": "NONE",  },
     "custom-scs-olccd":     {"wfn_type": "OCEPA",  "orb_opt": "TRUE",  "spin_scale_type": "CUSTOM",},
    }

    if name not in director:
        raise ValidationError(f"Invalid method {name} for OCC energy")

    for k, v in director[name].items():
        core.set_local_option("OCC", k.upper(), v)

    # Bypass the scf call if a reference wavefunction is given
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, **kwargs)  # C1 certified

    # Ensure IWL files have been written
    proc_util.check_iwl_file_from_scf_type(core.get_global_option('SCF_TYPE'), ref_wfn)

    if core.get_option('SCF', 'REFERENCE') == 'ROHF':
        ref_wfn.semicanonicalize()

    occ_wfn = core.occ(ref_wfn)

    # Shove variables into global space
    keep_custom_spin_scaling = core.has_option_changed("OCC", "SS_SCALE") or core.has_option_changed("OCC", "OS_SCALE")
    for k, v in occ_wfn.variables().items():
        # Custom spin component scaling variables are meaningless if custom scalings hasn't been set. Delete them.
        if k.startswith("CUSTOM SCS") and not keep_custom_spin_scaling:
            occ_wfn.del_variable(k)
        else:
            core.set_variable(k, v)

    optstash.restore()
    return occ_wfn


def run_occ_gradient(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a conventional integral (O)MPN computation
    """
    optstash = p4util.OptionsState(
        ['OCC', 'ORB_OPT'],
        ['OCC', 'WFN_TYPE'],
        ['OCC', 'DO_SCS'],
        ['OCC', 'DO_SOS'],
        ['GLOBALS', 'DERTYPE'])

    if core.get_global_option('SCF_TYPE') in ['CD', 'DF', 'MEM_DF', 'DISK_DF']:
        raise ValidationError('OCC gradients need conventional SCF reference.')

    director = {
           "mp2":   {"wfn_type": "OMP2",   "orb_opt": "FALSE",},
          "omp2":   {"wfn_type": "OMP2",   "orb_opt": "TRUE", },
     "conv-omp2":   {"wfn_type": "OMP2",   "orb_opt": "TRUE", },

           "mp2.5": {"wfn_type": "OMP2.5", "orb_opt": "FALSE",},
          "omp2.5": {"wfn_type": "OMP2.5", "orb_opt": "TRUE", },

           "mp3":   {"wfn_type": "OMP3",   "orb_opt": "FALSE",},
          "omp3":   {"wfn_type": "OMP3",   "orb_opt": "TRUE", },

        "oremp2":   {"wfn_type": "OREMP",  "orb_opt": "TRUE", },

          "lccd":   {"wfn_type": "OCEPA",  "orb_opt": "FALSE",},
         "olccd":   {"wfn_type": "OCEPA",  "orb_opt": "TRUE", },
    }

    if name not in director:
        raise ValidationError(f"Invalid method {name} for OCC gradient")

    for k, v in director[name].items():
        core.set_local_option("OCC", k.upper(), v)

    core.set_global_option('DERTYPE', 'FIRST')

    # locking out SCS through explicit keyword setting
    # * so that current energy must match call
    # * since grads not avail for scs
    core.set_local_option('OCC', 'SPIN_SCALE_TYPE', 'NONE')

    # Bypass the scf call if a reference wavefunction is given
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, **kwargs)  # C1 certified

    # Ensure IWL files have been written
    proc_util.check_iwl_file_from_scf_type(core.get_global_option('SCF_TYPE'), ref_wfn)

    if core.get_option('SCF', 'REFERENCE') == 'ROHF':
        ref_wfn.semicanonicalize()

    occ_wfn = core.occ(ref_wfn)

    derivobj = core.Deriv(occ_wfn)
    grad = derivobj.compute()

    occ_wfn.set_gradient(grad)
    occ_wfn.set_variable(f"{name.upper()} TOTAL GRADIENT", grad)

    # Shove variables into global space
    keep_custom_spin_scaling = core.has_option_changed("OCC", "SS_SCALE") or core.has_option_changed("OCC", "OS_SCALE")
    for k, v in occ_wfn.variables().items():
        # Custom spin component scaling variables are meaningless if custom scalings hasn't been set. Delete them.
        if k.startswith("CUSTOM SCS") and not keep_custom_spin_scaling:
            occ_wfn.del_variable(k)
        else:
            core.set_variable(k, v)

    optstash.restore()
    return occ_wfn


def run_scf(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a self-consistent-field theory (HF & DFT) calculation.
    """
    optstash_mp2 = p4util.OptionsState(
        ['DF_BASIS_MP2'],
        ['DFMP2', 'MP2_OS_SCALE'],
        ['DFMP2', 'MP2_SS_SCALE'])

    dft_func = False
    if "dft_functional" in kwargs:
        dft_func = True

    optstash_scf = proc_util.scf_set_reference_local(name, is_dft=dft_func)

    # See if we're doing TDSCF after, keep JK if so
    if sum(core.get_option("SCF", "TDSCF_STATES")) > 0:
        core.set_local_option("SCF", "SAVE_JK", True)

    # Alter default algorithm
    if not core.has_global_option_changed('SCF_TYPE'):
        core.set_global_option('SCF_TYPE', 'DF')

    scf_wfn = scf_helper(name, post_scf=False, **kwargs)
    returnvalue = scf_wfn.energy()

    ssuper = scf_wfn.functional()

    if ssuper.is_c_hybrid():

        # throw exception for CONV/CD MP2
        if (mp2_type := core.get_global_option("MP2_TYPE")) != "DF":
            raise ValidationError(f"Invalid MP2 type {mp2_type} for DF-DFT energy. See capabilities Table.")

        core.tstart()
        aux_basis = core.BasisSet.build(scf_wfn.molecule(), "DF_BASIS_MP2",
                                        core.get_option("DFMP2", "DF_BASIS_MP2"),
                                        "RIFIT", core.get_global_option('BASIS'),
                                        puream=-1)
        scf_wfn.set_basisset("DF_BASIS_MP2", aux_basis)
        if ssuper.is_c_scs_hybrid():
            core.set_local_option('DFMP2', 'MP2_OS_SCALE', ssuper.c_os_alpha())
            core.set_local_option('DFMP2', 'MP2_SS_SCALE', ssuper.c_ss_alpha())
            dfmp2_wfn = core.dfmp2(scf_wfn)
            dfmp2_wfn.compute_energy()

            vdh = dfmp2_wfn.variable('CUSTOM SCS-MP2 CORRELATION ENERGY')

        else:
            dfmp2_wfn = core.dfmp2(scf_wfn)
            dfmp2_wfn.compute_energy()
            vdh = ssuper.c_alpha() * dfmp2_wfn.variable('MP2 CORRELATION ENERGY')

        # remove misleading MP2 psivars computed with DFT, not HF, reference
        for var in dfmp2_wfn.variables():
            if var.startswith('MP2 ') and ssuper.name() not in ['MP2D']:
                scf_wfn.del_variable(var)

        scf_wfn.set_variable("DOUBLE-HYBRID CORRECTION ENERGY", vdh)  # P::e SCF
        scf_wfn.set_variable("{} DOUBLE-HYBRID CORRECTION ENERGY".format(ssuper.name()), vdh)
        returnvalue += vdh
        scf_wfn.set_variable("DFT TOTAL ENERGY", returnvalue)  # P::e SCF
        for pv, pvv in scf_wfn.variables().items():
            if pv.endswith('DISPERSION CORRECTION ENERGY') and pv.startswith(ssuper.name()):
                fctl_plus_disp_name = pv.split()[0]
                scf_wfn.set_variable(fctl_plus_disp_name + ' TOTAL ENERGY', returnvalue)
                break
        else:
            scf_wfn.set_variable('{} TOTAL ENERGY'.format(ssuper.name()), returnvalue)

        scf_wfn.set_variable('CURRENT ENERGY', returnvalue)
        scf_wfn.set_energy(returnvalue)
        core.print_out('\n\n')
        core.print_out('    %s Energy Summary\n' % (name.upper()))
        core.print_out('    ' + '-' * (15 + len(name)) + '\n')
        core.print_out('    DFT Reference Energy                  = %22.16lf\n' % (returnvalue - vdh))
        core.print_out('    Scaled MP2 Correlation                = %22.16lf\n' % (vdh))
        core.print_out('    @Final double-hybrid DFT total energy = %22.16lf\n\n' % (returnvalue))
        core.tstop()

        if ssuper.name() == 'MP2D':
            for pv, pvv in dfmp2_wfn.variables().items():
                scf_wfn.set_variable(pv, pvv)

            # Conversely, remove DFT qcvars from MP2D
            for var in scf_wfn.variables():
                if 'DFT ' in var or 'DOUBLE-HYBRID ' in var:
                    scf_wfn.del_variable(var)

            # DFT groups dispersion with SCF. Reshuffle so dispersion with MP2 for MP2D.
            for pv in ['SCF TOTAL ENERGY', 'SCF ITERATION ENERGY', 'MP2 TOTAL ENERGY']:
                scf_wfn.set_variable(pv, scf_wfn.variable(pv) - scf_wfn.variable('DISPERSION CORRECTION ENERGY'))

            scf_wfn.set_variable('MP2D CORRELATION ENERGY', scf_wfn.variable('MP2 CORRELATION ENERGY') + scf_wfn.variable('DISPERSION CORRECTION ENERGY'))
            scf_wfn.set_variable('MP2D TOTAL ENERGY', scf_wfn.variable('MP2D CORRELATION ENERGY') + scf_wfn.variable('HF TOTAL ENERGY'))
            scf_wfn.set_variable('CURRENT ENERGY', scf_wfn.variable('MP2D TOTAL ENERGY'))
            scf_wfn.set_variable('CURRENT CORRELATION ENERGY', scf_wfn.variable('MP2D CORRELATION ENERGY'))
            scf_wfn.set_variable('CURRENT REFERENCE ENERGY', scf_wfn.variable('SCF TOTAL ENERGY'))

    # Shove variables into global space
    for k, v in scf_wfn.variables().items():
        core.set_variable(k, v)

    optstash_scf.restore()
    optstash_mp2.restore()
    return scf_wfn


def run_scf_gradient(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a SCF gradient calculation.

    """

    dft_func = False
    if "dft_functional" in kwargs:
        dft_func = True

    optstash = proc_util.scf_set_reference_local(name, is_dft=dft_func)

    # Bypass the scf call if a reference wavefunction is given
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = run_scf(name, **kwargs)

    if core.get_option('SCF', 'REFERENCE') in ['ROHF', 'CUHF']:
        ref_wfn.semicanonicalize()

    if hasattr(ref_wfn, "_disp_functor"):
        disp_grad = ref_wfn._disp_functor.compute_gradient(ref_wfn.molecule(), ref_wfn)
        ref_wfn.set_variable("-D Gradient", disp_grad)

    grad = core.scfgrad(ref_wfn)

    ref_wfn.set_gradient(grad)

    ref_wfn.set_variable("SCF TOTAL GRADIENT", grad)  # P::e SCF
    if ref_wfn.functional().needs_xc():
        ref_wfn.set_variable("DFT TOTAL GRADIENT", grad)  # overwritten later for DH -- TODO when DH gradients  # P::e SCF
    else:
        ref_wfn.set_variable("HF TOTAL GRADIENT", grad)  # P::e SCF

    # Shove variables into global space
    for k, v in ref_wfn.variables().items():
        core.set_variable(k, v)

    optstash.restore()
    return ref_wfn


def run_scf_hessian(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    an SCF hessian calculation.

    """
    optstash = proc_util.scf_set_reference_local(name)

    # Bypass the scf call if a reference wavefunction is given
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = run_scf(name, **kwargs)

    badref = core.get_option('SCF', 'REFERENCE') in ['ROHF', 'CUHF', 'UKS']
    badint = core.get_global_option('SCF_TYPE') in [ 'CD', 'OUT_OF_CORE']
    if badref or badint:
        raise ValidationError("Only RHF/UHF/RKS Hessians are currently implemented. SCF_TYPE either CD or OUT_OF_CORE not supported")

    if hasattr(ref_wfn, "_disp_functor"):
        disp_hess = ref_wfn._disp_functor.compute_hessian(ref_wfn.molecule(), ref_wfn)
        ref_wfn.set_variable("-D Hessian", disp_hess)

    H = core.scfhess(ref_wfn)
    ref_wfn.set_hessian(H)

    ref_wfn.set_variable("SCF TOTAL HESSIAN", H)  # P::e SCF
    if ref_wfn.functional().needs_xc():
        ref_wfn.set_variable("DFT TOTAL HESSIAN", H)  # overwritten later for DH -- TODO when DH Hessians # P::e SCF
    else:
        ref_wfn.set_variable("HF TOTAL HESSIAN", H)  # P::e SCF

    # Shove variables into global space
    for k, v in ref_wfn.variables().items():
        core.set_variable(k, v)

    optstash.restore()
    return ref_wfn


def run_mcscf(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a multiconfigurational self-consistent-field calculation.

    """
    # Make sure the molecule the user provided is the active one
    mcscf_molecule = kwargs.get('molecule', core.get_active_molecule())
    mcscf_molecule.update_geometry()
    if 'ref_wfn' in kwargs:
        raise ValidationError("It is not possible to pass run_mcscf a reference wavefunction")
    new_wfn = core.Wavefunction.build(mcscf_molecule, core.get_global_option('BASIS'))

    return core.mcscf(new_wfn)


def run_dfmp2_gradient(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a DFMP2 gradient calculation.

    """
    optstash = p4util.OptionsState(
        ['DF_BASIS_SCF'],
        ['DF_BASIS_MP2'],
        ['SCF_TYPE'])

    # Alter default algorithm
    if not core.has_global_option_changed('SCF_TYPE'):
        core.set_global_option('SCF_TYPE', 'DF')
        core.print_out("""    SCF Algorithm Type (re)set to DF.\n""")

    if "DF" not in core.get_global_option('SCF_TYPE'):
        raise ValidationError('DF-MP2 gradients need DF-SCF reference.')

    # Bypass the scf call if a reference wavefunction is given
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, **kwargs)  # C1 certified

    if ref_wfn.basisset().has_ECP():
        raise ValidationError('DF-MP2 gradients with an ECP are not yet available.  Use dertype=0 to select numerical gradients.')

    core.tstart()
    core.print_out('\n')
    p4util.banner('DFMP2')
    core.print_out('\n')

    aux_basis = core.BasisSet.build(ref_wfn.molecule(), "DF_BASIS_MP2",
                                    core.get_option("DFMP2", "DF_BASIS_MP2"),
                                    "RIFIT", core.get_global_option('BASIS'))
    ref_wfn.set_basisset("DF_BASIS_MP2", aux_basis)

    dfmp2_wfn = core.dfmp2(ref_wfn)
    grad = dfmp2_wfn.compute_gradient()
    dfmp2_wfn.set_gradient(grad)

    # Shove variables into global space
    dfmp2_wfn.set_variable("MP2 TOTAL GRADIENT", grad)  # P::e DFMP2
    dfmp2_wfn.set_variable('CURRENT ENERGY', dfmp2_wfn.variable('MP2 TOTAL ENERGY'))
    dfmp2_wfn.set_variable('CURRENT CORRELATION ENERGY', dfmp2_wfn.variable('MP2 CORRELATION ENERGY'))
    for k, v in dfmp2_wfn.variables().items():
        core.set_variable(k, v)

    optstash.restore()
    core.tstop()
    return dfmp2_wfn


def run_dfmp2d_gradient(name, **kwargs):
    """Encode MP2-D method."""

    dfmp2_wfn = run_dfmp2_gradient('mp2', **kwargs)
    wfn_grad = dfmp2_wfn.gradient().clone()

    _, _disp_functor = build_disp_functor('MP2D', restricted=True)
    disp_grad = _disp_functor.compute_gradient(dfmp2_wfn.molecule(), dfmp2_wfn)
    wfn_grad.add(disp_grad)
    dfmp2_wfn.set_gradient(wfn_grad)

    dfmp2_wfn.set_variable('MP2D CORRELATION ENERGY', dfmp2_wfn.variable('MP2 CORRELATION ENERGY') + dfmp2_wfn.variable('DISPERSION CORRECTION ENERGY'))
    dfmp2_wfn.set_variable('MP2D TOTAL ENERGY', dfmp2_wfn.variable('MP2D CORRELATION ENERGY') + dfmp2_wfn.variable('HF TOTAL ENERGY'))
    dfmp2_wfn.set_variable('CURRENT ENERGY', dfmp2_wfn.variable('MP2D TOTAL ENERGY'))
    dfmp2_wfn.set_variable('CURRENT CORRELATION ENERGY', dfmp2_wfn.variable('MP2D CORRELATION ENERGY'))

    # Shove variables into global space
    for k, v in dfmp2_wfn.variables().items():
        core.set_variable(k, v)

    return dfmp2_wfn


def run_ccenergy(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a CCSD, CC2, and CC3 calculation.

    """
    optstash = p4util.OptionsState(
        ['TRANSQT2', 'WFN'],
        ['CCSORT', 'WFN'],
        ['CCENERGY', 'WFN'])

    if name == 'ccsd':
        core.set_local_option('TRANSQT2', 'WFN', 'CCSD')
        core.set_local_option('CCSORT', 'WFN', 'CCSD')
        core.set_local_option('CCTRANSORT', 'WFN', 'CCSD')
        core.set_local_option('CCENERGY', 'WFN', 'CCSD')
    elif name == 'ccsd(t)':
        core.set_local_option('TRANSQT2', 'WFN', 'CCSD_T')
        core.set_local_option('CCSORT', 'WFN', 'CCSD_T')
        core.set_local_option('CCTRANSORT', 'WFN', 'CCSD_T')
        core.set_local_option('CCENERGY', 'WFN', 'CCSD_T')
    elif name == 'a-ccsd(t)':
        core.set_local_option('TRANSQT2', 'WFN', 'CCSD_AT')
        core.set_local_option('CCSORT', 'WFN', 'CCSD_AT')
        core.set_local_option('CCTRANSORT', 'WFN', 'CCSD_AT')
        core.set_local_option('CCENERGY', 'WFN', 'CCSD_AT')
        core.set_local_option('CCHBAR', 'WFN', 'CCSD_AT')
        core.set_local_option('CCLAMBDA', 'WFN', 'CCSD_AT')
    elif name == 'cc2':
        core.set_local_option('TRANSQT2', 'WFN', 'CC2')
        core.set_local_option('CCSORT', 'WFN', 'CC2')
        core.set_local_option('CCTRANSORT', 'WFN', 'CC2')
        core.set_local_option('CCENERGY', 'WFN', 'CC2')
    elif name == 'cc3':
        core.set_local_option('TRANSQT2', 'WFN', 'CC3')
        core.set_local_option('CCSORT', 'WFN', 'CC3')
        core.set_local_option('CCTRANSORT', 'WFN', 'CC3')
        core.set_local_option('CCENERGY', 'WFN', 'CC3')
    elif name == 'eom-cc2':
        core.set_local_option('TRANSQT2', 'WFN', 'EOM_CC2')
        core.set_local_option('CCSORT', 'WFN', 'EOM_CC2')
        core.set_local_option('CCTRANSORT', 'WFN', 'EOM_CC2')
        core.set_local_option('CCENERGY', 'WFN', 'EOM_CC2')
    elif name == 'eom-ccsd':
        core.set_local_option('TRANSQT2', 'WFN', 'EOM_CCSD')
        core.set_local_option('CCSORT', 'WFN', 'EOM_CCSD')
        core.set_local_option('CCTRANSORT', 'WFN', 'EOM_CCSD')
        core.set_local_option('CCENERGY', 'WFN', 'EOM_CCSD')
    # Call a plain energy('ccenergy') and have full control over options, incl. wfn
    elif name == 'ccenergy':
        pass

    # Bypass routine scf if user did something special to get it to converge
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, **kwargs)  # C1 certified

    if core.get_global_option("CC_TYPE") == "DF":
        aux_basis = core.BasisSet.build(ref_wfn.molecule(), "DF_BASIS_CC",
                                            core.get_global_option("DF_BASIS_CC"),
                                            "RIFIT", core.get_global_option("BASIS"))
        ref_wfn.set_basisset("DF_BASIS_CC", aux_basis)

    # Ensure IWL files have been written
    proc_util.check_iwl_file_from_scf_type(core.get_global_option('SCF_TYPE'), ref_wfn)

    # Obtain semicanonical orbitals
    if ((core.get_option('SCF', 'REFERENCE') == 'ROHF')
        and ((name in ['ccsd(t)', 'a-ccsd(t)', 'cc2', 'cc3', 'eom-cc2', 'eom-cc3'])
            or core.get_option('CCTRANSORT', 'SEMICANONICAL'))):
        ref_wfn.semicanonicalize()

    if core.get_global_option('RUN_CCTRANSORT'):
        core.cctransort(ref_wfn)
    else:
        try:
            from psi4.driver.pasture import addins
            addins.ccsort_transqt2(ref_wfn)
        except Exception:
            raise PastureRequiredError("RUN_CCTRANSORT")


    ccwfn = core.ccenergy(ref_wfn)
    if core.get_global_option('PE'):
        ccwfn.pe_state = ref_wfn.pe_state

    if name == 'a-ccsd(t)':
        core.cchbar(ref_wfn)
        lambdawfn = core.cclambda(ref_wfn)
        for k, v in lambdawfn.variables().items():
            ccwfn.set_variable(k, v)

    optstash.restore()
    return ccwfn


def run_ccenergy_gradient(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a CCSD and CCSD(T) gradient calculation.

    """
    optstash = p4util.OptionsState(
        ['GLOBALS', 'DERTYPE'],
        ['CCLAMBDA', 'WFN'],
        ['CCDENSITY', 'WFN'])

    core.set_global_option('DERTYPE', 'FIRST')

    if core.get_global_option('FREEZE_CORE') not in ["FALSE", "0"]:
        raise ValidationError('Frozen core is not available for the CC gradients.')

    ccwfn = run_ccenergy(name, **kwargs)

    if name == 'cc2':
        core.set_local_option('CCHBAR', 'WFN', 'CC2')
        core.set_local_option('CCLAMBDA', 'WFN', 'CC2')
        core.set_local_option('CCDENSITY', 'WFN', 'CC2')
    if name == 'ccsd':
        core.set_local_option('CCLAMBDA', 'WFN', 'CCSD')
        core.set_local_option('CCDENSITY', 'WFN', 'CCSD')
    elif name == 'ccsd(t)':
        core.set_local_option('CCLAMBDA', 'WFN', 'CCSD_T')
        core.set_local_option('CCDENSITY', 'WFN', 'CCSD_T')

    core.cchbar(ccwfn)
    core.cclambda(ccwfn)
    core.ccdensity(ccwfn)

    derivobj = core.Deriv(ccwfn)
    grad = derivobj.compute()
    del derivobj

    ccwfn.set_gradient(grad)
    ccwfn.set_variable(f"{name.upper()} TOTAL GRADIENT", grad)
    core.set_variable(f"{name.upper()} TOTAL GRADIENT", grad)
    core.set_variable("CURRENT GRADIENT", grad)

    optstash.restore()
    return ccwfn


def run_bccd(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a Brueckner CCD calculation.

    """
    dtl = docs_table_link("dummy", "ccenergy")

    optstash = p4util.OptionsState(
        ['TRANSQT2', 'WFN'],
        ['CCSORT', 'WFN'],
        ['CCENERGY', 'WFN'])

    if name == 'bccd':
        core.set_local_option('TRANSQT2', 'WFN', 'BCCD')
        core.set_local_option('CCSORT', 'WFN', 'BCCD')
        core.set_local_option('CCTRANSORT', 'WFN', 'BCCD')
        core.set_local_option('CCENERGY', 'WFN', 'BCCD')

    elif name == 'bccd(t)':
        core.set_local_option('TRANSQT2', 'WFN', 'BCCD_T')
        core.set_local_option('CCSORT', 'WFN', 'BCCD_T')
        core.set_local_option('CCENERGY', 'WFN', 'BCCD_T')
        core.set_local_option('CCTRANSORT', 'WFN', 'BCCD_T')
        core.set_local_option('CCTRIPLES', 'WFN', 'BCCD_T')
    else:
        raise ValidationError("proc.py:run_bccd name %s not recognized" % name)

    if (corl_type := method_algorithm_type(name).now) != "CONV":
        raise ValidationError(f"Invalid type {corl_type} for CCENERGY energy through `run_bccd`. See Capabilities Table at {dtl}")

    # Bypass routine scf if user did something special to get it to converge
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, **kwargs)  # C1 certified

    # Needed for (T).
    if (core.get_option('SCF', 'REFERENCE') == 'ROHF'):
        ref_wfn.semicanonicalize()

    # Ensure IWL files have been written
    proc_util.check_iwl_file_from_scf_type(core.get_global_option('SCF_TYPE'), ref_wfn)

    core.set_local_option('CCTRANSORT', 'DELETE_TEI', 'false')

    bcc_iter_cnt = 0
    if (core.get_global_option("RUN_CCTRANSORT")):
        sort_func = core.cctransort
    else:
        try:
            from psi4.driver.pasture import addins
            core.set_local_option('TRANSQT2', 'DELETE_TEI', 'false')
            sort_func = addins.ccsort_transqt2
        except Exception:
            raise PastureRequiredError("RUN_CCTRANSORT")

    hold_qcvars = {
        "MP2 TOTAL ENERGY": None,
        "MP2 CORRELATION ENERGY": None,
        "MP2 SAME-SPIN CORRELATION ENERGY": None,
        "MP2 OPPOSITE-SPIN CORRELATION ENERGY": None,
        "MP2 SINGLES ENERGY": None,
        "MP2 DOUBLES ENERGY": None,
        "CCSD TOTAL ENERGY": None,
        "CCSD CORRELATION ENERGY": None,
        "CCSD SAME-SPIN CORRELATION ENERGY": None,
        "CCSD OPPOSITE-SPIN CORRELATION ENERGY": None,
        "CCSD SINGLES ENERGY": None,
        "CCSD DOUBLES ENERGY": None,
    }

    while True:
        sort_func(ref_wfn)

        ref_wfn = core.ccenergy(ref_wfn)
        core.print_out('Brueckner convergence check: %s\n' % bool(core.variable('BRUECKNER CONVERGED')))
        if core.variable('BRUECKNER CONVERGED'):
            break

        if bcc_iter_cnt >= core.get_option('CCENERGY', 'BCCD_MAXITER'):
            core.print_out("\n\nWarning! BCCD did not converge within the maximum number of iterations.")
            core.print_out("You can increase the number of BCCD iterations by changing BCCD_MAXITER.\n\n")
            break
        bcc_iter_cnt += 1

        if bcc_iter_cnt == 1:
            for pv in hold_qcvars:
                hold_qcvars[pv] = ref_wfn.variable(pv)

    ref_wfn.set_variable("BCCD TOTAL ENERGY", ref_wfn.variable("CCSD TOTAL ENERGY"))
    ref_wfn.set_variable("BCCD CORRELATION ENERGY", ref_wfn.variable("BCCD TOTAL ENERGY") - ref_wfn.variable("SCF TOTAL ENERGY"))
    ref_wfn.set_variable("CURRENT CORRELATION ENERGY", ref_wfn.variable("BCCD CORRELATION ENERGY"))

    # copy back canonical MP2 and CCSD from initial iteration
    for pv, v in hold_qcvars.items():
        if v is not None:
            ref_wfn.set_variable(pv, v)
            core.set_variable(pv, v)

    if name == 'bccd(t)':
        core.cctriples(ref_wfn)
        ref_wfn.set_variable("B(T) CORRECTION ENERGY", ref_wfn.variable("(T) CORRECTION ENERGY"))
        ref_wfn.set_variable("BCCD(T) TOTAL ENERGY", ref_wfn.variable("CCSD(T) TOTAL ENERGY"))
        ref_wfn.set_variable("BCCD(T) CORRELATION ENERGY", ref_wfn.variable("BCCD(T) TOTAL ENERGY") - ref_wfn.variable("SCF TOTAL ENERGY"))  # note != CCSD(T) CORRELATION ENERGY
        ref_wfn.set_variable("CURRENT CORRELATION ENERGY", ref_wfn.variable("BCCD(T) CORRELATION ENERGY"))

        for pv in ["(T) CORRECTION ENERGY", "CCSD(T) TOTAL ENERGY", "CCSD(T) CORRELATION ENERGY"]:
            ref_wfn.del_variable(pv)
            core.del_variable(pv)

    for pv in [
        "BCCD TOTAL ENERGY",
        "BCCD CORRELATION ENERGY",
        "B(T) CORRECTION ENERGY",
        "BCCD(T) TOTAL ENERGY",
        "BCCD(T) CORRELATION ENERGY",
        "CURRENT CORRELATION ENERGY",
    ]:
        if ref_wfn.has_variable(pv):
            core.set_variable(pv, ref_wfn.variable(pv))

    # Notes
    # * BCCD or BCCD(T) correlation energy is total energy of last Brueckner iteration minus HF energy of first Brueckner iteration

    optstash.restore()
    return ref_wfn


def run_tdscf_excitations(wfn,**kwargs):

    states = core.get_option("SCF","TDSCF_STATES")

    # some sanity checks
    if sum(states) == 0:
        raise ValidationError("TDSCF: No states requested in TDSCF_STATES")

    # unwrap 1-membered list of states, regardless of symmetry
    # we will apportion states per irrep later on
    if len(states) == 1:
        states = states[0]

    # Tie TDSCF_R_CONVERGENCE to D_CONVERGENCE in SCF reference
    if core.has_option_changed('SCF', 'TDSCF_R_CONVERGENCE'):
        r_convergence = core.get_option('SCF', 'TDSCF_R_CONVERGENCE')
    else:
        r_convergence = min(1.e-4, core.get_option('SCF', 'D_CONVERGENCE') * 1.e2)

    # "anonymous" return value, as we stash observables in the passed Wavefunction object internally
    _ = response.scf_response.tdscf_excitations(wfn,
                                                states=states,
                                                triplets=core.get_option("SCF", "TDSCF_TRIPLETS"),
                                                tda=core.get_option("SCF", "TDSCF_TDA"),
                                                r_convergence=r_convergence,
                                                maxiter=core.get_option("SCF", "TDSCF_MAXITER"),
                                                guess=core.get_option("SCF", "TDSCF_GUESS"),
                                                verbose=core.get_option("SCF", "TDSCF_PRINT"),
                                                coeff_cutoff=core.get_option("SCF", "TDSCF_COEFF_CUTOFF"),
                                                tdm_print=core.get_option("SCF", "TDSCF_TDM_PRINT"))

    # Shove variables into global space
    for k, v in wfn.variables().items():
        core.set_variable(k, v)

    return wfn


def run_tdscf_energy(name, **kwargs):

    # Get a wfn in case we aren't given one
    ref_wfn = kwargs.get('ref_wfn', None)

    if ref_wfn is None:
        if name is None:
            raise ValidationError("TDSCF: No reference wave function!")
        else:
            ref_wfn = run_scf(name.strip('td-'), **kwargs)

    return run_tdscf_excitations(ref_wfn, **kwargs)


def run_scf_property(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    SCF calculations. This is a simple alias to :py:func:`~proc.run_scf`
    since SCF properties all handled through oeprop.

    """

    core.tstart()
    optstash = proc_util.scf_set_reference_local(name)

    properties = kwargs.pop('properties')

    # What response do we need?
    response_list_vals = list(response.scf_response.property_dicts)
    oeprop_list_vals = core.OEProp.valid_methods

    oe_properties = []
    linear_response = []
    unknown_property = []
    for prop in properties:

        prop = prop.upper()
        if prop in response_list_vals:
            linear_response.append(prop)
        elif (prop in oeprop_list_vals) or ("MULTIPOLE(" in prop):
            oe_properties.append(prop)
        else:
            unknown_property.append(prop)

    if "DIPOLE" not in oe_properties:
        oe_properties.append("DIPOLE")

    # Throw if we dont know what something is
    if len(unknown_property):
        complete_options = oeprop_list_vals + response_list_vals
        alt_method_name = p4util.text.find_approximate_string_matches(unknown_property[0],
                                                         complete_options, 2)
        alternatives = ""
        if len(alt_method_name) > 0:
            alternatives = " Did you mean? %s" % (" ".join(alt_method_name))

        raise ValidationError("SCF Property: Feature '%s' is not recognized. %s" % (unknown_property[0], alternatives))

    # Validate OEProp
    if len(oe_properties):
        proc_util.oeprop_validator(oe_properties)

    if len(linear_response):
        optstash_jk = p4util.OptionsState(["SAVE_JK"])
        core.set_global_option("SAVE_JK", True)

    # Compute the Wavefunction
    scf_wfn = run_scf(name, scf_do_properties=False, do_timer=False, **kwargs)

    # Run OEProp
    oe = core.OEProp(scf_wfn)
    oe.set_title(name.upper())
    for prop in oe_properties:
        oe.add(prop.upper())
    oe.compute()
    scf_wfn.oeprop = oe

    core.set_variable("SCF DIPOLE", core.variable(name + " DIPOLE"))  # P::e SCF

    # Run Linear Respsonse
    if len(linear_response):
        core.prepare_options_for_module("SCF")
        ret = response.scf_response.cpscf_linear_response(scf_wfn, *linear_response,
                                                            conv_tol = core.get_global_option("SOLVER_CONVERGENCE"),
                                                            max_iter = core.get_global_option("SOLVER_MAXITER"),
                                                            print_lvl = (core.get_global_option("PRINT") + 1))
        optstash_jk.restore()

    core.tstop()
    optstash.restore()
    return scf_wfn


def run_cc_property(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    all CC property calculations.

    """
    optstash = p4util.OptionsState(
            ['WFN'],
            ['DERTYPE'],
            ['ONEPDM'],
            ['PROPERTY'],
            ['CCLAMBDA', 'R_CONVERGENCE'],
            ['CCEOM', 'R_CONVERGENCE'],
            ['CCEOM', 'E_CONVERGENCE']) # yapf:disable

    oneel_properties = core.OEProp.valid_methods
    twoel_properties = []
    response_properties = ['POLARIZABILITY', 'ROTATION', 'ROA', 'ROA_TENSOR']
    excited_properties = ['OSCILLATOR_STRENGTH', 'ROTATIONAL_STRENGTH']

    one = []
    two = []
    response = []
    excited = []
    invalid = []

    if 'properties' in kwargs:
        properties = kwargs['properties']

        for prop in properties:

            prop = prop.upper()
            if prop in oneel_properties:
                one.append(prop)
            elif prop in twoel_properties:
                two.append(prop)
            elif prop in response_properties:
                response.append(prop)
            elif prop in excited_properties:
                excited.append(prop)
            else:
                invalid.append(prop)
    else:
        raise ValidationError("""The "properties" keyword is required with the property() function.""")

    # People are used to requesting dipole/quadrupole and getting dipole,quadrupole,mulliken_charges and NO_occupations
    if ('DIPOLE' in one) or ('QUADRUPOLE' in one):
        one = list(set(one + ['DIPOLE', 'QUADRUPOLE', 'MULLIKEN_CHARGES', 'NO_OCCUPATIONS']))

    n_one = len(one)
    n_two = len(two)
    n_response = len(response)
    n_excited = len(excited)
    n_invalid = len(invalid)

    if n_invalid > 0:
        print("""The following properties are not currently supported: %s""" % invalid)

    if n_excited > 0 and (name not in ['eom-ccsd', 'eom-cc2']):
        raise ValidationError("""Excited state CC properties require EOM-CC2 or EOM-CCSD.""")

    if (name in ['eom-ccsd', 'eom-cc2']) and n_response > 0:
        raise ValidationError("""Cannot (yet) compute response properties for excited states.""")

    if 'roa' in response:
        # Perform distributed roa job
        run_roa(name, **kwargs)
        return  # Don't do anything further

    if (n_one > 0 or n_two > 0) and (n_response > 0):
        print("""Computing both density- and response-based properties.""")

    if n_response > 0:
        if ("ref_wfn" in kwargs and not kwargs["ref_wfn"].same_a_b_orbs()) or core.get_option('SCF', 'REFERENCE') != 'RHF':
            raise ValidationError(f"Non-RHF CC response properties are not implemented.")

    if name in ['ccsd', 'cc2', 'eom-ccsd', 'eom-cc2']:
        this_name = name.upper().replace('-', '_')
        core.set_global_option('WFN', this_name)
        ccwfn = run_ccenergy(name, **kwargs)
        core.set_global_option('WFN', this_name)
    else:
        raise ValidationError(f"CC property name {name.upper()} not recognized")

    # Need cchbar for everything
    core.cchbar(ccwfn)

    # Need ccdensity at this point only for density-based props
    if n_one > 0 or n_two > 0:
        if name == 'eom-ccsd':
            core.set_global_option('WFN', 'EOM_CCSD')
            core.set_global_option('DERTYPE', 'NONE')
            core.cceom(ccwfn)
        elif name == 'eom-cc2':
            core.set_global_option('WFN', 'EOM_CC2')
            core.set_global_option('DERTYPE', 'NONE')
            core.cceom(ccwfn)
        core.set_global_option('DERTYPE', 'NONE')
        if core.get_option('CCDENSITY', 'OPDM_RELAX') or n_two > 0:
            # WARNING!!! A one-particle property computed _with_ a two-particle property will differ
            # from a one-particle property computed by itself. There are no two-particle properties at
            # present, so we can kick the issue further down the road.
            core.set_global_option('OPDM_ONLY', 'FALSE')
        else:
            core.set_global_option('OPDM_ONLY', 'TRUE')
        core.cclambda(ccwfn)
        core.ccdensity(ccwfn)

    # Need ccresponse only for response-type props
    if n_response > 0:
        core.set_global_option('DERTYPE', 'RESPONSE')
        core.cclambda(ccwfn)
        for prop in response:
            core.set_global_option('PROPERTY', prop)
            core.ccresponse(ccwfn)

    # Excited-state transition properties
    if n_excited > 0:
        if name == 'eom-ccsd':
            core.set_global_option('WFN', 'EOM_CCSD')
        elif name == 'eom-cc2':
            core.set_global_option('WFN', 'EOM_CC2')
        else:
            raise ValidationError("""Unknown excited-state CC wave function.""")
        core.set_global_option('DERTYPE', 'NONE')
        if core.get_option('CCDENSITY', 'OPDM_RELAX'):
            core.set_global_option('OPDM_ONLY', 'FALSE')
        else:
            core.set_global_option('OPDM_ONLY', 'TRUE')
        # Tight convergence unnecessary for transition properties
        core.set_local_option('CCLAMBDA', 'R_CONVERGENCE', 1e-4)
        core.set_local_option('CCEOM', 'R_CONVERGENCE', 1e-4)
        core.set_local_option('CCEOM', 'E_CONVERGENCE', 1e-5)
        core.cceom(ccwfn)
        core.cclambda(ccwfn)
        core.ccdensity(ccwfn)

    # => Make OEProp calls <=
    if n_one > 0:
        # ==> Initialize OEProp  <==
        oe = core.OEProp(ccwfn)
        for oe_prop_name in one:
            oe.add(oe_prop_name.upper())
        # ==> OEProp for the ground state <==
        # TODO: When Psi is Py 3.9+, transition to the removeprefix version.
        title = name.upper().replace("EOM-", "")
        #title = name.upper().removeprefix("EOM-")
        oe.set_title(title)
        set_of_names = {title + " {}", "CC {}"}
        if name.startswith("eom"):
            gs_h = 0
            for h, i in enumerate(ccwfn.soccpi()):
                if i % 2:
                    gs_h = gs_h ^ h
            ct = ccwfn.molecule().point_group().char_table()
            total_h_lbl = ct.gamma(0).symbol()
            gs_h_lbl = ct.gamma(gs_h).symbol()
            set_of_names.update({title + " ROOT 0 {}", "CC ROOT 0 {}",
                                 f"{title} ROOT 0 {{}} - {total_h_lbl} TRANSITION",
                                      f"CC ROOT 0 {{}} - {total_h_lbl} TRANSITION",
                                 f"{title} ROOT 0 ({gs_h_lbl}) {{}}", f"CC ROOT 0 ({gs_h_lbl}) {{}}",
                                 f"{title} ROOT 0 (IN {gs_h_lbl}) {{}}", f"CC ROOT 0 (IN {gs_h_lbl}) {{}}"})
        oe.set_names(set_of_names)
        oe.compute()

        # ==> OEProp for Excited States <==
        if name.startswith('eom'):
            n_root_pi = core.get_global_option("ROOTS_PER_IRREP")
            for h in range(ccwfn.nirrep()):
                root_h_lbl = ct.gamma(h).symbol()
                trans_h_lbl = ct.gamma(h ^ gs_h).symbol()
                # Don't forget to count the ground state!
                for i in range(n_root_pi[h]):
                    if h == gs_h: i += 1
                    root_title = title + f" ROOT {i} (IN {root_h_lbl})"
                    oe.set_title(root_title)
                    total_idx = ccwfn.total_index(i, h)
                    set_of_names = {f"{title} ROOT {total_idx} {{}}", f"CC ROOT {total_idx} {{}}",
                                    f"{title} ROOT {total_idx} {{}} - {trans_h_lbl} TRANSITION",
                                         f"CC ROOT {total_idx} {{}} - {trans_h_lbl} TRANSITION",
                                    f"{title} ROOT {total_idx} ({root_h_lbl}) {{}}", f"CC ROOT {total_idx} ({root_h_lbl}) {{}}",
                                    f"{title} ROOT {i} (IN {root_h_lbl}) {{}}", f"CC ROOT {i} (IN {root_h_lbl}) {{}}"}
                    oe.set_names(set_of_names)
                    Da = ccwfn.get_density(root_title + " ALPHA")
                    oe.set_Da_so(Da)
                    if not ccwfn.same_a_b_dens():
                        Db = ccwfn.get_density(root_title + " BETA")
                        oe.set_Db_so(Db)
                    oe.compute()

    core.set_global_option('WFN', 'SCF')
    core.revoke_global_option_changed('WFN')
    core.set_global_option('DERTYPE', 'NONE')
    core.revoke_global_option_changed('DERTYPE')

    optstash.restore()
    return ccwfn


def run_dfmp2_property(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a DFMP2 property calculation.

    """
    optstash = p4util.OptionsState(
        ['DF_BASIS_SCF'],
        ['DF_BASIS_MP2'],
        ['ONEPDM'],
        ['OPDM_RELAX'],
        ['SCF_TYPE'])

    core.set_global_option('ONEPDM', 'TRUE')
    core.set_global_option('OPDM_RELAX', 'TRUE')

    # Alter default algorithm
    if not core.has_global_option_changed('SCF_TYPE'):
        core.set_global_option('SCF_TYPE', 'DF')  # local set insufficient b/c SCF option read in DFMP2
        core.print_out("""    SCF Algorithm Type (re)set to DF.\n""")

    if 'DF' not in core.get_global_option('SCF_TYPE'):
        raise ValidationError('DF-MP2 properties need DF-SCF reference.')

    properties = kwargs.pop('properties')
    proc_util.oeprop_validator(properties)

    # Bypass the scf call if a reference wavefunction is given
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, scf_do_properties=False, use_c1=True, **kwargs)  # C1 certified

    aux_basis = core.BasisSet.build(ref_wfn.molecule(), "DF_BASIS_MP2",
                                    core.get_option("DFMP2", "DF_BASIS_MP2"),
                                    "RIFIT", core.get_global_option('BASIS'))
    ref_wfn.set_basisset("DF_BASIS_MP2", aux_basis)

    core.tstart()
    core.print_out('\n')
    p4util.banner('DFMP2')
    core.print_out('\n')

    dfmp2_wfn = core.dfmp2(ref_wfn)
    grad = dfmp2_wfn.compute_gradient()

    if name == 'scs-mp2':
        dfmp2_wfn.set_variable('CURRENT ENERGY', dfmp2_wfn.variable('SCS-MP2 TOTAL ENERGY'))
        dfmp2_wfn.set_variable('CURRENT CORRELATION ENERGY', dfmp2_wfn.variable('SCS-MP2 CORRELATION ENERGY'))
    elif name == 'mp2':
        dfmp2_wfn.set_variable('CURRENT ENERGY', dfmp2_wfn.variable('MP2 TOTAL ENERGY'))
        dfmp2_wfn.set_variable('CURRENT CORRELATION ENERGY', dfmp2_wfn.variable('MP2 CORRELATION ENERGY'))

    # Run OEProp
    oe = core.OEProp(dfmp2_wfn)
    oe.set_title(name.upper())
    for prop in properties:
        oe.add(prop.upper())
    oe.compute()
    dfmp2_wfn.oeprop = oe

    # Shove variables into global space
    for k, v in dfmp2_wfn.variables().items():
        core.set_variable(k, v)

    optstash.restore()
    core.tstop()
    return dfmp2_wfn


def _clean_detci(keep: bool=True):
    psioh = core.IOManager.shared_object()
    psio = core.IO.shared_object()
    cifl = core.get_option("DETCI", "CI_FILE_START")
    for fl in range(cifl, cifl + 4):
        if psio.open_check(fl):
            psio.close(fl, keep)


def run_detci_property(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a configuration interaction calculation, namely FCI,
    CIn, MPn, and ZAPTn, computing properties.

    """
    optstash = p4util.OptionsState(
        ['OPDM'],
        ['TDM'])

    # Find valid properties
    valid_transition = ['TRANSITION_DIPOLE', 'TRANSITION_QUADRUPOLE']

    ci_prop = []
    ci_trans = []
    properties = kwargs.pop('properties')
    for prop in properties:
        if prop.upper() in valid_transition:
            ci_trans.append(prop)
        else:
            ci_prop.append(prop)

    proc_util.oeprop_validator(ci_prop)

    core.set_global_option('OPDM', 'TRUE')
    if len(ci_trans):
        core.set_global_option('TDM', 'TRUE')

    # Compute
    if name in ['mcscf', 'rasscf', 'casscf']:
        ciwfn = run_detcas(name, **kwargs)
    else:
        ciwfn = run_detci(name, **kwargs)

    # All property names are just CI
    if 'CI' in name.upper():
        name = 'CI'

    states = core.get_global_option('avg_states')
    nroots = core.get_global_option('num_roots')
    if len(states) != nroots:
        states = range(nroots)

    # Run OEProp
    oe = core.OEProp(ciwfn)
    oe.set_title(name.upper())
    for prop in ci_prop:
        oe.add(prop.upper())

    # Compute "the" CI density
    oe.compute()
    ciwfn.oeprop = oe

    # If we have more than one root, compute all data
    if nroots > 1:
        core.print_out("\n   ===> %s properties for all CI roots <=== \n\n" % name.upper())
        for root in states:
            oe.set_title("%s ROOT %d" % (name.upper(), root))
            if ciwfn.same_a_b_dens():
                oe.set_Da_mo(ciwfn.get_opdm(root, root, "A", True))
            else:
                oe.set_Da_mo(ciwfn.get_opdm(root, root, "A", True))
                oe.set_Db_mo(ciwfn.get_opdm(root, root, "B", True))
            oe.compute()

    # Transition density matrices
    if (nroots > 1) and len(ci_trans):
        oe.clear()
        for tprop in ci_trans:
            oe.add(tprop.upper())

        core.print_out("\n   ===> %s properties for all CI transition density matrices <=== \n\n" % name.upper())
        for root in states[1:]:
            oe.set_title("%s ROOT %d -> ROOT %d" % (name.upper(), 0, root))
            if ciwfn.same_a_b_dens():
                oe.set_Da_mo(ciwfn.get_opdm(0, root, "A", True))
            else:
                oe.set_Da_mo(ciwfn.get_opdm(0, root, "A", True))
                oe.set_Db_mo(ciwfn.get_opdm(0, root, "B", True))
            oe.compute()

    _clean_detci()
    optstash.restore()
    return ciwfn


def run_eom_cc(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    an EOM-CC calculation, namely EOM-CC2, EOM-CCSD, and EOM-CC3.

    """
    optstash = p4util.OptionsState(
        ['TRANSQT2', 'WFN'],
        ['CCSORT', 'WFN'],
        ['CCENERGY', 'WFN'],
        ['CCHBAR', 'WFN'],
        ['CCEOM', 'WFN'])

    if name == 'eom-ccsd':
        core.set_local_option('TRANSQT2', 'WFN', 'EOM_CCSD')
        core.set_local_option('CCSORT', 'WFN', 'EOM_CCSD')
        core.set_local_option('CCENERGY', 'WFN', 'EOM_CCSD')
        core.set_local_option('CCHBAR', 'WFN', 'EOM_CCSD')
        core.set_local_option('CCEOM', 'WFN', 'EOM_CCSD')
        ref_wfn = run_ccenergy('ccsd', **kwargs)

    elif name == 'eom-cc2':
        user_ref = core.get_option('CCENERGY', 'REFERENCE')
        if (user_ref != 'RHF') and (user_ref != 'UHF'):
            raise ValidationError('Reference %s for EOM-CC2 is not available.' % user_ref)
        core.set_local_option('TRANSQT2', 'WFN', 'EOM_CC2')
        core.set_local_option('CCSORT', 'WFN', 'EOM_CC2')
        core.set_local_option('CCENERGY', 'WFN', 'EOM_CC2')
        core.set_local_option('CCHBAR', 'WFN', 'EOM_CC2')
        core.set_local_option('CCEOM', 'WFN', 'EOM_CC2')
        ref_wfn = run_ccenergy('cc2', **kwargs)

    elif name == 'eom-cc3':
        core.set_local_option('TRANSQT2', 'WFN', 'EOM_CC3')
        core.set_local_option('CCSORT', 'WFN', 'EOM_CC3')
        core.set_local_option('CCENERGY', 'WFN', 'EOM_CC3')
        core.set_local_option('CCHBAR', 'WFN', 'EOM_CC3')
        core.set_local_option('CCEOM', 'WFN', 'EOM_CC3')
        ref_wfn = run_ccenergy('cc3', **kwargs)

    core.cchbar(ref_wfn)
    core.cceom(ref_wfn)

    optstash.restore()
    return ref_wfn
    # TODO ask if all these cc modules not actually changing wfn


def run_eom_cc_gradient(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    an EOM-CCSD gradient calculation.

    """
    optstash = p4util.OptionsState(
        ['CCDENSITY', 'XI'],
        ['CCDENSITY', 'ZETA'],
        ['CCLAMBDA', 'ZETA'],
        ['DERTYPE'],
        ['CCDENSITY', 'WFN'],
        ['CCLAMBDA', 'WFN'])

    core.set_global_option('DERTYPE', 'FIRST')

    if name == 'eom-ccsd':
        core.set_local_option('CCLAMBDA', 'WFN', 'EOM_CCSD')
        core.set_local_option('CCDENSITY', 'WFN', 'EOM_CCSD')
        ref_wfn = run_eom_cc(name, **kwargs)
    else:
        core.print_out('DGAS: proc.py:1599 hitting an undefined sequence')
        core.clean()
        raise ValueError('Hit a wall in proc.py:1599')

    core.set_local_option('CCLAMBDA', 'ZETA', 'FALSE')
    core.set_local_option('CCDENSITY', 'ZETA', 'FALSE')
    core.set_local_option('CCDENSITY', 'XI', 'TRUE')
    core.cclambda(ref_wfn)
    core.ccdensity(ref_wfn)
    core.set_local_option('CCLAMBDA', 'ZETA', 'TRUE')
    core.set_local_option('CCDENSITY', 'ZETA', 'TRUE')
    core.set_local_option('CCDENSITY', 'XI', 'FALSE')
    core.cclambda(ref_wfn)
    core.ccdensity(ref_wfn)

    derivobj = core.Deriv(ref_wfn)
    grad = derivobj.compute()

    ref_wfn.set_gradient(grad)

    optstash.restore()
    return ref_wfn


def run_adcc(name, **kwargs):
    """Prepare and run an ADC calculation in adcc, interpret the result and return
    as a wavefunction.

    """
    # TODO Maybe it would improve readability if this function was spilt
    #      up and the whole thing went to a separate file (like for sapt,
    #      interface_cfour.py, ...

    try:
        import adcc
        from adcc.exceptions import InvalidReference
    except ModuleNotFoundError:
        raise ValidationError("adcc extras qc_module not available. Try installing "
            "via 'pip install adcc' or 'conda install -c adcc adcc'.")


    if core.get_option('ADC', 'REFERENCE') not in ["RHF", "UHF"]:
        raise ValidationError('adcc requires reference RHF or UHF')

    # Bypass the scf call if a reference wavefunction is given
    ref_wfn = kwargs.pop('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, use_c1=True, **kwargs)

    # Start timer
    do_timer = kwargs.pop("do_timer", True)
    if do_timer:
        core.tstart()

    #
    # Build kwargs for adcc
    #
    kwargs.pop("molecule", None)

    if ref_wfn.frzcpi()[0] > 0:
        kwargs["frozen_core"] = ref_wfn.frzcpi()[0]
    if ref_wfn.frzvpi()[0] > 0:
        kwargs["frozen_virtual"] = ref_wfn.frzvpi()[0]
    if core.get_option("ADC", "NUM_CORE_ORBITALS"):
        kwargs["core_orbitals"] = core.get_option("ADC", "NUM_CORE_ORBITALS")

    scf_accuracy = max(core.get_option("SCF", "E_CONVERGENCE"),
                       core.get_option("SCF", "D_CONVERGENCE"))
    if core.get_option("ADC", "R_CONVERGENCE") < 0:
        kwargs["conv_tol"] = max(100 * scf_accuracy, 1e-6)
    else:
        kwargs["conv_tol"] = core.get_option("ADC", "R_CONVERGENCE")

    n_roots = core.get_option('ADC', 'ROOTS_PER_IRREP')
    if len(n_roots) > 1:
        raise ValidationError("adcc can only deal with a single irrep.")
    kwargs["n_states"] = n_roots[0]

    if core.get_option("ADC", "NUM_GUESSES") > 0:
        kwargs["n_guesses"] = core.get_option("ADC", "NUM_GUESSES")
    if core.get_option("ADC", "MAX_NUM_VECS") > 0:
        kwargs["max_subspace"] = core.get_option("ADC", "MAX_NUM_VECS")

    kind = core.get_option("ADC", "KIND").lower()
    if isinstance(ref_wfn, core.UHF):
        if not core.has_option_changed("ADC", "KIND"):
            kind = "any"
        elif kind not in ["any", "spin_flip"]:
            raise ValidationError("For UHF references the only valid values for 'KIND' are "
                                  "'SPIN_FLIP' or 'ANY' and not '{}.".format(kind.upper()))
    elif kind not in ["singlet", "triplet", "any"]:
        raise ValidationError("For RHF references the value '{}' for 'KIND' is "
                              "not supported.".format(kind.upper()))
    kwargs["kind"] = kind
    kwargs["max_iter"] = core.get_option("ADC", "MAXITER")

    #
    # Determine ADC function method from adcc to run ADC
    #
    adcrunner = {
        "cvs-adc(1)": adcc.cvs_adc1, "cvs-adc(2)": adcc.cvs_adc2,
        "cvs-adc(2)-x": adcc.cvs_adc2x, "cvs-adc(3)": adcc.cvs_adc3,
        "adc(1)": adcc.adc1, "adc(2)": adcc.adc2,
        "adc(2)-x": adcc.adc2x, "adc(3)": adcc.adc3,
    }
    if name not in adcrunner:
        raise ValidationError(f"Unsupported ADC method: {name}")
    if "cvs" in name and "core_orbitals" not in kwargs:
        raise ValidationError("If a CVS-ADC method is requested, the NUM_CORE_ORBITALS option "
                              "needs to be set.")
    if "core_orbitals" in kwargs and "cvs" not in name:
        raise ValidationError("The NUM_CORE_ORBITALS option needs to be set to '0' or absent "
                              "unless a CVS ADC method is requested.")
    if "cvs" in name and kwargs["kind"] in ["spin_flip"]:
        raise ValidationError("Spin-flip for CVS-ADC variants is not available.")

    #
    # Launch the rocket
    #
    # Copy thread setup from psi4
    adcc.set_n_threads(core.get_num_threads())

    # Hack to direct the stream-like interface adcc expects to the string interface of Psi4 core
    class CoreStream:
        def write(self, text):
            core.print_out(text)

    core.print_out("\n" + adcc.banner(colour=False) + "\n")
    try:
        state = adcrunner[name](ref_wfn, **kwargs, output=CoreStream())
    except InvalidReference as ex:
        raise ValidationError("Cannot run adcc because the passed reference wavefunction is "
                              "not supported in adcc. Check Psi4 SCF parameters. adcc reports: "
                              f"{ex}")
    except Exception as ex:
        raise ValidationError("Unknown exception occured while "
                              f"running adcc: '{ex}' ({type(ex).__name__})")
    core.print_out("\n")

    # TODO Should a non-converged calculation throw?

    #
    # Interpret results
    #
    # Note: This wavefunction is not consistent ... the density
    # is e.g. not the proper one (i.e. not the MP(n) one)
    adc_wfn = core.Wavefunction(ref_wfn.molecule(), ref_wfn.basisset())
    adc_wfn.shallow_copy(ref_wfn)
    adc_wfn.set_reference_wavefunction(ref_wfn)
    adc_wfn.set_name(name)
    adc_wfn.set_module("adcc")

    # MP(3) energy for CVS-ADC(3) calculations is still a missing feature in adcc
    # ... we store this variant here to be able to fall back to MP(2) energies.
    is_cvs_adc3 = state.method.level >= 3 and state.ground_state.has_core_occupied_space

    # Ground-state energies
    mp = state.ground_state
    mp_energy = mp.energy(state.method.level if not is_cvs_adc3 else 2)
    mp_corr = 0.0
    if state.method.level > 1:
        core.print_out("Ground state energy breakdown:\n")
        core.print_out("    Energy             SCF   {0:15.8g} [Eh]\n".format(ref_wfn.energy()))
        for level in range(2, state.method.level + 1):
            if level >= 3 and is_cvs_adc3:
                continue
            energy = mp.energy_correction(level)
            mp_corr += energy
            adc_wfn.set_variable(f"MP{level} CORRELATION ENERGY", energy)
            adc_wfn.set_variable(f"MP{level} TOTAL ENERGY", mp.energy(level))
            core.print_out(f"    Energy correlation MP{level}   {energy:15.8g} [Eh]\n")
        core.print_out("    Energy             total {0:15.8g} [Eh]\n".format(mp_energy))
    adc_wfn.set_variable("CURRENT CORRELATION ENERGY", mp_corr)  # P::e ADC
    adc_wfn.set_variable("CURRENT ENERGY", mp_energy)  # P::e ADC

    # Set results of excited-states computation
    # TODO Does not work: Can't use strings
    # adc_wfn.set_variable("excitation kind", state.kind)
    adc_wfn.set_variable("number of excited states", len(state.excitation_energy))
    adc_wfn.set_variable("ADC ITERATIONS", state.n_iter)  # P::e ADC
    methods = [name.upper(), 'ADC']
    for excitation in state.excitations:
        root_index = excitation.index + 1
        for method in methods:
            adc_wfn.set_variable(f"{method} ROOT 0 (A) -> ROOT {root_index} (A) EXCITATION ENERGY",
                                 excitation.excitation_energy)
            adc_wfn.set_variable(f"{method} ROOT 0 (IN A) -> ROOT {root_index} (IN A) EXCITATION ENERGY",
                                 excitation.excitation_energy)
            adc_wfn.set_variable(f"{method} ROOT 0 -> ROOT {root_index} EXCITATION ENERGY",
                                 excitation.excitation_energy)
            adc_wfn.set_variable(f"{method} ROOT 0 -> ROOT {root_index} EXCITATION ENERGY - A TRANSITION",
                                 excitation.excitation_energy)

    core.print_out("\n\n  ==> Excited states summary <==  \n")
    core.print_out("\n" + state.describe(oscillator_strengths=False) + "\n")

    # TODO Setting the excitation amplitude elements inside the wavefunction is a little
    #      challenging, since for each excitation vector one needs to extract the elements
    #      and map the indices from the adcc to the Psi4 convention. For this reason it
    #      is not yet done.

    core.print_out("\n  ==> Dominant amplitudes per state <==  \n\n")
    tol_ampl = core.get_option("ADC", "CUTOFF_AMPS_PRINT")
    core.print_out(state.describe_amplitudes(tolerance=tol_ampl) + "\n\n")

    # Shove variables into global space
    for k, v in adc_wfn.variables().items():
        core.set_variable(k, v)

    if do_timer:
        core.tstop()
    adc_wfn.adcc_state = state
    return adc_wfn


def run_adcc_property(name, **kwargs):
    """Run a ADC excited-states property calculation in adcc
    and return the resulting properties.

    """
    # TODO Things available in ADCC, but not yet implemented here:
    #      Export of difference and transition density matrices for all states

    properties = [prop.upper() for prop in kwargs.pop('properties')]
    valid_properties = ['DIPOLE', 'OSCILLATOR_STRENGTH', 'TRANSITION_DIPOLE',
                        'ROTATIONAL_STRENGTH']
    unknown_properties = [prop for prop in properties if prop not in valid_properties]

    if unknown_properties:
        alternatives = ""
        alt_method_name = p4util.text.find_approximate_string_matches(unknown_properties[0],
                                                                      valid_properties, 2)
        if alt_method_name:
            alternatives = " Did you mean? " + " ".join(alt_method_name)

        raise ValidationError("ADC property: Feature '{}' is not recognized. {}"
                              "".format(unknown_properties[0], alternatives))

    # Start timer
    do_timer = kwargs.pop("do_timer", True)
    if do_timer:
        core.tstart()
    adc_wfn = run_adcc(name, do_timer=False, **kwargs)
    state = adc_wfn.adcc_state
    hf = state.reference_state
    mp = state.ground_state

    # Formats and indention
    ind = "    "
    def format_vector(label, data):
        assert data.ndim == 1
        return f"{label:<40s} " + " ".join(f"{d:12.6g}" for d in data)

    if "DIPOLE" in properties:
        lines = ["\nGround state properties"]
        lines += [ind + "Hartree-Fock (HF)"]
        lines += [ind + ind + format_vector("Dipole moment (in a.u.)", hf.dipole_moment)]

        if state.method.level > 1:
            lines += [ind + "Møller Plesset 2nd order (MP2)"]
            lines += [ind + ind + format_vector("Dipole moment (in a.u.)", mp.dipole_moment(2))]
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                for i, cart in enumerate(["X", "Y", "Z"]):
                    # retire components at v1.5
                    adc_wfn.set_variable("MP2 dipole " + cart, mp.dipole_moment(2)[i])
                adc_wfn.set_variable("current dipole " + cart, mp.dipole_moment(2)[i])
            adc_wfn.set_variable("MP2 dipole", mp.dipole_moment(2))
            adc_wfn.set_variable("current dipole", mp.dipole_moment(2))
        lines += [""]
        core.print_out("\n".join(lines) + "\n")

    gauge = core.get_option("ADC", "GAUGE").lower()
    if gauge == "velocity":
        gauge_short = "VEL"
    elif gauge == "length":
        gauge_short = "LEN"
    else:
        raise ValidationError(f"Gauge {gauge} not recognised for ADC calculations.")

    computed = []
    methods = [name.upper(), 'ADC']
    for excitation in state.excitations:
        root_index = excitation.index + 1
        props = {}
        if any(prop in properties for prop in ("TRANSITION_DIPOLE", "OSCILLATOR_STRENGTH")):
            data = excitation.transition_dipole_moment
            props["Transition dipole moment (in a.u.)"] = data
            data_mat = data.reshape(1, 3)
            for method in methods:
                adc_wfn.set_variable(f"{method} ROOT 0 (A) -> ROOT {root_index} (A) "
                                    "ELECTRIC TRANSITION DIPOLE MOMENT (LEN)", data_mat)
                adc_wfn.set_variable(f"{method} ROOT 0 (IN A) -> ROOT {root_index} (IN A) "
                                    "ELECTRIC TRANSITION DIPOLE MOMENT (LEN)", data_mat)
                adc_wfn.set_variable(f"{method} ROOT 0 -> ROOT {root_index} "
                                    "ELECTRIC TRANSITION DIPOLE MOMENT (LEN)", data_mat)
                adc_wfn.set_variable(f"{method} ROOT 0 -> ROOT {root_index} "
                                    "ELECTRIC TRANSITION DIPOLE MOMENT (LEN) - A TRANSITION", data_mat)

        if "OSCILLATOR_STRENGTH" in properties:
            if gauge == "velocity":
                data = excitation.oscillator_strength_velocity
            else:
                data = excitation.oscillator_strength
            props[f"Oscillator strength ({gauge} gauge)"] = data.reshape(-1)
            for method in methods:
                adc_wfn.set_variable(f"{method} ROOT 0 (A) -> ROOT {root_index} (A) "
                                    f"OSCILLATOR STRENGTH ({gauge_short})", data)
                adc_wfn.set_variable(f"{method} ROOT 0 (IN A) -> ROOT {root_index} (IN A) "
                                    f"OSCILLATOR STRENGTH ({gauge_short})", data)
                adc_wfn.set_variable(f"{method} ROOT 0 -> ROOT {root_index} "
                                    f"OSCILLATOR STRENGTH ({gauge_short})", data)
                adc_wfn.set_variable(f"{method} ROOT 0 -> ROOT {root_index} "
                                    f"OSCILLATOR STRENGTH ({gauge_short}) - A TRANSITION", data)

        if "ROTATIONAL_STRENGTH" in properties:
            data = excitation.rotatory_strength
            props["Rotational strength (velocity gauge)"] = data.reshape(-1)
            for method in methods:
                adc_wfn.set_variable(f"{method} ROOT 0 (A) -> ROOT {root_index} (A) "
                                    "ROTATORY STRENGTH (VEL)", data)
                adc_wfn.set_variable(f"{method} ROOT 0 (IN A) -> ROOT {root_index} (IN A) "
                                    "ROTATORY STRENGTH (VEL)", data)
                adc_wfn.set_variable(f"{method} ROOT 0 -> ROOT {root_index} "
                                    "ROTATORY STRENGTH (VEL)", data)
                adc_wfn.set_variable(f"{method} ROOT 0 -> ROOT {root_index} "
                                    "ROTATORY STRENGTH (VEL) - A TRANSITION", data)

        if "DIPOLE" in properties:
            data = excitation.state_dipole_moment
            props["State dipole moment (in a.u.)"] = data
            data_mat = data.reshape(1, 3)
            for method in methods:
                adc_wfn.set_variable(f"{method} ROOT {root_index} DIPOLE MOMENT", data_mat)
                adc_wfn.set_variable(f"{method} ROOT {root_index} DIPOLE MOMENT - A TRANSITION", data_mat)
                adc_wfn.set_variable(f"{method} ROOT {root_index} (A) DIPOLE MOMENT", data_mat)
                adc_wfn.set_variable(f"{method} ROOT {root_index} (IN A) DIPOLE MOMENT", data_mat)

        computed.append(props)

        # for Psivar scraper
        # wfn.set_variable("ADC ROOT 0 -> ROOT n EXCITATION ENERGY")               # P::e ADC
        # wfn.set_variable("ADC ROOT 0 (IN h) -> ROOT n (IN i) EXCITATION ENERGY")       # P::e ADC
        # wfn.set_variable("ADC ROOT 0 (h) -> ROOT n (i) EXCITATION ENERGY")       # P::e ADC
        # wfn.set_variable("ADC ROOT 0 -> ROOT n EXCITATION ENERGY - h TRANSITION")  # P::e ADC
        # wfn.set_array_variable("ADC ROOT 0 -> ROOT n ELECTRIC TRANSITION DIPOLE MOMENT (LEN)")                # P::e ADC
        # wfn.set_array_variable("ADC ROOT 0 (IN h) -> ROOT n (IN i) ELECTRIC TRANSITION DIPOLE MOMENT (LEN)")        # P::e ADC
        # wfn.set_array_variable("ADC ROOT 0 (h) -> ROOT n (i) ELECTRIC TRANSITION DIPOLE MOMENT (LEN)")        # P::e ADC
        # wfn.set_array_variable("ADC ROOT 0 -> ROOT n ELECTRIC TRANSITION DIPOLE MOMENT (LEN) - h TRANSITION")  # P::e ADC
        # wfn.set_array_variable("ADC ROOT n DIPOLE MOMENT")                # P::e ADC
        # wfn.set_array_variable("ADC ROOT n (IN i) DIPOLE MOMENT")                # P::e ADC
        # wfn.set_array_variable("ADC ROOT n (i) DIPOLE MOMENT")                # P::e ADC
        # wfn.set_array_variable("ADC ROOT n DIPOLE MOMENT - h TRANSITION")                # P::e ADC
        # wfn.set_variable("ADC ROOT 0 -> ROOT n OSCILLATOR STRENGTH (LEN)")               # P::e ADC
        # wfn.set_variable("ADC ROOT 0 (IN h) -> ROOT n (IN i) OSCILLATOR STRENGTH (LEN)")       # P::e ADC
        # wfn.set_variable("ADC ROOT 0 (h) -> ROOT n (i) OSCILLATOR STRENGTH (LEN)")       # P::e ADC
        # wfn.set_variable("ADC ROOT 0 -> ROOT n OSCILLATOR STRENGTH (LEN) - h TRANSITION")  # P::e ADC
        # wfn.set_variable("ADC ROOT 0 -> ROOT n OSCILLATOR STRENGTH (VEL)")               # P::e ADC
        # wfn.set_variable("ADC ROOT 0 (IN h) -> ROOT n (IN i) OSCILLATOR STRENGTH (VEL)")       # P::e ADC
        # wfn.set_variable("ADC ROOT 0 (h) -> ROOT n (i) OSCILLATOR STRENGTH (VEL)")       # P::e ADC
        # wfn.set_variable("ADC ROOT 0 -> ROOT n OSCILLATOR STRENGTH (VEL) - h TRANSITION")  # P::e ADC
        # wfn.set_variable("ADC ROOT 0 -> ROOT n ROTATORY STRENGTH (VEL)")               # P::e ADC
        # wfn.set_variable("ADC ROOT 0 (IN h) -> ROOT n (IN i) ROTATORY STRENGTH (VEL)")       # P::e ADC
        # wfn.set_variable("ADC ROOT 0 (h) -> ROOT n (i) ROTATORY STRENGTH (VEL)")       # P::e ADC
        # wfn.set_variable("ADC ROOT 0 -> ROOT n ROTATORY STRENGTH (VEL) - h TRANSITION")  # P::e ADC

    core.print_out("\nExcited state properties:\n")
    for i, props in enumerate(computed):
        lines = [ind + f"Excited state  {i}"]
        for prop, data in sorted(props.items()):
            lines += [ind + ind + format_vector(prop, data)]
        core.print_out("\n".join(lines) + "\n")

    # Shove variables into global space
    for k, v in adc_wfn.variables().items():
        core.set_variable(k, v)

    if do_timer:
        core.tstop()
    return adc_wfn


def run_detci(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a configuration interaction calculation, namely FCI,
    CIn, MPn, and ZAPTn.

    """
    # dtl = docs_table_link("dummy", "detci")

    optstash = p4util.OptionsState(
        ['DETCI', 'WFN'],
        ['DETCI', 'MAX_NUM_VECS'],
        ['DETCI', 'MPN_ORDER_SAVE'],
        ['DETCI', 'MPN'],
        ['DETCI', 'FCI'],
        ['DETCI', 'EX_LEVEL'])

    # throw exception for UHF
    if core.get_option('DETCI', 'REFERENCE') not in ['RHF', 'ROHF']:
        raise ValidationError('Reference %s for DETCI is not available.' %
            core.get_option('DETCI', 'REFERENCE'))

    # throw exception for DF/CD. many of these pre-trapped by select_* functions but some escape, incl. zapt
    if (corl_type := method_algorithm_type(name).now) != "CONV":
        raise ValidationError(f"Invalid type {corl_type} for DETCI energy through `run_detci`.")  # See Capabilities Table")

    mtdlvl_mobj = re.match(r"""\A(?P<method>[a-z]+)(?P<level>\d+)\Z""", name.lower())

    if mtdlvl_mobj and mtdlvl_mobj.group("method") == "zapt":
        level = int(mtdlvl_mobj.group("level"))

        # throw exception for non-ROHF
        if (ref := core.get_option("SCF", "REFERENCE")) != "ROHF":
            raise UpgradeHelper(f"energy('zapt{level}')", f"energy('mp{level}')", 1.7,
                    " Replace method ZAPT with method MP for RHF reference. DETCI is orders-of-magnitude inefficient for perturbation theory.")

        core.set_local_option('DETCI', 'WFN', 'ZAPTN')
        maxnvect = int((level + 1) / 2) + (level + 1) % 2
        core.set_local_option('DETCI', 'MAX_NUM_VECS', maxnvect)
        if (level + 1) % 2:
            core.set_local_option('DETCI', 'MPN_ORDER_SAVE', 2)
        else:
            core.set_local_option('DETCI', 'MPN_ORDER_SAVE', 1)
    elif mtdlvl_mobj and mtdlvl_mobj.group("method") == "mp":
        core.set_local_option('DETCI', 'WFN', 'DETCI')
        core.set_local_option('DETCI', 'MPN', 'TRUE')
        level = int(mtdlvl_mobj.group("level"))
        maxnvect = int((level + 1) / 2) + (level + 1) % 2
        core.set_local_option('DETCI', 'MAX_NUM_VECS', maxnvect)
        if (level + 1) % 2:
            core.set_local_option('DETCI', 'MPN_ORDER_SAVE', 2)
        else:
            core.set_local_option('DETCI', 'MPN_ORDER_SAVE', 1)
    elif name == 'ccsd':
        # untested
        core.set_local_option('DETCI', 'WFN', 'DETCI')
        core.set_local_option('DETCI', 'CC', 'TRUE')
        core.set_local_option('DETCI', 'CC_EX_LEVEL', 2)
    elif name == 'fci':
        core.set_local_option('DETCI', 'WFN', 'DETCI')
        core.set_local_option('DETCI', 'FCI', 'TRUE')
    elif name == 'cisd':
        core.set_local_option('DETCI', 'WFN', 'DETCI')
        core.set_local_option('DETCI', 'EX_LEVEL', 2)
    elif name == 'cisdt':
        core.set_local_option('DETCI', 'WFN', 'DETCI')
        core.set_local_option('DETCI', 'EX_LEVEL', 3)
    elif name == 'cisdtq':
        core.set_local_option('DETCI', 'WFN', 'DETCI')
        core.set_local_option('DETCI', 'EX_LEVEL', 4)
    elif mtdlvl_mobj and mtdlvl_mobj.group("method") == "ci":
        core.set_local_option('DETCI', 'WFN', 'DETCI')
        level = int(mtdlvl_mobj.group("level"))
        core.set_local_option('DETCI', 'EX_LEVEL', level)
    elif name == 'detci':
        pass

    # Bypass the scf call if a reference wavefunction is given
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, **kwargs)  # C1 certified

    # Ensure IWL files have been written
    proc_util.check_iwl_file_from_scf_type(core.get_global_option('SCF_TYPE'), ref_wfn)

    ciwfn = core.detci(ref_wfn)

    # Shove variables into global space
    for k, v in ciwfn.variables().items():
        core.set_variable(k, v)

    print_nos = False
    if core.get_option("DETCI", "NAT_ORBS"):
        ciwfn.ci_nat_orbs()
        print_nos = True

    proc_util.print_ci_results(ciwfn, name.upper(), ciwfn.variable("HF TOTAL ENERGY"), ciwfn.variable("CURRENT ENERGY"), print_nos)

    core.print_out("\t\t \"A good bug is a dead bug\" \n\n")
    core.print_out("\t\t\t - Starship Troopers\n\n")
    core.print_out("\t\t \"I didn't write FORTRAN.  That's the problem.\"\n\n")
    core.print_out("\t\t\t - Edward Valeev\n")

    if core.get_global_option("DIPMOM") and ("mp" not in name.lower()):
        # We always would like to print a little dipole information
        oeprop = core.OEProp(ciwfn)
        oeprop.set_title(name.upper())
        oeprop.add("DIPOLE")
        oeprop.compute()
        ciwfn.oeprop = oeprop
        core.set_variable("CURRENT DIPOLE", core.variable(name.upper() + " DIPOLE"))

    ciwfn.cleanup_ci()
    ciwfn.cleanup_dpd()
    _clean_detci()

    for lvl in range(4, 11):
        if ciwfn.has_variable(f"MP{lvl} CORRELATION ENERGY") and ciwfn.has_variable(f"MP{lvl-1} CORRELATION ENERGY"):
            ciwfn.set_variable(f"MP{lvl} CORRECTION ENERGY", ciwfn.variable(f"MP{lvl} CORRELATION ENERGY") - ciwfn.variable(f"MP{lvl-1} CORRELATION ENERGY"))
            core.set_variable(f"MP{lvl} CORRECTION ENERGY", ciwfn.variable(f"MP{lvl} CORRELATION ENERGY") - ciwfn.variable(f"MP{lvl-1} CORRELATION ENERGY"))

    optstash.restore()
    return ciwfn


def run_dfmp2(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a density-fitted MP2 calculation.

    """
    optstash = p4util.OptionsState(
        ['DF_BASIS_MP2'],
        ['SCF_TYPE'])

    # Alter default algorithm
    if not core.has_global_option_changed('SCF_TYPE'):
        core.set_global_option('SCF_TYPE', 'DF')
        core.print_out("""    SCF Algorithm Type (re)set to DF.\n""")

    # Bypass the scf call if a reference wavefunction is given
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, **kwargs)  # C1 certified

    core.tstart()
    core.print_out('\n')
    p4util.banner('DFMP2')
    core.print_out('\n')

    if core.get_global_option('REFERENCE') == "ROHF":
        ref_wfn.semicanonicalize()

    aux_basis = core.BasisSet.build(ref_wfn.molecule(), "DF_BASIS_MP2",
                                    core.get_option("DFMP2", "DF_BASIS_MP2"),
                                    "RIFIT", core.get_global_option('BASIS'))
    ref_wfn.set_basisset("DF_BASIS_MP2", aux_basis)

    dfmp2_wfn = core.dfmp2(ref_wfn)
    dfmp2_wfn.compute_energy()

    if name == 'scs-mp2':
        dfmp2_wfn.set_variable('CURRENT ENERGY', dfmp2_wfn.variable('SCS-MP2 TOTAL ENERGY'))
        dfmp2_wfn.set_variable('CURRENT CORRELATION ENERGY', dfmp2_wfn.variable('SCS-MP2 CORRELATION ENERGY'))

    elif name == 'mp2':
        dfmp2_wfn.set_variable('CURRENT ENERGY', dfmp2_wfn.variable('MP2 TOTAL ENERGY'))
        dfmp2_wfn.set_variable('CURRENT CORRELATION ENERGY', dfmp2_wfn.variable('MP2 CORRELATION ENERGY'))

    # Shove variables into global space
    for k, v in dfmp2_wfn.variables().items():
        core.set_variable(k, v)

    optstash.restore()
    core.tstop()
    return dfmp2_wfn


def run_dfep2(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a density-fitted MP2 calculation.

    """
    core.tstart()
    optstash = p4util.OptionsState(
        ['DF_BASIS_MP2'],
        ['SCF_TYPE'])

    # Alter default algorithm
    if not core.has_global_option_changed('SCF_TYPE'):
        core.set_global_option('SCF_TYPE', 'DF')
        core.print_out("""    SCF Algorithm Type (re)set to DF.\n""")

    # Bypass the scf call if a reference wavefunction is given
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, **kwargs)  # C1 certified

    if core.get_global_option('REFERENCE') != "RHF":
        raise ValidationError("DF-EP2 is not available for %s references.",
                              core.get_global_option('REFERENCE'))


    # Build the wavefunction
    aux_basis = core.BasisSet.build(ref_wfn.molecule(), "DF_BASIS_EP2",
                                    core.get_option("DFEP2", "DF_BASIS_EP2"),
                                    "RIFIT", core.get_global_option('BASIS'))
    ref_wfn.set_basisset("DF_BASIS_EP2", aux_basis)

    dfep2_wfn = core.DFEP2Wavefunction(ref_wfn)

    # Figure out what were doing
    if core.has_option_changed('DFEP2', 'EP2_ORBITALS'):
        ep2_input = core.get_global_option("EP2_ORBITALS")

    else:
        n_ip = core.get_global_option("EP2_NUM_IP")
        n_ea = core.get_global_option("EP2_NUM_EA")

        eps = np.hstack(dfep2_wfn.epsilon_a().nph)
        irrep_map = np.hstack([np.ones_like(dfep2_wfn.epsilon_a().nph[x]) * x for x in range(dfep2_wfn.nirrep())])
        sort = np.argsort(eps)

        ip_map = sort[dfep2_wfn.nalpha() - n_ip:dfep2_wfn.nalpha()]
        ea_map = sort[dfep2_wfn.nalpha():dfep2_wfn.nalpha() + n_ea]

        ep2_input = [[] for x in range(dfep2_wfn.nirrep())]
        nalphapi = tuple(dfep2_wfn.nalphapi())

        # Add IP info
        ip_info = np.unique(irrep_map[ip_map], return_counts=True)
        for irrep, cnt in zip(*ip_info):
            irrep = int(irrep)
            ep2_input[irrep].extend(range(nalphapi[irrep] - cnt, nalphapi[irrep]))

        # Add EA info
        ea_info = np.unique(irrep_map[ea_map], return_counts=True)
        for irrep, cnt in zip(*ea_info):
            irrep = int(irrep)
            ep2_input[irrep].extend(range(nalphapi[irrep], nalphapi[irrep] + cnt))

    # Compute
    ret = dfep2_wfn.compute(ep2_input)

    # Resort it...
    ret_eps = []
    for h in range(dfep2_wfn.nirrep()):
        ep2_data = ret[h]
        inp_data = ep2_input[h]

        for i in range(len(ep2_data)):
            tmp = [h, ep2_data[i][0], ep2_data[i][1], dfep2_wfn.epsilon_a().get(h, inp_data[i]), inp_data[i]]
            ret_eps.append(tmp)

    ret_eps.sort(key=lambda x: x[3])

    h2ev = constants.hartree2ev
    irrep_labels = dfep2_wfn.molecule().irrep_labels()

    core.print_out("  ==> Results <==\n\n")
    core.print_out("   %8s  %12s %12s %8s\n" % ("Orbital", "Koopmans (eV)", "EP2 (eV)", "EP2 PS"))
    core.print_out("  ----------------------------------------------\n")
    for irrep, ep2, ep2_ps, kt, pos in ret_eps:
        label = str(pos + 1)  + irrep_labels[irrep]
        core.print_out("  %8s    % 12.3f  % 12.3f   % 6.3f\n" % (label, (kt * h2ev), (ep2 * h2ev), ep2_ps))
        core.set_variable("EP2 " + label.upper() + " ENERGY", ep2)
    core.print_out("  ----------------------------------------------\n\n")

    # Figure out the IP and EA
    sorted_vals = np.array([x[1] for x in ret_eps])
    ip_vals = sorted_vals[sorted_vals < 0]
    ea_vals = sorted_vals[sorted_vals > 0]

    ip_value = None
    ea_value = None
    if len(ip_vals):
        core.set_variable("EP2 IONIZATION POTENTIAL", ip_vals[-1])
        core.set_variable("CURRENT ENERGY", ip_vals[-1])
    if len(ea_vals):
        core.set_variable("EP2 ELECTRON AFFINITY", ea_vals[0])
        if core.variable("EP2 IONIZATION POTENTIAL") == 0.0:
            core.set_variable("CURRENT ENERGY", ea_vals[0])

    core.print_out("  EP2 has completed successfully!\n\n")

    core.tstop()
    return dfep2_wfn


def run_dlpnomp2(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a DLPNO-MP2 calculation.

    """
    optstash = p4util.OptionsState(
        ['DF_BASIS_MP2'],
        ['SCF_TYPE'])

    # Alter default algorithm
    if not core.has_global_option_changed('SCF_TYPE'):
        core.set_global_option('SCF_TYPE', 'DF')
        core.print_out("""    SCF Algorithm Type (re)set to DF.\n""")

    # DLPNO-MP2 is only DF
    if core.get_global_option('MP2_TYPE') != "DF":
        raise ValidationError("""  DLPNO-MP2 is only implemented with density fitting.\n"""
                              """  'mp2_type' must be set to 'DF'.\n""")

    # Bypass the scf call if a reference wavefunction is given
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, use_c1=True, **kwargs)  # C1 certified
    elif ref_wfn.molecule().schoenflies_symbol() != 'c1':
        raise ValidationError("""  DLPNO-MP2 does not make use of molecular symmetry: """
                              """reference wavefunction must be C1.\n""")

    if core.get_global_option('REFERENCE') != "RHF":
        raise ValidationError("DLPNO-MP2 is not available for %s references.",
                              core.get_global_option('REFERENCE'))

    core.tstart()
    core.print_out('\n')
    p4util.banner('DLPNO-MP2')
    core.print_out('\n')

    aux_basis = core.BasisSet.build(ref_wfn.molecule(), "DF_BASIS_MP2",
                                    core.get_option("DLPNO", "DF_BASIS_MP2"),
                                    "RIFIT", core.get_global_option('BASIS'))
    ref_wfn.set_basisset("DF_BASIS_MP2", aux_basis)

    dlpnomp2_wfn = core.dlpno(ref_wfn)
    dlpnomp2_wfn.compute_energy()

    if name == 'scs-dlpno-mp2':
        dlpnomp2_wfn.set_variable('CURRENT ENERGY', dlpnomp2_wfn.variable('SCS-MP2 TOTAL ENERGY'))
        dlpnomp2_wfn.set_variable('CURRENT CORRELATION ENERGY', dlpnomp2_wfn.variable('SCS-MP2 CORRELATION ENERGY'))

    elif name == 'dlpno-mp2':
        dlpnomp2_wfn.set_variable('CURRENT ENERGY', dlpnomp2_wfn.variable('MP2 TOTAL ENERGY'))
        dlpnomp2_wfn.set_variable('CURRENT CORRELATION ENERGY', dlpnomp2_wfn.variable('MP2 CORRELATION ENERGY'))

    # Shove variables into global space
    for k, v in dlpnomp2_wfn.variables().items():
        core.set_variable(k, v)

    optstash.restore()
    core.tstop()
    return dlpnomp2_wfn


def run_dmrgscf(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    an DMRG calculation.

    """
    optstash = p4util.OptionsState(
        ['SCF_TYPE'],
        ['DMRG', 'DMRG_CASPT2_CALC'])

    # Bypass the scf call if a reference wavefunction is given
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, **kwargs)

    # Ensure IWL files have been written
    proc_util.check_iwl_file_from_scf_type(core.get_global_option('SCF_TYPE'), ref_wfn)

    if 'CASPT2' in name.upper():
        core.set_local_option("DMRG", "DMRG_CASPT2_CALC", True)

    dmrg_wfn = core.dmrg(ref_wfn)
    optstash.restore()

    # Shove variables into global space
    for k, v in dmrg_wfn.variables().items():
        core.set_variable(k, v)

    return dmrg_wfn


def run_dmrgci(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    an DMRG calculation.

    """
    optstash = p4util.OptionsState(
        ['SCF_TYPE'],
        ['DMRG', 'DMRG_SCF_MAX_ITER'])

    # Bypass the scf call if a reference wavefunction is given
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, **kwargs)

    # Ensure IWL files have been written
    proc_util.check_iwl_file_from_scf_type(core.get_global_option('SCF_TYPE'), ref_wfn)

    core.set_local_option('DMRG', 'DMRG_SCF_MAX_ITER', 1)

    dmrg_wfn = core.dmrg(ref_wfn)
    optstash.restore()

    # Shove variables into global space
    for k, v in dmrg_wfn.variables().items():
        core.set_variable(k, v)

    return dmrg_wfn


def run_psimrcc(name, **kwargs):
    """Function encoding sequence of PSI module calls for a PSIMRCC computation
     using a reference from the MCSCF module

    """
    mcscf_wfn = run_mcscf(name, **kwargs)
    psimrcc_wfn = core.psimrcc(mcscf_wfn)

    # Shove variables into global space
    for k, v in psimrcc_wfn.variables().items():
        core.set_variable(k, v)

    return psimrcc_wfn


def run_psimrcc_scf(name, **kwargs):
    """Function encoding sequence of PSI module calls for a PSIMRCC computation
     using a reference from the SCF module

    """
    # Bypass the scf call if a reference wavefunction is given
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, **kwargs)

    # Ensure IWL files have been written
    proc_util.check_iwl_file_from_scf_type(core.get_global_option('SCF_TYPE'), ref_wfn)

    psimrcc_wfn = core.psimrcc(ref_wfn)

    # Shove variables into global space
    for k, v in psimrcc_wfn.variables().items():
        core.set_variable(k, v)

    return psimrcc_wfn


def run_sapt(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a SAPT calculation of any level.

    """
    optstash = p4util.OptionsState(['SCF_TYPE'])

    # Alter default algorithm
    if not core.has_global_option_changed('SCF_TYPE'):
        core.set_global_option('SCF_TYPE', 'DF')

    # Get the molecule of interest
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        sapt_dimer = kwargs.pop('molecule', core.get_active_molecule())
    else:
        core.print_out('Warning! SAPT argument "ref_wfn" is only able to use molecule information.')
        sapt_dimer = ref_wfn.molecule()

    sapt_basis = kwargs.pop('sapt_basis', 'dimer')

    sapt_dimer, monomerA, monomerB = proc_util.prepare_sapt_molecule(sapt_dimer, sapt_basis)

    # Need to ensure consistent orbital freezing
    # between monomer and dimer computations
    monomerA_basis = core.BasisSet.build(monomerA, "BASIS", core.get_global_option("BASIS"))
    monomerB_basis = core.BasisSet.build(monomerB, "BASIS", core.get_global_option("BASIS"))
    nfc_ab = monomerA_basis.n_frozen_core() + monomerB_basis.n_frozen_core()

    if (core.get_option('SCF', 'REFERENCE') != 'RHF') and (name.upper() != "SAPT0"):
        raise ValidationError('Only SAPT0 supports a reference different from \"reference rhf\".')

    do_delta_mp2 = True if name.endswith('dmp2') else False
    do_empirical_disp = True if '-d' in name.lower() else False

    if do_empirical_disp:
        ## Make sure we are turning SAPT0 dispersion off
        core.set_local_option('SAPT', 'SAPT0_E10', True)
        core.set_local_option('SAPT', 'SAPT0_E20IND', True)
        core.set_local_option('SAPT', 'SAPT0_E20Disp', False)

    ri = core.get_global_option('SCF_TYPE')
    df_ints_io = core.get_option('SCF', 'DF_INTS_IO')
    # inquire if above at all applies to dfmp2

    core.IO.set_default_namespace('dimer')
    core.print_out('\n')
    p4util.banner('Dimer HF')
    core.print_out('\n')

    # Compute dimer wavefunction

    if (sapt_basis == 'dimer') and (ri == 'DF'):
        core.set_global_option('DF_INTS_IO', 'SAVE')

    optstash2 = p4util.OptionsState(['NUM_FROZEN_DOCC'])
    core.set_global_option("NUM_FROZEN_DOCC", nfc_ab)
    core.timer_on("SAPT: Dimer SCF")
    dimer_wfn = scf_helper('RHF', molecule=sapt_dimer, **kwargs)
    core.timer_off("SAPT: Dimer SCF")

    if do_delta_mp2:
        select_mp2("mp2", ref_wfn=dimer_wfn, **kwargs)
        mp2_corl_interaction_e = core.variable('MP2 CORRELATION ENERGY')

    optstash2.restore()
    if (sapt_basis == 'dimer') and (ri == 'DF'):
        core.set_global_option('DF_INTS_IO', 'LOAD')

    # Compute Monomer A wavefunction
    if (sapt_basis == 'dimer') and (ri == 'DF'):
        core.IO.change_file_namespace(97, 'dimer', 'monomerA')

    core.IO.set_default_namespace('monomerA')
    core.print_out('\n')
    p4util.banner('Monomer A HF')
    core.print_out('\n')

    core.timer_on("SAPT: Monomer A SCF")
    monomerA_wfn = scf_helper('RHF', molecule=monomerA, **kwargs)
    core.timer_off("SAPT: Monomer A SCF")

    if do_delta_mp2:
        select_mp2("mp2", ref_wfn=monomerA_wfn, **kwargs)
        mp2_corl_interaction_e -= core.variable('MP2 CORRELATION ENERGY')

    # Compute Monomer B wavefunction
    if (sapt_basis == 'dimer') and (ri == 'DF'):
        core.IO.change_file_namespace(97, 'monomerA', 'monomerB')
    core.IO.set_default_namespace('monomerB')
    core.print_out('\n')
    p4util.banner('Monomer B HF')
    core.print_out('\n')

    core.timer_on("SAPT: Monomer B SCF")
    monomerB_wfn = scf_helper('RHF', molecule=monomerB, **kwargs)
    core.timer_off("SAPT: Monomer B SCF")

    # Delta MP2
    if do_delta_mp2:
        select_mp2("mp2", ref_wfn=monomerB_wfn, **kwargs)
        mp2_corl_interaction_e -= core.variable('MP2 CORRELATION ENERGY')
        core.set_variable("SAPT MP2 CORRELATION ENERGY", mp2_corl_interaction_e)  # P::e SAPT
    core.set_global_option('DF_INTS_IO', df_ints_io)

    if core.get_option('SCF', 'REFERENCE') == 'RHF':
        core.IO.change_file_namespace(psif.PSIF_SAPT_MONOMERA, 'monomerA', 'dimer')
        core.IO.change_file_namespace(psif.PSIF_SAPT_MONOMERB, 'monomerB', 'dimer')

    core.IO.set_default_namespace('dimer')
    core.set_local_option('SAPT', 'E_CONVERGENCE', 10e-10)
    core.set_local_option('SAPT', 'D_CONVERGENCE', 10e-10)
    if name in ['sapt0', 'ssapt0']:
        core.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT0')
    elif name == 'sapt2':
        core.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT2')
    elif name in ['sapt2+', 'sapt2+dmp2']:
        core.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT2+')
        core.set_local_option('SAPT', 'DO_CCD_DISP', False)
    elif name in ['sapt2+(3)', 'sapt2+(3)dmp2']:
        core.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT2+3')
        core.set_local_option('SAPT', 'DO_THIRD_ORDER', False)
        core.set_local_option('SAPT', 'DO_CCD_DISP', False)
    elif name in ['sapt2+3', 'sapt2+3dmp2']:
        core.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT2+3')
        core.set_local_option('SAPT', 'DO_THIRD_ORDER', True)
        core.set_local_option('SAPT', 'DO_CCD_DISP', False)
    elif name in ['sapt2+(ccd)', 'sapt2+(ccd)dmp2']:
        core.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT2+')
        core.set_local_option('SAPT', 'DO_CCD_DISP', True)
    elif name in ['sapt2+(3)(ccd)', 'sapt2+(3)(ccd)dmp2']:
        core.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT2+3')
        core.set_local_option('SAPT', 'DO_THIRD_ORDER', False)
        core.set_local_option('SAPT', 'DO_CCD_DISP', True)
    elif name in ['sapt2+3(ccd)', 'sapt2+3(ccd)dmp2']:
        core.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT2+3')
        core.set_local_option('SAPT', 'DO_THIRD_ORDER', True)
        core.set_local_option('SAPT', 'DO_CCD_DISP', True)

    # Make sure we are not going to run CPHF on ROHF, since its MO Hessian
    # is not SPD
    if core.get_option('SCF', 'REFERENCE') == 'ROHF':
        core.set_local_option('SAPT', 'COUPLED_INDUCTION', False)
        core.print_out('  Coupled induction not available for ROHF.\n')
        core.print_out('  Proceeding with uncoupled induction only.\n')

    core.print_out("  Constructing Basis Sets for SAPT...\n\n")
    aux_basis = core.BasisSet.build(dimer_wfn.molecule(), "DF_BASIS_SAPT", core.get_global_option("DF_BASIS_SAPT"),
                                    "RIFIT", core.get_global_option("BASIS"))
    dimer_wfn.set_basisset("DF_BASIS_SAPT", aux_basis)

    aux_basis = core.BasisSet.build(dimer_wfn.molecule(), "DF_BASIS_ELST", core.get_global_option("DF_BASIS_ELST"),
                                    "JKFIT", core.get_global_option("BASIS"))
    dimer_wfn.set_basisset("DF_BASIS_ELST", aux_basis)

    core.print_out('\n')
    p4util.banner(name.upper())
    core.print_out('\n')
    e_sapt = core.sapt(dimer_wfn, monomerA_wfn, monomerB_wfn)
    dimer_wfn.set_module("sapt")

    from psi4.driver.qcdb.psivardefs import sapt_psivars
    p4util.expand_psivars(sapt_psivars())
    optstash.restore()

    # Get the SAPT name right if doing empirical dispersion
    if do_empirical_disp:
        sapt_name = "sapt0"
    else:
        sapt_name = name

    # Make sure we got induction, otherwise replace it with uncoupled induction
    which_ind = 'IND'
    target_ind = 'IND'
    if not core.has_variable(' '.join((sapt_name.upper(), which_ind, 'ENERGY'))):
        which_ind = 'IND,U'

    for term in ['ELST', 'EXCH', 'DISP', 'TOTAL']:
        core.set_variable(' '.join(['SAPT', term, 'ENERGY']),
                          core.variable(' '.join([sapt_name.upper(), term, 'ENERGY'])))
    # Special induction case
    core.set_variable(' '.join(['SAPT', target_ind, 'ENERGY']),
                      core.variable(' '.join([sapt_name.upper(), which_ind, 'ENERGY'])))
    core.set_variable('CURRENT ENERGY', core.variable('SAPT TOTAL ENERGY'))

    # Empirical dispersion
    if do_empirical_disp:
        proc_util.sapt_empirical_dispersion(name, dimer_wfn)

    return dimer_wfn


def run_sapt_ct(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a charge-transfer SAPT calcuation of any level.

    """
    optstash = p4util.OptionsState(
        ['SCF_TYPE'])

    if 'ref_wfn' in kwargs:
        core.print_out('\nWarning! Argument ref_wfn is not valid for sapt computations\n')

    # Alter default algorithm
    if not core.has_global_option_changed('SCF_TYPE'):
        core.set_global_option('SCF_TYPE', 'DF')

    # Get the molecule of interest
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        sapt_dimer = kwargs.pop('molecule', core.get_active_molecule())
    else:
        core.print_out('Warning! SAPT argument "ref_wfn" is only able to use molecule information.')
        sapt_dimer = ref_wfn.molecule()

    sapt_dimer, monomerA, monomerB = proc_util.prepare_sapt_molecule(sapt_dimer, "dimer")
    monomerAm = sapt_dimer.extract_subsets(1)
    monomerAm.set_name('monomerAm')
    monomerBm = sapt_dimer.extract_subsets(2)
    monomerBm.set_name('monomerBm')

    if core.get_option('SCF', 'REFERENCE') != 'RHF':
        raise ValidationError('SAPT requires requires \"reference rhf\".')

    ri = core.get_global_option('SCF_TYPE')
    df_ints_io = core.get_option('SCF', 'DF_INTS_IO')
    # inquire if above at all applies to dfmp2

    core.IO.set_default_namespace('dimer')
    core.print_out('\n')
    p4util.banner('Dimer HF')
    core.print_out('\n')
    core.set_global_option('DF_INTS_IO', 'SAVE')
    dimer_wfn = scf_helper('RHF', molecule=sapt_dimer, **kwargs)
    core.set_global_option('DF_INTS_IO', 'LOAD')

    if (ri == 'DF'):
        core.IO.change_file_namespace(97, 'dimer', 'monomerA')
    core.IO.set_default_namespace('monomerA')
    core.print_out('\n')
    p4util.banner('Monomer A HF (Dimer Basis)')
    core.print_out('\n')
    monomerA_wfn = scf_helper('RHF', molecule=monomerA, **kwargs)

    if (ri == 'DF'):
        core.IO.change_file_namespace(97, 'monomerA', 'monomerB')
    core.IO.set_default_namespace('monomerB')
    core.print_out('\n')
    p4util.banner('Monomer B HF (Dimer Basis)')
    core.print_out('\n')
    monomerB_wfn = scf_helper('RHF', molecule=monomerB, **kwargs)
    core.set_global_option('DF_INTS_IO', df_ints_io)

    core.IO.set_default_namespace('monomerAm')
    core.print_out('\n')
    p4util.banner('Monomer A HF (Monomer Basis)')
    core.print_out('\n')
    monomerAm_wfn = scf_helper('RHF', molecule=monomerAm, **kwargs)

    core.IO.set_default_namespace('monomerBm')
    core.print_out('\n')
    p4util.banner('Monomer B HF (Monomer Basis)')
    core.print_out('\n')
    monomerBm_wfn = scf_helper('RHF', molecule=monomerBm, **kwargs)

    core.IO.set_default_namespace('dimer')
    core.set_local_option('SAPT', 'E_CONVERGENCE', 10e-10)
    core.set_local_option('SAPT', 'D_CONVERGENCE', 10e-10)
    if name == 'sapt0-ct':
        core.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT0')
    elif name == 'sapt2-ct':
        core.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT2')
    elif name == 'sapt2+-ct':
        core.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT2+')
    elif name == 'sapt2+(3)-ct':
        core.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT2+3')
        core.set_local_option('SAPT', 'DO_THIRD_ORDER', False)
    elif name == 'sapt2+3-ct':
        core.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT2+3')
        core.set_local_option('SAPT', 'DO_THIRD_ORDER', True)
    elif name == 'sapt2+(ccd)-ct':
        core.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT2+')
        core.set_local_option('SAPT', 'DO_CCD_DISP', True)
    elif name == 'sapt2+(3)(ccd)-ct':
        core.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT2+3')
        core.set_local_option('SAPT', 'DO_THIRD_ORDER', False)
        core.set_local_option('SAPT', 'DO_CCD_DISP', True)
    elif name == 'sapt2+3(ccd)-ct':
        core.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT2+3')
        core.set_local_option('SAPT', 'DO_THIRD_ORDER', True)
        core.set_local_option('SAPT', 'DO_CCD_DISP', True)

    core.print_out('\n')

    aux_basis = core.BasisSet.build(dimer_wfn.molecule(), "DF_BASIS_SAPT",
                                        core.get_global_option("DF_BASIS_SAPT"),
                                        "RIFIT", core.get_global_option("BASIS"))
    dimer_wfn.set_basisset("DF_BASIS_SAPT", aux_basis)
    if core.get_global_option("DF_BASIS_ELST") == "":
        dimer_wfn.set_basisset("DF_BASIS_ELST", aux_basis)
    else:
        aux_basis = core.BasisSet.build(dimer_wfn.molecule(), "DF_BASIS_ELST",
                                            core.get_global_option("DF_BASIS_ELST"),
                                            "RIFIT", core.get_global_option("BASIS"))
        dimer_wfn.set_basisset("DF_BASIS_ELST", aux_basis)


    core.print_out('\n')
    p4util.banner('SAPT Charge Transfer')
    core.print_out('\n')

    core.print_out('\n')
    p4util.banner('Dimer Basis SAPT')
    core.print_out('\n')
    core.IO.change_file_namespace(psif.PSIF_SAPT_MONOMERA, 'monomerA', 'dimer')
    core.IO.change_file_namespace(psif.PSIF_SAPT_MONOMERB, 'monomerB', 'dimer')
    e_sapt = core.sapt(dimer_wfn, monomerA_wfn, monomerB_wfn)
    CTd = core.variable('SAPT CT ENERGY')
    dimer_wfn.set_module("sapt")

    core.print_out('\n')
    p4util.banner('Monomer Basis SAPT')
    core.print_out('\n')
    core.IO.change_file_namespace(psif.PSIF_SAPT_MONOMERA, 'monomerAm', 'dimer')
    core.IO.change_file_namespace(psif.PSIF_SAPT_MONOMERB, 'monomerBm', 'dimer')
    e_sapt = core.sapt(dimer_wfn, monomerAm_wfn, monomerBm_wfn)
    CTm = core.variable('SAPT CT ENERGY')
    CT = CTd - CTm

    units = (1000.0, constants.hartree2kcalmol, constants.hartree2kJmol)
    core.print_out('\n\n')
    core.print_out('    SAPT Charge Transfer Analysis\n')
    core.print_out('  ------------------------------------------------------------------------------------------------\n')
    core.print_out('    SAPT Induction (Dimer Basis)  %12.4lf [mEh] %12.4lf [kcal/mol] %12.4lf [kJ/mol]\n' %
        tuple(CTd * u for u in units))
    core.print_out('    SAPT Induction (Monomer Basis)%12.4lf [mEh] %12.4lf [kcal/mol] %12.4lf [kJ/mol]\n' %
        tuple(CTm * u for u in units))
    core.print_out('    SAPT Charge Transfer          %12.4lf [mEh] %12.4lf [kcal/mol] %12.4lf [kJ/mol]\n\n' %
        tuple(CT * u for u in units))
    core.set_variable("SAPT CT ENERGY", CT)  # P::e SAPT

    optstash.restore()
    return dimer_wfn

def run_fisapt(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    an F/ISAPT0 computation

    """
    optstash = p4util.OptionsState(['SCF_TYPE'])

    # Alter default algorithm
    if not core.has_global_option_changed('SCF_TYPE'):
        core.set_global_option('SCF_TYPE', 'DF')

    # Get the molecule of interest
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        sapt_dimer = kwargs.pop('molecule', core.get_active_molecule())
    else:
        core.print_out('Warning! FISAPT argument "ref_wfn" is only able to use molecule information.')
        sapt_dimer = ref_wfn.molecule()
    sapt_dimer.update_geometry()  # make sure since mol from wfn, kwarg, or P::e

    # Shifting to C1 so we need to copy the active molecule
    if sapt_dimer.schoenflies_symbol() != 'c1':
        core.print_out('  FISAPT does not make use of molecular symmetry, further calculations in C1 point group.\n')
        sapt_dimer = sapt_dimer.clone()
        sapt_dimer.reset_point_group('c1')
        sapt_dimer.fix_orientation(True)
        sapt_dimer.fix_com(True)
        sapt_dimer.update_geometry()

    if core.get_option('SCF', 'REFERENCE') != 'RHF':
        raise ValidationError('FISAPT requires requires \"reference rhf\".')

    if ref_wfn is None:
        core.timer_on("FISAPT: Dimer SCF")
        ref_wfn = scf_helper('RHF', molecule=sapt_dimer, **kwargs)
        core.timer_off("FISAPT: Dimer SCF")

    core.print_out("  Constructing Basis Sets for FISAPT...\n\n")
    scf_aux_basis = core.BasisSet.build(ref_wfn.molecule(),
                                        "DF_BASIS_SCF",
                                        core.get_option("SCF", "DF_BASIS_SCF"),
                                        "JKFIT",
                                        core.get_global_option('BASIS'),
                                        puream=ref_wfn.basisset().has_puream())
    ref_wfn.set_basisset("DF_BASIS_SCF", scf_aux_basis)

    sapt_basis = core.BasisSet.build(ref_wfn.molecule(), "DF_BASIS_SAPT", core.get_global_option("DF_BASIS_SAPT"),
                                     "RIFIT", core.get_global_option("BASIS"),
                                     ref_wfn.basisset().has_puream())
    ref_wfn.set_basisset("DF_BASIS_SAPT", sapt_basis)

    minao = core.BasisSet.build(ref_wfn.molecule(), "BASIS", core.get_global_option("MINAO_BASIS"))
    ref_wfn.set_basisset("MINAO", minao)

    # Turn of dispersion for -d
    if "-d" in name.lower():
        core.set_local_option("FISAPT", "FISAPT_DO_FSAPT_DISP", False)

    fisapt_wfn = core.FISAPT(ref_wfn)
    from .sapt import fisapt_proc
    fisapt_wfn.compute_energy(external_potentials=kwargs.get("external_potentials", None))

    # Compute -D dispersion
    if "-d" in name.lower():
        proc_util.sapt_empirical_dispersion(name, ref_wfn)

    optstash.restore()
    return ref_wfn

def run_mrcc(name, **kwargs):
    """Function that prepares environment and input files
    for a calculation calling Kallay's MRCC code.

    """
    from .proc_data import mrcc_methods

    # level is a dictionary of settings to be passed to core.mrcc
    try:
        level = mrcc_methods[name.lower()]
    except KeyError:
        if name.lower() == "a-ccsd(t)":
            level = mrcc_methods["ccsd(t)_l"]
        else:
            raise ValidationError(f"""MRCC method '{name}' invalid.""")

    # Check to see if we really need to run the SCF code.
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, **kwargs)
    vscf = ref_wfn.variable('SCF TOTAL ENERGY')

    # Fullname is the string we need to search for in iface
    fullname = level['fullname']

    # User can provide 'keep' to the method.
    # When provided, do not delete the MRCC scratch directory.
    keep = False
    if 'keep' in kwargs:
        keep = kwargs['keep']

    # Save current directory location
    current_directory = os.getcwd()

    # Find environment by merging PSIPATH and PATH environment variables
    lenv = {
        'PATH': ':'.join([os.path.abspath(x) for x in os.environ.get('PSIPATH', '').split(':') if x != '']) +
                ':' + os.environ.get('PATH'),
        'LD_LIBRARY_PATH': os.environ.get('LD_LIBRARY_PATH')
        }
    #   Filter out None values as subprocess will fault on them
    lenv = {k: v for k, v in lenv.items() if v is not None}

    # Need to move to the scratch directory, perferrably into a separate directory in that location
    psi_io = core.IOManager.shared_object()
    os.chdir(psi_io.get_default_path())

    # Make new directory specifically for mrcc
    mrcc_tmpdir = 'mrcc_' + str(os.getpid())
    if 'path' in kwargs:
        mrcc_tmpdir = kwargs['path']

    # Check to see if directory already exists, if not, create.
    if os.path.exists(mrcc_tmpdir) is False:
        os.mkdir(mrcc_tmpdir)

    # Move into the new directory
    os.chdir(mrcc_tmpdir)

    # Generate integrals and input file (dumps files to the current directory)
    core.mrcc_generate_input(ref_wfn, level)
    ref_wfn.set_module("mrcc")

    # Load the fort.56 file
    # and dump a copy into the outfile
    core.print_out('\n===== Begin fort.56 input for MRCC ======\n')
    core.print_out(open('fort.56', 'r').read())
    core.print_out('===== End   fort.56 input for MRCC ======\n')

    # Modify the environment:
    #    PGI Fortan prints warning to screen if STOP is used
    lenv['NO_STOP_MESSAGE'] = '1'

    # Obtain the number of threads MRCC should use
    lenv['OMP_NUM_THREADS'] = str(core.get_num_threads())

    # If the user provided MRCC_OMP_NUM_THREADS set the environ to it
    if core.has_option_changed('MRCC', 'MRCC_OMP_NUM_THREADS'):
        lenv['OMP_NUM_THREADS'] = str(core.get_option('MRCC', 'MRCC_OMP_NUM_THREADS'))

    # Call dmrcc, directing all screen output to the output file
    external_exe = 'dmrcc'
    try:
        retcode = subprocess.Popen([external_exe], bufsize=0, stdout=subprocess.PIPE, env=lenv)
    except OSError as e:
        sys.stderr.write('Program %s not found in path or execution failed: %s\n' % (external_exe, e.strerror))
        core.print_out('Program %s not found in path or execution failed: %s\n' % (external_exe, e.strerror))
        message = ("Program %s not found in path or execution failed: %s\n" % (external_exe, e.strerror))
        raise ValidationError(message)

    c4out = ''
    while True:
        data = retcode.stdout.readline()
        if not data:
            break
        core.print_out(data.decode('utf-8'))
        c4out += data.decode('utf-8')

    # Scan iface file and grab the file energy.
    ene = 0.0
    for line in open('iface'):
        fields = line.split()
        m = fields[1]
        try:
            ene = float(fields[5])
            if m == "MP(2)":
                m = "MP2"
            elif m == "CCSD[T]":
                m = "CCSD+T(CCSD)"
            elif m == "CCSD(T)_L":
                m = "A-CCSD(T)"
            ref_wfn.set_variable(m + ' TOTAL ENERGY', ene)
            ref_wfn.set_variable(m + ' CORRELATION ENERGY', ene - vscf)
        except ValueError:
            continue

    # The last 'ene' in iface is the one the user requested.
    ref_wfn.set_variable('CURRENT ENERGY', ene)
    ref_wfn.set_variable('CURRENT CORRELATION ENERGY', ene - vscf)

    for subm in ["MP2", "CCSD"]:
        if ref_wfn.has_variable(f"{subm} TOTAL ENERGY") and core.get_option("SCF", "REFERENCE") in ["RHF", "UHF"]:
            ref_wfn.set_variable(f"{subm} SINGLES ENERGY", 0.0)
            ref_wfn.set_variable(f"{subm} DOUBLES ENERGY", ref_wfn.variable(f"{subm} CORRELATION ENERGY"))

    if ref_wfn.has_variable("CCSD(T) CORRELATION ENERGY") and ref_wfn.has_variable("CCSD CORRELATION ENERGY"):
        ref_wfn.set_variable("(T) CORRECTION ENERGY", ref_wfn.variable("CCSD(T) CORRELATION ENERGY") - ref_wfn.variable("CCSD CORRELATION ENERGY"))

    if ref_wfn.has_variable("A-CCSD(T) CORRELATION ENERGY") and ref_wfn.has_variable("CCSD CORRELATION ENERGY"):
        ref_wfn.set_variable("A-(T) CORRECTION ENERGY", ref_wfn.variable("A-CCSD(T) CORRELATION ENERGY") - ref_wfn.variable("CCSD CORRELATION ENERGY"))

    # Load the iface file
    iface = open('iface', 'r')
    iface_contents = iface.read()

    # Delete mrcc tempdir
    os.chdir('..')
    try:
        # Delete unless we're told not to
        if (keep is False and not('path' in kwargs)):
            shutil.rmtree(mrcc_tmpdir)
    except OSError as e:
        print('Unable to remove MRCC temporary directory %s' % e, file=sys.stderr)
        exit(1)

    # Return to submission directory
    os.chdir(current_directory)

    # If we're told to keep the files or the user provided a path, do nothing.
    if keep or ('path' in kwargs):
        core.print_out('\nMRCC scratch files have been kept.\n')
        core.print_out('They can be found in ' + mrcc_tmpdir)

    # Dump iface contents to output
    core.print_out('\n')
    p4util.banner('Full results from MRCC')
    core.print_out('\n')
    core.print_out(iface_contents)

    # Shove variables into global space
    for k, v in ref_wfn.variables().items():
        core.set_variable(k, v)

    core.print_variables()
    return ref_wfn


def run_fnodfcc(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a [FNO-](DF|CD)-CCSD[(T)] computation.

    >>> set cc_type df
    >>> energy('fno-ccsd(t)')

    """
    kwargs = p4util.kwargs_lower(kwargs)
    dtl = docs_table_link("dummy", "fnocc")

    # stash user options
    optstash = p4util.OptionsState(
        ['FNOCC', 'COMPUTE_TRIPLES'],
        ['FNOCC', 'DFCC'],
        ['FNOCC', 'NAT_ORBS'],
        ['FNOCC', 'RUN_CEPA'],
        ['FNOCC', 'DF_BASIS_CC'],
        ['SCF', 'DF_BASIS_SCF'],
        ['SCF', 'DF_INTS_IO'])

    def set_cholesky_from(type_val):
        if type_val == 'CD':
            core.set_local_option('FNOCC', 'DF_BASIS_CC', 'CHOLESKY')
            # Alter default algorithm
            if not core.has_global_option_changed('SCF_TYPE'):
                optstash.add_option(['SCF_TYPE'])
                core.set_global_option('SCF_TYPE', 'CD')
                core.print_out("""    SCF Algorithm Type (re)set to CD.\n""")

        elif type_val in ['DISK_DF', 'DF']:
            if core.get_option('FNOCC', 'DF_BASIS_CC') == 'CHOLESKY':
                core.set_local_option('FNOCC', 'DF_BASIS_CC', '')

            proc_util.check_disk_df(name.upper(), optstash)
        else:
            raise ValidationError(f"Invalid type {type_val} for FNOCC energy through `run_fnodfcc`. See Capabilities Table at {dtl}")

    director = {
                 # Note "nat_orbs" not set defensively False for non-"fno" calls
           "ccsd":     {                   "dfcc": True,  "run_cepa": False, "compute_triples": False,},
       "fno-ccsd":     {"nat_orbs": True,  "dfcc": True,  "run_cepa": False, "compute_triples": False,},

           "ccsd(t)":  {                   "dfcc": True,  "run_cepa": False, "compute_triples": True, },
       "fno-ccsd(t)":  {"nat_orbs": True,  "dfcc": True,  "run_cepa": False, "compute_triples": True, },
    }

    if name not in director:
        raise ValidationError(f"Invalid method {name} for FNOCC energy")

    # throw exception for open-shells
    if (ref := core.get_option("SCF", "REFERENCE")) != "RHF":
        raise ValidationError(f"Invalid reference type {ref} != RHF for FNOCC energy. See Capabilities Table at {dtl}.")

    # throw exception for CONV (approximately). after defaulting logic, throw exception for SCF_TYPE CONV (approximately)
    set_cholesky_from(method_algorithm_type(name).now)
    if (scf_type := core.get_global_option("SCF_TYPE")) not in ["CD", "DISK_DF"]:
        raise ValidationError(f"Invalid {scf_type=} for FNOCC energy through `run_fnodfcc`. See Capabilities Table at {dtl}")

    for k, v in director[name].items():
        core.set_local_option("FNOCC", k.upper(), v)

    # save DF or CD ints generated by SCF for use in CC
    core.set_local_option('SCF', 'DF_INTS_IO', 'SAVE')

    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, use_c1=True, **kwargs)  # C1 certified
    else:
        if ref_wfn.molecule().schoenflies_symbol() != 'c1':
            raise ValidationError("""  FNOCC does not make use of molecular symmetry: """
                                  """reference wavefunction must be C1.\n""")

    core.print_out("  Constructing Basis Sets for FNOCC...\n\n")
    scf_aux_basis = core.BasisSet.build(ref_wfn.molecule(), "DF_BASIS_SCF",
                                        core.get_option("SCF", "DF_BASIS_SCF"),
                                        "JKFIT", core.get_global_option('BASIS'),
                                        puream=ref_wfn.basisset().has_puream())
    ref_wfn.set_basisset("DF_BASIS_SCF", scf_aux_basis)

    aux_basis = core.BasisSet.build(ref_wfn.molecule(), "DF_BASIS_CC",
                                        core.get_global_option("DF_BASIS_CC"),
                                        "RIFIT", core.get_global_option("BASIS"))
    ref_wfn.set_basisset("DF_BASIS_CC", aux_basis)

    if core.get_global_option("RELATIVISTIC") in ["X2C", "DKH"]:
        rel_bas = core.BasisSet.build(ref_wfn.molecule(), "BASIS_RELATIVISTIC",
                                          core.get_option("SCF", "BASIS_RELATIVISTIC"),
                                          "DECON", core.get_global_option('BASIS'),
                                          puream=ref_wfn.basisset().has_puream())
        ref_wfn.set_basisset('BASIS_RELATIVISTIC',rel_bas)

    fnocc_wfn = core.fnocc(ref_wfn)

    # Shove variables into global space
    for k, v in fnocc_wfn.variables().items():
        core.set_variable(k, v)

    optstash.restore()
    return fnocc_wfn


def run_fnocc(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a QCISD(T), CCSD(T), MP2.5, MP3, and MP4 computation.

    >>> energy('fno-ccsd(t)')

    """
    kwargs = p4util.kwargs_lower(kwargs)
    dtl = docs_table_link("dummy", "fnocc")

    # stash user options:
    optstash = p4util.OptionsState(
        ['TRANSQT2', 'WFN'],
        ['FNOCC', 'RUN_MP2'],
        ['FNOCC', 'RUN_MP3'],
        ['FNOCC', 'RUN_MP4'],
        ['FNOCC', 'RUN_CCSD'],
        ['FNOCC', 'COMPUTE_TRIPLES'],
        ['FNOCC', 'COMPUTE_MP4_TRIPLES'],
        ['FNOCC', 'DFCC'],
        ['FNOCC', 'RUN_CEPA'],
        ['FNOCC', 'USE_DF_INTS'],
        ['FNOCC', 'NAT_ORBS'])

    # AED reinforces default
    core.set_local_option('FNOCC', 'USE_DF_INTS', False)

    if name in ["mp3", "fno-mp3"] and not core.has_global_option_changed("MP_TYPE"):
        core.print_out(f"    Information: {name.upper()} default algorithm changed to DF in August 2020. Use `set mp_type conv` for previous behavior.\n")

    director = {
                 # Note "nat_orbs" not set defensively False for non-"fno" calls
           "mp2":      {                  "dfcc": False, "run_cepa": False,                                              "run_mp2": True,                                                                   },

           "mp3":      {                  "dfcc": False, "run_cepa": False,                                                                "run_mp3": True,                                                 },
       "fno-mp3":      {"nat_orbs": True, "dfcc": False, "run_cepa": False,                                                                "run_mp3": True,                                                 },

           "mp4(sdq)": {                  "dfcc": False, "run_cepa": False,                    "compute_triples": False,                                     "run_mp4": True,  "compute_mp4_triples": False,},
       "fno-mp4(sdq)": {"nat_orbs": True, "dfcc": False, "run_cepa": False,                    "compute_triples": False,                                     "run_mp4": True,  "compute_mp4_triples": False,},

           "mp4":      {                  "dfcc": False, "run_cepa": False,                    "compute_triples": True,                                      "run_mp4": True,  "compute_mp4_triples": True, },
       "fno-mp4":      {"nat_orbs": True, "dfcc": False, "run_cepa": False,                    "compute_triples": True,                                      "run_mp4": True,  "compute_mp4_triples": True, },

          "qcisd":     {                  "dfcc": False, "run_cepa": False, "run_ccsd": False, "compute_triples": False,                                                                                    },
      "fno-qcisd":     {"nat_orbs": True, "dfcc": False, "run_cepa": False, "run_ccsd": False, "compute_triples": False,                                                                                    },

          "qcisd(t)":  {                  "dfcc": False, "run_cepa": False, "run_ccsd": False, "compute_triples": True,                                                                                     },
      "fno-qcisd(t)":  {"nat_orbs": True, "dfcc": False, "run_cepa": False, "run_ccsd": False, "compute_triples": True,                                                                                     },

           "ccsd":     {                  "dfcc": False, "run_cepa": False, "run_ccsd": True,  "compute_triples": False,                                                                                    },
       "fno-ccsd":     {"nat_orbs": True, "dfcc": False, "run_cepa": False, "run_ccsd": True,  "compute_triples": False,                                                                                    },

           "ccsd(t)":  {                  "dfcc": False, "run_cepa": False, "run_ccsd": True,  "compute_triples": True,                                                                                     },
       "fno-ccsd(t)":  {"nat_orbs": True, "dfcc": False, "run_cepa": False, "run_ccsd": True,  "compute_triples": True,                                                                                     },
    }

    if name not in director:
        raise ValidationError(f"Invalid method {name} for FNOCC energy")

    # throw exception for open-shells
    if (ref := core.get_option("SCF", "REFERENCE")) != "RHF":
        raise ValidationError(f"Invalid reference type {ref} != RHF for FNOCC energy. See Capabilities Table at {dtl}")

    # throw exception for DF/CD. most of these pre-trapped by select_* functions but some escape, incl. mp4(sdq) and qcisd variants
    if (corl_type := method_algorithm_type(name).now) != "CONV":
        raise ValidationError(f"Invalid type {corl_type} for FNOCC energy through `run_fnocc`. See Capabilities Table at {dtl}")

    for k, v in director[name].items():
        core.set_local_option("FNOCC", k.upper(), v)

    # Bypass the scf call if a reference wavefunction is given
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, **kwargs)  # C1 certified

    if not core.get_option('FNOCC', 'USE_DF_INTS'):
        # Ensure IWL files have been written
        proc_util.check_iwl_file_from_scf_type(core.get_global_option('SCF_TYPE'), ref_wfn)
    else:
        core.print_out("  Constructing Basis Sets for FNOCC...\n\n")
        scf_aux_basis = core.BasisSet.build(ref_wfn.molecule(), "DF_BASIS_SCF",
                                            core.get_option("SCF", "DF_BASIS_SCF"),
                                            "JKFIT", core.get_global_option('BASIS'),
                                            puream=ref_wfn.basisset().has_puream())
        ref_wfn.set_basisset("DF_BASIS_SCF", scf_aux_basis)

    if core.get_global_option("RELATIVISTIC") in ["X2C", "DKH"]:
        rel_bas = core.BasisSet.build(ref_wfn.molecule(), "BASIS_RELATIVISTIC",
                                          core.get_option("SCF", "BASIS_RELATIVISTIC"),
                                          "DECON", core.get_global_option('BASIS'),
                                          puream=ref_wfn.basisset().has_puream())
        ref_wfn.set_basisset('BASIS_RELATIVISTIC',rel_bas)

    fnocc_wfn = core.fnocc(ref_wfn)

    # set current correlation energy and total energy. only need to treat mpn here.
    if (lbl := name.replace("fno-", "")) in ["mp3", "mp4(sdq)", "mp4"]:
        fnocc_wfn.set_variable("CURRENT ENERGY", fnocc_wfn.variable(f"{lbl.upper()} TOTAL ENERGY"))
        fnocc_wfn.set_variable("CURRENT CORRELATION ENERGY", fnocc_wfn.variable(f"{lbl.upper()} CORRELATION ENERGY"))
        if lbl == "mp4":
            fnocc_wfn.set_variable("MP4 CORRECTION ENERGY", fnocc_wfn.variable("MP4 CORRELATION ENERGY") - fnocc_wfn.variable("MP3 CORRELATION ENERGY"))

    # Shove variables into global space
    for k, v in fnocc_wfn.variables().items():
        core.set_variable(k, v)

    optstash.restore()
    return fnocc_wfn


def run_cepa(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a cepa-like calculation.

    >>> energy('cepa(1)')

    """
    kwargs = p4util.kwargs_lower(kwargs)
    dtl = docs_table_link("dummy", "fnocc")

    # save user options
    optstash = p4util.OptionsState(
        ['TRANSQT2', 'WFN'],
        ['FNOCC', 'NAT_ORBS'],
        ['FNOCC', 'DFCC'],
        ['FNOCC', 'RUN_CEPA'],
        ['FNOCC', 'USE_DF_INTS'],
        ['FNOCC', 'CEPA_LEVEL'],
        ['FNOCC', 'CEPA_NO_SINGLES'])

    # AED reinforces default
    core.set_local_option('FNOCC', 'USE_DF_INTS', False)

    director = {
                 # Note "nat_orbs" not set defensively False for non-"fno" calls
           "lccd":     {                   "dfcc": False, "run_cepa": True,  "cepa_level": "cepa(0)", "cepa_no_singles": True, },
       "fno-lccd":     {"nat_orbs": True,  "dfcc": False, "run_cepa": True,  "cepa_level": "cepa(0)", "cepa_no_singles": True, },

           "lccsd":    {                   "dfcc": False, "run_cepa": True,  "cepa_level": "cepa(0)", "cepa_no_singles": False,},
       "fno-lccsd":    {"nat_orbs": True,  "dfcc": False, "run_cepa": True,  "cepa_level": "cepa(0)", "cepa_no_singles": False,},
            "cepa(0)": {                   "dfcc": False, "run_cepa": True,  "cepa_level": "cepa(0)", "cepa_no_singles": False,},
        "fno-cepa(0)": {"nat_orbs": True,  "dfcc": False, "run_cepa": True,  "cepa_level": "cepa(0)", "cepa_no_singles": False,},

            "cepa(1)": {                   "dfcc": False, "run_cepa": True,  "cepa_level": "cepa(1)", "cepa_no_singles": False,},
        "fno-cepa(1)": {"nat_orbs": True,  "dfcc": False, "run_cepa": True,  "cepa_level": "cepa(1)", "cepa_no_singles": False,},

            "cepa(3)": {                   "dfcc": False, "run_cepa": True,  "cepa_level": "cepa(3)", "cepa_no_singles": False,},
        "fno-cepa(3)": {"nat_orbs": True,  "dfcc": False, "run_cepa": True,  "cepa_level": "cepa(3)", "cepa_no_singles": False,},

            "acpf":    {                   "dfcc": False, "run_cepa": True,  "cepa_level": "acpf",    "cepa_no_singles": False,},
        "fno-acpf":    {"nat_orbs": True,  "dfcc": False, "run_cepa": True,  "cepa_level": "acpf",    "cepa_no_singles": False,},

          "aqcc":      {                   "dfcc": False, "run_cepa": True,  "cepa_level": "aqcc",    "cepa_no_singles": False,},
      "fno-aqcc":      {"nat_orbs": True,  "dfcc": False, "run_cepa": True,  "cepa_level": "aqcc",    "cepa_no_singles": False,},

          "cisd":      {                   "dfcc": False, "run_cepa": True,  "cepa_level": "cisd",    "cepa_no_singles": False,},
      "fno-cisd":      {"nat_orbs": True,  "dfcc": False, "run_cepa": True,  "cepa_level": "cisd",    "cepa_no_singles": False,},
    }

    if name not in director:
        raise ValidationError(f"Invalid method {name} for FNOCC energy")

    # throw exception for open-shells
    if (ref := core.get_option("SCF", "REFERENCE")) != "RHF":
        raise ValidationError(f"Invalid reference type {ref} != RHF for FNOCC energy. See Capabilities Table at {dtl}")

    # throw exception for DF/CD. some of these pre-trapped by select_* functions but others escape, incl. cepa variants
    if (corl_type := method_algorithm_type(name).now) != "CONV":
        raise ValidationError(f"Invalid type {corl_type} for FNOCC energy through `run_cepa`. See Capabilities Table at {dtl}")

    for k, v in director[name].items():
        core.set_local_option("FNOCC", k.upper(), v)

    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, **kwargs)  # C1 certified

    if not core.get_option('FNOCC', 'USE_DF_INTS'):
        # Ensure IWL files have been written
        proc_util.check_iwl_file_from_scf_type(core.get_global_option('SCF_TYPE'), ref_wfn)
    else:
        core.print_out("  Constructing Basis Sets for FISAPT...\n\n")
        scf_aux_basis = core.BasisSet.build(ref_wfn.molecule(), "DF_BASIS_SCF",
                                            core.get_option("SCF", "DF_BASIS_SCF"),
                                            "JKFIT", core.get_global_option('BASIS'),
                                            puream=ref_wfn.basisset().has_puream())
        ref_wfn.set_basisset("DF_BASIS_SCF", scf_aux_basis)

    fnocc_wfn = core.fnocc(ref_wfn)

    # one-electron properties
    cepa_level = director[name]["cepa_level"]
    if core.get_option('FNOCC', 'DIPMOM'):
        if cepa_level in ['cepa(1)', 'cepa(3)']:
            core.print_out("""\n    Error: one-electron properties not implemented for %s\n\n""" % name)
        elif core.get_option('FNOCC', 'NAT_ORBS'):
            core.print_out("""\n    Error: one-electron properties not implemented for %s\n\n""" % name)
        else:
            p4util.oeprop(fnocc_wfn, 'DIPOLE', 'QUADRUPOLE', 'MULLIKEN_CHARGES', 'NO_OCCUPATIONS', title=cepa_level.upper())

    # Shove variables into global space
    for k, v in fnocc_wfn.variables().items():
        core.set_variable(k, v)

    optstash.restore()
    return fnocc_wfn


def run_detcas(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    determinant-based multireference wavefuncations,
    namely CASSCF and RASSCF.
    """
    optstash = p4util.OptionsState(
        ['DETCI', 'WFN'],
        ['SCF_TYPE'],
        ['ONEPDM'],
        ['OPDM_RELAX']
        )

    user_ref = core.get_option('DETCI', 'REFERENCE')
    if user_ref not in ['RHF', 'ROHF']:
        raise ValidationError('Reference %s for DETCI is not available.' % user_ref)

    if name == 'rasscf':
        core.set_local_option('DETCI', 'WFN', 'RASSCF')
    elif name == 'casscf':
        core.set_local_option('DETCI', 'WFN', 'CASSCF')
    else:
        raise ValidationError("Run DETCAS: Name %s not understood" % name)

    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:

        ref_optstash = p4util.OptionsState(
            ['SCF_TYPE'],
            ['DF_BASIS_SCF'],
            ['DF_BASIS_MP2'],
            ['ONEPDM'],
            ['OPDM_RELAX']
            )

        # No real reason to do a conventional guess
        if not core.has_global_option_changed('SCF_TYPE'):
            core.set_global_option('SCF_TYPE', 'DF')

        # If RHF get MP2 NO's
        # Why doesnt this work for conv?
        if (('DF' in core.get_global_option('SCF_TYPE')) and (user_ref == 'RHF') and
                    (core.get_option('DETCI', 'MCSCF_TYPE') in ['DF', 'AO']) and
                    (core.get_option("DETCI", "MCSCF_GUESS") == "MP2")):
            core.set_global_option('ONEPDM', True)
            core.set_global_option('OPDM_RELAX', False)
            ref_wfn = run_dfmp2_gradient(name, **kwargs)
        else:
            ref_wfn = scf_helper(name, **kwargs)

        # Ensure IWL files have been written
        if (core.get_option('DETCI', 'MCSCF_TYPE') == 'CONV'):
            mints = core.MintsHelper(ref_wfn.basisset())
            mints.set_print(1)
            mints.integrals()

        ref_optstash.restore()

    # The DF case
    if core.get_option('DETCI', 'MCSCF_TYPE') == 'DF':
        if not core.has_global_option_changed('SCF_TYPE'):
            core.set_global_option('SCF_TYPE', 'DF')

        core.print_out("  Constructing Basis Sets for MCSCF...\n\n")
        scf_aux_basis = core.BasisSet.build(ref_wfn.molecule(), "DF_BASIS_SCF",
                                            core.get_option("SCF", "DF_BASIS_SCF"),
                                            "JKFIT", core.get_global_option('BASIS'),
                                            puream=ref_wfn.basisset().has_puream())
        ref_wfn.set_basisset("DF_BASIS_SCF", scf_aux_basis)

    # The AO case
    elif core.get_option('DETCI', 'MCSCF_TYPE') == 'AO':
        if not core.has_global_option_changed('SCF_TYPE'):
            core.set_global_option('SCF_TYPE', 'DIRECT')

    # The conventional case
    elif core.get_option('DETCI', 'MCSCF_TYPE') == 'CONV':
        if not core.has_global_option_changed('SCF_TYPE'):
            core.set_global_option('SCF_TYPE', 'PK')

        # Ensure IWL files have been written
        proc_util.check_iwl_file_from_scf_type(core.get_global_option('SCF_TYPE'), ref_wfn)
    else:
        raise ValidationError("Run DETCAS: MCSCF_TYPE %s not understood." % str(core.get_option('DETCI', 'MCSCF_TYPE')))


    # Second-order SCF requires non-symmetric density matrix support
    if core.get_option('DETCI', 'MCSCF_ALGORITHM') in ['AH', 'OS']:
        proc_util.check_non_symmetric_jk_density("Second-order MCSCF")

    ciwfn = mcscf.mcscf_solver(ref_wfn)

    # We always would like to print a little dipole information
    oeprop = core.OEProp(ciwfn)
    oeprop.set_title(name.upper())
    oeprop.add("DIPOLE")
    oeprop.compute()
    ciwfn.oeprop = oeprop
    core.set_variable("CURRENT DIPOLE", core.variable(name.upper() + " DIPOLE"))

    # Shove variables into global space
    for k, v in ciwfn.variables().items():
        core.set_variable(k, v)

    optstash.restore()
    return ciwfn

def run_efp(name, **kwargs):
    """Function encoding sequence of module calls for a pure EFP
    computation (ignore any QM atoms).

    """

    efp_molecule = kwargs.get('molecule', core.get_active_molecule())
    try:
        efpobj = efp_molecule.EFP
    except AttributeError:
        raise ValidationError("""Method 'efp' not available without EFP fragments in molecule""")

    # print efp geom in [A]
    core.print_out(efpobj.banner())
    core.print_out(efpobj.geometry_summary(units_to_bohr=constants.bohr2angstroms))

    # set options
    # * 'chtr', 'qm_exch', 'qm_disp', 'qm_chtr' may be enabled in a future libefp release
    efpopts = {}
    for opt in ['elst', 'exch', 'ind', 'disp',
                   'elst_damping', 'ind_damping', 'disp_damping']:
        psiopt = 'EFP_' + opt.upper()
        if core.has_option_changed('EFP', psiopt):
            efpopts[opt] = core.get_option('EFP', psiopt)
    efpopts['qm_elst'] = False
    efpopts['qm_ind'] = False
    efpobj.set_opts(efpopts, label='psi', append='psi')
    do_gradient = core.get_option('EFP', 'DERTYPE') == 'FIRST'

    # compute and report
    efpobj.compute(do_gradient=do_gradient)
    core.print_out(efpobj.energy_summary(label='psi'))

    ene = efpobj.get_energy(label='psi')
    core.set_variable('EFP ELST ENERGY', ene['electrostatic'] + ene['charge_penetration'] + ene['electrostatic_point_charges'])
    core.set_variable('EFP IND ENERGY', ene['polarization'])
    core.set_variable('EFP DISP ENERGY', ene['dispersion'])
    core.set_variable('EFP EXCH ENERGY', ene['exchange_repulsion'])
    core.set_variable('EFP TOTAL ENERGY', ene['total'])
    core.set_variable('CURRENT ENERGY', ene['total'])

    if do_gradient:
        core.print_out(efpobj.gradient_summary())

        torq = efpobj.get_gradient()
        torq = core.Matrix.from_array(np.asarray(torq).reshape(-1, 6))
        core.set_variable("EFP TORQUE", torq)  # P::e EFP

    return ene['total']
