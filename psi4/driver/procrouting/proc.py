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

"""Module with functions that encode the sequence of PSI module
calls for each of the *name* values of the energy(), optimize(),
response(), and frequency() function. *name* can be assumed lowercase by here.

"""
import os
import sys
import shutil
import subprocess
import warnings

import numpy as np
from qcelemental import constants

from psi4 import extras
from psi4 import core
from psi4.driver import p4util
from psi4.driver import qcdb
from psi4.driver import psifiles as psif
from psi4.driver.p4util.exceptions import ManagedMethodError, PastureRequiredError, ValidationError
#from psi4.driver.molutil import *
from psi4.driver.qcdb.basislist import corresponding_basis
# never import driver, wrappers, or aliases into this file

from .roa import run_roa
from . import proc_util
from . import empirical_dispersion
from . import dft
from . import mcscf
from . import response
from . import solvent


# ATTN NEW ADDITIONS!
# consult http://psicode.org/psi4manual/master/proc_py.html

def select_mp2(name, **kwargs):
    """Function selecting the algorithm for a MP2 energy call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    mtd_type = core.get_global_option('MP2_TYPE')
    module = core.get_global_option('QC_MODULE')
    # Considering only [df]occ/dfmp2/detci/fnocc

    # MP2_TYPE exists largely for py-side reasoning, so must manage it
    #   here rather than passing to c-side unprepared for validation

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
                func = run_detci
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
                       """byproduct of a CISD computation.\n  DETCI ROHF MP2 will produce non-standard results.\n""")

    if func is None:
        raise ManagedMethodError(['select_mp2', name, 'MP2_TYPE', mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_mp2_gradient(name, **kwargs):
    """Function selecting the algorithm for a MP2 gradient call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    mtd_type = core.get_global_option('MP2_TYPE')
    module = core.get_global_option('QC_MODULE')
    all_electron = (core.get_global_option('FREEZE_CORE') == "FALSE")
    # Considering only [df]occ/dfmp2

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
        raise ManagedMethodError(['select_mp2_gradient', name, 'MP2_TYPE', mtd_type, reference, module, all_electron])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_mp2_property(name, **kwargs):
    """Function selecting the algorithm for a MP2 property call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    mtd_type = core.get_global_option('MP2_TYPE')
    module = core.get_global_option('QC_MODULE')
    # Considering only dfmp2 for now

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
        raise ManagedMethodError(['select_mp2_property', name, 'MP2_TYPE', mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_omp2(name, **kwargs):
    """Function selecting the algorithm for an OMP2 energy call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    mtd_type = core.get_global_option('MP2_TYPE')
    module = core.get_global_option('QC_MODULE')
    # Considering only [df]occ

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
        raise ManagedMethodError(['select_omp2', name, 'MP2_TYPE', mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_omp2_gradient(name, **kwargs):
    """Function selecting the algorithm for an OMP2 gradient call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    mtd_type = core.get_global_option('MP2_TYPE')
    module = core.get_global_option('QC_MODULE')
    # Considering only [df]occ

    func = None
    if reference in ['RHF', 'UHF', 'ROHF', 'RKS', 'UKS']:
        if mtd_type == 'CONV':
            if module in ['', 'OCC']:
                func = run_occ_gradient
        elif mtd_type == 'DF':
            if module in ['', 'OCC']:
                func = run_dfocc_gradient

    if func is None:
        raise ManagedMethodError(['select_omp2_gradient', name, 'MP2_TYPE', mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_omp2_property(name, **kwargs):
    """Function selecting the algorithm for an OMP2 property call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    mtd_type = core.get_global_option('MP2_TYPE')
    module = core.get_global_option('QC_MODULE')
    # Considering only [df]occ

    func = None
    if reference in ['RHF', 'UHF', 'ROHF', 'RKS', 'UKS']:
        if mtd_type == 'DF':
            if module in ['', 'OCC']:
                func = run_dfocc_property

    if func is None:
        raise ManagedMethodError(['select_omp2_property', name, 'MP2_TYPE', mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_omp2p5_property(name, **kwargs):
    """Function selecting the algorithm for an OMP2.5 property call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    mtd_type = core.get_global_option('MP_TYPE')
    module = core.get_global_option('QC_MODULE')
    # Considering only [df]occ

    func = None
    if reference in ['RHF', 'UHF', 'ROHF', 'RKS', 'UKS']:
        if mtd_type == 'DF':
            if module in ['', 'OCC']:
                func = run_dfocc_property

    if func is None:
        raise ManagedMethodError(['select_omp2p5_property', name, 'MP_TYPE', mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_omp3_property(name, **kwargs):
    """Function selecting the algorithm for an OMP3 property call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    mtd_type = core.get_global_option('MP_TYPE')
    module = core.get_global_option('QC_MODULE')
    # Considering only [df]occ

    func = None
    if reference in ['RHF', 'UHF', 'ROHF', 'RKS', 'UKS']:
        if mtd_type == 'DF':
            if module in ['', 'OCC']:
                func = run_dfocc_property

    if func is None:
        raise ManagedMethodError(['select_omp3_property', name, 'MP_TYPE', mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_olccd_property(name, **kwargs):
    """Function selecting the algorithm for an OLCCD property call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    mtd_type = core.get_global_option('CC_TYPE')
    module = core.get_global_option('QC_MODULE')
    # Considering only [df]occ

    func = None
    if reference in ['RHF', 'UHF', 'ROHF', 'RKS', 'UKS']:
        if mtd_type == 'DF':
            if module in ['', 'OCC']:
                func = run_dfocc_property

    if func is None:
        raise ManagedMethodError(['select_olccd_property', name, 'CC_TYPE', mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_mp3(name, **kwargs):
    """Function selecting the algorithm for a MP3 energy call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    mtd_type = core.get_global_option('MP_TYPE') if core.has_global_option_changed("MP_TYPE") else "DF"
    module = core.get_global_option('QC_MODULE')
    # Considering only [df]occ/fnocc/detci

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
            if module == 'DETCI':  # no default for this case
                func = run_detci
            elif module in ['']:
                core.print_out("""\nThis method is available inefficiently as a """
                               """byproduct of a CISD computation.\n  Add "set """
                               """qc_module detci" to input to access this route.\n""")

    if func is None:
        raise ManagedMethodError(['select_mp3', name, 'MP_TYPE', mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_mp3_gradient(name, **kwargs):
    """Function selecting the algorithm for a MP3 gradient call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    mtd_type = core.get_global_option('MP_TYPE') if core.has_global_option_changed("MP_TYPE") else "DF"
    module = core.get_global_option('QC_MODULE')
    all_electron = (core.get_global_option('FREEZE_CORE') == "FALSE")
    # Considering only [df]occ

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
        raise ManagedMethodError(['select_mp3_gradient', name, 'MP_TYPE', mtd_type, reference, module, all_electron])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_omp3(name, **kwargs):
    """Function selecting the algorithm for an OMP3 energy call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    mtd_type = core.get_global_option('MP_TYPE')
    module = core.get_global_option('QC_MODULE')
    # Considering only [df]occ

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
        raise ManagedMethodError(['select_omp3', name, 'MP_TYPE', mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_omp3_gradient(name, **kwargs):
    """Function selecting the algorithm for an OMP3 gradient call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    mtd_type = core.get_global_option('MP_TYPE')
    module = core.get_global_option('QC_MODULE')
    # Considering only [df]occ

    func = None
    if reference in ['RHF', 'UHF', 'ROHF', 'RKS', 'UKS']:
        if mtd_type == 'CONV':
            if module in ['', 'OCC']:
                func = run_occ_gradient
        elif mtd_type == 'DF':
            if module in ['', 'OCC']:
                func = run_dfocc_gradient

    if func is None:
        raise ManagedMethodError(['select_omp3_gradient', name, 'MP_TYPE', mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_mp2p5(name, **kwargs):
    """Function selecting the algorithm for a MP2.5 energy call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    mtd_type = core.get_global_option('MP_TYPE') if core.has_global_option_changed("MP_TYPE") else "DF"
    module = core.get_global_option('QC_MODULE')
    # Considering only [df]occ

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
        raise ManagedMethodError(['select_mp2p5', name, 'MP_TYPE', mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_mp2p5_gradient(name, **kwargs):
    """Function selecting the algorithm for a MP2.5 gradient call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    mtd_type = core.get_global_option('MP_TYPE') if core.has_global_option_changed("MP_TYPE") else "DF"
    module = core.get_global_option('QC_MODULE')
    all_electron = (core.get_global_option('FREEZE_CORE') == "FALSE")
    # Considering only [df]occ

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
        raise ManagedMethodError(['select_mp2p5_gradient', name, 'MP_TYPE', mtd_type, reference, module, all_electron])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_omp2p5(name, **kwargs):
    """Function selecting the algorithm for an OMP2.5 energy call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    mtd_type = core.get_global_option('MP_TYPE')
    module = core.get_global_option('QC_MODULE')
    # Considering only [df]occ

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
        raise ManagedMethodError(['select_omp2p5', name, 'MP_TYPE', mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_omp2p5_gradient(name, **kwargs):
    """Function selecting the algorithm for an OMP2.5 gradient call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    mtd_type = core.get_global_option('MP_TYPE')
    module = core.get_global_option('QC_MODULE')
    # Considering only [df]occ

    func = None
    if reference in ['RHF', 'UHF', 'ROHF', 'RKS', 'UKS']:
        if mtd_type == 'CONV':
            if module in ['', 'OCC']:
                func = run_occ_gradient
        elif mtd_type == 'DF':
            if module in ['', 'OCC']:
                func = run_dfocc_gradient

    if func is None:
        raise ManagedMethodError(['select_omp2p5_gradient', name, 'MP_TYPE', mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_lccd(name, **kwargs):
    """Function selecting the algorithm for a LCCD energy call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    mtd_type = core.get_global_option('CC_TYPE')
    module = core.get_global_option('QC_MODULE')
    # Considering only [df]occ/fnocc

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
        raise ManagedMethodError(['select_lccd', name, 'CC_TYPE', mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)

def select_lccd_gradient(name, **kwargs):
    """Function selecting the algorithm for a LCCD gradient call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    mtd_type = core.get_global_option('CC_TYPE')
    module = core.get_global_option('QC_MODULE')
    all_electron = (core.get_global_option('FREEZE_CORE') == "FALSE")
    # Considering only [df]occ

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
        raise ManagedMethodError(['select_lccd_gradient', name, 'CC_TYPE', mtd_type, reference, module, all_electron])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_olccd(name, **kwargs):
    """Function selecting the algorithm for an OLCCD energy call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    mtd_type = core.get_global_option('CC_TYPE')
    module = core.get_global_option('QC_MODULE')
    # Considering only [df]occ

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
        raise ManagedMethodError(['select_olccd', name, 'CC_TYPE', mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_olccd_gradient(name, **kwargs):
    """Function selecting the algorithm for an OLCCD gradient call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    mtd_type = core.get_global_option('CC_TYPE')
    module = core.get_global_option('QC_MODULE')
    # Considering only [df]occ

    func = None
    if reference in ['RHF', 'UHF', 'ROHF', 'RKS', 'UKS']:
        if mtd_type == 'CONV':
            if module in ['', 'OCC']:
                func = run_occ_gradient
        elif mtd_type == 'DF':
            if module in ['', 'OCC']:
                func = run_dfocc_gradient

    if func is None:
        raise ManagedMethodError(['select_olccd_gradient', name, 'CC_TYPE', mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_fnoccsd(name, **kwargs):
    """Function selecting the algorithm for a FNO-CCSD energy call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    mtd_type = core.get_global_option('CC_TYPE')
    module = core.get_global_option('QC_MODULE')
    # Considering only fnocc

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
        raise ManagedMethodError(['select_fnoccsd', name, 'CC_TYPE', mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_ccsd(name, **kwargs):
    """Function selecting the algorithm for a CCSD energy call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    mtd_type = core.get_global_option('CC_TYPE')
    module = core.get_global_option('QC_MODULE')
    # Considering only [df]occ/ccenergy/detci/fnocc

    func = None
    if reference == 'RHF':
        if mtd_type == 'CONV':
            if module == 'FNOCC':
                func = run_fnocc
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
            if module in ['', 'CCENERGY']:
                func = run_ccenergy
    elif reference == 'ROHF':
        if mtd_type == 'CONV':
            if module == 'CCT3' and extras.addons("cct3"):
                import cct3
                func = cct3.run_cct3
            elif module in ['', 'CCENERGY']:
                func = run_ccenergy

    if func is None:
        raise ManagedMethodError(['select_ccsd', name, 'CC_TYPE', mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_ccsd_gradient(name, **kwargs):
    """Function selecting the algorithm for a CCSD gradient call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    mtd_type = core.get_global_option('CC_TYPE')
    module = core.get_global_option('QC_MODULE')
    # Considering only [df]occ/ccenergy

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
    elif reference == 'ROHF':
        if mtd_type == 'CONV':
            if module in ['', 'CCENERGY']:
                func = run_ccenergy_gradient

    if func is None:
        raise ManagedMethodError(['select_ccsd_gradient', name, 'CC_TYPE', mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_fnoccsd_t_(name, **kwargs):
    """Function selecting the algorithm for a FNO-CCSD(T) energy call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    mtd_type = core.get_global_option('CC_TYPE')
    module = core.get_global_option('QC_MODULE')
    # Considering only fnocc

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
        raise ManagedMethodError(['select_fnoccsd_t_', name, 'CC_TYPE', mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_ccsd_t_(name, **kwargs):
    """Function selecting the algorithm for a CCSD(T) energy call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    mtd_type = core.get_global_option('CC_TYPE')
    module = core.get_global_option('QC_MODULE')
    # Considering only [df]occ/ccenergy/fnocc

    func = None
    if reference == 'RHF':
        if mtd_type == 'CONV':
            if module == 'FNOCC':
                func = run_fnocc
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
    elif reference in ['UHF', 'ROHF']:
        if mtd_type == 'CONV':
            if module in ['', 'CCENERGY']:
                func = run_ccenergy

    if func is None:
        raise ManagedMethodError(['select_ccsd_t_', name, 'CC_TYPE', mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_ccsd_t__gradient(name, **kwargs):
    """Function selecting the algorithm for a CCSD(T) gradient call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    mtd_type = core.get_global_option('CC_TYPE')
    module = core.get_global_option('QC_MODULE')
    # Considering only ccenergy

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

    if func is None:
        raise ManagedMethodError(['select_ccsd_t__gradient', name, 'CC_TYPE', mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_ccsd_at_(name, **kwargs):
    """Function selecting the algorithm for a CCSD(AT) energy call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    mtd_type = core.get_global_option('CC_TYPE')
    module = core.get_global_option('QC_MODULE')
    # Considering only [df]occ/ccenergy

    func = None
    if reference == 'RHF':
        if mtd_type == 'CONV':
            if module in ['', 'CCENERGY']:
                func = run_ccenergy
        elif mtd_type == 'DF':
            if module in ['', 'OCC']:
                func = run_dfocc
        elif mtd_type == 'CD':
            if module in ['', 'OCC']:
                func = run_dfocc

    if func is None:
        raise ManagedMethodError(['select_ccsd_at_', name, 'CC_TYPE', mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_cisd(name, **kwargs):
    """Function selecting the algorithm for a CISD energy call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    mtd_type = core.get_global_option('CI_TYPE')
    module = core.get_global_option('QC_MODULE')
    # Considering only detci/fnocc

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
        raise ManagedMethodError(['select_cisd', name, 'CI_TYPE', mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_mp4(name, **kwargs):
    """Function selecting the algorithm for a MP4 energy call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    mtd_type = core.get_global_option('MP_TYPE')
    module = core.get_global_option('QC_MODULE')
    # Considering only detci/fnocc

    func = None
    if reference == 'RHF':
        if mtd_type == 'CONV':
            if module == 'DETCI':
                func = run_detci
            elif module in ['', 'FNOCC']:
                func = run_fnocc
    elif reference == 'ROHF':
        if mtd_type == 'CONV':
            if module == 'DETCI':  # no default for this case
                func = run_detci
            elif module in ['']:
                core.print_out("""\nThis method is available inefficiently as a """
                               """byproduct of a CISDT computation.\n  Add "set """
                               """qc_module detci" to input to access this route.\n""")

    if func is None:
        raise ManagedMethodError(['select_mp4', name, 'MP_TYPE', mtd_type, reference, module])

    if kwargs.pop('probe', False):
        return
    else:
        return func(name, **kwargs)


def select_adc2(name, **kwargs):
    """Function selecting the algorithm for ADC(2) excited state energy
    call and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    mtd_type = core.get_global_option('MP_TYPE')
    module = core.get_global_option('QC_MODULE')
    # Considering only adcc/adc

    # TODO Actually one should do selection on a couple of other options here
    #      as well, e.g. adcc supports frozen-core and frozen-virtual,
    #      spin-specific states or spin-flip methods.
    #      But as far as I (mfherbst) know the BUILTIN ADC routine only supports
    #      singlet states and without freezing some core or some virtual orbitals.

    func = None
    if reference == 'RHF':
        if mtd_type == 'CONV':
            if module == 'ADCC' and extras.addons("adcc"):
                func = run_adcc
            elif module in ['', 'BUILTIN']:
                func = run_adc

    if reference == 'UHF':
        if mtd_type == 'CONV':
            if module in ['ADCC', ''] and extras.addons("adcc"):
                func = run_adcc

    # Note: ROHF is theoretically available in adcc, but are not fully tested
    #       ... so will be added later.

    if func is None:
        raise ManagedMethodError(['select_adc2', name, 'MP_TYPE', mtd_type, reference, module])

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
    if (("DF" in core.get_global_option("SCF_TYPE")) or
            (core.get_option("SCF", "DF_SCF_GUESS") and (core.get_global_option("SCF_TYPE") == "DIRECT"))):
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

    # If EXTERN is set, then place that potential on the wfn
    if hasattr(core, "EXTERN"):
        wfn.set_potential_variable("C", core.EXTERN) # This is for the FSAPT procedure
        wfn.set_external_potential(core.EXTERN)

    elif 'external_potentials' in kwargs:
        # For FSAPT, we can take a dictionary of external potentials, e.g.,
        # external_potentials={'A': potA, 'B': potB, 'C': potC} (any optional)
        # For the dimer SAPT calculation, we need to account for the external potential
        # in all of the subsystems A, B, C. So we add them all in total_external_potential
        # and set the external potential to the dimer wave function
        total_external_potential = core.ExternalPotential()

        for frag in kwargs['external_potentials']:
            if frag.upper() in "ABC":
                wfn.set_potential_variable(frag.upper(), kwargs['external_potentials'][frag].extern)
                total_external_potential.appendCharges(kwargs['external_potentials'][frag].extern.getCharges())

            else:
                core.print_out("\n  Warning! Unknown key for the external_potentials argument: %s" %frag)

        wfn.set_external_potential(total_external_potential)

    return wfn


def scf_helper(name, post_scf=True, **kwargs):
    """Function serving as helper to SCF, choosing whether to cast
    up or just run SCF with a standard guess. This preserves
    previous SCF options set by other procedures (e.g., SAPT
    output file types for SCF).

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

    # PCM needs to be run w/o symmetry
    if core.get_option("SCF", "PCM"):
        c1_molecule = scf_molecule.clone()
        c1_molecule.reset_point_group('c1')
        c1_molecule.update_geometry()

        scf_molecule = c1_molecule
        core.print_out("""  PCM does not make use of molecular symmetry: """
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

    # sort out broken_symmetry settings.
    if 'brokensymmetry' in kwargs:
        multp = scf_molecule.multiplicity()
        if multp != 1:
            raise ValidationError('Broken symmetry is only for singlets.')
        if core.get_option('SCF', 'REFERENCE') not in ['UHF', 'UKS']:
            raise ValidationError("""You must specify 'set reference uhf' to use broken symmetry.""")
        do_broken = True
    else:
        do_broken = False

    if cast and read_orbitals:
        raise ValidationError("""Detected options to both cast and read orbitals""")

    if cast and do_broken:
        raise ValidationError("""Detected options to both cast and perform a broken symmetry computation""")

    if (core.get_option('SCF', 'STABILITY_ANALYSIS') == 'FOLLOW') and (core.get_option('SCF', 'REFERENCE') != 'UHF'):
        raise ValidationError("""Stability analysis root following is only available for UHF""")

    # broken set-up
    if do_broken:
        raise ValidationError("""Broken symmetry computations are not currently enabled.""")
        scf_molecule.set_multiplicity(3)
        core.print_out('\n')
        p4util.banner('  Computing high-spin triplet guess  ')
        core.print_out('\n')

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
    if cast or do_broken:
        # Cast or broken are special cases
        base_wfn = core.Wavefunction.build(scf_molecule, core.get_global_option('BASIS'))
        core.print_out("\n         ---------------------------------------------------------\n");
        if banner:
            core.print_out("         " + banner.center(58));
        if cast:
            core.print_out("         " + "SCF Castup computation".center(58));
        ref_wfn = scf_wavefunction_factory(name, base_wfn, core.get_option('SCF', 'REFERENCE'), **kwargs)
        core.set_legacy_wavefunction(ref_wfn)

        # Compute dftd3
        if hasattr(ref_wfn, "_disp_functor"):
            disp_energy = ref_wfn._disp_functor.compute_energy(ref_wfn.molecule())
            ref_wfn.set_variable("-D Energy", disp_energy)
        ref_wfn.compute_energy()

    # broken clean-up
    if do_broken:
        raise ValidationError("Broken Symmetry computations are temporarily disabled.")
        scf_molecule.set_multiplicity(1)
        core.set_local_option('SCF', 'GUESS', 'READ')
        core.print_out('\n')
        p4util.banner('  Computing broken symmetry solution from high-spin triplet guess  ')
        core.print_out('\n')

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
        core.print_out("\n         ---------------------------------------------------------\n");
        core.print_out("         " + banner.center(58));

    scf_wfn = scf_wavefunction_factory(name, base_wfn, core.get_option('SCF', 'REFERENCE'), **kwargs)
    core.set_legacy_wavefunction(scf_wfn)

    # The wfn from_file routine adds the npy suffix if needed, but we add it here so that
    # we can use os.path.isfile to query whether the file exists before attempting to read
    read_filename = scf_wfn.get_scratch_filename(180) + '.npy'

    if (core.get_option('SCF', 'GUESS') == 'READ') and os.path.isfile(read_filename):
        old_wfn = core.Wavefunction.from_file(read_filename)
        Ca_occ = old_wfn.Ca_subset("SO", "OCC")
        Cb_occ = old_wfn.Cb_subset("SO", "OCC")

        if old_wfn.molecule().schoenflies_symbol() != scf_molecule.schoenflies_symbol():
            raise ValidationError("Cannot compute projection of different symmetries.")

        if old_wfn.basisset().name() == scf_wfn.basisset().name():
            core.print_out("  Reading orbitals from file 180, no projection.\n\n")
            scf_wfn.guess_Ca(Ca_occ)
            scf_wfn.guess_Cb(Cb_occ)
        else:
            core.print_out("  Reading orbitals from file 180, projecting to new basis.\n\n")
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
        core.print_out("  Unable to find file 180, defaulting to SAD guess.\n")
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

    # Compute dftd3
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

        # Compute properties
        oeprop.compute()
        for obj in [core, scf_wfn]:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                # component qcvars can be retired at v1.5
                for xyz in 'XYZ':
                    obj.set_variable('CURRENT DIPOLE ' + xyz, obj.variable('SCF DIPOLE ' + xyz))
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

    # Write out orbitals and basis; Can be disabled, e.g., for findif displacements
    if kwargs.get('write_orbitals', True):
        write_filename = scf_wfn.get_scratch_filename(180)

        scf_wfn.to_file(write_filename)
        extras.register_numpy_file(write_filename)

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

    core.set_local_option('DCT', 'OPDM', 'true');
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
            raise ValidationError(f"""Invalid type '{corl_type}' for DFOCC""")

    if name in ['mp2', 'omp2']:
        core.set_local_option('DFOCC', 'WFN_TYPE', 'DF-OMP2')
        corl_type = core.get_global_option('MP2_TYPE')
    elif name in ['mp2.5']:
        core.set_local_option('DFOCC', 'WFN_TYPE', 'DF-OMP2.5')
        corl_type = core.get_global_option('MP_TYPE') if core.has_global_option_changed("MP_TYPE") else "DF"
    elif name in ['omp2.5']:
        core.set_local_option('DFOCC', 'WFN_TYPE', 'DF-OMP2.5')
        corl_type = core.get_global_option('MP_TYPE')
    elif name in ['mp3']:
        core.set_local_option('DFOCC', 'WFN_TYPE', 'DF-OMP3')
        corl_type = core.get_global_option('MP_TYPE') if core.has_global_option_changed("MP_TYPE") else "DF"
    elif name in ['omp3']:
        core.set_local_option('DFOCC', 'WFN_TYPE', 'DF-OMP3')
        corl_type = core.get_global_option('MP_TYPE')
    elif name in ['lccd', 'olccd']:
        core.set_local_option('DFOCC', 'WFN_TYPE', 'DF-OLCCD')
        corl_type = core.get_global_option('CC_TYPE')

    elif name == 'ccd':
        core.set_local_option('DFOCC', 'WFN_TYPE', 'DF-CCD')
        corl_type = core.get_global_option('CC_TYPE')
    elif name == 'ccsd':
        core.set_local_option('DFOCC', 'WFN_TYPE', 'DF-CCSD')
        corl_type = core.get_global_option('CC_TYPE')
    elif name == 'ccsd(t)':
        core.set_local_option('DFOCC', 'WFN_TYPE', 'DF-CCSD(T)')
        corl_type = core.get_global_option('CC_TYPE')
    elif name == 'ccsd(at)':
        core.set_local_option('DFOCC', 'CC_LAMBDA', 'TRUE')
        core.set_local_option('DFOCC', 'WFN_TYPE', 'DF-CCSD(AT)')
        corl_type = core.get_global_option('CC_TYPE')
    elif name == 'dfocc':
        pass
    else:
        raise ValidationError('Unidentified method %s' % (name))

    set_cholesky_from(corl_type)

    # conventional vs. optimized orbitals
    if name in ['mp2', 'mp2.5', 'mp3', 'lccd',
                     'ccd', 'ccsd', 'ccsd(t)', 'ccsd(at)']:
        core.set_local_option('DFOCC', 'ORB_OPT', 'FALSE')
    elif name in ['omp2', 'omp2.5', 'omp3', 'olccd']:
        core.set_local_option('DFOCC', 'ORB_OPT', 'TRUE')

    core.set_local_option('DFOCC', 'DO_SCS', 'FALSE')
    core.set_local_option('DFOCC', 'DO_SOS', 'FALSE')
    core.set_local_option('SCF', 'DF_INTS_IO', 'SAVE')

    if name in ["mp2.5", "mp3"] and not core.has_global_option_changed("MP_TYPE"):
        core.print_out(f"    Information: {name.upper()} default algorithm changed to DF in August 2020. Use `set mp_type conv` for previous behavior.\n")

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
    if name in ['mp2', 'omp2', 'mp2.5', 'mp3', 'lccd',]:
        for k, v in dfocc_wfn.variables().items():
            core.set_variable(k, v)

    optstash.restore()
    return dfocc_wfn


def run_dfocc_gradient(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a density-fitted (non-)orbital-optimized MPN or CC computation.

    """
    optstash = p4util.OptionsState(
        ['SCF', 'DF_INTS_IO'],
        ['REFERENCE'],
        ['DFOCC', 'WFN_TYPE'],
        ['DFOCC', 'ORB_OPT'],
        ['DFOCC', 'CC_LAMBDA'],
        ['GLOBALS', 'DERTYPE'])


    proc_util.check_disk_df(name.upper(), optstash)

    if core.get_global_option('SCF_TYPE') != 'DISK_DF':
        raise ValidationError('DFOCC gradients need DF-SCF reference.')

    if name in ['mp2', 'omp2']:
        core.set_local_option('DFOCC', 'WFN_TYPE', 'DF-OMP2')
    elif name in ['mp2.5', 'omp2.5']:
        core.set_local_option('DFOCC', 'WFN_TYPE', 'DF-OMP2.5')
    elif name in ['mp3', 'omp3']:
        core.set_local_option('DFOCC', 'WFN_TYPE', 'DF-OMP3')
    elif name in ['lccd', 'olccd']:
        core.set_local_option('DFOCC', 'WFN_TYPE', 'DF-OLCCD')
    elif name in ['ccd']:
        core.set_local_option('DFOCC', 'WFN_TYPE', 'DF-CCD')
        core.set_local_option('DFOCC', 'CC_LAMBDA', 'TRUE')
    elif name in ['ccsd']:
        core.set_local_option('DFOCC', 'WFN_TYPE', 'DF-CCSD')
        core.set_local_option('DFOCC', 'CC_LAMBDA', 'TRUE')
    elif name in ['ccsd(t)']:
        core.set_local_option('DFOCC', 'WFN_TYPE', 'DF-CCSD(T)')
        core.set_local_option('DFOCC', 'CC_LAMBDA', 'TRUE')
    else:
        raise ValidationError('Unidentified method %s' % (name))

    if name in ['mp2', 'mp2.5', 'mp3', 'lccd', 'ccd', 'ccsd', 'ccsd(t)']:
        core.set_local_option('DFOCC', 'ORB_OPT', 'FALSE')
    elif name in ['omp2', 'omp2.5', 'omp3', 'olccd']:
        core.set_local_option('DFOCC', 'ORB_OPT', 'TRUE')

    core.set_global_option('DERTYPE', 'FIRST')
    core.set_local_option('DFOCC', 'DO_SCS', 'FALSE')
    core.set_local_option('DFOCC', 'DO_SOS', 'FALSE')
    core.set_local_option('SCF', 'DF_INTS_IO', 'SAVE')

    if name in ["mp2.5", "mp3"] and not core.has_global_option_changed("MP_TYPE"):
        core.print_out(f"    Information: {name.upper()} default algorithm changed to DF in August 2020. Use `set mp_type conv` for previous behavior.\n")

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
    if name in ['mp2', 'mp2.5', 'mp3', 'lccd', 'ccsd', 'omp2']:
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
    an density-fitted orbital-optimized MP2 computation

    """
    optstash = p4util.OptionsState(
        ['SCF', 'DF_INTS_IO'],
        ['DF_BASIS_SCF'],
        ['DIE_IF_NOT_CONVERGED'],
        ['MAXITER'],
        ['DFOCC', 'ORB_OPT'],
        ['DFOCC', 'WFN_TYPE'],
        ['DFOCC', 'QCHF'],
        ['DFOCC', 'E_CONVERGENCE'])

    core.set_local_option('DFOCC', 'ORB_OPT', 'TRUE')
    core.set_local_option('DFOCC', 'WFN_TYPE', 'QCHF')
    core.set_local_option('DFOCC', 'QCHF', 'TRUE')
    core.set_local_option('DFOCC', 'E_CONVERGENCE', 8)

    core.set_local_option('SCF', 'DF_INTS_IO', 'SAVE')
    core.set_local_option('SCF', 'DIE_IF_NOT_CONVERGED', 'FALSE')
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

    if name == 'mp2':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP2')
        core.set_local_option('OCC', 'ORB_OPT', 'FALSE')
        core.set_local_option('OCC', 'SPIN_SCALE_TYPE', 'NONE')
    elif name == 'scs-mp2':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP2')
        core.set_local_option('OCC', 'ORB_OPT', 'FALSE')
        core.set_local_option('OCC', 'SPIN_SCALE_TYPE', 'SCS')
    elif name == 'scs(n)-mp2':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP2')
        core.set_local_option('OCC', 'ORB_OPT', 'FALSE')
        core.set_local_option('OCC', 'SPIN_SCALE_TYPE', 'SCSN')
    elif name == 'scs-mp2-vdw':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP2')
        core.set_local_option('OCC', 'ORB_OPT', 'FALSE')
        core.set_local_option('OCC', 'SPIN_SCALE_TYPE', 'SCSVDW')
    elif name == 'sos-mp2':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP2')
        core.set_local_option('OCC', 'ORB_OPT', 'FALSE')
        core.set_local_option('OCC', 'SPIN_SCALE_TYPE', 'SOS')
    elif name == 'sos-pi-mp2':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP2')
        core.set_local_option('OCC', 'ORB_OPT', 'FALSE')
        core.set_local_option('OCC', 'SPIN_SCALE_TYPE', 'SOSPI')
    elif name == 'custom-scs-mp2':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP2')
        core.set_local_option('OCC', 'ORB_OPT', 'FALSE')
        core.set_local_option('OCC', 'SPIN_SCALE_TYPE', 'CUSTOM')

    elif name == 'omp2':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP2')
        core.set_local_option('OCC', 'ORB_OPT', 'TRUE')
        core.set_local_option('OCC', 'SPIN_SCALE_TYPE', 'NONE')
    elif name == 'scs-omp2':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP2')
        core.set_local_option('OCC', 'ORB_OPT', 'TRUE')
        core.set_local_option('OCC', 'SPIN_SCALE_TYPE', 'SCS')
    elif name == 'sos-omp2':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP2')
        core.set_local_option('OCC', 'ORB_OPT', 'TRUE')
        core.set_local_option('OCC', 'SPIN_SCALE_TYPE', 'SOS')
    elif name == 'custom-scs-omp2':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP2')
        core.set_local_option('OCC', 'ORB_OPT', 'TRUE')
        core.set_local_option('OCC', 'SPIN_SCALE_TYPE', 'CUSTOM')

    elif name == 'mp2.5':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP2.5')
        core.set_local_option('OCC', 'ORB_OPT', 'FALSE')
        core.set_local_option('OCC', 'SPIN_SCALE_TYPE', 'NONE')
    elif name == 'custom-scs-mp2.5':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP2.5')
        core.set_local_option('OCC', 'ORB_OPT', 'FALSE')
        core.set_local_option('OCC', 'SPIN_SCALE_TYPE', 'CUSTOM')

    elif name == 'omp2.5':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP2.5')
        core.set_local_option('OCC', 'ORB_OPT', 'TRUE')
        core.set_local_option('OCC', 'SPIN_SCALE_TYPE', 'NONE')
    elif name == 'custom-scs-omp2.5':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP2.5')
        core.set_local_option('OCC', 'ORB_OPT', 'TRUE')
        core.set_local_option('OCC', 'SPIN_SCALE_TYPE', 'CUSTOM')

    elif name == 'mp3':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP3')
        core.set_local_option('OCC', 'ORB_OPT', 'FALSE')
        core.set_local_option('OCC', 'SPIN_SCALE_TYPE', 'NONE')
    elif name == 'scs-mp3':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP3')
        core.set_local_option('OCC', 'ORB_OPT', 'FALSE')
        core.set_local_option('OCC', 'SPIN_SCALE_TYPE', 'SCS')
    elif name == 'custom-scs-mp3':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP3')
        core.set_local_option('OCC', 'ORB_OPT', 'FALSE')
        core.set_local_option('OCC', 'SPIN_SCALE_TYPE', 'CUSTOM')

    elif name == 'omp3':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP3')
        core.set_local_option('OCC', 'ORB_OPT', 'TRUE')
        core.set_local_option('OCC', 'SPIN_SCALE_TYPE', 'NONE')
    elif name == 'scs-omp3':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP3')
        core.set_local_option('OCC', 'ORB_OPT', 'TRUE')
        core.set_local_option('OCC', 'SPIN_SCALE_TYPE', 'SCS')
    elif name == 'sos-omp3':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP3')
        core.set_local_option('OCC', 'ORB_OPT', 'TRUE')
        core.set_local_option('OCC', 'SPIN_SCALE_TYPE', 'SOS')
    elif name == 'custom-scs-omp3':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP3')
        core.set_local_option('OCC', 'ORB_OPT', 'TRUE')
        core.set_local_option('OCC', 'SPIN_SCALE_TYPE', 'CUSTOM')

    elif name == 'lccd':
        core.set_local_option('OCC', 'WFN_TYPE', 'OCEPA')
        core.set_local_option('OCC', 'ORB_OPT', 'FALSE')
        core.set_local_option('OCC', 'SPIN_SCALE_TYPE', 'NONE')
    elif name == 'custom-scs-lccd':
        core.set_local_option('OCC', 'WFN_TYPE', 'OCEPA')
        core.set_local_option('OCC', 'ORB_OPT', 'FALSE')
        core.set_local_option('OCC', 'SPIN_SCALE_TYPE', 'CUSTOM')

    elif name == 'olccd':
        core.set_local_option('OCC', 'WFN_TYPE', 'OCEPA')
        core.set_local_option('OCC', 'ORB_OPT', 'TRUE')
        core.set_local_option('OCC', 'SPIN_SCALE_TYPE', 'NONE')
    elif name == 'custom-scs-olccd':
        core.set_local_option('OCC', 'WFN_TYPE', 'OCEPA')
        core.set_local_option('OCC', 'ORB_OPT', 'TRUE')
        core.set_local_option('OCC', 'SPIN_SCALE_TYPE', 'CUSTOM')

    else:
        raise ValidationError("""Invalid method %s""" % name)

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

    if name == 'mp2':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP2')
        core.set_local_option('OCC', 'ORB_OPT', 'FALSE')
    elif name in ['omp2', 'conv-omp2']:
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP2')
        core.set_local_option('OCC', 'ORB_OPT', 'TRUE')

    elif name == 'mp2.5':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP2.5')
        core.set_local_option('OCC', 'ORB_OPT', 'FALSE')
    elif name == 'omp2.5':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP2.5')
        core.set_local_option('OCC', 'ORB_OPT', 'TRUE')

    elif name == 'mp3':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP3')
        core.set_local_option('OCC', 'ORB_OPT', 'FALSE')
    elif name == 'omp3':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP3')
        core.set_local_option('OCC', 'ORB_OPT', 'TRUE')

    elif name == 'lccd':
        core.set_local_option('OCC', 'WFN_TYPE', 'OCEPA')
        core.set_local_option('OCC', 'ORB_OPT', 'FALSE')
    elif name == 'olccd':
        core.set_local_option('OCC', 'WFN_TYPE', 'OCEPA')
        core.set_local_option('OCC', 'ORB_OPT', 'TRUE')
    else:
        raise ValidationError("""Invalid method %s""" % name)

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

    if ref_wfn.basisset().has_ECP():
        core.print_out("\n\n  ==> Adding ECP gradient terms (computed numerically) <==\n")
        # Build a map of atom->ECP number
        old_print = ref_wfn.get_print()
        ref_wfn.set_print(0)
        delta = 0.0001
        natom = ref_wfn.molecule().natom()
        mints = core.MintsHelper(ref_wfn)
        ecpgradmat = core.Matrix("ECP Gradient", natom, 3)
        ecpgradmat.zero()
        ecpgrad = np.asarray(ecpgradmat)
        Dmat = ref_wfn.Da_subset("AO")
        Dmat.add(ref_wfn.Db_subset("AO"))
        def displaced_energy(atom, displacement):
            mints.basisset().move_atom(atom, displacement)
            E = Dmat.vector_dot(mints.ao_ecp())
            mints.basisset().move_atom(atom, -1*displacement)
            return E

        for atom in range(natom):
            for xyz in range(3):
                transvec = core.Vector3(0.0)
                transvec[xyz] += delta
                # +1 displacement
                Ep1 = displaced_energy(atom,  1*transvec)
                # -1 displacement
                Em1 = displaced_energy(atom, -1*transvec)
                # +2 displacement
                Ep2 = displaced_energy(atom,  2*transvec)
                # -2 displacement
                Em2 = displaced_energy(atom, -2*transvec)
                # Evaluate
                ecpgrad[atom, xyz] = (Em2 + 8*Ep1 - 8*Em1 - Ep2) / (12*delta)
        ecpgradmat.symmetrize_gradient(ref_wfn.molecule())
        ecpgradmat.print_atom_vector()
        grad.add(ecpgradmat)
        grad.print_atom_vector()
        ref_wfn.set_print(old_print)

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
        raise ValidationError("Only RHF/UHF Hessians are currently implemented. SCF_TYPE either CD or OUT_OF_CORE not supported")

    if hasattr(ref_wfn, "_disp_functor"):
        disp_hess = ref_wfn._disp_functor.compute_hessian(ref_wfn.molecule(), ref_wfn)
        ref_wfn.set_variable("-D Hessian", disp_hess)

    H = core.scfhess(ref_wfn)
    ref_wfn.set_hessian(H)

    # Clearly, add some logic when the reach of this fn expands
    ref_wfn.set_variable("HF TOTAL HESSIAN", H)  # P::e SCF

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
    elif name == 'ccsd(at)':
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
    if (core.get_option('SCF', 'REFERENCE') == 'ROHF') and \
            ((name in ['ccsd(t)', 'ccsd(at)', 'cc2', 'cc3', 'eom-cc2', 'eom-cc3']) or
              core.get_option('CCTRANSORT', 'SEMICANONICAL')):
        ref_wfn.semicanonicalize()

    if core.get_global_option('RUN_CCTRANSORT'):
        core.cctransort(ref_wfn)
    else:
        try:
            from psi4.driver.pasture import addins
            addins.ccsort_transqt2(ref_wfn)
        except:
            raise PastureRequiredError("RUN_CCTRANSORT")


    ccwfn = core.ccenergy(ref_wfn)
    if core.get_global_option('PE'):
        ccwfn.pe_state = ref_wfn.pe_state

    if name == 'ccsd(at)':
        core.cchbar(ref_wfn)
        core.cclambda(ref_wfn)

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

    if core.get_global_option('FREEZE_CORE') == 'TRUE':
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
        except:
            raise PastureRequiredError("RUN_CCTRANSORT")

    while True:
        sort_func(ref_wfn)

        ref_wfn = core.ccenergy(ref_wfn)
        core.print_out('Brueckner convergence check: %s\n' % bool(core.variable('BRUECKNER CONVERGED')))
        if (core.variable('BRUECKNER CONVERGED') == True):
            break

        if bcc_iter_cnt >= core.get_option('CCENERGY', 'BCCD_MAXITER'):
            core.print_out("\n\nWarning! BCCD did not converge within the maximum number of iterations.")
            core.print_out("You can increase the number of BCCD iterations by changing BCCD_MAXITER.\n\n")
            break
        bcc_iter_cnt += 1

    if name == 'bccd(t)':
        core.cctriples(ref_wfn)

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
                                                verbose=core.get_option("SCF", "TDSCF_PRINT"))

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

    # Always must set SCF dipole (retire components at v1.5)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for cart in ["X", "Y", "Z"]:
            core.set_variable("SCF DIPOLE " + cart, core.variable(name + " DIPOLE " + cart))
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

    if name in ['ccsd', 'cc2', 'eom-ccsd', 'eom-cc2']:
        this_name = name.upper().replace('-', '_')
        core.set_global_option('WFN', this_name)
        ccwfn = run_ccenergy(name, **kwargs)
        core.set_global_option('WFN', this_name)
    else:
        raise ValidationError("""CC property name %s not recognized""" % name.upper())

    # Need cchbar for everything
    core.cchbar(ccwfn)

    # Need ccdensity at this point only for density-based props
    if n_one > 0 or n_two > 0:
        if name == 'eom-ccsd':
            core.set_global_option('WFN', 'EOM_CCSD')
            core.set_global_option('DERTYPE', 'NONE')
            core.set_global_option('ONEPDM', 'TRUE')
            core.cceom(ccwfn)
        elif name == 'eom-cc2':
            core.set_global_option('WFN', 'EOM_CC2')
            core.set_global_option('DERTYPE', 'NONE')
            core.set_global_option('ONEPDM', 'TRUE')
            core.cceom(ccwfn)
        core.set_global_option('DERTYPE', 'NONE')
        core.set_global_option('ONEPDM', 'TRUE')
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
        core.set_global_option('ONEPDM', 'TRUE')
        # Tight convergence unnecessary for transition properties
        core.set_local_option('CCLAMBDA', 'R_CONVERGENCE', 1e-4)
        core.set_local_option('CCEOM', 'R_CONVERGENCE', 1e-4)
        core.set_local_option('CCEOM', 'E_CONVERGENCE', 1e-5)
        core.cceom(ccwfn)
        core.cclambda(ccwfn)
        core.ccdensity(ccwfn)

    if n_one > 0:
        # call oe prop for GS density
        oe = core.OEProp(ccwfn)
        oe.set_title(name.upper())
        for oe_name in one:
            oe.add(oe_name.upper())
        oe.compute()
        # call oe prop for each ES density
        if name.startswith('eom'):
            # copy GS CC DIP/QUAD ... to CC ROOT 0 DIP/QUAD ... if we are doing multiple roots
            # retire components at v1.5
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                if 'dipole' in one:
                    core.set_variable("CC ROOT 0 DIPOLE X", core.variable("CC DIPOLE X"))
                    core.set_variable("CC ROOT 0 DIPOLE Y", core.variable("CC DIPOLE Y"))
                    core.set_variable("CC ROOT 0 DIPOLE Z", core.variable("CC DIPOLE Z"))
                if 'quadrupole' in one:
                    core.set_variable("CC ROOT 0 QUADRUPOLE XX", core.variable("CC QUADRUPOLE XX"))
                    core.set_variable("CC ROOT 0 QUADRUPOLE XY", core.variable("CC QUADRUPOLE XY"))
                    core.set_variable("CC ROOT 0 QUADRUPOLE XZ", core.variable("CC QUADRUPOLE XZ"))
                    core.set_variable("CC ROOT 0 QUADRUPOLE YY", core.variable("CC QUADRUPOLE YY"))
                    core.set_variable("CC ROOT 0 QUADRUPOLE YZ", core.variable("CC QUADRUPOLE YZ"))
                    core.set_variable("CC ROOT 0 QUADRUPOLE ZZ", core.variable("CC QUADRUPOLE ZZ"))
            if 'dipole' in one:
                core.set_variable("CC ROOT 0 DIPOLE", core.variable("CC DIPOLE"))
                # core.set_variable("CC ROOT n DIPOLE", core.variable("CC DIPOLE"))  # P::e CCENERGY
            if 'quadrupole' in one:
                core.set_variable("CC ROOT 0 QUADRUPOLE", core.variable("CC QUADRUPOLE"))
                # core.set_variable("CC ROOT n QUADRUPOLE", core.variable("CC QUADRUPOLE"))  # P::e CCENERGY

            n_root = sum(core.get_global_option("ROOTS_PER_IRREP"))
            for rn in range(n_root):
                oe.set_title("CC ROOT {}".format(rn + 1))
                Da = ccwfn.variable("CC ROOT {} Da".format(rn + 1))
                oe.set_Da_so(Da)
                if core.get_global_option("REFERENCE") == "UHF":
                    Db = ccwfn.variable("CC ROOT {} Db".format(rn + 1))
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

    if not 'DF' in core.get_global_option('SCF_TYPE'):
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


def run_adc_deprecated(*args, **kwargs):
     warnings.warn("The method 'adc' has been deprecated, please use 'adc2' instead."
                   "The method key 'adc' will be removed Psi4 1.6.", DeprecationWarning)
     return select_adc2(*args, **kwargs)


def run_adc(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    an algebraic diagrammatic construction calculation.

    .. caution:: Get rid of active molecule lines- should be handled in energy.

    """
    if core.get_option('ADC', 'REFERENCE') != 'RHF':
        raise ValidationError('ADC requires reference RHF')

    # Bypass the scf call if a reference wavefunction is given
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, **kwargs)

    # Ensure IWL files have been written
    proc_util.check_iwl_file_from_scf_type(core.get_global_option('SCF_TYPE'), ref_wfn)

    return core.adc(ref_wfn)


def run_adcc(name, **kwargs):
    """Prepare and run an ADC calculation in adcc, interpret the result and return
    as a wavefunction.

    """
    # TODO Maybe it would improve readability if this function was spilt
    #      up and the whole thing went to a separate file (like for sapt,
    #      interface_cfour.py, ...

    try:
        import adcc
        from adcc.backends import InvalidReference
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
        elif not kind in ["any", "spin_flip"]:
            raise ValidationError("For UHF references the only valid values for 'KIND' are "
                                  "'SPIN_FLIP' or 'ANY' and not '{}.".format(kind.upper()))
    elif not kind in ["singlet", "triplet", "any"]:
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
    if "core_orbitals" in kwargs and not "cvs" in name:
        raise ValidationError("The NUM_CORE_ORBITALS option needs to be set to '0' or absent "
                              "unless a CVS ADC method is requested.")
    if "cvs" in name and kwargs["kind"] in ["spin_flip"]:
        raise ValidationError("Spin-flip for CVS-ADC variants is not available.")

    #
    # Check for unsupported options
    #
    for option in ["PR", "NORM_TOLERANCE", "POLE_MAXITER", "SEM_MAXITER",
                   "NEWTON_CONVERGENCE", "MEMORY", "CACHELEVEL", "NUM_AMPS_PRINT"]:
        if core.has_option_changed("ADC", option):
            raise ValidationError(f"ADC backend adcc does not support option '{option}'")

    #
    # Launch the rocket
    #
    # Copy thread setup from psi4
    try:
        adcc.set_n_threads(core.get_num_threads())
    except AttributeError:
        # Before adcc 0.13.3:
        adcc.thread_pool.reinit(core.get_num_threads(), core.get_num_threads())

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
                              "{}".format(str(ex)))
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
    adc_wfn.set_variable("ADC ITERATIONS", state.n_iter)  # P::e ADC
    adc_wfn.set_variable(name + " excitation energies",
                         core.Matrix.from_array(state.excitation_energy.reshape(-1, 1)))
    adc_wfn.set_variable("number of excited states", len(state.excitation_energy))

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

    computed = {}
    if any(prop in properties for prop in ("TRANSITION_DIPOLE", "OSCILLATOR_STRENGTH")):
        data = state.transition_dipole_moment
        computed["Transition dipole moment (in a.u.)"] = data
        adc_wfn.set_variable(f"{name} transition dipoles", core.Matrix.from_array(data))

    if "OSCILLATOR_STRENGTH" in properties:
        if gauge == "velocity":
            data = state.oscillator_strength_velocity.reshape(-1, 1)
        else:
            data = state.oscillator_strength.reshape(-1, 1)
        computed[f"Oscillator strength ({gauge} gauge)"] = data
        adc_wfn.set_variable(f"{name} oscillator strengths ({gauge_short})",
                             core.Matrix.from_array(data))

    if "ROTATIONAL_STRENGTH" in properties:
        data = state.rotatory_strength.reshape(-1, 1)
        computed["Rotational strength (velocity gauge)"] = data
        adc_wfn.set_variable(f"{name} rotational strengths (VEL)",
                             core.Matrix.from_array(data))

    if "DIPOLE" in properties:
        data = state.state_dipole_moment
        computed["State dipole moment (in a.u.)"] = data
        adc_wfn.set_variable(f"{name} state dipoles", core.Matrix.from_array(data))

    core.print_out("\nExcited state properties:\n")
    n_states = adc_wfn.variable("number of excited states")
    for i in range(int(n_states)):
        lines = [ind + f"Excited state  {i}"]
        for prop, data in sorted(computed.items()):
            lines += [ind + ind + format_vector(prop, data[i])]
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
    optstash = p4util.OptionsState(
        ['DETCI', 'WFN'],
        ['DETCI', 'MAX_NUM_VECS'],
        ['DETCI', 'MPN_ORDER_SAVE'],
        ['DETCI', 'MPN'],
        ['DETCI', 'FCI'],
        ['DETCI', 'EX_LEVEL'])

    if core.get_option('DETCI', 'REFERENCE') not in ['RHF', 'ROHF']:
        raise ValidationError('Reference %s for DETCI is not available.' %
            core.get_option('DETCI', 'REFERENCE'))

    if name == 'zapt':
        core.set_local_option('DETCI', 'WFN', 'ZAPTN')
        level = kwargs['level']
        maxnvect = int((level + 1) / 2) + (level + 1) % 2
        core.set_local_option('DETCI', 'MAX_NUM_VECS', maxnvect)
        if (level + 1) % 2:
            core.set_local_option('DETCI', 'MPN_ORDER_SAVE', 2)
        else:
            core.set_local_option('DETCI', 'MPN_ORDER_SAVE', 1)
    elif name in ['mp', 'mp2', 'mp3', 'mp4']:
        core.set_local_option('DETCI', 'WFN', 'DETCI')
        core.set_local_option('DETCI', 'MPN', 'TRUE')
        if name == 'mp2':
            level = 2
        elif name == 'mp3':
            level = 3
        elif name == 'mp4':
            level = 4
        else:
            level = kwargs['level']
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
    elif name == 'ci':
        core.set_local_option('DETCI', 'WFN', 'DETCI')
        level = kwargs['level']
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

    core.print_out("\t\t \"A good bug is a dead bug\" \n\n");
    core.print_out("\t\t\t - Starship Troopers\n\n");
    core.print_out("\t\t \"I didn't write FORTRAN.  That's the problem.\"\n\n");
    core.print_out("\t\t\t - Edward Valeev\n");

    if core.get_global_option("DIPMOM") and ("mp" not in name.lower()):
        # We always would like to print a little dipole information
        oeprop = core.OEProp(ciwfn)
        oeprop.set_title(name.upper())
        oeprop.add("DIPOLE")
        oeprop.compute()
        ciwfn.oeprop = oeprop
        # retire components in v1.5
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            core.set_variable("CURRENT DIPOLE X", core.variable(name.upper() + " DIPOLE X"))
            core.set_variable("CURRENT DIPOLE Y", core.variable(name.upper() + " DIPOLE Y"))
            core.set_variable("CURRENT DIPOLE Z", core.variable(name.upper() + " DIPOLE Z"))
        core.set_variable("CURRENT DIPOLE", core.variable(name.upper() + " DIPOLE"))

    ciwfn.cleanup_ci()
    ciwfn.cleanup_dpd()
    _clean_detci()

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

    if (core.get_option('SCF', 'REFERENCE') != 'RHF') and (name.upper() != "SAPT0"):
        raise ValidationError('Only SAPT0 supports a reference different from \"reference rhf\".')

    do_delta_mp2 = True if name.endswith('dmp2') else False
    do_empirical_disp = True if '-d' in name.lower() else False

    if do_empirical_disp:
        ## Make sure we are turning SAPT0 dispersion off
        core.set_local_option('SAPT', 'SAPT0_E10', True)
        core.set_local_option('SAPT', 'SAPT0_E20IND', True)
        core.set_local_option('SAPT', 'SAPT0_E20Disp', False)

    # raise Exception("")

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

    core.timer_on("SAPT: Dimer SCF")
    dimer_wfn = scf_helper('RHF', molecule=sapt_dimer, **kwargs)
    core.timer_off("SAPT: Dimer SCF")

    if do_delta_mp2:
        select_mp2(name, ref_wfn=dimer_wfn, **kwargs)
        mp2_corl_interaction_e = core.variable('MP2 CORRELATION ENERGY')

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
        select_mp2(name, ref_wfn=monomerA_wfn, **kwargs)
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
        select_mp2(name, ref_wfn=monomerB_wfn, **kwargs)
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
    if core.get_global_option("DF_BASIS_ELST") == "":
        dimer_wfn.set_basisset("DF_BASIS_ELST", aux_basis)
    else:
        aux_basis = core.BasisSet.build(dimer_wfn.molecule(), "DF_BASIS_ELST", core.get_global_option("DF_BASIS_ELST"),
                                        "RIFIT", core.get_global_option("BASIS"))
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

    # Check to see if we really need to run the SCF code.
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, **kwargs)
    vscf = core.variable('SCF TOTAL ENERGY')

    # The parse_arbitrary_order method provides us the following information
    # We require that level be provided. level is a dictionary
    # of settings to be passed to core.mrcc
    if not('level' in kwargs):
        raise ValidationError('level parameter was not provided.')

    level = kwargs['level']

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
        'PATH': ':'.join([os.path.abspath(x) for x in os.environ.get('PSIPATH', '').split(':') if x != '']) + \
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
    if core.has_option_changed('MRCC', 'MRCC_OMP_NUM_THREADS') == True:
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
            core.set_variable(m + ' TOTAL ENERGY', ene)
            core.set_variable(m + ' CORRELATION ENERGY', ene - vscf)
        except ValueError:
            continue

    # The last 'ene' in iface is the one the user requested.
    core.set_variable('CURRENT ENERGY', ene)
    core.set_variable('CURRENT CORRELATION ENERGY', ene - vscf)

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
    if (keep != False or ('path' in kwargs)):
        core.print_out('\nMRCC scratch files have been kept.\n')
        core.print_out('They can be found in ' + mrcc_tmpdir)

    # Dump iface contents to output
    core.print_out('\n')
    p4util.banner('Full results from MRCC')
    core.print_out('\n')
    core.print_out(iface_contents)

    return ref_wfn


def run_fnodfcc(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a DF-CCSD(T) computation.

    >>> set cc_type df
    >>> energy('fno-ccsd(t)')

    """
    kwargs = p4util.kwargs_lower(kwargs)

    # stash user options
    optstash = p4util.OptionsState(
        ['FNOCC', 'COMPUTE_TRIPLES'],
        ['FNOCC', 'DFCC'],
        ['FNOCC', 'NAT_ORBS'],
        ['FNOCC', 'RUN_CEPA'],
        ['FNOCC', 'DF_BASIS_CC'],
        ['SCF', 'DF_BASIS_SCF'],
        ['SCF', 'DF_INTS_IO'])

    core.set_local_option('FNOCC', 'DFCC', True)
    core.set_local_option('FNOCC', 'RUN_CEPA', False)

    # throw an exception for open-shells
    if core.get_option('SCF', 'REFERENCE') != 'RHF':
        raise ValidationError(f"""Error: {name} requires 'reference rhf'.""")

    def set_cholesky_from(mtd_type):
        type_val = core.get_global_option(mtd_type)
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
            raise ValidationError("""Invalid type '%s' for DFCC""" % type_val)

    # triples?
    if name == 'ccsd':
        core.set_local_option('FNOCC', 'COMPUTE_TRIPLES', False)
        set_cholesky_from('CC_TYPE')
    elif name == 'ccsd(t)':
        core.set_local_option('FNOCC', 'COMPUTE_TRIPLES', True)
        set_cholesky_from('CC_TYPE')
    elif name == 'fno-ccsd':
        core.set_local_option('FNOCC', 'COMPUTE_TRIPLES', False)
        core.set_local_option('FNOCC', 'NAT_ORBS', True)
        set_cholesky_from('CC_TYPE')
    elif name == 'fno-ccsd(t)':
        core.set_local_option('FNOCC', 'COMPUTE_TRIPLES', True)
        core.set_local_option('FNOCC', 'NAT_ORBS', True)
        set_cholesky_from('CC_TYPE')

    if core.get_global_option('SCF_TYPE') not in ['CD', 'DISK_DF']:
        raise ValidationError("""Invalid scf_type for DFCC.""")

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
    level = kwargs.get('level', 0)

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

    core.set_local_option('FNOCC', 'DFCC', False)
    core.set_local_option('FNOCC', 'RUN_CEPA', False)
    core.set_local_option('FNOCC', 'USE_DF_INTS', False)

    # which method?
    if name == 'ccsd':
        core.set_local_option('FNOCC', 'COMPUTE_TRIPLES', False)
        core.set_local_option('FNOCC', 'RUN_CCSD', True)
    elif name == 'ccsd(t)':
        core.set_local_option('FNOCC', 'COMPUTE_TRIPLES', True)
        core.set_local_option('FNOCC', 'RUN_CCSD', True)
    elif name == 'fno-ccsd':
        core.set_local_option('FNOCC', 'COMPUTE_TRIPLES', False)
        core.set_local_option('FNOCC', 'RUN_CCSD', True)
        core.set_local_option('FNOCC', 'NAT_ORBS', True)
    elif name == 'fno-ccsd(t)':
        core.set_local_option('FNOCC', 'COMPUTE_TRIPLES', True)
        core.set_local_option('FNOCC', 'RUN_CCSD', True)
        core.set_local_option('FNOCC', 'NAT_ORBS', True)
    elif name == 'qcisd':
        core.set_local_option('FNOCC', 'COMPUTE_TRIPLES', False)
        core.set_local_option('FNOCC', 'RUN_CCSD', False)
    elif name == 'qcisd(t)':
        core.set_local_option('FNOCC', 'COMPUTE_TRIPLES', True)
        core.set_local_option('FNOCC', 'RUN_CCSD', False)
    elif name == 'fno-qcisd':
        core.set_local_option('FNOCC', 'COMPUTE_TRIPLES', False)
        core.set_local_option('FNOCC', 'RUN_CCSD', False)
        core.set_local_option('FNOCC', 'NAT_ORBS', True)
    elif name == 'fno-qcisd(t)':
        core.set_local_option('FNOCC', 'COMPUTE_TRIPLES', True)
        core.set_local_option('FNOCC', 'NAT_ORBS', True)
        core.set_local_option('FNOCC', 'RUN_CCSD', False)
    elif name == 'mp2':
        core.set_local_option('FNOCC', 'RUN_MP2', True)
    elif name == 'fno-mp3':
        core.set_local_option('FNOCC', 'RUN_MP3', True)
        core.set_local_option('FNOCC', 'NAT_ORBS', True)
    elif name == 'fno-mp4':
        core.set_local_option('FNOCC', 'RUN_MP4', True)
        core.set_local_option('FNOCC', 'COMPUTE_MP4_TRIPLES', True)
        core.set_local_option('FNOCC', 'COMPUTE_TRIPLES', True)
        core.set_local_option('FNOCC', 'NAT_ORBS', True)
    elif name == 'mp4(sdq)':
        core.set_local_option('FNOCC', 'RUN_MP4', True)
        core.set_local_option('FNOCC', 'COMPUTE_MP4_TRIPLES', False)
        core.set_local_option('FNOCC', 'COMPUTE_TRIPLES', False)
    elif name == 'fno-mp4(sdq)':
        core.set_local_option('FNOCC', 'RUN_MP4', True)
        core.set_local_option('FNOCC', 'COMPUTE_MP4_TRIPLES', False)
        core.set_local_option('FNOCC', 'COMPUTE_TRIPLES', False)
        core.set_local_option('FNOCC', 'NAT_ORBS', True)
    elif name == 'mp3':
        core.set_local_option('FNOCC', 'RUN_MP3', True)
    elif name == 'mp4':
        core.set_local_option('FNOCC', 'RUN_MP4', True)
        core.set_local_option('FNOCC', 'COMPUTE_MP4_TRIPLES', True)
        core.set_local_option('FNOCC', 'COMPUTE_TRIPLES', True)

    # throw an exception for open-shells
    if core.get_option('SCF', 'REFERENCE') != 'RHF':
        raise ValidationError(f"""Error: {name} requires 'reference rhf'.""")

    # Bypass the scf call if a reference wavefunction is given
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, **kwargs)  # C1 certified

    if core.get_option('FNOCC', 'USE_DF_INTS') == False:
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

    # set current correlation energy and total energy.  only need to treat mpn here.
    if name in ["mp3", "fno-mp3"]:
        fnocc_wfn.set_variable("CURRENT ENERGY", fnocc_wfn.variable("MP3 TOTAL ENERGY"))
        fnocc_wfn.set_variable("CURRENT CORRELATION ENERGY", fnocc_wfn.variable("MP3 CORRELATION ENERGY"))
    elif name in ["mp4(sdq)", "fno-mp4(sdq)"]:
        fnocc_wfn.set_variable("CURRENT ENERGY", fnocc_wfn.variable("MP4(SDQ) TOTAL ENERGY"))
        fnocc_wfn.set_variable("CURRENT CORRELATION ENERGY", fnocc_wfn.variable("MP4(SDQ) CORRELATION ENERGY"))
    elif name in ["mp4", "fno-mp4"]:
        fnocc_wfn.set_variable("CURRENT ENERGY", fnocc_wfn.variable("MP4 TOTAL ENERGY"))
        fnocc_wfn.set_variable("CURRENT CORRELATION ENERGY", fnocc_wfn.variable("MP4 CORRELATION ENERGY"))

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

    # save user options
    optstash = p4util.OptionsState(
        ['TRANSQT2', 'WFN'],
        ['FNOCC', 'NAT_ORBS'],
        ['FNOCC', 'RUN_CEPA'],
        ['FNOCC', 'USE_DF_INTS'],
        ['FNOCC', 'CEPA_NO_SINGLES'])

    core.set_local_option('FNOCC', 'RUN_CEPA', True)
    core.set_local_option('FNOCC', 'USE_DF_INTS', False)

    # what type of cepa?
    if name in ['lccd', 'fno-lccd']:
        cepa_level = 'cepa(0)'
        core.set_local_option('FNOCC', 'CEPA_NO_SINGLES', True)
    elif name in ['cepa(0)', 'fno-cepa(0)', 'lccsd', 'fno-lccsd']:
        cepa_level = 'cepa(0)'
        core.set_local_option('FNOCC', 'CEPA_NO_SINGLES', False)
    elif name in ['cepa(1)', 'fno-cepa(1)']:
        cepa_level = 'cepa(1)'
    elif name in ['cepa(3)', 'fno-cepa(3)']:
        cepa_level = 'cepa(3)'
    elif name in ['acpf', 'fno-acpf']:
        cepa_level = 'acpf'
    elif name in ['aqcc', 'fno-aqcc']:
        cepa_level = 'aqcc'
    elif name in ['cisd', 'fno-cisd']:
        cepa_level = 'cisd'
    else:
        raise ValidationError("""Error: %s not implemented\n""" % name)

    core.set_local_option('FNOCC', 'CEPA_LEVEL', cepa_level.upper())

    if name in ['fno-lccd', 'fno-lccsd', 'fno-cepa(0)', 'fno-cepa(1)', 'fno-cepa(3)',
                'fno-acpf', 'fno-aqcc', 'fno-cisd']:
        core.set_local_option('FNOCC', 'NAT_ORBS', True)

    # throw an exception for open-shells
    if core.get_option('SCF', 'REFERENCE') != 'RHF':
        raise ValidationError("""Error: %s requires 'reference rhf'.""" % name)

    reference = core.get_option('SCF', 'REFERENCE')
    if core.get_global_option('CC_TYPE') != "CONV":
        raise ValidationError("""CEPA methods from FNOCC module require 'cc_type conv'.""")

    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, **kwargs)  # C1 certified

    if core.get_option('FNOCC', 'USE_DF_INTS') == False:
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
    # retire components by v1.5
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        core.set_variable("CURRENT DIPOLE X", core.variable(name.upper() + " DIPOLE X"))
        core.set_variable("CURRENT DIPOLE Y", core.variable(name.upper() + " DIPOLE Y"))
        core.set_variable("CURRENT DIPOLE Z", core.variable(name.upper() + " DIPOLE Z"))
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

