#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2017 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

"""Module with functions that encode the sequence of PSI module
calls for each of the *name* values of the energy(), optimize(),
response(), and frequency() function. *name* can be assumed lowercase by here.

"""
from __future__ import print_function
from __future__ import absolute_import
import shutil
import os
import subprocess
import re
import numpy as np

from psi4 import extras
from psi4.driver import p4util
from psi4.driver import qcdb
from psi4.driver.p4util.exceptions import *
from psi4.driver.molutil import *

from .roa import *
from . import proc_util
from . import empirical_dispersion
from . import dft_functional
from . import mcscf

# never import driver, wrappers, or aliases into this file

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
    # Considering only [df]occ/dfmp2

    func = None
    if reference == 'RHF':
        if mtd_type == 'CONV':
            if module in ['', 'OCC']:
                func = run_occ_gradient
        elif mtd_type == 'DF':
            if module == 'OCC':
                func = run_dfocc_gradient
            elif module in ['', 'DFMP2']:
                func = run_dfmp2_gradient
    elif reference == 'UHF':
        if mtd_type == 'CONV':
            if module in ['', 'OCC']:
                func = run_occ_gradient
        elif mtd_type == 'DF':
            if module in ['', 'OCC']:
                func = run_dfocc_gradient

    if func is None:
        raise ManagedMethodError(['select_mp2_gradient', name, 'MP2_TYPE', mtd_type, reference, module])

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


def select_mp3(name, **kwargs):
    """Function selecting the algorithm for a MP3 energy call
    and directing to specified or best-performance default modules.

    """
    reference = core.get_option('SCF', 'REFERENCE')
    mtd_type = core.get_global_option('MP_TYPE')
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
    mtd_type = core.get_global_option('MP_TYPE')
    module = core.get_global_option('QC_MODULE')
    # Considering only [df]occ

    func = None
    if reference == 'RHF':
        if mtd_type == 'CONV':
            if module in ['', 'OCC']:
                func = run_occ_gradient
        elif mtd_type == 'DF':
            if module in ['', 'OCC']:
                func = run_dfocc_gradient
    elif reference == 'UHF':
        if mtd_type == 'CONV':
            if module in ['', 'OCC']:
                func = run_occ_gradient
        elif mtd_type == 'DF':
            if module in ['', 'OCC']:
                func = run_dfocc_gradient

    if func is None:
        raise ManagedMethodError(['select_mp3_gradient', name, 'MP_TYPE', mtd_type, reference, module])

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
    mtd_type = core.get_global_option('MP_TYPE')
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
    mtd_type = core.get_global_option('MP_TYPE')
    module = core.get_global_option('QC_MODULE')
    # Considering only [df]occ

    func = None
    if reference in ['RHF', 'UHF']:
        if mtd_type == 'CONV':
            if module in ['', 'OCC']:
                func = run_occ_gradient
        elif mtd_type == 'DF':
            if module in ['', 'OCC']:
                func = run_dfocc_gradient

    if func is None:
        raise ManagedMethodError(['select_mp2p5_gradient', name, 'MP_TYPE', mtd_type, reference, module])

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
    # Considering only [df]occ

    func = None
    if reference in ['RHF', 'UHF']:
        if mtd_type == 'CONV':
            if module in ['', 'OCC']:
                func = run_occ_gradient
        elif mtd_type == 'DF':
            if module in ['', 'OCC']:
                func = run_dfocc_gradient

    if func is None:
        raise ManagedMethodError(['select_lccd_gradient', name, 'CC_TYPE', mtd_type, reference, module])

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
            if module == 'DETCI':
                func = run_detci
            elif module == 'FNOCC':
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
    elif reference == 'UHF':
        if mtd_type == 'CONV':
            if module in ['', 'CCENERGY']:
                func = run_ccenergy
    elif reference == 'ROHF':
        if mtd_type == 'CONV':
            if module == 'DETCI':
                func = run_detci
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
    if reference in ['RHF', 'UHF']:
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


def scf_wavefunction_factory(reference, ref_wfn, functional=None):
    """Builds the correct wavefunction from the provided information
    """

    if core.has_option_changed("SCF", "DFT_DISPERSION_PARAMETERS"):
        modified_disp_params = core.get_option("SCF", "DFT_DISPERSION_PARAMETERS")
    else:
        modified_disp_params = None

    # Figure out functional
    if functional is None:
        superfunc, disp_type = dft_functional.build_superfunctional(core.get_option("SCF", "DFT_FUNCTIONAL"))
    elif isinstance(functional, core.SuperFunctional):
        superfunc = functional
        disp_type = False
    elif isinstance(functional, (str, unicode)):
        superfunc, disp_type = dft_functional.build_superfunctional(functional)
    else:
        raise ValidationError("Functional %s is not understood" % str(functional))

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

    if disp_type:
        wfn._disp_functor = empirical_dispersion.EmpericalDispersion(disp_type[0], disp_type[1],
                                                                              tuple_params = modified_disp_params)
        wfn._disp_functor.print_out()

    # Set the multitude of SAD basis sets
    if (core.get_option("SCF", "SCF_TYPE") == "DF") or \
       (core.get_option("SCF", "DF_SCF_GUESS") and (core.get_option("SCF", "SCF_TYPE") == "DIRECT")):
        aux_basis = core.BasisSet.build(wfn.molecule(), "DF_BASIS_SCF",
                                        core.get_option("SCF", "DF_BASIS_SCF"),
                                        "JKFIT", core.get_global_option('BASIS'),
                                        puream=wfn.basisset().has_puream())
        wfn.set_basisset("DF_BASIS_SCF", aux_basis)

    if core.get_global_option("RELATIVISTIC") in ["X2C", "DKH"]:
        decon_basis = core.BasisSet.build(wfn.molecule(), "BASIS_RELATIVISTIC",
                                        core.get_option("SCF", "BASIS_RELATIVISTIC"),
                                        "DECON", core.get_global_option('BASIS'),
                                        puream=wfn.basisset().has_puream())
        wfn.set_basisset("BASIS_RELATIVISTIC", decon_basis)

    if (core.get_option("SCF", "GUESS") == "SAD"):
        sad_basis_list = core.BasisSet.build(wfn.molecule(), "ORBITAL",
                                             core.get_global_option("BASIS"),
                                             puream=wfn.basisset().has_puream(),
                                             return_atomlist=True)
        wfn.set_sad_basissets(sad_basis_list)

        if (core.get_option("SCF", "SAD_SCF_TYPE") == "DF"):
            sad_fitting_list = core.BasisSet.build(wfn.molecule(), "DF_BASIS_SAD",
                                                   core.get_option("SCF", "DF_BASIS_SAD"),
                                                   puream=wfn.basisset().has_puream(),
                                                   return_atomlist=True)
            wfn.set_sad_fitting_basissets(sad_fitting_list)

    return wfn


def scf_helper(name, **kwargs):
    """Function serving as helper to SCF, choosing whether to cast
    up or just run SCF with a standard guess. This preserves
    previous SCF options set by other procedures (e.g., SAPT
    output file types for SCF).

    """

    optstash = p4util.OptionsState(
        ['PUREAM'],
        ['BASIS'],
        ['QMEFP'],
        ['DF_BASIS_SCF'],
        ['SCF', 'GUESS'],
        ['SCF', 'DF_INTS_IO'],
        ['SCF', 'SCF_TYPE']  # Hack: scope gets changed internally with the Andy trick
    )

    optstash2 = p4util.OptionsState(
        ['BASIS'],
        ['DF_BASIS_SCF'],
        ['SCF', 'SCF_TYPE'],
        ['SCF', 'DF_INTS_IO'])

    # Grab a few kwargs
    use_c1 = kwargs.get('use_c1', False)
    scf_molecule = kwargs.get('molecule', core.get_active_molecule())
    read_orbitals = core.get_option('SCF', 'GUESS') is "READ"
    do_timer = kwargs.pop("do_timer", True)
    ref_wfn = kwargs.pop('ref_wfn', None)
    if ref_wfn is not None:
        raise Exception("Cannot supply a SCF wavefunction a ref_wfn.")

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

        # A use can set "BASIS_GUESS" to True and we default to 3-21G
        if cast is True:
            guessbasis = '3-21G'
        else:
            guessbasis = cast
        core.set_global_option('BASIS', guessbasis)

        castdf = core.get_option('SCF', 'SCF_TYPE') == 'DF'

        if core.has_option_changed('SCF', 'DF_BASIS_GUESS'):
            castdf = core.get_option('SCF', 'DF_BASIS_GUESS')
            if p4util.yes.match(str(castdf)):
                castdf = True
            elif p4util.no.match(str(castdf)):
                castdf = False

        if castdf:
            core.set_local_option('SCF', 'SCF_TYPE', 'DF')
            core.set_local_option('SCF', 'DF_INTS_IO', 'none')

            # Figure out the fitting basis set
            if castdf is True:
                core.set_global_option('DF_BASIS_SCF', '')
            elif isinstance(castdf, (unicode, str)):
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

    # broken set-up
    if do_broken:
        raise ValidationError("""Broken symmetry computations are not currently enabled.""")
        scf_molecule.set_multiplicity(3)
        core.print_out('\n')
        p4util.banner('  Computing high-spin triplet guess  ')
        core.print_out('\n')


    # If we force c1 copy the active molecule
    if use_c1:
        scf_molecule.update_geometry()
        if scf_molecule.schoenflies_symbol() != 'c1':
            core.print_out("""  A requested method does not make use of molecular symmetry: """
                           """further calculations in C1 point group.\n""")
            scf_molecule = scf_molecule.clone()
            scf_molecule.reset_point_group('c1')
            scf_molecule.fix_orientation(True)
            scf_molecule.fix_com(True)
            scf_molecule.update_geometry()

    # If GUESS is auto guess what it should be
    if core.get_option('SCF', 'GUESS') == "AUTO":
        if (core.get_option('SCF', 'REFERENCE') in ['RHF', 'RKS']) and \
                ((scf_molecule.natom() > 1) or core.get_option('SCF', 'SAD_FRAC_OCC')):
            core.set_local_option('SCF', 'GUESS', 'SAD')
        elif core.get_option('SCF', 'REFERENCE') in ['ROHF', 'ROKS', 'UHF', 'UKS']:
            core.set_local_option('SCF', 'GUESS', 'GWH')
        else:
            core.set_local_option('SCF', 'GUESS', 'CORE')

    if core.get_global_option('BASIS') == '':
        if name == 'hf3c':
            core.set_global_option('BASIS', 'minix')
        elif name == 'pbeh3c':
            core.set_global_option('BASIS', 'def2-msvp')

    # the FIRST scf call
    if cast or do_broken:
        # Cast or broken are special cases
        base_wfn = core.Wavefunction.build(scf_molecule, core.get_global_option('BASIS'))
        ref_wfn = scf_wavefunction_factory(core.get_option('SCF', 'REFERENCE'), base_wfn)
        core.set_legacy_wavefunction(ref_wfn)

        # Compute dftd3
        if "_disp_functor" in dir(ref_wfn):
            disp_energy = ref_wfn._disp_functor.compute_energy(ref_wfn.molecule())
            ref_wfn.set_variables("-D Energy", disp_energy)
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

        # Set to read and project, and reset bases to final ones
        optstash2.restore()
        core.set_local_option('SCF', 'GUESS', 'READ')

        # Print the banner for the standard operation
        core.print_out('\n')
        p4util.banner(name.upper())
        core.print_out('\n')


    # EFP preparation
    efp = core.get_active_efp()
    if efp.nfragments() > 0:
        core.set_legacy_molecule(scf_molecule)
        core.set_global_option('QMEFP', True)  # apt to go haywire if set locally to efp
        core.efp_set_options()
        efp.set_qm_atoms()
        efp.print_out()

    # the SECOND scf call
    base_wfn = core.Wavefunction.build(scf_molecule, core.get_global_option('BASIS'))
    scf_wfn = scf_wavefunction_factory(core.get_option('SCF', 'REFERENCE'), base_wfn)
    core.set_legacy_wavefunction(scf_wfn)

    fname = os.path.split(os.path.abspath(core.get_writer_file_prefix(scf_molecule.name())))[1]
    read_filename = os.path.join(core.get_environment("PSI_SCRATCH"), fname + ".180.npz")

    if (core.get_option('SCF', 'GUESS') == 'READ') and os.path.isfile(read_filename):
        data = np.load(read_filename)
        Ca_occ = core.Matrix.np_read(data, "Ca_occ")
        Cb_occ = core.Matrix.np_read(data, "Cb_occ")
        symmetry = str(data["symmetry"])
        basis_name = str(data["BasisSet"])

        if symmetry != scf_molecule.schoenflies_symbol():
            raise ValidationError("Cannot compute projection of different symmetries.")

        if basis_name == scf_wfn.basisset().name():
            core.print_out("  Reading orbitals from file 180, no projection.\n\n")
            scf_wfn.guess_Ca(Ca_occ)
            scf_wfn.guess_Cb(Cb_occ)
        else:
            core.print_out("  Reading orbitals from file 180, projecting to new basis.\n\n")

            puream = int(data["BasisSet PUREAM"])

            if ".gbs" in basis_name:
                basis_name = basis_name.split('/')[-1].replace('.gbs', '')

            old_basis = core.BasisSet.build(scf_molecule, "ORBITAL", basis_name, puream=puream)
            core.print_out("  Computing basis projection from %s to %s\n\n" % (basis_name, base_wfn.basisset().name()))

            nalphapi = core.Dimension.from_list(data["nalphapi"])
            nbetapi = core.Dimension.from_list(data["nbetapi"])
            pCa = scf_wfn.basis_projection(Ca_occ, nalphapi, old_basis, base_wfn.basisset())
            pCb = scf_wfn.basis_projection(Cb_occ, nbetapi, old_basis, base_wfn.basisset())
            scf_wfn.guess_Ca(pCa)
            scf_wfn.guess_Cb(pCb)

        # Strip off headers to only get R, RO, U, CU
        old_ref = str(data["reference"]).replace("KS", "").replace("HF", "")
        new_ref = core.get_option('SCF', 'REFERENCE').replace("KS", "").replace("HF", "")
        if old_ref != new_ref:
            scf_wfn.reset_occ(True)


    elif (core.get_option('SCF', 'GUESS') == 'READ') and not os.path.isfile(read_filename):
        core.print_out("  Unable to find file 180, defaulting to SAD guess.\n")
        core.set_local_option('SCF', 'GUESS', 'SAD')
        sad_basis_list = core.BasisSet.build(scf_wfn.molecule(), "ORBITAL",
                                             core.get_global_option("BASIS"),
                                             puream=scf_wfn.basisset().has_puream(),
                                             return_atomlist=True)
        scf_wfn.set_sad_basissets(sad_basis_list)

        if (core.get_option("SCF", "SAD_SCF_TYPE") == "DF"):
            sad_fitting_list = core.BasisSet.build(scf_wfn.molecule(), "DF_BASIS_SAD",
                                                   core.get_option("SCF", "DF_BASIS_SAD"),
                                                   puream=scf_wfn.basisset().has_puream(),
                                                   return_atomlist=True)
            scf_wfn.set_sad_fitting_basissets(sad_fitting_list)


    if cast:
        core.print_out("\n  Computing basis projection from %s to %s\n\n" % (ref_wfn.basisset().name(), base_wfn.basisset().name()))
        pCa = ref_wfn.basis_projection(ref_wfn.Ca(), ref_wfn.nalphapi(), ref_wfn.basisset(), scf_wfn.basisset())
        pCb = ref_wfn.basis_projection(ref_wfn.Cb(), ref_wfn.nbetapi(), ref_wfn.basisset(), scf_wfn.basisset())
        scf_wfn.guess_Ca(pCa)
        scf_wfn.guess_Cb(pCb)


    # Print basis set info
    if core.get_option("SCF", "PRINT_BASIS"):
        scf_wfn.basisset().print_detail_out()

    # Compute dftd3
    if "_disp_functor" in dir(scf_wfn):
        disp_energy = scf_wfn._disp_functor.compute_energy(scf_wfn.molecule())
        scf_wfn.set_variable("-D Energy", disp_energy)

    e_scf = scf_wfn.compute_energy()
    core.set_variable("SCF TOTAL ENERGY", e_scf)
    core.set_variable("CURRENT ENERGY", e_scf)
    core.set_variable("CURRENT REFERENCE ENERGY", e_scf)

    # We always would like to print a little dipole information
    if kwargs.get('scf_do_dipole', True):
        oeprop = core.OEProp(scf_wfn)
        oeprop.set_title("SCF")
        oeprop.add("DIPOLE")
        oeprop.compute()
        core.set_variable("CURRENT DIPOLE X", core.get_variable("SCF DIPOLE X"))
        core.set_variable("CURRENT DIPOLE Y", core.get_variable("SCF DIPOLE Y"))
        core.set_variable("CURRENT DIPOLE Z", core.get_variable("SCF DIPOLE Z"))

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

    # Write out orbitals and basis
    fname = os.path.split(os.path.abspath(core.get_writer_file_prefix(scf_molecule.name())))[1]
    write_filename = os.path.join(core.get_environment("PSI_SCRATCH"), fname + ".180.npz")
    data = {}
    data.update(scf_wfn.Ca().np_write(None, prefix="Ca"))
    data.update(scf_wfn.Cb().np_write(None, prefix="Cb"))

    Ca_occ = scf_wfn.Ca_subset("SO", "OCC")
    data.update(Ca_occ.np_write(None, prefix="Ca_occ"))

    Cb_occ = scf_wfn.Cb_subset("SO", "OCC")
    data.update(Cb_occ.np_write(None, prefix="Cb_occ"))

    data["reference"] = core.get_option('SCF', 'REFERENCE')
    data["nsoccpi"] = scf_wfn.soccpi().to_tuple()
    data["ndoccpi"] = scf_wfn.doccpi().to_tuple()
    data["nalphapi"] = scf_wfn.nalphapi().to_tuple()
    data["nbetapi"] = scf_wfn.nbetapi().to_tuple()
    data["symmetry"] = scf_molecule.schoenflies_symbol()
    data["BasisSet"] = scf_wfn.basisset().name()
    data["BasisSet PUREAM"] = scf_wfn.basisset().has_puream()
    np.savez(write_filename, **data)
    extras.register_numpy_file(write_filename)

    if do_timer:
        core.tstop()

    optstash.restore()
    return scf_wfn


def run_dcft(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a density cumulant functional theory calculation.

    """

    if (core.get_global_option('FREEZE_CORE') == 'TRUE'):
        raise ValidationError('Frozen core is not available for DCFT.')

    # Bypass the scf call if a reference wavefunction is given
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, **kwargs)

    if (core.get_global_option("DCFT_TYPE") == "DF"):
        aux_basis = core.BasisSet.build(ref_wfn.molecule(), "DF_BASIS_DCFT",
                                        core.get_global_option("DF_BASIS_DCFT"),
                                        "RIFIT", core.get_global_option("BASIS"))
        ref_wfn.set_basisset("DF_BASIS_DCFT", aux_basis)

        scf_aux_basis = core.BasisSet.build(ref_wfn.molecule(), "DF_BASIS_SCF",
                                            core.get_option("SCF", "DF_BASIS_SCF"),
                                            "JKFIT", core.get_global_option('BASIS'),
                                            puream=ref_wfn.basisset().has_puream())
        ref_wfn.set_basisset("DF_BASIS_SCF", scf_aux_basis)

    # Ensure IWL files have been written
    proc_util.check_iwl_file_from_scf_type(core.get_option('SCF', 'SCF_TYPE'), ref_wfn)

    dcft_wfn = core.dcft(ref_wfn)
    return dcft_wfn


def run_dcft_gradient(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    DCFT gradient calculation.

    """
    optstash = p4util.OptionsState(
        ['GLOBALS', 'DERTYPE'])


    core.set_global_option('DERTYPE', 'FIRST')
    dcft_wfn = run_dcft(name, **kwargs)

    derivobj = core.Deriv(dcft_wfn)
    derivobj.set_tpdm_presorted(True)
    grad = derivobj.compute()

    dcft_wfn.set_gradient(grad)

    optstash.restore()
    return dcft_wfn


def run_dfocc(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a density-fitted or Cholesky-decomposed
    (non-)orbital-optimized MPN or CC computation.

    """
    optstash = p4util.OptionsState(
        ['SCF', 'SCF_TYPE'],
        ['SCF', 'DF_INTS_IO'],
        ['DFOCC', 'WFN_TYPE'],
        ['DFOCC', 'ORB_OPT'],
        ['DFOCC', 'DO_SCS'],
        ['DFOCC', 'DO_SOS'],
        ['DFOCC', 'READ_SCF_3INDEX'],
        ['DFOCC', 'CHOLESKY'],
        ['DFOCC', 'CC_LAMBDA'])

    def set_cholesky_from(mtd_type):
        type_val = core.get_global_option(mtd_type)
        if type_val == 'DF':
            core.set_local_option('DFOCC', 'CHOLESKY', 'FALSE')
            # Alter default algorithm
            if not core.has_option_changed('SCF', 'SCF_TYPE'):
                core.set_global_option('SCF_TYPE', 'DF')
                core.print_out("""    SCF Algorithm Type (re)set to DF.\n""")
        elif type_val == 'CD':
            core.set_local_option('DFOCC', 'CHOLESKY', 'TRUE')
            # Alter default algorithm
            if not core.has_option_changed('SCF', 'SCF_TYPE'):
                core.set_global_option('SCF_TYPE', 'CD')
                core.print_out("""    SCF Algorithm Type (re)set to CD.\n""")
            if core.get_option('SCF', 'SCF_TYPE') != 'CD':
                core.set_local_option('DFOCC', 'READ_SCF_3INDEX', 'FALSE')
        else:
            raise ValidationError("""Invalid type '%s' for DFOCC""" % type_val)

    if name in ['mp2', 'omp2']:
        core.set_local_option('DFOCC', 'WFN_TYPE', 'DF-OMP2')
        set_cholesky_from('MP2_TYPE')
    elif name in ['mp2.5', 'omp2.5']:
        core.set_local_option('DFOCC', 'WFN_TYPE', 'DF-OMP2.5')
        set_cholesky_from('MP_TYPE')
    elif name in ['mp3', 'omp3']:
        core.set_local_option('DFOCC', 'WFN_TYPE', 'DF-OMP3')
        set_cholesky_from('MP_TYPE')
    elif name in ['lccd', 'olccd']:
        core.set_local_option('DFOCC', 'WFN_TYPE', 'DF-OLCCD')
        set_cholesky_from('CC_TYPE')

    elif name == 'ccd':
        core.set_local_option('DFOCC', 'WFN_TYPE', 'DF-CCD')
        set_cholesky_from('CC_TYPE')
    elif name == 'ccsd':
        core.set_local_option('DFOCC', 'WFN_TYPE', 'DF-CCSD')
        set_cholesky_from('CC_TYPE')
    elif name == 'ccsd(t)':
        core.set_local_option('DFOCC', 'WFN_TYPE', 'DF-CCSD(T)')
        set_cholesky_from('CC_TYPE')
    elif name == 'ccsd(at)':
        core.set_local_option('DFOCC', 'CC_LAMBDA', 'TRUE')
        core.set_local_option('DFOCC', 'WFN_TYPE', 'DF-CCSD(AT)')
        set_cholesky_from('CC_TYPE')
    elif name == 'dfocc':
        pass
    else:
        raise ValidationError('Unidentified method %s' % (name))

    # conventional vs. optimized orbitals
    if name in ['mp2', 'mp2.5', 'mp3', 'lccd',
                     'ccd', 'ccsd', 'ccsd(t)', 'ccsd(at)']:
        core.set_local_option('DFOCC', 'ORB_OPT', 'FALSE')
    elif name in ['omp2', 'omp2.5', 'omp3', 'olccd']:
        core.set_local_option('DFOCC', 'ORB_OPT', 'TRUE')

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

    optstash.restore()
    return dfocc_wfn


def run_dfocc_gradient(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a density-fitted (non-)orbital-optimized MPN or CC computation.

    """
    optstash = p4util.OptionsState(
        ['SCF', 'SCF_TYPE'],
        ['SCF', 'DF_INTS_IO'],
        ['REFERENCE'],
        ['DFOCC', 'WFN_TYPE'],
        ['DFOCC', 'ORB_OPT'],
        ['DFOCC', 'CC_LAMBDA'],
        ['GLOBALS', 'DERTYPE'])

    # Alter default algorithm
    if not core.has_option_changed('SCF', 'SCF_TYPE'):
        core.set_global_option('SCF_TYPE', 'DF')
        core.print_out("""    SCF Algorithm Type (re)set to DF.\n""")

    if core.get_option('SCF', 'SCF_TYPE') != 'DF':
        raise ValidationError('DFOCC gradients need DF-HF reference, for now.')

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
    else:
        raise ValidationError('Unidentified method %s' % (name))

    if name in ['mp2', 'mp2.5', 'mp3', 'lccd', 'ccd', 'ccsd']:
        core.set_local_option('DFOCC', 'ORB_OPT', 'FALSE')
    elif name in ['omp2', 'omp2.5', 'omp3', 'olccd']:
        core.set_local_option('DFOCC', 'ORB_OPT', 'TRUE')

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

    optstash.restore()
    return dfocc_wfn


def run_dfocc_property(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a density-fitted (non-)orbital-optimized MPN or CC computation.

    """
    optstash = p4util.OptionsState(
        ['SCF', 'SCF_TYPE'],
        ['SCF', 'DF_INTS_IO'],
        ['DFOCC', 'WFN_TYPE'],
        ['DFOCC', 'ORB_OPT'],
        ['DFOCC', 'OEPROP'])

    if name in ['mp2', 'omp2']:
        core.set_local_option('DFOCC', 'WFN_TYPE', 'DF-OMP2')
    else:
        raise ValidationError('Unidentified method ' % (name))

    # Alter default algorithm
    if not core.has_option_changed('SCF', 'SCF_TYPE'):
        core.set_global_option('SCF_TYPE', 'DF')
        core.print_out("""    SCF Algorithm Type (re)set to DF.\n""")

    if core.get_option('SCF', 'SCF_TYPE') != 'DF':
        raise ValidationError('DFOCC gradients need DF-HF reference, for now.')

    if name in ['mp2']:
        core.set_local_option('DFOCC', 'ORB_OPT', 'FALSE')
    elif name in ['omp2']:
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
    optstash = p4util.OptionsState(
        ['OCC', 'SCS_TYPE'],
        ['OCC', 'DO_SCS'],
        ['OCC', 'SOS_TYPE'],
        ['OCC', 'DO_SOS'],
        ['OCC', 'ORB_OPT'],
        ['OCC', 'WFN_TYPE'])

    if name == 'mp2':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP2')
        core.set_local_option('OCC', 'ORB_OPT', 'FALSE')
        core.set_local_option('OCC', 'DO_SCS', 'FALSE')
        core.set_local_option('OCC', 'DO_SOS', 'FALSE')
    elif name == 'omp2':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP2')
        core.set_local_option('OCC', 'ORB_OPT', 'TRUE')
        core.set_local_option('OCC', 'DO_SCS', 'FALSE')
        core.set_local_option('OCC', 'DO_SOS', 'FALSE')
    elif name == 'scs-omp2':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP2')
        core.set_local_option('OCC', 'ORB_OPT', 'TRUE')
        core.set_local_option('OCC', 'DO_SCS', 'TRUE')
        core.set_local_option('OCC', 'SCS_TYPE', 'SCS')
    elif name == 'scs(n)-omp2':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP2')
        core.set_local_option('OCC', 'ORB_OPT', 'TRUE')
        core.set_local_option('OCC', 'DO_SCS', 'TRUE')
        core.set_local_option('OCC', 'SCS_TYPE', 'SCSN')
    elif name == 'scs-omp2-vdw':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP2')
        core.set_local_option('OCC', 'ORB_OPT', 'TRUE')
        core.set_local_option('OCC', 'DO_SCS', 'TRUE')
        core.set_local_option('OCC', 'SCS_TYPE', 'SCSVDW')
    elif name == 'sos-omp2':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP2')
        core.set_local_option('OCC', 'ORB_OPT', 'TRUE')
        core.set_local_option('OCC', 'DO_SOS', 'TRUE')
        core.set_local_option('OCC', 'SOS_TYPE', 'SOS')
    elif name == 'sos-pi-omp2':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP2')
        core.set_local_option('OCC', 'ORB_OPT', 'TRUE')
        core.set_local_option('OCC', 'DO_SOS', 'TRUE')
        core.set_local_option('OCC', 'SOS_TYPE', 'SOSPI')

    elif name == 'mp2.5':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP2.5')
        core.set_local_option('OCC', 'ORB_OPT', 'FALSE')
        core.set_local_option('OCC', 'DO_SCS', 'FALSE')
        core.set_local_option('OCC', 'DO_SOS', 'FALSE')
    elif name == 'omp2.5':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP2.5')
        core.set_local_option('OCC', 'ORB_OPT', 'TRUE')
        core.set_local_option('OCC', 'DO_SCS', 'FALSE')
        core.set_local_option('OCC', 'DO_SOS', 'FALSE')

    elif name == 'mp3':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP3')
        core.set_local_option('OCC', 'ORB_OPT', 'FALSE')
        core.set_local_option('OCC', 'DO_SCS', 'FALSE')
        core.set_local_option('OCC', 'DO_SOS', 'FALSE')
    elif name == 'omp3':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP3')
        core.set_local_option('OCC', 'ORB_OPT', 'TRUE')
        core.set_local_option('OCC', 'DO_SCS', 'FALSE')
        core.set_local_option('OCC', 'DO_SOS', 'FALSE')
    elif name == 'scs-omp3':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP3')
        core.set_local_option('OCC', 'ORB_OPT', 'TRUE')
        core.set_local_option('OCC', 'DO_SCS', 'TRUE')
        core.set_local_option('OCC', 'SCS_TYPE', 'SCS')
    elif name == 'scs(n)-omp3':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP3')
        core.set_local_option('OCC', 'ORB_OPT', 'TRUE')
        core.set_local_option('OCC', 'DO_SCS', 'TRUE')
        core.set_local_option('OCC', 'SCS_TYPE', 'SCSN')
    elif name == 'scs-omp3-vdw':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP3')
        core.set_local_option('OCC', 'ORB_OPT', 'TRUE')
        core.set_local_option('OCC', 'DO_SCS', 'TRUE')
        core.set_local_option('OCC', 'SCS_TYPE', 'SCSVDW')
    elif name == 'sos-omp3':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP3')
        core.set_local_option('OCC', 'ORB_OPT', 'TRUE')
        core.set_local_option('OCC', 'DO_SOS', 'TRUE')
        core.set_local_option('OCC', 'SOS_TYPE', 'SOS')
    elif name == 'sos-pi-omp3':
        core.set_local_option('OCC', 'WFN_TYPE', 'OMP3')
        core.set_local_option('OCC', 'ORB_OPT', 'TRUE')
        core.set_local_option('OCC', 'DO_SOS', 'TRUE')
        core.set_local_option('OCC', 'SOS_TYPE', 'SOSPI')

    elif name == 'lccd':
        core.set_local_option('OCC', 'WFN_TYPE', 'OCEPA')
        core.set_local_option('OCC', 'ORB_OPT', 'FALSE')
        core.set_local_option('OCC', 'DO_SCS', 'FALSE')
        core.set_local_option('OCC', 'DO_SOS', 'FALSE')
    elif name == 'olccd':
        core.set_local_option('OCC', 'WFN_TYPE', 'OCEPA')
        core.set_local_option('OCC', 'ORB_OPT', 'TRUE')
        core.set_local_option('OCC', 'DO_SCS', 'FALSE')
        core.set_local_option('OCC', 'DO_SOS', 'FALSE')
    else:
        raise ValidationError("""Invalid method %s""" % name)

    # Bypass the scf call if a reference wavefunction is given
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, **kwargs)  # C1 certified

    # Ensure IWL files have been written
    proc_util.check_iwl_file_from_scf_type(core.get_option('SCF', 'SCF_TYPE'), ref_wfn)

    if core.get_option('SCF', 'REFERENCE') == 'ROHF':
        ref_wfn.semicanonicalize()

    occ_wfn = core.occ(ref_wfn)

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
    core.set_local_option('OCC', 'DO_SCS', 'FALSE')
    core.set_local_option('OCC', 'DO_SOS', 'FALSE')

    # Bypass the scf call if a reference wavefunction is given
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, **kwargs)  # C1 certified

    # Ensure IWL files have been written
    proc_util.check_iwl_file_from_scf_type(core.get_option('SCF', 'SCF_TYPE'), ref_wfn)

    if core.get_option('SCF', 'REFERENCE') == 'ROHF':
        ref_wfn.semicanonicalize()

    occ_wfn = core.occ(ref_wfn)

    derivobj = core.Deriv(occ_wfn)
    grad = derivobj.compute()

    occ_wfn.set_gradient(grad)

    optstash.restore()
    return occ_wfn


def run_scf(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a self-consistent-field theory (HF & DFT) calculation.

    """

    optstash = proc_util.scf_set_reference_local(name)

    scf_wfn = scf_helper(name, **kwargs)

    optstash.restore()
    return scf_wfn


def run_scf_gradient(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a SCF gradient calculation.

    """
    optstash = proc_util.scf_set_reference_local(name)

    # Bypass the scf call if a reference wavefunction is given
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = run_scf(name, **kwargs)

    if core.get_option('SCF', 'REFERENCE') in ['ROHF', 'CUHF']:
        ref_wfn.semicanonicalize()

    if "_disp_functor" in dir(ref_wfn):
        disp_grad = ref_wfn._disp_functor.compute_gradient(ref_wfn.molecule())
        ref_wfn.set_array("-D Gradient", disp_grad)

    grad = core.scfgrad(ref_wfn)
    ref_wfn.set_gradient(grad)

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

    badref = core.get_option('SCF', 'REFERENCE') in ['UHF', 'ROHF', 'CUHF', 'RKS', 'UKS']
    badint = core.get_option('SCF', 'SCF_TYPE') in [ 'CD', 'OUT_OF_CORE']
    if badref or badint:
        raise ValidationError("Only RHF Hessians are currently implemented. SCF_TYPE either CD or OUT_OF_CORE not supported")
    H = core.scfhess(ref_wfn)
    ref_wfn.set_hessian(H)

    # Temporary freq code.  To be replaced with proper frequency analysis later...
    import numpy as np

    mol = ref_wfn.molecule()
    natoms = mol.natom()
    masses = np.zeros(natoms)

    for atom in range(natoms):
        masses[atom] = mol.mass(atom)

    m = np.repeat( np.divide(1.0, np.sqrt(masses)), 3)
    mwhess = np.einsum('i,ij,j->ij', m, H, m)

    # Are we linear?
    if mol.get_full_point_group() in [ "C_inf_v", "D_inf_h" ]:
        nexternal = 5
    else:
        nexternal = 6

    fcscale = constants.hartree2J / (constants.bohr2m * constants.bohr2m * constants.amu2kg);
    fc = fcscale * np.linalg.eigvalsh(mwhess)
    # Sort by magnitude of the force constants, to project out rot/vib
    ordering = np.argsort(np.abs(fc))
    projected = fc[ordering][nexternal:]
    freqs = np.sqrt(np.abs(projected))
    freqs *= 1.0 / (2.0 * np.pi * constants.c * 100.0)
    freqs[projected < 0] *= -1
    freqs.sort()

    freqvec = core.Vector.from_array(freqs)
    ref_wfn.set_frequencies(freqvec)
    # End of temporary freq hack.  Remove me later!

    # Write Hessian out.  This probably needs a more permanent home, too.
    # This is a drop-in replacement for the code that lives in findif
    if core.get_option('FINDIF', 'HESSIAN_WRITE'):
        molname = ref_wfn.molecule().name()
        prefix = core.get_writer_file_prefix(molname)
        with open(prefix+".hess", 'w') as fp:
            fp.write("%5d%5d\n" % (natoms, 6*natoms))
            for row in np.reshape(H, (-1, 3)):
                fp.write("%20.10f%20.10f%20.10f\n" % tuple(row))

    optstash.restore()
    return ref_wfn


def run_libfock(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a calculation through libfock, namely RCPHF,
    RCIS, RTDHF, RTDA, and RTDDFT.

    """
    if name == 'cphf':
        core.set_global_option('MODULE', 'RCPHF')
    if name == 'cis':
        core.set_global_option('MODULE', 'RCIS')
    if name == 'tdhf':
        core.set_global_option('MODULE', 'RTDHF')
    if name == 'cpks':
        core.set_global_option('MODULE', 'RCPKS')
    if name == 'tda':
        core.set_global_option('MODULE', 'RTDA')
    if name == 'tddft':
        core.set_global_option('MODULE', 'RTDDFT')

    # Bypass the scf call if a reference wavefunction is given
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, **kwargs)

    libfock_wfn = core.libfock(ref_wfn)
    libfock_wfn.compute_energy()
    return libfock_wfn


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
    core.tstart()
    optstash = p4util.OptionsState(
        ['DF_BASIS_SCF'],
        ['DF_BASIS_MP2'],
        ['SCF_TYPE'])  # yes, this really must be global, not local to SCF

    # Alter default algorithm
    if not core.has_option_changed('SCF', 'SCF_TYPE'):
        core.set_global_option('SCF_TYPE', 'DF')
        core.print_out("""    SCF Algorithm Type (re)set to DF.\n""")

    if core.get_option('SCF', 'SCF_TYPE') != 'DF':
        raise ValidationError('DF-MP2 gradients need DF-SCF reference.')

    # Bypass the scf call if a reference wavefunction is given
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, **kwargs)  # C1 certified

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

    optstash.restore()
    core.tstop()
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
        wfn.set_basisset("DF_BASIS_CC", aux_basis)

    # Ensure IWL files have been written
    proc_util.check_iwl_file_from_scf_type(core.get_option('SCF', 'SCF_TYPE'), ref_wfn)

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
    proc_util.check_iwl_file_from_scf_type(core.get_option('SCF', 'SCF_TYPE'), ref_wfn)

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
        core.print_out('Brueckner convergence check: %s\n' % bool(core.get_variable('BRUECKNER CONVERGED')))
        if (core.get_variable('BRUECKNER CONVERGED') == True):
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

def run_dft_property(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    DFT calculations. This is a simple alias to :py:func:`~proc.run_scf`
    since DFT properties all handled through oeprop.

    """

    core.tstart()
    optstash = proc_util.dft_set_reference_local(name)

    properties = kwargs.pop('properties')
    proc_util.oeprop_validator(properties)

    scf_wfn = run_scf(name, scf_do_dipole=False, do_timer=False, **kwargs)

    # Run OEProp
    oe = core.OEProp(scf_wfn)
    oe.set_title(name.upper())
    for prop in properties:
        oe.add(prop.upper())
    oe.compute()
    scf_wfn.set_oeprop(oe)

    core.tstop()
    optstash.restore()
    return scf_wfn


def run_scf_property(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    SCF calculations. This is a simple alias to :py:func:`~proc.run_scf`
    since SCF properties all handled through oeprop.

    """

    core.tstart()
    optstash = proc_util.scf_set_reference_local(name)

    properties = kwargs.pop('properties')
    proc_util.oeprop_validator(properties)

    scf_wfn = run_scf(name, scf_do_dipole=False, do_timer=False, **kwargs)

    # Run OEProp
    oe = core.OEProp(scf_wfn)
    oe.set_title(name.upper())
    for prop in properties:
        oe.add(prop.upper())
    oe.compute()
    scf_wfn.set_oeprop(oe)

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
        ['CCEOM', 'E_CONVERGENCE'])

    oneel_properties = ['dipole', 'quadrupole']
    twoel_properties = []
    response_properties = ['polarizability', 'rotation', 'roa', 'roa_tensor']
    excited_properties = ['oscillator_strength', 'rotational_strength']

    one = []
    two = []
    response = []
    excited = []
    invalid = []

    if 'properties' in kwargs:
        properties = kwargs['properties']

        for prop in properties:
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
        core.set_local_option('CCLAMBDA','R_CONVERGENCE',1e-4)
        core.set_local_option('CCEOM','R_CONVERGENCE',1e-4)
        core.set_local_option('CCEOM','E_CONVERGENCE',1e-5)
        core.cceom(ccwfn)
        core.cclambda(ccwfn)
        core.ccdensity(ccwfn)

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
    if not core.has_option_changed('SCF', 'SCF_TYPE'):
        core.set_global_option('SCF_TYPE', 'DF')  # local set insufficient b/c SCF option read in DFMP2
        core.print_out("""    SCF Algorithm Type (re)set to DF.\n""")

    if not core.get_option('SCF', 'SCF_TYPE') == 'DF':
        raise ValidationError('DF-MP2 properties need DF-SCF reference.')

    properties = kwargs.pop('properties')
    proc_util.oeprop_validator(properties)

    # Bypass the scf call if a reference wavefunction is given
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, scf_do_dipole=False, use_c1=True, **kwargs)  # C1 certified

    aux_basis = core.BasisSet.build(ref_wfn.molecule(), "DF_BASIS_MP2",
                                    core.get_option("DFMP2", "DF_BASIS_MP2"),
                                    "RIFIT", core.get_global_option('BASIS'))
    ref_wfn.set_basisset("DF_BASIS_MP2", aux_basis)

    core.print_out('\n')
    p4util.banner('DFMP2')
    core.print_out('\n')

    dfmp2_wfn = core.dfmp2(ref_wfn)
    grad = dfmp2_wfn.compute_gradient()

    if name == 'scs-mp2':
        core.set_variable('CURRENT ENERGY', core.get_variable('SCS-MP2 TOTAL ENERGY'))
        core.set_variable('CURRENT CORRELATION ENERGY', core.get_variable('SCS-MP2 CORRELATION ENERGY'))
    elif name == 'mp2':
        core.set_variable('CURRENT ENERGY', core.get_variable('MP2 TOTAL ENERGY'))
        core.set_variable('CURRENT CORRELATION ENERGY', core.get_variable('MP2 CORRELATION ENERGY'))

    # Run OEProp
    oe = core.OEProp(dfmp2_wfn)
    oe.set_title(name.upper())
    for prop in properties:
        oe.add(prop.upper())
    oe.compute()
    dfmp2_wfn.set_oeprop(oe)

    optstash.restore()
    return dfmp2_wfn


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
    ciwfn.set_oeprop(oe)

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
    proc_util.check_iwl_file_from_scf_type(core.get_option('SCF', 'SCF_TYPE'), ref_wfn)

    return core.adc(ref_wfn)


def run_dft(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a density-functional-theory calculation.

    """
    optstash = p4util.OptionsState(
        ['SCF', 'DFT_FUNCTIONAL'],
        ['SCF', 'REFERENCE'],
        ['SCF', 'SCF_TYPE'],
        ['DF_BASIS_MP2'],
        ['DFMP2', 'MP2_OS_SCALE'],
        ['DFMP2', 'MP2_SS_SCALE'])

    # Alter default algorithm
    if not core.has_option_changed('SCF', 'SCF_TYPE'):
        core.set_local_option('SCF', 'SCF_TYPE', 'DF')

    core.set_local_option('SCF', 'DFT_FUNCTIONAL', name)

    user_ref = core.get_option('SCF', 'REFERENCE')
    if (user_ref == 'RHF'):
        core.set_local_option('SCF', 'REFERENCE', 'RKS')
    elif (user_ref == 'UHF'):
        core.set_local_option('SCF', 'REFERENCE', 'UKS')
    elif (user_ref == 'ROHF'):
        raise ValidationError('ROHF reference for DFT is not available.')
    elif (user_ref == 'CUHF'):
        raise ValidationError('CUHF reference for DFT is not available.')

    scf_wfn = run_scf(name, **kwargs)
    returnvalue = core.get_variable('CURRENT ENERGY')

    for ssuper in dft_functional.superfunctional_list:
        if ssuper.name().lower() == name:
            dfun = ssuper

    if dfun.is_c_hybrid():
        core.tstart()
        aux_basis = core.BasisSet.build(scf_wfn.molecule(), "DF_BASIS_MP2",
                                        core.get_option("DFMP2", "DF_BASIS_MP2"),
                                        "RIFIT", core.get_global_option('BASIS'),
                                        puream=-1)
        scf_wfn.set_basisset("DF_BASIS_MP2", aux_basis)
        if dfun.is_c_scs_hybrid():
            core.set_local_option('DFMP2', 'MP2_OS_SCALE', dfun.c_os_alpha())
            core.set_local_option('DFMP2', 'MP2_SS_SCALE', dfun.c_ss_alpha())
            dfmp2_wfn = core.dfmp2(scf_wfn)
            dfmp2_wfn.compute_energy()

            vdh = dfun.c_alpha() * core.get_variable('SCS-MP2 CORRELATION ENERGY')

        else:
            dfmp2_wfn = core.dfmp2(scf_wfn)
            dfmp2_wfn.compute_energy()
            vdh = dfun.c_alpha() * core.get_variable('MP2 CORRELATION ENERGY')

        # TODO: delete these variables, since they don't mean what they look to mean?
        # 'MP2 TOTAL ENERGY',
        # 'MP2 CORRELATION ENERGY',
        # 'MP2 SAME-SPIN CORRELATION ENERGY']

        core.set_variable('DOUBLE-HYBRID CORRECTION ENERGY', vdh)
        returnvalue += vdh
        core.set_variable('DFT TOTAL ENERGY', returnvalue)
        core.set_variable('CURRENT ENERGY', returnvalue)
        core.print_out('\n\n')
        core.print_out('    %s Energy Summary\n' % (name.upper()))
        core.print_out('    -------------------------\n')
        core.print_out('    DFT Reference Energy                  = %22.16lf\n' % (returnvalue - vdh))
        core.print_out('    Scaled MP2 Correlation                = %22.16lf\n' % (vdh))
        core.print_out('    @Final double-hybrid DFT total energy = %22.16lf\n\n' % (returnvalue))
        core.tstop()

    optstash.restore()
    return scf_wfn


def run_dft_gradient(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a density-functional-theory gradient calculation.

    """
    optstash = p4util.OptionsState(
        ['SCF', 'DFT_FUNCTIONAL'],
        ['SCF', 'REFERENCE'],
        ['SCF', 'SCF_TYPE'])

    # Alter default algorithm
    if not core.has_option_changed('SCF', 'SCF_TYPE'):
        core.set_local_option('SCF', 'SCF_TYPE', 'DF')

    core.set_local_option('SCF', 'DFT_FUNCTIONAL', name.upper())

    user_ref = core.get_option('SCF', 'REFERENCE')
    if (user_ref == 'RHF'):
        core.set_local_option('SCF', 'REFERENCE', 'RKS')
    elif (user_ref == 'UHF'):
        core.set_local_option('SCF', 'REFERENCE', 'UKS')
    elif (user_ref == 'ROHF'):
        raise ValidationError('ROHF reference for DFT is not available.')
    elif (user_ref == 'CUHF'):
        raise ValidationError('CUHF reference for DFT is not available.')

    wfn = run_scf_gradient(name, **kwargs)

    optstash.restore()
    return wfn


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
    proc_util.check_iwl_file_from_scf_type(core.get_option('SCF', 'SCF_TYPE'), ref_wfn)

    ciwfn = core.detci(ref_wfn)

    print_nos = False
    if core.get_option("DETCI", "NAT_ORBS"):
        ciwfn.ci_nat_orbs()
        print_nos = True

    proc_util.print_ci_results(ciwfn, name.upper(), core.get_variable("HF TOTAL ENERGY"), core.get_variable("CURRENT ENERGY"), print_nos)

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
        ciwfn.set_oeprop(oeprop)
        core.set_variable("CURRENT DIPOLE X", core.get_variable(name.upper() + " DIPOLE X"))
        core.set_variable("CURRENT DIPOLE Y", core.get_variable(name.upper() + " DIPOLE Y"))
        core.set_variable("CURRENT DIPOLE Z", core.get_variable(name.upper() + " DIPOLE Z"))

    ciwfn.cleanup_ci()
    ciwfn.cleanup_dpd()

    optstash.restore()
    return ciwfn


def run_dfmp2(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a density-fitted MP2 calculation.

    """
    optstash = p4util.OptionsState(
        ['DF_BASIS_MP2'],
        ['SCF', 'SCF_TYPE'])

    # Alter default algorithm
    if not core.has_option_changed('SCF', 'SCF_TYPE'):
        core.set_local_option('SCF', 'SCF_TYPE', 'DF')
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
        core.set_variable('CURRENT ENERGY', core.get_variable('SCS-MP2 TOTAL ENERGY'))
        core.set_variable('CURRENT CORRELATION ENERGY', core.get_variable('SCS-MP2 CORRELATION ENERGY'))
    elif name == 'mp2':
        core.set_variable('CURRENT ENERGY', core.get_variable('MP2 TOTAL ENERGY'))
        core.set_variable('CURRENT CORRELATION ENERGY', core.get_variable('MP2 CORRELATION ENERGY'))

    optstash.restore()
    core.tstop()
    return dfmp2_wfn


def run_dmrgscf(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    an DMRG calculation.

    """
    optstash = p4util.OptionsState(
        ['SCF', 'SCF_TYPE'],
        ['DMRG', 'DMRG_CASPT2_CALC'])

    # Bypass the scf call if a reference wavefunction is given
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, **kwargs)

    # Ensure IWL files have been written
    proc_util.check_iwl_file_from_scf_type(core.get_option('SCF', 'SCF_TYPE'), ref_wfn)

    if 'CASPT2' in name.upper():
        core.set_local_option("DMRG", "DMRG_CASPT2_CALC", True)

    dmrg_wfn = core.dmrg(ref_wfn)
    optstash.restore()

    return dmrg_wfn


def run_dmrgci(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    an DMRG calculation.

    """
    optstash = p4util.OptionsState(
        ['SCF', 'SCF_TYPE'],
        ['DMRG', 'DMRG_SCF_MAX_ITER'])

    # Bypass the scf call if a reference wavefunction is given
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, **kwargs)

    # Ensure IWL files have been written
    proc_util.check_iwl_file_from_scf_type(core.get_option('SCF', 'SCF_TYPE'), ref_wfn)

    core.set_local_option('DMRG', 'DMRG_SCF_MAX_ITER', 1)

    dmrg_wfn = core.dmrg(ref_wfn)
    optstash.restore()

    return dmrg_wfn


def run_psimrcc(name, **kwargs):
    """Function encoding sequence of PSI module calls for a PSIMRCC computation
     using a reference from the MCSCF module

    """
    mcscf_wfn = run_mcscf(name, **kwargs)
    psimrcc_e = core.psimrcc(mcscf_wfn)

    return mcscf_wfn


def run_psimrcc_scf(name, **kwargs):
    """Function encoding sequence of PSI module calls for a PSIMRCC computation
     using a reference from the SCF module

    """
    # Bypass the scf call if a reference wavefunction is given
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, **kwargs)

    psimrcc_e = core.psimrcc(ref_wfn)

    return ref_wfn


def run_sapt(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a SAPT calculation of any level.

    """
    optstash = p4util.OptionsState(
        ['SCF', 'SCF_TYPE'])

    # Alter default algorithm
    if not core.has_option_changed('SCF', 'SCF_TYPE'):
        core.set_local_option('SCF', 'SCF_TYPE', 'DF')

    # Get the molecule of interest
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        sapt_dimer = kwargs.pop('molecule', core.get_active_molecule())
    else:
        core.print_out('Warning! SAPT argument "ref_wfn" is only able to use molecule information.')
        sapt_dimer = ref_wfn.molecule()
    sapt_dimer.update_geometry()  # make sure since mol from wfn, kwarg, or P::e
    sapt_dimer.fix_orientation(True)
    sapt_dimer.fix_com(True)

    # Shifting to C1 so we need to copy the active molecule
    if sapt_dimer.schoenflies_symbol() != 'c1':
        core.print_out('  SAPT does not make use of molecular symmetry, further calculations in C1 point group.\n')
        sapt_dimer = sapt_dimer.clone()
        sapt_dimer.reset_point_group('c1')
        sapt_dimer.fix_orientation(True)
        sapt_dimer.fix_com(True)
        sapt_dimer.update_geometry()

    if (core.get_option('SCF', 'REFERENCE') != 'RHF') and (name.upper() != "SAPT0"):
        raise ValidationError('Only SAPT0 supports a reference different from \"reference rhf\".')

    nfrag = sapt_dimer.nfragments()
    if nfrag != 2:
        raise ValidationError('SAPT requires active molecule to have 2 fragments, not %s.' % (nfrag))

    do_delta_mp2 = True if name.endswith('dmp2') else False

    sapt_basis = 'dimer'
    if 'sapt_basis' in kwargs:
        sapt_basis = kwargs.pop('sapt_basis')
    sapt_basis = sapt_basis.lower()

    if sapt_basis == 'dimer':
        monomerA = sapt_dimer.extract_subsets(1, 2)
        monomerA.set_name('monomerA')
        monomerB = sapt_dimer.extract_subsets(2, 1)
        monomerB.set_name('monomerB')
    elif sapt_basis == 'monomer':
        monomerA = sapt_dimer.extract_subsets(1)
        monomerA.set_name('monomerA')
        monomerB = sapt_dimer.extract_subsets(2)
        monomerB.set_name('monomerB')

    ri = core.get_option('SCF', 'SCF_TYPE')
    df_ints_io = core.get_option('SCF', 'DF_INTS_IO')
    # inquire if above at all applies to dfmp2

    core.IO.set_default_namespace('dimer')
    core.print_out('\n')
    p4util.banner('Dimer HF')
    core.print_out('\n')

    # Compute dimer wavefunction
    if (sapt_basis == 'dimer') and (ri == 'DF'):
        core.set_global_option('DF_INTS_IO', 'SAVE')

    dimer_wfn = scf_helper('RHF', molecule=sapt_dimer, **kwargs)
    if do_delta_mp2:
        select_mp2(name, ref_wfn=dimer_wfn, **kwargs)
        mp2_corl_interaction_e = core.get_variable('MP2 CORRELATION ENERGY')

    if (sapt_basis == 'dimer') and (ri == 'DF'):
        core.set_global_option('DF_INTS_IO', 'LOAD')

    # Compute Monomer A wavefunction
    if (sapt_basis == 'dimer') and (ri == 'DF'):
        core.IO.change_file_namespace(97, 'dimer', 'monomerA')

    core.IO.set_default_namespace('monomerA')
    core.print_out('\n')
    p4util.banner('Monomer A HF')
    core.print_out('\n')
    monomerA_wfn = scf_helper('RHF', molecule=monomerA, **kwargs)
    if do_delta_mp2:
        select_mp2(name, ref_wfn=monomerA_wfn, **kwargs)
        mp2_corl_interaction_e -= core.get_variable('MP2 CORRELATION ENERGY')

    # Compute Monomer B wavefunction
    if (sapt_basis == 'dimer') and (ri == 'DF'):
        core.IO.change_file_namespace(97, 'monomerA', 'monomerB')
    core.IO.set_default_namespace('monomerB')
    core.print_out('\n')
    p4util.banner('Monomer B HF')
    core.print_out('\n')
    monomerB_wfn = scf_helper('RHF', molecule=monomerB, **kwargs)

    # Delta MP2
    if do_delta_mp2:
        select_mp2(name, ref_wfn=monomerB_wfn, **kwargs)
        mp2_corl_interaction_e -= core.get_variable('MP2 CORRELATION ENERGY')
        core.set_variable('SAPT MP2 CORRELATION ENERGY', mp2_corl_interaction_e)
    core.set_global_option('DF_INTS_IO', df_ints_io)

    if core.get_option('SCF', 'REFERENCE') == 'RHF':
        core.IO.change_file_namespace(constants.PSIF_SAPT_MONOMERA, 'monomerA', 'dimer')
        core.IO.change_file_namespace(constants.PSIF_SAPT_MONOMERB, 'monomerB', 'dimer')

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
        core.set_local_option('SAPT','COUPLED_INDUCTION',False)
        core.print_out('  Coupled induction not available for ROHF.\n')
        core.print_out('  Proceeding with uncoupled induction only.\n')

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
    p4util.banner(name.upper())
    core.print_out('\n')
    e_sapt = core.sapt(dimer_wfn, monomerA_wfn, monomerB_wfn)

    from psi4.driver.qcdb.psivardefs import sapt_psivars
    p4util.expand_psivars(sapt_psivars())
    optstash.restore()
    for term in ['ELST', 'EXCH', 'IND', 'DISP', 'TOTAL']:
        core.set_variable(' '.join(['SAPT', term, 'ENERGY']),
            core.get_variable(' '.join([name.upper(), term, 'ENERGY'])))
    core.set_variable('CURRENT ENERGY', core.get_variable('SAPT TOTAL ENERGY'))

    return dimer_wfn


def run_sapt_ct(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a charge-transfer SAPT calcuation of any level.

    """
    optstash = p4util.OptionsState(
        ['SCF', 'SCF_TYPE'])

    if 'ref_wfn' in kwargs:
        core.print_out('\nWarning! Argument ref_wfn is not valid for sapt computations\n')

    # Alter default algorithm
    if not core.has_option_changed('SCF', 'SCF_TYPE'):
        core.set_local_option('SCF', 'SCF_TYPE', 'DF')

    # Get the molecule of interest
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        sapt_dimer = kwargs.pop('molecule', core.get_active_molecule())
    else:
        core.print_out('Warning! SAPT argument "ref_wfn" is only able to use molecule information.')
        sapt_dimer = ref_wfn.molecule()
    sapt_dimer.update_geometry()  # make sure since mol from wfn, kwarg, or P::e
    sapt_dimer.fix_orientation(True)
    sapt_dimer.fix_com(True)

    # Shifting to C1 so we need to copy the active molecule
    if sapt_dimer.schoenflies_symbol() != 'c1':
        core.print_out('  SAPT does not make use of molecular symmetry, further calculations in C1 point group.\n')
        sapt_dimer = sapt_dimer.clone()
        sapt_dimer.reset_point_group('c1')
        sapt_dimer.fix_orientation(True)
        sapt_dimer.fix_com(True)
        sapt_dimer.update_geometry()

    if core.get_option('SCF', 'REFERENCE') != 'RHF':
        raise ValidationError('SAPT requires requires \"reference rhf\".')

    nfrag = sapt_dimer.nfragments()
    if nfrag != 2:
        raise ValidationError('SAPT requires active molecule to have 2 fragments, not %s.' % (nfrag))

    monomerA = sapt_dimer.extract_subsets(1, 2)
    monomerA.set_name('monomerA')
    monomerB = sapt_dimer.extract_subsets(2, 1)
    monomerB.set_name('monomerB')
    sapt_dimer.update_geometry()
    monomerAm = sapt_dimer.extract_subsets(1)
    monomerAm.set_name('monomerAm')
    monomerBm = sapt_dimer.extract_subsets(2)
    monomerBm.set_name('monomerBm')

    ri = core.get_option('SCF', 'SCF_TYPE')
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
    core.IO.change_file_namespace(constants.PSIF_SAPT_MONOMERA, 'monomerA', 'dimer')
    core.IO.change_file_namespace(constants.PSIF_SAPT_MONOMERB, 'monomerB', 'dimer')
    e_sapt = core.sapt(dimer_wfn, monomerA_wfn, monomerB_wfn)
    CTd = core.get_variable('SAPT CT ENERGY')

    core.print_out('\n')
    p4util.banner('Monomer Basis SAPT')
    core.print_out('\n')
    core.IO.change_file_namespace(constants.PSIF_SAPT_MONOMERA, 'monomerAm', 'dimer')
    core.IO.change_file_namespace(constants.PSIF_SAPT_MONOMERB, 'monomerBm', 'dimer')
    e_sapt = core.sapt(dimer_wfn, monomerAm_wfn, monomerBm_wfn)
    CTm = core.get_variable('SAPT CT ENERGY')
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
    core.set_variable('SAPT CT ENERGY', CT)

    optstash.restore()
    return dimer_wfn


def run_fisapt(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    an F/ISAPT0 computation

    """
    optstash = p4util.OptionsState(
        ['SCF', 'SCF_TYPE'])

    # Alter default algorithm
    if not core.has_option_changed('SCF', 'SCF_TYPE'):
        core.set_local_option('SCF', 'SCF_TYPE', 'DF')

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
        ref_wfn = scf_helper('RHF', molecule=sapt_dimer, **kwargs)

    scf_aux_basis = core.BasisSet.build(ref_wfn.molecule(), "DF_BASIS_SCF",
                                        core.get_option("SCF", "DF_BASIS_SCF"),
                                        "JKFIT", core.get_global_option('BASIS'),
                                        puream=ref_wfn.basisset().has_puream())
    ref_wfn.set_basisset("DF_BASIS_SCF", scf_aux_basis)

    sapt_basis = core.BasisSet.build(ref_wfn.molecule(), "DF_BASIS_SAPT",
                                     core.get_global_option("DF_BASIS_SAPT"),
                                     "RIFIT", core.get_global_option("BASIS"),
                                     ref_wfn.basisset().has_puream())
    ref_wfn.set_basisset("DF_BASIS_SAPT", sapt_basis)

    minao = core.BasisSet.build(ref_wfn.molecule(), "BASIS",
                                core.get_global_option("MINAO_BASIS"))
    ref_wfn.set_basisset("MINAO", minao)


    fisapt_wfn = core.fisapt(ref_wfn)

    optstash.restore()
    return fisapt_wfn


def run_mrcc(name, **kwargs):
    """Function that prepares environment and input files
    for a calculation calling Kallay's MRCC code.

    """

    # Check to see if we really need to run the SCF code.
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, **kwargs)
    vscf = core.get_variable('SCF TOTAL ENERGY')

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
        ['SCF', 'DF_INTS_IO'],
        ['SCF', 'SCF_TYPE'])

    core.set_local_option('FNOCC', 'DFCC', True)
    core.set_local_option('FNOCC', 'RUN_CEPA', False)

    # throw an exception for open-shells
    if core.get_option('SCF', 'REFERENCE') != 'RHF':
        raise ValidationError("""Error: %s requires 'reference rhf'.""" % name)

    def set_cholesky_from(mtd_type):
        type_val = core.get_global_option(mtd_type)
        if type_val == 'CD':
            core.set_local_option('FNOCC', 'DF_BASIS_CC', 'CHOLESKY')
            # Alter default algorithm
            if not core.has_option_changed('SCF', 'SCF_TYPE'):
                core.set_global_option('SCF_TYPE', 'CD')
                core.print_out("""    SCF Algorithm Type (re)set to CD.\n""")
        elif type_val == 'DF':
            if core.get_option('FNOCC', 'DF_BASIS_CC') == 'CHOLESKY':
                core.set_local_option('FNOCC', 'DF_BASIS_CC', '')
            # Alter default algorithm
            if not core.has_option_changed('SCF', 'SCF_TYPE'):
                core.set_global_option('SCF_TYPE', 'DF')
                core.print_out("""    SCF Algorithm Type (re)set to DF.\n""")
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

    if core.get_option('SCF', 'SCF_TYPE') not in ['CD', 'DF']:
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

    scf_aux_basis = core.BasisSet.build(ref_wfn.molecule(), "DF_BASIS_SCF",
                                        core.get_option("SCF", "DF_BASIS_SCF"),
                                        "JKFIT", core.get_global_option('BASIS'),
                                        puream=ref_wfn.basisset().has_puream())
    ref_wfn.set_basisset("DF_BASIS_SCF", scf_aux_basis)

    aux_basis = core.BasisSet.build(ref_wfn.molecule(), "DF_BASIS_CC",
                                        core.get_global_option("DF_BASIS_CC"),
                                        "RIFIT", core.get_global_option("BASIS"))
    ref_wfn.set_basisset("DF_BASIS_CC", aux_basis)

    fnocc_wfn = core.fnocc(ref_wfn)

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
        raise ValidationError("""Error: %s requires 'reference rhf'.""" % name)

    # Bypass the scf call if a reference wavefunction is given
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, **kwargs)  # C1 certified

    if core.get_option('FNOCC', 'USE_DF_INTS') == False:
        # Ensure IWL files have been written
        proc_util.check_iwl_file_from_scf_type(core.get_option('SCF', 'SCF_TYPE'), ref_wfn)
    else:
        scf_aux_basis = core.BasisSet.build(ref_wfn.molecule(), "DF_BASIS_SCF",
                                            core.get_option("SCF", "DF_BASIS_SCF"),
                                            "JKFIT", core.get_global_option('BASIS'),
                                            puream=ref_wfn.basisset().has_puream())
        ref_wfn.set_basisset("DF_BASIS_SCF", scf_aux_basis)

    fnocc_wfn = core.fnocc(ref_wfn)

    # set current correlation energy and total energy.  only need to treat mpn here.
    if name == 'mp3':
        emp3 = core.get_variable("MP3 TOTAL ENERGY")
        cemp3 = core.get_variable("MP3 CORRELATION ENERGY")
        core.set_variable("CURRENT ENERGY", emp3)
        core.set_variable("CURRENT CORRELATION ENERGY", cemp3)
    elif name == 'fno-mp3':
        emp3 = core.get_variable("MP3 TOTAL ENERGY")
        cemp3 = core.get_variable("MP3 CORRELATION ENERGY")
        core.set_variable("CURRENT ENERGY", emp3)
        core.set_variable("CURRENT CORRELATION ENERGY", cemp3)
    elif name == 'mp4(sdq)':
        emp4sdq = core.get_variable("MP4(SDQ) TOTAL ENERGY")
        cemp4sdq = core.get_variable("MP4(SDQ) CORRELATION ENERGY")
        core.set_variable("CURRENT ENERGY", emp4sdq)
        core.set_variable("CURRENT CORRELATION ENERGY", cemp4sdq)
    elif name == 'fno-mp4(sdq)':
        emp4sdq = core.get_variable("MP4(SDQ) TOTAL ENERGY")
        cemp4sdq = core.get_variable("MP4(SDQ) CORRELATION ENERGY")
        core.set_variable("CURRENT ENERGY", emp4sdq)
        core.set_variable("CURRENT CORRELATION ENERGY", cemp4sdq)
    elif name == 'fno-mp4':
        emp4 = core.get_variable("MP4 TOTAL ENERGY")
        cemp4 = core.get_variable("MP4 CORRELATION ENERGY")
        core.set_variable("CURRENT ENERGY", emp4)
        core.set_variable("CURRENT CORRELATION ENERGY", cemp4)
    elif name == 'mp4':
        emp4 = core.get_variable("MP4 TOTAL ENERGY")
        cemp4 = core.get_variable("MP4 CORRELATION ENERGY")
        core.set_variable("CURRENT ENERGY", emp4)
        core.set_variable("CURRENT CORRELATION ENERGY", cemp4)

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

    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, **kwargs)  # C1 certified

    if core.get_option('FNOCC', 'USE_DF_INTS') == False:
        # Ensure IWL files have been written
        proc_util.check_iwl_file_from_scf_type(core.get_option('SCF', 'SCF_TYPE'), ref_wfn)
    else:
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

    optstash.restore()
    return fnocc_wfn


def run_detcas(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    determinant-based multireference wavefuncations,
    namely CASSCF and RASSCF.
    """

    optstash = p4util.OptionsState(
        ['DETCI', 'WFN'],
        ['SCF', 'SCF_TYPE']
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
        if not core.has_option_changed('SCF', 'SCF_TYPE'):
            core.set_global_option('SCF_TYPE', 'DF')

        # If RHF get MP2 NO's
        # Why doesnt this work for conv?
        if ((core.get_option('SCF', 'SCF_TYPE') == 'DF') and (user_ref == 'RHF') and
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
        if not core.has_option_changed('SCF', 'SCF_TYPE'):
            core.set_global_option('SCF_TYPE', 'DF')

        scf_aux_basis = core.BasisSet.build(ref_wfn.molecule(), "DF_BASIS_SCF",
                                            core.get_option("SCF", "DF_BASIS_SCF"),
                                            "JKFIT", core.get_global_option('BASIS'),
                                            puream=ref_wfn.basisset().has_puream())
        ref_wfn.set_basisset("DF_BASIS_SCF", scf_aux_basis)

    # The AO case
    elif core.get_option('DETCI', 'MCSCF_TYPE') == 'AO':
        if not core.has_option_changed('SCF', 'SCF_TYPE'):
            core.set_global_option('SCF_TYPE', 'DIRECT')

    # The conventional case
    elif core.get_option('DETCI', 'MCSCF_TYPE') == 'CONV':
        if not core.has_option_changed('SCF', 'SCF_TYPE'):
            core.set_global_option('SCF_TYPE', 'PK')

        # Ensure IWL files have been written
        proc_util.check_iwl_file_from_scf_type(core.get_option('SCF', 'SCF_TYPE'), ref_wfn)
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
    ciwfn.set_oeprop(oeprop)
    core.set_variable("CURRENT DIPOLE X", core.get_variable(name.upper() + " DIPOLE X"))
    core.set_variable("CURRENT DIPOLE Y", core.get_variable(name.upper() + " DIPOLE Y"))
    core.set_variable("CURRENT DIPOLE Z", core.get_variable(name.upper() + " DIPOLE Z"))

    optstash.restore()
    return ciwfn

def run_efp(name, **kwargs):
    """Function encoding sequence of module calls for a pure EFP
    computation (ignore any QM atoms).

    """
    # initialize library
    efp = core.get_active_efp()

    if efp.nfragments() == 0:
        raise ValidationError("""Method 'efp' not available without EFP fragments in molecule""")

    # set options
    core.set_global_option('QMEFP', False)  # apt to go haywire if set locally to efp
    core.efp_set_options()

    efp.print_out()
    returnvalue = efp.compute()
    return returnvalue
