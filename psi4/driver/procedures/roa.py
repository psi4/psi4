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
from __future__ import absolute_import
from __future__ import print_function

import collections
import shelve
import copy
import sys
import inspect
import os

from psi4.driver.constants import *
from psi4.driver.p4util import *
from psi4 import core
from . import findif_response_utils


def run_roa(name, **kwargs):
    """
        Main driver for managing Raman Optical activity computations with
        CC response theory.

        Uses distributed finite differences approach -->
            1. Sets up a database to keep track of running/finished/waiting
                computations.
            2. Generates separate input files for displaced geometries.
            3. When all displacements are run, collects the necessary information
                from each displaced computation, and computes final result.
    """

    # Get list of omega values -> Make sure we only have one wavelength
    # Catch this now before any real work gets done
    omega = core.get_option('CCRESPONSE', 'OMEGA')
    if len(omega) > 2:
        raise Exception('ROA scattering can only be performed for one wavelength.')
    else:
        pass

    core.print_out(
        'Running ROA computation. Subdirectories for each '
        'required displaced geometry have been created.\n\n')

    dbno = 0
    # Initialize database
    db = shelve.open('database', writeback=True)
    # Check if final result is in here
    # ->if we have already computed roa, back up the dict
    # ->copy it setting this flag to false and continue
    if ('roa_computed' in db) and ( db['roa_computed'] ):
        db2 = shelve.open('.database.bak{}'.format(dbno), writeback=True)
        dbno += 1
        for key,value in db.items():
            db2[key]=value

        db2.close()
        db['roa_computed'] = False
    else:
        db['roa_computed'] = False

    if 'inputs_generated' not in db:
        findif_response_utils.initialize_database(db,name,"roa", ["roa_tensor"])

    # Generate input files
    if not db['inputs_generated']:
        findif_response_utils.generate_inputs(db,name)
        # handled by helper db['inputs_generated'] = True

    # Check job status
    if db['inputs_generated'] and not db['jobs_complete']:
        print('Checking status')
        findif_response_utils.stat(db)
        for job, status in db['job_status'].items():
            print("{} --> {}".format(job, status))

    # Compute ROA Scattering
    if db['jobs_complete']:
        mygauge = core.get_option('CCRESPONSE', 'GAUGE')
        consider_gauge = {
            'LENGTH': ['Length Gauge'],
            'VELOCITY': ['Modified Velocity Gauge'],
            'BOTH': ['Length Gauge', 'Modified Velocity Gauge']
        }
        gauge_list = ["{} Results".format(x) for x in consider_gauge[mygauge]]
        # Gather data
        dip_polar_list = findif_response_utils.collect_displaced_matrix_data(
            db, 'Dipole Polarizability', 3)
        opt_rot_list = [
            x for x in (
                findif_response_utils.collect_displaced_matrix_data(
                    db,
                    "Optical Rotation Tensor ({})".format(gauge),
                    3
                )
                for gauge in consider_gauge[mygauge]
            )
        ]
        dip_quad_polar_list = findif_response_utils.collect_displaced_matrix_data(
            db, "Electric-Dipole/Quadrupole Polarizability", 9)
        # Compute Scattering
        # Run new function (src/bin/ccresponse/scatter.cc)
        core.print_out('Running scatter function')
        step = core.get_local_option('FINDIF', 'DISP_SIZE')
        for g_idx, gauge in enumerate(opt_rot_list):
            print('\n\n----------------------------------------------------------------------')
            print('\t%%%%%%%%%% {} %%%%%%%%%%'.format(gauge_list[g_idx]))
            print('----------------------------------------------------------------------\n\n')
            core.print_out('\n\n----------------------------------------------------------------------\n')
            core.print_out('\t%%%%%%%%%% {} %%%%%%%%%%\n'.format(gauge_list[g_idx]))
            core.print_out('----------------------------------------------------------------------\n\n')
            print('roa.py:85 I am not being passed a molecule, grabbing from global :(')
            core.scatter(core.get_active_molecule(), step, dip_polar_list, gauge, dip_quad_polar_list)

        db['roa_computed'] = True

    db.close()

#   SAVE this for when multiple wavelengths works
# # Get list of omega values
# omega = core.get_option('CCRESPONSE','OMEGA')
# if len(omega) > 1:
#     units = copy.copy(omega[-1])
#     omega.pop()
# else:
#     units = 'atomic'
# wavelength = copy.copy(omega[0])
# # Set up units for scatter.cc
# if units == 'NM':
#     wavelength = (constants.c * constants.h * 1*(10**-9))/(wavelength * constants.hartree2J)
# if units == 'HZ':
#     wavelength = wavelength * constants.h / constants.hartree2J
# if units == 'EV':
#     wavelength = wavelength / constants.hartree2ev
# if units == 'atomic':
#     pass
# ################################
# ###                          ###
# ###    DATABASE STRUCTURE    ###
# ###                          ###
# ################################

# Dict of dicts
# inputs_generated (boolean)
# job_status: (ordered Dict)
#    key-> {atom}_{cord}_{p/m}
#       val-> (not_started,running,finished)
#    job_list: (string)
#        status (string)
# jobs_complete (boolean)
# roa_computed (boolean)
# prop (string) = roa
#

# ?
# data: dipole_polarizability
#    : optical_rotation
#    : dipole_quadrupole_polarizability
# ?
# results:
