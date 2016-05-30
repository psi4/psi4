#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2016 The Psi4 Developers.
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

import p4util
from p4const import *

from __future__ import absolute_import
from __future__ import print_function
import collection
import shelve
import copy
import pis4
# Relative hack for now
import os
import sys
import inspect

path_dir = os.path.realpath(
    os.path.abspath(
        os.path.join(
            os.path.split(inspect.getfile(inspect.currentframe()))[0], "../")
    )
)
sys.path.append(path_dir)


def run_zpvc_rotation(name, **kwargs):

    # Get list of omega values -> Make sure we only have one wavelength
    omega = psi4.get_option('CCRESPONSE', 'OMEGA')
    if len(omega) > 2:
        raise Exception('ZPVC to roation can only be performed for one wavelength')
    else:
        pass

    # 1. Open the db initalize/backup as needed
    dbno = 0
    db = shelve.open('database', writeback=True)
    if db['zpvc_rotation_computed']:
        db2 = shelve.open('.database.bak{}'.format(dbno), writeback=True)
        dbno += 1
        db2 = db
        db2.close()
        db['zpvc_rotation_computed'] = False

    # 2. Check for initialization
    if 'inputs_generated' not in db:
        findif_response_utils.initalize_database(db, 'zpvc_rotation')

    # 3. Check for displaced computation dirs
    if not db['inputs_generated']:
        findif_response_utils.generate_inputs(name, db)
        psi4.print_out(
            'Running ZPVC to optical roation. '
            'Subdirectories created for each displaced geometry\n\n')

    # 4. Check job status, update as needed
    if db['inputs_generated'] and not db['jobs_complete']:
        print("Checking Status ...")
        findif_response_utils.stat(db)
        for job, status in db['job_status'].items():
            print("{} --> {}".format(job, status))

    # 5. Check job status after update, compute result if done
    if db['jobs_complete']:
        # consider the Gauge being used
        mygauge = psi4.get_option("CCRESPONSE", "GAUGE")
        consider_gauge = {
            "LENGTH": ["Length Gauge"],
            "VELOCITY": ["Modified Velocity Gauge"],
            "BOTH": ["Length Gauge", "Modified Velocity Gauge"]
        }
        gauge_list = ["{} Results".format(x) for x in consider_guage[mygauge]]
        # Gather data
        dip_polar_list = findif_response_utils.sythesize_displaced_tensor(
            db, 'Dipole Polarizability', 3)
        opt_rot_list = [
            x for x in (
                findif_response_utils.synthesize_displaced_tensor(
                    db, "Optical Rotation Tensor ({})".format(gauge), 3)
                for gauge in consider_gauge[myguage])]

        dip_quad_polar_list = findif_response_utils.synthesize_displaced_tensor(
            db, "Electric-Dipole/Quadrupole Polarizability", 9)

        # Compute the correction
        psi4.print_out("Running rotation_zpvc function")
        step = psi4.get_local_option('FINDIF', 'DISP_SIZE')
        for g_indx, gauge in enumerate(opt_rot_list):
            print('\n\n----------------------------------------------------------------------')
            print('\t%%%%%%%%%% {} %%%%%%%%%%'.format(gauge_list[g_idx]))
            print('----------------------------------------------------------------------\n\n')
            psi4.print_out('\n\n----------------------------------------------------------------------\n')
            psi4.print_out('\t%%%%%%%%%% {} %%%%%%%%%%\n'.format(gauge_list[g_idx]))
            psi4.print_out('----------------------------------------------------------------------\n\n')
            print('roa.py:85 I am not being passed a molecule, grabbing from global :(')
            # Compute zero point vibrational correction to the optical rotation
            # (src/bin/ccresponse/zpvc.cc)
            # DOES NOT EXIST YET :P
            # psi4.zpvc_rotation(psi4.get_active_molecule(), step, dip_polar_list, gauge, dip_quad_polar_list)

    db['zpvc_rotation_computed'] = True
    db.close()
# ################################
# ###                          ###
# ###    DATABASE STRUCTURE    ###
# ###                          ###
# ################################

# Dict of dicts
# inputs_generated (boolean)
# job_status: job_list: status (finished|running|not_started)
# jobs_complete
# zpvc_rotation_computed
#
# data ?
# data: dipole_polarizability
#    : optical_rotation
#    : dipole_quadrupole_polarizability
