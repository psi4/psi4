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

"""
Module of helper functions for ccresponse distributed property calculations.
Defines functions for interacting with the database created by the run_XXX
driver function.

Properties that are able to use this module should be added to
the registered_props dictionary.

"""
from __future__ import absolute_import
from __future__ import print_function
import collections
import shelve
import copy
import os
import psi4
import p4util
from p4const import *

"""
registered_props (dict)

The dictionary key is the top-level input properties array argument that indicates
the driver should be used. The value is the properties array argument that
should be used in each subdir computation.  Drivers using this module will be
added as needed.
"""
registered_props = {
    "roa": "roa_tensor"
}


def generate_inputs(name, db):
    """
        Generates the input files in each sub-directory of the
        distributed finite differences property calculation.

    name: ( string ) method name passed to calling driver,
    db:   (database) The database object associated with this property
          calculation. On exit this db['inputs_generated'] has been set True

    Returns: nothing
    Throws: Exception if the number of atomic displacements is not correct.
    """
    molecule = psi4.get_active_molecule()
    natom = molecule.natom()

    # get list of displacements
    displacement_geoms = psi4.atomic_displacements(molecule)

    # Sanity Check
    # there should be 3 cords * natoms *2 directions (+/-)
    if not (6 * natom) == len(displacement_geoms):
        raise Exception('The number of atomic displacements should be 6 times'
                        ' the number of atoms!')

    displacement_names = db['job_status'].keys()

    for n, entry in enumerate(displacement_names):
        if not os.path.exists(entry):
            os.makedirs(entry)

        # Setup up input file string
        inp_template = 'molecule {molname}_{disp}'
        inp_template += ' {{\n{molecule_info}\n}}\n{options}\n{jobspec}\n'
        molecule.set_geometry(displacement_geoms[n])
        molecule.fix_orientation(True)
        molecule.fix_com(True)
        inputfile = open('{0}/input.dat'.format(entry), 'w')
        inputfile.write("# This is a psi4 input file auto-generated for"
            "computing properties by finite differences.\n\n")
        inputfile.write(
            inp_template.format(
                molname=molecule.name(),
                disp=entry,
                molecule_info=molecule.create_psi4_string_from_molecule(),
                options=p4util.format_options_for_input(),
                jobspec="property('{0}', properties=['{1}'])".format(
                    name, db['prop_string'])))
        inputfile.close()
    db['inputs_generated'] = True

    # END generate_inputs


def initialize_database(database, prop):
    """
        Initialize the database for computation of some property
        using distributed finite differences driver

    database: (database) the database object passed from the caller
    prop:  (string) the property that is being computed

    Returns: nothing
    Throws: nothing
    """
    database['inputs_generated'] = False
    database['jobs_complete'] = False
    database['{}_computed'.format(prop)] = False
    database['prop_string'] = registered_props[prop]
    database['job_status'] = collections.OrderedDict()
    # Populate the job_status dict
    molecule = psi4.get_active_molecule()
    natom = molecule.natom()
    coordinates = ['x', 'y', 'z']
    step_direction = ['p', 'm']

    for atom in range(1, natom + 1):
        for coord in coordinates:
            for step in step_direction:
                job_name = '{}_{}_{}'.format(atom, coord, step)
                database['job_status'].update({job_name: 'not_started'})

    # END initialize_database()


def stat(db):
    """
        Checks displacement sub_directories for the status of each
        displacement computation

    db: (database) the database storing information for this distributed
        property calculation

    Returns: nothing
    Throws: nothing
    """
    n_finished = 0
    for job, status in db['job_status'].items():
        if status == 'finished':
            n_finished += 1
        elif status in ('not_started', 'running'):
            try:
                with open("{}/output.dat".format(job)) as outfile:
                    outfile.seek(-150, 2)
                    for line in outfile:
                        if 'Psi4 exiting successfully' in line:
                            db['job_status'][job] = 'finished'
                            n_finished += 1
                            break
                        else:
                            db['job_status'][job] = 'running'
            except:
                pass
    # check all jobs done?
    if n_finished == len(db['job_status'].keys()):
        db['jobs_complete'] = True

    # END stat()
