#
# @BEGIN LICENSE
#
# PSI4: an ab initio quantum chemistry software package
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

from __future__ import print_function
from __future__ import absolute_import
# Relative hack for now
import os, sys, inspect
path_dir = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"../")))
sys.path.append(path_dir)
from p4util.exceptions import *
import p4util
import psi4

def scf_set_reference_local(name):
    """
    Figures out the correct SCF reference to set locally
    """

    optstash = p4util.OptionsState(
        ['SCF', 'DFT_FUNCTIONAL'],
        ['SCF', 'SCF_TYPE'],
        ['SCF', 'REFERENCE'])

    # Alter default algorithm
    if not psi4.has_option_changed('SCF', 'SCF_TYPE'):
        psi4.set_local_option('SCF', 'SCF_TYPE', 'DF')

    if name == 'hf':
        if psi4.get_option('SCF','REFERENCE') == 'RKS':
            psi4.set_local_option('SCF','REFERENCE','RHF')
        elif psi4.get_option('SCF','REFERENCE') == 'UKS':
            psi4.set_local_option('SCF','REFERENCE','UHF')
    elif name == 'scf':
        if psi4.get_option('SCF','REFERENCE') == 'RKS':
            if (len(psi4.get_option('SCF', 'DFT_FUNCTIONAL')) > 0) or psi4.get_option('SCF', 'DFT_CUSTOM_FUNCTIONAL') is not None:
                pass
            else:
                psi4.set_local_option('SCF','REFERENCE','RHF')
        elif psi4.get_option('SCF','REFERENCE') == 'UKS':
            if (len(psi4.get_option('SCF', 'DFT_FUNCTIONAL')) > 0) or psi4.get_option('SCF', 'DFT_CUSTOM_FUNCTIONAL') is not None:
                pass
            else:
                psi4.set_local_option('SCF','REFERENCE','UHF')
    return optstash

def dft_set_reference_local(name):
    """
    Figures out the correct DFT reference to set locally
    """

    optstash = p4util.OptionsState(
        ['SCF', 'DFT_FUNCTIONAL'],
        ['SCF', 'REFERENCE'],
        ['SCF', 'SCF_TYPE'],
        ['DF_BASIS_MP2'],
        ['DFMP2', 'MP2_OS_SCALE'],
        ['DFMP2', 'MP2_SS_SCALE'])

    # Alter default algorithm
    if not psi4.has_option_changed('SCF', 'SCF_TYPE'):
        psi4.set_local_option('SCF', 'SCF_TYPE', 'DF')

    psi4.set_local_option('SCF', 'DFT_FUNCTIONAL', name)

    user_ref = psi4.get_option('SCF', 'REFERENCE')
    if (user_ref == 'RHF'):
        psi4.set_local_option('SCF', 'REFERENCE', 'RKS')
    elif (user_ref == 'UHF'):
        psi4.set_local_option('SCF', 'REFERENCE', 'UKS')
    elif (user_ref == 'ROHF'):
        raise ValidationError('ROHF reference for DFT is not available.')
    elif (user_ref == 'CUHF'):
        raise ValidationError('CUHF reference for DFT is not available.')

    return optstash

def oeprop_validator(prop_list):
    """
    Validations a list of OEProp computations. Throws if not found

    """
    oeprop_methods = ['DIPOLE', 'QUADRUPOLE', 'MULLIKEN_CHARGES', 'LOWDIN_CHARGES',
                      'WIBERG_LOWDIN_INDICES', 'MAYER_INDICES', 'MAYER_INDICES',
                      'MO_EXTENTS', 'GRID_FIELD', 'GRID_ESP', 'ESP_AT_NUCLEI',
                      'NO_OCCUPATIONS']

    if not len(prop_list):
        raise ValidationnError("OEProp: No properties specified!")


    for prop in prop_list:
        prop = prop.upper()

        if 'MULTIPOLE(' in prop: continue

        if prop not in oeprop_methods:
            alt_method_name = p4util.text.find_approximate_string_matches(prop,
                                                                     oeprop_methods, 2)
            alternatives = ""
            if len(alt_method_name) > 0:
                alternatives = " Did you mean? %s" % (" ".join(alt_method_name))

            raise ValidationError("OEProp: Feature '%s' is not recognized. %s" % (prop, alternatives))
    

def check_iwl_file_from_scf_type(scf_type, wfn):
    """
    Ensures that a IWL file has been written based on input SCF type.
    """

    if scf_type in ['DF', 'CD', 'PK', 'DIRECT']:
        mints = psi4.MintsHelper(wfn.basisset())
        mints.set_print(1)
        mints.integrals()

def check_non_symmetric_jk_density(name):
    """
    Ensure non-symmetric density matrices are supported for the selected JK routine.
    """
    scf_type = psi4.get_option('SCF', 'SCF_TYPE')
    supp_jk_type = ['DF', 'CD', 'PK', 'DIRECT', 'OUT_OF_CORE']
    supp_string = ', '.join(supp_jk_type[:-1]) + ', or ' + supp_jk_type[-1] + '.'

    if scf_type not in supp_jk_type:
        raise ValidationError("Method %s: Requires support for non-symmetric density matrices.\n"
                              "     Please set SCF_TYPE to %s" % (name, supp_string)) 



