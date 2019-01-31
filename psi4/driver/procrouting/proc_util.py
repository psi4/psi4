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

import numpy as np

from psi4 import core
from psi4.driver import p4util
from psi4.driver.p4util.exceptions import *
from psi4.driver.procrouting.dft import functionals, build_superfunctional_from_dictionary

def scf_set_reference_local(name, is_dft=False):
    """
    Figures out the correct SCF reference to set locally
    """

    optstash = p4util.OptionsState(
        ['SCF_TYPE'],
        ['SCF', 'REFERENCE'])

    # Alter default algorithm
    if not core.has_global_option_changed('SCF_TYPE'):
        core.set_global_option('SCF_TYPE', 'DF')

    # Alter reference name if needed
    user_ref = core.get_option('SCF', 'REFERENCE')

    sup = build_superfunctional_from_dictionary(functionals[name], 1, 1, True)[0]
    if sup.needs_xc() or is_dft:
        if (user_ref == 'RHF'):
            core.set_local_option('SCF', 'REFERENCE', 'RKS')
        elif (user_ref == 'UHF'):
            core.set_local_option('SCF', 'REFERENCE', 'UKS')
        elif (user_ref == 'ROHF'):
            raise ValidationError('ROHF reference for DFT is not available.')
        elif (user_ref == 'CUHF'):
            raise ValidationError('CUHF reference for DFT is not available.')
    # else we are doing HF and nothing needs to be overloaded

    return optstash


def oeprop_validator(prop_list):
    """
    Validations a list of OEProp computations. Throws if not found

    """
    oeprop_methods = core.OEProp.valid_methods

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


    if scf_type in ['DF', 'DISK_DF', 'MEM_DF', 'CD', 'PK', 'DIRECT']:
        mints = core.MintsHelper(wfn.basisset())
        if core.get_global_option("RELATIVISTIC") in ["X2C", "DKH"]:
            rel_bas = core.BasisSet.build(wfn.molecule(), "BASIS_RELATIVISTIC",
                                          core.get_option("SCF", "BASIS_RELATIVISTIC"),
                                          "DECON", core.get_global_option('BASIS'),
                                          puream=wfn.basisset().has_puream())
            mints.set_rel_basisset(rel_bas)

        mints.set_print(1)
        mints.integrals()

def check_non_symmetric_jk_density(name):
    """
    Ensure non-symmetric density matrices are supported for the selected JK routine.
    """
    scf_type = core.get_global_option('SCF_TYPE')
    supp_jk_type = ['DF', 'DISK_DF', 'MEM_DF', 'CD', 'PK', 'DIRECT', 'OUT_OF_CORE']
    supp_string = ', '.join(supp_jk_type[:-1]) + ', or ' + supp_jk_type[-1] + '.'

    if scf_type not in supp_jk_type:
        raise ValidationError("Method %s: Requires support for non-symmetric density matrices.\n"
                              "     Please set SCF_TYPE to %s" % (name, supp_string))

def check_disk_df(name, optstash):

    optstash.add_option(['SCF_TYPE'])

    # Alter default algorithm
    if not core.has_global_option_changed('SCF_TYPE'):
        core.set_global_option('SCF_TYPE', 'DISK_DF')
        core.print_out("""    Method '%s' requires SCF_TYPE = DISK_DF, setting.\n""" % name)
    elif core.get_global_option('SCF_TYPE') == "DF":
        core.set_global_option('SCF_TYPE', 'DISK_DF')
        core.print_out("""    Method '%s' requires SCF_TYPE = DISK_DF, setting.\n""" % name)
    else:
        if core.get_global_option('SCF_TYPE') != "DISK_DF":
            raise ValidationError("  %s requires SCF_TYPE = DISK_DF, please use SCF_TYPE = DF to automatically choose the correct DFJK implementation." % name)

def print_ci_results(ciwfn, rname, scf_e, ci_e, print_opdm_no=False):
    """
    Printing for all CI Wavefunctions
    """

    # Print out energetics
    core.print_out("\n   ==> Energetics <==\n\n")
    core.print_out("    SCF energy =         %20.15f\n" % scf_e)
    if "CI" in rname:
        core.print_out("    Total CI energy =    %20.15f\n" % ci_e)
    elif "MP" in rname:
        core.print_out("    Total MP energy =    %20.15f\n" % ci_e)
    elif "ZAPT" in rname:
        core.print_out("    Total ZAPT energy =  %20.15f\n" % ci_e)
    else:
        core.print_out("    Total MCSCF energy = %20.15f\n" % ci_e)

    # Nothing to be done for ZAPT or MP
    if ("MP" in rname) or ("ZAPT" in rname):
        core.print_out("\n")
        return

    # Initial info
    ci_nroots = core.get_option("DETCI", "NUM_ROOTS")
    irrep_labels = ciwfn.molecule().irrep_labels()

    # Grab the D-vector
    dvec = ciwfn.D_vector()
    dvec.init_io_files(True)

    for root in range(ci_nroots):
        core.print_out("\n   ==> %s root %d information <==\n\n" % (rname, root))

        # Print total energy
        root_e = core.variable("CI ROOT %d TOTAL ENERGY" % (root))
        core.print_out("    %s Root %d energy =  %20.15f\n" % (rname, root, root_e))

        # Print natural occupations
        if print_opdm_no:
            core.print_out("\n   Active Space Natural occupation numbers:\n\n")

            occs_list = []
            r_opdm = ciwfn.get_opdm(root, root, "SUM", False)
            for h in range(len(r_opdm.nph)):
                if 0 in r_opdm.nph[h].shape:
                    continue
                nocc, rot = np.linalg.eigh(r_opdm.nph[h])
                for e in nocc:
                    occs_list.append((e, irrep_labels[h]))

            occs_list.sort(key=lambda x: -x[0])

            cnt = 0
            for value, label in occs_list:
                value, label = occs_list[cnt]
                core.print_out("      %4s  % 8.6f" % (label, value))
                cnt += 1
                if (cnt % 3) == 0:
                    core.print_out("\n")

            if (cnt % 3):
                core.print_out("\n")

        # Print CIVector information
        ciwfn.print_vector(dvec, root)

    # True to keep the file
    dvec.close_io_files(True)


def prepare_sapt_molecule(sapt_dimer, sapt_basis):
    """
    Prepares a dimer molecule for a SAPT computations. Returns the dimer, monomerA, and monomerB.
    """

    # Shifting to C1 so we need to copy the active molecule
    sapt_dimer = sapt_dimer.clone()
    if sapt_dimer.schoenflies_symbol() != 'c1':
        core.print_out('  SAPT does not make use of molecular symmetry, further calculations in C1 point group.\n')
        sapt_dimer.reset_point_group('c1')
        sapt_dimer.fix_orientation(True)
        sapt_dimer.fix_com(True)
        sapt_dimer.update_geometry()
    else:
        sapt_dimer.update_geometry()  # make sure since mol from wfn, kwarg, or P::e
        sapt_dimer.fix_orientation(True)
        sapt_dimer.fix_com(True)

    nfrag = sapt_dimer.nfragments()

    if nfrag == 3:
        # Midbond case
        if sapt_basis == 'monomer':
            raise ValidationError("SAPT basis cannot both be monomer centered and have midbond functions.")

        midbond = sapt_dimer.extract_subsets(3)
        ztotal = 0
        for n in range(midbond.natom()):
            ztotal += midbond.Z(n)

        if ztotal > 0:
            raise ValidationError("SAPT third monomr must be a midbond function (all ghosts).")

        ghosts = ([2, 3], [1, 3])
    elif nfrag == 2:
        # Classical dimer case
        ghosts = (2, 1)
    else:
        raise ValidationError('SAPT requires active molecule to have 2 fragments, not %s.' % (nfrag))

    if sapt_basis == 'dimer':
        monomerA = sapt_dimer.extract_subsets(1, ghosts[0])
        monomerA.set_name('monomerA')
        monomerB = sapt_dimer.extract_subsets(2, ghosts[1])
        monomerB.set_name('monomerB')
    elif sapt_basis == 'monomer':
        monomerA = sapt_dimer.extract_subsets(1)
        monomerA.set_name('monomerA')
        monomerB = sapt_dimer.extract_subsets(2)
        monomerB.set_name('monomerB')
    else:
        raise ValidationError("SAPT basis %s not recognized" % sapt_basis)

    return (sapt_dimer, monomerA, monomerB)
