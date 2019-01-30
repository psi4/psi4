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
"""Module with a *procedures* dictionary specifying available quantum
chemical methods.
"""

from . import sapt
from . import proc
from . import interface_cfour

from psi4.driver.procrouting.dft import functionals, build_superfunctional_from_dictionary

# never import wrappers or aliases into this file

# Procedure lookup tables
procedures = {
    'energy': {
        'hf'            : proc.run_scf,
        'scf'           : proc.run_scf,
        'mcscf'         : proc.run_mcscf,
        'dcft'          : proc.run_dcft,
        'ep2'           : proc.run_dfep2,
        'mp3'           : proc.select_mp3,
        'mp2.5'         : proc.select_mp2p5,
        'mp2'           : proc.select_mp2,
        'omp2'          : proc.select_omp2,
        'scs-omp2'      : proc.run_occ,
        'scs(n)-omp2'   : proc.run_occ,
        'scs-omp2-vdw'  : proc.run_occ,
        'sos-omp2'      : proc.run_occ,
        'sos-pi-omp2'   : proc.run_occ,
        'omp3'          : proc.select_omp3,
        'scs-omp3'      : proc.run_occ,
        'scs(n)-omp3'   : proc.run_occ,
        'scs-omp3-vdw'  : proc.run_occ,
        'sos-omp3'      : proc.run_occ,
        'sos-pi-omp3'   : proc.run_occ,
        'olccd'         : proc.select_olccd,
        'omp2.5'        : proc.select_omp2p5,
        'dfocc'         : proc.run_dfocc,  # full control over dfocc
        'qchf'          : proc.run_qchf,
        'ccd'           : proc.run_dfocc,
        'sf-sapt'       : sapt.run_sf_sapt,
        'sapt(dft)'     : sapt.run_sapt_dft,
        'sapt0'         : proc.run_sapt,
        'ssapt0'        : proc.run_sapt,
        'sapt2'         : proc.run_sapt,
        'sapt2+'        : proc.run_sapt,
        'sapt2+(3)'     : proc.run_sapt,
        'sapt2+3'       : proc.run_sapt,
        'sapt2+(ccd)'   : proc.run_sapt,
        'sapt2+(3)(ccd)': proc.run_sapt,
        'sapt2+3(ccd)'  : proc.run_sapt,
        'sapt2+dmp2'    : proc.run_sapt,
        'sapt2+(3)dmp2' : proc.run_sapt,
        'sapt2+3dmp2'   : proc.run_sapt,
        'sapt2+(ccd)dmp2' : proc.run_sapt,
        'sapt2+(3)(ccd)dmp2' : proc.run_sapt,
        'sapt2+3(ccd)dmp2' : proc.run_sapt,
        'sapt0-ct'      : proc.run_sapt_ct,
        'sapt2-ct'      : proc.run_sapt_ct,
        'sapt2+-ct'     : proc.run_sapt_ct,
        'sapt2+(3)-ct'  : proc.run_sapt_ct,
        'sapt2+3-ct'    : proc.run_sapt_ct,
        'sapt2+(ccd)-ct'     : proc.run_sapt_ct,
        'sapt2+(3)(ccd)-ct'  : proc.run_sapt_ct,
        'sapt2+3(ccd)-ct'    : proc.run_sapt_ct,
        'fisapt0'       : proc.run_fisapt,
        'ccenergy'      : proc.run_ccenergy,  # full control over ccenergy
        'ccsd'          : proc.select_ccsd,
        'ccsd(t)'       : proc.select_ccsd_t_,
        'ccsd(at)'      : proc.select_ccsd_at_,
        'cc2'           : proc.run_ccenergy,
        'cc3'           : proc.run_ccenergy,
        'mrcc'          : proc.run_mrcc,  # interface to Kallay's MRCC program
        'bccd'          : proc.run_bccd,
        'bccd(t)'       : proc.run_bccd,
        'eom-ccsd'      : proc.run_eom_cc,
        'eom-cc2'       : proc.run_eom_cc,
        'eom-cc3'       : proc.run_eom_cc,
        'detci'         : proc.run_detci,  # full control over detci
        'mp'            : proc.run_detci,  # arbitrary order mp(n)
        'zapt'          : proc.run_detci,  # arbitrary order zapt(n)
        'cisd'          : proc.select_cisd,
        'cisdt'         : proc.run_detci,
        'cisdtq'        : proc.run_detci,
        'ci'            : proc.run_detci,  # arbitrary order ci(n)
        'fci'           : proc.run_detci,
        'casscf'        : proc.run_detcas,
        'rasscf'        : proc.run_detcas,
        'adc'           : proc.run_adc,
#        'cphf'          : proc.run_libfock,
#        'cis'           : proc.run_libfock,
#        'tdhf'          : proc.run_libfock,
#        'cpks'          : proc.run_libfock,
#        'tda'           : proc.run_libfock,
#        'tddft'         : proc.run_libfock,
        'psimrcc'       : proc.run_psimrcc,
        'psimrcc_scf'   : proc.run_psimrcc_scf,
        'qcisd'         : proc.run_fnocc,
        'qcisd(t)'      : proc.run_fnocc,
        'mp4'           : proc.select_mp4,
        'mp4(sdq)'      : proc.run_fnocc,
        'fno-ccsd'      : proc.select_fnoccsd,
        'fno-ccsd(t)'   : proc.select_fnoccsd_t_,
        'fno-qcisd'     : proc.run_fnocc,
        'fno-qcisd(t)'  : proc.run_fnocc,
        'fno-mp3'       : proc.run_fnocc,
        'fno-mp4(sdq)'  : proc.run_fnocc,
        'fno-mp4'       : proc.run_fnocc,
        'fno-lccd'      : proc.run_cepa,
        'fno-lccsd'     : proc.run_cepa,
        'fno-cepa(0)'   : proc.run_cepa,
        'fno-cepa(1)'   : proc.run_cepa,
        'fno-cepa(3)'   : proc.run_cepa,
        'fno-acpf'      : proc.run_cepa,
        'fno-aqcc'      : proc.run_cepa,
        'fno-cisd'      : proc.run_cepa,
        'lccd'          : proc.select_lccd,
        'lccsd'         : proc.run_cepa,
        'cepa(0)'       : proc.run_cepa,
        'cepa(1)'       : proc.run_cepa,
        'cepa(3)'       : proc.run_cepa,
        'acpf'          : proc.run_cepa,
        'aqcc'          : proc.run_cepa,
        'efp'           : proc.run_efp,
        'dmrg-scf'      : proc.run_dmrgscf,
        'dmrg-caspt2'   : proc.run_dmrgscf,
        'dmrg-ci'       : proc.run_dmrgci,
        # Upon adding a method to this list, add it to the docstring in energy() below
        # Aliases are discouraged. If you must add an alias to this list (e.g.,
        #    lccsd/cepa(0)), please search the whole driver to find uses of
        #    name in return values and psi variables and extend the logic to
        #    encompass the new alias.
    },
    'gradient' : {
        'hf'            : proc.run_scf_gradient,
        'scf'           : proc.run_scf_gradient,
        'cc2'           : proc.run_ccenergy_gradient,
        'ccsd'          : proc.select_ccsd_gradient,
        'ccsd(t)'       : proc.select_ccsd_t__gradient,
        'mp2'           : proc.select_mp2_gradient,
        'eom-ccsd'      : proc.run_eom_cc_gradient,
        'dcft'          : proc.run_dcft_gradient,
        'omp2'          : proc.select_omp2_gradient,
        'omp3'          : proc.select_omp3_gradient,
        'mp3'           : proc.select_mp3_gradient,
        'mp2.5'         : proc.select_mp2p5_gradient,
        'omp2.5'        : proc.select_omp2p5_gradient,
        'lccd'          : proc.select_lccd_gradient,
        'olccd'         : proc.select_olccd_gradient,
        'ccd'           : proc.run_dfocc_gradient,
        # Upon adding a method to this list, add it to the docstring in optimize() below
    },
    'hessian' : {
        # Upon adding a method to this list, add it to the docstring in frequency() below
        'hf'            : proc.run_scf_hessian,
        'scf'           : proc.run_scf_hessian,
    },
    'properties' : {
        'hf'       : proc.run_scf_property,
        'scf'      : proc.run_scf_property,
        'mp2'      : proc.select_mp2_property,
        'cc2'      : proc.run_cc_property,
        'ccsd'     : proc.run_cc_property,
        'eom-cc2'  : proc.run_cc_property,
        'eom-ccsd' : proc.run_cc_property,
        'detci'    : proc.run_detci_property,  # full control over detci
        'cisd'     : proc.run_detci_property,
        'cisdt'    : proc.run_detci_property,
        'cisdtq'   : proc.run_detci_property,
        'ci'       : proc.run_detci_property,  # arbitrary order ci(n)
        'fci'      : proc.run_detci_property,
        'rasscf'   : proc.run_detci_property,
        'casscf'   : proc.run_detci_property,
        # Upon adding a method to this list, add it to the docstring in property() below
    }} # yapf: disable

# Will only allow energy to be run for the following methods
energy_only_methods = [x for x in procedures['energy'].keys() if 'sapt' in x]
energy_only_methods += ['adc', 'efp', 'cphf', 'tdhf', 'cis']

# Integrate DFT with driver routines
for key in functionals:
    ssuper = build_superfunctional_from_dictionary(functionals[key], 1, 1, True)[0]

    procedures['energy'][key] = proc.run_scf

    # Properties
    if not ssuper.is_c_hybrid():
        procedures['properties'][key] = proc.run_scf_property

    # Gradients
    if not (ssuper.is_c_hybrid() or ssuper.is_c_lrc() or ssuper.needs_vv10()):
        procedures['gradient'][key] = proc.run_scf_gradient

    # Hessians
    if not ssuper.needs_xc():
        procedures['hessian'][key] = proc.run_scf_hessian

# Integrate CFOUR with driver routines
for ssuper in interface_cfour.cfour_list():
    procedures['energy'][ssuper.lower()] = interface_cfour.run_cfour

for ssuper in interface_cfour.cfour_gradient_list():
    procedures['gradient'][ssuper.lower()] = interface_cfour.run_cfour

# dictionary to register pre- and post-compute hooks for driver routines
hooks = dict((k1, dict((k2, []) for k2 in ['pre', 'post'])) for k1 in ['energy', 'optimize', 'frequency'])
