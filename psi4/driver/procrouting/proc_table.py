#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2025 The Psi4 Developers.
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

from qcelemental.util import which

from . import interface_cfour, proc, proc_data, sapt
from .dft import build_superfunctional_from_dictionary, functionals

# never import wrappers or aliases into this file


# ADVICE upon adding to the `procedures` dict:
# * (1) add entry to `procedures` below. See ADVICE in psi4/driver/procrouting/proc.py on run_ vs. select_
# * (2) add entry to `method_governing_type_keywords` in psi4/driver/procrouting/proc_data.py
# * (3) add entry to table in docstring of `def energy`, etc. in psi4/driver/driver.py
# * (4) add entry to capabilities table in doc/sphinxman/source/introduction.rst
# * aliases discouraged but allowed. See `lccsd` and `a-ccsd(t)` for examples
# * for `hessian` entries, program up and set DIPOLE GRADIENT, too, otherwise IR intensities logic will fail

# Procedure lookup tables
procedures = {
    'energy': {
        'hf'            : proc.run_scf,
        'scf'           : proc.run_scf,
        'mcscf'         : proc.run_mcscf,
        'dct'           : proc.run_dct,
        'ep2'           : proc.run_dfep2,
        'mp2'           : proc.select_mp2,
        'scs-mp2'       : proc.run_occ,
        'scs(n)-mp2'    : proc.run_occ,
        'scs-mp2-vdw'   : proc.run_occ,
        'sos-mp2'       : proc.run_occ,
        'sos-pi-mp2'    : proc.run_occ,
        'custom-scs-mp2': proc.run_occ,
        'omp2'          : proc.select_omp2,
        'scs-omp2'      : proc.run_occ,
        'sos-omp2'       : proc.run_occ,
        'custom-scs-omp2' : proc.run_occ,
        'dlpno-mp2'     : proc.run_dlpnomp2,
        'scs-dlpno-mp2' : proc.run_dlpnomp2,
        'mp2.5'         : proc.select_mp2p5,
        'custom-scs-mp2.5' : proc.run_occ,
        'omp2.5'        : proc.select_omp2p5,
        'custom-scs-omp2.5' : proc.run_occ,
        'mp3'           : proc.select_mp3,
        'scs-mp3'       : proc.run_occ,
        'custom-scs-mp3' : proc.run_occ,
        'omp3'          : proc.select_omp3,
        'scs-omp3'      : proc.run_occ,
        'sos-omp3'      : proc.run_occ,
        'custom-scs-omp3' : proc.run_occ,
        'lccd'          : proc.select_lccd,
        'custom-scs-lccd' : proc.run_occ,
        'olccd'         : proc.select_olccd,
        'custom-scs-olccd' : proc.run_occ,
        'remp2'         : proc.select_remp2,
        'oremp2'        : proc.select_olccd,
        # 'dfocc'         : proc.run_dfocc,  # full control over dfocc  # canceled Jul 2022 as Error raising and not useful
        'qchf'          : proc.run_qchf,
        'ccd'           : proc.select_ccd,
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
        'ccsd(at)'      : proc.select_ccsd_at_,  # alias
        'a-ccsd(t)'     : proc.select_ccsd_at_,
        'lambda-ccsd(t)': proc.select_ccsd_at_,  # alias
        'ccsd(t)_l'     : proc.select_ccsd_at_,  # alias
        'cc2'           : proc.select_cc2,
        'cc3'           : proc.select_cc3,
        'mrcc'          : proc.run_mrcc,  # interface to Kallay's MRCC program  # Aug 2022 deprecated
        'bccd'          : proc.run_bccd,
        'bccd(t)'       : proc.run_bccd,
        'eom-ccsd'      : proc.run_eom_cc,
        'eom-cc2'       : proc.run_eom_cc,
        'eom-cc3'       : proc.run_eom_cc,
        'detci'         : proc.run_detci,  # full control over detci
        # 'mp'          : proc.run_detci,  # arbitrary order mp(n)  # Aug 2022 reworked below to add levels directly
        # 'zapt'        : proc.run_detci,  # arbitrary order zapt(n)  # Aug 2022 reworked below to add levels directly
        'cisd'          : proc.select_cisd,
        'cisdt'         : proc.run_detci,
        'cisdtq'        : proc.run_detci,
        # 'ci'          : proc.run_detci,  # arbitrary order ci(n)  # Aug 2022 reworked below to add levels directly
        'fci'           : proc.run_detci,
        'casscf'        : proc.run_detcas,
        'rasscf'        : proc.run_detcas,
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
        'lccsd'         : proc.run_cepa,
        'cepa(0)'       : proc.run_cepa,  # alias
        'cepa(1)'       : proc.run_cepa,
        'cepa(3)'       : proc.run_cepa,
        'acpf'          : proc.run_cepa,
        'aqcc'          : proc.run_cepa,
        'efp'           : proc.run_efp,
        'dmrg-scf'      : proc.run_dmrgscf,
        'dmrg-caspt2'   : proc.run_dmrgscf,
        'dmrg-ci'       : proc.run_dmrgci,
        'adc(1)'        : proc.run_adcc,
        'adc(2)'        : proc.run_adcc,
        'adc(2)-x'      : proc.run_adcc,
        'adc(3)'        : proc.run_adcc,
        'cvs-adc(1)'    : proc.run_adcc,
        'cvs-adc(2)'    : proc.run_adcc,
        'cvs-adc(2)-x'  : proc.run_adcc,
        'cvs-adc(3)'    : proc.run_adcc,
    },
    'gradient' : {
        'hf'            : proc.select_scf_gradient,
        'scf'           : proc.select_scf_gradient,
        'cc2'           : proc.select_cc2_gradient,
        'ccsd'          : proc.select_ccsd_gradient,
        'ccsd(t)'       : proc.select_ccsd_t__gradient,
        'mp2'           : proc.select_mp2_gradient,
        'eom-ccsd'      : proc.run_eom_cc_gradient,
        'dct'           : proc.run_dct_gradient,
        'omp2'          : proc.select_omp2_gradient,
        'omp3'          : proc.select_omp3_gradient,
        'mp3'           : proc.select_mp3_gradient,
        'mp2.5'         : proc.select_mp2p5_gradient,
        'omp2.5'        : proc.select_omp2p5_gradient,
        'mp2-d'         : proc.run_dfmp2d_gradient,
        'mp2d'          : proc.run_dfmp2d_gradient,  # alias to match dft aliasing
        'lccd'          : proc.select_lccd_gradient,
        'olccd'         : proc.select_olccd_gradient,
        'oremp2'        : proc.select_olccd_gradient,
        'ccd'           : proc.select_ccd_gradient,
    },
    'hessian' : {
        'hf'            : proc.run_scf_hessian,
        'scf'           : proc.run_scf_hessian,
    },
    'properties' : {
        'hf'           : proc.run_scf_property,
        'scf'          : proc.run_scf_property,
        'mp2'          : proc.select_mp2_property,
        'cc2'          : proc.run_cc_property,
        'ccsd'         : proc.run_cc_property,
        'eom-cc2'      : proc.run_cc_property,
        'eom-ccsd'     : proc.run_cc_property,
        'dct'          : proc.run_dct_property,
        'detci'        : proc.run_detci_property,  # full control over detci
        'cisd'         : proc.run_detci_property,
        'cisdt'        : proc.run_detci_property,
        'cisdtq'       : proc.run_detci_property,
        # 'ci'         : proc.run_detci_property,  # arbitrary order ci(n)  # Aug 2022 reworked below to add levels directly
        'fci'          : proc.run_detci_property,
        'rasscf'       : proc.run_detci_property,
        'casscf'       : proc.run_detci_property,
        'omp2'         : proc.select_omp2_property,
        'omp2.5'       : proc.select_omp2p5_property,
        'omp3'         : proc.select_omp3_property,
        'olccd'        : proc.select_olccd_property,
        # TODO 'oremp2'       : proc.select_olccd_property,
        'adc(1)'       : proc.run_adcc_property,
        'adc(2)'       : proc.run_adcc_property,
        'adc(2)-x'     : proc.run_adcc_property,
        'adc(3)'       : proc.run_adcc_property,
        'cvs-adc(1)'   : proc.run_adcc_property,
        'cvs-adc(2)'   : proc.run_adcc_property,
        'cvs-adc(2)-x' : proc.run_adcc_property,
        'cvs-adc(3)'   : proc.run_adcc_property,
    }} # yapf: disable

# Will only allow energy to be run for the following methods
energy_only_methods = [x for x in procedures['energy'].keys() if 'sapt' in x]
energy_only_methods += ['efp', 'cphf', 'tdhf', 'cis']

# Catch all SAPT-D variants
for key in functionals:
    # Grab the available -Ds from HF, since that's what SAPT0-D calls
    if key.startswith('hf-d'):
        disp = key.split('-')[-1]
        procedures['energy']['sapt0-' + disp] = proc.run_sapt
        procedures['energy']['fisapt0-' + disp] = proc.run_fisapt

# Will complete modelchem spec with basis='(auto)' for following methods
integrated_basis_methods = [
    'g2', 'gaussian-2',
    'hf3c', 'hf-3c',
    'pbeh3c', 'pbeh-3c',
    'b973c', 'b97-3c',
    'r2scan3c', 'r2scan-3c',
    'wb97x3c', 'wb97x-3c',
    'sns-mp2',
]

# Integrate arbitrary order with driver routines
for lvl in range(2, 99):
    procedures['energy'][f'ci{lvl}'] = proc.run_detci
    procedures['energy'][f'zapt{lvl}'] = proc.run_detci
    if lvl >= 5:
        procedures['energy'][f'mp{lvl}'] = proc.run_detci
    procedures['properties'][f'ci{lvl}'] = proc.run_detci_property

# Integrate MRCC with driver routines
if which("dmrcc", return_bool=True):
    for key in proc_data.mrcc_methods:
        if key not in ["cc2", "ccsd", "ccsd(t)", "cc3", "ccsd(t)_l"]:  # covered by select_ routines above
            procedures["energy"][key] = proc.select_mrcc

# Integrate DFT with driver routines
for key in functionals:
    ssuper = build_superfunctional_from_dictionary(functionals[key], 1, 1, True)[0]

    # Energy
    procedures['energy'][key] = proc.run_scf

    if not (ssuper.is_c_hybrid() or ssuper.is_c_lrc() or ssuper.needs_vv10()):
        procedures['energy']['td-' + key] = proc.run_tdscf_energy

    # Properties
    if not ssuper.is_c_hybrid():
        procedures['properties'][key] = proc.run_scf_property

    # Gradients
    if not (ssuper.is_c_hybrid() or ssuper.is_c_lrc() or ssuper.needs_vv10()):
        procedures['gradient'][key] = proc.select_scf_gradient

    # Hessians
    if not ssuper.is_gga(): # N.B. this eliminates both GGA and m-GGA, as the latter contains GGA terms
        procedures['hessian'][key] = proc.run_scf_hessian

# Integrate CFOUR with driver routines
for ssuper in interface_cfour.cfour_list():
    procedures['energy'][ssuper.lower()] = interface_cfour.run_cfour

for ssuper in interface_cfour.cfour_gradient_list():
    procedures['gradient'][ssuper.lower()] = interface_cfour.run_cfour

for ssuper in interface_cfour.cfour_hessian_list():
    procedures['hessian'][ssuper.lower()] = interface_cfour.run_cfour

# dictionary to register pre- and post-compute hooks for driver routines
hooks = dict((k1, dict((k2, []) for k2 in ['pre', 'post'])) for k1 in ['energy', 'optimize', 'frequency'])
