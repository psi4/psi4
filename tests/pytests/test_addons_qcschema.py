import json
import pprint

import pytest
from addons import hardware_nvidia_gpu, uusing, using
import qcengine as qcng

import psi4

# Notes
# * options-setting NOT cummulative if a run_qcschema in between

# Generating
# * equivalent to test_psi4. copy over the job, then run below to generate atomicinput
#    atin = psi4.driver.p4util.state_to_atomicinput(driver="energy", method="ccsd", molecule=ethene_ethyne)
#    print(f'    jatin = """{atin.serialize("json")}"""')
#    assert 0
# * switch to json running
#    atres = psi4.schema_wrapper.run_qcschema(json.loads(jatin))
#    pprint.pprint(atres.dict())

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick, pytest.mark.smoke]


@uusing("gdma")
def test_gdma():
    """gdma1"""
    #! Water RHF/cc-pVTZ distributed multipole analysis

    ref_energy = -76.0571685433842219
    ref_dma_mat = psi4.core.Matrix.from_list([[
        -0.43406697290168, -0.18762673939633, 0.00000000000000, 0.00000000000000, 0.03206686487531, 0.00000000000000,
        -0.00000000000000, -0.53123477172696, 0.00000000000000
    ],
                                              [
                                                  0.21703348903257, -0.06422316619952, 0.00000000000000,
                                                  -0.11648289410022, 0.01844320206227, 0.00000000000000,
                                                  0.07409226544133, -0.07115302332866, 0.00000000000000
                                              ],
                                              [
                                                  0.21703348903257, -0.06422316619952, 0.00000000000000,
                                                  0.11648289410022, 0.01844320206227, 0.00000000000000,
                                                  -0.07409226544133, -0.07115302332866, 0.00000000000000
                                              ]])
    ref_tot_mat = psi4.core.Matrix(1, 9)
    ref_tot_mat.name = "Reference total values"
    ref_tot_arr = [
        0.00000000516346, -0.79665315928128, 0.00000000000000, 0.00000000000000, 0.10813259329390, 0.00000000000000,
        0.00000000000000, -2.01989585894142, 0.00000000000000
    ]
    for i in range(9):
        ref_tot_mat.set(0, i, ref_tot_arr[i])

    # noreorient/nocom are not needed, but are used here to guarantee that the
    #   GDMA origin placement defined below is at the O atom.

    # added protocols.wavefunction = orbitals_and_eigenvalues (needs to be all)
    jatin = """{"id": null, "schema_name": "qcschema_input", "schema_version": 1, "molecule": {"schema_name": "qcschema_molecule", "schema_version": 2, "validated": true, "symbols": ["O", "H", "H"], "geometry": [0.0, 0.0, 0.22143054847664648, 4.379423262771008e-17, -1.4304281906653031, -0.8857259733588368, -4.379423262771008e-17, 1.4304281906653031, -0.8857259733588368], "name": "H2O", "molecular_charge": 0.0, "molecular_multiplicity": 1, "masses": [15.99491461957, 1.00782503223, 1.00782503223], "real": [true, true, true], "atom_labels": ["", "", ""], "atomic_numbers": [8, 1, 1], "mass_numbers": [16, 1, 1], "fragments": [[0, 1, 2]], "fragment_charges": [0.0], "fragment_multiplicities": [1], "fix_com": true, "fix_orientation": true, "provenance": {"creator": "QCElemental", "version": "v0.17.0+7.gf55d5ac.dirty", "routine": "qcelemental.molparse.from_string"}}, "driver": "energy", "model": {"method": "scf", "basis": "CC-PVTZ"}, "keywords": {"d_convergence": 1e-10, "gdma_limit": 2, "gdma_origin": [0.0, 0.0, 0.117176], "gdma_radius": ["H", 0.65], "gdma_switch": 0.0, "scf_type": "PK"}, "protocols": {"wavefunction": "orbitals_and_eigenvalues"}, "extras": {"wfn_qcvars_only": true}, "provenance": {"creator": "Psi4", "version": "1.4a2.dev1089", "routine": "psi4.driver.p4util.procutil"}}"""

    atres = psi4.schema_wrapper.run_qcschema(json.loads(jatin))

    try:
        psi4.gdma(atres.wavefunction)
    except TypeError as e:
        if "Invoked with: WavefunctionProperties" in str(e):
            pytest.xfail("GDMA not processable from AtomicResult.wavefunction")

    dmavals = psi4.core.variable("DMA DISTRIBUTED MULTIPOLES")
    totvals = psi4.core.variable("DMA TOTAL MULTIPOLES")
    assert psi4.compare_values(ref_energy, energy, 8, "SCF Energy")
    assert psi4.compare_matrices(ref_dma_mat, dmavals, 6, "DMA Distributed Multipoles")
    assert psi4.compare_matrices(ref_tot_mat, totvals, 6, "DMA Total Multipoles")


@uusing("ipi")
def test_ipi_broker1():
    """ipi_broker1"""

    pytest.xfail("IPI doesn't use basic psi4 functions, so can't transmit as AtomicInput")


#    water = psi4.geometry("""
#      O -1.216  -0.015  -0.261
#      H -1.946   0.681  -0.378
#      H -1.332  -0.754   0.283
#      units angstrom
#      no_reorient
#      no_com
#    """)
#
#    psi4.set_options({
#        'basis': 'sto-3g',
#        'reference': 'rhf',
#    })
#
#    options = {}
#
#    #ipi_broker(serverdata="inet:localhost:21340", options=options)
#    b = psi4.ipi_broker("ccsd", serverdata=False, options=options)
#
#    refnuc   =   9.05843673637
#    refscf   = -74.9417588868628
#    refccsd  = -0.04895074370294
#    reftotal = -74.9907096305658
#
#    frc = [[ 0.08704801,  0.1067644 , -0.11170374],
#           [-0.02216499, -0.03279655,  0.03215871],
#           [-0.06488302, -0.07396785,  0.07954503]]
#
#    b.calculate_force()
#
##    atin = psi4.driver.p4util.state_to_atomicinput(driver="energy", method="ccsd", molecule=ethene_ethyne)
##    print(f'    jatin = """{atin.serialize("json")}"""')
##    assert 0
##
##    atres = psi4.schema_wrapper.run_qcschema(json.loads(jatin))
##    pprint.pprint(atres.dict())
#
#
#    assert psi4.compare_values(refnuc,   water.nuclear_repulsion_energy(),              3, "Nuclear repulsion energy")
#    assert psi4.compare_values(refscf,   psi4.core.variable("SCF total energy"),        5, "SCF energy")
#    assert psi4.compare_values(refccsd,  psi4.core.variable("CCSD correlation energy"), 4, "CCSD contribution")
#    assert psi4.compare_values(reftotal, psi4.core.variable("Current energy"),          7, "Total energy")
#    assert psi4.compare_values(reftotal, b._potential,                                  7, "Total energy (Broker)")
#    assert psi4.compare_arrays(frc,      b._force,                                      4, "Total force (Broker)")
#
#    water_mirror = psi4.geometry("""
#      O  1.216   0.015   0.261
#      H  1.946  -0.681   0.378
#      H  1.332   0.754  -0.283
#      units angstrom
#      no_reorient
#      no_com
#    """)
#
#    b.calculate_force()
#
##    atin = psi4.driver.p4util.state_to_atomicinput(driver="energy", method="ccsd", molecule=ethene_ethyne)
##    print(f'    jatin = """{atin.serialize("json")}"""')
##    assert 0
##
##    atres = psi4.schema_wrapper.run_qcschema(json.loads(jatin))
##    pprint.pprint(atres.dict())
#
#
#    assert psi4.compare_values(refnuc,   water_mirror.nuclear_repulsion_energy(),    3, "Nuclear repulsion energy")
#    assert psi4.compare_values(refscf,   psi4.core.variable("SCF total energy"),        5, "SCF energy")
#    assert psi4.compare_values(refccsd,  psi4.core.variable("CCSD correlation energy"), 4, "CCSD contribution")
#    assert psi4.compare_values(reftotal, psi4.core.variable("Current energy"),          7, "Total energy")
#    assert psi4.compare_values(reftotal, b._potential,                                  7, "Total energy (Broker)")
#    assert psi4.compare_arrays(frc,     -b._force,                                      4, "Total force (Broker)")

#@uusing("mrcc")
#def test_mrcc():
#    """mrcc/ccsdt"""
#    #! CCSDT cc-pVDZ energy for the H2O molecule using MRCC
#
#    h2o = psi4.geometry("""
#        o
#        h 1 1.0
#        h 1 1.0 2 104.5
#    """)
#
#    psi4.set_options({
#        'basis': 'cc-pvdz',
#        'freeze_core': 'true'})
#
#    psi4.energy('mrccsdt')
#
##    atin = psi4.driver.p4util.state_to_atomicinput(driver="energy", method="ccsd", molecule=ethene_ethyne)
##    print(f'    jatin = """{atin.serialize("json")}"""')
##    assert 0
##
##    atres = psi4.schema_wrapper.run_qcschema(json.loads(jatin))
##    pprint.pprint(atres.dict())
#
#
#    assert psi4.compare_values(  8.801465529972, psi4.variable("NUCLEAR REPULSION ENERGY"), 6, 'NRE')
#    assert psi4.compare_values(-76.021418445155, psi4.variable("SCF TOTAL ENERGY"), 6, 'SCF')
#    assert psi4.compare_values( -0.204692406830, psi4.variable("MP2 CORRELATION ENERGY") , 6, 'MP2 correlation')
#    assert psi4.compare_values( -0.217715210258, psi4.variable("CCSDT CORRELATION ENERGY"), 6, 'CCSDT correlation')
#    assert psi4.compare_values(-76.239133655413, psi4.variable("CURRENT ENERGY"), 6, 'CCSDT')


@uusing("chemps2")
def test_chemps2():
    """chemps2/scf-n2"""
    #! dmrg-scf on N2

    jatin = """{"id": null, "schema_name": "qcschema_input", "schema_version": 1, "molecule": {"schema_name": "qcschema_molecule", "schema_version": 2, "validated": true, "symbols": ["N", "N"], "geometry": [0.0, 0.0, -1.059, 0.0, 0.0, 1.059], "name": "N2", "molecular_charge": 0.0, "molecular_multiplicity": 1, "masses": [14.00307400443, 14.00307400443], "real": [true, true], "atom_labels": ["", ""], "atomic_numbers": [7, 7], "mass_numbers": [14, 14], "fragments": [[0, 1]], "fragment_charges": [0.0], "fragment_multiplicities": [1], "fix_com": false, "fix_orientation": false, "provenance": {"creator": "QCElemental", "version": "v0.17.0+7.gf55d5ac.dirty", "routine": "qcelemental.molparse.from_string"}}, "driver": "energy", "model": {"method": "dmrg-scf", "basis": "CC-PVDZ"}, "keywords": {"active": [2, 0, 1, 1, 0, 2, 1, 1], "dmrg_diis": 1, "dmrg_diis_write": 1, "dmrg_excitation": 0, "dmrg_irrep": 0, "dmrg_local_init": 0, "dmrg_mps_write": 0, "dmrg_multiplicity": 1, "dmrg_print_corr": 1, "dmrg_scf_active_space": "NO", "dmrg_scf_diis_thr": 0.01, "dmrg_scf_state_avg": 0, "dmrg_sweep_dvdson_rtol": [0.0001, 1e-06, 1e-08], "dmrg_sweep_energy_conv": [1e-10, 1e-10, 1e-10], "dmrg_sweep_max_sweeps": [5, 5, 10], "dmrg_sweep_noise_prefac": [0.05, 0.05, 0.0], "dmrg_sweep_states": [500, 1000, 1000], "dmrg_unitary_write": 1, "d_convergence": 1e-12, "e_convergence": 1e-12, "reference": "RHF", "restricted_docc": [1, 0, 0, 0, 0, 1, 0, 0]}, "protocols": {}, "extras": {"wfn_qcvars_only": true}, "provenance": {"creator": "Psi4", "version": "1.4a2.dev1089", "routine": "psi4.driver.p4util.procutil"}}"""

    atres = psi4.schema_wrapper.run_qcschema(json.loads(jatin))
    assert psi4.compare_values(-109.1035023353, atres.return_result, 6, "DMRG Energy")


@uusing("mp2d")
def test_mp2d():

    jatin = """{"id": null, "schema_name": "qcschema_input", "schema_version": 1, "molecule": {"schema_name": "qcschema_molecule", "schema_version": 2, "validated": true, "symbols": ["C", "C", "H", "H", "H", "H", "C", "C", "H", "H"], "geometry": [3.8623510442091865e-17, -1.261539587380886, -4.02163254721969, -3.862351044209187e-17, 1.261539587380886, -4.02163254721969, 1.745390733721485, -2.3286206872737854, -4.0245162692871395, -1.7453907337214847, -2.3286206872737854, -4.0245162692871395, -1.745390733721485, 2.3286206872737854, -4.0245162692871395, 1.7453907337214847, 2.3286206872737854, -4.0245162692871395, -5.7777898331617076e-34, 0.0, 5.47454736883822, -5.7777898331617076e-34, 0.0, 3.193150937439626, -5.7777898331617076e-34, 0.0, 1.1789145370276326, -5.7777898331617076e-34, 0.0, 7.484131263529336], "name": "C4H6", "molecular_charge": 0.0, "molecular_multiplicity": 1, "masses": [12.0, 12.0, 1.00782503223, 1.00782503223, 1.00782503223, 1.00782503223, 12.0, 12.0, 1.00782503223, 1.00782503223], "real": [true, true, true, true, true, true, true, true, true, true], "atom_labels": ["", "", "", "", "", "", "", "", "", ""], "atomic_numbers": [6, 6, 1, 1, 1, 1, 6, 6, 1, 1], "mass_numbers": [12, 12, 1, 1, 1, 1, 12, 12, 1, 1], "fragments": [[0, 1, 2, 3, 4, 5], [6, 7, 8, 9]], "fragment_charges": [0.0, 0.0], "fragment_multiplicities": [1, 1], "fix_com": false, "fix_orientation": false, "provenance": {"creator": "QCElemental", "version": "v0.17.0+7.gf55d5ac.dirty", "routine": "qcelemental.molparse.from_string"}}, "driver": "energy", "model": {"method": "mp2-d", "basis": "cc-pvdz"}, "keywords": {}, "protocols": {}, "extras": {"wfn_qcvars_only": true}, "provenance": {"creator": "Psi4", "version": "1.4a2.dev1089", "routine": "psi4.driver.p4util.procutil"}}"""

    atres = psi4.schema_wrapper.run_qcschema(json.loads(jatin))

    expected = 0.00632174635953
    assert psi4.compare_values(expected, atres.extras["qcvars"]["DISPERSION CORRECTION ENERGY"], 7, 'disp E')
    assert psi4.compare_values(expected, atres.properties.scf_dispersion_correction_energy, 7, 'mp2d disp E')


@uusing("dftd3")
def test_dftd3():
    """dftd3/energy"""

    ref_d2 = [-0.00390110, -0.00165271, -0.00058118]
    ref_d3zero = [-0.00285088, -0.00084340, -0.00031923]
    ref_d3bj = [-0.00784595, -0.00394347, -0.00226683]

    ref_pbe_d2 = [-0.00278650, -0.00118051, -0.00041513]
    ref_pbe_d3zero = [-0.00175474, -0.00045421, -0.00016839]
    ref_pbe_d3bj = [-0.00475937, -0.00235265, -0.00131239]

    # mA, b3lyp-d2, libdisp
    jatin = """{"id": null, "schema_name": "qcschema_input", "schema_version": 1, "molecule": {"schema_name": "qcschema_molecule", "schema_version": 2, "validated": true, "symbols": ["C", "C", "H", "H", "H", "H"], "geometry": [3.8623510442091865e-17, -1.261539587380886, 0.00041472029798115525, -3.862351044209187e-17, 1.261539587380886, 0.00041472029798115525, 1.745390733721485, -2.3286206872737854, -0.0024690017694669127, -1.7453907337214847, -2.3286206872737854, -0.0024690017694669127, -1.745390733721485, 2.3286206872737854, -0.0024690017694669127, 1.7453907337214847, 2.3286206872737854, -0.0024690017694669127], "name": "C2H4", "molecular_charge": 0.0, "molecular_multiplicity": 1, "masses": [12.0, 12.0, 1.00782503223, 1.00782503223, 1.00782503223, 1.00782503223], "real": [true, true, true, true, true, true], "atom_labels": ["", "", "", "", "", ""], "atomic_numbers": [6, 6, 1, 1, 1, 1], "mass_numbers": [12, 12, 1, 1, 1, 1], "fragments": [[0, 1, 2, 3, 4, 5]], "fragment_charges": [0.0], "fragment_multiplicities": [1], "fix_com": false, "fix_orientation": false, "provenance": {"creator": "QCElemental", "version": "v0.17.0+7.gf55d5ac.dirty", "routine": "qcelemental.molparse.from_string"}}, "driver": "energy", "model": {"method": "b3lyp-d2", "basis": "STO-3G"}, "keywords": {"dft_radial_points": 50, "dft_spherical_points": 110, "scf_type": "DF", "function_kwargs": {"engine": "libdisp"}}, "protocols": {}, "extras": {"wfn_qcvars_only": true}, "provenance": {"creator": "Psi4", "version": "1.4a2.dev1089", "routine": "psi4.driver.p4util.procutil"}}"""

    atres = psi4.schema_wrapper.run_qcschema(json.loads(jatin))
    assert psi4.compare_values(ref_d2[1], atres.extras["qcvars"]["DISPERSION CORRECTION ENERGY"], 7,
                               'Ethene -D2 (calling psi4 Disp class)')

    # mA, b3lyp-d3bj, dftd3
    jatin = """{"id": null, "schema_name": "qcschema_input", "schema_version": 1, "molecule": {"schema_name": "qcschema_molecule", "schema_version": 2, "validated": true, "symbols": ["C", "C", "H", "H", "H", "H"], "geometry": [3.8623510442091865e-17, -1.261539587380886, 0.00041472029798115525, -3.862351044209187e-17, 1.261539587380886, 0.00041472029798115525, 1.745390733721485, -2.3286206872737854, -0.0024690017694669127, -1.7453907337214847, -2.3286206872737854, -0.0024690017694669127, -1.745390733721485, 2.3286206872737854, -0.0024690017694669127, 1.7453907337214847, 2.3286206872737854, -0.0024690017694669127], "name": "C2H4", "molecular_charge": 0.0, "molecular_multiplicity": 1, "masses": [12.0, 12.0, 1.00782503223, 1.00782503223, 1.00782503223, 1.00782503223], "real": [true, true, true, true, true, true], "atom_labels": ["", "", "", "", "", ""], "atomic_numbers": [6, 6, 1, 1, 1, 1], "mass_numbers": [12, 12, 1, 1, 1, 1], "fragments": [[0, 1, 2, 3, 4, 5]], "fragment_charges": [0.0], "fragment_multiplicities": [1], "fix_com": false, "fix_orientation": false, "provenance": {"creator": "QCElemental", "version": "v0.17.0+7.gf55d5ac.dirty", "routine": "qcelemental.molparse.from_string"}}, "driver": "energy", "model": {"method": "B3LYP-D3BJ", "basis": "STO-3G"}, "keywords": {"dft_radial_points": 50, "dft_spherical_points": 110, "scf_type": "DF"}, "protocols": {}, "extras": {"wfn_qcvars_only": true}, "provenance": {"creator": "Psi4", "version": "1.4a2.dev1089", "routine": "psi4.driver.p4util.procutil"}}"""

    atres = psi4.schema_wrapper.run_qcschema(json.loads(jatin))
    assert psi4.compare_values(ref_d3bj[1], atres.extras["qcvars"]["DISPERSION CORRECTION ENERGY"], 7,
                               'Ethene -D3 (calling dftd3 -bj)')

    # mB, pbe-d3bj, custom parmeters
    jatin = """{"id": null, "schema_name": "qcschema_input", "schema_version": 1, "molecule": {"schema_name": "qcschema_molecule", "schema_version": 2, "validated": true, "symbols": ["C", "C", "H", "H"], "geometry": [3.7092061506874214e-68, -3.7092061506874214e-68, 1.1408784499704554, 3.7092061506874214e-68, -3.7092061506874214e-68, -1.140517981428138, 3.7092061506874214e-68, -3.7092061506874214e-68, -3.154754381840132, 3.7092061506874214e-68, -3.7092061506874214e-68, 3.1504623446615714], "name": "C2H2", "molecular_charge": 0.0, "molecular_multiplicity": 1, "masses": [12.0, 12.0, 1.00782503223, 1.00782503223], "real": [true, true, true, true], "atom_labels": ["", "", "", ""], "atomic_numbers": [6, 6, 1, 1], "mass_numbers": [12, 12, 1, 1], "fragments": [[0, 1, 2, 3]], "fragment_charges": [0.0], "fragment_multiplicities": [1], "fix_com": false, "fix_orientation": false, "provenance": {"creator": "QCElemental", "version": "v0.17.0+7.gf55d5ac.dirty", "routine": "qcelemental.molparse.from_string"}}, "driver": "energy", "model": {"method": "pbe-d3bj", "basis": "STO-3G"}, "keywords": {"dft_dispersion_parameters": [2.0, 0.7875, 0.4289, 4.4407], "dft_radial_points": 50, "dft_spherical_points": 110, "scf_type": "DF"}, "protocols": {}, "extras": {"wfn_qcvars_only": true}, "provenance": {"creator": "Psi4", "version": "1.4a2.dev1089", "routine": "psi4.driver.p4util.procutil"}}"""

    atres = psi4.schema_wrapper.run_qcschema(json.loads(jatin))
    assert psi4.compare_values(-0.002238400, atres.extras["qcvars"]["DISPERSION CORRECTION ENERGY"], 7,
                               'Ethene -D3 (calling dftd3 -bj)')


@uusing("libefp")
def test_libefp():
    """libefp/qchem-qmefp-sp"""
    #! EFP on mixed QM (water) and EFP (water + 2 * ammonia) system.
    #! An EFP-only calc performed first to test vales against q-chem.

    pytest.xfail("EFP not transmittable through QCSchema")

    qmefp = psi4.geometry("""
    # QM fragment
    0 1
    units bohr
    O1     0.000000000000     0.000000000000     0.224348285559
    H2    -1.423528800232     0.000000000000    -0.897393142237
    H3     1.423528800232     0.000000000000    -0.897393142237
    # EFP as EFP fragments
    --
    efp h2o -4.014110144291     2.316749370493    -1.801514729931 -2.902133 1.734999 -1.953647
    --
    efp NH3,1.972094713645,,3.599497221584 ,    5.447701074734 -1.105309 2.033306 -1.488582
    --
    efp NH3 -7.876296399270    -1.854372164887    -2.414804197762  2.526442 1.658262 -2.742084
    """)

    #  <<<  EFP calc  >>>
    psi4.set_options({'basis': '6-31g*', 'scf_type': 'pk', 'guess': 'core', 'df_scf_guess': False})

    #    psi4.energy('efp')
    #
    atin = psi4.driver.p4util.state_to_atomicinput(driver="energy", method="efp")
    print(f'    jatin = """{atin.serialize("json")}"""')
    assert 0


##    atres = psi4.schema_wrapper.run_qcschema(json.loads(jatin))
##    pprint.pprint(atres.dict())
#
#
#    assert psi4.compare_values( 9.1793879214, qmefp.nuclear_repulsion_energy(), 6, 'QM NRE')
#    assert psi4.compare_values(-0.0004901368, psi4.variable('efp elst energy'), 6, 'EFP-EFP Elst')  # from q-chem
#    assert psi4.compare_values(-0.0003168768, psi4.variable('efp ind energy'), 6, 'EFP-EFP Indc')
#    assert psi4.compare_values(-0.0021985285, psi4.variable('efp disp energy'), 6, 'EFP-EFP Disp')  # from q-chem
#    assert psi4.compare_values( 0.0056859871, psi4.variable('efp exch energy'), 6, 'EFP-EFP Exch')  # from q-chem
#    assert psi4.compare_values( 0.0026804450, psi4.variable('efp total energy'), 6, 'EFP-EFP Totl')
#    assert psi4.compare_values( 0.0026804450, psi4.variable('current energy'), 6, 'Current')
#    psi4.core.print_variables()
#
#    psi4.core.clean()
#    psi4.core.clean_variables()
#
#    #  <<<  QM + EFP calc  >>>
#    psi4.set_options({
#        'e_convergence': 12,
#        'd_convergence': 12})
#    psi4.energy('scf')
#
#    assert psi4.compare_values( 9.1793879214, qmefp.nuclear_repulsion_energy(), 6, 'QM NRE')
#    assert psi4.compare_values( 0.2622598847, psi4.variable('efp total energy') - psi4.variable('efp ind energy'), 6, 'EFP corr to SCF')  # from q-chem
#    assert psi4.compare_values(-0.0117694790, psi4.variable('efp ind energy'), 6, 'QM-EFP Indc')  # from q-chem
#    assert psi4.compare_values(-0.0021985285, psi4.variable('efp disp energy'), 6, 'EFP-EFP Disp')  # from q-chem
#    assert psi4.compare_values( 0.0056859871, psi4.variable('efp exch energy'), 6, 'EFP-EFP Exch')  # from q-chem
#    assert psi4.compare_values( 0.2504904057, psi4.variable('efp total energy'), 6, 'EFP-EFP Totl')  # from q-chem
#    assert psi4.compare_values(-76.0139362744, psi4.variable('scf total energy'), 6, 'SCF')  # from q-chem
#    psi4.core.print_variables()


@uusing("pcmsolver")
def test_pcmsolver():
    """pcmsolver/scf"""
    #! pcm

    nucenergy = 12.0367196636183458
    polenergy = -0.0053060443528559
    totalenergy = -55.4559426361734040

    #    pcm_string = """
    #       Units = Angstrom
    #       Medium {
    #       SolverType = IEFPCM
    #       Solvent = Water
    #       }
    #
    #       Cavity {
    #       RadiiSet = UFF
    #       Type = GePol
    #       Scaling = False
    #       Area = 0.3
    #       Mode = Implicit
    #       }
    #    """
    #    psi4.set_options({
    #      'basis': 'STO-3G',
    #      'scf_type': 'pk',
    #      'pcm': True,
    #      'pcm_scf_type': 'total',
    #    "pcm__input": pcm_string,
    #    })
    #pcmfile = "\\n".join(pcm_string.split("\n"))
    #print(f"{pcmfile=}")

    jatin = """{"id": null, "schema_name": "qcschema_input", "schema_version": 1, "molecule": {"schema_name": "qcschema_molecule", "schema_version": 2, "validated": true, "symbols": ["N", "H", "H", "H"], "geometry": [-1e-10, -0.1040380466, 0.0, -0.9015844116, 0.4818470201, -1.5615900098, -0.9015844116, 0.4818470201, 1.5615900098, 1.8031688251, 0.4818470204, 0.0], "name": "H3N", "molecular_charge": 0.0, "molecular_multiplicity": 1, "masses": [14.00307400443, 1.00782503223, 1.00782503223, 1.00782503223], "real": [true, true, true, true], "atom_labels": ["", "", "", ""], "atomic_numbers": [7, 1, 1, 1], "mass_numbers": [14, 1, 1, 1], "fragments": [[0, 1, 2, 3]], "fragment_charges": [0.0], "fragment_multiplicities": [1], "fix_com": false, "fix_orientation": false, "fix_symmetry": "c1", "provenance": {"creator": "QCElemental", "version": "v0.17.0+7.gf55d5ac.dirty", "routine": "qcelemental.molparse.from_string"}}, "driver": "energy", "model": {"method": "scf", "basis": "STO-3G"}, "keywords": {"pcm": 1, "pcm__input": "\\n       Units = Angstrom\\n       Medium {\\n       SolverType = IEFPCM\\n       Solvent = Water\\n       }\\n    \\n       Cavity {\\n       RadiiSet = UFF\\n       Type = GePol\\n       Scaling = False\\n       Area = 0.3\\n       Mode = Implicit\\n       }\\n    ", "pcm_scf_type": "TOTAL", "scf_type": "PK"}, "protocols": {}, "extras": {"wfn_qcvars_only": true}, "provenance": {"creator": "Psi4", "version": "1.4a2.dev1089", "routine": "psi4.driver.p4util.procutil"}}"""

    atres1 = psi4.schema_wrapper.run_qcschema(json.loads(jatin))
    pprint.pprint(atres1.dict())
    assert psi4.compare_values(nucenergy, atres1.properties.nuclear_repulsion_energy, 10,
                               "Nuclear repulsion energy (PCM, total algorithm)")
    assert psi4.compare_values(totalenergy, atres1.return_result, 10, "Total energy (PCM, total algorithm)")
    assert psi4.compare_values(polenergy, atres1.extras["qcvars"]["PCM POLARIZATION ENERGY"], 6,
                               "Polarization energy (PCM, total algorithm)")

    jatin = """{"id": null, "schema_name": "qcschema_input", "schema_version": 1, "molecule": {"schema_name": "qcschema_molecule", "schema_version": 2, "validated": true, "symbols": ["N", "H", "H", "H"], "geometry": [-1e-10, -0.1040380466, 0.0, -0.9015844116, 0.4818470201, -1.5615900098, -0.9015844116, 0.4818470201, 1.5615900098, 1.8031688251, 0.4818470204, 0.0], "name": "H3N", "molecular_charge": 0.0, "molecular_multiplicity": 1, "masses": [14.00307400443, 1.00782503223, 1.00782503223, 1.00782503223], "real": [true, true, true, true], "atom_labels": ["", "", "", ""], "atomic_numbers": [7, 1, 1, 1], "mass_numbers": [14, 1, 1, 1], "fragments": [[0, 1, 2, 3]], "fragment_charges": [0.0], "fragment_multiplicities": [1], "fix_com": false, "fix_orientation": false, "fix_symmetry": "c1", "provenance": {"creator": "QCElemental", "version": "v0.17.0+7.gf55d5ac.dirty", "routine": "qcelemental.molparse.from_string"}}, "driver": "energy", "model": {"method": "scf", "basis": "STO-3G"}, "keywords": {"pcm": 1, "pcm__input": "\\n       Units = Angstrom\\n       Medium {\\n       SolverType = IEFPCM\\n       Solvent = Water\\n       }\\n    \\n       Cavity {\\n       RadiiSet = UFF\\n       Type = GePol\\n       Scaling = False\\n       Area = 0.3\\n       Mode = Implicit\\n       }\\n    ", "pcm_scf_type": "SEPARATE", "scf_type": "PK"}, "protocols": {}, "extras": {"wfn_qcvars_only": true}, "provenance": {"creator": "Psi4", "version": "1.4a2.dev1089", "routine": "psi4.driver.p4util.procutil"}}"""

    atres2 = psi4.schema_wrapper.run_qcschema(json.loads(jatin))
    pprint.pprint(atres2.dict())
    assert psi4.compare_values(totalenergy, atres2.return_result, 10, "Total energy (PCM, separate algorithm)")
    assert psi4.compare_values(polenergy, atres2.extras["qcvars"]["PCM POLARIZATION ENERGY"], 6,
                               "Polarization energy (PCM, separate algorithm)")

    jatin = """{"id": null, "schema_name": "qcschema_input", "schema_version": 1, "molecule": {"schema_name": "qcschema_molecule", "schema_version": 2, "validated": true, "symbols": ["N", "H", "H", "H"], "geometry": [-1e-10, -0.1040380466, 0.0, -0.9015844116, 0.4818470201, -1.5615900098, -0.9015844116, 0.4818470201, 1.5615900098, 1.8031688251, 0.4818470204, 0.0], "name": "H3N", "molecular_charge": 0.0, "molecular_multiplicity": 1, "masses": [14.00307400443, 1.00782503223, 1.00782503223, 1.00782503223], "real": [true, true, true, true], "atom_labels": ["", "", "", ""], "atomic_numbers": [7, 1, 1, 1], "mass_numbers": [14, 1, 1, 1], "fragments": [[0, 1, 2, 3]], "fragment_charges": [0.0], "fragment_multiplicities": [1], "fix_com": false, "fix_orientation": false, "fix_symmetry": "c1", "provenance": {"creator": "QCElemental", "version": "v0.17.0+7.gf55d5ac.dirty", "routine": "qcelemental.molparse.from_string"}}, "driver": "energy", "model": {"method": "scf", "basis": "STO-3G"}, "keywords": {"pcm": 1, "pcm__input": "\\n       Units = Angstrom\\n       Medium {\\n       SolverType = IEFPCM\\n       Solvent = Water\\n       }\\n    \\n       Cavity {\\n       RadiiSet = UFF\\n       Type = GePol\\n       Scaling = False\\n       Area = 0.3\\n       Mode = Implicit\\n       }\\n    ", "pcm_scf_type": "TOTAL", "scf_type": "PK", "reference": "uhf"}, "protocols": {}, "extras": {"wfn_qcvars_only": true}, "provenance": {"creator": "Psi4", "version": "1.4a2.dev1089", "routine": "psi4.driver.p4util.procutil"}}"""

    atres3 = psi4.schema_wrapper.run_qcschema(json.loads(jatin))
    pprint.pprint(atres3.dict())
    assert psi4.compare_values(totalenergy, atres3.return_result, 10, "Total energy (PCM, separate algorithm)")
    assert psi4.compare_values(polenergy, atres3.extras["qcvars"]["PCM POLARIZATION ENERGY"], 6,
                               "Polarization energy (PCM, separate algorithm)")

    jatin = """{"id": null, "schema_name": "qcschema_input", "schema_version": 1, "molecule": {"schema_name": "qcschema_molecule", "schema_version": 2, "validated": true, "symbols": ["N", "H", "H", "H"], "geometry": [-1e-10, -0.1040380466, 0.0, -0.9015844116, 0.4818470201, -1.5615900098, -0.9015844116, 0.4818470201, 1.5615900098, 1.8031688251, 0.4818470204, 0.0], "name": "H3N", "molecular_charge": 0.0, "molecular_multiplicity": 1, "masses": [14.00307400443, 1.00782503223, 1.00782503223, 1.00782503223], "real": [true, true, true, true], "atom_labels": ["", "", "", ""], "atomic_numbers": [7, 1, 1, 1], "mass_numbers": [14, 1, 1, 1], "fragments": [[0, 1, 2, 3]], "fragment_charges": [0.0], "fragment_multiplicities": [1], "fix_com": false, "fix_orientation": false, "fix_symmetry": "c1", "provenance": {"creator": "QCElemental", "version": "v0.17.0+7.gf55d5ac.dirty", "routine": "qcelemental.molparse.from_string"}}, "driver": "energy", "model": {"method": "scf", "basis": "STO-3G"}, "keywords": {"pcm": 1, "pcm__input": "\\n       Units = Angstrom\\n       Medium {\\n       SolverType = IEFPCM\\n       Solvent = Water\\n       }\\n    \\n       Cavity {\\n       RadiiSet = UFF\\n       Type = GePol\\n       Scaling = False\\n       Area = 0.3\\n       Mode = Implicit\\n       }\\n    ", "pcm_scf_type": "TOTAL", "scf_type": "PK", "reference": "rohf"}, "protocols": {}, "extras": {"wfn_qcvars_only": true}, "provenance": {"creator": "Psi4", "version": "1.4a2.dev1089", "routine": "psi4.driver.p4util.procutil"}}"""

    atres4 = psi4.schema_wrapper.run_qcschema(json.loads(jatin))
    pprint.pprint(atres4.dict())
    assert psi4.compare_values(totalenergy, atres4.return_result, 10, "Total energy (PCM, separate algorithm)")
    assert psi4.compare_values(polenergy, atres4.extras["qcvars"]["PCM POLARIZATION ENERGY"], 6,
                               "Polarization energy (PCM, separate algorithm)")


@pytest.mark.parametrize("integral_package", [
    pytest.param("libint2"),
    pytest.param("simint", marks=using("simint")),
])
def test_integrals(integral_package):
    """scf5"""
    #! Test of all different algorithms and reference types for SCF, on singlet and triplet O2, using the cc-pVTZ basis set and using ERD integrals.

    Eref_sing_can = -149.58723684929720
    Eref_sing_df = -149.58715054487624
    Eref_uhf_can = -149.67135517240553
    Eref_uhf_df = -149.67125624291961
    Eref_rohf_can = -149.65170765757173
    Eref_rohf_df = -149.65160796208073

    jatin = """{"id": null, "schema_name": "qcschema_input", "schema_version": 1, "molecule": {"schema_name": "qcschema_molecule", "schema_version": 2, "validated": true, "symbols": ["O", "O"], "geometry": [0.0, 0.0, -1.0393493690018054, 0.0, 0.0, 1.0393493690018054], "name": "O2", "molecular_charge": 0.0, "molecular_multiplicity": 1, "masses": [15.99491461957, 15.99491461957], "real": [true, true], "atom_labels": ["", ""], "atomic_numbers": [8, 8], "mass_numbers": [16, 16], "fragments": [[0, 1]], "fragment_charges": [0.0], "fragment_multiplicities": [1], "fix_com": false, "fix_orientation": false, "provenance": {"creator": "QCElemental", "version": "v0.17.0+7.gf55d5ac.dirty", "routine": "qcelemental.molparse.from_string"}}, "driver": "energy", "model": {"method": "scf", "basis": "CC-PVTZ"}, "keywords": {"df_basis_scf": "CC-PVTZ-JKFIT", "print": 2, "scf__scf_type": "DIRECT"}, "protocols": {}, "extras": {"wfn_qcvars_only": true}, "provenance": {"creator": "Psi4", "version": "1.4a2.dev1089", "routine": "psi4.driver.p4util.procutil"}}"""

    datin = json.loads(jatin)
    datin["keywords"]["integral_package"] = integral_package
    atres = psi4.schema_wrapper.run_qcschema(datin)
    assert psi4.compare_values(Eref_sing_can, atres.return_result, 6, 'Singlet Direct RHF energy')

    jatin = """{"id": null, "schema_name": "qcschema_input", "schema_version": 1, "molecule": {"schema_name": "qcschema_molecule", "schema_version": 2, "validated": true, "symbols": ["O", "O"], "geometry": [0.0, 0.0, -1.0393493690018054, 0.0, 0.0, 1.0393493690018054], "name": "O2", "molecular_charge": 0.0, "molecular_multiplicity": 1, "masses": [15.99491461957, 15.99491461957], "real": [true, true], "atom_labels": ["", ""], "atomic_numbers": [8, 8], "mass_numbers": [16, 16], "fragments": [[0, 1]], "fragment_charges": [0.0], "fragment_multiplicities": [1], "fix_com": false, "fix_orientation": false, "provenance": {"creator": "QCElemental", "version": "v0.17.0+7.gf55d5ac.dirty", "routine": "qcelemental.molparse.from_string"}}, "driver": "energy", "model": {"method": "scf", "basis": "CC-PVTZ"}, "keywords": {"df_basis_scf": "CC-PVTZ-JKFIT", "print": 2, "scf_type": "df", "reference": "uhf"}, "protocols": {}, "extras": {"wfn_qcvars_only": true}, "provenance": {"creator": "Psi4", "version": "1.4a2.dev1089", "routine": "psi4.driver.p4util.procutil"}}"""

    datin = json.loads(jatin)
    datin["keywords"]["integral_package"] = integral_package
    atres = psi4.schema_wrapper.run_qcschema(datin)
    assert psi4.compare_values(Eref_sing_df, atres.return_result, 6, 'Singlet DF UHF energy')

    jatin = """{"id": null, "schema_name": "qcschema_input", "schema_version": 1, "molecule": {"schema_name": "qcschema_molecule", "schema_version": 2, "validated": true, "symbols": ["O", "O"], "geometry": [0.0, 0.0, -1.0393493690018054, 0.0, 0.0, 1.0393493690018054], "name": "O2", "molecular_charge": 0.0, "molecular_multiplicity": 3, "masses": [15.99491461957, 15.99491461957], "real": [true, true], "atom_labels": ["", ""], "atomic_numbers": [8, 8], "mass_numbers": [16, 16], "fragments": [[0, 1]], "fragment_charges": [0.0], "fragment_multiplicities": [1], "fix_com": false, "fix_orientation": false, "provenance": {"creator": "QCElemental", "version": "v0.17.0+7.gf55d5ac.dirty", "routine": "qcelemental.molparse.from_string"}}, "driver": "energy", "model": {"method": "scf", "basis": "CC-PVTZ"}, "keywords": {"df_basis_scf": "CC-PVTZ-JKFIT", "print": 2, "scf_type": "out_of_core", "reference": "rohf"}, "protocols": {}, "extras": {"wfn_qcvars_only": true}, "provenance": {"creator": "Psi4", "version": "1.4a2.dev1089", "routine": "psi4.driver.p4util.procutil"}}"""

    datin = json.loads(jatin)
    datin["keywords"]["integral_package"] = integral_package
    atres = psi4.schema_wrapper.run_qcschema(datin)
    assert psi4.compare_values(Eref_rohf_can, atres.return_result, 6, 'Triplet Disk ROHF energy')


def test_json():
    """json/energy"""

    import numpy as np

    # Generate JSON data
    json_input = {
        "schema_name": "qcschema_input",
        "schema_version": 1,
        "molecule": {
            "symbols": ["He", "He"],
            "geometry": [0, 0, -1, 0, 0, 1]
        },
        "driver": "gradient",
        "model": {
            "method": "SCF",
            "basis": "sto-3g"
        },
        "keywords": {}
    }

    json_ret = psi4.schema_wrapper.run_qcschema(json_input)
    json_ret = json_ret.dict()
    pprint.pprint(json_ret)

    assert psi4.compare_integers(True, json_ret["success"], "Success")
    assert psi4.compare_values(-5.474227786274896, json_ret["properties"]["return_energy"], 4, "SCF ENERGY")

    bench_gradient = np.array([[0.0, 0.0, 0.32746933], [0.0, 0.0, -0.32746933]])
    cgradient = np.array(json_ret["return_result"]).reshape(-1, 3)
    assert psi4.compare_arrays(bench_gradient, cgradient, 4, "SCF RETURN GRADIENT")


#@pytest.mark.smoke
#@uusing("cfour")
#def test_cfour():
#    """cfour/sp-rhf-ccsd_t_"""
#    #! single-point CCSD(T)/qz2p on water
#
#    print('        <<< Translation of ZMAT to Psi4 format to Cfour >>>')
#
#    psi4.geometry("""
#    O
#    H 1 R
#    H 1 R 2 A
#
#    R=0.958
#    A=104.5
#    """)
#
#    psi4.set_options({
#    'cfour_CALC_level': 'CCSD(T)',
#    'cfour_BASIS': 'qz2p',
#    'cfour_SCF_CONV': 12,
#    'cfour_CC_CONV': 12,
#    })
#
#    psi4.energy('cfour')
#
##    atin = psi4.driver.p4util.state_to_atomicinput(driver="energy", method="ccsd", molecule=ethene_ethyne)
##    print(f'    jatin = """{atin.serialize("json")}"""')
##    assert 0
##
##    atres = psi4.schema_wrapper.run_qcschema(json.loads(jatin))
##    pprint.pprint(atres.dict())
#
#
#    assert psi4.compare_values(-76.062748460117, psi4.variable('scf total energy'), 6, 'SCF')
#    assert psi4.compare_values(-76.332940127333, psi4.variable('mp2 total energy'), 6, 'MP2')
#    assert psi4.compare_values(-76.338453951890, psi4.variable('ccsd total energy'), 6, 'CCSD')
#    assert psi4.compare_values(-0.275705491773, psi4.variable('ccsd correlation energy'), 6, 'CCSD corl')
#    assert psi4.compare_values(-76.345717549886, psi4.variable('ccsd(t) total energy'), 6, 'CCSD(T)')
#    assert psi4.compare_values(-0.282969089769, psi4.variable('ccsd(t) correlation energy'), 6, 'CCSD(T) corl')


@uusing("v2rdm_casscf")
def test_v2rdm_casscf():
    """v2rdm_casscf/tests/v2rdm1"""
    #! cc-pvdz N2 (6,6) active space Test DQG

    print('        N2 / cc-pVDZ / DQG(6,6), scf_type = CD / 1e-12, rNN = 0.5 A')

    interloper = psi4.geometry("""
    0 1
    O
    H 1 1.0
    H 1 1.0 2 90.0
    """)

    # NOTES
    # * add plugin keywords to AtomicInput by hand, not set_module_options, since the module list doesn't know about v2rdm

    jatin = """{"id": null, "schema_name": "qcschema_input", "schema_version": 1, "molecule": {"schema_name": "qcschema_molecule", "schema_version": 2, "validated": true, "symbols": ["N", "N"], "geometry": [0.0, 0.0, -0.4724315332214108, 0.0, 0.0, 0.4724315332214108], "name": "N2", "molecular_charge": 0.0, "molecular_multiplicity": 1, "masses": [14.00307400443, 14.00307400443], "real": [true, true], "atom_labels": ["", ""], "atomic_numbers": [7, 7], "mass_numbers": [14, 14], "fragments": [[0, 1]], "fragment_charges": [0.0], "fragment_multiplicities": [1], "fix_com": false, "fix_orientation": false, "provenance": {"creator": "QCElemental", "version": "v0.17.0+7.gf55d5ac.dirty", "routine": "qcelemental.molparse.from_string"}}, "driver": "energy", "model": {"method": "v2rdm-casscf", "basis": "CC-PVDZ"}, "keywords": {"active": [1, 0, 1, 1, 0, 1, 1, 1], "cholesky_tolerance": 1e-12, "d_convergence": 1e-10, "maxiter": 500, "restricted_docc": [2, 0, 0, 0, 0, 2, 0, 0], "scf_type": "CD", "v2rdm_casscf__r_convergence": 1e-5, "v2rdm_casscf__e_convergence": 1e-6, "v2rdm_casscf__maxiter": 20000}, "protocols": {}, "extras": {"wfn_qcvars_only": true}, "provenance": {"creator": "Psi4", "version": "1.4a2.dev1089", "routine": "psi4.driver.p4util.procutil"}}"""

    atres = psi4.schema_wrapper.run_qcschema(json.loads(jatin))
    pprint.pprint(atres.dict())

    assert psi4.compare_values(-103.04337420425350, atres.extras["qcvars"]["SCF TOTAL ENERGY"], 8, "SCF total energy")
    assert psi4.compare_values(-103.086205379481, atres.return_result, 5, "v2RDM-CASSCF total energy")


#@hardware_nvidia_gpu
#@uusing("gpu_dfcc")
#def test_gpu_dfcc():
#    """gpu_dfcc/tests/gpu_dfcc1"""
#    #! cc-pvdz (H2O)2 Test DF-CCSD vs GPU-DF-CCSD
#
#    import gpu_dfcc
#
#    H20 = psi4.geometry("""
#               O          0.000000000000     0.000000000000    -0.068516219310
#               H          0.000000000000    -0.790689573744     0.543701060724
#               H          0.000000000000     0.790689573744     0.543701060724
#    """)
#
#    psi4.set_memory(32000000000)
#    psi4.set_options({
#      'cc_timings': False,
#      'num_gpus': 1,
#      'cc_type': 'df',
#      'df_basis_cc':  'aug-cc-pvdz-ri',
#      'df_basis_scf': 'aug-cc-pvdz-jkfit',
#      'basis':        'aug-cc-pvdz',
#      'freeze_core': 'true',
#      'e_convergence': 1e-8,
#      'd_convergence': 1e-8,
#      'r_convergence': 1e-8,
#      'scf_type': 'df',
#      'maxiter': 30})
#    psi4.set_num_threads(2)
#    en_dfcc     = psi4.energy('ccsd', molecule=H20)
#    en_gpu_dfcc = psi4.energy('gpu-df-ccsd', molecule=H20)
#
##    atin = psi4.driver.p4util.state_to_atomicinput(driver="energy", method="ccsd", molecule=ethene_ethyne)
##    print(f'    jatin = """{atin.serialize("json")}"""')
##    assert 0
##
##    atres = psi4.schema_wrapper.run_qcschema(json.loads(jatin))
##    pprint.pprint(atres.dict())
#
##    atin = psi4.driver.p4util.state_to_atomicinput(driver="energy", method="ccsd", molecule=ethene_ethyne)
##    print(f'    jatin = """{atin.serialize("json")}"""')
##    assert 0
##
##    atres = psi4.schema_wrapper.run_qcschema(json.loads(jatin))
##    pprint.pprint(atres.dict())
#
#
#    assert psi4.compare_values(en_gpu_dfcc, en_dfcc, 8, "CCSD total energy")


@pytest.mark.nbody
@uusing("dftd3")
@uusing("gcp")
def test_grimme_3c():

    # NOTES
    # * add `model.basis = "(auto)"`

    jatin = """{"id": null, "schema_name": "qcschema_input", "schema_version": 1, "molecule": {"schema_name": "qcschema_molecule", "schema_version": 2, "validated": true, "symbols": ["C", "C", "H", "H", "H", "H", "C", "C", "H", "H"], "geometry": [2.676440890378796e-18, -1.261539587380886, -4.02163254721969, 2.676440890378796e-18, 1.261539587380886, -4.02163254721969, 1.7453907337214847, -2.3286206872737854, -4.0245162692871395, -1.7453907337214847, -2.3286206872737854, -4.0245162692871395, -1.7453907337214847, 2.3286206872737854, -4.0245162692871395, 1.7453907337214847, 2.3286206872737854, -4.0245162692871395, 2.676440890378796e-18, -1.8740279466317074e-17, 5.47454736883822, 2.676440890378796e-18, -1.8740279466317074e-17, 3.193150937439626, 2.676440890378796e-18, -1.8740279466317074e-17, 1.1789145370276326, 2.676440890378796e-18, -1.8740279466317074e-17, 7.484131263529336], "name": "C4H6", "molecular_charge": 0.0, "molecular_multiplicity": 1, "masses": [12.0, 12.0, 1.00782503223, 1.00782503223, 1.00782503223, 1.00782503223, 12.0, 12.0, 1.00782503223, 1.00782503223], "real": [true, true, true, true, true, true, true, true, true, true], "atom_labels": ["", "", "", "", "", "", "", "", "", ""], "atomic_numbers": [6, 6, 1, 1, 1, 1, 6, 6, 1, 1], "mass_numbers": [12, 12, 1, 1, 1, 1, 12, 12, 1, 1], "fragments": [[0, 1, 2, 3, 4, 5], [6, 7, 8, 9]], "fragment_charges": [0.0, 0.0], "fragment_multiplicities": [1, 1], "fix_com": false, "fix_orientation": false, "fix_symmetry": "c1", "provenance": {"creator": "QCElemental", "version": "v0.17.0+7.gf55d5ac.dirty", "routine": "qcelemental.molparse.from_string"}}, "driver": "energy", "model": {"method": "pbeh3c", "basis": "(auto)"}, "keywords": {"function_kwargs": {"bsse_type": "nocp"}}, "protocols": {}, "extras": {"wfn_qcvars_only": true}, "provenance": {"creator": "Psi4", "version": "1.4a2.dev1089", "routine": "psi4.driver.p4util.procutil"}}"""

    atres = psi4.schema_wrapper.run_qcschema(json.loads(jatin))
    assert psi4.compare_values(-2.153, atres.return_result * psi4.constants.hartree2kcalmol, 0.03,
                               'S22-16 PBEh-3c/def2-mSVP')

    jatin = """{"id": null, "schema_name": "qcschema_input", "schema_version": 1, "molecule": {"schema_name": "qcschema_molecule", "schema_version": 2, "validated": true, "symbols": ["C", "C", "H", "H", "H", "H", "C", "C", "H", "H"], "geometry": [2.676440890378796e-18, -1.261539587380886, -4.02163254721969, 2.676440890378796e-18, 1.261539587380886, -4.02163254721969, 1.7453907337214847, -2.3286206872737854, -4.0245162692871395, -1.7453907337214847, -2.3286206872737854, -4.0245162692871395, -1.7453907337214847, 2.3286206872737854, -4.0245162692871395, 1.7453907337214847, 2.3286206872737854, -4.0245162692871395, 2.676440890378796e-18, -1.8740279466317074e-17, 5.47454736883822, 2.676440890378796e-18, -1.8740279466317074e-17, 3.193150937439626, 2.676440890378796e-18, -1.8740279466317074e-17, 1.1789145370276326, 2.676440890378796e-18, -1.8740279466317074e-17, 7.484131263529336], "name": "C4H6", "molecular_charge": 0.0, "molecular_multiplicity": 1, "masses": [12.0, 12.0, 1.00782503223, 1.00782503223, 1.00782503223, 1.00782503223, 12.0, 12.0, 1.00782503223, 1.00782503223], "real": [true, true, true, true, true, true, true, true, true, true], "atom_labels": ["", "", "", "", "", "", "", "", "", ""], "atomic_numbers": [6, 6, 1, 1, 1, 1, 6, 6, 1, 1], "mass_numbers": [12, 12, 1, 1, 1, 1, 12, 12, 1, 1], "fragments": [[0, 1, 2, 3, 4, 5], [6, 7, 8, 9]], "fragment_charges": [0.0, 0.0], "fragment_multiplicities": [1, 1], "fix_com": false, "fix_orientation": false, "fix_symmetry": "c1", "provenance": {"creator": "QCElemental", "version": "v0.17.0+7.gf55d5ac.dirty", "routine": "qcelemental.molparse.from_string"}}, "driver": "energy", "model": {"method": "hf3c", "basis": ""}, "keywords": {"scf_type": "pk", "function_kwargs": {"bsse_type": "nocp"}}, "protocols": {}, "extras": {"wfn_qcvars_only": true}, "provenance": {"creator": "Psi4", "version": "1.4a2.dev1089", "routine": "psi4.driver.p4util.procutil"}}"""

    atres = psi4.schema_wrapper.run_qcschema(json.loads(jatin))
    assert psi4.compare_values(-0.00240232, atres.return_result, 6, 'S22-16 HF-3c/minix')


@uusing("dkh")
def test_dkh():
    """dkh/molpro-2order"""

    jatin = """{"id": null, "schema_name": "qcschema_input", "schema_version": 1, "molecule": {"schema_name": "qcschema_molecule", "schema_version": 2, "validated": true, "symbols": ["Ne"], "geometry": [0.0, 0.0, 0.0], "name": "Ne", "molecular_charge": 0.0, "molecular_multiplicity": 1, "masses": [19.9924401762], "real": [true], "atom_labels": [""], "atomic_numbers": [10], "mass_numbers": [20], "fragments": [[0]], "fragment_charges": [0.0], "fragment_multiplicities": [1], "fix_com": false, "fix_orientation": false, "provenance": {"creator": "QCElemental", "version": "v0.17.0+7.gf55d5ac.dirty", "routine": "qcelemental.molparse.from_string"}}, "driver": "energy", "model": {"method": "scf", "basis": "CC-PVTZ-DK"}, "keywords": {"dkh_order": 2, "print": 2, "reference": "RHF", "relativistic": "DKH", "scf_type": "PK"}, "protocols": {}, "extras": {"wfn_qcvars_only": true}, "provenance": {"creator": "Psi4", "version": "1.4a2.dev1089", "routine": "psi4.driver.p4util.procutil"}}"""

    atres = psi4.schema_wrapper.run_qcschema(json.loads(jatin))

    assert psi4.compare_values(-128.66891610, atres.return_result, 6, '2nd order vs Molpro')


#@uusing("ambit")
#@uusing("forte")
#def disabled_test_forte():
#    """aci-10: Perform aci on benzyne"""
#
#    import forte
#
#    refscf    = -229.20378006852584
#    refaci    = -229.359450812283
#    refacipt2 = -229.360444943286
#
#    mbenzyne = psi4.geometry("""
#      0 1
#       C   0.0000000000  -2.5451795941   0.0000000000
#       C   0.0000000000   2.5451795941   0.0000000000
#       C  -2.2828001669  -1.3508352528   0.0000000000
#       C   2.2828001669  -1.3508352528   0.0000000000
#       C   2.2828001669   1.3508352528   0.0000000000
#       C  -2.2828001669   1.3508352528   0.0000000000
#       H  -4.0782187459  -2.3208602146   0.0000000000
#       H   4.0782187459  -2.3208602146   0.0000000000
#       H   4.0782187459   2.3208602146   0.0000000000
#       H  -4.0782187459   2.3208602146   0.0000000000
#
#      units bohr
#    """)
#
#    psi4.set_options({
#       'basis': 'DZ',
#       'df_basis_mp2': 'cc-pvdz-ri',
#       'reference': 'uhf',
#       'scf_type': 'pk',
#       'd_convergence': 10,
#       'e_convergence': 12,
#       'guess': 'gwh',
#    })
#
#    psi4.set_module_options("FORTE", {
#      'root_sym': 0,
#      'frozen_docc':     [2,1,0,0,0,0,2,1],
#      'restricted_docc': [3,2,0,0,0,0,2,3],
#      'active':          [1,0,1,2,1,2,1,0],
#      'multiplicity': 1,
#      'aci_nroot': 1,
#      'job_type': 'aci',
#      'sigma': 0.001,
#      'aci_select_type': 'aimed_energy',
#      'aci_spin_projection': 1,
#      'aci_enforce_spin_complete': True,
#      'aci_add_aimed_degenerate': False,
#      'aci_project_out_spin_contaminants': False,
#      'diag_algorithm': 'full',
#      'aci_quiet_mode': True,
#    })
#
##    atin = psi4.driver.p4util.state_to_atomicinput(driver="energy", method="ccsd", molecule=ethene_ethyne)
##    print(f'    jatin = """{atin.serialize("json")}"""')
##    assert 0
##
##    atres = psi4.schema_wrapper.run_qcschema(json.loads(jatin))
##    pprint.pprint(atres.dict())
#
#
#    scf = psi4.energy('scf')
#    assert psi4.compare_values(refscf, scf,10,"SCF Energy")
#
#    psi4.energy('forte')
#    assert psi4.compare_values(refaci, psi4.variable("ACI ENERGY"),10,"ACI energy")
#    assert psi4.compare_values(refacipt2, psi4.variable("ACI+PT2 ENERGY"),8,"ACI+PT2 energy")


@uusing("snsmp2")
def test_snsmp2():
    """snsmp2/he-he"""

    # NOTES
    # * add `model.basis = "(auto)"`

    jatin = """{"id": null, "schema_name": "qcschema_input", "schema_version": 1, "molecule": {"schema_name": "qcschema_molecule", "schema_version": 2, "validated": true, "symbols": ["He", "He"], "geometry": [-1.8897261254578286, -2.892808813508824e-17, 0.0, 1.8897261254578286, 2.892808813508824e-17, 0.0], "name": "He2", "molecular_charge": 0.0, "molecular_multiplicity": 1, "masses": [4.00260325413, 4.00260325413], "real": [true, true], "atom_labels": ["", ""], "atomic_numbers": [2, 2], "mass_numbers": [4, 4], "fragments": [[0], [1]], "fragment_charges": [0.0, 0.0], "fragment_multiplicities": [1, 1], "fix_com": false, "fix_orientation": false, "provenance": {"creator": "QCElemental", "version": "v0.17.0+7.gf55d5ac.dirty", "routine": "qcelemental.molparse.from_string"}}, "driver": "energy", "model": {"method": "sns-mp2", "basis": "(auto)"}, "keywords": {}, "protocols": {}, "extras": {"wfn_qcvars_only": true}, "provenance": {"creator": "Psi4", "version": "1.4a2.dev1089", "routine": "psi4.driver.p4util.procutil"}}"""

    atres = psi4.schema_wrapper.run_qcschema(json.loads(jatin))
    pprint.pprint(atres.dict())

    assert psi4.compare_values(0.00176708227, atres.return_result, 5, "SNS-MP2 IE [Eh]")


@uusing("resp")
def test_resp():
    """resp/tests/test_resp_1"""

    pytest.xfail("RESP calls Psi4, not the reverse, so no AtomicInput possible")


@uusing("fockci")
def test_psi4fockci():
    """psi4fockci/n2"""

    jatin = """{"id": null, "schema_name": "qcschema_input", "schema_version": 1, "molecule": {"schema_name": "qcschema_molecule", "schema_version": 2, "validated": true, "symbols": ["N", "N"], "geometry": [0.0, 0.0, -2.3621576568222853, 0.0, 0.0, 2.3621576568222853], "name": "N2", "molecular_charge": 0.0, "molecular_multiplicity": 7, "masses": [14.00307400443, 14.00307400443], "real": [true, true], "atom_labels": ["", ""], "atomic_numbers": [7, 7], "mass_numbers": [14, 14], "fragments": [[0, 1]], "fragment_charges": [0.0], "fragment_multiplicities": [7], "fix_com": false, "fix_orientation": false, "fix_symmetry": "c1", "provenance": {"creator": "QCElemental", "version": "v0.17.0+7.gf55d5ac.dirty", "routine": "qcelemental.molparse.from_string"}}, "driver": "energy", "model": {"method": "psi4fockci", "basis": "cc-pvdz"}, "keywords": {"function_kwargs": {"new_charge": 0, "new_multiplicity": 1}}, "protocols": {}, "extras": {"wfn_qcvars_only": true}, "provenance": {"creator": "Psi4", "version": "1.4a2.dev1089", "routine": "psi4.driver.p4util.procutil"}}"""

    atres = psi4.schema_wrapper.run_qcschema(json.loads(jatin))
    assert psi4.compare_values(-108.776024394853295, atres.extras["qcvars"]["CI ROOT 0 TOTAL ENERGY"], 7, "3SF Energy")

    jatin = """{"id": null, "schema_name": "qcschema_input", "schema_version": 1, "molecule": {"schema_name": "qcschema_molecule", "schema_version": 2, "validated": true, "symbols": ["N", "N"], "geometry": [0.0, 0.0, -2.3621576568222853, 0.0, 0.0, 2.3621576568222853], "name": "N2", "molecular_charge": 0.0, "molecular_multiplicity": 7, "masses": [14.00307400443, 14.00307400443], "real": [true, true], "atom_labels": ["", ""], "atomic_numbers": [7, 7], "mass_numbers": [14, 14], "fragments": [[0, 1]], "fragment_charges": [0.0], "fragment_multiplicities": [7], "fix_com": false, "fix_orientation": false, "fix_symmetry": "c1", "provenance": {"creator": "QCElemental", "version": "v0.17.0+7.gf55d5ac.dirty", "routine": "qcelemental.molparse.from_string"}}, "driver": "energy", "model": {"method": "psi4fockci", "basis": "cc-pvdz"}, "keywords": {"function_kwargs": {"new_charge": 1, "new_multiplicity": 2}}, "protocols": {}, "extras": {"wfn_qcvars_only": true}, "provenance": {"creator": "Psi4", "version": "1.4a2.dev1089", "routine": "psi4.driver.p4util.procutil"}}"""

    atres = psi4.schema_wrapper.run_qcschema(json.loads(jatin))
    assert psi4.compare_values(-108.250639579451, atres.extras["qcvars"]["CI ROOT 0 TOTAL ENERGY"], 7, "2SF-IP Energy")

    jatin = """{"id": null, "schema_name": "qcschema_input", "schema_version": 1, "molecule": {"schema_name": "qcschema_molecule", "schema_version": 2, "validated": true, "symbols": ["N", "N"], "geometry": [0.0, 0.0, -2.3621576568222853, 0.0, 0.0, 2.3621576568222853], "name": "N2", "molecular_charge": 0.0, "molecular_multiplicity": 7, "masses": [14.00307400443, 14.00307400443], "real": [true, true], "atom_labels": ["", ""], "atomic_numbers": [7, 7], "mass_numbers": [14, 14], "fragments": [[0, 1]], "fragment_charges": [0.0], "fragment_multiplicities": [7], "fix_com": false, "fix_orientation": false, "fix_symmetry": "c1", "provenance": {"creator": "QCElemental", "version": "v0.17.0+7.gf55d5ac.dirty", "routine": "qcelemental.molparse.from_string"}}, "driver": "energy", "model": {"method": "psi4fockci", "basis": "cc-pvdz"}, "keywords": {"function_kwargs": {"new_charge": -1, "new_multiplicity": 2}}, "protocols": {}, "extras": {"wfn_qcvars_only": true}, "provenance": {"creator": "Psi4", "version": "1.4a2.dev1089", "routine": "psi4.driver.p4util.procutil"}}"""

    atres = psi4.schema_wrapper.run_qcschema(json.loads(jatin))
    assert psi4.compare_values(-108.600832070267, atres.extras["qcvars"]["CI ROOT 0 TOTAL ENERGY"], 7, "2SF-EA Energy")


@uusing("cppe")
def test_cppe():
    #! PE-SCF of PNA in presence of 6 water molecules
    #! Reference data from Q-Chem calculation

    #    potfile = """
    #! Generated by PyFraME 0.1.0
    #@COORDINATES
    #18
    #AA
    #O      9.37100000     2.95300000    -6.07800000       1
    #H      8.87200000     2.13400000    -6.04900000       2
    #...
    #@MULTIPOLES
    #ORDER 0
    #...
    #ORDER 1
    #...
    #ORDER 2
    #...
    #@POLARIZABILITIES
    #ORDER 1 1
    #...
    #EXCLISTS
    #...
    #17   16  18
    #18   16  17
    #"""
    #    potfile = "\\n".join(potfile.split("\n"))
    #
    #    psi4.set_options({
    #     'pe': True,
    #     ...
    #     'pe__potfile': potfile,
    #    })

    jatin = """{"id": null, "schema_name": "qcschema_input", "schema_version": 1, "molecule": {"schema_name": "qcschema_molecule", "schema_version": 2, "validated": true, "symbols": ["C", "C", "C", "C", "C", "C", "H", "H", "H", "H", "N", "N", "O", "O", "H", "H"], "geometry": [16.3423515329593, 2.031455584867165, -3.233321400658344, 17.918383121591127, 0.8125822339468661, -1.5268987093699253, 17.755866674801755, 1.417294594093371, 1.0166726554963117, 16.028656996133297, 3.2352111267838017, 1.880277494830539, 14.462074038128758, 4.431407764198608, 0.10393493690018055, 14.611362402039928, 3.838033760804849, -2.441526154091514, 19.243081135537064, -0.5839253727664689, -2.1996412100329117, 18.980409204098425, 0.4762109836153727, 2.356488478445912, 13.118478762928243, 5.837364001539231, 0.7351034628030951, 13.411386312374207, 4.752661205526438, -3.813467321173897, 15.875589179971215, 3.826695404052102, 4.393613241689451, 16.50486797974867, 1.4002870589642507, -5.912953046557544, 15.08001448115347, 2.515225472984369, -7.371821615410987, 18.058222854875005, -0.2078698738003611, -6.5497907508368325, 14.64348774617271, 5.123047526116172, 5.01155368471416, 16.99052759399133, 2.9763186475960794, 5.659729745746196], "name": "C6H6N2O2", "molecular_charge": 0.0, "molecular_multiplicity": 1, "masses": [12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 1.00782503223, 1.00782503223, 1.00782503223, 1.00782503223, 14.00307400443, 14.00307400443, 15.99491461957, 15.99491461957, 1.00782503223, 1.00782503223], "real": [true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true], "atom_labels": ["", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""], "atomic_numbers": [6, 6, 6, 6, 6, 6, 1, 1, 1, 1, 7, 7, 8, 8, 1, 1], "mass_numbers": [12, 12, 12, 12, 12, 12, 1, 1, 1, 1, 14, 14, 16, 16, 1, 1], "fragments": [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]], "fragment_charges": [0.0], "fragment_multiplicities": [1], "fix_com": false, "fix_orientation": false, "provenance": {"creator": "QCElemental", "version": "v0.17.0+7.gf55d5ac.dirty", "routine": "qcelemental.molparse.from_string"}}, "driver": "energy", "model": {"method": "scf", "basis": "STO-3G"}, "keywords": {"d_convergence": 1e-10, "e_convergence": 1e-10, "pe": 1, "pe__potfile": "\\n! Generated by PyFraME 0.1.0\\n@COORDINATES\\n18\\nAA\\nO      9.37100000     2.95300000    -6.07800000       1\\nH      8.87200000     2.13400000    -6.04900000       2\\nH      9.87300000     2.94000000    -5.26100000       3\\nO      7.72000000     5.12000000    -5.51900000       4\\nH      7.63800000     5.70400000    -6.27600000       5\\nH      8.29100000     4.41700000    -5.83600000       6\\nO     10.45300000     3.07700000    -3.43400000       7\\nH      9.94500000     3.80300000    -3.06600000       8\\nH     11.35900000     3.30000000    -3.21200000       9\\nO      6.15200000     4.88500000    -1.44700000      10\\nH      5.50700000     5.59100000    -1.36900000      11\\nH      5.89100000     4.42800000    -2.24900000      12\\nO      5.82300000     3.53700000    -3.94100000      13\\nH      6.31400000     2.71500000    -4.01100000      14\\nH      6.27500000     4.12300000    -4.55200000      15\\nO      8.86000000     5.34600000    -2.74900000      16\\nH      8.46500000     5.48700000    -3.61200000      17\\nH      8.10300000     5.25400000    -2.16700000      18\\n@MULTIPOLES\\nORDER 0\\n18\\n1      -0.67072060\\n2       0.33528566\\n3       0.33543494\\n4      -0.67055041\\n5       0.33526795\\n6       0.33528246\\n7      -0.67071744\\n8       0.33530196\\n9       0.33541547\\n10     -0.67067328\\n11      0.33530711\\n12      0.33536617\\n13     -0.67033801\\n14      0.33511794\\n15      0.33522007\\n16     -0.67061076\\n17      0.33527067\\n18      0.33534009\\nORDER 1\\n18\\n1       0.00080560    -0.20614866     0.20971724\\n2       0.12061233     0.19534842    -0.00437020\\n3      -0.12142692     0.00052018    -0.19496774\\n4       0.12128208    -0.02952128    -0.26636264\\n5       0.02123859    -0.14133483     0.17957852\\n6      -0.13641247     0.16937441     0.07336077\\n7       0.09874019     0.23525681     0.14626863\\n8       0.12395337    -0.17250412    -0.08710686\\n9      -0.21781638    -0.05098712    -0.05185290\\n10     -0.22463142     0.06171102    -0.17953533\\n11      0.15306787    -0.16979678    -0.02103891\\n12      0.06031266     0.11119740     0.19160536\\n13      0.23374507    -0.05843888    -0.16882627\\n14     -0.11567937     0.19764719     0.01487010\\n15     -0.10631260    -0.14219349     0.14548634\\n16     -0.28549425     0.01213605    -0.06959331\\n17      0.09187853    -0.03392224     0.20768211\\n18      0.17941729     0.02239518    -0.14158334\\nORDER 2\\n18\\n1      -3.82946639     0.38325366     0.37670941    -3.74967413     0.05052330    -3.75455511\\n2      -0.30950705     0.24207772    -0.02491219    -0.03438601    -0.02932713    -0.45703642\\n3      -0.30728365    -0.02024246     0.24345483    -0.45782760    -0.02099977    -0.03545985\\n4      -4.01843560    -0.42071343    -0.05131494    -3.51375485    -0.22366647    -3.80188906\\n5      -0.45058853    -0.01355906     0.03221603    -0.26454826    -0.26763342    -0.08600138\\n6      -0.24883865    -0.23689663    -0.12134441    -0.16776157     0.15242467    -0.38450342\\n7      -3.29787344    -0.20331433    -0.01414521    -3.86392383     0.23538575    -4.17197450\\n8      -0.32338453    -0.21790162    -0.11343377    -0.11675396     0.16862725    -0.36084577\\n9       0.03250781     0.14275944     0.13176702    -0.41814694     0.03110778    -0.41507187\\n10     -3.94275214    -0.29331092     0.07399890    -3.64071323     0.42271528    -3.75038841\\n11     -0.18497857    -0.27928790    -0.02349247    -0.15879913     0.01388996    -0.45721403\\n12     -0.40467967     0.08377330     0.14054460    -0.34162736     0.21069480    -0.05454010\\n13     -3.98817938    -0.10562980    -0.21981136    -3.34423732    -0.30490142    -4.00186338\\n14     -0.29272797    -0.25436911    -0.02365191    -0.05984813     0.05198185    -0.44875356\\n15     -0.31577317     0.16791558    -0.17641064    -0.26949266    -0.21058649    -0.21581787\\n16     -3.76934440     0.01997274    -0.13330175    -4.28009915    -0.16523022    -3.28439494\\n17     -0.34779244    -0.03700706     0.22659090    -0.43533186    -0.07018391    -0.01780670\\n18     -0.08414919     0.04219414    -0.26720999    -0.44242465    -0.02713904    -0.27419185\\n@POLARIZABILITIES\\nORDER 1 1\\n18\\n1       2.30791521     0.59643991     0.58658837     2.61100398    -0.10257978     2.60785108\\n2       1.30711897     0.90889808    -0.14203759     2.23138041     0.03064426     0.56363899\\n3       1.31473738    -0.12335800     0.91322978     0.56426081     0.06343491     2.22049295\\n4       2.07587445    -0.67088756    -0.21506644     2.80527167    -0.31578724     2.64886958\\n5       0.72279353     0.00751649     0.24603845     1.44539448    -1.05520399     1.93525324\\n6       1.51667446    -0.87135876    -0.35833460     1.82612167     0.59796749     0.76044651\\n7       3.17597352    -0.21979725     0.03827106     2.48641276     0.51074529     1.86433465\\n8       1.15711080    -0.91546924    -0.49534521     1.87545624     0.51559151     1.06960779\\n9       2.55558469     0.50338563     0.46922356     0.68876452    -0.02563589     0.85563095\\n10      2.34307583    -0.51517500     0.28388438     2.61854341     0.61181317     2.56593966\\n11      1.63329191    -1.01651663    -0.24252266     1.86507920     0.04484738     0.60387225\\n12      0.76841032     0.41201391     0.40278140     1.14519478     0.81893980     2.18754212\\n13      2.29130369    -0.22262257    -0.50824047     3.08292101    -0.43438854     2.16011928\\n14      1.20133160    -0.94148299     0.07075699     2.22478033     0.20477680     0.68049498\\n15      1.10983723     0.72098321    -0.53040485     1.39825592    -0.82860371     1.59671491\\n16      2.74369667     0.01646994    -0.12367515     1.60703710    -0.26089660     3.17839863\\n17      0.86347479    -0.13575237     0.83717068     0.86374986    -0.25419187     2.37549183\\n18      1.89994613     0.17584219    -1.10556667     0.83549707    -0.08475173     1.36594908\\nEXCLISTS\\n18 3\\n1     2   3\\n2     1   3\\n3     1   2\\n4     5   6\\n5     4   6\\n6     4   5\\n7     8   9\\n8     7   9\\n9     7   8\\n10   11  12\\n11   10  12\\n12   10  11\\n13   14  15\\n14   13  15\\n15   13  14\\n16   17  18\\n17   16  18\\n18   16  17\\n", "scf_type": "PK"}, "protocols": {}, "extras": {"wfn_qcvars_only": true}, "provenance": {"creator": "Psi4", "version": "1.4a2.dev1090", "routine": "psi4.driver.p4util.procutil"}}"""

    atres = psi4.schema_wrapper.run_qcschema(json.loads(jatin))
    assert psi4.compare_values(-0.03424830892844, atres.extras["qcvars"]["PE ENERGY"], 6, "PE Energy contribution")
    assert psi4.compare_values(-482.9411084900, atres.return_result, 6, "Total PE-SCF Energy")


#@pytest.mark.smoke
#@uusing("cct3")
#def test_cct3():
#    import cct3
#
#    psi4.geometry("""
#        units bohr
#        h -2.514213562373  -1.000000000000   0.000000000000
#        h -2.514213562373   1.000000000000   0.000000000000
#        h  2.514213562373  -1.000000000000   0.000000000000
#        h  2.514213562373   1.000000000000   0.000000000000
#        h -1.000000000000  -2.414213562373   0.000000000000
#        h -1.000000000000   2.414213562373   0.000000000000
#        h  1.000000000000  -2.414213562373   0.000000000000
#        h  1.000000000000   2.414213562373   0.000000000000
#        symmetry d2h
#    """)
#
#    def basisspec_psi4_yo__anonymous1234(mol, role):
#        bas = """
#            cartesian
#            ****
#            H   0
#            S   3  1.0000
#                  4.50038     0.0704800
#                  0.681277    0.407890
#                  0.151374    0.647670
#            ****
#        """
#        mol.set_basis_all_atoms("mbs_my", role=role)
#        return {"mbs_my": bas}
#
#    psi4.driver.qcdb.libmintsbasisset.basishorde["ANONYMOUS1234"] = basisspec_psi4_yo__anonymous1234
#
#    psi4.set_options({
#        "cct3__froz": 0,
#        "cct3__act_occ": 1,
#        "cct3__act_unocc": 1,
#        "cct3__etol": 16,
#        "cct3__calc_type": "cct3",
#        "basis": "anonymous1234",
#    })
#
##    atin = psi4.driver.p4util.state_to_atomicinput(driver="energy", method="ccsd", molecule=ethene_ethyne)
##    print(f'    jatin = """{atin.serialize("json")}"""')
##    assert 0
##
##    atres = psi4.schema_wrapper.run_qcschema(json.loads(jatin))
##    pprint.pprint(atres.dict())
#
#
#    ene = psi4.energy("cct3")
#    assert psi4.compare_values(-4.220587742726, ene, 10, "cc(t;3) energy")
#
