import pytest
import psi4

from .utils import compare_values

methods = pytest.mark.parametrize("inp", [
    pytest.param({"method": 'PNO', 'options': {"local_cutoff": 0.0}}, id='PNO 0'),
    pytest.param({"method": 'PNO', 'options': {"local_cutoff": 1e-5}}, id='PNO 1e-5'),
    pytest.param({"method": 'PNO++', 'options': {"local_cutoff": 0.0, "local_pert": "DIPOLE"}}, id='PNO++ 0'),
    pytest.param({"method": 'PNO++', 'options': {"local_cutoff": 1e-5, "local_pert": "DIPOLE"}}, id='PNO++ 1e-5'),
    pytest.param({"method": 'CPNO++', 'options': {"local_cutoff": 0.0, "local_pert": "DIPOLE", "unpert_cutoff": 0.0}}, id='CPNO++ 0'),
    pytest.param({"method": 'CPNO++', 'options': {"local_cutoff": 1e-5, "local_pert": "DIPOLE", "unpert_cutoff": 1e-6}}, id='CPNO++ 1e-5'),
    ]
)

@methods
def test_en(inp):
    ccsd_e_dict = {"PNO": -0.3832736390, "PNO++": -0.2726823209, "CPNO++": -0.3948947549}
    h2o2 = psi4.geometry("""
         O     -0.028962160801    -0.694396279686    -0.049338350190
         O      0.028962160801     0.694396279686    -0.049338350190
         H      0.350498145881    -0.910645626300     0.783035421467
         H     -0.350498145881     0.910645626300     0.783035421467
         symmetry c1
    """)
    psi4.set_module_options("SCF", {"E_CONVERGENCE": 1e-12})
    psi4.set_module_options("SCF", {"D_CONVERGENCE": 1e-12})
    psi4.set_options({"basis": "cc-pVDZ", "scf_type": "pk",
              "freeze_core": "false", "e_convergence": 1e-10,
              "d_convergence": 1e-10, "r_convergence": 1e-10, 
              "local": "true", "local_filter_singles": "true",
              "restart": 0, "maxiter": 100,
              "local_method": inp["method"], "omega": ["0.077357"]})
    psi4.set_options(inp['options'])

    # Do the local simulation energy calculation
    psi4.energy("ccsd")

    ccsd_e = psi4.core.variable("CCSD correlation energy")

    if inp["options"]["local_cutoff"] == 0.0:
        # Do the regular energy calculation
        psi4.set_options({"local": "false", 
            "local_filter_singles": "false"})
        psi4.energy("ccsd")

        conv_ccsd_e = psi4.core.variable("CCSD correlation energy")

        assert compare_values(ccsd_e, conv_ccsd_e, 6, "CCSD correlation energy at Tcut 0")
    else:
        assert compare_values(ccsd_e, ccsd_e_dict[inp["method"]], 6, "CCSD correlation energy at Tcut 1e-5")

@methods
def test_properties(inp):
    polar_dict = {"PNO": 6.7388409598, "PNO++": 7.6884792679, "CPNO++": 8.3693956048}
    lg_dict = {"PNO": -6.12398, "PNO++": -50.00546, "CPNO++": -68.78573}
    mvg_dict = {"PNO": -17.50370, "PNO++": -140.39388, "CPNO++": -291.91706}

    h2o2 = psi4.geometry("""
         O     -0.028962160801    -0.694396279686    -0.049338350190
         O      0.028962160801     0.694396279686    -0.049338350190
         H      0.350498145881    -0.910645626300     0.783035421467
         H     -0.350498145881     0.910645626300     0.783035421467
         symmetry c1
    """)
    psi4.set_module_options("SCF", {"E_CONVERGENCE": 1e-12})
    psi4.set_module_options("SCF", {"D_CONVERGENCE": 1e-12})
    psi4.set_options({"basis": "cc-pVDZ", "scf_type": "pk",
              "freeze_core": "false", "e_convergence": 1e-10,
              "d_convergence": 1e-10, "r_convergence": 1e-10, 
              "local": "true", "local_filter_singles": "true",
              "local_method": inp["method"], "restart": 0,
              "maxiter": 100,
              "gauge": "both", "omega": ["0.077357"]})
    psi4.set_options(inp['options'])
    # Do the local simulation linear response calculation
    psi4.properties("ccsd", properties=["polarizability", "rotation"])

    polar = psi4.core.variable("CCSD DIPOLE POLARIZABILITY @ 589NM")
    optrot_lg = psi4.core.variable("CCSD SPECIFIC ROTATION (LEN) @ 589NM")
    optrot_mvg = psi4.core.variable("CCSD SPECIFIC ROTATION (MVG) @ 589NM")

    if inp["options"]["local_cutoff"] == 0.0:
        # Do the regular linear response calculation
        psi4.set_options({"local": "false", 
            "local_filter_singles": "false"})
        psi4.properties("ccsd", properties=["polarizability", "rotation"])

        conv_polar = psi4.core.variable("CCSD DIPOLE POLARIZABILITY @ 589NM")
        conv_optrot_lg = psi4.core.variable("CCSD SPECIFIC ROTATION (LEN) @ 589NM")
        conv_optrot_mvg = psi4.core.variable("CCSD SPECIFIC ROTATION (MVG) @ 589NM")

        assert compare_values(polar, conv_polar, 4, "CCSD dipole polarizability at Tcut 0")
        assert compare_values(optrot_lg, conv_optrot_lg, 4, "CCSD optical rotation (LG) at Tcut 0")
        assert compare_values(optrot_mvg, conv_optrot_mvg, 4, "CCSD optical rotation (MVG) at Tcut 0")
    else:
        assert compare_values(polar, polar_dict[inp["method"]], 4, "CCSD dipole polarizability at Tcut 1e-5")
        assert compare_values(optrot_lg, lg_dict[inp["method"]], 4, "CCSD optical rotation (LG) at Tcut 1e-5")
        assert compare_values(optrot_mvg, mvg_dict[inp["method"]], 4, "CCSD optical rotation (MVG) at Tcut 1e-5")
