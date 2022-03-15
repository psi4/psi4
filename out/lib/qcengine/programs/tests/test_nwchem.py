"""Tests for NWChem functionality"""
import numpy as np
import pytest
import qcelemental as qcel
from qcelemental.testing import compare_values

import qcengine as qcng
from qcengine.testing import using

# Molecule where autoz fails
_auto_z_problem = xyz = """C                    15.204188380000    -3.519180270000   -10.798726560000
C                    15.097645630000    -2.650246400000    -8.505033680000
C                    14.976892130000    -1.867030510000    -6.378827230000
C                    14.827093540000    -1.105281910000    -4.239222620000
C                    14.635399720000    -0.362921900000    -2.101035300000
N                    14.357893440000     0.315376390000    -0.000359050000
C                    16.031227020000     2.397930170000     0.801905280000
C                    15.337300690000     2.925144860000     3.538890120000
C                    14.742074760000     3.377413020000     5.886685850000
C                    14.048167320000     3.904627710000     8.623670690000
N                    12.688887320000     2.073331920000     9.820037400000
C                    14.474640720000     6.115210430000     9.591909660000
C                    14.946826580000     8.308766730000    10.603064320000
C                    17.506555110000     9.181593430000    10.935939570000
C                    18.224273090000    10.758664370000    12.778309160000
O                    20.684526430000    11.555996510000    13.018002020000
C                    20.952489590000    13.130270650000    14.964495520000
C                    21.186059740000    14.590046290000    16.687812370000
C                    21.422747940000    16.379824800000    18.663483240000
C                    22.929993500000    15.902177630000    20.539319870000
C                    24.426014080000    15.474645990000    22.431805000000
N                    24.775991360000    13.082933010000    23.323528960000
O                    18.709082330000     1.890803270000     0.670399240000
O                    20.008439110000     4.184817400000     1.512272230000
H                    16.847872170000    -3.183338140000   -11.974003930000
H                    13.649907540000    -4.604695650000   -11.575422900000
H                    15.665451630000     4.107792170000    -0.318135390000
H                    13.044004660000     0.318248780000     9.124920540000
H                    12.835662350000     2.164511200000    11.730569410000
H                    13.401843200000     9.476598570000    11.281381510000
H                    18.900965120000     8.515559460000     9.580741380000
H                    16.945759980000    11.469295880000    14.211420760000
H                    20.389540180000    18.145263640000    18.570527610000
H                    25.288598470000    16.991453560000    23.511896860000
H                    24.436823320000    11.711936710000    22.016121940000
H                    26.442597520000    12.834282850000    24.246149950000
H                    20.850425490000     3.414376060000     2.960577230000"""


@pytest.fixture
def nh2():
    smol = """
 # R=1.008 #A=105.0
 0 2
 N   0.000000000000000   0.000000000000000  -0.145912918634892
 H   0.000000000000000  -1.511214298139000   1.013682596946108
 H   0.000000000000000   1.511214298139000   1.013682596946108
 units au
 symmetry c1
"""
    return qcel.models.Molecule.from_data(smol)


@using("nwchem")
def test_b3lyp(nh2):
    # Run NH2
    resi = {"molecule": nh2, "driver": "energy", "model": {"method": "b3lyp", "basis": "3-21g"}}
    res = qcng.compute(resi, "nwchem", raise_error=True, return_dict=True)

    # Make sure the calculation completed successfully
    assert compare_values(-55.554037, res["return_result"], atol=1e-3)
    assert res["driver"] == "energy"
    assert "provenance" in res
    assert res["success"] is True

    # Check the other status information
    assert res["extras"]["qcvars"]["N ALPHA ELECTRONS"] == "5"
    assert res["extras"]["qcvars"]["N ATOMS"] == "3"
    assert res["extras"]["qcvars"]["N BASIS FUNCTIONS"] == "13"

    # Make sure the properties parsed correctly
    assert compare_values(-55.554037, res["properties"]["return_energy"], atol=1e-3)
    assert res["properties"]["calcinfo_natom"] == 3
    assert res["properties"]["calcinfo_nalpha"] == 5
    assert res["properties"]["calcinfo_nbasis"] == 13


@using("nwchem")
def test_hess(nh2):
    resi = {"molecule": nh2, "driver": "hessian", "model": {"method": "b3lyp", "basis": "3-21g"}}
    res = qcng.compute(resi, "nwchem", raise_error=True, return_dict=False)
    assert compare_values(-3.5980754370e-02, res.return_result[0, 0], atol=1e-3)
    assert compare_values(0, res.return_result[1, 0], atol=1e-3)
    assert compare_values(0.018208307756, res.return_result[3, 0], atol=1e-3)
    assert np.allclose(res.return_result, res.return_result.T, atol=1e-8)  # Should be symmetric about diagonal

    # Test that the Hessian changes with rotation, but that its determinants remain the same
    shifted_nh2, _ = nh2.scramble(do_shift=False, do_mirror=False, do_rotate=True, do_resort=False)

    resi["molecule"] = shifted_nh2
    res_shifted = qcng.compute(resi, "nwchem", raise_error=True, return_dict=False)
    assert not np.allclose(res.return_result, res_shifted.return_result, atol=1e-8)
    assert np.isclose(np.linalg.det(res.return_result), np.linalg.det(res_shifted.return_result))


@using("nwchem")
def test_gradient(nh2):
    resi = {
        "molecule": nh2,
        "driver": "gradient",
        "model": {"method": "b3lyp", "basis": "3-21g"},
        "keywords": {"dft__convergence__gradient": "1e-6"},
    }
    res = qcng.compute(resi, "nwchem", raise_error=True, return_dict=True)
    assert compare_values(4.22418267e-2, res["return_result"][2], atol=1e-7)  # Beyond accuracy of NWChem stdout

    # Rotate the molecule and verify that the gradient changes
    shifted_nh2, _ = nh2.scramble(do_shift=False, do_mirror=False, do_rotate=True, do_resort=False)

    resi["molecule"] = shifted_nh2
    res_shifted = qcng.compute(resi, "nwchem", raise_error=True, return_dict=True)

    assert not compare_values(4.22418267e-2, res_shifted["return_result"][2], atol=1e-7)

    # Make sure the two matrices still have the same determinant and norms, as they are just rotations of each other
    #  I am leveraging the fact that the gradients are square, 3x3 matrices just by happenstance of the
    #  test molecule having 3 atoms
    orig_grads = np.reshape(res["return_result"], (-1, 3))
    shif_grads = np.reshape(res_shifted["return_result"], (-1, 3))

    # Test that the magnitude of forces are the same
    assert np.allclose(np.linalg.norm(orig_grads, ord=2, axis=1), np.linalg.norm(shif_grads, ord=2, axis=1))

    # Test that the determinants are the same
    orig_det = np.linalg.det(orig_grads)
    shif_det = np.linalg.det(shif_grads)

    assert np.allclose(orig_det, shif_det)


@pytest.fixture
def h20():
    water = """
-1 2
O 0 0 0
H 0 0 1
H 0 1 0
    """
    return qcel.models.Molecule.from_data(water)


@using("nwchem")
def test_dipole(h20):
    # Run NH2
    resi = {
        "molecule": h20,
        "driver": "properties",
        "model": {"method": "dft", "basis": "3-21g"},
        "keywords": {"dft__xc": "b3lyp", "property__dipole": True},
    }
    res = qcng.compute(resi, "nwchem", raise_error=True, return_dict=True)

    # Make sure the calculation completed successfully
    assert compare_values(-75.764944, res["return_result"], atol=1e-3)
    assert res["driver"] == "properties"
    assert "provenance" in res
    assert res["success"] is True

    # Check the other status information
    assert res["extras"]["qcvars"]["N ALPHA ELECTRONS"] == "6"
    assert res["extras"]["qcvars"]["N ATOMS"] == "3"
    assert res["extras"]["qcvars"]["N BASIS FUNCTIONS"] == "13"

    # Make sure the properties parsed correctly
    assert compare_values(-75.764944, res["properties"]["return_energy"], atol=1e-3)
    assert res["properties"]["calcinfo_natom"] == 3
    assert res["properties"]["calcinfo_nalpha"] == 6
    assert res["properties"]["calcinfo_nbasis"] == 13
    # Make sure Dipole Moment and center of charge parsed correctly
    assert compare_values(0.272949872, float(res["extras"]["qcvars"]["TOTAL DIPOLE MOMENT"]), atol=1e-5)
    assert compare_values(-0.00, float(res["extras"]["qcvars"]["DIPOLE MOMENT"][0]), atol=1e-3)
    assert compare_values(-0.00, float(res["extras"]["qcvars"]["DIPOLE MOMENT"][1]), atol=1e-3)
    assert compare_values(-0.272949872, float(res["extras"]["qcvars"]["DIPOLE MOMENT"][2]), atol=1e-5)


@pytest.fixture
def h20v2():
    water = """
O 0 0 0
H 0 0 1
H 0 1 0
    """
    return qcel.models.Molecule.from_data(water)


@using("nwchem")
def test_homo_lumo(h20v2):
    # Run NH2
    resi = {
        "molecule": h20v2,
        "driver": "energy",
        "model": {"method": "dft", "basis": "3-21g"},
        "keywords": {"dft__xc": "b3lyp"},
    }
    res = qcng.compute(resi, "nwchem", raise_error=True, return_dict=True)

    # Make sure the calculation completed successfully
    assert compare_values(-75.968095, res["return_result"], atol=1e-3)
    assert res["driver"] == "energy"
    assert "provenance" in res
    assert res["success"] is True

    # Check the other status information
    assert res["extras"]["qcvars"]["N ALPHA ELECTRONS"] == "5"
    assert res["extras"]["qcvars"]["N ATOMS"] == "3"
    assert res["extras"]["qcvars"]["N BASIS FUNCTIONS"] == "13"

    # Make sure the properties parsed correctly
    assert compare_values(-75.968095, res["properties"]["return_energy"], atol=1e-3)
    assert res["properties"]["calcinfo_natom"] == 3
    assert res["properties"]["calcinfo_nalpha"] == 5
    assert res["properties"]["calcinfo_nbasis"] == 13
    # Make sure Dipole Moment and center of charge parsed correctly
    assert compare_values(-0.2636515, float(res["extras"]["qcvars"]["HOMO"][0]), atol=1e-5)
    assert compare_values(0.08207131, float(res["extras"]["qcvars"]["LUMO"][0]), atol=1e-5)


@using("nwchem")
def test_geometry_bug():
    """Make sure that the harvester does not crash if NWChem's autosym moves atoms too far"""

    # Example molecule that has an RMSD of 2e-4 after NWChem symmetrizes the coordinates
    xyz = """6
Properties=species:S:1:pos:R:3 unique_id=2da214efcf4e4c277fabc5b2b6ca6f32 pbc="F F F"
C      -0.00828817       1.39046978      -0.00560069
O      -0.00797038      -0.02504537       0.02030606
H       1.00658338       1.81556366       0.00348335
H      -0.54657475       1.79916975      -0.87390126
H      -0.52288871       1.72555240       0.89907326
H       0.44142019      -0.33354425      -0.77152059"""
    mol = qcel.models.Molecule.from_data(xyz)
    qcng.compute(
        {"molecule": mol, "model": {"method": "b3lyp", "basis": "6-31g"}, "driver": "gradient"},
        "nwchem",
        raise_error=True,
    )


@using("nwchem")
def test_autoz_error():
    """Test ability to turn off autoz"""
    # Large molecule that leads to an AutoZ error
    mol = qcel.models.Molecule.from_data(_auto_z_problem)
    result = qcng.compute(
        {
            "molecule": mol,
            "model": {"method": "hf", "basis": "sto-3g"},
            "driver": "energy",
            "protocols": {"error_correction": {"default_policy": False}},
        },  # Turn off error correction
        "nwchem",
        raise_error=False,
    )

    assert not result.success
    assert "Error when generating redundant atomic coordinates" in result.error.error_message

    # Turn off autoz
    result = qcng.compute(
        {
            "molecule": mol,
            "model": {"method": "hf", "basis": "sto-3g"},
            "driver": "energy",
            "keywords": {"geometry__noautoz": True},
        },
        "nwchem",
        raise_error=False,
    )

    # Ok if it crashes for other reasons
    assert "Error when generating redundant atomic coordinates" not in result.error.error_message


@using("nwchem")
def test_autoz_error_correction():
    """See if error correction for autoz works"""

    # Large molecule that leads to an AutoZ error
    mol = qcel.models.Molecule.from_data(_auto_z_problem)
    result = qcng.compute(
        {
            "molecule": mol,
            "model": {"method": "hf", "basis": "sto-3g"},
            "driver": "energy",
            "keywords": {"scf__maxiter": 250, "scf__thresh": 1e-1},
        },
        "nwchem",
        raise_error=True,
    )

    assert result.success
    assert "geom_binvr" in result.extras["observed_errors"]
    assert result.extras["observed_errors"]["geom_binvr"]["keyword_updates"] == {"geometry__noautoz": True}


@pytest.mark.parametrize(
    "method, keyword, init_iters, use_tce",
    [
        ["b3lyp", "dft__maxiter", 4, False],
        ["b3lyp", "dft__iterations", 4, False],
        ["hf", "scf__maxiter", 2, False],
        ["mp2", "scf__maxiter", 2, False],
        ["mp2", "scf__maxiter", 2, True],
        ["ccsd", "tce__maxiter", 8, True],
        ["ccsd", "ccsd__maxiter", 4, False],
    ],
)
@using("nwchem")
def test_conv_threshold(h20v2, method, keyword, init_iters, use_tce):
    result = qcng.compute(
        {
            "molecule": h20v2,
            "model": {"method": method, "basis": "sto-3g"},
            "driver": "energy",
            "keywords": {
                keyword: init_iters,
                "qc_module": use_tce,
                "scf__uhf": True,  # UHF needed for SCF test
            },
        },
        "nwchem",
        raise_error=True,
    )

    assert result.success
    assert "convergence_failed" in result.extras["observed_errors"]
    assert result.extras["observed_errors"]["convergence_failed"]["keyword_updates"] == {keyword: init_iters * 4}
