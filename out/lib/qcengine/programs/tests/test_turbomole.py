import numpy as np
import pytest
import qcelemental
from qcelemental.testing import compare_values

import qcengine as qcng
from qcengine.testing import using


@pytest.fixture
def h2o():
    mol = qcelemental.models.Molecule.from_data(
        """
            O 0.000000000000     0.000000000000    -0.068516245955
            H 0.000000000000    -0.790689888800     0.543701278274
            H 0.000000000000     0.790689888800     0.543701278274
    """
    )
    return mol


@pytest.fixture
def h2o_ricc2_def2svp():
    """NumForce calls only make sense for stationary points. So this
    geometry was optimized at the ricc2/def2-svp level of theory and
    can be used to run NumForce with ricc2."""

    mol = qcelemental.models.Molecule.from_data(
        """
        O     0.0000000    0.0000000   -0.0835835
        H     0.7501772    0.0000000    0.5210589
        H    -0.7501772    0.0000000    0.5210589
    """
    )
    return mol


@pytest.mark.parametrize(
    "method, keywords, ref_energy",
    [
        pytest.param("hf", {}, -75.95536954370, marks=using("turbomole")),
        pytest.param("pbe0", {"grid": "m5"}, -76.27371135900, marks=using("turbomole")),
        pytest.param("ricc2", {}, -76.1603807755, marks=using("turbomole")),
        pytest.param("rimp2", {}, -76.1593614075, marks=using("turbomole")),
    ],
)
def test_turbomole_energy(method, keywords, ref_energy, h2o):
    resi = {"molecule": h2o, "driver": "energy", "model": {"method": method, "basis": "def2-SVP"}, "keywords": keywords}

    res = qcng.compute(resi, "turbomole", raise_error=True, return_dict=True)

    assert res["driver"] == "energy"
    assert res["success"] is True

    assert compare_values(ref_energy, res["return_result"])


@pytest.mark.parametrize(
    "method, keywords, ref_norm",
    [
        pytest.param("hf", {}, 0.099340, marks=using("turbomole")),
        pytest.param("pbe0", {"grid": "m5"}, 0.0606266, marks=using("turbomole")),
        pytest.param("ricc2", {}, 0.059378, marks=using("turbomole")),
        pytest.param("rimp2", {}, 0.061576, marks=using("turbomole")),
    ],
)
def test_turbomole_gradient(method, keywords, ref_norm, h2o):
    resi = {
        "molecule": h2o,
        "driver": "gradient",
        "model": {"method": method, "basis": "def2-SVP"},
        "keywords": keywords,
    }

    res = qcng.compute(resi, "turbomole", raise_error=True)

    assert res.driver == "gradient"
    assert res.success is True
    assert res.properties.return_energy

    grad = res.return_result
    grad_norm = np.linalg.norm(grad)
    assert compare_values(ref_norm, grad_norm)


@using("turbomole")
def test_turbomole_ri_dsp(h2o):
    resi = {
        "molecule": h2o,
        "driver": "energy",
        "model": {"method": "b-p", "basis": "def2-SVP"},
        "keywords": {"ri": True, "d3bj": True},
    }

    res = qcng.compute(resi, "turbomole", raise_error=True)

    assert res.driver == "energy"
    assert res.success is True

    energy = res.return_result
    ref_energy = -76.36275642866
    assert compare_values(ref_energy, energy)


def assert_hessian(H, ref_eigvals, ref_size):
    w, v = np.linalg.eigh(H)
    last_eigvals = w[-3:]
    # Hessian must be symmetric
    np.testing.assert_allclose(H, H.T)
    # Check eigenvalues
    np.testing.assert_allclose(last_eigvals, ref_eigvals)
    # Hessian must be of shape (3N x 3N)
    assert H.shape == (ref_size, ref_size)


@using("turbomole")
@pytest.mark.parametrize(
    "method, keywords, ref_eigvals",
    [
        ("hf", {}, (2.00771683e-01, 7.77977644e-01, 9.91091318e-01)),
        ("pbe0", {"grid": "m5"}, (1.72092719e-01, 7.38603449e-01, 9.73783598e-01)),
        ("b-p", {"grid": "m5", "ri": True}, (1.59729409e-01, 7.21364827e-01, 9.63399519e-01)),
    ],
)
def test_turbomole_hessian(method, keywords, ref_eigvals, h2o):
    resi = {
        "molecule": h2o,
        "driver": "hessian",
        "model": {
            "method": method,
            "basis": "def2-SVP",
        },
        "keywords": keywords,
    }

    res = qcng.compute(resi, "turbomole", raise_error=True)
    H = res.return_result
    size = h2o.geometry.size

    assert res.driver == "hessian"
    assert res.success is True
    assert res.properties.return_energy
    assert_hessian(H, ref_eigvals, size)


@using("turbomole")
@pytest.mark.parametrize(
    "method, keywords, ref_eigvals",
    [
        ("ricc2", {}, (1.65405531e-01, 9.63690706e-01, 1.24676634e00)),
    ],
)
def test_turbomole_num_hessian(method, keywords, ref_eigvals, h2o_ricc2_def2svp):
    resi = {
        "molecule": h2o_ricc2_def2svp,
        "driver": "hessian",
        "model": {
            "method": method,
            "basis": "def2-SVP",
        },
        "keywords": keywords,
    }

    res = qcng.compute(resi, "turbomole", raise_error=True)
    H = res.return_result

    size = h2o_ricc2_def2svp.geometry.size

    assert res.driver == "hessian"
    assert res.success is True
    assert res.properties.return_energy
    assert_hessian(H, ref_eigvals, size)
