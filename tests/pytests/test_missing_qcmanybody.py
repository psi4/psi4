import sys

import pytest
from qcelemental.util import which_import

import psi4
from psi4.driver.task_planner import task_planner

pytestmark = [
    pytest.mark.quick,
    pytest.mark.psi,
    pytest.mark.api,
    pytest.mark.skipif(
        which_import("qcmanybody", return_bool=True),
        reason="qcmanybody installed; tests target missing-addon behavior",
    ),
]


_MSG = "Python module qcmanybody not found"


@pytest.fixture
def dimer():
    return psi4.geometry(
        """
        He 0 0 -2
        --
        He 0 0 2
        units au
        """
    )


@pytest.mark.parametrize(
    "driver,method,kwargs",
    [
        (psi4.energy, "hf", {"bsse_type": "cp"}),
        (psi4.gradient, "hf", {"bsse_type": "cp"}),
        (psi4.hessian, "hf", {"bsse_type": "cp", "dertype": "energy"}),
    ],
)
def test_missing_qcmanybody_via_driver(driver, method, kwargs, dimer):
    with pytest.raises(ModuleNotFoundError, match=_MSG):
        driver(method, molecule=dimer, **kwargs)


def test_missing_qcmanybody_via_return_plan(dimer):
    with pytest.raises(ModuleNotFoundError, match=_MSG):
        psi4.gradient("HF/cc-pV[D,T]Z", molecule=dimer, bsse_type="cp", return_plan=True)


@pytest.mark.parametrize("driver,method", [("energy", "MP2/cc-pVDZ"), ("gradient", "MP2/cc-pV[D,T]Z")])
def test_missing_qcmanybody_via_task_planner(driver, method, dimer):
    with pytest.raises(ModuleNotFoundError, match=_MSG):
        task_planner(driver, method, dimer, bsse_type="cp")


@pytest.mark.parametrize("schver", [1, 2])
def test_missing_qcmanybody_via_qcschema(schver, dimer):
    if schver == 1:
        qcs_input = {
            "schema_name": "qcschema_input",
            "schema_version": 1,
            "molecule": dimer.to_schema(dtype=2),
            "driver": "energy",
            "model": {"method": "hf", "basis": "sto-3g"},
            "keywords": {"function_kwargs": {"bsse_type": ["cp", "nocp"], "return_total_data": True}},
        }
    else:
        qcs_input = {
            "schema_name": "qcschema_atomic_input",
            "schema_version": 2,
            "molecule": dimer.to_schema(dtype=3),
            "specification": {
                "driver": "energy",
                "model": {"method": "hf", "basis": "sto-3g"},
                "keywords": {"function_kwargs": {"bsse_type": ["cp", "nocp"], "return_total_data": True}},
            },
        }

    return_version = 2 if schver == 1 and sys.version_info >= (3, 14) else -1
    ret = psi4.schema_wrapper.run_qcschema(qcs_input, return_version=return_version)

    assert not ret.success
    assert ret.error.error_type == "ModuleNotFoundError"
    assert _MSG in ret.error.error_message


def test_missing_qcmanybody_via_qcschema_levels(dimer):
    qcs_input = {
        "schema_name": "qcschema_atomic_input",
        "schema_version": 2,
        "molecule": dimer.to_schema(dtype=3),
        "specification": {
            "driver": "energy",
            "model": {"method": "", "basis": "(auto)"},
            "keywords": {
                "function_kwargs": {
                    "bsse_type": "nocp",
                    "levels": {1: "mp2/sto-3g", "supersystem": "scf/sto-3g"},
                }
            },
        },
    }

    ret = psi4.schema_wrapper.run_qcschema(qcs_input)

    assert not ret.success
    assert ret.error.error_type == "ModuleNotFoundError"
    assert _MSG in ret.error.error_message
