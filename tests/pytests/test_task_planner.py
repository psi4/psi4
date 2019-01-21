"""
Tests the Psi4 task driver
"""

import pytest
from utils import *
from addons import *
import time

import numpy as np

import psi4
from psi4.driver.task_planner import task_planner
from psi4.driver.task_base import SingleResult
from psi4.driver.driver_cbs import CBSComputer
from psi4.driver.driver_nbody import NBodyComputer
from psi4.driver.driver_findif import FinDifComputer

pytestmark = pytest.mark.quick

def test_single_result():
    mol = psi4.geometry("He")
    plan = task_planner("energy", "HF", mol, basis="sto-3G")

    assert isinstance(plan, SingleResult)
    assert plan.basis == "sto-3g"
    assert plan.method == "hf"
    assert plan.driver == "energy"
    assert plan.keywords == {}

    # Local options currently taking ~3 seconds per call!
    # psi4.set_options({"basis": "cc-pVQZ"})
    # plan = task_planner("energy", "HF", mol)
    # assert plan.basis == "cc-pvqz"


def test_single_result_cbs_unpack():
    mol = psi4.geometry("He")
    plan = task_planner("energy", "MP2/cc-pVDZ", mol)

    assert isinstance(plan, SingleResult)
    assert plan.basis == "cc-pvdz"
    assert plan.method == "mp2"
    assert plan.driver == "energy"


def test_cbs_extrapolation():
    mol = psi4.geometry("He")
    plan = task_planner("energy", "MP2/cc-pV[D,T]Z", mol)

    assert isinstance(plan, CBSComputer)
    assert len(plan.task_list) == 2

    assert isinstance(plan.task_list[0], SingleResult)
    assert plan.task_list[0].basis == "cc-pvdz"
    assert plan.task_list[0].driver == "energy"

    assert isinstance(plan.task_list[1], SingleResult)
    assert plan.task_list[1].basis == "cc-pvtz"


def test_cbs_extrapolation_grad():
    mol = psi4.geometry("He")
    plan = task_planner("gradient", "MP2/cc-pV[D,T]Z", mol)

    assert isinstance(plan, CBSComputer)
    assert len(plan.task_list) == 4

    assert isinstance(plan.task_list[0], SingleResult)
    assert plan.task_list[0].basis == "cc-pvtz"
    assert plan.task_list[0].method == "hf"
    assert plan.task_list[0].driver == "gradient"

    assert isinstance(plan.task_list[1], SingleResult)
    assert plan.task_list[1].basis == "cc-pvdz"
    assert plan.task_list[1].method == "mp2"
    assert plan.task_list[1].driver == "gradient"

    assert isinstance(plan.task_list[2], SingleResult)
    assert plan.task_list[2].basis == "cc-pvtz"
    assert plan.task_list[2].method == "mp2"
    assert plan.task_list[2].driver == "gradient"

    assert isinstance(plan.task_list[3], SingleResult)
    assert plan.task_list[3].basis == "cc-pvdz"
    assert plan.task_list[3].method == "hf"
    assert plan.task_list[3].driver == "gradient"


def test_nbody_dimer():
    mol = psi4.geometry("""
                        He 0 0 -5
                        --
                        He 0 0 5""")
    plan = task_planner("energy", "MP2/cc-pVDZ", mol, bsse_type="cp")

    assert isinstance(plan, NBodyComputer)
    # DGAS Note: This is wrong, for some reason cp is adding monomers in monomer basis.
    # See the CP builder in nbody:build_nbody_compute_list
    assert len(plan.task_list) == 5

    for key, result in plan.task_list.items():
        assert isinstance(result, SingleResult)


def test_nbody_dimer_cbs():
    mol = psi4.geometry("""
                        He 0 0 -5
                        --
                        He 0 0 5""")
    plan = task_planner("energy", "MP2/cc-pV[D,T]Z", mol, bsse_type="cp")

    assert isinstance(plan, NBodyComputer)
    # DGAS Note: This is wrong, for some reason cp is adding monomers in monomer basis.
    # See the CP builder in nbody:build_nbody_compute_list
    assert len(plan.task_list) == 5

    for key, result in plan.task_list.items():
        assert isinstance(result, CBSComputer)


def test_cbs_extrapolation_delta():
    mol = psi4.geometry("He")
    plan = task_planner("energy", "MP2/cc-pV[D,T]Z + D:ccsd(t)/cc-pv[tq]z", mol)

    assert isinstance(plan, CBSComputer)
    assert len(plan.task_list) == 3

    assert isinstance(plan.task_list[0], SingleResult)
    assert plan.task_list[0].basis == "cc-pvdz"
    assert plan.task_list[0].method == "mp2"

    assert isinstance(plan.task_list[1], SingleResult)
    assert plan.task_list[1].basis == "cc-pvtz"
    assert plan.task_list[1].method == "ccsd(t)"

    assert isinstance(plan.task_list[2], SingleResult)
    assert plan.task_list[2].basis == "cc-pvqz"
    assert plan.task_list[2].method == "ccsd(t)"


def test_findif_1_1():
    mol = psi4.geometry("H\nH 1 2.0\nunits au")
    plan = task_planner("gradient", "MP2/cc-pVDZ", mol)

    assert isinstance(plan, SingleResult)
    assert plan.basis == "cc-pvdz"
    assert plan.method == "mp2"
    assert plan.driver == "gradient"


def test_findif_1_0():
    mol = psi4.geometry("H\nH 1 2.0\nunits au")
    plan = task_planner("gradient", "MP2/cc-pVDZ", mol, dertype=0, findif_stencil_size=3, findif_step_size=0.005/math.sqrt(2/1.00782503223))

    displacements = {
        '0: -1': np.array([[ 0.    ,  0.    , -1.0025], [ 0.    ,  0.    ,  1.0025]]),
        '0: 1':  np.array([[ 0.    ,  0.    , -0.9975], [ 0.    ,  0.    ,  0.9975]]),
        'reference':  np.array([[ 0.    ,  0.    , -1.0], [ 0.    ,  0.    ,  1.0]]),
    }

    assert isinstance(plan, FinDifComputer)
    assert len(plan.task_list) == 3

    key = 'reference'
    assert isinstance(plan.task_list[key], SingleResult)
    assert plan.task_list[key].basis == "cc-pvdz"
    assert plan.task_list[key].method == "mp2"
    assert plan.task_list[key].driver == "energy"
    assert np.allclose(plan.task_list[key].molecule.geometry().np, displacements[key])

    key = '0: -1'
    assert isinstance(plan.task_list[key], SingleResult)
    assert plan.task_list[key].basis == "cc-pvdz"
    assert plan.task_list[key].method == "mp2"
    assert plan.task_list[key].driver == "energy"
    assert np.allclose(plan.task_list[key].molecule.geometry().np, displacements[key])

    key = '0: 1'
    assert isinstance(plan.task_list[key], SingleResult)
    assert plan.task_list[key].basis == "cc-pvdz"
    assert plan.task_list[key].method == "mp2"
    assert plan.task_list[key].driver == "energy"
    assert np.allclose(plan.task_list[key].molecule.geometry().np, displacements[key])


