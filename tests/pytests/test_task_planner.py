"""
Tests the Psi4 task driver
"""

import pytest
from .utils import *
from .addons import *

import math

import numpy as np

import psi4
from psi4.driver.task_planner import task_planner
from psi4.driver.task_base import AtomicComputer
from psi4.driver.driver_cbs import CompositeComputer
from psi4.driver.driver_nbody import ManyBodyComputer
from psi4.driver.driver_findif import FiniteDifferenceComputer

pytestmark = pytest.mark.quick


def test_single_result():
    mol = psi4.geometry("He")
    plan = task_planner("energy", "HF", mol, basis="sto-3G")

    assert isinstance(plan, AtomicComputer)
    assert plan.basis == "sto-3g"
    assert plan.method == "hf"
    assert plan.driver == "energy"
    assert plan.keywords == {'function_kwargs': {}}

    # Local options currently taking ~3 seconds per call!
    # psi4.set_options({"basis": "cc-pVQZ"})
    # plan = task_planner("energy", "HF", mol)
    # assert plan.basis == "cc-pvqz"


def test_single_result_cbs_unpack():
    mol = psi4.geometry("He")
    plan = task_planner("energy", "MP2/cc-pVDZ", mol)

    assert isinstance(plan, AtomicComputer)
    assert plan.basis == "cc-pvdz"
    assert plan.method == "mp2"
    assert plan.driver == "energy"


def test_cbs_extrapolation():
    mol = psi4.geometry("He")
    plan = task_planner("energy", "MP2/cc-pV[D,T]Z", mol)

    assert isinstance(plan, CompositeComputer)
    assert len(plan.task_list) == 2

    for ires, plan2 in enumerate(plan.task_list):
        assert isinstance(plan2, AtomicComputer)
        assert plan2.basis == ["cc-pvdz", "cc-pvtz"][ires]
        assert plan2.method == "mp2"
        assert plan2.driver == "energy"


def test_cbs_extrapolation_gradient():
    mol = psi4.geometry("He")
    plan = task_planner("gradient", "MP2/cc-pV[D,T]Z", mol)

    assert isinstance(plan, CompositeComputer)
    assert len(plan.task_list) == 4

    for ires, plan2 in enumerate(plan.task_list):
        assert isinstance(plan2, AtomicComputer)
        assert plan2.basis == ["cc-pvtz", "cc-pvdz", "cc-pvtz", "cc-pvdz"][ires]
        assert plan2.method == ["hf", "mp2", "mp2", "hf"][ires]
        assert plan2.driver == "gradient"


@pytest.mark.parametrize("mtd, kw", [
    ('mp2', {'dertype': 0}),
    ('mp5', {}),
])
def test_cbs_extrapolation_gradient_1_0(mtd, kw):
    mol = psi4.geometry("H\nH 1 2.0\nunits au")
    plan = task_planner("gradient", f"{mtd}/cc-pV[D,T]Z", mol, **kw, findif_stencil_size=3, findif_step_size=0.005/math.sqrt(2/1.00782503223), findif_verbose=4)

    assert isinstance(plan, FiniteDifferenceComputer)
    assert len(plan.task_list) == 3

    for key, plan2 in plan.task_list.items():
        assert isinstance(plan2, CompositeComputer)
        assert len(plan2.task_list) == 2

        for ires, plan3 in enumerate(plan2.task_list):
            assert isinstance(plan3, AtomicComputer)
            assert plan3.basis == ["cc-pvdz", "cc-pvtz"][ires]
            assert plan3.method == mtd
            assert plan3.driver == "energy"


def test_nbody_dimer():
    mol = psi4.geometry("""
                        He 0 0 -5
                        --
                        He 0 0 5
                        units au""")
    plan = task_planner("energy", "MP2/cc-pVDZ", mol, bsse_type="cp")

    ghostiness = {
        "1_((2,), (2,))": (['He'], [True]),
        "1_((1,), (1,))": (['He'], [True]),
        "1_((1, 2), (1, 2))": (['He', 'He'], [True, True]),
        "1_((1,), (1, 2))": (['He', 'He'], [True, False]),
        "1_((2,), (1, 2))": (['He', 'He'], [False, True]),
    }

    assert isinstance(plan, ManyBodyComputer)
    # DGAS Note: This is wrong, for some reason cp is adding monomers in monomer basis.
    # See the CP builder in nbody:build_nbody_compute_list
    assert len(plan.task_list) == 5

    for k2, plan2 in plan.task_list.items():
        assert isinstance(plan2, AtomicComputer)
        assert plan2.basis == "cc-pvdz"
        assert plan2.method == "mp2"
        assert plan2.driver == "energy"

        kmol = plan2.molecule.to_schema(dtype=2)
        assert kmol['symbols'] == ghostiness[k2][0]
        assert kmol['real'] == ghostiness[k2][1]


def test_nbody_dimer_gradient():
    mol = psi4.geometry("""
                        He 0 0 -5
                        --
                        He 0 0 5
                        units au""")
    plan = task_planner("gradient", "MP2/cc-pVDZ", mol, bsse_type="cp")

    ghostiness = {
        "1_((2,), (2,))": (['He'], [True]),
        "1_((1,), (1,))": (['He'], [True]),
        "1_((1, 2), (1, 2))": (['He', 'He'], [True, True]),
        "1_((1,), (1, 2))": (['He', 'He'], [True, False]),
        "1_((2,), (1, 2))": (['He', 'He'], [False, True]),
    }

    assert isinstance(plan, ManyBodyComputer)
    # DGAS Note: This is wrong, for some reason cp is adding monomers in monomer basis.
    # See the CP builder in nbody:build_nbody_compute_list
    assert len(plan.task_list) == 5

    for k2, plan2 in plan.task_list.items():
        assert isinstance(plan2, AtomicComputer)
        assert plan2.basis == "cc-pvdz"
        assert plan2.method == "mp2"
        assert plan2.driver == "gradient"

        kmol = plan2.molecule.to_schema(dtype=2)
        assert kmol['symbols'] == ghostiness[k2][0]
        assert kmol['real'] == ghostiness[k2][1]


@pytest.mark.parametrize("mtd, kw", [
    ('mp2', {'dertype': 0}),
    ('mp5', {}),
])
def test_nbody_dimer_gradient_1_0(mtd, kw):
    mol = psi4.geometry("""
                        He 0 0 -5
                        --
                        He 0 0 5
                        units au
                        """)
    plan = task_planner("gradient", f"{mtd}/cc-pVDZ", mol, **kw, bsse_type="cp", findif_stencil_size=3, findif_step_size=0.005/math.sqrt(2/4.00260325413))

    displacements = {
        '0: -1': np.array([[ 0.    ,  0.    , -5.0025], [ 0.    ,  0.    ,  5.0025]]),
        '0: 1':  np.array([[ 0.    ,  0.    , -4.9975], [ 0.    ,  0.    ,  4.9975]]),
        'reference':  np.array([[ 0.    ,  0.    , -5.0], [ 0.    ,  0.    ,  5.0]]),
    }

    nbody_displacements = {
        "1_((2,), (2,))": {k: v[1] for k, v in displacements.items()},
        "1_((1,), (1,))": {k: v[0] for k, v in displacements.items()},
        "1_((1, 2), (1, 2))": displacements,
        "1_((1,), (1, 2))": displacements,
        "1_((2,), (1, 2))": displacements,
    }

    assert isinstance(plan, ManyBodyComputer)
    # DGAS Note: This is wrong, for some reason cp is adding monomers in monomer basis.
    # See the CP builder in nbody:build_nbody_compute_list
    assert len(plan.task_list) == 5

    ghostiness = {
        "1_((2,), (2,))": (['He'], [True]),
        "1_((1,), (1,))": (['He'], [True]),
        "1_((1, 2), (1, 2))": (['He', 'He'], [True, True]),
        "1_((1,), (1, 2))": (['He', 'He'], [True, False]),
        "1_((2,), (1, 2))": (['He', 'He'], [False, True]),
    }

    for k2, plan2 in plan.task_list.items():
        assert isinstance(plan2, FiniteDifferenceComputer)
        assert plan2.basis == "cc-pvdz"
        assert plan2.method == f"{mtd}"
        assert plan2.driver == "gradient"

        kmol = plan2.molecule.to_schema(dtype=2)
        assert kmol['symbols'] == ghostiness[k2][0]
        assert kmol['real'] == ghostiness[k2][1]

        for k3, plan3 in plan2.task_list.items():
            assert isinstance(plan3, AtomicComputer)
            assert plan3.basis == "cc-pvdz"
            assert plan3.method == f"{mtd}"
            assert plan3.driver == "energy"
            assert np.allclose(plan3.molecule.geometry().np, nbody_displacements[k2][k3])
            #assert plan3.keywords['SCF__E_CONVERGENCE'] == 1.e-6
            #assert plan3.keywords['SCF__D_CONVERGENCE'] == 1.e-11
            #assert plan3.keywords['E_CONVERGENCE'] == 1.e-10


def test_nbody_dimer_cbs():
    mol = psi4.geometry("""
                        He 0 0 -5
                        --
                        He 0 0 5
                        units au""")
    plan = task_planner("energy", "MP2/cc-pV[D,T]Z", mol, bsse_type="cp")

    assert isinstance(plan, ManyBodyComputer)
    # DGAS Note: This is wrong, for some reason cp is adding monomers in monomer basis.
    # See the CP builder in nbody:build_nbody_compute_list
    assert len(plan.task_list) == 5

    for k2, plan2 in plan.task_list.items():
        assert isinstance(plan2, CompositeComputer)
        assert len(plan2.task_list) == 2

        for i3, plan3 in enumerate(plan2.task_list):
            assert isinstance(plan3, AtomicComputer)
            assert plan3.basis == ["cc-pvdz", "cc-pvtz"][i3]
            assert plan3.method == "mp2"
            assert plan3.driver == "energy"


def test_nbody_dimer_cbs_gradient():
    mol = psi4.geometry("""
                        He 0 0 -5
                        --
                        He 0 0 5
                        units au""")
    plan = task_planner("gradient", "MP2/cc-pV[D,T]Z", mol, bsse_type="cp")

    assert isinstance(plan, ManyBodyComputer)
    # DGAS Note: This is wrong, for some reason cp is adding monomers in monomer basis.
    # See the CP builder in nbody:build_nbody_compute_list
    assert len(plan.task_list) == 5

    for k2, plan2 in plan.task_list.items():
        assert isinstance(plan2, CompositeComputer)
        assert len(plan2.task_list) == 4

        for i3, plan3 in enumerate(plan2.task_list):
            assert isinstance(plan3, AtomicComputer)
            assert plan3.basis == ["cc-pvtz", "cc-pvdz", "cc-pvtz", "cc-pvdz"][i3]
            assert plan3.method == ["hf", "mp2", "mp2", "hf"][i3]
            assert plan3.driver == "gradient"


@pytest.mark.parametrize("mtd, kw", [
    ('mp2', {'dertype': 0}),
    ('mp5', {}),
])
def test_nbody_dimer_cbs_gradient_1_0(mtd, kw):
    mol = psi4.geometry("""
                        He 0 0 -5
                        --
                        He 0 0 5
                        units au""")
    plan = task_planner("gradient", f"{mtd}/cc-pV[D,T]Z", mol, **kw, bsse_type="cp", findif_stencil_size=3, findif_step_size=0.005/math.sqrt(2/4.00260325413))

    displacements = {
        '0: -1': np.array([[ 0.    ,  0.    , -5.0025], [ 0.    ,  0.    ,  5.0025]]),
        '0: 1':  np.array([[ 0.    ,  0.    , -4.9975], [ 0.    ,  0.    ,  4.9975]]),
        'reference':  np.array([[ 0.    ,  0.    , -5.0], [ 0.    ,  0.    ,  5.0]]),
    }

    nbody_displacements = {
        "1_((2,), (2,))": {k: v[1] for k, v in displacements.items()},
        "1_((1,), (1,))": {k: v[0] for k, v in displacements.items()},
        "1_((1, 2), (1, 2))": displacements,
        "1_((1,), (1, 2))": displacements,
        "1_((2,), (1, 2))": displacements,
    }

    ghostiness = {
        "1_((2,), (2,))": (['He'], [True]),
        "1_((1,), (1,))": (['He'], [True]),
        "1_((1, 2), (1, 2))": (['He', 'He'], [True, True]),
        "1_((1,), (1, 2))": (['He', 'He'], [True, False]),
        "1_((2,), (1, 2))": (['He', 'He'], [False, True]),
    }

    assert isinstance(plan, ManyBodyComputer)
    # DGAS Note: This is wrong, for some reason cp is adding monomers in monomer basis.
    # See the CP builder in nbody:build_nbody_compute_list
    assert len(plan.task_list) == 5

    for k2, plan2 in plan.task_list.items():
        assert isinstance(plan2, FiniteDifferenceComputer)
        assert plan2.driver == "gradient"

        kmol = plan2.molecule.to_schema(dtype=2)
        assert kmol['symbols'] == ghostiness[k2][0]
        assert kmol['real'] == ghostiness[k2][1]

        for k3, plan3 in plan2.task_list.items():
            assert isinstance(plan3, CompositeComputer)
            assert len(plan3.task_list) == 2
            assert np.allclose(plan3.molecule.geometry().np, nbody_displacements[k2][k3])

            for i4, plan4 in enumerate(plan3.task_list):
                assert isinstance(plan4, AtomicComputer)
                assert plan4.basis == ["cc-pvdz", "cc-pvtz"][i4]
                assert plan4.method == mtd
                assert plan4.driver == "energy"


def test_cbs_extrapolation_delta():
    mol = psi4.geometry("He")
    plan = task_planner("energy", "MP2/cc-pV[D,T]Z + D:ccsd(t)/cc-pv[tq]z", mol)

    assert isinstance(plan, CompositeComputer)
    assert len(plan.task_list) == 3

    for i2, plan2 in enumerate(plan.task_list):
        assert isinstance(plan2, AtomicComputer)
        assert plan2.basis == ["cc-pvdz", "cc-pvtz", "cc-pvqz"][i2]
        assert plan2.method == ["mp2", "ccsd(t)", "ccsd(t)"][i2]
        assert plan2.driver == "energy"


def test_findif_1_1():
    mol = psi4.geometry("H\nH 1 2.0\nunits au")
    plan = task_planner("gradient", "MP2/cc-pVDZ", mol)

    assert isinstance(plan, AtomicComputer)
    assert plan.basis == "cc-pvdz"
    assert plan.method == "mp2"
    assert plan.driver == "gradient"
    # below are now back to optstash
    #assert plan.task_list[key].keywords['E_CONVERGENCE'] == 10
    #assert plan.task_list[key].keywords['SCF__E_CONVERGENCE'] == 10
    #assert plan.task_list[key].keywords['SCF__D_CONVERGENCE'] == 10


@pytest.mark.parametrize("mtd, kw", [
    ('mp2', {'dertype': 0}),
    ('mp5', {}),
])
def test_findif_1_0(mtd, kw):
    mol = psi4.geometry("H\nH 1 2.0\nunits au")
    plan = task_planner("gradient", f"{mtd}/cc-pVDZ", mol, **kw, findif_stencil_size=3, findif_step_size=0.005/math.sqrt(2/1.00782503223))

    displacements = {
        '0: -1': np.array([[ 0.    ,  0.    , -1.0025], [ 0.    ,  0.    ,  1.0025]]),
        '0: 1':  np.array([[ 0.    ,  0.    , -0.9975], [ 0.    ,  0.    ,  0.9975]]),
        'reference':  np.array([[ 0.    ,  0.    , -1.0], [ 0.    ,  0.    ,  1.0]]),
    }

    assert isinstance(plan, FiniteDifferenceComputer)
    assert len(plan.task_list) == 3

    for k2, plan2 in plan.task_list.items():
        assert isinstance(plan2, AtomicComputer)
        assert plan2.basis == "cc-pvdz"
        assert plan2.method == mtd
        assert plan2.driver == "energy"
        assert np.allclose(plan2.molecule.geometry().np, displacements[k2])
        assert plan2.keywords['SCF__E_CONVERGENCE'] == 1.e-10
        assert plan2.keywords['SCF__D_CONVERGENCE'] == 1.e-10
        assert plan2.keywords['E_CONVERGENCE'] == 1.e-8


def test_findif_2_1():
    mol = psi4.geometry("H\nH 1 2.0\nunits au")
    psi4.set_options({"E_CONVERGENCE": 6})
    plan = task_planner("hessian", "MP2/cc-pVDZ", mol, dertype=1, findif_stencil_size=3, findif_step_size=0.005/math.sqrt(2/1.00782503223))

    displacements = {
        '0: -1': np.array([[ 0.    ,  0.    , -1.0025], [ 0.    ,  0.    ,  1.0025]]),
        '0: 1':  np.array([[ 0.    ,  0.    , -0.9975], [ 0.    ,  0.    ,  0.9975]]),
        'reference':  np.array([[ 0.    ,  0.    , -1.0], [ 0.    ,  0.    ,  1.0]]),
    }

    assert isinstance(plan, FiniteDifferenceComputer)
    assert len(plan.task_list) == 3

    for k2, plan2 in plan.task_list.items():
        assert isinstance(plan2, AtomicComputer)
        assert plan2.basis == "cc-pvdz"
        assert plan2.method == "mp2"
        assert plan2.driver == "gradient"
        assert np.allclose(plan2.molecule.geometry().np, displacements[k2])
        assert plan2.keywords['SCF__D_CONVERGENCE'] == 1.e-10
        assert plan2.keywords['E_CONVERGENCE'] == 1.e-6


def test_findif_2_0():
    mol = psi4.geometry("H\nH 1 2.0\nunits au")
    psi4.set_options({"scf__E_CONVERGENCE": 6})
    plan = task_planner("hessian", "MP2/cc-pVDZ", mol, dertype=0, findif_stencil_size=5, findif_step_size=0.005/math.sqrt(2/1.00782503223))

    displacements = {
        '0: -2': np.array([[ 0.    ,  0.    , -1.0050], [ 0.    ,  0.    ,  1.0050]]),
        '0: 2':  np.array([[ 0.    ,  0.    , -0.9950], [ 0.    ,  0.    ,  0.9950]]),
        '0: -1': np.array([[ 0.    ,  0.    , -1.0025], [ 0.    ,  0.    ,  1.0025]]),
        '0: 1':  np.array([[ 0.    ,  0.    , -0.9975], [ 0.    ,  0.    ,  0.9975]]),
        'reference':  np.array([[ 0.    ,  0.    , -1.0], [ 0.    ,  0.    ,  1.0]]),
    }

    assert isinstance(plan, FiniteDifferenceComputer)
    assert len(plan.task_list) == 5

    for k2, plan2 in plan.task_list.items():
        assert isinstance(plan2, AtomicComputer)
        assert plan2.basis == "cc-pvdz"
        assert plan2.method == "mp2"
        assert plan2.driver == "energy"
        assert np.allclose(plan2.molecule.geometry().np, displacements[k2])
        assert plan2.keywords['SCF__E_CONVERGENCE'] == 1.e-6
        assert plan2.keywords['SCF__D_CONVERGENCE'] == 1.e-11
        assert plan2.keywords['E_CONVERGENCE'] == 1.e-10
