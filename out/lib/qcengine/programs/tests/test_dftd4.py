"""
Testing for the DFT-D4 program harness.

Most of the tests use mindless molecules for the diversity of element species
to test as much different interactions as possible.
"""

import numpy as np
import pytest
import qcelemental as qcel

import qcengine as qcng
from qcengine.testing import using


@using("dftd4")
def test_dftd4_task_b97m_m01():

    thr = 1.0e-8

    return_result = -0.025024986301735823

    atomic_input = qcel.models.AtomicInput(
        molecule=qcng.get_molecule("mindless-01"),
        model={"method": "b97m"},
        driver="energy",
    )

    atomic_result = qcng.compute(atomic_input, "dftd4")

    print(atomic_result.return_result)
    assert atomic_result.success
    assert pytest.approx(atomic_result.return_result, abs=thr) == return_result


@using("dftd4")
def test_dftd4_task_tpss_m02():

    thr = 1.0e-8

    return_result = np.array(
        [
            [-1.1006300442501966e-04, -1.9882677409412426e-04, +9.7969330873203453e-05],
            [-9.0564826517969396e-04, -1.1221656220960802e-04, +7.4420608164684399e-04],
            [-8.7360785425552928e-05, +2.8146769568946061e-04, +5.8528720048366145e-04],
            [-1.7602264871464991e-04, -3.0981396137083546e-04, -1.1315671174001885e-03],
            [+6.2854950429890028e-05, +3.1662886959557454e-04, +1.8192682295895164e-04],
            [-1.2033574062747143e-04, -5.8117271931259600e-05, -2.4103707001727414e-04],
            [-9.6727467530178130e-05, -1.2937750080676020e-04, -3.0565838265275271e-05],
            [-1.2139762550383215e-05, +2.2334694600929773e-05, +6.7681766980205114e-05],
            [+2.9131236664176958e-04, -2.3159856865643386e-04, +3.9334087528658355e-04],
            [-8.9512907514913525e-05, +9.3857773376367446e-06, -2.6351331811570630e-04],
            [+3.9621947274011518e-04, -5.8926113817768316e-04, +1.8422719633870202e-04],
            [-1.6466329096860149e-04, +1.3720763821983788e-04, -8.3839852301748649e-05],
            [-3.2173735474065368e-04, +7.2068238657310593e-04, -3.1811972689724591e-04],
            [+2.1298820577529375e-04, -9.6810383911168402e-05, +1.0273702281858908e-04],
            [-4.4599095547508355e-04, -1.8320769460464104e-04, +1.8653129262487302e-04],
            [+1.5668271875651332e-03, +4.2152279374596904e-04, -4.7526466701417510e-04],
        ]
    )

    atomic_input = qcel.models.AtomicInput(
        molecule=qcng.get_molecule("mindless-02"),
        model={"method": ""},
        keywords={
            "params_tweaks": {
                "s8": 1.76596355,
                "a1": 0.42822303,
                "a2": 4.54257102,
            },
        },
        driver="gradient",
    )

    atomic_result = qcng.compute(atomic_input, "dftd4")

    assert atomic_result.success
    assert pytest.approx(atomic_result.return_result, abs=thr) == return_result


@using("dftd4")
def test_dftd4_task_r2scan_m03():

    thr = 1.0e-8

    return_result = np.array(
        [
            [-1.0184389290110640e-05, -1.2902131714706869e-05, -4.7991389860088184e-05],
            [+3.9008038585852290e-05, -8.6543337237836941e-05, +1.0086268708372256e-04],
            [-2.1792111331028049e-05, +2.8414264003578678e-05, +1.8916279332120321e-05],
            [-3.8847190487495055e-06, +6.4069296782535270e-06, -1.7362611928711208e-05],
            [-1.5358378150451566e-05, -5.7611513737809471e-06, -1.9871764856476196e-05],
            [-3.7225957329162548e-05, -1.5687552833166010e-06, -3.4717349142275990e-05],
            [+3.2408291983641028e-07, -6.7683941908276048e-05, -6.1602680430568183e-05],
            [-8.1819149671629329e-06, -9.5470867672313039e-06, -1.7007294079721094e-05],
            [-1.7245516321582448e-05, +1.6416848978506148e-05, +1.6005589814916458e-05],
            [-3.0116095523698349e-06, +7.3382310356380993e-07, -1.3652036487357768e-05],
            [+1.0378663330288130e-04, +4.4180555120489236e-05, +3.2576724474326750e-05],
            [-2.8021453980080289e-05, +3.9008636389849655e-05, -1.2083469461741150e-04],
            [-3.0452136766248191e-05, +7.0712047105794630e-05, +8.3682948356528971e-05],
            [+1.6902554427499383e-06, -2.9147110970219634e-05, +3.6402009185785855e-05],
            [+4.6966237842630462e-05, -1.2554523547742754e-05, +3.0600183501647347e-05],
            [-1.6417061357004406e-05, +1.9834934423075429e-05, +1.3993399653561848e-05],
        ]
    )

    atomic_input = qcel.models.AtomicInput(
        molecule=qcng.get_molecule("mindless-03"),
        keywords={"level_hint": "D4"},
        driver="gradient",
        model={"method": "r2scan"},
    )

    atomic_result = qcng.compute(atomic_input, "dftd4")

    print(atomic_result.return_result)
    assert atomic_result.success
    assert pytest.approx(atomic_result.return_result, abs=thr) == return_result


@using("dftd4")
def test_dftd4_task_unknown_method():

    atomic_input = qcel.models.AtomicInput(
        molecule=qcng.get_molecule("water"),
        keywords={"level_hint": "D4"},
        model={"method": "non-existent-method"},
        driver="energy",
    )
    error = qcel.models.ComputeError(
        error_type="input error", error_message="Functional 'non-existent-method' not known"
    )

    atomic_result = qcng.compute(atomic_input, "dftd4")

    print(atomic_result.error)
    assert not atomic_result.success
    assert atomic_result.error == error


@using("dftd4")
def test_dftd4_task_cold_fusion():

    atomic_input = qcel.models.AtomicInput(
        molecule={
            "symbols": ["Li", "Li", "Li", "Li"],
            "geometry": [
                [-1.58746019997201, +1.58746019997201, +1.58746019997201],
                [-1.58746019997201, +1.58746019997201, +1.58746019997201],
                [-1.58746019997201, -1.58746019997201, -1.58746019997201],
                [+1.58746019997201, +1.58746019997201, -1.58746019997201],
            ],
            "validated": True,  # Force a nuclear fusion input, to make dftd4 fail
        },
        keywords={"level_hint": "D4"},
        model={"method": "pbe"},
        driver="energy",
    )
    error = qcel.models.ComputeError(
        error_type="input error",
        error_message="Too close interatomic distances found",
    )

    atomic_result = qcng.compute(atomic_input, "dftd4")

    print(atomic_result.error)
    assert not atomic_result.success
    assert atomic_result.error == error
