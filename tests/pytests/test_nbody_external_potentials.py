import numpy as np
import pytest

import psi4
from psi4.driver.driver_nbody import ManyBodyComputer
from psi4.driver.p4util.exceptions import ValidationError
from psi4.driver.task_planner import task_planner
from qcmanybody.utils import delabeler

from addons import uusing


pytestmark = [pytest.mark.psi, pytest.mark.api]


@pytest.fixture(scope="module")
def water_cluster():
    trimer = psi4.geometry(
        """
        0 1
        H  0.0290 -1.1199 -1.5243
        O  0.9481 -1.3990 -1.3587
        H  1.4371 -0.5588 -1.3099
        --
        0 1
        H  1.0088 -1.5240  0.5086
        O  1.0209 -1.1732  1.4270
        H  1.5864 -0.3901  1.3101
        --
        0 1
        H -1.0231  1.6243 -0.8743
        O -0.5806  2.0297 -0.1111
        H -0.9480  1.5096  0.6281
        symmetry c1
        no_reorient
        no_com
        """
    )
    dimer = trimer.extract_subsets([1, 2])
    # fmt: off
    b_external_potential = [
        [0.417, (np.array([-2.5628, -0.8269, -1.6696]) / psi4.constants.bohr2angstroms).tolist()],
        [-0.834, (np.array([-1.7899, -0.4027, -1.2768]) / psi4.constants.bohr2angstroms).tolist()],
        [0.417, (np.array([-1.8988, -0.4993, -0.3072]) / psi4.constants.bohr2angstroms).tolist()],
    ]
    # fmt: on
    return trimer, dimer, b_external_potential


@uusing("qcmanybody")
def test_nbody_external_potentials_nocp_dimer(water_cluster):
    """Verify fragment potentials yield the manually assembled non-CP interaction energy."""
    _, dimer, b_external_potential = water_cluster
    psi4.set_options(
        {
            "basis": "sto-3g",
            "scf_type": "pk",
            "e_convergence": 10,
            "d_convergence": 10,
        }
    )

    mbe_ie = psi4.energy(
        "hf",
        molecule=dimer,
        bsse_type="nocp",
        external_potentials={2: b_external_potential},
    )

    monomer_a = psi4.core.Molecule.from_schema(
        dimer.extract_subsets(1).to_schema(dtype=2)
    )
    monomer_b = psi4.core.Molecule.from_schema(
        dimer.extract_subsets(2).to_schema(dtype=2)
    )
    e_a = psi4.energy("hf", molecule=monomer_a)
    e_b = psi4.energy(
        "hf", molecule=monomer_b, external_potentials=b_external_potential
    )
    supersystem = psi4.core.Molecule.from_schema(dimer.to_schema(dtype=2))
    e_ab = psi4.energy(
        "hf", molecule=supersystem, external_potentials=b_external_potential
    )

    assert mbe_ie == pytest.approx(e_ab - e_a - e_b, abs=1.0e-10)


@uusing("qcmanybody")
def test_nbody_external_potentials_cp_dimer(water_cluster):
    """Verify fragment potentials propagate through counterpoise-corrected calculations."""
    _, dimer, b_external_potential = water_cluster
    psi4.set_options(
        {
            "basis": "sto-3g",
            "scf_type": "pk",
            "e_convergence": 10,
            "d_convergence": 10,
        }
    )

    cp_ie = psi4.energy(
        "hf",
        molecule=dimer,
        bsse_type="cp",
        external_potentials={2: b_external_potential},
    )

    assert cp_ie == pytest.approx(-0.005286696523, abs=1.0e-10)


@uusing("qcmanybody")
def test_nbody_external_potentials_trimer_plan(water_cluster):
    """Verify planning assigns potentials only to tasks containing the target fragment."""
    trimer, _, b_external_potential = water_cluster
    plan = task_planner(
        "energy",
        "hf/sto-3g",
        trimer,
        bsse_type=["nocp", "cp"],
        external_potentials={2: b_external_potential},
    )

    assert isinstance(plan, ManyBodyComputer)
    real_fragments_seen = set()
    for label, component in plan.task_list.items():
        _, real_fragments, _ = delabeler(label)
        real_fragments = tuple(real_fragments)
        real_fragments_seen.add(real_fragments)
        component_potentials = component.keywords["function_kwargs"].get(
            "external_potentials"
        )
        if 2 in real_fragments:
            assert component_potentials == b_external_potential
        else:
            assert component_potentials is None

    assert real_fragments_seen == {(1,), (2,), (3,), (1, 2), (1, 3), (2, 3), (1, 2, 3)}

    with pytest.raises(ValidationError, match="1-indexed fragment integers"):
        task_planner(
            "energy",
            "hf/sto-3g",
            trimer,
            bsse_type="nocp",
            external_potentials={0: b_external_potential},
        )
