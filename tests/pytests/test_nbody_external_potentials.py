import numpy as np
import pytest

import psi4
from psi4.driver.driver_nbody import ManyBodyComputer
from psi4.driver.p4util.exceptions import ValidationError
from psi4.driver.task_planner import task_planner

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

    assert cp_ie == pytest.approx(-0.005286696523, abs=1.0e-6)


@uusing("qcmanybody")
def test_nbody_external_potentials_trimer_plan(water_cluster):
    """Verify planning assigns potentials only to tasks containing the target fragment."""
    from qcmanybody.utils import delabeler

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


@uusing("qcmanybody")
@pytest.mark.parametrize(
    "spelling",
    [
        pytest.param("points", id="mode-keys"),
        pytest.param("C", id="sapt-fragment-keys"),
    ],
)
def test_nbody_global_external_potential_dict_plan(water_cluster, spelling):
    """Verify whole-molecule dict spellings are not mistaken for a fragment mapping."""
    trimer, _, b_external_potential = water_cluster
    global_potential = {spelling: b_external_potential}

    plan = task_planner(
        "energy",
        "hf/sto-3g",
        trimer,
        bsse_type="nocp",
        external_potentials=global_potential,
    )

    assert isinstance(plan, ManyBodyComputer)
    assert plan.task_list
    for component in plan.task_list.values():
        assert component.keywords["function_kwargs"]["external_potentials"] == global_potential


@uusing("qcmanybody")
def test_nbody_global_external_potential_dict_energy(water_cluster):
    """Verify a whole-molecule dict spelling gives the same MBE energy as the equivalent flat list."""
    _, dimer, b_external_potential = water_cluster
    psi4.set_options(
        {
            "basis": "sto-3g",
            "scf_type": "pk",
            "e_convergence": 10,
            "d_convergence": 10,
        }
    )

    dict_ie = psi4.energy(
        "hf",
        molecule=dimer,
        bsse_type="nocp",
        external_potentials={"points": b_external_potential},
    )
    list_ie = psi4.energy(
        "hf",
        molecule=dimer,
        bsse_type="nocp",
        external_potentials=b_external_potential,
    )

    assert dict_ie == pytest.approx(list_ie, abs=1.0e-10)


@uusing("qcmanybody")
def test_nbody_external_potentials_with_embedding_charges(water_cluster, monkeypatch):
    """Verify embedding charges add to, rather than replace, fragment-scoped potentials."""
    from qcmanybody.utils import delabeler

    monkeypatch.setenv("QCMANYBODY_EMBEDDING_CHARGES", "1")

    _, dimer, b_external_potential = water_cluster
    plan = task_planner(
        "energy",
        "hf/sto-3g",
        dimer,
        bsse_type="nocp",
        return_total_data=True,
        embedding_charges={1: [0.417, -0.834, 0.417], 2: [0.417, -0.834, 0.417]},
        external_potentials={2: b_external_potential},
    )

    assert isinstance(plan, ManyBodyComputer)
    flat_b = [[charge, *position] for charge, position in b_external_potential]
    real_fragments_seen = set()
    for label, component in plan.task_list.items():
        _, real_fragments, _ = delabeler(label)
        real_fragments = tuple(real_fragments)
        real_fragments_seen.add(real_fragments)
        component_potentials = component.keywords["function_kwargs"]["external_potentials"]

        # whatever spelling the planner chose, check the point charges the component will actually receive
        validated = psi4.driver.p4util.validate_external_potential(component_potentials)
        assert set(validated) == {"C"}
        assert set(validated["C"]) == {"points"}
        points = [list(point) for point in validated["C"]["points"]]

        # one embedding charge per atom of each fragment not in the component's basis
        nembedded = 3 * (2 - len(real_fragments))
        if 2 in real_fragments:
            assert len(points) == nembedded + len(flat_b)
            assert points[nembedded:] == flat_b
        else:
            assert len(points) == nembedded
            for potential in flat_b:
                assert potential not in points

    assert real_fragments_seen == {(1,), (2,), (1, 2)}


@uusing("qcmanybody")
@pytest.mark.parametrize(
    "c_partition",
    [
        pytest.param(None, id="absent-c"),
        pytest.param("points-list", id="existing-c-points-list"),
        pytest.param("mode-dict", id="existing-c-points-and-diffuse"),
    ],
)
def test_nbody_sapt_dict_external_potentials_with_embedding_charges(
    water_cluster, monkeypatch, c_partition
):
    """Verify embedding charges join the C partition without disturbing A, B, or C diffuse content."""
    from qcmanybody.utils import delabeler

    monkeypatch.setenv("QCMANYBODY_EMBEDDING_CHARGES", "1")

    _, dimer, b_external_potential = water_cluster
    pot_a = [[0.2, [1.0, 2.0, 3.0]]]
    pot_b = b_external_potential
    c_points = [[-0.3, [4.0, 5.0, 6.0]]]
    c_diffuse = [[0.1, 7.0, 8.0, 9.0, 0.5]]

    global_potential = {"A": pot_a, "B": pot_b}
    if c_partition == "points-list":
        global_potential["C"] = c_points
    elif c_partition == "mode-dict":
        global_potential["C"] = {"points": c_points, "diffuse": c_diffuse}

    plan = task_planner(
        "energy",
        "hf/sto-3g",
        dimer,
        bsse_type="nocp",
        return_total_data=True,
        embedding_charges={1: [0.417, -0.834, 0.417], 2: [0.417, -0.834, 0.417]},
        external_potentials=global_potential,
    )

    assert isinstance(plan, ManyBodyComputer)
    expected_c_points = [] if c_partition is None else c_points
    real_fragments_seen = set()
    for label, component in plan.task_list.items():
        _, real_fragments, _ = delabeler(label)
        real_fragments_seen.add(tuple(real_fragments))
        component_potentials = component.keywords["function_kwargs"]["external_potentials"]

        assert component_potentials["A"] == pot_a
        assert component_potentials["B"] == pot_b

        # one embedding charge per atom of each fragment not in the component's basis
        nembedded = 3 * (2 - len(real_fragments))
        if nembedded == 0:
            assert component_potentials == global_potential
            continue

        assert set(component_potentials) == {"A", "B", "C"}
        c_spec = component_potentials["C"]
        assert set(c_spec) == ({"points", "diffuse"} if c_partition == "mode-dict" else {"points"})
        assert len(c_spec["points"]) == nembedded + len(expected_c_points)
        assert all(len(point) == 4 for point in c_spec["points"])
        if c_partition == "mode-dict":
            assert [list(row) for row in c_spec["diffuse"]] == c_diffuse

        # the merged structure remains acceptable to the downstream consumer
        assert psi4.driver.p4util.validate_external_potential(component_potentials)

    assert real_fragments_seen == {(1,), (2,), (1, 2)}


@uusing("qcmanybody")
@pytest.mark.parametrize(
    "spelling",
    [
        pytest.param("matrix-only", id="positional-matrix-only"),
        pytest.param("diffuse-only", id="positional-diffuse-only"),
        pytest.param("points-diffuse", id="positional-points-diffuse"),
        pytest.param("all-three", id="positional-points-diffuse-matrix"),
        pytest.param("nested-points", id="positional-points-only"),
        pytest.param("nested-single-row", id="positional-single-row-points"),
        pytest.param("nested-single-row-xyz", id="positional-single-row-points-nested-xyz"),
        pytest.param("nested-single-row-array", id="positional-single-row-points-ndarray"),
        pytest.param("points-trailing-none", id="positional-points-trailing-none"),
    ],
)
def test_nbody_positional_external_potentials_with_embedding_charges(
    water_cluster, monkeypatch, spelling
):
    """Verify positional [points, diffuse, matrix] potentials survive the embedding-charge merge."""
    from qcmanybody.utils import delabeler

    monkeypatch.setenv("QCMANYBODY_EMBEDDING_CHARGES", "1")

    _, dimer, _ = water_cluster
    points = [[0.5, 0.0, 0.0, 1.0], [-0.5, 0.0, 0.0, -1.0]]
    diffuse = [[0.3, 1.0, 1.0, 1.0, 0.5]]
    matrix = np.eye(5).tolist()

    single_row = [[-0.5, 1.0, 1.0, 0.0]]
    single_row_xyz = [[-0.5, [1.0, 1.0, 0.0]]]

    global_potential, expected_points, expected_modes = {
        "matrix-only": ([None, None, matrix], [], {"points", "matrix"}),
        "diffuse-only": ([None, diffuse], [], {"points", "diffuse"}),
        "points-diffuse": ([points, diffuse], points, {"points", "diffuse"}),
        "all-three": ([points, diffuse, matrix], points, {"points", "diffuse", "matrix"}),
        "nested-points": ([points], points, {"points"}),
        "nested-single-row": ([single_row], single_row, {"points"}),
        "nested-single-row-xyz": ([single_row_xyz], single_row, {"points"}),
        "nested-single-row-array": ([np.array(single_row)], single_row, {"points"}),
        "points-trailing-none": ([points, None], points, {"points"}),
    }[spelling]

    plan = task_planner(
        "energy",
        "hf/sto-3g",
        dimer,
        bsse_type="nocp",
        return_total_data=True,
        embedding_charges={1: [0.417, -0.834, 0.417], 2: [0.417, -0.834, 0.417]},
        external_potentials=global_potential,
    )

    assert isinstance(plan, ManyBodyComputer)
    real_fragments_seen = set()
    for label, component in plan.task_list.items():
        _, real_fragments, _ = delabeler(label)
        real_fragments_seen.add(tuple(real_fragments))
        component_potentials = component.keywords["function_kwargs"]["external_potentials"]

        # one embedding charge per atom of each fragment not in the component's basis
        nembedded = 3 * (2 - len(real_fragments))
        if nembedded == 0:
            assert psi4.driver.p4util.validate_external_potential(
                component_potentials
            ) == psi4.driver.p4util.validate_external_potential(global_potential)
            continue

        assert set(component_potentials) == expected_modes
        assert len(component_potentials["points"]) == nembedded + len(expected_points)
        assert all(len(point) == 4 for point in component_potentials["points"])
        assert [list(point) for point in component_potentials["points"][nembedded:]] == expected_points
        if "diffuse" in expected_modes:
            assert [list(row) for row in component_potentials["diffuse"]] == diffuse
        if "matrix" in expected_modes:
            assert component_potentials["matrix"] == matrix

        # the merged structure remains acceptable to the downstream consumer
        assert psi4.driver.p4util.validate_external_potential(component_potentials)

    assert real_fragments_seen == {(1,), (2,), (1, 2)}


@uusing("qcmanybody")
@pytest.mark.parametrize(
    "value",
    [
        pytest.param([None, [[0.3, 1.0, 1.0, 1.0, 0.5]]], id="positional-diffuse"),
        pytest.param([None, None, np.eye(5).tolist()], id="positional-matrix"),
        pytest.param({"points": [[0.4, [1.0, 2.0, 3.0]]]}, id="mode-dict"),
        pytest.param({"diffuse": [[0.3, 1.0, 1.0, 1.0, 0.5]]}, id="mode-dict-diffuse"),
        pytest.param(np.eye(5).tolist(), id="bare-matrix"),
        pytest.param(0.4, id="scalar"),
    ],
)
def test_nbody_fragment_potential_rejects_non_point_charge_values(water_cluster, value):
    """Verify a per-fragment value that is not a flat point-charge list is refused, not garbled."""
    _, dimer, b_external_potential = water_cluster

    with pytest.raises(ValidationError, match="flat point-charge list"):
        task_planner(
            "energy",
            "hf/sto-3g",
            dimer,
            bsse_type="nocp",
            external_potentials={1: b_external_potential, 2: value},
        )


@uusing("qcmanybody")
def test_nbody_fragment_potentials_combine_across_fragments(water_cluster):
    """Verify point charges from several real fragments concatenate into one valid potential."""
    from qcmanybody.utils import delabeler

    trimer, _, b_external_potential = water_cluster
    pot_1 = [[0.4, [1.0, 2.0, 3.0]]]
    pot_2 = np.array([[0.5, 4.0, 5.0, 6.0], [-0.5, 7.0, 8.0, 9.0]])

    plan = task_planner(
        "energy",
        "hf/sto-3g",
        trimer,
        bsse_type="nocp",
        external_potentials={1: pot_1, 2: pot_2},
    )

    per_fragment = {1: 1, 2: 2, 3: 0}
    for label, component in plan.task_list.items():
        _, real_fragments, _ = delabeler(label)
        expected = sum(per_fragment[ifr] for ifr in real_fragments)
        component_potentials = component.keywords["function_kwargs"].get("external_potentials")
        if expected == 0:
            assert component_potentials is None
            continue

        validated = psi4.driver.p4util.validate_external_potential(component_potentials)
        assert set(validated) == {"C"}
        assert set(validated["C"]) == {"points"}
        assert len(validated["C"]["points"]) == expected


@uusing("qcmanybody")
def test_fragment_keyed_external_potentials_require_bsse_type(water_cluster):
    """Verify an integer-keyed mapping without bsse_type reports the missing many-body request."""
    trimer, _, b_external_potential = water_cluster
    psi4.set_options({"basis": "sto-3g", "scf_type": "pk"})

    with pytest.raises(ValidationError, match="1-indexed fragment mapping"):
        psi4.energy("hf", molecule=trimer, external_potentials={2: b_external_potential})

    with pytest.raises(ValidationError, match="1-indexed fragment mapping"):
        psi4.driver.p4util.validate_external_potential({2: b_external_potential})
