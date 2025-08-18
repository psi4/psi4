import pytest
import numpy as np
import psi4

from qcelemental.testing import compare_recursive

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.extern]

@pytest.mark.parametrize(
    "frame",
    # * this parameter alters the input molecule by fix_com and fix_orientation to imitate user signaling frame matters or not.
    [
        pytest.param("fixed"),  # fix_=True (no_com/no_reorient)
        pytest.param("free"),  # fix_=False (def)
    ],
)
@pytest.mark.parametrize("deriv", ["0_0", "1_1",
    # pytest.param("1_0", marks=pytest.mark.long),  # not really long but 6x rest of test so no need to run routinely
])
@pytest.mark.parametrize("molmode", [ "geometry", "from_arrays", "from_string", "from_dict", "from_schema_2",])
def test_nopotential(frame, deriv, molmode):
    #! extern5 test with broader checking:
    #! External potential sanity check with 0 charge far away
    #! Checks if all units behave the same and energy is same as no potential

    driver, dertype = deriv.split("_")
    driver = {"0": psi4.energy, "1": psi4.gradient}[driver]

    b2a=0.529177249
    smol_bohr = """
        O  -1.47172427  0.          2.1404605
        H  -1.2598463   1.44393774  3.22442245
        H  -1.2598463  -1.44393774  3.22442056
        units au
        """
    smol_ang = """
        O   -0.778803000000  0.000000000000  1.132683000000
        H   -0.666682000000  0.764099000000  1.706291000000
        H   -0.666682000000  -0.764099000000  1.706290000000
        units ang
        """

    if frame == "fixed":
        fixlines = """
        symmetry c1
        nocom
        noreorient
        """
        smol_bohr += fixlines
        smol_ang += fixlines

    if molmode == "geometry":
        molecule_bohr = psi4.geometry(smol_bohr)
        molecule_ang = psi4.geometry(smol_ang)

    elif molmode == "from_arrays":
        elements = ["O","H","H"]
        # Coordinates added in angstrom
        coords = np.array([[ -0.778803000000 , 0.000000000000,  1.132683000000],
         [ -0.666682000000,  0.764099000000,  1.706291000000],
         [ -0.666682000000,  -0.764099000000 , 1.706290000000]])

        fixargs = {"fix_symmetry": "c1", "fix_com": True, "fix_orientation": True} if frame == "fixed" else {}
        molecule_bohr = psi4.core.Molecule.from_arrays(geom=coords/b2a, elem=elements, **fixargs, units="Bohr")
        molecule_ang  = psi4.core.Molecule.from_arrays(geom=coords, elem=elements, **fixargs)

    elif molmode == "from_string":
        molecule_bohr = psi4.core.Molecule.from_string(smol_bohr)
        molecule_ang = psi4.core.Molecule.from_string(smol_ang)

    elif molmode  == "from_dict":
        dict_bohr = {'provenance': {'creator': 'QCElemental', 'routine': 'qcelemental.molparse.from_string', 'version': '0.29.0'}, 'units': 'Bohr', 'geom': [-1.47172427,  0.        ,  2.1404605 , -1.2598463 ,  1.44393774,
        3.22442245, -1.2598463 , -1.44393774,  3.22442056], 'elea': [16,  1,  1], 'elez': [8, 1, 1], 'elem': ['O', 'H', 'H'], 'mass': [15.99491462,  1.00782503,  1.00782503], 'real': [ True,  True,  True], 'elbl': ['', '', ''], 'fragment_separators': [], 'molecular_charge': 0.0, 'fragment_charges': [0.0], 'molecular_multiplicity': 1, 'fragment_multiplicities': [1]}
        dict_ang = {'provenance': {'creator': 'QCElemental', 'routine': 'qcelemental.molparse.from_string', 'version': '0.29.0'}, 'units': 'Angstrom', 'geom': [-0.778803,  0.      ,  1.132683, -0.666682,  0.764099,  1.706291,
       -0.666682, -0.764099,  1.70629 ], 'elea': [16,  1,  1], 'elez': [8, 1, 1], 'elem': ['O', 'H', 'H'], 'mass': [15.99491462,  1.00782503,  1.00782503], 'real': [ True,  True,  True], 'elbl': ['', '', ''], 'fragment_separators': [], 'molecular_charge': 0.0, 'fragment_charges': [0.0], 'molecular_multiplicity': 1, 'fragment_multiplicities': [1]}

        if frame == "fixed":
            fixargs = {"fix_symmetry": "c1", "fix_com": True, "fix_orientation": True}
        else:
            fixargs = {"fix_com": False, "fix_orientation": False}
        molecule_bohr = psi4.core.Molecule.from_dict({**dict_bohr, **fixargs})
        molecule_ang = psi4.core.Molecule.from_dict({**dict_ang, **fixargs})

    elif molmode == "from_schema_2":
        qcsk_bohr = {'symbols': ['O', 'H', 'H'], 'geometry': [-1.47172427, 0.0, 2.1404605, -1.2598463, 1.44393774, 3.22442245, -1.2598463, -1.44393774, 3.22442056], 'masses': [15.99491461957, 1.00782503223, 1.00782503223], 'atomic_numbers': [8, 1, 1], 'mass_numbers': [16, 1, 1], 'atom_labels': ['', '', ''], 'name': 'H2O', 'molecular_charge': 0.0, 'molecular_multiplicity': 1, 'real': [True, True, True], 'fragments': [[0, 1, 2]], 'fragment_charges': [0.0], 'fragment_multiplicities': [1], 'provenance': {'creator': 'QCElemental', 'routine': 'qcelemental.molparse.from_string', 'version': '0.29.0'}, 'schema_name': 'qcschema_molecule', 'schema_version': 2}
        # note ang same as bohr in QCSchema v2
        qcsk_ang = {'validated': True, 'symbols': ['O', 'H', 'H'], 'geometry': [-1.471724375684933, 0.0, 2.1404606569619493, -1.2598463927724757, 1.443937842736201, 3.224422680333563, -1.2598463927724757, -1.443937842736201, 3.2244207906074376], 'masses': [15.99491461957, 1.00782503223, 1.00782503223], 'atomic_numbers': [8, 1, 1], 'mass_numbers': [16, 1, 1], 'atom_labels': ['', '', ''], 'name': 'H2O', 'molecular_charge': 0.0, 'molecular_multiplicity': 1, 'real': [True, True, True], 'fragments': [[0, 1, 2]], 'fragment_charges': [0.0], 'fragment_multiplicities': [1], 'provenance': {'creator': 'QCElemental', 'routine': 'qcelemental.molparse.from_string', 'version': '0.29.0'}, 'schema_name': 'qcschema_molecule', 'schema_version': 2}

        if frame == "fixed":
            fixargs = {"fix_symmetry": "c1", "fix_com": True, "fix_orientation": True}
        else:
            fixargs = {"fix_com": False, "fix_orientation": False}
        molecule_bohr = psi4.core.Molecule.from_schema({**qcsk_bohr, **fixargs})
        molecule_ang = psi4.core.Molecule.from_schema({**qcsk_ang, **fixargs})

    nre = 9.14756
    psi4.compare_values(nre, molecule_bohr.nuclear_repulsion_energy(), 5, f"{molmode}: [1b] Bohr geometry NRE")
    psi4.compare_values(nre, molecule_ang.nuclear_repulsion_energy(), 5, f"{molmode}: [1a] Ang geometry NRE")

    external_potentials = [[0.00, np.array([10.0,10.0,10.0]) / b2a]]

    refs = {psi4.energy: -74.96341862235225, psi4.gradient: np.array([[ 1.09256884e-02, -7.29151482e-07,  5.58954770e-02],
 [-5.46279071e-03, -1.75246287e-02, -2.79474763e-02],
 [-5.46289770e-03,  1.75253578e-02, -2.79480007e-02]])}

    psi4.set_options( {
        "scf_type": "df",
        "d_convergence": 12,
        "basis": "STO-3G",
        "print": 0,
        "debug": 0,
    })
    if psi4.core.get_option("scf", "orbital_optimizer_package") != "INTERNAL":
        psi4.set_options({"e_convergence": 9, "d_convergence": 3e-8})

    atol = 3.e-6 if deriv == "1_0" else 1.e-6

    ret_bohr_pure = driver('scf', molecule=molecule_bohr, dertype=dertype)
    ene_bohr_pure = psi4.variable("CURRENT ENERGY")
    ret_ang_pure = driver('scf', molecule=molecule_ang, dertype=dertype)
    ene_ang_pure = psi4.variable("CURRENT ENERGY")

    # check energies (always available)
    psi4.compare_values(refs[psi4.energy], ene_bohr_pure, atol, f"[2] {molmode} {deriv}: Bohr geometry, no charges vs reference equality")
    psi4.compare_values(ene_ang_pure, ene_bohr_pure, atol, f"[3] {molmode} {deriv}: No charges, Bohr vs Angstrom geometry energy equality")

    if molmode == "from_dict" and deriv in ["0_0", "1_1"]:
        with pytest.raises(psi4.ValidationError) as err:
            driver('scf', molecule=molecule_bohr, external_potentials=external_potentials, dertype=dertype)
            driver('scf', molecule=molecule_ang, external_potentials=external_potentials, dertype=dertype)
        assert "Set no_com/no_reorient/symmetry c1 by hand" in str(err.value)
        return

    ret_bohr_charges = driver('scf', molecule=molecule_bohr, external_potentials=external_potentials, dertype=dertype)
    ene_bohr_charges = psi4.variable("CURRENT ENERGY")
    ret_ang_charges = driver('scf', molecule=molecule_ang, external_potentials=external_potentials, dertype=dertype)
    ene_ang_charges = psi4.variable("CURRENT ENERGY")

    # check energies continued
    psi4.compare_values(ene_bohr_charges, ene_bohr_pure, atol, f"[4] {molmode} {deriv}: Bohr geometry, charges vs no charges energy equality")
    psi4.compare_values(ene_ang_charges, ene_ang_pure, atol, f"[5] {molmode} {deriv}: Angstrom geometry, charges vs no charges energy equality")

    # check returns (redundant for driver=energy; new for gradients)
    # * driver should be imposing fix_*=True when EP specified
    charges_and_pure_frames_should_match = (driver == psi4.energy) or (frame == "fixed")
    if charges_and_pure_frames_should_match:
        psi4.compare_values(refs[driver], ret_bohr_pure, atol, f"[6] {molmode} {deriv}: Bohr geometry, no charges vs reference equality")
    psi4.compare_values(ret_ang_pure, ret_bohr_pure, atol, f"[7] {molmode} {deriv}: No charges, Bohr vs Angstrom geometry return equality")
    psi4.compare_values(ret_ang_charges, ret_bohr_charges, atol, f"[8] {molmode} {deriv}: With charges, Bohr vs Angstrom geometry return equality")
    if charges_and_pure_frames_should_match:
        psi4.compare_values(ret_bohr_charges, ret_bohr_pure, atol, f"[9] {molmode} {deriv}: Bohr geometry, charges vs no charges return equality")
        psi4.compare_values(ret_ang_charges, ret_ang_pure, atol, f"[10] {molmode} {deriv}: Angstrom geometry, charges vs no charges return equality")


_qxyz1a = [ [-0.5, 1.0, 1.0, 0.0] ]
_qxyz1b = [ [-0.5, [1.0, 1.0, 0.0]] ]
_ans1 = {"C": {"points": _qxyz1a}}
_ans2 = {"B": {"points": _qxyz1a}}
_qxyz4a = [ [-0.5, 1.0, 1.0, 0.0],
            [ 0.5, 0.0, 1.0, 0.0],
            [ 7.7, 0.0, 0.0, 1.0],
            [-7.7, 0.0, 1.0, 1.0] ]
_qxyz4b = [ [-0.5, [1.0, 1.0, 0.0]],
            [ 0.5, [0.0, 1.0, 0.0]],
            [ 7.7, [0.0, 0.0, 1.0]],
            [-7.7, [0.0, 1.0, 1.0]] ]
_ans3 = {"C": {"points": _qxyz4a}}
_ans4 = {"B": {"points": _qxyz4a}}
_qxyzw4a = [ [-0.5, 1.0, 1.0, 0.0, 0.8],
             [ 0.5, 0.0, 1.0, 0.0, 0.7],
             [ 7.7, 0.0, 0.0, 1.0, 0.6],
             [-7.7, 0.0, 1.0, 1.0, 0.5] ]
_qxyzw4b = [ [-0.5, [1.0, 1.0, 0.0], 0.8],
             [ 0.5, [0.0, 1.0, 0.0], 0.7],
             [ 7.7, [0.0, 0.0, 1.0], 0.6],
             [-7.7, [0.0, 1.0, 1.0], 0.5] ]
_ans5 = {"C": {"diffuse": _qxyzw4a}}
_ans6 = {"B": {"diffuse": _qxyzw4a}}
_ans7 = {"C": {"points": _qxyz4a, "diffuse": _qxyzw4a}}
_ans8 = {"B": {"points": _qxyz4a, "diffuse": _qxyzw4a}}
_mat5b = np.identity(5)
_mat5c = psi4.core.Matrix.from_array(_mat5b)
_mat5a = _mat5b.tolist()
_ans9 = {"C": {"matrix": _mat5a}}
_ans10 = {"B": {"matrix": _mat5a}}
_ans11 = {"C": {"points": _qxyz4a, "diffuse": _qxyzw4a, "matrix": _mat5a}}
_ans12 = {"B": {"points": _qxyz4a, "diffuse": _qxyzw4a}, "A": {"matrix": _mat5a}}

@pytest.mark.parametrize("ep,ans", [
    # lone points, 1-atom
    (_qxyz1a, _ans1),
    (_qxyz1b, _ans1),
    (np.array(_qxyz1a), _ans1),
    ([_qxyz1a], _ans1),
    ([_qxyz1b], _ans1),
    ([np.array(_qxyz1a)], _ans1),
    ({"points": _qxyz1a}, _ans1),
    ({"points": _qxyz1b}, _ans1),
    ({"points": np.array(_qxyz1a)}, _ans1),
    ({"b": _qxyz1a}, _ans2),
    ({"B": [_qxyz1b]}, _ans2),
    ({"b": {"points": _qxyz1a}}, _ans2),
    # lone points, 4-atom
    (_qxyz4a, _ans3),
    (_qxyz4b, _ans3),
    (np.array(_qxyz4a), _ans3),
    ([_qxyz4a], _ans3),
    ([_qxyz4b], _ans3),
    ([np.array(_qxyz4a)], _ans3),
    ({"points": _qxyz4a}, _ans3),
    ({"points": _qxyz4b}, _ans3),
    ({"points": np.array(_qxyz4a)}, _ans3),
    ({"b": _qxyz4a}, _ans4),
    ({"B": [_qxyz4b]}, _ans4),
    ({"b": {"points": _qxyz4a}}, _ans4),
    # lone diffuse
    ([None, _qxyzw4a], _ans5),
    ([None, _qxyzw4b], _ans5),
    ([None, np.array(_qxyzw4a)], _ans5),
    ([None, np.array(_qxyzw4a), None], _ans5),
    ({"diffuse": _qxyzw4a}, _ans5),
    ({"diffuse": _qxyzw4b}, _ans5),
    ({"diffuse": np.array(_qxyzw4a)}, _ans5),
    ({"b": {"diffuse": _qxyzw4a}}, _ans6),
    ({"B": {"diffuse": _qxyzw4b}}, _ans6),
    ({"b": [None, _qxyzw4b]}, _ans6),
    # point and diffuse
    ([_qxyz4a, _qxyzw4a], _ans7),
    ([_qxyz4b, _qxyzw4b], _ans7),
    ([_qxyz4a, np.array(_qxyzw4a), None], _ans7),
    ({"points": _qxyz4a, "diffuse": _qxyzw4a, "matrix": None}, _ans7),
    ({"points": _qxyz4b, "diffuse": _qxyzw4b}, _ans7),
    ({"points": _qxyz4a, "diffuse": np.array(_qxyzw4a)}, _ans7),
    ({"b": {"points": _qxyz4a, "diffuse": _qxyzw4a}}, _ans8),
    ({"B": {"points": _qxyz4b, "diffuse": _qxyzw4b}}, _ans8),
    ({"b": [_qxyz4b, _qxyzw4b]}, _ans8),
    # lone matrix
    ([None, None, _mat5a], _ans9),
    ([None, None, _mat5b], _ans9),
    ([None, None, _mat5c], _ans9),
    ({"matrix": _mat5a}, _ans9),
    ({"b": {"matrix": _mat5b}}, _ans10),
    ({"b": [None, None, _mat5c]}, _ans10),
    # altogether now
    ([_qxyz4a, _qxyzw4a, _mat5a], _ans11),
    ([_qxyz4b, _qxyzw4b, _mat5b], _ans11),
    ([_qxyz4a, np.array(_qxyzw4a), _mat5c], _ans11),
    ({"points": _qxyz4a, "diffuse": _qxyzw4a, "matrix": _mat5c}, _ans11),
    ({"points": _qxyz4b, "diffuse": _qxyzw4b, "matrix": _mat5a}, _ans11),
    ({"points": _qxyz4a, "diffuse": np.array(_qxyzw4a), "matrix": _mat5b}, _ans11),
    ({"b": {"points": _qxyz4a, "diffuse": _qxyzw4a}, "A": [None, None, _mat5b]}, _ans12),
    ({"B": {"points": _qxyz4b, "diffuse": _qxyzw4b}, "A": [None, None, _mat5c]}, _ans12),
    ({"b": [_qxyz4b, _qxyzw4b], "a": {"matrix": _mat5a}}, _ans12),
])
def test_extern_parsing(ep, ans):
    cptd = psi4.p4util.validate_external_potential(ep)
    assert compare_recursive(ans, cptd)


@pytest.mark.parametrize("ep", [
    # alone, unwrapped
    (_qxyz1a[0]),
    (_qxyz1b[0]),
    (_qxyz4a[0]),
    (_qxyz4b[0]),
    (np.array(_qxyz1a[0])),
    (np.array(_qxyz4a[0])),
    # lone diffuse, unlabeled/unlisted
    (_qxyzw4a),
    (_qxyzw4b),
    (np.array(_qxyzw4a)),
    ([_qxyzw4a]),
    ([_qxyzw4b]),
    ([np.array(_qxyzw4a)]),
    ({"b": _qxyzw4a}),
    ({"b": [_qxyzw4b]}),
    # lone matrix unlabeled/unlisted
    (_mat5a),
    (_mat5b),
    (_mat5c),
    ([_mat5a]),
    ([_mat5b]),
    ([_mat5c]),
    ({"b": [_mat5a]}),
    ({"b": [_mat5b]}),
    ({"b": [_mat5c]}),
    ({"C": _mat5a}),
    ({"C": _mat5b}),
    ({"C": _mat5c}),
    # point and diffuse
    ([_qxyzw4a, _qxyz4a]),
    ([_qxyzw4b, _qxyz4b]),
    ([np.array(_qxyzw4a), np.array(_qxyz4a)]),
    ({"points": [_qxyz1a]}),
    ({"points": [_qxyz1b]}),
    ({"points": [np.array(_qxyz1a)]}),
    ({"points": [_qxyz4a]}),
    ({"points": [_qxyz4b]}),
    ({"points": [np.array(_qxyz4a)]}),
    ({"b": {"points": [_qxyz1b]}}),
    ({"b": {"points": [_qxyz4b]}}),
    ({"diffuse": [_qxyzw4a]}),
    ({"diffuse": [_qxyzw4b]}),
    ({"diffuse": [np.array(_qxyzw4a)]}),
    ({"b": {"diffuse": [_qxyzw4b]}}),
    ({"points": _qxyz4a, "diffuse": [_qxyzw4a]}),
    ({"points": _qxyz4b, "diffuse": [_qxyzw4b]}),
    ({"points": _qxyz4a, "diffuse": [np.array(_qxyzw4a)]}),
    ({"b": {"points": _qxyz4b, "diffuse": [_qxyzw4b]}}),
    # bad keys
    ({"pointses": np.array(_qxyz1a)}),
    ({"d": _qxyz1a}),
    # bad dims
    ([None, None, np.zeros((6, 5))]),
    ([None, None, np.zeros((6, 6, 4))]),
    # empty
    # ([None, None, None]),
])
def test_extern_parsing_error(ep):
    with pytest.raises((psi4.ValidationError, TypeError)):
        psi4.p4util.validate_external_potential(ep)
