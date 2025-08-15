#! Multilevel computation of water trimer energy (geometry from J. Chem. Theory Comput. 11, 2126-2136 (2015))
import copy
import pprint

import pytest

import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.nbody]

@pytest.fixture
def base_schema():
    h2o_trimer = psi4.geometry("""
    O      -2.76373224  -1.24377706  -0.15444566
    H      -1.12357791  -2.06227970  -0.05243799
    H      -3.80792362  -2.08705525   1.06090407
    --
    O       2.46924614  -1.75437739  -0.17092884
    H       3.76368260  -2.21425403   1.00846104
    H       2.30598330   0.07098445  -0.03942473
    --
    O       0.29127930   3.00875625   0.20308515
    H      -1.21253048   1.95820900   0.10303324
    H       0.10002049   4.24958115  -1.10222079
    units bohr
    """)

    _base_schema = {
        'schema_name': 'qcschema_input',
        'schema_version': 1,
        'molecule': h2o_trimer.to_schema(dtype=2),
        'keywords': {
            'function_kwargs': {
            },
            'cc_type': 'df',
        },
        'model': {
            'basis': '(auto)',
        },
        'driver': 'energy',
    }

    return _base_schema


@pytest.mark.parametrize("inp,expected", [
    # Compute 1-body contribution with ccsd(t) and 2-body contribution with mp2
    pytest.param({'method': '', 'kfk': {'bsse_type': ['nocp', 'cp', 'vmfc'], 'return_total_data': True, 'levels': {1: 'mp2/sto-3g', 2: 'scf/sto-3g'}}},
                 #{'2NOCP': -225.019408434635, '2CP': -225.000173661598, '2VMFC': -224.998744381484},
                 {'NOCP-CORRECTED TOTAL ENERGY THROUGH 2-BODY': -225.019408434635,
                  'CP-CORRECTED TOTAL ENERGY THROUGH 2-BODY': -225.000173661598,
                  'VMFC-CORRECTED TOTAL ENERGY THROUGH 2-BODY': -224.998744381484},
                 id='nbody-multilevel'),
    # Compute 1-body contribution with ccsd(t) and estimate all higher order contributions with scf
    pytest.param({'method': '', 'kfk': {'bsse_type': 'nocp', 'return_total_data': True, 'levels': {1: 'mp2/sto-3g', 'supersystem': 'scf/sto-3g'}}},
                 #{'1NOCP': -224.998373505116, '3NOCP': -225.023509855159},
                 {'NOCP-CORRECTED TOTAL ENERGY THROUGH 1-BODY': -224.998373505116,
                  'NOCP-CORRECTED TOTAL ENERGY THROUGH 3-BODY': -225.023509855159},
                 id='nbody-multilevel-supersys'),
    # Compute electrostatically embedded  many-body expansion energy.with TIP3P charges
    pytest.param({'method': 'scf/sto-3g', 'kfk': {'bsse_type': 'vmfc', 'return_total_data': True, 'levels': None, 'max_nbody': 2,
                                                  'embedding_charges': {i: [j for j in [-0.834, 0.417, 0.417]] for i in range(1, 4)}}},
                 #{'1': -224.940138148882, '2': -224.943882712817},
                 {'VMFC-CORRECTED TOTAL ENERGY THROUGH 1-BODY': -224.940138148882,
                  'VMFC-CORRECTED TOTAL ENERGY THROUGH 2-BODY': -224.943882712817},
                 id='nbody-embedded', marks=pytest.mark.extern),
])
def test_nbody_levels(inp, expected, base_schema, monkeypatch):
    monkeypatch.setenv("QCMANYBODY_EMBEDDING_CHARGES", "1")
    # reference for nbody-multilevel generated with this larger fitting basis for sto-3g. fails otherwise by 3.e-5
    basfams = psi4.driver.qcdb.basislist.load_basis_families()
    for fam in basfams:
        if fam.ornate == "STO-3G":
            fam.add_rifit("def2-qzvpp-ri")

    jin = copy.deepcopy(base_schema)
    jin['model']['method'] = inp['method']
    jin['keywords']['function_kwargs'] = inp['kfk']

    if psi4.core.get_option("scf", "orbital_optimizer_package") != "INTERNAL":
        jin["keywords"].update({"e_convergence": 9, "d_convergence": 5e-8})

    otp = psi4.schema_wrapper.run_qcschema(jin)
    pprint.pprint(otp)

    for b, v in expected.items():
        assert psi4.compare_values(v, otp.extras["qcvars"][b], 6, b)

