import pytest
from utils import *
from addons import *
import numpy as np

import psi4
from psi4.driver.procrouting.response.scf_response import tdscf_excitations
from psi4.driver.procrouting.response.scf_products import TDRSCFEngine, TDUSCFEngine

@pytest.fixture
def tddft_systems():
    psi4.core.clean()
    ch2 = psi4.geometry("""
    0 3
    C        0.000000    0.000000    0.159693
    H       -0.000000    0.895527   -0.479080
    H       -0.000000   -0.895527   -0.479080
    symmetry c1
    no_reorient
    no_com
    """)
    h2o = psi4.geometry("""0 1
    O          0.000000    0.000000    0.135446
    H         -0.000000    0.866812   -0.541782
    H         -0.000000   -0.866812   -0.541782
    symmetry c1
    no_reorient
    no_com
    """)
    psi4.set_options({
        'scf_type': 'pk',
        'e_convergence': 8,
        'd_convergence': 8,
        'save_jk': True
        })

    return {'UHF': ch2, 'RHF': h2o}

def name_tdscf_run(reference, func, tda, triplets, basis, expected):
    if reference == 'RHF':
        if triplets:
            name = 'RHF [Triplet]'
        else:
            name = 'RHF [Singlet]'
    else:
        name     = 'UHF          '
    name += "/{}".format(func)
    if tda:
        name += "/TDA"
    else:
        name += "/RPA"
    name += "/{}".format(basis)
    return name

@pytest.mark.tdscf
@pytest.mark.parametrize("reference,func,tda,triplets,basis,expected", [
    pytest.param(*args, id=name_tdscf_run(*args)) for args in [
    # ('RHF','HF', False, False,'cc-pvdz',[[0.37554133], [0.33365445], [0.27878771], [0.42481114]]), # G09 rev E.01
    # ('RHF','HF',  True, False,'cc-pvdz',[[0.37988510], [0.33726691], [0.28225317], [0.43006997]]), # G09 rev E.01
    # ('UHF','HF', False, False,'cc-pvdz',[[], [0.28785744, 0.35473019], [0.31791104], [0.24457042]]), # G09 rev E.01
    # ('UHF','HF',  True, False,'cc-pvdz',[[], [0.28968756, 0.35741289], [0.32054964], [0.24880026]]), # G09 rev E.01
    # ('RHF','HF', False,  True,'cc-pvdz',[[0.26804588], [0.30135179], [0.23573588], [0.30942930]]), # G09 rev E.01
    # ('RHF','HF',  True,  True,'cc-pvdz',[[0.29549763], [0.30883763], [0.24288362], [0.33788430]]), # G09 rev E.01
    # ('RHF',    'HF', False, False,'sto-3g',[[0.50011418], [0.41531879], [0.35477796], [0.55140153]]), # G09 rev E.01
    # ('RHF',    'HF',  True, False,'sto-3g',[[0.50564128], [0.41607215], [0.35646108], [0.55521611]]), # G09 rev E.01
    # ('UHF',    'HF', False, False,'sto-3g',[[0.53924486], [0.39001321, 0.56554635], [], [0.32329849]]), # G09 rev E.01
    # ('UHF',    'HF',  True, False,'sto-3g',[[0.55163673], [0.39029985, 0.56582932], [], [0.32673823]]), # G09 rev E.01
    # ('RHF',    'HF', False,  True,'sto-3g',[[0.29976422], [0.36513024], [0.28516372], [0.35263180]]), # G09 rev E.01
    # ('RHF',    'HF',  True,  True,'sto-3g',[[0.34444405], [0.36599018], [0.28725475], [0.39451500]]), # G09 rev E.01
    ('RHF','HCTH93', False, False,'cc-pvdz',[[0.32890644], [0.28941929], [0.22782742], [0.38849179]]), # G09 rev E.01
    ('RHF','HCTH93',  True, False,'cc-pvdz',[[0.33356625], [0.28974636], [0.22904383], [0.39102014]]), # G09 rev E.01
    # ('UHF','HCTH93', False, False,'cc-pvdz',[[], [0.24517678, 0.30895155], [0.26726680], [0.18563553]]), # G09 rev E.01
    # ('UHF','HCTH93',  True, False,'cc-pvdz',[[], [0.24640421, 0.30919042], [0.26811571], [0.18766409]]), # G09 rev E.01
    # ('RHF','HCTH93', False,  True,'cc-pvdz',[[0.27850107], [0.26766737], [0.20068070], [0.34208107]]), # G09 rev E.01
    # ('RHF','HCTH93',  True,  True,'cc-pvdz',[[0.28172766], [0.26861182], [0.20200735], [0.34495487]]), # G09 rev E.01
        ]])
def test_tdscf(reference, func, tda, triplets, expected, basis, tddft_systems):
    mol = tddft_systems[reference]
    psi4.set_options({'reference': reference, 'basis': basis})



    # gather tddft driver args
    tdscf_kwargs = {
            'tda': tda,
            'states_per_irrep': [4],
            'r_tol': 1.0e-5,
            'e_tol': 1.0e-4,
            'triplets': triplets
            }
    label_fmt = "{ref}-{func}-{count}-{sym}"
    if triplets:
        label_fmt+= "-Triplet"
    e, wfn = psi4.energy(func, molecule=mol, return_wfn=True)
    e_out = tdscf_excitations(wfn, **tdscf_kwargs)
    energies = []
    for ee in e_out:
        energies.extend(ee)
    exp = []
    for ee in expected:
        exp.extend(ee)
    expected = np.array(exp)
    expected = expected[np.argsort(expected)]
    energies = np.array(energies)
    energies = energies[np.argsort(energies)]
    for i, (ref, me) in enumerate(zip(expected, energies[:len(expected)])):
        assert compare_values(ref, me, 3, label_fmt.format(ref=reference, func=func, count = str(i+1), sym="??"))



