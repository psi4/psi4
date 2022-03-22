import pytest
import psi4
import os
from shutil import copytree
from psi4.driver.p4util.testing import compare_strings, compare_values, compare_integers
from psi4.driver.p4util.exceptions import ValidationError

pytestmark = [pytest.mark.psi, pytest.mark.api]

# Checks for
# Molden files are the same as reference

@pytest.fixture
def datadir(tmpdir, request):
    """
    from: https://stackoverflow.com/a/29631801
    Fixture responsible for searching a folder with the same name of test
    module and, if available, moving all contents to a temporary directory so
    tests can use them freely.
    """
    filename = request.module.__file__
    test_dir, _ = os.path.splitext(filename)

    if os.path.isdir(test_dir):
        copytree(test_dir, str(tmpdir), dirs_exist_ok=True)

    return tmpdir


@pytest.mark.parametrize('inp_h2o', [
    pytest.param({'name': 'h2o_normal', 'energy': 'scf', 'do_virtual':True, 'use_natural': False, 'options': {'e_convergence': 10}}, id='h2o_normal'),
    pytest.param({'name': 'dovirt_false', 'energy': 'scf', 'do_virtual':False, 'use_natural': False, 'options': {'e_convergence': 10}}, id='dovirt_false'),
    pytest.param({'name': 'orbso_detci', 'energy': 'cisd', 'do_virtual':True, 'use_natural': True, 'options': {'e_convergence': 10, 'qc_module':'detci', 'opdm':True}}, id='orbso_detci')
    ])
def test_H2O_molden(inp_h2o, datadir):
    mol = psi4.geometry("""
            0 1
            O
            H   1   0.951342
            H   1   0.951342   2   112.505645
            """)
    psi4.set_options({
        'basis': 'dz',
        'scf_type': 'pk',
        })
    psi4.set_options(inp_h2o['options'])
    molden_file = f"{inp_h2o['name']}.molden"
    ref = datadir.join(f"{inp_h2o['name']}.ref")
    e, wfn = psi4.energy(inp_h2o['energy'], return_wfn=True, molecule=mol)
    wfn.write_molden(molden_file, do_virtual=inp_h2o['do_virtual'], use_natural=inp_h2o['use_natural'])
    assert psi4.compare_moldenfiles(ref, molden_file)

@pytest.mark.parametrize('inp_h2o_density', [
    pytest.param({'name': 'orbso_density', 'energy': 'ccsd', 'do_virtual':True, 'use_natural': True}, id='orbso_density'),
    ])
def test_H2O_density_molden(inp_h2o_density, datadir):
    mol = psi4.geometry("""
            0 1
            O
            H   1   0.951342
            H   1   0.951342   2   112.505645
            """)
    psi4.set_options({
        'basis': 'dz',
        'scf_type': 'pk',
        'e_convergence': 10
        })
    molden_file = f"{inp_h2o_density['name']}.molden"
    ref = datadir.join(f"{inp_h2o_density['name']}.ref")
    e, wfn = psi4.properties(inp_h2o_density['energy'], return_wfn=True, molecule=mol)
    wfn.write_molden(molden_file, do_virtual=inp_h2o_density['do_virtual'], use_natural=inp_h2o_density['use_natural'])
    assert psi4.compare_moldenfiles(ref, molden_file)

@pytest.mark.parametrize('inp_oh', [
    pytest.param({'name': 'ref_uhf', 'ref':'uhf'}, id='ref_uhf'),
    pytest.param({'name': 'ref_rohf', 'ref':'rohf'}, id='ref_rohf')
    ])
def test_OH_molden(inp_oh, datadir):
    mol = psi4.geometry("""
            0 2
            O
            H   1   0.970369
            symmetry c1
            """)
    psi4.set_options({
        'basis': 'dz',
        'scf_type': 'pk',
        'e_convergence': 11,
        'reference':inp_oh['ref']
        })
    molden_file = f"{inp_oh['name']}.molden"
    ref = datadir.join(f"{inp_oh['name']}.ref")
    e, wfn = psi4.energy('scf', return_wfn=True, molecule=mol)
    wfn.write_molden(molden_file, do_virtual=True, use_natural=False)
    assert psi4.compare_moldenfiles(ref, molden_file)

@pytest.mark.parametrize('inp_h2s', [
    pytest.param({'name': 'dorbs_cartesian', 'options': {'basis': '6-31g'}}, id='dorbs_cartesian'),
    pytest.param({'name': 'dorbs_spherical', 'options': {'basis': 'dz'}}, id='dorbs_spherical'),
    ])
def test_H2S_molden(inp_h2s, datadir):
    mol = psi4.geometry("""
            0 1
            S
            H   1   1.350490
            H   1   1.350490   2   96.061977
            """)
    psi4.set_options({
        'scf_type': 'pk',
        'e_convergence': 10
        })
    psi4.set_options(inp_h2s['options'])
    molden_file = f"{inp_h2s['name']}.molden"
    ref = datadir.join(f"{inp_h2s['name']}.ref")
    e, wfn = psi4.energy('scf', return_wfn=True, molecule=mol)
    wfn.write_molden(molden_file, do_virtual=True, use_natural=False)
    assert psi4.compare_moldenfiles(ref, molden_file)

@pytest.mark.parametrize('inp_clfhcoh', [
    pytest.param({'name': 'sym_trivial'}, id='sym_trivial'),
    ])
def test_ClFHCOH_molden(inp_clfhcoh, datadir):
    mol = psi4.geometry("""
            0 1
            C
            F   1   1.395520
            Cl  1   1.853978   2   106.297922
            H   1   1.066516   2   109.322008   3  -116.352650
            O   1   1.363622   2   110.838591   3   122.400775
            H   5   0.955096   1   116.200547   3    59.282816
            """)
    psi4.set_options({
        'basis': 'dz',
        'scf_type': 'pk',
        'e_convergence': 11
        })
    molden_file = f"{inp_clfhcoh['name']}.molden"
    ref = datadir.join(f"{inp_clfhcoh['name']}.ref")
    e, wfn = psi4.energy('scf', return_wfn=True, molecule=mol)
    wfn.write_molden(molden_file, do_virtual=True, use_natural=False)
    assert psi4.compare_moldenfiles(ref, molden_file)
