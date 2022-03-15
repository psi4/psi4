import pytest
import psi4
import os
from distutils import dir_util
import re
from psi4.driver.p4util.testing import compare_strings, compare_values, compare_integers
from psi4.driver.p4util.exceptions import ValidationError


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
        dir_util.copy_tree(test_dir, str(tmpdir))

    return tmpdir


@pytest.mark.parametrize('inp_h2o', [
    pytest.param({'name': 'h2o_normal', 'energy': 'scf', 'do_virtual':True, 'use_natural': False, 'options': {'r_convergence':12}}, id='h2o_normal'),
    pytest.param({'name': 'dovirt_false', 'energy': 'scf', 'do_virtual':False, 'use_natural': False, 'options': {'r_convergence':12}}, id='dovirt_false'),
    pytest.param({'name': 'orbso_detci', 'energy': 'cisd', 'do_virtual':True, 'use_natural': True, 'options': {'r_convergence':11, 'qc_module':'detci', 'opdm':True}}, id='orbso_detci')
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
    assert compare_moldenfiles(ref, molden_file)

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
        'r_convergence': 12
        })
    molden_file = f"{inp_h2o_density['name']}.molden"
    ref = datadir.join(f"{inp_h2o_density['name']}.ref")
    e, wfn = psi4.properties(inp_h2o_density['energy'], return_wfn=True, molecule=mol)
    wfn.write_molden(molden_file, do_virtual=inp_h2o_density['do_virtual'], use_natural=inp_h2o_density['use_natural'])
    assert compare_moldenfiles(ref, molden_file)

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
        'd_convergence': 8,
        'reference':inp_oh['ref']
        })
    molden_file = f"{inp_oh['name']}.molden"
    ref = datadir.join(f"{inp_oh['name']}.ref")
    e, wfn = psi4.energy('scf', return_wfn=True, molecule=mol)
    wfn.write_molden(molden_file, do_virtual=True, use_natural=False)
    assert compare_moldenfiles(ref, molden_file)

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
        'r_convergence': 12
        })
    psi4.set_options(inp_h2s['options'])
    molden_file = f"{inp_h2s['name']}.molden"
    ref = datadir.join(f"{inp_h2s['name']}.ref")
    e, wfn = psi4.energy('scf', return_wfn=True, molecule=mol)
    wfn.write_molden(molden_file, do_virtual=True, use_natural=False)
    assert compare_moldenfiles(ref, molden_file)

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
        'r_convergence': 12
        })
    molden_file = f"{inp_clfhcoh['name']}.molden"
    ref = datadir.join(f"{inp_clfhcoh['name']}.ref")
    e, wfn = psi4.energy('scf', return_wfn=True, molecule=mol)
    wfn.write_molden(molden_file, do_virtual=True, use_natural=False)
    assert compare_moldenfiles(ref, molden_file)

def moldenfile_to_string(fname):
    with open(fname, 'r') as fn:
        molden_string = fn.read()
    return molden_string
    
def compare_moldenfiles(expected, computed, digits=7, label='Compare Molden'):
    ref = moldenfile_to_string(expected).splitlines()
    calc = moldenfile_to_string(computed).splitlines()
    if len(ref) != len(calc):
        raise ValidationError(f"These two molden files have different lengths...\n")
    
    high_accuracy = digits
    index = 0
    max_len = len(calc)
    tests = []
    section = 0

    geom_re = re.compile(r'^\s*(\w*)\s+(\d+)\s+(\d+)\s+(-?\d+.\d+)\s+(-?\d+.\d+)\s+(-?\d+.\d+)\s*$')
    basis_header_re = re.compile(r'^\s*([s,p,d,f,g])\s*(\d*)\s*(\d*.\d*)\s*$')
    s1_re = re.compile(r'^\s*(\d+.?\d*)\s+(\d+.?\d*)$')
    s2_re = re.compile(r'^\s*(\d+)\s+(-?\d+.\d+[e,E][\+,-]\d+)\s*$')
    sym_re = re.compile(r'^\s*Sym\s*=\s*(\w*)\s*$')
    energy_re = re.compile(r'^\s*Ene\s*=\s*(-?\d*.?\d*[e,E]?\+?-?\d*)\s*$')
    spin_re = re.compile(r'^\s*Spin\s*=\s*(\w*)\s*$')
    occ_re = re.compile(r'^\s*Occup\s*=\s*(-?\d*.\d*[e,E]?-?\+?\d*)\s*$')

    for i in range(max_len):
        line = calc[i]
        
        if geom_re.match(line):
            c1, c2, c3, c4, c5, c6 = geom_re.match(line).groups()
            r1, r2, r3, r4, r5, r6 = geom_re.match(line).groups()
            test = compare_strings(r1, c1) and compare_integers(r2, c2) and compare_integers(r3, c3) and compare_values(r4, c4, high_accuracy) and compare_values(r5, c5, high_accuracy) and compare_values(r6, c6, high_accuracy)

        elif basis_header_re.match(line):
            c1, c2, c3 = basis_header_re.match(line).groups()
            r1, r2, r3 = basis_header_re.match(ref[i]).groups()
            test = compare_strings(r1,c1) and compare_integers(r2,c2) and compare_values(r3,c3,3)

        elif s1_re.match(line):
            c1, c2 = s1_re.match(line).groups()
            r1, r2 = s1_re.match(ref[i]).groups()
            test = compare_values(r1, c1, high_accuracy) and compare_values(r2, c2, high_accuracy)

        elif sym_re.match(line):
            c = sym_re.match(line).group(1)
            r = sym_re.match(ref[i]).group(1)
            test = compare_strings(r, c, f'text line: {line}')
        
        elif energy_re.match(line):
            c = energy_re.match(line).group(1)
            r = energy_re.match(ref[i]).group(1)
            test = compare_values(r, c, high_accuracy, f'float value: {line}')

        elif spin_re.match(line):
            c = spin_re.match(line).group(1)
            r = spin_re.match(ref[i]).group(1)
            test = compare_strings(r, c, f'text line: {line}')
        
        elif occ_re.match(line):
            c = occ_re.match(line).group(1)
            r = occ_re.match(ref[i]).group(1)
            test = compare_values(r, c, high_accuracy, f'float value: {line}')
        
        elif s2_re.match(line):
            c1, c2 = s2_re.match(line).groups()
            r1, r2 = s2_re.match(line).groups()
            test = compare_integers(r1, c1, f'int value: {line}') and compare_values(r2, c2, high_accuracy, f'float value: {line}')

        else:
            test = compare_strings(line, ref[i])

        tests.append(test)

    return compare_integers(True, all(tests), label)


