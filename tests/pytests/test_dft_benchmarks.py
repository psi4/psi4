import numpy as np

import pytest
from utils import *
from addons import using

import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api]


@pytest.fixture
def dft_bench_systems():
    ang = np.array([
      [ -1.551007,  -0.114520,   0.000000],
      [ -1.934259,   0.762503,   0.000000],
      [ -0.599677,   0.040712,   0.000000]])
    oldang = ang * 0.52917721067 / 0.52917720859
    oldmass = [15.99491461956, 1.00782503207, 1.00782503207]

    h2o = psi4.core.Molecule.from_arrays(geom=oldang, elez=[8, 1, 1], units='Angstrom', mass=oldmass)
    h2o_plus = psi4.core.Molecule.from_arrays(geom=oldang, elez=[8, 1, 1], units='Angstrom', mass=oldmass,
                                              molecular_charge=1, molecular_multiplicity=2)

    h2o_dimer = psi4.geometry("""
        0 1
        O  -1.551007  -0.114520   0.000000
        H  -1.934259   0.762503   0.000000
        H  -0.599677   0.040712   0.000000
        --
        0 1
        O   1.350625   0.111469   0.000000
        H   1.680398  -0.373741  -0.758561
        H   1.680398  -0.373741   0.758561
        no_reorient
        """)

    psi4.set_options({
       'dft_radial_points': 200,
       'dft_spherical_points': 590,
       'guess': 'sad',
       'e_convergence': 9,
       'd_convergence': 9,
    })

    return {'h2o': h2o, 'h2o_plus': h2o_plus, 'h2o_dimer': h2o_dimer}


def name_dft_test(val):
    if isinstance(val, float):
        return ' '  # NOT empty string or it'll just default back
    if val in ['6-31G', 'cc-pVDZ']:
        return ' '


@pytest.mark.scf
@pytest.mark.dft
@pytest.mark.parametrize("func,expected,basis", [
    pytest.param(          'B1B95', 0.4535328884, '6-31G'),  # Q-Chem
    pytest.param(          'B1LYP', 0.4497408502, '6-31G'),  # Q-Chem
    pytest.param(         'B1PW91', 0.4538643591, '6-31G'),  # Q-Chem
    pytest.param(          'B3LYP', 0.4572870674, '6-31G'),  # Q-Chem
    pytest.param(         'B3LYP5', 0.4535575461, '6-31G'),  # Q-Chem
    pytest.param(          'B3P86', 0.4638089921, '6-31G'),  # Q-Chem
    pytest.param(         'B3PW91', 0.4568637908, '6-31G'),  # Q-Chem
    pytest.param(         'B3TLAP', 0.4473764358, '6-31G'),  # Q-Chem
    pytest.param(       'B5050LYP', 0.4508657425, '6-31G'),  # Q-Chem
    pytest.param(          'B97-0', 0.4555394941, '6-31G'),  # Q-Chem
    pytest.param(          'B97-1', 0.4551137,    '6-31G'),  # Q-Chem
    pytest.param(          'B97-2', 0.4562077948, '6-31G'),  # Q-Chem
    pytest.param(          'B97-3', 0.4573836831, '6-31G'),  # Q-Chem
    pytest.param(          'B97-D', 0.4562801577, '6-31G'),  # Q-Chem
    pytest.param(         'B97-D3', 0.4562801577, '6-31G', marks=using("dftd3")),  # Q-Chem
    pytest.param(          'B97-K', 0.4498949296, '6-31G'),  # Q-Chem
    pytest.param(         'B97M-V', 0.4561660762, '6-31G'),  # Q-Chem
    pytest.param(           'BB1K', 0.4523318654, '6-31G'),  # Q-Chem
    pytest.param(         'BHHLYP', 0.4474902387, '6-31G'),  # Q-Chem
    pytest.param(           'BLYP', 0.452062348,  '6-31G'),  # Q-Chem
    pytest.param(            'BMK', 0.4586754862, '6-31G'),  # Q-Chem  # not in LibXC 3.0.0
    pytest.param(            'BOP', 0.4518711518, '6-31G'),  # Q-Chem
    pytest.param(           'BP86', 0.4599766012, '6-31G'),  # Q-Chem
    pytest.param(        'BP86VWN', 0.4588991901, '6-31G'),  # Q-Chem
    pytest.param(      'CAM-B3LYP', 0.4568003604, '6-31G'),  # Q-Chem
    pytest.param(       'CAM-LDA0', 0.4688834734, '6-31G'),  # ERKALE
    pytest.param(           'dlDF', 0.4620507089, '6-31G'),  # Q-Chem
    pytest.param('DSD-PBEB95-D3BJ', 0.44615418,   '6-31G', marks=using("dftd3")),  # Q-Chem
    pytest.param('DSD-PBEP86-D3BJ', 0.44599836,   '6-31G', marks=using("dftd3")),  # Q-Chem
    pytest.param('DSD-PBEPBE-D3BJ', 0.44608512,   '6-31G', marks=[*using("dftd3"), pytest.mark.quick]),  # Q-Chem
    pytest.param(           'EDF1', 0.4557270241, '6-31G'),  # Q-Chem
    pytest.param(           'EDF2', 0.457161542,  '6-31G'),  # Q-Chem
    pytest.param(            'GAM', 0.4509109774, '6-31G'),  # Q-Chem
    pytest.param(        'HCTH120', 0.4609912502, '6-31G'),  # Q-Chem
    pytest.param(        'HCTH147', 0.462383627,  '6-31G'),  # Q-Chem
    pytest.param(        'HCTH407', 0.4627707723, '6-31G'),  # Q-Chem
    pytest.param(         'HCTH93', 0.4587611969, '6-31G'),  # Q-Chem
    pytest.param(        'LC-VV10', 0.4556872545, '6-31G'),  # Q-Chem
    pytest.param(           'LDA0', 0.4506129625, '6-31G'),  # ERKALE
    pytest.param(        'LRC-BOP', 0.4595497115, '6-31G'),  # Q-Chem
    pytest.param(       'LRC-wPBE', 0.4580992958, '6-31G'),  # Q-Chem
    pytest.param(      'LRC-wPBEh', 0.4549011451, '6-31G', marks=pytest.mark.quick),  # Q-Chem
    pytest.param(            'M05', 0.457437892,  '6-31G'),  # Q-Chem
    pytest.param(         'M05-2X', 0.4583363493, '6-31G'),  # Q-Chem
    pytest.param(            'M06', 0.4610467104, '6-31G'),  # Q-Chem
    pytest.param(         'M06-2X', 0.4584074697, '6-31G'),  # Q-Chem
    pytest.param(         'M06-HF', 0.4582368218, '6-31G'),  # Q-Chem
    pytest.param(          'M06-L', 0.4568710032, '6-31G'),  # Q-Chem
    pytest.param(         'M08-HX', 0.4616204211, '6-31G'),  # Q-Chem
    pytest.param(         'M08-SO', 0.4656227382, '6-31G'),  # Q-Chem
    pytest.param(            'M11', 0.4596599711, '6-31G'),  # Q-Chem
    pytest.param(          'M11-L', 0.4598708424, '6-31G'),  # Q-Chem
    pytest.param(       'MGGA_MS0', 0.4465224137, '6-31G'),  # Q-Chem
    pytest.param(       'MGGA_MS1', 0.4464339073, '6-31G'),  # Q-Chem
    pytest.param(       'MGGA_MS2', 0.446697433,  '6-31G'),  # Q-Chem
    pytest.param(      'MGGA_MS2h', 0.4467710458, '6-31G'),  # Q-Chem
    pytest.param(       'MGGA_MVS', 0.4646978551, '6-31G'),  # Q-Chem
    pytest.param(      'MGGA_MVSh', 0.4606380903, '6-31G'),  # Q-Chem
    pytest.param(         'MN12-L', 0.4541560093, '6-31G'),  # Q-Chem
    pytest.param(        'MN12-SX', 0.461815166,  '6-31G'),  # Q-Chem
    pytest.param(         'MN15-L', 0.457945977,  '6-31G'),  # Q-Chem
    pytest.param(        'MPW1B95', 0.4537920011, '6-31G'),  # Q-Chem
    pytest.param(          'MPW1K', 0.4527968482, '6-31G'),  # Q-Chem
    pytest.param(        'MPW1LYP', 0.4503072868, '6-31G'),  # Q-Chem
    pytest.param(        'MPW1PBE', 0.4536406557, '6-31G'),  # Q-Chem
    pytest.param(       'MPW1PW91', 0.454429221,  '6-31G'),  # Q-Chem
    pytest.param(         'MPWB1K', 0.4525753564, '6-31G'),  # Q-Chem
    pytest.param(            'N12', 0.4479283792, '6-31G'),  # Q-Chem
    pytest.param(         'N12-SX', 0.4522295631, '6-31G'),  # Q-Chem
    pytest.param(          'O3LYP', 0.4543928605, '6-31G'),  # Q-Chem
    pytest.param(            'PBE', 0.4552076767, '6-31G'),  # Q-Chem
    pytest.param(           'PBE0', 0.4530422829, '6-31G'),  # Q-Chem
    pytest.param(          'PBE50', 0.4509652861, '6-31G'),  # Q-Chem
    pytest.param(          'PBEOP', 0.4515908966, '6-31G'),  # Q-Chem
    pytest.param(         'PBEsol', 0.455663707,  '6-31G'),  # Q-Chem
    pytest.param(           'PKZB', 0.4468952399, '6-31G'),  # Q-Chem
    pytest.param(         'PW6B95', 0.4566580191, '6-31G'),  # Q-Chem
    pytest.param(           'PW91', 0.4573407123, '6-31G'),  # Q-Chem
    pytest.param(          'PWB6K', 0.4535664415, '6-31G', marks=pytest.mark.quick),  # Q-Chem
    pytest.param(         'revPBE', 0.4527402025, '6-31G'),  # Q-Chem
    pytest.param(        'revPBE0', 0.4512137933, '6-31G'),  # Q-Chem
    pytest.param(        'revTPSS', 0.4499706673, '6-31G'),  # Q-Chem
    pytest.param(       'revTPSSh', 0.4497777532, '6-31G'),  # Q-Chem
    pytest.param(           'SCAN', 0.4517181611, '6-31G'),  # Q-Chem
    pytest.param(          'SCAN0', 0.4496859365, '6-31G'),  # Q-Chem
    pytest.param(          'SOGGA', 0.451387499,  '6-31G'),  # Q-Chem
    pytest.param(        'SOGGA11', 0.4529071593, '6-31G'),  # Q-Chem
    pytest.param(      'SOGGA11-X', 0.4601258852, '6-31G'),  # Q-Chem
    pytest.param(         't-HCTH', 0.4623451143, '6-31G'),  # Q-Chem
    pytest.param(        't-HCTHh', 0.4601544314, '6-31G'),  # Q-Chem
    pytest.param(           'TPSS', 0.4510368445, '6-31G'),  # Q-Chem
    pytest.param(          'TPSSh', 0.4505861795, '6-31G'),  # Q-Chem
    pytest.param(           'VSXC', 0.4547894146, '6-31G'),  # Q-Chem
    pytest.param(           'VV10', 0.4551366594, '6-31G'),  # Q-Chem
    pytest.param(           'wB97', 0.4561211941, '6-31G'),  # Q-Chem
    pytest.param(        'wB97M-V', 0.4544676075, '6-31G'),  # Q-Chem
    pytest.param(          'wB97X', 0.4564711283, '6-31G'),  # Q-Chem
    pytest.param(        'wB97X-D', 0.4575912358, '6-31G', marks=pytest.mark.quick),  # Q-Chem
    pytest.param(       'wB97X-D3', 0.4570744381, '6-31G', marks=using("dftd3")),   # Q-Chem
    pytest.param(        'wB97X-V', 0.455302602,  '6-31G'),  # Q-Chem
    pytest.param(         'wM05-D', 0.4560790902, '6-31G'),  # Q-Chem  # https://gitlab.com/libxc/libxc/-/issues/180
    pytest.param(        'wM06-D3', 0.4563459267, '6-31G', marks=using("dftd3")),  # https://gitlab.com/libxc/libxc/-/issues/180
    pytest.param(          'X3LYP', 0.4549534082, '6-31G'),  # Q-Chem
    pytest.param(          'KMLYP', 0.4665060738, '6-31G'),  # ERKALE
], ids=name_dft_test)

def test_dft_bench_ionization(func, expected, basis, dft_bench_systems, request):
    """functionals ionization energies vs. Q-Chem"""
    if func.lower() in psi4.driver.procedures['energy']:
        mols = dft_bench_systems
        psi4.set_options({'basis': basis, 'reference': 'uks'})
        cation  = psi4.energy(func, molecule=mols['h2o_plus'])
        psi4.set_options({'reference': 'rks'})
        neutral = psi4.energy(func, molecule=mols['h2o'])
        assert compare_values(expected, cation - neutral, 4, request.node.name), (cation - neutral - expected)
    else:
        pytest.skip("{0:s} not in Psi4.".format(func))



@pytest.mark.nbody
@pytest.mark.scf
@pytest.mark.dft
@pytest.mark.long
@pytest.mark.parametrize("func,expected,basis", [
    pytest.param(     'B1B95', -0.0132525420,  '6-31G'),  #   Q-Chem
    pytest.param(     'B1LYP', -0.0144226937,  '6-31G'),  #   Q-Chem
    pytest.param(    'B1PW91', -0.0130688293,  '6-31G'),  #   Q-Chem
    pytest.param(     'B3LYP', -0.0145306919,  '6-31G'),  #   Q-Chem
    pytest.param(    'B3LYP5', -0.0144957247,  '6-31G'),  #   Q-Chem
    pytest.param(     'B3P86', -0.0142045385,  '6-31G'),  #   Q-Chem
    pytest.param(    'B3PW91', -0.0133871214,  '6-31G'),  #   Q-Chem
    pytest.param(    'B3TLAP', -0.0155221815,  '6-31G'),  #   Q-Chem
    pytest.param(  'B5050LYP', -0.0149474908,  '6-31G'),  #   Q-Chem
    pytest.param(     'B97-0', -0.0139591073,  '6-31G'),  #   Q-Chem
    pytest.param(     'B97-1', -0.0146773380,  '6-31G'),  #   Q-Chem
    pytest.param(     'B97-2', -0.0132026212,  '6-31G'),  #   Q-Chem
    pytest.param(     'B97-3', -0.0131076955,  '6-31G'),  #   Q-Chem
    pytest.param(     'B97-D', -0.0140422551,  '6-31G'),  #   Q-Chem
    pytest.param(    'B97-D3', -0.0144570179,  '6-31G', marks=using("dftd3")), #   Q-Chem
    pytest.param(     'B97-K', -0.0143346234,  '6-31G'),  #   Q-Chem
    pytest.param( 'B97M-D3BJ', -0.01242543212, 'cc-pVDZ', marks=using("dftd3")), #   Orca
    pytest.param(    'B97M-V', -0.0140939544,  '6-31G'),  #   Q-Chem
    pytest.param(      'BB1K', -0.0134816841,  '6-31G'),  #   Q-Chem
    pytest.param(    'BHHLYP', -0.0148261610,  '6-31G'),  #   Q-Chem
    pytest.param(      'BLYP', -0.0142958198,  '6-31G'),  #   Q-Chem
    pytest.param(       'BMK', -0.0134597648,  '6-31G'),  #   Q-Chem
    pytest.param(       'BOP', -0.0117582062,  '6-31G'),  #   Q-Chem
    pytest.param(      'BP86', -0.0138889593,  '6-31G'),  #   Q-Chem
    pytest.param(   'BP86VWN', -0.0138884252,  '6-31G'),  #   Q-Chem
    pytest.param( 'CAM-B3LYP', -0.0159494700,  '6-31G'),  #   Q-Chem
    pytest.param(  'CAM-LDA0', -0.0172600346,  '6-31G'),  #   ERKALE
    pytest.param(      'dlDF', -0.0100880118,  '6-31G'),  #   Q-Chem
    pytest.param(      'EDF1', -0.0112660306,  '6-31G'),  #   Q-Chem
    pytest.param(      'EDF2', -0.0155177142,  '6-31G'),  #   Q-Chem
    pytest.param(       'GAM', -0.0149107345,  '6-31G'),  #   Q-Chem
    pytest.param(   'HCTH120', -0.0140214969,  '6-31G'),  #   Q-Chem
    pytest.param(   'HCTH147', -0.0132522028,  '6-31G'),  #   Q-Chem
    pytest.param(   'HCTH407', -0.0135769693,  '6-31G'),  #   Q-Chem
    pytest.param(    'HCTH93', -0.0110370744,  '6-31G'),  #   Q-Chem
    pytest.param(   'LC-VV10', -0.0156463704,  '6-31G'),  #   Q-Chem
    pytest.param(      'LDA0', -0.0178827299,  '6-31G'),  #   ERKALE
    pytest.param(   'LRC-BOP', -0.0157298937,  '6-31G'),  #   Q-Chem
    pytest.param(  'LRC-wPBE', -0.0147152576,  '6-31G'),  #   Q-Chem
    pytest.param( 'LRC-wPBEh', -0.0147072182,  '6-31G'),  #   Q-Chem
    pytest.param(       'M05', -0.0156581803,  '6-31G'),  #   Q-Chem
    pytest.param(    'M05-2X', -0.0151408376,  '6-31G'),  #   Q-Chem
    pytest.param(       'M06', -0.0143374034,  '6-31G'),  #   Q-Chem
    pytest.param(    'M06-2X', -0.0149081161,  '6-31G'),  #   Q-Chem
    pytest.param(    'M06-HF', -0.0153050734,  '6-31G'),  #   Q-Chem
    pytest.param(     'M06-L', -0.0131734851,  '6-31G'),  #   Q-Chem
    pytest.param(    'M08-HX', -0.0150796776,  '6-31G'),  #   Q-Chem
    pytest.param(    'M08-SO', -0.0154578319,  '6-31G'),  #   Q-Chem
    pytest.param(       'M11', -0.0152106441,  '6-31G'),  #   Q-Chem
    pytest.param(     'M11-L', -0.0120196336,  '6-31G'),  #   Q-Chem
    pytest.param(  'MGGA_MS0', -0.0141108035,  '6-31G'),  #   Q-Chem
    pytest.param(  'MGGA_MS1', -0.0135190604,  '6-31G'),  #   Q-Chem
    pytest.param(  'MGGA_MS2', -0.0138641589,  '6-31G'),  #   Q-Chem
    pytest.param( 'MGGA_MS2h', -0.0138846318,  '6-31G'),  #   Q-Chem
    pytest.param(  'MGGA_MVS', -0.0145231638,  '6-31G'),  #   Q-Chem
    pytest.param( 'MGGA_MVSh', -0.0145138947,  '6-31G'),  #   Q-Chem
    pytest.param(    'MN12-L', -0.0122476159,  '6-31G'),  #   Q-Chem
    pytest.param(   'MN12-SX', -0.0130140323,  '6-31G'),  #   Q-Chem
    pytest.param(    'MN15-L', -0.0122572417,  '6-31G'),  #   Q-Chem
    pytest.param(   'MPW1B95', -0.0143828237,  '6-31G'),  #   Q-Chem
    pytest.param(     'MPW1K', -0.0142603959,  '6-31G'),  #   Q-Chem
    pytest.param(   'MPW1LYP', -0.0156088011,  '6-31G'),  #   Q-Chem
    pytest.param(   'MPW1PBE', -0.0142137779,  '6-31G'),  #   Q-Chem
    pytest.param(  'MPW1PW91', -0.0142534651,  '6-31G'),  #   Q-Chem
    pytest.param(    'MPWB1K', -0.0143991465,  '6-31G'),  #   Q-Chem
    pytest.param(       'N12', -0.0153151670,  '6-31G'),  #   Q-Chem
    pytest.param(    'N12-SX', -0.015641602,   '6-31G'),  #   Q-Chem
    pytest.param(     'O3LYP', -0.0124587623,  '6-31G'),  #   Q-Chem
    pytest.param(       'PBE', -0.0154580371,  '6-31G'),  #   Q-Chem
    pytest.param(      'PBE0', -0.0149591069,  '6-31G'),  #   Q-Chem
    pytest.param(     'PBE50', -0.01474869,    '6-31G'),  #   Q-Chem
    pytest.param(     'PBEOP', -0.0142184544,  '6-31G'),  #   Q-Chem
    pytest.param(    'PBEsol', -0.0169733502,  '6-31G'),  #   Q-Chem
    pytest.param(      'PKZB', -0.0107309601,  '6-31G'),  #   Q-Chem
    pytest.param(    'PW6B95', -0.0145444936,  '6-31G'),  #   Q-Chem
    pytest.param(      'PW91', -0.0160147465,  '6-31G'),  #   Q-Chem
    pytest.param(     'PWB6K', -0.0150275431,  '6-31G'),  #   Q-Chem
    pytest.param(    'revPBE', -0.0129363604,  '6-31G'),  #   Q-Chem
    pytest.param(   'revPBE0', -0.0131251963,  '6-31G'),  #   Q-Chem
    pytest.param(   'revTPSS', -0.0140137038,  '6-31G'),  #   Q-Chem
    pytest.param(  'revTPSSh', -0.0139563411,  '6-31G'),  #   Q-Chem
    pytest.param(      'SCAN', -0.0151696769,  '6-31G'),  #   Q-Chem
    pytest.param(     'SCAN0', -0.0151046408,  '6-31G'),  #   Q-Chem
    pytest.param(     'SOGGA', -0.0171800268,  '6-31G'),  #   Q-Chem
    pytest.param(   'SOGGA11', -0.0133998476,  '6-31G'),  #   Q-Chem
    pytest.param( 'SOGGA11-X', -0.0141398574,  '6-31G'),  #   Q-Chem
    pytest.param(    't-HCTH', -0.0139348053,  '6-31G'),  #   Q-Chem
    pytest.param(   't-HCTHh', -0.0143195773,  '6-31G'),  #   Q-Chem
    pytest.param(      'TPSS', -0.0141667663,  '6-31G'),  #   Q-Chem
    pytest.param(     'TPSSh', -0.0140808774,  '6-31G'),  #   Q-Chem
    pytest.param(      'VSXC', -0.0138314231,  '6-31G'),  #   Q-Chem
    pytest.param(      'VV10', -0.0160812898,  '6-31G'),  #   Q-Chem
    pytest.param(      'wB97', -0.0160148410,  '6-31G'),  #   Q-Chem
    pytest.param('wB97M-D3BJ', -0.01351626999, 'cc-pVDZ', marks=using("dftd3")), #   Orca
    pytest.param(   'wB97M-V', -0.0152233456,  '6-31G'),  #   Q-Chem
    pytest.param(     'wB97X', -0.0156289024,  '6-31G'),  #   Q-Chem
    pytest.param(   'wB97X-D', -0.0146246032,  '6-31G'),  #   Q-Chem
    pytest.param(  'wB97X-D3', -0.0148666307,  '6-31G', marks=using("dftd3")), #  Q-Chem
    pytest.param('wB97X-D3BJ', -0.01295664452, 'cc-pVDZ', marks=using("dftd3")), #   Orca
    pytest.param(   'wB97X-V', -0.0151102751,  '6-31G'),  #   Q-Chem
    pytest.param(    'wM05-D', -0.0147496512,  '6-31G'),  #   Q-Chem  # https://gitlab.com/libxc/libxc/-/issues/180
    pytest.param(   'wM06-D3', -0.0151611219,  '6-31G'),  #   Q-Chem  # https://gitlab.com/libxc/libxc/-/issues/180
    pytest.param(     'X3LYP', -0.0151870467,  '6-31G'),  #   Q-Chem
    pytest.param(     'KMLYP', -0.0173175768,  '6-31G'),  #   ERKALE
], ids=name_dft_test)

def test_dft_bench_interaction(func, expected, basis, dft_bench_systems, request):
    """functionals interaction energies vs. Q-Chem & Orca"""
    if func.lower() in psi4.driver.procedures['energy']:
        mols = dft_bench_systems
        psi4.set_options({'basis': basis})
        psi4_ie = psi4.energy(func, molecule=mols['h2o_dimer'], bsse_type='nocp')
        assert compare_values(expected, psi4_ie, 4, request.node.name), (psi4_ie - expected)
    else:
        pytest.skip("{0:s} not in Psi4.".format(func))

# Current version of Psi4 does not match Q-Chem for these tests
#expected_fail_qchem = ['B97-D', 'wB97X-D3'] #TEST


# ionization energy references from an older version of Psi4 (~April 2017, SHA: 53e752c)
old_psi4_data = {
     'BP86': 0.45997666061779796,
     'BLYP': 0.4520631769723451,
   'B97-D3': 0.4562813415641358,      # @using_dftd3
   'B3LYP5': 0.4535570976370309,
      'PBE': 0.45520846640530976,
   'M05-2X': 0.4583328495011756,
  'wB97X-D': 0.4576016432866794,
     'wB97': 0.45615296910915504,
  'HCTH120': 0.46099230571738303,
     'dlDF': 0.46204786302111245,
     'PBE0': 0.4530414801763101,
    'B97-1': 0.4551134986781449,
     'PW91': 0.45735016732695044,
    'B97-D': 0.4562813415641216,
      'M05': 0.4574377079341474,
    'B3LYP': 0.4572866176143293,
  'HCTH407': 0.46277176698060885,
  'HCTH147': 0.46238465771986625,
    'SOGGA': 0.45733831961008775,     # this pre-LibXC value known to not match current Psi4
    'B97-2': 0.4562076272818558,
    'wB97X': 0.4588992888233463,      # this pre-LibXC value known to not match current Psi4
}


# interaction energy references from an older version of Psi4 (~April 2017, SHA: 53e752c)
old_psi4_data = {
      'BP86': -0.013889878815405154,
      'BLYP': -0.014294874492634335,
    'B97-D3': -0.014455902448531788,  # @using_dftd3
    'B3LYP5': -0.014494768950584103,
       'PBE': -0.015457134998626998,
    'M05-2X': -0.015145056492713138,
   'wB97X-D': -0.014613361929264101,
      'wB97': -0.015991295341081013,
   'HCTH120': -0.014020533321001949,
      'dlDF': -0.010087591449917,
      'PBE0': -0.014958205517643819,
     'B97-1': -0.014676434498966273,
      'PW91': -0.01602902561316455,
       'M05': -0.015654414236166758,
     'B3LYP': -0.014529738435527406,
   'HCTH407': -0.013575929357500627,
   'HCTH147': -0.01325130624383064,
     'SOGGA': -0.01806810274382542,
     'B97-2': -0.013201680603032173,
     'B97-D': -0.014041140231540794,  # this pre-LibXC value known to not match current Psi4
     'wB97X': -0.01638428793378921,   # this pre-LibXC value known to not match current Psi4
}
