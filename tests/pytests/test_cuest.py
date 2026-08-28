import subprocess
import sys

import psi4
import pytest
import numpy as np

from psi4 import core

from utils import *
from addons import using, uusing

pytestmark = [pytest.mark.psi, pytest.mark.api]

do_fd = False

__geoms = {
    "methylamine":
    """
    units bohr
    N   -1.2443656662    1.5116296128    1.1401094834
    C    1.2072325506    0.2976885350    0.9288948097
    H    2.1975245576    0.3690634856    2.7411987245
    H    2.3661826021    1.2428916800   -0.4976593391
    H    0.9782544525   -1.6841238010    0.3901528276
    H   -2.1240520112    1.4799955998   -0.5728893307
    H   -1.0022350753    3.3701184308    1.5836092757
    """,
    "methylamine_cation":
    """
    units bohr
    1 2
    N   -1.2443656662    1.5116296128    1.1401094834
    C    1.2072325506    0.2976885350    0.9288948097
    H    2.1975245576    0.3690634856    2.7411987245
    H    2.3661826021    1.2428916800   -0.4976593391
    H    0.9782544525   -1.6841238010    0.3901528276
    H   -2.1240520112    1.4799955998   -0.5728893307
    H   -1.0022350753    3.3701184308    1.5836092757
    """,
    }

# Reference values below are taken from regular Psi4 unless otherwise stated.
# HF gradients are from analytic calculations, while DFT gradients are from
# 5 point finite differences due to that fact that the grid weight derivatives
# are missing from the regular Psi4 analytic gradients.

@pytest.mark.smoke
@pytest.mark.quick
@uusing("cuest")
@uusing("cuda_cc8")
@pytest.mark.parametrize("inp", [

    # This does not work - the aux basis must be spherical for cuEST.  The
    # orbital basis can be Cartesian, however.
    #pytest.param({
    #    "geom": __geoms["methylamine"],
    #    "methodname": "hf",
    #    "args": {
    #        "basis": "def2-svp",
    #        "puream": False,
    #        },
    #    "ref_energy":  -95.14919233891949 ,
    #    "ref_gradient":  [[ 0.0048455767, -0.009222819 ,  0.0055741773],
    #                      [ 0.0047066585, -0.0063561858,  0.0031168156],
    #                      [ 0.0007370865,  0.0009772166,  0.0015993905],
    #                      [-0.0084520345,  0.0044652126,  0.0004904717],
    #                      [-0.0007077581, -0.0014647869, -0.0012034852],
    #                      [-0.0045594848, -0.000927345 , -0.0124510086],
    #                      [ 0.0034299557,  0.0125287075,  0.0028736387]]
    #}, id='methylamine_rhf_vacuum_cartesian'),]

    pytest.param({
        "geom": __geoms["methylamine"],
        "methodname": "hf",
        "args": {
            "basis": "def2-svp",
            },
        "ref_energy":  -95.14430698801766 ,
        "ref_gradient":  [[ 0.0040058528,  -0.008936694 ,   0.0057602451],
                          [ 0.0048293759,  -0.0064352169,   0.0031220897],
                          [ 0.0010464924,   0.0009827615,   0.0020279838],
                          [-0.0081563885,   0.0046904509,   0.0001404389],
                          [-0.0007015473,  -0.0019697556,  -0.001359591 ],
                          [-0.0044886294,  -0.0008634602,  -0.0124732382],
                          [ 0.0034648442,   0.0125319143,   0.0027820718]]
    }, id='methylamine_rhf_vacuum_spherical'),

    pytest.param({
        "geom": __geoms["methylamine"],
        "methodname": "svwn",
        "args": {
            "basis": "def2-svp",
            "dft_spherical_points": 590,
            "dft_radial_points": 100,
            },
        "ref_energy":  -95.27059071129463,
        "ref_gradient":  [[-0.0042079931,  0.0110837226, -0.0075103912],
                          [ 0.021872753 , -0.015601343 ,  0.0022902887],
                          [-0.0053019608,  0.0011495745, -0.0084300329],
                          [-0.0168052084, -0.000320317 ,  0.0090010852],
                          [-0.0001021064,  0.0098862187,  0.0015760428],
                          [ 0.0039998826, -0.000186438 ,  0.0048714784],
                          [ 0.0005446332, -0.0060114179, -0.0017984709]]
    }, id='methylamine_rsvwn_100_590_vacuum_spherical'),

    pytest.param({
        "geom": __geoms["methylamine_cation"],
        "methodname": "svwn",
        "args": {
            "basis": "def2-svp",
            "dft_spherical_points": 590,
            "dft_radial_points": 100,
            "reference": "uhf",
            },
        "ref_energy":  -94.90671106428948,
        "ref_gradient":  [[-0.0747781353,  0.0012338152,  0.0377768587],
                          [ 0.0577686097,  0.0118455182, -0.0401835875],
                          [-0.0104786522, -0.0064496597, -0.0091122221],
                          [-0.0125291386, -0.0137666696,  0.0184633873],
                          [-0.000168829 ,  0.0108593863,  0.01069133  ],
                          [ 0.028840043 ,  0.012876872 ,  0.0080063547],
                          [ 0.0113461023, -0.0165992621, -0.0256421211]]
    }, id='methylamine_usvwn_100_590_vacuum_spherical'),

    pytest.param({
        "geom": __geoms["methylamine"],
        "methodname": "blyp",
        "args": {
            "basis": "def2-svp",
            "dft_spherical_points": 590,
            "dft_radial_points": 100,
            },
        "ref_energy":  -95.72487707870928,
        "ref_gradient":  [[ 0.0127674565,  0.0057335608, -0.0116499233],
                          [ 0.0055262857, -0.0073053824,  0.0035248171],
                          [-0.005638603 ,  0.0008815796, -0.0085831624],
                          [-0.0175365195, -0.0002660706,  0.0093324223],
                          [-0.0001518177,  0.0101017388,  0.0019746095],
                          [ 0.0044823102, -0.0012587237,  0.0064942249],
                          [ 0.0005508877, -0.0078867025, -0.0010929882]]
    }, id='methylamine_rblyp_vacuum_spherical'),

    pytest.param({
        "geom": __geoms["methylamine_cation"],
        "methodname": "blyp",
        "args": {
            "basis": "def2-svp",
            "dft_spherical_points": 590,
            "dft_radial_points": 100,
            "reference": "uhf",
            },
        "ref_energy":  -95.3888474306933,
        "ref_gradient":  [[-0.0552469372, -0.0106376991,  0.0380240886],
                          [ 0.0321608685,  0.0249975328, -0.038409722 ],
                          [-0.0091110922, -0.0056762757, -0.0080018817],
                          [-0.0101799653, -0.0141288954,  0.0175679031],
                          [-0.0000751885,  0.0094926445,  0.0093506282],
                          [ 0.029649432 ,  0.0121691824,  0.0069369835],
                          [ 0.0128028826, -0.0162164895, -0.0254679997]]
    }, id='methylamine_ublyp_vacuum_spherical'),

    pytest.param({
        "geom": __geoms["methylamine"],
        "methodname": "b3lyp",
        "args": {
            "basis": "def2-svp",
            "dft_spherical_points": 590,
            "dft_radial_points": 100,
            },
        "ref_energy":  -95.78351424952231,
        "ref_gradient":  [[ 0.0076817599,  0.0020692916, -0.0057947484],
                          [ 0.0070551461, -0.0080576472,  0.0033875956],
                          [-0.0030799965,  0.0009953079, -0.0044523463],
                          [-0.0139792364,  0.0017055423,  0.005766114 ],
                          [-0.0004310461,  0.0054397525,  0.0006354435],
                          [ 0.0013863736, -0.0010576975,  0.0002627017],
                          [ 0.0013669995, -0.0010945497,  0.0001952399]]
    }, id='methylamine_rb3lyp_vacuum_spherical'),

    pytest.param({
        "geom": __geoms["methylamine_cation"],
        "methodname": "b3lyp",
        "args": {
            "basis": "def2-svp",
            "dft_spherical_points": 590,
            "dft_radial_points": 100,
            "reference": "uhf",
            },
        "ref_energy":  -95.43805064558163,
        "ref_gradient":  [[-0.0596547165, -0.0157778335,  0.0448214444],
                          [ 0.0276504454,  0.0264348671, -0.0373476431],
                          [-0.0055176907, -0.0047137131, -0.0037513197],
                          [-0.0046882016, -0.0125150164,  0.013317778 ],
                          [-0.0000338852,  0.0044822062,  0.0067611333],
                          [ 0.0280176243,  0.0126632811,  0.0013651831],
                          [ 0.0142264242, -0.0105737913, -0.0251665761]]
    }, id='methylamine_ub3lyp_vacuum_spherical'),

    pytest.param({
        "geom": __geoms["methylamine"],
        "methodname": "m06",
        "args": {
            "basis": "def2-svp",
            "dft_spherical_points": 590,
            "dft_radial_points": 100,
            },
        "ref_energy":  -95.70796917156567,
        "ref_gradient":  [[ 0.0022244308,  0.0034835216, -0.00419669  ],
                          [ 0.0128967341, -0.0110354976,  0.0029576458],
                          [-0.0030004485,  0.0012345365, -0.0039862551],
                          [-0.0135249132,  0.0019323333,  0.0053328599],
                          [-0.0007170777,  0.0050620552,  0.000397552 ],
                          [ 0.0007445523, -0.0008692979, -0.0008460016],
                          [ 0.0013767221,  0.0001923487,  0.0003408889]]
    }, id='methylamine_rm06_vacuum_spherical'),

    pytest.param({
        "geom": __geoms["methylamine_cation"],
        "methodname": "m06",
        "args": {
            "basis": "def2-svp",
            "dft_spherical_points": 590,
            "dft_radial_points": 100,
            "reference": "uhf",
            },
        "ref_energy":  -95.35812516766704,
        "ref_gradient":  [[-0.0601192783, -0.0136144738,  0.0431677739],
                          [ 0.0279783072,  0.0227113845, -0.0342771294],
                          [-0.0047710394, -0.0039605914, -0.0031205887],
                          [-0.0033875789, -0.0110793165,  0.0113954085],
                          [-0.0001443984,  0.0037976122,  0.0057474575],
                          [ 0.026788119 ,  0.0121357601,  0.0011760326],
                          [ 0.0136558688, -0.0099903752, -0.0240889545]]
    }, id='methylamine_um06_vacuum_spherical'),

    pytest.param({
        "geom": __geoms["methylamine"],
        "methodname": "pw6b95",
        "args": {
            "basis": "def2-svp",
            "dft_spherical_points": 590,
            "dft_radial_points": 100,
            },
        "ref_energy": -95.891979193935711,
        "ref_gradient": [[ 0.0035623300,  0.0011535385, -0.0028500117],
                         [ 0.0092168605, -0.0094680146,  0.0034988231],
                         [-0.0011688988,  0.0010378132, -0.0014311974],
                         [-0.0116414462,  0.0032048803,  0.0032445416],
                         [-0.0005877578,  0.0020011832, -0.0003293703],
                         [-0.0009796034, -0.0011343016, -0.0035287499],
                         [ 0.0015980964,  0.0032045977,  0.0013964163]]
    }, id='methylamine_rpw6b95_vacuum_spherical'),

    pytest.param({
        "geom": __geoms["methylamine_cation"],
        "methodname": "pw6b95",
        "args": {
            "basis": "def2-svp",
            "dft_spherical_points": 590,
            "dft_radial_points": 100,
            "reference": "uhf",
            },
        "ref_energy":  -95.545516905743852,
        "ref_gradient":  [[-0.0641789973, -0.0162144750,  0.0475569199],
                          [ 0.0298746048,  0.0250976880, -0.0373331997],
                          [-0.0036131281, -0.0046448509, -0.0008717052],
                          [-0.0024133635, -0.0111916599,  0.0109828317],
                          [-0.0001348544,  0.0011744476,  0.0057748601],
                          [ 0.0259943907,  0.0125974876, -0.0019683850],
                          [ 0.0144707663, -0.0068188042, -0.0241408920]]
    }, id='methylamine_upw6b95_vacuum_spherical'),

    pytest.param({
        "geom": __geoms["methylamine"],
        "methodname": "b97m-v",
        "args": {
            "basis": "def2-svp",
            "reference": "rhf",
            "dft_spherical_points": 590,
            "dft_radial_points": 100,
            },
        "ref_energy":  -95.76449121751773,
        "ref_gradient":  [[ 0.0055922095,  0.0024444584, -0.0050355658],
                          [ 0.0097374657, -0.0103601528,  0.0040118796],
                          [-0.0010394312,  0.0010197292, -0.0010268594],
                          [-0.0124738259,  0.0035498608,  0.0033752089],
                          [-0.0006911452,  0.0015901425, -0.0003744161],
                          [-0.0021707994, -0.0018294385, -0.0035510444],
                          [ 0.0010455264,  0.0035854005,  0.0026007971]]
    }, id='methylamine_ub97m-v_vacuum_spherical'),

    pytest.param({
        "geom": __geoms["methylamine_cation"],
        "methodname": "b97m-v",
        "args": {
            "basis": "def2-svp",
            "dft_spherical_points": 590,
            "dft_radial_points": 100,
            "reference": "uhf",
            },
        "ref_energy":  -95.4171292692646,
        "ref_gradient":  [[-0.0630133701, -0.0180904545,  0.0485969564],
                          [ 0.0286855434,  0.0253307189, -0.0369193625],
                          [-0.0032461931, -0.0045893328,  0.0002275563],
                          [-0.0018300779, -0.0102145822,  0.009825065 ],
                          [-0.0004492454,  0.0000803324,  0.0055558165],
                          [ 0.0249685327,  0.0122370494, -0.0039408211],
                          [ 0.0148848105, -0.0047537312, -0.0233452105]]
    }, id='methylamine_ub97m-v_vacuum_spherical'),

    pytest.param({
        "geom": __geoms["methylamine"],
        "methodname": "hf",
        "args": {
            "basis": "def2-tzvp",
            "reference": "rhf",
            "cuest_pcm": True,
            "cuest_pcm_dielectric": 80.0,
            },
        "ref_energy":  -9.5260345723745388E+01,
        "ref_gradient": [[ 0.01061845, -0.01408918,  0.00682433],
                         [-0.00613936, -0.00073797,  0.00383922],
                         [ 0.003623  ,  0.00051632,  0.00620779],
                         [-0.00312115,  0.00625346, -0.00383818],
                         [-0.0007377 , -0.00683027, -0.002212  ],
                         [-0.00722152, -0.00112063, -0.01518268],
                         [ 0.00297828,  0.01600828,  0.00436151]]
    }, id='methylamine_hf_pcm_spherical'),

    pytest.param({
        "geom": __geoms["methylamine"],
        "methodname": "hf",
        "args": {
            "basis": "def2-tzvp",
            "reference": "uhf",
            "cuest_pcm": True,
            "cuest_pcm_dielectric": 80.0,
            },
        "ref_energy":  -9.5260345723745388E+01,
        "ref_gradient": [[ 0.01061845, -0.01408918,  0.00682433],
                         [-0.00613936, -0.00073797,  0.00383922],
                         [ 0.003623  ,  0.00051632,  0.00620779],
                         [-0.00312115,  0.00625346, -0.00383818],
                         [-0.0007377 , -0.00683027, -0.002212  ],
                         [-0.00722152, -0.00112063, -0.01518268],
                         [ 0.00297828,  0.01600828,  0.00436151]]
    }, id='methylamine_hf_pcm_spherical'),
])

def test_cuest_scf(inp, request):
    psi4.core.set_num_threads(4)

    molecule = psi4.geometry(inp["geom"])

    # Default generic options
    psi4.set_options({
        # FD gradient options
        'points': 5,
        'disp_size': 0.005,
        'fd_project' : False,
        'scf_type': 'df',
        'dft_nuclear_scheme': 'stratmann', # To get cuEST and Psi4 to agree exactly for DFT
        'df_basis_scf': 'def2-universal-JKFIT',
        'maxiter': 300,
        'e_convergence': 10, 
        'd_convergence': 9,
        "puream": True,
        "reference": "rhf",
        'use_cuest': True,
        'cuest_mixed_precision': False,
    })

    # Override with the test-specific options from above
    psi4.set_options(inp['args'])

    G, wfn = psi4.gradient(inp['methodname'], return_wfn=True)
    assert psi4.compare_values(inp["ref_energy"], wfn.get_energies("Total Energy"), 5e-6, f'{request.node.callspec.id} energy')
    assert psi4.compare_values(inp["ref_gradient"], G, 5e-5, f'{request.node.callspec.id} gradient')


# ===========================================================================
# Molecular symmetry (non-C1) in the cuEST DFT XC path
# ===========================================================================
#
# The methylamine geometry used above is effectively C1, so none of those cases
# exercise a symmetry-blocked wavefunction.  cuEST's XC potential code
# (RV/UV::compute_V) works in the dense AO basis, while Va_/Vb_ are SO-basis,
# symmetry-blocked matrices; sizing the GPU buffers off ret[i]->size() picked up
# only the first irrep block, and nothing back-transformed AO -> SO.  Any
# nirrep() > 1 DFT wavefunction therefore died inside cuestXCPotentialCompute.
#
# That was worked around in the driver by forcing C1 for cuEST DFT, and is now
# fixed properly in v.cc, so the workaround is gone.  These tests pin both
# halves of that: cuEST must run in the *native* point group (not be silently
# demoted to C1), and must give the same answer as the same geometry in C1.

__symm_geoms = {
    "water": """
    0 1
    O
    H 1 0.9584
    H 1 0.9584 2 104.45
    """,
    "water_cation": """
    1 2
    O
    H 1 0.9584
    H 1 0.9584 2 104.45
    """,
}

# Reference energies from regular (non-cuEST) Psi4 in the native C2v point
# group, same options as the test below.
@pytest.mark.quick
@uusing("cuest")
@uusing("cuda_cc8")
@pytest.mark.parametrize("inp", [

    pytest.param({
        "geom": "water",
        "methodname": "hf",
        "reference": "rhf",
        "ref_energy": -75.9608932798,
    }, id='water_c2v_rhf'),

    pytest.param({
        "geom": "water",
        "methodname": "b3lyp",
        "reference": "rhf",
        "ref_energy": -76.3581721562,
    }, id='water_c2v_rb3lyp'),

    # VV10 functional: exercises the nonlocal-correlation branch of compute_V,
    # which folds its contribution into the same AO -> SO round trip.
    pytest.param({
        "geom": "water",
        "methodname": "b97m-v",
        "reference": "rhf",
        "ref_energy": -76.3299936452,
    }, id='water_c2v_rb97m-v'),

    pytest.param({
        "geom": "water_cation",
        "methodname": "hf",
        "reference": "uhf",
        "ref_energy": -75.5622929120,
    }, id='water_cation_c2v_uhf'),

    # Open shell: UV::compute_V is a separate code path from RV::compute_V and
    # was broken in the same way (two Vxc buffers instead of one).
    pytest.param({
        "geom": "water_cation",
        "methodname": "b3lyp",
        "reference": "uhf",
        "ref_energy": -75.9020172059,
    }, id='water_cation_c2v_ub3lyp'),

    pytest.param({
        "geom": "water_cation",
        "methodname": "m06",
        "reference": "uhf",
        "ref_energy": -75.8633697758,
    }, id='water_cation_c2v_um06'),

    pytest.param({
        "geom": "water_cation",
        "methodname": "b97m-v",
        "reference": "uhf",
        "ref_energy": -75.8758535607,
    }, id='water_cation_c2v_ub97m-v'),
])
def test_cuest_symmetry(inp, request):
    """cuEST must support non-C1 wavefunctions and be invariant to symmetry."""
    psi4.core.set_num_threads(4)

    psi4.set_options({
        'scf_type': 'df',
        'df_basis_scf': 'def2-universal-JKFIT',
        'basis': 'def2-svp',
        'puream': True,
        'dft_nuclear_scheme': 'stratmann',  # To get cuEST and Psi4 to agree exactly for DFT
        'dft_radial_points': 100,
        'dft_spherical_points': 590,
        'd_convergence': 9,
        'maxiter': 300,
        'use_cuest': True,
        'reference': inp['reference'],
    })

    mol_symm = psi4.geometry(__symm_geoms[inp['geom']])
    mol_symm.update_geometry()

    # Identical geometry with symmetry switched off. This is exactly what the
    # removed scf_helper() workaround used to do behind the user's back, so the
    # two gradients are directly comparable element by element.
    mol_c1 = mol_symm.clone()
    mol_c1.reset_point_group('c1')
    mol_c1.update_geometry()

    G_symm, wfn_symm = psi4.gradient(inp['methodname'], molecule=mol_symm, return_wfn=True)
    G_c1, wfn_c1 = psi4.gradient(inp['methodname'], molecule=mol_c1, return_wfn=True)

    tid = request.node.callspec.id

    # (1) cuEST really ran in the native point group. Before the v.cc fix, DFT
    # here was silently demoted to C1 (nirrep 1) by the driver workaround, so
    # this is the assertion that pins the workaround's removal.
    assert compare_integers(4, wfn_symm.nirrep(), f'{tid} nirrep')
    assert compare_strings('c2v', mol_symm.point_group().symbol(), f'{tid} point group')
    assert compare_integers(1, wfn_c1.nirrep(), f'{tid} C1 reference nirrep')

    # (2) Absolute energy against a non-cuEST Psi4 reference in native symmetry.
    assert compare_values(inp['ref_energy'], wfn_symm.energy(), 5e-6, f'{tid} energy')

    # (3) Symmetry invariance. These are the same physics, so agreement should
    # be near machine precision; observed worst case is ~1e-9 (E) and ~3e-10
    # (gradient) across these cases.
    assert compare_values(wfn_c1.energy(), wfn_symm.energy(), 1e-7, f'{tid} energy symmetry invariance')
    assert compare_values(np.array(G_c1), np.array(G_symm), 1e-7, f'{tid} gradient symmetry invariance')


# ===========================================================================
# DFT code paths that cuEST does not implement must fail loudly
# ===========================================================================
#
# Under USE_CUEST the quadrature grid is built on the GPU and the CPU-side grid
# post-processing is skipped entirely (cubature.cc: "if (is_cuest) return;"
# before postProcess()), so grid_->blocks() is empty and functional_workers_ is
# never populated. compute_V()/compute_gradient() have cuEST branches and are
# fine, but compute_Vx(), compute_fock_derivatives(), SAP::compute_V() and
# set_grac_shift() do not -- they used to fall through to the CPU code and index
# an empty std::vector, terminating the process with SIGSEGV.
#
# Note that simply populating functional_workers_ is NOT a fix: with no grid
# blocks the CPU routines integrate over nothing and silently return a *zero* XC
# contribution. These paths must therefore raise, not "work".
#
# Each case runs in a subprocess on purpose. If the guard regresses, the failure
# mode is a segfault, which would otherwise take down the whole pytest session
# and report nothing useful; as a subprocess it shows up as a negative return
# code that we can assert on and describe.

_UNSUPPORTED_ROUTE_SCRIPT = """
import sys
import psi4

route, outfile = sys.argv[1], sys.argv[2]
psi4.core.set_output_file(outfile, False)
psi4.set_options({
    'scf_type': 'df',
    'df_basis_scf': 'def2-universal-JKFIT',
    'basis': 'def2-svp',
    'puream': True,
    'reference': 'rhf',
    'dft_nuclear_scheme': 'stratmann',
    'd_convergence': 8,
    'save_jk': True,
    'use_cuest': True,
})
mol = psi4.geometry('''
0 1
O
H 1 0.9584
H 1 0.9584 2 104.45
symmetry c1
''')

try:
    if route == 'tdscf':
        _, wfn = psi4.energy('b3lyp', molecule=mol, return_wfn=True)
        psi4.set_options({'tdscf_states': 2})
        psi4.tdscf(wfn, states=2)
    elif route == 'polarizability':
        psi4.properties('b3lyp', properties=['DIPOLE_POLARIZABILITIES'], molecule=mol)
    elif route == 'sap_guess':
        psi4.set_options({'guess': 'sap'})
        psi4.energy('b3lyp', molecule=mol)
    elif route == 'grac':
        psi4.set_options({'dft_grac_shift': 0.5})
        psi4.energy('b3lyp', molecule=mol)
    else:
        raise SystemExit('unknown route ' + route)
except Exception as exc:
    print('CAUGHT ' + type(exc).__name__ + ': ' + ' '.join(str(exc).split()))
else:
    print('COMPLETED WITHOUT ERROR')
"""


@uusing("cuest")
@uusing("cuda_cc8")
@pytest.mark.parametrize("route,expected", [
    # TDDFT and CPHF response properties both go through RV::compute_Vx.
    pytest.param('tdscf', 'XC response kernel', id='tdscf'),
    pytest.param('polarizability', 'XC response kernel', id='polarizability'),
    pytest.param('sap_guess', 'SAP guess', id='sap_guess'),
    pytest.param('grac', 'GRAC asymptotic correction', id='grac'),
])
def test_cuest_unsupported_route(route, expected, tmp_path):
    """Unimplemented cuEST DFT paths raise a clear error instead of segfaulting."""
    script = tmp_path / "cuest_unsupported_route.py"
    script.write_text(_UNSUPPORTED_ROUTE_SCRIPT)

    proc = subprocess.run(
        [sys.executable, str(script), route, str(tmp_path / "psi.out")],
        capture_output=True, text=True, timeout=1800,
    )

    # A negative return code means killed by signal; -11 (SIGSEGV) is precisely
    # the pre-guard behaviour this test exists to prevent.
    assert proc.returncode == 0, (
        f"cuEST '{route}' subprocess exited with returncode {proc.returncode} "
        f"(negative means killed by a signal; -11 is SIGSEGV).\n"
        f"stdout:\n{proc.stdout}\nstderr tail:\n{proc.stderr[-2000:]}"
    )

    # Must actually have raised, rather than quietly returning XC-free numbers.
    assert 'CAUGHT' in proc.stdout, (
        f"cuEST '{route}' did not raise. A silent success here means the CPU "
        f"fallback ran with an empty grid and returned a zero XC contribution.\n"
        f"stdout:\n{proc.stdout}"
    )

    # And the error must be the informative cuEST guard, naming the operation.
    assert 'has no cuEST implementation' in proc.stdout, f"unexpected error:\n{proc.stdout}"
    assert expected in proc.stdout, (
        f"cuEST '{route}' raised, but not from the expected guard "
        f"(wanted {expected!r}).\nstdout:\n{proc.stdout}"
    )
