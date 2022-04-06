import pytest

from utils import *

import psi4
import numpy as np
from addons import uusing, using

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]


# This will remove the 32 file in each iteration of the for loop 
#psi4_io.set_specific_retention(32,False)

ref = {
    'nobas': {
        'NUCLEAR REPULSION ENERGY': 0.587974678522222,
        'MP2-D DISPERSION CORRECTION ENERGY': -0.00053454696076,
        'MP2-D DISPERSION CORRECTION GRADIENT': np.array([
             0.000, 0.000,  0.00264290716623,
             0.000, 0.000, -0.00264290716623]).reshape((-1, 3)),
             # original MP2D had bad gradients: 0.001513831369
    },
    'cc-pVDZ': {
        'HF TOTAL ENERGY': -1.11639202566179,
        'MP2 TOTAL ENERGY': -1.1441651458266191,
        'MP2 TOTAL GRADIENT': np.array([
             0.000, 0.000, -6.96221095e-02,
             0.000, 0.000,  6.96221095e-02]).reshape((-1, 3)),
    },
    'cc-pVTZ': {
        'HF TOTAL ENERGY': -1.11872756109405,
        'MP2 TOTAL ENERGY': -1.1511940290102121,
        'MP2 TOTAL GRADIENT': np.array([
             0.000, 0.000, -7.43721409e-02,
             0.000, 0.000,  7.43721409e-02]).reshape((-1, 3)),
    },
}
for bas in ['cc-pVDZ', 'cc-pVTZ']:
    ref[bas]['MP2D TOTAL ENERGY'] = ref[bas]['MP2 TOTAL ENERGY'] + ref['nobas']['MP2-D DISPERSION CORRECTION ENERGY']
    ref[bas]['MP2 CORRELATION ENERGY'] = ref[bas]['MP2 TOTAL ENERGY'] - ref[bas]['HF TOTAL ENERGY']
    ref[bas]['MP2D TOTAL GRADIENT'] = ref[bas]['MP2 TOTAL GRADIENT'] + ref['nobas']['MP2-D DISPERSION CORRECTION GRADIENT']


@pytest.mark.parametrize("inp", [
    pytest.param({'driver': 'energy', 'name': 'Mp2', 'pv': 'MP2', 'options': {}}, id='mp2-energy'),
    pytest.param({'driver': 'energy', 'name': 'MP2-d', 'pv': 'MP2D', 'options': {}}, id='mp2d-energy', marks=using('mp2d')),
    pytest.param({'driver': 'gradient', 'name': 'Mp2', 'pv': 'MP2', 'options': {}}, id='mp2-gradient'),
    pytest.param({'driver': 'gradient', 'name': 'Mp2', 'pv': 'MP2', 'dertype': 0, 'options': {}}, id='mp2-gradient-findif'),
    pytest.param({'driver': 'gradient', 'name': 'MP2-d', 'pv': 'MP2D', 'options': {}}, id='mp2d-gradient', marks=using('mp2d')),
    pytest.param({'driver': 'gradient', 'name': 'MP2-d', 'pv': 'MP2D', 'dertype': 0, 'options': {}}, id='mp2d-gradient-findif', marks=using('mp2d')),
#    ('mp2mp2'),
#    ('mp2-dmp2'),
])
def test_dft_mp2(inp):
    """Adapted from dfmp2-2"""

    h2 = psi4.geometry("""
        0 1
        H
        H 1 1.7007535129120455
        units bohr
        """)

    basisset = "cc-pVDZ"
    psi4.set_options({
        'scf_type': 'df',
        'guess': 'sad',
        'd_convergence': 12,
        'e_convergence': 12,
        'points': 5,
        'basis': basisset,
    })
    psi4.set_options(inp['options'])
    kwargs = {'return_wfn': True}
    if inp.get('dertype') is not None:
        kwargs.update({'dertype': inp['dertype']})

    if inp['driver'] == 'energy':
        ene, wfn = psi4.energy(inp['name'], **kwargs)
    elif inp['driver'] == 'gradient':
        grad, wfn = psi4.gradient(inp['name'], **kwargs)
        ene = wfn.energy()  # simplifies testing logic

    #for pv, pvv in psi4.core.variables().items():
    #    print('CORE -- ', pv, pvv)
    #for pv, pvv in wfn.variables().items():
    #    print('WFN  -- ', pv, pvv)
    assert compare_values(ref['nobas']['NUCLEAR REPULSION ENERGY'], h2.nuclear_repulsion_energy(), 10, "Nuclear Repulsion Energy")

    for obj in [psi4.core, wfn]:
        for qcvar in ['MP2 TOTAL ENERGY',
                      'MP2 CORRELATION ENERGY']:
            assert compare_values(ref[basisset][qcvar], obj.variable(qcvar), 10, basisset + " " + qcvar)

        for qcvar in [inp['pv'] + ' TOTAL ENERGY',
                      'CURRENT ENERGY']:
            assert compare_values(ref[basisset][inp['pv'] + ' TOTAL ENERGY'], obj.variable(qcvar), 10, basisset + " " + qcvar)

        ## apparently not widespread
        #if inp['driver'] == 'gradient':
        #    for qcvar in [inp['pv'] + ' TOTAL GRADIENT',
        #                  'CURRENT GRADIENT']:
        #        assert compare_values(ref[basisset][inp['pv'] + ' TOTAL GRADIENT'], obj.variable(qcvar), 10, basisset + " " + qcvar)

        for qcvar in ['HF TOTAL ENERGY', 
                      'SCF TOTAL ENERGY', 
                      'SCF ITERATION ENERGY', 
                      'CURRENT REFERENCE ENERGY']:
            assert compare_values(ref[basisset]['HF TOTAL ENERGY'], obj.variable(qcvar), 10, basisset + " " + qcvar)

        if inp['pv'] == 'MP2D':
            for qcvar in ['DISPERSION CORRECTION ENERGY',
                          'MP2D DISPERSION CORRECTION ENERGY']:
                assert compare_values(ref['nobas']['MP2-D DISPERSION CORRECTION ENERGY'], obj.variable(qcvar), 10, basisset + " " + qcvar)

            if inp['driver'] == 'gradient' and 'dertype' not in inp:
                for qcvar in ['DISPERSION CORRECTION GRADIENT',
                              'MP2D DISPERSION CORRECTION GRADIENT']:
                    assert compare_arrays(ref['nobas']['MP2-D DISPERSION CORRECTION GRADIENT'], np.asarray(obj.variable(qcvar)), 10, basisset + " disp grad")

    for retrn in [ene,
                  wfn.energy()]:
        assert compare_values(ref[basisset][inp['pv'] + ' TOTAL ENERGY'], retrn, 10, basisset + " tot ene")

    if inp['driver'] == 'gradient':
        for retrn in [grad,
                      wfn.gradient()]:
            atol = 2.e-8 if 'dertype' in inp else 1.e-10
            assert compare_values(ref[basisset][inp['pv'] + ' TOTAL GRADIENT'], np.asarray(retrn), basisset + " tot grad", atol=atol)


#TABLE 14259 -1.155358302362078 0.7013114524160179


@uusing("mp2d")
def test_mp2d_opt():

    h2 = psi4.geometry("""
        0 1
        H
        H 1 R
        units bohr
        R = 1.7007535129120455
        """)

    psi4.set_options({
        'scf_type': 'df',
        'd_convergence': 12,
        'e_convergence': 12,
        'g_convergence': 'gau_verytight',
    })

    ene, wfn = psi4.optimize('mp2d/cc-pvdz', return_wfn=True, molecule=h2)

    assert compare_values(1.4259, h2.R, 'h2 bond length', atol=1.e-3)
