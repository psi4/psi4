import pytest

from utils import *
from addons import *

import psi4
import numpy as np
from qcengine.testing import using_mp2d

pytestmark = [pytest.mark.quick]

#! Density fitted MP2 energy of H2, using density fitted reference and automatic looping over cc-pVDZ and cc-pVTZ basis sets.
#! Results are tabulated using the built in table functions by using the default options and by specifiying the format.

# This will remove the 32 file in each iteration of the for loop 
#psi4_io.set_specific_retention(32,False)

ref = {
    'nobas': {
        'NUCLEAR REPULSION ENERGY': 0.587974678522222,
        'MP2-D DISPERSION CORRECTION ENERGY': -0.00053454696076,
        'MP2-D DISPERSION CORRECTION GRADIENT': np.array([
             0.000, 0.000,  0.001513831369,
             0.000, 0.000, -0.001513831369]).reshape((-1, 3)),
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
    pytest.param({'driver': 'energy', 'name': 'Mp2', 'pv': 'MP2', 'options': {}}, id='mp2 energy'),
    pytest.param({'driver': 'energy', 'name': 'MP2-d', 'pv': 'MP2D', 'options': {}}, id='mp2d energy', marks=using_mp2d),
    pytest.param({'driver': 'gradient', 'name': 'Mp2', 'pv': 'MP2', 'options': {}}, id='mp2 gradient'),
#    pytest.param({'driver': 'gradient', 'name': 'MP2-d', 'pv': 'MP2D', 'options': {}}, id='mp2d gradient' marks=using_mp2d),  # not working even in findif
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

    if inp['driver'] == 'energy':
        ene, wfn = psi4.energy(inp['name'], return_wfn=True)
    elif inp['driver'] == 'gradient':
        grad, wfn = psi4.gradient(inp['name'], return_wfn=True) #, dertype=0)
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

        for qcvar in ['HF TOTAL ENERGY', 
                      'SCF TOTAL ENERGY', 
                      'SCF ITERATION ENERGY', 
                      'CURRENT REFERENCE ENERGY']:
            assert compare_values(ref[basisset]['HF TOTAL ENERGY'], obj.variable(qcvar), 10, basisset + " " + qcvar)

        for qcvar in ['DISPERSION CORRECTION ENERGY']: #, 
                      #'MP2D DISPERSION CORRECTION ENERGY']:
            if inp['pv'] == 'MP2D':
                assert compare_values(ref['nobas']['MP2-D DISPERSION CORRECTION ENERGY'], obj.variable(qcvar), 10, basisset + " " + qcvar)

        if inp['pv'] == 'MP2D':
            assert compare_values(ref['nobas']['MP2-D DISPERSION CORRECTION ENERGY'], obj.variable('DISPERSION CORRECTION ENERGY'), 10, basisset + " disp ene")
            if inp['driver'] == 'gradient':
                #print(ref['nobas']['MP2-D DISPERSION CORRECTION GRADIENT'].shape)
                #print(np.asarray(obj.variable('DISPERSION CORRECTION GRADIENT')).shape)
                assert compare_arrays(ref['nobas']['MP2-D DISPERSION CORRECTION GRADIENT'], np.asarray(obj.variable('DISPERSION CORRECTION GRADIENT')), 10, basisset + " disp grad")
                pass

    for retrn in [ene,
                  wfn.energy()]:
        assert compare_values(ref[basisset][inp['pv'] + ' TOTAL ENERGY'], retrn, 10, basisset + " tot ene")

    if inp['driver'] == 'gradient':
        for retrn in [grad,
                      wfn.gradient()]:
            assert compare_arrays(ref[basisset][inp['pv'] + ' TOTAL GRADIENT'], np.asarray(retrn), 10, basisset + " tot grad")

#    assert 0

#TABLE 14259 -1.155358302362078 0.7013114524160179

