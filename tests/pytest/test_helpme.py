import pytest
from addons import *

import numpy as np

@pytest.mark.smoke
@using_helpme
@pytest.mark.parametrize("tol,dtype", [
    pytest.param(1.e-8, np.float64, id='double'),
    pytest.param(1.e-4, np.float32, id='float'),
])
def test_helpme(tol, dtype):

    import helpme as pme

    def print_results(label, e, f, v):
        np.set_printoptions(precision=10, linewidth=100)
        print(label)
        print("Energy = {:16.10f}".format(e))
        print("Forces:")
        print(f)
        print("Virial:")
        print(v)
        print()

    expectedEnergy = 5.864957414;
    expectedForces = np.array([[-1.20630693, -1.49522843, 12.65589187],
                               [ 1.00695882,  0.88956328, -5.08428301],
                               [ 0.69297661,  1.09547848, -5.22771480],
                               [-2.28988057, -2.10832506, 10.18914165],
                               [ 0.81915340,  0.92013663, -6.43738026],
                               [ 0.97696467,  0.69833887, -6.09492437]]);
    expectedVirial = np.array([[0.65613058, 0.49091167, 0.61109732,
                                2.26906257, 2.31925449, -10.04901641]]);
    expectedPotential = np.array([[ 1.18119329, -0.72320559, -0.89641992, 7.58746515],
                                 [ 7.69247982, -1.20738468, -1.06662264, 6.09626260],
                                 [ 8.73449635, -0.83090721, -1.31352336, 6.26824317],
                                 [-9.98483179, -1.37283008, -1.26398385, 6.10859811],
                                 [-3.50591589, -0.98219832, -1.10328133, 7.71868137],
                                 [-2.39904512, -1.17142047, -0.83733677, 7.30806279]]);
    
    # Instantiate single or double precision PME object
    coords = np.array([
        [ 2.00000,  2.00000, 2.00000],
        [ 2.50000,  2.00000, 3.00000],
        [ 1.50000,  2.00000, 3.00000],
        [ 0.00000,  0.00000, 0.00000],
        [ 0.50000,  0.00000, 1.00000],
        [-0.50000,  0.00000, 1.00000]
    ], dtype=dtype)
    charges = np.array([[-0.834, 0.417, 0.417, -0.834, 0.417, 0.417]], dtype=dtype).T

    energy = 0
    forces = np.zeros((6,3), dtype=dtype)
    virial = np.zeros((1,6), dtype=dtype)
    potentialAndGradient = np.zeros((6,4), dtype=dtype)

    if dtype == np.float64:
        pmei = pme.PMEInstanceD()
        mat = pme.MatrixD
    elif dtype == np.float32:
        pmei = pme.PMEInstanceF()
        mat = pme.MatrixF
    else:
        raise KeyError('undef helPME mode')

    pmei.setup(1, 0.3, 5, 32, 32, 32, 332.0716, 1)
    pmei.set_lattice_vectors(20, 20, 20, 90, 90, 90, pmei.LatticeType.XAligned)
    # Compute just the energy
    print_results("Before pmei.compute_E_rec", energy, forces, virial)
    energy = pmei.compute_E_rec(0, mat(charges), mat(coords))
    print_results("After pmei.compute_E_rec", energy, forces, virial)
    # Compute the energy and forces
    energy = pmei.compute_EF_rec(0, mat(charges), mat(coords), mat(forces))
    print_results("After pmei.compute_EF_rec", energy, forces, virial)
    # Compute the energy, forces and virial
    energy = pmei.compute_EFV_rec(0, mat(charges), mat(coords), mat(forces), mat(virial))
    print_results("After pmei.compute_EFV_rec", energy, forces, virial)
    # Compute the reciprocal space potential and its gradient
    pmei.compute_P_rec(0, mat(charges), mat(coords), mat(coords), 1, mat(potentialAndGradient))
    print("Potential and its gradient:")
    print(potentialAndGradient, "\n")

    assert np.allclose([expectedEnergy], [energy], atol=tol)
    assert np.allclose(expectedForces, forces, atol=tol)
    assert np.allclose(expectedVirial, virial, atol=tol)
    assert np.allclose(expectedPotential, potentialAndGradient, atol=tol)
