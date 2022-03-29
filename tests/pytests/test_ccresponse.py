import numpy as np
import pytest
import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api]

def test_cc_polaroptrot():
    # cc29 + polarizabilities + tensors
    mol = psi4.geometry("""
        O     -0.028962160801    -0.694396279686    -0.049338350190
        O      0.028962160801     0.694396279686    -0.049338350190
        H      0.350498145881    -0.910645626300     0.783035421467
        H     -0.350498145881     0.910645626300     0.783035421467
        symmetry c1
        noreorient
        """)

    mw = sum([mol.mass(i) for i in range(0,4)])
    c = psi4.constants.c
    h = psi4.constants.h
    h2j = psi4.constants.hartree2J
    Na = psi4.constants.na
    me = psi4.constants.me
    hbar = psi4.constants.hbar 
    prefactor = -72E6 * hbar**2 * Na / c**2 / me**2 / mw / 3.0
    au_589 = c * h * 1E9 / h2j / 589
    au_355 = c * h * 1E9 / h2j / 355
    
    psi4.set_options({"basis":"cc-pVDZ"})
    omega = [589, 355, 'nm']
    psi4.set_options({'gauge': 'both', 'omega': omega, 'freeze_core': True, 'r_convergence': 10})

    e,wfn = psi4.properties('CCSD',properties=["polarizability","rotation"],return_wfn=True)

    pol_589      =      8.92296  #TEST
    pol_355      =      9.16134  #TEST
    rotlen_589   =    -76.93388  #TEST
    rotvel_589   =    850.24712  #TEST
    rotmvg_589   =   -179.08820  #TEST
    rotlen_355   =   -214.73273  #TEST
    rotvel_355   =    384.88786  #TEST
    rotmvg_355   =   -644.44746  #TEST

    polr_589     =  np.trace(wfn.variable("CCSD DIPOLE POLARIZABILITY TENSOR @ 589NM").to_array()) / 3.0                 #RESULT 
    polr_355     =  np.trace(wfn.variable("CCSD DIPOLE POLARIZABILITY TENSOR @ 355NM").to_array()) / 3.0                 #RESULT
    rotlenr_589  =  np.trace(wfn.variable("CCSD OPTICAL ROTATION TENSOR (LEN) @ 589NM").to_array()) * prefactor * au_589 #RESULT 
    rotvelr_589  =  np.trace(wfn.variable("CCSD OPTICAL ROTATION TENSOR (VEL) @ 589NM").to_array()) * prefactor          #RESULT 
    rotmvgr_589  =  np.trace(wfn.variable("CCSD OPTICAL ROTATION TENSOR (MVG) @ 589NM").to_array()) * prefactor          #RESULT 
    rotlenr_355  =  np.trace(wfn.variable("CCSD OPTICAL ROTATION TENSOR (LEN) @ 355NM").to_array()) * prefactor * au_355 #RESULT 
    rotvelr_355  =  np.trace(wfn.variable("CCSD OPTICAL ROTATION TENSOR (VEL) @ 355NM").to_array()) * prefactor          #RESULT 
    rotmvgr_355  =  np.trace(wfn.variable("CCSD OPTICAL ROTATION TENSOR (MVG) @ 355NM").to_array()) * prefactor          #RESULT 

    polt_589 = np.array([[ 5.44184081e+00, -8.30270852e-01,  9.39878004e-14],
     [-8.30270852e-01,  1.27518991e+01, -6.94720506e-14],
     [ 9.39878004e-14, -6.94720506e-14,  8.57516214e+00]]) #TEST
    polt_355 = np.array([[ 5.57913165e+00, -8.66951957e-01,  1.02575190e-13],
     [-8.66951957e-01,  1.31256109e+01, -7.71806997e-14],
     [ 1.02575190e-13, -7.71806997e-14,  8.77928702e+00]]) #TEST
    rottvel_0 = np.array([[ 3.63702792e-01, -3.62754886e-01, -3.37227243e-14],
     [-1.16320923e-01,  1.20784891e-02,  1.76483272e-14],
     [ 1.90062651e-14, -2.33369227e-15, -3.92022202e-01]]) #TEST
    rottlen_589 = np.array([[ 1.41679366e-01, -7.07140738e-02, -2.26060131e-14],
     [-2.72710108e-01,  1.28928522e-03,  1.01578936e-14],
     [ 3.74236371e-14,  3.88522300e-15, -1.27276912e-01]]) #TEST
    rottvel_589 = np.array([[ 3.75624000e-01, -3.74245928e-01, -3.61489832e-14],
     [-1.39409385e-01,  1.26867931e-02,  1.94630079e-14],
     [ 2.33993608e-14, -1.61361648e-15, -4.01726048e-01]]) #TEST
    rottmvg_589 = np.array([[ 1.19212077e-02, -1.14910412e-02, -2.42625886e-15],
     [-2.30884623e-02,  6.08303937e-04,  1.81468070e-15],
     [ 4.39309571e-15,  7.20075789e-16, -9.70384612e-03]]) #TEST

    assert psi4.compare_values(pol_589, polr_589,  3, "CCSD polarizability @ 589nm") #TEST
    assert psi4.compare_values(pol_355, polr_355,  3, "CCSD polarizability @ 355nm") #TEST
    assert psi4.compare_values(rotlen_589, rotlenr_589,  3, "CCSD rotation @ 589nm in length gauge") #TEST
    assert psi4.compare_values(rotvel_589, rotvelr_589,  3, "CCSD rotation @ 589nm in velocity gauge") #TEST
    assert psi4.compare_values(rotmvg_589, rotmvgr_589,  3, "CCSD rotation @ 589nm in modified velocity gauge") #TEST
    assert psi4.compare_values(rotlen_355, rotlenr_355,  3, "CCSD rotation @ 355nm in length gauge") #TEST
    assert psi4.compare_values(rotvel_355, rotvelr_355,  3, "CCSD rotation @ 355nm in velocity gauge") #TEST
    assert psi4.compare_values(rotmvg_355, rotmvgr_355,  3, "CCSD rotation @ 355nm in modified velocity gauge") #TEST

    assert psi4.compare_values(polt_589, wfn.variable("CCSD DIPOLE POLARIZABILITY TENSOR @ 589NM"),  3, "CCSD polarizability tensor @ 589nm") #TEST
    assert psi4.compare_values(polt_355, wfn.variable("CCSD DIPOLE POLARIZABILITY TENSOR @ 355NM"),  3, "CCSD polarizability tensor @ 355nm") #TEST
    assert psi4.compare_arrays(rottvel_0, wfn.variable("CCSD OPTICAL ROTATION TENSOR (VEL) @ 0NM"), 3, "CCSD optrot tensor @ 0NM in velocity guage") #TEST
    assert psi4.compare_arrays(rottlen_589, wfn.variable("CCSD OPTICAL ROTATION TENSOR (LEN) @ 589NM"), 3, "CCSD optrot tensor @ 589NM in length guage") #TEST 
    assert psi4.compare_arrays(rottvel_589, wfn.variable("CCSD OPTICAL ROTATION TENSOR (VEL) @ 589NM"), 3, "CCSD optrot tensor @ 589NM in velocity guage") #TEST
    assert psi4.compare_arrays(rottmvg_589, wfn.variable("CCSD OPTICAL ROTATION TENSOR (MVG) @ 589NM"), 3, "CCSD optrot tensor @ 589NM in modified velocity guage") #TEST 


def test_cc_roa():
    # three tensors for ROA (polar, optrot, and dip-quad)
    mol = psi4.geometry("""
        O     -0.028962160801    -0.694396279686    -0.049338350190
        O      0.028962160801     0.694396279686    -0.049338350190
        H      0.350498145881    -0.910645626300     0.783035421467
        H     -0.350498145881     0.910645626300     0.783035421467
        symmetry c1
        noreorient
        """)
    
    psi4.set_options({"basis":"cc-pVDZ"})
    omega = [589, 'nm']
    psi4.set_options({'gauge': 'both', 'omega': omega, 'freeze_core': True, 'r_convergence': 10})

    e,wfn = psi4.properties('CCSD',properties=["roa_tensor"],return_wfn=True)

    polt_589 = np.array([[ 5.44184081e+00, -8.30270852e-01,  9.39878004e-14],
     [-8.30270852e-01,  1.27518991e+01, -6.94720506e-14],
     [ 9.39878004e-14, -6.94720506e-14,  8.57516214e+00]])
    rottvel_0 = np.array([[ 3.63702792e-01, -3.62754886e-01, -3.37227243e-14],
     [-1.16320923e-01,  1.20784891e-02,  1.76483272e-14],
     [ 1.90062651e-14, -2.33369227e-15, -3.92022202e-01]])
    rottlen_589 = np.array([[ 1.41679366e-01, -7.07140738e-02, -2.26060131e-14],
     [-2.72710108e-01,  1.28928522e-03,  1.01578936e-14],
     [ 3.74236371e-14,  3.88522300e-15, -1.27276912e-01]])
    rottvel_589 = np.array([[ 3.75624000e-01, -3.74245928e-01, -3.61489832e-14],
     [-1.39409385e-01,  1.26867931e-02,  1.94630079e-14],
     [ 2.33993608e-14, -1.61361648e-15, -4.01726048e-01]])
    rottmvg_589 = np.array([[ 1.19212077e-02, -1.14910412e-02, -2.42625886e-15],
     [-2.30884623e-02,  6.08303937e-04,  1.81468070e-15],
     [ 4.39309571e-15,  7.20075789e-16, -9.70384612e-03]])
    quad0 = np.array(
    [
    [[-2.52229282e-14, -4.33846468e-13,  3.90176525e+00],
     [-4.33846468e-13, -6.81906994e-15, -5.83659009e+00],
     [ 3.90176525e+00, -5.83659009e+00,  3.18177152e-14]]
    ,
    [[ 1.20834634e-13,  6.99638702e-14, -2.37995874e+00],
     [ 6.99638702e-14, -1.17603861e-13,  6.85651590e+00],
     [-2.37995874e+00,  6.85651590e+00, -3.99609068e-15]]
    ,
    [[-4.55024365e+00, -5.30335671e+00,  1.65990108e-13],
     [-5.30335671e+00, -1.35415335e+00, -8.75522529e-13],
     [ 1.65990108e-13, -8.75522529e-13,  5.90439700e+00]]
    ])

    assert psi4.compare_values(polt_589, wfn.variable("CCSD DIPOLE POLARIZABILITY TENSOR @ 589NM"),  3, "CCSD polarizability tensor @ 589nm") #TEST
    assert psi4.compare_arrays(rottvel_0, wfn.variable("CCSD OPTICAL ROTATION TENSOR (VEL) @ 0NM"), 3, "CCSD optrot tensor @ 0NM in velocity guage") #TEST
    assert psi4.compare_arrays(rottlen_589, wfn.variable("CCSD OPTICAL ROTATION TENSOR (LEN) @ 589NM"), 3, "CCSD optrot tensor @ 589NM in length guage") #TEST 
    assert psi4.compare_arrays(rottvel_589, wfn.variable("CCSD OPTICAL ROTATION TENSOR (VEL) @ 589NM"), 3, "CCSD optrot tensor @ 589NM in velocity guage") #TEST
    assert psi4.compare_arrays(rottmvg_589, wfn.variable("CCSD OPTICAL ROTATION TENSOR (MVG) @ 589NM"), 3, "CCSD optrot tensor @ 589NM in modified velocity guage") #TEST 
    assert psi4.compare_values(quad0, wfn.variable("CCSD QUADRUPOLE POLARIZABILITY TENSOR @ 589NM"),  3, "CCSD quadrupole tensor @ 589nm") #TEST
