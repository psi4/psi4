import psi4

def test_chemps2():
    """chemps2/scf-n2"""
    #! dmrg-scf on N2

    N2 = psi4.geometry("""
      N       0.0000   0.0000   0.0000
      N       0.0000   0.0000   2.1180
    units au
    """)

    psi4.set_options({
    'basis': 'cc-pVDZ',
    'reference': 'rhf',
    'e_convergence': 1e-12,
    'd_convergence': 1e-12,

    'dmrg_irrep': 0,
    'dmrg_multiplicity': 1,
    'restricted_docc': [ 1 , 0 , 0 , 0 , 0 , 1 , 0 , 0 ],
    'active': [ 2 , 0 , 1 , 1 , 0 , 2 , 1 , 1 ],

    'dmrg_sweep_states': [   500,  1000,  1000 ],
    'dmrg_sweep_energy_conv': [ 1e-10, 1e-10, 1e-10 ],
    'dmrg_sweep_dvdson_rtol': [  1e-4,  1e-6,  1e-8 ],
    'dmrg_sweep_max_sweeps': [     5,     5,    10 ],
    'dmrg_sweep_noise_prefac': [  0.05,  0.05,   0.0 ],
    'dmrg_print_corr': True,
    'dmrg_mps_write': False,

    'dmrg_unitary_write': True,
    'dmrg_diis': True,
    'dmrg_scf_diis_thr': 1e-2,
    'dmrg_diis_write': True,

    'dmrg_excitation': 0,   # Ground state
    'dmrg_scf_state_avg': False,
    'dmrg_scf_active_space': 'NO',  # INPUT; NO; LOC
    'dmrg_local_init': False,
    })

    psi4.energy("dmrg-scf")

    ref_energy = -109.1035023353
    ans = psi4.compare_values(ref_energy, psi4.get_variable("CURRENT ENERGY"), 6, "DMRG Energy")
    assert ans
