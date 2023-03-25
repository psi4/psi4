#!/usr/bin/env python3
import psi4

#! UHF PCM TDSCF 
# comparison against other implementations:
#CPCM solvation model excitation energies of NH2 radical [eV]:
#(Same radii, epsilon, refraction index, probe radius in each program)
#
#state  PSI4    ORCA    G09
#1      2.43618 2.457   2.4411
#2      7.39834 7.203   7.3993
#3      9.11054 8.991   9.1777
#4      10.1063 9.848   10.1017
#5      10.2178 10.365  10.3305

exc_energies_uhf = [ #TEST
0.08972598844884663, #TEST
0.2719972189552112, #TEST
0.33525624703045037, #TEST
0.3713900898382711, #TEST
0.3762084431466903, #TEST
] #TEST

osc_strengths_uhf = [ #TEST
0.00414272407997719, #TEST
4.590977316768732e-28, #TEST
0.005715258102198367, #TEST
0.018432750255865125, #TEST
0.006434117688452513,#TEST
] #TEST

scf_ref = -55.5218518303635165 #TEST


mol = psi4.geometry("""
    0, 2
    N  0.000000000000     0.000000000000    -0.079859033927
    H  0.000000000000    -0.803611003426     0.554794694632
    H  0.000000000000     0.803611003426     0.554794694632
    symmetry c1
    units angstrom
    """)

psi4.set_options({
    "reference":      "uhf",
    "scf_type":       "pk",
    "basis":          "def2-SVP",
    "e_convergence":  10,
    "ddx":            True,
    "maxiter":        50,
    "tdscf_states":   5,
    #
    # "ddx_lmax": 15,
    # "ddx_n_lebedev": 590,
    # "ddx_solute_spherical_points": 590,
    # "ddx_solute_radial_points": 90,
    # "ddx_solvation_convergence": 1e-11,
    # "ddx_fmm": False,
    "ddx_model": "pcm",
    "ddx_solvent": "water",
    "ddx_radii_set": "uff",
    'ddx_radii_scaling': 1.2,
})

scf_energy = psi4.energy("TD-SCF")

r_calc = []
e_calc = []
for i in range(5):
    e_calc.append(variable(f'TD-HF ROOT 0 -> ROOT {i+1} EXCITATION ENERGY - A TRANSITION'))
    r_calc.append(variable(f'TD-HF ROOT 0 -> ROOT {i+1} OSCILLATOR STRENGTH (LEN) - A TRANSITION'))

print("SCF        ", scf_energy)
print("e_calc     ", e_calc)
print("r_calc     ", r_calc)
print("e_calc (eV)", e_calc * 27.211)
