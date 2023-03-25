#!/usr/bin/env python3
import psi4

mol = psi4.geometry("""
    C   -0.502049    0.594262    0.000000
    H   -0.145376    1.098660    0.873652
    H   -1.572049    0.594275    0.000000
    Cl   0.084628    1.423927   -1.437034
    F   -0.052065   -0.678535    0.000000
    symmetry c1
    units angstrom
    """)

psi4.set_options({
    "ddx": True,
    "basis": "cc-pvdz",
    #
    "ddx_lmax": 10,
    "ddx_n_lebedev": 302,
    "ddx_solute_spherical_points": 302,
    "ddx_solute_radial_points": 75,
    "ddx_solvation_convergence": 1e-10,
    "ddx_model": "pcm",
    "ddx_solvent": "water",
    "ddx_radii_set": "uff",
})

def rundd(model, solvent_kappa=0.0):
    psi4.set_options({"ddx_model": model, "ddx_solvent_kappa": solvent_kappa})

    scf_e, wfn = psi4.energy('SCF', return_wfn=True, molecule=mol)
    dd_e = wfn.variable("DD SOLVATION ENERGY")
    return scf_e, dd_e


results = dict()
results["PCM    "] = rundd("pcm")
results["LPB0   "] = rundd("lpb", solvent_kappa=1e-6)
results["LPB    "] = rundd("lpb", solvent_kappa=0.1)
results["LPB8   "] = rundd("lpb", solvent_kappa=18)
results["COSMO  "] = rundd("cosmo")

for k in results:
    print(f"{k} {results[k][0]:15.8f} {results[k][1]:15.8f}")
