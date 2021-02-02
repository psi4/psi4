#! API access to CCSD amplitudes

import numpy as np

import psi4

Ne = psi4.geometry("""
0 1
Ne 0.0 0.0 0.0
symmetry c1
""")
psi4.set_options({'basis': 'cc-pvdz', 'freeze_core': 'false'})

_, wfn = psi4.energy('ccsd', return_wfn=True, molecule=Ne)
amps = wfn.get_amplitudes()

TIjAb = amps['tIjAb'].to_array()
TIA = amps['tIA'].to_array()
tau_IjAb = TIjAb + np.einsum("ia,jb->ijab", TIA, TIA)

mints = psi4.core.MintsHelper(wfn.basisset())
D = mints.mo_eri(
    wfn.Ca_subset("AO", "OCC"), wfn.Ca_subset("AO", "VIR"), wfn.Ca_subset("AO", "OCC"),
    wfn.Ca_subset("AO", "VIR")).to_array()
D = D.swapaxes(1, 2)

RHF_ccsd_corr_e = 2 * np.einsum("ijab,ijab->", tau_IjAb, D) - np.einsum("ijab,ijba->", tau_IjAb, D)
psi4.compare_values(psi4.variable('CCSD CORRELATION ENERGY'), RHF_ccsd_corr_e, 8, "RHF CCSD CORRELATION ENERGY")

# END RHF

psi4.core.clean()

# UHF

psi4.set_options({
    'basis': 'cc-pvdz',
    'freeze_core': 'false',
    'reference': 'UHF',
})

_, wfn = psi4.energy('ccsd', return_wfn=True, molecule=Ne)
amps = wfn.get_amplitudes()
TIJAB = amps['tIJAB'].to_array()
Tijab = amps['tijab'].to_array()
TIjAb = amps['tIjAb'].to_array()
Tia = amps['tia'].to_array()
TIA = amps['tIA'].to_array()

tauIJAB = TIJAB + np.einsum('IA,JB->IJAB', TIA, TIA) - np.einsum("IA,JB->IJBA", TIA, TIA)
tauijab = Tijab + np.einsum('ia,jb->ijab', Tia, Tia) - np.einsum("ia,jb->ijba", Tia, Tia)
tauIjAb = TIjAb + np.einsum('IA,jb->IjAb', TIA, Tia)

CO = wfn.Ca_subset("AO", "OCC")
Co = wfn.Cb_subset("AO", "OCC")
CV = wfn.Ca_subset("AO", "VIR")
Cv = wfn.Cb_subset("AO", "VIR")

mints = psi4.core.MintsHelper(wfn.basisset())

D_IJAB = mints.mo_eri(CO, CV, CO, CV).to_array().swapaxes(1, 2)
D_ijab = mints.mo_eri(Co, Cv, Co, Cv).to_array().swapaxes(1, 2)
D_IjAb = mints.mo_eri(CO, CV, Co, Cv).to_array().swapaxes(1, 2)

E2AA = 0.5 * np.einsum("IJAB,IJAB->", tauIJAB, D_IJAB)
E2BB = 0.5 * np.einsum("ijab,ijab->", tauijab, D_ijab)
E2AB = np.einsum("IjAb,IjAb->", tauIjAb, D_IjAb)

UHF_ccsd_corr_e = E2AA + E2BB + E2AB
psi4.compare_values(psi4.variable("CCSD CORRELATION ENERGY"), UHF_ccsd_corr_e, 8, "UHF CCSD CORRELATION ENERGY")

# END UHF

psi4.core.clean()

# ROHF-semicanonical

CN = psi4.geometry("""
  units bohr
  0 2
  C  0.000000000000      0.000000000000      1.195736583480
  N  0.000000000000      0.000000000000     -1.024692078304
  symmetry c1
""")

psi4.set_options({
    'basis': 'cc-pvdz',
    'freeze_core': 'false',
    'reference': 'ROHF',
    'semicanonical': 'True',
})

_, wfn = psi4.energy('ccsd', return_wfn=True, molecule=CN)
amps = wfn.get_amplitudes()
TIJAB = amps['tIJAB'].to_array()
Tijab = amps['tijab'].to_array()
TIjAb = amps['tIjAb'].to_array()
Tia = amps['tia'].to_array()
TIA = amps['tIA'].to_array()

tauIJAB = TIJAB + np.einsum('IA,JB->IJAB', TIA, TIA) - np.einsum("IA,JB->IJBA", TIA, TIA)
tauijab = Tijab + np.einsum('ia,jb->ijab', Tia, Tia) - np.einsum("ia,jb->ijba", Tia, Tia)
tauIjAb = TIjAb + np.einsum('IA,jb->IjAb', TIA, Tia)

CO = wfn.Ca_subset("AO", "OCC")
Co = wfn.Cb_subset("AO", "OCC")
CV = wfn.Ca_subset("AO", "VIR")
Cv = wfn.Cb_subset("AO", "VIR")
fIA = psi4.core.Matrix.triplet(CO, wfn.Fa(), CV, True, False, False).to_array()
fia = psi4.core.Matrix.triplet(Co, wfn.Fb(), Cv, True, False, False).to_array()

mints = psi4.core.MintsHelper(wfn.basisset())

D_IJAB = mints.mo_eri(CO, CV, CO, CV).to_array().swapaxes(1, 2)
D_ijab = mints.mo_eri(Co, Cv, Co, Cv).to_array().swapaxes(1, 2)
D_IjAb = mints.mo_eri(CO, CV, Co, Cv).to_array().swapaxes(1, 2)

E1_A = np.einsum("IA,IA->", TIA, fIA)
E1_B = np.einsum("ia,ia->", Tia, fia)

E2AA = 0.5 * np.einsum("IJAB,IJAB->", tauIJAB, D_IJAB)
E2BB = 0.5 * np.einsum("ijab,ijab->", tauijab, D_ijab)
E2AB = np.einsum("IjAb,IjAb->", tauIjAb, D_IjAb)

ROHF_ccsd_corr_e = E2AA + E2BB + E2AB + E1_A + E1_B
psi4.compare_values(psi4.variable("CCSD CORRELATION ENERGY"), ROHF_ccsd_corr_e, 8, "ROHF CCSD CORRELATION ENERGY")
# END ROHF
