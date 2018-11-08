import psi4

Ne = psi4.geometry("""
0 1
Ne 0.0 0.0 0.0
symmetry c1
""")
psi4.set_options({'basis': 'cc-pvdz', 'freeze_core': 'false'})

_, wfn = psi4.energy('ccsd', return_wfn=True, molecule=Ne)
amps = wfn.get_amplitudes()
for k, v in amps.items():
    print(k, type(v).__name__)

TIjAb = amps['tIjAb'].to_array()
TIA = amps['tIA'].to_array()
tau_IjAb = TIjAb + np.einsum("ia,jb->ijab", TIA, TIA)

mints = psi4.core.MintsHelper(wfn.basisset())
D = mints.mo_eri(
    wfn.Ca_subset("AO", "OCC"), wfn.Ca_subset("AO", "VIR"), wfn.Ca_subset("AO", "OCC"),
    wfn.Ca_subset("AO", "VIR")).to_array()
D = D.swapaxes(1, 2)

RHF_ccsd_corr_e = 2 * np.einsum("ijab,ijab->", tau_IjAb, D) - np.einsum("ijab,ijba->", tau_IjAb, D)
psi4.compare_values(RHF_ccsd_corr_e, psi4.get_variable('CCSD CORRELATION ENERGY'), 8, "RHF CCSD CORRELATION ENERGY")

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

E2AA = 0.5*np.einsum("IJAB,IJAB->", tauIJAB, D_IJAB)
E2BB = 0.5*np.einsum("ijab,ijab->", tauijab, D_ijab)
E2AB = np.einsum("IjAb,IjAb->", tauIjAb, D_IjAb)

UHF_ccsd_corr_e = E2AA + E2BB + E2AB
psi4.compare_values(UHF_ccsd_corr_e, psi4.get_variable("CCSD CORRELATION ENERGY"), 8, "UHF CCSD CORRELATION ENERGY")

# END UHF

psi4.core.clean()
