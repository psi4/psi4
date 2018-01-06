import psi4

psi4.set_output_file("output.dat", False)

ne = psi4.geometry("""
   Ne  0  0  0
   """)

psi4.set_options({'basis': '6-311g', 'scf_type': 'pk'})

scf_e, scf_wfn = psi4.energy('SCF', return_wfn=True)
psi4.fcidump(scf_wfn, fname='Ne.6311G.INTDUMP', oe_ints=['EIGENVALUES'])

psi4.compare_fcidumps('Ne.6311G.INTDUMP.ref', 'Ne.6311G.INTDUMP', 'Ne, 6-311g, integrals in FCIDUMP format')  # TEST
