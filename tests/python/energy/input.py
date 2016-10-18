import psi4

psi4.set_output_file("output.dat", False)

geom = psi4.geometry("""
He
He 1 5
""")

psi4.set_global_option("SCF_TYPE", "DIRECT")
psi4.set_global_option("BASIS", "STO-3G")
scf_e = psi4.energy('SCF')

psi4.compare_values(-5.6155679150795779, scf_e, 6, 'SCF Energy')
