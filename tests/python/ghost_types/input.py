#! Testing mass() and Z() for both ghost types: A (Gh(x) or @(x)) and B (extract_subsets())
import psi4

m_He = 4.00260325415 #TEST
m_Gh = 0.0           #TEST
Z_He = 2.0           #TEST
Z_Gh = 0.0           #TEST

frag = psi4.geometry("""
    He
    --
    He 1 1
""")

A_frag = psi4.geometry(""" 
    He
    --
    Gh(He) 1 1
""")

B_frag = frag.extract_subsets(1,2)

mix_frag = A_frag.extract_subsets(2,1) # two ghost atoms: 1 is type B, 2 is type A

def mol_mass(mol, opt=False): # for testing total mass
    mass = 0
    for i in range(0, mol.natom()):
        mass += mol.mass(i,opt)
    return mass

# mass defaults to zero_ghost = false
assert psi4.compare_values(m_He, A_frag.mass(1), 9, 'A_frag mass default')
assert psi4.compare_values(m_He, A_frag.mass(1, False), 9, 'A_frag mass False')
assert psi4.compare_values(m_Gh, A_frag.mass(1, True), 9, 'A_frag mass True')
assert psi4.compare_values(m_He, B_frag.mass(1), 9, 'B_frag mass default')
assert psi4.compare_values(m_He, B_frag.mass(1, False), 9, 'B_frag mass False')
assert psi4.compare_values(m_Gh, B_frag.mass(1, True), 9, 'B_frag mass True')

# z defaults to zero_ghost = true
assert psi4.compare_values(Z_Gh, A_frag.Z(1), 9, 'A_frag Z default')
assert psi4.compare_values(Z_Gh, A_frag.Z(1, True), 9, 'A_frag Z True')
assert psi4.compare_values(Z_He, A_frag.Z(1, False), 9, 'A_frag Z False')
assert psi4.compare_values(Z_Gh, B_frag.Z(1), 9, 'B_frag Z default')
assert psi4.compare_values(Z_Gh, B_frag.Z(1, True), 9, 'B_frag Z True')
assert psi4.compare_values(Z_He, B_frag.Z(1, False), 9, 'B_frag Z False')
