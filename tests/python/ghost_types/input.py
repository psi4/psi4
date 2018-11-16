#! Testing mass() and Z() for both ghost types: A (Gh(x) or @(x)) and B (extract_subsets())
import psi4

m_H = 1.00782503207      #TEST
m_He = 4.00260325415     #TEST
m_Gh = 0.0               #TEST
Z_H = 1.0                #TEST
Z_He = 2.0               #TEST
Z_Gh = 0.0               #TEST
NRE_real = 1.05835442134 #TEST
NRE_Gh = 0.0             #TEST

frag = psi4.geometry("""
    H
    --
    He 1 1
""")

A_frag = psi4.geometry(""" 
    H
    --
    Gh(He) 1 1
""")

B_frag = frag.extract_subsets(1,2)

mix_frag = A_frag.extract_subsets(2,1) # two ghost atoms: H is type B, He is type A

def mol_mass(mol,opt): # for testing total mass
    mass = 0
    for i in range(0, mol.natom()):
        mass += mol.mass(i,opt)
    return mass

# mass defaults to zero_ghost = false
# check mass of He for frag, A_frag, and B_frag
assert psi4.compare_values(m_He, A_frag.mass(1), 9, 'A_frag mass default')
assert psi4.compare_values(m_He, A_frag.mass(1, False), 9, 'A_frag mass False')
assert psi4.compare_values(m_Gh, A_frag.mass(1, True), 9, 'A_frag mass True')
assert psi4.compare_values(m_He, B_frag.mass(1), 9, 'B_frag mass default')
assert psi4.compare_values(m_He, B_frag.mass(1, False), 9, 'B_frag mass False')
assert psi4.compare_values(m_Gh, B_frag.mass(1, True), 9, 'B_frag mass True')

# check masses for mix_frag
assert psi4.compare_values(m_Gh, mix_frag.mass(0, True), 9, 'mix_frag H mass True')
assert psi4.compare_values(m_H, mix_frag.mass(0, False), 9, 'mix_frag H mass False')
assert psi4.compare_values(m_Gh, mix_frag.mass(1, True), 9, 'mix_frag He mass True')
assert psi4.compare_values(m_He, mix_frag.mass(1, False), 9, 'mix_frag He mass False')

# z defaults to zero_ghost = true
# check Z of He for frag, A_frag, and B_frag
assert psi4.compare_values(Z_Gh, A_frag.Z(1), 9, 'A_frag Z default')
assert psi4.compare_values(Z_Gh, A_frag.Z(1, True), 9, 'A_frag Z True')
assert psi4.compare_values(Z_He, A_frag.Z(1, False), 9, 'A_frag Z False')
assert psi4.compare_values(Z_Gh, B_frag.Z(1), 9, 'B_frag Z default')
assert psi4.compare_values(Z_Gh, B_frag.Z(1, True), 9, 'B_frag Z True')
assert psi4.compare_values(Z_He, B_frag.Z(1, False), 9, 'B_frag Z False')

# check Z for mix_frag
assert psi4.compare_values(Z_Gh, mix_frag.Z(0, True), 9, 'mix_frag H Z True')
assert psi4.compare_values(Z_H, mix_frag.Z(0, False), 9, 'mix_frag H Z False')
assert psi4.compare_values(Z_Gh, mix_frag.Z(1, True), 9, 'mix_frag He Z True')
assert psi4.compare_values(Z_He, mix_frag.Z(1, False), 9, 'mix_frag He Z False')

# check molecular mass of all frag types + mix_frag
# both atoms return non-zero for False
assert psi4.compare_values((m_H + m_He), mol_mass(frag, False), 9, 'frag mol_mass False')
assert psi4.compare_values((m_H + m_He), mol_mass(A_frag, False), 9, 'A_frag mol_mass False') 
assert psi4.compare_values((m_H + m_He), mol_mass(B_frag, False), 9, 'B_frag mol_mass False')
assert psi4.compare_values((m_H + m_He), mol_mass(mix_frag, False), 9, 'mix_frag mol_mass False')
# ghost returns zero mass for True
assert psi4.compare_values((m_H + m_He), mol_mass(frag, True), 9, 'frag mol_mass True')
assert psi4.compare_values(m_H, mol_mass(A_frag, True), 9, 'A_frag mol_mass True')
assert psi4.compare_values(m_H, mol_mass(B_frag, True), 9, 'B_frag mol_mass True')
assert psi4.compare_values(m_Gh, mol_mass(mix_frag, True), 9, 'mix_frag mol_mass True')

# check nuclear_repulsion_energy of all frag types
# any dimer NRE w/ ghost atom(s) is 0
assert psi4.compare_values(NRE_real, frag.nuclear_repulsion_energy(), 9, 'frag NRE')
assert psi4.compare_values(NRE_Gh, A_frag.nuclear_repulsion_energy(), 9, 'A_frag NRE')
assert psi4.compare_values(NRE_Gh, B_frag.nuclear_repulsion_energy(), 9, 'B_frag NRE')
assert psi4.compare_values(NRE_Gh, mix_frag.nuclear_repulsion_energy(), 9, 'mix_frag NRE')
