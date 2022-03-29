import pytest
from addons import uusing

pytestmark = [pytest.mark.psi, pytest.mark.api]


@uusing("openfermionpsi4")
def test_openfermionpsi4_demo():
    import sys
    import math
    from openfermion.chem import MolecularData
    from openfermionpsi4 import run_psi4
    
    # Set molecule parameters.
    basis = 'sto-3g'
    multiplicity = 1
    bond_length_interval = 0.2
    n_points = 10
    
    # Set calculation parameters.
    run_scf = 1
    run_mp2 = 1
    run_cisd = 0
    run_ccsd = 0
    run_fci = 1
    delete_input = True
    delete_output = True
    
    # Generate molecule at different bond lengths.
    hf_energies = []
    fci_energies = []
    bond_lengths = []
    for point in range(1, n_points + 1):
        bond_length = bond_length_interval * float(point)
        bond_lengths += [bond_length]
        geometry = [('H', (0., 0., 0.)), ('H', (0., 0., bond_length))]
        molecule = MolecularData(
            geometry, basis, multiplicity,
            description=str(round(bond_length, 2)))
        
        # Run Psi4.
        molecule = run_psi4(molecule,
                            run_scf=run_scf,
                            run_mp2=run_mp2,
                            run_cisd=run_cisd,
                            run_ccsd=run_ccsd,
                            run_fci=run_fci)
    
        # Print out some results of calculation.
        print('\nAt bond length of {} angstrom, molecular hydrogen has:'.format(
            bond_length))
        print('Hartree-Fock energy of {} Hartree.'.format(molecule.hf_energy))
        print('MP2 energy of {} Hartree.'.format(molecule.mp2_energy))
        print('FCI energy of {} Hartree.'.format(molecule.fci_energy))
        print('Nuclear repulsion energy between protons is {} Hartree.'.format(
            molecule.nuclear_repulsion))
        for orbital in range(molecule.n_orbitals):
            print('Spatial orbital {} has energy of {} Hartree.'.format(
                orbital, molecule.orbital_energies[orbital]))
        hf_energies += [molecule.hf_energy]
        fci_energies += [molecule.fci_energy]
    
        if point == 5:
            # At bond length of 1.0 angstrom, molecular hydrogen has:
            assert math.isclose(molecule.hf_energy, -1.06610864808, abs_tol=1e-5)
            assert math.isclose(molecule.mp2_energy, -1.08666480357, abs_tol=1e-5)
            assert math.isclose(molecule.fci_energy, -1.1011503293, abs_tol=1e-5)
            assert math.isclose(molecule.nuclear_repulsion, 0.52917720859, abs_tol=1e-4)

