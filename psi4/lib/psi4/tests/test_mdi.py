import pytest
from utils import *
from addons import uusing

from qcelemental.testing import compare_values

import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick, pytest.mark.mdi]

@pytest.mark.smoke
@uusing("mdi")
def test_mdi_water():

    psi4.geometry("""
    0 1
    O
    H 1 1.0
    H 1 1.0 2 104.5
    """)

    psi4.set_options({
        'reference': 'uhf',
        'scf_type': 'direct',
        'e_convergence': 1e-12,
    })

    scf_method = "scf/cc-pVDZ"

    psi4.mdi_engine.mdi_init("-role DRIVER -name Psi4 -method TEST")
    engine = psi4.mdi_engine.MDIEngine(scf_method)

    # Test the <NATOMS command
    natom = engine.send_natoms()
    assert natom == 3

    # Test the <COORDS command
    coords = engine.send_coords()
    expected = [
        0.0, 0.0, -0.12947688966627896, 0.0, -1.4941867446308548, 1.027446098871551, 0.0, 1.4941867446308548,
        1.027446098871551
    ]
    assert compare_values(expected, coords, atol=1.e-7)

    # Test the >COORDS command
    expected = [
        0.1, 0.0, -0.12947688966627896, 0.0, -1.4341867446308548, 1.027446098871551, 0.0, 1.4941867446308548,
        1.027446098871551
    ]
    engine.recv_coords(expected)
    coords = engine.send_coords()
    assert compare_values(expected, coords, atol=1.e-7)

    # Test the <CHARGES command
    charges = engine.send_charges()
    expected = [8.0, 1.0, 1.0]
    assert compare_values(expected, charges, atol=1.e-7)

    # Test the <ELEMENTS command
    elements = engine.send_elements()
    expected = [8.0, 1.0, 1.0]
    assert compare_values(expected, elements, atol=1.e-7)

    # Test the <MASSES command
    masses = engine.send_masses()
    expected = [15.99491461957, 1.00782503223, 1.00782503223]
    assert compare_values(expected, masses, atol=1.e-6)

    # Test the SCF command
    engine.run_scf()

    # Test the <ENERGY command
    energy = engine.send_energy()
    expected = -76.02320201768676
    assert compare_values(expected, energy, atol=1.e-6)

    # Test the <FORCES command
    forces = engine.send_forces()
    expected = [
        0.004473733542292732, -0.01775852359379196, -0.051757651796320636, -0.0016762687719661835,
        -0.024093270269019917, 0.019393138772564877, -0.0027974647702799747, 0.04185179386282589, 0.03236451302360538
    ]
    assert compare_values(expected, forces, atol=1.e-6)

    # Test the >MASSES command
    expected = [15.99491461957, 3.0, 2.0]
    engine.recv_masses(expected)
    masses = engine.send_masses()
    assert compare_values(expected, masses, atol=1.e-7)

    # Test the <DIMENSIONS command
    expected = [1, 1, 1]
    dimensions = engine.send_dimensions()
    assert compare_values(expected, dimensions, atol=1.e-7)

    # Test the <TOTCHARGE command
    totcharge = engine.send_total_charge()
    expected = 0.0
    assert compare_values(expected, totcharge, atol=1.e-7)

    # Test the >TOTCHARGE command
    expected = 1.0
    engine.recv_total_charge(expected)
    totcharge = engine.send_total_charge()
    assert compare_values(expected, totcharge, atol=1.e-7)

    # Test the <ELEC_MULT command
    multiplicity = engine.send_multiplicity()
    expected = 1
    assert compare_values(expected, multiplicity, atol=1.e-7)

    # Test the >ELEC_MULT command
    expected = 2
    engine.recv_multiplicity(expected)
    multiplicity = engine.send_multiplicity()
    assert compare_values(expected, multiplicity, atol=1.e-7)

    # Test the >NLATTICE, >CLATTICE, and >LATTICE commands
    nlattice = 2
    clattice = [-3.0, 1.0, 0.0, 4.0, 3.0, 0.0]
    lattice = [-1.0, -0.5]
    engine.recv_nlattice(nlattice)
    engine.recv_clattice(clattice)
    engine.recv_lattice(lattice)

    # Test the final energy
    engine.run_scf()
    energy = engine.send_energy()
    expected = -76.03284774075954
    assert compare_values(expected, energy, atol=1.e-6)
    forces = engine.send_forces()
    expected = [
        0.032656171616, -0.027255620629, -0.02211206073, 0.015827924001, -0.009687198705,
        0.011685024025, 0.026352703866, 0.008313469614, 0.033873286169,
    ]
    assert compare_values(expected, forces, atol=1.e-6)

    # Test the EXIT command
    engine.exit()

    return None
