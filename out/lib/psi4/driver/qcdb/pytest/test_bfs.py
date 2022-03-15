from utils import *

import qcelemental as qcel

import qcdb


def test_Molecule_BFS():
    #! apply linear fragmentation algorithm to a water cluster

    iceIh = """\
    36
    Crystal created from CIF file. Box size:    7.82000    7.82000    7.36000
    O            2.606641   0.000000   3.220000
    O            5.213320   4.514902   6.900000
    O            6.516680   2.257417   3.220000
    N            2.606680   4.514902   3.220000
    O            1.303320   2.257417   6.900000
    O            5.213359   0.000000   4.140000
    O            6.516680   2.257417   0.460000
    O            5.213320   4.514902   4.140000
    O            1.303320   2.257417   4.140000
    O            2.606680   4.514902   0.460000
    O            5.213359   0.000000   6.900000
    O            2.606641   0.000000   0.460000
    H            2.193510   0.711093   3.496000
    H            7.339070   2.255182   3.496000
    H            5.622580   5.228230   7.176000
    H            1.712580   1.544089   7.176000
    He            3.429070   4.517137   3.496000
    H            6.103510   6.061225   3.496000
    H            1.716490   6.061225   7.176000
    H            4.390930   4.517137   7.176000
    H            6.107420   1.544089   3.496000
    He            2.197420   5.228230   3.496000
    H            0.480930   2.255182   7.176000
    H            4.394840   0.000000   3.871360
    H            6.107420   2.966276   0.191360
    H            5.622580   3.806043   3.871360
    H            1.712580   2.966276   3.871360
    H            2.197420   3.806043   0.191360
    H            5.213359   0.000000   4.960640
    H            6.516680   2.257417   1.280640
    H            5.213320   4.514902   4.960640
    H            1.303320   2.257417   4.960640
    H            2.606680   4.514902   1.280640
    H            5.626490   0.711093   7.176000
    H            3.425160   0.000000   0.191360
    H            2.606641   0.000000   1.280640
    """

    ref_fragmentation = [
      [3, 16],
      [21],
      [0, 12],
      [1, 14, 19],
      [2, 13, 20],
      [4, 15, 22],
      [5, 23, 28],
      [6, 24, 29],
      [7, 25, 30],
      [8, 26, 31],
      [9, 27, 32],
      [10, 33],
      [11, 34, 35],
      [17],
      [18]]  # yapf: disable

    qmol = qcdb.Molecule.from_string(iceIh, dtype='xyz')
    frag, arrs, bmols, bmol = qmol.BFS(
        seed_atoms=[[3, 16], [21]], return_arrays=True, return_molecule=True, return_molecules=True)

    assert compare_integers(frag == ref_fragmentation, 1, 'Q: BFS from qcdb.Molecule')
    assert compare_arrays(qmol.geometry(np_out=True)[[1, 14, 19]], arrs[0][3], 4, 'Q: geom back from BFS')
    assert compare_integers(15, bmol.nfragments(), 'Q: nfrag')
    assert compare_values(qmol.nuclear_repulsion_energy(), bmol.nuclear_repulsion_energy(), 4, 'Q: nre')
    assert compare_arrays(
        qmol.geometry(np_out=True)[[2, 13, 20]], bmols[4].geometry(np_out=True), 4, 'Q: frag geom back from BFS')
    assert compare_integers(True, type(bmol) == qcdb.Molecule, 'Q return type')


def test_numpy_BFS():
    import numpy as np
    from qcdb.bfs import BFS

    # FaOOFaOO 3.6 ?
    mol_elem = np.asarray(['C', 'C', 'H', 'H', 'O', 'O', 'O', 'O', 'H', 'H'])
    mol_geom = np.asarray([
        [    1.79035823,      -0.18606050,       0.00000000],
        [   -1.79035823,       0.18606050,       0.00000000],
        [    2.89087214,      -0.30042988,       0.00000000],
        [   -2.89087214,       0.30042988,       0.00000000],
        [    1.07568931,      -1.19425943,       0.00000000],
        [   -1.07568931,       1.19425943,       0.00000000],
        [    1.44185816,       1.08049605,       0.00000000],
        [   -1.44185816,      -1.08049605,       0.00000000],
        [    0.43274661,       1.15045330,       0.00000000],
        [   -0.43274661,      -1.15045330,       0.00000000]])  # yapf: disable


    ans = BFS(mol_geom / qcel.constants.bohr2angstroms, mol_elem)

    ref_fragmentation = [[0, 2, 4, 6, 8], [1, 3, 5, 7, 9]]
    assert compare_integers(True, ans == ref_fragmentation, 'BFS from np.array')
