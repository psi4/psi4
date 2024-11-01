import pytest
import psi4
import numpy as np

pytestmark = [pytest.mark.psi, pytest.mark.api]


def test_pdms():
    h2o = psi4.geometry("""
        H
        H 1 0.7
    """)

    psi4.set_options({"opdm": True, "num_roots": 3})
    wfn = psi4.energy("fci/cc-pvdz", return_wfn=True)[1]
    # Okay, I need three matrix comparisons...

    ref_opdm = psi4.core.Matrix.from_array([np.array([[ 9.84317436e-01,  4.69564137e-03,  3.68310242e-03],
        [ 4.69564137e-03,  2.71595295e-03, -8.89586171e-04],
        [ 3.68310242e-03, -8.89586171e-04,  4.24202798e-04]]),
        np.zeros((0,0)),
        np.array([[7.6893825e-05]]),
        np.array([[7.6893825e-05]]),
        np.zeros((0,0)),
        np.array([[0.00416892, 0.00440128, 0.00073837],
        [0.00440128, 0.00473526, 0.00083457],
        [0.00073837, 0.00083457, 0.00017527]]),
        np.array([[0.00165458]]),
        np.array([[0.00165458]])
        ])

    # Test OPDM
    assert psi4.compare_matrices(ref_opdm, wfn.get_opdm(-1, -1, "A", True), 6, "OPDM")

    ref_tdm = psi4.core.Matrix.from_array([np.array([[ 1.11543287e-02,  6.66837005e-02, -9.01185592e-03],
        [-3.95119873e-02,  2.57584925e-02,  4.69185262e-02],
        [-1.44471817e-03,  7.47375607e-05,  4.47155049e-04]]),
        np.zeros((0,0)),
        np.array([[-1.17823927e-05]]),
        np.array([[-1.17823927e-05]]),
        np.zeros((0,0)),
        np.array([[-0.0436292 , -0.01418167, -0.00124078],
        [ 0.04298612,  0.00671799,  0.00043999],
        [ 0.01168194,  0.0019797 ,  0.0001121 ]]),
        np.array([[-0.00026865]]),
        np.array([[-0.00026865]])
        ])

    # Test transition matrix
    assert psi4.compare_matrices(ref_tdm, wfn.get_opdm(1, 2, "A", True), 6, "TDM", equal_phase=True)
    # Test the swapping the bra and the key produces the transpose matrix.
    assert psi4.compare_matrices(wfn.get_opdm(2, 1, "A", True), wfn.get_opdm(1, 2, "A", True).transpose(), 6, "TDM SWAP")

