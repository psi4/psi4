import pytest

from psi4 import compare_values

import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.dft, pytest.mark.scf]


@pytest.mark.parametrize("db,refie", [
    (2, -5.388),
    (4, -16.352),
])
def test_xyg3_s22_ie(db, refie):
    """S22 unCP interaction energy at XYG3/aug-cc-pVDZ."""

    dh = 0.3211
    labels = ["§(auto)_(1, 2)@(1, 2)", "§(auto)_(1)@(1)", "§(auto)_(2)@(2)"]
    def ie(d, ma, mb):
        return (d - ma - mb) * psi4.constants.hartree2kcalmol

    ref = {
        2: {
            "DFT ORBITAL TOTAL ENERGY": [-152.896556997805, -76.444577153690, -76.444574136778],  # NWChem: Total DFT energy =
            "DFT FUNCTIONAL TOTAL ENERGY": [-152.584365930251, -76.288376390838, -76.288510494955],  # NWChem: DFT energy
            "DOUBLE-HYBRID CORRECTION ENERGY": [dh * -0.628008042867, dh * -0.312449744115, dh * -0.312114201506],  # NWChem: Unscaled MP2 energy
            "DFT TOTAL ENERGY": [-152.786019312815, -76.388704003673, -76.388730365058]  # NWChem: Total DFT+MP2 energy
        },
        4: {
            "DFT ORBITAL TOTAL ENERGY": [-339.872394738091, -169.924556335410, -169.924556335478],
            "DFT FUNCTIONAL TOTAL ENERGY": [-339.125711995348, -169.550887686232, -169.550887686235],
            "DOUBLE-HYBRID CORRECTION ENERGY": [dh * -1.542536206757, dh * -0.767970003064, dh * -0.767970003069],
            "DFT TOTAL ENERGY": [-339.621020371338, -169.797482854216, -169.797482854220],
        },
    }
    for item in list(ref.keys()):
        for pv in list(ref[item].keys()):
            ref[item][pv] = dict(zip(labels, ref[item][pv]))

    sdimer = {
        2: """
O  -1.551007  -0.114520   0.000000
H  -1.934259   0.762503   0.000000
H  -0.599677   0.040712   0.000000
--
0 1
O   1.350625   0.111469   0.000000
H   1.680398  -0.373741  -0.758561
H   1.680398  -0.373741   0.758561
        """,
        4: """
0 1
C  -2.018649   0.052883   0.000000
O  -1.452200   1.143634   0.000000
N  -1.407770  -1.142484   0.000000
H  -1.964596  -1.977036   0.000000
H  -0.387244  -1.207782   0.000000
H  -3.117061  -0.013701   0.000000
--
0 1
C   2.018649  -0.052883   0.000000
O   1.452200  -1.143634   0.000000
N   1.407770   1.142484   0.000000
H   1.964596   1.977036   0.000000
H   0.387244   1.207782   0.000000
H   3.117061   0.013701   0.000000
        """,
    }[db]
    dimer = psi4.geometry(sdimer)

    psi4.set_options(
        {
            "basis": "aug-cc-pvdz",
            #"puream": False,
            #"scf_type": "pk",
            "scf_type": "df",
            "mp2_type": "df",
            #"freeze_core": True,
            "dft_radial_points": 99,
            "dft_spherical_points": 590,
            "e_convergence": 8,
            "d_convergence": 8,
        }
    )

    ie_au, wfn = psi4.energy("xyg3", molecule=dimer, bsse_type="nocp", return_wfn=True)

    assert compare_values(refie, ie_au * psi4.constants.hartree2kcalmol, 2, f"S22-{db} IE")

    if db == 2:
        tottol, ietol = 0.0001, 0.001
    elif db == 4:
        tottol, ietol = 0.0005, 0.008

    import qcmanybody as qcmb
    cluster_qcvars = {qcmb.labeler(*qcmb.delabeler(k), opaque=False): v.extras["qcvars"] for k, v in wfn.qcschema.component_results.items()}
    for pv in ref[db].keys():
        for cluster in ref[db][pv].keys():
            assert compare_values(ref[db][pv][cluster], cluster_qcvars[cluster][pv], tottol, f"{cluster}: {pv}")
        assert compare_values(ie(ref[db][pv][labels[0]], ref[db][pv][labels[1]], ref[db][pv][labels[2]]),
                              ie(cluster_qcvars[labels[0]][pv], cluster_qcvars[labels[1]][pv], cluster_qcvars[labels[2]][pv]), ietol, "IE")


if __name__ == "__main__":

    test_xyg3_s22_ie(2, -5.388)
    #test_xyg3_s22_ie(4, -16.352)

    #test_xyg3_s22_ie(2, -5.48)  # -152.74379273     -76.36749396     -76.36756211"  IE= -5.48225  # paper with Alvaro
