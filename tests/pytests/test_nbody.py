import re
import pytest
import numpy as np
from qcelemental.testing import compare, compare_values
import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.nbody]

_tot_cp_ene = -155.40761029
_ie_cp_ene = -0.00172611
_tot_uncp_ene = -155.40887161
_ie_uncp_ene = -0.00298743

_tot_cp_grad = np.array(
    [[-0.            ,  0.003989481299,  0.000257671876],
     [ 0.            , -0.003989481299,  0.000257671876],
     [-0.007206301768,  0.004333142371,  0.000029433226],
     [ 0.007206301768,  0.004333142371,  0.000029433226],
     [ 0.007206301768, -0.004333142371,  0.000029433226],
     [-0.007206301768, -0.004333142371,  0.000029433226],
     [-0.            , -0.            , -0.035555548002],
     [-0.            ,  0.            ,  0.035547129205],
     [ 0.            , -0.            ,  0.00840988492 ],
     [ 0.            ,  0.            , -0.009034542779]])
_ie_cp_grad = np.array(
    [[-0.           ,   0.000831523566,  0.000095334214],
     [ 0.           ,  -0.000831523566,  0.000095334214],
     [-0.00011127085,   0.000096261386,  0.000110602057],
     [ 0.00011127085,   0.000096261386,  0.000110602057],
     [ 0.00011127085,  -0.000096261386,  0.000110602057],
     [-0.00011127085,  -0.000096261386,  0.000110602057],
     [-0.           ,   0.            , -0.00139615575 ],
     [-0.           ,   0.            , -0.000684195246],
     [ 0.           ,   0.            ,  0.001393980177],
     [ 0.           ,   0.            ,  0.000053294161]])
_tot_uncp_grad = np.array(
    [[-0.            ,  0.004688102439,  0.000019589414],
     [ 0.            , -0.004688102439,  0.000019589414],
     [-0.007232900028,  0.004303778033,  0.000005580357],
     [ 0.007232900028,  0.004303778033,  0.000005580357],
     [ 0.007232900028, -0.004303778033,  0.000005580357],
     [-0.007232900028, -0.004303778033,  0.000005580357],
     [-0.            ,  0.            , -0.035419642667],
     [-0.            ,  0.            ,  0.03585801241 ],
     [ 0.            ,  0.            ,  0.008543047852],
     [ 0.            ,  0.            , -0.009042917851]])
_ie_uncp_grad = np.array(
    [[-0.           ,   0.001530144706, -0.000142748248],
     [ 0.           ,  -0.001530144706, -0.000142748248],
     [-0.00013786911,   0.000066897048,  0.000086749188],
     [ 0.00013786911,   0.000066897048,  0.000086749188],
     [ 0.00013786911,  -0.000066897048,  0.000086749188],
     [-0.00013786911,  -0.000066897048,  0.000086749188],
     [-0.           ,   0.            , -0.001260250416],
     [-0.           ,  -0.            , -0.000373312041],
     [ 0.           ,   0.            ,  0.00152714311 ],
     [ 0.           ,   0.            ,  0.000044919089]])


_stdouts = {
    "cp_T": r"""
\s*   ==> N-Body: Counterpoise Corrected \(CP\) energies <==
\s*        n-Body     Total Energy            Interaction Energy                          N-body Contribution to Interaction Energy
\s*                   \[Eh\]                    \[Eh\]                  \[kcal/mol\]            \[Eh\]                  \[kcal/mol\]
\s*             1     -155.4058841\d\d\d\d\d        0.000000000000        0.000000000000        0.000000000000        0.000000000000
\s*  FULL/RTN   2     -155.4076102\d\d\d\d\d       -0.0017261\d\d\d\d\d       -1.08315\d\d\d\d\d\d\d       -0.0017261\d\d\d\d\d       -1.08315\d\d\d\d\d\d\d
""",
    "cp_F": r"""
\s*   ==> N-Body: Counterpoise Corrected \(CP\) energies <==
\s*        n-Body     Total Energy            Interaction Energy                          N-body Contribution to Interaction Energy
\s*                   \[Eh\]                    \[Eh\]                  \[kcal/mol\]            \[Eh\]                  \[kcal/mol\]
\s*             1        N/A                   0.000000000000        0.000000000000        0.000000000000        0.000000000000
\s*  FULL/RTN   2        N/A                  -0.0017261\d\d\d\d\d       -1.08315\d\d\d\d\d\d\d       -0.0017261\d\d\d\d\d       -1.08315\d\d\d\d\d\d\d
""",
    "uncp": r"""
\s*   ==> N-Body: Non-Counterpoise Corrected \(NoCP\) energies <==
\s*        n-Body     Total Energy            Interaction Energy                          N-body Contribution to Interaction Energy
\s*                   \[Eh\]                    \[Eh\]                  \[kcal/mol\]            \[Eh\]                  \[kcal/mol\]
\s*             1     -155.4058841\d\d\d\d\d        0.000000000000        0.000000000000        0.000000000000        0.000000000000
\s*  FULL/RTN   2     -155.4088716\d\d\d\d\d       -0.0029874\d\d\d\d\d       -1.87464\d\d\d\d\d\d\d       -0.0029874\d\d\d\d\d       -1.87464\d\d\d\d\d\d\d
""",
}
_stdouts["cpuncp"] = _stdouts["cp_T"] + _stdouts["uncp"]
_stdouts["uncpcp"] = _stdouts["uncp"] + _stdouts["cp_T"]


@pytest.mark.parametrize("driver,bsse_type,return_total_data,nbody_number,return_result,stdoutkey", [
    ("energy",   ["cp"],         True , 5, _tot_cp_ene,   "cp_T"),     # return CP tot         5
    ("energy",   ["cp"],         False, 3, _ie_cp_ene,    "cp_F"),     # return CP IE          3
    ("energy",   ["cp"],         None , 3, _ie_cp_ene,    "cp_F"),     # return CP IE          3
    ("energy",   ["nocp"],       True , 3, _tot_uncp_ene, "uncp"),     # return tot            3
    ("energy",   ["nocp"],       False, 3, _ie_uncp_ene,  "uncp"),     # return IE             3
    ("energy",   ["nocp"],       None , 3, _ie_uncp_ene,  "uncp"),     # return IE             3
    ("energy",   ["cp", "nocp"], True , 5, _tot_cp_ene,   "cpuncp"),   # return CP tot         5
    ("energy",   ["cp", "nocp"], False, 5, _ie_cp_ene,    "cpuncp"),   # return CP IE          5
    ("energy",   ["cp", "nocp"], None , 5, _ie_cp_ene,    "cpuncp"),   # return CP IE          5
    ("energy",   ["nocp", "cp"], True , 5, _tot_uncp_ene, "uncpcp"),   # return tot            5
    ("energy",   ["nocp", "cp"], False, 5, _ie_uncp_ene,  "uncpcp"),   # return IE             5
    ("energy",   ["nocp", "cp"], None , 5, _ie_uncp_ene,  "uncpcp"),   # return IE             5
    ("energy",   ["ssfc"],       None , 3, _ie_cp_ene,    "cp_F"),     # return CP IE          3
    ("gradient", ["cp"],         True , 5, _tot_cp_grad,   "cp_T"),    # return CP tot G       5
    ("gradient", ["cp"],         False, 3, _ie_cp_grad,    "cp_F"),    # return CP IE G        3
    ("gradient", ["cp"],         None , 5, _tot_cp_grad,   "cp_T"),    # return CP tot G       5
    ("gradient", ["nocp"],       True , 3, _tot_uncp_grad, "uncp"),    # return tot G          3
    ("gradient", ["nocp"],       False, 3, _ie_uncp_grad,  "uncp"),    # return IE G           3
    ("gradient", ["nocp"],       None , 3, _tot_uncp_grad, "uncp"),    # return tot G          3
    ("gradient", ["cp", "nocp"], True , 5, _tot_cp_grad,   "cpuncp"),  # return CP tot G       5
    ("gradient", ["cp", "nocp"], False, 5, _ie_cp_grad,    "cpuncp"),  # return CP IE G        5
    ("gradient", ["cp", "nocp"], None , 5, _tot_cp_grad,   "cpuncp"),  # return CP tot G       5
    ("gradient", ["nocp", "cp"], True , 5, _tot_uncp_grad, "uncpcp"),  # return tot G          5
    ("gradient", ["nocp", "cp"], False, 5, _ie_uncp_grad,  "uncpcp"),  # return IE G           5
    ("gradient", ["nocp", "cp"], None , 5, _tot_uncp_grad, "uncpcp"),  # return tot G          5
    ("gradient", ["nocp", "ssfc"], None , 5, _tot_uncp_grad, "uncpcp"),  # return tot G          5
])
def test_nbody_number(driver, bsse_type, return_total_data, nbody_number, return_result, stdoutkey):

    eneyne = psi4.geometry("""
C   0.000000  -0.667578  -2.124659
C   0.000000   0.667578  -2.124659
H   0.923621  -1.232253  -2.126185
H  -0.923621  -1.232253  -2.126185
H  -0.923621   1.232253  -2.126185
H   0.923621   1.232253  -2.126185
--
C   0.000000   0.000000   2.900503
C   0.000000   0.000000   1.693240
H   0.000000   0.000000   0.627352
H   0.000000   0.000000   3.963929
""")

    atin = {
        "driver": driver,
        "model": {
            "method": "mp2",
            "basis": "cc-pvdz",
        },
        "molecule": eneyne.to_schema(dtype=2),
        "keywords": {
            "function_kwargs": {
                "bsse_type": bsse_type,
                "return_total_data": return_total_data,
            },
        },
    }

    ret = psi4.schema_wrapper.run_qcschema(atin)

    assert compare_values(return_result, ret.return_result, atol=1.e-6, label="manybody")
    assert compare(nbody_number, ret.extras["qcvars"]["NBODY NUMBER"], label="nbody number")
    assert re.search(_stdouts[stdoutkey], ret.stdout, re.MULTILINE), f"N-Body pattern not found: {_stdouts[stdoutkey]}"

