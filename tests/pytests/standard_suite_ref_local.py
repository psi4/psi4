from typing import Tuple

import numpy as np
from qcengine.programs.tests.standard_suite_ref import answer_hash, compute_derived_qcvars, _std_suite, _std_generics


# in-repo extensions for _std_suite above
# * ideally empty. PR to QCEngine ASAP and empty this after QCEngine release.
_std_suite_psi4_extension = [
    # <<<  CONV-AE-CONV  >>>
    {
        "meta": {
            "system": "hf",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "rhf",
            "fcae": "ae",
            "corl_type": "conv",
            "sdsc": "sd"
        },
        "data": {
        },
    },
    {
        "meta": {
            "system": "h2o",
            "basis": "aug-cc-pvdz",
            "scf_type": "pk",
            "reference": "rhf",
            "fcae": "ae",
            "corl_type": "conv",
            "sdsc": "sd"
        },
        "data": {
            "HF-CABS TOTAL ENERGY": -76.049013663448,
            "MP2-F12 CORRELATION ENERGY": -0.31082166907,
            "MP2-F12 SINGLES ENERGY": 0.0,
            "MP2-F12 SAME-SPIN CORRELATION ENERGY": -0.07400595,
        },
    },
    {
        "meta": {
            "system": "h2o",
            "basis": "cfour-qz2p",
            "scf_type": "pk",
            "reference": "rhf",
            "fcae": "ae",
            "corl_type": "conv",
            "sdsc": "sd"
        },
        "data": {
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "uhf",
            "fcae": "ae",
            "corl_type": "conv",
            "sdsc": "sd",
        },
        "data": {
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "aug-cc-pvdz",
            "scf_type": "pk",
            "reference": "uhf",
            "fcae": "ae",
            "corl_type": "conv",
            "sdsc": "sd"
        },
        "data": {
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "cfour-qz2p",
            "scf_type": "pk",
            "reference": "uhf",
            "fcae": "ae",
            "corl_type": "conv",
            "sdsc": "sd"
        },
        "data": {
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "rohf",
            "fcae": "ae",
            "corl_type": "conv",
            "sdsc": "sd",
        },
        "data": {
            "CISD CORRELATION ENERGY": -0.08142433,  # detci != cfour's vcc ???  # locally, replacing the rohf cisd vcc=tce value (stored in qcng) by the detci=guga value. correct sdsc label unclear.
            "FCI CORRELATION ENERGY": -0.084637876308811,  # detci
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "aug-cc-pvdz",
            "scf_type": "pk",
            "reference": "rohf",
            "fcae": "ae",
            "corl_type": "conv",
            "sdsc": "sd",
        },
        "data": {
            "CISD CORRELATION ENERGY": -0.1723668643052676,  # detci != vcc ???
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "cfour-qz2p",
            "scf_type": "pk",
            "reference": "rohf",
            "fcae": "ae",
            "corl_type": "conv",
            "sdsc": "sd",
        },
        "data": {
            "CISD CORRELATION ENERGY": -0.21038651,  # detci != vcc ???
        },
    },

    # <<<  CONV-FC-CONV  >>>
    {
        "meta": {
            "system": "h2o",
            "basis": "aug-cc-pvdz",
            "scf_type": "pk",
            "reference": "rhf",
            "fcae": "fc",
            "corl_type": "conv",
            "sdsc": "sd"
        },
        "data": {
            # molpro adz, default fitting. p4 agrees ~e-3 w/ cabs_basis=cc-pvdz-jkfit
            # "HF-CABS TOTAL ENERGY": -76.062354120090,
            # "MP2-F12 CORRELATION ENERGY": -0.292865859255,
            # "MP2-F12 SAME-SPIN CORRELATION ENERGY": -0.067098997832,
            "HF-CABS TOTAL ENERGY": -76.049013663448,
            "MP2-F12 CORRELATION ENERGY": -0.292676554889,
            "MP2-F12 SINGLES ENERGY": 0.0,
            "MP2-F12 SAME-SPIN CORRELATION ENERGY": -0.07262946,
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "rohf",
            "fcae": "fc",
            "corl_type": "conv",
            "sdsc": "sd",
        },
        "data": {
            "CISD CORRELATION ENERGY": -0.08045048714872,  # detci only != vcc ???
            "FCI CORRELATION ENERGY": -0.083612606639434,  # detci
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "aug-cc-pvdz",
            "scf_type": "pk",
            "reference": "rohf",
            "fcae": "fc",
            "corl_type": "conv",
            "sdsc": "sd",
        },
        "data": {
            "CISD CORRELATION ENERGY": -0.170209639586457,  # detci only != vcc ???
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "cfour-qz2p",
            "scf_type": "pk",
            "reference": "rohf",
            "fcae": "fc",
            "corl_type": "conv",
            "sdsc": "sd",
        },
        "data": {
            "CISD CORRELATION ENERGY": -0.186640254417867,  # detci only != vcc ???
        },
    },

    # <<<  DF-AE-DF  >>>
    {
        "meta": {
            "system": "hf",
            "basis": "cc-pvdz",
            "scf_type": "df",
            "reference": "rhf",
            "fcae": "ae",
            "corl_type": "df",
            "sdsc": "sd"
        },
        "data": {
        },
    },
    {
        "meta": {
            "system": "h2o",
            "basis": "aug-cc-pvdz",
            "scf_type": "df",
            "reference": "rhf",
            "fcae": "ae",
            "corl_type": "df",
            "sdsc": "sd"
        },
        "data": {
            "HF-CABS TOTAL ENERGY": -76.05009136215748,
            "MP2-F12 CORRELATION ENERGY": -0.310633975041,
            "MP2-F12 SINGLES ENERGY": 0.0,
            "MP2-F12 SAME-SPIN CORRELATION ENERGY": -0.07403290,
        },
    },
    {
        "meta": {
            "system": "h2o",
            "basis": "cfour-qz2p",
            "scf_type": "df",
            "reference": "rhf",
            "fcae": "ae",
            "corl_type": "df",
            "sdsc": "sd"
        },
        "data": {
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "df",
            "reference": "uhf",
            "fcae": "ae",
            "corl_type": "df",
            "sdsc": "sd",
        },
        "data": {
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "aug-cc-pvdz",
            "scf_type": "df",
            "reference": "uhf",
            "fcae": "ae",
            "corl_type": "df",
            "sdsc": "sd"
        },
        "data": {
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "cfour-qz2p",
            "scf_type": "df",
            "reference": "uhf",
            "fcae": "ae",
            "corl_type": "df",
            "sdsc": "sd"
        },
        "data": {
        },
    },
    # <<<  DF-FC-DF  >>>
    {
        "meta": {
            "system": "hf",
            "basis": "cc-pvdz",
            "scf_type": "df",
            "reference": "rhf",
            "fcae": "fc",
            "corl_type": "df",
            "sdsc": "sd"
        },
        "data": {
        },
    },
    {
        "meta": {
            "system": "h2o",
            "basis": "aug-cc-pvdz",
            "scf_type": "df",
            "reference": "rhf",
            "fcae": "fc",
            "corl_type": "df",
            "sdsc": "sd"
        },
        "data": {
            # molpro adz, default fitting. p4 agrees ~e-3
            # "HF-CABS TOTAL ENERGY": -76.062354120090,
            # "MP2-F12 CORRELATION ENERGY": -0.292858198131,
            # "MP-F12 SAME-SPIN CORRELATION ENERGY": -0.067147423336667,
            # molpro https://github.com/cdsgroup/qcdb/tree/data/data/f12dilabiousemefiles diff H2O geom
            # adz raw        -76.04117083
            # adz f12 raw    -76.06210955
            # adz mp2 corl    -0.21962112
            # adz mp2 trip    -0.05593210
            # adz mp2f12 corl -0.29304999
            # adz mp2f12 trip -0.06719202
            "HF-CABS TOTAL ENERGY": -76.05009136215748,
            "MP2-F12 CORRELATION ENERGY": -0.292535867049,
            "MP2-F12 SINGLES ENERGY": 0.0,
            "MP2-F12 SAME-SPIN CORRELATION ENERGY": -0.07265948,
        },
    },
    {
        "meta": {
            "system": "h2o",
            "basis": "cfour-qz2p",
            "scf_type": "df",
            "reference": "rhf",
            "fcae": "fc",
            "corl_type": "df",
            "sdsc": "sd"
        },
        "data": {
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "df",
            "reference": "uhf",
            "fcae": "fc",
            "corl_type": "df",
            "sdsc": "sd",
        },
        "data": {
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "aug-cc-pvdz",
            "scf_type": "df",
            "reference": "uhf",
            "fcae": "fc",
            "corl_type": "df",
            "sdsc": "sd"
        },
        "data": {
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "cfour-qz2p",
            "scf_type": "df",
            "reference": "uhf",
            "fcae": "fc",
            "corl_type": "df",
            "sdsc": "sd"
        },
        "data": {
        },
    },
    # <<<  CD-AE-CD  >>>
    # <<<  CD-FC-CD  >>>
]


for calc1 in _std_suite_psi4_extension:
    metahash1 = answer_hash(**calc1["meta"])
    for calc0 in _std_suite:
        metahash0 = answer_hash(**calc0["meta"])
        if metahash0 == metahash1:
            calc0["data"].update(calc1["data"])
            break

compute_derived_qcvars(_std_suite)
std_suite = {answer_hash(**calc["meta"]): calc["data"] for calc in _std_suite}


def contractual_mp2f12(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    f"""Of the list of QCVariables an ideal MP2-F12 should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {{_contractual_docstring}}
    """
    contractual_qcvars = [
        "HF-CABS TOTAL ENERGY",
        "MP2-F12 CORRELATION ENERGY",
        "MP2-F12 TOTAL ENERGY",
        "MP2-F12 SAME-SPIN CORRELATION ENERGY",
        "MP2-F12 SINGLES ENERGY",
        "MP2-F12 DOUBLES ENERGY",
        "MP2-F12 OPPOSITE-SPIN CORRELATION ENERGY",
    ]
    if driver == "gradient" and method == "mp2-f12":
        contractual_qcvars.append("MP2-F12 TOTAL GRADIENT")
    elif driver == "hessian" and method == "mp2-f12":
        # contractual_qcvars.append("MP2-F12 TOTAL GRADIENT")
        contractual_qcvars.append("MP2-F12 TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True
        if (
        ):
            expected = False

        yield (pv, pv, expected)
