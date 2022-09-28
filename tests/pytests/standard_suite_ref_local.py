import numpy as np
from qcengine.programs.tests.standard_suite_ref import answer_hash, compute_derived_qcvars, _std_suite, _std_generics


# in-repo extensions for _std_suite above
# * ideally empty. PR to QCEngine ASAP and empty this after QCEngine release.
_std_suite_psi4_extension = [
    # <<<  CONV-AE-CONV  >>>
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
            "CISD CORRELATION ENERGY": -0.08142433,  # detci != vcc ???  # locally, replacing the rohf cisd vcc=tce value by the detci=guga value. correct sdsc label unclear.
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
