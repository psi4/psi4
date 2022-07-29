import numpy as np
from qcengine.programs.tests.standard_suite_ref import answer_hash, compute_derived_qcvars, _std_suite, _std_generics


# in-repo extensions for _std_suite above
# * ideally empty. PR to QCEngine ASAP and empty this after QCEngine release.
_std_suite_psi4_extension = [
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
