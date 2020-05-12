import pytest
#from qcengine.programs.tests.standard_suite_contracts import (
#    contractual_ccsd,
#    contractual_ccsd_prt_pr,
#    contractual_current,
#    contractual_mp2,
#    query_has_qcvar,
#    query_qcvar,
#)  # skip is temporary until references in place at qcng
#from qcengine.programs.tests.standard_suite_ref import answer_hash, std_suite  # skip is temporary until references in place at qcng

import psi4

from .utils import *


def runner_asserter(inp, subject, method, basis, tnm):

    qc_module = "-".join(["psi4", inp["keywords"].get("qc_module", "")]).strip("-")
    driver = inp["driver"]
    reference = inp["keywords"]["reference"]
    fcae = {"true": "fc", "false": "ae"}[inp["keywords"]["freeze_core"]]

    # <<<  Reference Values  >>>

    # ? precedence on next two
    mp2_type = inp.get("corl_type", inp["keywords"].get("mp2_type", "df"))  # hard-code of read_options.cc MP2_TYPE
    cc_type = inp.get("corl_type", inp["keywords"].get("cc_type", "conv"))  # hard-code of read_options.cc CC_TYPE
    corl_natural_values = {"mp2": mp2_type, "ccsd": cc_type, "ccsd(t)": cc_type}
    corl_type = corl_natural_values[method]

    natural_ref = {"conv": "pk", "df": "df", "cd": "cd"}
    scf_type = inp["keywords"].get("scf_type", natural_ref[corl_type])
    natural_values = {"pk": "pk", "direct": "pk", "df": "df", "mem_df": "df", "disk_df": "df", "cd": "cd"}
    scf_type = natural_values[scf_type]

    atol = 1.0e-6
    chash = answer_hash(
        system=subject.name(), basis=basis, fcae=fcae, scf_type=scf_type, reference=reference, corl_type=corl_type,
    )

    # check all calcs against conventional reference to looser tolerance
    atol_conv = 3.0e-4  # for df-ccsd. mp2 ok with 1.e-4
    chash_conv = answer_hash(
        system=subject.name(), basis=basis, fcae=fcae, reference=reference, corl_type="conv", scf_type="pk",
    )
    ref_block_conv = std_suite[chash_conv]

    # <<<  Prepare Calculation and Call API  >>>

    driver_call = {"energy": psi4.energy, "gradient": psi4.gradient}

    psi4.set_options(
        {
            # "guess": "sad",
            # "e_convergence": 8,
            # "d_convergence": 7,
            # "r_convergence": 7,

            # "e_convergence": 10,
            # "d_convergence": 9,
            # "r_convergence": 9,
            # "pcg_convergence": 9,
            "pcg_convergence": 7,
            "points": 5,
        }
    )
    extra_kwargs = inp["keywords"].pop("function_kwargs", {})
    psi4.set_options(inp["keywords"])

    if "error" in inp:
        errtype, errmsg = inp["error"]
        with pytest.raises(errtype) as e:
            driver_call[driver](inp["call"], molecule=subject, **extra_kwargs)

        assert errmsg in str(e.value)
        return

    ret, wfn = driver_call[driver](inp["call"], molecule=subject, return_wfn=True, **extra_kwargs)

    # <<<  Comparison Tests  >>>

    ref_block = std_suite[chash]

    # qcvars
    contractual_args = [
        qc_module,
        driver,
        reference,
        method,
        corl_type,
        fcae,
    ]
    asserter_args = [
        [psi4.core, wfn],
        ref_block,
        atol,
        ref_block_conv,
        atol_conv,
        tnm,
    ]

    def qcvar_assertions():
        print("BLOCK", chash, contractual_args)
        if method == "mp2":
            _asserter(asserter_args, contractual_args, contractual_mp2)
        elif method == "ccsd":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_ccsd)
        elif method == "ccsd(t)":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_ccsd)
            _asserter(asserter_args, contractual_args, contractual_ccsd_prt_pr)

    if "wrong" in inp:
        errmsg, reason = inp["wrong"]
        with pytest.raises(AssertionError) as e:
            qcvar_assertions()

        # print("WRONG", errmsg, reason, str(e.value), "ENDW")
        assert errmsg in str(e.value)
        pytest.xfail(reason)

    qcvar_assertions()

    # aliases
    _asserter(asserter_args, contractual_args, contractual_current)

    # returns
    tf, errmsg = compare_values(ref_block[f"{method.upper()} TOTAL ENERGY"], wfn.energy(), tnm + " wfn", atol=atol, return_message=True, quiet=True)
    assert compare_values(ref_block[f"{method.upper()} TOTAL ENERGY"], wfn.energy(), tnm + " wfn", atol=atol), errmsg

    if driver == "energy":
        assert compare_values(ref_block[f"{method.upper()} TOTAL ENERGY"], ret, tnm + " return")

    elif driver == "gradient":
        assert compare_values(
            ref_block[f"{method.upper()} TOTAL GRADIENT"], wfn.gradient().np, tnm + " grad wfn", atol=atol
        )
        assert compare_values(ref_block[f"{method.upper()} TOTAL GRADIENT"], ret.np, tnm + " grad return", atol=atol)


def _asserter(asserter_args, contractual_args, contractual_fn):
    """For expectations in `contractual_fn`, check that the QCVars are present in P::e.globals and wfn and match expected ref_block."""

    qcvar_stores, ref_block, atol, ref_block_conv, atol_conv, tnm = asserter_args

    for obj in qcvar_stores:
        for rpv, pv, present in contractual_fn(*contractual_args):
            label = tnm + " " + pv

            if present:
                # verify exact match to method (may be df) and near match to conventional (non-df) method
                tf, errmsg = compare_values(
                    ref_block[rpv], query_qcvar(obj, pv), label, atol=atol, return_message=True, quiet=True
                )
                assert compare_values(ref_block[rpv], query_qcvar(obj, pv), label, atol=atol), errmsg
                tf, errmsg = compare_values(
                    ref_block_conv[rpv], query_qcvar(obj, pv), label, atol=atol_conv, return_message=True, quiet=True
                )
                assert compare_values(ref_block_conv[rpv], query_qcvar(obj, pv), label, atol=atol_conv), errmsg

                # Note that the double compare_values lines are to collect the errmsg in the first for assertion in the second.
                #   If the errmsg isn't present in the assert, the string isn't accessible through `e.value`.
                #   If a plain bool is compared in the assert, the printed message will show booleans and not numbers.
            else:
                # verify and forgive known contract violations
                assert compare(False, query_has_qcvar(obj, pv), label + " SKIP")
