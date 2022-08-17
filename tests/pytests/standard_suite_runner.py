import re
from typing import Callable, List

import pytest
from qcengine.programs.tests.standard_suite_contracts import *
from qcengine.programs.tests.standard_suite_ref import answer_hash
from standard_suite_ref_local import std_suite

import psi4

from utils import *


def runner_asserter(inp, subject, method, basis, tnm):

    qc_module_in = "-".join(["psi4", inp["keywords"].get("qc_module", "")]).strip("-")  # returns "psi4"|"psi4-<module>"  # input-specified routing
    qc_module_xptd = (
        ("psi4-" + inp["xptd"]["qc_module"]) if inp.get("xptd", {}).get("qc_module", None) else None
    )  # expected routing
    driver = inp["driver"]
    reference = inp["keywords"]["reference"]
    fcae = {"true": "fc", "false": "ae"}[inp["keywords"]["freeze_core"]]

    if qc_module_in == "psi4-detci" and basis != "cc-pvdz":
        pytest.skip(f"basis {basis} too big for {qc_module_in}")
    if qc_module_in == "psi4-ccenergy" and basis != "cc-pvdz" and method == "ccsd(t)" and reference == "uhf" and driver == "gradient" and inp["keywords"]["function_kwargs"]["dertype"] == 1:
        pytest.skip(f"ccenergy uhf analytic gradients add 10m")

    # <<<  Reference Values  >>>

    # ? precedence on next two
    scf_type = inp.get("corl_type", inp["keywords"].get("scf_type", "df"))  # hard-code of read_options.cc SCF_TYPE
    mp2_type = inp.get("corl_type", inp["keywords"].get("mp2_type", "df"))  # hard-code of read_options.cc MP2_TYPE
    if method in ["mp2.5", "mp3"]:
        mp_type = inp.get("corl_type", inp["keywords"].get("mp_type", "df"))  # hard-code of proc.py run_dfocc MP_TYPE
    else:
        mp_type = inp.get("corl_type", inp["keywords"].get("mp_type", "conv"))  # hard-code of read_options.cc MP_TYPE
    cc_type = inp.get("corl_type", inp["keywords"].get("cc_type", "conv"))  # hard-code of read_options.cc CC_TYPE
    corl_natural_values = {
        "hf": "df",  # dummy to assure df/cd/conv scf_type refs available
        "mp2": mp2_type,
        "mp2.5": mp_type,
        "mp3": mp_type,
        "remp2": cc_type,
        "lccd": cc_type,
        "lccsd": cc_type,
        "ccsd": cc_type,
        "ccsd(t)": cc_type,
        "a-ccsd(t)": cc_type,
        "omp2": mp2_type,
        "omp2.5": mp_type,
        "omp3": mp_type,
        "oremp2": cc_type,
        "olccd": cc_type,
    }
    corl_type = corl_natural_values[method]

    natural_ref = {"conv": "pk", "df": "df", "cd": "cd"}
    scf_type = inp["keywords"].get("scf_type", natural_ref[corl_type])
    natural_values = {"pk": "pk", "direct": "pk", "df": "df", "mem_df": "df", "disk_df": "df", "cd": "cd"}
    scf_type = natural_values[scf_type]

    atol = 1.0e-6
    if driver == "hessian":
        atol = 1.0e-5  # todo implement more elaborate e/g/h atol like at qcdb & qcng: https://github.com/qcdb/qcdb/blob/master/qcdb/tests/standard_suite_runner.py#L144-L152
    if driver == "gradient" and inp["keywords"]["function_kwargs"]["dertype"] == 0:
        atol = 2.0e-6
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

    driver_call = {"energy": psi4.energy, "gradient": psi4.gradient, "hessian": psi4.hessian}

    psi4.set_options(
        # ADVICE: adding new tests, as soon as you start editing any standard_suite_ref* file, uncomment
        #         "reference generation" block below, so any numbers you may copy are of best quality.
        {
            # reference generation conv crit
            # "guess": "sad",
            # "e_convergence": 10,
            # "d_convergence": 9,
            # "r_convergence": 9,
            # "pcg_convergence": 9,
            # "mo_maxiter": 150,

            # runtime conv crit
            "points": 5,
            "fd_project": False,
        }
    )
    extra_kwargs = inp["keywords"].pop("function_kwargs", {})
    psi4.set_options(inp["keywords"])

    if "error" in inp:
        errtype, errmatch, reason = inp["error"]
        with pytest.raises(errtype) as e:
            driver_call[driver](inp["call"], molecule=subject, **extra_kwargs)

        assert re.search(errmatch, str(e.value)), f"Not found: {errtype} '{errmatch}' in {e.value}"
        # _recorder(qcprog, qc_module_in, driver, method, reference, fcae, scf_type, corl_type, "error", "nyi: " + reason)
        return

    psi4.set_output_file("asdf")

    ret, wfn = driver_call[driver](inp["call"], molecule=subject, return_wfn=True, **extra_kwargs)
    qc_module_out = "psi4-" + ("occ" if wfn.module() == "dfocc" else wfn.module())  # returns "psi4-<module>"

    # <<<  Comparison Tests  >>>

    if qc_module_in != "psi4":
        assert qc_module_out == qc_module_in, f"QC_MODULE used ({qc_module_out}) != requested ({qc_module_in})"
    if qc_module_xptd:
        assert qc_module_out == qc_module_xptd, f"QC_MODULE used ({qc_module_out}) != expected ({qc_module_xptd})"

    ref_block = std_suite[chash]

    # qcvars
    contractual_args = [
        qc_module_out,
        driver,
        reference,
        method,
        corl_type,
        fcae,
    ]
    # ADVICE: adding new tests, when a variable claims "not set", and you're sure it was, investigate core vs wfn.
    #         best practice is to set in C++ code on wfn, then copy wfn to P::e at end of proc.py routine.
    asserter_args = [
        {
            "P::e": psi4.core,
            "wfn": wfn,
        },
        ref_block,
        atol,
        ref_block_conv,
        atol_conv,
        tnm,
    ]

    # ADVICE: adding new tests, comment out lesser methods to focus first on target method
    #         add lesser back later, fulfilled by either adding qcvars to psi4 code or excusing them in qcng contracts
    def qcvar_assertions():
        print("BLOCK", chash, contractual_args)
        if method == "hf":
            _asserter(asserter_args, contractual_args, contractual_hf)
        elif method == "mp2":
            _asserter(asserter_args, contractual_args, contractual_mp2)
        elif method == "mp2.5":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_mp3)
            _asserter(asserter_args, contractual_args, contractual_mp2p5)
        elif method == "mp3":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_mp2p5)
            _asserter(asserter_args, contractual_args, contractual_mp3)
        elif method == "remp2":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_remp2)
        elif method == "lccd":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_lccd)
        elif method == "lccsd":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_lccsd)
        elif method == "ccsd":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_ccsd)
        elif method == "ccsd(t)":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_ccsd)
            _asserter(asserter_args, contractual_args, contractual_ccsd_prt_pr)
        elif method == "a-ccsd(t)":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_ccsd)
            _asserter(asserter_args, contractual_args, contractual_accsd_prt_pr)
        elif method == "omp2":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_omp2)
        elif method == "omp2.5":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_mp3)
            _asserter(asserter_args, contractual_args, contractual_mp2p5)
            _asserter(asserter_args, contractual_args, contractual_omp2p5)
        elif method == "omp3":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_mp3)
            _asserter(asserter_args, contractual_args, contractual_omp3)
        elif method == "oremp2":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_remp2)  # assert skipped
            _asserter(asserter_args, contractual_args, contractual_oremp2)
        elif method == "olccd":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_olccd)

    if "wrong" in inp:
        errmatch, reason = inp["wrong"]
        with pytest.raises(AssertionError) as e:
            qcvar_assertions()

        assert errmatch in str(e.value), f"Not found: AssertionError '{errmatch}' for '{reason}' in {e.value}"
        # _recorder(qcprog, qc_module_out, driver, method, reference, fcae, scf_type, corl_type, "wrong", reason + f" First wrong at `{errmatch}`.")
        pytest.xfail(reason)

    qcvar_assertions()

    # aliases
    _asserter(asserter_args, contractual_args, contractual_current)

    # returns
    tf, errmsg = compare_values(
        ref_block[f"{method.upper()} TOTAL ENERGY"],
        wfn.energy(),
        tnm + " wfn",
        atol=atol,
        return_message=True,
        quiet=True,
    )
    assert compare_values(ref_block[f"{method.upper()} TOTAL ENERGY"], wfn.energy(), tnm + " wfn", atol=atol), errmsg

    if driver == "energy":
        assert compare_values(ref_block[f"{method.upper()} TOTAL ENERGY"], ret, tnm + " return")

    elif driver == "gradient":
        assert compare_values(
            ref_block[f"{method.upper()} TOTAL GRADIENT"], wfn.gradient().np, tnm + " grad wfn", atol=atol
        )
        assert compare_values(ref_block[f"{method.upper()} TOTAL GRADIENT"], ret.np, tnm + " grad return", atol=atol)

    elif driver == "hessian":
        tf, errmsg = compare_values(ref_block[f"{method.upper()} TOTAL HESSIAN"], wfn.hessian().np, tnm + " hess wfn", atol=atol, return_message=True, quiet=True)
        assert compare_values(ref_block[f"{method.upper()} TOTAL HESSIAN"], wfn.hessian().np, tnm + " hess wfn", atol=atol), errmsg
        assert compare_values(ref_block[f"{method.upper()} TOTAL HESSIAN"], wfn.hessian().np, tnm + " hess wfn", atol=atol)
        assert compare_values(
            ref_block[f"{method.upper()} TOTAL GRADIENT"], wfn.gradient().np, tnm + " grad wfn", atol=atol
        )
        assert compare_values(ref_block[f"{method.upper()} TOTAL HESSIAN"], ret.np, tnm + " hess return", atol=atol)

    # generics
    # yapf: disable
    assert compare(ref_block["N BASIS FUNCTIONS"], wfn.nso(), tnm + " nbasis wfn"), f"nbasis {wfn.nso()} != {ref_block['N BASIS FUNCTIONS']}"
    assert compare(ref_block["N MOLECULAR ORBITALS"], wfn.nmo(), tnm + " nmo wfn"), f"nmo {wfn.nmo()} != {ref_block['N MOLECULAR ORBITALS']}"
    assert compare(ref_block["N ALPHA ELECTRONS"], wfn.nalpha(), tnm + " nalpha wfn"), f"nalpha {wfn.nalpha()} != {ref_block['N ALPHA ELECTRONS']}"
    assert compare(ref_block["N BETA ELECTRONS"], wfn.nbeta(), tnm + " nbeta wfn"), f"nbeta {wfn.nbeta()} != {ref_block['N BETA ELECTRONS']}"
    # yapf: enable


def _asserter(asserter_args: List, contractual_args: List[str], contractual_fn: Callable) -> None:
    """For expectations in `contractual_fn`, check that the QCVars are present in P::e.globals and wfn and match expected ref_block."""

    qcvar_stores, ref_block, atol, ref_block_conv, atol_conv, tnm = asserter_args

    for objlbl, obj in qcvar_stores.items():
        for rpv, pv, present in contractual_fn(*contractual_args):
            label = tnm + " " + objlbl + " " + pv

            if present:
                # verify exact match to method (may be df) and near match to conventional (non-df) method
                tf, errmsg = compare_values(
                    ref_block[rpv], query_qcvar(obj, pv), label, atol=atol, return_message=True, quiet=True
                )
                assert compare_values(ref_block[rpv], query_qcvar(obj, pv), label, atol=atol), errmsg

                if isinstance(ref_block_conv[rpv], str) and ref_block_conv[rpv] == "KnownMissing":
                    assert compare(False, False, label + " v. CONV SKIP"), f"{label} wrongly present in CONV"
                else:
                    # ADVICE: adding new tests, when a conventional reference is not present and isn't a priority, comment next two statements.
                    #         To long-term forgive, as last resort, add `_knownmissing` to CONV in qcng _std_suite.
                    tf, errmsg = compare_values(
                        ref_block_conv[rpv], query_qcvar(obj, pv), label + " v. CONV", atol=atol_conv, return_message=True, quiet=True
                    )
                    assert compare_values(ref_block_conv[rpv], query_qcvar(obj, pv), label + " v. CONV", atol=atol_conv), errmsg

                # Note that the double compare_values lines are to collect the errmsg in the first for assertion in the second.
                #   If the errmsg isn't present in the assert, the string isn't accessible through `e.value`.
                #   If a plain bool is compared in the assert, the printed message will show booleans and not numbers.
            else:
                # verify and forgive known contract violations
                assert compare(False, query_has_qcvar(obj, pv), label + " SKIP"), f"{label} wrongly present"
