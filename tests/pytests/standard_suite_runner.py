import re
from typing import Callable, List

import pytest
from qcengine.programs.tests.standard_suite_contracts import *
from qcengine.programs.tests.standard_suite_ref import answer_hash
from standard_suite_ref_local import std_suite, contractual_mp2f12

import psi4

from utils import *


def runner_asserter(inp, subject, method, basis, tnm):

    qcprog = "psi4"
    qc_module_in = "-".join([qcprog, inp["keywords"].get("qc_module", "")]).strip("-")  # returns "psi4"|"psi4-<module>"  # input-specified routing
    qc_module_xptd = (
        (qcprog + "-" + inp["xptd"]["qc_module"]) if inp.get("xptd", {}).get("qc_module", None) else None
    )  # expected routing
    driver = inp["driver"]
    reference = inp["keywords"]["reference"]
    fcae = {"true": "fc", "false": "ae"}[inp["keywords"]["freeze_core"]]
    sdsc = inp.get("xptd", {}).get("sdsc", False) or ("sc" if reference == "rohf" else "sd")

    if (qc_module_in == "psi4-detci" or qc_module_xptd == "psi4-detci") and basis != "cc-pvdz":
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
    ci_type = inp.get("corl_type", inp["keywords"].get("ci_type", "conv"))  # hard-code of read_options.cc CI_TYPE
    cc_type = inp.get("corl_type", inp["keywords"].get("cc_type", "conv"))  # hard-code of read_options.cc CC_TYPE
    corl_natural_values = {
        "hf": "df",  # dummy to assure df/cd/conv scf_type refs available
        "mp2": mp2_type,
        "mp2.5": mp_type,
        "mp3": mp_type,
        "mp4(sdq)": mp_type,
        "mp4": mp_type,
        "zapt2": mp_type,
        "cisd": ci_type,
        "qcisd": ci_type,
        "qcisd(t)": ci_type,
        "fci": ci_type,
        "remp2": cc_type,
        "lccd": cc_type,
        "lccsd": cc_type,
        "cepa(1)": cc_type,
        "cepa(3)": cc_type,
        "acpf": cc_type,
        "aqcc": cc_type,
        "ccd": cc_type,
        "bccd": cc_type,
        "cc2": cc_type,
        "ccsd": cc_type,
        "ccsd(t)": cc_type,
        "a-ccsd(t)": cc_type,
        "bccd(t)": cc_type,
        "cc3": cc_type,
        "omp2": mp2_type,
        "omp2.5": mp_type,
        "omp3": mp_type,
        "oremp2": cc_type,
        "olccd": cc_type,
        "mp2-f12": mp2_type,
        "svwn": "df",  # dummy
        "pbe": "df",  # dummy
        "b3lyp": "df",  # dummy
        "wb97x": "df",  # dummy
        "b2plyp": "df",  # dummy
    }
    corl_type = corl_natural_values[method]

    natural_ref = {"conv": "pk", "df": "df", "cd": "cd"}
    natural_ref_inv = {v: k for k, v in natural_ref.items()}
    scf_type = inp["keywords"].get("scf_type", natural_ref[corl_type])
    natural_values = {"pk": "pk", "direct": "pk", "df": "df", "mem_df": "df", "disk_df": "df", "cd": "cd"}
    scf_type = natural_values[scf_type]

    is_dft = method in ["svwn", "pbe", "b3lyp", "wb97x", "b3lyp5", "b2plyp"]
    is_dhdft = method in ["b2plyp"]

    if is_dhdft:
        corl_type = inp["keywords"].get("mp2_type", corl_natural_values["mp2"])
    elif is_dft:
        # Hartree--Fock doesn't need this b/c it does value sharing by variable so CONV-AE-CONV also stored at CONV-AE-DF
        corl_type = inp["keywords"].get("mp2_type", natural_ref_inv[scf_type])

    using_fd = (driver == "gradient" and inp["keywords"]["function_kwargs"]["dertype"] == 0 or driver == "hessian" and inp["keywords"]["function_kwargs"]["dertype"] in [0, 1])  # T/F: notate fd vs. analytic for docs table
    defaultness = "defaultdefault" if ("default" in tnm and all([not k.endswith("_type") for k in inp["keywords"]])) else ("default" if "default" in tnm else "")

    # * absolute and relative tolerances function approx as `or` operation. see https://numpy.org/doc/stable/reference/generated/numpy.allclose.html
    atol_e, rtol_e = 1.0e-6, 1.0e-16  # qcdb: 2.0e-7, 1.0e-16
    atol_g, rtol_g = 1.0e-6, 1.0e-16  # qcdb: 5.0e-7, 2.0e-5
    atol_h, rtol_h = 2.0e-6, 1.0e-16  # qcdb: 1.0e-5, 2.0e-5
    if method == "wb97x":
        # lrc very sensitive to grid, so ref values are computed at 500,1202 at which tightness grd1==grd0. so, at runtime default grid, need to relax tol to match ref
        atol_e = 2.5e-6
        atol_g = 1.0e-5

    if driver == "gradient" and inp["keywords"]["function_kwargs"]["dertype"] == 0:
        # relax to this pre-e/g/h-separated value if necessary: atol = 2.0e-6
        pass
    if driver == "hessian" and inp["keywords"]["function_kwargs"]["dertype"] in [0, 1]:
        atol_h = 1.0e-5
        if is_dft:
            atol_h = 7.0e-5


    chash = answer_hash(
        system=subject.name(), basis=basis, fcae=fcae, scf_type=scf_type, reference=reference, corl_type=corl_type, sdsc=sdsc
    )

    # check all calcs against conventional reference to looser tolerance
    atol_conv = 1.5e-4
    rtol_conv = 1.0e-16  # qcdb 1.0e-3
    if (method in ["ccsd", "ccsd(t)"] and corl_type in ["df", "cd"]) or (method == "wb97x"):
        # fnocc non-conv ccsd opposite-spin df error is larger than corl df error. TODO: see if other methods/module/mols show this or if suspicious.
        atol_conv = 3.0e-4
    if method in ["mp2-f12"] and corl_type in ["df"]:
        # hf-cabs is sensitive
        atol_conv = 2.0e-3
    chash_conv = answer_hash(
        system=subject.name(), basis=basis, fcae=fcae, reference=reference, corl_type="conv", scf_type="pk", sdsc=sdsc
    )
    ref_block_conv = std_suite[chash_conv]

    # <<<  Prepare Calculation and Call API  >>>

    driver_call = {"energy": psi4.energy, "gradient": psi4.gradient, "hessian": psi4.hessian}

    # psi4.set_memory("3 GiB")
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
            # "brueckner_orbs_r_convergence": 9,

            # runtime conv crit
            "points": 5,
            "fd_project": False,
        }
    )
    extra_kwargs = inp["keywords"].pop("function_kwargs", {})

    if "error" in inp:
        errtype, errmatch, reason = inp["error"]
        with pytest.raises(errtype) as e:
            psi4.set_options(inp["keywords"])
            driver_call[driver](inp["call"], molecule=subject, **extra_kwargs)

        assert re.search(errmatch, str(e.value)), f"Not found: {errtype} '{errmatch}' in {e.value}"
        if "scftype" not in tnm:
            _recorder(qcprog, qc_module_in, driver, method, reference, fcae, sdsc, scf_type, corl_type, "error", "nyi: " + reason)
        return

    psi4.set_output_file("asdf")  # easy name to find output files. TODO: why doesn't .out remain w/o this?

    psi4.set_options(inp["keywords"])
    ret, wfn = driver_call[driver](inp["call"], molecule=subject, return_wfn=True, **extra_kwargs)
    qc_module_out = qcprog + "-" + ("occ" if wfn.module() == "dfocc" else wfn.module())  # returns "psi4-<module>"

    # <<<  Comparison Tests  >>>

    # routing checks
    if qc_module_in != qcprog:
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
        sdsc,
    ]
    # ADVICE: adding new tests, when a variable claims "not set", and you're sure it was, investigate core vs wfn.
    #         best practice is to set in C++ code on wfn, then copy wfn to P::e at end of proc.py routine.
    asserter_args = [
        {
            "P::e": psi4.core,
            "wfn": wfn,
        },
        ref_block,
        [atol_e, atol_g, atol_h],
        [rtol_e, rtol_g, rtol_h],
        ref_block_conv,
        atol_conv,
        rtol_conv,
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
        elif method == "mp4(sdq)":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_mp2p5)
            _asserter(asserter_args, contractual_args, contractual_mp3)
            _asserter(asserter_args, contractual_args, contractual_mp4_prsdq_pr)
        elif method == "mp4":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_mp2p5)
            _asserter(asserter_args, contractual_args, contractual_mp3)
            _asserter(asserter_args, contractual_args, contractual_mp4_prsdq_pr)
            _asserter(asserter_args, contractual_args, contractual_mp4)
        elif method == "zapt2":
            _asserter(asserter_args, contractual_args, contractual_zapt2)
        elif method == "cisd":
            _asserter(asserter_args, contractual_args, contractual_cisd)
        elif method == "qcisd":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_qcisd)
        elif method == "qcisd(t)":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_qcisd)
            _asserter(asserter_args, contractual_args, contractual_qcisd_prt_pr)
        elif method == "fci":
            _asserter(asserter_args, contractual_args, contractual_fci)
        elif method == "remp2":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_remp2)
        elif method == "lccd":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_lccd)
        elif method == "lccsd":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_lccsd)
        elif method == "cepa(1)":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_cepa_pr1_pr)
        elif method == "cepa(3)":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_cepa_pr3_pr)
        elif method == "acpf":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_acpf)
        elif method == "aqcc":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_aqcc)
        elif method == "ccd":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_ccd)
        elif method == "bccd":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_ccsd)
            _asserter(asserter_args, contractual_args, contractual_bccd)
        elif method == "cc2":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_cc2)
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
        elif method == "bccd(t)":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_ccsd)
            _asserter(asserter_args, contractual_args, contractual_bccd)
            _asserter(asserter_args, contractual_args, contractual_ccsd_prt_pr)  # assert skipped
            _asserter(asserter_args, contractual_args, contractual_bccd_prt_pr)
        elif method == "cc3":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_cc3)
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
            _asserter(asserter_args, contractual_args, contractual_lccd)  # assert skipped
            _asserter(asserter_args, contractual_args, contractual_olccd)
        elif method == "mp2-f12":
            # TODO resolve df not matching ref, may need exclusions in qcng
            # _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_mp2f12)

    if "wrong" in inp:
        errmatch, reason = inp["wrong"]
        with pytest.raises(AssertionError) as e:
            qcvar_assertions()

        assert errmatch in str(e.value), f"Not found: AssertionError '{errmatch}' for '{reason}' in {e.value}"
        if "scftype" not in tnm:
            _recorder(qcprog, qc_module_out, driver, method, reference, fcae, sdsc, scf_type, corl_type, "wrong", reason + f" First wrong at `{errmatch}`.")
        pytest.xfail(reason)

    # primary label checks
    qcvar_assertions()

    # aliases checks
    if is_dhdft:
        _asserter(asserter_args, contractual_args, contractual_dhdft_current)
    elif is_dft:
        _asserter(asserter_args, contractual_args, contractual_dft_current)
    else:
        _asserter(asserter_args, contractual_args, contractual_current)

    # returns checks
    tf, errmsg = compare_values(
        ref_block[f"{method.upper()} TOTAL ENERGY"],
        wfn.energy(),
        tnm + " wfn",
        atol=atol_e,
        rtol=rtol_e,
        return_message=True,
        quiet=True,
    )
    assert compare_values(ref_block[f"{method.upper()} TOTAL ENERGY"], wfn.energy(), tnm + " wfn", atol=atol_e, rtol=rtol_e), errmsg

    if driver == "energy":
        assert compare_values(ref_block[f"{method.upper()} TOTAL ENERGY"], ret, tnm + " return", atol=atol_e, rtol=rtol_e)

    elif driver == "gradient":
        assert compare_values(
            ref_block[f"{method.upper()} TOTAL GRADIENT"], wfn.gradient().np, tnm + " grad wfn", atol=atol_g, rtol=rtol_g
        )
        assert compare_values(ref_block[f"{method.upper()} TOTAL GRADIENT"], ret.np, tnm + " grad return", atol=atol_g, rtol=rtol_g)

    elif driver == "hessian":
        tf, errmsg = compare_values(ref_block[f"{method.upper()} TOTAL HESSIAN"], wfn.hessian().np, tnm + " hess wfn", atol=atol_h, rtol=rtol_h, return_message=True, quiet=True)
        assert compare_values(ref_block[f"{method.upper()} TOTAL HESSIAN"], wfn.hessian().np, tnm + " hess wfn", atol=atol_h, rtol=rtol_h), errmsg
        assert compare_values(
            ref_block[f"{method.upper()} TOTAL GRADIENT"], wfn.gradient().np, tnm + " grad wfn", atol=atol_g, rtol=rtol_g
        )
        assert compare_values(ref_block[f"{method.upper()} TOTAL HESSIAN"], ret.np, tnm + " hess return", atol=atol_h, rtol=rtol_h)

    # generics checks
    # yapf: disable
    assert compare(ref_block["N BASIS FUNCTIONS"], wfn.nso(), tnm + " nbasis wfn"), f"nbasis {wfn.nso()} != {ref_block['N BASIS FUNCTIONS']}"
    assert compare(ref_block["N MOLECULAR ORBITALS"], wfn.nmo(), tnm + " nmo wfn"), f"nmo {wfn.nmo()} != {ref_block['N MOLECULAR ORBITALS']}"
    assert compare(ref_block["N ALPHA ELECTRONS"], wfn.nalpha(), tnm + " nalpha wfn"), f"nalpha {wfn.nalpha()} != {ref_block['N ALPHA ELECTRONS']}"
    assert compare(ref_block["N BETA ELECTRONS"], wfn.nbeta(), tnm + " nbeta wfn"), f"nbeta {wfn.nbeta()} != {ref_block['N BETA ELECTRONS']}"
    # yapf: enable

    if "scftype" not in tnm:
        _recorder(qcprog, qc_module_out, driver, method, reference, fcae, sdsc, scf_type, corl_type, "fd" if using_fd else "pass", defaultness)
        if "default" in tnm:
            # add'l entry for "default" line
            _recorder(qcprog, "aaaa-", driver, method, reference, fcae, sdsc, scf_type, corl_type, "fd" if using_fd else "pass", defaultness if defaultness == "defaultdefault" else "")

    # assert 0


def _asserter(asserter_args: List, contractual_args: List[str], contractual_fn: Callable) -> None:
    """For expectations in `contractual_fn`, check that the QCVars are present in P::e.globals and wfn and match expected ref_block."""

    qcvar_stores, ref_block, atol_egh, rtol_egh, ref_block_conv, atol_conv, rtol_conv, tnm = asserter_args

    for objlbl, obj in qcvar_stores.items():
        for rpv, pv, present in contractual_fn(*contractual_args):
            label = tnm + " " + objlbl + " " + pv
            atol = atol_egh["EGH".index(rpv.split()[-1][0])]
            rtol = rtol_egh["EGH".index(rpv.split()[-1][0])]

            if present:
                # verify exact match to method (may be df) and near match to conventional (non-df) method
                tf, errmsg = compare_values(
                    ref_block[rpv], query_qcvar(obj, pv), label, atol=atol, rtol=rtol, return_message=True, quiet=True
                )
                assert compare_values(ref_block[rpv], query_qcvar(obj, pv), label, atol=atol, rtol=rtol), errmsg

                if isinstance(ref_block_conv[rpv], str) and ref_block_conv[rpv] == "KnownMissing":
                    assert compare(False, False, label + " v. CONV SKIP"), f"{label} wrongly present in CONV"
                else:
                    # ADVICE: adding new tests, when a conventional reference is not present and isn't a priority, comment next two statements.
                    #         To long-term forgive, as last resort, add `_knownmissing` to CONV in qcng _std_suite.
                    tf, errmsg = compare_values(
                        ref_block_conv[rpv],
                        query_qcvar(obj, pv),
                        label + " v. CONV",
                        atol=atol_conv,
                        rtol=rtol_conv,
                        return_message=True,
                        quiet=True,
                    )
                    assert compare_values(ref_block_conv[rpv], query_qcvar(obj, pv), label + " v. CONV", atol=atol_conv, rtol=rtol_conv), errmsg

                # Note that the double compare_values lines are to collect the errmsg in the first for assertion in the second.
                #   If the errmsg isn't present in the assert, the string isn't accessible through `e.value`.
                #   If a plain bool is compared in the assert, the printed message will show booleans and not numbers.
            else:
                # verify and forgive known contract violations
                assert compare(False, query_has_qcvar(obj, pv), label + " SKIP"), f"{label} wrongly present"


def _recorder(engine, module, driver, method, reference, fcae, sdsc, scf_type, corl_type, status, note):
    """Accumulates record of stdsuite tests circumstances and status.
    See ../../psi4/share/psi4/scripts/merge_stdsuite.py for discussion on viewing results."""

    with open("stdsuite_psi4.txt", "a") as fp:
        stuff = {
            "module": module,
            "driver": driver,
            "method": method,
            "reference": reference,
            "fcae": fcae,
            "sdsc": sdsc,
            "scf_type": scf_type,
            "corl_type": corl_type,
            "status": status,
            "note": note,
        }
        fp.write(f"{stuff!r}\n")
