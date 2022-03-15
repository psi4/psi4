import pprint
import re
from typing import Any, Dict

import numpy as np
import pytest
from qcelemental.models import AtomicInput
from qcelemental.molutil import compute_scramble
from qcelemental.testing import compare, compare_values

import qcengine as qcng
from qcengine.programs.util import mill_qcvars

from .standard_suite_ref import answer_hash, std_suite

from .standard_suite_contracts import (  # isort:skip
    contractual_hf,
    contractual_mp2,
    contractual_mp2p5,
    contractual_mp3,
    contractual_mp4_prsdq_pr,
    contractual_mp4,
    contractual_cisd,
    contractual_qcisd,
    contractual_qcisd_prt_pr,
    contractual_lccd,
    contractual_lccsd,
    contractual_ccd,
    contractual_ccsd,
    contractual_ccsdpt_prccsd_pr,
    contractual_ccsd_prt_pr,
    contractual_accsd_prt_pr,
    contractual_ccsdt1a,
    contractual_ccsdt1b,
    contractual_ccsdt2,
    contractual_ccsdt3,
    contractual_ccsdt,
    contractual_ccsdt_prq_pr,
    contractual_ccsdtq,
    contractual_dft_current,
    contractual_current,
    query_has_qcvar,
    query_qcvar,
)

pp = pprint.PrettyPrinter(width=120)


def runner_asserter(inp, ref_subject, method, basis, tnm, scramble, frame):

    qcprog = inp["call"]
    qc_module_in = inp["qc_module"]  # returns "<qcprog>"|"<qcprog>-<module>"  # input-specified routing
    qc_module_xptd = (
        (qcprog + "-" + inp["xptd"]["qc_module"]) if inp.get("xptd", {}).get("qc_module", None) else None
    )  # expected routing
    driver = inp["driver"]
    reference = inp["reference"]
    fcae = inp["fcae"]

    if basis == "cfour-qz2p" and qcprog in ["gamess", "nwchem", "qchem"]:
        pytest.skip(f"basis {basis} not available in {qcprog} library")

    # <<<  Molecule  >>>

    # 1. ref mol: `ref_subject` nicely oriented mol taken from standard_suite_ref.py
    min_nonzero_coords = np.count_nonzero(np.abs(ref_subject.geometry) > 1.0e-10)

    if scramble is None:
        subject = ref_subject
        ref2in_mill = compute_scramble(
            len(subject.symbols), do_resort=False, do_shift=False, do_rotate=False, do_mirror=False
        )  # identity AlignmentMill

    else:
        subject, data = ref_subject.scramble(**scramble, do_test=False)
        ref2in_mill = data["mill"]

    # 2. input mol: `subject` now ready for `atin.molecule`. may have been scrambled away from nice ref orientation

    # <<<  Reference Values  >>>

    # ? precedence on these types -- not really applicable to qcng
    mp2_type = inp.get("corl_type", inp["keywords"].get("mp2_type", "df"))  # hard-code of read_options.cc MP2_TYPE
    mp_type = inp.get("corl_type", inp["keywords"].get("mp_type", "conv"))  # hard-code of read_options.cc MP_TYPE
    ci_type = inp.get("corl_type", inp["keywords"].get("ci_type", "conv"))  # hard-code of read_options.cc CI_TYPE
    cc_type = inp.get("corl_type", inp["keywords"].get("cc_type", "conv"))  # hard-code of read_options.cc CC_TYPE
    corl_natural_values = {
        "hf": "conv",  # dummy to assure df/cd/conv scf_type refs available
        "mp2": mp2_type,
        "mp3": mp_type,
        "mp4(sdq)": mp_type,
        "mp4": mp_type,
        "cisd": ci_type,
        "qcisd": ci_type,
        "qcisd(t)": ci_type,
        "lccd": cc_type,
        "lccsd": cc_type,
        "ccd": cc_type,
        "ccsd": cc_type,
        "ccsd+t(ccsd)": cc_type,
        "ccsd(t)": cc_type,
        "a-ccsd(t)": cc_type,
        "ccsdt-1a": cc_type,
        "ccsdt-1b": cc_type,
        "ccsdt-2": cc_type,
        "ccsdt-3": cc_type,
        "ccsdt": cc_type,
        "ccsdt(q)": cc_type,
        "ccsdtq": cc_type,
        "pbe": "conv",
        "b3lyp": "conv",
        "b3lyp5": "conv",
    }
    corl_type = corl_natural_values[method]

    natural_ref = {"conv": "pk", "df": "df", "cd": "cd"}
    scf_type = inp["keywords"].get("scf_type", natural_ref[corl_type])
    natural_values = {"pk": "pk", "direct": "pk", "df": "df", "mem_df": "df", "disk_df": "df", "cd": "cd"}
    scf_type = natural_values[scf_type]

    is_dft = method in ["pbe", "b3lyp", "b3lyp5"]

    atol_e, rtol_e = 2.0e-7, 1.0e-16
    atol_g, rtol_g = 5.0e-7, 2.0e-5
    atol_h, rtol_h = 1.0e-5, 2.0e-5
    if is_dft:
        atol_g = 6.0e-6
    chash = answer_hash(
        system=subject.name,
        basis=basis,
        fcae=fcae,
        scf_type=scf_type,
        reference=reference,
        corl_type=corl_type,
    )
    ref_block = std_suite[chash]

    # check all calcs against conventional reference to looser tolerance
    atol_conv = 1.0e-4
    rtol_conv = 1.0e-3
    chash_conv = answer_hash(
        system=subject.name,
        basis=basis,
        fcae=fcae,
        reference=reference,
        corl_type="conv",
        scf_type="pk",
    )
    ref_block_conv = std_suite[chash_conv]

    # <<<  Prepare Calculation and Call API  >>>

    atin = AtomicInput(
        **{
            "molecule": subject,
            "driver": driver,
            "model": {
                "method": method,
                "basis": inp.get("basis", "(auto)"),
            },
            "keywords": inp["keywords"],
        }
    )

    local_options = {}
    local_options = {"nnodes": 1, "ncores": 2}  # debug. temp fix low for testing to avoid too many files open w/gamess

    if "error" in inp:
        errtype, errmatch, reason = inp["error"]
        with pytest.raises(errtype) as e:
            qcng.compute(atin, qcprog, raise_error=True, return_dict=True, local_options=local_options)

        assert re.search(errmatch, str(e.value)), f"Not found: {errtype} '{errmatch}' in {e.value}"
        # _recorder(qcprog, qc_module_in, driver, method, reference, fcae, scf_type, corl_type, "error", "nyi: " + reason)
        return

    wfn = qcng.compute(atin, qcprog, raise_error=True, local_options=local_options)

    print("WFN")
    pp.pprint(wfn.dict())

    qc_module_out = wfn.provenance.creator.lower()
    if hasattr(wfn.provenance, "module"):
        qc_module_out += "-" + wfn.provenance.module  # returns "<qcprog>-<module>"
    # assert 0, f"{qc_module_xptd=} {qc_module_in=} {qc_module_out=}"  # debug

    # 3. output mol: `wfn.molecule` after calc. orientation for nonscalar quantities may be different from `subject` if fix_=False

    _, data = ref_subject.align(wfn.molecule, atoms_map=False, mols_align=True, verbose=0)
    ref2out_mill = data["mill"]

    if subject.fix_com and subject.fix_orientation:
        with np.printoptions(precision=3, suppress=True):
            assert compare_values(
                subject.geometry, wfn.molecule.geometry, atol=5.0e-8
            ), f"coords: atres ({wfn.molecule.geometry}) != atin ({subject.geometry})"  # 10 too much
        assert (
            ref_subject.fix_com
            and ref_subject.fix_orientation
            and subject.fix_com
            and subject.fix_orientation
            and wfn.molecule.fix_com
            and wfn.molecule.fix_orientation
        ), f"fixed, so all T: {ref_subject.fix_com} {ref_subject.fix_orientation} {subject.fix_com} {subject.fix_orientation} {wfn.molecule.fix_com} {wfn.molecule.fix_orientation}"

        ref_block = mill_qcvars(ref2in_mill, ref_block)
        ref_block_conv = mill_qcvars(ref2in_mill, ref_block_conv)

    else:
        # this check assumes the qcprog will adjust an ugly Cartesian geometry into a pretty one (with more symmetry for computational efficiency).
        # if qcprog doesn't have that behavior, it will need to be excused from this check.
        with np.printoptions(precision=3, suppress=True):
            assert compare(
                min_nonzero_coords, np.count_nonzero(np.abs(wfn.molecule.geometry) > 1.0e-10), tnm + " !0 coords wfn"
            ), f"count !0 coords {wfn.molecule.geometry} != {min_nonzero_coords}"
        assert (
            (not ref_subject.fix_com)
            and (not ref_subject.fix_orientation)
            and (not subject.fix_com)
            and (not subject.fix_orientation)
            and (not wfn.molecule.fix_com)
            and (not wfn.molecule.fix_orientation)
        ), f"free, so all F: {ref_subject.fix_com} {ref_subject.fix_orientation} {subject.fix_com} {subject.fix_orientation} {wfn.molecule.fix_com} {wfn.molecule.fix_orientation}"

        ref_block = mill_qcvars(ref2out_mill, ref_block)
        ref_block_conv = mill_qcvars(ref2out_mill, ref_block_conv)

    # <<<  Comparison Tests  >>>

    assert wfn.success is True
    if qc_module_in != qcprog:
        assert qc_module_out == qc_module_in, f"QC_MODULE used ({qc_module_out}) != requested ({qc_module_in})"
    if qc_module_xptd:
        assert qc_module_out == qc_module_xptd, f"QC_MODULE used ({qc_module_out}) != expected ({qc_module_xptd})"

    # qcvars
    contractual_args = [
        qc_module_out,
        driver,
        reference,
        method,
        corl_type,
        fcae,
    ]
    asserter_args = [
        [wfn.extras["qcvars"], wfn.properties],
        ref_block,
        [atol_e, atol_g, atol_h],
        [rtol_e, rtol_g, rtol_h],
        ref_block_conv,
        atol_conv,
        rtol_conv,
        tnm,
    ]

    def qcvar_assertions():
        print("BLOCK", chash, contractual_args)
        if method == "hf":
            _asserter(asserter_args, contractual_args, contractual_hf)
        elif method == "mp2":
            _asserter(asserter_args, contractual_args, contractual_mp2)
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
        elif method == "cisd":
            _asserter(asserter_args, contractual_args, contractual_cisd)
        elif method == "qcisd":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_qcisd)
        elif method == "qcisd(t)":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_qcisd)
            _asserter(asserter_args, contractual_args, contractual_qcisd_prt_pr)
        elif method == "lccd":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_lccd)
        elif method == "lccsd":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_lccsd)
        elif method == "ccd":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_ccd)
        elif method == "ccsd":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_ccsd)
        elif method == "ccsd+t(ccsd)":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_ccsd)
            _asserter(asserter_args, contractual_args, contractual_ccsdpt_prccsd_pr)
        elif method == "ccsd(t)":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_ccsd)
            _asserter(asserter_args, contractual_args, contractual_ccsd_prt_pr)
        elif method == "a-ccsd(t)":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_ccsd)
            _asserter(asserter_args, contractual_args, contractual_accsd_prt_pr)
        elif method == "ccsdt-1a":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_ccsdt1a)
        elif method == "ccsdt-1b":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_ccsdt1b)
        elif method == "ccsdt-2":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_ccsdt2)
        elif method == "ccsdt-3":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_ccsdt3)
        elif method == "ccsdt":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_ccsdt)
        elif method == "ccsdt(q)":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_ccsdt)
            _asserter(asserter_args, contractual_args, contractual_ccsdt_prq_pr)
        elif method == "ccsdtq":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_ccsdtq)
        # separations here for DFT appropriate when qcvars are labeled by functional

    if "wrong" in inp:
        errmatch, reason = inp["wrong"]
        with pytest.raises(AssertionError) as e:
            qcvar_assertions()

        assert errmatch in str(e.value), f"Not found: AssertionError '{errmatch}' for '{reason}' in {e.value}"
        # _recorder(qcprog, qc_module_out, driver, method, reference, fcae, scf_type, corl_type, "wrong", reason + f" First wrong at `{errmatch}`.")
        pytest.xfail(reason)

    # primary label checks
    qcvar_assertions()

    # aliases checks
    asserter_args[0].pop()  # checks not appropriate for properties
    if is_dft:
        _asserter(asserter_args, contractual_args, contractual_dft_current)
    else:
        _asserter(asserter_args, contractual_args, contractual_current)

    # returns checks
    if driver == "energy":
        tf, errmsg = compare_values(
            ref_block[f"{method.upper()} TOTAL ENERGY"],
            wfn.return_result,
            tnm + " wfn",
            atol=atol_e,
            rtol=rtol_e,
            return_message=True,
            quiet=True,
        )
        assert compare_values(
            ref_block[f"{method.upper()} TOTAL ENERGY"], wfn.return_result, tnm + " wfn", atol=atol_e, rtol=rtol_e
        ), errmsg

    elif driver == "gradient":
        tf, errmsg = compare_values(
            ref_block[f"{method.upper()} TOTAL GRADIENT"],
            wfn.return_result,
            tnm + " grad wfn",
            atol=atol_g,
            rtol=rtol_g,
            return_message=True,
            quiet=True,
        )
        assert compare_values(
            ref_block[f"{method.upper()} TOTAL GRADIENT"],
            wfn.return_result,
            tnm + " grad wfn",
            atol=atol_g,
            rtol=rtol_g,
        ), errmsg

    elif driver == "hessian":
        tf, errmsg = compare_values(
            ref_block[f"{method.upper()} TOTAL HESSIAN"],
            wfn.return_result,
            tnm + " hess wfn",
            atol=atol_h,
            rtol=rtol_h,
            return_message=True,
            quiet=True,
        )
        assert compare_values(
            ref_block[f"{method.upper()} TOTAL HESSIAN"], wfn.return_result, tnm + " hess wfn", atol=atol_h, rtol=rtol_h
        ), errmsg

    if driver in ["energy", "gradient", "hessian"]:
        tf, errmsg = compare_values(
            ref_block[f"{method.upper()} TOTAL ENERGY"],
            wfn.properties.return_energy,
            tnm + " prop",
            atol=atol_e,
            rtol=rtol_e,
            return_message=True,
            quiet=True,
        )
        assert compare_values(
            ref_block[f"{method.upper()} TOTAL ENERGY"],
            wfn.properties.return_energy,
            tnm + " prop",
            atol=atol_e,
            rtol=rtol_e,
        ), errmsg

    if driver in ["gradient"]:  # , "hessian"]:
        tf, errmsg = compare_values(
            ref_block[f"{method.upper()} TOTAL GRADIENT"],
            wfn.properties.return_gradient,
            tnm + " grad prop",
            atol=atol_g,
            rtol=rtol_g,
            return_message=True,
            quiet=True,
        )
        assert compare_values(
            ref_block[f"{method.upper()} TOTAL GRADIENT"],
            wfn.properties.return_gradient,
            tnm + " grad prop",
            atol=atol_g,
            rtol=rtol_g,
        ), errmsg

    if driver == "hessian":
        tf, errmsg = compare_values(
            ref_block[f"{method.upper()} TOTAL HESSIAN"],
            wfn.properties.return_hessian,
            tnm + " hess prop",
            atol=atol_h,
            rtol=rtol_h,
            return_message=True,
            quiet=True,
        )
        assert compare_values(
            ref_block[f"{method.upper()} TOTAL HESSIAN"],
            wfn.properties.return_hessian,
            tnm + " hess prop",
            atol=atol_h,
            rtol=rtol_h,
        ), errmsg

    # generics
    # yapf: disable
    assert compare(ref_block["N BASIS FUNCTIONS"], wfn.properties.calcinfo_nbasis, tnm + " nbasis wfn"), f"nbasis {wfn.properties.calcinfo_nbasis} != {ref_block['N BASIS FUNCTIONS']}"
    assert compare(ref_block["N MOLECULAR ORBITALS"], wfn.properties.calcinfo_nmo, tnm + " nmo wfn"), f"nmo {wfn.properties.calcinfo_nmo} != {ref_block['N MOLECULAR ORBITALS']}"
    assert compare(ref_block["N ALPHA ELECTRONS"], wfn.properties.calcinfo_nalpha, tnm + " nalpha wfn"), f"nalpha {wfn.properties.calcinfo_nalpha} != {ref_block['N ALPHA ELECTRONS']}"
    assert compare(ref_block["N BETA ELECTRONS"], wfn.properties.calcinfo_nbeta, tnm + " nbeta wfn"), f"nbeta {wfn.properties.calcinfo_nbeta} != {ref_block['N BETA ELECTRONS']}"
    # yapf: enable

    # record
    # _recorder(qcprog, qc_module_out, driver, method, reference, fcae, scf_type, corl_type, "fd" if using_fd else "pass", "")

    # assert 0


def _asserter(asserter_args, contractual_args, contractual_fn):
    """For expectations in `contractual_fn`, check that the QCVars are present in P::e.globals and wfn and match expected ref_block."""

    qcvar_stores, ref_block, atol_egh, rtol_egh, ref_block_conv, atol_conv, rtol_conv, tnm = asserter_args

    for obj in qcvar_stores:
        for rpv, pv, present in contractual_fn(*contractual_args):
            label = tnm + " " + pv
            atol = atol_egh["EGH".index(rpv.split()[-1][0])]
            rtol = rtol_egh["EGH".index(rpv.split()[-1][0])]

            if present:
                # verify exact match to method (may be df) and near match to conventional (non-df) method
                tf, errmsg = compare_values(
                    ref_block[rpv], query_qcvar(obj, pv), label, atol=atol, rtol=rtol, return_message=True, quiet=True
                )
                assert compare_values(ref_block[rpv], query_qcvar(obj, pv), label, atol=atol, rtol=rtol), errmsg
                tf, errmsg = compare_values(
                    ref_block_conv[rpv],
                    query_qcvar(obj, pv),
                    label,
                    atol=atol_conv,
                    rtol=rtol_conv,
                    return_message=True,
                    quiet=True,
                )
                assert compare_values(
                    ref_block_conv[rpv], query_qcvar(obj, pv), label, atol=atol_conv, rtol=rtol_conv
                ), errmsg

                # Note that the double compare_values lines are to collect the errmsg in the first for assertion in the second.
                #   If the errmsg isn't present in the assert, the string isn't accessible through `e.value`.
                #   If a plain bool is compared in the assert, the printed message will show booleans and not numbers.
            else:
                # verify and forgive known contract violations
                assert compare(False, query_has_qcvar(obj, pv), label + " SKIP")
