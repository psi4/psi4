# This file consolidates the expectations of what properties/QCVariables
#  a method should produce and what in practice a QC program does produce.

from typing import Any, Tuple

from qcengine.programs.qcvar_identities_resources import qcvars_to_atomicproperties


def query_qcvar(obj: Any, pv: str) -> Any:
    """Uniform interface to value of variable `pv` or its QCSchema alias in QCVariable store `obj`."""

    try:
        # psi4.core (module -- P::e.globals)
        # psi4.Wavefunction
        # qcdb (module)
        vval = obj.variable(pv)
    except AttributeError:
        try:
            # qcdb jobrec["qcvars"] qcel.Datum
            vval = obj[pv].data
        except (AttributeError, KeyError):
            # qcel.AtomicResult.extras["qcvars"]
            vval = obj.get(pv)
            if vval is None:
                # qcel.AtomicResult.properties dict
                vval = obj.get(qcvars_to_atomicproperties[pv])
        except TypeError:
            # qcel.AtomicResult.properties object
            vval = getattr(obj, qcvars_to_atomicproperties[pv])

    return vval


def query_has_qcvar(obj: Any, pv: str) -> bool:
    """Uniform interface to whether variable `pv` or its QCSchema alias in QCVariable store `obj`."""

    try:
        bval = obj.has_variable(pv)
    except AttributeError:
        bval = pv in obj

        if not bval:
            try:
                bval = qcvars_to_atomicproperties[pv] in obj
            except KeyError:
                bval = False

    return bval


_contractual_docstring = """
    Parameters
    ----------
    qc_module
        The program or subprogram running the job (e.g., "cfour" or "cfour-ecc").
    driver
        {"energy", "gradient", "hessian"}
        The derivative level that should be expected.
    reference
        {"rhf", "uhf", "rohf"}
        The SCF reference since programs often output differently based on it.
    method
        The target AtomicInput.model.method since "free" methods may not always be
        output (e.g., MP2 available when target is MP2 but not when target is CCSD).
    corl_type
        {"conv", "df", "cd"}
        The algorithm for the target method since programs often output differently
        based on it.
    fcae
        {"ae", "fc"}
        The all-electron vs. frozen-orbital aspect.

    Returns
    -------
    (rpv, pv, expected)
        Of all the QCVariables `pv` that should be available, returns tuple of
        whether `expected` and what key `rpv` in the reference `pv` should match.

"""


def contractual_current(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str
) -> Tuple[str, str, bool]:
    """Given the target method, returns the CURRENT QCVariables that should be produced.

    {_contractual_docstring}
    """
    contractual_qcvars = [
        ("HF TOTAL ENERGY", "SCF TOTAL ENERGY"),
        ("HF TOTAL ENERGY", "CURRENT REFERENCE ENERGY"),
        (f"{method.upper()} CORRELATION ENERGY", "CURRENT CORRELATION ENERGY"),
        (f"{method.upper()} TOTAL ENERGY", "CURRENT ENERGY"),
    ]
    if driver == "gradient":
        contractual_qcvars.append((f"{method.upper()} TOTAL GRADIENT", "CURRENT GRADIENT"))
    elif driver == "hessian":
        # contractual_qcvars.append((f"{method.upper()} TOTAL GRADIENT", "CURRENT GRADIENT"))
        contractual_qcvars.append((f"{method.upper()} TOTAL HESSIAN", "CURRENT HESSIAN"))

    for rpv, pv in contractual_qcvars:
        expected = True
        if method == "hf" and rpv == f"{method.upper()} CORRELATION ENERGY":
            expected = False

        yield (rpv, pv, expected)


def contractual_hf(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal HF should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """

    contractual_qcvars = [
        ("HF TOTAL ENERGY", "HF TOTAL ENERGY"),
        ("HF TOTAL ENERGY", "SCF TOTAL ENERGY"),
    ]
    if driver == "gradient" and method == "hf":
        contractual_qcvars.append(("HF TOTAL GRADIENT", "HF TOTAL GRADIENT"))
        # contractual_qcvars.append(("HF TOTAL GRADIENT", "SCF TOTAL GRADIENT"))
    elif driver == "hessian" and method == "hf":
        # contractual_qcvars.append(("HF TOTAL GRADIENT", "HF TOTAL GRADIENT"))
        contractual_qcvars.append(("HF TOTAL HESSIAN", "HF TOTAL HESSIAN"))
        # contractual_qcvars.append(("HF TOTAL GRADIENT", "SCF TOTAL GRADIENT"))

    for rpv, pv in contractual_qcvars:
        yield (rpv, pv, True)


def contractual_mp2(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal MP2 should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "MP2 CORRELATION ENERGY",
        "MP2 TOTAL ENERGY",
        "MP2 SAME-SPIN CORRELATION ENERGY",
        "MP2 SINGLES ENERGY",
        "MP2 DOUBLES ENERGY",
        "MP2 OPPOSITE-SPIN CORRELATION ENERGY",
    ]
    if driver == "gradient" and method == "mp2":
        contractual_qcvars.append("MP2 TOTAL GRADIENT")
    elif driver == "hessian" and method == "mp2":
        # contractual_qcvars.append("MP2 TOTAL GRADIENT")
        contractual_qcvars.append("MP2 TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True
        if (
            (
                (
                    (qc_module == "cfour" and reference == "rohf" and method == "mp2" and driver == "hessian")
                    or (
                        qc_module in ["gamess-serial", "gamess-ddi"]
                        and reference in ["uhf", "rohf"]
                        and method == "mp2"
                    )
                    or (
                        qc_module in ["gamess-serial", "gamess-ims"]
                        and reference == "rhf"
                        and method == "mp2"
                        and driver in ["gradient", "hessian"]
                    )
                    or (
                        qc_module == "gamess"
                        and reference in ["rhf"]
                        and method in ["lccd", "ccd", "ccsd", "ccsd+t(ccsd)", "ccsd(t)"]
                    )
                    or (qc_module == "nwchem-tce" and method in ["mp2", "mp3", "mp4"])
                    or (
                        qc_module == "nwchem-cc"
                        and reference in ["rhf"]
                        and method in ["ccsd", "ccsd+t(ccsd)", "ccsd(t)"]
                    )
                    or (qc_module == "nwchem-directmp2" and reference == "rhf" and method == "mp2")
                    or (
                        qc_module == "psi4-occ"
                        and reference == "rhf"
                        and corl_type in ["df", "cd"]
                        and method in ["mp2", "mp2.5", "mp3", "lccd", "ccsd", "ccsd(t)"]
                    )
                )
                and pv in ["MP2 SAME-SPIN CORRELATION ENERGY", "MP2 OPPOSITE-SPIN CORRELATION ENERGY"]
            )
            or (
                (
                    (qc_module == "psi4-detci" and method in ["mp2", "mp3"])
                    or (
                        qc_module == "qchem" and method == "mp2"
                    )  # for structured -- can probably get these from parsing
                )
                and pv
                in [
                    "MP2 SAME-SPIN CORRELATION ENERGY",
                    "MP2 OPPOSITE-SPIN CORRELATION ENERGY",
                    "MP2 SINGLES ENERGY",
                    "MP2 DOUBLES ENERGY",
                ]
            )
            or (
                ((qc_module == "psi4-occ" and reference == "rohf" and method in ["olccd"]))
                and pv
                in [
                    "MP2 CORRELATION ENERGY",
                    "MP2 TOTAL ENERGY",
                    "MP2 SINGLES ENERGY",
                ]
            )
            or (
                (
                    (qc_module == "psi4-ccenergy" and reference == "rohf" and method == "ccsd")
                    or (
                        qc_module == "nwchem-tce"
                        and method in ["qcisd", "lccd", "lccsd", "ccd", "ccsd", "ccsd+t(ccsd)", "ccsd(t)", "ccsdt"]
                    )
                    or (qc_module == "gamess" and reference == "rohf" and method == "ccsd")
                    or (
                        qc_module.startswith("cfour")
                        and reference == "rohf"
                        and fcae == "fc"
                        and method in ["lccsd", "ccd", "ccsd", "ccsd(t)", "ccsdt"]
                    )  # this is a cop out as c4 perfectly able to produce good rohf mp2 but not with same orbitals as ref definition on ccsd
                )
                and pv
                in [
                    "MP2 CORRELATION ENERGY",
                    "MP2 TOTAL ENERGY",
                    "MP2 SAME-SPIN CORRELATION ENERGY",
                    "MP2 OPPOSITE-SPIN CORRELATION ENERGY",
                    "MP2 SINGLES ENERGY",
                    "MP2 DOUBLES ENERGY",
                ]
            )
        ):
            expected = False

        yield (pv, pv, expected)


#        # TODO check CUSTOM SCS-MP2 _absent_


def contractual_mp2p5(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str
) -> Tuple[str, str, bool]:
    f"""Of the list of QCVariables an ideal MP2.5 should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "MP2.5 CORRELATION ENERGY",
        "MP2.5 TOTAL ENERGY",
        "MP2.5 SAME-SPIN CORRELATION ENERGY",
        "MP2.5 SINGLES ENERGY",
        "MP2.5 DOUBLES ENERGY",
        "MP2.5 OPPOSITE-SPIN CORRELATION ENERGY",
    ]
    if driver == "gradient" and method == "mp2.5":
        contractual_qcvars.append("MP2.5 TOTAL GRADIENT")
    elif driver == "hessian" and method == "mp2.5":
        # contractual_qcvars.append("MP2.5 TOTAL GRADIENT")
        contractual_qcvars.append("MP2.5 TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True
        if (
            (
                (qc_module.startswith("cfour") and method in ["mp3", "mp4(sdq)", "mp4"])
                or (
                    qc_module == "psi4-occ"
                    and reference == "rhf"
                    and corl_type in ["df", "cd"]
                    and method in ["mp2.5", "mp3"]
                )
                or (qc_module.startswith("nwchem") and method in ["mp3", "mp4"])
            )
            and pv in ["MP2.5 SAME-SPIN CORRELATION ENERGY", "MP2.5 OPPOSITE-SPIN CORRELATION ENERGY"]
        ) or (
            (qc_module == "psi4-detci" and method in ["mp3"])
            and pv
            in [
                "MP2.5 CORRELATION ENERGY",
                "MP2.5 TOTAL ENERGY",
                "MP2.5 SAME-SPIN CORRELATION ENERGY",
                "MP2.5 OPPOSITE-SPIN CORRELATION ENERGY",
                "MP2.5 SINGLES ENERGY",
                "MP2.5 DOUBLES ENERGY",
            ]
        ):
            expected = False

        yield (pv, pv, expected)


def contractual_mp3(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal MP3 should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "MP3 CORRELATION ENERGY",
        "MP3 TOTAL ENERGY",
        "MP3 SAME-SPIN CORRELATION ENERGY",
        "MP3 SINGLES ENERGY",
        "MP3 DOUBLES ENERGY",
        "MP3 OPPOSITE-SPIN CORRELATION ENERGY",
    ]
    if driver == "gradient" and method == "mp3":
        contractual_qcvars.append("MP3 TOTAL GRADIENT")
    elif driver == "hessian" and method == "mp3":
        # contractual_qcvars.append("MP3 TOTAL GRADIENT")
        contractual_qcvars.append("MP3 TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True
        if (
            (
                (qc_module.startswith("cfour") and method in ["mp3", "mp4(sdq)", "mp4"])
                or (
                    qc_module == "psi4-occ"
                    and reference == "rhf"
                    and corl_type in ["df", "cd"]
                    and method in ["mp2.5", "mp3"]
                )
                or (qc_module.startswith("nwchem") and method in ["mp3", "mp4"])
            )
            and pv in ["MP3 SAME-SPIN CORRELATION ENERGY", "MP3 OPPOSITE-SPIN CORRELATION ENERGY"]
        ) or (
            ((qc_module == "psi4-detci" and method == "mp3"))
            and pv
            in [
                "MP3 SAME-SPIN CORRELATION ENERGY",
                "MP3 OPPOSITE-SPIN CORRELATION ENERGY",
                "MP3 SINGLES ENERGY",
                "MP3 DOUBLES ENERGY",
            ]
        ):
            expected = False

        yield (pv, pv, expected)


def contractual_mp4_prsdq_pr(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal MP4(SDQ) should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "MP4(SDQ) CORRELATION ENERGY",
        "MP4(SDQ) TOTAL ENERGY",
    ]
    if driver == "gradient" and method == "mp4(sdq)":
        contractual_qcvars.append("MP4(SDQ) TOTAL GRADIENT")
    elif driver == "hessian" and method == "mp4(sdq)":
        # contractual_qcvars.append("MP4(SDQ) TOTAL GRADIENT")
        contractual_qcvars.append("MP4(SDQ) TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True
        if (qc_module.startswith("nwchem") and method == "mp4") and pv in [
            "MP4(SDQ) TOTAL ENERGY",
            "MP4(SDQ) CORRELATION ENERGY",
        ]:
            expected = False

        yield (pv, pv, expected)


def contractual_mp4(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal MP4 should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "MP4(T) CORRECTION ENERGY",
        "MP4 CORRELATION ENERGY",
        "MP4 CORRECTION ENERGY",
        "MP4 TOTAL ENERGY",
    ]
    if driver == "gradient" and method == "mp4":
        contractual_qcvars.append("MP4 TOTAL GRADIENT")
    elif driver == "hessian" and method == "mp4":
        # contractual_qcvars.append("MP4 TOTAL GRADIENT")
        contractual_qcvars.append("MP4 TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True
        if (qc_module.startswith("nwchem") and method == "mp4") and pv in ["MP4(T) CORRECTION ENERGY"]:
            expected = False

        yield (pv, pv, expected)


def contractual_cisd(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str
) -> Tuple[str, str, bool]:
    f"""Of the list of QCVariables an ideal CISD should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "CISD CORRELATION ENERGY",
        "CISD TOTAL ENERGY",
    ]
    if driver == "gradient" and method == "cisd":
        contractual_qcvars.append("CISD TOTAL GRADIENT")
    elif driver == "hessian" and method == "cisd":
        # contractual_qcvars.append("CISD TOTAL GRADIENT")
        contractual_qcvars.append("CISD TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True

        yield (pv, pv, expected)


def contractual_qcisd(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str
) -> Tuple[str, str, bool]:
    f"""Of the list of QCVariables an ideal QCISD should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "QCISD CORRELATION ENERGY",
        "QCISD TOTAL ENERGY",
    ]
    if driver == "gradient" and method == "qcisd":
        contractual_qcvars.append("QCISD TOTAL GRADIENT")
    elif driver == "hessian" and method == "qcisd":
        # contractual_qcvars.append("QCISD TOTAL GRADIENT")
        contractual_qcvars.append("QCISD TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True

        yield (pv, pv, expected)


def contractual_qcisd_prt_pr(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str
) -> Tuple[str, str, bool]:
    f"""Of the list of QCVariables an ideal QCISD(T) should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "QCISD(T) CORRECTION ENERGY",
        "QCISD(T) CORRELATION ENERGY",
        "QCISD(T) TOTAL ENERGY",
    ]
    if driver == "gradient" and method == "qcisd(t)":
        contractual_qcvars.append("QCISD(T) TOTAL GRADIENT")
    elif driver == "hessian" and method == "qcisd(t)":
        # contractual_qcvars.append("QCISD(T) TOTAL GRADIENT")
        contractual_qcvars.append("QCISD(T) TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True

        yield (pv, pv, expected)


def contractual_fci(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str
) -> Tuple[str, str, bool]:
    f"""Of the list of QCVariables an ideal FCI should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "FCI CORRELATION ENERGY",
        "FCI TOTAL ENERGY",
    ]
    if driver == "gradient" and method == "fci":
        contractual_qcvars.append("FCI TOTAL GRADIENT")
    elif driver == "hessian" and method == "fci":
        # contractual_qcvars.append("FCI TOTAL GRADIENT")
        contractual_qcvars.append("FCI TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True

        yield (pv, pv, expected)


def contractual_lccd(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str
) -> Tuple[str, str, bool]:
    f"""Of the list of QCVariables an ideal LCCD should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "LCCD CORRELATION ENERGY",
        "LCCD TOTAL ENERGY",
        "LCCD SAME-SPIN CORRELATION ENERGY",
        "LCCD SINGLES ENERGY",
        "LCCD DOUBLES ENERGY",
        "LCCD OPPOSITE-SPIN CORRELATION ENERGY",
    ]
    if driver == "gradient" and method == "lccd":
        contractual_qcvars.append("LCCD TOTAL GRADIENT")
    elif driver == "hessian" and method == "lccd":
        # contractual_qcvars.append("LCCD TOTAL GRADIENT")
        contractual_qcvars.append("LCCD TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True
        if (
            (
                (qc_module == "psi4-occ" and reference == "rhf" and corl_type in ["df", "cd"] and method == "lccd")
                or (qc_module == "cfour-ncc" and reference in ["rhf"] and method == "lccd")
                or (qc_module == "nwchem-tce" and reference in ["rhf", "uhf"] and method == "lccd")
                or (qc_module == "gamess" and reference in ["rhf"] and method == "lccd")
            )
            and pv in ["LCCD SAME-SPIN CORRELATION ENERGY", "LCCD OPPOSITE-SPIN CORRELATION ENERGY"]
        ) or (
            (qc_module == "nwchem-tce" and reference in ["rohf"] and method in ["lccd"])
            and pv
            in [
                "LCCD SAME-SPIN CORRELATION ENERGY",
                "LCCD OPPOSITE-SPIN CORRELATION ENERGY",
                "LCCD SINGLES ENERGY",
                "LCCD DOUBLES ENERGY",
            ]
        ):
            expected = False

        yield (pv, pv, expected)


def contractual_lccsd(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str
) -> Tuple[str, str, bool]:
    f"""Of the list of QCVariables an ideal LCCSD should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "LCCSD CORRELATION ENERGY",
        "LCCSD TOTAL ENERGY",
        "LCCSD SAME-SPIN CORRELATION ENERGY",
        "LCCSD SINGLES ENERGY",
        "LCCSD DOUBLES ENERGY",
        "LCCSD OPPOSITE-SPIN CORRELATION ENERGY",
    ]
    if driver == "gradient" and method == "lccsd":
        contractual_qcvars.append("LCCSD TOTAL GRADIENT")
    elif driver == "hessian" and method == "lccsd":
        # contractual_qcvars.append("LCCSD TOTAL GRADIENT")
        contractual_qcvars.append("LCCSD TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True
        if (
            (qc_module == "cfour-ncc" and reference in ["rhf"] and method == "lccsd")
            or (qc_module == "nwchem-tce" and reference in ["rhf", "uhf"] and method == "lccsd")
        ) and pv in ["LCCSD SAME-SPIN CORRELATION ENERGY", "LCCSD OPPOSITE-SPIN CORRELATION ENERGY"]:
            expected = False

        yield (pv, pv, expected)


def contractual_ccd(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal CCD should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "CCD CORRELATION ENERGY",
        "CCD TOTAL ENERGY",
        "CCD SAME-SPIN CORRELATION ENERGY",
        "CCD SINGLES ENERGY",
        "CCD DOUBLES ENERGY",
        "CCD OPPOSITE-SPIN CORRELATION ENERGY",
    ]
    if driver == "gradient" and method == "ccd":
        contractual_qcvars.append("CCD TOTAL GRADIENT")
    elif driver == "hessian" and method == "ccd":
        # contractual_qcvars.append("CCD TOTAL GRADIENT")
        contractual_qcvars.append("CCD TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True
        if (
            (
                (qc_module == "cfour-ecc" and reference in ["rhf"] and method == "ccd")
                or (qc_module == "cfour-ncc" and reference in ["rhf"] and method == "ccd")
                or (qc_module == "nwchem-tce" and reference in ["rhf", "uhf"] and method == "ccd")
                or (qc_module == "gamess" and reference in ["rhf"] and method == "ccd")
            )
            and pv in ["CCD SAME-SPIN CORRELATION ENERGY", "CCD OPPOSITE-SPIN CORRELATION ENERGY"]
        ) or (
            (
                (qc_module == "cfour-vcc" and reference in ["rohf"] and method in ["ccd"])
                or (qc_module == "nwchem-tce" and reference in ["rohf"] and method in ["ccd"])
            )
            and pv
            in [
                "CCD SAME-SPIN CORRELATION ENERGY",
                "CCD SINGLES ENERGY",
                "CCD DOUBLES ENERGY",
                "CCD OPPOSITE-SPIN CORRELATION ENERGY",
            ]
        ):
            expected = False

        yield (pv, pv, expected)


def contractual_ccsd(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal CCSD should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "CCSD CORRELATION ENERGY",
        "CCSD TOTAL ENERGY",
        "CCSD SAME-SPIN CORRELATION ENERGY",
        "CCSD SINGLES ENERGY",
        "CCSD DOUBLES ENERGY",
        "CCSD OPPOSITE-SPIN CORRELATION ENERGY",
    ]
    if driver == "gradient" and method == "ccsd":
        contractual_qcvars.append("CCSD TOTAL GRADIENT")
    elif driver == "hessian" and method == "ccsd":
        # contractual_qcvars.append("CCSD TOTAL GRADIENT")
        contractual_qcvars.append("CCSD TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True
        if (
            (
                (
                    (qc_module == "gamess" and reference == "rhf" and method in ["ccsd", "ccsd+t(ccsd)", "ccsd(t)"])
                    or (
                        qc_module == "nwchem-tce"
                        and reference in ["rhf", "uhf"]
                        and method in ["ccsd", "ccsd+t(ccsd)", "ccsd(t)"]
                    )
                    or (
                        qc_module in ["cfour-ncc", "cfour-ecc"]
                        and reference in ["rhf"]
                        and method in ["ccsd", "ccsd+t(ccsd)", "ccsd(t)", "a-ccsd(t)"]
                    )
                    or (
                        qc_module == "psi4-occ"
                        and reference == "rhf"
                        and corl_type in ["df", "cd"]
                        and method in ["ccsd", "ccsd(t)"]
                    )
                    or (qc_module in ["cfour-vcc"] and reference in ["rhf", "uhf"] and method in ["ccsd+t(ccsd)"])
                )
                and pv in ["CCSD SAME-SPIN CORRELATION ENERGY", "CCSD OPPOSITE-SPIN CORRELATION ENERGY"]
            )
            or (
                (qc_module == "cfour-vcc" and reference in ["rohf"] and method in ["ccsd", "ccsd(t)"])
                and pv
                in [
                    "CCSD SAME-SPIN CORRELATION ENERGY",
                    "CCSD SINGLES ENERGY",
                    "CCSD DOUBLES ENERGY",
                ]
            )
            or (
                (qc_module == "cfour-ecc" and reference in ["rohf"] and method in ["ccsd", "ccsd(t)"])
                and pv
                in [
                    "CCSD OPPOSITE-SPIN CORRELATION ENERGY",
                    "CCSD SINGLES ENERGY",
                    "CCSD DOUBLES ENERGY",
                ]
            )
            or (
                (
                    (qc_module == "gamess" and reference in ["rohf"] and method == "ccsd")
                    or (qc_module == "nwchem-tce" and reference in ["rohf"] and method in ["ccsd", "ccsd(t)"])
                )
                and pv
                in [
                    "CCSD SAME-SPIN CORRELATION ENERGY",
                    "CCSD OPPOSITE-SPIN CORRELATION ENERGY",
                    "CCSD SINGLES ENERGY",
                    "CCSD DOUBLES ENERGY",
                ]
            )
            or (
                (False)
                and pv
                in [
                    "CCSD CORRELATION ENERGY",
                    "CCSD TOTAL ENERGY",
                    "CCSD SAME-SPIN CORRELATION ENERGY",
                    "CCSD OPPOSITE-SPIN CORRELATION ENERGY",
                    "CCSD SINGLES ENERGY",
                    "CCSD DOUBLES ENERGY",
                ]
            )
        ):
            expected = False

        yield (pv, pv, expected)

    # TODO check CUSTOM SCS-CCSD _absent_


def contractual_ccsdpt_prccsd_pr(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal CCSD+T(CCSD) (aka CCSD[T]) should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "T(CCSD) CORRECTION ENERGY",
        "CCSD+T(CCSD) CORRELATION ENERGY",
        "CCSD+T(CCSD) TOTAL ENERGY",
    ]
    if driver == "gradient" and method == "ccsd+t(ccsd)":
        contractual_qcvars.append("CCSD+T(CCSD) TOTAL GRADIENT")
    elif driver == "hessian" and method == "ccsd+t(ccsd)":
        # contractual_qcvars.append("CCSD+T(CCSD) TOTAL GRADIENT")
        contractual_qcvars.append("CCSD+T(CCSD) TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True
        if False:
            expected = False

        yield (pv, pv, expected)


def contractual_ccsd_prt_pr(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal CCSD(T) should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "(T) CORRECTION ENERGY",
        "CCSD(T) CORRELATION ENERGY",
        "CCSD(T) TOTAL ENERGY",
    ]
    if driver == "gradient" and method == "ccsd(t)":
        contractual_qcvars.append("CCSD(T) TOTAL GRADIENT")
    elif driver == "hessian" and method == "ccsd(t)":
        # contractual_qcvars.append("CCSD(T) TOTAL GRADIENT")
        contractual_qcvars.append("CCSD(T) TOTAL HESSIAN")

    for pv in contractual_qcvars:
        # print("WW", qc_module, driver, reference, method, corl_type, fcae, pv)
        expected = True
        if False:
            expected = False

        yield (pv, pv, expected)


def contractual_accsd_prt_pr(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal A-CCSD(T) (aka Lambda-CCSD(T), aka CCSD(aT) should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "A-(T) CORRECTION ENERGY",
        "A-CCSD(T) CORRELATION ENERGY",
        "A-CCSD(T) TOTAL ENERGY",
    ]
    if driver == "gradient" and method == "a-ccsd(t)":
        contractual_qcvars.append("A-CCSD(T) TOTAL GRADIENT")
    elif driver == "hessian" and method == "a-ccsd(t)":
        # contractual_qcvars.append("A-CCSD(T) TOTAL GRADIENT")
        contractual_qcvars.append("A-CCSD(T) TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True
        if False:
            expected = False

        yield (pv, pv, expected)


def contractual_ccsdt(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str
) -> Tuple[str, str, bool]:
    f"""Of the list of QCVariables an ideal CCSDT should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "CCSDT CORRELATION ENERGY",
        "CCSDT TOTAL ENERGY",
    ]
    if driver == "gradient" and method == "ccsdt":
        contractual_qcvars.append("CCSDT TOTAL GRADIENT")
    elif driver == "hessian" and method == "ccsdt":
        # contractual_qcvars.append("CCSDT TOTAL GRADIENT")
        contractual_qcvars.append("CCSDT TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True

        yield (pv, pv, expected)


def contractual_ccsdt1a(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str
) -> Tuple[str, str, bool]:
    f"""Of the list of QCVariables an ideal CCSDT-1A should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "CCSDT-1A CORRELATION ENERGY",
        "CCSDT-1A TOTAL ENERGY",
    ]
    if driver == "gradient" and method == "ccsdt-1a":
        contractual_qcvars.append("CCSDT-1A TOTAL GRADIENT")
    elif driver == "hessian" and method == "ccsdt-1a":
        # contractual_qcvars.append("CCSDT-1A TOTAL GRADIENT")
        contractual_qcvars.append("CCSDT-1A TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True

        yield (pv, pv, expected)


def contractual_ccsdt1b(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str
) -> Tuple[str, str, bool]:
    f"""Of the list of QCVariables an ideal CCSDT-1B should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "CCSDT-1B CORRELATION ENERGY",
        "CCSDT-1B TOTAL ENERGY",
    ]
    if driver == "gradient" and method == "ccsdt-1b":
        contractual_qcvars.append("CCSDT-1B TOTAL GRADIENT")
    elif driver == "hessian" and method == "ccsdt-1b":
        # contractual_qcvars.append("CCSDT-1B TOTAL GRADIENT")
        contractual_qcvars.append("CCSDT-1B TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True

        yield (pv, pv, expected)


def contractual_ccsdt2(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str
) -> Tuple[str, str, bool]:
    f"""Of the list of QCVariables an ideal CCSDT-2 should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "CCSDT-2 CORRELATION ENERGY",
        "CCSDT-2 TOTAL ENERGY",
    ]
    if driver == "gradient" and method == "ccsdt-2":
        contractual_qcvars.append("CCSDT-2 TOTAL GRADIENT")
    elif driver == "hessian" and method == "ccsdt-2":
        # contractual_qcvars.append("CCSDT-2 TOTAL GRADIENT")
        contractual_qcvars.append("CCSDT-2 TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True

        yield (pv, pv, expected)


def contractual_ccsdt3(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str
) -> Tuple[str, str, bool]:
    f"""Of the list of QCVariables an ideal CCSDT-3 should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "CCSDT-3 CORRELATION ENERGY",
        "CCSDT-3 TOTAL ENERGY",
    ]
    if driver == "gradient" and method == "ccsdt-3":
        contractual_qcvars.append("CCSDT-3 TOTAL GRADIENT")
    elif driver == "hessian" and method == "ccsdt-3":
        # contractual_qcvars.append("CCSDT-3 TOTAL GRADIENT")
        contractual_qcvars.append("CCSDT-3 TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True

        yield (pv, pv, expected)


def contractual_ccsdt_prq_pr(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal CCSDT(Q) should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "(Q) CORRECTION ENERGY",
        "CCSDT(Q) CORRELATION ENERGY",
        "CCSDT(Q) TOTAL ENERGY",
    ]
    if driver == "gradient" and method == "ccsdt(q)":
        contractual_qcvars.append("CCSDT(Q) TOTAL GRADIENT")
    elif driver == "hessian" and method == "ccsdt(q)":
        # contractual_qcvars.append("CCSDT(Q) TOTAL GRADIENT")
        contractual_qcvars.append("CCSDT(Q) TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True

        yield (pv, pv, expected)


def contractual_ccsdtq(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal CCSDTQ should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "CCSDTQ CORRELATION ENERGY",
        "CCSDTQ TOTAL ENERGY",
    ]
    if driver == "gradient" and method == "ccsdtq":
        contractual_qcvars.append("CCSDTQ TOTAL GRADIENT")
    elif driver == "hessian" and method == "ccsdtq":
        # contractual_qcvars.append("CCSDTQ TOTAL GRADIENT")
        contractual_qcvars.append("CCSDTQ TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True

        yield (pv, pv, expected)


def contractual_olccd(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal OLCCD should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "OLCCD CORRELATION ENERGY",
        "OLCCD TOTAL ENERGY",
        "OLCCD REFERENCE CORRECTION ENERGY",
        "OLCCD SAME-SPIN CORRELATION ENERGY",
        "OLCCD OPPOSITE-SPIN CORRELATION ENERGY",
    ]
    if driver == "gradient" and method == "olccd":
        contractual_qcvars.append("OLCCD TOTAL GRADIENT")
    elif driver == "hessian" and method == "olccd":
        # contractual_qcvars.append("OLCCD TOTAL GRADIENT")
        contractual_qcvars.append("OLCCD TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True

        yield (pv, pv, expected)


def contractual_dft_current(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str
) -> Tuple[str, str, bool]:
    """Given the target DFT method, returns the CURRENT QCVariables that should be produced.

    {_contractual_docstring}
    """
    contractual_qcvars = [
        (f"{method.upper()} TOTAL ENERGY", "DFT TOTAL ENERGY"),
        (f"{method.upper()} TOTAL ENERGY", "CURRENT REFERENCE ENERGY"),
        (f"{method.upper()} TOTAL ENERGY", "CURRENT ENERGY"),
    ]
    if driver == "gradient":
        contractual_qcvars.append((f"{method.upper()} TOTAL GRADIENT", "CURRENT GRADIENT"))
    elif driver == "hessian":
        # contractual_qcvars.append((f"{method.upper()} TOTAL GRADIENT", "CURRENT GRADIENT"))
        contractual_qcvars.append((f"{method.upper()} TOTAL HESSIAN", "CURRENT HESSIAN"))

    for rpv, pv in contractual_qcvars:
        expected = True
        if method == "hf" and rpv == f"{method.upper()} CORRELATION ENERGY":
            expected = False

        yield (rpv, pv, expected)
