import numpy as np
from qcengine.programs.tests.standard_suite_ref import answer_hash, _std_suite, _std_generics


_std_suite_psi4_extension = [
    # <<<  CD-AE-CD  >>>
    {
        "meta": {
            "system": "hf",
            "basis": "cc-pvdz",
            "scf_type": "cd",
            "reference": "rhf",
            "fcae": "ae",
            "corl_type": "cd",
        },
        "data": {
            "MP3 TOTAL GRADIENT": np.array(
                # dfocc findif-5
                [ 0., 0., -0.000926981449, 0., 0., 0.000926981449]
            ).reshape((-1, 3)),
            "LCCD TOTAL GRADIENT": np.array(
                # dfocc findif-5
                [ 0., 0., 0.002193849073, 0., 0., -0.002193849073]
            ).reshape((-1, 3)),
        },
    },
    # <<<  CD-FC-CD  >>>
    {
        "meta": {
            "system": "hf", 
            "basis": "cc-pvdz",
            "scf_type": "cd", 
            "reference": "rhf",
            "fcae": "fc", 
            "corl_type": "cd", 
        },    
        "data": {
            "MP3 TOTAL GRADIENT": np.array(
                # dfocc findif-5 fc cd+cd
                [ 0., 0., -0.000588974421, 0., 0., 0.000588974421]
            ).reshape((-1, 3)),
            "LCCD TOTAL GRADIENT": np.array(
                # dfocc findif-5 fc cd+cd
                [ 0., 0., 0.002525704147, 0., 0., -0.002525704147]
            ).reshape((-1, 3)),
        },
    },
]


for calc1 in _std_suite_psi4_extension:
    metahash1 = answer_hash(**calc1["meta"])
    for calc0 in _std_suite:
        metahash0 = answer_hash(**calc0["meta"])
        if metahash0 == metahash1:
            calc0["data"].update(calc1["data"])
            break


# TODO: this should become a fn in qcng so doesn't have to be copied here
for calc in _std_suite:
    if calc["data"]:
        if "MP2 CORRELATION ENERGY" in calc["data"]:
            calc["data"]["MP2 TOTAL ENERGY"] = calc["data"]["MP2 CORRELATION ENERGY"] + calc["data"]["HF TOTAL ENERGY"]
            if "MP2 SINGLES ENERGY" in calc["data"]:
                calc["data"]["MP2 DOUBLES ENERGY"] = (
                    calc["data"]["MP2 CORRELATION ENERGY"] - calc["data"]["MP2 SINGLES ENERGY"]
                )
                if "MP2 SAME-SPIN CORRELATION ENERGY" in calc["data"]:
                    calc["data"]["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = (
                        calc["data"]["MP2 CORRELATION ENERGY"]
                        - calc["data"]["MP2 SAME-SPIN CORRELATION ENERGY"]
                        - calc["data"]["MP2 SINGLES ENERGY"]
                    )
                    calc["data"]["SCS-MP2 CORRELATION ENERGY"] = (
                        (1 / 3) * calc["data"]["MP2 SAME-SPIN CORRELATION ENERGY"]
                        + (6 / 5) * calc["data"]["MP2 OPPOSITE-SPIN CORRELATION ENERGY"]
                        + calc["data"]["MP2 SINGLES ENERGY"]
                    )
                    calc["data"]["SCS-MP2 TOTAL ENERGY"] = (
                        calc["data"]["SCS-MP2 CORRELATION ENERGY"] + calc["data"]["HF TOTAL ENERGY"]
                    )

        if "MP3 CORRELATION ENERGY" in calc["data"]:
            calc["data"]["MP3 TOTAL ENERGY"] = calc["data"]["MP3 CORRELATION ENERGY"] + calc["data"]["HF TOTAL ENERGY"]
            if "MP2 CORRELATION ENERGY" in calc["data"]:
                calc["data"]["MP2.5 CORRELATION ENERGY"] = 0.5 * (
                    calc["data"]["MP3 CORRELATION ENERGY"] + calc["data"]["MP2 CORRELATION ENERGY"]
                )
                calc["data"]["MP2.5 TOTAL ENERGY"] = (
                    calc["data"]["MP2.5 CORRELATION ENERGY"] + calc["data"]["HF TOTAL ENERGY"]
                )

            if "MP3 SINGLES ENERGY" in calc["data"]:
                calc["data"]["MP3 DOUBLES ENERGY"] = (
                    calc["data"]["MP3 CORRELATION ENERGY"] - calc["data"]["MP3 SINGLES ENERGY"]
                )
                if "MP2 SINGLES ENERGY" in calc["data"]:
                    calc["data"]["MP2.5 SINGLES ENERGY"] = 0.5 * (
                        calc["data"]["MP3 SINGLES ENERGY"] + calc["data"]["MP2 SINGLES ENERGY"]
                    )
                    calc["data"]["MP2.5 DOUBLES ENERGY"] = (
                        calc["data"]["MP2.5 CORRELATION ENERGY"] - calc["data"]["MP2.5 SINGLES ENERGY"]
                    )
                if "MP3 SAME-SPIN CORRELATION ENERGY" in calc["data"]:
                    calc["data"]["MP3 OPPOSITE-SPIN CORRELATION ENERGY"] = (
                        calc["data"]["MP3 CORRELATION ENERGY"]
                        - calc["data"]["MP3 SAME-SPIN CORRELATION ENERGY"]
                        - calc["data"]["MP3 SINGLES ENERGY"]
                    )
                    if "MP2 SAME-SPIN CORRELATION ENERGY" in calc["data"]:
                        calc["data"]["MP2.5 SAME-SPIN CORRELATION ENERGY"] = 0.5 * (
                            calc["data"]["MP3 SAME-SPIN CORRELATION ENERGY"]
                            + calc["data"]["MP2 SAME-SPIN CORRELATION ENERGY"]
                        )
                        calc["data"]["MP2.5 OPPOSITE-SPIN CORRELATION ENERGY"] = (
                            calc["data"]["MP2.5 CORRELATION ENERGY"]
                            - calc["data"]["MP2.5 SAME-SPIN CORRELATION ENERGY"]
                            - calc["data"]["MP2.5 SINGLES ENERGY"]
                        )

        if (
            "MP3 TOTAL GRADIENT" in calc["data"]
            and "MP2 TOTAL GRADIENT" in calc["data"]
            and "HF TOTAL GRADIENT" in calc["data"]
        ):
            calc["data"]["MP2.5 TOTAL GRADIENT"] = 0.5 * (
                calc["data"]["MP3 TOTAL GRADIENT"] + calc["data"]["MP2 TOTAL GRADIENT"]
            )

        if "MP4(SDQ) CORRELATION ENERGY" in calc["data"]:
            calc["data"]["MP4(SDQ) TOTAL ENERGY"] = (
                calc["data"]["MP4(SDQ) CORRELATION ENERGY"] + calc["data"]["HF TOTAL ENERGY"]
            )
            if "MP4(T) CORRECTION ENERGY" in calc["data"]:
                calc["data"]["MP4 CORRELATION ENERGY"] = (
                    calc["data"]["MP4(SDQ) CORRELATION ENERGY"] + calc["data"]["MP4(T) CORRECTION ENERGY"]
                )
                calc["data"]["MP4 TOTAL ENERGY"] = (
                    calc["data"]["MP4 CORRELATION ENERGY"] + calc["data"]["HF TOTAL ENERGY"]
                )
                if "MP3 CORRELATION ENERGY" in calc["data"]:
                    calc["data"]["MP4 CORRECTION ENERGY"] = (
                        calc["data"]["MP4 CORRELATION ENERGY"] - calc["data"]["MP3 CORRELATION ENERGY"]
                    )

        if "CISD CORRELATION ENERGY" in calc["data"]:
            calc["data"]["CISD TOTAL ENERGY"] = (
                calc["data"]["CISD CORRELATION ENERGY"] + calc["data"]["HF TOTAL ENERGY"]
            )

        if "QCISD CORRELATION ENERGY" in calc["data"]:
            calc["data"]["QCISD TOTAL ENERGY"] = (
                calc["data"]["QCISD CORRELATION ENERGY"] + calc["data"]["HF TOTAL ENERGY"]
            )
            if "QCISD(T) CORRECTION ENERGY" in calc["data"]:
                calc["data"]["QCISD(T) CORRELATION ENERGY"] = (
                    calc["data"]["QCISD CORRELATION ENERGY"] + calc["data"]["QCISD(T) CORRECTION ENERGY"]
                )
                calc["data"]["QCISD(T) TOTAL ENERGY"] = (
                    calc["data"]["QCISD(T) CORRELATION ENERGY"] + calc["data"]["HF TOTAL ENERGY"]
                )

        if "FCI CORRELATION ENERGY" in calc["data"]:
            calc["data"]["FCI TOTAL ENERGY"] = calc["data"]["FCI CORRELATION ENERGY"] + calc["data"]["HF TOTAL ENERGY"]

        if "LCCD CORRELATION ENERGY" in calc["data"]:
            calc["data"]["LCCD TOTAL ENERGY"] = (
                calc["data"]["LCCD CORRELATION ENERGY"] + calc["data"]["HF TOTAL ENERGY"]
            )
            if "LCCD SINGLES ENERGY" in calc["data"]:
                calc["data"]["LCCD DOUBLES ENERGY"] = (
                    calc["data"]["LCCD CORRELATION ENERGY"] - calc["data"]["LCCD SINGLES ENERGY"]
                )
                if "LCCD SAME-SPIN CORRELATION ENERGY" in calc["data"]:
                    calc["data"]["LCCD OPPOSITE-SPIN CORRELATION ENERGY"] = (
                        calc["data"]["LCCD CORRELATION ENERGY"]
                        - calc["data"]["LCCD SAME-SPIN CORRELATION ENERGY"]
                        - calc["data"]["LCCD SINGLES ENERGY"]
                    )

        if "LCCSD CORRELATION ENERGY" in calc["data"]:
            calc["data"]["LCCSD TOTAL ENERGY"] = (
                calc["data"]["LCCSD CORRELATION ENERGY"] + calc["data"]["HF TOTAL ENERGY"]
            )
            if "LCCSD SINGLES ENERGY" in calc["data"]:
                calc["data"]["LCCSD DOUBLES ENERGY"] = (
                    calc["data"]["LCCSD CORRELATION ENERGY"] - calc["data"]["LCCSD SINGLES ENERGY"]
                )
                if "LCCSD SAME-SPIN CORRELATION ENERGY" in calc["data"]:
                    calc["data"]["LCCSD OPPOSITE-SPIN CORRELATION ENERGY"] = (
                        calc["data"]["LCCSD CORRELATION ENERGY"]
                        - calc["data"]["LCCSD SAME-SPIN CORRELATION ENERGY"]
                        - calc["data"]["LCCSD SINGLES ENERGY"]
                    )

        if "CCD CORRELATION ENERGY" in calc["data"]:
            calc["data"]["CCD TOTAL ENERGY"] = calc["data"]["CCD CORRELATION ENERGY"] + calc["data"]["HF TOTAL ENERGY"]
            if "CCD SINGLES ENERGY" in calc["data"]:
                calc["data"]["CCD DOUBLES ENERGY"] = (
                    calc["data"]["CCD CORRELATION ENERGY"] - calc["data"]["CCD SINGLES ENERGY"]
                )
                if "CCD SAME-SPIN CORRELATION ENERGY" in calc["data"]:
                    calc["data"]["CCD OPPOSITE-SPIN CORRELATION ENERGY"] = (
                        calc["data"]["CCD CORRELATION ENERGY"]
                        - calc["data"]["CCD SAME-SPIN CORRELATION ENERGY"]
                        - calc["data"]["CCD SINGLES ENERGY"]
                    )

        if "CCSD CORRELATION ENERGY" in calc["data"]:
            calc["data"]["CCSD TOTAL ENERGY"] = (
                calc["data"]["CCSD CORRELATION ENERGY"] + calc["data"]["HF TOTAL ENERGY"]
            )
            if "CCSD SINGLES ENERGY" in calc["data"]:
                calc["data"]["CCSD DOUBLES ENERGY"] = (
                    calc["data"]["CCSD CORRELATION ENERGY"] - calc["data"]["CCSD SINGLES ENERGY"]
                )
                if "CCSD SAME-SPIN CORRELATION ENERGY" in calc["data"]:
                    calc["data"]["CCSD OPPOSITE-SPIN CORRELATION ENERGY"] = (
                        calc["data"]["CCSD CORRELATION ENERGY"]
                        - calc["data"]["CCSD SAME-SPIN CORRELATION ENERGY"]
                        - calc["data"]["CCSD SINGLES ENERGY"]
                    )

        if "T(CCSD) CORRECTION ENERGY" in calc["data"]:
            calc["data"]["CCSD+T(CCSD) CORRELATION ENERGY"] = (
                calc["data"]["CCSD CORRELATION ENERGY"] + calc["data"]["T(CCSD) CORRECTION ENERGY"]
            )
            calc["data"]["CCSD+T(CCSD) TOTAL ENERGY"] = (
                calc["data"]["CCSD+T(CCSD) CORRELATION ENERGY"] + calc["data"]["HF TOTAL ENERGY"]
            )

        if "(T) CORRECTION ENERGY" in calc["data"]:
            calc["data"]["CCSD(T) CORRELATION ENERGY"] = (
                calc["data"]["CCSD CORRELATION ENERGY"] + calc["data"]["(T) CORRECTION ENERGY"]
            )
            calc["data"]["CCSD(T) TOTAL ENERGY"] = (
                calc["data"]["CCSD(T) CORRELATION ENERGY"] + calc["data"]["HF TOTAL ENERGY"]
            )

        if "A-(T) CORRECTION ENERGY" in calc["data"]:
            calc["data"]["A-CCSD(T) CORRELATION ENERGY"] = (
                calc["data"]["CCSD CORRELATION ENERGY"] + calc["data"]["A-(T) CORRECTION ENERGY"]
            )
            calc["data"]["A-CCSD(T) TOTAL ENERGY"] = (
                calc["data"]["A-CCSD(T) CORRELATION ENERGY"] + calc["data"]["HF TOTAL ENERGY"]
            )

        if "CCSDT-1A CORRELATION ENERGY" in calc["data"]:
            calc["data"]["CCSDT-1A TOTAL ENERGY"] = (
                calc["data"]["CCSDT-1A CORRELATION ENERGY"] + calc["data"]["HF TOTAL ENERGY"]
            )

        if "CCSDT-1B CORRELATION ENERGY" in calc["data"]:
            calc["data"]["CCSDT-1B TOTAL ENERGY"] = (
                calc["data"]["CCSDT-1B CORRELATION ENERGY"] + calc["data"]["HF TOTAL ENERGY"]
            )

        if "CCSDT-2 CORRELATION ENERGY" in calc["data"]:
            calc["data"]["CCSDT-2 TOTAL ENERGY"] = (
                calc["data"]["CCSDT-2 CORRELATION ENERGY"] + calc["data"]["HF TOTAL ENERGY"]
            )

        if "CCSDT-3 CORRELATION ENERGY" in calc["data"]:
            calc["data"]["CCSDT-3 TOTAL ENERGY"] = (
                calc["data"]["CCSDT-3 CORRELATION ENERGY"] + calc["data"]["HF TOTAL ENERGY"]
            )

        if "CCSDT CORRELATION ENERGY" in calc["data"]:
            calc["data"]["CCSDT TOTAL ENERGY"] = (
                calc["data"]["CCSDT CORRELATION ENERGY"] + calc["data"]["HF TOTAL ENERGY"]
            )

        if "(Q) CORRECTION ENERGY" in calc["data"]:
            calc["data"]["CCSDT(Q) CORRELATION ENERGY"] = (
                calc["data"]["CCSDT CORRELATION ENERGY"] + calc["data"]["(Q) CORRECTION ENERGY"]
            )
            calc["data"]["CCSDT(Q) TOTAL ENERGY"] = (
                calc["data"]["CCSDT(Q) CORRELATION ENERGY"] + calc["data"]["HF TOTAL ENERGY"]
            )

        if "CCSDTQ CORRELATION ENERGY" in calc["data"]:
            calc["data"]["CCSDTQ TOTAL ENERGY"] = (
                calc["data"]["CCSDTQ CORRELATION ENERGY"] + calc["data"]["HF TOTAL ENERGY"]
            )

        if "OLCCD CORRELATION ENERGY" in calc["data"]:
            calc["data"]["OLCCD TOTAL ENERGY"] = (
                calc["data"]["OLCCD CORRELATION ENERGY"] + calc["data"]["HF TOTAL ENERGY"]
            )
            calc["data"]["OLCCD OPPOSITE-SPIN CORRELATION ENERGY"] = (
                calc["data"]["OLCCD CORRELATION ENERGY"]
                - calc["data"]["OLCCD REFERENCE CORRECTION ENERGY"]
                - calc["data"]["OLCCD SAME-SPIN CORRELATION ENERGY"]
            )

    calc["data"].update(_std_generics[f"{calc['meta']['system']}_{calc['meta']['basis']}_{calc['meta']['fcae']}"])


std_suite = {answer_hash(**calc["meta"]): calc["data"] for calc in _std_suite}
