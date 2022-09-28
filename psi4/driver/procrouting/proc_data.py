from collections import namedtuple
from typing import Tuple


# A negative order indicates perturbative method
mrcc_methods = {
            'ccsd'          : { 'method': 1, 'order':  2, 'fullname': 'CCSD'         },
            'ccsdt'         : { 'method': 1, 'order':  3, 'fullname': 'CCSDT'        },
            'ccsdtq'        : { 'method': 1, 'order':  4, 'fullname': 'CCSDTQ'       },
            'ccsdtqp'       : { 'method': 1, 'order':  5, 'fullname': 'CCSDTQP'      },
            'ccsdtqph'      : { 'method': 1, 'order':  6, 'fullname': 'CCSDTQPH'     },
            'ccsd(t)'       : { 'method': 3, 'order': -3, 'fullname': 'CCSD(T)'      },
            'ccsdt(q)'      : { 'method': 3, 'order': -4, 'fullname': 'CCSDT(Q)'     },
            'ccsdtq(p)'     : { 'method': 3, 'order': -5, 'fullname': 'CCSDTQ(P)'    },
            'ccsdtqp(h)'    : { 'method': 3, 'order': -6, 'fullname': 'CCSDTQP(H)'   },
            'ccsd(t)_l'     : { 'method': 4, 'order': -3, 'fullname': 'CCSD(T)_L'    },
            'ccsdt(q)_l'    : { 'method': 4, 'order': -4, 'fullname': 'CCSDT(Q)_L'   },
            'ccsdtq(p)_l'   : { 'method': 4, 'order': -5, 'fullname': 'CCSDTQ(P)_L'  },
            'ccsdtqp(h)_l'  : { 'method': 4, 'order': -6, 'fullname': 'CCSDTQP(H)_L' },
            'ccsdt-1a'      : { 'method': 5, 'order':  3, 'fullname': 'CCSDT-1a'     },
            'ccsdtq-1a'     : { 'method': 5, 'order':  4, 'fullname': 'CCSDTQ-1a'    },
            'ccsdtqp-1a'    : { 'method': 5, 'order':  5, 'fullname': 'CCSDTQP-1a'   },
            'ccsdtqph-1a'   : { 'method': 5, 'order':  6, 'fullname': 'CCSDTQPH-1a'  },
            'ccsdt-1b'      : { 'method': 6, 'order':  3, 'fullname': 'CCSDT-1b'     },
            'ccsdtq-1b'     : { 'method': 6, 'order':  4, 'fullname': 'CCSDTQ-1b'    },
            'ccsdtqp-1b'    : { 'method': 6, 'order':  5, 'fullname': 'CCSDTQP-1b'   },
            'ccsdtqph-1b'   : { 'method': 6, 'order':  6, 'fullname': 'CCSDTQPH-1b'  },
            'cc2'           : { 'method': 7, 'order':  2, 'fullname': 'CC2'          },
            'cc3'           : { 'method': 7, 'order':  3, 'fullname': 'CC3'          },
            'cc4'           : { 'method': 7, 'order':  4, 'fullname': 'CC4'          },
            'cc5'           : { 'method': 7, 'order':  5, 'fullname': 'CC5'          },
            'cc6'           : { 'method': 7, 'order':  6, 'fullname': 'CC6'          },
            'ccsdt-3'       : { 'method': 8, 'order':  3, 'fullname': 'CCSDT-3'      },
            'ccsdtq-3'      : { 'method': 8, 'order':  4, 'fullname': 'CCSDTQ-3'     },
            'ccsdtqp-3'     : { 'method': 8, 'order':  5, 'fullname': 'CCSDTQP-3'    },
            'ccsdtqph-3'    : { 'method': 8, 'order':  6, 'fullname': 'CCSDTQPH-3'   }
}  # yapf: disable


method_governing_type_keywords = {
        "hf"           : "scf_type",
        "scf"          : "scf_type",
        "qchf"         : "scf_type",
        "dct"          : "dct_type",
        "mp2"          : "mp2_type",
        "mp2.5"        : "mp_type",  # default more complex than read_options
        "mp3"          : "mp_type",  # default more complex than read_options
        "mp4(sdq)"     : "mp_type",
        "mp4"          : "mp_type",
        "cisd"         : "ci_type",
        "cisdt"        : "ci_type",
        "cisdtq"       : "ci_type",
        "fci"          : "ci_type",
        "detci"        : "ci_type",  # imprecise, but only used to check CONV
        "qcisd"        : "ci_type",
        "qcisd(t)"     : "ci_type",
        "remp2"        : "cc_type",
        "lccd"         : "cc_type",
        "lccsd"        : "cc_type",
        "cepa(0)"      : "cc_type",
        "cepa(1)"      : "cc_type",
        "cepa(3)"      : "cc_type",
        "acpf"         : "cc_type",
        "aqcc"         : "cc_type",
        "ccd"          : "cc_type",
        "bccd"         : "cc_type",
        "cc2"          : "cc_type",
        "ccsd"         : "cc_type",
        "ccsd(t)"      : "cc_type",
        "a-ccsd(t)"    : "cc_type",
        "bccd(t)"      : "cc_type",
        "cc3"          : "cc_type",
        "ccenergy"     : "cc_type",

        "omp2"         : "mp2_type",
        "omp2.5"       : "mp_type",
        "omp3"         : "mp_type",
        "oremp2"       : "cc_type",
        "olccd"        : "cc_type",

        "fno-mp3"      : "mp_type",  # default more complex than read_options
        "fno-mp4(sdq)" : "mp_type",
        "fno-mp4"      : "mp_type",
        "fno-cisd"     : "ci_type",
        "fno-qcisd"    : "ci_type",
        "fno-qcisd(t)" : "ci_type",
        "fno-lccd"     : "cc_type",
        "fno-lccsd"    : "cc_type",
        "fno-cepa(0)"  : "cc_type",
        "fno-cepa(1)"  : "cc_type",
        "fno-cepa(3)"  : "cc_type",
        "fno-acpf"     : "cc_type",
        "fno-aqcc"     : "cc_type",
        "fno-ccsd"     : "cc_type",
        "fno-ccsd(t)"  : "cc_type",

        "dlpno-mp2"    : "mp2_type",

        "ep2"          : "mp2_type",
        "eom-cc2"      : "cc_type",
        "eom-ccsd"     : "cc_type",
        "eom-cc3"      : "cc_type",

        "scs-mp2"           : "mp2_type",
        "scs(n)-mp2"        : "mp2_type",
        "scs-mp2-vdw"       : "mp2_type",
        "sos-mp2"           : "mp2_type",
        "sos-pi-mp2"        : "mp2_type",
        "custom-scs-mp2"    : "mp2_type",
        "scs-dlpno-mp2"     : "mp2_type",
        "scs-omp2"          : "mp2_type",
        "sos-omp2"          : "mp2_type",
        "custom-scs-omp2"   : "mp2_type",
        "custom-scs-mp2.5"  : "mp_type",  # default more complex than read_options
        "custom-scs-omp2.5" : "mp_type",
        "scs-mp3"           : "mp_type",  # default more complex than read_options
        "custom-scs-mp3"    : "mp_type",  # default more complex than read_options
        "scs-omp3"          : "mp_type",
        "sos-omp3"          : "mp_type",
        "custom-scs-omp3"   : "mp_type",
        "custom-scs-lccd"   : "cc_type",
        "custom-scs-olccd"  : "cc_type",
}  # yapf: disable

# not adding to lookup until needed:
# * sapt, until LAB understands its types better
# * adc, until builtin is removed and proc edited
# * efp, no type in its own right
# * multiref/cas: mcscf/psimrcc/psimrcc_scf, casscf/rasscf, dmrg-scf/dmrg-ci/dmrg-caspt2

for m in mrcc_methods:
    method_governing_type_keywords[m] = "cc_type"

for lvl in range(2, 99):
    method_governing_type_keywords[f"ci{lvl}"] = "ci_type"
    method_governing_type_keywords[f"zapt{lvl}"] = "mp_type"
    if lvl >= 5:
        method_governing_type_keywords[f"mp{lvl}"] = "mp_type"


def method_algorithm_type(name: str) -> Tuple[str, str, str]:
    """For method, return useful information about algorithm keyword.

    Parameters
    ----------
    name
        Method argument, case insensitive, e.g. `mp4`.

    Returns
    -------
    keyword : str
        Name of read_options keyword governing the method, e.g., `MP_TYPE` for MP4.
    default : str
        Default value of keyword, e.g., `CONV` for `MP_TYPE`.
    now : str
        Currently set value of keyword, e.g., `DF` after `set mp_type df`.

    """
    from psi4 import core

    # not single source-of-truth, since "truth" lives in un-probe-able read_options.cc defaults and in proc.py per-method coding
    type_keywords_defaults = {
        "SCF_TYPE": "DF",  # actually "PK" by read_options, but "DF" as practical default when SCF is target method
        "DCT_TYPE": "CONV",
        "MP2_TYPE": "DF",
        "MP_TYPE": "CONV",
        "CC_TYPE": "CONV",
        "CI_TYPE": "CONV",
    }

    lowername = name.lower()
    keyword = method_governing_type_keywords[lowername].upper()
    default = type_keywords_defaults[keyword]

    if lowername in ["mp2.5", "fno-mp2.5", "mp3", "fno-mp3",
            # below are newly treated like their parent above c. Sept 2022
            "custom-scs-mp2.5", "scs-mp3", "custom-scs-mp3"]:
        now = core.get_global_option(keyword) if core.has_global_option_changed(keyword) else "DF"
    elif keyword == "DCT_TYPE":
        now = core.get_option("DCT", keyword)
    else:
        now = core.get_global_option(keyword)

    MethodType = namedtuple("MethodType", "keyword default now")
    return MethodType(keyword, default, now)
