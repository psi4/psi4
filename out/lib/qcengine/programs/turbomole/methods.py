# Available methods

METHODS = {
    # Dummy for HF
    "hf": ["hf"],
    "ricc2": ["rimp2", "rimp3", "rimp4", "ricc2"],
    # Hardcoded XC-functionals that can be selected from the dft submenu
    # of define.
    "dft_hardcoded": [
        # Hardcoded in V7.3
        "s-vwn",
        "s-vwn_Gaussian",
        "pwlda",
        "b-lyp",
        "b-vwn",
        "b-p",
        "pbe",
        "tpss",
        "bh-lyp",
        "b3-lyp",
        "b3-lyp_Gaussian",
        "pbe0",
        "tpssh",
        "pw6b95",
        "m06",
        "m06-l",
        "m06-2x",
        "lhf",
        "oep",
        "b97-d",
        "pbeh-3c",
        "b97-3c",
        "lh07t-svwn",
        "lh07s-svwn",
        "lh12ct-ssirpw92",
        "lh12ct-ssifpw92",
        "lh14t-calpbe",
        # Hardcoded in V7.4
        "cam-b3lyp",
        # B2PLYP his is not easily supported right now as we would need an
        # additional MP2 calculation from rimp2/ricc2.
        # "b2-plyp",
    ],
    # Shorctus for XC functionals in V7.4 using LibXC
    "dft_libxc": [
        "wb97",
        "wb97x",
        "sogga11",
        "sogga-11x",
        "mn12-l",
        "mn12-sx",
        "mn15",
        "mn15-l",
        "m06-libxc",
        "cam-b3lyp-libxc",
        "hse06-libxc",
    ],
}

# Available keywords
KEYWORDS = {
    # Resolution of identity
    "ri": ["rijk", "ri", "marij"],
    # Dispersion correction
    "dsp": ["d3", "d3bj"],
}
