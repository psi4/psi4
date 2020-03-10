import numpy as np

#   "bh_h2p" "cc-pvdz"
#   "hf" "cc-pvdz"
# "h2o" "aug-cc-pvdz"
# "h2o" "qz2p"
# "nh2" "aug-cc-pvdz"
# "nh2" "qz2p"

_scf_hf_dz_df_rhf = -100.019400605629
_scf_bh3p_dz_df_uhf = -25.945130559147
_scf_bh3p_dz_df_rohf = -25.943606522029
_scf_hf_dz_pk_rhf = -100.01941126902270
_scf_bh3p_dz_pk_uhf = -25.94513842869638
_scf_bh3p_dz_pk_rohf = -25.943614318546

_std_suite = [
    # <<<  scf DF, fc  >>>
    {
        "meta": {
            "system": "hf",
            "basis": "cc-pvdz",
            "scf_type": "df",
            "reference": "rhf",
            "fc": "fc",
            "mp2_type": "conv",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_hf_dz_df_rhf,
            "MP2 CORRELATION ENERGY": -0.201612517228,
            #            "MP2 TOTAL ENERGY": -100.221013122857,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.05348507322421174,
        },
    },
    {
        "meta": {
            "system": "hf",
            "basis": "cc-pvdz",
            "scf_type": "df",
            "reference": "rhf",
            "fc": "fc",
            "mp2_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_hf_dz_df_rhf,
            "MP2 CORRELATION ENERGY": -0.201610660387,
            #            "MP2 TOTAL ENERGY": -100.221011266016,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.0535212487451535,
            "MP2 TOTAL GRADIENT": np.array(
                [  # dfmp2 findif-5 fc df+df
                    0.00000000000000,
                    0.00000000000000,
                    0.00314716362539,
                    0.00000000000000,
                    0.00000000000000,
                    -0.00314716362539,
                ]
            ).reshape((-1, 3)),
        },
    },
    {
        "meta": {
            "system": "hf",
            "basis": "cc-pvdz",
            "scf_type": "df",
            "reference": "rhf",
            "fc": "fc",
            "mp2_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_hf_dz_df_rhf,
            "MP2 CORRELATION ENERGY": -0.201609396752,
            #            "MP2 TOTAL ENERGY": -100.221010002381,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.4,
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "df",
            "reference": "uhf",
            "fc": "fc",
            "mp2_type": "conv",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_df_uhf,
            "MP2 CORRELATION ENERGY": -0.058421122206,
            #            "MP2 TOTAL ENERGY": -26.003551681354,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.5,
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "df",
            "reference": "uhf",
            "fc": "fc",
            "mp2_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_df_uhf,
            "MP2 CORRELATION ENERGY": -0.058390006825,
            #            "MP2 TOTAL ENERGY": -26.003520565972,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.001768919072594215,
            "MP2 TOTAL GRADIENT": np.array(
                [  # dfmp2 findif-5 fc df+df
                    0.00000000000000,
                    0.00000000000000,
                    0.01231996225662,
                    0.00000000000000,
                    0.00000000000000,
                    -0.01186374280678,
                    0.00000000000000,
                    0.01031743020277,
                    -0.00022810972492,
                    0.00000000000000,
                    -0.01031743020277,
                    -0.00022810972492,
                ]
            ).reshape((-1, 3)),
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "df",
            "reference": "uhf",
            "fc": "fc",
            "mp2_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_df_uhf,
            "MP2 CORRELATION ENERGY": -0.058409837177,
            #            "MP2 TOTAL ENERGY": -26.003540396324,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.7,
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "df",
            "reference": "rohf",
            "fc": "fc",
            "mp2_type": "conv",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_df_rohf,
            "MP2 CORRELATION ENERGY": -0.060939211739,
            #            "MP2 TOTAL ENERGY": -26.004545733768,
            "MP2 SINGLES ENERGY": 1.1,
            "MP2 SAME-SPIN CORRELATION ENERGY": -1.1,
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "df",
            "reference": "rohf",
            "fc": "fc",
            "mp2_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_df_rohf,
            "MP2 CORRELATION ENERGY": -0.059372748391,
            #            "MP2 TOTAL ENERGY": -26.002979270420,
            "MP2 SINGLES ENERGY": -0.000688391888527046,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.0018535174789756292,
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "df",
            "reference": "rohf",
            "fc": "fc",
            "mp2_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_df_rohf,
            "MP2 CORRELATION ENERGY": -0.059393510962,
            #            "MP2 TOTAL ENERGY": -26.003000032991,
            "MP2 SINGLES ENERGY": 1.3,
            "MP2 SAME-SPIN CORRELATION ENERGY": -1.3,
        },
    },
    # <<<  scf DF, ae  >>>
    {
        "meta": {
            "system": "hf",
            "basis": "cc-pvdz",
            "scf_type": "df",
            "reference": "rhf",
            "fc": "ae",
            "mp2_type": "conv",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_hf_dz_df_rhf,
            "MP2 CORRELATION ENERGY": -0.2037668844651997,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.05427252944164894,
        },
    },
    {
        "meta": {
            "system": "hf",
            "basis": "cc-pvdz",
            "scf_type": "df",
            "reference": "rhf",
            "fc": "ae",
            "mp2_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_hf_dz_df_rhf,
            "MP2 CORRELATION ENERGY": -0.2037649370559149,
            #            "MP2 TOTAL ENERGY": -100.2231655426856776,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.05430875283333263,
            "MP2 TOTAL GRADIENT": np.array(
                [  # dfmp2 findif-5 ae df+df
                    0.00000000000000,
                    0.00000000000000,
                    0.00279211492833,
                    0.00000000000000,
                    0.00000000000000,
                    -0.00279211492833,
                ]
            ).reshape((-1, 3)),
        },
    },
    {
        "meta": {
            "system": "hf",
            "basis": "cc-pvdz",
            "scf_type": "df",
            "reference": "rhf",
            "fc": "ae",
            "mp2_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_hf_dz_df_rhf,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -2.3,
            "MP2 CORRELATION ENERGY": -2.3,
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "df",
            "reference": "uhf",
            "fc": "ae",
            "mp2_type": "conv",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_df_uhf,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -2.4,
            "MP2 CORRELATION ENERGY": -2.4,
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "df",
            "reference": "uhf",
            "fc": "ae",
            "mp2_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_df_uhf,
            "MP2 CORRELATION ENERGY": -0.0594557966607590,
            #            "MP2 TOTAL ENERGY": -26.0045863558097601,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.001920220330437888,
            "MP2 TOTAL GRADIENT": np.array(
                [  # dfmp2 findif-5 ae df+df
                    0.00000000000000,
                    0.00000000000000,
                    0.01252024755551,
                    0.00000000000000,
                    0.00000000000000,
                    -0.01207773525598,
                    0.00000000000000,
                    0.01032204616770,
                    -0.00022125614977,
                    0.00000000000000,
                    -0.01032204616770,
                    -0.00022125614977,
                ]
            ).reshape((-1, 3)),
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "df",
            "reference": "uhf",
            "fc": "ae",
            "mp2_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_df_uhf,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -2.5,
            "MP2 CORRELATION ENERGY": -2.5,
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "df",
            "reference": "rohf",
            "fc": "ae",
            "mp2_type": "conv",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_df_rohf,
            "MP2 SINGLES ENERGY": -2.7,
            "MP2 SAME-SPIN CORRELATION ENERGY": -2.7,
            "MP2 CORRELATION ENERGY": -2.7,
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "df",
            "reference": "rohf",
            "fc": "ae",
            "mp2_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_df_rohf,
            "MP2 CORRELATION ENERGY": -0.0604436327328384,
            #            "MP2 TOTAL ENERGY": -26.0040501547626377,
            "MP2 SINGLES ENERGY": -0.0006940750313001934,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.0020065676314032863,
            "MP2 TOTAL GRADIENT": np.array(
                [  # occ findif-5 ae df+df
                    0.00000000000000,
                    0.00000000000000,
                    0.01361287313486,
                    0.00000000000000,
                    0.00000000000000,
                    -0.01314329502424,
                    0.00000000000000,
                    0.01029838165151,
                    -0.00023478905531,
                    0.00000000000000,
                    -0.01029838165151,
                    -0.00023478905531,
                ]
            ).reshape((-1, 3)),
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "df",
            "reference": "rohf",
            "fc": "ae",
            "mp2_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_df_rohf,
            "MP2 SINGLES ENERGY": -2.8,
            "MP2 SAME-SPIN CORRELATION ENERGY": -2.8,
            "MP2 CORRELATION ENERGY": -2.8,
        },
    },
    # <<<  scf CONV, fc  >>>
    {
        "meta": {
            "system": "hf",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "rhf",
            "fc": "fc",
            "mp2_type": "conv",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_hf_dz_pk_rhf,
            "MP2 CORRELATION ENERGY": -0.201627516796,
            #            "MP2 TOTAL ENERGY": -100.221038785818,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.0534895025483176,
            "MP2 TOTAL GRADIENT": np.array(
                [  # fnocc findif-5 fc pk+conv
                    0.00000000000000,
                    0.00000000000000,
                    0.00317450456474,
                    0.00000000000000,
                    0.00000000000000,
                    -0.00317450456474,
                ]
            ).reshape((-1, 3)),
        },
    },
    {
        "meta": {
            "system": "hf",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "rhf",
            "fc": "fc",
            "mp2_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_hf_dz_pk_rhf,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.05352569481658172,
            "MP2 CORRELATION ENERGY": -0.20162566806258586,
        },
    },
    {
        "meta": {
            "system": "hf",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "rhf",
            "fc": "fc",
            "mp2_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_hf_dz_pk_rhf,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -2.9,
            "MP2 CORRELATION ENERGY": -2.9,
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "uhf",
            "fc": "fc",
            "mp2_type": "conv",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_pk_uhf,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.001767468898,
            "MP2 CORRELATION ENERGY": -0.058423513790,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 TOTAL GRADIENT": np.array(
                [
                    0.000000000000000,
                    0.000000000000000,
                    -0.012305278627642,
                    0.000000000000000,
                    0.000000000000000,
                    0.011851332672482,
                    0.000000000000000,
                    -0.010327045553422,
                    0.000226972977580,
                    0.000000000000000,
                    0.010327045553422,
                    0.000226972977580,
                ]
            ).reshape((-1, 3)),
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "uhf",
            "fc": "fc",
            "mp2_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_pk_uhf,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.0017690135626491292,
            "MP2 CORRELATION ENERGY": -0.058392397606538686,
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "uhf",
            "fc": "fc",
            "mp2_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_pk_uhf,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -2.33,
            "MP2 CORRELATION ENERGY": -2.33,
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "rohf",
            "fc": "fc",
            "mp2_type": "conv",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_pk_rohf,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.001851937488,
            "MP2 SINGLES ENERGY": -0.000688368657,
            "MP2 CORRELATION ENERGY": -0.058718885599,
            "MP2 TOTAL GRADIENT": np.array(
                [
                    0.000000000000000,
                    0.000000000000000,
                    -0.013388410166131,
                    0.000000000000000,
                    0.000000000000000,
                    0.012907368096590,
                    0.000000000000000,
                    -0.010303507439169,
                    0.000240521034770,
                    0.000000000000000,
                    0.010303507439169,
                    0.000240521034770,
                ]
            ).reshape((-1, 3)),
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "rohf",
            "fc": "fc",
            "mp2_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_pk_rohf,
            "MP2 SINGLES ENERGY": -0.0006883686516107368,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.0018536363586657242,
            "MP2 CORRELATION ENERGY": -0.05937514348825628,
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "rohf",
            "fc": "fc",
            "mp2_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_pk_rohf,
            "MP2 SINGLES ENERGY": -2.2,
            "MP2 SAME-SPIN CORRELATION ENERGY": -2.2,
            "MP2 CORRELATION ENERGY": -2.2,
        },
    },
    # <<<  scf CONV, ae  >>>
    {
        "meta": {
            "system": "hf",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "rhf",
            "fc": "ae",
            "mp2_type": "conv",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_hf_dz_pk_rhf,
            "MP2 CORRELATION ENERGY": -0.203781911950,
            #            "MP2 TOTAL ENERGY": -100.223193180973,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.05427697023782003,
            "MP2 TOTAL GRADIENT": np.array(
                [  # fnocc findif-5 ae pk+conv
                    0.0000000000,
                    0.0000000000,
                    0.0028193375,
                    0.0000000000,
                    0.0000000000,
                    -0.0028193375,
                ]
            ).reshape((-1, 3)),
        },
    },
    {
        "meta": {
            "system": "hf",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "rhf",
            "fc": "ae",
            "mp2_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_hf_dz_pk_rhf,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 CORRELATION ENERGY": -0.20377997248921056,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.05431321036920538,
        },
    },
    {
        "meta": {
            "system": "hf",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "rhf",
            "fc": "ae",
            "mp2_type": "cd",
        },
        "data": {},
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "uhf",
            "fc": "ae",
            "mp2_type": "conv",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_pk_uhf,
            "MP2 CORRELATION ENERGY": -0.05948928003552,
            #            "MP2 TOTAL ENERGY": -26.00462770873190,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.001918693775,
            "MP2 TOTAL GRADIENT": np.array(
                [  # occ findif-5 ae pk+conv
                    0.00000000000000,
                    0.00000000000000,
                    0.01250561195911,
                    0.00000000000000,
                    0.00000000000000,
                    -0.01206536529299,
                    0.00000000000000,
                    0.01033165380573,
                    -0.00022012333306,
                    0.00000000000000,
                    -0.01033165380573,
                    -0.00022012333306,
                ]
            ).reshape((-1, 3)),
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "uhf",
            "fc": "ae",
            "mp2_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_pk_uhf,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 CORRELATION ENERGY": -0.05945820694747983,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.0019203155958724552,
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "uhf",
            "fc": "ae",
            "mp2_type": "cd",
        },
        "data": {},
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "rohf",
            "fc": "ae",
            "mp2_type": "conv",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_pk_rohf,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.002004909679,
            "MP2 SINGLES ENERGY": -0.000694049865,
            "MP2 CORRELATION ENERGY": -0.059784065292,
            "MP2 TOTAL GRADIENT": np.array(
                [
                    0.000000000000000,
                    0.000000000000000,
                    -0.013594741747853,
                    0.000000000000000,
                    0.000000000000000,
                    0.013127629532095,
                    0.000000000000000,
                    -0.010308255599051,
                    0.000233556107879,
                    0.000000000000000,
                    0.010308255599051,
                    0.000233556107879,
                ]
            ).reshape((-1, 3)),
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "rohf",
            "fc": "ae",
            "mp2_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_pk_rohf,
            "MP2 SINGLES ENERGY": -0.0006940498589629459,
            "MP2 CORRELATION ENERGY": -0.0604460449537298,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.0020066877639503184,
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "rohf",
            "fc": "ae",
            "mp2_type": "cd",
        },
        "data": {},
    },
]


for calc in _std_suite:
    if calc["data"]:
        calc["data"]["MP2 TOTAL ENERGY"] = calc["data"]["MP2 CORRELATION ENERGY"] + calc["data"]["HF TOTAL ENERGY"]
        calc["data"]["MP2 DOUBLES ENERGY"] = calc["data"]["MP2 CORRELATION ENERGY"] - calc["data"]["MP2 SINGLES ENERGY"]
        calc["data"]["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = (
            calc["data"]["MP2 CORRELATION ENERGY"]
            - calc["data"]["MP2 SAME-SPIN CORRELATION ENERGY"]
            - calc["data"]["MP2 SINGLES ENERGY"]
        )


def answer_hash(**kwargs):
    system = kwargs.pop("system")
    basis = kwargs.pop("basis")
    scf_type = kwargs.pop("scf_type")
    reference = kwargs.pop("reference")
    fc = kwargs.pop("fc")
    mp2_type = kwargs.pop("mp2_type")

    return "_".join([system, basis, scf_type, reference, fc, mp2_type])


std_suite = {answer_hash(**calc["meta"]): calc["data"] for calc in _std_suite}
