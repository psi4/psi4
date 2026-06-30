#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2026 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#
"""
List of XC functionals
"""

funcs = []

# yapf: disable
funcs.append({"name": "TETER93"        , "xc_functionals": {"LDA_XC_TETER93"            : {}}})
funcs.append({"name": "ZLP"            , "xc_functionals": {"LDA_XC_ZLP"                : {}}})
funcs.append({"name": "KSDT"           , "xc_functionals": {"LDA_XC_KSDT"               : {}}})
funcs.append({"name": "LDA0"           , "xc_functionals": {"HYB_LDA_XC_LDA0"           : {}}})
funcs.append({"name": "CAM-LDA0"       , "xc_functionals": {"HYB_LDA_XC_CAM_LDA0"       : {}}})
funcs.append({"name": "OPBE-D"         , "xc_functionals": {"GGA_XC_OPBE_D"             : {}}, "dispersion": {"type": "d2", "params": {'s6': 1.0,  'alpha6': 20.0, 'sr6': 1.15}, "citation": '    L. Goerigk, S. Grimme, J. Chem. Theory. Comput. 6, 107-126, 2010\n'}, 'alias': ['OPBE-D2']})
funcs.append({"name": "OPWLYP-D"       , "xc_functionals": {"GGA_XC_OPWLYP_D"           : {}}, "dispersion": {"type": "d2", "params": {'s6': 1.0,  'alpha6': 20.0, 'sr6': 1.15}, "citation": '    L. Goerigk, S. Grimme, J. Chem. Theory. Comput. 6, 107-126, 2010\n'}, 'alias': ['OPWLYP-D2']})
funcs.append({"name": "OBLYP-D"        , "xc_functionals": {"GGA_XC_OBLYP_D"            : {}}, "dispersion": {"type": "d2", "params": {'s6': 1.0,  'alpha6': 20.0, 'sr6': 1.15}, "citation": '    L. Goerigk, S. Grimme, J. Chem. Theory. Comput. 6, 107-126, 2010\n'}, 'alias': ['OBLYP-D2']})
funcs.append({"name": "HCTH407P"       , "xc_functionals": {"GGA_XC_HCTH_407P"          : {}}})
funcs.append({"name": "HCTHP76"        , "xc_functionals": {"GGA_XC_HCTH_P76"           : {}}})
funcs.append({"name": "HCTHP14"        , "xc_functionals": {"GGA_XC_HCTH_P14"           : {}}})
funcs.append({"name": "B97-GGA1"       , "xc_functionals": {"GGA_XC_B97_GGA1"           : {}}})
funcs.append({"name": "KT2"            , "xc_functionals": {"GGA_XC_KT2"                : {}}})
funcs.append({"name": "TH1"            , "xc_functionals": {"GGA_XC_TH1"                : {}}})
funcs.append({"name": "TH2"            , "xc_functionals": {"GGA_XC_TH2"                : {}}})
funcs.append({"name": "TH3"            , "xc_functionals": {"GGA_XC_TH3"                : {}}})
funcs.append({"name": "TH4"            , "xc_functionals": {"GGA_XC_TH4"                : {}}})
funcs.append({"name": "HCTH93"         , "xc_functionals": {"GGA_XC_HCTH_93"            : {}}, "alias": ["HCTH"]})
funcs.append({"name": "HCTH120"        , "xc_functionals": {"GGA_XC_HCTH_120"           : {}}})
funcs.append({"name": "HCTH147"        , "xc_functionals": {"GGA_XC_HCTH_147"           : {}}})
funcs.append({"name": "HCTH407"        , "xc_functionals": {"GGA_XC_HCTH_407"           : {}}})
funcs.append({"name": "EDF1"           , "xc_functionals": {"GGA_XC_EDF1"               : {}}})
funcs.append({"name": "XLYP"           , "xc_functionals": {"GGA_XC_XLYP"               : {}}})
funcs.append({"name": "B97-D"          , "xc_functionals": {"GGA_XC_B97_D"              : {}}, "dispersion": {"type": "d2", "params": {'s6': 1.25, 'alpha6': 20.0, 'sr6': 1.1}}, "alias": ["B97-D2"]})
funcs.append({"name": "PBE1W"          , "xc_functionals": {"GGA_XC_PBE1W"              : {}}})
funcs.append({"name": "mPWLYP1W"       , "xc_functionals": {"GGA_XC_MPWLYP1W"           : {}}})
funcs.append({"name": "PBELYP1W"       , "xc_functionals": {"GGA_XC_PBELYP1W"           : {}}})
funcs.append({"name": "MOHLYP"         , "xc_functionals": {"GGA_XC_MOHLYP"             : {}}})
funcs.append({"name": "MOHLYP2"        , "xc_functionals": {"GGA_XC_MOHLYP2"            : {}}})
funcs.append({"name": "TH-FL"          , "xc_functionals": {"GGA_XC_TH_FL"              : {}}})
funcs.append({"name": "TH-FC"          , "xc_functionals": {"GGA_XC_TH_FC"              : {}}})
funcs.append({"name": "TH-FCFO"        , "xc_functionals": {"GGA_XC_TH_FCFO"            : {}}})
funcs.append({"name": "TH-FCO"         , "xc_functionals": {"GGA_XC_TH_FCO"             : {}}})
funcs.append({"name": "VV10"           , "xc_functionals": {"GGA_XC_VV10"               : {}}})
funcs.append({"name": "B97-1p"         , "xc_functionals": {"HYB_GGA_XC_B97_1p"         : {}}})
funcs.append({"name": "B3PW91"         , "xc_functionals": {"HYB_GGA_XC_B3PW91"         : {}}})
funcs.append({"name": "B3LYP"          , "xc_functionals": {"HYB_GGA_XC_B3LYP"          : {}}})
funcs.append({"name": "B3P86"          , "xc_functionals": {"HYB_GGA_XC_B3P86"          : {}}})
funcs.append({"name": "O3LYP"          , "xc_functionals": {"HYB_GGA_XC_O3LYP"          : {}}})
funcs.append({"name": "mPW1K"          , "xc_functionals": {"HYB_GGA_XC_mPW1K"          : {}}})
funcs.append({"name": "PBE0"           , "xc_functionals": {"HYB_GGA_XC_PBEH"           : {}}, "alias": ["PBEH"]})
funcs.append({"name": "B97-0"          , "xc_functionals": {"HYB_GGA_XC_B97"            : {}}})
funcs.append({"name": "B97-1"          , "xc_functionals": {"HYB_GGA_XC_B97_1"          : {}}})
funcs.append({"name": "B97-2"          , "xc_functionals": {"HYB_GGA_XC_B97_2"          : {}}})
funcs.append({"name": "X3LYP"          , "xc_functionals": {"HYB_GGA_XC_X3LYP"          : {}}})
funcs.append({"name": "B1WC"           , "xc_functionals": {"HYB_GGA_XC_B1WC"           : {}}})
funcs.append({"name": "B97-K"          , "xc_functionals": {"HYB_GGA_XC_B97_K"          : {}}})
funcs.append({"name": "B97-3"          , "xc_functionals": {"HYB_GGA_XC_B97_3"          : {}}})
funcs.append({"name": "mPW3PW"         , "xc_functionals": {"HYB_GGA_XC_MPW3PW"         : {}}})
funcs.append({"name": "B1LYP"          , "xc_functionals": {"HYB_GGA_XC_B1LYP"          : {}}})
funcs.append({"name": "B1PW91"         , "xc_functionals": {"HYB_GGA_XC_B1PW91"         : {}}})
funcs.append({"name": "mPW1PW"         , "xc_functionals": {"HYB_GGA_XC_mPW1PW"         : {}}, "alias": ["MPW1PW91"]})
funcs.append({"name": "mPW3LYP"        , "xc_functionals": {"HYB_GGA_XC_MPW3LYP"        : {}}})
funcs.append({"name": "SB98-1a"        , "xc_functionals": {"HYB_GGA_XC_SB98_1a"        : {}}})
funcs.append({"name": "SB98-1b"        , "xc_functionals": {"HYB_GGA_XC_SB98_1b"        : {}}})
funcs.append({"name": "SB98-1c"        , "xc_functionals": {"HYB_GGA_XC_SB98_1c"        : {}}})
funcs.append({"name": "SB98-2a"        , "xc_functionals": {"HYB_GGA_XC_SB98_2a"        : {}}})
funcs.append({"name": "SB98-2b"        , "xc_functionals": {"HYB_GGA_XC_SB98_2b"        : {}}})
funcs.append({"name": "SB98-2c"        , "xc_functionals": {"HYB_GGA_XC_SB98_2c"        : {}}})
funcs.append({"name": "HSE03"          , "xc_functionals": {"HYB_GGA_XC_HSE03"          : {}}})
funcs.append({"name": "HSE06"          , "xc_functionals": {"HYB_GGA_XC_HSE06"          : {}}})
funcs.append({"name": "HJS-PBE"        , "xc_functionals": {"HYB_GGA_XC_HJS_PBE"        : {}}})
funcs.append({"name": "HJS-PBE-SOL"    , "xc_functionals": {"HYB_GGA_XC_HJS_PBE_SOL"    : {}}, "alias": ["HJS-PBESOL"]})
funcs.append({"name": "HJS-B88"        , "xc_functionals": {"HYB_GGA_XC_HJS_B88"        : {}}})
funcs.append({"name": "HJS-B97X"       , "xc_functionals": {"HYB_GGA_XC_HJS_B97X"       : {}}})
funcs.append({"name": "CAM-B3LYP"      , "xc_functionals": {"HYB_GGA_XC_CAM_B3LYP"      : {}}})
funcs.append({"name": "TUNED-CAM-B3LYP", "xc_functionals": {"HYB_GGA_XC_TUNED_CAM_B3LYP": {}}})
funcs.append({"name": "BHandH"         , "xc_functionals": {"HYB_GGA_XC_BHANDH"         : {}}})
funcs.append({"name": "BHandHLYP"      , "xc_functionals": {"HYB_GGA_XC_BHANDHLYP"      : {}}, "alias": ["BHHLYP"]})
funcs.append({"name": "MB3LYP-RC04"    , "xc_functionals": {"HYB_GGA_XC_MB3LYP_RC04"    : {}}})
funcs.append({"name": "mPWLYP1M"       , "xc_functionals": {"HYB_GGA_XC_MPWLYP1M"       : {}}})
funcs.append({"name": "revB3LYP"       , "xc_functionals": {"HYB_GGA_XC_REVB3LYP"       : {}}})
#funcs.append({"name": "CAMY-BLYP"      , "xc_functionals": {"HYB_GGA_XC_CAMY_BLYP"      : {}}}) # CAMY range separation not supported by psi4
funcs.append({"name": "PBE0-13"        , "xc_functionals": {"HYB_GGA_XC_PBE0_13"        : {}}})
funcs.append({"name": "B3LYPs"         , "xc_functionals": {"HYB_GGA_XC_B3LYPs"         : {}}})
funcs.append({"name": "wB97"           , "xc_functionals": {"HYB_GGA_XC_WB97"           : {}}})
funcs.append({"name": "wB97X"          , "xc_functionals": {"HYB_GGA_XC_WB97X"          : {}}})
funcs.append({"name": "wB97X-D"        , "xc_functionals": {"HYB_GGA_XC_WB97X_D"        : {}}, "dispersion": {"type": "chg", "params": {"s6": 1.0}}})
funcs.append({"name": "wB97X-D3"       , "xc_functionals": {"HYB_GGA_XC_WB97X_D3"       : {}}, "dispersion": {"type": "d3zero2b", "params": {'s6': 1.0,'s8':1.0,'sr6': 1.281,'sr8': 1.094,'alpha6': 14.0}}})
funcs.append({"name": "LRC-wPBEh"      , "xc_functionals": {"HYB_GGA_XC_LRC_WPBEH"      : {}}})
funcs.append({"name": "wB97X-V"        , "xc_functionals": {"HYB_GGA_XC_WB97X_V"        : {}}, "alias": ["WB97XV"]})
#funcs.append({"name": "LCY-PBE"        , "xc_functionals": {"HYB_GGA_XC_LCY_PBE"        : {}}}) # LCY range separation not supported by psi4
#funcs.append({"name": "LCY-BLYP"       , "xc_functionals": {"HYB_GGA_XC_LCY_BLYP"       : {}}}) # LCY range separation not supported by psi4
funcs.append({"name": "LC-VV10"        , "xc_functionals": {"HYB_GGA_XC_LC_VV10"        : {}}})
#funcs.append({"name": "CAMY-B3LYP"     , "xc_functionals": {"HYB_GGA_XC_CAMY_B3LYP"     : {}}}) # CAMY range separation not supported by psi4
funcs.append({"name": "HPBEINT"        , "xc_functionals": {"HYB_GGA_XC_HPBEINT"        : {}}})
funcs.append({"name": "LRC-WPBE"       , "xc_functionals": {"HYB_GGA_XC_LRC_WPBE"       : {}}})
funcs.append({"name": "B3LYP5"         , "xc_functionals": {"HYB_GGA_XC_B3LYP5"         : {}}})
funcs.append({"name": "EDF2"           , "xc_functionals": {"HYB_GGA_XC_EDF2"           : {}}})
funcs.append({"name": "CAP0"           , "xc_functionals": {"HYB_GGA_XC_CAP0"           : {}}})
funcs.append({"name": "B88B95"         , "xc_functionals": {"HYB_MGGA_XC_B88B95"        : {}}, "alias": ["B1B95"]})
funcs.append({"name": "B86B95"         , "xc_functionals": {"HYB_MGGA_XC_B86B95"        : {}}})
funcs.append({"name": "PW86B95"        , "xc_functionals": {"HYB_MGGA_XC_PW86B95"       : {}}})
funcs.append({"name": "BB1K"           , "xc_functionals": {"HYB_MGGA_XC_BB1K"          : {}}})
funcs.append({"name": "mPW1B95"        , "xc_functionals": {"HYB_MGGA_XC_MPW1B95"       : {}}})
funcs.append({"name": "mPWB1K"         , "xc_functionals": {"HYB_MGGA_XC_MPWB1K"        : {}}})
funcs.append({"name": "X1B95"          , "xc_functionals": {"HYB_MGGA_XC_X1B95"         : {}}})
funcs.append({"name": "XB1K"           , "xc_functionals": {"HYB_MGGA_XC_XB1K"          : {}}})
#funcs.append({"name": "PW6B95"         , "xc_functionals": {"HYB_MGGA_XC_PW6B95"        : {}}})  # duplicate of dict-based in hyb_functionals.py
funcs.append({"name": "PWB6K"          , "xc_functionals": {"HYB_MGGA_XC_PWB6K"         : {}}})
funcs.append({"name": "TPSSh"          , "xc_functionals": {"HYB_MGGA_XC_TPSSH"         : {}}, "alias": ["TPSS0"]})
funcs.append({"name": "revTPSSh"       , "xc_functionals": {"HYB_MGGA_XC_REVTPSSH"      : {}}})
funcs.append({"name": "wB97M-V"        , "xc_functionals": {"HYB_MGGA_XC_WB97M_V"       : {}}, "alias": ["WB97MV"]})
funcs.append({"name": "ZLP"            , "xc_functionals": {"MGGA_XC_ZLP"               : {}}})
funcs.append({"name": "OTPSS-D"        , "xc_functionals": {"MGGA_XC_OTPSS_D"           : {}}, "dispersion": {"type": "d2", "params": {'s6': 1.0,  'alpha6': 20.0, 'sr6': 1.15}, "citation": '    L. Goerigk, S. Grimme, J. Chem. Theory. Comput. 6, 107-126, 2010\n'}, 'alias': ['OTPSS-D2']})
funcs.append({"name": "TPSSLYP1W"      , "xc_functionals": {"MGGA_XC_TPSSLYP1W"         : {}}})
funcs.append({"name": "B97M-V"         , "xc_functionals": {"MGGA_XC_B97M_V"            : {}}})
funcs.append({"name": "B5050LYP"       , "xc_functionals": {"HYB_GGA_XC_B5050LYP"       : {}}})
funcs.append({"name": "LC-BOP"         , "xc_functionals": {"HYB_GGA_XC_LC_BOP"         : {}}, "alias": ["LRC-BOP"]})
funcs.append({"name": "KMLYP"          , "xc_functionals": {"HYB_GGA_XC_KMLYP"          : {}}})
funcs.append({"name": "R2SCAN0"        , "xc_functionals": {"HYB_MGGA_XC_R2SCAN0"       : {}}})
funcs.append({"name": "R2SCANh"        , "xc_functionals": {"HYB_MGGA_XC_R2SCANH"       : {}}})
funcs.append({"name": "R2SCAN50"       , "xc_functionals": {"HYB_MGGA_XC_R2SCAN50"      : {}}})
funcs.append({"name": "KT1"             , "xc_functionals": {"GGA_XC_KT1"                 : {}}})
funcs.append({"name": "NCAP"            , "xc_functionals": {"GGA_XC_NCAP"                : {}}})
funcs.append({"name": "PBE-MOL0"        , "xc_functionals": {"HYB_GGA_XC_PBE_MOL0"        : {}}})
funcs.append({"name": "PBE-SOL0"        , "xc_functionals": {"HYB_GGA_XC_PBE_SOL0"        : {}}})
funcs.append({"name": "PBEB0"           , "xc_functionals": {"HYB_GGA_XC_PBEB0"           : {}}})
funcs.append({"name": "PBE-MOLB0"       , "xc_functionals": {"HYB_GGA_XC_PBE_MOLB0"       : {}}})
funcs.append({"name": "BEEFVDW"         , "xc_functionals": {"GGA_XC_BEEFVDW"             : {}}})
funcs.append({"name": "HLE17"           , "xc_functionals": {"MGGA_XC_HLE17"              : {}}})
funcs.append({"name": "HFLYP"           , "xc_functionals": {"HYB_GGA_XC_HFLYP"           : {}}})
funcs.append({"name": "B3P86-NWCHEM"    , "xc_functionals": {"HYB_GGA_XC_B3P86_NWCHEM"    : {}}})
funcs.append({"name": "RELPBE0"         , "xc_functionals": {"HYB_GGA_XC_RELPBE0"         : {}}})
funcs.append({"name": "OPB3LYP"         , "xc_functionals": {"HYB_GGA_XC_OPB3LYP"         : {}}})
funcs.append({"name": "CASE21"          , "xc_functionals": {"HYB_GGA_XC_CASE21"          : {}}})
funcs.append({"name": "PBE-2X"          , "xc_functionals": {"HYB_GGA_XC_PBE_2X"          : {}}})
funcs.append({"name": "PBE38"           , "xc_functionals": {"HYB_GGA_XC_PBE38"           : {}}})
funcs.append({"name": "B3LYP3"          , "xc_functionals": {"HYB_GGA_XC_B3LYP3"          : {}}})
funcs.append({"name": "CAM-O3LYP"       , "xc_functionals": {"HYB_GGA_XC_CAM_O3LYP"       : {}}})
funcs.append({"name": "LC-BLYP"         , "xc_functionals": {"HYB_GGA_XC_LC_BLYP"         : {}}})
funcs.append({"name": "APF"             , "xc_functionals": {"HYB_GGA_XC_APF"             : {}}})
funcs.append({"name": "QTP17"           , "xc_functionals": {"HYB_GGA_XC_QTP17"           : {}}})
funcs.append({"name": "B3LYP-MCM1"      , "xc_functionals": {"HYB_GGA_XC_B3LYP_MCM1"      : {}}})
funcs.append({"name": "B3LYP-MCM2"      , "xc_functionals": {"HYB_GGA_XC_B3LYP_MCM2"      : {}}})
funcs.append({"name": "HSE12"           , "xc_functionals": {"HYB_GGA_XC_HSE12"           : {}}})
funcs.append({"name": "HSE12S"          , "xc_functionals": {"HYB_GGA_XC_HSE12S"          : {}}})
funcs.append({"name": "HSE-SOL"         , "xc_functionals": {"HYB_GGA_XC_HSE_SOL"         : {}}})
funcs.append({"name": "CAM-QTP-01"      , "xc_functionals": {"HYB_GGA_XC_CAM_QTP_01"      : {}}})
funcs.append({"name": "LC-WPBE-WHS"     , "xc_functionals": {"HYB_GGA_XC_LC_WPBE_WHS"     : {}}})
funcs.append({"name": "LC-WPBEH-WHS"    , "xc_functionals": {"HYB_GGA_XC_LC_WPBEH_WHS"    : {}}})
funcs.append({"name": "LC-WPBE08-WHS"   , "xc_functionals": {"HYB_GGA_XC_LC_WPBE08_WHS"   : {}}})
funcs.append({"name": "LC-WPBESOL-WHS"  , "xc_functionals": {"HYB_GGA_XC_LC_WPBESOL_WHS"  : {}}})
funcs.append({"name": "CAM-QTP-00"      , "xc_functionals": {"HYB_GGA_XC_CAM_QTP_00"      : {}}})
funcs.append({"name": "CAM-QTP-02"      , "xc_functionals": {"HYB_GGA_XC_CAM_QTP_02"      : {}}})
funcs.append({"name": "LC-QTP"          , "xc_functionals": {"HYB_GGA_XC_LC_QTP"          : {}}})
funcs.append({"name": "BLYP35"          , "xc_functionals": {"HYB_GGA_XC_BLYP35"          : {}}})
funcs.append({"name": "1D-EHWLRG-1"     , "xc_functionals": {"LDA_XC_1D_EHWLRG_1"         : {}}})
funcs.append({"name": "1D-EHWLRG-2"     , "xc_functionals": {"LDA_XC_1D_EHWLRG_2"         : {}}})
funcs.append({"name": "1D-EHWLRG-3"     , "xc_functionals": {"LDA_XC_1D_EHWLRG_3"         : {}}})
funcs.append({"name": "HLE16"           , "xc_functionals": {"GGA_XC_HLE16"               : {}}})
funcs.append({"name": "LP-A"            , "xc_functionals": {"LDA_XC_LP_A"                : {}}})
funcs.append({"name": "LP-B"            , "xc_functionals": {"LDA_XC_LP_B"                : {}}})
funcs.append({"name": "B0KCIS"          , "xc_functionals": {"HYB_MGGA_XC_B0KCIS"         : {}}})
funcs.append({"name": "MPW1KCIS"        , "xc_functionals": {"HYB_MGGA_XC_MPW1KCIS"       : {}}})
funcs.append({"name": "MPWKCIS1K"       , "xc_functionals": {"HYB_MGGA_XC_MPWKCIS1K"      : {}}})
funcs.append({"name": "PBE1KCIS"        , "xc_functionals": {"HYB_MGGA_XC_PBE1KCIS"       : {}}})
funcs.append({"name": "TPSS1KCIS"       , "xc_functionals": {"HYB_MGGA_XC_TPSS1KCIS"      : {}}})
funcs.append({"name": "KT3"             , "xc_functionals": {"GGA_XC_KT3"                 : {}}})
funcs.append({"name": "LB07"            , "xc_functionals": {"HYB_GGA_XC_LB07"            : {}}})
funcs.append({"name": "TIH"             , "xc_functionals": {"LDA_XC_TIH"                 : {}}})
funcs.append({"name": "APBE0"           , "xc_functionals": {"HYB_GGA_XC_APBE0"           : {}}})
funcs.append({"name": "HAPBE"           , "xc_functionals": {"HYB_GGA_XC_HAPBE"           : {}}})
funcs.append({"name": "RCAM-B3LYP"      , "xc_functionals": {"HYB_GGA_XC_RCAM_B3LYP"      : {}}})
funcs.append({"name": "WC04"            , "xc_functionals": {"HYB_GGA_XC_WC04"            : {}}})
funcs.append({"name": "WP04"            , "xc_functionals": {"HYB_GGA_XC_WP04"            : {}}})
funcs.append({"name": "CAMH-B3LYP"      , "xc_functionals": {"HYB_GGA_XC_CAMH_B3LYP"      : {}}})
funcs.append({"name": "WHPBE0"          , "xc_functionals": {"HYB_GGA_XC_WHPBE0"          : {}}})
funcs.append({"name": "LC-BLYP-EA"      , "xc_functionals": {"HYB_GGA_XC_LC_BLYP_EA"      : {}}})
funcs.append({"name": "LC-PBEOP"        , "xc_functionals": {"HYB_GGA_XC_LC_PBEOP"        : {}}})
funcs.append({"name": "LC-BLYPR"        , "xc_functionals": {"HYB_GGA_XC_LC_BLYPR"        : {}}})
funcs.append({"name": "MCAM-B3LYP"      , "xc_functionals": {"HYB_GGA_XC_MCAM_B3LYP"      : {}}})
funcs.append({"name": "VCML-RVV10"      , "xc_functionals": {"MGGA_XC_VCML_RVV10"         : {}}})
funcs.append({"name": "GAS22"           , "xc_functionals": {"HYB_MGGA_XC_GAS22"          : {}}})
funcs.append({"name": "CAM-PBEH"        , "xc_functionals": {"HYB_GGA_XC_CAM_PBEH"        : {}}})
funcs.append({"name": "LC-TMLYP"        , "xc_functionals": {"HYB_MGGA_XC_LC_TMLYP"       : {}}})
# yapf: enable

functional_list = {}
for functional in funcs:
    functional_list[functional["name"].lower()] = functional
