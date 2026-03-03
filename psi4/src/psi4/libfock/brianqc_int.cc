/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifdef USING_brianxc

#include "brianqc_int.h"

#include "psi4/libfunctional/LibXCfunctional.h"
#include "psi4/libfunctional/superfunctional.h"

#include "psi4/libmints/basisset.h"
#include "psi4/libmints/molecule.h"

#include "psi4/liboptions/liboptions.h"

#include <use_brian_wrapper.h>
#include <brian_macros.h>
#include <brian_common.h>
#include <brian_scf.h>

extern void checkBrian();
extern BrianCookie brianCookie;
extern bool brianEnable;
extern bool brianEnableDFT;

struct BrianGrid {
    std::vector<brianInt> atomBlockCounts;
    std::vector<brianInt> atomBlockOffsets;

    std::vector<brianInt> blockRadialCounts;
    std::vector<brianInt> blockRadialOffsets;
    std::vector<double> radialCoordinates;
    std::vector<double> radialWeights;

    std::vector<brianInt> blockAngularCounts;
    std::vector<brianInt> blockAngularOffsets;
    std::vector<double> angularCoordinates;
    std::vector<double> angularWeights;

    std::vector<double> atomRotationMatrices;
};

extern bool brianBuildingNLCGrid;

namespace psi {

void BrianQCBase::initialize() {
        // TODO: Brian's grid initialization needs to be separated from Psi4's and moved in here.
        static const std::map<std::string, brianInt> functionalIDMap = {
            {"XC_GGA_X_AIRY", BRIAN_FUNCTIONAL_GGA_AIRY_X},
            {"XC_GGA_X_AK13", BRIAN_FUNCTIONAL_GGA_AK13_X},
            {"XC_GGA_C_AM05", BRIAN_FUNCTIONAL_GGA_AM05_C},
            {"XC_GGA_X_AM05", BRIAN_FUNCTIONAL_GGA_AM05_X},
            {"XC_GGA_C_APBE", BRIAN_FUNCTIONAL_GGA_APBE_C},
            {"XC_GGA_X_APBE", BRIAN_FUNCTIONAL_GGA_APBE_X},
            {"XC_GGA_X_B86_MGC", BRIAN_FUNCTIONAL_GGA_B86_MGC_X},
            {"XC_GGA_X_B86_R", BRIAN_FUNCTIONAL_GGA_B86_R_X},
            {"XC_GGA_X_B86", BRIAN_FUNCTIONAL_GGA_B86_X},
            {"XC_GGA_X_B88M", BRIAN_FUNCTIONAL_GGA_B88M_X},
            {"XC_GGA_X_B88", BRIAN_FUNCTIONAL_GGA_B88_X},
            {"XC_GGA_XC_B97_D", BRIAN_FUNCTIONAL_GGA_B97_D3_XC},
            {"XC_GGA_XC_B97_D", BRIAN_FUNCTIONAL_GGA_B97_D_XC},
            {"XC_GGA_XC_B97_GGA1", BRIAN_FUNCTIONAL_GGA_B97_GGA1_XC},
            {"XC_GGA_X_BAYESIAN", BRIAN_FUNCTIONAL_GGA_BAYESIAN_X},
            {"XC_GGA_C_BCGP", BRIAN_FUNCTIONAL_GGA_BCGP_C},
            {"XC_GGA_X_BCGP", BRIAN_FUNCTIONAL_GGA_BCGP_X},
            {"XC_GGA_X_BEEFVDW", BRIAN_FUNCTIONAL_GGA_BEEFVDW_X},
            {"XC_GGA_XC_BEEFVDW", BRIAN_FUNCTIONAL_GGA_BEEFVDW_XC},
            {"XC_GGA_C_BMK", BRIAN_FUNCTIONAL_GGA_BMK_C},
            {"XC_GGA_X_BPCCAC", BRIAN_FUNCTIONAL_GGA_BPCCAC_X},
            {"XC_GGA_X_C09X", BRIAN_FUNCTIONAL_GGA_C09X_X},
            {"XC_GGA_X_CAP", BRIAN_FUNCTIONAL_GGA_CAP_X},
            {"XC_GGA_X_CHACHIYO", BRIAN_FUNCTIONAL_GGA_CHACHIYO_X},
            {"XC_GGA_C_CS1", BRIAN_FUNCTIONAL_GGA_CS1_C},
            {"XC_GGA_X_DK87_R1", BRIAN_FUNCTIONAL_GGA_DK87_R1_X},
            {"XC_GGA_X_DK87_R2", BRIAN_FUNCTIONAL_GGA_DK87_R2_X},
            {"XC_GGA_X_EB88", BRIAN_FUNCTIONAL_GGA_EB88_X},
            {"XC_GGA_XC_EDF1", BRIAN_FUNCTIONAL_GGA_EDF1_XC},
            {"XC_GGA_X_EV93", BRIAN_FUNCTIONAL_GGA_EV93_X},
            {"XC_GGA_X_FT97_A", BRIAN_FUNCTIONAL_GGA_FT97_A_X},
            {"XC_GGA_X_FT97_B", BRIAN_FUNCTIONAL_GGA_FT97_B_X},
            {"XC_GGA_C_FT97", BRIAN_FUNCTIONAL_GGA_FT97_C},
            {"XC_GGA_X_G96", BRIAN_FUNCTIONAL_GGA_G96_X},
            {"XC_GGA_X_GAM", BRIAN_FUNCTIONAL_GGA_GAM_X},
            {"XC_GGA_C_GAPC", BRIAN_FUNCTIONAL_GGA_GAPC_C},
            {"XC_GGA_C_GAPLOC", BRIAN_FUNCTIONAL_GGA_GAPLOC_C},
            {"XC_GGA_X_GG99", BRIAN_FUNCTIONAL_GGA_GG99_X},
            {"XC_GGA_XC_HCTH_120", BRIAN_FUNCTIONAL_GGA_HCTH_120_XC},
            {"XC_GGA_XC_HCTH_147", BRIAN_FUNCTIONAL_GGA_HCTH_147_XC},
            {"XC_GGA_XC_HCTH_407P", BRIAN_FUNCTIONAL_GGA_HCTH_407P_XC},
            {"XC_GGA_XC_HCTH_407", BRIAN_FUNCTIONAL_GGA_HCTH_407_XC},
            {"XC_GGA_XC_HCTH_93", BRIAN_FUNCTIONAL_GGA_HCTH_93_XC},
            {"XC_GGA_C_HCTH_A", BRIAN_FUNCTIONAL_GGA_HCTH_A_C},
            {"XC_GGA_X_HCTH_A", BRIAN_FUNCTIONAL_GGA_HCTH_A_X},
            {"XC_GGA_XC_HCTH_P14", BRIAN_FUNCTIONAL_GGA_HCTH_P14_XC},
            {"XC_GGA_XC_HCTH_P76", BRIAN_FUNCTIONAL_GGA_HCTH_P76_XC},
            {"XC_GGA_X_HERMAN", BRIAN_FUNCTIONAL_GGA_HERMAN_X},
            {"XC_GGA_X_HJS_B88_V2", BRIAN_FUNCTIONAL_GGA_HJS_B88_V2_X},
            {"XC_GGA_X_HJS_B88", BRIAN_FUNCTIONAL_GGA_HJS_B88_X},
            {"XC_GGA_X_HJS_B97X", BRIAN_FUNCTIONAL_GGA_HJS_B97X_X},
            {"XC_GGA_X_HJS_PBE_SOL", BRIAN_FUNCTIONAL_GGA_HJS_PBE_SOL_X},
            {"XC_GGA_X_HJS_PBE", BRIAN_FUNCTIONAL_GGA_HJS_PBE_X},
            {"XC_GGA_XC_HLE16", BRIAN_FUNCTIONAL_GGA_HLE16_XC},
            {"XC_GGA_X_HTBS", BRIAN_FUNCTIONAL_GGA_HTBS_X},
            {"XC_GGA_C_HYB_TAU_HCTH", BRIAN_FUNCTIONAL_GGA_HYB_TAU_HCTH_C},
            {"XC_GGA_X_ITYH", BRIAN_FUNCTIONAL_GGA_ITYH_X},
            {"XC_GGA_X_KGG99", BRIAN_FUNCTIONAL_GGA_KGG99_X},
            {"XC_GGA_X_KT1", BRIAN_FUNCTIONAL_GGA_KT1_X},
            {"XC_GGA_XC_KT1", BRIAN_FUNCTIONAL_GGA_KT1_XC},
            {"XC_GGA_XC_KT2", BRIAN_FUNCTIONAL_GGA_KT2_XC},
            {"XC_GGA_X_LAG", BRIAN_FUNCTIONAL_GGA_LAG_X},
            {"XC_GGA_X_LAMBDA_CH_N", BRIAN_FUNCTIONAL_GGA_LAMBDA_CH_N_X},
            {"XC_GGA_X_LAMBDA_LO_N", BRIAN_FUNCTIONAL_GGA_LAMBDA_LO_N_X},
            {"XC_GGA_X_LAMBDA_OC2_N", BRIAN_FUNCTIONAL_GGA_LAMBDA_OC2_N_X},
            {"XC_GGA_X_LBM", BRIAN_FUNCTIONAL_GGA_LBM_X},
            {"XC_GGA_X_LB", BRIAN_FUNCTIONAL_GGA_LB_X},
            {"XC_GGA_X_LG93", BRIAN_FUNCTIONAL_GGA_LG93_X},
            {"XC_GGA_C_LM", BRIAN_FUNCTIONAL_GGA_LM_C},
            {"XC_GGA_X_LV_RPW86", BRIAN_FUNCTIONAL_GGA_LV_RPW86_X},
            {"XC_GGA_C_LYP", BRIAN_FUNCTIONAL_GGA_LYP_C},
            {"XC_GGA_X_MB88", BRIAN_FUNCTIONAL_GGA_MB88_X},
            {"XC_GGA_XC_MOHLYP2", BRIAN_FUNCTIONAL_GGA_MOHLYP2_XC},
            {"XC_GGA_XC_MOHLYP", BRIAN_FUNCTIONAL_GGA_MOHLYP_XC},
            {"XC_GGA_X_MPBE", BRIAN_FUNCTIONAL_GGA_MPBE_X},
            {"XC_GGA_X_MPW91", BRIAN_FUNCTIONAL_GGA_MPW91_X},
            {"XC_GGA_XC_MPWLYP1W", BRIAN_FUNCTIONAL_GGA_MPWLYP1W_XC},
            {"XC_GGA_C_N12", BRIAN_FUNCTIONAL_GGA_N12_C},
            {"XC_GGA_C_N12_SX", BRIAN_FUNCTIONAL_GGA_N12_SX_C},
            {"XC_GGA_X_N12", BRIAN_FUNCTIONAL_GGA_N12_X},
            {"XC_GGA_XC_OBLYP_D", BRIAN_FUNCTIONAL_GGA_OBLYP_D_XC},
            {"XC_GGA_X_OL2", BRIAN_FUNCTIONAL_GGA_OL2_X},
            {"XC_GGA_XC_OPBE_D", BRIAN_FUNCTIONAL_GGA_OPBE_D_XC},
            {"XC_GGA_X_OPTB88_VDW", BRIAN_FUNCTIONAL_GGA_OPTB88_VDW_X},
            {"XC_GGA_C_OPTC", BRIAN_FUNCTIONAL_GGA_OPTC_C},
            {"XC_GGA_X_OPTPBE_VDW", BRIAN_FUNCTIONAL_GGA_OPTPBE_VDW_X},
            {"XC_GGA_X_OPTX", BRIAN_FUNCTIONAL_GGA_OPTX_X},
            {"XC_GGA_XC_OPWLYP_D", BRIAN_FUNCTIONAL_GGA_OPWLYP_D_XC},
            {"XC_GGA_C_OP_B88", BRIAN_FUNCTIONAL_GGA_OP_B88_C},
            {"XC_GGA_C_OP_G96", BRIAN_FUNCTIONAL_GGA_OP_G96_C},
            {"XC_GGA_C_OP_PBE", BRIAN_FUNCTIONAL_GGA_OP_PBE_C},
            {"XC_GGA_C_OP_PW91", BRIAN_FUNCTIONAL_GGA_OP_PW91_C},
            {"XC_GGA_C_OP_XALPHA", BRIAN_FUNCTIONAL_GGA_OP_XALPHA_C},
            {"XC_GGA_C_P86", BRIAN_FUNCTIONAL_GGA_P86_C},
            {"XC_GGA_XC_PBE1W", BRIAN_FUNCTIONAL_GGA_PBE1W_XC},
            {"XC_GGA_X_PBEA", BRIAN_FUNCTIONAL_GGA_PBEA_X},
            {"XC_GGA_C_PBEFE", BRIAN_FUNCTIONAL_GGA_PBEFE_C},
            {"XC_GGA_X_PBEFE", BRIAN_FUNCTIONAL_GGA_PBEFE_X},
            {"XC_GGA_C_PBEINT", BRIAN_FUNCTIONAL_GGA_PBEINT_C},
            {"XC_GGA_X_PBEINT", BRIAN_FUNCTIONAL_GGA_PBEINT_X},
            {"XC_GGA_X_PBEK1_VDW", BRIAN_FUNCTIONAL_GGA_PBEK1_VDW_X},
            {"XC_GGA_C_PBELOC", BRIAN_FUNCTIONAL_GGA_PBELOC_C},
            {"XC_GGA_XC_PBELYP1W", BRIAN_FUNCTIONAL_GGA_PBELYP1W_XC},
            {"XC_GGA_X_PBEPOW", BRIAN_FUNCTIONAL_GGA_PBEPOW_X},
            {"XC_GGA_X_PBETRANS", BRIAN_FUNCTIONAL_GGA_PBETRANS_X},
            {"XC_GGA_C_PBE", BRIAN_FUNCTIONAL_GGA_PBE_C},
            {"XC_GGA_C_PBE_JRGX", BRIAN_FUNCTIONAL_GGA_PBE_JRGX_C},
            {"XC_GGA_X_PBE_JSJR", BRIAN_FUNCTIONAL_GGA_PBE_JSJR_X},
            {"XC_GGA_C_PBE_MOL", BRIAN_FUNCTIONAL_GGA_PBE_MOL_C},
            {"XC_GGA_X_PBE_MOL", BRIAN_FUNCTIONAL_GGA_PBE_MOL_X},
            {"XC_GGA_X_PBE_R", BRIAN_FUNCTIONAL_GGA_PBE_R_X},
            {"XC_GGA_C_PBE_SOL", BRIAN_FUNCTIONAL_GGA_PBE_SOL_C},
            {"XC_GGA_X_PBE_SOL", BRIAN_FUNCTIONAL_GGA_PBE_SOL_X},
            {"XC_GGA_X_PBE_TCA", BRIAN_FUNCTIONAL_GGA_PBE_TCA_X},
            {"XC_GGA_X_PBE", BRIAN_FUNCTIONAL_GGA_PBE_X},
            {"XC_GGA_X_PW86", BRIAN_FUNCTIONAL_GGA_PW86_X},
            {"XC_GGA_C_PW91", BRIAN_FUNCTIONAL_GGA_PW91_C},
            {"XC_GGA_X_PW91", BRIAN_FUNCTIONAL_GGA_PW91_X},
            {"XC_GGA_C_Q2D", BRIAN_FUNCTIONAL_GGA_Q2D_C},
            {"XC_GGA_X_Q2D", BRIAN_FUNCTIONAL_GGA_Q2D_X},
            {"XC_GGA_C_REGTPSS", BRIAN_FUNCTIONAL_GGA_REGTPSS_C},
            {"XC_GGA_C_REVTCA", BRIAN_FUNCTIONAL_GGA_REVTCA_C},
            {"XC_GGA_C_RGE2", BRIAN_FUNCTIONAL_GGA_RGE2_C},
            {"XC_GGA_X_RGE2", BRIAN_FUNCTIONAL_GGA_RGE2_X},
            {"XC_GGA_X_RPBE", BRIAN_FUNCTIONAL_GGA_RPBE_X},
            {"XC_GGA_X_RPW86", BRIAN_FUNCTIONAL_GGA_RPW86_X},
            {"XC_GGA_C_SCAN_E0", BRIAN_FUNCTIONAL_GGA_SCAN_E0_C},
            {"XC_GGA_X_SFAT", BRIAN_FUNCTIONAL_GGA_SFAT_X},
            {"XC_GGA_C_SG4", BRIAN_FUNCTIONAL_GGA_SG4_C},
            {"XC_GGA_X_SG4", BRIAN_FUNCTIONAL_GGA_SG4_X},
            {"XC_GGA_C_SOGGA11", BRIAN_FUNCTIONAL_GGA_SOGGA11_C},
            {"XC_GGA_X_SOGGA11", BRIAN_FUNCTIONAL_GGA_SOGGA11_X},
            {"XC_GGA_C_SOGGA11_X", BRIAN_FUNCTIONAL_GGA_SOGGA11_X_C},
            {"XC_GGA_X_SOGGA", BRIAN_FUNCTIONAL_GGA_SOGGA_X},
            {"XC_GGA_C_SPBE", BRIAN_FUNCTIONAL_GGA_SPBE_C},
            {"XC_GGA_X_SSB_D", BRIAN_FUNCTIONAL_GGA_SSB_D_X},
            {"XC_GGA_X_SSB_SW", BRIAN_FUNCTIONAL_GGA_SSB_SW_X},
            {"XC_GGA_X_SSB", BRIAN_FUNCTIONAL_GGA_SSB_X},
            {"XC_GGA_C_TAU_HCTH", BRIAN_FUNCTIONAL_GGA_TAU_HCTH_C},
            {"XC_GGA_C_TCA", BRIAN_FUNCTIONAL_GGA_TCA_C},
            {"XC_GGA_XC_TH1", BRIAN_FUNCTIONAL_GGA_TH1_XC},
            {"XC_GGA_XC_TH2", BRIAN_FUNCTIONAL_GGA_TH2_XC},
            {"XC_GGA_XC_TH3", BRIAN_FUNCTIONAL_GGA_TH3_XC},
            {"XC_GGA_XC_TH4", BRIAN_FUNCTIONAL_GGA_TH4_XC},
            {"XC_GGA_XC_TH_FCFO", BRIAN_FUNCTIONAL_GGA_TH_FCFO_XC},
            {"XC_GGA_XC_TH_FCO", BRIAN_FUNCTIONAL_GGA_TH_FCO_XC},
            {"XC_GGA_XC_TH_FC", BRIAN_FUNCTIONAL_GGA_TH_FC_XC},
            {"XC_GGA_XC_TH_FL", BRIAN_FUNCTIONAL_GGA_TH_FL_XC},
            {"XC_GGA_C_TM_LYP", BRIAN_FUNCTIONAL_GGA_TM_LYP_C},
            {"XC_GGA_C_TM_PBE", BRIAN_FUNCTIONAL_GGA_TM_PBE_C},
            {"XC_GGA_X_VMT84_GE", BRIAN_FUNCTIONAL_GGA_VMT84_GE_X},
            {"XC_GGA_X_VMT84_PBE", BRIAN_FUNCTIONAL_GGA_VMT84_PBE_X},
            {"XC_GGA_X_VMT_GE", BRIAN_FUNCTIONAL_GGA_VMT_GE_X},
            {"XC_GGA_X_VMT_PBE", BRIAN_FUNCTIONAL_GGA_VMT_PBE_X},
            {"XC_GGA_XC_VV10", BRIAN_FUNCTIONAL_GGA_VV10_XC},
            {"XC_GGA_C_W94", BRIAN_FUNCTIONAL_GGA_W94_C},
            {"XC_GGA_X_WC", BRIAN_FUNCTIONAL_GGA_WC_X},
            {"XC_GGA_C_WI0", BRIAN_FUNCTIONAL_GGA_WI0_C},
            {"XC_GGA_C_WI", BRIAN_FUNCTIONAL_GGA_WI_C},
            {"XC_GGA_C_WL", BRIAN_FUNCTIONAL_GGA_WL_C},
            {"XC_GGA_X_WPBEH", BRIAN_FUNCTIONAL_GGA_WPBEH_X},
            {"XC_GGA_XC_XLYP", BRIAN_FUNCTIONAL_GGA_XLYP_XC},
            {"XC_GGA_C_XPBE", BRIAN_FUNCTIONAL_GGA_XPBE_C},
            {"XC_GGA_X_XPBE", BRIAN_FUNCTIONAL_GGA_XPBE_X},
            {"XC_GGA_C_ZPBEINT", BRIAN_FUNCTIONAL_GGA_ZPBEINT_C},
            {"XC_GGA_C_ZPBESOL", BRIAN_FUNCTIONAL_GGA_ZPBESOL_C},
            {"XC_GGA_C_ZVPBEINT", BRIAN_FUNCTIONAL_GGA_ZVPBEINT_C},
            {"XC_GGA_C_ZVPBESOL", BRIAN_FUNCTIONAL_GGA_ZVPBESOL_C},
            {"XC_HYB_GGA_XC_B1LYP", BRIAN_FUNCTIONAL_HGGA_B1LYP_XC},
            {"XC_HYB_GGA_XC_B1PW91", BRIAN_FUNCTIONAL_HGGA_B1PW91_XC},
            {"XC_HYB_GGA_XC_B1WC", BRIAN_FUNCTIONAL_HGGA_B1WC_XC},
            {"XC_HYB_GGA_XC_B3LYP5", BRIAN_FUNCTIONAL_HGGA_B3LYP5_XC},
            {"XC_HYB_GGA_XC_B3LYPS", BRIAN_FUNCTIONAL_HGGA_B3LYPS_XC},
            {"XC_HYB_GGA_XC_B3LYP", BRIAN_FUNCTIONAL_HGGA_B3LYP_XC},
            {"XC_HYB_GGA_XC_B3P86", BRIAN_FUNCTIONAL_HGGA_B3P86_XC},
            {"XC_HYB_GGA_XC_B3PW91", BRIAN_FUNCTIONAL_HGGA_B3PW91_XC},
            {"XC_HYB_GGA_XC_B5050LYP", BRIAN_FUNCTIONAL_HGGA_B5050LYP_XC},
            {"XC_HYB_GGA_XC_B97_1P", BRIAN_FUNCTIONAL_HGGA_B97_1P_XC},
            {"XC_HYB_GGA_XC_B97_1", BRIAN_FUNCTIONAL_HGGA_B97_1_XC},
            {"XC_HYB_GGA_XC_B97_2", BRIAN_FUNCTIONAL_HGGA_B97_2_XC},
            {"XC_HYB_GGA_XC_B97_3", BRIAN_FUNCTIONAL_HGGA_B97_3_XC},
            {"XC_HYB_GGA_XC_B97_K", BRIAN_FUNCTIONAL_HGGA_B97_K_XC},
            {"XC_HYB_GGA_XC_B97", BRIAN_FUNCTIONAL_HGGA_B97_XC},
            {"XC_HYB_GGA_XC_BHANDHLYP", BRIAN_FUNCTIONAL_HGGA_BHANDHLYP_XC},
            {"XC_HYB_GGA_XC_BHANDH", BRIAN_FUNCTIONAL_HGGA_BHANDH_XC},
            {"XC_HYB_GGA_XC_CAMY_B3LYP", BRIAN_FUNCTIONAL_HGGA_CAMY_B3LYP_XC},
            {"XC_HYB_GGA_XC_CAMY_BLYP", BRIAN_FUNCTIONAL_HGGA_CAMY_BLYP_XC},
            {"XC_HYB_GGA_XC_CAM_B3LYP", BRIAN_FUNCTIONAL_HGGA_CAM_B3LYP_XC},
            {"XC_HYB_GGA_XC_CAM_QTP_01", BRIAN_FUNCTIONAL_HGGA_CAM_QTP_01_XC},
            {"XC_HYB_GGA_XC_CAP0", BRIAN_FUNCTIONAL_HGGA_CAP0_XC},
            {"XC_HYB_GGA_XC_EDF2", BRIAN_FUNCTIONAL_HGGA_EDF2_XC},
            {"XC_HYB_GGA_XC_HJS_B88", BRIAN_FUNCTIONAL_HGGA_HJS_B88_XC},
            {"XC_HYB_GGA_XC_HJS_B97X", BRIAN_FUNCTIONAL_HGGA_HJS_B97X_XC},
            {"XC_HYB_GGA_XC_HJS_PBE_SOL", BRIAN_FUNCTIONAL_HGGA_HJS_PBE_SOL_XC},
            {"XC_HYB_GGA_XC_HJS_PBE", BRIAN_FUNCTIONAL_HGGA_HJS_PBE_XC},
            {"XC_HYB_GGA_XC_HSE03", BRIAN_FUNCTIONAL_HGGA_HSE03_XC},
            {"XC_HYB_GGA_XC_HSE06", BRIAN_FUNCTIONAL_HGGA_HSE06_XC},
            {"XC_HYB_GGA_XC_HSE12S", BRIAN_FUNCTIONAL_HGGA_HSE12S_XC},
            {"XC_HYB_GGA_XC_HSE12", BRIAN_FUNCTIONAL_HGGA_HSE12_XC},
            {"XC_HYB_GGA_XC_HSE_SOL", BRIAN_FUNCTIONAL_HGGA_HSE_SOL_XC},
            {"XC_HYB_GGA_XC_KMLYP", BRIAN_FUNCTIONAL_HGGA_KMLYP_XC},
            {"XC_HYB_GGA_XC_LCY_BLYP", BRIAN_FUNCTIONAL_HGGA_LCY_BLYP_XC},
            {"XC_HYB_GGA_XC_LCY_PBE", BRIAN_FUNCTIONAL_HGGA_LCY_PBE_XC},
            {"XC_HYB_GGA_XC_LC_VV10", BRIAN_FUNCTIONAL_HGGA_LC_VV10_XC},
            {"XC_HYB_GGA_XC_LC_WPBE", BRIAN_FUNCTIONAL_HGGA_LC_WPBE_XC},
            {"XC_HYB_GGA_XC_LRC_WPBEH", BRIAN_FUNCTIONAL_HGGA_LRC_WPBEH_XC},
            {"XC_HYB_GGA_XC_LRC_WPBE", BRIAN_FUNCTIONAL_HGGA_LRC_WPBE_XC},
            {"XC_HYB_GGA_XC_MB3LYP_RC04", BRIAN_FUNCTIONAL_HGGA_MB3LYP_RC04_XC},
            {"XC_HYB_GGA_XC_MPW1K", BRIAN_FUNCTIONAL_HGGA_MPW1K_XC},
            {"XC_HYB_GGA_XC_MPW1LYP", BRIAN_FUNCTIONAL_HGGA_MPW1LYP_XC},
            {"XC_HYB_GGA_XC_MPW1PBE", BRIAN_FUNCTIONAL_HGGA_MPW1PBE_XC},
            {"XC_HYB_GGA_XC_MPW1PW", BRIAN_FUNCTIONAL_HGGA_MPW1PW_XC},
            {"XC_HYB_GGA_XC_MPW3LYP", BRIAN_FUNCTIONAL_HGGA_MPW3LYP_XC},
            {"XC_HYB_GGA_XC_MPW3PW", BRIAN_FUNCTIONAL_HGGA_MPW3PW_XC},
            {"XC_HYB_GGA_XC_MPWLYP1M", BRIAN_FUNCTIONAL_HGGA_MPWLYP1M_XC},
            {"XC_HYB_GGA_X_N12_SX", BRIAN_FUNCTIONAL_HGGA_N12_SX_X},
            {"XC_HYB_GGA_XC_O3LYP", BRIAN_FUNCTIONAL_HGGA_O3LYP_XC},
            {"XC_HYB_GGA_XC_PBE0_13", BRIAN_FUNCTIONAL_HGGA_PBE0_13_XC},
            {"XC_HYB_GGA_XC_PBE50", BRIAN_FUNCTIONAL_HGGA_PBE50_XC},
            {"XC_HYB_GGA_XC_PBEB0", BRIAN_FUNCTIONAL_HGGA_PBEB0_XC},
            {"XC_HYB_GGA_XC_PBEH", BRIAN_FUNCTIONAL_HGGA_PBEH_XC},
            {"XC_HYB_GGA_XC_PBE_MOL0", BRIAN_FUNCTIONAL_HGGA_PBE_MOL0_XC},
            {"XC_HYB_GGA_XC_PBE_MOLB0", BRIAN_FUNCTIONAL_HGGA_PBE_MOLB0_XC},
            {"XC_HYB_GGA_XC_PBE_SOL0", BRIAN_FUNCTIONAL_HGGA_PBE_SOL0_XC},
            {"XC_HYB_GGA_XC_REVB3LYP", BRIAN_FUNCTIONAL_HGGA_REVB3LYP_XC},
            {"XC_HYB_GGA_XC_SB98_1A", BRIAN_FUNCTIONAL_HGGA_SB98_1A_XC},
            {"XC_HYB_GGA_XC_SB98_1B", BRIAN_FUNCTIONAL_HGGA_SB98_1B_XC},
            {"XC_HYB_GGA_XC_SB98_1C", BRIAN_FUNCTIONAL_HGGA_SB98_1C_XC},
            {"XC_HYB_GGA_XC_SB98_2A", BRIAN_FUNCTIONAL_HGGA_SB98_2A_XC},
            {"XC_HYB_GGA_XC_SB98_2B", BRIAN_FUNCTIONAL_HGGA_SB98_2B_XC},
            {"XC_HYB_GGA_XC_SB98_2C", BRIAN_FUNCTIONAL_HGGA_SB98_2C_XC},
            {"XC_HYB_GGA_X_SOGGA11_X", BRIAN_FUNCTIONAL_HGGA_SOGGA11_X_X},
            {"XC_HYB_GGA_XC_TUNED_CAM_B3LYP", BRIAN_FUNCTIONAL_HGGA_TUNED_CAM_B3LYP_XC},
            {"XC_HYB_GGA_XC_WB97X_D", BRIAN_FUNCTIONAL_HGGA_WB97X_D_XC},
            {"XC_HYB_GGA_XC_WB97X_V", BRIAN_FUNCTIONAL_HGGA_WB97X_V_XC},
            {"XC_HYB_GGA_XC_WB97X", BRIAN_FUNCTIONAL_HGGA_WB97X_XC},
            {"XC_HYB_GGA_XC_WB97", BRIAN_FUNCTIONAL_HGGA_WB97_XC},
            {"XC_HYB_GGA_XC_X3LYP", BRIAN_FUNCTIONAL_HGGA_X3LYP_XC},
            {"XC_HYB_MGGA_XC_B86B95", BRIAN_FUNCTIONAL_HMGGA_B86B95_XC},
            {"XC_HYB_MGGA_XC_B88B95", BRIAN_FUNCTIONAL_HMGGA_B88B95_XC},
            {"XC_HYB_MGGA_XC_BB1K", BRIAN_FUNCTIONAL_HMGGA_BB1K_XC},
            {"XC_HYB_MGGA_X_BMK", BRIAN_FUNCTIONAL_HMGGA_BMK_X},
            {"XC_HYB_MGGA_X_DLDF", BRIAN_FUNCTIONAL_HMGGA_DLDF_X},
            {"XC_HYB_MGGA_X_M05_2X", BRIAN_FUNCTIONAL_HMGGA_M05_2X_X},
            {"XC_HYB_MGGA_X_M05", BRIAN_FUNCTIONAL_HMGGA_M05_X},
            {"XC_HYB_MGGA_X_M06_2X", BRIAN_FUNCTIONAL_HMGGA_M06_2X_X},
            {"XC_HYB_MGGA_X_M06_HF", BRIAN_FUNCTIONAL_HMGGA_M06_HF_X},
            {"XC_HYB_MGGA_X_M06", BRIAN_FUNCTIONAL_HMGGA_M06_X},
            {"XC_HYB_MGGA_X_M08_HX", BRIAN_FUNCTIONAL_HMGGA_M08_HX_X},
            {"XC_HYB_MGGA_X_M08_SO", BRIAN_FUNCTIONAL_HMGGA_M08_SO_X},
            {"XC_HYB_MGGA_X_M11", BRIAN_FUNCTIONAL_HMGGA_M11_X},
            {"XC_HYB_MGGA_X_MN12_SX", BRIAN_FUNCTIONAL_HMGGA_MN12_SX_X},
            {"XC_HYB_MGGA_X_MN15", BRIAN_FUNCTIONAL_HMGGA_MN15_X},
            {"XC_HYB_MGGA_XC_MPW1B95", BRIAN_FUNCTIONAL_HMGGA_MPW1B95_XC},
            {"XC_HYB_MGGA_XC_MPWB1K", BRIAN_FUNCTIONAL_HMGGA_MPWB1K_XC},
            {"XC_HYB_MGGA_X_MS2H", BRIAN_FUNCTIONAL_HMGGA_MS2H_X},
            {"XC_HYB_MGGA_X_MVSH", BRIAN_FUNCTIONAL_HMGGA_MVSH_X},
            {"XC_HYB_MGGA_XC_PW6B95", BRIAN_FUNCTIONAL_HMGGA_PW6B95_XC},
            {"XC_HYB_MGGA_XC_PW86B95", BRIAN_FUNCTIONAL_HMGGA_PW86B95_XC},
            {"XC_HYB_MGGA_XC_PWB6K", BRIAN_FUNCTIONAL_HMGGA_PWB6K_XC},
            {"XC_HYB_MGGA_X_REVSCAN0", BRIAN_FUNCTIONAL_HMGGA_REVSCAN0_X},
            {"XC_HYB_MGGA_XC_REVTPSSH", BRIAN_FUNCTIONAL_HMGGA_REVTPSSH_XC},
            {"XC_HYB_MGGA_X_SCAN0", BRIAN_FUNCTIONAL_HMGGA_SCAN0_X},
            {"XC_HYB_MGGA_X_TAU_HCTH", BRIAN_FUNCTIONAL_HMGGA_TAU_HCTH_X},
            {"XC_HYB_MGGA_XC_TPSSH", BRIAN_FUNCTIONAL_HMGGA_TPSSH_XC},
            {"XC_HYB_MGGA_XC_WB97M_V", BRIAN_FUNCTIONAL_HMGGA_WB97M_V_XC},
            {"XC_HYB_MGGA_XC_X1B95", BRIAN_FUNCTIONAL_HMGGA_X1B95_XC},
            {"XC_HYB_MGGA_XC_XB1K", BRIAN_FUNCTIONAL_HMGGA_XB1K_XC},
            {"XC_LDA_C_BR78", BRIAN_FUNCTIONAL_LDA_BR78_C},
            {"XC_LDA_C_CHACHIYO", BRIAN_FUNCTIONAL_LDA_CHACHIYO_C},
            {"XC_LDA_X_ERF", BRIAN_FUNCTIONAL_LDA_ERF_X},
            {"XC_LDA_XC_GDSMFB", BRIAN_FUNCTIONAL_LDA_GDSMFB_XC},
            {"XC_LDA_C_GK72", BRIAN_FUNCTIONAL_LDA_GK72_C},
            {"XC_LDA_C_GL", BRIAN_FUNCTIONAL_LDA_GL_C},
            {"XC_LDA_C_GOMBAS", BRIAN_FUNCTIONAL_LDA_GOMBAS_C},
            {"XC_LDA_C_HL", BRIAN_FUNCTIONAL_LDA_HL_C},
            {"XC_LDA_C_KARASIEV", BRIAN_FUNCTIONAL_LDA_KARASIEV_C},
            {"XC_LDA_XC_KSDT", BRIAN_FUNCTIONAL_LDA_KSDT_XC},
            {"XC_LDA_C_LP96", BRIAN_FUNCTIONAL_LDA_LP96_C},
            {"XC_LDA_XC_LP_A", BRIAN_FUNCTIONAL_LDA_LP_A_XC},
            {"XC_LDA_XC_LP_B", BRIAN_FUNCTIONAL_LDA_LP_B_XC},
            {"XC_LDA_C_MCWEENY", BRIAN_FUNCTIONAL_LDA_MCWEENY_C},
            {"XC_LDA_C_ML1", BRIAN_FUNCTIONAL_LDA_ML1_C},
            {"XC_LDA_C_ML2", BRIAN_FUNCTIONAL_LDA_ML2_C},
            {"XC_LDA_C_OB_PW", BRIAN_FUNCTIONAL_LDA_OB_PW_C},
            {"XC_LDA_C_OB_PZ", BRIAN_FUNCTIONAL_LDA_OB_PZ_C},
            {"XC_LDA_C_OW", BRIAN_FUNCTIONAL_LDA_OW_C},
            {"XC_LDA_C_OW_LYP", BRIAN_FUNCTIONAL_LDA_OW_LYP_C},
            {"XC_LDA_C_PK09", BRIAN_FUNCTIONAL_LDA_PK09_C},
            {"XC_LDA_C_PW", BRIAN_FUNCTIONAL_LDA_PW_C},
            {"XC_LDA_C_PW_MOD", BRIAN_FUNCTIONAL_LDA_PW_MOD_C},
            {"XC_LDA_C_PW_RPA", BRIAN_FUNCTIONAL_LDA_PW_RPA_C},
            {"XC_LDA_C_PZ", BRIAN_FUNCTIONAL_LDA_PZ_C},
            {"XC_LDA_C_PZ_MOD", BRIAN_FUNCTIONAL_LDA_PZ_MOD_C},
            {"XC_LDA_X_RAE", BRIAN_FUNCTIONAL_LDA_RAE_X},
            {"XC_LDA_C_RC04", BRIAN_FUNCTIONAL_LDA_RC04_C},
            {"XC_LDA_X_REL", BRIAN_FUNCTIONAL_LDA_REL_X},
            {"XC_LDA_C_RPA", BRIAN_FUNCTIONAL_LDA_RPA_C},
            {"XC_LDA_X", BRIAN_FUNCTIONAL_LDA_SLATER_X},
            {"XC_LDA_XC_TETER93", BRIAN_FUNCTIONAL_LDA_TETER93_XC},
            {"XC_LDA_C_VBH", BRIAN_FUNCTIONAL_LDA_VBH_C},
            {"XC_LDA_C_VWN_RPA", BRIAN_FUNCTIONAL_LDA_VWN1RPA_C},
            {"XC_LDA_C_VWN_1", BRIAN_FUNCTIONAL_LDA_VWN1_C},
            {"XC_LDA_C_VWN_RPA", BRIAN_FUNCTIONAL_LDA_VWN5RPA_C},
            {"XC_LDA_C_VWN", BRIAN_FUNCTIONAL_LDA_VWN5_C},
            {"XC_LDA_C_VWN_1", BRIAN_FUNCTIONAL_LDA_VWN_1_C},
            {"XC_LDA_C_VWN_2", BRIAN_FUNCTIONAL_LDA_VWN_2_C},
            {"XC_LDA_C_VWN_3", BRIAN_FUNCTIONAL_LDA_VWN_3_C},
            {"XC_LDA_C_VWN_4", BRIAN_FUNCTIONAL_LDA_VWN_4_C},
            {"XC_LDA_C_VWN", BRIAN_FUNCTIONAL_LDA_VWN_C},
            {"XC_LDA_C_VWN_RPA", BRIAN_FUNCTIONAL_LDA_VWN_RPA_C},
            {"XC_LDA_C_WIGNER", BRIAN_FUNCTIONAL_LDA_WIGNER_C},
            {"XC_LDA_C_XALPHA", BRIAN_FUNCTIONAL_LDA_XALPHA_C},
            {"XC_LDA_XC_ZLP", BRIAN_FUNCTIONAL_LDA_ZLP_XC},
            {"XC_MGGA_C_B88", BRIAN_FUNCTIONAL_MGGA_B88_C},
            {"XC_MGGA_XC_B97M_V", BRIAN_FUNCTIONAL_MGGA_B97M_V_XC},
            {"XC_MGGA_C_BC95", BRIAN_FUNCTIONAL_MGGA_BC95_C},
            {"XC_MGGA_X_BLOC", BRIAN_FUNCTIONAL_MGGA_BLOC_X},
            {"XC_MGGA_C_DLDF", BRIAN_FUNCTIONAL_MGGA_DLDF_C},
            {"XC_MGGA_X_GVT4", BRIAN_FUNCTIONAL_MGGA_GVT4_X},
            {"XC_MGGA_X_GX", BRIAN_FUNCTIONAL_MGGA_GX_X},
            {"XC_MGGA_XC_HLE17", BRIAN_FUNCTIONAL_MGGA_HLE17_XC},
            {"XC_MGGA_X_LTA", BRIAN_FUNCTIONAL_MGGA_LTA_X},
            {"XC_MGGA_C_M05_2X", BRIAN_FUNCTIONAL_MGGA_M05_2X_C},
            {"XC_MGGA_C_M05", BRIAN_FUNCTIONAL_MGGA_M05_C},
            {"XC_MGGA_C_M06_2X", BRIAN_FUNCTIONAL_MGGA_M06_2X_C},
            {"XC_MGGA_C_M06", BRIAN_FUNCTIONAL_MGGA_M06_C},
            {"XC_MGGA_C_M06_HF", BRIAN_FUNCTIONAL_MGGA_M06_HF_C},
            {"XC_MGGA_C_M06_L", BRIAN_FUNCTIONAL_MGGA_M06_L_C},
            {"XC_MGGA_X_M06_L", BRIAN_FUNCTIONAL_MGGA_M06_L_X},
            {"XC_MGGA_C_M08_HX", BRIAN_FUNCTIONAL_MGGA_M08_HX_C},
            {"XC_MGGA_C_M08_SO", BRIAN_FUNCTIONAL_MGGA_M08_SO_C},
            {"XC_MGGA_C_M11", BRIAN_FUNCTIONAL_MGGA_M11_C},
            {"XC_MGGA_C_M11_L", BRIAN_FUNCTIONAL_MGGA_M11_L_C},
            {"XC_MGGA_X_M11_L", BRIAN_FUNCTIONAL_MGGA_M11_L_X},
            {"XC_MGGA_X_MBEEFVDW", BRIAN_FUNCTIONAL_MGGA_MBEEFVDW_X},
            {"XC_MGGA_X_MBEEF", BRIAN_FUNCTIONAL_MGGA_MBEEF_X},
            {"XC_MGGA_C_MN12_L", BRIAN_FUNCTIONAL_MGGA_MN12_L_C},
            {"XC_MGGA_X_MN12_L", BRIAN_FUNCTIONAL_MGGA_MN12_L_X},
            {"XC_MGGA_C_MN12_SX", BRIAN_FUNCTIONAL_MGGA_MN12_SX_C},
            {"XC_MGGA_C_MN15", BRIAN_FUNCTIONAL_MGGA_MN15_C},
            {"XC_MGGA_C_MN15_L", BRIAN_FUNCTIONAL_MGGA_MN15_L_C},
            {"XC_MGGA_X_MN15_L", BRIAN_FUNCTIONAL_MGGA_MN15_L_X},
            {"XC_MGGA_X_MODTPSS", BRIAN_FUNCTIONAL_MGGA_MODTPSS_X},
            {"XC_MGGA_X_MS0", BRIAN_FUNCTIONAL_MGGA_MS0_X},
            {"XC_MGGA_X_MS1", BRIAN_FUNCTIONAL_MGGA_MS1_X},
            {"XC_MGGA_X_MS2", BRIAN_FUNCTIONAL_MGGA_MS2_X},
            {"XC_MGGA_X_MVS", BRIAN_FUNCTIONAL_MGGA_MVS_X},
            {"XC_MGGA_XC_OTPSS_D", BRIAN_FUNCTIONAL_MGGA_OTPSS_D_XC},
            {"XC_MGGA_X_PBE_GX", BRIAN_FUNCTIONAL_MGGA_PBE_GX_X},
            {"XC_MGGA_C_PKZB", BRIAN_FUNCTIONAL_MGGA_PKZB_C},
            {"XC_MGGA_X_PKZB", BRIAN_FUNCTIONAL_MGGA_PKZB_X},
            {"XC_MGGA_C_REVM06_L", BRIAN_FUNCTIONAL_MGGA_REVM06_L_C},
            {"XC_MGGA_X_REVM06_L", BRIAN_FUNCTIONAL_MGGA_REVM06_L_X},
            {"XC_MGGA_C_REVSCAN", BRIAN_FUNCTIONAL_MGGA_REVSCAN_C},
            {"XC_MGGA_C_REVSCAN_VV10", BRIAN_FUNCTIONAL_MGGA_REVSCAN_VV10_C},
            {"XC_MGGA_X_REVSCAN", BRIAN_FUNCTIONAL_MGGA_REVSCAN_X},
            {"XC_MGGA_C_REVTPSS", BRIAN_FUNCTIONAL_MGGA_REVTPSS_C},
            {"XC_MGGA_X_REVTPSS", BRIAN_FUNCTIONAL_MGGA_REVTPSS_X},
            {"XC_MGGA_X_SA_TPSS", BRIAN_FUNCTIONAL_MGGA_SA_TPSS_X},
            {"XC_MGGA_C_SCAN", BRIAN_FUNCTIONAL_MGGA_SCAN_C},
            {"XC_MGGA_C_SCAN_RVV10", BRIAN_FUNCTIONAL_MGGA_SCAN_RVV10_C},
            {"XC_MGGA_C_SCAN_VV10", BRIAN_FUNCTIONAL_MGGA_SCAN_VV10_C},
            {"XC_MGGA_X_SCAN", BRIAN_FUNCTIONAL_MGGA_SCAN_X},
            {"XC_MGGA_X_TAU_HCTH", BRIAN_FUNCTIONAL_MGGA_TAU_HCTH_X},
            {"XC_MGGA_X_TM", BRIAN_FUNCTIONAL_MGGA_TM_X},
            {"XC_MGGA_C_TPSSLOC", BRIAN_FUNCTIONAL_MGGA_TPSSLOC_C},
            {"XC_MGGA_XC_TPSSLYP1W", BRIAN_FUNCTIONAL_MGGA_TPSSLYP1W_XC},
            {"XC_MGGA_C_TPSS", BRIAN_FUNCTIONAL_MGGA_TPSS_C},
            {"XC_MGGA_X_TPSS", BRIAN_FUNCTIONAL_MGGA_TPSS_X},
            {"XC_MGGA_C_VSXC", BRIAN_FUNCTIONAL_MGGA_VSXC_C},
            {"XC_MGGA_X_VT84", BRIAN_FUNCTIONAL_MGGA_VT84_X},
        };
        
        std::vector<brianInt> functionalIDs;
        std::vector<double> functionalWeights;
        for (std::shared_ptr<Functional> functionalComponent: functional_->x_functionals()) {
            if (functionalIDMap.count(functionalComponent->name()) == 0) {
                throw PSIEXCEPTION("This DFT functional cannot be handled by BrianQC");
            }
            functionalIDs.push_back(functionalIDMap.at(functionalComponent->name()));
            functionalWeights.push_back(functionalComponent->alpha());
        }
        for (std::shared_ptr<Functional> functionalComponent: functional_->c_functionals()) {
            if (functionalIDMap.count(functionalComponent->name()) == 0) {
                throw PSIEXCEPTION("This DFT functional cannot be handled by BrianQC");
            }
            functionalIDs.push_back(functionalIDMap.at(functionalComponent->name()));
            functionalWeights.push_back(functionalComponent->alpha());
        }
        
        static const std::map<std::string, brianInt> functionalParameterIDMap = {
            // TODO currently, Brian doesn't handle any functional parameters
        };
        
        std::map<brianInt, double> functionalParameterMap;
        for (std::shared_ptr<Functional> functionalComponent: functional_->x_functionals()) {
            for (const std::pair<std::string, double>& parameter: functionalComponent->parameters()) {
                if (functionalParameterIDMap.count(parameter.first) == 0) {
                    throw PSIEXCEPTION("This DFT functional parameter cannot be handled by BrianQC");
                }
                
                brianInt functionalParameterID = functionalParameterIDMap.at(parameter.first);
                
                if (functionalParameterMap.count(functionalParameterID)) {
                    if (functionalParameterMap.at(functionalParameterID) != parameter.second) {
                        throw PSIEXCEPTION("BrianQC cannot handle different values of the same parameter for different DFT functional components");
                    }
                }
                else {
                    functionalParameterMap.insert({functionalParameterID, parameter.second});
                }
            }
        }
        
        if (options_.exists("DFT_OMEGA") and options_.get_double("DFT_OMEGA") > 0) {
            functionalParameterMap.insert({BRIAN_FUNCTIONAL_PARAMETER_RANGE, options_.get_double("DFT_OMEGA")});
        }
        
        std::vector<brianInt> functionalParameterIDs;
        std::vector<double> functionalParameterValues;
        for (const std::pair<brianInt, double>& parameter: functionalParameterMap) {
            functionalParameterIDs.push_back(parameter.first);
            functionalParameterValues.push_back(parameter.second);
        }
        
        brianInt functionalCount = functionalIDs.size();
        brianInt functionalParameterCount = functionalParameterIDs.size();
        brianCOMSetDFTFunctional(&brianCookie,
            &functionalCount,
            functionalIDs.data(),
            functionalWeights.data(),
            &functionalParameterCount,
            functionalParameterIDs.data(),
            functionalParameterValues.data()
        );
        checkBrian();
        
        if (functional_->needs_vv10()) {
            // Psi4 would only generate the NLC grid when actually computing the NLC term
            // (and that code path is not even called when BrianQC is active),
            // but we need the NLC grid before initializing, so we replicate the grid building here
            {
                std::map<std::string, std::string> opt_map;
                opt_map["DFT_PRUNING_SCHEME"] = "FLAT";
                
                std::map<std::string, int> opt_int_map;
                opt_int_map["DFT_RADIAL_POINTS"] = options_.get_int("DFT_VV10_RADIAL_POINTS");
                opt_int_map["DFT_SPHERICAL_POINTS"] = options_.get_int("DFT_VV10_SPHERICAL_POINTS");
                
                brianBuildingNLCGrid = true;
                DFTGrid nlgrid = DFTGrid(primary_->molecule(), primary_, opt_int_map, opt_map, options_);
                brianBuildingNLCGrid = false;
            }
            
            std::vector<brianInt> NLCParameterIDs;
            std::vector<double> NLCParameterValues;
            if (options_.exists("DFT_VV10_B") and options_.get_double("DFT_VV10_B") > 0) {
                NLCParameterIDs.push_back(BRIAN_NLC_PARAMETER_VV_B);
                NLCParameterValues.push_back(options_.get_double("DFT_VV10_B"));
            }
            if (options_.exists("DFT_VV10_C") and options_.get_double("DFT_VV10_C") > 0) {
                NLCParameterIDs.push_back(BRIAN_NLC_PARAMETER_VV_C);
                NLCParameterValues.push_back(options_.get_double("DFT_VV10_C"));
            }
            
            if (not NLCParameterIDs.empty()) {
                brianInt NLCID = BRIAN_NLC_VV10;
                double NLCWeight = 1.0;
                brianInt NLCParameterCount = NLCParameterIDs.size();
                brianCOMSetNLC(&brianCookie, &NLCID, &NLCWeight, &NLCParameterCount, NLCParameterIDs.data(), NLCParameterValues.data());
            }
        }
        
        // BrianQC's DFT takes the required precision into account
        double SCFConvergenceThreshold = std::min(options_.get_double("E_CONVERGENCE"), options_.get_double("D_CONVERGENCE"));
        brianSCFSetThresholds(&brianCookie, &SCFConvergenceThreshold);
        
        brianCOMInitDFT(&brianCookie);
        checkBrian();
}

std::map<std::string, double> BrianRV::compute_V(std::vector<SharedMatrix> ret) {
    double DFTEnergy;
    
    brianSCFBuildFockDFT(&brianCookie,
        D_AO_[0]->get_pointer(0),
        nullptr,
        ret[0]->get_pointer(0),
        nullptr,
        &DFTEnergy
    );
    checkBrian();
    
    std::map<std::string, double> quad_values;
    quad_values["VV10"] = 0.0; // NOTE: BrianQC doesn't compute the VV10 term separately, it just includes it in the DFT energy term
    quad_values["FUNCTIONAL"] = DFTEnergy;
    quad_values["RHO_A"] = 0.0;
    quad_values["RHO_AX"] = 0.0;
    quad_values["RHO_AY"] = 0.0;
    quad_values["RHO_AZ"] = 0.0;
    quad_values["RHO_B"] = 0.0;
    quad_values["RHO_BX"] = 0.0;
    quad_values["RHO_BY"] = 0.0;
    quad_values["RHO_BZ"] = 0.0;
    return quad_values;
}
std::map<std::string, double> BrianUV::compute_V(std::vector<SharedMatrix> ret) {
    double DFTEnergy;
    brianSCFBuildFockDFT(&brianCookie,
        D_AO_[0]->get_pointer(0),
        D_AO_[1]->get_pointer(0),
        ret[0]->get_pointer(0),
        ret[1]->get_pointer(0),
        &DFTEnergy
    );
    checkBrian();
        
    std::map<std::string, double> quad_values;
    quad_values["VV10"] = 0.0; // NOTE: BrianQC doesn't compute the VV10 term separately, it just includes it in the DFT energy term
    quad_values["FUNCTIONAL"] = DFTEnergy;
    quad_values["RHO_A"] = 0.0;
    quad_values["RHO_AX"] = 0.0;
    quad_values["RHO_AY"] = 0.0;
    quad_values["RHO_AZ"] = 0.0;
    quad_values["RHO_B"] = 0.0;
    quad_values["RHO_BX"] = 0.0;
    quad_values["RHO_BY"] = 0.0;
    quad_values["RHO_BZ"] = 0.0;
    return quad_values;
}


} // namespace psi

#endif
