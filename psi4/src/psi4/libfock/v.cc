/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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

#include "v.h"
#include "cubature.h"
#include "points.h"
#include "dft_integrators.h"
#include "sap.h"

#include "psi4/libfunctional/LibXCfunctional.h"
#include "psi4/libfunctional/functional.h"
#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libqt/qt.h"
#include "psi4/psi4-dec.h"

#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/petitelist.h"
#include "psi4/libmints/vector.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"

#include <cstdlib>
#include <numeric>
#include <sstream>
#include <string>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef USING_BrianQC

#include <use_brian_wrapper.h>
#include <brian_macros.h>
#include <brian_common.h>
#include <brian_scf.h>

extern void checkBrian();
extern BrianCookie brianCookie;
extern bool brianEnable;
extern bool brianEnableDFT;

struct brianGrid {
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

#endif

namespace psi {

VBase::VBase(std::shared_ptr<SuperFunctional> functional, std::shared_ptr<BasisSet> primary, Options& options)
    : options_(options), primary_(primary), functional_(functional) {
    common_init();
}
VBase::~VBase() {}
void VBase::common_init() {
    print_ = options_.get_int("PRINT");
    debug_ = options_.get_int("DEBUG");
    v2_rho_cutoff_ = options_.get_double("DFT_V2_RHO_CUTOFF");
    vv10_rho_cutoff_ = options_.get_double("DFT_VV10_RHO_CUTOFF");
    grac_initialized_ = false;
    cache_map_deriv_ = -1;
    num_threads_ = 1;
#ifdef _OPENMP
    num_threads_ = omp_get_max_threads();
#endif
}
std::shared_ptr<VBase> VBase::build_V(std::shared_ptr<BasisSet> primary, std::shared_ptr<SuperFunctional> functional,
                                      Options& options, const std::string& type) {
    std::shared_ptr<VBase> v;
    if (type == "RV") {
        if (!functional->is_unpolarized()) {
            throw PSIEXCEPTION("Passed in functional was polarized for RV reference.");
        }
        v = std::make_shared<RV>(functional, primary, options);
    } else if (type == "UV") {
        if (functional->is_unpolarized()) {
            throw PSIEXCEPTION("Passed in functional was unpolarized for UV reference.");
        }
        v = std::make_shared<UV>(functional, primary, options);
    } else if (type == "SAP") {
        v = std::make_shared<SAP>(functional, primary, options);
    } else {
        throw PSIEXCEPTION("V: V type is not recognized");
    }

    return v;
}
void VBase::set_D(std::vector<SharedMatrix> Dvec) {
    if (Dvec.size() > 2) {
        throw PSIEXCEPTION("VBase::set_D: Can only set up to two D vectors.");
    }

    // Build AO2USO matrix, if needed
    if (!AO2USO_ && (Dvec[0]->nirrep() != 1)) {
        auto integral = std::make_shared<IntegralFactory>(primary_);
        PetiteList pet(primary_, integral);
        AO2USO_ = SharedMatrix(pet.aotoso());
        USO2AO_ = AO2USO_->transpose();
    }

    if (AO2USO_) {
        nbf_ = AO2USO_->rowspi()[0];
    } else {
        nbf_ = Dvec[0]->rowspi()[0];
    }

    // Allocate the densities
    if (D_AO_.size() != Dvec.size()) {
        D_AO_.clear();
        for (size_t i = 0; i < Dvec.size(); i++) {
            D_AO_.push_back(std::make_shared<Matrix>("D AO temp", nbf_, nbf_));
        }
    }

    // Copy over the AO
    for (size_t i = 0; i < Dvec.size(); i++) {
        if (Dvec[i]->nirrep() != 1) {
            D_AO_[i]->remove_symmetry(Dvec[i], USO2AO_);
        } else {
            D_AO_[i]->copy(Dvec[i]);
        }
    }
}
void VBase::initialize() {
    timer_on("V: Grid");
    grid_ = std::make_shared<DFTGrid>(primary_->molecule(), primary_, options_);
    timer_off("V: Grid");

    for (size_t i = 0; i < num_threads_; i++) {
        // Need a functional worker per thread
        functional_workers_.push_back(functional_->build_worker());
    }
    
#ifdef USING_BrianQC
    if (brianEnable and brianEnableDFT)
    {
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
#endif
}
SharedMatrix VBase::compute_gradient() { throw PSIEXCEPTION("VBase: gradient not implemented for this V instance."); }
SharedMatrix VBase::compute_hessian() { throw PSIEXCEPTION("VBase: hessian not implemented for this V instance."); }
void VBase::compute_V(std::vector<SharedMatrix> ret) {
    throw PSIEXCEPTION("VBase: deriv not implemented for this V instance.");
}
void VBase::compute_Vx(std::vector<SharedMatrix> Dx, std::vector<SharedMatrix> ret) {
    throw PSIEXCEPTION("VBase: deriv not implemented for this Vx instance.");
}
std::vector<SharedMatrix> VBase::compute_fock_derivatives() {
    throw PSIEXCEPTION("VBase: compute_fock_derivatives not implemented for this Vx instance.");
}
void VBase::set_grac_shift(double grac_shift) {
    // Well this is a flaw in my plan
    if (!grac_initialized_) {
        double grac_alpha = options_.get_double("DFT_GRAC_ALPHA");
        double grac_beta = options_.get_double("DFT_GRAC_BETA");
        auto grac_x_func = std::make_shared<LibXCFunctional>(options_.get_str("DFT_GRAC_X_FUNC"), functional_->is_unpolarized());
        auto grac_c_func = std::make_shared<LibXCFunctional>(options_.get_str("DFT_GRAC_C_FUNC"), functional_->is_unpolarized());

        // Special case for LRC, needs to be this way due to defaults.
        if (functional_->is_x_lrc()) {
            double lr_exch = functional_->x_alpha() + functional_->x_beta();
            grac_x_func->set_alpha(1.0 - lr_exch);
        } else {
            grac_x_func->set_alpha(1.0 - functional_->x_alpha());
        }

        functional_->set_lock(false);
        functional_->set_grac_alpha(grac_alpha);
        functional_->set_grac_beta(grac_beta);
        functional_->set_grac_x_functional(grac_x_func);
        functional_->set_grac_c_functional(grac_c_func);
        functional_->allocate();
        functional_->set_lock(true);
        for (size_t i = 0; i < num_threads_; i++) {
            functional_workers_[i]->set_lock(false);
            functional_workers_[i]->set_grac_alpha(grac_alpha);
            functional_workers_[i]->set_grac_beta(grac_beta);
            functional_workers_[i]->set_grac_x_functional(grac_x_func->build_worker());
            functional_workers_[i]->set_grac_c_functional(grac_c_func->build_worker());
            functional_workers_[i]->allocate();
            functional_workers_[i]->set_lock(true);
        }
        grac_initialized_ = true;
    }

    functional_->set_lock(false);
    functional_->set_grac_shift(grac_shift);
    functional_->set_lock(true);
    for (size_t i = 0; i < num_threads_; i++) {
        functional_workers_[i]->set_lock(false);
        functional_workers_[i]->set_grac_shift(grac_shift);
        functional_workers_[i]->set_lock(true);
    }
}
void VBase::print_header() const {
    outfile->Printf("  ==> DFT Potential <==\n\n");
    functional_->print("outfile", print_);
    grid_->print("outfile", print_);
    if (print_ > 2) grid_->print_details("outfile", print_);
}
std::shared_ptr<BlockOPoints> VBase::get_block(int block) { return grid_->blocks()[block]; }
size_t VBase::nblocks() { return grid_->blocks().size(); }
void VBase::finalize() { grid_.reset(); }
void VBase::build_collocation_cache(size_t memory) {
    // Figure out many blocks to skip

    size_t collocation_size = grid_->collocation_size();
    if (functional_->ansatz() == 1) {
        collocation_size *= 4;  // For gradients
    }
    if (functional_->ansatz() == 2) {
        collocation_size *= 10;  // For gradients and Hessians
    }

    // Figure out stride as closest whole number to amount we need
    size_t stride = (size_t)(std::ceil(collocation_size / (double)memory));

    // More memory than needed
    if (stride == 0) {
        stride = 1;
    }
    cache_map_.clear();

    // Effectively zero blocks saved.
    if (stride > grid_->blocks().size()) {
        return;
    }

    cache_map_deriv_ = point_workers_[0]->deriv();
    auto saved_size_rank = std::vector<size_t>(num_threads_, 0);
    auto ncomputed_rank = std::vector<size_t>(num_threads_, 0);

// Loop over the blocks
#pragma omp parallel for schedule(guided) num_threads(num_threads_)
    for (size_t Q = 0; Q < grid_->blocks().size(); Q += stride) {
        // Get thread info
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif

        // Compute a collocation block
        std::shared_ptr<BlockOPoints> block = grid_->blocks()[Q];
        std::shared_ptr<PointFunctions> pworker = point_workers_[rank];
        pworker->compute_functions(block);

        // Build temps
        size_t nrows = block->npoints();
        size_t ncols = block->local_nbf();
        std::map<std::string, SharedMatrix> collocation_map;

        // Loop over components PHI, PHI_X, PHI_Y, ...
        for (auto& kv : pworker->basis_values()) {
            auto coll = std::make_shared<Matrix>(kv.second->name(), nrows, ncols);

            double** sourcep = kv.second->pointer();
            double** collp = coll->pointer();

            // Matrices are packed in a upper left rectangle, cannot use pure DCOPY
            for (size_t i = 0; i < nrows; i++) {
                C_DCOPY(ncols, sourcep[i], 1, collp[i], 1);
            }
            collocation_map[kv.first] = coll;

            saved_size_rank[rank] += nrows * ncols;
        }
        ncomputed_rank[rank]++;
#pragma omp critical
        cache_map_[block->index()] = collocation_map;
    }

    size_t saved_size = std::accumulate(saved_size_rank.begin(), saved_size_rank.end(), 0.0);
    size_t ncomputed = std::accumulate(ncomputed_rank.begin(), ncomputed_rank.end(), 0.0);

    double gib_saved = 8.0 * (double)saved_size / 1024.0 / 1024.0 / 1024.0;
    double fraction = (double)ncomputed / grid_->blocks().size() * 100;
    if (print_) {
        outfile->Printf("  Cached %.1lf%% of DFT collocation blocks in %.3lf [GiB].\n\n", fraction, gib_saved);
    }
}
void VBase::prepare_vv10_cache(DFTGrid& nlgrid, SharedMatrix D,
                               std::vector<std::map<std::string, SharedVector>>& vv10_cache,
                               std::vector<std::shared_ptr<PointFunctions>>& nl_point_workers, int ansatz) {
    // Densities should be set by the calling functional
    int rank = 0;

    // Build local points workers as they max_points/max_funcs may differ
    const int max_points = nlgrid.max_points();
    const int max_functions = nlgrid.max_functions();

    for (size_t i = 0; i < num_threads_; i++) {
        // Need a points worker per thread, only need RKS-like terms
        auto point_tmp = std::make_shared<RKSFunctions>(primary_, max_points, max_functions);
        point_tmp->set_ansatz(ansatz);
        point_tmp->set_pointers(D);
        nl_point_workers.push_back(point_tmp);
    }

    // => Make the return and "interior" cache <=
    std::vector<std::map<std::string, SharedVector>> vv10_tmp_cache;
    vv10_tmp_cache.resize(nlgrid.blocks().size());

#pragma omp parallel for private(rank) schedule(guided) num_threads(num_threads_)
    for (size_t Q = 0; Q < nlgrid.blocks().size(); Q++) {
// Get thread info
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif

        // Get workers and compute data
        std::shared_ptr<SuperFunctional> fworker = functional_workers_[rank];
        std::shared_ptr<PointFunctions> pworker = nl_point_workers[rank];
        std::shared_ptr<BlockOPoints> block = nlgrid.blocks()[Q];
        // printf("Block %zu\n", Q);

        pworker->compute_points(block);
        vv10_tmp_cache[Q] =
            fworker->compute_vv10_cache(pworker->point_values(), block, vv10_rho_cutoff_, block->npoints(), false);
    }

    // Stitch the cache together to make a single contiguous cache
    size_t total_size = 0;
    for (auto cache : vv10_tmp_cache) {
        total_size += cache["W"]->dimpi()[0];
    }

    // printf("VV10 NL Total size %zu\n", total_size);

    // Leave this as a vector of maps in case we ever revisit the on-the fly manipulation
    vv10_cache.resize(1);
    vv10_cache[0]["W"] = std::make_shared<Vector>("W Grid points", total_size);
    vv10_cache[0]["X"] = std::make_shared<Vector>("X Grid points", total_size);
    vv10_cache[0]["Y"] = std::make_shared<Vector>("Y Grid points", total_size);
    vv10_cache[0]["Z"] = std::make_shared<Vector>("Z Grid points", total_size);
    vv10_cache[0]["RHO"] = std::make_shared<Vector>("RHO Grid points", total_size);
    vv10_cache[0]["W0"] = std::make_shared<Vector>("W0 Grid points", total_size);
    vv10_cache[0]["KAPPA"] = std::make_shared<Vector>("KAPPA Grid points", total_size);

    double* w_vecp = vv10_cache[0]["W"]->pointer();
    double* x_vecp = vv10_cache[0]["X"]->pointer();
    double* y_vecp = vv10_cache[0]["Y"]->pointer();
    double* z_vecp = vv10_cache[0]["Z"]->pointer();
    double* rho_vecp = vv10_cache[0]["RHO"]->pointer();
    double* w0_vecp = vv10_cache[0]["W0"]->pointer();
    double* kappa_vecp = vv10_cache[0]["KAPPA"]->pointer();

    size_t offset = 0;
    for (auto cache : vv10_tmp_cache) {
        size_t csize = cache["W"]->dimpi()[0];
        C_DCOPY(csize, cache["W"]->pointer(), 1, (w_vecp + offset), 1);
        C_DCOPY(csize, cache["X"]->pointer(), 1, (x_vecp + offset), 1);
        C_DCOPY(csize, cache["Y"]->pointer(), 1, (y_vecp + offset), 1);
        C_DCOPY(csize, cache["Z"]->pointer(), 1, (z_vecp + offset), 1);
        C_DCOPY(csize, cache["RHO"]->pointer(), 1, (rho_vecp + offset), 1);
        C_DCOPY(csize, cache["W0"]->pointer(), 1, (w0_vecp + offset), 1);
        C_DCOPY(csize, cache["KAPPA"]->pointer(), 1, (kappa_vecp + offset), 1);

        offset += csize;
    }
}
double VBase::vv10_nlc(SharedMatrix D, SharedMatrix ret) {
    timer_on("V: VV10");
    timer_on("Setup");

    // => VV10 Grid and Cache <=
    std::map<std::string, std::string> opt_map;
    opt_map["DFT_PRUNING_SCHEME"] = "FLAT";

    std::map<std::string, int> opt_int_map;
    opt_int_map["DFT_RADIAL_POINTS"] = options_.get_int("DFT_VV10_RADIAL_POINTS");
    opt_int_map["DFT_SPHERICAL_POINTS"] = options_.get_int("DFT_VV10_SPHERICAL_POINTS");

    DFTGrid nlgrid = DFTGrid(primary_->molecule(), primary_, opt_int_map, opt_map, options_);
    std::vector<std::map<std::string, SharedVector>> vv10_cache;
    std::vector<std::shared_ptr<PointFunctions>> nl_point_workers;
    prepare_vv10_cache(nlgrid, D, vv10_cache, nl_point_workers);

    timer_off("Setup");

    // => Setup info <=
    int rank = 0;
    const int max_functions = nlgrid.max_functions();
    double** Vp = ret->pointer();

    // VV10 temps
    std::vector<double> vv10_exc(num_threads_);

    // Build local points workers as they max_points/max_funcs may differ
    std::vector<SharedMatrix> V_local;
    for (size_t i = 0; i < num_threads_; i++) {
        V_local.push_back(std::make_shared<Matrix>("V Temp", max_functions, max_functions));
    }

// => Compute the kernel <=
// -11.948063
#pragma omp parallel for private(rank) schedule(guided) num_threads(num_threads_)
    for (size_t Q = 0; Q < nlgrid.blocks().size(); Q++) {
// Get thread info
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif

        // Get per rank-workers
        std::shared_ptr<BlockOPoints> block = nlgrid.blocks()[Q];
        std::shared_ptr<SuperFunctional> fworker = functional_workers_[rank];
        std::shared_ptr<PointFunctions> pworker = nl_point_workers[rank];

        // Compute Rho, Phi, etc
        pworker->compute_points(block);

        // Updates the vals map and returns the energy
        std::map<std::string, SharedVector> vals = fworker->values();

        parallel_timer_on("Kernel", rank);
        vv10_exc[rank] += fworker->compute_vv10_kernel(pworker->point_values(), vv10_cache, block);
        parallel_timer_off("Kernel", rank);

        parallel_timer_on("VV10 Fock", rank);

        // => LSDA and GGA contribution (symmetrized) <= //
        dft_integrators::rks_integrator(block, fworker, pworker, V_local[rank], 1);

        // => Unpacking <= //
        const std::vector<int>& function_map = block->functions_local_to_global();
        int nlocal = function_map.size();
        double** V2p = V_local[rank]->pointer();

        for (int ml = 0; ml < nlocal; ml++) {
            int mg = function_map[ml];
            for (int nl = 0; nl < ml; nl++) {
                int ng = function_map[nl];
#pragma omp atomic update
                Vp[mg][ng] += V2p[ml][nl];
#pragma omp atomic update
                Vp[ng][mg] += V2p[ml][nl];
            }
#pragma omp atomic update
            Vp[mg][mg] += V2p[ml][ml];
        }
        parallel_timer_off("VV10 Fock", rank);
    }

    double vv10_e = std::accumulate(vv10_exc.begin(), vv10_exc.end(), 0.0);
    timer_off("V: VV10");
    return vv10_e;
}
SharedMatrix VBase::vv10_nlc_gradient(SharedMatrix D) {
    /* Not yet finished, missing several components*/
    throw PSIEXCEPTION("V: Cannot compute VV10 gradient contribution.");

    timer_on("V: VV10");
    timer_on("Setup");

    // => VV10 Grid and Cache <=
    std::map<std::string, std::string> opt_map;
    opt_map["DFT_PRUNING_SCHEME"] = "FLAT";
    // opt_map["DFT_NUCLEAR_SCHEME"] = "BECKE";

    std::map<std::string, int> opt_int_map;
    opt_int_map["DFT_RADIAL_POINTS"] = options_.get_int("DFT_VV10_RADIAL_POINTS");
    opt_int_map["DFT_SPHERICAL_POINTS"] = options_.get_int("DFT_VV10_SPHERICAL_POINTS");

    DFTGrid nlgrid = DFTGrid(primary_->molecule(), primary_, opt_int_map, opt_map, options_);
    std::vector<std::map<std::string, SharedVector>> vv10_cache;
    std::vector<std::shared_ptr<PointFunctions>> nl_point_workers;
    prepare_vv10_cache(nlgrid, D, vv10_cache, nl_point_workers, 2);

    timer_off("Setup");

    // => Setup info <=
    int rank = 0;
    const int max_functions = nlgrid.max_functions();
    const int max_points = nlgrid.max_points();
    const int natom = primary_->molecule()->natom();

    // VV10 temps
    std::vector<double> vv10_exc(num_threads_);

    // Per thread temporaries
    std::vector<SharedMatrix> G_local, U_local;
    for (size_t i = 0; i < num_threads_; i++) {
        G_local.push_back(std::make_shared<Matrix>("G Temp", natom, 3));
        U_local.push_back(std::make_shared<Matrix>("U Temp", max_points, max_functions));
    }

// => Compute the kernel <=
#pragma omp parallel for private(rank) schedule(guided) num_threads(num_threads_)
    for (size_t Q = 0; Q < nlgrid.blocks().size(); Q++) {
// Get thread info
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif

        // Get per rank-workers
        std::shared_ptr<BlockOPoints> block = nlgrid.blocks()[Q];
        std::shared_ptr<SuperFunctional> fworker = functional_workers_[rank];
        std::shared_ptr<PointFunctions> pworker = nl_point_workers[rank];
        const std::vector<int>& function_map = block->functions_local_to_global();
        const int nlocal = function_map.size();
        const int npoints = block->npoints();
        double** Tp = pworker->scratch()[0]->pointer();

        // Compute Rho, Phi, etc
        pworker->compute_points(block);

        // Updates the vals map and returns the energy
        std::map<std::string, SharedVector> vals = fworker->values();

        parallel_timer_on("Kernel", rank);
        vv10_exc[rank] += fworker->compute_vv10_kernel(pworker->point_values(), vv10_cache, block, npoints, true);
        parallel_timer_off("Kernel", rank);

        parallel_timer_on("V_xc gradient", rank);

        // => LSDA and GGA gradient contributions <= //
        dft_integrators::rks_gradient_integrator(primary_, block, fworker, pworker, G_local[rank], U_local[rank]);

        // => Grid gradient contributions <= //
        double** Gp = G_local[rank]->pointer();
        const double* x_grid = fworker->vv_value("GRID_WX")->pointer();
        const double* y_grid = fworker->vv_value("GRID_WY")->pointer();
        const double* z_grid = fworker->vv_value("GRID_WZ")->pointer();
        double** phi = pworker->basis_value("PHI")->pointer();
        double** phi_x = pworker->basis_value("PHI_X")->pointer();
        double** phi_y = pworker->basis_value("PHI_Y")->pointer();
        double** phi_z = pworker->basis_value("PHI_Z")->pointer();

        // These terms are incorrect until they are able to isolate blocks on a single atom due to
        // the requirement of the sum to not include blocks on the same atom
        for (int P = 0; P < npoints; P++) {
            std::fill(Tp[P], Tp[P] + nlocal, 0.0);
            C_DAXPY(nlocal, z_grid[P], phi[P], 1, Tp[P], 1);
        }
        for (int ml = 0; ml < nlocal; ml++) {
            int A = primary_->function_to_center(function_map[ml]);
            // Gp[A][0] += C_DDOT(npoints, &Tp[0][ml], max_functions, &phi_x[0][ml], max_functions);
            // Gp[A][1] += C_DDOT(npoints, &Tp[0][ml], max_functions, &phi_y[0][ml], max_functions);
            Gp[A][2] += C_DDOT(npoints, &Tp[0][ml], max_functions, &phi_z[0][ml], max_functions);
            // printf("Value %d %16.15lf\n", A, C_DDOT(npoints, &Tp[0][ml], max_functions, &phi_z[0][ml],
            // max_functions));
        }

        // printf("--\n");

        parallel_timer_off("V_xc gradient", rank);
    }

    // Sum up the matrix
    auto G = std::make_shared<Matrix>("XC Gradient", natom, 3);
    for (auto const& val : G_local) {
        G->add(val);
    }
    G->print();
    G->zero();

    double vv10_e = std::accumulate(vv10_exc.begin(), vv10_exc.end(), 0.0);
    timer_off("V: VV10");
    return G;
}

SAP::SAP(std::shared_ptr<SuperFunctional> functional, std::shared_ptr<BasisSet> primary, Options& options)
    : VBase(functional, primary, options) {}
SAP::~SAP() {}
void SAP::initialize() {
    VBase::initialize();
    int max_points = grid_->max_points();
    int max_functions = grid_->max_functions();
    for (size_t i = 0; i < num_threads_; i++) {
        // Need a points worker per thread
        auto point_tmp = std::make_shared<SAPFunctions>(primary_, max_points, max_functions);
        // This is like LDA
        point_tmp->set_ansatz(0);
        point_tmp->set_cache_map(&cache_map_);
        point_workers_.push_back(point_tmp);
    }

    // Initialize symmetry
    auto integral = std::make_shared<IntegralFactory>(primary_);
    PetiteList pet(primary_, integral);
    AO2USO_ = SharedMatrix(pet.aotoso());
    USO2AO_ = AO2USO_->transpose();
    nbf_ = AO2USO_->rowspi()[0];
}
void SAP::finalize() { VBase::finalize(); }
void SAP::print_header() const {
    outfile->Printf("  ==> SAP guess <==\n\n");
    grid_->print("outfile", print_);
    if (print_ > 2) grid_->print_details("outfile", print_);
}
void SAP::compute_V(std::vector<SharedMatrix> ret) {
    timer_on("SAP: Form V");

    if (ret.size() != 1) {
        throw PSIEXCEPTION("SAP outputs only one V Matrix");
    }

    // Thread info
    int rank = 0;

    // How many functions are there (for lda in Vtemp, T)
    int max_functions = grid_->max_functions();
    int max_points = grid_->max_points();

    // Per thread temporaries
    std::vector<SharedMatrix> V_local;
    for (size_t i = 0; i < num_threads_; i++) {
        V_local.push_back(std::make_shared<Matrix>("V Temp", max_functions, max_functions));
    }

    auto V_AO = std::make_shared<Matrix>("V AO Temp", nbf_, nbf_);
    double** Vp = V_AO->pointer();

    // Nuclear coordinates
    std::vector<double> nucx, nucy, nucz, nucZ;
    nucx.resize(primary_->molecule()->natom());
    nucy.resize(primary_->molecule()->natom());
    nucz.resize(primary_->molecule()->natom());
    nucZ.resize(primary_->molecule()->natom());
    for (size_t iatom = 0; iatom < nucx.size(); iatom++) {
        nucx[iatom] = primary_->molecule()->x(iatom);
        nucy[iatom] = primary_->molecule()->y(iatom);
        nucz[iatom] = primary_->molecule()->z(iatom);
        nucZ[iatom] = primary_->molecule()->Z(iatom);
    }

// Traverse the blocks of points
#pragma omp parallel for private(rank) schedule(guided) num_threads(num_threads_)
    for (size_t Q = 0; Q < grid_->blocks().size(); Q++) {
// Get thread info
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif

        // Get per-rank workers
        std::shared_ptr<BlockOPoints> block = grid_->blocks()[Q];
        std::shared_ptr<SuperFunctional> fworker = functional_workers_[rank];
        std::shared_ptr<PointFunctions> pworker = point_workers_[rank];

        // Compute Rho, Phi, etc
        parallel_timer_on("Properties", rank);
        pworker->compute_points(block, false);
        parallel_timer_off("Properties", rank);

        // Compute the SAP potential
        parallel_timer_on("Functional", rank);
        SharedVector sap_potential = std::make_shared<Vector>("sappot", block->npoints());
        for (int ip = 0; ip < block->npoints(); ip++) {
            // Coordinates of the point
            double xi = block->x()[ip];
            double yi = block->y()[ip];
            double zi = block->z()[ip];
            // and its potential
            double V = 0.0;

            // Loop over nuclei
            for (size_t iatom = 0; iatom < nucx.size(); iatom++) {
                // Distance to nucleus is
                double dx = xi - nucx[iatom];
                double dy = yi - nucy[iatom];
                double dz = zi - nucz[iatom];
                double r = sqrt(dx * dx + dy * dy + dz * dz);
                // and the SAP potential at this point is
                V -= ::sap_effective_charge(nucZ[iatom], r) / r;
            }

            // Store
            (*sap_potential)[ip] = V;
        }

        parallel_timer_off("Functional", rank);

        if (debug_ > 4) {
            block->print("outfile", debug_);
            pworker->print("outfile", debug_);
        }

        parallel_timer_on("V_xc", rank);

        // => LSDA contribution (symmetrized) <= //
        dft_integrators::sap_integrator(block, sap_potential, pworker, V_local[rank]);

        // => Unpacking <= //
        double** V2p = V_local[rank]->pointer();
        const std::vector<int>& function_map = block->functions_local_to_global();
        int nlocal = function_map.size();

        for (int ml = 0; ml < nlocal; ml++) {
            int mg = function_map[ml];
            for (int nl = 0; nl < ml; nl++) {
                int ng = function_map[nl];
#pragma omp atomic update
                Vp[mg][ng] += V2p[ml][nl];
#pragma omp atomic update
                Vp[ng][mg] += V2p[ml][nl];
            }
#pragma omp atomic update
            Vp[mg][mg] += V2p[ml][ml];
        }
        parallel_timer_off("V_xc", rank);
    }

    // Set the result
    if (AO2USO_) {
        ret[0]->apply_symmetry(V_AO, AO2USO_);
    } else {
        ret[0]->copy(V_AO);
    }
    timer_off("SAP: Form V");
}

RV::RV(std::shared_ptr<SuperFunctional> functional, std::shared_ptr<BasisSet> primary, Options& options)
    : VBase(functional, primary, options) {}
RV::~RV() {}
void RV::initialize() {
    VBase::initialize();
    int max_points = grid_->max_points();
    int max_functions = grid_->max_functions();
    for (size_t i = 0; i < num_threads_; i++) {
        // Need a points worker per thread
        auto point_tmp = std::make_shared<RKSFunctions>(primary_, max_points, max_functions);
        point_tmp->set_ansatz(functional_->ansatz());
        point_tmp->set_cache_map(&cache_map_);
        point_workers_.push_back(point_tmp);
    }
}
void RV::finalize() { VBase::finalize(); }
void RV::print_header() const { VBase::print_header(); }
void RV::compute_V(std::vector<SharedMatrix> ret) {
    timer_on("RV: Form V");
    
    if ((D_AO_.size() != 1) || (ret.size() != 1)) {
        throw PSIEXCEPTION("V: RKS should have only one D/V Matrix");
    }
    
#ifdef USING_BrianQC
    if (brianEnable and brianEnableDFT) {
        double DFTEnergy;
        
        brianSCFBuildFockDFT(&brianCookie,
            D_AO_[0]->get_pointer(0),
            nullptr,
            ret[0]->get_pointer(0),
            nullptr,
            &DFTEnergy
        );
        checkBrian();
        
        quad_values_["VV10"] = 0.0; // NOTE: BrianQC doesn't compute the VV10 term separately, it just includes it in the DFT energy term
        quad_values_["FUNCTIONAL"] = DFTEnergy;
        quad_values_["RHO_A"] = 0.0;
        quad_values_["RHO_AX"] = 0.0;
        quad_values_["RHO_AY"] = 0.0;
        quad_values_["RHO_AZ"] = 0.0;
        quad_values_["RHO_B"] = 0.0;
        quad_values_["RHO_BX"] = 0.0;
        quad_values_["RHO_BY"] = 0.0;
        quad_values_["RHO_BZ"] = 0.0;
        
        return;
    }
#endif

    // Thread info
    int rank = 0;

    // What local XC ansatz are we in?
    int ansatz = functional_->ansatz();

    // How many functions are there (for lda in Vtemp, T)
    int max_functions = grid_->max_functions();
    int max_points = grid_->max_points();

    // Setup the pointers
    for (size_t i = 0; i < num_threads_; i++) {
        point_workers_[i]->set_pointers(D_AO_[0]);
    }

    // Per thread temporaries
    std::vector<SharedMatrix> V_local;
    for (size_t i = 0; i < num_threads_; i++) {
        V_local.push_back(std::make_shared<Matrix>("V Temp", max_functions, max_functions));
    }

    auto V_AO = std::make_shared<Matrix>("V AO Temp", nbf_, nbf_);
    auto Vp = V_AO->pointer();

    std::vector<double> functionalq(num_threads_);
    std::vector<double> rhoaq(num_threads_);
    std::vector<double> rhoaxq(num_threads_);
    std::vector<double> rhoayq(num_threads_);
    std::vector<double> rhoazq(num_threads_);

// VV10 kernel data if requested

// Traverse the blocks of points
#pragma omp parallel for private(rank) schedule(guided) num_threads(num_threads_)
    for (size_t Q = 0; Q < grid_->blocks().size(); Q++) {
// Get thread info
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif

        // Get per-rank workers
        auto block = grid_->blocks()[Q];
        auto fworker = functional_workers_[rank];
        auto pworker = point_workers_[rank];

        // Compute Rho, Phi, etc
        parallel_timer_on("Properties", rank);
        pworker->compute_points(block, false);
        parallel_timer_off("Properties", rank);

        // Compute functional values
        parallel_timer_on("Functional", rank);
        fworker->compute_functional(pworker->point_values());
        parallel_timer_off("Functional", rank);

        if (debug_ > 4) {
            block->print("outfile", debug_);
            pworker->print("outfile", debug_);
        }

        parallel_timer_on("V_xc", rank);

        // => Compute quadrature <= //
        auto qvals = dft_integrators::rks_quadrature_integrate(block, fworker, pworker);
        functionalq[rank] += qvals[0];
        rhoaq[rank] += qvals[1];
        rhoaxq[rank] += qvals[2];
        rhoayq[rank] += qvals[3];
        rhoazq[rank] += qvals[4];

        // => LSDA, GGA, and meta contribution (symmetrized) <= //
        dft_integrators::rks_integrator(block, fworker, pworker, V_local[rank]);

        // => Unpacking <= //
        auto V2p = V_local[rank]->pointer();
        const auto& function_map = block->functions_local_to_global();
        int nlocal = function_map.size();

        for (int ml = 0; ml < nlocal; ml++) {
            int mg = function_map[ml];
            for (int nl = 0; nl < ml; nl++) {
                int ng = function_map[nl];
#pragma omp atomic update
                Vp[mg][ng] += V2p[ml][nl];
#pragma omp atomic update
                Vp[ng][mg] += V2p[ml][nl];
            }
#pragma omp atomic update
            Vp[mg][mg] += V2p[ml][ml];
        }
        parallel_timer_off("V_xc", rank);
    }

    // Do we need VV10?
    double vv10_e = 0.0;
    if (functional_->needs_vv10()) {
        vv10_e = vv10_nlc(D_AO_[0], V_AO);
    }

    // Set the result
    if (AO2USO_) {
        ret[0]->apply_symmetry(V_AO, AO2USO_);
    } else {
        ret[0]->copy(V_AO);
    }

    quad_values_["VV10"] = vv10_e;
    quad_values_["FUNCTIONAL"] = std::accumulate(functionalq.begin(), functionalq.end(), 0.0);
    quad_values_["RHO_A"] = std::accumulate(rhoaq.begin(), rhoaq.end(), 0.0);
    quad_values_["RHO_AX"] = std::accumulate(rhoaxq.begin(), rhoaxq.end(), 0.0);
    quad_values_["RHO_AY"] = std::accumulate(rhoayq.begin(), rhoayq.end(), 0.0);
    quad_values_["RHO_AZ"] = std::accumulate(rhoazq.begin(), rhoazq.end(), 0.0);
    quad_values_["RHO_B"] = quad_values_["RHO_A"];
    quad_values_["RHO_BX"] = quad_values_["RHO_AX"];
    quad_values_["RHO_BY"] = quad_values_["RHO_AY"];
    quad_values_["RHO_BZ"] = quad_values_["RHO_AZ"];

    if (std::isnan(quad_values_["FUNCTIONAL"])) {
        throw PSIEXCEPTION("V: Integrated DFT functional to get NaN. The functional is not numerically stable. Pick a different one.");
    }

    if (debug_) {
        outfile->Printf("   => Numerical Integrals <=\n\n");
        outfile->Printf("    VV10 Value:         %24.16E\n", quad_values_["VV10"]);
        outfile->Printf("    Functional Value:   %24.16E\n", quad_values_["FUNCTIONAL"]);
        outfile->Printf("    <\\rho_a>        :  %24.16E\n", quad_values_["RHO_A"]);
        outfile->Printf("    <\\rho_b>        :  %24.16E\n", quad_values_["RHO_B"]);
        outfile->Printf("    <\\vec r\\rho_a> : <%24.16E,%24.16E,%24.16E>\n", quad_values_["RHO_AX"],
                        quad_values_["RHO_AY"], quad_values_["RHO_AZ"]);
        outfile->Printf("    <\\vec r\\rho_b> : <%24.16E,%24.16E,%24.16E>\n\n", quad_values_["RHO_BX"],
                        quad_values_["RHO_BY"], quad_values_["RHO_BZ"]);
    }
    timer_off("RV: Form V");
}

std::vector<SharedMatrix> RV::compute_fock_derivatives() {
    timer_on("RV: Form Fx");

    int natoms = primary_->molecule()->natom();
    std::vector<SharedMatrix> Vx(3*natoms);
    for(int n = 0; n < 3*natoms; ++n)
        Vx[n] = std::make_shared<Matrix>("Vx for Perturbation " + std::to_string(n), nbf_, nbf_);
    if (D_AO_.size() != 1) {
        throw PSIEXCEPTION("DFT Hessian: RKS should have only one D Matrix");
    }

    if (functional_->needs_vv10()) {
        throw PSIEXCEPTION("DFT Hessian: RKS cannot compute VV10 Fx contribution.");
    }

    // Thread info
    int rank = 0;

    // What local XC ansatz are we in?
    int ansatz = functional_->ansatz();
    if (ansatz >= 1) {
        throw PSIEXCEPTION("DFT Hessian: RKS does not support GGAs or MGGAs yet");
    }

    int old_func_deriv = functional_->deriv();

    // How many functions are there (for lda in Vtemp, T)
    int max_functions = grid_->max_functions();
    int max_points = grid_->max_points();

    // Set pointers to SCF density
    for (size_t i = 0; i < num_threads_; i++) {
        point_workers_[i]->set_pointers(D_AO_[0]);
        point_workers_[i]->set_deriv(1);
    }

    // Per [R]ank quantities
    std::vector<std::shared_ptr<Vector>> R_rho_x, R_rho_y, R_rho_z;
    std::vector<SharedMatrix> R_Vx_local;
    for (size_t i = 0; i < num_threads_; i++) {
        R_Vx_local.push_back(std::make_shared<Matrix>("Vx Temp", max_functions, max_functions));

        functional_workers_[i]->set_deriv(2);
        functional_workers_[i]->allocate();
    }
    // Output quantities
    std::vector<SharedMatrix> Vx_AO;
    for (size_t i = 0; i < 3*natoms; i++) {
        Vx_AO.push_back(std::make_shared<Matrix>("Vx AO Temp", nbf_, nbf_));
    }

// Traverse the blocks of points
#pragma omp parallel for private(rank) schedule(guided) num_threads(num_threads_)
    for (size_t Q = 0; Q < grid_->blocks().size(); Q++) {
// Get thread info
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif

        // => Setup <= //
        std::shared_ptr<SuperFunctional> fworker = functional_workers_[rank];
        std::shared_ptr<PointFunctions> pworker = point_workers_[rank];
        double **Vx_localp = R_Vx_local[rank]->pointer();

        // => Compute blocks <= //
        double** Tp = pworker->scratch()[0]->pointer();

        std::shared_ptr<BlockOPoints> block = grid_->blocks()[Q];
        int npoints = block->npoints();
        double* w = block->w();
        const std::vector<int>& function_map = block->functions_local_to_global();
        double** Dp = pworker->D_scratch()[0]->pointer();
        int nlocal = function_map.size();

        // Compute Rho, Phi, etc
        parallel_timer_on("Properties", rank);
        pworker->compute_points(block);
        parallel_timer_off("Properties", rank);

        // Compute functional values

        parallel_timer_on("Functional", rank);
        std::map<std::string, SharedVector>& vals = fworker->compute_functional(pworker->point_values(), npoints);
        parallel_timer_off("Functional", rank);

        // => Grab quantities <= //
        // LDA
        double** phi = pworker->basis_value("PHI")->pointer();
        double** phi_x = pworker->basis_value("PHI_X")->pointer();
        double** phi_y = pworker->basis_value("PHI_Y")->pointer();
        double** phi_z = pworker->basis_value("PHI_Z")->pointer();
        double* rho_a = pworker->point_value("RHO_A")->pointer();
        double* v_rho_a = vals["V_RHO_A"]->pointer();
        double* v_rho_aa = vals["V_RHO_A_RHO_A"]->pointer();
        for (int P = 0; P < npoints; P++) {
            if (std::fabs(rho_a[P]) < v2_rho_cutoff_) {
                v_rho_a[P] = 0.0;
                v_rho_aa[P] = 0.0;
            }
        }
        size_t coll_funcs = pworker->basis_value("PHI")->ncol();
        for(int atom = 0; atom < primary_->molecule()->natom(); ++atom){
            // Find first and last basis functions on this atom, from the subset of bfs being handled by this block of points
            auto first_func_iter = std::find_if(function_map.begin(), function_map.end(), [&](int i) {return primary_->function_to_center(i) == atom;});
            if(first_func_iter == function_map.end()) continue;
            auto last_func_riter = std::find_if(function_map.rbegin(), function_map.rend(), [&](int i) {return primary_->function_to_center(i) == atom;});
            if(last_func_riter == function_map.rend()) continue;
            auto last_func_iter = last_func_riter.base(); // convert to forward iterator

            int first_func_addr = std::distance(function_map.begin(), first_func_iter);
            int nfuncs = std::distance(first_func_iter, last_func_iter);

            /*
             * X derivatives
             */
            // T = ɸ D, remembering that only the bfs centered on the current atom of interest contribute to the derivative
            C_DGEMM('N', 'N', npoints, nfuncs, nlocal, 1.0, phi[0], coll_funcs, &Dp[0][first_func_addr], max_functions, 0.0, Tp[0], max_functions);
            for (int P = 0; P < npoints; P++) {
                // ρ_x  = T ɸ_x^t
                double rho_xP = C_DDOT(nfuncs, Tp[P], 1, &phi_x[P][first_func_addr], 1);
                // Now redefine the intermediate T:
                //
                //       /  | ∂^2 F
                // T <- | ɸ | ----- ρ_x
                //       \  | ∂ ρ^2
                std::fill(Tp[P], Tp[P] + nlocal, 0);
                C_DAXPY(nlocal, -v_rho_aa[P] * w[P] * rho_xP, phi[P], 1, Tp[P], 1);
                //
                //       /    | ∂ F
                // T <- | ɸ_x | ---
                //       \    | ∂ ρ
                C_DAXPY(nfuncs, -0.5 * v_rho_a[P] * w[P], &phi_x[P][first_func_addr], 1, &Tp[P][first_func_addr], 1);
            }
            //         |  \
            // Vx <- T | ɸ |
            //         |  /
            C_DGEMM('T', 'N', nlocal, nlocal, npoints, 1.0, Tp[0], max_functions, phi[0], coll_funcs, 0.0, Vx_localp[0], max_functions);
            // => Accumulate the result <= //
            double **Vxp = Vx[3*atom + 0]->pointer();
            for (int ml = 0; ml < nlocal; ml++) {
                int mg = function_map[ml];
                for (int nl = 0; nl < nlocal; nl++) {
                    int ng = function_map[nl];
                    double result = Vx_localp[ml][nl] + Vx_localp[nl][ml];
#pragma omp atomic update
                     Vxp[mg][ng] += result;
#pragma omp atomic update
                     Vxp[ng][mg] += result;
                }
            }

            /*
             * Y derivatives
             */
            // T = ɸ D, remembering that only the bfs centered on the current atom of interest contribute to the derivative
            C_DGEMM('N', 'N', npoints, nfuncs, nlocal, 1.0, phi[0], coll_funcs, &Dp[0][first_func_addr], max_functions, 0.0, Tp[0], max_functions);
            for (int P = 0; P < npoints; P++) {
                // ρ_y  = T ɸ_y^t
                double rho_yP = C_DDOT(nfuncs, Tp[P], 1, &phi_y[P][first_func_addr], 1);
                // Now redefine the intermediate T:
                //
                //       /  | ∂^2 F
                // T <- | ɸ | ----- ρ_y
                //       \  | ∂ ρ^2
                std::fill(Tp[P], Tp[P] + nlocal, 0);
                C_DAXPY(nlocal, -v_rho_aa[P] * w[P] * rho_yP, phi[P], 1, Tp[P], 1);
                //
                //       /    | ∂ F
                // T <- | ɸ_y | ---
                //       \    | ∂ ρ
                C_DAXPY(nfuncs, -0.5*v_rho_a[P] * w[P], &phi_y[P][first_func_addr], 1, &Tp[P][first_func_addr], 1);
            }
            //         |  \
            // Vx <- T | ɸ |
            //         |  /
            C_DGEMM('T', 'N', nlocal, nlocal, npoints, 1.0, Tp[0], max_functions, phi[0], coll_funcs, 0.0, Vx_localp[0], max_functions);
            // => Accumulate the result <= //
            double **Vyp = Vx[3*atom + 1]->pointer();
            for (int ml = 0; ml < nlocal; ml++) {
                int mg = function_map[ml];
                for (int nl = 0; nl < nlocal; nl++) {
                    int ng = function_map[nl];
                    double result = Vx_localp[ml][nl] + Vx_localp[nl][ml];
#pragma omp atomic update
                     Vyp[mg][ng] += result;
#pragma omp atomic update
                     Vyp[ng][mg] += result;
                }
            }

            /*
             * Z derivatives
             */
            // T = ɸ D, remembering that only the bfs centered on the current atom of interest contribute to the derivative
            C_DGEMM('N', 'N', npoints, nfuncs, nlocal, 1.0, phi[0], coll_funcs, &Dp[0][first_func_addr], max_functions, 0.0, Tp[0], max_functions);
            for (int P = 0; P < npoints; P++) {
                // ρ_z  = T ɸ_z^t
                double rho_zP = C_DDOT(nfuncs, Tp[P], 1, &phi_z[P][first_func_addr], 1);
                // Now redefine the intermediate T:
                //
                //       /  | ∂^2 F
                // T <- | ɸ | ----- ρ_z
                //       \  | ∂ ρ^2
                std::fill(Tp[P], Tp[P] + nlocal, 0);
                C_DAXPY(nlocal, -v_rho_aa[P] * w[P] * rho_zP, phi[P], 1, Tp[P], 1);
                //
                //       /    | ∂ F
                // T <- | ɸ_z | ---
                //       \    | ∂ ρ
                C_DAXPY(nfuncs, -0.5*v_rho_a[P] * w[P], &phi_z[P][first_func_addr], 1, &Tp[P][first_func_addr], 1);
            }
            //         |  \
            // Vx <- T | ɸ |
            //         |  /
            C_DGEMM('T', 'N', nlocal, nlocal, npoints, 1.0, Tp[0], max_functions, phi[0], coll_funcs, 0.0, Vx_localp[0], max_functions);
            // => Accumulate the result <= //
            double **Vzp = Vx[3*atom + 2]->pointer();
            for (int ml = 0; ml < nlocal; ml++) {
                int mg = function_map[ml];
                for (int nl = 0; nl < nlocal; nl++) {
                    int ng = function_map[nl];
                    double result = Vx_localp[ml][nl] + Vx_localp[nl][ml];
#pragma omp atomic update
                     Vzp[mg][ng] += result;
#pragma omp atomic update
                     Vzp[ng][mg] += result;
                }
            }
        }
    }

    // Reset the workers
    for (size_t i = 0; i < num_threads_; i++) {
        functional_workers_[i]->set_deriv(old_func_deriv);
        functional_workers_[i]->allocate();
    }
    timer_off("RV: Form Fx");
    return Vx;
}

void RV::compute_Vx(std::vector<SharedMatrix> Dx, std::vector<SharedMatrix> ret) {
    timer_on("RV: Form Vx");

    if (D_AO_.size() != 1) {
        throw PSIEXCEPTION("Vx: RKS should have only one D Matrix");
    }
    if ((Dx.size() != ret.size()) || (Dx.size() == 0)) {
        throw PSIEXCEPTION("Vx: RKS input and output sizes should be the same.");
    }

    if (functional_->needs_vv10()) {
        throw PSIEXCEPTION("Vx: RKS cannot compute VV10 Vx contribution.");
    }

    // Thread info
    int rank = 0;

    // What local XC ansatz are we in?
    int ansatz = functional_->ansatz();
    if (ansatz >= 2) {
        throw PSIEXCEPTION("Vx: RKS does not support rotated V builds for MGGA's");
    }

    int old_point_deriv = point_workers_[0]->deriv();
    int old_func_deriv = functional_->deriv();

    // How many functions are there (for lda in Vtemp, T)
    int max_functions = grid_->max_functions();
    int max_points = grid_->max_points();

    // Set pointers to SCF density
    for (size_t i = 0; i < num_threads_; i++) {
        point_workers_[i]->set_pointers(D_AO_[0]);
    }

    std::vector<SharedMatrix> Dx_vec;
    for (size_t i = 0; i < Dx.size(); i++) {
        if (Dx[i]->nirrep() != 1) {
            auto Dx_mat = std::make_shared<Matrix>("D AO temp", nbf_, nbf_);
            Dx_mat->remove_symmetry(Dx[i], USO2AO_);
            Dx_vec.push_back(Dx_mat);
        } else {
            Dx_vec.push_back(Dx[i]);
        }
    }

    // Per [R]ank quantities
    std::vector<SharedMatrix> R_Vx_local, R_Dx_local;
    std::vector<std::shared_ptr<Vector>> R_rho_k, R_rho_k_x, R_rho_k_y, R_rho_k_z, R_gamma_k;
    for (size_t i = 0; i < num_threads_; i++) {
        R_Vx_local.push_back(std::make_shared<Matrix>("Vx Temp", max_functions, max_functions));
        R_Dx_local.push_back(std::make_shared<Matrix>("Dk Temp", max_functions, max_functions));

        R_rho_k.push_back(std::make_shared<Vector>("Rho K Temp", max_points));

        if (ansatz >= 1) {
            R_rho_k_x.push_back(std::make_shared<Vector>("RHO K X Temp", max_points));
            R_rho_k_y.push_back(std::make_shared<Vector>("RHO K Y Temp", max_points));
            R_rho_k_z.push_back(std::make_shared<Vector>("Rho K Z Temp", max_points));
            R_gamma_k.push_back(std::make_shared<Vector>("Gamma K Temp", max_points));
        }

        functional_workers_[i]->set_deriv(2);
        functional_workers_[i]->allocate();
    }

    // Output quantities
    std::vector<SharedMatrix> Vx_AO;
    for (size_t i = 0; i < Dx.size(); i++) {
        Vx_AO.push_back(std::make_shared<Matrix>("Vx AO Temp", nbf_, nbf_));
    }

// Traverse the blocks of points
#pragma omp parallel for private(rank) schedule(guided) num_threads(num_threads_)
    for (size_t Q = 0; Q < grid_->blocks().size(); Q++) {
// Get thread info
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif

        // => Setup <= //
        std::shared_ptr<SuperFunctional> fworker = functional_workers_[rank];
        std::shared_ptr<PointFunctions> pworker = point_workers_[rank];
        double** Vx_localp = R_Vx_local[rank]->pointer();
        double** Dx_localp = R_Dx_local[rank]->pointer();

        // => Compute blocks <= //
        double** Tp = pworker->scratch()[0]->pointer();

        std::shared_ptr<BlockOPoints> block = grid_->blocks()[Q];
        int npoints = block->npoints();
        double* w = block->w();
        const std::vector<int>& function_map = block->functions_local_to_global();
        int nlocal = function_map.size();

        // Compute Rho, Phi, etc
        parallel_timer_on("Properties", rank);
        pworker->compute_points(block);
        parallel_timer_off("Properties", rank);

        // Compute functional values

        parallel_timer_on("Functional", rank);
        std::map<std::string, SharedVector>& vals = fworker->compute_functional(pworker->point_values(), npoints);
        parallel_timer_off("Functional", rank);

        // => Grab quantities <= //
        // LDA
        double** phi = pworker->basis_value("PHI")->pointer();
        double* rho_a = pworker->point_value("RHO_A")->pointer();
        double* v2_rho2 = vals["V_RHO_A_RHO_A"]->pointer();
        double* rho_k = R_rho_k[rank]->pointer();
        size_t coll_funcs = pworker->basis_value("PHI")->ncol();

        // GGA
        double* rho_k_x;
        double* rho_k_y;
        double* rho_k_z;
        double* gamma_k;
        double** phi_x;
        double** phi_y;
        double** phi_z;
        double* rho_x;
        double* rho_y;
        double* rho_z;
        if (ansatz >= 1) {
            rho_k_x = R_rho_k_x[rank]->pointer();
            rho_k_y = R_rho_k_y[rank]->pointer();
            rho_k_z = R_rho_k_z[rank]->pointer();
            gamma_k = R_gamma_k[rank]->pointer();
            phi_x = pworker->basis_value("PHI_X")->pointer();
            phi_y = pworker->basis_value("PHI_Y")->pointer();
            phi_z = pworker->basis_value("PHI_Z")->pointer();
            rho_x = pworker->point_value("RHO_AX")->pointer();
            rho_y = pworker->point_value("RHO_AY")->pointer();
            rho_z = pworker->point_value("RHO_AZ")->pointer();
        }

        // Meta
        // Forget that!

        // Loop over perturbation tensors
        for (size_t dindex = 0; dindex < Dx_vec.size(); dindex++) {
            double** Dxp = Dx_vec[dindex]->pointer();

            // => Build Rotated Densities <= //
            for (int ml = 0; ml < nlocal; ml++) {
                int mg = function_map[ml];
                for (int nl = 0; nl < nlocal; nl++) {
                    int ng = function_map[nl];
                    Dx_localp[ml][nl] = Dxp[mg][ng];
                }
            }

            parallel_timer_on("Derivative Properties", rank);
            // Rho_a = D^k_xy phi_xa phi_ya
            C_DGEMM('N', 'N', npoints, nlocal, nlocal, 1.0, phi[0], coll_funcs, Dx_localp[0], max_functions, 0.0, Tp[0],
                    max_functions);
            C_DGEMM('N', 'T', npoints, nlocal, nlocal, 1.0, phi[0], coll_funcs, Dx_localp[0], max_functions, 1.0, Tp[0],
                    max_functions);

            for (int P = 0; P < npoints; P++) {
                rho_k[P] = 0.5 * C_DDOT(nlocal, phi[P], 1, Tp[P], 1);
            }

            // Rho^d_k and gamma_k
            if (ansatz >= 1) {
                for (int P = 0; P < npoints; P++) {
                    rho_k_x[P] = C_DDOT(nlocal, phi_x[P], 1, Tp[P], 1);
                    rho_k_y[P] = C_DDOT(nlocal, phi_y[P], 1, Tp[P], 1);
                    rho_k_z[P] = C_DDOT(nlocal, phi_z[P], 1, Tp[P], 1);
                    gamma_k[P] = rho_k_x[P] * rho_x[P];
                    gamma_k[P] += rho_k_y[P] * rho_y[P];
                    gamma_k[P] += rho_k_z[P] * rho_z[P];
                    gamma_k[P] *= 2;
                }
            }
            parallel_timer_off("Derivative Properties", rank);

            // => LSDA contribution (symmetrized) <= //
            parallel_timer_on("V_XCd", rank);
            // parallel_timer_on("LSDA", rank);
            for (int P = 0; P < npoints; P++) {
                std::fill(Tp[P], Tp[P] + nlocal, 0.0);
                if (rho_a[P] < v2_rho_cutoff_) continue;
                C_DAXPY(nlocal, 0.5 * v2_rho2[P] * w[P] * rho_k[P], phi[P], 1, Tp[P], 1);
            }
            // parallel_timer_off("LSDA", rank);

            // => GGA contribution <= //
            // parallel_timer_on("GGA", rank);
            if (ansatz >= 1) {
                double* v_gamma = vals["V_GAMMA_AA"]->pointer();
                double* v2_gamma_gamma = vals["V_GAMMA_AA_GAMMA_AA"]->pointer();
                double* v2_rho_gamma = vals["V_RHO_A_GAMMA_AA"]->pointer();
                double tmp_val = 0.0, v2_val = 0.0;

                for (int P = 0; P < npoints; P++) {
                    if (rho_a[P] < v2_rho_cutoff_) continue;

                    // V contributions
                    C_DAXPY(nlocal, (0.5 * w[P] * v2_rho_gamma[P] * gamma_k[P]), phi[P], 1, Tp[P], 1);

                    // W contributions
                    v2_val = (v2_rho_gamma[P] * rho_k[P] + v2_gamma_gamma[P] * gamma_k[P]);

                    tmp_val = 2.0 * w[P] * (v_gamma[P] * rho_k_x[P] + v2_val * rho_x[P]);
                    C_DAXPY(nlocal, tmp_val, phi_x[P], 1, Tp[P], 1);

                    tmp_val = 2.0 * w[P] * (v_gamma[P] * rho_k_y[P] + v2_val * rho_y[P]);
                    C_DAXPY(nlocal, tmp_val, phi_y[P], 1, Tp[P], 1);

                    tmp_val = 2.0 * w[P] * (v_gamma[P] * rho_k_z[P] + v2_val * rho_z[P]);
                    C_DAXPY(nlocal, tmp_val, phi_z[P], 1, Tp[P], 1);
                }
            }
            // parallel_timer_off("GGA", rank);

            // Put it all together
            // parallel_timer_on("LSDA", rank);
            C_DGEMM('T', 'N', nlocal, nlocal, npoints, 1.0, phi[0], coll_funcs, Tp[0], max_functions, 0.0, Vx_localp[0],
                    max_functions);
            // parallel_timer_off("LSDA", rank);

            // Symmetrization (V is *always* Hermitian)
            for (int m = 0; m < nlocal; m++) {
                for (int n = 0; n <= m; n++) {
                    Vx_localp[m][n] = Vx_localp[n][m] = Vx_localp[m][n] + Vx_localp[n][m];
                }
            }

            // => Unpacking <= //
            double** Vxp = Vx_AO[dindex]->pointer();
            for (int ml = 0; ml < nlocal; ml++) {
                int mg = function_map[ml];
                for (int nl = 0; nl < ml; nl++) {
                    int ng = function_map[nl];
#pragma omp atomic update
                    Vxp[mg][ng] += Vx_localp[ml][nl];
#pragma omp atomic update
                    Vxp[ng][mg] += Vx_localp[ml][nl];
                }
#pragma omp atomic update
                Vxp[mg][mg] += Vx_localp[ml][ml];
            }
            parallel_timer_off("V_XCd", rank);
        }
    }

    // Set the result
    for (size_t i = 0; i < Dx.size(); i++) {
        if (Dx[i]->nirrep() != 1) {
            ret[i]->apply_symmetry(Vx_AO[i], AO2USO_);
        } else {
            ret[i]->copy(Vx_AO[i]);
        }
    }

    // Reset the workers
    for (size_t i = 0; i < num_threads_; i++) {
        functional_workers_[i]->set_deriv(old_func_deriv);
        functional_workers_[i]->allocate();
    }
    timer_off("RV: Form Vx");
}
SharedMatrix RV::compute_gradient() {
    if ((D_AO_.size() != 1)) throw PSIEXCEPTION("V: RKS should have only one D Matrix");

    if (functional_->needs_vv10()) {
        throw PSIEXCEPTION("V: RKS cannot compute VV10 gradient contribution.");
    }

    // How many atoms?
    int natom = primary_->molecule()->natom();

    // Set Hessian derivative level in properties
    int old_deriv = point_workers_[0]->deriv();

    // Thread info
    int rank = 0;

    // What local XC ansatz are we in?
    int ansatz = functional_->ansatz();

    // How many functions are there (for lda in Vtemp, T)
    int max_functions = grid_->max_functions();
    int max_points = grid_->max_points();

    // Setup the pointers
    for (size_t i = 0; i < num_threads_; i++) {
        point_workers_[i]->set_pointers(D_AO_[0]);
        point_workers_[i]->set_deriv((functional_->is_gga() || functional_->is_meta() ? 2 : 1));
    }

    // Per thread temporaries
    std::vector<SharedMatrix> G_local, U_local;
    for (size_t i = 0; i < num_threads_; i++) {
        G_local.push_back(std::make_shared<Matrix>("G Temp", natom, 3));
        U_local.push_back(std::make_shared<Matrix>("U Temp", max_points, max_functions));
    }

    std::vector<double> functionalq(num_threads_);
    std::vector<double> rhoaq(num_threads_);
    std::vector<double> rhoaxq(num_threads_);
    std::vector<double> rhoayq(num_threads_);
    std::vector<double> rhoazq(num_threads_);

// Traverse the blocks of points
#pragma omp parallel for private(rank) schedule(dynamic) num_threads(num_threads_)
    for (size_t Q = 0; Q < grid_->blocks().size(); Q++) {
// Get thread info
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif

        // Get per-rank workers
        auto block = grid_->blocks()[Q];
        auto fworker = functional_workers_[rank];
        auto pworker = point_workers_[rank];

        parallel_timer_on("Properties", rank);
        pworker->compute_points(block);
        parallel_timer_off("Properties", rank);

        parallel_timer_on("Functional", rank);
        auto& vals = fworker->compute_functional(pworker->point_values());
        parallel_timer_off("Functional", rank);

        parallel_timer_on("V_xc gradient", rank);

        // => Compute quadrature <= //
        auto qvals = dft_integrators::rks_quadrature_integrate(block, fworker, pworker);
        functionalq[rank] += qvals[0];
        rhoaq[rank] += qvals[1];
        rhoaxq[rank] += qvals[2];
        rhoayq[rank] += qvals[3];
        rhoazq[rank] += qvals[4];

        // => Integrate all contributions into G <= //
        dft_integrators::rks_gradient_integrator(primary_, block, fworker, pworker, G_local[rank], U_local[rank]);

        parallel_timer_off("V_xc gradient", rank);
    }

    // Sum up the matrix
    auto G = std::make_shared<Matrix>("XC Gradient", natom, 3);
    for (auto const& val : G_local) {
        G->add(val);
    }

    quad_values_["FUNCTIONAL"] = std::accumulate(functionalq.begin(), functionalq.end(), 0.0);
    quad_values_["RHO_A"] = std::accumulate(rhoaq.begin(), rhoaq.end(), 0.0);
    quad_values_["RHO_AX"] = std::accumulate(rhoaxq.begin(), rhoaxq.end(), 0.0);
    quad_values_["RHO_AY"] = std::accumulate(rhoayq.begin(), rhoayq.end(), 0.0);
    quad_values_["RHO_AZ"] = std::accumulate(rhoazq.begin(), rhoazq.end(), 0.0);
    quad_values_["RHO_B"] = quad_values_["RHO_A"];
    quad_values_["RHO_BX"] = quad_values_["RHO_AX"];
    quad_values_["RHO_BY"] = quad_values_["RHO_AY"];
    quad_values_["RHO_BZ"] = quad_values_["RHO_AZ"];

    if (std::isnan(quad_values_["FUNCTIONAL"])) {
        throw PSIEXCEPTION("V: Integrated DFT functional to get NaN. The functional is not numerically stable. Pick a different one.");
    }

    if (debug_) {
        outfile->Printf("   => XC Gradient: Numerical Integrals <=\n\n");
        outfile->Printf("    Functional Value:  %24.16E\n", quad_values_["FUNCTIONAL"]);
        outfile->Printf("    <\\rho_a>        :  %24.16E\n", quad_values_["RHO_A"]);
        outfile->Printf("    <\\rho_b>        :  %24.16E\n", quad_values_["RHO_B"]);
        outfile->Printf("    <\\vec r\\rho_a>  : <%24.16E,%24.16E,%24.16E>\n", quad_values_["RHO_AX"],
                        quad_values_["RHO_AY"], quad_values_["RHO_AZ"]);
        outfile->Printf("    <\\vec r\\rho_b>  : <%24.16E,%24.16E,%24.16E>\n\n", quad_values_["RHO_BX"],
                        quad_values_["RHO_BY"], quad_values_["RHO_BZ"]);
    }

    for (size_t i = 0; i < num_threads_; i++) {
        point_workers_[i]->set_deriv(old_deriv);
    }
    if (functional_->needs_vv10()) {
        G->add(vv10_nlc_gradient(D_AO_[0]));
    }

    // RKS
    G->scale(2.0);

    return G;
}

SharedMatrix RV::compute_hessian() {
    if (functional_->is_gga() || functional_->is_meta())
        throw PSIEXCEPTION("Hessians for GGA and meta GGA functionals are not yet implemented.");

    if ((D_AO_.size() != 1)) throw PSIEXCEPTION("V: RKS should have only one D Matrix");

    if (functional_->needs_vv10()) {
        throw PSIEXCEPTION("V: RKS cannot compute VV10 Hessian contribution.");
    }

    // Build the target Hessian Matrix
    int natom = primary_->molecule()->natom();
    auto H = std::make_shared<Matrix>("XC Hessian", 3 * natom, 3 * natom);
    double** Hp = H->pointer();

    // Thread info
    int rank = 0;

    // Set Hessian derivative level in properties
    int old_deriv = point_workers_[0]->deriv();
    int old_func_deriv = functional_->deriv();

    // How many functions are there (for lda in Vtemp, T)
    int max_functions = grid_->max_functions();
    int max_points = grid_->max_points();

    int derivlev = (functional_->is_gga() || functional_->is_meta()) ? 3 : 2;
    functional_->set_deriv(derivlev);

    // Setup the pointers
    for (size_t i = 0; i < num_threads_; i++) {
        point_workers_[i]->set_pointers(D_AO_[0]);
        point_workers_[i]->set_deriv(derivlev);
        functional_workers_[i]->set_deriv(derivlev);
        functional_workers_[i]->allocate();
    }

    // Per thread temporaries
    std::vector<SharedMatrix> V_local;
    std::vector<std::shared_ptr<Vector>> Q_temp;
    for (size_t i = 0; i < num_threads_; i++) {
        V_local.push_back(std::make_shared<Matrix>("V Temp", max_functions, max_functions));
        Q_temp.push_back(std::make_shared<Vector>("Quadrature Tempt", max_points));
    }

    auto QT = std::make_shared<Vector>("Quadrature Temp", max_points);
    double* QTp = QT->pointer();
    const std::vector<std::shared_ptr<BlockOPoints>>& blocks = grid_->blocks();

    for (size_t Q = 0; Q < blocks.size(); Q++) {
// Get thread info
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif

        std::shared_ptr<SuperFunctional> fworker = functional_workers_[rank];
        std::shared_ptr<PointFunctions> pworker = point_workers_[rank];
        double** V2p = V_local[rank]->pointer();
        double* QTp = Q_temp[rank]->pointer();
        double** Dp = pworker->D_scratch()[0]->pointer();
        SharedMatrix tmpHXX = pworker->D_scratch()[0]->clone();
        SharedMatrix tmpHXY = pworker->D_scratch()[0]->clone();
        SharedMatrix tmpHXZ = pworker->D_scratch()[0]->clone();
        SharedMatrix tmpHYX = pworker->D_scratch()[0]->clone();
        SharedMatrix tmpHYY = pworker->D_scratch()[0]->clone();
        SharedMatrix tmpHYZ = pworker->D_scratch()[0]->clone();
        SharedMatrix tmpHZX = pworker->D_scratch()[0]->clone();
        SharedMatrix tmpHZY = pworker->D_scratch()[0]->clone();
        SharedMatrix tmpHZZ = pworker->D_scratch()[0]->clone();
        double **pHXX = tmpHXX->pointer();
        double **pHXY = tmpHXY->pointer();
        double **pHXZ = tmpHXZ->pointer();
        double **pHYX = tmpHYX->pointer();
        double **pHYY = tmpHYY->pointer();
        double **pHYZ = tmpHYZ->pointer();
        double **pHZX = tmpHZX->pointer();
        double **pHZY = tmpHZY->pointer();
        double **pHZZ = tmpHZZ->pointer();

        // Scratch
        double** Tp = pworker->scratch()[0]->pointer();
        SharedMatrix U_local(pworker->scratch()[0]->clone());
        double** Up = U_local->pointer();

        // ACS TODO: these need to be threaded eventually, to fit in with the new infrastructure
        auto Tmpx = std::make_shared<Vector>("Tx", max_functions);
        auto Tmpy = std::make_shared<Vector>("Ty", max_functions);
        auto Tmpz = std::make_shared<Vector>("Tz", max_functions);
        double* pTx = Tmpx->pointer();
        double* pTy = Tmpy->pointer();
        double* pTz = Tmpz->pointer();
        SharedMatrix Tx(U_local->clone());
        SharedMatrix Ty(U_local->clone());
        SharedMatrix Tz(U_local->clone());
        SharedMatrix T2(U_local->clone());
        double** pTx2 = Tx->pointer();
        double** pTy2 = Ty->pointer();
        double** pTz2 = Tz->pointer();
        double** pT2 = T2->pointer();

        std::shared_ptr<BlockOPoints> block = blocks[Q];
        int npoints = block->npoints();
        double* x = block->x();
        double* y = block->y();
        double* z = block->z();
        double* w = block->w();
        const std::vector<int>& function_map = block->functions_local_to_global();
        int nlocal = function_map.size();

        pworker->compute_points(block);
        std::map<std::string, SharedVector>& vals = fworker->compute_functional(pworker->point_values(), npoints);

        double** phi = pworker->basis_value("PHI")->pointer();
        double** phi_x = pworker->basis_value("PHI_X")->pointer();
        double** phi_y = pworker->basis_value("PHI_Y")->pointer();
        double** phi_z = pworker->basis_value("PHI_Z")->pointer();
        double** phi_xx = pworker->basis_value("PHI_XX")->pointer();
        double** phi_xy = pworker->basis_value("PHI_XY")->pointer();
        double** phi_xz = pworker->basis_value("PHI_XZ")->pointer();
        double** phi_yy = pworker->basis_value("PHI_YY")->pointer();
        double** phi_yz = pworker->basis_value("PHI_YZ")->pointer();
        double** phi_zz = pworker->basis_value("PHI_ZZ")->pointer();
        double* rho_a = pworker->point_value("RHO_A")->pointer();
        double* v_rho_a = vals["V_RHO_A"]->pointer();
        double* v_rho_aa = vals["V_RHO_A_RHO_A"]->pointer();
        size_t coll_funcs = pworker->basis_value("PHI")->ncol();

        // => LSDA Contribution <= //

        for (int ml = 0; ml < nlocal; ml++) {
            std::fill(pHXX[ml], pHXX[ml] + nlocal, 0);
            std::fill(pHXY[ml], pHXY[ml] + nlocal, 0);
            std::fill(pHXZ[ml], pHXZ[ml] + nlocal, 0);
            std::fill(pHYX[ml], pHYX[ml] + nlocal, 0);
            std::fill(pHYY[ml], pHYY[ml] + nlocal, 0);
            std::fill(pHYZ[ml], pHYZ[ml] + nlocal, 0);
            std::fill(pHZX[ml], pHZX[ml] + nlocal, 0);
            std::fill(pHZY[ml], pHZY[ml] + nlocal, 0);
            std::fill(pHZZ[ml], pHZZ[ml] + nlocal, 0);
        }

        /*
         *                        mn  ∂ F
         *  H_mn <- 2 D_ab ɸ_a ɸ_b    ---
         *                            ∂ ρ
         */
        // T = ɸ D
        C_DGEMM('N', 'N', npoints, nlocal, nlocal, 1.0, phi[0], coll_funcs, Dp[0], max_functions, 0.0, Tp[0],
                max_functions);
        for (int P = 0; P < npoints; P++) {
            std::fill(Up[P], Up[P] + nlocal, 0.0);
            if (std::fabs(rho_a[P]) > v2_rho_cutoff_) {
                C_DAXPY(nlocal, 2.0 * w[P] * v_rho_a[P], Tp[P], 1, Up[P], 1);
            }
        }
        for (int ml = 0; ml < nlocal; ml++) {
            int A = primary_->function_to_center(function_map[ml]);
            double Txx = C_DDOT(npoints, &Up[0][ml], max_functions, &phi_xx[0][ml], coll_funcs);
            double Txy = C_DDOT(npoints, &Up[0][ml], max_functions, &phi_xy[0][ml], coll_funcs);
            double Txz = C_DDOT(npoints, &Up[0][ml], max_functions, &phi_xz[0][ml], coll_funcs);
            double Tyy = C_DDOT(npoints, &Up[0][ml], max_functions, &phi_yy[0][ml], coll_funcs);
            double Tyz = C_DDOT(npoints, &Up[0][ml], max_functions, &phi_yz[0][ml], coll_funcs);
            double Tzz = C_DDOT(npoints, &Up[0][ml], max_functions, &phi_zz[0][ml], coll_funcs);
            Hp[3 * A + 0][3 * A + 0] += Txx;
            Hp[3 * A + 0][3 * A + 1] += Txy;
            Hp[3 * A + 0][3 * A + 2] += Txz;
            Hp[3 * A + 1][3 * A + 0] += Txy;
            Hp[3 * A + 1][3 * A + 1] += Tyy;
            Hp[3 * A + 1][3 * A + 2] += Tyz;
            Hp[3 * A + 2][3 * A + 0] += Txz;
            Hp[3 * A + 2][3 * A + 1] += Tyz;
            Hp[3 * A + 2][3 * A + 2] += Tzz;
        }

        /*
         *                        m             n  ∂^2 F
         *  H_mn <- 4 D_ab ɸ_a ɸ_b  D_cd ɸ_c ɸ_d   ------
         *                                         ∂ ρ^2
         *  RHF prefactor gets multiplied by 4 to account for occupancy in both densities.
         *  A factor of two is applied at the end, so just double this contribution.
         */

        for (int P = 0; P < npoints; P++) {
            std::fill(Up[P], Up[P] + nlocal, 0.0);
            if (std::fabs(rho_a[P]) > v2_rho_cutoff_) {
                C_DAXPY(nlocal, 8.0 * w[P] * v_rho_aa[P], Tp[P], 1, Up[P], 1);
            }
        }
        for (int P = 0; P < npoints; P++) {
            for (int ml = 0; ml < nlocal; ml++) {
                pTx2[P][ml] = Tp[P][ml] * phi_x[P][ml];
                pTy2[P][ml] = Tp[P][ml] * phi_y[P][ml];
                pTz2[P][ml] = Tp[P][ml] * phi_z[P][ml];
            }
        }

        // x derivatives
        for (int P = 0; P < npoints; P++) {
            for (int ml = 0; ml < nlocal; ml++) {
                Tp[P][ml] = Up[P][ml] * phi_x[P][ml];
            }
        }
        C_DGEMM('T', 'N', nlocal, nlocal, npoints, 1.0, Tp[0], coll_funcs, pTx2[0], max_functions, 1.0, pHXX[0], max_functions);
        C_DGEMM('T', 'N', nlocal, nlocal, npoints, 1.0, Tp[0], coll_funcs, pTy2[0], max_functions, 1.0, pHXY[0], max_functions);
        C_DGEMM('T', 'N', nlocal, nlocal, npoints, 1.0, Tp[0], coll_funcs, pTz2[0], max_functions, 1.0, pHXZ[0], max_functions);

        // y derivatives
        for (int P = 0; P < npoints; P++) {
            for (int ml = 0; ml < nlocal; ml++) {
                Tp[P][ml] = Up[P][ml] * phi_y[P][ml];
            }
        }
        C_DGEMM('T', 'N', nlocal, nlocal, npoints, 1.0, Tp[0], coll_funcs, pTx2[0], max_functions, 1.0, pHYX[0], max_functions);
        C_DGEMM('T', 'N', nlocal, nlocal, npoints, 1.0, Tp[0], coll_funcs, pTy2[0], max_functions, 1.0, pHYY[0], max_functions);
        C_DGEMM('T', 'N', nlocal, nlocal, npoints, 1.0, Tp[0], coll_funcs, pTz2[0], max_functions, 1.0, pHYZ[0], max_functions);

        // z derivatives
        for (int P = 0; P < npoints; P++) {
            for (int ml = 0; ml < nlocal; ml++) {
                Tp[P][ml] = Up[P][ml] * phi_z[P][ml];
            }
        }
        C_DGEMM('T', 'N', nlocal, nlocal, npoints, 1.0, Tp[0], coll_funcs, pTx2[0], max_functions, 1.0, pHZX[0], max_functions);
        C_DGEMM('T', 'N', nlocal, nlocal, npoints, 1.0, Tp[0], coll_funcs, pTy2[0], max_functions, 1.0, pHZY[0], max_functions);
        C_DGEMM('T', 'N', nlocal, nlocal, npoints, 1.0, Tp[0], coll_funcs, pTz2[0], max_functions, 1.0, pHZZ[0], max_functions);

        /*
         *                    m    n  ∂ F
         *  H_mn <- 2 D_ab ɸ_a  ɸ_b   ---
         *                            ∂ ρ
         */
        // x derivatives
        for (int P = 0; P < npoints; P++) {
            std::fill(Up[P], Up[P] + nlocal, 0.0);
            if (std::fabs(rho_a[P]) > v2_rho_cutoff_) {
                C_DAXPY(nlocal, 2.0 * w[P] * v_rho_a[P], phi_x[P], 1, Up[P], 1);
            }
        }
        C_DGEMM('T', 'N', nlocal, nlocal, npoints, 1.0, phi_x[0], coll_funcs, Up[0], max_functions, 0.0, pTx2[0], max_functions);
        C_DGEMM('T', 'N', nlocal, nlocal, npoints, 1.0, phi_y[0], coll_funcs, Up[0], max_functions, 0.0, pTy2[0], max_functions);
        C_DGEMM('T', 'N', nlocal, nlocal, npoints, 1.0, phi_z[0], coll_funcs, Up[0], max_functions, 0.0, pTz2[0], max_functions);
        for (int ml = 0; ml < nlocal; ml++) {
            for (int nl = 0; nl < nlocal; nl++) {
                double D = Dp[ml][nl];
                pHXX[ml][nl] += pTx2[ml][nl] * D;
                pHYX[ml][nl] += pTy2[ml][nl] * D;
                pHZX[ml][nl] += pTz2[ml][nl] * D;
            }
        }

        // y derivatives
        for (int P = 0; P < npoints; P++) {
            std::fill(Up[P], Up[P] + nlocal, 0.0);
            if (std::fabs(rho_a[P]) > v2_rho_cutoff_) {
                C_DAXPY(nlocal, 2.0 * w[P] * v_rho_a[P], phi_y[P], 1, Up[P], 1);
            }
        }
        C_DGEMM('T', 'N', nlocal, nlocal, npoints, 1.0, phi_x[0], coll_funcs, Up[0], max_functions, 0.0, pTx2[0], max_functions);
        C_DGEMM('T', 'N', nlocal, nlocal, npoints, 1.0, phi_y[0], coll_funcs, Up[0], max_functions, 0.0, pTy2[0], max_functions);
        C_DGEMM('T', 'N', nlocal, nlocal, npoints, 1.0, phi_z[0], coll_funcs, Up[0], max_functions, 0.0, pTz2[0], max_functions);
        for (int ml = 0; ml < nlocal; ml++) {
            for (int nl = 0; nl < nlocal; nl++) {
                double D = Dp[ml][nl];
                pHXY[ml][nl] += pTx2[ml][nl] * D;
                pHYY[ml][nl] += pTy2[ml][nl] * D;
                pHZY[ml][nl] += pTz2[ml][nl] * D;
            }
        }

        // z derivatives
        for (int P = 0; P < npoints; P++) {
            std::fill(Up[P], Up[P] + nlocal, 0.0);
            if (std::fabs(rho_a[P]) > v2_rho_cutoff_) {
                C_DAXPY(nlocal, 2.0 * w[P] * v_rho_a[P], phi_z[P], 1, Up[P], 1);
            }
        }
        C_DGEMM('T', 'N', nlocal, nlocal, npoints, 1.0, phi_x[0], coll_funcs, Up[0], max_functions, 0.0, pTx2[0], max_functions);
        C_DGEMM('T', 'N', nlocal, nlocal, npoints, 1.0, phi_y[0], coll_funcs, Up[0], max_functions, 0.0, pTy2[0], max_functions);
        C_DGEMM('T', 'N', nlocal, nlocal, npoints, 1.0, phi_z[0], coll_funcs, Up[0], max_functions, 0.0, pTz2[0], max_functions);
        for (int ml = 0; ml < nlocal; ml++) {
            for (int nl = 0; nl < nlocal; nl++) {
                double D = Dp[ml][nl];
                pHXZ[ml][nl] += pTx2[ml][nl] * D;
                pHYZ[ml][nl] += pTy2[ml][nl] * D;
                pHZZ[ml][nl] += pTz2[ml][nl] * D;
            }
        }

        // Accumulate contributions to the full Hessian: N.B. these terms are not symmetric!
        for (int ml = 0; ml < nlocal; ml++) {
            int A = primary_->function_to_center(function_map[ml]);
            for (int nl = 0; nl < nlocal; nl++) {
                int B = primary_->function_to_center(function_map[nl]);
                Hp[3 * A + 0][3 * B + 0] += pHXX[ml][nl];
                Hp[3 * A + 1][3 * B + 0] += pHYX[ml][nl];
                Hp[3 * A + 2][3 * B + 0] += pHZX[ml][nl];
                Hp[3 * A + 0][3 * B + 1] += pHXY[ml][nl];
                Hp[3 * A + 1][3 * B + 1] += pHYY[ml][nl];
                Hp[3 * A + 2][3 * B + 1] += pHZY[ml][nl];
                Hp[3 * A + 0][3 * B + 2] += pHXZ[ml][nl];
                Hp[3 * A + 1][3 * B + 2] += pHYZ[ml][nl];
                Hp[3 * A + 2][3 * B + 2] += pHZZ[ml][nl];
            }
        }
    }

    if (std::isnan(quad_values_["FUNCTIONAL"])) {
        throw PSIEXCEPTION("V: Integrated DFT functional to get NaN. The functional is not numerically stable. Pick a different one.");
    }

    if (debug_) {
        outfile->Printf("   => XC Hessian: Numerical Integrals <=\n\n");
        outfile->Printf("    Functional Value:  %24.16E\n", quad_values_["FUNCTIONAL"]);
        outfile->Printf("    <\\rho_a>        :  %24.16E\n", quad_values_["RHO_A"]);
        outfile->Printf("    <\\rho_b>        :  %24.16E\n", quad_values_["RHO_B"]);
        outfile->Printf("    <\\vec r\\rho_a>  : <%24.16E,%24.16E,%24.16E>\n", quad_values_["RHO_AX"],
                        quad_values_["RHO_AY"], quad_values_["RHO_AZ"]);
        outfile->Printf("    <\\vec r\\rho_b>  : <%24.16E,%24.16E,%24.16E>\n\n", quad_values_["RHO_BX"],
                        quad_values_["RHO_BY"], quad_values_["RHO_BZ"]);
    }

    for (size_t i = 0; i < num_threads_; i++) {
        point_workers_[i]->set_deriv(old_deriv);
    }
    functional_->set_deriv(old_func_deriv);

    // RKS
    H->scale(2.0);
    H->hermitivitize();
//    H->print_out();

    return H;
}

UV::UV(std::shared_ptr<SuperFunctional> functional, std::shared_ptr<BasisSet> primary, Options& options)
    : VBase(functional, primary, options) {}
UV::~UV() {}
void UV::initialize() {
    VBase::initialize();
    int max_points = grid_->max_points();
    int max_functions = grid_->max_functions();
    for (size_t i = 0; i < num_threads_; i++) {
        // Need a points worker per thread
        std::shared_ptr<PointFunctions> point_tmp = std::make_shared<UKSFunctions>(primary_, max_points, max_functions);
        point_tmp->set_ansatz(functional_->ansatz());
        point_tmp->set_cache_map(&cache_map_);
        point_workers_.push_back(point_tmp);
    }
}
void UV::finalize() { VBase::finalize(); }
void UV::print_header() const { VBase::print_header(); }
void UV::compute_V(std::vector<SharedMatrix> ret) {
    timer_on("UV: Form V");
    if ((D_AO_.size() != 2) || (ret.size() != 2)) {
        throw PSIEXCEPTION("V: UKS should have two D/V Matrices");
    }
    
#ifdef USING_BrianQC
    if (brianEnable and brianEnableDFT) {
        double DFTEnergy;
        brianSCFBuildFockDFT(&brianCookie,
            D_AO_[0]->get_pointer(0),
            D_AO_[1]->get_pointer(0),
            ret[0]->get_pointer(0),
            ret[1]->get_pointer(0),
            &DFTEnergy
        );
        checkBrian();
        
        quad_values_["VV10"] = 0.0; // NOTE: BrianQC doesn't compute the VV10 term separately, it just includes it in the DFT energy term
        quad_values_["FUNCTIONAL"] = DFTEnergy;
        quad_values_["RHO_A"] = 0.0;
        quad_values_["RHO_AX"] = 0.0;
        quad_values_["RHO_AY"] = 0.0;
        quad_values_["RHO_AZ"] = 0.0;
        quad_values_["RHO_B"] = 0.0;
        quad_values_["RHO_BX"] = 0.0;
        quad_values_["RHO_BY"] = 0.0;
        quad_values_["RHO_BZ"] = 0.0;
        
        return;
    }
#endif

    if (functional_->needs_grac()) {
        throw PSIEXCEPTION("V: UKS cannot compute GRAC corrections.");
    }

    // Thread info
    int rank = 0;

    // What local XC ansatz are we in?
    int ansatz = functional_->ansatz();

    // How many functions are there (for lda in Vtemp, T)
    int max_functions = grid_->max_functions();
    int max_points = grid_->max_points();

    // Setup the pointers
    for (size_t i = 0; i < num_threads_; i++) {
        point_workers_[i]->set_pointers(D_AO_[0], D_AO_[1]);
    }

    // Per thread temporaries
    std::vector<SharedMatrix> Va_local, Vb_local;
    std::vector<std::shared_ptr<Vector>> Qa_temp, Qb_temp;
    for (size_t i = 0; i < num_threads_; i++) {
        Va_local.push_back(std::make_shared<Matrix>("Va Temp", max_functions, max_functions));
        Vb_local.push_back(std::make_shared<Matrix>("Vb Temp", max_functions, max_functions));
        Qa_temp.push_back(std::make_shared<Vector>("Quadrature A Temp", max_points));
        Qb_temp.push_back(std::make_shared<Vector>("Quadrature B Temp", max_points));
    }

    auto Va_AO = std::make_shared<Matrix>("Va Temp", nbf_, nbf_);
    auto Vb_AO = std::make_shared<Matrix>("Vb Temp", nbf_, nbf_);
    double** Vap = Va_AO->pointer();
    double** Vbp = Vb_AO->pointer();

    std::vector<double> functionalq(num_threads_);
    std::vector<double> rhoaq(num_threads_);
    std::vector<double> rhoaxq(num_threads_);
    std::vector<double> rhoayq(num_threads_);
    std::vector<double> rhoazq(num_threads_);
    std::vector<double> rhobq(num_threads_);
    std::vector<double> rhobxq(num_threads_);
    std::vector<double> rhobyq(num_threads_);
    std::vector<double> rhobzq(num_threads_);

    // Loop over grid
    for (size_t Q = 0; Q < grid_->blocks().size(); Q++) {
// Get thread info
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif

        std::shared_ptr<SuperFunctional> fworker = functional_workers_[rank];
        std::shared_ptr<PointFunctions> pworker = point_workers_[rank];
        double** Va2p = Va_local[rank]->pointer();
        double** Vb2p = Vb_local[rank]->pointer();
        double* QTap = Qa_temp[rank]->pointer();
        double* QTbp = Qb_temp[rank]->pointer();

        // Scratch
        double** Tap = pworker->scratch()[0]->pointer();
        double** Tbp = pworker->scratch()[1]->pointer();

        std::shared_ptr<BlockOPoints> block = grid_->blocks()[Q];
        int npoints = block->npoints();
        double* x = block->x();
        double* y = block->y();
        double* z = block->z();
        double* w = block->w();
        const std::vector<int>& function_map = block->functions_local_to_global();
        int nlocal = function_map.size();

        parallel_timer_on("Properties", rank);
        pworker->compute_points(block, false);
        parallel_timer_off("Properties", rank);

        parallel_timer_on("Functional", rank);
        std::map<std::string, SharedVector>& vals = fworker->compute_functional(pworker->point_values(), npoints);
        parallel_timer_off("Functional", rank);

        if (debug_ > 3) {
            block->print("outfile", debug_);
            pworker->print("outfile", debug_);
        }

        parallel_timer_on("V_xc", rank);
        double** phi = pworker->basis_value("PHI")->pointer();
        double* rho_a = pworker->point_value("RHO_A")->pointer();
        double* rho_b = pworker->point_value("RHO_B")->pointer();
        double* zk = vals["V"]->pointer();
        double* v_rho_a = vals["V_RHO_A"]->pointer();
        double* v_rho_b = vals["V_RHO_B"]->pointer();
        size_t coll_funcs = pworker->basis_value("PHI")->ncol();

        // => Quadrature values <= //
        functionalq[rank] += C_DDOT(npoints, w, 1, zk, 1);
        for (int P = 0; P < npoints; P++) {
            QTap[P] = w[P] * rho_a[P];
            QTbp[P] = w[P] * rho_b[P];
        }
        rhoaq[rank] += C_DDOT(npoints, w, 1, rho_a, 1);
        rhoaxq[rank] += C_DDOT(npoints, QTap, 1, x, 1);
        rhoayq[rank] += C_DDOT(npoints, QTap, 1, y, 1);
        rhoazq[rank] += C_DDOT(npoints, QTap, 1, z, 1);
        rhobq[rank] += C_DDOT(npoints, w, 1, rho_b, 1);
        rhobxq[rank] += C_DDOT(npoints, QTbp, 1, x, 1);
        rhobyq[rank] += C_DDOT(npoints, QTbp, 1, y, 1);
        rhobzq[rank] += C_DDOT(npoints, QTbp, 1, z, 1);

        // => LSDA contribution (symmetrized) <= //
        // timer_on("V: LSDA");
        for (int P = 0; P < npoints; P++) {
            std::fill(Tap[P], Tap[P] + nlocal, 0.0);
            std::fill(Tbp[P], Tbp[P] + nlocal, 0.0);
            C_DAXPY(nlocal, 0.5 * v_rho_a[P] * w[P], phi[P], 1, Tap[P], 1);
            C_DAXPY(nlocal, 0.5 * v_rho_b[P] * w[P], phi[P], 1, Tbp[P], 1);
        }
        // timer_off("V: LSDA");

        // => GGA contribution (symmetrized) <= //
        if (ansatz >= 1) {
            // timer_on("V: GGA");
            double** phix = pworker->basis_value("PHI_X")->pointer();
            double** phiy = pworker->basis_value("PHI_Y")->pointer();
            double** phiz = pworker->basis_value("PHI_Z")->pointer();
            double* rho_ax = pworker->point_value("RHO_AX")->pointer();
            double* rho_ay = pworker->point_value("RHO_AY")->pointer();
            double* rho_az = pworker->point_value("RHO_AZ")->pointer();
            double* rho_bx = pworker->point_value("RHO_BX")->pointer();
            double* rho_by = pworker->point_value("RHO_BY")->pointer();
            double* rho_bz = pworker->point_value("RHO_BZ")->pointer();
            double* v_sigma_aa = vals["V_GAMMA_AA"]->pointer();
            double* v_sigma_ab = vals["V_GAMMA_AB"]->pointer();
            double* v_sigma_bb = vals["V_GAMMA_BB"]->pointer();

            for (int P = 0; P < npoints; P++) {
                C_DAXPY(nlocal, w[P] * (2.0 * v_sigma_aa[P] * rho_ax[P] + v_sigma_ab[P] * rho_bx[P]), phix[P], 1,
                        Tap[P], 1);
                C_DAXPY(nlocal, w[P] * (2.0 * v_sigma_aa[P] * rho_ay[P] + v_sigma_ab[P] * rho_by[P]), phiy[P], 1,
                        Tap[P], 1);
                C_DAXPY(nlocal, w[P] * (2.0 * v_sigma_aa[P] * rho_az[P] + v_sigma_ab[P] * rho_bz[P]), phiz[P], 1,
                        Tap[P], 1);
                C_DAXPY(nlocal, w[P] * (2.0 * v_sigma_bb[P] * rho_bx[P] + v_sigma_ab[P] * rho_ax[P]), phix[P], 1,
                        Tbp[P], 1);
                C_DAXPY(nlocal, w[P] * (2.0 * v_sigma_bb[P] * rho_by[P] + v_sigma_ab[P] * rho_ay[P]), phiy[P], 1,
                        Tbp[P], 1);
                C_DAXPY(nlocal, w[P] * (2.0 * v_sigma_bb[P] * rho_bz[P] + v_sigma_ab[P] * rho_az[P]), phiz[P], 1,
                        Tbp[P], 1);
            }
            // timer_off("V: GGA");
        }

        // timer_on("V: LSDA");
        // Single GEMM slams GGA+LSDA together (man but GEM's hot!)
        C_DGEMM('T', 'N', nlocal, nlocal, npoints, 1.0, phi[0], coll_funcs, Tap[0], max_functions, 0.0, Va2p[0],
                max_functions);
        C_DGEMM('T', 'N', nlocal, nlocal, npoints, 1.0, phi[0], coll_funcs, Tbp[0], max_functions, 0.0, Vb2p[0],
                max_functions);

        // Symmetrization (V is Hermitian)
        for (int m = 0; m < nlocal; m++) {
            for (int n = 0; n <= m; n++) {
                Va2p[m][n] = Va2p[n][m] = Va2p[m][n] + Va2p[n][m];
                Vb2p[m][n] = Vb2p[n][m] = Vb2p[m][n] + Vb2p[n][m];
            }
        }
        // timer_off("V: LSDA");

        // => Meta contribution <= //
        if (ansatz >= 2) {
            // timer_on("V: Meta");
            double** phix = pworker->basis_value("PHI_X")->pointer();
            double** phiy = pworker->basis_value("PHI_Y")->pointer();
            double** phiz = pworker->basis_value("PHI_Z")->pointer();
            double* v_tau_a = vals["V_TAU_A"]->pointer();
            double* v_tau_b = vals["V_TAU_B"]->pointer();

            double** phi[3];
            phi[0] = phix;
            phi[1] = phiy;
            phi[2] = phiz;

            double* v_tau[2];
            v_tau[0] = v_tau_a;
            v_tau[1] = v_tau_b;

            double** V_val[2];
            V_val[0] = Va2p;
            V_val[1] = Vb2p;

            for (int s = 0; s < 2; s++) {
                double** V2p = V_val[s];
                double* v_taup = v_tau[s];
                for (int i = 0; i < 3; i++) {
                    double** phiw = phi[i];
                    for (int P = 0; P < npoints; P++) {
                        std::fill(Tap[P], Tap[P] + nlocal, 0.0);
                        C_DAXPY(nlocal, v_taup[P] * w[P], phiw[P], 1, Tap[P], 1);
                    }
                    C_DGEMM('T', 'N', nlocal, nlocal, npoints, 1.0, phiw[0], coll_funcs, Tap[0], max_functions, 1.0,
                            V2p[0], max_functions);
                }
            }

            // timer_off("V: Meta");
        }

        // => Unpacking <= //
        for (int ml = 0; ml < nlocal; ml++) {
            int mg = function_map[ml];
            for (int nl = 0; nl < ml; nl++) {
                int ng = function_map[nl];
#pragma omp atomic update
                Vap[mg][ng] += Va2p[ml][nl];
#pragma omp atomic update
                Vap[ng][mg] += Va2p[ml][nl];
#pragma omp atomic update
                Vbp[mg][ng] += Vb2p[ml][nl];
#pragma omp atomic update
                Vbp[ng][mg] += Vb2p[ml][nl];
            }
#pragma omp atomic update
            Vap[mg][mg] += Va2p[ml][ml];
#pragma omp atomic update
            Vbp[mg][mg] += Vb2p[ml][ml];
        }
        parallel_timer_off("V_xc", rank);
    }

    // Do we need VV10?
    double vv10_e = 0.0;
    if (functional_->needs_vv10()) {
        SharedMatrix Ds = D_AO_[0]->clone();
        Ds->axpy(1.0, D_AO_[1]);
        Ds->scale(0.5);  // Will be scaled by a factor of 2 in vv10_nlc

        SharedMatrix V_vv10 = Ds->clone();
        V_vv10->zero();

        vv10_e = vv10_nlc(Ds, V_vv10);

        Va_AO->add(V_vv10);
        Vb_AO->add(V_vv10);
    }

    // Set the result
    if (AO2USO_) {
        ret[0]->apply_symmetry(Va_AO, AO2USO_);
        ret[1]->apply_symmetry(Vb_AO, AO2USO_);
    } else {
        ret[0]->copy(Va_AO);
        ret[1]->copy(Vb_AO);
    }

    quad_values_["VV10"] = vv10_e;
    quad_values_["FUNCTIONAL"] = std::accumulate(functionalq.begin(), functionalq.end(), 0.0);
    quad_values_["RHO_A"] = std::accumulate(rhoaq.begin(), rhoaq.end(), 0.0);
    quad_values_["RHO_AX"] = std::accumulate(rhoaxq.begin(), rhoaxq.end(), 0.0);
    quad_values_["RHO_AY"] = std::accumulate(rhoayq.begin(), rhoayq.end(), 0.0);
    quad_values_["RHO_AZ"] = std::accumulate(rhoazq.begin(), rhoazq.end(), 0.0);
    quad_values_["RHO_B"] = std::accumulate(rhobq.begin(), rhobq.end(), 0.0);
    quad_values_["RHO_BX"] = std::accumulate(rhobxq.begin(), rhobxq.end(), 0.0);
    quad_values_["RHO_BY"] = std::accumulate(rhobyq.begin(), rhobyq.end(), 0.0);
    quad_values_["RHO_BZ"] = std::accumulate(rhobzq.begin(), rhobzq.end(), 0.0);

    if (std::isnan(quad_values_["FUNCTIONAL"])) {
        throw PSIEXCEPTION("V: Integrated DFT functional to get NaN. The functional is not numerically stable. Pick a different one.");
    }

    if (debug_) {
        outfile->Printf("   => Numerical Integrals <=\n\n");
        outfile->Printf("    Functional Value:  %24.16E\n", quad_values_["FUNCTIONAL"]);
        outfile->Printf("    <\\rho_a>        :  %24.16E\n", quad_values_["RHO_A"]);
        outfile->Printf("    <\\rho_b>        :  %24.16E\n", quad_values_["RHO_B"]);
        outfile->Printf("    <\\vec r\\rho_a>  : <%24.16E,%24.16E,%24.16E>\n", quad_values_["RHO_AX"],
                        quad_values_["RHO_AY"], quad_values_["RHO_AZ"]);
        outfile->Printf("    <\\vec r\\rho_b>  : <%24.16E,%24.16E,%24.16E>\n\n", quad_values_["RHO_BX"],
                        quad_values_["RHO_BY"], quad_values_["RHO_BZ"]);
    }
    timer_off("UV: Form V");
}
void UV::compute_Vx(std::vector<SharedMatrix> Dx, std::vector<SharedMatrix> ret) {
    timer_on("UV: Form Vx");
    if (D_AO_.size() != 2) {
        throw PSIEXCEPTION("Vx: UKS should have two D matrices.");
    }
    if ((Dx.size() != ret.size()) || (Dx.size() == 0)) {
        throw PSIEXCEPTION("Vx: UKS input and output sizes should be the same.");
    }
    if ((Dx.size() % 2) != 0) {
        throw PSIEXCEPTION("Vx: UKS input and output sizes should be the same.");
    }

    if (functional_->needs_vv10()) {
        throw PSIEXCEPTION("V: UKS cannot compute VV10 Vx contribution.");
    }

    // Thread info
    int rank = 0;

    // What local XC ansatz are we in?
    int ansatz = functional_->ansatz();
    if (ansatz >= 2) {
        throw PSIEXCEPTION("Vx: UKS does not support rotated V builds for MGGA's");
    }

    int old_point_deriv = point_workers_[0]->deriv();
    int old_func_deriv = functional_->deriv();

    // How many functions are there (for lda in Vtemp, T)
    int max_functions = grid_->max_functions();
    int max_points = grid_->max_points();

    // Set pointers to SCF density
    for (size_t i = 0; i < num_threads_; i++) {
        point_workers_[i]->set_pointers(D_AO_[0], D_AO_[1]);
    }

    std::vector<SharedMatrix> Dx_vec;
    for (size_t i = 0; i < Dx.size(); i++) {
        if (Dx[i]->nirrep() != 1) {
            auto Dx_mat = std::make_shared<Matrix>("D AO temp", nbf_, nbf_);
            Dx_mat->remove_symmetry(Dx[i], USO2AO_);
            Dx_vec.push_back(Dx_mat);
        } else {
            Dx_vec.push_back(Dx[i]);
        }
    }

    // Per [R]ank quantities
    std::vector<SharedMatrix> R_Vax_local, R_Dax_local;
    std::vector<SharedMatrix> R_Vbx_local, R_Dbx_local;
    std::vector<std::shared_ptr<Vector>> R_rho_ak, R_rho_ak_x, R_rho_ak_y, R_rho_ak_z, R_gamma_ak;
    std::vector<std::shared_ptr<Vector>> R_rho_bk, R_rho_bk_x, R_rho_bk_y, R_rho_bk_z, R_gamma_bk;
    std::vector<std::shared_ptr<Vector>> R_gamma_abk;
    for (size_t i = 0; i < num_threads_; i++) {
        R_Vax_local.push_back(std::make_shared<Matrix>("Vax Temp", max_functions, max_functions));
        R_Vbx_local.push_back(std::make_shared<Matrix>("Vbx Temp", max_functions, max_functions));
        R_Dax_local.push_back(std::make_shared<Matrix>("Dak Temp", max_functions, max_functions));
        R_Dbx_local.push_back(std::make_shared<Matrix>("Dbk Temp", max_functions, max_functions));

        R_rho_ak.push_back(std::make_shared<Vector>("Rho aK Temp", max_points));
        R_rho_bk.push_back(std::make_shared<Vector>("Rho bK Temp", max_points));

        if (ansatz >= 1) {
            R_rho_ak_x.push_back(std::make_shared<Vector>("RHO K X Temp", max_points));
            R_rho_ak_y.push_back(std::make_shared<Vector>("RHO K Y Temp", max_points));
            R_rho_ak_z.push_back(std::make_shared<Vector>("Rho K Z Temp", max_points));
            R_gamma_ak.push_back(std::make_shared<Vector>("Gamma K Temp", max_points));

            R_rho_bk_x.push_back(std::make_shared<Vector>("RHO K X Temp", max_points));
            R_rho_bk_y.push_back(std::make_shared<Vector>("RHO K Y Temp", max_points));
            R_rho_bk_z.push_back(std::make_shared<Vector>("Rho K Z Temp", max_points));
            R_gamma_bk.push_back(std::make_shared<Vector>("Gamma K Temp", max_points));

            R_gamma_abk.push_back(std::make_shared<Vector>("Gamma K Temp", max_points));
        }

        functional_workers_[i]->set_deriv(2);
        functional_workers_[i]->allocate();
    }

    // Output quantities
    std::vector<SharedMatrix> Vax_AO;
    std::vector<SharedMatrix> Vbx_AO;
    for (size_t i = 0; i < Dx.size(); i++) {
        Vbx_AO.push_back(std::make_shared<Matrix>("Vax AO Temp", nbf_, nbf_));
        Vax_AO.push_back(std::make_shared<Matrix>("Vbx AO Temp", nbf_, nbf_));
    }

// Traverse the blocks of points
#pragma omp parallel for private(rank) schedule(guided) num_threads(num_threads_)
    for (size_t Q = 0; Q < grid_->blocks().size(); Q++) {
// Get thread info
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif

        // => Setup <= //
        std::shared_ptr<SuperFunctional> fworker = functional_workers_[rank];
        std::shared_ptr<PointFunctions> pworker = point_workers_[rank];
        double** Vax_localp = R_Vax_local[rank]->pointer();
        double** Vbx_localp = R_Vbx_local[rank]->pointer();
        double** Dax_localp = R_Dax_local[rank]->pointer();
        double** Dbx_localp = R_Dbx_local[rank]->pointer();

        // => Compute blocks <= //
        double** Tap = pworker->scratch()[0]->pointer();
        double** Tbp = pworker->scratch()[1]->pointer();

        std::shared_ptr<BlockOPoints> block = grid_->blocks()[Q];
        int npoints = block->npoints();
        double* w = block->w();
        const std::vector<int>& function_map = block->functions_local_to_global();
        int nlocal = function_map.size();

        // Compute Rho, Phi, etc
        parallel_timer_on("Properties", rank);
        pworker->compute_points(block);
        parallel_timer_off("Properties", rank);

        // Compute functional values
        parallel_timer_on("Functional", rank);
        std::map<std::string, SharedVector>& vals = fworker->compute_functional(pworker->point_values(), npoints);
        parallel_timer_off("Functional", rank);

        // => Grab quantities <= //
        // LDA
        double** phi = pworker->basis_value("PHI")->pointer();
        double* rho_a = pworker->point_value("RHO_A")->pointer();
        double* rho_b = pworker->point_value("RHO_B")->pointer();
        double* v2_rho2_aa = vals["V_RHO_A_RHO_A"]->pointer();
        double* v2_rho2_ab = vals["V_RHO_A_RHO_B"]->pointer();
        double* v2_rho2_bb = vals["V_RHO_B_RHO_B"]->pointer();
        size_t coll_funcs = pworker->basis_value("PHI")->ncol();

        double* rho_ak = R_rho_ak[rank]->pointer();
        double* rho_bk = R_rho_bk[rank]->pointer();

        // GGA
        double** phi_x;
        double** phi_y;
        double** phi_z;

        double *rho_ak_x, *rho_bk_x;
        double *rho_ak_y, *rho_bk_y;
        double *rho_ak_z, *rho_bk_z;
        double *gamma_aak, *gamma_bbk;
        double *rho_ax, *rho_bx;
        double *rho_ay, *rho_by;
        double *rho_az, *rho_bz;
        double* gamma_abk;
        if (ansatz >= 1) {
            // Phi
            phi_x = pworker->basis_value("PHI_X")->pointer();
            phi_y = pworker->basis_value("PHI_Y")->pointer();
            phi_z = pworker->basis_value("PHI_Z")->pointer();

            // Alpha
            rho_ak_x = R_rho_ak_x[rank]->pointer();
            rho_ak_y = R_rho_ak_y[rank]->pointer();
            rho_ak_z = R_rho_ak_z[rank]->pointer();
            gamma_aak = R_gamma_ak[rank]->pointer();
            rho_ax = pworker->point_value("RHO_AX")->pointer();
            rho_ay = pworker->point_value("RHO_AY")->pointer();
            rho_az = pworker->point_value("RHO_AZ")->pointer();

            // Beta
            rho_bk_x = R_rho_bk_x[rank]->pointer();
            rho_bk_y = R_rho_bk_y[rank]->pointer();
            rho_bk_z = R_rho_bk_z[rank]->pointer();
            gamma_bbk = R_gamma_bk[rank]->pointer();
            rho_bx = pworker->point_value("RHO_AX")->pointer();
            rho_by = pworker->point_value("RHO_AY")->pointer();
            rho_bz = pworker->point_value("RHO_AZ")->pointer();

            gamma_abk = R_gamma_abk[rank]->pointer();
        }

        // Meta
        // Forget that!

        // Loop over perturbation tensors
        for (size_t dindex = 0; dindex < (Dx_vec.size() / 2); dindex++) {
            double** Daxp = Dx_vec[2 * dindex]->pointer();
            double** Dbxp = Dx_vec[2 * dindex + 1]->pointer();

            // => Build Rotated Densities <= //
            for (int ml = 0; ml < nlocal; ml++) {
                int mg = function_map[ml];
                for (int nl = 0; nl < nlocal; nl++) {
                    int ng = function_map[nl];
                    Dax_localp[ml][nl] = Daxp[mg][ng];
                    Dbx_localp[ml][nl] = Dbxp[mg][ng];
                    // printf("%d %d: %16.8f %16.8f\n", ml, nl, Dax_localp[ml][nl], Dbx_localp[ml][nl]);
                }
            }

            // Rho_a = D^k_xy phi_xa phi_ya
            // Alpha
            parallel_timer_on("Derivative Properties", rank);
            C_DGEMM('N', 'N', npoints, nlocal, nlocal, 1.0, phi[0], coll_funcs, Dax_localp[0], max_functions, 0.0,
                    Tap[0], max_functions);
            C_DGEMM('N', 'T', npoints, nlocal, nlocal, 1.0, phi[0], coll_funcs, Dax_localp[0], max_functions, 1.0,
                    Tap[0], max_functions);

            // Beta
            C_DGEMM('N', 'N', npoints, nlocal, nlocal, 1.0, phi[0], coll_funcs, Dbx_localp[0], max_functions, 0.0,
                    Tbp[0], max_functions);
            C_DGEMM('N', 'T', npoints, nlocal, nlocal, 1.0, phi[0], coll_funcs, Dbx_localp[0], max_functions, 1.0,
                    Tbp[0], max_functions);

            for (int P = 0; P < npoints; P++) {
                rho_ak[P] = 0.5 * C_DDOT(nlocal, phi[P], 1, Tap[P], 1);
                rho_bk[P] = 0.5 * C_DDOT(nlocal, phi[P], 1, Tbp[P], 1);
            }

            // Rho^d_k and gamma_k
            if (ansatz >= 1) {
                for (int P = 0; P < npoints; P++) {
                    // Alpha
                    rho_ak_x[P] = C_DDOT(nlocal, phi_x[P], 1, Tap[P], 1);
                    rho_ak_y[P] = C_DDOT(nlocal, phi_y[P], 1, Tap[P], 1);
                    rho_ak_z[P] = C_DDOT(nlocal, phi_z[P], 1, Tap[P], 1);
                    gamma_aak[P] = rho_ak_x[P] * rho_ax[P];
                    gamma_aak[P] += rho_ak_y[P] * rho_ay[P];
                    gamma_aak[P] += rho_ak_z[P] * rho_az[P];
                    gamma_aak[P] *= 2.0;

                    // Beta
                    rho_bk_x[P] = C_DDOT(nlocal, phi_x[P], 1, Tbp[P], 1);
                    rho_bk_y[P] = C_DDOT(nlocal, phi_y[P], 1, Tbp[P], 1);
                    rho_bk_z[P] = C_DDOT(nlocal, phi_z[P], 1, Tbp[P], 1);
                    gamma_bbk[P] = rho_bk_x[P] * rho_bx[P];
                    gamma_bbk[P] += rho_bk_y[P] * rho_by[P];
                    gamma_bbk[P] += rho_ak_z[P] * rho_bz[P];
                    gamma_bbk[P] *= 2.0;

                    // Alpha-Beta
                    gamma_abk[P] = rho_ak_x[P] * rho_bx[P] + rho_bk_x[P] * rho_ax[P];
                    gamma_abk[P] += rho_ak_y[P] * rho_by[P] + rho_bk_y[P] * rho_ay[P];
                    gamma_abk[P] += rho_ak_z[P] * rho_bz[P] + rho_bk_z[P] * rho_az[P];
                }
            }
            parallel_timer_off("Derivative Properties", rank);

            parallel_timer_on("V_XCd", rank);
            // => LSDA contribution (symmetrized) <= //
            double tmp_val = 0.0, tmp_ab_val = 0.0;
            for (int P = 0; P < npoints; P++) {
                std::fill(Tap[P], Tap[P] + nlocal, 0.0);
                std::fill(Tbp[P], Tbp[P] + nlocal, 0.0);

                if (rho_a[P] > v2_rho_cutoff_) {
                    tmp_val = v2_rho2_aa[P] * rho_ak[P];
                    tmp_val += v2_rho2_ab[P] * rho_bk[P];
                    tmp_val *= 0.5 * w[P];
                    C_DAXPY(nlocal, tmp_val, phi[P], 1, Tap[P], 1);
                }

                if (rho_b[P] > v2_rho_cutoff_) {
                    tmp_val = v2_rho2_bb[P] * rho_bk[P];
                    tmp_val += v2_rho2_ab[P] * rho_ak[P];
                    tmp_val *= 0.5 * w[P];
                    C_DAXPY(nlocal, tmp_val, phi[P], 1, Tbp[P], 1);
                }
            }

            // // => GGA contribution <= //
            if (ansatz >= 1) {
                double* gamma_aa = pworker->point_value("GAMMA_AA")->pointer();
                double* gamma_ab = pworker->point_value("GAMMA_AB")->pointer();
                double* gamma_bb = pworker->point_value("GAMMA_BB")->pointer();

                double* v_gamma_aa = vals["V_GAMMA_AA"]->pointer();
                double* v_gamma_ab = vals["V_GAMMA_AB"]->pointer();
                double* v_gamma_bb = vals["V_GAMMA_BB"]->pointer();

                double* v2_gamma_aa_gamma_aa = vals["V_GAMMA_AA_GAMMA_AA"]->pointer();
                double* v2_gamma_aa_gamma_ab = vals["V_GAMMA_AA_GAMMA_AB"]->pointer();
                double* v2_gamma_aa_gamma_bb = vals["V_GAMMA_AA_GAMMA_BB"]->pointer();
                double* v2_gamma_ab_gamma_ab = vals["V_GAMMA_AB_GAMMA_AB"]->pointer();
                double* v2_gamma_ab_gamma_bb = vals["V_GAMMA_AB_GAMMA_BB"]->pointer();
                double* v2_gamma_bb_gamma_bb = vals["V_GAMMA_BB_GAMMA_BB"]->pointer();

                double* v2_rho_a_gamma_aa = vals["V_RHO_A_GAMMA_AA"]->pointer();
                double* v2_rho_a_gamma_ab = vals["V_RHO_A_GAMMA_AB"]->pointer();
                double* v2_rho_a_gamma_bb = vals["V_RHO_A_GAMMA_BB"]->pointer();
                double* v2_rho_b_gamma_aa = vals["V_RHO_B_GAMMA_AA"]->pointer();
                double* v2_rho_b_gamma_ab = vals["V_RHO_B_GAMMA_AB"]->pointer();
                double* v2_rho_b_gamma_bb = vals["V_RHO_B_GAMMA_BB"]->pointer();

                double tmp_val = 0.0, v2_val_aa = 0.0, v2_val_ab = 0.0, v2_val_bb = 0.0;

                // This one is a doozy
                for (int P = 0; P < npoints; P++) {
                    // V alpha contributions
                    if (rho_a[P] > v2_rho_cutoff_) {
                        tmp_val = v2_rho_a_gamma_aa[P] * gamma_aak[P];
                        tmp_val += v2_rho_a_gamma_ab[P] * gamma_abk[P];
                        tmp_val += v2_rho_a_gamma_bb[P] * gamma_bbk[P];
                        C_DAXPY(nlocal, (0.5 * w[P] * tmp_val), phi[P], 1, Tap[P], 1);
                    }

                    // V beta contributions
                    if (rho_b[P] > v2_rho_cutoff_) {
                        tmp_val = v2_rho_b_gamma_aa[P] * gamma_aak[P];
                        tmp_val += v2_rho_b_gamma_ab[P] * gamma_abk[P];
                        tmp_val += v2_rho_b_gamma_bb[P] * gamma_bbk[P];
                        C_DAXPY(nlocal, (0.5 * w[P] * tmp_val), phi[P], 1, Tbp[P], 1);
                    }

                    // => Alpha W terms <= //
                    if ((rho_a[P] < v2_rho_cutoff_) || (rho_b[P] < v2_rho_cutoff_)) continue;

                    // rho_ak
                    v2_val_aa = v2_rho_a_gamma_aa[P] * rho_ak[P];
                    v2_val_ab = v2_rho_a_gamma_ab[P] * rho_ak[P];

                    // rho_bk
                    v2_val_aa += v2_rho_b_gamma_aa[P] * rho_bk[P];
                    v2_val_ab += v2_rho_b_gamma_ab[P] * rho_bk[P];

                    // gamma_aak
                    v2_val_aa += v2_gamma_aa_gamma_aa[P] * gamma_aak[P];
                    v2_val_ab += v2_gamma_aa_gamma_ab[P] * gamma_aak[P];

                    // gamma_abk
                    v2_val_aa += v2_gamma_aa_gamma_ab[P] * gamma_abk[P];
                    v2_val_ab += v2_gamma_ab_gamma_ab[P] * gamma_abk[P];

                    // gamma_bbk
                    v2_val_aa += v2_gamma_aa_gamma_bb[P] * gamma_bbk[P];
                    v2_val_ab += v2_gamma_ab_gamma_bb[P] * gamma_bbk[P];

                    // Wx
                    tmp_val = 2.0 * v_gamma_aa[P] * rho_ak_x[P];
                    tmp_val += v_gamma_ab[P] * rho_bk_x[P];
                    tmp_val += 2.0 * v2_val_aa * rho_ax[P];
                    tmp_val += v2_val_ab * rho_bx[P];
                    tmp_val *= w[P];

                    C_DAXPY(nlocal, tmp_val, phi_x[P], 1, Tap[P], 1);

                    // Wy
                    tmp_val = 2.0 * v_gamma_aa[P] * rho_ak_y[P];
                    tmp_val += v_gamma_ab[P] * rho_bk_y[P];
                    tmp_val += 2.0 * v2_val_aa * rho_ay[P];
                    tmp_val += v2_val_ab * rho_by[P];
                    tmp_val *= w[P];

                    C_DAXPY(nlocal, tmp_val, phi_y[P], 1, Tap[P], 1);

                    // Wz
                    tmp_val = 2.0 * v_gamma_aa[P] * rho_ak_z[P];
                    tmp_val += v_gamma_ab[P] * rho_bk_z[P];
                    tmp_val += 2.0 * v2_val_aa * rho_az[P];
                    tmp_val += v2_val_ab * rho_bz[P];
                    tmp_val *= w[P];

                    C_DAXPY(nlocal, tmp_val, phi_z[P], 1, Tap[P], 1);

                    // => Beta W terms <= //

                    // rho_ak
                    v2_val_bb = v2_rho_a_gamma_bb[P] * rho_ak[P];
                    v2_val_ab = v2_rho_a_gamma_ab[P] * rho_ak[P];

                    // rho_bk
                    v2_val_bb += v2_rho_b_gamma_bb[P] * rho_bk[P];
                    v2_val_ab += v2_rho_b_gamma_ab[P] * rho_bk[P];

                    // gamma_bbk
                    v2_val_bb += v2_gamma_bb_gamma_bb[P] * gamma_bbk[P];
                    v2_val_ab += v2_gamma_ab_gamma_bb[P] * gamma_bbk[P];

                    // gamma_abk
                    v2_val_bb += v2_gamma_ab_gamma_bb[P] * gamma_abk[P];
                    v2_val_ab += v2_gamma_ab_gamma_ab[P] * gamma_abk[P];

                    // gamma_aak
                    v2_val_bb += v2_gamma_aa_gamma_bb[P] * gamma_aak[P];
                    v2_val_ab += v2_gamma_aa_gamma_ab[P] * gamma_aak[P];

                    // Wx
                    tmp_val = 2.0 * v_gamma_bb[P] * rho_bk_x[P];
                    tmp_val += v_gamma_ab[P] * rho_ak_x[P];
                    tmp_val += 2.0 * v2_val_bb * rho_bx[P];
                    tmp_val += v2_val_ab * rho_ax[P];
                    tmp_val *= w[P];

                    C_DAXPY(nlocal, tmp_val, phi_x[P], 1, Tbp[P], 1);

                    // Wy
                    tmp_val = 2.0 * v_gamma_bb[P] * rho_bk_y[P];
                    tmp_val += v_gamma_ab[P] * rho_ak_y[P];
                    tmp_val += 2.0 * v2_val_bb * rho_by[P];
                    tmp_val += v2_val_ab * rho_ay[P];
                    tmp_val *= w[P];

                    C_DAXPY(nlocal, tmp_val, phi_y[P], 1, Tbp[P], 1);

                    // Wz
                    tmp_val = 2.0 * v_gamma_bb[P] * rho_bk_z[P];
                    tmp_val += v_gamma_ab[P] * rho_ak_z[P];
                    tmp_val += 2.0 * v2_val_bb * rho_bz[P];
                    tmp_val += v2_val_ab * rho_az[P];
                    tmp_val *= w[P];

                    C_DAXPY(nlocal, tmp_val, phi_z[P], 1, Tbp[P], 1);
                }
            }

            // Put it all together
            C_DGEMM('T', 'N', nlocal, nlocal, npoints, 1.0, phi[0], coll_funcs, Tap[0], max_functions, 0.0,
                    Vax_localp[0], max_functions);
            C_DGEMM('T', 'N', nlocal, nlocal, npoints, 1.0, phi[0], coll_funcs, Tbp[0], max_functions, 0.0,
                    Vbx_localp[0], max_functions);

            // Symmetrization (V is *always* Hermitian)
            for (int m = 0; m < nlocal; m++) {
                for (int n = 0; n <= m; n++) {
                    Vax_localp[m][n] = Vax_localp[n][m] = Vax_localp[m][n] + Vax_localp[n][m];
                    Vbx_localp[m][n] = Vbx_localp[n][m] = Vbx_localp[m][n] + Vbx_localp[n][m];
                }
            }
            // R_Vax_local[rank]->print();
            // R_Vbx_local[rank]->print();

            // => Unpacking <= //
            double** Vaxp = Vax_AO[dindex]->pointer();
            double** Vbxp = Vbx_AO[dindex]->pointer();
            for (int ml = 0; ml < nlocal; ml++) {
                int mg = function_map[ml];
                for (int nl = 0; nl < ml; nl++) {
                    int ng = function_map[nl];

#pragma omp atomic update
                    Vaxp[mg][ng] += Vax_localp[ml][nl];
#pragma omp atomic update
                    Vaxp[ng][mg] += Vax_localp[ml][nl];

#pragma omp atomic update
                    Vbxp[mg][ng] += Vbx_localp[ml][nl];
#pragma omp atomic update
                    Vbxp[ng][mg] += Vbx_localp[ml][nl];
                }
#pragma omp atomic update
                Vaxp[mg][mg] += Vax_localp[ml][ml];
#pragma omp atomic update
                Vbxp[mg][mg] += Vbx_localp[ml][ml];
            }
            parallel_timer_off("V_XCd", rank);
        }
    }

    // Set the result
    for (size_t i = 0; i < (Dx.size() / 2); i++) {
        if (Dx[i]->nirrep() != 1) {
            ret[2 * i]->apply_symmetry(Vax_AO[i], AO2USO_);
            ret[2 * i + 1]->apply_symmetry(Vbx_AO[i], AO2USO_);
        } else {
            ret[2 * i]->copy(Vax_AO[i]);
            ret[2 * i + 1]->copy(Vbx_AO[i]);
        }
    }

    // Reset the workers
    for (size_t i = 0; i < num_threads_; i++) {
        functional_workers_[i]->set_deriv(old_func_deriv);
        functional_workers_[i]->allocate();
    }

    timer_off("UV: Form Vx");
}
SharedMatrix UV::compute_gradient() {
    if ((D_AO_.size() != 2)) throw PSIEXCEPTION("V: UKS should have two D Matrices");

    if (functional_->needs_vv10()) {
        throw PSIEXCEPTION("V: UKS cannot compute VV10 gradient contribution.");
    }

    int rank = 0;

    // Build the target gradient Matrix
    int natom = primary_->molecule()->natom();
    auto G = std::make_shared<Matrix>("XC Gradient", natom, 3);
    double** Gp = G->pointer();

    // What local XC ansatz are we in?
    int ansatz = functional_->ansatz();

    // How many functions are there (for lda in Vtemp, T)
    int max_functions = grid_->max_functions();
    int max_points = grid_->max_points();

    // Set Hessian derivative level in properties
    int old_deriv = point_workers_[0]->deriv();

    // Setup the pointers
    for (size_t i = 0; i < num_threads_; i++) {
        point_workers_[i]->set_pointers(D_AO_[0], D_AO_[1]);
        point_workers_[i]->set_deriv((functional_->is_gga() || functional_->is_meta() ? 2 : 1));
    }

    // Thread scratch
    std::vector<std::shared_ptr<Vector>> Q_temp;
    for (size_t i = 0; i < num_threads_; i++) {
        Q_temp.push_back(std::make_shared<Vector>("Quadrature Temp", max_points));
    }

    std::vector<double> functionalq(num_threads_);
    std::vector<double> rhoaq(num_threads_);
    std::vector<double> rhoaxq(num_threads_);
    std::vector<double> rhoayq(num_threads_);
    std::vector<double> rhoazq(num_threads_);
    std::vector<double> rhobq(num_threads_);
    std::vector<double> rhobxq(num_threads_);
    std::vector<double> rhobyq(num_threads_);
    std::vector<double> rhobzq(num_threads_);

    // timer_off("V: V_XC");
    for (size_t Q = 0; Q < grid_->blocks().size(); Q++) {
// Get thread info
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif

        std::shared_ptr<SuperFunctional> fworker = functional_workers_[rank];
        std::shared_ptr<PointFunctions> pworker = point_workers_[rank];
        double* QTp = Q_temp[rank]->pointer();

        double** Tap = pworker->scratch()[0]->pointer();
        double** Tbp = pworker->scratch()[1]->pointer();
        double** Dap = pworker->D_scratch()[0]->pointer();
        double** Dbp = pworker->D_scratch()[1]->pointer();

        SharedMatrix Ua_local(pworker->scratch()[0]->clone());
        double** Uap = Ua_local->pointer();
        SharedMatrix Ub_local(pworker->scratch()[1]->clone());
        double** Ubp = Ub_local->pointer();

        // Grid info
        std::shared_ptr<BlockOPoints> block = grid_->blocks()[Q];
        int npoints = block->npoints();
        double* x = block->x();
        double* y = block->y();
        double* z = block->z();
        double* w = block->w();
        const std::vector<int>& function_map = block->functions_local_to_global();
        int nlocal = function_map.size();

        // Compute grid and functional
        parallel_timer_on("Properties", rank);
        pworker->compute_points(block);
        parallel_timer_off("Properties", rank);

        parallel_timer_on("Functional", rank);
        std::map<std::string, SharedVector>& vals = fworker->compute_functional(pworker->point_values(), npoints);
        parallel_timer_off("Functional", rank);

        // More pointers
        parallel_timer_on("V_xc gradient", rank);
        double** phi = pworker->basis_value("PHI")->pointer();
        double** phi_x = pworker->basis_value("PHI_X")->pointer();
        double** phi_y = pworker->basis_value("PHI_Y")->pointer();
        double** phi_z = pworker->basis_value("PHI_Z")->pointer();
        double* rho_a = pworker->point_value("RHO_A")->pointer();
        double* rho_b = pworker->point_value("RHO_B")->pointer();
        double* zk = vals["V"]->pointer();
        double* v_rho_a = vals["V_RHO_A"]->pointer();
        double* v_rho_b = vals["V_RHO_B"]->pointer();
        size_t coll_funcs = pworker->basis_value("PHI")->ncol();

        // => Quadrature values <= //
        functionalq[rank] += C_DDOT(npoints, w, 1, zk, 1);
        for (int P = 0; P < npoints; P++) {
            QTp[P] = w[P] * rho_a[P];
        }
        rhoaq[rank] += C_DDOT(npoints, w, 1, rho_a, 1);
        rhoaxq[rank] += C_DDOT(npoints, QTp, 1, x, 1);
        rhoayq[rank] += C_DDOT(npoints, QTp, 1, y, 1);
        rhoazq[rank] += C_DDOT(npoints, QTp, 1, z, 1);
        for (int P = 0; P < npoints; P++) {
            QTp[P] = w[P] * rho_b[P];
        }
        rhobq[rank] += C_DDOT(npoints, w, 1, rho_b, 1);
        rhobxq[rank] += C_DDOT(npoints, QTp, 1, x, 1);
        rhobyq[rank] += C_DDOT(npoints, QTp, 1, y, 1);
        rhobzq[rank] += C_DDOT(npoints, QTp, 1, z, 1);

        // => LSDA Contribution <= //
        for (int P = 0; P < npoints; P++) {
            std::fill(Tap[P], Tap[P] + nlocal, 0.0);
            std::fill(Tbp[P], Tbp[P] + nlocal, 0.0);
            C_DAXPY(nlocal, -2.0 * w[P] * v_rho_a[P], phi[P], 1, Tap[P], 1);
            C_DAXPY(nlocal, -2.0 * w[P] * v_rho_b[P], phi[P], 1, Tbp[P], 1);
        }

        // => GGA Contribution (Term 1) <= //
        if (fworker->is_gga()) {
            double* rho_ax = pworker->point_value("RHO_AX")->pointer();
            double* rho_ay = pworker->point_value("RHO_AY")->pointer();
            double* rho_az = pworker->point_value("RHO_AZ")->pointer();
            double* rho_bx = pworker->point_value("RHO_BX")->pointer();
            double* rho_by = pworker->point_value("RHO_BY")->pointer();
            double* rho_bz = pworker->point_value("RHO_BZ")->pointer();
            double* v_gamma_aa = vals["V_GAMMA_AA"]->pointer();
            double* v_gamma_ab = vals["V_GAMMA_AB"]->pointer();
            double* v_gamma_bb = vals["V_GAMMA_BB"]->pointer();

            for (int P = 0; P < npoints; P++) {
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_aa[P] * rho_ax[P] + v_gamma_ab[P] * rho_bx[P]), phi_x[P],
                        1, Tap[P], 1);
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_aa[P] * rho_ay[P] + v_gamma_ab[P] * rho_by[P]), phi_y[P],
                        1, Tap[P], 1);
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_aa[P] * rho_az[P] + v_gamma_ab[P] * rho_bz[P]), phi_z[P],
                        1, Tap[P], 1);
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_bb[P] * rho_bx[P] + v_gamma_ab[P] * rho_ax[P]), phi_x[P],
                        1, Tbp[P], 1);
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_bb[P] * rho_by[P] + v_gamma_ab[P] * rho_ay[P]), phi_y[P],
                        1, Tbp[P], 1);
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_bb[P] * rho_bz[P] + v_gamma_ab[P] * rho_az[P]), phi_z[P],
                        1, Tbp[P], 1);
            }
        }

        // => Synthesis <= //
        C_DGEMM('N', 'N', npoints, nlocal, nlocal, 1.0, Tap[0], max_functions, Dap[0], max_functions, 0.0, Uap[0],
                max_functions);
        C_DGEMM('N', 'N', npoints, nlocal, nlocal, 1.0, Tbp[0], max_functions, Dbp[0], max_functions, 0.0, Ubp[0],
                max_functions);

        for (int ml = 0; ml < nlocal; ml++) {
            int A = primary_->function_to_center(function_map[ml]);
            Gp[A][0] += C_DDOT(npoints, &Uap[0][ml], max_functions, &phi_x[0][ml], coll_funcs);
            Gp[A][1] += C_DDOT(npoints, &Uap[0][ml], max_functions, &phi_y[0][ml], coll_funcs);
            Gp[A][2] += C_DDOT(npoints, &Uap[0][ml], max_functions, &phi_z[0][ml], coll_funcs);
            Gp[A][0] += C_DDOT(npoints, &Ubp[0][ml], max_functions, &phi_x[0][ml], coll_funcs);
            Gp[A][1] += C_DDOT(npoints, &Ubp[0][ml], max_functions, &phi_y[0][ml], coll_funcs);
            Gp[A][2] += C_DDOT(npoints, &Ubp[0][ml], max_functions, &phi_z[0][ml], coll_funcs);
        }

        // => GGA Contribution (Term 2) <= //
        if (fworker->is_gga()) {
            double** phi_xx = pworker->basis_value("PHI_XX")->pointer();
            double** phi_xy = pworker->basis_value("PHI_XY")->pointer();
            double** phi_xz = pworker->basis_value("PHI_XZ")->pointer();
            double** phi_yy = pworker->basis_value("PHI_YY")->pointer();
            double** phi_yz = pworker->basis_value("PHI_YZ")->pointer();
            double** phi_zz = pworker->basis_value("PHI_ZZ")->pointer();
            double* rho_ax = pworker->point_value("RHO_AX")->pointer();
            double* rho_ay = pworker->point_value("RHO_AY")->pointer();
            double* rho_az = pworker->point_value("RHO_AZ")->pointer();
            double* rho_bx = pworker->point_value("RHO_BX")->pointer();
            double* rho_by = pworker->point_value("RHO_BY")->pointer();
            double* rho_bz = pworker->point_value("RHO_BZ")->pointer();
            double* v_gamma_aa = vals["V_GAMMA_AA"]->pointer();
            double* v_gamma_ab = vals["V_GAMMA_AB"]->pointer();
            double* v_gamma_bb = vals["V_GAMMA_BB"]->pointer();

            C_DGEMM('N', 'N', npoints, nlocal, nlocal, 1.0, phi[0], coll_funcs, Dap[0], max_functions, 0.0, Uap[0],
                    max_functions);
            C_DGEMM('N', 'N', npoints, nlocal, nlocal, 1.0, phi[0], coll_funcs, Dbp[0], max_functions, 0.0, Ubp[0],
                    max_functions);

            // x
            for (int P = 0; P < npoints; P++) {
                std::fill(Tap[P], Tap[P] + nlocal, 0.0);
                std::fill(Tbp[P], Tbp[P] + nlocal, 0.0);
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_aa[P] * rho_ax[P] + v_gamma_ab[P] * rho_bx[P]), Uap[P], 1,
                        Tap[P], 1);
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_bb[P] * rho_bx[P] + v_gamma_ab[P] * rho_ax[P]), Ubp[P], 1,
                        Tbp[P], 1);
            }
            for (int ml = 0; ml < nlocal; ml++) {
                int A = primary_->function_to_center(function_map[ml]);
                Gp[A][0] += C_DDOT(npoints, &Tap[0][ml], max_functions, &phi_xx[0][ml], coll_funcs);
                Gp[A][1] += C_DDOT(npoints, &Tap[0][ml], max_functions, &phi_xy[0][ml], coll_funcs);
                Gp[A][2] += C_DDOT(npoints, &Tap[0][ml], max_functions, &phi_xz[0][ml], coll_funcs);
                Gp[A][0] += C_DDOT(npoints, &Tbp[0][ml], max_functions, &phi_xx[0][ml], coll_funcs);
                Gp[A][1] += C_DDOT(npoints, &Tbp[0][ml], max_functions, &phi_xy[0][ml], coll_funcs);
                Gp[A][2] += C_DDOT(npoints, &Tbp[0][ml], max_functions, &phi_xz[0][ml], coll_funcs);
            }

            // y
            for (int P = 0; P < npoints; P++) {
                std::fill(Tap[P], Tap[P] + nlocal, 0.0);
                std::fill(Tbp[P], Tbp[P] + nlocal, 0.0);
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_aa[P] * rho_ay[P] + v_gamma_ab[P] * rho_by[P]), Uap[P], 1,
                        Tap[P], 1);
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_bb[P] * rho_by[P] + v_gamma_ab[P] * rho_ay[P]), Ubp[P], 1,
                        Tbp[P], 1);
            }
            for (int ml = 0; ml < nlocal; ml++) {
                int A = primary_->function_to_center(function_map[ml]);
                Gp[A][0] += C_DDOT(npoints, &Tap[0][ml], max_functions, &phi_xy[0][ml], coll_funcs);
                Gp[A][1] += C_DDOT(npoints, &Tap[0][ml], max_functions, &phi_yy[0][ml], coll_funcs);
                Gp[A][2] += C_DDOT(npoints, &Tap[0][ml], max_functions, &phi_yz[0][ml], coll_funcs);
                Gp[A][0] += C_DDOT(npoints, &Tbp[0][ml], max_functions, &phi_xy[0][ml], coll_funcs);
                Gp[A][1] += C_DDOT(npoints, &Tbp[0][ml], max_functions, &phi_yy[0][ml], coll_funcs);
                Gp[A][2] += C_DDOT(npoints, &Tbp[0][ml], max_functions, &phi_yz[0][ml], coll_funcs);
            }

            // z
            for (int P = 0; P < npoints; P++) {
                std::fill(Tap[P], Tap[P] + nlocal, 0.0);
                std::fill(Tbp[P], Tbp[P] + nlocal, 0.0);
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_aa[P] * rho_az[P] + v_gamma_ab[P] * rho_bz[P]), Uap[P], 1,
                        Tap[P], 1);
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_bb[P] * rho_bz[P] + v_gamma_ab[P] * rho_az[P]), Ubp[P], 1,
                        Tbp[P], 1);
            }
            for (int ml = 0; ml < nlocal; ml++) {
                int A = primary_->function_to_center(function_map[ml]);
                Gp[A][0] += C_DDOT(npoints, &Tap[0][ml], max_functions, &phi_xz[0][ml], coll_funcs);
                Gp[A][1] += C_DDOT(npoints, &Tap[0][ml], max_functions, &phi_yz[0][ml], coll_funcs);
                Gp[A][2] += C_DDOT(npoints, &Tap[0][ml], max_functions, &phi_zz[0][ml], coll_funcs);
                Gp[A][0] += C_DDOT(npoints, &Tbp[0][ml], max_functions, &phi_xz[0][ml], coll_funcs);
                Gp[A][1] += C_DDOT(npoints, &Tbp[0][ml], max_functions, &phi_yz[0][ml], coll_funcs);
                Gp[A][2] += C_DDOT(npoints, &Tbp[0][ml], max_functions, &phi_zz[0][ml], coll_funcs);
            }
        }

        // => Meta Contribution <= //
        if (fworker->is_meta()) {
            double** phi_xx = pworker->basis_value("PHI_XX")->pointer();
            double** phi_xy = pworker->basis_value("PHI_XY")->pointer();
            double** phi_xz = pworker->basis_value("PHI_XZ")->pointer();
            double** phi_yy = pworker->basis_value("PHI_YY")->pointer();
            double** phi_yz = pworker->basis_value("PHI_YZ")->pointer();
            double** phi_zz = pworker->basis_value("PHI_ZZ")->pointer();
            double* v_tau_a = vals["V_TAU_A"]->pointer();
            double* v_tau_b = vals["V_TAU_B"]->pointer();

            double** phi_i[3];
            phi_i[0] = phi_x;
            phi_i[1] = phi_y;
            phi_i[2] = phi_z;

            double** phi_ij[3][3];
            phi_ij[0][0] = phi_xx;
            phi_ij[0][1] = phi_xy;
            phi_ij[0][2] = phi_xz;
            phi_ij[1][0] = phi_xy;
            phi_ij[1][1] = phi_yy;
            phi_ij[1][2] = phi_yz;
            phi_ij[2][0] = phi_xz;
            phi_ij[2][1] = phi_yz;
            phi_ij[2][2] = phi_zz;

            double** Ds[2];
            Ds[0] = Dap;
            Ds[1] = Dbp;

            double* v_tau_s[2];
            v_tau_s[0] = v_tau_a;
            v_tau_s[1] = v_tau_b;

            for (int s = 0; s < 2; s++) {
                double* v_tau = v_tau_s[s];
                for (int i = 0; i < 3; i++) {
                    double*** phi_j = phi_ij[i];
                    C_DGEMM('N', 'N', npoints, nlocal, nlocal, 1.0, phi_i[i][0], coll_funcs, Ds[s][0], max_functions,
                            0.0, Uap[0], max_functions);
                    for (int P = 0; P < npoints; P++) {
                        std::fill(Tap[P], Tap[P] + nlocal, 0.0);
                        C_DAXPY(nlocal, -2.0 * w[P] * (v_tau[P]), Uap[P], 1, Tap[P], 1);
                    }
                    for (int ml = 0; ml < nlocal; ml++) {
                        int A = primary_->function_to_center(function_map[ml]);
                        Gp[A][0] += C_DDOT(npoints, &Tap[0][ml], max_functions, &phi_j[0][0][ml], coll_funcs);
                        Gp[A][1] += C_DDOT(npoints, &Tap[0][ml], max_functions, &phi_j[1][0][ml], coll_funcs);
                        Gp[A][2] += C_DDOT(npoints, &Tap[0][ml], max_functions, &phi_j[2][0][ml], coll_funcs);
                    }
                }
            }
        }
        Ua_local.reset();
        Ub_local.reset();
        parallel_timer_off("V_xc gradient", rank);
    }
    // timer_off("V: V_XC");

    quad_values_["FUNCTIONAL"] = std::accumulate(functionalq.begin(), functionalq.end(), 0.0);
    quad_values_["RHO_A"] = std::accumulate(rhoaq.begin(), rhoaq.end(), 0.0);
    quad_values_["RHO_AX"] = std::accumulate(rhoaxq.begin(), rhoaxq.end(), 0.0);
    quad_values_["RHO_AY"] = std::accumulate(rhoayq.begin(), rhoayq.end(), 0.0);
    quad_values_["RHO_AZ"] = std::accumulate(rhoazq.begin(), rhoazq.end(), 0.0);
    quad_values_["RHO_B"] = std::accumulate(rhobq.begin(), rhobq.end(), 0.0);
    quad_values_["RHO_BX"] = std::accumulate(rhobxq.begin(), rhobxq.end(), 0.0);
    quad_values_["RHO_BY"] = std::accumulate(rhobyq.begin(), rhobyq.end(), 0.0);
    quad_values_["RHO_BZ"] = std::accumulate(rhobzq.begin(), rhobzq.end(), 0.0);

    if (std::isnan(quad_values_["FUNCTIONAL"])) {
        throw PSIEXCEPTION("V: Integrated DFT functional to get NaN. The functional is not numerically stable. Pick a different one.");
    }

    if (debug_) {
        outfile->Printf("   => XC Gradient: Numerical Integrals <=\n\n");
        outfile->Printf("    Functional Value:  %24.16E\n", quad_values_["FUNCTIONAL"]);
        outfile->Printf("    <\\rho_a>        :  %24.16E\n", quad_values_["RHO_A"]);
        outfile->Printf("    <\\rho_b>        :  %24.16E\n", quad_values_["RHO_B"]);
        outfile->Printf("    <\\vec r\\rho_a>  : <%24.16E,%24.16E,%24.16E>\n", quad_values_["RHO_AX"],
                        quad_values_["RHO_AY"], quad_values_["RHO_AZ"]);
        outfile->Printf("    <\\vec r\\rho_b>  : <%24.16E,%24.16E,%24.16E>\n\n", quad_values_["RHO_BX"],
                        quad_values_["RHO_BY"], quad_values_["RHO_BZ"]);
    }

    for (size_t i = 0; i < num_threads_; i++) {
        point_workers_[i]->set_deriv(old_deriv);
    }

    return G;
}
}  // namespace psi
