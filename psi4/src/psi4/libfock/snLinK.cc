/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
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

#include "jk.h"
#include "SplitJK.h"
#include "psi4/libqt/qt.h"
#include "psi4/libfock/cubature.h"
#include "psi4/libfock/points.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/electrostatic.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/integral.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/lib3index/dftensor.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include <unordered_set>
#include <vector>
#include <map>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace psi;

namespace psi {

// converts a Psi4::Molecule object to a GauXC::Molecule object
GauXC::Molecule snLinK::psi4_to_gauxc_molecule(std::shared_ptr<Molecule> psi4_molecule) {
    GauXC::Molecule gauxc_molecule;
 
    // TODO: Check if Bohr/Angstrom conversion is needed
    for (size_t iatom = 0; iatom != psi4_molecule->natom(); ++iatom) {
        gauxc_molecule.emplace_back(
            GauXC::AtomicNumber(psi4_molecule->true_atomic_number(iatom)), 
            psi4_molecule->x(iatom),
            psi4_molecule->y(iatom),
            psi4_molecule->z(iatom)
        );
    }

    return gauxc_molecule;
}

// converts a Psi4::BasisSet object to a GauXC::BasisSet object
template <typename T>
GauXC::BasisSet<T> snLinK::psi4_to_gauxc_basisset(std::shared_ptr<BasisSet> psi4_basisset) {
    using prim_array = typename GauXC::Shell<T>::prim_array;
    using cart_array = typename GauXC::Shell<T>::cart_array;

    GauXC::BasisSet<T> gauxc_basisset;
  
    for (size_t ishell = 0; ishell != psi4_basisset->nshell(); ++ishell) {
        auto psi4_shell = psi4_basisset->shell(ishell);
       
        const auto nprim = GauXC::PrimSize(psi4_shell.nprimitive());
        prim_array alpha; 
        prim_array coeff;

        //TODO: Ensure normalization is okay
        for (size_t iprim = 0; iprim != psi4_shell.nprimitive(); ++iprim) {
            alpha.at(iprim) = psi4_shell.exp(iprim);
            coeff.at(iprim) = psi4_shell.coef(iprim);
        }

        auto psi4_shell_center = psi4_shell.center();
        cart_array center = { psi4_shell_center[0], psi4_shell_center[1], psi4_shell_center[2] };

        gauxc_basisset.emplace_back(
            nprim, 
            GauXC::AngularMomentum(psi4_shell.am()), //TODO: check if am() is 0-indexed
            GauXC::SphericalType(!(psi4_shell.is_cartesian())),
            alpha,
            coeff,
            center
        );
    }
    
    for (auto& sh : gauxc_basisset) {
        sh.set_shell_tolerance(std::numeric_limits<double>::epsilon()); 
    }

    return gauxc_basisset;
}

// create GauXC grid using Psi4 option parameters
// TODO: This! 
//GauXC::MolGrid snLinK::psi4_to_gauxc_grid() {
void snLinK::psi4_to_gauxc_grid() {
    using atomic_grid_map = std::unordered_map< GauXC::AtomicNumber, GauXC::Grid >;

    // create Psi4 grid  
    std::map<std::string, std::string> grid_str_options = {
        {"DFT_PRUNING_SCHEME", options_.get_str("COSX_PRUNING_SCHEME")},
        {"DFT_RADIAL_SCHEME",  "TREUTLER"},
        {"DFT_NUCLEAR_SCHEME", "TREUTLER"},
        {"DFT_GRID_NAME",      ""},
        {"DFT_BLOCK_SCHEME",   "OCTREE"},
    };

    std::map<std::string, int> grid_int_options = {
        {"DFT_SPHERICAL_POINTS", options_.get_int("COSX_SPHERICAL_POINTS_FINAL")},
        {"DFT_RADIAL_POINTS",    options_.get_int("COSX_RADIAL_POINTS_FINAL")},
        {"DFT_BLOCK_MIN_POINTS", 100},
        {"DFT_BLOCK_MAX_POINTS", 256},
    };

    std::map<std::string, double> grid_float_options = {
        {"DFT_BASIS_TOLERANCE",   options_.get_double("COSX_BASIS_TOLERANCE")},
        {"DFT_BS_RADIUS_ALPHA",   1.0},
        {"DFT_PRUNING_ALPHA",     1.0},
        {"DFT_BLOCK_MAX_RADIUS",  3.0},
        {"DFT_WEIGHTS_TOLERANCE", 1e-15},
    };

    psi4_grid_ = std::make_shared<DFTGrid>(primary_->molecule(), primary_, grid_int_options, grid_str_options, grid_float_options, options_);

    // interpret Psi4 options for GauXC
    atomic_grid_map atom_grids; 
    
    //std::unique_ptr<GauXC::MolGrid> gauxc_molgrid;
    return;
}

// converts a Psi4::Matrix to an Eigen::MatrixXd
Eigen::MatrixXd snLinK::psi4_to_eigen_matrix(SharedMatrix psi4_matrix) {
    // Sanity checks
    if (psi4_matrix->nirrep() != 1) {
        throw PSIEXCEPTION("Psi4::Matrix must be in C1 symmetry to be transformed into Eigen::MatrixXd!");
    }

    Eigen::MatrixXd eigen_matrix(psi4_matrix->nrow(), psi4_matrix->ncol());

    // equate matrices
    auto psi4_matrix_p = psi4_matrix->pointer(); 
    for (size_t irow = 0; irow != psi4_matrix->nrow(); ++irow) {
        for (size_t icol = 0; icol != psi4_matrix->ncol(); ++icol) {
            eigen_matrix(irow, icol) = psi4_matrix_p[irow][icol]; //TODO: check that Psi4::Matrix is column-major
        }
    }

    return eigen_matrix;
}

// converts an Eigen MatrixXd object to a Psi4::Matrix 
void snLinK::eigen_to_psi4_matrix(SharedMatrix psi4_matrix, const Eigen::MatrixXd& eigen_matrix) {
    // Sanity checks
    if (psi4_matrix->nirrep() != 1) {
        throw PSIEXCEPTION("Psi4::Matrix must be in C1 symmetry to be transformed into Eigen::MatrixXd!");
    }

    if (psi4_matrix->nrow() != eigen_matrix.rows()) {
        throw PSIEXCEPTION("Input matrix row counts don't match in snLinK::eigen_to_psi4_matrix!");
    } else if (psi4_matrix->ncol() != eigen_matrix.cols()) {
        throw PSIEXCEPTION("Input matrix column counts don't match in snLinK::eigen_to_psi4_matrix!");
    }
    
    // equate matrices
    auto psi4_matrix_p = psi4_matrix->pointer(); 
    for (size_t irow = 0; irow != psi4_matrix->nrow(); ++irow) {
        for (size_t icol = 0; icol != psi4_matrix->ncol(); ++icol) {
            psi4_matrix_p[irow][icol] = eigen_matrix(irow, icol); //TODO: check that Psi4::Matrix is column-major
        }
    }
}

snLinK::snLinK(std::shared_ptr<BasisSet> primary, Options& options) : SplitJK(primary, options) {

    // always executing on host (i.e., CPU) for now
    ex_ = std::make_unique<GauXC::ExecutionSpace>(GauXC::ExecutionSpace::Host); 
    rt_ = std::make_unique<GauXC::RuntimeEnvironment>(GAUXC_MPI_CODE(MPI_COMM_WORLD));

    // convert Psi4 fundamental quantities to GauXC 
    gauxc_mol_ = psi4_to_gauxc_molecule(primary_->molecule());
    gauxc_primary_ = psi4_to_gauxc_basisset<double>(primary_);
    
    // create snLinK grid for GauXC
    //gauxc_grid_ = psi4_to_gauxc_grid();  
    gauxc_grid_ = std::make_unique<GauXC::MolGrid>(
        GauXC::MolGridFactory::create_default_molgrid(
            gauxc_mol_, 
            GauXC::PruningScheme::Robust,
            GauXC::BatchSize(512), 
            GauXC::RadialQuad::MuraKnowles, 
            GauXC::AtomicGridSizeDefault::UltraFineGrid
        )
    );
  
    // construct load balancer
    const std::string load_balancer_kernel = "Default";
    const size_t quad_pad_value = 1;

    gauxc_load_balancer_factory_ = std::make_unique<GauXC::LoadBalancerFactory>(*ex_, load_balancer_kernel);
    GauXC::LoadBalancer gauxc_load_balancer = gauxc_load_balancer_factory_->get_instance(*rt_, gauxc_mol_, *gauxc_grid_, gauxc_primary_, quad_pad_value);
    
    // construct weights module
    const std::string mol_weights_kernel = "Default";
    gauxc_mol_weights_factory_ = std::make_unique<GauXC::MolecularWeightsFactory>(*ex_, mol_weights_kernel, GauXC::MolecularWeightsSettings{});
    GauXC::MolecularWeights gauxc_mol_weights = gauxc_mol_weights_factory_->get_instance();

    // apply partition weights
    gauxc_mol_weights.modify_weights(gauxc_load_balancer);
    
    // set integrator options 
    const std::string integrator_input_type = "Replicated";
    const std::string integrator_kernel = "Default";
    const std::string reduction_kernel = "Default";
    const std::string lwd_kernel = "Default";
    
    // construct dummy functional
    const std::string dummy_func_str = "B3LYP";
    auto spin = options_.get_str("REFERENCE") == "UHF" ? ExchCXX::Spin::Polarized : ExchCXX::Spin::Unpolarized;

    ExchCXX::Functional dummy_func_key = ExchCXX::functional_map.value(dummy_func_str);
    ExchCXX::XCFunctional dummy_func = { ExchCXX::Backend::builtin, dummy_func_key, spin };  
   
    // actually construct integrator 
    integrator_factory_ = std::make_unique<GauXC::XCIntegratorFactory<matrix_type> >(*ex_, integrator_input_type, integrator_kernel, lwd_kernel, reduction_kernel);
    integrator_ = integrator_factory_->get_shared_instance(dummy_func, gauxc_load_balancer); 
}

snLinK::~snLinK() {}

void snLinK::print_header() const {
    if (print_) {
        outfile->Printf("\n");
        outfile->Printf("  ==> snLinK: Semi-Numerical Linear Exchange K <==\n\n");

        //outfile->Printf("    K Screening Cutoff: %11.0E\n", kscreen_);
        //outfile->Printf("    K Density Cutoff:   %11.0E\n", dscreen_);
        //outfile->Printf("    K Basis Cutoff:     %11.0E\n", basis_tol_);
        //outfile->Printf("    K Overlap Fitting:  %11s\n", (overlap_fitted_ ? "Yes" : "No"));
    }
}

// build the K matrix using Neeses's Chain-of-Spheres Exchange algorithm
// algorithm is originally proposed in https://doi.org/10.1016/j.chemphys.2008.10.036
// overlap fitting is discussed in https://doi.org/10.1063/1.3646921
void snLinK::build_G_component(std::vector<std::shared_ptr<Matrix>>& D, std::vector<std::shared_ptr<Matrix>>& K,
    std::vector<std::shared_ptr<TwoBodyAOInt> >& eri_computers) {
    
    for (int iD = 0; iD != D.size(); ++iD) {
        // compute K for density Di using GauXC
        auto D_eigen = psi4_to_eigen_matrix(D[iD]);    
        auto K_eigen = integrator_->eval_exx(D_eigen);

        // convert result to Psi4 matrix
        eigen_to_psi4_matrix(K[iD], K_eigen);
    }
    
    return;
}

}  // namespace psi
