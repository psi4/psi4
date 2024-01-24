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

#include <algorithm>
#include <map>
#include <unordered_set>
#include <variant>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace psi;

namespace psi {

// creates maps for translating Psi4 option inputs to GauXC enums
void snLinK::generate_enum_mappings() {
    // generate map for grid pruning schemes 
    pruning_scheme_map_["ROBUST"] = GauXC::PruningScheme::Robust;
    pruning_scheme_map_["TREUTLER"] = GauXC::PruningScheme::Treutler;
    pruning_scheme_map_["NONE"] = GauXC::PruningScheme::Unpruned;

    // generate map for radial quadrature schemes 
    radial_scheme_map_["TREUTLER"] = GauXC::RadialQuad::TreutlerAldrichs;
    radial_scheme_map_["MURA"] = GauXC::RadialQuad::MuraKnowles;
    // TODO: confirm the correctness of this specific mapping
    // Answer: Yep, it is! The Murray, Handy, Laming literature reference
    // is mentioned in cubature.cc
    radial_scheme_map_["EM"] = GauXC::RadialQuad::MurrayHandyLaming; 
}

// converts a Psi4::Molecule object to a GauXC::Molecule object
GauXC::Molecule snLinK::psi4_to_gauxc_molecule(std::shared_ptr<Molecule> psi4_molecule) {
    GauXC::Molecule gauxc_molecule;
 
    outfile->Printf("snLinK::psi4_to_gauxc_molecule\n");
    outfile->Printf("------------------------------\n");
    // TODO: Check if Bohr/Angstrom conversion is needed 
    // Answer: Nope! GauXC accepts input in Bohrs, and Psi4 coords are in Bohr
    // at this point. 
    for (size_t iatom = 0; iatom != psi4_molecule->natom(); ++iatom) {
        auto atomic_number = psi4_molecule->true_atomic_number(iatom);
        auto x_coord = psi4_molecule->x(iatom);
        auto y_coord = psi4_molecule->y(iatom);
        auto z_coord = psi4_molecule->z(iatom);
        outfile->Printf("  Atom #%i: %f, %f, %f\n", atomic_number, x_coord, y_coord, z_coord); 
        
        gauxc_molecule.emplace_back(GauXC::AtomicNumber(atomic_number), x_coord, y_coord, z_coord);
    }

    return std::move(gauxc_molecule);
}

// converts a Psi4::BasisSet object to a GauXC::BasisSet object
template <typename T>
GauXC::BasisSet<T> snLinK::psi4_to_gauxc_basisset(std::shared_ptr<BasisSet> psi4_basisset) {
    using prim_array = typename GauXC::Shell<T>::prim_array;
    using cart_array = typename GauXC::Shell<T>::cart_array;

    GauXC::BasisSet<T> gauxc_basisset;
 
    outfile->Printf("snLinK::psi4_to_gauxc_basisset\n");
    outfile->Printf("------------------------------\n");
    for (size_t ishell = 0; ishell != psi4_basisset->nshell(); ++ishell) {
        auto psi4_shell = psi4_basisset->shell(ishell);
       
        const auto nprim = GauXC::PrimSize(psi4_shell.nprimitive());
        prim_array alpha; 
        prim_array coeff;

        outfile->Printf("  ");
        outfile->Printf("%s", (psi4_shell.is_cartesian()) ? "Cartesian" : "Spherical");
        outfile->Printf(" Shell #%i (AM %i): %i primitives\n", ishell, psi4_shell.am(), psi4_shell.nprimitive());
        //TODO: Ensure normalization is okay
        // It seems so! We need to turn explicit normalization off for the Psi4-to-GauXC interface, since Psi4 
        // basis set object contains normalized coefficients already
        for (size_t iprim = 0; iprim != psi4_shell.nprimitive(); ++iprim) {
            alpha.at(iprim) = psi4_shell.exp(iprim);
            coeff.at(iprim) = psi4_shell.coef(iprim);
            outfile->Printf("    Primitive #%i: %f, %f\n", iprim, psi4_shell.exp(iprim), psi4_shell.coef(iprim));
        }

        auto psi4_shell_center = psi4_shell.center();
        cart_array center = { psi4_shell_center[0], psi4_shell_center[1], psi4_shell_center[2] };

        gauxc_basisset.emplace_back(
            nprim,
            // TODO: check if am() is 0-indexed
            // Answer: It is! See tests/basisset_test.cxx in GauXC
            GauXC::AngularMomentum(psi4_shell.am()), 
            GauXC::SphericalType(!(psi4_shell.is_cartesian())),
            alpha,
            coeff,
            center,
            false // do not normalize shell via GauXC; it is normalized via Psi4
        );
    }
    
    for (auto& sh : gauxc_basisset) {
        sh.set_shell_tolerance(basis_tol_); 
    }

    return std::move(gauxc_basisset);
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
   
    // set Psi4-specific parameters
    radial_points_ = options_.get_int("SNLINK_RADIAL_POINTS");
    spherical_points_ = options_.get_int("SNLINK_SPHERICAL_POINTS");
    basis_tol_ = options.get_double("SNLINK_BASIS_TOLERANCE");

    generate_enum_mappings();
    
    pruning_scheme_ = options_.get_str("SNLINK_PRUNING_SCHEME");
    
    // define and sanity-check radial quadrature scheme
    // sanity-checking needed because generalist DFT_RADIAL_SCHEME is used here
    radial_scheme_ = options_.get_str("DFT_RADIAL_SCHEME");
    
    std::array<std::string, 3> valid_radial_schemes = { "TREUTLER", "MURA", "EM" };
    bool is_valid_radial_scheme = std::any_of(
      valid_radial_schemes.cbegin(),
      valid_radial_schemes.cend(),
      [&](std::string valid_radial_scheme) { return radial_scheme_.find(valid_radial_scheme) != std::string::npos; }
    );

    if (!is_valid_radial_scheme) {
        throw PSIEXCEPTION("Invalid radial quadrature scheme selected for snLinK! Please set DFT_RADIAL_SCHEME to either TREUTLER, MURA, or EM.");
    }

    // define execution space
    use_gpu_ = options_.get_bool("SNLINK_USE_GPU"); 
    auto ex = use_gpu_ ? GauXC::ExecutionSpace::Device : GauXC::ExecutionSpace::Host;  
    // only running on host (i.e., GPU disabled) for the moment
    if (ex == GauXC::ExecutionSpace::Device) {
        throw PSIEXCEPTION("GauXC GPU support has not yet been implemented in Psi4!");
    }
   
    // always executing on host (i.e., CPU) for now
    rt_ = std::make_unique<GauXC::RuntimeEnvironment>(GAUXC_MPI_CODE(MPI_COMM_WORLD));

    // convert Psi4 fundamental quantities to GauXC 
    auto gauxc_mol_ = psi4_to_gauxc_molecule(primary_->molecule());
    auto gauxc_primary_ = psi4_to_gauxc_basisset<double>(primary_);
    
    // create snLinK grid for GauXC
    //gauxc_grid_ = psi4_to_gauxc_grid();  
    
    gauxc_grid_ = std::make_unique<GauXC::MolGrid>(
        GauXC::MolGridFactory::create_default_molgrid(
            gauxc_mol_, 
            pruning_scheme_map_[pruning_scheme_],
            GauXC::BatchSize(512), 
            radial_scheme_map_[radial_scheme_], 
            GauXC::RadialSize(radial_points_),
            GauXC::AngularSize(spherical_points_)
        )
    );
  
    // construct load balancer
    const std::string load_balancer_kernel = "Default";
    const size_t quad_pad_value = 1;

    gauxc_load_balancer_factory_ = std::make_unique<GauXC::LoadBalancerFactory>(ex, load_balancer_kernel);
    std::shared_ptr<GauXC::LoadBalancer> gauxc_load_balancer = gauxc_load_balancer_factory_->get_shared_instance(*rt_, gauxc_mol_, *gauxc_grid_, gauxc_primary_, quad_pad_value);
    
    // construct weights module
    const std::string mol_weights_kernel = "Default";
    gauxc_mol_weights_factory_ = std::make_unique<GauXC::MolecularWeightsFactory>(ex, mol_weights_kernel, GauXC::MolecularWeightsSettings{});
    std::shared_ptr<GauXC::MolecularWeights> gauxc_mol_weights = gauxc_mol_weights_factory_->get_shared_instance();

    // apply partition weights
    gauxc_mol_weights->modify_weights(*gauxc_load_balancer);
    
    // set integrator options 
    bool do_no_screen = options.get_str("SCREENING") == "NONE";
    integrator_settings_.screen_ek = (do_no_screen || cutoff_ == 0.0) ? false : true;
    
    // TODO: Check correctness of these mappings
    auto ints_tol = options.get_double("SNLINK_INTS_TOLERANCE");
    integrator_settings_.energy_tol = ints_tol; 
    integrator_settings_.k_tol = ints_tol;

    const std::string integrator_input_type = "Replicated";
    const std::string integrator_kernel = "Default";
    const std::string reduction_kernel = "Default";
    const std::string lwd_kernel = "Default";
    
    // construct dummy functional
    // note that the snLinK used by the eventually-called eval_exx is agnostic to the input functional!
    const std::string dummy_func_str = "B3LYP";
    auto spin = options_.get_str("REFERENCE") == "UHF" ? ExchCXX::Spin::Polarized : ExchCXX::Spin::Unpolarized;

    ExchCXX::Functional dummy_func_key = ExchCXX::functional_map.value(dummy_func_str);
    ExchCXX::XCFunctional dummy_func = { ExchCXX::Backend::builtin, dummy_func_key, spin };  
   
    // actually construct integrator 
    integrator_factory_ = std::make_unique<GauXC::XCIntegratorFactory<matrix_type> >(ex, integrator_input_type, integrator_kernel, lwd_kernel, reduction_kernel);
    integrator_ = integrator_factory_->get_shared_instance(dummy_func, *gauxc_load_balancer); 

    // sanity-check integrator
     
}

snLinK::~snLinK() {}

void snLinK::print_header() const {
    if (print_) {
        outfile->Printf("\n");
        outfile->Printf("  ==> snLinK: GauXC Semi-Numerical Linear Exchange K <==\n\n");
        
        outfile->Printf("    K Grid Radial Points: %i\n", radial_points_);
        outfile->Printf("    K Grid Spherical/Angular Points: %i\n", spherical_points_);
        outfile->Printf("    K Screening?:     %s\n", (integrator_settings_.screen_ek) ? "Yes" : "No");
        outfile->Printf("    K Basis Cutoff:     %11.0E\n", basis_tol_);
        outfile->Printf("    K Ints Cutoff:     %11.0E\n", integrator_settings_.energy_tol);
        if (debug_) {
          outfile->Printf("    (Debug) K Grid Pruning Scheme:     %s\n", pruning_scheme_.c_str());
          outfile->Printf("    (Debug) K Radial Quadrature Scheme:     %s\n", radial_scheme_.c_str());
        }
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
        auto K_eigen = integrator_->eval_exx(D_eigen, integrator_settings_);

        //Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");
       
        //auto K_eigen_block = K_eigen.block<20, 20>(0, 0); 
        //std::cout << K_eigen_block.format(HeavyFmt) << std::endl;
 
        // convert result to Psi4 matrix
        if (incfock_iter_) { 
            auto delta_K = std::make_shared<Matrix>(K[iD]->nrow(), K[iD]->ncol());
            eigen_to_psi4_matrix(delta_K, K_eigen);

            K[iD]->add(delta_K);
        } else {
            eigen_to_psi4_matrix(K[iD], K_eigen);
        }
    }
    
    return;
}

}  // namespace psi
