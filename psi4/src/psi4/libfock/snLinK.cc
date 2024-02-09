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
#include "psi4/libmints/petitelist.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/lib3index/dftensor.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include <algorithm>
#include <iostream>
#include <map>
#include <tuple>
#include <unordered_set>
#include <variant>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace psi;

namespace psi {

// creates maps for translating Psi4 option inputs to GauXC enums
std::tuple<
    std::unordered_map<std::string, GauXC::PruningScheme>, 
    std::unordered_map<std::string, GauXC::RadialQuad> 
> snLinK::generate_enum_mappings() {
    // generate map for grid pruning schemes 
    std::unordered_map<std::string, GauXC::PruningScheme> pruning_scheme_map; 
    pruning_scheme_map["ROBUST"] = GauXC::PruningScheme::Robust;
    pruning_scheme_map["TREUTLER"] = GauXC::PruningScheme::Treutler;
    pruning_scheme_map["NONE"] = GauXC::PruningScheme::Unpruned;

    // generate map for radial quadrature schemes 
    std::unordered_map<std::string, GauXC::RadialQuad> radial_scheme_map; 
    radial_scheme_map["TREUTLER"] = GauXC::RadialQuad::TreutlerAldrichs;
    radial_scheme_map["MURA"] = GauXC::RadialQuad::MuraKnowles;
    // TODO: confirm the correctness of this specific mapping
    // Answer: Yep, it is! The Murray, Handy, Laming literature reference
    // is mentioned in cubature.cc
    radial_scheme_map["EM"] = GauXC::RadialQuad::MurrayHandyLaming; 

    // we are done
    //return std::move(pruning_scheme_map), std::move(radial_scheme_map);
    return std::move(std::make_tuple(pruning_scheme_map, radial_scheme_map));
}

template <typename T>
Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> snLinK::generate_permutation_matrix(const GauXC::BasisSet<T>& gauxc_basisset) {
  std::array<std::vector<int>, 5> am_to_mapping; // index i represents the shell of AM i
  am_to_mapping[0] = { 0 }; // s
  am_to_mapping[1] = { 0, 1, -1 }; // p
  am_to_mapping[2] = { 0, 1, 2, -2, -1 }; // d

  std::array<std::vector<int>, 7> am_to_mapping_auto; // index i represents the shell of AM i
  for (size_t am = 0; am != am_to_mapping_auto.size(); ++am) {
    am_to_mapping_auto[am] = std::vector<int>(2*am + 1, 0); // s
    for (size_t l = 1; l < am_to_mapping_auto[am].size(); l += 2) {
      am_to_mapping_auto[am][l] = l;
      am_to_mapping_auto[am][l+1] = -l;
    }
  }

  Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> permutation_matrix(primary_->nbf());

  //permutation_matrix.setIdentity(nbf);
  for (int ish = 0, ibf = 0; ish != gauxc_basisset.size(); ++ish) {
      auto& sh = gauxc_basisset[ish];

      auto sh_am_mapping = am_to_mapping[sh.l()];

      auto ibf_base = ibf;
      for (int ishbf = 0; ishbf != sh_am_mapping.size(); ++ishbf, ++ibf) {
        permutation_matrix.indices()[ibf] = ibf_base + sh_am_mapping[ishbf] + sh.l(); 
      }
  }
  
  return std::move(permutation_matrix);
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
 
    outfile->Printf("snLinK::psi4_to_gauxc_basisset (Psi-side)\n");
    outfile->Printf("-----------------------------------------\n");
    for (size_t ishell = 0; ishell != psi4_basisset->nshell(); ++ishell) {
        auto psi4_shell = psi4_basisset->shell(ishell);
       
        const auto nprim = GauXC::PrimSize(psi4_shell.nprimitive());
        prim_array alpha; 
        prim_array coeff;

        outfile->Printf("  ");
        outfile->Printf("%s", (force_cartesian_ || psi4_shell.is_cartesian()) ? "Cartesian" : "Spherical");
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
            (force_cartesian_ ? GauXC::SphericalType(false) : GauXC::SphericalType(!(psi4_shell.is_cartesian())) ),
            alpha,
            coeff,
            center,
            false // do not normalize shell via GauXC; it is normalized via Psi4
        );
    }
    
    for (auto& sh : gauxc_basisset) {
        sh.set_shell_tolerance(basis_tol_); 
    }

    outfile->Printf("snLinK::psi4_to_gauxc_basisset (GauXC-side)\n");
    outfile->Printf("-------------------------------------------\n");
    for (int ish = 0; ish != gauxc_basisset.size(); ++ish) {
      auto& sh = gauxc_basisset[ish];

      auto alpha = sh.alpha();
      auto coeff = sh.coeff();
      assert(alpha.size() == coeff.size());
      assert(alpha.size() == sh.nprim());

      outfile->Printf("  ");
      outfile->Printf("  %s", (!sh.pure() ? "Cartesian" : "Spherical"));
      outfile->Printf(" Shell #%i (AM %i): %i primitives\n", ish, sh.l(), sh.nprim());
      for (size_t iprim = 0; iprim != sh.nprim(); ++iprim) {
        outfile->Printf("      Primitive #%i: %f, %f\n", iprim, alpha.at(iprim), coeff.at(iprim));
      }
    }
 
    return std::move(gauxc_basisset);
}

// converts a Psi4::Matrix to an Eigen::MatrixXd map
Eigen::Map<Eigen::MatrixXd> snLinK::psi4_to_eigen_map(SharedMatrix psi4_matrix) {
    // Sanity checks
    // only works for C1 symmetry at the moment
    // could be improved if symmetry is utilized in JK builds in the future
    if (psi4_matrix->nirrep() != 1) {
        throw PSIEXCEPTION("Psi4::Matrix must be in C1 symmetry to be transformed into Eigen::MatrixXd!");
    }

    // create Eigen matrix "map" using Psi4 matrix data array directly
    return std::move(
        Eigen::Map<Eigen::MatrixXd>(
            psi4_matrix->pointer()[0], // pointer to first element of underlying Psi4::Matrix data array
            psi4_matrix->nrow(), psi4_matrix->ncol()
        )
    );
}

snLinK::snLinK(std::shared_ptr<BasisSet> primary, Options& options) : SplitJK(primary, options) {
   
    // set Psi4-specific parameters
    radial_points_ = options_.get_int("SNLINK_RADIAL_POINTS");
    spherical_points_ = options_.get_int("SNLINK_SPHERICAL_POINTS");
    basis_tol_ = options.get_double("SNLINK_BASIS_TOLERANCE");

    pruning_scheme_ = options_.get_str("SNLINK_PRUNING_SCHEME");
    radial_scheme_ = options_.get_str("DFT_RADIAL_SCHEME");
 
    format_ = Eigen::IOFormat(4, 0, ", ", "\n", "[", "]");
    
    // sanity-checking of radial scheme needed because generalist DFT_RADIAL_SCHEME option is used 
    std::array<std::string, 3> valid_radial_schemes = { "TREUTLER", "MURA", "EM" };
    bool is_valid_radial_scheme = std::any_of(
      valid_radial_schemes.cbegin(),
      valid_radial_schemes.cend(),
      [&](std::string valid_radial_scheme) { return radial_scheme_.find(valid_radial_scheme) != std::string::npos; }
    );

    if (!is_valid_radial_scheme) {
        throw PSIEXCEPTION("Invalid radial quadrature scheme selected for snLinK! Please set DFT_RADIAL_SCHEME to either TREUTLER, MURA, or EM.");
    }

    // create mappings for GauXC pruning and radial scheme enums
    auto [ pruning_scheme_map, radial_scheme_map ] = generate_enum_mappings();

    // define runtime environment and execution space
    use_gpu_ = options_.get_bool("SNLINK_USE_GPU"); 
    auto ex = use_gpu_ ? GauXC::ExecutionSpace::Device : GauXC::ExecutionSpace::Host;  

    std::unique_ptr<GauXC::RuntimeEnvironment> rt = nullptr; 
#ifdef USING_gauxc_GPU
    if (use_gpu_) {
        // 0.9 indicates to use maximum 90% of maximum GPU memory, I think?
        rt = std::make_unique<GauXC::DeviceRuntimeEnvironment>( GAUXC_MPI_CODE(MPI_COMM_WORLD,) 0.9 );
    } else { 
        rt = std::make_unique<GauXC::RuntimeEnvironment>( GAUXC_MPI_CODE(MPI_COMM_WORLD) );
    }
#else
        rt = std::make_unique<GauXC::RuntimeEnvironment>( GAUXC_MPI_CODE(MPI_COMM_WORLD) );
#endif

    // set whether we force use of cartesian coordinates or not
    // this is required for GPU execution when using spherical harmonic basis sets
    force_cartesian_ = options_.get_bool("SNLINK_FORCE_CARTESIAN"); 
    if (use_gpu_ && !force_cartesian_ && primary_->has_puream()) {
        throw PSIEXCEPTION("GPU snLinK must be executed with SNLINK_FORCE_CARTESIAN=true when using spherical harmonic basis sets!");  
    }

    const auto factory = std::make_shared<IntegralFactory>(primary_, primary_, primary_, primary);
    PetiteList petite(primary_, factory, true);
    sph_to_cart_matrix_ = petite.aotoso()->transpose();

    // convert Psi4 fundamental quantities to GauXC 
    auto gauxc_mol = psi4_to_gauxc_molecule(primary_->molecule());
    auto gauxc_primary = psi4_to_gauxc_basisset<double>(primary_);

    // create permutation matrix to handle integral ordering
    permutation_matrix_ = generate_permutation_matrix(gauxc_primary);
    force_permute_ = options_.get_bool("SNLINK_FORCE_PERMUTE"); 
 
    //std::cout << "Permutation Matrix (" << permutation_matrix_.rows() << ", " << permutation_matrix_.cols() << "):" << std::endl;
    //std::cout << "----------------------------  " << std::endl;
    //std::cout << permutation_matrix_.indices() << std::endl << std::endl;

    // create snLinK grid for GauXC
    auto grid_batch_size = options_.get_int("SNLINK_GRID_BATCH_SIZE");
    auto gauxc_grid = GauXC::MolGridFactory::create_default_molgrid(
        gauxc_mol, 
        pruning_scheme_map[pruning_scheme_],
        GauXC::BatchSize(grid_batch_size), 
        radial_scheme_map[radial_scheme_], 
        //GauXC::RadialSize(radial_points_),
        //GauXC::AngularSize(spherical_points_)
        GauXC::AtomicGridSizeDefault::UltraFineGrid
    );

    // construct load balancer
    const std::string load_balancer_kernel = options_.get_str("SNLINK_LOAD_BALANCER_KERNEL"); 
    const size_t quad_pad_value = 1;

    gauxc_load_balancer_factory_ = std::make_unique<GauXC::LoadBalancerFactory>(ex, load_balancer_kernel);
    auto gauxc_load_balancer = gauxc_load_balancer_factory_->get_instance(*rt, gauxc_mol, gauxc_grid, gauxc_primary, quad_pad_value);

    // construct weights module
    const std::string mol_weights_kernel = options_.get_str("SNLINK_MOL_WEIGHTS_KERNEL");
    gauxc_mol_weights_factory_ = std::make_unique<GauXC::MolecularWeightsFactory>(ex, mol_weights_kernel, GauXC::MolecularWeightsSettings{});
    auto gauxc_mol_weights = gauxc_mol_weights_factory_->get_instance();

    // apply partition weights
    gauxc_mol_weights.modify_weights(gauxc_load_balancer);

    // set integrator options 
    bool do_no_screen = options.get_str("SCREENING") == "NONE";
    integrator_settings_.screen_ek = (do_no_screen || cutoff_ == 0.0) ? false : true;

    // TODO: Check correctness of these mappings
    auto ints_tol = options.get_double("SNLINK_INTS_TOLERANCE");
    integrator_settings_.energy_tol = ints_tol; 
    integrator_settings_.k_tol = ints_tol;

    const std::string integrator_input_type = "Replicated";
    const std::string integrator_kernel = options_.get_str("SNLINK_INTEGRATOR_KERNEL");  
    const std::string reduction_kernel = options_.get_str("SNLINK_REDUCTION_KERNEL");  
    const std::string lwd_kernel = options_.get_str("SNLINK_LWD_KERNEL");  

    // construct dummy functional
    // note that the snLinK used by the eventually-called eval_exx function is agnostic to the input functional!
    const std::string dummy_func_str = "B3LYP";
    auto spin = options_.get_str("REFERENCE") == "UHF" ? ExchCXX::Spin::Polarized : ExchCXX::Spin::Unpolarized;

    ExchCXX::Functional dummy_func_key = ExchCXX::functional_map.value(dummy_func_str);
    ExchCXX::XCFunctional dummy_func = { ExchCXX::Backend::builtin, dummy_func_key, spin };  

    // actually construct integrator 
    integrator_factory_ = std::make_unique<GauXC::XCIntegratorFactory<matrix_type> >(ex, integrator_input_type, integrator_kernel, lwd_kernel, reduction_kernel);
    integrator_ = integrator_factory_->get_shared_instance(dummy_func, gauxc_load_balancer); 
}

snLinK::~snLinK() {}

void snLinK::print_header() const {
    if (print_) {
        outfile->Printf("\n");
        outfile->Printf("  ==> snLinK: GauXC Semi-Numerical Linear Exchange K <==\n\n");
        
        outfile->Printf("    K Execution Space: %s\n", (use_gpu_) ? "Device" : "Host");
        outfile->Printf("    K Grid Radial Points: %i\n", radial_points_);
        outfile->Printf("    K Grid Spherical/Angular Points: %i\n", spherical_points_);
        outfile->Printf("    K Screening?:     %s\n", (integrator_settings_.screen_ek) ? "Yes" : "No");
        outfile->Printf("    K Basis Cutoff:     %11.0E\n", basis_tol_);
        outfile->Printf("    K Ints Cutoff:     %11.0E\n", integrator_settings_.energy_tol);
        //if (debug_) {
        if (true) {
          outfile->Printf("\n");    
          outfile->Printf("    (Debug) K Grid Pruning Scheme:     %s\n", pruning_scheme_.c_str());
          outfile->Printf("    (Debug) K Radial Quadrature Scheme:     %s\n", radial_scheme_.c_str());
          outfile->Printf("    (Debug) K Grid Batch Size:     %i\n\n", options_.get_int("SNLINK_GRID_BATCH_SIZE"));
          
          outfile->Printf("    (Debug) K Load Balancer Kernel:     %s\n", options_.get_str("SNLINK_LOAD_BALANCER_KERNEL").c_str());
          outfile->Printf("    (Debug) K Mol. Weights Kernel:     %s\n\n", options_.get_str("SNLINK_MOL_WEIGHTS_KERNEL").c_str());
          
          outfile->Printf("    (Debug) K Integrator Type:     %s\n", "Replicated");
          outfile->Printf("    (Debug) K Integrator Kernel:     %s\n", options_.get_str("SNLINK_INTEGRATOR_KERNEL").c_str());
          outfile->Printf("    (Debug) K Reduction Kernel:     %s\n", options_.get_str("SNLINK_REDUCTION_KERNEL").c_str());
          outfile->Printf("    (Debug) K LWD Kernel:     %s\n\n", options_.get_str("SNLINK_LWD_KERNEL").c_str());

        }
    }
}

// build the K matrix using Ochsenfeld's sn-LinK algorithm
// Implementation is provided by the external GauXC library 
void snLinK::build_G_component(std::vector<std::shared_ptr<Matrix>>& D, std::vector<std::shared_ptr<Matrix>>& K,
    std::vector<std::shared_ptr<TwoBodyAOInt> >& eri_computers) {

    auto is_spherical_basis = primary_->has_puream();
    // compute K for density Di using GauXC
    for (int iD = 0; iD != D.size(); ++iD) {
        // map Psi4 density matrix to Eigen matrix map
        // also includes spherical->cartesian transformation if specified
        SharedMatrix D_temp = nullptr; 
        if (force_cartesian_ && is_spherical_basis) {
            D_temp = std::make_shared<Matrix>(sph_to_cart_matrix_->nrow(), sph_to_cart_matrix_->nrow());
            D_temp->transform(D[iD], sph_to_cart_matrix_);
        } else {
            D_temp = D[iD];
        }
        auto D_eigen = psi4_to_eigen_map(D_temp);    
        assert(D_eigen.rows() == permutation_matrix_.rows());
        assert(D_eigen.cols() == permutation_matrix_.cols());
        //std::cout << "D_eigen pre-permute(" << D_eigen.rows() << ", " << D_eigen.cols() << "): " << std::endl;
        //std::cout << "----------------------------  " << std::endl;
        //std::cout << D_eigen << std::endl << std::endl;
        if (force_permute_ && is_spherical_basis) {
            D_eigen = permutation_matrix_ * D_eigen; 
            D_eigen = D_eigen * permutation_matrix_.transpose();
        }
        //std::cout << "D_eigen post-permute(" << D_eigen.rows() << ", " << D_eigen.cols() << "): " << std::endl;
        //std::cout << "----------------------------  " << std::endl;
        //std::cout << D_eigen << std::endl << std::endl;
 
        // map Psi4 exchange matrix buffer to Eigen matrix map
        // buffer can be either K itself or a cartesian representation of K
        SharedMatrix K_temp = nullptr; 
        if (force_cartesian_ && is_spherical_basis) {
            K_temp = std::make_shared<Matrix>(sph_to_cart_matrix_->ncol(), sph_to_cart_matrix_->ncol());
        } else {
            K_temp = K[iD];
        }
        auto K_eigen = psi4_to_eigen_map(K_temp); 
        assert(K_eigen.rows() == permutation_matrix_.rows());
        assert(K_eigen.cols() == permutation_matrix_.cols());
        //std::cout << "K_eigen pre-permute: " << std::endl;
        //std::cout << "-------------------- " << std::endl;
        //std::cout << K_eigen << std::endl << std::endl;
        if (force_permute_ && is_spherical_basis) {
            K_eigen = permutation_matrix_ * K_eigen;
            K_eigen = K_eigen * permutation_matrix_.transpose();
        }
        //std::cout << "K_eigen post-permute: " << std::endl;
        //std::cout << "-------------------- " << std::endl;
        //std::cout << K_eigen << std::endl << std::endl;
 
        // compute delta K if incfock iteration... 
        if (incfock_iter_) { 
            // if cartesian transformation is forced, K_eigen is delta K and must be added to Psi4 K separately...
            if (force_cartesian_ && is_spherical_basis) {
                K_eigen = integrator_->eval_exx(D_eigen, integrator_settings_);
               
                K_temp->back_transform(sph_to_cart_matrix_);
                K[iD]->add(K_temp);
            // ... otherwise the computation and addition can be bundled together 
            } else {
                K_eigen += integrator_->eval_exx(D_eigen, integrator_settings_);

           }
        
        // ... else compute full K 
        } else {
            K_eigen = integrator_->eval_exx(D_eigen, integrator_settings_);        

            if (force_cartesian_ && is_spherical_basis) {
                K[iD]->back_transform(K_temp, sph_to_cart_matrix_);          
            }
        }

        if (force_permute_ && is_spherical_basis) {
            D_eigen = permutation_matrix_.transpose() * D_eigen;
            D_eigen *= permutation_matrix_; 

            K_eigen = permutation_matrix_.transpose() * K_eigen;
            K_eigen = K_eigen * permutation_matrix_; 
        }

        // symmetrize K if applicable
        //if (lr_symmetric_) {
        //    K[iD]->hermitivitize();
        //}

        //outfile->Printf("K[%i] norm: %f \n", iD, K[iD]->norm());   
    }
   
    return;
}

}  // namespace psi
