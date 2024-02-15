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

// constructs a permutation matrix for converting matrices to and from GauXC's integral ordering standard 
//template <typename T>
Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> snLinK::generate_permutation_matrix(const std::shared_ptr<BasisSet> psi4_basisset) {
  
  Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> permutation_matrix(psi4_basisset->nbf());
  
  // general array for how to reorder integrals 
  constexpr int max_am = 7;
  std::array<std::vector<int>, max_am> cca_integral_order; 

  // s shell, easy
  cca_integral_order[0] = { 0 }; 

  // p shell, actually convoluted
  if constexpr (is_cca_) {
#ifndef USING_gauxc_CCA
    throw PSIEXCEPTION("Shouldn't be here!");
#endif
    // dont even know what this is, but it works I guess?
    // GauXC internally represents p shells as cartesians even in spherical
    // harmonic basis sets
    // perhaps we just disagree on cartesian ordering of p shells???
    cca_integral_order[1] = { 1, -1, 0 }; 
    // forced cartesian is usual CCA ordering
    if (force_cartesian_) cca_integral_order[1] = { -1, 0, 1 }; 
  } else {
#ifndef USING_gauxc_GAUSSIAN
    throw PSIEXCEPTION("Shouldn't be here!");
#endif

    cca_integral_order[1] = { 1, -1, 0 };
    //cca_integral_order[1] = { 0, 1, -1 };
  }
 
  // d shells or larger
  for (size_t l = 2; l != max_am; ++l) {
    cca_integral_order[l] = std::vector<int>(2*l + 1, 0);
    for (size_t idx = 1, val = 1; idx < cca_integral_order[l].size(); idx += 2, ++val) {
      cca_integral_order.at(l)[idx] = val;
      cca_integral_order.at(l)[idx + 1] = -val;
    }
  }

  //if constexpr (is_cca_) {
  //  for (int ibf = 0; ibf != gauxc_basisset.nbf(); ++ibf) {
  //    permutation_matrix.indices()[ibf] = ibf;
  //  }
  //} else {
  for (int ish = 0, ibf = 0; ish != psi4_basisset->nshell(); ++ish) {
    auto& sh = psi4_basisset->shell(ish);
    auto am = sh.am();

    auto ibf_base = ibf;
    for (int ishbf = 0; ishbf != 2*am + 1; ++ishbf, ++ibf) {
      permutation_matrix.indices()[ibf] = ibf_base + cca_integral_order[am][ishbf] + am;
    }
  }
  
  std::cout << "Permutation Matrix (" << permutation_matrix.rows() << ", " << permutation_matrix_.cols() << "):" << std::endl;
  std::cout << "----------------------------  " << std::endl;
  std::cout << "    final_index_ref = [ " << std::endl;
  constexpr int num_col = 13;
  for (int irow = 0; irow != permutation_matrix.rows()/num_col; ++irow) {
    std::cout << "      ";
    for (int icol = 0; icol != num_col; ++icol) {
      std::cout << permutation_matrix.indices()[num_col * irow + icol] << ", ";
      if (icol == num_col - 1) std::cout << std::endl; 
    }
  }
  std::cout << "  ]" << std::endl;
  return std::move(permutation_matrix);
}

/*
// constructs a permutation matrix for converting matrices to and from GauXC's integral ordering standard 
template <typename T>
Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> snLinK::generate_permutation_matrix(const GauXC::BasisSet<T>& gauxc_basisset) {
  
  Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> permutation_matrix(gauxc_basisset.nbf());
  
  // general array for how to reorder integrals 
  constexpr int max_am = 7;
  std::array<std::vector<int>, max_am> cca_integral_order; 

  // s shell, easy
  cca_integral_order[0] = { 0 }; 

  // p shell, actually convoluted
  if constexpr (is_cca_) {
#ifndef USING_gauxc_CCA
    throw PSIEXCEPTION("Shouldn't be here!");
#endif
    // dont even know what this is, but it works I guess?
    // GauXC internally represents p shells as cartesians even in spherical
    // harmonic basis sets
    // perhaps we just disagree on cartesian ordering of p shells???
    cca_integral_order[1] = { 1, -1, 0 }; 
    // forced cartesian is usual CCA ordering
    if (force_cartesian_) cca_integral_order[1] = { -1, 0, 1 }; 
  } else {
#ifndef USING_gauxc_GAUSSIAN
    throw PSIEXCEPTION("Shouldn't be here!");
#endif

    cca_integral_order[1] = { 1, -1, 0 };
  }
 
  // d shells or larger
  for (size_t l = 2; l != max_am; ++l) {
    cca_integral_order[l] = std::vector<int>(2*l + 1, 0);
    for (size_t idx = 1, val = 1; idx < cca_integral_order[l].size(); idx += 2, ++val) {
      cca_integral_order.at(l)[idx] = val;
      cca_integral_order.at(l)[idx + 1] = -val;
    }
  }

  //if constexpr (is_cca_) {
  //  for (int ibf = 0; ibf != gauxc_basisset.nbf(); ++ibf) {
  //    permutation_matrix.indices()[ibf] = ibf;
  //  }
  //} else {
  for (int ish = 0, ibf = 0; ish != gauxc_basisset.size(); ++ish) {
    auto& sh = gauxc_basisset[ish];
    auto am = sh.l();

    auto ibf_base = ibf;
    for (int ishbf = 0; ishbf != 2*am + 1; ++ishbf, ++ibf) {
      permutation_matrix.indices()[ibf] = ibf_base + cca_integral_order[am][ishbf] + am;
    }
  }
  
  std::cout << "Permutation Matrix (" << permutation_matrix.rows() << ", " << permutation_matrix_.cols() << "):" << std::endl;
  std::cout << "----------------------------  " << std::endl;
  std::cout << permutation_matrix.indices() << std::endl << std::endl;
 
  return std::move(permutation_matrix);
}
*/

// converts a Psi4::Molecule object to a GauXC::Molecule object
GauXC::Molecule snLinK::psi4_to_gauxc_molecule(std::shared_ptr<Molecule> psi4_molecule) {
    GauXC::Molecule gauxc_molecule;
 
    outfile->Printf("  snLinK::psi4_to_gauxc_molecule (GauXC-side)\n");
    outfile->Printf("  -------------------------------------------\n");
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
GauXC::BasisSet<T> snLinK::psi4_to_gauxc_basisset(std::shared_ptr<BasisSet> psi4_basisset, bool force_cartesian) {
    using prim_array = typename GauXC::Shell<T>::prim_array;
    using cart_array = typename GauXC::Shell<T>::cart_array;

    //GauXC::BasisSet<T> gauxc_basisset;
    GauXC::BasisSet<T> gauxc_basisset(psi4_basisset->nshell());
 
    outfile->Printf("  snLinK::psi4_to_gauxc_basisset (Psi-side)\n");
    outfile->Printf("  -----------------------------------------\n");
    for (size_t ishell = 0; ishell != psi4_basisset->nshell(); ++ishell) {
        auto psi4_shell = psi4_basisset->shell(ishell);
       
        const auto nprim = GauXC::PrimSize(psi4_shell.nprimitive());
        prim_array alpha; 
        prim_array coeff;

        outfile->Printf("  ");
        outfile->Printf("%s", (force_cartesian || psi4_shell.is_cartesian()) ? "Cartesian" : "Spherical");
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

        //gauxc_basisset[permutation_matrix_.indices()[ishell]] = GauXC::Shell(
        gauxc_basisset[ishell] = GauXC::Shell(
            nprim,
            // TODO: check if am() is 0-indexed
            // Answer: It is! See tests/basisset_test.cxx in GauXC
            GauXC::AngularMomentum(psi4_shell.am()), 
            //(force_cartesian ? GauXC::SphericalType(false) : GauXC::SphericalType( !(psi4_shell.is_cartesian()) && psi4_shell.am() >= 2) ),
            (force_cartesian ? GauXC::SphericalType(false) : GauXC::SphericalType( !(psi4_shell.is_cartesian()) ) ),
            alpha,
            coeff,
            center,
            false // do not normalize shell via GauXC; it is normalized via Psi4
        );
    }
    
    for (auto& sh : gauxc_basisset) {
        sh.set_shell_tolerance(basis_tol_); 
    }

    outfile->Printf("  snLinK::psi4_to_gauxc_basisset (GauXC-side)\n");
    outfile->Printf("  -------------------------------------------\n");
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
    incfock_iter_ = false;
  
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

    // create matrix for spherical-to-cartesian matrix transformations
    const auto factory = std::make_shared<IntegralFactory>(primary_, primary_, primary_, primary);
    PetiteList petite(primary_, factory, true);
    sph_to_cart_matrix_ = petite.aotoso()->clone()->transpose();

    // create permutation matrix to handle integral ordering
    permutation_matrix_ = generate_permutation_matrix(primary_);
    force_permute_ = options_.get_bool("SNLINK_FORCE_PERMUTE"); 

    //if (force_cartesian_ && !is_cca_) {
    //{
    //    auto sph_to_cart_order_shift_buffer = psi4_to_eigen_map(sph_to_cart_matrix_);
    //    sph_to_cart_order_shift_buffer = permutation_matrix_ * sph_to_cart_order_shift_buffer; 
    //}

    // convert Psi4 fundamental quantities to GauXC 
    auto gauxc_mol = psi4_to_gauxc_molecule(primary_->molecule());
    auto gauxc_primary = psi4_to_gauxc_basisset<double>(primary_, force_cartesian_);
    auto gauxc_primary_spherical = psi4_to_gauxc_basisset<double>(primary_, false);

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
    
    auto gauxc_load_balancer_spherical = gauxc_load_balancer_factory_->get_instance(*rt, gauxc_mol, gauxc_grid, gauxc_primary_spherical, quad_pad_value);

    // construct weights module
    const std::string mol_weights_kernel = options_.get_str("SNLINK_MOL_WEIGHTS_KERNEL");
    gauxc_mol_weights_factory_ = std::make_unique<GauXC::MolecularWeightsFactory>(ex, mol_weights_kernel, GauXC::MolecularWeightsSettings{});
    auto gauxc_mol_weights = gauxc_mol_weights_factory_->get_instance();
    auto gauxc_mol_weights_spherical = gauxc_mol_weights_factory_->get_instance();

    // apply partition weights
    gauxc_mol_weights.modify_weights(gauxc_load_balancer);
    gauxc_mol_weights_spherical.modify_weights(gauxc_load_balancer_spherical);

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

    //integrator_spherical_ = nullptr;
    integrator_spherical_ = integrator_factory_->get_shared_instance(dummy_func, gauxc_load_balancer_spherical); 
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

    // we need to know fi we are using a spherical harmonic basis
    // much of the behavior here is influenced by this
    auto is_spherical_basis = primary_->has_puream();
   
    constexpr bool print_D = false;
    constexpr bool print_K = true;
 
    // compute K for density Di using GauXC
    for (int iD = 0; iD != D.size(); ++iD) {
        if (integrator_spherical_ != nullptr) {
        auto D_eigen_pre = psi4_to_eigen_map(D[iD]);    
        
        if constexpr(print_D) {
        std::cout << "D_eigen (pre-reorder): " << std::endl;
        std::cout << "------------------------- " << std::endl;
        std::cout << D_eigen_pre.block<6, 6>(0,0) << std::endl;
        }

        auto n_el_pre_d = integrator_spherical_->integrate_den(D_eigen_pre);
        outfile->Printf("  Integrator_den (pre-reorder): %f\n", n_el_pre_d);

        auto Db = D[iD]->clone();
        auto DSb = linalg::doublet(Db, S_, false, false);
        auto n_el_pre_b = DSb->trace();
        outfile->Printf("  trace(DS) (pre-reorder): %f\n", n_el_pre_b);
        }

        // map Psi4 density matrix to Eigen matrix map
        SharedMatrix D_eigen_buffer = nullptr; 
        if (is_spherical_basis) {
            // need to reorder Psi4 density matrix to CCA ordering if in spherical harmonics
            if constexpr (!is_cca_) {
                auto D_order_shift_buffer = psi4_to_eigen_map(D[iD]);
                D_order_shift_buffer = permutation_matrix_ * D_order_shift_buffer * permutation_matrix_.transpose();

                if (integrator_spherical_ != nullptr) {
                if constexpr (print_D) {
                std::cout << "D_eigen (post-reorder): " << std::endl;
                std::cout << "------------------------- " << std::endl;
                std::cout << D_order_shift_buffer.block<6, 6>(0,0) << std::endl;
                }

                auto n_el_post_a = integrator_spherical_->integrate_den(D_order_shift_buffer);
                outfile->Printf("  Integrator_den (post-reorder): %f\n", n_el_post_a);
       
                auto Dc = D[iD]->clone();
                auto DSc = linalg::doublet(Dc, S_, false, false);
                auto n_el_post_c = DSc->trace();
                outfile->Printf("  trace(DS) (post-reorder): %f\n", n_el_post_c);
                }
            }
            // also need to transform D to cartesian coordinates if requested/required
            if (force_cartesian_) {
                D_eigen_buffer = std::make_shared<Matrix>(sph_to_cart_matrix_->nrow(), sph_to_cart_matrix_->nrow());
                D_eigen_buffer->transform(D[iD], sph_to_cart_matrix_);
            } else {
                D_eigen_buffer = D[iD];
            }
        // natively-cartesian bases dont require the reordering/transforming above
        } else {
            D_eigen_buffer = D[iD];
        }
        auto D_eigen = psi4_to_eigen_map(D_eigen_buffer);    
        if constexpr (print_D) {
        std::cout << "D_eigen (post-transform): " << std::endl;
        std::cout << "------------------------- " << std::endl;
        std::cout << D_eigen.block<6, 6>(0,0) << std::endl;
        }

        auto n_el_post_a = integrator_->integrate_den(D_eigen);
        outfile->Printf("  Integrator_den (post-transform): %f\n", n_el_post_a);

        auto Dc = D[iD]->clone();
        auto DSc = linalg::doublet(Dc, S_, false, false);
        auto n_el_post_c = DSc->trace();
        outfile->Printf("  trace(DS) (post-transform): %f\n", n_el_post_c);
        
        // map Psi4 exchange matrix buffer to Eigen matrix map
        // buffer can be either K itself or a cartesian representation of K
        SharedMatrix K_eigen_buffer = nullptr; 
        if (is_spherical_basis) {
            // need to reorder Psi4 exchange matrix to CCA ordering if in spherical harmonics
            if constexpr (!is_cca_) {
                auto K_order_shift_buffer = psi4_to_eigen_map(K[iD]);
                K_order_shift_buffer = permutation_matrix_ * K_order_shift_buffer * permutation_matrix_.transpose();
            }

            // also need to transform D to cartesian coordinates if requested/required
            if (force_cartesian_) { 
                K_eigen_buffer = std::make_shared<Matrix>(sph_to_cart_matrix_->ncol(), sph_to_cart_matrix_->ncol());
                K_eigen_buffer->transform(K[iD], sph_to_cart_matrix_);
            } else {
                K_eigen_buffer = K[iD];
            }
        // natively-cartesian bases dont require the reordering/transforming above
        } else {
            K_eigen_buffer = K[iD];
        }
        auto K_eigen = psi4_to_eigen_map(K_eigen_buffer); 
        
        outfile->Printf("  Done.\n"); 
 
        // compute delta K if incfock iteration... 
        if (incfock_iter_) { 
            //throw PSIEXCEPTION("NOPE!");
            // if cartesian transformation is forced, K_eigen is delta K and must be added to Psi4 K separately...
            if (force_cartesian_ && is_spherical_basis) {
                K_eigen = integrator_->eval_exx(D_eigen, integrator_settings_);
                K_eigen_buffer->back_transform(sph_to_cart_matrix_);
                K[iD]->add(K_eigen_buffer);
            // ... otherwise the computation and addition can be bundled together 
            } else {
                K_eigen += integrator_->eval_exx(D_eigen, integrator_settings_);
            }

        // ... else compute full K 
        } else {
            K_eigen = integrator_->eval_exx(D_eigen, integrator_settings_);        
       
            if constexpr (print_K) { 
            std::cout << "K_eigen (pre-transform): " << std::endl;
            std::cout << "------------------------- " << std::endl;
            std::cout << K_eigen.block<6, 6>(0,0) << std::endl;
            std::cout << std::endl;
            std::cout << K_eigen.block<5, 5>(17,0) << std::endl;
            }

            if (force_cartesian_ && is_spherical_basis) {
                K[iD]->back_transform(K_eigen_buffer, sph_to_cart_matrix_);          
            }

            if constexpr (print_K) {
            std::cout << "K_eigen (post-transform): " << std::endl;
            std::cout << "------------------------- " << std::endl;
            std::cout << K_eigen.block<6, 6>(0,0) << std::endl;
            std::cout << std::endl;
            std::cout << K_eigen.block<5, 5>(17,0) << std::endl;
            }
        }

        // now we need to reverse the CCA reordering previously performed
        if (is_spherical_basis) {
            //for (int i = 0; i != 5; ++i) { for (int j = 0; j != 5; ++j) { std::cout << K[iD]->get(i,j) << ", "; }; }; std::cout << std::endl << std::endl;
            if constexpr (!is_cca_) {
                auto D_order_shift_buffer = psi4_to_eigen_map(D[iD]);
                if constexpr (print_D) {
                std::cout << "D_eigen (pre-back-reorder): " << std::endl;
                std::cout << "------------------------- " << std::endl;
                std::cout << D_order_shift_buffer.block<6, 6>(0,0) << std::endl;
                }
 
                D_order_shift_buffer = permutation_matrix_.transpose() * D_order_shift_buffer * permutation_matrix_; 

                if constexpr (print_D) {
                std::cout << "D_eigen (post-back-reorder): " << std::endl;
                std::cout << "------------------------- " << std::endl;
                std::cout << D_order_shift_buffer.block<6, 6>(0,0) << std::endl;
                }

                auto K_order_shift_buffer = psi4_to_eigen_map(K[iD]);
                if constexpr (print_K) {
                std::cout << "K_eigen (pre-back-reorder): " << std::endl;
                std::cout << "------------------------- " << std::endl;
                std::cout << K_order_shift_buffer.block<6, 6>(0,0) << std::endl;
                std::cout << std::endl;
                std::cout << K_eigen.block<5, 5>(17,0) << std::endl;
                }
 
                //K_order_shift_buffer = permutation_matrix_ * K_order_shift_buffer * permutation_matrix_.transpose(); 
                K_order_shift_buffer = permutation_matrix_.transpose() * K_order_shift_buffer * permutation_matrix_; 
                if constexpr (print_K) {
                std::cout << "K_eigen (post-back-reorder): " << std::endl;
                std::cout << "------------------------- " << std::endl;
                std::cout << K_order_shift_buffer.block<6, 6>(0,0) << std::endl;
                std::cout << std::endl;
                std::cout << K_eigen.block<5, 5>(17,0) << std::endl;
                }
            }

            //if (force_cartesian_) {
            //    ;
                //D_order_shift_buffer = permutation_matrix_.transpose() * D_order_shift_buffer * permutation_matrix_; 
            //} else {
                //D_order_shift_buffer = permutation_matrix_.transpose() * D_order_shift_buffer * permutation_matrix_; 
            //}
            
//            if (force_cartesian_ || is_cca_) {
//                ;
                //K_order_shift_buffer = permutation_matrix_.transpose() * K_order_shift_buffer * permutation_matrix_; 
//            } else {
//                auto K_order_shift_buffer = psi4_to_eigen_map(K[iD]);
//                K_order_shift_buffer = permutation_matrix_.transpose() * K_order_shift_buffer * permutation_matrix_; 
//            }
            //for (int i = 0; i != 5; ++i) { for (int j = 0; j != 5; ++j) { std::cout << K[iD]->get(i,j) << ", "; }; }; std::cout << std::endl << std::endl;
        }
      
        // symmetrize K if applicable
        if (lr_symmetric_) {
            K[iD]->hermitivitize();
        }

        //outfile->Printf("K[%i] norm: %f \n", iD, K[iD]->norm());   
    }
   
    return;
}

}  // namespace psi
