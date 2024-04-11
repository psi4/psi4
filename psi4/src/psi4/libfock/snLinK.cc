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
    // The Murray, Handy, Laming literature reference is mentioned in cubature.cc
    // with association to this keyword
    radial_scheme_map["EM"] = GauXC::RadialQuad::MurrayHandyLaming; 

    // we are done
    return std::make_tuple(pruning_scheme_map, radial_scheme_map);
}

// constructs a permutation matrix for converting matrices to and from GauXC's integral ordering standard 
Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> snLinK::generate_permutation_matrix(const std::shared_ptr<BasisSet> psi4_basisset) {
  
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> permutation_matrix(psi4_basisset->nbf());
    
    // general array for how to reorder integrals 
    constexpr int max_am = 10; 
    std::array<int, 2*max_am + 1> cca_integral_order; 
  
    // s shell, easy
    cca_integral_order[0] = { 0 }; 
   
    // p shells or larger
    for (size_t l = 1, idx = 1; l != max_am; idx += 2, ++l) {
        cca_integral_order[idx] = l;
        cca_integral_order[idx + 1] = -l;
    }
  
    // actually construct permutation matrix
    for (int ish = 0, ibf = 0; ish != psi4_basisset->nshell(); ++ish) {
        auto& sh = psi4_basisset->shell(ish);
        auto am = sh.am();
  
        auto ibf_base = ibf;
        for (int ishbf = 0; ishbf != 2*am + 1; ++ishbf, ++ibf) {
            permutation_matrix.indices()[ibf] = ibf_base + cca_integral_order[ishbf] + am;
        }
    }
   
    // we are done
    return permutation_matrix;
}

// converts a Psi4::Molecule object to a GauXC::Molecule object
GauXC::Molecule snLinK::psi4_to_gauxc_molecule(std::shared_ptr<Molecule> psi4_molecule) {
    GauXC::Molecule gauxc_molecule;

    for (size_t iatom = 0; iatom != psi4_molecule->natom(); ++iatom) {
        auto atomic_number = psi4_molecule->true_atomic_number(iatom);
        auto x_coord = psi4_molecule->x(iatom);
        auto y_coord = psi4_molecule->y(iatom);
        auto z_coord = psi4_molecule->z(iatom);
        
        gauxc_molecule.emplace_back(GauXC::AtomicNumber(atomic_number), x_coord, y_coord, z_coord);
    }

    return gauxc_molecule;
}

// converts a Psi4::BasisSet object to a GauXC::BasisSet object
template <typename T>
GauXC::BasisSet<T> snLinK::psi4_to_gauxc_basisset(std::shared_ptr<BasisSet> psi4_basisset, double basis_tol, bool force_cartesian) {
    using prim_array = typename GauXC::Shell<T>::prim_array;
    using cart_array = typename GauXC::Shell<T>::cart_array;

    GauXC::BasisSet<T> gauxc_basisset(psi4_basisset->nshell());
 
    for (size_t ishell = 0; ishell != psi4_basisset->nshell(); ++ishell) {
        auto psi4_shell = psi4_basisset->shell(ishell);
       
        const auto nprim = GauXC::PrimSize(psi4_shell.nprimitive());
        prim_array alpha; 
        prim_array coeff;

        for (size_t iprim = 0; iprim != psi4_shell.nprimitive(); ++iprim) {
            alpha.at(iprim) = psi4_shell.exp(iprim);
            coeff.at(iprim) = psi4_shell.coef(iprim);
        }

        auto psi4_shell_center = psi4_shell.center();
        cart_array center = { psi4_shell_center[0], psi4_shell_center[1], psi4_shell_center[2] };

        gauxc_basisset[ishell] = GauXC::Shell(
            nprim,
            GauXC::AngularMomentum(psi4_shell.am()), 
            (force_cartesian ? GauXC::SphericalType(false) : GauXC::SphericalType( !(psi4_shell.is_cartesian()) ) ),
            alpha,
            coeff,
            center,
            false // do not normalize shell via GauXC; it is normalized via Psi4
        );
    }
    
    for (auto& sh : gauxc_basisset) {
        sh.set_shell_tolerance(basis_tol); 
    }

    return gauxc_basisset;
}

snLinK::snLinK(std::shared_ptr<BasisSet> primary, Options& options) : SplitJK(primary, options) {
    timer_on("snLinK: Setup");
    
    // => Part #1: Options Processing <= //
    timer_on("snLinK: Options Processing");

    // set a few Psi4-specific parameters
    incfock_iter_ = false;
  
    basis_tol_ = options.get_double("SNLINK_BASIS_TOLERANCE");

    pruning_scheme_ = options_.get_str("SNLINK_PRUNING_SCHEME");
    radial_scheme_ = options_.get_str("SNLINK_RADIAL_SCHEME");
 
//#if psi4_SHGSHELL_ORDERING == LIBINT_SHGSHELL_ORDERING_STANDARD
#if psi4_SHGSHELL_ORDERING == 1 
    is_cca_ = true;
//#elif psi4_SHGSHELL_ORDERING == LIBINT_SHGSHELL_ORDERING_GAUSSIAN
#elif psi4_SHGSHELL_ORDERING == 2 
    is_cca_ = false;
#else
    #error "unknown value of macro psi4_SHGSHELL_ORDERING"
#endif

    format_ = Eigen::IOFormat(4, 0, ", ", "\n", "[", "]");
    
    // create mappings for GauXC pruning and radial scheme enums
    auto [ pruning_scheme_map, radial_scheme_map ] = generate_enum_mappings();

    // define runtime environment and execution space
    use_gpu_ = options_.get_bool("SNLINK_USE_GPU"); 
    auto ex = use_gpu_ ? GauXC::ExecutionSpace::Device : GauXC::ExecutionSpace::Host;  

    std::unique_ptr<GauXC::RuntimeEnvironment> rt = nullptr; 
#ifdef GAUXC_HAS_DEVICE 
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
    auto force_cartesian = options_.get_bool("SNLINK_FORCE_CARTESIAN"); 
    if (use_gpu_ && !force_cartesian && primary_->has_puream()) {
        outfile->Printf("    INFO: GPU snLinK must be executed with SNLINK_FORCE_CARTESIAN=true when using spherical harmonic basis sets!\n");  
        outfile->Printf("    Enabling forced spherical-to-Cartesian transform...\n\n");  
        force_cartesian = true; 
    } else if (force_cartesian && !(primary_->has_puream())) {
        outfile->Printf("    INFO: SNLINK_FORCE_CARTESIAN=true has no effect when using Cartesian basis sets!\n");  
        //outfile->Printf("    Enabling forced spherical-to-Cartesian transform...\n\n");  
        force_cartesian = false; 
    }

    // create ERI ordering permutation matrix to handle integral ordering
    if (!is_cca_ && primary_->has_puream()) {
        permutation_matrix_ = std::make_optional<Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> >(
            generate_permutation_matrix(primary_)
        );
    } else {
        permutation_matrix_ = std::nullopt;
    }

    // create matrix for spherical-to-cartesian matrix transformations
    if (force_cartesian && primary_->has_puream()) {
        const auto factory = std::make_shared<IntegralFactory>(primary, primary, primary, primary);
        PetiteList petite(primary, factory, true);
        sph_to_cart_matrix_ = petite.sotoao(); 

        // SNLINK_FORCE_CARTESIAN only works with C1 symmtry currently
        // TODO: Fix this!
        if (sph_to_cart_matrix_->nirrep() != 1) {
            auto point_group = primary->molecule()->point_group();
            
            std::string message = "SNLINK_FORCE_CARTESIAN only works with C1 symmetry! ";
            message += "Current molecular point group is ";
            message += point_group->symbol();   

            throw PSIEXCEPTION(message);
        }

        // if needed, we also need to transform the Spherical part of the permutation matrix itself
        // For whatever reason, psi4_to_eigen_map doesnt seem to work specifically in the constructor
        // this is the work-around 
        if (permutation_matrix_.has_value()) { 
            auto permutation_dense = permutation_matrix_.value().toDenseMatrix();
            auto sph_to_cart_permute = std::make_shared<Matrix>(permutation_dense.rows(), permutation_dense.cols());
            
            auto ptr = sph_to_cart_permute->pointer();
            for (int irow = 0; irow != permutation_dense.rows(); ++irow) {
                for (int icol = 0; icol != permutation_dense.cols(); ++icol) {  
                    ptr[irow][icol] = permutation_dense(irow, icol);
                } 
            }
            sph_to_cart_matrix_ = linalg::doublet(sph_to_cart_permute, sph_to_cart_matrix_->to_block_sharedmatrix());
        }
    } else {
        sph_to_cart_matrix_ = nullptr;
    }
    timer_off("snLinK: Options Processing");
    
    // => Part #2: Construct GauXC Factories <= //
    timer_on("snLinK: Construct GauXC Factories");
    
    // construct load balancer and molecular weights factory factories
    const std::string load_balancer_kernel = options_.get_str("SNLINK_LOAD_BALANCER_KERNEL"); 
    gauxc_load_balancer_factory_ = std::make_unique<GauXC::LoadBalancerFactory>(ex, load_balancer_kernel);

    const std::string mol_weights_kernel = options_.get_str("SNLINK_MOL_WEIGHTS_KERNEL");
    gauxc_mol_weights_factory_ = std::make_unique<GauXC::MolecularWeightsFactory>(ex, mol_weights_kernel, GauXC::MolecularWeightsSettings{});

    // set integrator ERI screening options 
    bool do_no_screen = options.get_str("SCREENING") == "NONE";
    integrator_settings_.screen_ek = (do_no_screen || cutoff_ == 0.0) ? false : true;

    auto ints_tol = options.get_double("SNLINK_INTS_TOLERANCE");
    integrator_settings_.energy_tol = ints_tol; 
    integrator_settings_.k_tol = ints_tol;

    // set other integrator options
    const std::string integrator_input_type = "Replicated";
    const std::string integrator_kernel = options_.get_str("SNLINK_INTEGRATOR_KERNEL");  
    const std::string reduction_kernel = options_.get_str("SNLINK_REDUCTION_KERNEL");  
    const std::string lwd_kernel = options_.get_str("SNLINK_LWD_KERNEL");  

    // construct integrator factory 
    integrator_factory_ = std::make_unique<GauXC::XCIntegratorFactory<matrix_type> >(ex, integrator_input_type, integrator_kernel, lwd_kernel, reduction_kernel);
 
    timer_off("snLinK: Construct GauXC Factories");
    
    // => Part #3: Construct GauXC Integrator <= //
    timer_on("snLinK: Construct GauXC Integrator");

    radial_points_ = options_.get_int("SNLINK_RADIAL_POINTS");
    spherical_points_ = options_.get_int("SNLINK_SPHERICAL_POINTS");

    // convert Psi4 fundamental quantities to GauXC 
    auto gauxc_mol = psi4_to_gauxc_molecule(primary_->molecule());
    auto gauxc_primary = psi4_to_gauxc_basisset<double>(primary_, basis_tol_, force_cartesian);

    // create snLinK grid for GauXC
    auto grid_batch_size = options_.get_int("SNLINK_GRID_BATCH_SIZE");
    auto use_debug_grid = options_.get_bool("SNLINK_USE_DEBUG_GRID");

    auto gauxc_grid = !use_debug_grid ? GauXC::MolGridFactory::create_default_molgrid(
        gauxc_mol, 
        pruning_scheme_map[pruning_scheme_],
        GauXC::BatchSize(grid_batch_size), 
        radial_scheme_map[radial_scheme_], 
        GauXC::RadialSize(radial_points_),
        GauXC::AngularSize(spherical_points_)
    ) :
    GauXC::MolGridFactory::create_default_molgrid(
        gauxc_mol, 
        pruning_scheme_map[pruning_scheme_],
        GauXC::BatchSize(grid_batch_size), 
        radial_scheme_map[radial_scheme_], 
        GauXC::AtomicGridSizeDefault::UltraFineGrid
    );
   
    // construct load balancer
    //const size_t quad_pad_value = 1;
    //auto gauxc_load_balancer = gauxc_load_balancer_factory_->get_instance(*rt, gauxc_mol, gauxc_grid, gauxc_primary, quad_pad_value);
    auto gauxc_load_balancer = gauxc_load_balancer_factory_->get_instance(*rt, gauxc_mol, gauxc_grid, gauxc_primary);
    
    // construct weights module
    auto gauxc_mol_weights = gauxc_mol_weights_factory_->get_instance();

    // apply partition weights
    gauxc_mol_weights.modify_weights(gauxc_load_balancer);

    // construct dummy functional
    // note that the snLinK used by the eventually-called eval_exx function is agnostic to the input functional!
    const std::string dummy_func_str = "B3LYP";
    auto spin = options_.get_str("REFERENCE") == "UHF" ? ExchCXX::Spin::Polarized : ExchCXX::Spin::Unpolarized;
    
    ExchCXX::Functional dummy_func_key = ExchCXX::functional_map.value(dummy_func_str);
    ExchCXX::XCFunctional dummy_func = { ExchCXX::Backend::builtin, dummy_func_key, spin };  

    // construct integrator
    integrator_ = integrator_factory_->get_shared_instance(dummy_func, gauxc_load_balancer); 
 
    timer_off("snLinK: Construct GauXC Integrator");
    
    timer_off("snLinK: Setup");
}

snLinK::~snLinK() {}

void snLinK::print_header() const {
    if (print_) {
        outfile->Printf("\n");
        outfile->Printf("  ==> snLinK: GauXC Semi-Numerical Linear Exchange K <==\n\n");
        
        outfile->Printf("    K Execution Space: %s\n", (use_gpu_) ? "Device" : "Host");
        outfile->Printf("    K Grid Radial Points: %i\n", radial_points_);
        outfile->Printf("    K Grid Spherical/Angular Points: %i\n\n", spherical_points_);
        outfile->Printf("    K Screening?:     %s\n", (integrator_settings_.screen_ek) ? "Yes" : "No");
        outfile->Printf("    K Ints Cutoff:     %11.0E\n", integrator_settings_.energy_tol);
        outfile->Printf("    K Basis Cutoff:     %11.0E\n", basis_tol_);
        if (debug_) {
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

    // we need to know if we are using a spherical harmonic basis
    // much of the behavior here is influenced by this
    auto is_spherical_basis = primary_->has_puream();
     
    // compute K for density Di using GauXC
    for (int iD = 0; iD != D.size(); ++iD) {
        timer_on("snLinK: Transform D");
        outfile->Printf("snLinK: Transform D\n");

        auto Did_eigen = D[iD]->eigen_map();    
        auto Kid_eigen = K[iD]->eigen_map();    
       
        // map Psi4 density matrix to Eigen matrix map
        SharedMatrix D_buffer = nullptr; 
        if (is_spherical_basis) {
            // need to reorder Psi4 density matrix to CCA ordering if in spherical harmonics
            if (permutation_matrix_.has_value()) {
                auto permutation_matrix_val = permutation_matrix_.value();
                auto D_eigen_permute = D[iD]->eigen_map();
                D_eigen_permute = permutation_matrix_val * D_eigen_permute * permutation_matrix_val.transpose();
            }
            
            // also need to transform D to cartesian coordinates if requested/required
            if (sph_to_cart_matrix_ != nullptr) {
                D_buffer = std::make_shared<Matrix>();
                D_buffer->transform(D[iD], sph_to_cart_matrix_);
            } else {
                D_buffer = D[iD];
            }
        // natively-cartesian bases dont require the reordering/transforming above
        } else {
            D_buffer = D[iD];
        }
        auto D_buffer_eigen = D_buffer->eigen_map();    
        
        timer_off("snLinK: Transform D");
        
        // map Psi4 exchange matrix buffer to Eigen matrix map
        timer_on("snLinK: Transform K");
        outfile->Printf("snLinK: Transform K\n");
        
        SharedMatrix K_buffer = nullptr; 
        if (is_spherical_basis) {
            // need to reorder Psi4 exchange matrix to CCA ordering if in spherical harmonics
            if (permutation_matrix_.has_value()) { 
                auto permutation_matrix_val = permutation_matrix_.value();
                auto K_eigen_permute = K[iD]->eigen_map();
                K_eigen_permute = permutation_matrix_val * K_eigen_permute * permutation_matrix_val.transpose();
            }
            
            // also need to transform D to cartesian coordinates if requested/required
            if (sph_to_cart_matrix_ != nullptr) { 
                K_buffer = std::make_shared<Matrix>();
                K_buffer->transform(K[iD], sph_to_cart_matrix_);
            } else {
                K_buffer = K[iD];
            }
        // natively-cartesian bases dont require the reordering/transforming above
        } else {
            K_buffer = K[iD];
        }
        auto K_buffer_eigen = K_buffer->eigen_map(); 
        timer_off("snLinK: Transform K");
        
        timer_on("snLinK: Execute integrator");
        outfile->Printf("snLinK: Execute integrator\n");
        // compute delta K if incfock iteration... 
        if (incfock_iter_) { 
            // K buffer (i.e., delta K) must be back-transformed and added to Psi4 K separately if cartesian transformation is forced...
            if ((sph_to_cart_matrix_ != nullptr) && is_spherical_basis) {
                K_buffer_eigen = integrator_->eval_exx(D_buffer_eigen, integrator_settings_);
                //K_buffer_eigen = K_buffer_eigen; 
                K_buffer->back_transform(sph_to_cart_matrix_);
                K[iD]->add(K_buffer);
            // ... otherwise the computation and addition can be bundled together 
            } else {
                K_buffer_eigen += integrator_->eval_exx(D_buffer_eigen, integrator_settings_);
                //K_buffer_eigen += K_buffer_eigen; 
            }

        // ... else compute full K 
        } else {
            K_buffer_eigen = integrator_->eval_exx(D_buffer_eigen, integrator_settings_);        
            //K_buffer_eigen = K_buffer_eigen; 
            if ((sph_to_cart_matrix_ != nullptr) && is_spherical_basis) {
                K[iD]->back_transform(K_buffer, sph_to_cart_matrix_);          
            }
        }
        timer_off("snLinK: Execute integrator");

        // now we need to reverse the CCA reordering previously performed
        timer_on("snLinK: Back-transform D and K");
        outfile->Printf("snLinK: Back-transform D and K\n");
        if (permutation_matrix_.has_value()) {
            auto permutation_matrix_val = permutation_matrix_ .value();
            auto D_eigen_permute = D[iD]->eigen_map();
            D_eigen_permute = permutation_matrix_val.transpose() * D_eigen_permute * permutation_matrix_val; 

            auto K_eigen_permute = K[iD]->eigen_map();
            K_eigen_permute = permutation_matrix_val.transpose() * K_eigen_permute * permutation_matrix_val; 
        }
        timer_off("snLinK: Back-transform D and K");
    }

    return;
}

}  // namespace psi
