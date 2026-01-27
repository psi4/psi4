#include "cghf.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <functional>
#include <tuple>
#include <utility>
#include <vector>
#include <any>

#include "psi4/physconst.h"

#include "psi4/libciomr/libciomr.h"
#include "psi4/libfock/jk.h"
#include "psi4/libfock/v.h"
#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/factory.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/pointgrp.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#define PSI4_NO_TIMER
#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libqt/qt.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libmints/sobasis.h"
#include "psi4/libmints/pointgrp.h"

//#include "Einsums/Tensor.hpp"
#include "Einsums/LinearAlgebra.hpp"
#include "Einsums/TensorAlgebra.hpp"
#include <Einsums/Config.hpp>


namespace psi {
namespace scf {

//using std::vector;
//using std::complex;

CGHF::CGHF(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> func)
    : HF(ref_wfn, func, Process::environment.options, PSIO::shared_object()) {
    common_init();
}

CGHF::CGHF(SharedWavefunction ref_wfn,
	   std::shared_ptr<SuperFunctional> func,
	   Options& options,
           std::shared_ptr<PSIO> psio)
    : HF(ref_wfn, func, options, psio) {
    common_init();
}

CGHF::~CGHF() {}

/*
 * 
 */



// Define global einsums::BlockTensors and variables needed throughout the CGHF class
void CGHF::common_init() {
    // => HF REQUIREMENTS from hf.cc <=
    Ca_ = SharedMatrix(factory_->create_matrix("alpha MO coefficients (C)"));
    Cb_ = SharedMatrix(factory_->create_matrix("beta MO coefficients (C)"));
    Fa_ = SharedMatrix(factory_->create_matrix("F alpha"));
    Fb_ = SharedMatrix(factory_->create_matrix("F beta"));
    Da_ = SharedMatrix(factory_->create_matrix("SCF alpha density"));
    Db_ = SharedMatrix(factory_->create_matrix("SCF beta density"));
    Va_ = SharedMatrix(factory_->create_matrix("G alpha"));
    Vb_ = SharedMatrix(factory_->create_matrix("G beta"));

    epsilon_a_ = SharedVector(factory_->create_vector());
    epsilon_a_->set_name("alpha orbital energies");
    epsilon_b_ = SharedVector(factory_->create_vector());
    epsilon_b_->set_name("beta orbital energies");

    same_a_b_dens_ = false;
    same_a_b_orbs_ = false;
    nuclearrep_ = molecule_->nuclear_repulsion_energy({0.0, 0.0, 0.0});

    // => DIIS variables <=
    err_vecs = std::deque<einsums::BlockTensor<std::complex<double>, 2>>(0);
    Fdiis = std::deque<einsums::BlockTensor<std::complex<double>, 2>>(0);
    diis_coeffs = std::vector<std::complex<double>>(0);
    error_doubles = std::vector<std::complex<double>>(0);

    find_occupation(); // Fills size of epsilon_a_ and epsilon_b_ so we can compute nalphapi_ + nbetapi_ = nelecpi_

    form_H(); // Fills a single spin block (AA or BB) with the core Hamiltonian 
    subclass_init(); // This appears to set up DFT stuff and external potentials

    // Combine nalphapi_ and nbetapi_ to give nelecpi_ for forming the density matrix
    for (int h = 0; h < nirrep_; h++) { nelecpi_.push_back(nalphapi_[h] + nbetapi_[h]); };

    // => GLOBAL VARS <=
    // irrep_sizes_: number of spin orbitals per irrep. Cannot use nsopi_ here since Einsums requires a vector
    // G_mat: two-electron ERI as a SharedMatrix
    // nelecpi_: nalphapi_ + nbetapi_
    G_mat = mintshelper()->ao_eri();

    // Double size of nsopi_ for GHF
    for (int h = 0; h < nirrep_; h++) { irrep_sizes_.push_back(2*nsopi_[h]); };

    // Combine nalphapi_ and nbetapi_ to give nelecpi_ for forming the density matrix
    
    // => BLOCKTENSORS <= will be (2nsopi_ x 2nsopi_) unless listed otherwise
    // Note that all of these are complex except for EINS_, EINX_, and the Fock evals
    F0_ = einsums::BlockTensor<std::complex<double>, 2>("F0", irrep_sizes_); // Spin blocked core Hamiltonian
    F_ = einsums::BlockTensor<std::complex<double>, 2>("F", irrep_sizes_);  // Fock matrix
    EINS_ = einsums::BlockTensor<double, 2>("Ovlp", irrep_sizes_); // Spin blocked overlap matrix. Needed as metric to compute gradient with FDS - SDF
    EINX_ = einsums::BlockTensor<std::complex<double>, 2>("EINX_", irrep_sizes_); // Spin blocked orthogonalization matrix
    Fp_ = einsums::BlockTensor<std::complex<double>, 2>("EINX_", irrep_sizes_); // Orthogonalized Fock matrix
    C_ = einsums::BlockTensor<std::complex<double>, 2>("C", irrep_sizes_);  // Coefficient matrix
    FDSmSDF_ = einsums::BlockTensor<std::complex<double>, 2>("C", irrep_sizes_);  // Gradient FDS - SDF

    Fevecs_ = einsums::BlockTensor<std::complex<double>, 2>("Fock eigenvectors", irrep_sizes_);  // Eigenvectors upon diagonalizing Fock matrix
    //Fevals_ = einsums::BlockTensor<double, 1>("Fock evals", irrep_sizes_); // Eigenvalues upon diagonalizing Fock matrix
    Fevals_ = einsums::BlockTensor<std::complex<double>, 1>("Fock evals", irrep_sizes_); // Eigenvalues upon diagonalizing Fock matrix

    Cocc_ = einsums::BlockTensor<std::complex<double>, 2>("Cocc_", irrep_sizes_); // Occupied coefficients matrix
    cCocc_ = einsums::BlockTensor<std::complex<double>, 2>("Cocc_", irrep_sizes_); // Complex conjugate of occupied coefficients matrix
    D_ = einsums::BlockTensor<std::complex<double>, 2>("D", irrep_sizes_);  // Density matrix
    J_ = einsums::BlockTensor<std::complex<double>, 2>("J", irrep_sizes_);  // Coulombic matrix
    K_ = einsums::BlockTensor<std::complex<double>, 2>("K", irrep_sizes_);  // Exchange matrix (might not need to be decomposed for efficiency reasons)
    temp1_ = einsums::BlockTensor<std::complex<double>, 2>("temp1", irrep_sizes_);  // Temporary storage container
    temp2_ = einsums::BlockTensor<std::complex<double>, 2>("temp2", irrep_sizes_);  // Temporary storage container

    F0_.zero();
    F_.zero();
    EINS_.zero();
    EINX_.zero();
    C_.zero();
    Fevecs_.zero();
    Fevals_.zero();
    Cocc_.zero();
    cCocc_.zero();
    D_.zero();
    J_.zero();
    K_.zero();
    temp1_.zero();
    temp2_.zero();
    FDSmSDF_.zero();

    preiterations(); // CGHF preiterations() routine
}


// Needed for initializing the Einsums spin-blocked overlap matrix EINS_
// and subsequently the orthogonalization matrix EINX_
// Also builds the spin-blocked Hamiltonian as F0_
void CGHF::preiterations() {
    // => FILL SPIN BLOCKED OVERLAP AND HAMILTONIAN MATRICES <=
    //
    // Note that nsopi_[h] will be HALF of irrep_sizes_[h]
    // EINS_ and F0_ filled within the same loop
    for (int h = 0; h < nirrep_; h++) 
    for (int p = 0; p < nsopi_[h]; p++)
    for (int q = 0; q < nsopi_[h]; q++) {
	// Offset by nsopi_[h] to get the beta index for p and q
        int beta_p = p+nsopi_[h];
        int beta_q = q+nsopi_[h];

        EINS_[h].subscript(p, q) = S_->get(h, p, q);           // overlap AA spin block
        EINS_[h].subscript(beta_p, beta_q) = S_->get(h, p, q); // overlap BB spin block is same as AA since same spatial orbitals
        
	F0_[h].subscript(p, q) = H_->get(h, p, q);           // core Hamiltonian AA spin block
	F0_[h].subscript(beta_p, beta_q) = H_->get(h, p, q); // core Hamiltonian BB spin block

    }

    // Constructs orthogonalization matrix X_ using Lowdin symmetric orthogonalization
    // 
    // X = S^{-1/2}
    // 
    // where eigenvalues less than the threshold 1e-7 are set to 0, and the eigenvectors removed
    auto X_real = einsums::linear_algebra::pow(EINS_, -0.5, 1e-7);

    // Cannot do the pow operation on a complex matrix for whatever reason, so 
    // simply transfer the real part to the complex matrix
    //
    // Apparently EINX_ could be real, but real/complex matrix-matrix multiplication
    // doesn't appear to work for me unless they're both the same type
    for (int i = 0; i < nirrep_; i++)
    for (int j = 0; j < irrep_sizes_[i]; j++)
    for (int k = 0; k < irrep_sizes_[i]; k++) {
            EINX_[i].subscript(j, k) = std::complex<double> (X_real[i](j, k), 0.0);
        }

    rotated_sad_guess(); // TODO modify the existing Psi4 SAD scheme (UHF, RHF, etc) since it causes errors otherwise
			 // e.g. the only way for this to work is by using a CORE guess

    //if (options_.get_str("GUESS") == "SAD") { rotated_sad_guess(); }; // Overwrites the D_ density matrix with the SAD initial guess
}

void CGHF::rotated_sad_guess() {
    // Experimental
    //
    // For each atom (including non-unique), takes the computed atomic density, and rotates about the y-axis with a given theta in degrees
    // yielding a spin symmetry-broken atomic density
    //
    // By itself, this would be completely meaningless (and extra costly), but combining them in a SAD fashion yields a nice non-collinear initial guess for the molecule
    //
    //
    int total_atoms = molecule_->natom();

    int start_idx = 0; // Where we will grab values from the molecular 1e Hamiltonian
    for (int atom_idx = 0; atom_idx < total_atoms; ++atom_idx) {
        int nshells = basisset_->nshell_on_center(atom_idx);
        int atom_nbf = 0;

        // Get the nelec for the atom by itself
	auto atom_nelec = molecule_->fZ(atom_idx); // Should subtract charge from this

        // Iterate through all shells on this center
        for (int shell_idx = 0; shell_idx < nshells; ++shell_idx) {
            // Get the specific GaussianShell object
            auto shell = basisset_->shell(atom_idx, shell_idx);
            atom_nbf += shell.nfunction(); // Add number of functions in this shell
        }


	// Now to effectively do a GHF calculation for each individual atom
        auto atom_F0_ = einsums::Tensor<std::complex<double>, 2>("F0", 2*atom_nbf, 2*atom_nbf);
        auto atom_F_ = einsums::Tensor<std::complex<double>, 2>("F", 2*atom_nbf, 2*atom_nbf);  
        auto atom_S_ = einsums::Tensor<double, 2>("Ovlp", 2*atom_nbf, 2*atom_nbf); 
        auto atom_EINX_ = einsums::Tensor<std::complex<double>, 2>("EINX_", 2*atom_nbf, 2*atom_nbf);
        auto atom_Fp_ = einsums::Tensor<std::complex<double>, 2>("EINX_", 2*atom_nbf, 2*atom_nbf); 
        auto atom_C_ = einsums::Tensor<std::complex<double>, 2>("C", 2*atom_nbf, 2*atom_nbf);
        auto atom_FDSmSDF_ = einsums::Tensor<std::complex<double>, 2>("C", 2*atom_nbf, 2*atom_nbf);  

        auto atom_Fevecs_ = einsums::Tensor<std::complex<double>, 2>("Fock eigenvectors", 2*atom_nbf, 2*atom_nbf);
        auto atom_Fevals_ = einsums::Tensor<std::complex<double>, 1>("Fock evals", 2*atom_nbf);

        auto atom_Cocc_ = einsums::Tensor<std::complex<double>, 2>("Cocc_", 2*atom_nbf, 2*atom_nbf); 
        auto atom_cCocc_ = einsums::Tensor<std::complex<double>, 2>("Cocc_", 2*atom_nbf, 2*atom_nbf);
        auto atom_D_ = einsums::Tensor<std::complex<double>, 2>("D", 2*atom_nbf, 2*atom_nbf);
        auto atom_J_ = einsums::Tensor<std::complex<double>, 2>("J", 2*atom_nbf, 2*atom_nbf);
        auto atom_K_ = einsums::Tensor<std::complex<double>, 2>("K", 2*atom_nbf, 2*atom_nbf); 
        auto atom_temp1_ = einsums::Tensor<std::complex<double>, 2>("temp1", 2*atom_nbf, 2*atom_nbf);
        auto atom_temp2_ = einsums::Tensor<std::complex<double>, 2>("temp2", 2*atom_nbf, 2*atom_nbf);

        atom_F0_.zero();
	atom_F_.zero();
	atom_S_.zero();
	atom_EINX_.zero();
	atom_Fp_.zero();
	atom_C_.zero();
	atom_FDSmSDF_.zero();
	atom_Fevecs_.zero();
	atom_Fevals_.zero();
	atom_Cocc_.zero();
	atom_cCocc_.zero();
	atom_D_.zero();
       
	atom_temp1_.zero();
	atom_temp2_.zero();

	auto F0_mol = F0_[0];
	auto S_mol = EINS_[0];

	// NOTE: we are stealing values from the core Hamiltonian for the entire molecule
	// I believe this should be legal since it's 1e, non-interacting, and in the AO basis

	// Fill AA and BB spin blocks of atom_F0_
	for (int p = 0; p < atom_nbf; p++)
	for (int q = 0; q < atom_nbf; q++) {
	    auto pshift = p + start_idx; // Shifted to map properly to the molecular 1e Hamiltonian
	    auto qshift = q + start_idx; // These shifts ONLY WORK for C1 symmetry

	    auto F0_AA = F0_mol.subscript(pshift, qshift);
            auto F0_BB = F0_mol.subscript(pshift+nsopi_[0], qshift+nsopi_[0]);
            auto F0_BA = F0_mol.subscript(pshift+nsopi_[0], qshift);
            auto F0_AB = F0_mol.subscript(pshift, qshift+nsopi_[0]);

            auto S_AA = S_mol.subscript(pshift, qshift);
            auto S_BB = S_mol.subscript(pshift+nsopi_[0], qshift+nsopi_[0]);
            auto S_BA = S_mol.subscript(pshift+nsopi_[0], qshift);
            auto S_AB = S_mol.subscript(pshift, qshift+nsopi_[0]);

	    atom_F0_.subscript(p, q) = F0_AA;
            atom_F0_.subscript(p+atom_nbf, q+atom_nbf) = F0_BB;
            atom_F0_.subscript(p+atom_nbf, q) = F0_BA;
            atom_F0_.subscript(p, q+atom_nbf) = F0_AB;

            atom_S_.subscript(p, q) = S_AA;
            atom_S_.subscript(p+atom_nbf, q+atom_nbf) = S_BB;
            atom_S_.subscript(p+atom_nbf, q) = S_BA;
            atom_S_.subscript(p, q+atom_nbf) = S_AB;

	}
	start_idx += atom_nbf;

    	auto X_real = einsums::linear_algebra::pow(atom_S_, -0.5, 1e-7);

        for (int j = 0; j < 2*atom_nbf; j++)
        for (int k = 0; k < 2*atom_nbf; k++) {
                atom_EINX_.subscript(j, k) = std::complex<double> (X_real(j, k), 0.0);
            }


	atom_F_ = atom_F0_; // Use Core Hamiltonian as first Fock matrix

	// Orthogonalize the Fock matrix
        einsums::linear_algebra::gemm<false, false>(std::complex<double>{1.0}, atom_F_, atom_EINX_, std::complex<double>{0.0},  &atom_temp1_);
        einsums::linear_algebra::gemm<true, false>(std::complex<double>{1.0}, atom_EINX_, atom_temp1_, std::complex<double>{0.0}, &atom_Fp_);

	// Diagonalize Fock matrix and sort eigenpairs	
        einsums::linear_algebra::geev(&atom_Fp_, &atom_Fevals_, &atom_temp1_, &atom_temp2_); // atom_temp1_ and atom_temp2_ are left and right eigenvectors, respectively
        sort_eigenpairs(atom_Fevals_, atom_temp2_, 1e-10);
        atom_Fp_ = atom_temp2_;

	// Back transform to get the coefficients
        einsums::linear_algebra::gemm<false, false>(std::complex<double>{1.0}, atom_EINX_, atom_Fp_, std::complex<double>{0.0}, &atom_C_);

        //Do alpha-beta mixing on the occupied coefficients BEFORE using to make a density matrix
	auto eps = 0.01; // Mixing parameter

        println(atom_C_);	
	for (int m = 0; m < atom_nbf; m++)
	for (int i = 0; i < atom_nelec; i++) {
	    auto alpha_coeff = atom_C_.subscript(m, i);
	    auto beta_coeff = atom_C_.subscript(m+atom_nbf, i);

	    atom_C_.subscript(m+atom_nbf, i) = beta_coeff + (eps * alpha_coeff); // Mixing happens here
	}

	// Loop through and grab the OCCUPIED coefficients and its conjugate matrix
        for (int j = 0; j < 2*atom_nbf; j++)
        for (int k = 0; k < atom_nelec; k++) {
                atom_Cocc_(j,k) = atom_C_(j, k);
                atom_cCocc_(j,k) = std::conj(atom_C_(j, k));
            }

        // Performs einsums contraction ui,vi->uv with (Cocc_, cCocc_ -> D_)
        einsums::tensor_algebra::einsum(
            einsums::Indices{einsums::index::u,
                             einsums::index::v}, // D_ with uv
            &atom_D_,
            einsums::Indices{einsums::index::u,
                             einsums::index::i}, // Cocc_ with ui
            atom_Cocc_,
            einsums::Indices{einsums::index::v,
                             einsums::index::i}, // cCocc_ with vi
            atom_cCocc_);

	//// Get Sx, Sy, and Sz of atomic density matrix so we can know how much to rotate
	//auto D_AA = einsums::Tensor<std::complex<double>, 2>("D_AA", atom_nbf, atom_nbf);
        //auto D_BB = einsums::Tensor<std::complex<double>, 2>("D_BB", atom_nbf, atom_nbf);
        //auto D_AB = einsums::Tensor<std::complex<double>, 2>("D_AB", atom_nbf, atom_nbf);
        //auto D_BA = einsums::Tensor<std::complex<double>, 2>("D_BA", atom_nbf, atom_nbf);
        //auto S_block = einsums::Tensor<std::complex<double>, 2>("Overlap matrix spin block", atom_nbf, atom_nbf);

        //auto Cont = einsums::Tensor<std::complex<double>, 2>("S container", atom_nbf, atom_nbf);

        //for (int p = 0; p < atom_nbf; p++)
	//for (int q = 0; q < atom_nbf; q++) {
	//    D_AA.subscript(p, q) = atom_D_.subscript(p, q);
        //    D_BB.subscript(p, q) = atom_D_.subscript(p+atom_nbf, q+atom_nbf);
        //    D_AB.subscript(p, q) = atom_D_.subscript(p, q+atom_nbf);
        //    D_BA.subscript(p, q) = atom_D_.subscript(p+atom_nbf, q);
	//    Cont.subscript(p, q) = atom_S_.subscript(p, q); // Since all spin-blocks will use the same overlap block
	//}	

	//// (D_AA - D_BB) @ S
        //auto AABB = D_AA;
	//AABB -= D_BB;

        //einsums::linear_algebra::gemm<false, false>(std::complex<double>{1.0}, AABB, S_block, std::complex<double>{0.0}, &Cont);

	//// Sz = 0.5 * Tr( (AA - BB) @ S).real
	//
	//// Get traces of AABB first
	//double AABB_trace = 0.0;

	//for (int p = 0; p < atom_nbf; p++) {
	//    AABB_trace += Cont.subscript(p, p).real();
	//}

	//auto Sz = 0.5 * AABB_trace;


	//// Sx = 0.5 * Tr( (AB - BA) @ S).real
	//
        //auto ABBA = D_AB;
        //ABBA -= D_BA;
        //
	//Cont.zero();
        //einsums::linear_algebra::gemm<false, false>(std::complex<double>{1.0}, ABBA, S_block, std::complex<double>{0.0}, &Cont);	

        //double ABBA_trace = 0.0;
        //for (int p = 0; p < atom_nbf; p++) {
        //    ABBA_trace += Cont.subscript(p, p).real();
        //}

	//auto Sx = 0.5 * ABBA_trace;

        //// Sy = 0.5 * Tr( (-1j * AB + 1j * BA) @ S).real
	//std::complex<double> IMAG = std::complex<double>(0.0, 1.0);
	//auto cAB = D_AB;
	//cAB *= -IMAG;

	//auto cBA = D_BA;
	//cBA *= IMAG;

	//auto cABBA = cAB;
	//cABBA -= cBA;

        //Cont.zero();
        //einsums::linear_algebra::gemm<false, false>(std::complex<double>{1.0}, cABBA, S_block, std::complex<double>{0.0}, &Cont);

	//double cABBA_trace = 0.0;
	//for (int p = 0; p < atom_nbf; p++) {
	//    cABBA_trace += Cont.subscript(p, p).real();
	//}

	//auto Sy = 0.5 * cABBA_trace;

	//std::cout << "Spin magnitudes\n";
	//std::cout << "Sz = " << Sz << "\n";
        //std::cout << "Sy = " << Sy << "\n";
        //std::cout << "Sx = " << Sx << "\n";

    }

}

void CGHF::finalize() {
}

void CGHF::save_density_and_energy() {
}

void CGHF::form_V() {
}


/*
* Does an explicit 4-index for loop to compute JKwK_
*
* The density matrix is decomposed into its 4 spin-blocks and evaluated separately.
* According to the Slater-Condon rules, the off-diagonal elements (alpha-beta, beta-alpha)
* have zero Coulomb energy contributions. Therefore, only J_aa and J_bb need to be evaluated
*
* There is only one symmetry to take advantage of: K_ab = -K_ba† (negative adjoints)
* Thus, 5 terms are computed here:
*
* J_aa/J_bb =      D_sr, G_pqrs -> J_pq
*
* K_aa/K_bb/K_ab = D_sr, G_psrq -> K_pq
*
* which are of course computed using their respective density spin block (e.g. K_aa is D_aa)
*
*/

void CGHF::form_G() {
    J_.zero();
    K_.zero();

    int nso = nso_;
    int dim = 2*nso;

    for (int h = 0; h < nirrep_; h++)
    for (int p = 0; p < nsopi_[h]; p++)
    for (int q = p; q < nsopi_[h]; q++) { //Note the q = p
        for (int s = 0; s < nsopi_[h]; s++) {
            for (int r = 0; r < nsopi_[h]; r++) {
                auto g_pqrs = G_mat->get(h, p*nsopi_[h] + q, r*nsopi_[h] + s);
                auto g_psrq = G_mat->get(h, p*nsopi_[h] + s, r*nsopi_[h] + q);
                //int idx_pqrs  = ((p*nso + q)*nso + r)*nso + s;
                //int idx_psrq  = ((p*nso + s)*nso + r)*nso + q;
                //auto g_pqrs = Gp[idx_pqrs];
                //auto g_psrq = Gp[idx_psrq];
                if (std::abs(g_pqrs) < 1e-12 && std::abs(g_psrq) < 1e-12) { continue; };

                auto D_AA = D_[h].subscript(s, r);
                auto D_BB = D_[h].subscript(s+nsopi_[h], r+nsopi_[h]);
		auto D_BA = D_[h].subscript(s+nsopi_[h], r);

                auto sum_D = D_AA + D_BB;

                J_[h].subscript(p, q) += g_pqrs * sum_D - g_psrq * D_AA;
                J_[h].subscript(p+nsopi_[h], q+nsopi_[h]) += g_pqrs * sum_D - g_psrq * D_BB;
            }
        }
    }

    //Then fill in lower triangle
    for (int h = 0; h < nirrep_; h++)
    for (int i = 0; i < irrep_sizes_[h]; i++) {
        for (int j = 0; j < i; j++) {
            J_[h].subscript(i, j) = std::conj(J_[h].subscript(j, i));
        }
    }



}



// Combines F = H + J - K
void CGHF::form_F() {
    F_.zero();
    F_ += F0_;
    F_ += J_;
    F_ -= K_;
}

void CGHF::sort_eigenpairs(
    einsums::Tensor<std::complex<double>, 1>& evals,
    einsums::Tensor<std::complex<double>, 2>& evecs,
    double tol)
{
    // Assumes evecs is square: (dim, dim)
    int dim = evecs.dims()[0];

    // Index map: unsorted -> sorted
    std::vector<int> idx(dim);
    for (int i = 0; i < dim; i++) idx[i] = i;

    // Extract eigenvalues for sorting
    std::vector<double> vevals(dim);
    for (int i = 0; i < dim; i++)
        vevals[i] = std::real(evals(i));

    std::sort(idx.begin(), idx.end(),
              [&](int a, int b) { return vevals[a] < vevals[b]; });

    // Output containers
    einsums::Tensor<std::complex<double>, 1> sorted_evals("Sorted evals", dim);
    einsums::Tensor<std::complex<double>, 2> sorted_evecs("Sorted evecs", dim, dim);
    sorted_evecs.zero();

    int i = 0;
    while (i < dim) {

        // Identify degenerate block
        int j = i;
        while (j + 1 < dim &&
               std::abs(vevals[idx[j + 1]] - vevals[idx[i]]) < tol)
        {
            j++;
        }

        int block_size = j - i + 1;

        // Copy eigenvalues + eigenvectors for this block
        for (int k = 0; k < block_size; k++) {
            int src_col = idx[i + k];

            sorted_evals(i + k) = evals(src_col);

            for (int row = 0; row < dim; row++)
                sorted_evecs(row, i + k) = evecs(row, src_col);
        }

        // (Gram–Schmidt will go here later)

        i += block_size;
    }

    // Copy back
    evals = sorted_evals;
    evecs = sorted_evecs;
}

/*
 * Orthogonalizes the Fock matrix via Fp_ = X_dagger * F * X_
 * Then diagonalizes Fp_ to give a set of eigenpairs in the MO basis
 * e.g. the eigenvalues are the MO energies
 * Finally, the eigenvectors are back-transformed back to the AO basis
 *
 * Ignores the shift for now
 *
*/
void CGHF::form_C(double shift) {
    // => ORTHOGONALIZE FOCK <=
    // F_ @ EINX_ = temp1_
    //

    if (options_.get_bool("DIIS") && iteration_ > 1) { 
        auto error_trace = do_diis();
    } else {
        einsums::linear_algebra::gemm<false, false>(std::complex<double>{1.0}, F_, EINX_, std::complex<double>{0.0},  &temp1_);
   
        //EINX_ @ temp1 = Fp_
        einsums::linear_algebra::gemm<true, false>(std::complex<double>{1.0}, EINX_, temp1_, std::complex<double>{0.0}, &Fp_);
    }
    // => DIAGONALIZE FOCK <=
    //Fevecs_ = Fp_;

    temp1_.zero();
    temp2_.zero();
    // Must be done in a loop over irreps -- cannot just do it with the BlockTensor
    for (int h = 0; h < nirrep_; h++) {
	if (nsopi_[h] == 0) { continue ; }; // Cannot diagonalize a 0x0 matrix
        //einsums::linear_algebra::heev<true>(&Fp_[h], &Fevals_[h]); // Hermitian eigensolver SHOULD be legitimate over general eigensolver
	einsums::linear_algebra::geev(&Fp_[h], &Fevals_[h], &temp1_[h], &temp2_[h]); // Where temp1 and temp2 are the left and right eigenvectors, respectively
        sort_eigenpairs(Fevals_[h], temp2_[h], 1e-10); // Overwrites Fevals_ and temp2_ with the sorted eigenpairs
	Fp_[h] = temp2_[h];
    }

    // => BACK TRANSFORM <=
    // EINX_ @ Fevecs_ = C_
    einsums::linear_algebra::gemm<false, false>(std::complex<double>{1.0}, EINX_, Fp_, std::complex<double>{0.0}, &C_);
}

/*
 * Fills the occupied coefficient matrix Cocc_ as well as the complex conjugate cCocc_
 * Then constructs the density matrix D_ with D_uv = C_ui * C_vi_{conj}
 */
void CGHF::form_D() {
    Cocc_.zero();
    cCocc_.zero();

    // Simply fills Cocc_ and cCocc_ with the (2*nsopi_ x nelecpi_) matrix
    // Note that Cocc_ and cCocc_ are both (2*nsopi_ x 2*nsopi_), but the 'remainder'
    // are all zeros and contribute nothing
    // This is just a minor inefficiency, since these matrices are larger than they should be

    for (int i = 0; i < nirrep_; i++)
    for (int j = 0; j < irrep_sizes_[i]; j++)
    for (int k = 0; k < nelecpi_[i]; k++) {
            Cocc_[i](j,k) = C_[i](j, k);
            cCocc_[i](j,k) = std::conj(C_[i](j, k));
        }

    // Performs einsums contraction ui,vi->uv with (Cocc_, cCocc_, D_)
    einsums::tensor_algebra::einsum(
        einsums::Indices{einsums::index::u,
                         einsums::index::v}, // D_ with uv
        &D_,
        einsums::Indices{einsums::index::u,
                         einsums::index::i}, // Cocc_ with ui
        Cocc_,
        einsums::Indices{einsums::index::v,
                         einsums::index::i}, // cCocc_ with vi
        cCocc_);

}

void CGHF::damping_update(double damping_percentage) {
}


/*
 * E = E_1e + E_2e + E_nuc
 *
 * E_1e = trace(D_ • F0_)
 * E_2e = trace(D_ • 0.5*JK_)
 */
double CGHF::compute_E() {
    double kinetic_E = 0.0;
    double one_electron_E = 0.0;
    double two_E = 0.0;
    double exchange_E = 0.0;

    // Because Psi4 likes these energies (me too) decomposed, it seems both easier and more
    // efficient to just knock them out in the same for loop
    for (int h = 0; h < nirrep_; h++) {
	if (nsopi_[h] == 0) { continue; };
        for (int i = 0; i < irrep_sizes_[h]; i++)
        for (int j = 0; j < irrep_sizes_[h]; j++) {
            auto Dji = D_[h].subscript(j, i);

            //kinetic_E += (T_->get(h, i, j) * Dji).real();            // T_ij * Dji
            one_electron_E += (F0_[h].subscript(i, j) * Dji).real(); // F0_ij * D_ji
            two_E += (J_[h].subscript(i, j) * Dji).real();       // J_ij * D_ji
            //exchange_E += (K_[h].subscript(i, j) * Dji).real();      // K_ij * D_ji
        }
    }

    // Add these to the global energies struct
    // These are seen in the output.dat file
    energies_["Nuclear"] = nuclearrep_;
    energies_["Kinetic"] = kinetic_E;
    energies_["One-Electron"] = one_electron_E;
    energies_["Two-Electron"] = 0.5 * two_E;

    double Etotal = 0.0;
    Etotal += nuclearrep_;
    Etotal += one_electron_E;
    Etotal += 0.5 * two_E;
    return Etotal;
}


/*
 * Simply computes the 1e energy 
 */
double CGHF::compute_initial_E() {
    double one_electron_E = 0.0;

    for (int h = 0; h < nirrep_; h++)
    for (int i = 0; i < irrep_sizes_[h]; i++)
    for (int j = 0; j < irrep_sizes_[h]; j++) {
        auto Dji = D_[h].subscript(j, i);
        one_electron_E += (F0_[h].subscript(i, j) * Dji).real(); // F0_ij * D_ji

    }
    
    return one_electron_E;
}

std::complex<double> CGHF::do_diis() {
    int diis_max = 8;
    int diis_count = 0;

    //FDS-SDF
    form_FDSmSDF();

    auto ortho_error = einsums::BlockTensor<std::complex<double>, 2>("Orthogonalized FDSmSDF", irrep_sizes_); //eorth_it
    auto temp1 = einsums::BlockTensor<std::complex<double>, 2>("Temp FDSmSDF", irrep_sizes_);
    auto ecurr = einsums::BlockTensor<std::complex<double>, 2>("Current error", irrep_sizes_);

    temp1_.zero();
    ortho_error.zero();
    ecurr.zero();
    temp2_.zero();

    einsums::linear_algebra::gemm<false, false>(std::complex<double>{1.0}, F_, EINX_, std::complex<double>{0.0},  &temp2_);
    einsums::linear_algebra::gemm<true, false>(std::complex<double>{1.0}, EINX_, temp2_, std::complex<double>{0.0}, &Fp_);

    einsums::linear_algebra::gemm<false, false>(std::complex<double>{1.0}, FDSmSDF_, EINX_, std::complex<double>{0.0},  &temp1_);
    einsums::linear_algebra::gemm<true, false>(std::complex<double>{1.0}, EINX_, temp1_, std::complex<double>{0.0},  &ortho_error);

    einsums::linear_algebra::gemm<true, false>(std::complex<double>{1.0}, ortho_error, ortho_error, std::complex<double>{0.0},  &ecurr);

    auto error_trace = std::complex<double> {0.0, 0.0};

    //Get collective trace of all irreps
    for (int i = 0; i < nirrep_; i++) {
        for (int j = 0; j < irrep_sizes_[i]; j++) {
                error_trace = ecurr[i](j, j);
        }
    }

    auto abs_trace = std::abs(error_trace);

    //If we have enough error vectors in storage, check for the worst one to remove
    if (err_vecs.size() == diis_max) {
        //Find the one with the max error to eliminate from storage
        //Give it some ludicrously low error so that anything is larger
        //C++ isn't like Python where we can assign the var in the loop
        //at least not at first
        auto max_error = std::complex<double> {-10000000000, -10000000000};
        int max_error_ind = 0;

        for (int i = 0; i < diis_max; i++) {
            auto curr_error = error_doubles.at(i);

            //You can't compare two complex numbers directly, so I compare the real parts
            //Doing it with imaginary is very tricky since these are very very close to 0
            //When I checked imag as well, we were converging at 17 iterations instead of 14
            if (curr_error.real() > max_error.real()) {
               max_error = curr_error;
               max_error_ind = i;
            }
        }
        /*
         * Use this to replace the Fock matrix with the most error
         * Seems to be more applicable to what I'm doing
         */
        //Then make the replacement based off the index
        Fdiis.at(max_error_ind) = Fp_;
        err_vecs.at(max_error_ind) = FDSmSDF_;

        //std::cout << "DIIS vector being replaced at index " << max_error_ind << "\n";
        //std::cout << "Current DIIS subspace size: " << Fdiis.size() << "\n";
        //Norm for the replacement error vector, which is why the index is different here compared to below
        auto norm = einsums::linear_algebra::dot(err_vecs.at(max_error_ind), err_vecs.at(max_error_ind));
        error_doubles.at(max_error_ind) = norm; //This is set to real for now, not sure what to do about this

        /*
         * Use this to reset the subspace at every nth iteration
         * Might be more useful for generally converging
         * energies rather than longer geometries
         *
        Fdiis.clear();
        err_vecs.clear();
        Fdiis.push_back(Fp_);
        err_vecs.push_back(error_);
        */
    }

    //If not, add the Fock matrix and error matrix to memory
    else {
        err_vecs.push_back(FDSmSDF_);
        Fdiis.push_back(Fp_);
        //Calculate the norm for the last error vector
        auto norm = einsums::linear_algebra::dot(err_vecs.at(err_vecs.size()-1), err_vecs.at(err_vecs.size()-1));
        error_doubles.push_back(norm); //See above comment with norm

    }

    //Next form the 'base' B matrix
    auto B_ = einsums::Tensor<std::complex<double>, 2>("Overlap matrix", err_vecs.size() + 1, err_vecs.size() + 1);
    B_.zero();

    //I do it in a weird way by filling the bottom row and
    //the right column with the 1.0+0i
    for (int i = 0; i < err_vecs.size()+1; i++) {
        B_(i, err_vecs.size()) = std::complex<double> {1.0, 0.0};
        B_(err_vecs.size(), i) = std::complex<double> {1.0, 0.0};
    }

    //Then simply zero out the last element in the matrix
    B_(err_vecs.size(), err_vecs.size()) = std::complex<double> {0.0, 0.0};

    //Actually populate it with the overlaps

    //for (int i = 0; i < err_vecs.size(); i++) {
    //    for (int j = 0; j <= i; j++) { //Had to consult Connor's UHF for this little trick as I tried 80 million things
    //        B_(i, j) = einsums::linear_algebra::dot(err_vecs[i], err_vecs[j]);
    //        B_(j, i) = B_(i, j); //Which makes perfect sense actually
    //    }
    //}
    for (int i = 0; i < err_vecs.size(); i++) {
        for (int j = 0; j < err_vecs.size(); j++) {
            B_(i, j) = {0.0,0.0};
            for (int p = 0; p < irrep_sizes_[0]; p++){
                for (int q = 0; q < irrep_sizes_[0]; q++){
                    B_(i, j) += std::conj(err_vecs[i](q,p)) * err_vecs[j](p,q);
                }
            }
        }
    }

    //Container for the coefficients after gesv
    //This is a weird artifact of gesv. This should be a column vector
    //But it doesn't accept the type? IDK
    auto C_temp = einsums::Tensor<std::complex<double>, 2>("Temp coefficient matrix", 1, err_vecs.size() + 1);

    C_temp.zero();

    //C_temp(0, err_vecs.size()+1) = std::complex<double> {1.0, 0.0};

    for (int i = 0; i < err_vecs.size()+1; i++) {
        C_temp(0, i) = std::complex<double> {1.0, 0.0};
    }

    einsums::linear_algebra::gesv(&B_, &C_temp);


    //Increase the size of the subspace
    diis_coeffs.resize(err_vecs.size());

    for (int i = 0; i < err_vecs.size(); i++) {
        diis_coeffs.at(i) = C_temp(0, i);
    }

    //Then we can extrapolate the Fock matrix!

    //println(Fdiis[0]);
    //
    Fp_.zero();

    if (err_vecs.size() == diis_max) {
        for (int i = 0; i < diis_coeffs.size(); i++) {
            einsums::linear_algebra::axpy(diis_coeffs[i], Fdiis[i], &Fp_);
            //std::cout << diis_coeffs[i] << "\n";
        }
    }
    else {
        temp2_.zero();
        einsums::linear_algebra::gemm<false, false>(std::complex<double>{1.0}, F_, EINX_, std::complex<double>{0.0},  &temp2_);
        einsums::linear_algebra::gemm<true, false>(std::complex<double>{1.0}, EINX_, temp2_, std::complex<double>{0.0}, &Fp_);
    }

    //Fp_ = temp2;

    //F0_ = temp2;
    //

    return error_trace;

}

void CGHF::form_FDSmSDF() {
    auto SComplex = einsums::BlockTensor<std::complex<double>, 2>("Complex overlap matrix", irrep_sizes_);
    SComplex.zero();

    for (int i = 0; i < nirrep_; i++) {
        for (int j = 0; j < irrep_sizes_[i]; j++) {
            for (int k = 0; k < irrep_sizes_[i]; k++) {
                SComplex[i].subscript(j, k) = EINS_[i].subscript(j, k);
            }
        }
    }

    for (int i = 0; i < nirrep_; i++) {
        for (int j = 0; j < irrep_sizes_[i]; j++) {
            for (int k = 0; k < irrep_sizes_[i]; k++) {
                FDSmSDF_[i](j,k) = 0.;
                for (int p = 0; p < irrep_sizes_[i]; p++) {
                    for (int q = 0; q < irrep_sizes_[i]; q++) {
                        FDSmSDF_[i](j,k) += F_[i](j,p) * D_[i](p,q) * SComplex[i](q,k);
                        FDSmSDF_[i](j,k) -= std::conj(F_[i](k,p) * D_[i](p,q) * SComplex[i](q,j));
                    }
                }
            }
        }
    }



}

double CGHF::compute_Dnorm() {
    double dnorm = einsums::linear_algebra::norm(einsums::linear_algebra::Norm::Frobenius, FDSmSDF_[0]);

    return dnorm;
}

void CGHF::setup_potential() {
}

void CGHF::openorbital_scf() {
}

int CGHF::soscf_update(double soscf_conv, int soscf_min_iter, int soscf_max_iter, int soscf_print) {
}

std::vector<SharedMatrix> CGHF::onel_Hx(std::vector<SharedMatrix> x) {
}

std::vector<SharedMatrix> CGHF::twoel_Hx(std::vector<SharedMatrix> x, bool combine,
                                   std::string return_basis) {
}

std::vector<SharedMatrix> CGHF::cphf_Hx(std::vector<SharedMatrix> x) {
}

std::vector<SharedMatrix> CGHF::cphf_solve(std::vector<SharedMatrix> x_vec, double conv_tol, int max_iter,
                                         int print_lvl) {
}

std::shared_ptr<CGHF> CGHF::c1_deep_copy(std::shared_ptr<BasisSet> basis) {
    auto wfn = Wavefunction::c1_deep_copy(basis);
    //auto hf_wfn = std::make_shared<CGHF>(wfn, functional_, wfn->options(), wfn->psio());
    
    auto hf_wfn = std::shared_ptr<CGHF>(new CGHF(wfn, functional_, wfn->options(), wfn->psio()));

    return hf_wfn;
}

}
}





