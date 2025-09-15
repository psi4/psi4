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
#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libqt/qt.h"
#include "psi4/libtrans/integraltransform.h"

#include "Einsums/Tensor.hpp"
#include "Einsums/LinearAlgebra.hpp"
#include "Einsums/TensorAlgebra.hpp"

namespace psi {
namespace scf {

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

auto GetRealVector(auto A, auto dim) {
    auto real = einsums::Tensor<double, 1>("R", dim);

    for (int i = 0; i < dim; i++) {
        real(i) = A(i).real();
    }

    return real;

}

void CGHF::common_init() {
    G_mat = mintshelper()->ao_eri(); 

    nelecpi_ = nalphapi_ + nbetapi_;
    for (int h = 0; h < nirrep_; h++) {
        irrep_sizes_.push_back(2*nsopi_[h]);
    }

    F0_ = einsums::BlockTensor<std::complex<double>, 2>("Core Fock", irrep_sizes_);
    EINT_ = einsums::BlockTensor<std::complex<double>, 2>("Core Fock", irrep_sizes_);
    F_ = einsums::BlockTensor<std::complex<double>, 2>("Core Fock", irrep_sizes_);
    Fp_ = einsums::BlockTensor<std::complex<double>, 2>("Core Fock", irrep_sizes_);
    EINX_ = einsums::BlockTensor<std::complex<double>, 2>("Core Fock", irrep_sizes_);
    EINS_ = einsums::BlockTensor<double, 2>("Core Fock", irrep_sizes_);
    C_ = einsums::BlockTensor<std::complex<double>, 2>("Core Fock", irrep_sizes_);
    Cocc_ = einsums::BlockTensor<std::complex<double>, 2>("Core Fock", irrep_sizes_);
    cCocc_ = einsums::BlockTensor<std::complex<double>, 2>("Core Fock", irrep_sizes_);
    FDSmSDF_ = einsums::BlockTensor<std::complex<double>, 2>("Orbital gradient", irrep_sizes_);

    D_ = einsums::BlockTensor<std::complex<double>, 2>("Core Fock", irrep_sizes_);
    Fevecs_ = einsums::BlockTensor<std::complex<double>, 2>("Core Fock", irrep_sizes_);
    Fevals_ = einsums::BlockTensor<std::complex<double>, 1>("Core Fock", irrep_sizes_);
    RealEvals_ = einsums::BlockTensor<double, 1>("Core Fock", irrep_sizes_);

    JKwK_ = einsums::BlockTensor<std::complex<double>, 2>("Core Fock", irrep_sizes_);
    temp = einsums::BlockTensor<std::complex<double>, 2>("Core Fock", irrep_sizes_);
    temp2 = einsums::BlockTensor<std::complex<double>, 2>("Core Fock", irrep_sizes_);
 
    F0_.zero();
    EINT_.zero();
    F_.zero();
    Fp_.zero();
    EINX_.zero();
    EINS_.zero();
    C_.zero();
    Cocc_.zero();
    cCocc_.zero();
    FDSmSDF_.zero();
    D_.zero();
    Fevecs_.zero();
    Fevals_.zero();
    RealEvals_.zero();
    JKwK_.zero();
    temp.zero();
    temp2.zero();

    //subclass_init();
    form_H();
    form_init_F();
    form_S();
    form_X();

    F_ += F0_;
    
    //Useless things I can't circumvent in the HF class
    //Perhaps I can overwrite some parent routines from HF
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
    //End useless things
    
    err_vecs = std::deque<einsums::BlockTensor<std::complex<double>, 2>>(0);
    Fdiis = std::deque<einsums::BlockTensor<std::complex<double>, 2>>(0);
    diis_coeffs = std::vector<std::complex<double>>(0);
    error_doubles = std::vector<std::complex<double>>(0);

    //preiterations();
}

void CGHF::form_S() {
    for (int i = 0; i < nirrep_; i++) {
        for (int j = 0; j < S_->rowdim(i); j++) {
                for (int k = 0; k < S_->coldim(i); k++) {
                        EINS_[i].subscript(j, k) = S_->get(i, j, k);
                        EINS_[i].subscript(j+S_->rowdim(i), k+S_->coldim(i)) = S_->get(i, j, k);
                }
        }

    }
}


void CGHF::form_X() {
    auto X_real = einsums::linear_algebra::pow(EINS_, -0.5, 1e-14);


    for (int i = 0; i < nirrep_; i++) {
        for (int j = 0; j < irrep_sizes_[i]; j++) {
            for (int k = 0; k < irrep_sizes_[i]; k++) {
                 EINX_[i].subscript(j, k) = X_real[i](j, k);
            }
        }
        //X_.push_block(X_block);
    }

}

void CGHF::form_init_F() {
    for (int h = 0; h < nirrep_; h++) {
        for (int j = 0; j < nsopi_[h]; j++) {
            for (int k = 0; k < nsopi_[h]; k++) {
                F0_[h].subscript(j, k) = H_->get(h, j, k);
                F0_[h].subscript(j+nsopi_[h], k+nsopi_[h]) = H_->get(h, j, k);

                EINT_[h].subscript(j, k) = T_->get(h, j, k);
	        EINT_[h].subscript(j+nsopi_[h], k+nsopi_[h]) = T_->get(h, j, k);	

            }
        }
    }
}

void CGHF::preiterations() {
    form_C(0.0);
    form_D();
    form_G();
    double E0 = compute_E();
}

void CGHF::finalize() {
}

void CGHF::save_density_and_energy() {
}

void CGHF::form_V() {
}

void CGHF::form_G() {
    JKwK_.zero();
    //This is for circumventing using the Psi4 JK object. Purely for proof that the JK hack works
    //Not recommended to use since no algorithms are used to boost efficiency (e.x. permutational symmetry)

    for (int h = 0; h < nirrep_; h++) {
        for (int p = 0; p < nsopi_[h]; p++) {
            for (int q = 0; q < nsopi_[h]; q++) {
                int pq = p*nsopi_[h] + q;
                for (int r = 0; r < nsopi_[h]; r++) {
                    for (int s = 0; s < nsopi_[h]; s++) {
                        int rs = r*nsopi_[h] + s;
                        int ps = p*nsopi_[h] + s;
                        int rq = r*nsopi_[h] + q;

			std::complex<double> D_aa  = D_[h].subscript(s, r);
                        std::complex<double> D_bb = D_[h].subscript(s+nsopi_[h], r+nsopi_[h]);
                        std::complex<double> D_ab = D_[h].subscript(s, r+nsopi_[h]);
                        //auto D_ba = D_[h].subscript(s+nsopi_[h], r);

                        double pqrs = G_mat->get(h, pq, rs);
                        double psrq = G_mat->get(h, ps, rq);
                        
			//Jaa gaaaa Daa
			auto J_val = pqrs*D_aa + pqrs*D_bb;

                        JKwK_[h].subscript(p, q) += J_val;
			JKwK_[h].subscript(p+nsopi_[h], q+nsopi_[h]) += J_val;

                        JKwK_[h].subscript(p, q) -= psrq*D_aa;

                        auto K_contract_ab = psrq*D_ab;
                        JKwK_[h].subscript(p, q+nsopi_[h]) -= K_contract_ab;
			JKwK_[h].subscript(p+nsopi_[h], p) += std::conj(K_contract_ab);

			JKwK_[h].subscript(p+nsopi_[h], q+nsopi_[h]) -= psrq*D_bb;
                    }
                }
            }
        }

    }
}

void CGHF::form_F() {
    F_.zero();
    F_ += F0_;
    F_ += JKwK_;


}

auto bubbleSort(std::map<int, double>& map) {
    // Convert map to vector of pairs for sorting
    std::vector<std::pair<int, double>> vec(map.begin(), map.end());

    int n = vec.size();
    for (int i = 0; i < n - 1; i++) {
        bool swapped = false;
        for (int j = 0; j < n - i - 1; j++) {
            if (vec[j].second > vec[j + 1].second) {
                std::swap(vec[j], vec[j + 1]);
                swapped = true;
            }
        }
        if (!swapped)
            break;
    }

   return vec;
}

void CGHF::sort_real_evals(){
    timer_on("GHF sort_real_evals");
    // Sorted Evecs
    auto Eval_blocks = RealEvals_.vector_data();
    auto C_unsorted = C_;

    for (int i = 0; i < nirrep_; i++) {
            std::map<int, double> eval_map_temp;
            int orderpi [irrep_sizes_[i]];
            for (int j = 0; j < 2*nsopi_[i]; j++) {
                double eval = Eval_blocks[i](j);
                eval_map_temp.insert(std::pair(j, eval));
            }
            auto eval_map = bubbleSort(eval_map_temp);

            int counter = 0;
            for (const auto& pair : eval_map) {
                orderpi[counter] = pair.first;
                counter += 1;
            }
            for (int j = 0; j < irrep_sizes_[i]; j++) {
                for (int k = 0; k < irrep_sizes_[i]; k++) {
                    C_[i](j, k) = C_unsorted[i](j, orderpi[k]);
                    RealEvals_[i](k) = Eval_blocks[i](orderpi[k]);
                }
                }
            }

    timer_off("GHF sort_real_evals");

}


void CGHF::form_C(double shift) {
    temp.zero();
    Fp_.zero();
    //Orthogonalize Fock matrix:
    //F' = X'FX

    if (options_.get_bool("DIIS") && iteration_ > 1) {
        auto error_trace = do_diis();
    } else {
        einsums::linear_algebra::gemm<false, false>(std::complex<double>{1.0}, F_, EINX_, std::complex<double>{0.0},  &temp);
        einsums::linear_algebra::gemm<true, false>(std::complex<double>{1.0}, EINX_, temp, std::complex<double>{0.0}, &Fp_);
    }

    //Diagonalize Fock matrix
    for (int i = 0; i < nirrep_; i++) {
         if (irrep_sizes_[i] > 0) {
             einsums::linear_algebra::geev(&Fp_[i], &Fevals_[i], &temp[i], &Fevecs_[i]);
             //einsums::linear_algebra::heev(&Fp_[i], &RealEvals_[i]);
         }
  
      }
  
    for (int i = 0; i < nirrep_; i++) {
        auto real_eval_irrep = GetRealVector(Fevals_[i], irrep_sizes_[i]);
        RealEvals_[i] = real_eval_irrep;
    }

    //sort_real_evals();

    einsums::linear_algebra::gemm<false, false>(std::complex<double>{1.0}, EINX_, Fevecs_, std::complex<double>{0.0},  &C_);
    sort_real_evals();
}

std::tuple<SharedMatrix, SharedMatrix> CGHF::einsums_to_numpy(std::string mat_str) {
    std::shared_ptr<Matrix> real_mat;
    std::shared_ptr<Matrix> imag_mat;
    einsums::BlockTensor<std::complex<double>, 2> ein_mat;

    if (mat_str == "D") {
        ein_mat = D_;
    } else if (mat_str == "F") {
        ein_mat = F_;
    } else if (mat_str == "C") {
        ein_mat = C_;
    } else if (mat_str == "JK") {
        ein_mat = JKwK_;
    }
    
    real_mat = std::make_shared<Matrix>(irrep_sizes_, irrep_sizes_);
    imag_mat = std::make_shared<Matrix>(irrep_sizes_, irrep_sizes_);

    for (int i = 0; i < nirrep_; i++) {
	auto dim1 = ein_mat[i].dim(0);
	auto dim2 = ein_mat[i].dim(0);

        for (int j = 0; j < dim1; j++) {
            for (int k = 0; k < dim2; k++) {
                //std::cout << j << " " << D_[i](j, k).real() << "\n";

                real_mat->set(j, k, ein_mat[i](j, k).real());
                imag_mat->set(j, k, ein_mat[i](j, k).imag());

            }
        }
    }

    return {real_mat, imag_mat};
}

void CGHF::form_numpy_D() {
    auto Dr_mat = std::make_shared<Matrix>(irrep_sizes_, irrep_sizes_);
    auto Di_mat = std::make_shared<Matrix>(irrep_sizes_, irrep_sizes_);

    for (int i = 0; i < nirrep_; i++) {
        for (int j = 0; j < irrep_sizes_[i]; j++) {
            for (int k = 0; k < irrep_sizes_[i]; k++) {
                //std::cout << j << " " << D_[i](j, k).real() << "\n";

                Dr_mat->set(j, k, D_[i](j, k).real());
                Di_mat->set(j, k, D_[i](j, k).imag());

            }
        }
    }
    set_array_variable("real_cghf_D", Dr_mat);
    set_array_variable("imag_cghf_D", Di_mat);

}

void CGHF::form_D() {
    D_.zero();
    Cocc_.zero();
    cCocc_.zero();
    auto nelec_ = nalpha_ + nbeta_;
    //std::cout << nelec_ << "\n";

    for (int i = 0; i < nirrep_; i++) {
        //Cocc.zero();
        //cCocc.zero();
        for (int j = 0; j < irrep_sizes_[i]; j++) {
                for (int k = 0; k < nelec_; k++) {
                        Cocc_[i](j,k) = C_[i](j, k);
                        cCocc_[i](j,k) = std::conj(C_[i](j, k));

                }
        }
    }

    einsums::tensor_algebra::einsum(
        einsums::Indices{einsums::index::i,
                                         einsums::index::j},
        &D_,
        einsums::Indices{einsums::index::i,
                                         einsums::index::m},
        Cocc_,
        einsums::Indices{einsums::index::j,
                                         einsums::index::m},
        cCocc_);

    //println(C_);
    //println(Cocc_);
    //println(cCocc_);
    //println(D_);
}

void CGHF::damping_update(double damping_percentage) {
}

double CGHF::compute_E() {
    auto TensorE = einsums::Tensor<std::complex<double>, 0>("E");
    auto Energy_1e = einsums::Tensor<std::complex<double>, 0>("JK E");
    auto Energy_T = einsums::Tensor<std::complex<double>, 0>("JK E");
    auto Energy_2e = einsums::Tensor<std::complex<double>, 0>("JK E");

    //auto temp2 = einsums::BlockTensor<std::complex<double>, 2>("temp", irrep_sizes_);
    //auto JKwK_ = (*JKwK_);

    temp2 = F0_;
    for (int i = 0; i < nirrep_; i++) {
        auto JK_block = JKwK_[i];
        //auto K_block = K_[i];

        JK_block *= 0.5;

        temp2[i] += JK_block;
    }

    einsums::tensor_algebra::einsum(
        0.0, einsums::Indices{}, &TensorE, 1.0,
        einsums::Indices{einsums::index::i,
                                         einsums::index::j},
        D_,
        einsums::Indices{einsums::index::i,
                                         einsums::index::j},
        temp2);

    std::complex<double> E_complex = (std::complex<double>)TensorE;

    einsums::tensor_algebra::einsum(
        0.0, einsums::Indices{}, &Energy_T, 1.0,
        einsums::Indices{einsums::index::i,
                                         einsums::index::j},
        D_,
        einsums::Indices{einsums::index::i,
                                         einsums::index::j},
        EINT_);

    std::complex<double> KE_complex = (std::complex<double>)Energy_T;

    einsums::tensor_algebra::einsum(
        0.0, einsums::Indices{}, &Energy_1e, 1.0,
        einsums::Indices{einsums::index::i,
                                         einsums::index::j},
        D_,
        einsums::Indices{einsums::index::i,
                                         einsums::index::j},
        F0_);

    std::complex<double> OneE_complex = (std::complex<double>)Energy_1e;

    einsums::tensor_algebra::einsum(
        0.0, einsums::Indices{}, &Energy_2e, 1.0,
        einsums::Indices{einsums::index::i,
                                         einsums::index::j},
        D_,
        einsums::Indices{einsums::index::i,
                                         einsums::index::j},
        JKwK_);

    std::complex<double> TwoE_complex = (std::complex<double>)Energy_2e;


    energies_["Nuclear"] = nuclearrep_;
    energies_["Kinetic"] = KE_complex.real();
    energies_["One-Electron"] = OneE_complex.real();
    energies_["Two-Electron"] = 0.5*TwoE_complex.real();

    double tot_energy = 0.0;

    tot_energy += nuclearrep_;
    //tot_energy += KE_complex.real();
    tot_energy += OneE_complex.real();
    tot_energy += 0.5*TwoE_complex.real();

    return tot_energy;
}

double CGHF::compute_initial_E() {
    double E0 = compute_E();
    return 0.0;
}

std::complex<double> CGHF::do_diis() {
    int diis_max = 8;
    int diis_count = 0;

    //FDS-SDF
    form_FDSmSDF();

    auto ortho_error = einsums::BlockTensor<std::complex<double>, 2>("Orthogonalized FDSmSDF", irrep_sizes_); //eorth_it
    auto temp1 = einsums::BlockTensor<std::complex<double>, 2>("Temp FDSmSDF", irrep_sizes_);
    auto ecurr = einsums::BlockTensor<std::complex<double>, 2>("Current error", irrep_sizes_);

    temp1.zero();
    ortho_error.zero();
    ecurr.zero();
    temp.zero();
    
    einsums::linear_algebra::gemm<false, false>(std::complex<double>{1.0}, F_, EINX_, std::complex<double>{0.0},  &temp);
    einsums::linear_algebra::gemm<true, false>(std::complex<double>{1.0}, EINX_, temp, std::complex<double>{0.0}, &Fp_);

    einsums::linear_algebra::gemm<false, false>(std::complex<double>{1.0}, FDSmSDF_, EINX_, std::complex<double>{0.0},  &temp1);
    einsums::linear_algebra::gemm<true, false>(std::complex<double>{1.0}, EINX_, FDSmSDF_, std::complex<double>{0.0},  &ortho_error);

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

    for (int i = 0; i < err_vecs.size(); i++) {
        for (int j = 0; j <= i; j++) { //Had to consult Connor's UHF for this little trick as I tried 80 million things
            B_(i, j) = einsums::linear_algebra::dot(err_vecs[i], err_vecs[j]);
            B_(j, i) = B_(i, j); //Which makes perfect sense actually
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

    auto temp2 = einsums::BlockTensor<std::complex<double>, 2>("temp2", irrep_sizes_);
    temp2.zero();
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
	temp.zero();
        einsums::linear_algebra::gemm<false, false>(std::complex<double>{1.0}, F_, EINX_, std::complex<double>{0.0},  &temp);
        einsums::linear_algebra::gemm<true, false>(std::complex<double>{1.0}, EINX_, temp, std::complex<double>{0.0}, &Fp_);
    }

    //Fp_ = temp2;

    //F0_ = temp2;
    //

    return error_trace;

}

void CGHF::form_FDSmSDF() {
    //FDS-SDF for DIIS
    auto SComplex = einsums::BlockTensor<std::complex<double>, 2>("Complex overlap matrix", irrep_sizes_);
    SComplex.zero();

    for (int i = 0; i < nirrep_; i++) {
        for (int j = 0; j < irrep_sizes_[i]; j++) {
            for (int k = 0; k < irrep_sizes_[i]; k++) {
                SComplex[i].subscript(j, k) = EINS_[i].subscript(j, k);
            }
        }
    }

    //auto FDS = einsums::BlockTensor<std::complex<double>, 2>("FDS", irrep_sizes_);
    //auto SDF = einsums::BlockTensor<std::complex<double>, 2>("SDF", irrep_sizes_);

    //auto FDS1 = einsums::BlockTensor<std::complex<double>, 2>("FDS Buffer", irrep_sizes_);
    //auto SDF1 = einsums::BlockTensor<std::complex<double>, 2>("SDF Buffer", irrep_sizes_);

    //FDS.zero();
    //SDF.zero();
    //FDS1.zero();
    //SDF1.zero();
    temp.zero();
    temp2.zero();

    //FDS
    einsums::linear_algebra::gemm<false, false>(std::complex<double>{1.0}, D_, SComplex, std::complex<double>{0.0},  &temp);
    einsums::linear_algebra::gemm<false, false>(std::complex<double>{1.0}, F_, temp, std::complex<double>{0.0}, &temp2);

    temp.zero();
    FDSmSDF_ = temp2;
    temp2.zero();

    //SDF
    einsums::linear_algebra::gemm<false, false>(std::complex<double>{1.0}, D_, F_, std::complex<double>{0.0}, &temp);
    einsums::linear_algebra::gemm<false, false>(std::complex<double>{1.0}, SComplex, temp, std::complex<double>{0.0}, &temp2);

    //FDSmSDF = FDS;
    FDSmSDF_ -= temp2;
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




