/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "direct_screening.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libmints/integral.h"
using namespace psi;

DirectScreening::DirectScreening(std::shared_ptr<BasisSet> basis_in,
           std::vector<SharedMatrix>& density_in)
:
basis_(basis_in),
do_J_(true),
do_K_(true),
D_(density_in)
{

  factory_ = std::shared_ptr<IntegralFactory>(new IntegralFactory(basis_,
                                                                    basis_,
                                                                    basis_,
                                                                    basis_));

  Options& options = Process::environment.options;
  double cutoff = options.get_double("SCHWARZ_CUTOFF");
  sieve_ = std::shared_ptr<ERISieve>(new ERISieve(basis_, cutoff));
  //std::cout << "cutoff: " << cutoff << "\n";

  eri_.clear();
  //for (int thread = 0; thread < omp_nthread_; thread++) {
  eri_.push_back(std::shared_ptr<TwoBodyAOInt>(factory_->eri()));
  //}

  for (size_t N = 0; N < D_.size() && do_J_; ++N) {
    std::stringstream s;
    s << "J " << N << " (AO)";
    J_.push_back(SharedMatrix(new Matrix(s.str(),D_[N]->nirrep(),
                                         D_[N]->rowspi(), D_[N]->rowspi(),
                                         D_[N]->symmetry())));
    J_[N]->zero();
  }

  // start K_
  for (size_t N = 0; N < D_.size() && do_K_; ++N) {
    std::stringstream s;
    s << "K " << N << " (AO)";
    K_.push_back(SharedMatrix(new Matrix(s.str(),D_[N]->nirrep(),
                                         D_[N]->rowspi(), D_[N]->rowspi(),
                                         D_[N]->symmetry())));
    K_[N]->zero();
  }

} // Constructor

DirectScreening::~DirectScreening() {}

void DirectScreening::Compute()
{

  std::cout << "Doing DirectScreening.\n";

  /*
  double most_negative_integral = 0.0;
  double most_positive_integral = 0.0;

  int most_negative_mu = -1;
  int most_negative_nu = -1;
  int most_negative_rho = -1;
  int most_negative_sig = -1;

  int most_positive_mu = -1;
  int most_positive_nu = -1;
  int most_positive_rho = -1;
  int most_positive_sig = -1;
   */
  const double* buffer = eri_[0]->buffer();

  // Index the significant MN shell pairs (triangular M,N)
  const std::vector<std::pair<int,int> >& MN = sieve_->shell_pairs();

  for (size_t mu_nu_ind = 0L; mu_nu_ind < MN.size(); ++mu_nu_ind) {

    int mu_ind = MN[mu_nu_ind].first;
    int nu_ind = MN[mu_nu_ind].second;

    int num_funs_mu = basis_->shell(mu_ind).nfunction();
    int num_funs_nu = basis_->shell(nu_ind).nfunction();

    int mu_fun_start = basis_->shell(mu_ind).function_index();
    int nu_fun_start = basis_->shell(nu_ind).function_index();


    for (size_t rho_sig_ind = 0L; rho_sig_ind < MN.size(); ++rho_sig_ind) {

      int rho_ind = MN[rho_sig_ind].first;
      int sig_ind = MN[rho_sig_ind].second;

      if (sieve_->shell_significant(mu_ind, nu_ind, rho_ind, sig_ind)) {

        //std::cout << "computing shell (" << mu_ind << " " << nu_ind << "|";
        //std::cout << rho_ind << " " << sig_ind << ")\n";
        eri_[0]->compute_shell(mu_ind, nu_ind, rho_ind, sig_ind);

        int num_funs_rho = basis_->shell(rho_ind).nfunction();
        int num_funs_sig = basis_->shell(sig_ind).nfunction();

        int rho_fun_start = basis_->shell(rho_ind).function_index();
        int sig_fun_start = basis_->shell(sig_ind).function_index();

        for (int mu_fun_ind = 0, index = 0; mu_fun_ind < num_funs_mu; mu_fun_ind++) {
          for (int nu_fun_ind = 0; nu_fun_ind < num_funs_nu; nu_fun_ind++) {
            for (int rho_fun_ind = 0; rho_fun_ind < num_funs_rho; rho_fun_ind++) {
              for (int sig_fun_ind = 0; sig_fun_ind < num_funs_sig; sig_fun_ind++, index++) {

                double val = buffer[index];
                /*
                if (val < most_negative_integral) {

                  most_negative_integral = val;

                  most_negative_mu = mu_ind;
                  most_negative_nu = nu_ind;
                  most_negative_rho = rho_ind;
                  most_negative_sig = sig_ind;

                }

                if (val > most_positive_integral) {

                  most_positive_integral = val;

                  most_positive_mu = mu_ind;
                  most_positive_nu = nu_ind;
                  most_positive_rho = rho_ind;
                  most_positive_sig = sig_ind;

                }
                */
                int m = mu_fun_ind + mu_fun_start;
                int n = nu_fun_ind + nu_fun_start;
                int r = rho_fun_ind + rho_fun_start;
                int s = sig_fun_ind + sig_fun_start;

                if (do_J_)
                {

                  double ref_sym_val = rho_ind == sig_ind ? 1.0 : 2.0;

                  for (size_t N = 0; N < J_.size(); N++) {

                    double coulomb_val = D_[N]->get(0,r,s) * val * ref_sym_val;

                    //std::cout << "Adding value to J: " << coulomb_val << "\n";
                    J_[N]->add(0,m,n, coulomb_val);
                    if (mu_ind != nu_ind) {
                      J_[N]->add(0,n,m, coulomb_val);
                    }

                  }

                } // contract to coulomb matrix

                if (do_K_) {
                  for (size_t N = 0; N < K_.size(); N++) {

                    K_[N]->add(0,m,r, D_[N]->get(0,n,s)*val);

                    if (mu_ind != nu_ind) {
                      K_[N]->add(0,n,r, D_[N]->get(0,m,s)*val);
                    }

                    if (rho_ind != sig_ind) {
                      K_[N]->add(0,m,s, D_[N]->get(0,n,r) * val);
                    }

                    if (rho_ind != sig_ind && mu_ind != nu_ind) {
                      K_[N]->add(0,n,s, D_[N]->get(0,m,r) * val);
                    }

                  } // loop over N
                }
              }}}} // functions in the shell

      } // if the shell is significant

    } // loop over ket shell pairs

  } // loop over bra shell pairs
  /*
  std::cout << "Most negative integral: " << most_negative_integral << "\n";
  std::cout << "Indices: " << most_negative_mu << ", ";
  std::cout << most_negative_nu << ", " << most_negative_rho;
  std::cout << ", " << most_negative_sig << "\n";
  std::cout << "Negative shells: \n";
  basis_->shell(most_negative_mu).print(stdout);
  basis_->shell(most_negative_nu).print(stdout);
  basis_->shell(most_negative_rho).print(stdout);
  basis_->shell(most_negative_sig).print(stdout);
  std::cout << "\n\n";

  std::cout << "Most positive integral: " << most_positive_integral << "\n";
  std::cout << "Indices: " << most_positive_mu << ", ";
  std::cout << most_positive_nu << ", " << most_positive_rho;
  std::cout << ", " << most_positive_sig << "\n\n";
  std::cout << "Positive shells: \n";
  basis_->shell(most_positive_mu).print(stdout);
  basis_->shell(most_positive_nu).print(stdout);
  basis_->shell(most_positive_rho).print(stdout);
  basis_->shell(most_positive_sig).print(stdout);
  std::cout << "\n\n";
  */
} // Compute

// Assuming this always gets called before Compute()
void DirectScreening::Update(const std::vector<SharedMatrix>& D_new)
{

  D_ = D_new;

  J_.clear();
  for (size_t N = 0; N < D_.size() && do_J_; ++N) {
    std::stringstream s;
    s << "J " << N << " (AO)";
    J_.push_back(SharedMatrix(new Matrix(s.str(),D_[N]->nirrep(),
                                         D_[N]->rowspi(), D_[N]->rowspi(),
                                         D_[N]->symmetry())));
    J_[N]->zero();

  }

  K_.clear();
  for (size_t N = 0; N < D_.size() && do_K_; ++N) {
    std::stringstream s;
    s << "K " << N << " (AO)";
    K_.push_back(SharedMatrix(new Matrix(s.str(),D_[N]->nirrep(),
                                         D_[N]->rowspi(), D_[N]->rowspi(),
                                         D_[N]->symmetry())));
    K_[N]->zero();

  }

}

std::vector<SharedMatrix>& DirectScreening::J()
{
  return J_;
}

std::vector<SharedMatrix>& DirectScreening::K() {
  return K_;
}

void DirectScreening::set_do_J(bool do_it)
{
  do_J_ = do_it;
}

void DirectScreening::set_do_K(bool do_it)
{
  do_K_ = do_it;
}

void DirectScreening::print_header() const
{

  outfile->Printf( "  ==> Direct Screening Exchange Matrix Calculation <==\n\n");

  outfile->Printf( "    J tasked:          %11s\n", (do_J_ ? "Yes" : "No"));
  outfile->Printf( "    K tasked:          %11s\n", (do_K_ ? "Yes" : "No"));
  outfile->Printf( "    wK tasked:         %11s\n", (do_wK_ ? "Yes" : "No"));

  //outfile->Printf( "    OpenMP threads:    %11d\n", omp_nthread_);
  //outfile->Printf( "    Memory (MB):       %11ld\n", (memory_ *8L) / (1024L * 1024L));
  //outfile->Printf( "    Schwarz Cutoff:    %11.0E\n\n", cutoff_);

}


