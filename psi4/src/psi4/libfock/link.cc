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

#include "link.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
using namespace psi;

LinK::LinK(std::shared_ptr<BasisSet> basis_in,
           std::vector<SharedMatrix>& density_in)
:
basis_(basis_in),
do_J_(false),
do_K_(true),
D_(density_in),
nu_for_mu_indices_(basis_->nshell()),
shell_to_shell_(basis_->nshell()),
shell_max_q_sqr_(basis_->nshell(), -DBL_MAX),
shell_pair_threshold_sqr_(0.0),
integral_threshold_sqr_(0.0),
num_integrals_(0),
total_num_integrals_(0)
{

  factory_ = std::shared_ptr<IntegralFactory>(new IntegralFactory(basis_,
                                                                    basis_,
                                                                    basis_,
                                                                    basis_));

  Options& options = Process::environment.options;
  double cutoff = options.get_double("SCHWARZ_CUTOFF");
  sieve_ = std::shared_ptr<ERISieve>(new ERISieve(basis_, cutoff));
  //std::cout << "cutoff: " << cutoff << "\n";

  shell_pair_threshold_sqr_ = cutoff * cutoff;
  integral_threshold_sqr_ = cutoff * cutoff;

  eri_.clear();
  //for (int thread = 0; thread < omp_nthread_; thread++) {
  eri_.push_back(std::shared_ptr<TwoBodyAOInt>(factory_->eri()));
  //}

  // start K_
  for (size_t N = 0; N < D_.size() && do_K_; ++N) {
    std::stringstream s;
    s << "K " << N << " (AO)";
    K_.push_back(SharedMatrix(new Matrix(s.str(),D_[N]->nirrep(),
                                         D_[N]->rowspi(), D_[N]->rowspi(),
                                         D_[N]->symmetry())));
    K_[N]->zero();
  }



  // fill in the max Q factor for each shell
  const std::vector<std::vector<int> >& shell_to_shell = sieve_->shell_to_shell();

  for (int mu_ind = 0; mu_ind < basis_->nshell(); mu_ind++) {

    for (std::vector<int>::const_iterator nu_ind = shell_to_shell[mu_ind].begin();
         nu_ind != shell_to_shell[mu_ind].end(); nu_ind++) {

      double this_val = sieve_->shell_pair_value(mu_ind, *nu_ind);

      shell_max_q_sqr_[mu_ind] = std::max(this_val, shell_max_q_sqr_[mu_ind]);

    } // loop over nu

  } // loop over mu

  // fill in shell_to_shell_
  const std::vector<std::vector<int> >& shell_to_shell_index = sieve_->shell_to_shell();

  for (int mu_ind = 0; mu_ind < basis_->nshell(); mu_ind++) {

    const std::vector<int>& mu_shells = shell_to_shell_index[mu_ind];

    for (std::vector<int>::const_iterator nu_ind = mu_shells.begin();
         nu_ind != mu_shells.end(); nu_ind++) {

      // This needs to be triangular for symmetry later
      if (*nu_ind >= mu_ind) {
        shell_to_shell_[mu_ind].push_back(std::make_pair(sieve_->shell_pair_value(mu_ind, *nu_ind),
                                                         *nu_ind));
      }

    } // loop over nu

    std::sort(shell_to_shell_[mu_ind].begin(), shell_to_shell_[mu_ind].end(),
              std::greater<std::pair<double, int> >());

  } // loop over mu

} // Constructor

LinK::~LinK() {
  std::cout << "LinK total integrals: " << total_num_integrals_ << "\n";
}

void LinK::Compute()
{

  std::cout << "Doing LINK for K matrix.\n";
  //std::cout << "Doing LinK computation.\n";
  // for each shell mu, need a list of significant nu
  // sort these lists (for each mu) when finished

  // for each bra shell pair, form list of ket shell pairs
  FormSignificantShellPairList_();

  std::cout << "Num integrals: " << num_integrals_ << "\n\n";

}

// Assuming this always gets called before Compute()
void LinK::Update(const std::vector<SharedMatrix>& D_new)
{

  D_ = D_new;

  K_.clear();
  for (size_t N = 0; N < D_.size() && do_K_; ++N) {
    std::stringstream s;
    s << "K " << N << " (AO)";
    K_.push_back(SharedMatrix(new Matrix(s.str(),D_[N]->nirrep(),
                                         D_[N]->rowspi(), D_[N]->rowspi(),
                                         D_[N]->symmetry())));
  }

  // clear them first so we can do symmetry here
  for (int i = 0; i < basis_->nshell(); i++) {
    nu_for_mu_indices_[i].clear();
  }

  for (int mu_ind = 0; mu_ind < basis_->nshell(); mu_ind++) {
    FindSignificantNuforMu_(mu_ind);
  }

  num_integrals_ = 0;

}

void LinK::FindSignificantNuforMu_(int mu_ind)
{

  // this tries to do symmetry right now -- is it correct?
  for (int nu_ind = mu_ind; nu_ind < basis_->nshell(); nu_ind++) {

    double density_mu_nu_sqr = 0.0;

    // we need to pick out the largest density matrix entry
    for (int m = basis_->shell(mu_ind).function_index();
         m < basis_->shell(mu_ind).function_index() + basis_->shell(mu_ind).nfunction();
         m++) {
      for (int n = basis_->shell(nu_ind).function_index();
           n < basis_->shell(nu_ind).function_index() + basis_->shell(nu_ind).nfunction(); n++) {

        // NOTE: this assumes that I'm both single threaded and C1 symmetry
        double this_density = D_[0]->get(0,m,n) * D_[0]->get(0,m,n);
        density_mu_nu_sqr = std::max(density_mu_nu_sqr, this_density);

      } // loop over nu functions
    } // loop over mu functions

    //std::cout << mu_ind << ", " << nu_ind << " density: " << density_mu_nu_sqr << "\n";

    if (density_mu_nu_sqr * shell_max_q_sqr_[mu_ind] * shell_max_q_sqr_[nu_ind] > shell_pair_threshold_sqr_) {

      nu_for_mu_indices_[mu_ind].push_back(std::make_pair(density_mu_nu_sqr * shell_max_q_sqr_[nu_ind], nu_ind));
      // can't do symmetry because they'll be sorted by density as well

      // do I really want to handle symmetry here
      if (nu_ind != mu_ind) {
        nu_for_mu_indices_[nu_ind].push_back(std::make_pair(density_mu_nu_sqr * shell_max_q_sqr_[mu_ind], mu_ind));
      }
    }
    /*
    else {
      std::cout << "Pruned a nu for mu: mu_ind: " << mu_ind << ", ";
      std::cout << "nu_ind: " << nu_ind << "\n";
    }
     */

  } // loop over nu

  // now, sort the nu's for each mu
  // sort by d_munu * max_nu
  // assuming this sorts by the first element of the pair
  std::sort(nu_for_mu_indices_[mu_ind].begin(),
            nu_for_mu_indices_[mu_ind].end(),
            std::greater<std::pair<double, int> >());

}

// for each bra shell pair, form list of ket shell pairs
void LinK::FormSignificantShellPairList_()
{

  // for all shell pairs
  // these need to be sorted by their value of Q
  const std::vector<std::pair<int,int> >& shell_pairs = sieve_->shell_pairs();

  //const std::vector<long int>& shell_pairs_reverse = sieve_->shell_pairs_reverse();

  for (size_t i = 0; i < shell_pairs.size(); i++) {

    int mu_ind = shell_pairs[i].first;
    int lambda_ind = shell_pairs[i].second;

    //std::cout << "Shell pair: " << mu_ind << ", " << lambda_ind << "\n";

    double mu_lambda_val_sqr = sieve_->shell_pair_value(mu_ind, lambda_ind);

    std::vector<std::pair<int,int> > ml_integrals;

    // Loop for mu
    std::vector<std::pair<int,int> > ml_integrals_mu;

    // these are already sorted
    for (std::vector<std::pair<double, int> >::iterator nu_ind = nu_for_mu_indices_[mu_ind].begin();
         nu_ind != nu_for_mu_indices_[mu_ind].end(); nu_ind++) {

      int num_added_sigmas = 0;

      double density_mu_nu_sqr = nu_ind->first;

      // loop over sigmas that count for nu
      // this needs to be triangular - it is now because shell_to_shell_ is
      for (std::vector<std::pair<double, int> >::iterator sigma_ind = shell_to_shell_[nu_ind->second].begin();
           sigma_ind != shell_to_shell_[nu_ind->second].end(); sigma_ind++) {

        // TODO: I think I can get this out of the stored values with nu
        double nu_sigma_val_sqr = sigma_ind->first;

        if (density_mu_nu_sqr * mu_lambda_val_sqr * nu_sigma_val_sqr > integral_threshold_sqr_) {

          // we found an integral we need to compute
          ml_integrals_mu.push_back(std::make_pair(nu_ind->second, sigma_ind->second));
          num_added_sigmas++;

        }
        else {
          // now that sigmas are sorted, this is ok
          break;
        }

      } // loop over sigma

      // if number of sigmas we found is zero, leave, since the array of
      // nu's is sorted
      if (num_added_sigmas == 0) {
        break;
      }

    } // loop over nu (for mu)

    std::sort(ml_integrals_mu.begin(), ml_integrals_mu.end());

    // Loop for lambda
    std::vector<std::pair<int,int> > ml_integrals_lambda;

    if (mu_ind != lambda_ind) {

      for (std::vector<std::pair<double, int> >::iterator nu_ind = nu_for_mu_indices_[lambda_ind].begin();
           nu_ind != nu_for_mu_indices_[lambda_ind].end(); nu_ind++) {

        double density_lambda_nu_sqr = nu_ind->first;

        int num_added_sigmas = 0;

        // loop over sigmas that count for nu
        for (std::vector<std::pair<double, int> >::iterator sigma_ind = shell_to_shell_[nu_ind->second].begin();
             sigma_ind != shell_to_shell_[nu_ind->second].end(); sigma_ind++) {

          double nu_sigma_val_sqr = sigma_ind->first;

          if (density_lambda_nu_sqr * mu_lambda_val_sqr * nu_sigma_val_sqr > integral_threshold_sqr_) {

            // we found an integral we need to compute
            ml_integrals_lambda.push_back(std::make_pair(nu_ind->second, sigma_ind->second));
            num_added_sigmas++;

          }
          else {
            break;
          }

        } // loop over sigma

        // if number of sigmas we found is zero, leave, since the array of
        // nu's is sorted
        if (num_added_sigmas == 0) {
          break;
        }

      } // loop over nu (for lambda)

      // the ordering doesn't matter here, just need them to be sorted for union
      std::sort(ml_integrals_lambda.begin(), ml_integrals_lambda.end());

    } // don't do anything if mu and lambda are equal

    // this is the largest one possible
    ml_integrals.resize(ml_integrals_mu.size() + ml_integrals_lambda.size());

    // now, merge the two lists
    std::vector<std::pair<int, int> >::iterator ml_int_end =
        std::set_union(ml_integrals_mu.begin(), ml_integrals_mu.end(),
                       ml_integrals_lambda.begin(), ml_integrals_lambda.end(),
                       ml_integrals.begin());

    int num_ml_ints = ml_int_end - ml_integrals.begin();
    //std::cout << "num ml_ints: " << num_ml_ints << "\n";
    ml_integrals.resize(num_ml_ints);

    // now, compute the integrals and contract into the result matrix
    ContractIntegrals_(mu_ind, lambda_ind, ml_integrals);

  } // loop over all bra shell pairs

} // FormSignificantShellPairList_()


void LinK::ContractIntegrals_(int mu_ind, int lambda_ind,
                              std::vector<std::pair<int,int> >& ml_integrals)
{

  //std::cout << "Contracting for " << mu_ind << ", " << lambda_ind;
  //std::cout << " and " << ml_integrals.size() << " others.\n";

  const double* buffer = eri_[0]->buffer();

  int num_funs_mu = basis_->shell(mu_ind).nfunction();
  int num_funs_lambda = basis_->shell(lambda_ind).nfunction();

  int mu_fun_start = basis_->shell(mu_ind).function_index();
  int lambda_fun_start = basis_->shell(lambda_ind).function_index();

  for (std::vector<std::pair<int,int> >::iterator ns_ind = ml_integrals.begin();
       ns_ind != ml_integrals.end(); ns_ind++) {

    int nu_ind = ns_ind->first;
    int sigma_ind = ns_ind->second;

    int num_funs_nu = basis_->shell(nu_ind).nfunction();
    int num_funs_sigma = basis_->shell(sigma_ind).nfunction();

    int nu_fun_start = basis_->shell(nu_ind).function_index();
    int sigma_fun_start = basis_->shell(sigma_ind).function_index();

    //std::cout << "LinK computing integral (" << mu_ind << " ";
    //std::cout << lambda_ind << "|" << nu_ind << " " << sigma_ind << ")\n";
    eri_[0]->compute_shell(mu_ind, lambda_ind, nu_ind, sigma_ind);
    num_integrals_++;
    total_num_integrals_++;

    for (int mu_fun_ind = 0, index = 0; mu_fun_ind < num_funs_mu; mu_fun_ind++) {
      for (int lambda_fun_ind = 0; lambda_fun_ind < num_funs_lambda; lambda_fun_ind++) {
        for (int nu_fun_ind = 0; nu_fun_ind < num_funs_nu; nu_fun_ind++) {
          for (int sigma_fun_ind = 0; sigma_fun_ind < num_funs_sigma; sigma_fun_ind++, index++) {

                    double val = buffer[index];

                    int m = mu_fun_ind + mu_fun_start;
                    int l = lambda_fun_ind + lambda_fun_start;
                    int n = nu_fun_ind + nu_fun_start;
                    int s = sigma_fun_ind + sigma_fun_start;

                    for (size_t N = 0; N < K_.size(); N++) {

                      K_[N]->add(0,l,s, D_[N]->get(0,m,n)*val);

                      if (nu_ind != sigma_ind) {
                        K_[N]->add(0,l,n, D_[N]->get(0,m,s)*val);
                      }

                      if (mu_ind != lambda_ind) {
                        K_[N]->add(0,m,s, D_[N]->get(0,l,n) * val);
                      }

                      if (nu_ind != sigma_ind && mu_ind != lambda_ind) {
                        K_[N]->add(0,m,n, D_[N]->get(0,l,s) * val);
                      }

                    } // loop over N

          }}}} // functions in the shell
  } // loop over ket shell pairs

} //contract integrals


std::vector<SharedMatrix>& LinK::J()
{
  throw PSIEXCEPTION("LinK can't do J matrix.");
}

std::vector<SharedMatrix>& LinK::K() {
  return K_;
}

void LinK::set_do_J(bool do_it)
{
  if (do_it)
    throw PSIEXCEPTION("LinK can't do J matrix.");
  else
    do_J_ = false;
}

void LinK::set_do_K(bool do_it)
{
  do_K_ = do_it;
}

void LinK::print_header() const
{

  outfile->Printf( "  ==> LinK Exchange Matrix Calculation <==\n\n");

  outfile->Printf( "    J tasked:          %11s\n", (do_J_ ? "Yes" : "No"));
  outfile->Printf( "    K tasked:          %11s\n", (do_K_ ? "Yes" : "No"));
  outfile->Printf( "    wK tasked:         %11s\n", (do_wK_ ? "Yes" : "No"));

  //outfile->Printf( "    OpenMP threads:    %11d\n", omp_nthread_);
  //outfile->Printf( "    Memory (MB):       %11ld\n", (memory_ *8L) / (1024L * 1024L));
  //outfile->Printf( "    Schwarz Cutoff:    %11.0E\n\n", cutoff_);

}


