/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

#include <libmints/mints.h>
#include <libqt/qt.h>
#include <psi4-dec.h>
#include "sieve.h"

using namespace std;
using namespace psi;

namespace psi {

ERISieve::ERISieve(boost::shared_ptr<BasisSet> primary, double sieve) :
    primary_(primary), sieve_(sieve)
{
    common_init();
}
ERISieve::~ERISieve()
{
    delete[] function_pair_values_;
    delete[] shell_pair_values_;
}
void ERISieve::common_init()
{
  // if sieve_ is 0, then erfc_inv is infinite
  // the boost function just throws an error in this case
  if (sieve_ > 0.0)
    erfc_thresh_ = boost::math::erfc_inv(sieve_);
  else
    erfc_thresh_ = DBL_MAX;
  
  Options& options = Process::environment.options;
  do_qqr_ = options.get_bool("DO_QQR_SIEVE");

    debug_ = 0;

    integrals();
    set_sieve(sieve_);
}
void ERISieve::set_sieve(double sieve)
{
    sieve_ = sieve;
    sieve2_ = sieve_ * sieve;
    sieve_over_max_ = sieve_ / max_;
    sieve2_over_max_ = sieve2_ / max_;

    shell_pairs_.clear();
    function_pairs_.clear();
    shell_pairs_reverse_.resize(nshell_ * (nshell_ + 1L) / 2L);
    function_pairs_reverse_.resize(nbf_ * (nbf_ + 1L) / 2L);

    long int offset = 0L;
    unsigned long int MUNU = 0L;
    for (int MU = 0; MU < nshell_; MU++) {
        for (int NU = 0; NU <= MU; NU++, MUNU++) {
            if (shell_pair_values_[MU * (unsigned long int) nshell_ + NU] >= sieve2_over_max_) {
                shell_pairs_.push_back(make_pair(MU,NU));
                shell_pairs_reverse_[MUNU] = offset;
                offset++;
            } else {
                shell_pairs_reverse_[MUNU] = -1L;
            }
        }
    }

    offset = 0L;
    unsigned long int munu = 0L;
    for (int mu = 0; mu < nbf_; mu++) {
        for (int nu = 0; nu <= mu; nu++, munu++) {
            if (function_pair_values_[mu * (unsigned long int) nbf_ + nu] >= sieve2_over_max_) {
                function_pairs_.push_back(make_pair(mu,nu));
                function_pairs_reverse_[munu] = offset;
                offset++;
            } else {
                function_pairs_reverse_[munu] = -1L;
            }
        }
    }

    shell_to_shell_.clear();
    function_to_function_.clear();
    shell_to_shell_.resize(nshell_);
    function_to_function_.resize(nbf_);

    for (int MU = 0; MU < nshell_; MU++) {
        for (int NU = 0; NU < nshell_; NU++) {
            if (shell_pair_values_[MU * (unsigned long int) nshell_ + NU] >= sieve2_over_max_) {
                shell_to_shell_[MU].push_back(NU);
            }
        }
    }

    for (int mu = 0; mu < nbf_; mu++) {
        for (int nu = 0; nu < nbf_; nu++) {
            if (function_pair_values_[mu * (unsigned long int) nbf_ + nu] >= sieve2_over_max_) {
                function_to_function_[mu].push_back(nu);
            }
        }
    }

    if (debug_) {
        outfile->Printf( "  ==> ERISieve Debug <==\n\n");
        outfile->Printf( "    Sieve Cutoff = %11.3E\n", sieve_);
        outfile->Printf( "    Sieve^2      = %11.3E\n", sieve2_);
        outfile->Printf( "    Max          = %11.3E\n", max_);
        outfile->Printf( "    Sieve/Max    = %11.3E\n", sieve_over_max_);
        outfile->Printf( "    Sieve^2/Max  = %11.3E\n\n", sieve2_over_max_);

        primary_->print_by_level("outfile",3);

        outfile->Printf( "   => Shell Pair Values <=\n\n");
        for (int M = 0; M < nshell_; M++) {
            for (int N = 0; N < nshell_; N++) {
                outfile->Printf( "    (%3d, %3d| = %11.3E\n", M, N, shell_pair_values_[M * nshell_ + N]);
            }
        }
        outfile->Printf( "\n");

        outfile->Printf( "   => Function Pair Values <=\n\n");
        for (int M = 0; M < nbf_; M++) {
            for (int N = 0; N < nbf_; N++) {
                outfile->Printf( "    (%3d, %3d| = %11.3E\n", M, N, function_pair_values_[M * nbf_ + N]);
            }
        }
        outfile->Printf( "\n");

        outfile->Printf( "   => Significant Shell Pairs <=\n\n");
        for (int MN = 0; MN < shell_pairs_.size(); MN++) {
            outfile->Printf( "    %6d = (%3d,%3d|\n", MN, shell_pairs_[MN].first,shell_pairs_[MN].second);
        }
        outfile->Printf( "\n");

        outfile->Printf( "   => Significant Function Pairs <=\n\n");
        for (int MN = 0; MN < function_pairs_.size(); MN++) {
            outfile->Printf( "    %6d = (%3d,%3d|\n", MN, function_pairs_[MN].first,function_pairs_[MN].second);
        }
        outfile->Printf( "\n");

        outfile->Printf( "   => Significant Shell Pairs Reverse <=\n\n");
        for (int M = 0; M < nshell_; M++) {
            for (int N = 0; N <= M; N++) {
                outfile->Printf( "    %6ld = (%3d,%3d|\n", shell_pairs_reverse_[M * (M + 1) / 2 + N], M, N);
            }
        }
        outfile->Printf( "\n");

        outfile->Printf( "   => Significant Function Pairs Reverse <=\n\n");
        for (int M = 0; M < nbf_; M++) {
            for (int N = 0; N <= M; N++) {
                outfile->Printf( "    %6ld = (%3d,%3d|\n", function_pairs_reverse_[M * (M + 1) / 2 + N], M, N);
            }
        }
        outfile->Printf( "\n");

        outfile->Printf( "   => Shell to Shell <=\n\n");
        for (int M = 0; M < nshell_; M++) {
            for (int N = 0; N < shell_to_shell_[M].size(); N++) {
                outfile->Printf( "    (%3d, %3d|\n", M, shell_to_shell_[M][N]);
            }
        }
        outfile->Printf( "\n");

        outfile->Printf( "   => Function to Function <=\n\n");
        for (int M = 0; M < nbf_; M++) {
            for (int N = 0; N < function_to_function_[M].size(); N++) {
                outfile->Printf( "    (%3d, %3d|\n", M, function_to_function_[M][N]);
            }
        }
        outfile->Printf( "\n");

    }

}
void ERISieve::integrals()
{
    int nshell = primary_->nshell();
    int nbf = primary_->nbf();

    nbf_ = nbf;
    nshell_ = nshell;

    function_pair_values_ = new double[nbf * (unsigned long int)nbf];
    shell_pair_values_ = new double[nshell * (unsigned long int)nshell];
    ::memset((void*) function_pair_values_, '\0', sizeof(double) * nbf * nbf);
    ::memset((void*) shell_pair_values_, '\0', sizeof(double) * nshell * nshell);
    max_ = 0.0;

    IntegralFactory schwarzfactory(primary_,primary_,primary_,primary_);
    boost::shared_ptr<TwoBodyAOInt> eri = boost::shared_ptr<TwoBodyAOInt>(schwarzfactory.eri());
    const double *buffer = eri->buffer();

    int MU, NU, mu, nu,omu,onu, nummu, numnu, index;
    for (MU=0; MU < nshell; ++MU) {
        nummu = primary_->shell(MU).nfunction();
        for (NU=0; NU <= MU; ++NU) {
            numnu = primary_->shell(NU).nfunction();
            eri->compute_shell(MU,NU,MU,NU);
            for (mu=0; mu < nummu; ++mu) {
                omu = primary_->shell(MU).function_index() + mu;
                for (nu=0; nu < numnu; ++nu) {
                    onu = primary_->shell(NU).function_index() + nu;

                    if (omu>=onu) {
                        index = mu*(numnu*nummu*numnu+numnu)+nu*(nummu*numnu+1);
                        if (max_ < fabs(buffer[index]))
                            max_ = fabs(buffer[index]);
                        if (shell_pair_values_[MU * (unsigned long int) nshell + NU] < fabs(buffer[index])) {
                            shell_pair_values_[MU * (unsigned long int) nshell + NU] = fabs(buffer[index]);
                            shell_pair_values_[NU * (unsigned long int) nshell + MU] = fabs(buffer[index]);
                        }
                        if (function_pair_values_[omu * (unsigned long int) nbf + onu] < fabs(buffer[index])) {
                            function_pair_values_[omu * (unsigned long int) nbf + onu] = fabs(buffer[index]);
                            function_pair_values_[onu * (unsigned long int) nbf + omu] = fabs(buffer[index]);
                        }
                    }
                }
            }
        }
    }

    if (do_qqr_) {
            
            //std::cout << "Doing QQR precomputations.\n";
            
            const GaussianShell& mu_shell = primary_->shell(MU);
            const GaussianShell& nu_shell = primary_->shell(NU);
            
            double coef_denominator = 0.0;

            std::vector<Vector3> primitive_centers;
            std::vector<double> primitive_extents;
            
            // IMPORTANT: need to compute this first, outside these loops
            Vector3 contracted_center(0.0, 0.0, 0.0);
            
            for (int p_ind = 0; p_ind < mu_shell.nprimitive(); p_ind++)
            {
              
              double p_exp = mu_shell.exp(p_ind);
              double p_coef = fabs(mu_shell.coef(p_ind));
              
              Vector3 p_center = mu_shell.center();
              p_center *= p_exp;
              
              for (int r_ind = 0; r_ind < nu_shell.nprimitive(); r_ind++)
              {
                
                double r_exp = nu_shell.exp(r_ind);
                double r_coef = fabs(nu_shell.coef(r_ind));
                
                Vector3 r_center = nu_shell.center();
                r_center *= r_exp;
                
                Vector3 primitive_center = (p_center + r_center) / (p_exp + r_exp);
                
                primitive_centers.push_back(primitive_center);
                
                coef_denominator += p_coef * r_coef;
                contracted_center += primitive_center;
                
                double this_extent = sqrt(2 / (p_exp + r_exp)) * erfc_thresh_;
                
                primitive_extents.push_back(this_extent);
                
              } // loop over primitives in NU
              
            } // loop over primitives in MU
            
            contracted_center /= coef_denominator;
            unsigned long int this_ind = MU * (unsigned long int) nshell + NU;
            unsigned long int sym_ind = NU * (unsigned long int) nshell + MU;
            
            for (int prim_it = 0; prim_it < primitive_extents.size(); prim_it++)
            {
              
              double this_dist = contracted_center.distance(primitive_centers[prim_it]);

              extents_[this_ind] = std::max(extents_[this_ind],
                                            primitive_extents[prim_it]
                                              + this_dist);
              extents_[sym_ind] = extents_[this_ind];
              
            } // now loop over the pairs of primitives
              
            
     } // doing QQR
}
bool ERISieve::shell_significant_qqr(int M, int N, int R, int S)
{
  
  double Q_mn = shell_pair_values_[N * (unsigned long int) nshell_ + M];
  double Q_rs = shell_pair_values_[R * (unsigned long int) nshell_ + S];
  
  double dist = contracted_centers_[N * (unsigned long int) nshell_ + M].distance(contracted_centers_[R * (unsigned long int) nshell_ + S]);
  
  double denom = dist - extents_[N * (unsigned long int) nshell_ + M]
                      - extents_[R * (unsigned long int) nshell_ + S];
  
  // this does the near field estimate if that's the only valid one
  // values of Q are squared
  double est = Q_mn * Q_rs / (denom > 0.0 ? denom * denom : 1.0);
  if (denom > 0.0) {
    std::cout << "Q_mn: " << Q_mn << ", ";
    std::cout << "Q_rs: " << Q_rs << ", ";
    std::cout << "dist: " << dist << ", ";
    std::cout << "denom: " << denom << ", ";
    std::cout << "est: " << est << ", ";
    std::cout << "sieve2: " << sieve2_ << "\n";
  }
  return est >= sieve2_;
  
}


double ERISieve::shell_pair_value(int m, int n) const
{
  
  return shell_pair_values_[m * nshell_ + n];
  
}
}
