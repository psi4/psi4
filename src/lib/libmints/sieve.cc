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
        fprintf(outfile, "  ==> ERISieve Debug <==\n\n");
        fprintf(outfile, "    Sieve Cutoff = %11.3E\n", sieve_);
        fprintf(outfile, "    Sieve^2      = %11.3E\n", sieve2_);
        fprintf(outfile, "    Max          = %11.3E\n", max_);
        fprintf(outfile, "    Sieve/Max    = %11.3E\n", sieve_over_max_);
        fprintf(outfile, "    Sieve^2/Max  = %11.3E\n\n", sieve2_over_max_);

        primary_->print_by_level(outfile,3);

        fprintf(outfile, "   => Shell Pair Values <=\n\n");
        for (int M = 0; M < nshell_; M++) {
            for (int N = 0; N < nshell_; N++) {
                fprintf(outfile, "    (%3d, %3d| = %11.3E\n", M, N, shell_pair_values_[M * nshell_ + N]);
            }
        }
        fprintf(outfile, "\n");

        fprintf(outfile, "   => Function Pair Values <=\n\n");
        for (int M = 0; M < nbf_; M++) {
            for (int N = 0; N < nbf_; N++) {
                fprintf(outfile, "    (%3d, %3d| = %11.3E\n", M, N, function_pair_values_[M * nbf_ + N]);
            }
        }
        fprintf(outfile, "\n");

        fprintf(outfile, "   => Significant Shell Pairs <=\n\n");
        for (int MN = 0; MN < shell_pairs_.size(); MN++) {
            fprintf(outfile, "    %6d = (%3d,%3d|\n", MN, shell_pairs_[MN].first,shell_pairs_[MN].second);
        }
        fprintf(outfile, "\n");

        fprintf(outfile, "   => Significant Function Pairs <=\n\n");
        for (int MN = 0; MN < function_pairs_.size(); MN++) {
            fprintf(outfile, "    %6d = (%3d,%3d|\n", MN, function_pairs_[MN].first,function_pairs_[MN].second);
        }
        fprintf(outfile, "\n");

        fprintf(outfile, "   => Significant Shell Pairs Reverse <=\n\n");
        for (int M = 0; M < nshell_; M++) {
            for (int N = 0; N <= M; N++) {
                fprintf(outfile, "    %6ld = (%3d,%3d|\n", shell_pairs_reverse_[M * (M + 1) / 2 + N], M, N);
            }
        }
        fprintf(outfile, "\n");

        fprintf(outfile, "   => Significant Function Pairs Reverse <=\n\n");
        for (int M = 0; M < nbf_; M++) {
            for (int N = 0; N <= M; N++) {
                fprintf(outfile, "    %6ld = (%3d,%3d|\n", function_pairs_reverse_[M * (M + 1) / 2 + N], M, N);
            }
        }
        fprintf(outfile, "\n");

        fprintf(outfile, "   => Shell to Shell <=\n\n");
        for (int M = 0; M < nshell_; M++) {
            for (int N = 0; N < shell_to_shell_[M].size(); N++) {
                fprintf(outfile, "    (%3d, %3d|\n", M, shell_to_shell_[M][N]);
            }
        }
        fprintf(outfile, "\n");

        fprintf(outfile, "   => Function to Function <=\n\n");
        for (int M = 0; M < nbf_; M++) {
            for (int N = 0; N < function_to_function_[M].size(); N++) {
                fprintf(outfile, "    (%3d, %3d|\n", M, function_to_function_[M][N]);
            }
        }
        fprintf(outfile, "\n");

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
}

}
