#include "composite.h"

#include "psi4/libmints/integral.h"
#include "psi4/libmints/vector.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/lib3index/dftensor.h"
#include "psi4/libqt/qt.h"

#include <vector>
#include <algorithm>
#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#include "psi4/libpsi4util/process.h"
#endif

using namespace psi;

namespace psi {

DirectDFJ::DirectDFJ(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary, Options& options)
                        : SplitJKBase(primary, options), auxiliary_(auxiliary) {
    cutoff_ = options_.get_double("INTS_TOLERANCE");
    build_metric();
    build_ints();
}

void DirectDFJ::print_header() {
    if (print_) {
        outfile->Printf("  ==> Direct Density-Fitted J <==\n\n");

        outfile->Printf("    Primary Basis: %11s\n", primary_->name().c_str());
        outfile->Printf("    Auxiliary Basis: %11s\n", auxiliary_->name().c_str());
        outfile->Printf("    ERI Screening Cutoff: %11.0E\n", cutoff_);
        outfile->Printf("\n");
    }
}

void DirectDFJ::build_metric() {
    timer_on("DirectDFJ: Build Coulomb Metric");

    FittingMetric J_metric_obj(auxiliary_, true);
    J_metric_obj.form_fitting_metric();
    J_metric_ = J_metric_obj.get_metric();

    timer_off("DirectDFJ: Build Coulomb Metric");
}

void DirectDFJ::build_ints() {
    timer_on("DirectDFJ: Build Ints");
    ints_.resize(nthread_);
    auto zero = BasisSet::zero_ao_basis_set();
    IntegralFactory rifactory(auxiliary_, zero, primary_, primary_);
    ints_[0] = std::shared_ptr<TwoBodyAOInt>(rifactory.eri());
    for (int rank = 1; rank < nthread_; rank++) {
        ints_[rank] = std::shared_ptr<TwoBodyAOInt>(ints_[0]->clone());
    }

    timer_off("DirectDFJ: Build Ints");
}

void DirectDFJ::build_G_component(const std::vector<SharedMatrix>& D, std::vector<SharedMatrix>& J) {

    timer_on("DirectDFJ: J");

    timer_on("DirectDFJ: Setup");

    // => Zeroing the J matrix <= //
    for (const auto& Jmat : J) {
        Jmat->zero();
    }

    // => Sizing <= //
    int njk = D.size();
    int nbf = primary_->nbf();
    int nshell = primary_->nshell();
    int nbf_aux = auxiliary_->nbf();
    int nshell_aux = auxiliary_->nshell();

    // benchmarking 
    size_t nshellpair = ints_[0]->shell_pairs().size();
    size_t nshelltriplet = nshell_aux * nshellpair;
    size_t computed_triplets1 = 0, computed_triplets2 = 0;

    // screening threshold
    double thresh2 = options_.get_double("INTS_TOLERANCE") * options_.get_double("INTS_TOLERANCE");

    // per-thread G Vector buffers (for accumulating thread contributions to G)
    // G is the contraction of the density matrix with the 3-index ERIs
    std::vector<std::vector<SharedVector>> GT(njk, std::vector<SharedVector>(nthread_));

    // H is the contraction of G with the inverse coulomb metric
    std::vector<SharedVector> H(njk);

    // per-thread J Matrix buffers (for accumulating thread contributions to J)
    std::vector<std::vector<SharedMatrix>> JT(njk, std::vector<SharedMatrix>(nthread_));

    // initialize per-thread objects
    for(size_t jki = 0; jki < njk; jki++) {
        for(size_t thread = 0; thread < nthread_; thread++) {
            JT[jki][thread] = std::make_shared<Matrix>(nbf, nbf);
            GT[jki][thread] = std::make_shared<Vector>(nbf_aux);
        }
        H[jki] = std::make_shared<Vector>(nbf_aux);
    }

    // diagonal shell maxima of J_metric_ for screening
    std::vector<double> J_metric_shell_diag(nshell_aux, 0.0);
    for (size_t s = 0; s < nshell_aux; s++) {
        int bf_start = auxiliary_->shell(s).function_index();
        int bf_end = bf_start + auxiliary_->shell(s).nfunction();
        for (size_t bf = bf_start; bf < bf_end; bf++) {
            J_metric_shell_diag[s] = std::max(J_metric_shell_diag[s], J_metric_->get(bf, bf));
        }
    }

    // shell maxima of D for screening
    Matrix Dshell(nshell, nshell);
    auto Dshellp = Dshell.pointer();

    for(size_t M = 0; M < nshell; M++) {
        int nm = primary_->shell(M).nfunction();
        int mstart = primary_->shell(M).function_index();
        for(size_t N = 0; N < nshell; N++) {
            int nn = primary_->shell(N).nfunction();
            int nstart = primary_->shell(N).function_index();
            for(size_t jki = 0; jki < njk; jki++) {
                auto Dp = D[jki]->pointer();
                for(size_t m = mstart; m < mstart + nm; m++) {
                    for(size_t n = nstart; n < nstart + nn; n++) {
                        Dshellp[M][N] = std::max(Dshellp[M][N], std::abs(Dp[m][n]));
                    }
                }
            }
        }
    }

    timer_off("DirectDFJ: Setup");

    //  => First Contraction <= //

    // contract D with three-index DF ERIs to get G:
    // G_{p} = D_{mn} * (mn|p)

    timer_on("DirectDFJ: ERI1");
#pragma omp parallel for schedule(guided) num_threads(nthread_) reduction(+ : computed_triplets1)
    for (size_t MNP = 0; MNP < nshelltriplet; MNP++) {

        size_t MN = MNP % nshellpair;
        size_t P = MNP / nshellpair;
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        auto bra = ints_[rank]->shell_pairs()[MN];
        size_t M = bra.first;
        size_t N = bra.second;
        if(Dshellp[M][N] * Dshellp[M][N] * J_metric_shell_diag[P] * ints_[rank]->shell_pair_value(M,N) < thresh2) {
            continue;
        }
        computed_triplets1++;
        int np = auxiliary_->shell(P).nfunction();
        int pstart = auxiliary_->shell(P).function_index();
        int nm = primary_->shell(M).nfunction();
        int mstart = primary_->shell(M).function_index();
        int nn = primary_->shell(N).nfunction();
        int nstart = primary_->shell(N).function_index();
        ints_[rank]->compute_shell(P, 0, M, N);
        const auto & buffer = ints_[rank]->buffers()[0];

        for(size_t jki = 0; jki < njk; jki++) {

            auto GTp = GT[jki][rank]->pointer();
            auto Dp = D[jki]->pointer();

            for (int p = pstart, index = 0; p < pstart + np; p++) {
                for (int m = mstart; m < mstart + nm; m++) {
                    for (int n = nstart; n < nstart + nn; n++, index++) {
                        GTp[p] += buffer[index] * Dp[m][n];
                        if (N != M) GTp[p] += buffer[index] * Dp[n][m];
                    }
                }
            }

        }

    }

    timer_off("DirectDFJ: ERI1");

    //  => Second Contraction <= //

    //  linear solve for H:
    //  G_{p} = H_{q} (q|p)

    timer_on("DirectDFJ: Metric Solve");

    std::vector<int> ipiv(nbf_aux);

    for(size_t jki = 0; jki < njk; jki++) {
        for(size_t thread = 0; thread < nthread_; thread++) {
            H[jki]->add(GT[jki][thread]);
        }
        C_DGESV(nbf_aux, 1, J_metric_->clone()->pointer()[0], nbf_aux, ipiv.data(), H[jki]->pointer(), nbf_aux);
    }


    // I believe C_DSYSV should be faster than C_GESV, but I've found the opposite to be true.
    // This performance issue should be investigated, but is not consequential here.
    // The cost of either linear solve is dwarfed by the actual integral computation.
    //
    //std::vector<double> work(3 * nbf_aux);
    //int errcode = C_DSYSV('U', nbf_aux, 1, J_metric_->clone()->pointer()[0], nbf_aux, ipiv.data(), H[jki]->pointer(), nbf_aux, work.data(), 3 * nbf_aux);

    // shell maxima of H for screening
    Vector H_shell_max(nshell_aux);
    auto H_shell_maxp = H_shell_max.pointer();

    for(size_t jki = 0; jki < njk; jki++) {

        auto Hp = H[jki]->pointer();

        for (int P = 0; P < nshell_aux; P++) {
            int np = auxiliary_->shell(P).nfunction();
            int pstart = auxiliary_->shell(P).function_index();
            for (int p = pstart; p < pstart + np; p++) {
                H_shell_maxp[P] = std::max(H_shell_maxp[P], std::abs(Hp[p]));
            }
        }

    }

    timer_off("DirectDFJ: Metric Solve");

    //  => Third Contraction <= //

    // contract H with three-index DF ERIs to get J
    // J_{mn} = H_{p} (mn|p)

    timer_on("DirectDFJ: ERI2");

#pragma omp parallel for schedule(guided) num_threads(nthread_) reduction(+ : computed_triplets2)
    for (size_t MNP = 0; MNP < nshelltriplet; MNP++) {

        size_t MN = MNP % nshellpair;
        size_t P = MNP / nshellpair;
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        auto bra = ints_[rank]->shell_pairs()[MN];
        size_t M = bra.first;
        size_t N = bra.second;
        if(H_shell_maxp[P] * H_shell_maxp[P] * J_metric_shell_diag[P] * ints_[rank]->shell_pair_value(M,N) < thresh2) {
            continue;
        }
        computed_triplets2++;
        int np = auxiliary_->shell(P).nfunction();
        int pstart = auxiliary_->shell(P).function_index();
        int nm = primary_->shell(M).nfunction();
        int mstart = primary_->shell(M).function_index();
        int nn = primary_->shell(N).nfunction();
        int nstart = primary_->shell(N).function_index();

        ints_[rank]->compute_shell(P, 0, M, N);
        const auto & buffer = ints_[rank]->buffers()[0];

        for(size_t jki = 0; jki < njk; jki++) {

            auto JTp = JT[jki][rank]->pointer();
            auto Hp = H[jki]->pointer();

            for (int p = pstart, index = 0; p < pstart + np; p++) {
                for (int m = mstart; m < mstart + nm; m++) {
                    for (int n = nstart; n < nstart + nn; n++, index++) {
                        JTp[m][n] += buffer[index] * Hp[p];
                        if (N != M) JTp[n][m] += buffer[index] * Hp[p];
                    }
                }
            }

        }

    }

    timer_off("DirectDFJ: ERI2");

    if (bench_) {
        auto mode = std::ostream::app;
        PsiOutStream printer("bench.dat", mode);
        printer.Printf("DFJ ERI Shells: %zu,%zu,%zu\n", computed_triplets1, computed_triplets2, nshelltriplet);
    }

    for(size_t jki = 0; jki < njk; jki++) {
        for (size_t thread = 0; thread < nthread_; thread++) {
            J[jki]->add(JT[jki][thread]);
        }
        J[jki]->hermitivitize();
    }

    timer_off("DirectDFJ: J");
}

} // namespace psi