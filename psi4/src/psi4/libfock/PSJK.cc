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


#include "psi4/lib3index/3index.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/aiohandler.h"
#include "psi4/libqt/qt.h"
#include "psi4/psi4-dec.h"
#include "psi4/psifiles.h"
#include "psi4/libmints/sieve.h"
#include "psi4/libiwl/iwl.hpp"
#include "jk.h"
#include "jk_independent.h"
#include "link.h"
#include "direct_screening.h"
#include "cubature.h"
#include "points.h"

#include "psi4/lib3index/cholesky.h"

#include <sstream>
#include "psi4/libparallel/ParallelPrinter.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace psi;

namespace psi {
#if 0
PSJK::PSJK(std::shared_ptr<BasisSet> primary,
    Options& options) :
    JK(primary), options_(options)
{
    common_init();
}
PSJK::~PSJK()
{
}
void PSJK::common_init()
{
    df_ints_num_threads_ = 1;
    #ifdef _OPENMP
        df_ints_num_threads_ = Process::environment.get_n_threads();
    #endif
    unit_ = PSIF_DFSCF_BJ;
    theta_ = 0.3;
    psio_ = PSIO::shared_object();
    dealiasing_ = "QUADRATURE";
}
void PSJK::print_header() const
{
    if (print_) {
        outfile->Printf( "  ==> PSJK: Pseudospectral J/K Matrices <==\n\n");

        outfile->Printf( "    J tasked:          %11s\n", (do_J_ ? "Yes" : "No"));
        outfile->Printf( "    K tasked:          %11s\n", (do_K_ ? "Yes" : "No"));
        outfile->Printf( "    wK tasked:         %11s\n", (do_wK_ ? "Yes" : "No"));
        if (do_wK_)
            outfile->Printf( "    Omega:             %11.3E\n", omega_);
        outfile->Printf( "    OpenMP threads:    %11d\n", omp_nthread_);
        outfile->Printf( "    Integrals threads: %11d\n", df_ints_num_threads_);
        outfile->Printf( "    Memory (MB):       %11ld\n", (memory_ *8L) / (1024L * 1024L));
        outfile->Printf( "    Schwarz Cutoff:    %11.0E\n", cutoff_);
        outfile->Printf( "    Theta:             %11.3E\n", theta_);
        outfile->Printf( "    Dealiasing:        %11s\n", dealiasing_.c_str());
        outfile->Printf( "\n");

        outfile->Printf( "   => Quadrature Grid <=\n\n");
        outfile->Printf( "    Total Points:      %11d\n", grid_->rowspi()[0]);
        outfile->Printf( "\n");
        // TODO print grid algorithm details

        if (dealiasing_ == "DEALIAS") {
            outfile->Printf( "   => Dealias Basis Set <=\n\n");
            dealias_->print_by_level("outfile",print_);
        }
    }
}
void PSJK::preiterations()
{
    // PS requires constant sieve, must be static througout object life
    if (!sieve_) {
        sieve_ = std::shared_ptr<ERISieve>(new ERISieve(primary_, cutoff_));
    }

    build_QR();

    if (do_J_ || do_K_) {
        build_Amn_disk(theta_,"(A|mn) JK");
    }

    if (do_wK_) {
        throw FeatureNotImplemented("PSJK", "wK", __FILE__, __LINE__);
    }
}
void PSJK::postiterations()
{
    Q_.reset();
    R_.reset();
    grid_.reset();
}
void PSJK::compute_JK()
{
    // Short Range
    build_JK_SR();

    // Long Range
    build_JK_LR();


    // TODO wK
}
void PSJK::build_QR()
{
    std::shared_ptr<PseudospectralGrid> grid(new PseudospectralGrid(primary_->molecule(),
        primary_, options_));
    int npoints = grid->npoints();
    int nbf = primary_->nbf();
    int max_points = grid->max_points();
    int max_functions = grid->max_functions();

    // Grid
    grid_ = SharedMatrix(new Matrix("xyzw", npoints, 4));
    double** gridp = grid_->pointer();
    double* x = grid->x();
    double* y = grid->y();
    double* z = grid->z();
    double* w = grid->w();
    for (int P = 0; P < npoints; P++) {
        gridp[P][0] = x[P];
        gridp[P][1] = y[P];
        gridp[P][2] = z[P];
        gridp[P][3] = w[P];
    }

    // R (Collocation)
    R_ = SharedMatrix(new Matrix("R", nbf, npoints));
    double** Rp = R_->pointer();
    std::shared_ptr<BasisFunctions> points(new BasisFunctions(primary_, max_points, max_functions));
    const std::vector<std::shared_ptr<BlockOPoints> >& blocks = grid->blocks();
    int offset = 0;
    for (int index = 0; index < blocks.size(); index++) {
        points->compute_functions(blocks[index]);
        SharedMatrix phi = points->basis_value("PHI");
        double** phip = phi->pointer();
        const std::vector<int>& funmap = blocks[index]->functions_local_to_global();
        int nP = blocks[index]->npoints();

        for (int i = 0; i < funmap.size(); i++) {
            int iglobal = funmap[i];
            C_DCOPY(nP,&phip[0][i],max_functions,&Rp[iglobal][offset],1);
        }

        offset += nP;
    }

    points.reset();
    grid.reset();

    // Q (Quadrature Rule)
    if (dealiasing_ == "QUADRATURE") {
        for (int P = 0; P < npoints; P++) {
            C_DSCAL(nbf,sqrt(gridp[P][3]),&Rp[0][P],npoints);
        }
        Q_ = R_;
    } else if (dealiasing_ == "RENORMALIZED") {

        for (int P = 0; P < npoints; P++) {
            C_DSCAL(nbf,sqrt(gridp[P][3]),&Rp[0][P],npoints);
        }

        bool warning;
        SharedMatrix Rplus = R_->pseudoinverse(1.0E-10, &warning);
        if (warning) {
            outfile->Printf( "    Warning, Renormalization had to be conditioned.\n\n");
        }

        std::shared_ptr<IntegralFactory> factory(new IntegralFactory(primary_));
        std::shared_ptr<OneBodyAOInt> ints(factory->ao_overlap());
        SharedMatrix S(new Matrix("S", nbf, nbf));
        ints->compute(S);
        ints.reset();

        Q_ = SharedMatrix(R_->clone());
        Q_->set_name("Q");

        double** Qp = Q_->pointer();
        double** Rp = Rplus->pointer();
        double** Sp = S->pointer();

        C_DGEMM('N','N',nbf,npoints,nbf,1.0,Sp[0],nbf,Rp[0],npoints,0.0,Qp[0],npoints);

    } else if (dealiasing_ == "DEALIAS") {
        // TODO
        throw FeatureNotImplemented("PSJK", "Dealiasing", __FILE__, __LINE__);
    }
}
void PSJK::build_Amn_disk(double theta, const std::string& entry)
{
    std::shared_ptr<IntegralFactory> factory(new IntegralFactory(primary_));
    std::vector<std::shared_ptr<PseudospectralInt> > ints;
    for (int i = 0; i < df_ints_num_threads_; i++) {
        ints.push_back(std::shared_ptr<PseudospectralInt>(static_cast<PseudospectralInt*>(factory->ao_pseudospectral())));
        ints[i]->set_omega(theta);
    }

    int maxrows = max_rows();
    int ntri  = sieve_->function_pairs().size();
    int naux  = grid_->rowspi()[0];

    const std::vector<long int>& function_pairs_r = sieve_->function_pairs_reverse();
    const std::vector<std::pair<int,int> >& shell_pairs = sieve_->shell_pairs();

    SharedMatrix Amn(new Matrix("Amn", maxrows, ntri));
    double** Amnp = Amn->pointer();
    double** Gp = grid_->pointer();

    psio_->open(unit_,PSIO_OPEN_OLD);
    psio_address addr = PSIO_ZERO;

    for (int Pstart = 0; Pstart < naux; Pstart += maxrows) {
        int nrows =  (Pstart + maxrows >= naux ? naux - Pstart : maxrows);

        #pragma omp parallel for num_threads(df_ints_num_threads_)
        for (int row = 0; row < (unsigned int)nrows; row++) {
            int thread = 0;
            #ifdef _OPENMP
                thread = omp_get_thread_num();
            #endif
            std::shared_ptr<PseudospectralInt> eri = ints[thread];
            const double* buffer = eri->buffer();
            eri->set_point(Gp[row + Pstart][0], Gp[row + Pstart][1], Gp[row + Pstart][2]);
            double* Amp = Amnp[row];

            for (int ind = 0; ind < shell_pairs.size(); ind++) {
                int M = shell_pairs[ind].first;
                int N = shell_pairs[ind].second;

                eri->compute_shell(M,N);

                int nM = primary_->shell(M).nfunction();
                int nN = primary_->shell(N).nfunction();
                int oM = primary_->shell(M).function_index();
                int oN = primary_->shell(N).function_index();

                for (int m = 0 ; m < nM; m++) {
                    for (int n = 0; n < nN; n++) {
                        int am = m + oM;
                        int an = n + oN;
                        if (am >= an) {
                            int mn = (am *(am + 1) >> 1) + an;
                            long int mn_local = function_pairs_r[mn];
                            if (mn >= 0) {
                                Amp[mn_local] = buffer[m * nN + n];
                            }
                        }
                    }
                }
            }
        }
        psio_->write(unit_, entry.c_str(), (char*) Amnp[0], sizeof(double) * ntri * nrows, addr, &addr);
    }

    psio_->close(unit_, 1);
}
int PSJK::max_rows()
{
    int ntri = sieve_->function_pairs().size();
    int naux = grid_->rowspi()[0];
    unsigned long int effective_memory = memory_ - memory_overhead();
    unsigned long int rows = effective_memory / ntri;
    rows = (rows > naux ? naux : rows);
    rows = (rows < 1 ? 1 : rows);
    return (int) rows;
}
void PSJK::block_J(double** Amnp, int Pstart, int nP, const std::vector<SharedMatrix>& J)
{
    int nbf = primary_->nbf();
    int npoints = grid_->rowspi()[0];
    const std::vector<std::pair<int,int> >& funmap = sieve_->function_pairs();
    int ntri = funmap.size();

    double** Qp = Q_->pointer();
    double** Rp = R_->pointer();
    double** Vp = V_->pointer();
    double*  dp = d_->pointer();
    double* J2p = J_temp_->pointer();

    for (int N = 0; N < D_ao_.size(); N++) {

        double** Dp = D_ao_[N]->pointer();
        double** Jp = J[N]->pointer();

        C_DGEMM('T','N',nP,nbf,nbf,1.0,&Qp[0][Pstart],npoints,Dp[0],nbf,0.0,Vp[0],nbf);
        for (int P = 0; P < nP; P++) {
            dp[P] = C_DDOT(nbf,&Rp[0][Pstart + P],npoints,Vp[P],1);
        }

        C_DGEMV('T',nP,ntri,1.0,Amnp[0],ntri,dp,1,0.0,J2p,1);

        for (long int mn = 0; mn < ntri; mn++) {
            int m = funmap[mn].first;
            int n = funmap[mn].second;
            Jp[m][n] += J2p[mn];
            Jp[n][m] += (m == n ? 0.0 : J2p[mn]);
        }
    }
}
void PSJK::block_K(double** Amnp, int Pstart, int nP, const std::vector<SharedMatrix>& K)
{
    int nbf = primary_->nbf();
    int npoints = grid_->rowspi()[0];
    const std::vector<std::pair<int,int> >& funmap = sieve_->function_pairs();
    int ntri = funmap.size();

    double** Qp = Q_->pointer();
    double** Rp = R_->pointer();
    double** Vp = V_->pointer();
    double** Wp = W_->pointer();

    for (int N = 0; N < D_ao_.size(); N++) {

        double** Dp = D_ao_[N]->pointer();
        double** Kp = K[N]->pointer();

        C_DGEMM('T','N',nP,nbf,nbf,1.0,&Qp[0][Pstart],npoints,Dp[0],nbf,0.0,Vp[0],nbf);

        W_->zero();
        #pragma omp parallel for
        for (int P = 0; P < nP; P++) {
            double* Arp = Amnp[P];
            double* Wrp = Wp[P];
            double* Vrp = Vp[P];
            for (long int mn = 0; mn < ntri; mn++) {
                int m = funmap[mn].first;
                int n = funmap[mn].second;
                double Aval = Arp[mn];
                Wrp[m] += Aval * Vrp[n];
                Wrp[n] += (m == n ? 0.0 : Aval * Vrp[m]);
            }
        }

        C_DGEMM('N','N',nbf,nbf,nP,1.0,&Rp[0][Pstart],npoints,Wp[0],nbf,1.0,Kp[0],nbf);
    }
}
void PSJK::build_JK_SR()
{
    std::shared_ptr<IntegralFactory> factory(new IntegralFactory(primary_));
    std::shared_ptr<TwoBodyAOInt> eri(factory->erf_complement_eri(theta_));
    const double* buffer = eri->buffer();

    double cutoff = sieve_->sieve();
    const std::vector<std::pair<int,int> >& shellmap = sieve_->shell_pairs();
    int nTRI = shellmap.size();
    for (int MN = 0; MN < nTRI; MN++) {
        int M = shellmap[MN].first;
        int N = shellmap[MN].second;
        int nM = primary_->shell(M).nfunction();
        int nN = primary_->shell(N).nfunction();
        int oM = primary_->shell(M).function_index();
        int oN = primary_->shell(N).function_index();
        for (int LS = 0; LS < nTRI; LS++) {
            int L = shellmap[LS].first;
            int S = shellmap[LS].second;
            if (!sieve_->shell_significant(M,N,L,S)) continue;
            // TODO: More sieving
            int nL = primary_->shell(L).nfunction();
            int nS = primary_->shell(S).nfunction();
            int oL = primary_->shell(L).function_index();
            int oS = primary_->shell(S).function_index();

            eri->compute_shell(M,N,L,S);

            for (int rM = 0, index = 0; rM < nM; rM++) {
                for (int rN = 0; rN < nN; rN++) {
                    for (int rL = 0; rL < nL; rL++) {
                        for (int rS = 0; rS < nS; rS++, index++) {
                            int m = rM + oM;
                            int n = rN + oN;
                            int l = rL + oL;
                            int s = rS + oS;
                            int mn = (m * (m + 1) >> 1) + n;
                            int ls = (l * (l + 1) >> 1) + s;
                            if (n > m || s > l || ls > mn) continue;
                            double val = buffer[index];

                            // J
                            if (do_J_) {
                                for (int A = 0; A < D_.size(); A++) {
                                    double** Jp = J_ao_[A]->pointer();
                                    double** Dp = D_ao_[A]->pointer();
                                    if (mn == ls && m == n) {
                                        // (mm|mm)
                                        Jp[m][m] += val * Dp[m][m]; // (mm|mm)
                                    } else if (m == n && l == s) {
                                        // (mm|ll)
                                        Jp[m][m] += val * Dp[l][l]; // (mm|ll)
                                        Jp[l][l] += val * Dp[m][m]; // (ll|mm)
                                    } else if (m == n) {
                                        // (mm|ls)
                                        Jp[m][m] += val * Dp[l][s]; // (mm|ls)
                                        Jp[m][m] += val * Dp[s][l]; // (mm|sl)
                                        Jp[l][s] += val * Dp[m][m]; // (ls|mm)
                                        Jp[s][l] += val * Dp[m][m]; // (sl|mm)
                                    } else if (l == s) {
                                        // (mn|ll)
                                        Jp[m][n] += val * Dp[l][l]; // (mn|ll)
                                        Jp[n][m] += val * Dp[l][l]; // (nm|ll)
                                        Jp[l][l] += val * Dp[m][n]; // (ll|mn)
                                        Jp[l][l] += val * Dp[n][m]; // (ll|nm)
                                    } else if (mn == ls) {
                                        // (mn|mn)
                                        Jp[m][n] += val * Dp[m][n]; // (mn|mn)
                                        Jp[m][n] += val * Dp[n][m]; // (mn|nm)
                                        Jp[n][m] += val * Dp[m][n]; // (nm|mn)
                                        Jp[n][m] += val * Dp[n][m]; // (nm|nm)
                                    } else {
                                        // (mn|ls)
                                        Jp[m][n] += val * Dp[l][s]; // (mn|ls)
                                        Jp[m][n] += val * Dp[s][l]; // (mn|sl)
                                        Jp[n][m] += val * Dp[l][s]; // (nm|ls)
                                        Jp[n][m] += val * Dp[s][l]; // (nm|sl)
                                        Jp[l][s] += val * Dp[m][n]; // (ls|mn)
                                        Jp[l][s] += val * Dp[n][m]; // (ls|nm)
                                        Jp[s][l] += val * Dp[m][n]; // (sl|mn)
                                        Jp[s][l] += val * Dp[n][m]; // (sl|nm)
                                    }
                                }
                            }
                            // K
                            if (do_K_) {
                                for (int A = 0; A < D_.size(); A++) {
                                    double** Kp = K_ao_[A]->pointer();
                                    double** Dp = D_ao_[A]->pointer();
                                    if (mn == ls && m == n) {
                                        // (mm|mm)
                                        Kp[m][m] += val * Dp[m][m]; // (mm|mm)
                                    } else if (m == n && l == s) {
                                        // (mm|ll)
                                        Kp[m][l] += val * Dp[m][l]; // (mm|ll)
                                        Kp[l][m] += val * Dp[l][m]; // (ll|mm)
                                    } else if (m == n) {
                                        // (mm|ls)
                                        Kp[m][s] += val * Dp[m][l]; // (mm|ls)
                                        Kp[m][l] += val * Dp[m][s]; // (mm|sl)
                                        Kp[l][m] += val * Dp[s][m]; // (ls|mm)
                                        Kp[s][m] += val * Dp[l][m]; // (sl|mm)
                                    } else if (l == s) {
                                        // (mn|ll)
                                        Kp[m][l] += val * Dp[n][l]; // (mn|ll)
                                        Kp[n][l] += val * Dp[m][l]; // (nm|ll)
                                        Kp[l][n] += val * Dp[l][m]; // (ll|mn)
                                        Kp[l][m] += val * Dp[l][n]; // (ll|nm)
                                    } else if (mn == ls) {
                                        // (mn|mn)
                                        Kp[m][n] += val * Dp[n][m]; // (mn|mn)
                                        Kp[m][m] += val * Dp[n][n]; // (mn|nm)
                                        Kp[n][n] += val * Dp[m][m]; // (nm|mn)
                                        Kp[n][m] += val * Dp[m][n]; // (nm|nm)
                                    } else {
                                        // (mn|ls)
                                        Kp[m][s] += val * Dp[n][l]; // (mn|ls)
                                        Kp[m][l] += val * Dp[n][s]; // (mn|sl)
                                        Kp[n][s] += val * Dp[m][l]; // (nm|ls)
                                        Kp[n][l] += val * Dp[m][s]; // (nm|sl)
                                        Kp[l][n] += val * Dp[s][m]; // (ls|mn)
                                        Kp[l][m] += val * Dp[s][n]; // (ls|nm)
                                        Kp[s][n] += val * Dp[l][m]; // (sl|mn)
                                        Kp[s][m] += val * Dp[l][n]; // (sl|nm)
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
void PSJK::build_JK_LR()
{
    int nbf  = primary_->nbf();
    int naux = grid_->rowspi()[0];
    int ntri = sieve_->function_pairs().size();
    int maxrows = max_rows();

    SharedMatrix Amn(new Matrix("Amn",maxrows,ntri));
    double** Amnp = Amn->pointer();

    V_ = SharedMatrix(new Matrix("V", maxrows,nbf));
    W_ = SharedMatrix(new Matrix("W", maxrows,nbf));
    d_ = SharedVector(new Vector("d", maxrows));;
    J_temp_ = SharedVector(new Vector("J temp", ntri));

    if (do_J_ || do_K_) {
        psio_->open(unit_, PSIO_OPEN_OLD);
        psio_address addr = PSIO_ZERO;
        for (int Pstart = 0; Pstart < naux; Pstart += maxrows) {
            int nrows = (naux - Pstart <= maxrows ? naux - Pstart : maxrows);
            psio_->read(unit_,"(A|mn) JK", (char*) Amnp[0], sizeof(double)*nrows*ntri,addr,&addr);
            if (do_J_) {
                block_J(Amnp,Pstart,nrows,J_ao_);
            }
            if (do_K_) {
                block_K(Amnp,Pstart,nrows,K_ao_);
            }
        }
        psio_->close(unit_, 1);
        if (lr_symmetric_ && do_K_) {
            for (int N = 0; N < D_.size(); N++) {
                K_ao_[N]->hermitivitize();
            }
        }
    }

    Amn.reset();
    V_.reset();
    W_.reset();
    d_.reset();
    J_temp_.reset();
}
void PSJK::build_JK_debug(const std::string& op, double theta)
{
    if (do_J_) {
        for (int A = 0; A < D_.size(); A++) {
            J_ao_[A]->zero();
        }
    }
    if (do_K_) {
        for (int A = 0; A < D_.size(); A++) {
            K_ao_[A]->zero();
        }
    }

    std::shared_ptr<IntegralFactory> factory(new IntegralFactory(primary_));
    std::shared_ptr<TwoBodyAOInt> eri;
    if (op == "") {
        eri = std::shared_ptr<TwoBodyAOInt>(factory->eri());
    } else if (op == "SR") {
        eri = std::shared_ptr<TwoBodyAOInt>(factory->erf_complement_eri(theta));
    } else if (op == "LR") {
        eri = std::shared_ptr<TwoBodyAOInt>(factory->erf_eri(theta));
    } else {
        throw PSIEXCEPTION("What is this?");
    }
    const double* buffer = eri->buffer();

    double cutoff = sieve_->sieve();
    const std::vector<std::pair<int,int> >& shellmap = sieve_->shell_pairs();
    int nTRI = shellmap.size();
    for (int MN = 0; MN < nTRI; MN++) {
        int M = shellmap[MN].first;
        int N = shellmap[MN].second;
        int nM = primary_->shell(M).nfunction();
        int nN = primary_->shell(N).nfunction();
        int oM = primary_->shell(M).function_index();
        int oN = primary_->shell(N).function_index();
        for (int LS = 0; LS < nTRI; LS++) {
            int L = shellmap[LS].first;
            int S = shellmap[LS].second;
            if (!sieve_->shell_significant(M,N,L,S)) continue;
            // TODO: More sieving
            int nL = primary_->shell(L).nfunction();
            int nS = primary_->shell(S).nfunction();
            int oL = primary_->shell(L).function_index();
            int oS = primary_->shell(S).function_index();

            eri->compute_shell(M,N,L,S);

            for (int rM = 0, index = 0; rM < nM; rM++) {
                for (int rN = 0; rN < nN; rN++) {
                    for (int rL = 0; rL < nL; rL++) {
                        for (int rS = 0; rS < nS; rS++, index++) {
                            int m = rM + oM;
                            int n = rN + oN;
                            int l = rL + oL;
                            int s = rS + oS;
                            int mn = (m * (m + 1) >> 1) + n;
                            int ls = (l * (l + 1) >> 1) + s;
                            if (n > m || s > l || ls > mn) continue;
                            double val = buffer[index];

                            // J
                            if (do_J_) {
                                for (int A = 0; A < D_.size(); A++) {
                                    double** Jp = J_ao_[A]->pointer();
                                    double** Dp = D_ao_[A]->pointer();
                                    if (mn == ls && m == n) {
                                        // (mm|mm)
                                        Jp[m][m] += val * Dp[m][m]; // (mm|mm)
                                    } else if (m == n && l == s) {
                                        // (mm|ll)
                                        Jp[m][m] += val * Dp[l][l]; // (mm|ll)
                                        Jp[l][l] += val * Dp[m][m]; // (ll|mm)
                                    } else if (m == n) {
                                        // (mm|ls)
                                        Jp[m][m] += val * Dp[l][s]; // (mm|ls)
                                        Jp[m][m] += val * Dp[s][l]; // (mm|sl)
                                        Jp[l][s] += val * Dp[m][m]; // (ls|mm)
                                        Jp[s][l] += val * Dp[m][m]; // (sl|mm)
                                    } else if (l == s) {
                                        // (mn|ll)
                                        Jp[m][n] += val * Dp[l][l]; // (mn|ll)
                                        Jp[n][m] += val * Dp[l][l]; // (nm|ll)
                                        Jp[l][l] += val * Dp[m][n]; // (ll|mn)
                                        Jp[l][l] += val * Dp[n][m]; // (ll|nm)
                                    } else if (mn == ls) {
                                        // (mn|mn)
                                        Jp[m][n] += val * Dp[m][n]; // (mn|mn)
                                        Jp[m][n] += val * Dp[n][m]; // (mn|nm)
                                        Jp[n][m] += val * Dp[m][n]; // (nm|mn)
                                        Jp[n][m] += val * Dp[n][m]; // (nm|nm)
                                    } else {
                                        // (mn|ls)
                                        Jp[m][n] += val * Dp[l][s]; // (mn|ls)
                                        Jp[m][n] += val * Dp[s][l]; // (mn|sl)
                                        Jp[n][m] += val * Dp[l][s]; // (nm|ls)
                                        Jp[n][m] += val * Dp[s][l]; // (nm|sl)
                                        Jp[l][s] += val * Dp[m][n]; // (ls|mn)
                                        Jp[l][s] += val * Dp[n][m]; // (ls|nm)
                                        Jp[s][l] += val * Dp[m][n]; // (sl|mn)
                                        Jp[s][l] += val * Dp[n][m]; // (sl|nm)
                                    }
                                }
                            }
                            // K
                            if (do_K_) {
                                for (int A = 0; A < D_.size(); A++) {
                                    double** Kp = K_ao_[A]->pointer();
                                    double** Dp = D_ao_[A]->pointer();
                                    if (mn == ls && m == n) {
                                        // (mm|mm)
                                        Kp[m][m] += val * Dp[m][m]; // (mm|mm)
                                    } else if (m == n && l == s) {
                                        // (mm|ll)
                                        Kp[m][l] += val * Dp[m][l]; // (mm|ll)
                                        Kp[l][m] += val * Dp[l][m]; // (ll|mm)
                                    } else if (m == n) {
                                        // (mm|ls)
                                        Kp[m][s] += val * Dp[m][l]; // (mm|ls)
                                        Kp[m][l] += val * Dp[m][s]; // (mm|sl)
                                        Kp[l][m] += val * Dp[s][m]; // (ls|mm)
                                        Kp[s][m] += val * Dp[l][m]; // (sl|mm)
                                    } else if (l == s) {
                                        // (mn|ll)
                                        Kp[m][l] += val * Dp[n][l]; // (mn|ll)
                                        Kp[n][l] += val * Dp[m][l]; // (nm|ll)
                                        Kp[l][n] += val * Dp[l][m]; // (ll|mn)
                                        Kp[l][m] += val * Dp[l][n]; // (ll|nm)
                                    } else if (mn == ls) {
                                        // (mn|mn)
                                        Kp[m][n] += val * Dp[n][m]; // (mn|mn)
                                        Kp[m][m] += val * Dp[n][n]; // (mn|nm)
                                        Kp[n][n] += val * Dp[m][m]; // (nm|mn)
                                        Kp[n][m] += val * Dp[m][n]; // (nm|nm)
                                    } else {
                                        // (mn|ls)
                                        Kp[m][s] += val * Dp[n][l]; // (mn|ls)
                                        Kp[m][l] += val * Dp[n][s]; // (mn|sl)
                                        Kp[n][s] += val * Dp[m][l]; // (nm|ls)
                                        Kp[n][l] += val * Dp[m][s]; // (nm|sl)
                                        Kp[l][n] += val * Dp[s][m]; // (ls|mn)
                                        Kp[l][m] += val * Dp[s][n]; // (ls|nm)
                                        Kp[s][n] += val * Dp[l][m]; // (sl|mn)
                                        Kp[s][m] += val * Dp[l][n]; // (sl|nm)
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    outfile->Printf( "  ==> JK Debug %s, theta = %11.3E <==\n\n", op.c_str(), theta);
    for (int A = 0; A < D_.size(); A++) {
        if (do_J_) {
            J_ao_[A]->print();
        }
        if (do_K_) {
            K_ao_[A]->print();
        }
    }
}
#endif
}
