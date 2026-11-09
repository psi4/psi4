/* 
 *  @BEGIN LICENSE
 * 
 *  Hilbert: a space for quantum chemistry plugins to Psi4 
 * 
 *  Copyright (c) 2020 by its authors (LICENSE).
 * 
 *  The copyrights for code used from other parties are included in
 *  the corresponding files.
 * 
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 * 
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program.  If not, see http://www.gnu.org/licenses/.
 * 
 *  @END LICENSE
 */

#include <psi4/psi4-dec.h>
#include <psi4/physconst.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libpsio/psio.hpp>

#include <psi4/libmints/molecule.h>
#include <psi4/libmints/wavefunction.h>
#include <psi4/libmints/mintshelper.h>
#include <psi4/libmints/matrix.h>
#include <psi4/libmints/factory.h>
#include <psi4/libmints/vector.h>
#include <psi4/libmints/basisset.h>
#include <psi4/lib3index/dftensor.h>
#include <psi4/libqt/qt.h>

// jk object
#include <psi4/libfock/jk.h>

// for dft
#include "psi4/libfock/v.h"
#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libscf_solver/hf.h"

#include <psi4/psifiles.h>
#include <psi4/libtrans/integraltransform.h>

#include <algorithm>

// plugin headers
#include "rks_response.h"
#include "diis.h"

namespace psi { namespace scf_response{ 

RKSResponseSolver::RKSResponseSolver(std::shared_ptr<Wavefunction> reference_wavefunction, Options& options_):
    KSResponseSolver(reference_wavefunction, options_) {
    reference_wavefunction_ = reference_wavefunction;

    // apparently compute_Vx wants me to set the density, and it should 
    // be done here to respect the difference between rks and uks
    potential_->set_D({Da_});
}

RKSResponseSolver::~RKSResponseSolver() {
}

std::vector<std::vector<double>> RKSResponseSolver::first_order_response(std::vector<std::shared_ptr<Matrix>> op_a, std::vector<std::shared_ptr<Matrix>> op_b, double omega) {

    double d_convergence = options_.get_double("D_CONVERGENCE");
    int maxiter = options_.get_int("MAXITER");

    outfile->Printf("\n");
    outfile->Printf("\n");
    outfile->Printf( "        *******************************************\n");
    outfile->Printf( "        *                                         *\n");
    outfile->Printf( "        *                                         *\n");
    outfile->Printf( "        *    Unrestricted TDDFT Response          *\n");
    outfile->Printf( "        *                                         *\n");
    outfile->Printf( "        *                                         *\n");
    outfile->Printf( "        *******************************************\n");
    outfile->Printf("\n");

    outfile->Printf("\n");
    outfile->Printf("    No. basis functions:            %5i\n",nso_);
    outfile->Printf("    No. alpha electrons:            %5i\n",nalpha_);
    outfile->Printf("    No. beta electrons:             %5i\n",nbeta_);
    outfile->Printf("    Omega:           %20.12lf\n",omega);
    outfile->Printf("    d_convergence:             %10.3le\n",d_convergence);
    outfile->Printf("    maxiter:                        %5i\n",maxiter);
    outfile->Printf("\n");

    // dimension of the problem
    int oa = nalpha_;
    int va = nso_ - oa;

    int N = oa*va;

    int X_dim = 3*N;
    int Y_dim = 3*N;

    // are we doing static properties? if so, no Y necessary
    bool have_Y = true;
    if ( fabs(omega) < 1e-14 ) {
        have_Y = false;
        Y_dim = 0;
    }

    double * X_and_Y = (double*)malloc((X_dim + Y_dim)*sizeof(double));
    double * X = X_and_Y;
    double * Y = X_and_Y + have_Y * X_dim;

    memset((void*)X, '\0', X_dim * sizeof(double));
    memset((void*)Y, '\0', Y_dim * sizeof(double));

    double * old_X_and_Y = (double*)malloc((X_dim + Y_dim)*sizeof(double));
    double * old_X = old_X_and_Y;
    double * old_Y = old_X_and_Y + have_Y * X_dim;

    memset((void*)old_X, '\0', X_dim * sizeof(double));
    memset((void*)old_Y, '\0', Y_dim * sizeof(double));

    double * X_and_Y_error = (double*)malloc((X_dim + Y_dim)*sizeof(double));
    double * X_error = X_and_Y_error;
    double * Y_error = X_and_Y_error + have_Y * X_dim;

    memset((void*)X_error, '\0', X_dim * sizeof(double));
    memset((void*)Y_error, '\0', Y_dim * sizeof(double));

    double * AB_X = (double*)malloc(3*nmo_*nmo_*sizeof(double));
    memset((void*)AB_X, '\0', 3*nmo_*nmo_ * sizeof(double));
    double * AB_Y = AB_X;
    if ( have_Y ) {
        AB_Y = (double*)malloc(3*nmo_*nmo_*sizeof(double));
        memset((void*)AB_Y, '\0', 3*nmo_*nmo_ * sizeof(double));
    }

    double * ea = epsilon_a_->pointer();
    double * eb = epsilon_b_->pointer();

    std::shared_ptr<DIIS> diis (new DIIS(X_dim + Y_dim));

    outfile->Printf("\n");
    outfile->Printf("    ==> First-Order Response Equations <==\n");
    outfile->Printf("\n");

    int iter = 0;
    double err = 0.0;
    outfile->Printf("    ");
    outfile->Printf(" Iter ");
    outfile->Printf("                dX/Y ");
    outfile->Printf("\n");
    do {
        double damp = 0.5;
        if (iter > 10) {
            damp = 0.0;
        }
        C_DCOPY(3*N, X, 1, old_X, 1);
        if ( have_Y ) {
            C_DCOPY(3*N, Y, 1, old_Y, 1);
        }
        for (int p = 0; p < 3; p++) {

            build_Au_Bu(N, 1, &X[p*N], &AB_X[p*nmo_*nmo_]);

            if ( have_Y ) {
                build_Au_Bu(N, 1, &Y[p*N], &AB_Y[p*nmo_*nmo_]);
            }

            // X, alpha
            for (int i = 0; i < oa; i++) {
                for (int a = 0; a < va; a++) {
                    double precon = 1.0 / (ea[a+oa] - ea[i] + omega);
                    int ia = i * va + a;
                    X[p*N+ia] = damp * X[p*N+ia] 
                              + (1.0 - damp) * precon * (-op_a[p]->pointer()[i][a+oa] 
                                                       - AB_X[p*nmo_*nmo_ + i*nmo_+(a+oa)] 
                                                       - AB_Y[p*nmo_*nmo_ + (a+oa)*nmo_+i] );
                }
            }
            if ( have_Y ) {
                // Y, alpha
                for (int i = 0; i < oa; i++) {
                    for (int a = 0; a < va; a++) {
                        double precon = 1.0 / (ea[a+oa] - ea[i] - omega);
                        int ia = i * va + a;
                        Y[p*N+ia] = damp * Y[p*N+ia] 
                                  + (1.0 - damp) * precon * (-op_a[p]->pointer()[i][a+oa] 
                                                           - AB_Y[p*nmo_*nmo_ + i*nmo_+(a+oa)] 
                                                           - AB_X[p*nmo_*nmo_ + (a+oa)*nmo_+i] );
                    }
                }
            }

            C_DCOPY(N, &X[p*N], 1, &X_error[p*N], 1);
            C_DAXPY(N, -1.0, &old_X[p*N], 1, &X_error[p*N], 1);

            if ( have_Y ) {
                C_DCOPY(N, &Y[p*N], 1, &Y_error[p*N], 1);
                C_DAXPY(N, -1.0, &old_Y[p*N], 1, &Y_error[p*N], 1);
            }
        }

        err = C_DDOT(N, X_error, 1, X_error, 1);
        if ( have_Y ) {
            err += C_DDOT(N, Y_error, 1, Y_error, 1);
        }
        err = sqrt(err);

        // DIIS extrapolation
        diis->WriteVector(X_and_Y);
        diis->WriteErrorVector(X_and_Y_error);
        diis->Extrapolate(X_and_Y);

        outfile->Printf("    %5i %20.12lf\n", iter, err);
        if ( err < d_convergence ) {
            break;
        }

        iter++;
        if ( iter > maxiter ) break;

    }while(err > options_.get_double("D_CONVERGENCE"));

    outfile->Printf("\n");
    if ( iter > maxiter && options_.get_bool("FAIL_ON_MAXITER") ){
        throw PsiException("Maximum number of iterations exceeded!",__FILE__,__LINE__);
    } else if ( iter > maxiter ) {
        outfile->Printf("    First-order response equations did not converge!\n");
    } else {
        outfile->Printf("    First-order response equations converged!\n");
    }
    outfile->Printf("\n");
    
    for (int p = 0; p < 3; p++) {
        build_Au_Bu(N, 1, &X[p*N], &AB_X[p*nmo_*nmo_]);
        if ( have_Y ) {
            build_Au_Bu(N, 1, &Y[p*N], &AB_Y[p*nmo_*nmo_]);
        }
    }

    free(old_X_and_Y);
    free(X_and_Y_error);

    std::vector<double> ret_X(X, X + X_dim); 
    std::vector<double> ret_Y(Y, Y + X_dim); // note this is not Y_dim because Y_dim could be zero

    std::vector<std::vector<double>> ret_X_and_Y;
    ret_X_and_Y.push_back(ret_X);
    ret_X_and_Y.push_back(ret_Y);

    free(X_and_Y);
    free(AB_X);
    if ( have_Y ) {
        free(AB_Y);
    }

    return ret_X_and_Y;
}

void RKSResponseSolver::compute_polarizability(std::vector<double>X, std::vector<double>Y, double omega){

    outfile->Printf("\n");
    outfile->Printf("    ==> Polarizability <==\n");
    outfile->Printf("\n");
    outfile->Printf("    omega: %12.8lf\n", omega);
    outfile->Printf("\n");
    outfile->Printf("\n");

    // dimension of the problem
    int oa = nalpha_;
    int va = nso_ - oa;

    int N = oa*va;

    std::vector<std::string> dir {"X", "Y", "Z"};
    for (int p = 0; p < 3; p++) {
        for (int q = p; q < 3; q++) {
            double alpha = 0.0;
            for (int i = 0; i < oa; i++) {
                for (int a = 0; a < va; a++) {
                    int ia = i * va + a;
                    alpha += mua_[p]->pointer()[i][a+oa] * X[q*N+ia];
                    alpha += mua_[p]->pointer()[i][a+oa] * Y[q*N+ia];
                }
            }
            alpha *= 2.0;
            outfile->Printf("    ALPHA(%s%s) %20.12lf\n", dir[p].c_str(), dir[q].c_str(), alpha);

            // add polarizabilities to psi variables
            std::string label = "ALPHA(";
            label += dir[p] + dir[q] + ")";
            std::transform(label.begin(), label.end(), label.begin(),
                [](unsigned char c) { return std::toupper(c); });
            Process::environment.globals[label] = alpha;
        }
    }
}

void RKSResponseSolver::compute_hyperpolarizability(

    std::vector<std::vector<double>>amps_wx, 
    std::vector<std::vector<double>>amps_wy, 
    std::vector<std::vector<double>>amps_wz,
    std::string type, double omega) {

    if ( !is_hf_ ) {
        throw PsiException("DFT hyperpolarizbalities are not supported", __FILE__, __LINE__);
    }

    // dimension of the problem
    int oa = nalpha_;
    int va = nso_ - oa;

    int N = oa*va;

    outfile->Printf("\n");
    outfile->Printf("    ==>  Hyperpolarizability <==\n");
    outfile->Printf("\n");
    outfile->Printf("    type:  %12s\n", type.c_str());
    outfile->Printf("    omega: %12.8lf\n", omega);
    outfile->Printf("\n");

    int L = 1;

    double * AX_wx = (double*)malloc(3*2*nmo_*nmo_*sizeof(double));
    double * AY_wx = (double*)malloc(3*2*nmo_*nmo_*sizeof(double));
    double * AX_wy = (double*)malloc(3*2*nmo_*nmo_*sizeof(double));
    double * AY_wy = (double*)malloc(3*2*nmo_*nmo_*sizeof(double));
    double * AX_wz = (double*)malloc(3*2*nmo_*nmo_*sizeof(double));
    double * AY_wz = (double*)malloc(3*2*nmo_*nmo_*sizeof(double));

    memset((void*)AX_wx, '\0', 3*2*nmo_*nmo_ * sizeof(double));
    memset((void*)AY_wx, '\0', 3*2*nmo_*nmo_ * sizeof(double));
    memset((void*)AX_wy, '\0', 3*2*nmo_*nmo_ * sizeof(double));
    memset((void*)AY_wy, '\0', 3*2*nmo_*nmo_ * sizeof(double));
    memset((void*)AX_wz, '\0', 3*2*nmo_*nmo_ * sizeof(double));
    memset((void*)AY_wz, '\0', 3*2*nmo_*nmo_ * sizeof(double));

    for (int p = 0; p < 3; p++) {

        build_Au_Bu(N, 1, &amps_wx[0].data()[p*N], &AX_wx[p*2*nmo_*nmo_]);
        build_Au_Bu(N, 1, &amps_wx[1].data()[p*N], &AY_wx[p*2*nmo_*nmo_]);

        build_Au_Bu(N, 1, &amps_wy[0].data()[p*N], &AX_wy[p*2*nmo_*nmo_]);
        build_Au_Bu(N, 1, &amps_wy[1].data()[p*N], &AY_wy[p*2*nmo_*nmo_]);

        build_Au_Bu(N, 1, &amps_wz[0].data()[p*N], &AX_wz[p*2*nmo_*nmo_]);
        build_Au_Bu(N, 1, &amps_wz[1].data()[p*N], &AY_wz[p*2*nmo_*nmo_]);
    }

    std::vector<std::string> dir {"X", "Y", "Z"};
    for (int p = 0; p < 3; p++) {
        for (int q = 0; q < 3; q++) {
            for (int r = 0; r < 3; r++) {
                double beta = 0.0;

                // alpha-spin
                for (int i = 0; i < oa; i++) {
                    for (int j = 0; j < oa; j++) {
                        double XY_pq = 0.0;
                        double YX_pq = 0.0;
                        double XY_pr = 0.0;
                        double YX_pr = 0.0;
                        double XY_qr = 0.0;
                        double YX_qr = 0.0;
                        for (int a = 0; a < va; a++) {
                            int ia = i * va + a;
                            int ja = j * va + a;

                            XY_pq += amps_wx[0][p*N + ia] * amps_wy[1][q*N + ja];
                            YX_pq += amps_wx[1][p*N + ia] * amps_wy[0][q*N + ja];

                            XY_pr += amps_wx[0][p*N + ia] * amps_wz[1][r*N + ja];
                            YX_pr += amps_wx[1][p*N + ia] * amps_wz[0][r*N + ja];

                            XY_qr += amps_wy[0][q*N + ia] * amps_wz[1][r*N + ja];
                            YX_qr += amps_wy[1][q*N + ia] * amps_wz[0][r*N + ja];
                        }
                        int ij = i * nmo_ + j;
                        int ji = j * nmo_ + i;

                        beta -= AX_wx[p*2*nmo_*nmo_ + ji] * XY_qr;
                        beta -= AX_wx[p*2*nmo_*nmo_ + ij] * YX_qr;
                        beta -= AY_wx[p*2*nmo_*nmo_ + ij] * XY_qr;
                        beta -= AY_wx[p*2*nmo_*nmo_ + ji] * YX_qr;

                        beta -= AX_wy[q*2*nmo_*nmo_ + ji] * XY_pr;
                        beta -= AX_wy[q*2*nmo_*nmo_ + ij] * YX_pr;
                        beta -= AY_wy[q*2*nmo_*nmo_ + ij] * XY_pr;
                        beta -= AY_wy[q*2*nmo_*nmo_ + ji] * YX_pr;

                        beta -= AX_wz[r*2*nmo_*nmo_ + ji] * XY_pq;
                        beta -= AX_wz[r*2*nmo_*nmo_ + ij] * YX_pq;
                        beta -= AY_wz[r*2*nmo_*nmo_ + ij] * XY_pq;
                        beta -= AY_wz[r*2*nmo_*nmo_ + ji] * YX_pq;
                    }
                }

                for (int a = 0; a < va; a++) {
                    for (int b = 0; b < va; b++) {
                        double XY_pq = 0.0;
                        double YX_pq = 0.0;
                        double XY_pr = 0.0;
                        double YX_pr = 0.0;
                        double XY_qr = 0.0;
                        double YX_qr = 0.0;
                        for (int i = 0; i < oa; i++) {
                            int ia = i * va + a;
                            int ib = i * va + b;

                            XY_pq += amps_wx[0][p*N + ia] * amps_wy[1][q*N + ib];
                            YX_pq += amps_wx[1][p*N + ia] * amps_wy[0][q*N + ib];

                            XY_pr += amps_wx[0][p*N + ia] * amps_wz[1][r*N + ib];
                            YX_pr += amps_wx[1][p*N + ia] * amps_wz[0][r*N + ib];

                            XY_qr += amps_wy[0][q*N + ia] * amps_wz[1][r*N + ib];
                            YX_qr += amps_wy[1][q*N + ia] * amps_wz[0][r*N + ib];
                        }
                        int ab = (a+oa) * nmo_ + (b+oa);
                        int ba = (b+oa) * nmo_ + (a+oa);

                        beta += AX_wx[p*2*nmo_*nmo_ + ab] * XY_qr;
                        beta += AX_wx[p*2*nmo_*nmo_ + ba] * YX_qr;
                        beta += AY_wx[p*2*nmo_*nmo_ + ba] * XY_qr;
                        beta += AY_wx[p*2*nmo_*nmo_ + ab] * YX_qr;

                        beta += AX_wy[q*2*nmo_*nmo_ + ab] * XY_pr;
                        beta += AX_wy[q*2*nmo_*nmo_ + ba] * YX_pr;
                        beta += AY_wy[q*2*nmo_*nmo_ + ba] * XY_pr;
                        beta += AY_wy[q*2*nmo_*nmo_ + ab] * YX_pr;

                        beta += AX_wz[r*2*nmo_*nmo_ + ab] * XY_pq;
                        beta += AX_wz[r*2*nmo_*nmo_ + ba] * YX_pq;
                        beta += AY_wz[r*2*nmo_*nmo_ + ba] * XY_pq;
                        beta += AY_wz[r*2*nmo_*nmo_ + ab] * YX_pq;
                    }
                }

                // derivative wrt field plus second derivative wrt wfn parameters
                //beta += -2.00 * einsum('ji,ai,aj', f[o, o], t1, t1, optimize=['einsum_path', (0, 1), (0, 1)])
                //beta +=  2.00 * einsum('ab,ai,bi', f[v, v], t1, t1, optimize=['einsum_path', (0, 1), (0, 1)])

                // also, third derivative wrt wfn parameters has similar structure
                //beta += -2.00 * einsum('ji,ai,aj,', dipole[o, o], t1, t1, t0_1p, optimize=['einsum_path', (0, 1), (0, 2), (0, 1)])
                //beta +=  2.00 * einsum('ab,ai,bi,', dipole[v, v], t1, t1, t0_1p, optimize=['einsum_path', (0, 1), (0, 2), (0, 1)])

                // alpha-spin
                for (int i = 0; i < oa; i++) {
                    for (int j = 0; j < oa; j++) {
                        double XY_pq = 0.0;
                        double YX_pq = 0.0;
                        double XY_pr = 0.0;
                        double YX_pr = 0.0;
                        double XY_qr = 0.0;
                        double YX_qr = 0.0;
                        for (int a = 0; a < va; a++) {
                            int ia = i * va + a;
                            int ja = j * va + a;

                            XY_pq += 0.5 * amps_wx[0][p*N + ja] * amps_wy[1][q*N + ia];
                            YX_pq += 0.5 * amps_wx[1][p*N + ja] * amps_wy[0][q*N + ia];

                            XY_pr += 0.5 * amps_wx[0][p*N + ja] * amps_wz[1][r*N + ia];
                            YX_pr += 0.5 * amps_wx[1][p*N + ja] * amps_wz[0][r*N + ia];

                            XY_qr += 0.5 * amps_wy[0][q*N + ja] * amps_wz[1][r*N + ia];
                            YX_qr += 0.5 * amps_wy[1][q*N + ja] * amps_wz[0][r*N + ia];
                        }
                        beta += 2 * (XY_qr + YX_qr) * mua_[p]->pointer()[j][i];
                        beta += 2 * (XY_pr + YX_pr) * mua_[q]->pointer()[j][i];
                        beta += 2 * (XY_pq + YX_pq) * mua_[r]->pointer()[j][i];
                    }
                }

                for (int a = 0; a < va; a++) {
                    for (int b = 0; b < va; b++) {
                        double XY_pq = 0.0;
                        double YX_pq = 0.0;
                        double XY_pr = 0.0;
                        double YX_pr = 0.0;
                        double XY_qr = 0.0;
                        double YX_qr = 0.0;
                        for (int i = 0; i < oa; i++) {
                            int ia = i * va + a;
                            int ib = i * va + b;

                            XY_pq += 0.5 * amps_wx[0][p*N + ia] * amps_wy[1][q*N + ib];
                            YX_pq += 0.5 * amps_wx[1][p*N + ia] * amps_wy[0][q*N + ib];

                            XY_pr += 0.5 * amps_wx[0][p*N + ia] * amps_wz[1][r*N + ib];
                            YX_pr += 0.5 * amps_wx[1][p*N + ia] * amps_wz[0][r*N + ib];

                            XY_qr += 0.5 * amps_wy[0][q*N + ia] * amps_wz[1][r*N + ib];
                            YX_qr += 0.5 * amps_wy[1][q*N + ia] * amps_wz[0][r*N + ib];
                        }
                        beta -= 2 * (XY_qr + YX_qr) * mua_[p]->pointer()[a+oa][b+oa];
                        beta -= 2 * (XY_pr + YX_pr) * mua_[q]->pointer()[a+oa][b+oa];
                        beta -= 2 * (XY_pq + YX_pq) * mua_[r]->pointer()[a+oa][b+oa];
                    }
                }

                beta *= -2;
                outfile->Printf("    BETA(%s%s%s) %20.12lf\n", dir[p].c_str(), dir[q].c_str(), dir[r].c_str(), beta);

                // add hyperpolarizabilities to psi variables
                std::string label = "BETA(";
                label += dir[p] + dir[q] + dir[r] + ")";
                std::transform(label.begin(), label.end(), label.begin(),
                    [](unsigned char c) { return std::toupper(c); });
                Process::environment.globals[label.c_str()] = beta;
            }
        }
    }

    free(AX_wx);
    free(AY_wx);
    free(AX_wy);
    free(AY_wy);
    free(AX_wz);
    free(AY_wz);

    return;
}

void RKSResponseSolver::build_Au_Bu(int N, int L, double *u, double *ABu){

    int oa = nalpha_;
    int va = nso_ - oa;

    double ** cap = Ca_->pointer();

    std::vector<SharedMatrix>& C_left  = jk_->C_left();
    std::vector<SharedMatrix>& C_right = jk_->C_right();
    C_left.clear();
    C_right.clear();

    std::vector<std::shared_ptr<Matrix> > Vx;
    std::vector<std::shared_ptr<Matrix> > Dx;

    // J/K-like contributions
    for (int I = 0; I < L; I++) {

        // point to current vector
        double * c = &u[I*N];

        // push density matrices on JK object
        // we need density matrices for alpha

        // Da(mu,nu) = ca(j,b) Ca(mu,j) Ca(nu,b) = Ca(mu,j) Ca'(nu,j)
        // Ca'(nu,j) = ca(j,b) Ca(nu,b)

        // etc.

        // singles

        // alpha
        std::shared_ptr<Matrix> cra (new Matrix(Ca_) );
        std::shared_ptr<Matrix> cla (new Matrix(Ca_) );

        double ** clap = cla->pointer();
        double ** crap = cra->pointer();

        cra->zero();
        cla->zero();

        for (int mu = 0; mu < nso_; mu++) {
            for (int i = 0; i < oa; i++) {

                // left is plain orbitals
                clap[mu][i] = cap[mu][i];

                // right is modified orbitals
                double dum = 0.0;
                for (int a = 0; a < va; a++) {
                    int ia = i * va + a;
                    dum += cap[mu][a+oa] * c[ia];
                }
                crap[mu][i] = dum;

            }
        }

        // push alpha orbitals onto JK object
        C_left.push_back(cla);
        C_right.push_back(cra);

        // build pseudo densities for xc contribution
        if (needs_xc_) {

            auto Dx_a = linalg::doublet(cla, cra, false, true);

            Vx.push_back(std::make_shared<Matrix>("Vax temp", Dx_a->rowspi(), Dx_a->colspi(), Dx_a->symmetry()));

            Dx.push_back(Dx_a);
        }
    }

    // form J/K
    jk_->compute();

    // form xc contributions
    if ( needs_xc_ ) {
        potential_->compute_Vx(Dx, Vx);
    }

    // now, accumulate sigma vectors

    // offset for j/k matrices
    int count = 0;

    double * ea = epsilon_a_->pointer();

    for (int I = 0; I < L; I++) {

        // point to current vector
        double * c = &u[I*N];

        // a <- a+b coulomb
        std::shared_ptr<Matrix> sa (new Matrix(jk_->J()[count]));
        sa->scale(2.0);

        // xc?
        // a <- a
        if ( needs_xc_ ) {
            sa->axpy(2.0, Vx[count]);
        }

        // exact exchange?
        // a <- a
        if (is_x_hybrid_) {
            sa->axpy(-x_alpha_,jk_->K()[count]);
        }

        // LRC functional?
        // a <- a
        if (is_x_lrc_) {
            double beta = 1.0 - x_alpha_;
            sa->axpy(-beta,jk_->wK()[count]);
        }

        // update counter
        count += 1;

        // transform jk, e.g., j(a,i) = j(mu,nu) c(mu,a) c(nu,i)

        sa->transform(Ca_);

        double ** sap = sa->pointer();

        // alpha
        C_DCOPY(nmo_*nmo_, &sap[0][0], 1, &ABu[I*2*nmo_*nmo_], 1);
    }
}

}} // End namespaces

