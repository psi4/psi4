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


#include "soscf.h"
#include "jk.h"

namespace psi {


/// SORHF class

SORHF::SORHF(boost::shared_ptr<JK> jk)
{
    jk_ = jk;
}

SORHF::~SORHF()
{
}

void SORHF::update(SharedMatrix Cocc, SharedMatrix Cvir, SharedMatrix Fock)
{
    // Build C matrices
    nocc_ = Cocc->ncol();
    nvir_ = Cvir->ncol();

    matrices_["Cocc"] = Cocc;
    matrices_["Cvir"] = Cvir;
    std::vector<boost::shared_ptr<Matrix> > fullC;
    fullC.push_back(Cocc);
    fullC.push_back(Cvir);
    matrices_["C"] = Matrix::horzcat(fullC);

    nao_ = matrices_["C"]->ncol();

    // => MO Fock Matrix <= //
    matrices_["IFock"] = Matrix::triplet(matrices_["C"], Fock, matrices_["C"], true, false, false);

    // => Gradient and Diagonal denominator <= //
    matrices_["Gradient"] = SharedMatrix(new Matrix("Gradient", nocc_, nvir_));
    matrices_["ia_denom"] = SharedMatrix(new Matrix("ia_denom", nocc_, nvir_));

    double* gp = matrices_["Gradient"]->pointer()[0];
    double* denomp = matrices_["ia_denom"]->pointer()[0];
    double** fp = matrices_["IFock"]->pointer();
    for (size_t i=0, target=0; i<nocc_; i++){
        double o_eps = fp[i][i];
        for (size_t j=nocc_; j < nao_; j++)
        {
            gp[target] = -4.0 * fp[i][j];
            denomp[target++] = 4.0 * (fp[j][j] - o_eps);
        }
    }
}

SharedMatrix SORHF::Hx(SharedMatrix x)
{

    // => Effective one electron part <= //

    double** IFp = matrices_["IFock"]->pointer();

    SharedMatrix Focc(new Matrix("Focc", nocc_, nocc_));
    double* Foccp = Focc->pointer()[0];

    SharedMatrix Fvir(new Matrix("Fvir", nvir_, nvir_));
    double* Fvirp = Fvir->pointer()[0];

    for (size_t i=0, target=0; i<nocc_; i++){
        for (size_t j=0; j < nocc_; j++){
            Foccp[target++] = IFp[i][j];
        }
    }

    for (size_t i=nocc_, target=0; i<nao_; i++){
        for (size_t j=nocc_; j < nao_; j++){
            Fvirp[target++] = IFp[i][j];
        }
    }

    SharedMatrix ret = Matrix::doublet(Focc, x);
    ret->gemm(false, false, -1.0, x, Fvir, 1.0);


    // => Two electron part <= //

    std::vector<SharedMatrix>& Cl = jk_->C_left();
    std::vector<SharedMatrix>& Cr = jk_->C_right();
    Cl.clear();
    Cr.clear();

    Cl.push_back(matrices_["Cocc"]);

    SharedMatrix R = Matrix::doublet(matrices_["Cvir"], x, false, true);
    R->scale(-1.0);
    Cr.push_back(R);

    jk_->compute();

    const std::vector<SharedMatrix>& J = jk_->J();
    const std::vector<SharedMatrix>& K = jk_->K();

    J[0]->scale(4.0);
    J[0]->subtract(K[0]);
    J[0]->subtract(K[0]->transpose());
    SharedMatrix M = Matrix::triplet(matrices_["Cocc"], J[0], matrices_["Cvir"], true, false, false);

    ret->add(M);
    ret->scale(-4.0);
    return ret;
}

SharedMatrix SORHF::solve(int max_iter, double conv, bool print)
{
    if (print){
        outfile->Printf("    ==> SORHF Iterations <==\n");
        outfile->Printf("    Maxiter     = %11d\n", max_iter);
        outfile->Printf("    Convergence = %11.3E\n", conv);
        outfile->Printf("    ---------------------------------------\n");
        outfile->Printf("    %-4s   %11s     %10s\n", "Iter", "Residual RMS", "Time [s]");
        outfile->Printf("    ---------------------------------------\n");
    }

    time_t start;
    time_t stop;
    start = time(NULL);

    // Initial guess
    SharedMatrix x = matrices_["Gradient"]->clone();
    x->apply_denominator(matrices_["ia_denom"]);

    SharedMatrix r = matrices_["Gradient"]->clone();
    r->subtract(Hx(x));

    SharedMatrix z = r->clone();
    z->apply_denominator(matrices_["ia_denom"]);

    SharedMatrix p = z->clone();

    size_t npairs = x->nrow() * x->ncol();
    for (size_t iter = 0; iter < max_iter; iter++) {
        SharedMatrix Ap = Hx(p);

        double rzpre = r->vector_dot(z);
        double alpha = rzpre / p->vector_dot(Ap);

        double* xp = x->pointer()[0];
        double* rp = r->pointer()[0];
        double* pp = p->pointer()[0];
        double* App = Ap->pointer()[0];
        C_DAXPY(npairs,  alpha, pp, 1, xp, 1);
        C_DAXPY(npairs, -alpha, App, 1, rp, 1);


        double rconv = r->rms();
        stop = time(NULL);
        if (print){
            outfile->Printf("    %-4d %11.3E %10ld\n", iter+1, rconv, stop-start);
        }

        if (rconv < conv){
            break;
        }

        z->copy(r);
        z->apply_denominator(matrices_["ia_denom"]);

        double beta = r->vector_dot(z) / rzpre;

        p->scale(beta);
        p->add(z);

    }
    if (print){
        outfile->Printf("    ---------------------------------------\n");
        outfile->Printf("\n");
    }
    return x;
}

/// SOMCSCF class

} // Namespace psi
