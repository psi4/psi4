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
    // => Build C matrices <= //
    nocc_ = Cocc->ncol();
    nvir_ = Cvir->ncol();

    nirrep_ = Cocc->nirrep();
    noccpi_ = Cocc->colspi();
    nvirpi_ = Cvir->colspi();

    matrices_["Cocc"] = Cocc;
    matrices_["Cvir"] = Cvir;
    std::vector<boost::shared_ptr<Matrix> > fullC;
    fullC.push_back(Cocc);
    fullC.push_back(Cvir);
    matrices_["C"] = Matrix::horzcat(fullC);

    nso_ = matrices_["C"]->ncol();
    nsopi_ = matrices_["C"]->colspi();

    // => MO Fock Matrix <= //
    matrices_["IFock"] = Matrix::triplet(matrices_["C"], Fock, matrices_["C"], true, false, false);

    // => Gradient and Diagonal denominator <= //
    matrices_["Gradient"] = SharedMatrix(new Matrix("Gradient", nirrep_, noccpi_, nvirpi_));
    matrices_["ia_denom"] = SharedMatrix(new Matrix("ia_denom", nirrep_, noccpi_, nvirpi_));

    for (size_t h=0; h<nirrep_; h++){

        if (!nsopi_[h] || !noccpi_[h] || !nvirpi_[h]) continue;
        double* gp = matrices_["Gradient"]->pointer(h)[0];
        double* denomp = matrices_["ia_denom"]->pointer(h)[0];
        double** fp = matrices_["IFock"]->pointer(h);

        for (size_t i=0, target=0; i<noccpi_[h]; i++){
            double o_eps = fp[i][i];
            for (size_t j=noccpi_[h]; j < nsopi_[h]; j++)
            {
                gp[target] = -4.0 * fp[i][j];
                denomp[target++] = 4.0 * (fp[j][j] - o_eps);
            }
        }
    }
}

SharedMatrix SORHF::Ck(SharedMatrix x)
{

    SharedMatrix tmp(new Matrix("Ck", nirrep_, nsopi_, nsopi_));

    // Form full antisymmetric matrix
    for (size_t h=0; h<nirrep_; h++){

        if (!nsopi_[h] || !noccpi_[h] || !nvirpi_[h]) continue;
        double** tp = tmp->pointer(h);
        double*  xp = x->pointer(h)[0];

        // Matrix::schmidt orthogonalizes rows not columns so we need to transpose
        for (size_t i=0, target=0; i<noccpi_[h]; i++){
            for (size_t a=noccpi_[h]; a < nsopi_[h]; a++){
                tp[i][a] = -1.0 * xp[target];
                tp[a][i] = xp[target++];
            }
        }
    }

    // Build exp(U) = 1 + U + 0.5 U U
    SharedMatrix U = tmp->clone();
    for (size_t h=0; h<nirrep_; h++){
        double** up = U->pointer(h);
        for (size_t i=0; i<U->rowspi()[h]; i++){
            up[i][i] += 1.0;
        }
    }
    U->gemm(false, false, 0.5, tmp, tmp, 1.0);

    // We did not fully exponentiate the matrix, need to orthogonalize
    U->schmidt();

    // C' = C U
    tmp->gemm(false, false, 1.0, matrices_["C"], U, 0.0);

    return tmp;
}

SharedMatrix SORHF::Hk(SharedMatrix x)
{

    // => Effective one electron part <= //
    SharedMatrix Focc(new Matrix("Focc", nirrep_, noccpi_, noccpi_));
    SharedMatrix Fvir(new Matrix("Fvir", nirrep_, nvirpi_, nvirpi_));

    for (size_t h=0; h<nirrep_; h++){
        if (!nsopi_[h] || !noccpi_[h] || !nvirpi_[h]) continue;
        double** IFp = matrices_["IFock"]->pointer(h);
        double* Foccp = Focc->pointer(h)[0];
        double* Fvirp = Fvir->pointer(h)[0];

        for (size_t i=0, target=0; i<noccpi_[h]; i++){
            for (size_t j=0; j < noccpi_[h]; j++){
                Foccp[target++] = IFp[i][j];
            }
        }

        for (size_t i=noccpi_[h], target=0; i<nsopi_[h]; i++){
            for (size_t j=noccpi_[h]; j < nsopi_[h]; j++){
                Fvirp[target++] = IFp[i][j];
            }
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
    Cl.clear();
    Cr.clear();

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
        outfile->Printf("\n");
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

    // Calc hessian vector product, find residual and conditioned residual
    SharedMatrix r = matrices_["Gradient"]->clone();
    SharedMatrix Ap = Hk(x);
    r->subtract(Ap);

    SharedMatrix z = r->clone();
    z->apply_denominator(matrices_["ia_denom"]);

    SharedMatrix p = z->clone();

    for (size_t iter = 0; iter < max_iter; iter++) {

        // Calc hessian vector product
        Ap = Hk(p);

        // Find factors and scale
        double rzpre = r->vector_dot(z);
        double alpha = rzpre / p->vector_dot(Ap);

        for (size_t h=0; h<nirrep_; h++){
            if (!nsopi_[h] || !noccpi_[h] || !nvirpi_[h]) continue;

            size_t npairs = noccpi_[h] * nvirpi_[h];
            double* xp = x->pointer(h)[0];
            double* rp = r->pointer(h)[0];
            double* pp = p->pointer(h)[0];
            double* App = Ap->pointer(h)[0];

            C_DAXPY(npairs,  alpha, pp, 1, xp, 1);
            C_DAXPY(npairs, -alpha, App, 1, rp, 1);
        }

        // Get residual
        double rconv = r->rms();
        stop = time(NULL);
        if (print){
            outfile->Printf("    %-4d %11.3E %10ld\n", iter+1, rconv, stop-start);
        }

        // Check convergence
        if (rconv < conv){
            break;
        }

        // Update p and z
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

/// End SORHF class

/// SOMCSCF class

} // Namespace psi
