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
#include <lib3index/3index.h>
#include <libpsio/psio.hpp>
#include <libpsio/psio.h>
#include <libpsio/aiohandler.h>
#include <libqt/qt.h>
#include <psi4-dec.h>
#include <psifiles.h>
#include <libmints/sieve.h>
#include <libiwl/iwl.hpp>
#include "jk.h"
#include "jk_independent.h"
#include "link.h"
#include "direct_screening.h"
#include "cubature.h"
#include "points.h"

#include<lib3index/cholesky.h>

#include <sstream>
#include "libparallel/ParallelPrinter.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace psi;

namespace psi {


CDJK::CDJK(boost::shared_ptr<BasisSet> primary, double cholesky_tolerance):
    DFJK(primary,primary), cholesky_tolerance_(cholesky_tolerance)
{
}
CDJK::~CDJK()
{
}
void CDJK::initialize_JK_disk()
{
    //throw PsiException("Disk algorithm for CD JK not implemented.",__FILE__,__LINE__);
}

void CDJK::initialize_JK_core()
{
    // generate CD integrals with lib3index
    timer_on("CD: cholesky decomposition");
    boost::shared_ptr<IntegralFactory> integral (new IntegralFactory(primary_,primary_,primary_,primary_));
    boost::shared_ptr<CholeskyERI> Ch (new CholeskyERI(boost::shared_ptr<TwoBodyAOInt>(integral->eri()),0.0,cholesky_tolerance_,memory_));
    if(df_ints_io_ == "LOAD")
    {
        Ch->read_previous_cholesky_vector();
    }
    Ch->choleskify();
    ncholesky_  = Ch->Q();

    boost::shared_ptr<Matrix> L = Ch->L();
    double ** Lp = L->pointer();
    timer_off("CD: cholesky decomposition");

    int ntri = sieve_->function_pairs().size();
    ULI three_memory = ncholesky_ * ntri;
    ULI nbf = primary_->nbf();

    if ( memory_  < ((ULI)sizeof(double) * three_memory + (ULI)sizeof(double)* ncholesky_ * nbf * nbf))
        throw PsiException("Not enough memory for CD.",__FILE__,__LINE__);

    Qmn_ = SharedMatrix(new Matrix("Qmn (CD Integrals)", ncholesky_ , ntri));

    double** Qmnp = Qmn_->pointer();

    const std::vector<long int>& schwarz_fun_pairs = sieve_->function_pairs_reverse();

    timer_on("CD: schwarz");
    for (size_t mu = 0; mu < nbf; mu++) {
        for (size_t nu = mu; nu < nbf; nu++) {
            if ( schwarz_fun_pairs[nu*(nu+1)/2+mu] < 0 ) continue;
            for (long int P = 0; P < ncholesky_; P++) {
                Qmnp[P][schwarz_fun_pairs[nu*(nu+1)/2+mu]] = Lp[P][mu * nbf + nu];
            }
        }
    }
    timer_off("CD: schwarz");

    if (df_ints_io_ == "SAVE") {
        // stick ncholesky in process environment for other codes that may use the integrals
        psio_->open(unit_,PSIO_OPEN_NEW);
        psio_->write_entry(unit_, "length", (char*)&ncholesky_, sizeof(long int));
        psio_->write_entry(unit_, "(Q|mn) Integrals", (char*) Lp[0], sizeof(double) * nbf * nbf * ncholesky_);
        psio_->close(unit_,1);
    }
}
void CDJK::manage_JK_core()
{
    for (int Q = 0 ; Q < ncholesky_; Q += max_rows_) {
        int naux = (ncholesky_ - Q <= max_rows_ ? ncholesky_ - Q : max_rows_);
        if (do_J_) {
            timer_on("JK: J");
            block_J(&Qmn_->pointer()[Q],naux);
            timer_off("JK: J");
        }
        if (do_K_) {
            timer_on("JK: K");
            block_K(&Qmn_->pointer()[Q],naux);
            timer_off("JK: K");
        }
    }
}
void CDJK::print_header() const
{
    if (print_) {
        outfile->Printf( "  ==> CDJK: Cholesky-decomposed J/K Matrices <==\n\n");

        outfile->Printf( "    J tasked:             %11s\n", (do_J_ ? "Yes" : "No"));
        outfile->Printf( "    K tasked:             %11s\n", (do_K_ ? "Yes" : "No"));
        outfile->Printf( "    wK tasked:            %11s\n", (do_wK_ ? "Yes" : "No"));
        if (do_wK_) {
            throw PsiException("no wk for scf_type cd.",__FILE__,__LINE__);
            //outfile->Printf( "    Omega:                %11.3E\n", omega_);
        }
        outfile->Printf( "    OpenMP threads:       %11d\n", omp_nthread_);
        outfile->Printf( "    Integrals threads:    %11d\n", df_ints_num_threads_);
        outfile->Printf( "    Memory (MB):          %11ld\n", (memory_ *8L) / (1024L * 1024L));
        outfile->Printf( "    Algorithm:            %11s\n",  (is_core_ ? "Core" : "Disk"));
        outfile->Printf( "    Integral Cache:       %11s\n",  df_ints_io_.c_str());
        outfile->Printf( "    Schwarz Cutoff:       %11.0E\n", cutoff_);
        outfile->Printf( "    Cholesky tolerance:   %11.2E\n", cholesky_tolerance_);
        outfile->Printf( "    No. Cholesky vectors: %11li\n\n", ncholesky_);
    }
}
}
