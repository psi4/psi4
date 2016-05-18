/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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
#include <libpsio/psio.hpp>
#include <libpsio/psio.h>
#include <psi4-dec.h>
#include <psifiles.h>
#include <libqt/qt.h>
#include "jk.h"
#include "PKmanagers.h"

#include <sstream>
#include "libparallel/ParallelPrinter.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace psi;

namespace psi {


PKJK::PKJK(boost::shared_ptr<BasisSet> primary, Options& options) :
    JK(primary), options_(options)
{
    common_init();
}

PKJK::~PKJK()
{
}


void PKJK::common_init()
{
    pk_file_ = PSIF_SO_PK;
    nthreads_ = 1;
#ifdef _OPENMP
    nthreads_ = omp_get_max_threads();
#endif
}

bool PKJK::C1() const {
    return true;
}

void PKJK::print_header() const
{
    if (print_) {
        outfile->Printf( "  ==> DiskJK: Disk-Based J/K Matrices <==\n\n");

        outfile->Printf( "    J tasked:          %11s\n", (do_J_ ? "Yes" : "No"));
        outfile->Printf( "    K tasked:          %11s\n", (do_K_ ? "Yes" : "No"));
        outfile->Printf( "    wK tasked:         %11s\n", (do_wK_ ? "Yes" : "No"));
        if (do_wK_)
            outfile->Printf( "    Omega:             %11.3E\n", omega_);
        outfile->Printf( "    Memory (MB):       %11ld\n", (memory_ *8L) / (1024L * 1024L));
        outfile->Printf( "    Schwarz Cutoff:    %11.0E\n\n", cutoff_);
        outfile->Printf( "    OpenMP threads:    %11d\n", nthreads_);
    }
}

void PKJK::preiterations()
{

    //Build PKManager to get proper algorithm set up
    Options& options = Process::environment.options;

    psio_ = _default_psio_lib_;

    timer_on("Total PK formation time");
    // We compute the integrals so that we can directly write the
    // PK file to disk. Also, do everything in the AO basis
    // like the modern JK algos, for adding sieving later

    PKmanager_ = pk::PKManager::build_PKManager(psio_,primary_,memory_,options,do_wK_,omega_);

    outfile->Printf(" Computing reordered integrals for PK\n\n");

    PKmanager_->initialize();

    PKmanager_->form_PK();

    // PK files are written at this point. We are done.
    timer_off("Total PK formation time");

    // If range-separated K needed, we redo all the above steps
    if(do_wK_) {
        outfile->Printf("  Computing range-separated integrals for PK\n");

        PKmanager_->initialize_wK();

        PKmanager_->form_PK_wK();

    }

    /*
     * For omega, we only need exchange and it's done separately from conventional terms, so we can
     * use fewer batches in principle.  For now we just use the batching scheme that's already been
     * computed for the conventional J/K combo, above
     */
/*    if(do_wK_){
        for(int batch = 0; batch < nbatches; ++batch){
            size_t min_index   = batch_index_min_[batch];
            size_t max_index   = batch_index_max_[batch];
            size_t batch_size = max_index - min_index;
            double *wk_block = new double[batch_size];
            ::memset(wk_block, '\0', batch_size * sizeof(double));

            IWL *iwl = new IWL(psio_.get(), PSIF_SO_ERF_TEI, cutoff_, 1, 1);
            Label *lblptr = iwl->labels();
            Value *valptr = iwl->values();
            int labelIndex, pabs, qabs, rabs, sabs, prel, qrel, rrel, srel, psym, qsym, rsym, ssym;
            size_t bra, ket, braket;
            double value;
            bool last_buffer;
            do{
                last_buffer = iwl->last_buffer();
                for(int index = 0; index < iwl->buffer_count(); ++index){
                    labelIndex = 4*index;
                    pabs  = abs((int) lblptr[labelIndex++]);
                    qabs  = (int) lblptr[labelIndex++];
                    rabs  = (int) lblptr[labelIndex++];
                    sabs  = (int) lblptr[labelIndex++];
                    prel  = so2index_[pabs];
                    qrel  = so2index_[qabs];
                    rrel  = so2index_[rabs];
                    srel  = so2index_[sabs];
                    psym  = so2symblk_[pabs];
                    qsym  = so2symblk_[qabs];
                    rsym  = so2symblk_[rabs];
                    ssym  = so2symblk_[sabs];
                    value = (double) valptr[index];

                    if ((psym == qsym) && (rsym == ssym)) {
                        // K (2nd sort)
                        if (build_k && (prel != qrel) && (rrel != srel)) {
                            if ((psym == ssym) && (qsym == rsym)) {
                                bra = INDEX2(prel, srel);
                                ket = INDEX2(qrel, rrel);
                                braket = INDEX2(bra + pk_symoffset[psym], ket + pk_symoffset[qsym]);
                                if((braket >= min_index) && (braket < max_index)){
                                    if ((prel == srel) || (qrel == rrel)) {
                                        wk_block[braket - min_index] += value;
                                    } else {
                                        wk_block[braket - min_index] += 0.5 * value;
                                    }
                                }
                            }
                        }
                    }

                    // K (1st sort)
                    if (build_k && (psym == rsym) && (qsym == ssym)) {
                        bra = INDEX2(prel, rrel);
                        ket = INDEX2(qrel, srel);
                        braket = INDEX2(bra + pk_symoffset[psym], ket + pk_symoffset[qsym]);
                        if((braket >= min_index) && (braket < max_index)){
                            if ((prel == rrel) || (qrel == srel)) {
                                wk_block[braket - min_index] += value;
                            } else {
                                wk_block[braket - min_index] += 0.5 * value;
                            }
                        }
                    }
                }
                if (!last_buffer) iwl->fetch();
            } while (!last_buffer);
            delete iwl;

            // Halve the diagonal elements held in core
            for(size_t pq = batch_pq_min_[batch]; pq < batch_pq_max_[batch]; ++pq){
                size_t address = INDEX2(pq, pq) - min_index;
                wk_block[address] *= 0.5;
            }

            char *label = new char[100];
            sprintf(label, "wK Block (Batch %d)", batch);
            psio_->write_entry(pk_file_, label, (char*) wk_block, batch_size * sizeof(double));
            delete [] label;

            delete [] wk_block;
        } // End of loop over batches
    }

    delete [] orb_offset;

//    if(!file_was_open);
        psio_->close(pk_file_, 1);
*/
}

void PKJK::compute_JK()
{
    // We form the vector containing the density matrix triangular elements
    PKmanager_->prepare_JK(D_ao_,C_left_ao_,C_right_ao_);

    if(J_ao_.size()) {
        // We can safely pass K here since its size is checked within
        // the routine
        PKmanager_->form_J(J_ao_,"",K_ao_);
    }
    if(K_ao_.size()) {
        PKmanager_->form_K(K_ao_);
    }
    if(wK_ao_.size()) {
        PKmanager_->form_wK(wK_ao_);
    }

    PKmanager_->finalize_JK();

}

void PKJK::postiterations()
{
}

}
