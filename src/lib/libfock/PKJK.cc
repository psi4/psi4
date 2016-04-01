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
#include "yoshimine.h"

#include<lib3index/cholesky.h>

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
        //outfile->Printf( "    OpenMP threads:    %11d\n", omp_nthread_);
    }
}

void PKJK::preiterations()
{
    // We need to access the option object to just get a few values

    Options& options = Process::environment.options;
    psio_ = _default_psio_lib_;

    algo_ = options.get_str("PK_ALGO");


    if (algo_ == "REORDER") {
    // We compute the integrals so that we can directly write the
    // PK file to disk. Also, do everything in the AO basis
    // like the modern JK algos, for adding sieving later
      integrals_reorder();
    // PK files are written at this point. We are done.
      return;
    }
    // ==> Section below only executed if algo != REORDER <==
    // #########################################################


    // Start by generating conventional integrals on disk
    integrals();
//    boost::shared_ptr<MintsHelper> mints(new MintsHelper());
//    mints->integrals();
//    if(do_wK_)
//        mints->integrals_erf(omega_);
//    mints.reset();

    nsopi_ = mints->sobasisset()->dimension();
    nso_ = nsopi_.sum();
    nirrep_ = nsopi_.n();

    //TODO: Delete the following ?
    so2symblk_ = new int[nso];
    so2index_  = new int[nso];
    int so_count = 0;
    int offset = 0;
    for (int h = 0; h < nirreps; ++h) {
        for (int i = 0; i < sopi[h]; ++i) {
            so2symblk_[so_count] = h;
            so2index_[so_count] = so_count-offset;
            ++so_count;
        }
        offset += nsopi_[h];
    }

    // The first SO in each irrep
    int* orb_offset = new int[nirrep_];
    orb_offset[0] = 0;
    for(int h = 1; h < nirrep_; ++h)
        orb_offset[h] = orb_offset[h-1] + nsopi_[h-1];

    // Compute PK symmetry mapping
    int *pk_symoffset = new int[nirrep_];
    pk_size_ = 0;  // Size of PK matrix triangle, i.e. ncol*(ncol - 1)/2 + ncol
    pk_pairs_ = 0; // Number of pk pairs, i.e. ncol of PK matrix
    for(int h = 0; h < nirrep_; ++h){
        pk_symoffset[h] = pk_pairs_;
        // Add up possible pair combinations that yield A1 symmetry
        pk_pairs_ += nsopi_[h]*(nsopi_[h] + 1)/2;
    }
    // Compute the number of pairs in PK
    pk_size_ = INDEX2(pk_pairs_-1, pk_pairs_-1) + 1;

    // I still have to figure out why this piece of code exists.
    // Count the pairs
    size_t npairs = 0;
    size_t *pairpi = new size_t[nirrep_];
    ::memset(pairpi, '\0', nirrep_*sizeof(size_t));
    for(int pq_sym = 0; pq_sym < nirrep_; ++pq_sym){
        for(int p_sym = 0; p_sym < nirrep_; ++p_sym){
            int q_sym = pq_sym ^ p_sym;
            if(p_sym >= q_sym){
                for(int p = 0; p < nsopi_[p_sym]; ++p){
                    for(int q = 0; q < nsopi_[q_sym]; ++q){
                        int p_abs = p + orb_offset[p_sym];
                        int q_abs = q + orb_offset[q_sym];
                        if(p_abs >= q_abs){
                            pairpi[pq_sym]++;
                            npairs++;
                        }
                    }
                }
            }
        }
    }
    delete [] orb_offset;

    /* Now we create the Yoshimine buffers for sorting. The PK supermatrix
     containing integrals (ij|kl) has all indices ij on the rows and indices
     kl on the columns. We need both the supermatrix for building the Coulomb
     matrix J and the exchange matrix K. For J, integrals must be sorted so
     that (ij|kl) goes to position (ij, kl) of the supermatrix. However, for K
     we want (ij|kl) to go in both (ik, jl) and (il, jk) in the supermatrix (not
     considering special cases for diagonal elements). As a result, we need
     to sort integrals in different files for J and K. For this purpose, we are going
     to use two Yoshimine buffers. They should each be initialized with half of the
     available memory, because they will be used simultaneously for the pre-sorting.
     */

    // Use a Yoshimine sorting, with objects adapted from TRANSQT implementation
    // of Yoshimine.

    int max_buckets = options.get_int("MAX_BUCKETS");

    // Initialize each buffer with only half of the memory. memory_ is in doubles
    // and already has a safety factor.
    long pk_sort_mem = memory_;
    long pk_presort_mem = memory_ * sizeof(double) / 2;

    // First temporary file for pre-sorting.
    int first_tmp_file_J = options.get_int("FIRST_TMP_FILE");

    // Integral tolerance
    double tolerance = options.get_double("INTS_TOLERANCE");

//DEBUG    outfile->Printf("The pk_sort_mem is %li doubles, there are %i max_buckets\n", pk_sort_mem, max_buckets);
//DEBUG    outfile->Printf("The first tmp file is %i and the integral tolerance %f\n", first_tmp_file_J, tolerance);

    Yosh YBuffJ(pk_pairs_, pk_presort_mem, pk_sort_mem,
                max_buckets, first_tmp_file_J, tolerance, psio_.get());
 //Oldstyle   YoshLight YBuffJ(pk_pairs_, pk_presort_mem, pk_sort_mem,
 //Oldstyle               max_buckets, first_tmp_file_J, tolerance, psio_.get());

    int first_tmp_file_K = first_tmp_file_J + YBuffJ.nbuckets();
//DEBUG    outfile->Printf("The first tmp K file is %i \n", first_tmp_file_K);

    Yosh YBuffK(pk_pairs_, pk_presort_mem, pk_sort_mem,
                max_buckets, first_tmp_file_K, tolerance, psio_.get());
 //Oldstyle   YoshLight YBuffK(pk_pairs_, pk_presort_mem, pk_sort_mem,
 //Oldstyle               max_buckets, first_tmp_file_K, tolerance, psio_.get());

//DEBUG    transqt::yosh_print(&YBuffJ, "outfile");
//DEBUG    transqt::yosh_print(&YBuffK, "outfile");

    // We copy the bucket info into our local arrays for later retrieval.
    // Low and high pq indices and number of buckets should be the same
    // between J and K

    int nbatches      = YBuffJ.nbuckets();
    size_t totally_symmetric_pairs = pairpi[0];
    delete [] pairpi;
    batch_pq_min_.clear();
    batch_pq_max_.clear();
    batch_index_min_.clear();
    batch_index_max_.clear();
    for(int i = 0; i < YBuffJ.nbuckets(); ++i) {
        int lowpq = (YBuffJ.buckets())[i]->lo();
        int hipq = (YBuffJ.buckets())[i]->hi();
        batch_pq_min_.push_back(lowpq);
        batch_pq_max_.push_back(++hipq);
        batch_index_min_.push_back( (size_t)(lowpq + 1L) * (size_t)lowpq / 2L);
        batch_index_max_.push_back((size_t)hipq * (size_t)(hipq + 1L) / 2L);
    }

    for(int batch = 0; batch < nbatches; ++batch){
        outfile->Printf("\tBatch %3d pq = [%8zu,%8zu] index = [%14zu,%zu]\n",
                batch + 1,
                batch_pq_min_[batch],batch_pq_max_[batch],
                batch_index_min_[batch],batch_index_max_[batch]);
    }

    //Now, what we are doing depends on the selected algorithm and the
    //number of buckets needed.

    if(algo_ == "DIRECT") {
        if(nbatches != 1) {
            // lol direct cannot work with 2 batches lol
            throw PSIEXCEPTION("Direct ? With 2 batches ? lol\n");
        } else {
            throw PSIEXCEPTION("Not implemented yet.\n");
        }
    } else if (algo_ == "NEW") {
       // Proceed with Yoshimine sorting.

    // Initiate buckets: allocate array memory and open temporary files.

    timer_on("Total PK formation time");
    YBuffJ.init_buckets();
    YBuffK.init_buckets();

//    transqt::yosh_print(&YBuffJ, "outfile");
//    transqt::yosh_print(&YBuffK, "outfile");
   //Do the pre-sorting step in temporary files

    //TODO: Solve the double sorting of K integrals that are not in the same
    // bucket, we are writing more integrals to disk than anticipated.
    YoshBase::rdtwo_pk(&YBuffJ,&YBuffK,PSIF_SO_TEI, 0, nirreps, so2index_, so2symblk_,
                           pk_symoffset, (debug_ > 5));

    // Close the buckets, but keep the temp. files.
    YBuffJ.close_buckets(0);
    YBuffK.close_buckets(0);

    // Really proceed with the sorting and writing of the PK files

    psio_->open(pk_file_, PSIO_OPEN_NEW);
//    psio_->open(pk_file_, PSIO_OPEN_OLD);

//DEBUG    outfile->Printf("Just before sorting\n");
    YBuffJ.sort_pk(0, pk_file_, 0, so2index_, so2symblk_, pk_symoffset, (debug_ > 5));
//    throw PSIEXCEPTION("Integral counting");
    YBuffK.sort_pk(1, pk_file_, 0, so2index_, so2symblk_, pk_symoffset, (debug_ > 5));

    timer_off("Total PK formation time");
 //   tbench.start();
 //   IWL *iwl = new IWL(psio_.get(), PSIF_SO_TEI, cutoff_, 1, 1);
 //   bool last_buffer;
 //   do{
 //       last_buffer = iwl->last_buffer();
 //       if (!last_buffer) {
 //           iwl->fetch();
 //       }
 //   } while (!last_buffer);
 //   delete iwl;
 //   tbench.stop("Dummy Read original IWL file");
    psio_->close(pk_file_, 1);

    YBuffJ.done();
    YBuffK.done();

    } // End of condition for Yoshimine sorting, algo NEW
    else if(algo_ == "JET") {
        // Here we go for the old algorithm

    // We can dismiss the Yoshimine buckets

    YBuffJ.done();
    YBuffK.done();

    // We might want to only build p in future...
    bool build_k = true;

    timer_on("Total PK formation time");
    psio_->open(pk_file_, PSIO_OPEN_NEW);
    for(int batch = 0; batch < nbatches; ++batch){
        size_t min_index   = batch_index_min_[batch];
        size_t max_index   = batch_index_max_[batch];
        size_t batch_size = max_index - min_index;
        double *j_block = new double[batch_size];
        double *k_block = new double[batch_size];
        ::memset(j_block, '\0', batch_size * sizeof(double));
        ::memset(k_block, '\0', batch_size * sizeof(double));

        std::stringstream readtag;
        readtag << "Read num. ";
        readtag << batch;
        readtag << " of IWL file.";
        timer_on(readtag.str().c_str());
        IWL *iwl = new IWL(psio_.get(), PSIF_SO_TEI, cutoff_, 1, 1);
        timer_off(readtag.str().c_str());
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
//                outfile->Printf("IWL <%d %d|%d %d>\n", pabs, qabs, rabs, sabs);

                // J
                if ((psym == qsym) && (rsym == ssym)) {
                    bra = INDEX2(prel, qrel);
                    ket = INDEX2(rrel, srel);
                    // pk_symoffset_ corrects for the symmetry offset in the pk_ vector
                    braket = INDEX2(bra + pk_symoffset[psym], ket + pk_symoffset[rsym]);
                    if((braket >= min_index) && (braket < max_index)){
                        j_block[braket - min_index] += value;
                    }

                    // K (2nd sort)
                    if (build_k && (prel != qrel) && (rrel != srel)) {
                        if ((psym == ssym) && (qsym == rsym)) {
                            bra = INDEX2(prel, srel);
                            ket = INDEX2(qrel, rrel);
                            braket = INDEX2(bra + pk_symoffset[psym], ket + pk_symoffset[qsym]);
                            if((braket >= min_index) && (braket < max_index)){
                                if ((prel == srel) || (qrel == rrel)) {
                                    k_block[braket - min_index] += value;
                                } else {
                                    k_block[braket - min_index] += 0.5 * value;
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
                            k_block[braket - min_index] += value;
                        } else {
                            k_block[braket - min_index] += 0.5 * value;
                        }
                    }
                }
            }
            if (!last_buffer) {
                timer_on(readtag.str().c_str());
                iwl->fetch();
                timer_off(readtag.str().c_str());
            }
        } while (!last_buffer);
        timer_on(readtag.str().c_str());
        delete iwl;
        timer_off(readtag.str().c_str());

        // Halve the diagonal elements held in core
        for(size_t pq = batch_pq_min_[batch]; pq < batch_pq_max_[batch]; ++pq){
            size_t address = INDEX2(pq, pq) - min_index;
            j_block[address] *= 0.5;
            k_block[address] *= 0.5;
        }

        char *label = new char[100];
        timer_on("Write J batches");
        sprintf(label, "J Block (Batch %d)", batch);
        psio_->write_entry(pk_file_, label, (char*) j_block, batch_size * sizeof(double));
        timer_off("Write J batches");
        timer_on("Write K batches");
        sprintf(label, "K Block (Batch %d)", batch);
        psio_->write_entry(pk_file_, label, (char*) k_block, batch_size * sizeof(double));
        timer_off("Write K batches");
        delete [] label;

        delete [] j_block;
        delete [] k_block;
    } // End of loop over batches
    timer_off("Total PK formation time");

  //  tread.start();
  //  IWL *iwl = new IWL(psio_.get(), PSIF_SO_TEI, cutoff_, 1, 1);
  //  bool last_buffer;
  //  do{
  //      last_buffer = iwl->last_buffer();
  //      if (!last_buffer) {
  //          iwl->fetch();
  //      }
  //  } while (!last_buffer);
  //  delete iwl;
  //  tread.stop("Dummy Read original IWL file");
    psio_->close(pk_file_, 1);

    } // End of algo conditional.
    // #########################################################
    // ==> End of section executed only for algo != REORDER <==

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
    if (algo_ == "REORDER") {
        // We are using AO integrals and two files, one for J and one for K
        // For now we are not handling asymmetric density matrices
        PKmanager_.open_files();

        // We form the vector containing the density matrix triangular elements
        for(int N = 0; N < D_ao_.size(); ++N) {
            if(C_left_[N] != C_right_[N]) {
                outfile->Printf("ERROR: REORDER algorithm does not support this computation\n");
                throw PSIEXCEPTION("Only symmetric density matrices implemented for now\n");
            }
        }
        PKmanager_.form_D_vec(D_ao_);

        // Now we actually compute J and K as needed
        if(J_.size()) {
            PKmanager_.form_J(J_);
        }
        if(K_.size()) {
            PKmanager_.form_K(K_);
        }
        PKmanager_.finalize_D();
        // Everything done, we close the files
        PKmanager_.close_files();

        return;

    }

    int nirreps = Process::environment.wavefunction()->nirrep();
    int *sopi   = Process::environment.wavefunction()->nsopi();

//    bool file_was_open = psio_->open_check(pk_file_);
//    if(!file_was_open);
        psio_->open(pk_file_, PSIO_OPEN_OLD);


    int nbatches = batch_pq_min_.size();
    std::vector<double*> J_vectors;
    std::vector<double*> K_vectors;
    std::vector<double*> D_vectors;

    // We do some more preliminary work in case a density matrix is asymmetric

    int *qind = NULL;
    int *pind = NULL;
    bool asym = false;
    bool asym_only = true;
    size_t square_dim = 0;

    for(size_t N = 0; N < D_.size(); ++N) {
        if(C_left_[N] == C_right_[N]) {
            asym_only = false;
        }
    }

    for(size_t N = 0; N < D_.size(); ++N) {
        if(C_left_[N] != C_right_[N]) {
            asym = true;
            // Limited support for asymmetric density matrices: only C1 symmetry.
            if ( nirrep_ != 1 ) {
                throw PSIEXCEPTION("PK integrals can be used for this type of calculation only in C1 symmetry");
            }
            square_dim = nsopi_[0] * nsopi_[0];
            // So we are going to do something dumb and store
            // explicitly all indices first.
            pind = new int[pk_pairs_];
            qind = new int[pk_pairs_];
            int counter = 0;
            for(int p = 0; p < nsopi_[0]; ++p) {
                for(int q = 0; q <= p; ++q) {
                    pind[counter] = p;
                    qind[counter] = q;
                    ++counter;
                }
            }
            break;
        }
    }

    /*
     * The J terms
     */
    for(size_t N = 0; N < J_.size(); ++N){
        if(D_[N]->symmetry())
            throw PSIEXCEPTION("PK integrals cannot be used for this type of calculation.");
        if(C_left_[N] != C_right_[N]) {
            // Placeholder value to preserve proper ordering.
            J_vectors.push_back(NULL);
            // We need to halve diagonal elements of the density matrix
            double *D_vector = new double[square_dim];
            ::memset(D_vector, 0, square_dim * sizeof(double));
            D_vectors.push_back(D_vector);
            for(int p = 0; p < nsopi_[0]; ++p) {
                for(int q = 0; q < nsopi_[0]; ++q) {
                    if(p != q) {
                        D_vector[p*nsopi_[0] + q] = D_[N]->get(p,q);
                    } else {
                        D_vector[p*nsopi_[0] + q] = 0.5 * D_[N]->get(p,q);
                    }
                }
            }
        } else {
            double *J_vector = new double[pk_pairs_];
            ::memset(J_vector,  0, pk_pairs_ * sizeof(double));
            J_vectors.push_back(J_vector);
            double *D_vector = new double[pk_pairs_];
            ::memset(D_vector,  0, pk_pairs_ * sizeof(double));
            D_vectors.push_back(D_vector);
            // The off-diagonal terms need to be doubled here
            size_t pqval = 0;
            for (int h = 0; h < nirrep_; ++h) {
                for (int p = 0; p < nsopi_[h]; ++p) {
                    for (int q = 0; q <= p; ++q) {
                        if (p != q) {
                            D_vector[pqval] = 2.0 * D_[N]->get(h, p, q);
                        }else{
                            D_vector[pqval] = D_[N]->get(h, p, q);
                        }
                        ++pqval;
                    }
                }
            }
        }
        J_[N]->zero();
        if(K_.size()) K_[N]->zero();
    }

    if(J_.size() || (K_.size() && asym)) {
        for(int batch = 0; batch < nbatches; ++batch){
            size_t min_pq      = batch_pq_min_[batch];
            size_t max_pq      = batch_pq_max_[batch];
            size_t min_index   = batch_index_min_[batch];
            size_t max_index   = batch_index_max_[batch];
            size_t batch_size = max_index - min_index;
            double *j_block = new double[batch_size];

            char *label = new char[100];
            sprintf(label, "J Block (Batch %d)", batch);
            psio_->read_entry(pk_file_, label, (char*) j_block, batch_size * sizeof(double));

            int nvectors = J_.size();
            for(int N = 0; N < nvectors; ++N){
                if (C_left_[N] != C_right_[N]) {
                    double **J_vector = J_[N]->pointer();
                    double *D_vec = D_vectors[N];
                    double *j_ptr = j_block;
                    int p, q, r, s;
                    for (size_t pq = min_pq; pq < max_pq; ++pq) {
                        for (size_t rs = 0; rs <= pq; ++rs) {
                            p = pind[pq];
                            q = qind[pq];
                            r = pind[rs];
                            s = qind[rs];
                            J_vector[p][q] += *j_ptr * (D_vec[r * nsopi_[0] + s] + D_vec[s * nsopi_[0] + r]);
                            J_vector[q][p] += *j_ptr * (D_vec[r * nsopi_[0] + s] + D_vec[s * nsopi_[0] + r]);
                            J_vector[r][s] += *j_ptr * (D_vec[p * nsopi_[0] + q] + D_vec[q * nsopi_[0] + p]);
                            J_vector[s][r] += *j_ptr * (D_vec[p * nsopi_[0] + q] + D_vec[q * nsopi_[0] + p]);
                            ++j_ptr;
                        }
                    }
                    if(K_.size()) {
                        double **K_vector = K_[N]->pointer();
                        j_ptr = j_block;
                        double fac;
                        // Each integral has a symmetry factor. Ugly in-loop condition
                        // is ugly. This is only a quickfix to get asymmetric density matrices.
                        // Better code to come later. -JFG
                        for (size_t pq = min_pq; pq < max_pq; ++pq) {
                            for (size_t rs = 0; rs <= pq; ++rs) {
                                fac = 1.0;
                                p = pind[pq];
                                q = qind[pq];
                                r = pind[rs];
                                s = qind[rs];
                                if (p == q && r == s && p == r) {
                                    fac = 0.25; // Divide only be 4, PK stores integral with a
                                    // factor 0.5 on the (pq|pq) diagonal.
                                } else if ( (p == q && q == r) || (q == r && r == s) ) {
                                    fac = 0.5;
                                } else if ( p == q && r == s) {
                                    fac = 0.25;
                                } else if (p == q || r == s) {
                                    fac = 0.5;
                                }
                                K_vector[p][r] += fac * (*j_ptr) * D_[N]->get(q, s);
                                K_vector[r][p] += fac * (*j_ptr) * D_[N]->get(s, q);
                                K_vector[q][r] += fac * (*j_ptr) * D_[N]->get(p, s);
                                K_vector[p][s] += fac * (*j_ptr) * D_[N]->get(q, r);
                                K_vector[s][p] += fac * (*j_ptr) * D_[N]->get(r, q);
                                K_vector[r][q] += fac * (*j_ptr) * D_[N]->get(s, p);
                                K_vector[s][q] += fac * (*j_ptr) * D_[N]->get(r, p);
                                K_vector[q][s] += fac * (*j_ptr) * D_[N]->get(p, r);
                                ++j_ptr;
                            }
                        }
                    }

                } else if (J_.size()) {
                    double *D_vector = D_vectors[N];
                    double *J_vector = J_vectors[N];
                    double *j_ptr = j_block;
                    for (size_t pq = min_pq; pq < max_pq; ++pq) {
                        double D_pq = D_vector[pq];
                        double *D_rs = D_vector;
                        double J_pq = 0.0;
                        double *J_rs = J_vector;
                        for (size_t rs = 0; rs <= pq; ++rs) {
//                            outfile->Printf("PK int (%lu|%lu) = %20.16f\n", pq, rs, *j_ptr);
                            J_pq  += *j_ptr * (*D_rs);
                            *J_rs += *j_ptr * D_pq;
                            ++D_rs;
                            ++J_rs;
                            ++j_ptr;
                        }
                        J_vector[pq] += J_pq;
                    }
                }
            }
            delete[] label;
            delete[] j_block;
        }
    }
//    throw PSIEXCEPTION("Reading PK file\n");

    for(size_t N = 0; N < J_.size(); ++N){
        if(C_left_[N] != C_right_[N]) {
            for(int p = 0; p < nsopi_[0]; ++p) {
                J_[N]->set(p, p, J_[N]->get(p,p) * 0.5);
            }
            delete [] D_vectors[N];
        } else {
            // Copy the results from the vector to the buffer
            double *J = J_vectors[N];
            for (int h = 0; h < nirrep_; ++h) {
                for (int p = 0; p < nsopi_[h]; ++p) {
                    for (int q = 0; q <= p; ++q) {
                        J_[N]->set(h, p, q, *J++);
                    }
                }
            }
            J_[N]->copy_lower_to_upper();
            delete [] D_vectors[N];
            delete [] J_vectors[N];
        }
    }

    /*
     * The K terms
     */
    D_vectors.clear();
    for(size_t N = 0; N < K_.size(); ++N){
        if(C_left_[N] != C_right_[N]) {
            // Fill in dummy vectors so that the ordering is not messed up.
            K_vectors.push_back(NULL);
            D_vectors.push_back(NULL);
        } else {
            double *K_vector = new double[pk_pairs_];
            ::memset(K_vector,  0, pk_pairs_ * sizeof(double));
            K_vectors.push_back(K_vector);
            double *D_vector = new double[pk_pairs_];
            ::memset(D_vector,  0, pk_pairs_ * sizeof(double));
            D_vectors.push_back(D_vector);
            // The off-diagonal terms need to be doubled here
            size_t pqval = 0;
            for (int h = 0; h < nirrep_; ++h) {
                for (int p = 0; p < nsopi_[h]; ++p) {
                    for (int q = 0; q <= p; ++q) {
                        if (p != q) {
                            D_vector[pqval] = 2.0 * D_[N]->get(h, p, q);
                        }else{
                            D_vector[pqval] = D_[N]->get(h, p, q);
                        }
                        ++pqval;
                    }
                }
            }
        }
    }

    if(K_.size() && !asym_only){
        for(int batch = 0; batch < nbatches; ++batch){
            size_t min_pq      = batch_pq_min_[batch];
            size_t max_pq      = batch_pq_max_[batch];
            size_t min_index   = batch_index_min_[batch];
            size_t max_index   = batch_index_max_[batch];
            size_t batch_size = max_index - min_index;
            double *k_block = new double[batch_size];

            char *label = new char[100];
            sprintf(label, "K Block (Batch %d)", batch);
            psio_->read_entry(pk_file_, label, (char*) k_block, batch_size * sizeof(double));

            int nvectors = K_.size();
            for(int N = 0; N < nvectors; ++N){
                if (C_left_[N] == C_right_[N]) {
                    double *D_vector = D_vectors[N];
                    double *K_vector = K_vectors[N];
                    double *k_ptr = k_block;
                    for (size_t pq = min_pq; pq < max_pq; ++pq) {
                        double D_pq = D_vector[pq];
                        double *D_rs = D_vector;
                        double K_pq = 0.0;
                        double *K_rs = K_vector;
                        for (size_t rs = 0; rs <= pq; ++rs) {
//                            outfile->Printf("PK int (%lu|%lu) = %20.16f\n", pq, rs, *k_ptr);
                            K_pq  += *k_ptr * (*D_rs);
                            *K_rs += *k_ptr * D_pq;
                            ++D_rs;
                            ++K_rs;
                            ++k_ptr;
                        }
                        K_vector[pq] += K_pq;
                    }
                }
            }
            delete[] label;
            delete[] k_block;
        }
//        throw PSIEXCEPTION("Reading K blocks\n");
    }

    for(size_t N = 0; N < K_.size(); ++N){
        if(C_left_[N] == C_right_[N]) {
            // Copy the results from the vector to the buffer
            double *K = K_vectors[N];
            for (int h = 0; h < nirrep_; ++h) {
                for (int p = 0; p < nsopi_[h]; ++p) {
                    for (int q = 0; q <= p; ++q) {
                        K_[N]->set(h, p, q, *K++);
                    }
                }
            }
            K_[N]->copy_lower_to_upper();
            delete [] D_vectors[N];
            delete [] K_vectors[N];
        }
    }

    /*
     * The wK terms
     */
    std::vector<double*> wK_vectors;
    D_vectors.clear();
    for(size_t N = 0; N < wK_.size(); ++N){
        if(D_[N]->symmetry()) {
            throw PSIEXCEPTION("PK integrals cannot be used for this type of calculation.");
        }
        double *K_vector = new double[pk_pairs_];
        ::memset(K_vector,  0, pk_pairs_ * sizeof(double));
        wK_vectors.push_back(K_vector);
        double *D_vector = new double[pk_pairs_];
        ::memset(D_vector,  0, pk_pairs_ * sizeof(double));
        D_vectors.push_back(D_vector);
        // The off-diagonal terms need to be doubled here
        size_t pqval = 0;
        for (int h = 0; h < nirrep_; ++h) {
            for (int p = 0; p < nsopi_[h]; ++p) {
                for (int q = 0; q <= p; ++q) {
                    if (p != q) {
                        D_vector[pqval] = 2.0 * D_[N]->get(h, p, q);
                    }else{
                        D_vector[pqval] = D_[N]->get(h, p, q);
                    }
                    ++pqval;
                }
            }
        }
    }

    if(wK_.size()){
        for(int batch = 0; batch < nbatches; ++batch){
            size_t min_pq      = batch_pq_min_[batch];
            size_t max_pq      = batch_pq_max_[batch];
            size_t min_index   = batch_index_min_[batch];
            size_t max_index   = batch_index_max_[batch];
            size_t batch_size = max_index - min_index;
            double *k_block = new double[batch_size];

            char *label = new char[100];
            sprintf(label, "wK Block (Batch %d)", batch);
            psio_->read_entry(pk_file_, label, (char*) k_block, batch_size * sizeof(double));

            int nvectors = wK_.size();
            for(int N = 0; N < nvectors; ++N){
                double *D_vector = D_vectors[N];
                double *K_vector = wK_vectors[N];
                double *k_ptr = k_block;
                for (size_t pq = min_pq; pq < max_pq; ++pq) {
                    double D_pq = D_vector[pq];
                    double *D_rs = D_vector;
                    double K_pq = 0.0;
                    double *K_rs = K_vector;
                    for (size_t rs = 0; rs <= pq; ++rs) {
                        K_pq  += *k_ptr * (*D_rs);
                        *K_rs += *k_ptr * D_pq;
                        ++D_rs;
                        ++K_rs;
                        ++k_ptr;
                    }
                    K_vector[pq] += K_pq;
                }
            }
            delete[] label;
            delete[] k_block;
        }
    }

    for(size_t N = 0; N < wK_.size(); ++N){
        // Copy the results from the vector to the buffer
        double *K = wK_vectors[N];
        for (int h = 0; h < nirrep_; ++h) {
            for (int p = 0; p < nsopi_[h]; ++p) {
                for (int q = 0; q <= p; ++q) {
                    wK_[N]->set(h, p, q, *K++);
                }
            }
        }
        wK_[N]->copy_lower_to_upper();

        delete [] D_vectors[N];
        delete [] wK_vectors[N];
    }

//    if(!file_was_open);
        psio_->close(pk_file_, 1);
}


//TODO: We may not need this anymore.
void PKJK::postiterations()
{
    if(algo_ != REORDER) {
        delete[] so2symblk_;
        delete[] so2index_;
    }
}
}
