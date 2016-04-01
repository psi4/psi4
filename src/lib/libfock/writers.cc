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

#include <cstdio>

#include <libpsio/psio.h>
#include <libqt/qt.h>
#include "writers.h"
#include <libmints/integral.h>
#include <psi4-dec.h>
#include <psifiles.h>
#include <libmints/mints.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// Defining the static maps here

std::map<psi::Value*, bool> psi::IWLAsync::busy_;
std::map<int, unsigned long int> psi::IWLAsync::bytes_written_;

namespace psi {

IWLAIOWriter::IWLAIOWriter(IWLAsync &writeto) : writeto_(writeto),
    count_(0) {
}

void IWLAIOWriter::operator ()(int i, int j, int k, int l, int, int, int, int,
                               int, int, int, int, double value) {
    // Fill the values, the internal integral counter of the buffer is
    // automatically incremented
    writeto_.fill_values(i, j, k, l, value);
    // Increment global counter
    count_++;
    // If buffer is full, dump it to disk

    if(writeto_.index() == writeto_.ints_per_buffer()) {
        writeto_.buffer_count() = writeto_.ints_per_buffer();
        writeto_.last_buffer() = 0;
        writeto_.put();
        writeto_.index() = 0;
    }
}

IWLAsync::IWLAsync(boost::shared_ptr<PSIO> psio, boost::shared_ptr<AIOHandler> aio, int itap)
{

    /*! Set up buffer info */
    itap_ = itap;
    bufpos_ = PSIO_ZERO;
    AIO_ = aio;
    ints_per_buf_ = IWL_INTS_PER_BUF;
    bufszc_ = 2 * sizeof(int) + 2 * ints_per_buf_ * 4 * sizeof(Label) +
            2 * ints_per_buf_ * sizeof(Value);
    lastbuf_[0] = 0;
    lastbuf_[1] = 0;
    inbuf_[0] = 0;
    inbuf_[1] = 0;
    idx_ = 0;
    whichbuf_ = 0;
    psio_ = psio;
    keep_ = true;

    /*! Allocate necessary space for two buffers */
    labels_[0] = new Label[4 * ints_per_buf_];
    values_[0] = new Value[ints_per_buf_];
    labels_[1] = new Label[4 * ints_per_buf_];
    values_[1] = new Value[ints_per_buf_];
    IWLAsync::busy_[values_[0]] = false;
    IWLAsync::busy_[values_[1]] = false;
    IWLAsync::bytes_written_[itap_] = 0;
}

IWLAsync::~IWLAsync() {
//    if(psio_->open_check(itap_))
//        psio_->close(itap_, keep_);
    if(labels_[0]) delete[] labels_[0];
    if(values_[0]) delete[] values_[0];
    if(labels_[1]) delete[] labels_[1];
    if(values_[1]) delete[] values_[1];
    labels_[0] = NULL;
    values_[0] = NULL;
    labels_[1] = NULL;
    values_[1] = NULL;
}

void IWLAsync::open_file(int oldfile) {
    psio_->open(itap_, oldfile ? PSIO_OPEN_OLD : PSIO_OPEN_NEW);
    if (oldfile && (psio_->tocscan(itap_, IWL_KEY_BUF) == NULL)) {
        outfile->Printf("IWLAsync::open_file: Can't find IWL buffers in file %d\n", itap_);
        psio_->close(itap_, 0);
        return;
    }
}

/*! Functions setting and getting values */
Label IWLAsync::get_p() {
    return labels_[whichbuf_][4 * idx_];
}
Label IWLAsync::get_q() {
    return labels_[whichbuf_][4 * idx_ + 1];
}
Label IWLAsync::get_r() {
    return labels_[whichbuf_][4 * idx_ + 2];
}
Label IWLAsync::get_s() {
    return labels_[whichbuf_][4 * idx_ + 3];
}
Value IWLAsync::get_val() {
    return values_[whichbuf_][idx_];
}

void IWLAsync::set_p(Label input) {
    labels_[whichbuf_][4 * idx_] = input;
}
void IWLAsync::set_q(Label input) {
    labels_[whichbuf_][4 * idx_ + 1] = input;
}
void IWLAsync::set_r(Label input) {
    labels_[whichbuf_][4 * idx_ + 2] = input;
}
void IWLAsync::set_s(Label input) {
    labels_[whichbuf_][4 * idx_ + 3] = input;
}
void IWLAsync::set_val(Value input) {
    values_[whichbuf_][idx_] = input;
}

/*! Sets all values and increment number of integrals in buffer */
void IWLAsync::fill_values(Label p, Label q, Label r, Label s, Value val) {
    if(idx_ >= ints_per_buf_) {
        outfile->Printf("Error: idx_ is %d\n", idx_);
    }
    set_p(p);
    set_q(q);
    set_r(r);
    set_s(s);
    set_val(val);
    ++idx_;
}

/*! Function that actually performs the asynchronous I/O and the
 * management of the buffers */

void IWLAsync::put() {

    // We probably need to put everything in there in a critical section

    // Let's start the asynchronous writing of the current buffer
    // We need to explicitly compute all starting addresses since we
    // will not be done writing by the time we place the next write in
    // the queue

    psio_address start1;
    psio_address start2;
    psio_address start3;
    psio_address start4;

    unsigned long int total_size = 2 * sizeof(int) + ints_per_buf_ * (4 * sizeof(Label) + sizeof(Value));

    int next_buf;

    next_buf = (whichbuf_ == 0) ? 1 : 0;

    // Now we make sure that the buffer in which we are going to write is not
    // currently being written
    bool busy;
    int thread = 0;
    #ifdef _OPENMP
      thread = omp_get_thread_num();
    #endif
    //printf("Thread number %d entering critical\n", thread);
    #pragma omp critical(IWLASYNC)
    // Could do that through getting a Job ID from AIOHandler
    busy = IWLAsync::busy_[values_[next_buf]];
    #pragma omp critical(IWLASYNC)
    if(busy) {
        //printf("Current buffer busy, synchronizing\n");
        //#pragma omp critical(IWLASYNC)
        AIO_->synchronize();
        //printf("Synch done, noone is busy any more\n");
        std::map<Value*,bool>::iterator it;
//#pragma omp critical(IWLASYNC)
        {
        for(it = IWLAsync::busy_.begin(); it != IWLAsync::busy_.end(); ++it) {
          it->second = false;
        }
        }
    }

#pragma omp critical(IWLWRITING)
    {

    //printf("Bytes written is now %lu\n",bytes_written_[itap_]);
    start1 = psio_get_address(PSIO_ZERO, bytes_written_[itap_]);
//#pragma omp atomic
    bytes_written_[itap_] += total_size;

    start2 = psio_get_address(start1, sizeof(int));
    start3 = psio_get_address(start2, sizeof(int));
    start4 = psio_get_address(start3, ints_per_buf_ * 4 * sizeof(Label));
    //printf("Writing lastbuf now\n");
//#pragma omp critical(IWLWRITING)
//    {
    AIO_->write(itap_, IWL_KEY_BUF, (char *) &(lastbuf_[whichbuf_]),
                sizeof(int), start1, &bufpos_);
    //bytes_written_[itap_] += sizeof(int);

    //printf("Writing inbuf now\n");
    //printf("Bytes written is now %lu\n",bytes_written_[itap_]);
    //start = psio_get_address(PSIO_ZERO, bytes_written_[itap_]);
    AIO_->write(itap_, IWL_KEY_BUF, (char*) &(inbuf_[whichbuf_]),
                sizeof(int), start2, &bufpos_);
    //bytes_written_[itap_] += sizeof(int);

    //printf("Writing labels now\n");
    //printf("Bytes written is now %lu\n",bytes_written_[itap_]);
    //start = psio_get_address(PSIO_ZERO, bytes_written_[itap_]);
    AIO_->write(itap_,IWL_KEY_BUF, (char *)(labels_[whichbuf_]),
                ints_per_buf_ * 4 * sizeof(Label), start3, &bufpos_);
    //bytes_written_[itap_] += ints_per_buf_ * 4 * sizeof(Label);

    //printf("Writing values now\n");
    //printf("Bytes written is now %lu\n",bytes_written_[itap_]);
    //start = psio_get_address(PSIO_ZERO, bytes_written_[itap_]);
    AIO_->write(itap_, IWL_KEY_BUF, (char *)(values_[whichbuf_]),
                ints_per_buf_ * sizeof(Value), start4, &bufpos_);
    //bytes_written_[itap_] += ints_per_buf_ * sizeof(Value);
//    }

//#pragma omp critical(IWLASYNC)
    IWLAsync::busy_[values_[whichbuf_]] = true;
    //printf("Thread number %d exiting critical\n", thread);
    }
    // Now we change the value to write in the other buffer
    whichbuf_ = next_buf;

}

void IWLAsync::flush(int lastbuf) {
    while(idx_ < ints_per_buf_) {
        fill_values(0, 0, 0, 0, 0.0);
    }

    if (lastbuf) last_buffer() = 1;
    else last_buffer() = 0;

    inbuf_[whichbuf_] = idx_;
    put();
    idx_ = 0;
}

void ijklBasisIterator::first() {
    i_ = 0;
    j_ = 0;
    k_ = 0;
    l_ = 0;
}

void ijklBasisIterator::next() {
    ++l_;
    if (l_ > j_ && k_ == i_ ) {
        l_ = 0;
        ++k_;
    }
    if (l_ > k_) {
        l_ = 0;
        ++k_;
    }
    if (k_ > i_) {
        k_ = 0;
        ++j_;
        if (j_ > i_) {
            j_ = 0;
            ++i_;
            if (i_ >= nbas_) {
                done_=true;
            }
        }
    }
}


PK_integrals::PK_integrals(boost::shared_ptr<BasisSet> primary, boost::shared_ptr<PSIO> psio,
                           int max_batches, size_t memory) {
    primary_ = primary;
    nbf_ = primary_->nbf();
    max_batches_ = max_batches;
    memory_ = memory;

    pk_pairs_ = nbf_ * (nbf_ + 1) / 2;
    pk_size_ = pk_pairs_ * (pk_pairs_ + 1) / 2;

    if (memory_ < pk_pairs_) {
        throw PSIEXCEPTION("Not enough memory for PK algorithm");
    }

    J_buf_[0] = NULL;
    J_buf_[1] = NULL;
    K_buf_[0] = NULL;
    K_buf_[1] = NULL;

    nbuffers_ = 0;

    psio_ = psio;
    AIO_ = boost::shared_ptr<AIOHandler>(new AIOHandler(psio_));

    itap_J_ = PSIF_SO_PKSUPER1;
    itap_K_ = PSIF_SO_PKSUPER2;

}

PK_integrals::~PK_integrals() {
    // In case we forgot something before

    if (J_buf_[0] != NULL) {
        delete [] J_buf_[0];
        J_buf_[0] = NULL;
        outfile->Printf("Clean up your code!!!");
    }
    if (J_buf_[1] != NULL) {
        delete [] J_buf_[1];
        J_buf_[1] = NULL;
        outfile->Printf("Clean up your code!!!");
    }
    if (K_buf_[0] != NULL) {
        delete [] K_buf_[0];
        K_buf_[0] = NULL;
        outfile->Printf("Clean up your code!!!");
    }
    if (K_buf_[1] != NULL) {
        delete [] K_buf_[1];
        K_buf_[1] = NULL;
        outfile->Printf("Clean up your code!!!");
    }
    if (label_J_[0].size() != 0) {
        outfile->Printf("Labels: Clean up your code!!!");
    }
    if (label_J_[1].size() != 0) {
        outfile->Printf("Labels: Clean up your code!!!");
    }
    if (label_K_[0].size() != 0) {
        outfile->Printf("Labels: Clean up your code!!!");
    }
    if (label_K_[1].size() != 0) {
        outfile->Printf("Labels: Clean up your code!!!");
    }
}

void PK_integrals::batch_sizing() {

    double batch_thresh = 0.1;

    ijklBasisIterator AOintsiter(nbf_);

    size_t old_pq = 0;
    size_t old_max = 0;
    size_t nintpq = 0;
    size_t nintbatch = 0;
    size_t min_pqrs = 0;
    size_t min_pq = 0;
    size_t pq = 0;
    size_t pb, qb, rb, sb;

    batch_index_min_.push_back(0);
    batch_pq_min_.push_back(0);
    for (AOintsiter.first(); AOintsiter.is_done() == false; AOintsiter.next()) {
        pb = AOintsiter.i();
        qb = AOintsiter.j();
        rb = AOintsiter.k();
        sb = AOintsiter.l();

        pq = INDEX2(pb, qb);

        if (old_pq == pq) {
            ++nintpq;
        } else {
            size_t pqrs = INDEX2(pq, INDEX2(rb, sb));
            nintbatch += nintpq;
            if (nintbatch > memory_) {
                batch_index_max_.push_back(old_max);
                batch_pq_max_.push_back(old_pq);
                batch_index_min_.push_back(old_max);
                batch_pq_min_.push_back(old_pq);
                nintbatch = nintpq;
            }
            nintpq = 1;
            old_pq = pq;
            old_max = pqrs;
        }
    }
    batch_index_max_.push_back(INDEX4(pb,qb,rb,sb) + 1);
    batch_pq_max_.push_back(INDEX2(pb,qb) + 1);

    // A little check here: if the last batch is less than 10% full,
    // we just transfer it to the previous batch.

    int lastb = batch_index_max_.size() - 1;
    if (lastb > 0) {
    size_t size_lastb = batch_index_max_[lastb] - batch_index_min_[lastb];
        if (((double) size_lastb / memory_) < batch_thresh) {
            batch_index_max_[lastb - 1] = batch_index_max_[lastb];
            batch_pq_max_[lastb - 1] = batch_pq_max_[lastb];
            batch_pq_max_.pop_back();
            batch_index_max_.pop_back();
            batch_pq_min_.pop_back();
            batch_index_min_.pop_back();
        }
    }

    int nbatches = batch_pq_min_.size();
    if (nbatches > max_batches_) {
      outfile->Printf( "PKJK: maximum number of batches exceeded\n") ;
      outfile->Printf( "   requested %d batches\n", nbatches) ;
      exit(PSI_RETURN_FAILURE) ;
    }

}

void PK_integrals::print_batches() {
    // Print batches for the user and for control
    for(int batch = 0; batch < batch_pq_min_.size(); ++batch){
        outfile->Printf("\tBatch %3d pq = [%8zu,%8zu] index = [%14zu,%zu] size = %12zu\n",
                batch + 1,
                batch_pq_min_[batch],batch_pq_max_[batch],
                batch_index_min_[batch],batch_index_max_[batch],
                batch_index_max_[batch] - batch_index_min_[batch]);
    }
}

void PK_integrals::allocate_buffers() {
    buffer_size_ = memory_ / 4;
    nbuffers_ = (pk_size_ - 1) / buffer_size_ + 1;
    outfile->Printf("We need %3d buffers during integral computation\n",nbuffers_);
    // Now we proceed to actually allocate buffer memory
    J_buf_[0] = new double[buffer_size_];
    J_buf_[1] = new double[buffer_size_];
    K_buf_[0] = new double[buffer_size_];
    K_buf_[1] = new double[buffer_size_];
    ::memset((void*) J_buf_[0], '\0', buffer_size_ * sizeof(double));
    ::memset((void*) J_buf_[1], '\0', buffer_size_ * sizeof(double));
    ::memset((void*) K_buf_[0], '\0', buffer_size_ * sizeof(double));
    ::memset((void*) K_buf_[1], '\0', buffer_size_ * sizeof(double));

    bufidx_ = 0;
    offset_ = 0;
    buf_ = 0;

    jobid_J_[0] = 0;
    jobid_J_[1] = 0;
    jobid_K_[0] = 0;
    jobid_K_[1] = 0;
}

void PK_integrals::deallocate_buffers() {

    // We need to make sure writing is over before deallocating!
    timer_on("AIO synchronize");
    AIO_->synchronize();
    timer_off("AIO synchronize");

    delete [] J_buf_[0];
    delete [] J_buf_[1];
    delete [] K_buf_[0];
    delete [] K_buf_[1];

    J_buf_[0] = NULL;
    J_buf_[1] = NULL;
    K_buf_[0] = NULL;
    K_buf_[1] = NULL;

    for(int i = 0; i < 2; ++i) {
        for(int j = 0; j < label_J_[i].size(); ++j) {
            delete[] label_J_[i][j];
        }
        for(int j = 0; j < label_K_[i].size(); ++j) {
            delete[] label_K_[i][j];
        }
        label_J_[i].clear();
        label_K_[i].clear();
    }

}

size_t PK_integrals::task_quartets() {
    // Need iterator over AO shells here
    AOShellCombinationsIterator shelliter(primary_,primary_,primary_,primary_);
    size_t nsh_task = 0;
    for(shelliter.first(); shelliter.is_done() == false; shelliter.next()) {
        int P = shelliter.p();
        int Q = shelliter.q();
        int R = shelliter.r();
        int S = shelliter.s();
        // Get the lowest basis function indices in the shells

        size_t lowi = primary_->shell_to_basis_function(P);
        size_t lowj = primary_->shell_to_basis_function(Q);
        size_t lowk = primary_->shell_to_basis_function(R);
        size_t lowl = primary_->shell_to_basis_function(S);

        size_t low_ijkl = INDEX4(lowi, lowj, lowk, lowl);
        size_t low_ikjl = INDEX4(lowi, lowk, lowj, lowl);
        size_t low_iljk = INDEX4(lowi, lowl, lowj, lowk);

        size_t nbJ_lo = low_ijkl / buffer_size_;
        size_t nbK1_lo = low_ikjl / buffer_size_;
        size_t nbK2_lo = low_iljk / buffer_size_;

        // Can we filter out this shell because all of its basis function
        // indices are too high ?

        if (nbJ_lo > bufidx_ && nbK1_lo > bufidx_ && nbK2_lo > bufidx_) {
            continue;
        }

        int ni = primary_->shell(P).nfunction();
        int nj = primary_->shell(Q).nfunction();
        int nk = primary_->shell(R).nfunction();
        int nl = primary_->shell(S).nfunction();

        size_t hii = lowi + ni - 1;
        size_t hij = lowj + nj - 1;
        size_t hik = lowk + nk - 1;
        size_t hil = lowl + nl - 1;

        size_t hi_ijkl = INDEX4(hii, hij, hik, hil);
        size_t hi_ikjl = INDEX4(hii, hik, hij, hil);
        size_t hi_iljk = INDEX4(hii, hil, hij, hik);

        size_t nbJ_hi = hi_ijkl / buffer_size_;
        size_t nbK1_hi = hi_ikjl / buffer_size_;
        size_t nbK2_hi = hi_iljk / buffer_size_;

        // Can we filter out this shell because all of its basis function
        // indices are too low ?

        if (nbJ_hi < bufidx_ && nbK1_hi < bufidx_ && nbK2_hi < bufidx_) {
            continue;
        }

        // Now we loop over unique basis function quartets in the shell quartet
        AOIntegralsIterator bfiter = shelliter.integrals_iterator();
        for(bfiter.first(); bfiter.is_done() == false; bfiter.next()) {
            size_t i = bfiter.i();
            size_t j = bfiter.j();
            size_t k = bfiter.k();
            size_t l = bfiter.l();

            size_t ijkl = INDEX4(i,j,k,l);
            size_t ikjl = INDEX4(i,k,j,l);
            size_t iljk = INDEX4(i,l,j,k);

            size_t nbJ = ijkl / buffer_size_;
            size_t nbK1 = ikjl / buffer_size_;
            size_t nbK2 = iljk / buffer_size_;

            if (nbJ == bufidx_ || nbK1 == bufidx_ || nbK2 == bufidx_) {
                buf_P.push_back(P);
                buf_Q.push_back(Q);
                buf_R.push_back(R);
                buf_S.push_back(S);
                ++nsh_task;
                break;
            }
        }
    }
    return nsh_task;
}

// This function is thread-safe
void PK_integrals::integrals_buffering(const double *buffer, int P, int Q, int R, int S) {

    AOIntegralsIterator bfiter(primary_->shell(P), primary_->shell(Q), primary_->shell(R), primary_->shell(S));
    int i0 = primary_->shell_to_basis_function(P);
    int j0 = primary_->shell_to_basis_function(Q);
    int k0 = primary_->shell_to_basis_function(R);
    int l0 = primary_->shell_to_basis_function(S);

    int ni = primary_->shell(P).nfunction();
    int nj = primary_->shell(Q).nfunction();
    int nk = primary_->shell(R).nfunction();
    int nl = primary_->shell(S).nfunction();

    for (bfiter.first(); bfiter.is_done() == false; bfiter.next()) {
        int i = bfiter.i();
        int j = bfiter.j();
        int k = bfiter.k();
        int l = bfiter.l();
        size_t idx = bfiter.index();

        double val = buffer[idx];
        fill_values(val, i, j, k, l);
    }
}

// This function is thread-safe
void PK_integrals::fill_values(double val, size_t i, size_t j, size_t k, size_t l) {
    size_t ijkl = INDEX4(i, j, k, l);
    if (ijkl / buffer_size_ == bufidx_) {
#pragma omp atomic
        J_buf_[buf_][ijkl - offset_] += val;
    }

    size_t ikjl = INDEX4(i, k, j, l);
    if(ikjl / buffer_size_ == bufidx_) {
        if (i == k || j == l) {
#pragma omp atomic
            K_buf_[buf_][ikjl - offset_] += val;
        } else {
#pragma omp atomic
            K_buf_[buf_][ikjl - offset_] += 0.5 * val;
        }
    }

    if(i != j && k != l) {
        size_t iljk = INDEX4(i, l, j, k);
        if (iljk / buffer_size_ == bufidx_) {
            if ( i == l || j == k) {
#pragma omp atomic
                K_buf_[buf_][iljk - offset_] += val;
            } else {
#pragma omp atomic
                K_buf_[buf_][iljk - offset_] += 0.5 * val;
            }
        }
    }
}

void PK_integrals::open_files(bool old) {
    psio_->open(itap_J_, old ? PSIO_OPEN_OLD : PSIO_OPEN_NEW);
    psio_->open(itap_K_, old ? PSIO_OPEN_OLD : PSIO_OPEN_NEW);
}

// This function is not designed thread-safe and is supposed to execute outside parallel
// environment for now.
void PK_integrals::write() {
    // Clear up the buffers holding P, Q, R, S since the task is over
    buf_P.clear();
    buf_Q.clear();
    buf_R.clear();
    buf_S.clear();
    // Compute initial ijkl index for current buffer
    size_t ijkl0 = buffer_size_ * bufidx_;
    size_t ijkl_end = ijkl0 + buffer_size_ - 1;

    std::vector<int> target_batches;  // Vector of batch number to which we want to write

    for (int i = 0; i < batch_index_min_.size(); ++i) {
        if (ijkl0 >= batch_index_min_[i] && ijkl0 < batch_index_max_[i]) {
            target_batches.push_back(i);
            continue;
        }
        if (ijkl_end >= batch_index_min_[i] && ijkl_end < batch_index_max_[i]) {
            target_batches.push_back(i);
            continue;
        }
        if (ijkl0 < batch_index_min_[i] && ijkl_end >= batch_index_max_[i]) {
            target_batches.push_back(i);
            continue;
        }
    }

    // Now that all buffers are full, and before we write integrals, we need
    // to divide diagonal elements by 2.

    for (size_t pq = 0; pq < pk_pairs_; ++pq) {
        size_t pqpq = INDEX2(pq, pq);
        if (pqpq >= ijkl0 && pqpq <= ijkl_end) {
            J_buf_[buf_][pqpq - offset_] *= 0.5;
            K_buf_[buf_][pqpq - offset_] *= 0.5;
        }
    }

    // And now we write to the files in the appropriate entries
    for (int i = 0; i < target_batches.size(); ++i) {
        char* label_J = new char[100];
        sprintf(label_J,"J Block (Batch %d)", target_batches[i]);
        label_J_[buf_].push_back(label_J);
        size_t start = std::max(ijkl0, batch_index_min_[target_batches[i]]);
        size_t stop = std::min(ijkl_end + 1, batch_index_max_[target_batches[i]]);
        psio_address adr = psio_get_address(PSIO_ZERO, (start - batch_index_min_[target_batches[i]]) * sizeof(double));
        size_t nints = stop - start;
        jobid_J_[buf_] = AIO_->write(itap_J_, label_J_[buf_][i], (char *)(&J_buf_[buf_][start - offset_]),
                    nints * sizeof(double), adr, &dummy_);
        char* label_K = new char[100];
        sprintf(label_K,"K Block (Batch %d)",target_batches[i]);
        label_K_[buf_].push_back(label_K);
        jobid_K_[buf_] = AIO_->write(itap_K_, label_K_[buf_][i], (char *)(&K_buf_[buf_][start - offset_]),
                    nints * sizeof(double), adr, &dummy_);
    }

    // Update the buffer being written into
    buf_ = buf_ == 0 ? 1 : 0;
    // Make sure the buffer has been written to disk and we can erase it
    AIO_->wait_for_job(jobid_J_[buf_]);
    AIO_->wait_for_job(jobid_K_[buf_]);
    // We can delete the labels for these buffers
    for(int i = 0; i < label_J_[buf_].size(); ++i) {
        delete [] label_J_[buf_][i];
    }
    for(int i = 0; i < label_K_[buf_].size(); ++i) {
        delete [] label_K_[buf_][i];
    }
    label_J_[buf_].clear();
    label_K_[buf_].clear();
    // Make sure the buffer in which we are going to write is set to zero
    ::memset((void *) J_buf_[buf_], '\0', buffer_size_ * sizeof(double));
    ::memset((void *) K_buf_[buf_], '\0', buffer_size_ * sizeof(double));
    // Update the offset for integral storage
    offset_ += buffer_size_;
    // And we are now writing to the next set of buffers
    ++bufidx_;

}

void PK_integrals::close_files() {
    psio_->close(itap_J_, 1);
    psio_->close(itap_K_, 1);
}

// We should really use safer pointers right here
void PK_integrals::form_D_vec(std::vector<SharedMatrix> D_ao) {
    // Assume symmetric density matrix for now.
    for (int N = 0; N < D_ao.size(); ++N) {
        double* D_vec = new double[pk_pairs_];
        ::memset((void *)D_vec,'\0',pk_pairs_ * sizeof(double));
        D_vec_.push_back(D_vec);
        size_t pqval = 0;
        for(int p = 0; p < nbf_; ++p) {
            for(int q = 0; q <= p; ++q) {
                if(p != q) {
                    D_vec[pqval] = 2.0 * D_ao[N]->get(0,p,q);
                } else {
                    D_vec[pqval] = D_ao[N]->get(0,p,q);
                }
                ++pqval;
            }
        }
    }
}

void PK_integrals::finalize_D() {
    for(int N = 0; N < D_vec_.size(); ++N) {
        delete [] D_vec_[N];
    }
    D_vec_.clear();
}

void PK_integrals::form_J(std::vector<SharedMatrix> J, bool exch) {
    std::vector<double*> Jvecs;
    // Begin by allocating vector for triangular J
    for(int N = 0; N < J.size(); ++N) {
        double* J_vec = new double[pk_pairs_];
        ::memset((void*)J_vec,'\0',pk_pairs_ * sizeof(double));
        Jvecs.push_back(J_vec);
    }

    // Now loop over batches
    for(int batch = 0; batch < batch_pq_min_.size(); ++batch) {
        size_t min_index = batch_index_min_[batch];
        size_t max_index = batch_index_max_[batch];
        size_t batch_size = max_index - min_index;
        size_t min_pq = batch_pq_min_[batch];
        size_t max_pq = batch_pq_max_[batch];
        double* j_block = new double[batch_size];

        int filenum;

        char* label = new char[100];
        if (exch) {
            sprintf(label,"K Block (Batch %d)", batch);
            filenum = itap_K_;
        } else {
            sprintf(label,"J Block (Batch %d)", batch);
            filenum = itap_J_;
        }
        psio_->read_entry(filenum, label, (char *) j_block, batch_size * sizeof(double));

        // Read one entry, use it for all density matrices
        for(int N = 0; N < J.size(); ++N) {
            double* D_vec = D_vec_[N];
            double* J_vec = Jvecs[N];
            double* j_ptr = j_block;
            for(size_t pq = min_pq; pq < max_pq; ++pq) {
                double D_pq = D_vec[pq];
                double *D_rs = D_vec;
                double J_pq = 0.0;
                double *J_rs = J_vec;
                for(size_t rs = 0; rs <= pq; ++rs) {
                    J_pq += *j_ptr * (*D_rs);
                    *J_rs += *j_ptr * D_pq;
                    ++D_rs;
                    ++J_rs;
                    ++j_ptr;
                }
                J_vec[pq] += J_pq;
            }
        }

        delete [] label;
        delete [] j_block;
    }

    // Now, directly transfer data to resulting matrices
    for(int N = 0; N < J.size(); ++N) {
        double *Jp = Jvecs[N];
        for(int p = 0; p < nbf_; ++p) {
            for(int q = 0; q <= p; ++q) {
                J[N]->set(0,p,q,*Jp++);
            }
        }
        J[N]->copy_lower_to_upper();
        delete [] Jvecs[N];
    }
}

void PK_integrals::form_K(std::vector<SharedMatrix> K) {
    // For asymmetric densities, we do exactly the same than for J
    // but we read another entry
    form_J(K, true);
}

}
