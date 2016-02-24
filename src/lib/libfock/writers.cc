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
#include "writers.h"
#include <psi4-dec.h>
#include <psifiles.h>
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

    psio_address start;

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
#pragma omp critical(IWLWRITING)
    {
    //printf("Thread number %d entering critical\n", thread);
//    #pragma omp critical(IWLASYNC)
    // Could do that through aio->get_thread()->joinable()
    busy = IWLAsync::busy_[values_[next_buf]];
    if(busy) {
        //printf("Current buffer busy, synchronizing\n");
        AIO_->synchronize();
        //printf("Synch done, noone is busy any more\n");
        std::map<Value*,bool>::iterator it;
        for(it = IWLAsync::busy_.begin(); it != IWLAsync::busy_.end(); ++it) {
//#pragma omp critical(IWLASYNC)
          it->second = false;
        }
    }


    //printf("Bytes written is now %lu\n",bytes_written_[itap_]);
    start = psio_get_address(PSIO_ZERO, bytes_written_[itap_]);
//#pragma omp atomic
    //bytes_written_[itap_] += total_size;

    //printf("Writing lastbuf now\n");
    AIO_->write(itap_, IWL_KEY_BUF, (char *) &(lastbuf_[whichbuf_]),
                sizeof(int), start, &bufpos_);
    bytes_written_[itap_] += sizeof(int);

    //printf("Writing inbuf now\n");
    //start = psio_get_address(start, sizeof(int));
    //printf("Bytes written is now %lu\n",bytes_written_[itap_]);
    start = psio_get_address(PSIO_ZERO, bytes_written_[itap_]);
    AIO_->write(itap_, IWL_KEY_BUF, (char*) &(inbuf_[whichbuf_]),
                sizeof(int), start, &bufpos_);
    bytes_written_[itap_] += sizeof(int);

    //printf("Writing labels now\n");
    //start = psio_get_address(start, sizeof(int));
    //printf("Bytes written is now %lu\n",bytes_written_[itap_]);
    start = psio_get_address(PSIO_ZERO, bytes_written_[itap_]);
    AIO_->write(itap_,IWL_KEY_BUF, (char *)(labels_[whichbuf_]),
                ints_per_buf_ * 4 * sizeof(Label), start, &bufpos_);
    bytes_written_[itap_] += ints_per_buf_ * 4 * sizeof(Label);

    //printf("Writing values now\n");
    //start = psio_get_address(start, ints_per_buf_ * 4 * sizeof(Label));
    //printf("Bytes written is now %lu\n",bytes_written_[itap_]);
    start = psio_get_address(PSIO_ZERO, bytes_written_[itap_]);
    AIO_->write(itap_, IWL_KEY_BUF, (char *)(values_[whichbuf_]),
                ints_per_buf_ * sizeof(Value), start, &bufpos_);
    bytes_written_[itap_] += ints_per_buf_ * sizeof(Value);

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

}
