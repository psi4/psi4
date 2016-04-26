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

#include <libmints/integral.h>
#include <libpsio/psio.h>
#include <libpsio/aiohandler.h>
#include <libiwl/config.h>
#include "PKmanagers.h"
#include "PK_workers.h"

namespace psi {

namespace pk {

//TODO Finish this constructor, of course
PKWorker::PKWorker(boost::shared_ptr<BasisSet> primary, boost::shared_ptr<AIOHandler> AIO,
                   int target_file, size_t buf_size) {
    AIO_ = AIO;
    target_file_ = target_file;
    primary_ = primary;
    buf_size_ = buf_size;
    bufidx_ = 0;
    offset_ = 0;

}

void PKWorker::first_quartet(size_t i) {
    shelliter_ = AOShellCombinationsIterator(primary_,primary_,primary_,primary_);
    bufidx_ = i;
    offset_ = bufidx_ * buf_size_;
    shells_left_ = false;
    for(shelliter_.first(); (shells_left_ || shelliter_.is_done()) == false; shelliter_.next()) {
        P_ = shelliter_.p();
        Q_ = shelliter_.q();
        R_ = shelliter_.r();
        S_ = shelliter_.s();

        shells_left_ = is_shell_relevant();
    }
    max_idx_ = get_max_idx();

}

bool PKWorker::is_shell_relevant() {
    // May implement the sieve here

    size_t lowi = primary_->shell_to_basis_function(P_);
    size_t lowj = primary_->shell_to_basis_function(Q_);
    size_t lowk = primary_->shell_to_basis_function(R_);
    size_t lowl = primary_->shell_to_basis_function(S_);

    size_t low_ijkl = INDEX4(lowi, lowj, lowk, lowl);
    size_t low_ikjl = INDEX4(lowi, lowk, lowj, lowl);
    size_t low_iljk = INDEX4(lowi, lowl, lowj, lowk);

    // Can we filter out this shell because all of its basis function
    // indices are too high ?

    if (low_ijkl > max_idx_ && low_ikjl > max_idx_ && low_iljk > max_idx_) {
        return false;
    }

    int ni = primary_->shell(P_).nfunction();
    int nj = primary_->shell(Q_).nfunction();
    int nk = primary_->shell(R_).nfunction();
    int nl = primary_->shell(S_).nfunction();

    size_t hii = lowi + ni - 1;
    size_t hij = lowj + nj - 1;
    size_t hik = lowk + nk - 1;
    size_t hil = lowl + nl - 1;

    size_t hi_ijkl = INDEX4(hii, hij, hik, hil);
    size_t hi_ikjl = INDEX4(hii, hik, hij, hil);
    size_t hi_iljk = INDEX4(hii, hil, hij, hik);

    // Use the min ijkl that can be in buffer, which is the offset_

    // Can we filter out this shell because all of its basis function
    // indices are too low ?

    if (hi_ijkl < offset_ && hi_ikjl < offset_ && hi_iljk < offset_) {
        return false;
    }

    // Now we loop over unique basis function quartets in the shell quartet
    AOIntegralsIterator bfiter = shelliter_.integrals_iterator();
    for(bfiter.first(); bfiter.is_done() == false; bfiter.next()) {
        size_t i = bfiter.i();
        size_t j = bfiter.j();
        size_t k = bfiter.k();
        size_t l = bfiter.l();

        size_t ijkl = INDEX4(i,j,k,l);
        size_t ikjl = INDEX4(i,k,j,l);
        size_t iljk = INDEX4(i,l,j,k);

        bool bJ = ijkl >= offset_ && ijkl <= max_idx_;
        bool bK1 = ikjl >= offset_ && ikjl <= max_idx_;
        bool bK2 = iljk >= offset_ && iljk <= max_idx_;

        if (bJ || bK1 || bK2) {
            // This shell should be computed by the present thread.
            return true;
        }
    }

    return false;

}

void PKWorker::next_quartet() {
    if(shelliter_.is_done()) {
        shells_left_ = false;
        return;
    }
    bool shell_found = false;
    for(; (shell_found || shelliter_.is_done()) == false; shelliter_.next()) {
        P_ = shelliter_.p();
        Q_ = shelliter_.q();
        R_ = shelliter_.r();
        S_ = shelliter_.s();

        shell_found = is_shell_relevant();
    }

    shells_left_ = shell_found;
}

PKWrkrReord::PKWrkrReord(boost::shared_ptr<BasisSet> primary, boost::shared_ptr<AIOHandler> AIO,
                         int target_file, size_t buf_size, unsigned int nbuf) :
    PKWorker(primary,AIO,target_file,buf_size) {

    nbuf_ = nbuf;
    buf_ = 0;

    double * mem;
    // Allocate everything we need
    for(int i = 0; i < nbuf_; ++i) {
        mem = new double[buf_size_];
        J_bufs_.push_back(mem);
        mem = new double[buf_size_];
        K_bufs_.push_back(mem);
        std::vector< char* > labelJ;
        std::vector< char* > labelK;
        labels_J_.push_back(labelJ);
        labels_K_.push_back(labelK);
        std::vector< size_t > idJ;
        std::vector< size_t > idK;
        jobID_J_.push_back(idJ);
        jobID_K_.push_back(idK);
    }
    // We only need to set to zero the first buffers.
    // The others are taken care of in the write function
    ::memset((void*) J_bufs_[0], '\0', buf_size_ * sizeof(double));
    ::memset((void*) K_bufs_[0], '\0', buf_size_ * sizeof(double));
}

PKWrkrReord::~PKWrkrReord() {
    // Deallocate all memory cleanly
    std::vector<double*>::iterator it;
    for(it = J_bufs_.begin(); it != J_bufs_.end(); ++it) {
        delete [] *it;
    }
    J_bufs_.clear();
    for(it = K_bufs_.begin(); it != K_bufs_.end(); ++it) {
        delete [] *it;
    }
    K_bufs_.clear();
    for(int i = 0; i < labels_J_.size(); ++i) {
        for(int j = 0; j < labels_J_[i].size(); ++j) {
            delete [] labels_J_[i][j];
        }
        labels_J_[i].clear();
    }
    labels_J_.clear();
    for(int i = 0; i < labels_K_.size(); ++i) {
        for(int j = 0; j < labels_K_[i].size(); ++j) {
            delete [] labels_K_[i][j];
        }
        labels_K_[i].clear();
    }
    labels_K_.clear();

}

size_t PKWrkrReord::get_max_idx() {
    max_idx_ = buf_size_ * (bufidx_ + 1) - 1;
}

void PKWrkrReord::fill_values(double val, size_t i, size_t j, size_t k, size_t l) {

    size_t ijkl = INDEX4(i,j,k,l);
    if (ijkl >= offset_ && ijkl <= max_idx_) {
        J_bufs_[buf_][ijkl - offset_] += val;
    }
    size_t ikjl = INDEX4(i, k, j, l);
    if(ikjl >= offset_ && ikjl <= max_idx_) {
        if (i == k || j == l) {
            K_bufs_[buf_][ikjl - offset_] += val;
        } else {
            K_bufs_[buf_][ikjl - offset_] += 0.5 * val;
        }
    }

    if(i != j && k != l) {
        size_t iljk = INDEX4(i, l, j, k);
        if (iljk >= offset_ && iljk <= max_idx_) {
            if ( i == l || j == k) {
                K_bufs_[buf_][iljk - offset_] += val;
            } else {
                K_bufs_[buf_][iljk - offset_] += 0.5 * val;
            }
        }
    }
}

void PKWrkrReord::write(std::vector<size_t> min_ind, std::vector<size_t> max_ind, size_t pk_pairs) {
    // Compute initial ijkl index for current buffer

    // Vector of batch number to which we want to write
    std::vector<unsigned int> target_batches;

    for (unsigned int i = 0; i < min_ind.size(); ++i) {
        if (offset_ >= min_ind[i] && offset_ < max_ind[i]) {
            target_batches.push_back(i);
            continue;
        }
        if (max_idx_ >= min_ind[i] && max_idx_ < max_ind[i]) {
            target_batches.push_back(i);
            continue;
        }
        if (offset_ < min_ind[i] && max_idx_ >= max_ind[i]) {
            target_batches.push_back(i);
            continue;
        }
    }

    // Now that all buffers are full, and before we write integrals, we need
    // to divide diagonal elements by 2.

    for (size_t pq = 0; pq < pk_pairs; ++pq) {
        size_t pqpq = INDEX2(pq, pq);
        if (pqpq >= offset_ && pqpq <= max_idx_) {
            J_bufs_[buf_][pqpq - offset_] *= 0.5;
            K_bufs_[buf_][pqpq - offset_] *= 0.5;
        }
    }

    // And now we write to the file in the appropriate entries
    for (int i = 0; i < target_batches.size(); ++i) {
        unsigned int b = target_batches[i];
        labels_J_[buf_].push_back(PKManager::get_label_J(b));
        size_t start = std::max(offset_, min_ind[b]);
        size_t stop = std::min(max_idx_ + 1, max_ind[b]);
        psio_address adr = psio_get_address(PSIO_ZERO, (start - min_ind[b]) * sizeof(double));
        size_t nints = stop - start;
        jobID_J_[buf_].push_back(AIO_->write(target_file_, labels_J_[buf_][i], (char *)(&J_bufs_[buf_][start - offset_]),
                    nints * sizeof(double), adr, &dummy_));
        labels_K_[buf_].push_back(PKManager::get_label_K(b));
        jobID_K_[buf_].push_back(AIO_->write(target_file_, labels_K_[buf_][i], (char *)(&K_bufs_[buf_][start - offset_]),
                    nints * sizeof(double), adr, &dummy_));
    }

    // Update the buffer being written into
    ++buf_;
    if (buf_ >= nbuf_) buf_ = 0;
    // Make sure the buffer has been written to disk and we can erase it
    for(int i = 0; i < jobID_J_[buf_].size(); ++i) {
      AIO_->wait_for_job(jobID_J_[buf_][i]);
    }
    jobID_J_[buf_].clear();
    for(int i = 0; i < jobID_K_[buf_].size(); ++i) {
      AIO_->wait_for_job(jobID_K_[buf_][i]);
    }
    jobID_K_[buf_].clear();
    // We can delete the labels for these buffers
    for(int i = 0; i < labels_J_[buf_].size(); ++i) {
        delete [] labels_J_[buf_][i];
    }
    for(int i = 0; i < labels_K_[buf_].size(); ++i) {
        delete [] labels_K_[buf_][i];
    }
    label_J_[buf_].clear();
    label_K_[buf_].clear();
    // Make sure the buffer in which we are going to write is set to zero
    ::memset((void *) J_bufs_[buf_], '\0', buf_size_ * sizeof(double));
    ::memset((void *) K_bufs_[buf_], '\0', buf_size_ * sizeof(double));

}

PKWrkrInCore::PKWrkrInCore(boost::shared_ptr<BasisSet> primary, size_t buf_size,
                           size_t lastbuf, double *Jbuf, double *Kbuf) :
    PKWorker(primary,NULL,0,buf_size) {

    last_buf_ = lastbuf;
    J_bufp_ = Jbuf;
    K_bufp_ = Kbuf;
}

size_t PKWrkrInCore::get_max_idx() {
    max_idx_ = (buf_size_ * (bufidx_ + 1) + last_buf_) - 1;
}

void PKWrkrInCore::fill_values(double val, size_t i, size_t j, size_t k, size_t l) {
    size_t ijkl = INDEX4(i,j,k,l);
    size_t ikjl = INDEX4(i, k, j, l);

    if (ijkl >= offset_ && ijkl <= max_idx) {
        J_bufp_[ijkl - offset_] += val;
    }
    if(ikjl >= offset_ && ikjl <= max_idx) {
        if (i == k || j == l) {
            K_bufp_[ikjl - offset_] += val;
        } else {
            K_bufp_[ikjl - offset_] += 0.5 * val;
        }
    }

    if(i != j && k != l) {
        size_t iljk = INDEX4(i, l, j, k);
        if (iljk >= offset_ && iljk <= max_idx) {
            if ( i == l || j == k) {
                K_bufp_[iljk - offset_] += val;
            } else {
                K_bufp_[iljk - offset_] += 0.5 * val;
            }
        }
    }

}

void PKWrkrInCore::write(std::vector<size_t> min_ind, std::vector<size_t> max_ind,
                         size_t pk_pairs) {

    for (size_t pq = 0; pq < pk_pairs; ++pq) {
        size_t pqpq = INDEX2(pq,pq);
        if(pqpq >= offset_ && pqpq <= max_idx) {
            J_bufp_[pqpq - offset_] *= 0.5;
            K_bufp_[pqpq - offset_] *= 0.5;
        }
    }

}

PKWrkrIWL::PKWrkrIWL(boost::shared_ptr<BasisSet> primary, boost::shared_ptr<AIOHandler> AIO,
                     int target_file, int K_file, size_t buf_size, std::vector<int> &bufforpq) :
    PKWorker(primary,AIO,target_file,buf_size) {
    K_file_ = K_file;
    buf_for_pq_ = bufforpq;
    size_t lastpq = buf_for_pq_.size() - 1;
    set_nbuf(buf_for_pq_[lastpq]);

    // Constructing the IWL buffers needed
    for(int i = 0; i < nbuf_; ++i) {
        IWL_J_.push_back(new IWLAsync_PK(&addresses_[2 * i], AIO_, target_file_));
        IWL_K_.push_back(new IWLAsync_PK(&addresses_[2 * i + 1], AIO_, target_file_));
    }
}

PKWrkrIWL::~PKWorker() {
    for(int i = 0; i < nbuf_; ++i) {
        delete IWL_J_[i];
        delete IWL_K_[i];
    }
}

void PKWrkrIWL::fill_values(double val, size_t i, size_t j, size_t k, size_t l) {
    // Pre-sorting for J
    size_t pq = INDEX2(i,j);
    IWLAsync_PK* buf = bufs_J_[buf_for_pq_[pq]];
    buf->fill_values(val,i,j,k,l);
    if(buf->nints() == buf->maxints()) {
        buf->write();
    }

    // Pre-sorting for K
    pq = INDEX2(i,k);
    int bufK1 = buf_for_pq_[pq];
    buf = bufs_K_[bufK1];
    buf->fill_values(val,i,j,k,l);
    if(buf->nints() == buf->maxints()) {
        buf->write();
    }
    // Second pre-sorting for K
    if (i != j && k != l) {
        pq = std::max(INDEX2(i,l), INDEX2(j,k));
        int bufK2 = buf_for_pq_[pq];
        if (bufK2 != bufK1) {
            buf = bufs_K_[bufK2];
            buf->fill_values(val,i,j,k,l);
            if(buf->nints() == buf->maxints()) {
                buf->write();
            }
        }
    }
}

bool PKWrkrIWL::pop_value(unsigned int bufid, double &val, size_t &i, size_t &j, size_t &k, size_t &l) {
    IWLAsync_PK* buf;
    if(bufid < nbuf()) {
        buf = bufs_J_[bufid];
    } else {
        buf = bufs_K_[bufid - nbuf()];
    }
    if(buf->nints() == 0) {
        return false;
    }
    buf->pop_value(val,i,j,k,l);
    return true;
}

void PKWrkrIWL::insert_value(unsigned int bufid, double val, size_t i, size_t j, size_t k, size_t l) {
    IWLAsync_PK* buf;
    if(bufid < nbuf()) {
        buf = bufs_J_[bufid];
    } else {
        buf = bufs_K_[bufid - nbuf()];
    }
    buf->fill_values(val,i,j,k,l);

}

void PKWrkrIWL::flush() {
    IWLAsync_PK* buf;
    for(int bufid = 0; bufid < nbuf(); ++bufid) {
        buf = bufs_J_[bufid];
        buf->flush();
        buf = bufs_K_[bufid];
        buf->flush();
    }
}

IWLAsync_PK::IWLAsync_PK(size_t *address, boost::shared_ptr<AIOHandler> AIO, int itap) {
    itap_ = itap;
    address_ = address;
    AIO_ = AIO;
    ints_per_buf_ = IWL_INTS_PER_BUF;
    nints_ = 0;
    idx_ = 0;
    labels_[0] = new Label[4 * ints_per_buf_];
    labels_[1] = new Label[4 * ints_per_buf_];
    values_[0] = new Value[ints_per_buf_];
    values_[1] = new Value[ints_per_buf_];
    JobID_[0] = 0;
    JobID_[1] = 0;
    lastbuf_ = 0;

}

IWLAsync_PK::~IWLAsync_PK() {
    delete [] labels_[0];
    delete [] labels_[1];
    delete [] values_[0];
    delete [] values_[1];
}

void IWLAsync_PK::fill_values(double val, size_t i, size_t j, size_t k, size_t l) {
    size_t id = 4 * nints_;
    labels_[idx_][id++] = i;
    labels_[idx_][id++] = j;
    labels_[idx_][id++] = k;
    labels_[idx_][id] = l;
    values_[idx_][nints_] = val;
    ++nints_;
}

/// Buffer is full, write it using AIO to disk
void IWLAsync_PK::write() {
    size_t lab_size = 4 * ints_per_buf_ * sizeof(Label);
    size_t val_size = ints_per_buf_ * sizeof(Value);
    // We need a special function in AIO to take care
    // of IWL buffer writing, since these contain four parts.
    JobID_[idx_] = AIO_->write_iwl(itap_,IWL_KEY_BUF,nints_,lastbuf_,(char *) labels_[idx_],
                    (char*) values_[idx_], lab_size, val_size, address_);

    // Now we need to switch the internal buffer to which we are writing.
    idx_ = idx_ == 0 ? 1 : 0;
    nints_ = 0;
    AIO_->wait_for_job(JobID_[idx_]);
}

// Pop a value from the buffer storage
void IWLAsync_PK::pop_value(double &val, size_t &i, size_t &j, size_t &k, size_t &l) {
    if(nints_ == 0) {
        throw PSIEXCEPTION("Cannot pop value from empty buffer\n");
    }
    --nints_;
    size_t id = 4 * nints_;
    i = labels_[idx_][id++];
    j = labels_[idx_][id++];
    k = labels_[idx_][id++];
    l = labels_[idx_][id];
    val = values_[idx_][nints_];
}

void IWLAsync_PK::flush() {
    unsigned int nints = nints_;
    while(nints_ != ints_per_buf_) {
        fill_values(0.0,0,0,0,0);
    }
    nints_ = nints;
    lastbuf_ = 1;
    write();
}

}  // End namespace pk
}  // End namespace psi
