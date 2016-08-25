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

#include "psi4/libmints/integral.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/gshell.h"
#include "psi4/libmints/sieve.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsio/aiohandler.h"
#include "psi4/libiwl/config.h"
#include "PK_workers.h"

namespace psi {

namespace pk {

AOShellSieveIterator::AOShellSieveIterator(std::shared_ptr<BasisSet> prim,
                                           SharedSieve sieve_input) :
shell_pairs_(sieve_input->shell_pairs()) {
    bs_ = prim;
    sieve_ = sieve_input;
    npairs_ = shell_pairs_.size();
    PQ_ = 0;
    RS_ = 0;
    done_ = false;
}

void AOShellSieveIterator::populate_indices() {
    P_ = shell_pairs_[PQ_].first;
    Q_ = shell_pairs_[PQ_].second;
    R_ = shell_pairs_[RS_].first;
    S_ = shell_pairs_[RS_].second;
}

void AOShellSieveIterator::first() {

    PQ_ = 0;
    RS_ = 0;
    populate_indices();
    while(!sieve_->shell_significant(P_,Q_,R_,S_)) {
        // We do not use a function to increment so we can directly
        // return when needed
        ++RS_;
        if(RS_ > PQ_) {
            RS_ = 0;
            ++PQ_;
            if(PQ_ >= npairs_) {
                done_ = true;
                return;
            }
        }
        populate_indices();
    }

}

void AOShellSieveIterator::next() {
    ++RS_;
    if(RS_ > PQ_) {
        RS_ = 0;
        ++PQ_;
        if(PQ_ >= npairs_) {
            done_ = true;
            return;
        }
    }
    populate_indices();
    while(!sieve_->shell_significant(P_,Q_,R_,S_)) {
        // We do not use a function to increment so we can directly
        // return when needed
        ++RS_;
        if(RS_ > PQ_) {
            RS_ = 0;
            ++PQ_;
            if(PQ_ >= npairs_) {
                done_ = true;
                return;
            }
        }
        populate_indices();
    }
}

AOFctSieveIterator AOShellSieveIterator::integrals_iterator() {
    return AOFctSieveIterator(bs_->shell(P_),bs_->shell(Q_),
                              bs_->shell(R_),bs_->shell(S_),sieve_);

}

AOFctSieveIterator::AOFctSieveIterator(const GaussianShell &s1, const GaussianShell &s2,
                    const GaussianShell &s3, const GaussianShell &s4,
                    std::shared_ptr<ERISieve> siev) : usi_(s1), usj_(s2), usk_(s3), usl_(s4) {
    sieve_ = siev;
    done_ = false;

    ni_ = usi_.nfunction();
    nj_ = usj_.nfunction();
    nk_ = usk_.nfunction();
    nl_ = usl_.nfunction();

    fi_ = usi_.function_index();
    fj_ = usj_.function_index();
    fk_ = usk_.function_index();
    fl_ = usl_.function_index();

    sh_aaaa_ = (&usi_ == &usj_ && &usk_ == &usl_ && &usi_ == &usk_);
    sh_abab_ = (&usi_ == &usk_ && &usj_ == &usl_);
    maxi_ = ni_ - 1;
}

void AOFctSieveIterator::populate_indices() {
    i_ = irel_ + fi_;
    j_ = jrel_ + fj_;
    k_ = krel_ + fk_;
    l_ = lrel_ + fl_;
}

void AOFctSieveIterator::increment_bra() {
    lrel_ = 0;
    krel_ = 0;
    ++jrel_;
    if(jrel_ > maxj_) {
        jrel_ = 0;
        ++irel_;
        if(irel_ > maxi_) {
            done_ = true;
        }
        if(sh_aaaa_) {
            maxj_ = irel_;
        } else if(!sh_abab_) {
            maxj_ = (&usi_ == &usj_) ? irel_ : nj_ - 1;
        }
    }

    if(sh_aaaa_) {
        maxk_ = irel_;
        maxl_ = (krel_ == irel_) ? jrel_ : krel_;
    } else if(sh_abab_) {
        maxk_ = irel_;
        maxl_ = (krel_ == irel_) ? jrel_ : nl_ - 1;
    } else {
        maxl_ = (&usk_ == &usl_) ? krel_ : nl_ - 1;
    }
    populate_indices();
}

void AOFctSieveIterator::increment_ket() {
    // Here we adopt the original implementation, might be more efficient
    if (sh_aaaa_) {
        ++lrel_;
        if(lrel_ > maxl_){
            ++krel_;
            lrel_ = 0;
            if(krel_ > maxk_){
                krel_ = 0;
                ++jrel_;
                if(jrel_ > maxj_){
                    jrel_ = 0;
                    ++irel_;
                    if(irel_ > maxi_){
                        done_ = true;
                    }
                    maxj_ = irel_;
                }
                maxk_ = irel_;

            }
            maxl_ = (krel_==irel_) ? jrel_ : krel_;
        }
    } else if(sh_abab_){
        ++lrel_;
        if(lrel_ > maxl_){
            ++krel_;
            lrel_ = 0;
            if(krel_ > maxk_){
                krel_ = 0;
                ++jrel_;
                if(jrel_ > maxj_){
                    jrel_ = 0;
                    ++irel_;
                    if(irel_ > maxi_){
                        done_ = true;
                    }
                }
                maxk_ = irel_;
            }
            maxl_ = (krel_ == irel_) ? jrel_ : nl_ - 1;
        }
    } else {
        ++lrel_;
        if(lrel_ > maxl_){
            ++krel_;
            lrel_ = 0;
            if(krel_ > maxk_){
                krel_ = 0;
                ++jrel_;
                if(jrel_ > maxj_){
                    jrel_ = 0;
                    ++irel_;
                    if(irel_ > maxi_){
                        done_ = true;
                    }
                    maxj_ = (&usi_ == &usj_) ? irel_ : nj_ - 1;
                }
            }
            maxl_ = (&usk_ == &usl_) ? krel_ : nl_ - 1;
        }
    }
    populate_indices();

}

void AOFctSieveIterator::reorder_inds() {
    if(!sh_aaaa_) {
        if(sh_abab_) {
            if(i_ < j_) {
                std::swap(i_,j_);
                std::swap(k_,l_);
            }
            if(i_ < k_) {
                std::swap(i_,k_);
                std::swap(j_,l_);
            }
        } else {
            if(i_ < j_) {
                std::swap(i_,j_);
            }
            if(k_ < l_) {
                std::swap(k_,l_);
            }
            if((i_ < k_) || (i_ == k_ && j_ < l_)) {
                std::swap(i_,k_);
                std::swap(j_,l_);

            }
        }
    }

}

// Could be more efficient to include increment steps
// explicitly to be able to return directly and avoid
// the if(done_) conditional (much more ugly though)
void AOFctSieveIterator::first() {
    if (sh_aaaa_) {
        maxk_ = 0;
        maxl_ = 0;
        maxj_ = 0;
    } else if (sh_abab_) {
        maxk_ = 0;
        maxl_ = 0;
        maxj_ = nj_ - 1;
    } else {
        maxk_ = nk_ - 1;
        maxj_ = (&usi_ == &usj_) ? 0 : nj_ - 1;
        maxl_ = (&usk_ == &usl_) ? 0 : nl_ - 1;
    }

    irel_ = 0;
    jrel_ = 0;
    krel_ = 0;
    lrel_ = 0;
    populate_indices();
    // find a significant ij
    while(!sieve_->function_pair_significant(i_,j_)) {
        increment_bra();
        if(done_) return;
    }
    // find a significant integral
    while(!sieve_->function_significant(i_,j_,k_,l_)) {
        increment_ket();
        if(done_) return;
        // If ij changes, find next significant pair
        while(!sieve_->function_pair_significant(i_,j_)) {
            increment_bra();
            if(done_) return;
        }
    }

    reorder_inds();
}

void AOFctSieveIterator::next() {
    increment_ket();
    if(done_) return;
    while(!sieve_->function_pair_significant(i_,j_)) {
        increment_bra();
        if(done_) return;
    }
    // find a significant integral
    while(!sieve_->function_significant(i_,j_,k_,l_)) {
        increment_ket();
        if(done_) return;
        // If ij changes, find next significant pair
        while(!sieve_->function_pair_significant(i_,j_)) {
            increment_bra();
            if(done_) return;
        }
    }

    reorder_inds();
}



PKWorker::PKWorker(std::shared_ptr<BasisSet> primary, SharedSieve sieve,
                   std::shared_ptr<AIOHandler> AIO, int target_file,
                   size_t buf_size) {
    AIO_ = AIO;
    sieve_ = sieve;
    target_file_ = target_file;
    primary_ = primary;
    buf_size_ = buf_size;
    bufidx_ = 0;
    offset_ = 0;
    do_wK_ = false;

}

char* PKWorker::get_label_J(const int batch) {
    char* label = new char[100];
    sprintf(label, "J Block (Batch %d)", batch);
    return label;
}

char* PKWorker::get_label_K(const int batch) {
    char* label = new char[100];
    sprintf(label, "K Block (Batch %d)", batch);
    return label;
}

char* PKWorker::get_label_wK(const int batch) {
    char* label = new char[100];
    sprintf(label, "wK Block (Batch %d)", batch);
    return label;
}

void PKWorker::first_quartet(size_t i) {
    shelliter_ = UniqueAOShellIt(new AOShellSieveIterator(primary_,sieve_)) ;
    bufidx_ = i;
    offset_ = bufidx_ * buf_size_;
    initialize_task();
//DEBUG    #pragma omp critical
//DEBUG    outfile->Printf("thread %d, offset is %lu and max_idx is %lu\n",omp_get_thread_num(),offset_,max_idx_);
//DEBUG    std::cout << "thread" << omp_get_thread_num() << ", offset is " << offset_ << " and max_idx is " << max_idx_ << std::endl;
    shells_left_ = false;
    for(shelliter_->first(); (shells_left_ || shelliter_->is_done()) == false; shelliter_->next()) {
        P_ = shelliter_->p();
        Q_ = shelliter_->q();
        R_ = shelliter_->r();
        S_ = shelliter_->s();

        shells_left_ = is_shell_relevant();
    }

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
//DEBUG        outfile->Printf("Rejecting1 shell <%d %d|%d %d>\n",P_,Q_,R_,S_);
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
//DEBUG        outfile->Printf("Rejecting2 shell <%d %d|%d %d>\n",P_,Q_,R_,S_);
        return false;
    }

    // Now we loop over unique basis function quartets in the shell quartet
    AOFctSieveIterator bfiter = shelliter_->integrals_iterator();
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
//DEBUG        outfile->Printf("Accepting shell <%d %d|%d %d>\n",P_,Q_,R_,S_);
            return true;
        }
    }

//DEBUG    outfile->Printf("Rejecting shell3 <%d %d|%d %d>\n",P_,Q_,R_,S_);
    return false;

}

void PKWorker::next_quartet() {
    if(shelliter_->is_done()) {
        shells_left_ = false;
        return;
    }
    bool shell_found = false;
    for(; (shell_found || shelliter_->is_done()) == false; shelliter_->next()) {
        P_ = shelliter_->p();
        Q_ = shelliter_->q();
        R_ = shelliter_->r();
        S_ = shelliter_->s();

        shell_found = is_shell_relevant();
    }

    shells_left_ = shell_found;
}

PKWrkrReord::PKWrkrReord(std::shared_ptr<BasisSet> primary, SharedSieve sieve, std::shared_ptr<AIOHandler> AIO,
                         int target_file, size_t buffer_size, unsigned int nbuffer) :
    PKWorker(primary,sieve,AIO,target_file,buffer_size) {

    set_nbuf(nbuffer);
    buf_ = 0;

    double * mem;
    // Allocate everything we need
    for(int i = 0; i < nbuf(); ++i) {
        mem = new double[buf_size()];
        J_bufs_.push_back(mem);
        mem = new double[buf_size()];
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
    ::memset((void*) J_bufs_[0], '\0', buf_size() * sizeof(double));
    ::memset((void*) K_bufs_[0], '\0', buf_size() * sizeof(double));
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
    for(it = wK_bufs_.begin(); it != wK_bufs_.end(); ++it) {
        delete [] *it;
    }
    wK_bufs_.clear();
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
    for(int i = 0; i < labels_wK_.size(); ++i) {
        for(int j = 0; j < labels_wK_[i].size(); ++j) {
            delete [] labels_wK_[i][j];
        }
        labels_wK_[i].clear();
    }
    labels_wK_.clear();

}

void PKWrkrReord::allocate_wK(size_t bufsize, unsigned int buf_per_thread) {
    // We want to deallocate the previous J and K buffers
    // Make sure the corresponding integrals are written to disk
    for(int i = 0; i < nbuf(); ++i) {
        for(int j = 0; j < jobID_J_[i].size(); ++j) {
            AIO()->wait_for_job(jobID_J_[i][j]);
            delete [] labels_J_[i][j];
        }
        labels_J_[i].clear();
        jobID_J_[i].clear();
        delete [] J_bufs_[i];
        for(int j = 0; j < jobID_K_[i].size(); ++j) {
            AIO()->wait_for_job(jobID_K_[i][j]);
            delete [] labels_K_[i][j];
        }
        labels_K_[i].clear();
        jobID_K_[i].clear();
        delete [] K_bufs_[i];
    }
    labels_J_.clear();
    labels_K_.clear();
    jobID_J_.clear();
    jobID_K_.clear();
    K_bufs_.clear();
    J_bufs_.clear();

    // wK gets buffers with updated size
    set_nbuf(buf_per_thread);
    set_bufsize(bufsize);
    buf_ = 0;

    for(int i = 0; i < nbuf(); ++i) {
        wK_bufs_.push_back(new double[buf_size()]);
        std::vector< char* > labelwK;
        labels_wK_.push_back(labelwK);
        std::vector< size_t > idwK;
        jobID_wK_.push_back(idwK);
    }
    // Only need to set to zero the first buffer
    ::memset((void*) wK_bufs_[0],'\0', buf_size() * sizeof(double));
}

void PKWrkrReord::initialize_task() {
    set_max_idx(buf_size() * (bufidx() + 1) - 1);
}

void PKWrkrReord::fill_values(double val, size_t i, size_t j, size_t k, size_t l) {

    size_t ijkl = INDEX4(i,j,k,l);
    if (ijkl >= offset() && ijkl <= max_idx()) {
        J_bufs_[buf_][ijkl - offset()] += val;
    }
    size_t ikjl = INDEX4(i, k, j, l);
    if(ikjl >= offset() && ikjl <= max_idx()) {
        if (i == k || j == l) {
            K_bufs_[buf_][ikjl - offset()] += val;
        } else {
            K_bufs_[buf_][ikjl - offset()] += 0.5 * val;
        }
    }

    if(i != j && k != l) {
        size_t iljk = INDEX4(i, l, j, k);
        if (iljk >= offset() && iljk <= max_idx()) {
            if ( i == l || j == k) {
                K_bufs_[buf_][iljk - offset()] += val;
            } else {
                K_bufs_[buf_][iljk - offset()] += 0.5 * val;
            }
        }
    }
}

void PKWrkrReord::fill_values_wK(double val, size_t i, size_t j, size_t k, size_t l) {
    size_t ijkl = INDEX4(i, j, k, l);
    if(ijkl >= offset() && ijkl <= max_idx()) {
        wK_bufs_[buf_][ijkl - offset()] += val;
    }
}

void PKWrkrReord::write(std::vector<size_t> min_ind, std::vector<size_t> max_ind, size_t pk_pairs) {
    // Compute initial ijkl index for current buffer

    // Vector of batch number to which we want to write
    std::vector<unsigned int> target_batches;

    for (unsigned int i = 0; i < min_ind.size(); ++i) {
        if (offset() >= min_ind[i] && offset() < max_ind[i]) {
            target_batches.push_back(i);
            continue;
        }
        if (max_idx() >= min_ind[i] && max_idx() < max_ind[i]) {
            target_batches.push_back(i);
            continue;
        }
        if (offset() < min_ind[i] && max_idx() >= max_ind[i]) {
            target_batches.push_back(i);
            continue;
        }
    }

    // Now that all buffers are full, and before we write integrals, we need
    // to divide diagonal elements by 2.

    for (size_t pq = 0; pq < pk_pairs; ++pq) {
        size_t pqpq = INDEX2(pq, pq);
        if (pqpq >= offset() && pqpq <= max_idx()) {
            J_bufs_[buf_][pqpq - offset()] *= 0.5;
            K_bufs_[buf_][pqpq - offset()] *= 0.5;
        }
    }

    // And now we write to the file in the appropriate entries
    for (int i = 0; i < target_batches.size(); ++i) {
        unsigned int b = target_batches[i];
        labels_J_[buf_].push_back(get_label_J(b));
        size_t start = std::max(offset(), min_ind[b]);
        size_t stop = std::min(max_idx() + 1, max_ind[b]);
        psio_address adr = psio_get_address(PSIO_ZERO, (start - min_ind[b]) * sizeof(double));
        size_t nints = stop - start;
        jobID_J_[buf_].push_back(AIO()->write(target_file(), labels_J_[buf_][i], (char *)(&J_bufs_[buf_][start - offset()]),
                    nints * sizeof(double), adr, &dummy_));
        labels_K_[buf_].push_back(get_label_K(b));
        jobID_K_[buf_].push_back(AIO()->write(target_file(), labels_K_[buf_][i], (char *)(&K_bufs_[buf_][start - offset()]),
                    nints * sizeof(double), adr, &dummy_));
    }

    // Update the buffer being written into
    ++buf_;
    if (buf_ >= nbuf()) buf_ = 0;
    // Make sure the buffer has been written to disk and we can erase it
    for(int i = 0; i < jobID_J_[buf_].size(); ++i) {
      AIO()->wait_for_job(jobID_J_[buf_][i]);
    }
    jobID_J_[buf_].clear();
    for(int i = 0; i < jobID_K_[buf_].size(); ++i) {
      AIO()->wait_for_job(jobID_K_[buf_][i]);
    }
    jobID_K_[buf_].clear();
    // We can delete the labels for these buffers
    for(int i = 0; i < labels_J_[buf_].size(); ++i) {
        delete [] labels_J_[buf_][i];
    }
    for(int i = 0; i < labels_K_[buf_].size(); ++i) {
        delete [] labels_K_[buf_][i];
    }
    labels_J_[buf_].clear();
    labels_K_[buf_].clear();
    // Make sure the buffer in which we are going to write is set to zero
    ::memset((void *) J_bufs_[buf_], '\0', buf_size() * sizeof(double));
    ::memset((void *) K_bufs_[buf_], '\0', buf_size() * sizeof(double));

}

void PKWrkrReord::write_wK(std::vector<size_t> min_ind, std::vector<size_t> max_ind, size_t pk_pairs) {
    // Compute initial ijkl index for current buffer

    // Vector of batch number to which we want to write
    std::vector<unsigned int> target_batches;

    for (unsigned int i = 0; i < min_ind.size(); ++i) {
        if (offset() >= min_ind[i] && offset() < max_ind[i]) {
            target_batches.push_back(i);
            continue;
        }
        if (max_idx() >= min_ind[i] && max_idx() < max_ind[i]) {
            target_batches.push_back(i);
            continue;
        }
        if (offset() < min_ind[i] && max_idx() >= max_ind[i]) {
            target_batches.push_back(i);
            continue;
        }
    }

    // Now that all buffers are full, and before we write integrals, we need
    // to divide diagonal elements by 2.

    for (size_t pq = 0; pq < pk_pairs; ++pq) {
        size_t pqpq = INDEX2(pq, pq);
        if (pqpq >= offset() && pqpq <= max_idx()) {
            wK_bufs_[buf_][pqpq - offset()] *= 0.5;
        }
    }

    // And now we write to the file in the appropriate entries
    for (int i = 0; i < target_batches.size(); ++i) {
        unsigned int b = target_batches[i];
        labels_wK_[buf_].push_back(get_label_wK(b));
        size_t start = std::max(offset(), min_ind[b]);
        size_t stop = std::min(max_idx() + 1, max_ind[b]);
        psio_address adr = psio_get_address(PSIO_ZERO, (start - min_ind[b]) * sizeof(double));
        size_t nints = stop - start;
        jobID_wK_[buf_].push_back(AIO()->write(target_file(), labels_wK_[buf_][i], (char *)(&wK_bufs_[buf_][start - offset()]),
                    nints * sizeof(double), adr, &dummy_));
    }

    // Update the buffer being written into
    ++buf_;
    if (buf_ >= nbuf()) buf_ = 0;
    // Make sure the next buffer has been written to disk and we can erase it
    for(int i = 0; i < jobID_wK_[buf_].size(); ++i) {
      AIO()->wait_for_job(jobID_wK_[buf_][i]);
    }
    jobID_wK_[buf_].clear();
    // We can delete the labels for these buffers
    for(int i = 0; i < labels_wK_[buf_].size(); ++i) {
        delete [] labels_wK_[buf_][i];
    }
    labels_wK_[buf_].clear();
    // Make sure the buffer in which we are going to write is set to zero
    ::memset((void *) wK_bufs_[buf_], '\0', buf_size() * sizeof(double));

}

PKWrkrInCore::PKWrkrInCore(std::shared_ptr<BasisSet> primary, SharedSieve sieve, size_t buf_size,
                           size_t lastbuf, double *Jbuf, double *Kbuf, double *wKbuf, int nworkers) :
    PKWorker(primary,sieve,std::shared_ptr<AIOHandler>(),0,buf_size) {

    nworkers_ = nworkers;
    last_buf_ = lastbuf;
    J_buf0_ = Jbuf;
    K_buf0_ = Kbuf;
    // wKbuf is NULL if we don't compute wK
    wK_buf0_ = wKbuf;

    J_bufp_ = NULL;
    K_bufp_ = NULL;
    wK_bufp_ = NULL;
}

void PKWrkrInCore::initialize_task() {
    size_t maxid = buf_size() * (bufidx() + 1);
    // If we are at the last worker, we extend the buffer to
    // include the last integrals
    if (bufidx() == nworkers_ - 1) {
        maxid += last_buf_;
    }
    set_max_idx(maxid - 1);
    //We set the pointers to the beginning of the attributed buffer section
    if(do_wK()) {
        wK_bufp_ = wK_buf0_ + offset();
    } else {
        J_bufp_ = J_buf0_ + offset();
        K_bufp_ = K_buf0_ + offset();
    }
}

void PKWrkrInCore::fill_values(double val, size_t i, size_t j, size_t k, size_t l) {
    size_t ijkl = INDEX4(i,j,k,l);
    size_t ikjl = INDEX4(i, k, j, l);

//DEBUG    size_t pqd = INDEX2(i,j);
//DEBUG    size_t rsd = INDEX2(k,l);
//DEBUG    if(rsd == 0) {
//DEBUG#pragma omp critical
//DEBUG      outfile->Printf("PK int (%lu|%lu) = %20.16f\n",pqd,rsd,val);
//DEBUG    }
    if (ijkl >= offset() && ijkl <= max_idx()) {
        J_bufp_[ijkl - offset()] += val;
//DEBUG    size_t pqd = INDEX2(i,j);
//DEBUG    size_t rsd = INDEX2(k,l);
//DEBUG    if(rsd == 0) {
//DEBUG#pragma omp critical
//DEBUG{
//DEBUG      outfile->Printf("For thread %d, ijkl = %lu, offset = %lu, max_idx = %lu\n", omp_get_thread_num(), ijkl, offset(), max_idx());
//DEBUG      outfile->Printf("PK int (%lu|%lu) = %20.16f\n",pqd,rsd,J_bufp_[ijkl-offset()]);
//DEBUG}
//DEBUG    }
//DEBUG#pragma omp critical
//DEBUG        outfile->Printf("Value at 0 for J: %16.12f\n",J_bufp_[0]);
    }
    if(ikjl >= offset() && ikjl <= max_idx()) {
        if (i == k || j == l) {
            K_bufp_[ikjl - offset()] += val;
        } else {
            K_bufp_[ikjl - offset()] += 0.5 * val;
        }
    }

    if(i != j && k != l) {
        size_t iljk = INDEX4(i, l, j, k);
        if (iljk >= offset() && iljk <= max_idx()) {
            if ( i == l || j == k) {
                K_bufp_[iljk - offset()] += val;
            } else {
                K_bufp_[iljk - offset()] += 0.5 * val;
            }
        }
    }

}

void PKWrkrInCore::fill_values_wK(double val, size_t i, size_t j, size_t k, size_t l) {
    size_t ijkl = INDEX4(i, j, k, l);

    if(ijkl >= offset() && ijkl <= max_idx()) {
//DEBUG        size_t pqd = INDEX2(i,j);
//DEBUG        size_t rsd = INDEX2(k,l);
//DEBUG        if(rsd == 0) {
//DEBUG#pragma omp critical
//DEBUG          outfile->Printf("PK int (%lu|%lu) = %20.16f from thread %d\n",pqd,rsd,val,omp_get_thread_num());
//DEBUG        }
        wK_bufp_[ijkl - offset()] += val;
    }

}

void PKWrkrInCore::finalize_ints(size_t pk_pairs) {

    for (size_t pq = 0; pq < pk_pairs; ++pq) {
        size_t pqpq = INDEX2(pq,pq);
        if(pqpq >= offset() && pqpq <= max_idx()) {
            J_bufp_[pqpq - offset()] *= 0.5;
            K_bufp_[pqpq - offset()] *= 0.5;
        }
    }

}

void PKWrkrInCore::finalize_ints_wK(size_t pk_pairs) {

    for (size_t pq = 0; pq < pk_pairs; ++pq) {
        size_t pqpq = INDEX2(pq,pq);
        if(pqpq >= offset() && pqpq <= max_idx()) {
            wK_bufp_[pqpq - offset()] *= 0.5;
        }
    }
}

PKWrkrIWL::PKWrkrIWL(std::shared_ptr<BasisSet> primary, SharedSieve sieve, std::shared_ptr<AIOHandler> AIOp,
                     int targetfile, int K_file, size_t buf_size, std::vector<int> &bufforpq,
                     std::shared_ptr<std::vector<size_t>> pos) :
    PKWorker(primary,sieve,AIOp,targetfile,buf_size) {
    K_file_ = K_file;
    buf_for_pq_ = bufforpq;
    size_t lastpq = buf_for_pq_.size() - 1;
    set_nbuf(buf_for_pq_[lastpq] + 1);
    addresses_ = pos;

    // Constructing the IWL buffers needed
    for(int i = 0; i < nbuf(); ++i) {
        IWL_J_.push_back(new IWLAsync_PK(&((*addresses_)[2 * i]), AIO(), target_file()));
        IWL_K_.push_back(new IWLAsync_PK(&((*addresses_)[2 * i + 1]), AIO(), K_file_));
    }
}

PKWrkrIWL::~PKWrkrIWL() {
    for(int i = 0; i < nbuf(); ++i) {
        delete IWL_J_[i];
        delete IWL_K_[i];
    }
    for(int i = 0; i < IWL_wK_.size(); ++i) {
        delete IWL_wK_[i];
    }
}

void PKWrkrIWL::allocate_wK(std::shared_ptr<std::vector<size_t>> pos, int wKfile) {
    wK_file_ = wKfile;
    addresses_wK_ = pos;

    // Constructing the IWL buffers needed for wK
    for(int i = 0; i < nbuf(); ++i) {
        IWL_wK_.push_back(new IWLAsync_PK(&((*addresses_wK_)[i]), AIO(), wK_file_));
    }
}

void PKWrkrIWL::fill_values(double val, size_t i, size_t j, size_t k, size_t l) {
    // Pre-sorting for J
    size_t pq = INDEX2(i,j);
    IWLAsync_PK* buf = IWL_J_[buf_for_pq_[pq]];
//DEBUG    outfile->Printf("Adding value to J\n");
//DEBUG    if(INDEX2(k,l) == 0) {
//DEBUG      outfile->Printf("Buffering int <%d|%d> = %16.12f\n",pq,INDEX2(k,l),val);
//DEBUG    }
    buf->fill_values(val,i,j,k,l);
    if(buf->nints() == buf->maxints()) {
        buf->write();
    }

    // Pre-sorting for K
    pq = INDEX2(i,k);
    int bufK1 = buf_for_pq_[pq];
    buf = IWL_K_[bufK1];
//DEBUG    outfile->Printf("Adding value to K1\n");
    buf->fill_values(val,i,j,k,l);
    if(buf->nints() == buf->maxints()) {
        buf->write();
    }
    // Second pre-sorting for K
    if (i != j && k != l) {
        pq = std::max(INDEX2(i,l), INDEX2(j,k));
        int bufK2 = buf_for_pq_[pq];
        if (bufK2 != bufK1) {
            buf = IWL_K_[bufK2];
//DEBUG    outfile->Printf("Adding value to K2\n");
            buf->fill_values(val,i,j,k,l);
            if(buf->nints() == buf->maxints()) {
                buf->write();
            }
        }
    }
}

void PKWrkrIWL::fill_values_wK(double val, size_t i, size_t j, size_t k, size_t l) {
    size_t pq = INDEX2(i,j);
    IWLAsync_PK* buf = IWL_wK_[buf_for_pq_[pq]];
    buf->fill_values(val,i,j,k,l);
    if(buf->nints() == buf->maxints()) {
        buf->write();
    }
}

bool PKWrkrIWL::pop_value(unsigned int bufid, double &val, size_t &i, size_t &j, size_t &k, size_t &l) {
    IWLAsync_PK* buf;
    if(bufid < nbuf()) {
        buf = IWL_J_[bufid];
    } else {
        buf = IWL_K_[bufid - nbuf()];
    }
    if(buf->nints() == 0) {
        return false;
    }
    buf->pop_value(val,i,j,k,l);
    return true;
}

bool PKWrkrIWL::pop_value_wK(unsigned int bufid, double &val, size_t &i, size_t &j, size_t &k, size_t &l) {
    IWLAsync_PK* buf = IWL_wK_[bufid];
    if(buf->nints() == 0) {
        return false;
    }
    buf->pop_value(val,i,j,k,l);
    return true;
}

void PKWrkrIWL::insert_value(unsigned int bufid, double val, size_t i, size_t j, size_t k, size_t l) {
    IWLAsync_PK* buf;
    if(bufid < nbuf()) {
        buf = IWL_J_[bufid];
    } else {
        buf = IWL_K_[bufid - nbuf()];
    }
    buf->fill_values(val,i,j,k,l);
    if(buf->nints() == buf->maxints()) {
        buf->write();
    }
}

void PKWrkrIWL::insert_value_wK(unsigned int bufid, double val, size_t i, size_t j, size_t k, size_t l) {
    IWLAsync_PK* buf = IWL_wK_[bufid];
    buf->fill_values(val,i,j,k,l);
    if(buf->nints() == buf->maxints()) {
        buf->write();
    }
}

void PKWrkrIWL::flush() {
    IWLAsync_PK* buf;
    for(int bufid = 0; bufid < nbuf(); ++bufid) {
        buf = IWL_J_[bufid];
        buf->flush();
        buf = IWL_K_[bufid];
        buf->flush();
    }
}

void PKWrkrIWL::flush_wK() {
    IWLAsync_PK* buf;
    for(int bufid = 0; bufid < nbuf(); ++bufid) {
        buf = IWL_wK_[bufid];
        buf->flush();
    }
}

IWLAsync_PK::IWLAsync_PK(size_t *address, std::shared_ptr<AIOHandler> AIO, int itap) {
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
    values_[idx_][nints_++] = val;
}

/// Buffer is full, write it using AIO to disk
void IWLAsync_PK::write() {
    size_t lab_size = 4 * ints_per_buf_ * sizeof(Label);
    size_t val_size = ints_per_buf_ * sizeof(Value);
    // We need a special function in AIO to take care
    // of IWL buffer writing, since these contain four parts.
//DEBUG    outfile->Printf("Writing to %d.\n",itap_);
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
    while(nints_ < ints_per_buf_) {
        fill_values(0.0,0,0,0,0);
    }
    nints_ = nints;
    lastbuf_ = 1;
    write();
}

}  // End namespace pk
}  // End namespace psi
