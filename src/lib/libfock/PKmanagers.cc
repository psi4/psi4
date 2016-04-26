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

#include <psi4-dec.h>
#include "psifiles.h"
#include <libpsio/psio.hpp>
#include <libiwl/iwl.hpp>
#include "PKmanagers.h"
#include "PK_workers.h"
#include <liboptions/liboptions.h>
#include <libmints/integral.h>
#include <libmints/twobody.h>
#include <libpsio/aiohandler.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace psi {

namespace pk {

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


std::shared_ptr<PKManager> PKManager::build_PKManager(boost::shared_ptr<PSIO> psio,
                  boost::shared_ptr<BasisSet> primary, size_t memory, Options &options) {

    std::string algo = options.get_str("PK_ALGO");
    bool noincore = options.get_bool("PK_NO_INCORE");

    size_t nbf = primary->nbf();
    size_t pk_size = nbf * (nbf + 1) / 2;
    pk_size = pk_size * (pk_size + 1) / 2;

    if(2 * pk_size < memory && !noincore) {
        PKMgrInCore* pkmgr = new PKMgrInCore(primary,memory,options);
        return std::shared_ptr<PKManager>(pkmgr);
    } else if(algo == "REORDER") {
        PKMgrReorder* pkmgr = new PKMgrReorder(psio,primary,memory,options);
        return std::shared_ptr<PKManager>(pkmgr);
    } else if (algo == "YOSHIMINE") {
        PKMgrYoshimine* pkmgr = new PKMgrYoshimine(psio,primary,memory,options);
        return std::shared_ptr<PKManager>(pkmgr);
    }

}

PKManager::PKManager(boost::shared_ptr<BasisSet> primary, size_t memory, Options &options) {

    options_ = options;

    primary_ = primary;
    nbf_ = primary_->nbf();
    memory_ = memory;

    pk_pairs_ = (size_t) nbf_ * ((size_t) nbf_ + 1) / 2;
    pk_size_ = pk_pairs_ * (pk_pairs_ + 1) / 2;
    cutoff_ = options.get_double("INTS_TOLERANCE");
    ntasks_ = 0;

    if(memory_ < pk_pairs_) {
        throw PSIEXCEPTION("Not enough memory for PK algorithm\n");
    }

    // Get number of threads
    nthreads_ = 1;
#ifdef _OPENMP
    nthreads_ = omp_get_max_threads();
#endif
}

SharedPKWrkr PKManager::get_buffer() {
    int thread = 0;
#ifdef _OPENMP
    thread = omp_get_thread_num();
#endif
    return iobuffers_[thread];
}

void PKManager::print_batches() {
    outfile->Printf( "   Calculation information:\n");
    outfile->Printf( "      Number of atoms:                %4d\n", primary_->molecule()->natom());
    outfile->Printf( "      Number of AO shells:            %4d\n", primary_->nshell());
    outfile->Printf( "      Number of primitives:           %4d\n", primary_->nprimitive());
    outfile->Printf( "      Number of atomic orbitals:      %4d\n", primary_->nao());
    outfile->Printf( "      Number of basis functions:      %4d\n\n", primary_->nbf());
    outfile->Printf( "      Integral cutoff                 %4.2e\n", cutoff_);
    outfile->Printf( "      Number of threads:              %4d\n", nthreads_);
    outfile->Printf( "\n");
}


// Integral computation based on several tasks
// computing each an interval of canonical indices ijkl
void PKManager::compute_integrals() {

    // Get an AO integral factory
    boost::shared_ptr<IntegralFactory> intfact(new IntegralFactory(primary_));

    // Get ERI object, one per thread
    std::vector<boost::shared_ptr<TwoBodyAOInt> > tb;
    for(int i = 0; i < nthreads_; ++i) {
        tb.push_back(boost::shared_ptr<TwoBodyAOInt>(intfact->erd_eri()));
    }


    size_t nshqu = 0;
#pragma omp parallel for schedule(dynamic) reduction(+:nshqu)
    for(size_t i = 0; i < ntasks_; ++i) {
        // We need to get the list of shell quartets for each task
        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif
        SharedPKWrkr buf = get_buffer();
        for(buf->first_quartet(i); buf->more_work(); buf->next_quartet()) {
            unsigned int P = buf->P();
            unsigned int Q = buf->Q();
            unsigned int R = buf->R();
            unsigned int S = buf->S();
            tb[thread]->compute_shell(P,Q,R,S);
            integrals_buffering(tb[thread]->buffer(),P,Q,R,S);
            ++nshqu;
        }
        write();
    }
    outfile->Printf("  We computed %lu shell quartets total.\n",nshqu);
    size_t nsh = primary_->nshell();
    size_t nsh_u = nsh * (nsh + 1) / 2;
    nsh_u = nsh_u * (nsh_u + 1) / 2;
    outfile->Printf("  Whereas there are %lu unique shell quartets.\n",nsh_u);
    outfile->Printf("  %7.2f percent of shell quartets recomputed.\n", (nshqu - nsh_u) / float(nsh_u) * 100);

}

void PKManager::integrals_buffering(const double *buffer, unsigned int P, unsigned int Q,
                                    unsigned int R, unsigned int S) {
    int thread = 0;
#ifdef _OPENMP
    thread = omp_get_thread_num();
#endif
    AOIntegralsIterator bfiter(primary_->shell(P), primary_->shell(Q), primary_->shell(R), primary_->shell(S));

    for (bfiter.first(); bfiter.is_done() == false; bfiter.next()) {
        int i = bfiter.i();
        int j = bfiter.j();
        int k = bfiter.k();
        int l = bfiter.l();
        size_t idx = bfiter.index();

        double val = buffer[idx];
        if(val > cutoff_) {
            iobuffers_[thread]->fill_values(val, i, j, k, l);
        }
    }

}

void PKManager::write() {
    int thread = 0;
#ifdef _OPENMP
    thread = omp_get_thread_num();
#endif
    iobuffers_[thread]->write(batch_index_min_,batch_index_max_,pk_pairs_);
}

void PKManager::form_D_vec(std::vector<SharedMatrix> D) {

    // Assume symmetric density matrix for now.
    for (int N = 0; N < D.size(); ++N) {
        double* Dvec = new double[pk_pairs_];
        ::memset((void *)Dvec,'\0',pk_pairs_ * sizeof(double));
        D_vec_.push_back(Dvec);
        size_t pqval = 0;
        for(int p = 0; p < nbf_; ++p) {
            for(int q = 0; q <= p; ++q) {
                if(p != q) {
                    Dvec[pqval] = 2.0 * D[N]->get(0,p,q);
                } else {
                    Dvec[pqval] = D[N]->get(0,p,q);
                }
                ++pqval;
            }
        }
    }

}

void PKManager::make_J_vec(std::vector<SharedMatrix> J) {
    for(int N = 0; N < J.size(); ++N) {
        double* Jvec = new double[pk_pairs_];
        ::memset((void*)Jvec,'\0',pk_pairs_ * sizeof(double));
        JK_vec_.push_back(Jvec);
    }
}

void PKManager::get_results(std::vector<SharedMatrix> J) {
    // Now, directly transfer data to resulting matrices
    for(int N = 0; N < J.size(); ++N) {
        double *Jp = JK_vec_[N];
        for(int p = 0; p < nbf_; ++p) {
            for(int q = 0; q <= p; ++q) {
                J[N]->set(0,p,q,*Jp++);
            }
        }
        J[N]->copy_lower_to_upper();
        delete [] JK_vec_[N];
    }

}

void PKManager::form_K(std::vector<SharedMatrix> K) {
    // For symmetric densities, we do exactly the same than for J
    // but we read another PK entry
    form_J(K,true);
}

void PKManager::finalize_D() {
    for(int N = 0; N < D_vec_.size(); ++N) {
        delete [] D_vec_[N];
    }
    D_vec_.clear();
}

char* PKManager::get_label_J(const int batch) {
    char* label = new char[100];
    sprintf(label, "J Block (Batch %d)", i);
    return label;
}

char* PKManager::get_label_K(const int batch) {
    char* label = new char[100];
    sprintf(label, "K Block (Batch %d)", i);
    return label;
}

PKMgrDisk::PKMgrDisk(boost::shared_ptr<PSIO> psio, boost::shared_ptr<BasisSet> primary,
          size_t memory, Options &options) : PKManager(primary,memory,options) {
    psio_ = psio;
    AIO_ = new AIOHandler(psio_);
    max_batches_ = options.get_int("MAX_BUCKETS");
    pk_file_ = PSIF_SO_PK;

    // No current writing since we are constructing
    writing_ = false;
}


void PKMgrDisk::initialize() {
    batch_sizing();
    prestripe_files();
    print_batches();
    allocate_buffers();
}

void PKMgrDisk::batch_sizing() {

    double batch_thresh = 0.1;

    ijklBasisIterator AOintsiter(nbf_,sieve_);

    size_t old_pq = 0;
    size_t old_max = 0;
    size_t nintpq = 0;
    size_t nintbatch = 0;
    size_t pq = 0;
    size_t pb, qb, rb, sb;
    int batch = 0;

    batch_index_min_.push_back(0);
    batch_pq_min_.push_back(0);
    batch_for_pq_.push_back(0);
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
                batch_for_pq_.pop_back();
                ++batch;
                batch_for_pq_.push_back(batch);
                batch_index_min_.push_back(old_max);
                batch_pq_min_.push_back(old_pq);
                nintbatch = nintpq;
            }
            nintpq = 1;
            old_pq = pq;
            batch_for_pq_.push_back(batch);
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
      outfile->Printf( "  PKJK: maximum number of batches exceeded\n") ;
      throw PSIEXCEPTION( "  PK computation needs %d batches, max. number: %d\n", nbatches, max_batches_);
    }


}

void PKMgrDisk::print_batches() {
    PKManager::print_batches();
    // Print batches for the user and for control
    for(int batch = 0; batch < batch_pq_min_.size(); ++batch){
        outfile->Printf("\tBatch %3d pq = [%8zu,%8zu] index = [%14zu,%zu] size = %12zu\n",
                batch + 1,
                batch_pq_min_[batch],batch_pq_max_[batch],
                batch_index_min_[batch],batch_index_max_[batch],
                batch_index_max_[batch] - batch_index_min_[batch]);
    }

}

void PKMgrDisk::open_PK_file() {
    psio_->open(pk_file_, PSIO_OPEN_OLD);
}

void PKMgrDisk::close_PK_file(bool keep) {
    psio_->close(pk_file_,keep);
}

void PKMgrDisk::prepare_JK(std::vector<SharedMatrix> D) {
    if (writing()) {
        finalize_PK();
        set_writing(false);
    } else {
        open_PK_file(true);
    }
    form_D_vec(D);
}


void PKMgrDisk::finalize_PK() {

    timer_on("AIO synchronize");
    AIO_->synchronize();
    timer_off("AIO synchronize");

//    // Can get rid of the pre-striping labels
//      for(int i = 0; i < label_J_[0].size(); ++i ) {
//          delete [] label_J_[0][i];
//      }
//      label_J_[0].clear();
//      for(int i = 0; i < label_K_[0].size(); ++i ) {
//          delete [] label_K_[0][i];
//      }
//      label_K_[0].clear();
    for(int i = 0; i < nthreads_; ++i) {
        iobuffers_[i].reset();
    }

}

void PKMgrDisk::form_J(std::vector<SharedMatrix> J, bool exch) {
    make_J_vec(J);

    // Now loop over batches
    for(int batch = 0; batch < batch_pq_min_.size(); ++batch) {
        size_t min_index = batch_index_min_[batch];
        size_t max_index = batch_index_max_[batch];
        size_t batch_size = max_index - min_index;
        size_t min_pq = batch_pq_min_[batch];
        size_t max_pq = batch_pq_max_[batch];
        double* j_block = new double[batch_size];

        char* label;
        if (exch) {
            label = get_label_K(batch);
        } else {
            label = get_label_J(batch);
        }
        psio_->read_entry(pk_file_, label, (char *) j_block, batch_size * sizeof(double));

        // Read one entry, use it for all density matrices
        for(int N = 0; N < J.size(); ++N) {
            double* D_vec = D_vec_[N];
            double* J_vec = JK_vec_[N];
            double* j_ptr = j_block;
        //TODO Could consider parallelizing this loop
            for(size_t pq = min_pq; pq < max_pq; ++pq) {
                double D_pq = D_vec[pq];
                double *D_rs = D_vec;
                double J_pq = 0.0;
                double *J_rs = J_vec;
                for(size_t rs = 0; rs <= pq; ++rs) {
//DEBUG                    if(exch) { // && rs == 0) {
//DEBUG                      outfile->Printf("PK int (%lu|%lu) = %20.16f\n",pq,rs,*j_ptr);
//DEBUG                    }
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
    get_results(J);
}

void PKMgrDisk::finalize_JK() {
    finalize_D();
    close_PK_file(true);
}

PKMgrReorder::PKMgrReorder(boost::shared_ptr<PSIO> psio, boost::shared_ptr<BasisSet> primary,
                           size_t memory, Options &options) :
    PKMgrDisk(psio,primary,memory,options)  {
    max_mem_buf_ = options.get_int("MAX_MEM_BUF");
}

// Pre-striping the PK file
void PKMgrReorder::prestripe_files() {
    psio_->open(pk_file_, PSIO_OPEN_NEW);
    // Create the files ? Pre-stripe them.
    for(int batch = 0; batch < batch_index_min_.size(); ++batch) {
        size_t batch_size = batch_index_max_[batch] - batch_index_min_[batch];
        // We need to keep the labels around in a vector
        label_J_.push_back(get_label_J(batch));
        AIO_->zero_disk(pk_file_,label_J_[batch],1,batch_size);
        label_K_.push_back(get_label_K(batch));
        AIO_->zero_disk(pk_file_,label_K_[batch],1,batch_size);
    }

}

void PKMgrReorder::allocate_buffers() {
    // Factor 2 because we need memory for J and K buffers
    size_t mem_per_thread = memory_ / (2 * nthreads_);
    // Factor 2 because we need 2 buffers for asynchronous I/O
    size_t buf_size = mem_per_thread / 2;
    if (max_mem_buf_ != 0) buf_size = std::min(buf_size, max_mem_buf_);

    // Number of tasks needed with this buffer size
    ntasks_ = pk_size_ / buf_size + 1;
    // Is there less tasks than threads ?
    if(ntasks_ < nthreads_) {
        ntasks_ = ntasks_ * nthreads_;
        size_t tmp_size = pk_size_ / ntasks_ + 1;
        buf_size = tmp_size;
        ntasks_ = pk_size_ / buf_size + 1;
    }
    size_t buf_per_thread = std::min(mem_per_thread / buf_size, ntasks_ / nthreads_);
    // Some printing for us
    outfile->Printf("  Task number: %lu\n",ntasks_);
    outfile->Printf("  Buffer size: %lu\n",buf_size);
    outfile->Printf("  Buffer per thread: %lu\n",buf_per_thread);

    // Ok, now we have the size of a buffer and how many buffers
    // we want for each thread. We can allocate IO buffers.
    for(int i = 0; i < nthreads_; ++i) {
        iobuffers_.push_back(new PKWrkrReord(primary,AIO_,pk_file_,buf_size,buf_per_thread));
    }

}

void PKMgrReorder::form_PK() {
    compute_integrals();
    set_writing(true);
}

void PKMgrReorder::finalize_PK() {
    timer_on("AIO synchronize");
    AIO_->synchronize();
    timer_off("AIO synchronize");

    // Can get rid of the pre-striping labels
    for(int i = 0; i < label_J_.size(); ++i ) {
        delete [] label_J_[i];
    }
    label_J_.clear();
    for(int i = 0; i < label_K_.size(); ++i ) {
        delete [] label_K_[i];
    }
    label_K_.clear();
    for(int i = 0; i < nthreads_; ++i) {
        iobuffers_[i].reset();
    }
}

PKMgrYoshimine::PKMgrYoshimine(boost::shared_ptr<PSIO> psio, boost::shared_ptr<BasisSet> primary,
                               size_t memory, Options &options) :
    PKMgrDisk(psio,primary,memory,options) {
    iwl_file_J_ = PSIF_SO_PKSUPER1;
    iwl_file_K_ = PSIF_SO_PKSUPER2;
    ints_per_buf_ = IWL_INTS_PER_BUF;
    iwl_int_size_ = ints_per_buf_ * (4L * sizeof(Label) + sizeof(Value)) + 2L * sizeof(int);
}

void PKMgrYoshimine::prestripe_files() {

    psio_->open(iwl_file_J_, PSIO_OPEN_NEW);
    // Pre-stripe the new file
    // Number of IWL buffers necessary to write all integrals
    size_t num_iwlbuf = pk_size_ / ints_per_buf_ + 1;
    // The last buffer of each batch may be partially filled only
    num_iwlbuf += batch_index_min_.size();
    size_t iwlsize_bytes = num_iwlbuf * iwlintsize_;
    size_t iwlsize = iwlsize_bytes / sizeof(double) + 1;
    AIO_->zero_disk(iwl_file_J_,IWL_KEY_BUF,1,iwlsize);

    // And we do the same for the file containing K buckets
    // We need at most twice as much IWL buffers since integrals can be
    // pre-sorted in two different buckets for K
    psio_->open(iwl_file_K_, PSIO_OPEN_NEW);
    AIO_->zero_disk(iwl_file_K_, IWL_KEY_BUF, 1, 2*iwlsize);
}

void PKMgrYoshimine::allocate_buffers() {
    int buf_per_thread = batch_index_min_.size();
    // Factor 2 because we need an address for J and an address for K
    // J and K addresses for the same buckets are stored contiguously, i.e.
    // element [0] is J address for first bucket and element [1] the
    // K address for the first bucket
    size_t* current_pos = new size_t[2 * buf_per_thread];

    current_pos[0] = 0;
    current_pos[1] = 0;
    //Current position of the bucket in the IWL file
    // For K, each integral can potentially go into two buckets, so we
    // give each bucket twice the J bucket's size. It's wasting quite a bit
    // of space, other solution is to explicitly count integrals.
    for(int i = 1; i < buf_per_thread; ++i) {
        size_t batchsize = batch_index_max_[i - 1] - batch_index_min_[i - 1];
        size_t iwlperbatch = batchsize / ints_per_buf_ + 1;
        current_pos[2 * i] = iwlperbatch * iwlintsize_ + current_pos[2 * i - 2];
        current_pos[2 * i + 1] = 2 * (iwlperbatch * iwlintsize_) + current_pos[2 * i - 1];
    }

    // Ok, now we have the size of a buffer and how many buffers
    // we want for each thread. We can allocate IO buffers.
    for(int i = 0; i < nthreads_; ++i) {
        iobuffers_.push_back(new PKWrkrIWL(primary,AIO_,iwl_file_J_,iwl_file_K_,ints_per_buf_,batch_for_pq_));
    }

}

void PKMgrYoshimine::form_PK() {
    compute_integrals();

    // Now we need to actually read the integral file and sort
    // the buckets to write them to the PK file
    // We deallocate buffers and synchronize AIO writing in there.

    sort_ints();

}

void PKMgrYoshimine::compute_integrals() {
#pragma omp parallel for schedule(dynamic) num_threads(nthreads_)
    for(int i = 0; i < primary_->nshell(); ++i) {
      int num_uniq_pk;
      int P_arr[3];
      int Q_arr[3];
      int R_arr[3];
      int S_arr[3];
      int P, Q, R, S;
      int thread = 0;
      #ifdef _OPENMP
        thread = omp_get_thread_num();
      #endif
      for(int j = 0; j <= i; ++j) {
        for (int k = 0; k <= j; ++k) {
          for(int l = 0; l <= k; ++l) {
            P_arr[0] = i;
            Q_arr[0] = j;
            R_arr[0] = k;
            S_arr[0] = l;
            // Now we take care of all symmetry cases
            if((i == j && i == k) || (j == k && j == l)) {
                num_uniq_pk = 1;
            } else if (i == k || j == l) {
                num_uniq_pk = 2;
                P_arr[1] = i;
                Q_arr[1] = k;
                R_arr[1] = j;
                S_arr[1] = l;
            } else if (j == k) {
                num_uniq_pk = 2;
                P_arr[1] = i;
                Q_arr[1] = l;
                R_arr[1] = j;
                S_arr[1] = k;
            } else if (i == j || k == l) {
                num_uniq_pk = 2;
                P_arr[1] = i;
                Q_arr[1] = k;
                R_arr[1] = j;
                S_arr[1] = l;
            } else {
                num_uniq_pk = 3;
                P_arr[1] = i;
                Q_arr[1] = k;
                R_arr[1] = j;
                S_arr[1] = l;

                P_arr[2] = i;
                Q_arr[2] = l;
                R_arr[2] = j;
                S_arr[2] = k;
            }

            for(int npk = 0; npk < num_uniq_pk; ++npk) {
              P = P_arr[npk];
              Q = Q_arr[npk];
              R = R_arr[npk];
              S = S_arr[npk];
              //Sort shells based on AM to save ERI some work doing permutation resorting
              if (primary_->shell(P).am() < primary_->shell(Q).am()) {
                  std::swap(P,Q);
              }
              if (primary_->shell(R).am() < primary_->shell(S).am()) {
                  std::swap(R,S);
              }
              if (primary_->shell(P).am() + primary_->shell(Q).am() >
                      primary_->shell(R).am() + primary_->shell(S).am()) {
                  std::swap(P, R);
                  std::swap(Q, S);
              }
    //          printf("Computing integral <%i %i|%i %i>\n", p, q, r, s);
              tb[thread]->compute_shell(P, Q, R, S);
              integrals_buffering(tb[thread]->buffer(),P,Q,R,S);
            }

          }
        }
      }
    } // end of parallelized loop

    // We write all remaining buffers to disk.
    write();

}

void PKMgrYoshimine::write() {
    // Each thread buffer has several buckets that may be partially full.
    // Concatenate all internal buckets to the ones of thread buffer 0.
    double val;
    size_t i, j, k, l;
    PKWorker* buf0 = iobuffers_[0];
    PKWorker* buftarget;
    //TODO Man if that loop actually works I'll be lucky
    for(int t = 1; t < nthreads_; ++t) {
        buftarget = iobuffers_[t];
        unsigned int nbufs = buftarget->nbuf();
        // Factor 2 to get buffers for J and K
        for(int buf = 0; buf < 2 * nbufs; ++buf) {
            while(buftarget->pop_value(buf,val,i,j,k,l)) {
                buf0->insert_value(buf,val,i,j,k,l);
            }
        }

    }
    // By now we should have transferred all remaining integrals to
    // buffer 0, we now flush it.
    buf0->flush();

}

void PKMgrYoshimine::sort_ints() {
    // compute max batch size
    size_t max_size = 0;
    size_t batch_size;
    int nbatches = batch_index_min_.size();

    for(int i = 0; i < nbatches; ++i) {
        batch_size = batch_index_max_[i] - batch_index_min_[i];
        if(batch_size > max_size) max_size = batch_size;
    }

    //TODO: other possibility here: have an array half the size to do AIO
    // Allocate a big array to contain the batch.
    double* twoel_ints = new double[max_size];
    ::memset((void*)twoel_ints,'\0',max_size * sizeof(double));

    // At this point we need to close the IWL file to create
    // an IWL object which will open it. Dumb but minor inconvenience (hopefully)
    // We also open the PK file
    psio_->open(pk_file_, PSIO_OPEN_NEW);
    //TODO: this function may need renaming
    finalize_PK();
    set_writing(false);
    close_iwl_buckets();

    generate_J_PK(twoel_ints,max_size);
    // Need to reset to zero the two-el integral array
    ::memset((void*)twoel_ints,'\0',max_size * sizeof(double));
    generate_K_PK(twoel_ints,max_size);

    //delete two-el int array
    delete [] twoel_ints;

    psio_->close(pk_file_,1);

}

void PKMgrYoshimine::close_iwl_buckets() {
    psio_->close(iwl_file_J_, 1);
    psio_->close(iwl_file_K_, 1);
}

void PKMgrYoshimine::generate_J_PK(double *twoel_ints, size_t max_size) {

    IWL inbuf(psio_.get(),iwl_file_J_,0.0, 1, 0);

    int idx, id;
    size_t p, q, r, s;
    Label* lblptr = inbuf.labels();
    Value* valptr = inbuf.values();
    int lastbuf;

    size_t offset, maxind, nintegrals;
    size_t pqrs;

    int batch = 0;
    int nbatches = batch_index_min_.size();
    while(batch < nbatches) {
        inbuf.fetch();
        offset = batch_index_min_[batch];
        maxind = batch_index_max_[batch];
        nintegrals = batch_index_max_[batch] - batch_index_min_[batch];

        for(idx = 0; idx < inbuf.buffer_count(); ++idx) {
            id = idx;
            p = lblptr[id++];
            q = lblptr[id++];
            r = lblptr[id++];
            s = lblptr[id];

            pqrs = INDEX4(p,q,r,s);
            //TODO: remove the check, for debug only
            if(pqrs < offset || pqrs > maxind) {
                throw PSIEXCEPTION("Integrals not written properly \n");
            }
            twoel_ints[pqrs - offset] += valptr[idx];
        }

        lastbuf = inbuf.last_buffer();
        if(lastbuf) {
            // We just completed a batch, write it to disk
            char* label = get_label_J(batch);
            // Divide diagonal elements by two
            for(size_t pq = batch_pq_min_[batch]; pq < batch_pq_max_[batch]; ++pq) {
                pqrs = INDEX2(pq,pq);
                twoel_ints[pqrs - offset] *= 0.5;
            }
            psio_->write_entry(pk_file_, label, (char*)twoel_ints, nintegrals * sizeof(double));
            delete [] label;
            ++batch;
            if(batch < nbatches) {
                ::memset((void*)twoel_ints,'\0',max_size * sizeof(double));
            }
        }

    }

    inbuf.set_keep_flag(false);

}

void PKMgrYoshimine::generate_K_PK(double *twoel_ints, size_t max_size) {

    IWL inbuf(psio_.get(),iwl_file_K_,0.0, 1, 0);

    int idx,id;
    size_t p, q, r, s;
    Label* lblptr = inbuf.labels();
    Value* valptr = inbuf.values();
    int lastbuf;

    size_t offset, maxind, nintegrals;
    size_t pqrs;

    int batch = 0;
    int nbatches = batch_index_min_.size();
    while(batch < nbatches) {
        inbuf.fetch();
        offset = batch_index_min_[batch];
        maxind = batch_index_max_[batch];
        nintegrals = batch_index_max_[batch] - batch_index_min_[batch];

        for(idx = 0; idx < inbuf.buffer_count(); ++idx) {
            id = idx;
            p = lblptr[id++];
            q = lblptr[id++];
            r = lblptr[id++];
            s = lblptr[id];


            // K first sort
            pqrs = INDEX4(p,r,q,s);
//DEBUG            size_t pqd = INDEX2(p,r);
//DEBUG            size_t rsd = INDEX2(q,s);
            if(pqrs <= maxind && pqrs >= offset) {
                if(p == r || q == s) {
                    twoel_ints[pqrs - offset] += valptr[idx];
                } else {
                    twoel_ints[pqrs - offset] += 0.5 * valptr[idx];
                }
//DEBUG                if(pqd == 0 && rsd == 0) {
//DEBUG                    outfile->Printf("Int is %20.16f\n",twoel_ints[pqrs - offset]);
//DEBUG                }
            }

            // K second sort
            if(p != q && r != s) {
                pqrs = INDEX4(p,s,q,r);
                if(pqrs <= maxind && pqrs >= offset) {
                    if(p == s || q == r) {
                        twoel_ints[pqrs - offset] += valptr[idx];
                    } else {
                        twoel_ints[pqrs - offset] += 0.5 * valptr[idx];
                    }
//DEBUG                if(pqd == 0 && rsd == 0) {
//DEBUG                    outfile->Printf("Int is %20.16f\n",twoel_ints[pqrs - offset]);
//DEBUG                }
                }

            }
        }

        lastbuf = inbuf.last_buffer();
        if(lastbuf) {
            // We just completed a batch, write it to disk
            char* label = get_label_K(batch);
            // Divide by two diagonal elements
            for(size_t pq = batch_pq_min_[batch]; pq < batch_pq_max_[batch]; ++pq) {
                pqrs = INDEX2(pq,pq);
                twoel_ints[pqrs - offset] *= 0.5;
            }
            psio_->write_entry(pk_file_, label, (char*)twoel_ints, nintegrals * sizeof(double));
            delete [] label;
            ++batch;
            if(batch < nbatches) {
                ::memset((void*)twoel_ints,'\0',max_size * sizeof(double));
            }
        }

    }

    inbuf.set_keep_flag(false);

}

void PKMgrInCore::initialize() {
    print_batches();
    allocate_buffers();
}

void PKMgrInCore::print_batches() {
    PKManager::print_batches();
    outfile->Printf("  Performing in-core PK\n");
    outfile->Printf("  Using %lu doubles for integral storage.\n", 2 * pk_size());
}

void PKMgrInCore::allocate_buffers() {
    // Need to allocate two big arrays
    J_ints_ = new double[pk_size_];
    K_ints_ = new double[pk_size_];
    ::memset((void*) J_ints_.get(), '\0', pk_size_ * sizeof(double));
    ::memset((void*) K_ints_.get(), '\0', pk_size_ * sizeof(double));

    // Now we allocate a derived class of IOBuffer_PK that takes care of
    // giving the tasks to the threads and storing results properly.
    // Other possibility: make all accesses to J_ints_ and K_ints_ atomic,
    // probably more expensive

    //TODO We probably need more buffer for more tasks for better load balancing
    size_t buffer_size = pk_size_ / nthreads_;
    size_t lastbuf = pk_size_ % nthreads_;
    size_t start = 0;

    for(size_t i = 0; i < nthreads_; ++i) {
        start = i * buffer_size;
//DEBUG        outfile->Printf("start is %lu\n",start);
        SharedPKWrkr buf;
        if(i < nthreads_ - 1) {
            buf = new PKWrkrInCore(primary_,buffer_size,0,&J_ints_[start],&K_ints_[start]);
        } else {
            buf = new PKWrkrInCore(primary_,buffer_size,lastbuf,&J_ints_[start],&K_ints_[start]);
        }
        iobuffers_.push_back(buf);
        ntasks_ = nthreads_;
    }
}

void PKMgrInCore::form_PK() {
    compute_integrals();
    finalize_PK();
}

void PKMgrInCore::finalize_PK() {
    for(int i = 0; i < nthreads_; ++i) {
        iobuffers_[i].reset();
    }
}

void PKMgrInCore::prepare_JK(std::vector<SharedMatrix> D) {
    form_D_vec(D);
}

void PKMgrInCore::form_J(std::vector<SharedMatrix> J, bool exch) {

    make_J_vec(J);

    for(int N = 0; N < J.size(); ++N) {
        double* D_vec = D_vec_[N];
        double* J_vec = JK_vec_[N];
        double* j_ptr;
        if(exch) {
            j_ptr = K_ints_.get();
        } else {
            j_ptr = J_ints_.get();
        }
//TODO We should totally parallelize this loop now.
        for(size_t pq = 0; pq < pk_pairs_; ++pq) {
            double D_pq = D_vec[pq];
            double *D_rs = D_vec;
            double J_pq = 0.0;
            double *J_rs = J_vec;
            for(size_t rs = 0; rs <= pq; ++rs) {
//DEBUG                    if(exch) { // && rs == 0) {
//DEBUG                      outfile->Printf("PK int (%lu|%lu) = %20.16f\n",pq,rs,*j_ptr);
//DEBUG
                J_pq += *j_ptr * (*D_rs);
                *J_rs += *j_ptr * D_pq;
                ++D_rs;
                ++J_rs;
                ++j_ptr;
            }
            J_vec[pq] += J_pq;
        }
    }

    get_results(J);

}

void PKMgrInCore::finalize_JK() {
    finalize_D();
}

}  // End namespace pk
}  // End namespace psi
