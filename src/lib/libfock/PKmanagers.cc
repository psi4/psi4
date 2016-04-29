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
#include <libmints/basisset.h>
#include <libmints/typedefs.h>
#include <libmints/matrix.h>
#include <libmints/sieve.h>
#include <libqt/qt.h>
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

PKManager::PKManager(boost::shared_ptr<BasisSet> primary, size_t memory, Options& options) :
primary_(primary), memory_(memory), options_(options) {

    nbf_ = primary_->nbf();

    pk_pairs_ = (size_t) nbf_ * ((size_t) nbf_ + 1) / 2;
    pk_size_ = pk_pairs_ * (pk_pairs_ + 1) / 2;
    cutoff_ = 1.0e-12;
    if(options["INTS_TOLERANCE"].has_changed()) {
        cutoff_ = options.get_double("INTS_TOLERANCE");
    }
    ntasks_ = 0;
    sieve_ = std::shared_ptr<ERISieve>(new ERISieve(primary_, cutoff_));

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
//DEBUG        outfile->Printf("Starting task %d\n", i);
        for(buf->first_quartet(i); buf->more_work(); buf->next_quartet()) {
            unsigned int P = buf->P();
            unsigned int Q = buf->Q();
            unsigned int R = buf->R();
            unsigned int S = buf->S();
            //Sort shells based on AM to save ERI some work doing permutation resorting
            if (primary()->shell(P).am() < primary()->shell(Q).am()) {
                std::swap(P,Q);
            }
            if (primary()->shell(R).am() < primary()->shell(S).am()) {
                std::swap(R,S);
            }
            if (primary()->shell(P).am() + primary()->shell(Q).am() >
                    primary()->shell(R).am() + primary()->shell(S).am()) {
                std::swap(P, R);
                std::swap(Q, S);
            }
//DEBUG#pragma omp critical
//DEBUG            outfile->Printf("Computing shell <%d %d|%d %d>\n",P,Q,R,S);
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
    if (nshqu > nsh_u) {
        outfile->Printf("  %7.2f percent of shell quartets recomputed by reordering.\n", (nshqu - nsh_u) / float(nsh_u) * 100);
    }

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
        if(fabs(val) > cutoff_) {
//DEBUG#pragma omp critical
//DEBUG{
//DEBUG            if(INDEX2(k,l) == 0) {
//DEBUG              outfile->Printf("Integral <%d %d|%d %d>,<%d|%d> = %16.12f\n",i,j,k,l,INDEX2(i,j),INDEX2(k,l),val);
//DEBUG              outfile->Printf("Integral <%d %d|%d %d> = %16.12f\n",i,j,k,l,val);
//DEBUG            }
//DEBUG}
            iobuffers_[thread]->fill_values(val, i, j, k, l);
        }
    }

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
    JK_vec_.clear();

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

PKMgrDisk::PKMgrDisk(boost::shared_ptr<PSIO> psio, boost::shared_ptr<BasisSet> primary,
          size_t memory, Options &options) : PKManager(primary,memory,options) {
    psio_ = psio;
    AIO_ = std::shared_ptr<AIOHandler>(new AIOHandler(psio_));
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

    ijklBasisIterator AOintsiter(nbf(),sieve());

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
            if (nintbatch > memory()) {
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
        if (((double) size_lastb / memory()) < batch_thresh) {
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
      outfile->Printf( "  PK computation needs %d batches, max. number: %d\n", nbatches, max_batches_);
      throw PSIEXCEPTION("  PK Failure: max batches exceeded\n");

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

void PKMgrDisk::write() {
    get_buffer()->write(batch_index_min_,batch_index_max_,pk_pairs());
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
        open_PK_file();
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
    for(int i = 0; i < nthreads(); ++i) {
        buffer(i).reset();
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
            label = PKWorker::get_label_K(batch);
        } else {
            label = PKWorker::get_label_J(batch);
        }
        psio_->read_entry(pk_file_, label, (char *) j_block, batch_size * sizeof(double));

        // Read one entry, use it for all density matrices
        for(int N = 0; N < J.size(); ++N) {
            double* D_vec = D_glob_vecs(N);
            double* J_vec = JK_glob_vecs(N);
            double* j_ptr = j_block;
        //TODO Could consider parallelizing this loop
            for(size_t pq = min_pq; pq < max_pq; ++pq) {
                double D_pq = D_vec[pq];
                double *D_rs = D_vec;
                double J_pq = 0.0;
                double *J_rs = J_vec;
                for(size_t rs = 0; rs <= pq; ++rs) {
//DEBUG                    if(!exch && rs == 0) {
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
    psio()->open(pk_file(), PSIO_OPEN_NEW);
    // Create the files ? Pre-stripe them.
    for(int batch = 0; batch < batch_ind_min().size(); ++batch) {
        size_t batch_size = batch_ind_max()[batch] - batch_ind_min()[batch];
        // We need to keep the labels around in a vector
        label_J_.push_back(PKWorker::get_label_J(batch));
        AIO()->zero_disk(pk_file(),label_J_[batch],1,batch_size);
        label_K_.push_back(PKWorker::get_label_K(batch));
        AIO()->zero_disk(pk_file(),label_K_[batch],1,batch_size);
    }

}

void PKMgrReorder::allocate_buffers() {
    // Factor 2 because we need memory for J and K buffers
    size_t mem_per_thread = memory() / (2 * nthreads());
    // Factor 2 because we need 2 buffers for asynchronous I/O
    size_t buf_size = mem_per_thread / 2;
    if (max_mem_buf_ != 0) buf_size = std::min(buf_size, max_mem_buf_);

    // Number of tasks needed with this buffer size
    set_ntasks(pk_size() / buf_size + 1);
    // Is there less tasks than threads ?
    if(ntasks() < nthreads()) {
        set_ntasks(ntasks() * nthreads());
        size_t tmp_size = pk_size() / ntasks() + 1;
        buf_size = tmp_size;
        set_ntasks(pk_size() / buf_size + 1);
    }
    size_t buf_per_thread = std::min(mem_per_thread / buf_size, ntasks() / nthreads());
    // Some printing for us
    outfile->Printf("  Task number: %lu\n",ntasks());
    outfile->Printf("  Buffer size: %lu\n",buf_size);
    outfile->Printf("  Buffer per thread: %lu\n",buf_per_thread);

    // Ok, now we have the size of a buffer and how many buffers
    // we want for each thread. We can allocate IO buffers.
    for(int i = 0; i < nthreads(); ++i) {
        fill_buffer(SharedPKWrkr(new PKWrkrReord(primary(),sieve(),AIO(),pk_file(),buf_size,buf_per_thread)));
    }

}

void PKMgrReorder::form_PK() {
    compute_integrals();
    set_writing(true);
}

void PKMgrReorder::finalize_PK() {
    timer_on("AIO synchronize");
    AIO()->synchronize();
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
    for(int i = 0; i < nthreads(); ++i) {
        buffer(i).reset();
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

    psio()->open(iwl_file_J_, PSIO_OPEN_NEW);
    // Pre-stripe the new file
    // Number of IWL buffers necessary to write all integrals
    size_t num_iwlbuf = pk_size() / ints_per_buf_ + 1;
    // The last buffer of each batch may be partially filled only
    num_iwlbuf += batch_ind_min().size();
    size_t iwlsize_bytes = num_iwlbuf * iwl_int_size_;
    size_t iwlsize = iwlsize_bytes / sizeof(double) + 1;
    AIO()->zero_disk(iwl_file_J_,IWL_KEY_BUF,1,iwlsize);

    // And we do the same for the file containing K buckets
    // We need at most twice as much IWL buffers since integrals can be
    // pre-sorted in two different buckets for K
    psio()->open(iwl_file_K_, PSIO_OPEN_NEW);
    AIO()->zero_disk(iwl_file_K_, IWL_KEY_BUF, 1, 2*iwlsize);
}

void PKMgrYoshimine::allocate_buffers() {
    int buf_per_thread = batch_ind_min().size();
    // Factor 2 because we need an address for J and an address for K
    // J and K addresses for the same buckets are stored contiguously, i.e.
    // element [0] is J address for first bucket and element [1] the
    // K address for the first bucket
    // we actually need the boost shared_ptr, the std::shared_ptr does
    // not support array syntax
    boost::shared_ptr<size_t []> current_pos(new size_t[2 * buf_per_thread]);

    current_pos[0] = 0;
    current_pos[1] = 0;
    //Current position of the bucket in the IWL file
    // For K, each integral can potentially go into two buckets, so we
    // give each bucket twice the J bucket's size. It's wasting quite a bit
    // of space, other solution is to explicitly count integrals.
    for(int i = 1; i < buf_per_thread; ++i) {
        size_t batchsize = batch_ind_max()[i - 1] - batch_ind_min()[i - 1];
        size_t iwlperbatch = batchsize / ints_per_buf_ + 1;
        current_pos[2 * i] = iwlperbatch * iwl_int_size_ + current_pos[2 * i - 2];
        current_pos[2 * i + 1] = 2 * (iwlperbatch * iwl_int_size_) + current_pos[2 * i - 1];
    }

    // Ok, now we have the size of a buffer and how many buffers
    // we want for each thread. We can allocate IO buffers.
    for(int i = 0; i < nthreads(); ++i) {
        fill_buffer(SharedPKWrkr(new PKWrkrIWL(primary(),sieve(),AIO(),iwl_file_J_,iwl_file_K_,ints_per_buf_,batch_for_pq(),current_pos)));
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

    // Get an AO integral factory
    boost::shared_ptr<IntegralFactory> intfact(new IntegralFactory(primary()));

    // Get ERI object, one per thread
    std::vector<boost::shared_ptr<TwoBodyAOInt> > tb;
    for(int i = 0; i < nthreads(); ++i) {
        tb.push_back(boost::shared_ptr<TwoBodyAOInt>(intfact->erd_eri()));
    }

    // Loop over significant shell pairs from ERISieve
    const std::vector< std::pair<int, int> >& sh_pairs = sieve()->shell_pairs();
    size_t npairs = sh_pairs.size();
#pragma omp parallel for schedule(dynamic) num_threads(nthreads())
    for(size_t i = 0; i < npairs; ++i) {
        int PP = sh_pairs[i].first;
        int QQ = sh_pairs[i].second;
        int thread = 0;
      #ifdef _OPENMP
        thread = omp_get_thread_num();
      #endif

        for(size_t j = 0; j <= i; ++j) {
            int RR = sh_pairs[j].first;
            int SS = sh_pairs[j].second;
            if(sieve()->shell_significant(PP,QQ,RR,SS)) {
                int P = PP;
                int Q = QQ;
                int R = RR;
                int S = SS;
                //Sort shells based on AM to save ERI some work doing permutation resorting
                if (primary()->shell(P).am() < primary()->shell(Q).am()) {
                    std::swap(P,Q);
                }
                if (primary()->shell(R).am() < primary()->shell(S).am()) {
                    std::swap(R,S);
                }
                if (primary()->shell(P).am() + primary()->shell(Q).am() >
                        primary()->shell(R).am() + primary()->shell(S).am()) {
                    std::swap(P, R);
                    std::swap(Q, S);
                }
//DEBUG                outfile->Printf("Computing shell <%i %i|%i %i>\n", P, Q, R, S);
                tb[thread]->compute_shell(P, Q, R, S);
                integrals_buffering(tb[thread]->buffer(),P,Q,R,S);
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
    SharedPKWrkr buf0 = buffer(0);
    SharedPKWrkr buftarget;
    for(int t = 1; t < nthreads(); ++t) {
        buftarget = buffer(t);
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
    int nbatches = batch_ind_min().size();

    for(int i = 0; i < nbatches; ++i) {
        batch_size = batch_ind_max()[i] - batch_ind_min()[i];
        if(batch_size > max_size) max_size = batch_size;
    }

    //TODO: other possibility here: have an array half the size to do AIO
    // Allocate a big array to contain the batch.
    double* twoel_ints = new double[max_size];
    ::memset((void*)twoel_ints,'\0',max_size * sizeof(double));

    // At this point we need to close the IWL file to create
    // an IWL object which will open it. Dumb but minor inconvenience (hopefully)
    // We also open the PK file
    psio()->open(pk_file(), PSIO_OPEN_NEW);
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

    psio()->close(pk_file(),1);

}

void PKMgrYoshimine::close_iwl_buckets() {
    psio()->close(iwl_file_J_, 1);
    psio()->close(iwl_file_K_, 1);
}

void PKMgrYoshimine::generate_J_PK(double *twoel_ints, size_t max_size) {

    IWL inbuf(psio().get(),iwl_file_J_,0.0, 1, 0);

    int idx, id;
    size_t p, q, r, s;
    Label* lblptr = inbuf.labels();
    Value* valptr = inbuf.values();
    int lastbuf;

    size_t offset, maxind, nintegrals;
    size_t pqrs;

    int batch = 0;
    int nbatches = batch_ind_min().size();
    while(batch < nbatches) {
        inbuf.fetch();
        offset = batch_ind_min()[batch];
        maxind = batch_ind_max()[batch];
        nintegrals = batch_ind_max()[batch] - batch_ind_min()[batch];

        for(idx = 0; idx < inbuf.buffer_count(); ++idx) {
            id = 4 * idx;
            p = lblptr[id++];
            q = lblptr[id++];
            r = lblptr[id++];
            s = lblptr[id];

//DEBUG            size_t pqd = INDEX2(p,q);
//DEBUG            size_t rsd = INDEX2(r,s);
//DEBUG            if (rsd == 0) {
//DEBUG              outfile->Printf("Integral <%d |%d> = %16.12f added\n",pqd,rsd,valptr[idx]);
//DEBUG            }
            pqrs = INDEX4(p,q,r,s);
            twoel_ints[pqrs - offset] += valptr[idx];
        }

        lastbuf = inbuf.last_buffer();
        if(lastbuf) {
            // We just completed a batch, write it to disk
            char* label = PKWorker::get_label_J(batch);
            // Divide diagonal elements by two
            for(size_t pq = batch_pq_min()[batch]; pq < batch_pq_max()[batch]; ++pq) {
                pqrs = INDEX2(pq,pq);
                twoel_ints[pqrs - offset] *= 0.5;
            }
            psio()->write_entry(pk_file(), label, (char*)twoel_ints, nintegrals * sizeof(double));
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

    IWL inbuf(psio().get(),iwl_file_K_,0.0, 1, 0);

    int idx,id;
    size_t p, q, r, s;
    Label* lblptr = inbuf.labels();
    Value* valptr = inbuf.values();
    int lastbuf;

    size_t offset, maxind, nintegrals;
    size_t pqrs;

    int batch = 0;
    int nbatches = batch_ind_min().size();
    while(batch < nbatches) {
//DEBUG        outfile->Printf("fetching for batch %d\n",batch);
        inbuf.fetch();
        offset = batch_ind_min()[batch];
        maxind = batch_ind_max()[batch];
        nintegrals = batch_ind_max()[batch] - batch_ind_min()[batch];

        for(idx = 0; idx < inbuf.buffer_count(); ++idx) {
            id = 4 * idx;
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
            char* label = PKWorker::get_label_K(batch);
            // Divide by two diagonal elements
            for(size_t pq = batch_pq_min()[batch]; pq < batch_pq_max()[batch]; ++pq) {
                pqrs = INDEX2(pq,pq);
                twoel_ints[pqrs - offset] *= 0.5;
            }
            psio()->write_entry(pk_file(), label, (char*)twoel_ints, nintegrals * sizeof(double));
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
    J_ints_ = std::unique_ptr< double [] >(new double[pk_size()]);
    K_ints_ = std::unique_ptr< double [] >(new double[pk_size()]);
    ::memset((void*) J_ints_.get(), '\0', pk_size() * sizeof(double));
    ::memset((void*) K_ints_.get(), '\0', pk_size() * sizeof(double));

    // Now we allocate a derived class of IOBuffer_PK that takes care of
    // giving the tasks to the threads and storing results properly.
    // Other possibility: make all accesses to J_ints_ and K_ints_ atomic,
    // probably more expensive

    //TODO We probably need more buffer for more tasks for better load balancing
    size_t buffer_size = pk_size() / nthreads();
    size_t lastbuf = pk_size() % nthreads();

    for(size_t i = 0; i < nthreads(); ++i) {
//DEBUG        outfile->Printf("start is %lu\n",start);
        SharedPKWrkr buf;
        if(i < nthreads() - 1) {
            buf = SharedPKWrkr(new PKWrkrInCore(primary(),sieve(),buffer_size,0,J_ints_.get(),K_ints_.get()));
        } else {
            buf = SharedPKWrkr(new PKWrkrInCore(primary(),sieve(),buffer_size,lastbuf,J_ints_.get(),K_ints_.get()));
        }
        fill_buffer(buf);
        set_ntasks(nthreads());
    }
}

void PKMgrInCore::form_PK() {
    compute_integrals();
    finalize_PK();
}

void PKMgrInCore::write() {
    get_buffer()->finalize_ints(pk_pairs());
}

void PKMgrInCore::finalize_PK() {
    for(int i = 0; i < nthreads(); ++i) {
        buffer(i).reset();
    }
}

void PKMgrInCore::prepare_JK(std::vector<SharedMatrix> D) {
    form_D_vec(D);
}

void PKMgrInCore::form_J(std::vector<SharedMatrix> J, bool exch) {

    make_J_vec(J);

    for(int N = 0; N < J.size(); ++N) {
        double* D_vec = D_glob_vecs(N);
        double* J_vec = JK_glob_vecs(N);
        double* j_ptr;
        if(exch) {
            j_ptr = K_ints_.get();
        } else {
            j_ptr = J_ints_.get();
        }
//TODO We should totally parallelize this loop now.
        for(size_t pq = 0; pq < pk_pairs(); ++pq) {
            double D_pq = D_vec[pq];
            double *D_rs = D_vec;
            double J_pq = 0.0;
            double *J_rs = J_vec;
            for(size_t rs = 0; rs <= pq; ++rs) {
//DEBUG                    if(!exch && rs == 0) {
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

    get_results(J);

}

void PKMgrInCore::finalize_JK() {
    finalize_D();
}

}  // End namespace pk
}  // End namespace psi
