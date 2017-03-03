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

#include "psi4/psi4-dec.h"
#include "psi4/psifiles.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libiwl/iwl.hpp"
#include "PKmanagers.h"
#include "PK_workers.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/typedefs.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/sieve.h"
#include "psi4/libqt/qt.h"
#include "psi4/libpsio/aiohandler.h"

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

std::shared_ptr<PKManager> PKManager::build_PKManager(std::shared_ptr<PSIO> psio,
                  std::shared_ptr<BasisSet> primary, size_t memory, Options &options,
                  bool dowK, double omega_in) {

    std::string algo = options.get_str("PK_ALGO");
    bool noincore = options.get_bool("PK_NO_INCORE");

    // We introduce another safety factor in the memory, otherwise
    // we are apparently prone to being killed by the OS.
    //TODO: Check for memory leaks ? Trace memory usage ?
    memory = memory * 9 / 10;

    // Approximate number of batches beyond which the Yoshimine algorithm should be preferred
    // Estimated from a single example on a nucleic base, there may be a better number
    int algo_factor = 40;

    size_t nbf = primary->nbf();
    size_t pk_size = nbf * (nbf + 1) / 2;
    pk_size = pk_size * (pk_size + 1) / 2;

    int ncorebuf = 2;
    if(dowK) {
        ncorebuf = 3;
    }

    bool do_reord = false;
    bool do_yosh = false;
    bool do_incore = false;
    if(options["PK_ALGO"].has_changed()) {
      if(algo == "REORDER") {
        do_reord = true;
      } else if (algo == "YOSHIMINE") {
        do_yosh = true;
      }
    } else {
      if ( algo_factor * memory > pk_size) {
        do_reord = true;
      } else {
        do_yosh = true;
      }
    }

    if(ncorebuf * pk_size < memory && !noincore) do_incore = true;

    std::shared_ptr<PKManager> pkmgr;

    if(do_incore) {
        outfile->Printf("  Using in-core PK algorithm.\n");
        pkmgr = std::shared_ptr<PKManager>(new PKMgrInCore(primary,memory,options));
    // Estimate that we'll need less than 40 buffers: do integral reorder
    } else if (do_reord) {
        outfile->Printf("  Using integral reordering PK algorithm.\n");
        pkmgr = std::shared_ptr<PKManager>(new PKMgrReorder(psio,primary,memory,options));
    // Low memory case: Yoshimine should be faster
    } else if (do_yosh) {
        outfile->Printf("  Using Yoshimine PK algorithm.\n");
        pkmgr = std::shared_ptr<PKManager>(new PKMgrYoshimine(psio,primary,memory,options));
    } else {
      throw PSIEXCEPTION("PK algorithm selection error.\n");
    }

    // Set do_wK and value of omega if necessary
    // Might need to get these functions public
    pkmgr->set_wK(dowK);
    pkmgr->set_omega(omega_in);

    return pkmgr;

}

PKManager::PKManager(std::shared_ptr<BasisSet> primary, size_t memory, Options& options) :
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
    nthreads_ = Process::environment.get_n_threads();
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
void PKManager::compute_integrals(bool wK) {

    // Get an AO integral factory
    std::shared_ptr<IntegralFactory> intfact(new IntegralFactory(primary_));

    // Get ERI object, one per thread
    std::vector<std::shared_ptr<TwoBodyAOInt> > tb;

    if(wK) {
        for(int i = 0; i < nthreads_; ++i) {
            tb.push_back(std::shared_ptr<TwoBodyAOInt>(intfact->erf_eri(omega())));
        }
    } else {
        for(int i = 0; i < nthreads_; ++i) {
            tb.push_back(std::shared_ptr<TwoBodyAOInt>(intfact->erd_eri()));
        }
    }


    size_t nshqu = 0;
#pragma omp parallel for num_threads(nthreads_) schedule(dynamic) reduction(+:nshqu)
    for(size_t i = 0; i < ntasks_; ++i) {
        // We need to get the list of shell quartets for each task
        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif
        SharedPKWrkr buf = get_buffer();
//DEBUG        outfile->Printf("Starting task %d\n", i);
        if(!wK) {  // Computing usual integrals
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
//DEBUG#    pragma omp critical
//DEBUG            outfile->Printf("Computing shell <%d %d|%d %d>\n",P,Q,R,S);
                tb[thread]->compute_shell(P,Q,R,S);
                integrals_buffering(tb[thread]->buffer(),P,Q,R,S);
//DEBUG#pragma omp critical
//DEBUG              {
//DEBUG                outfile->Printf("After buffering\n");
//DEBUG                debug_wrt();
//DEBUG              }
                ++nshqu;
            }
        } else {  // Computing range-separated integrals
            buf->set_do_wK(true);
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
//DEBUG#    pragma omp critical
//DEBUG            outfile->Printf("Computing shell <%d %d|%d %d>\n",P,Q,R,S);
                tb[thread]->compute_shell(P,Q,R,S);
                integrals_buffering_wK(tb[thread]->buffer(),P,Q,R,S);
                ++nshqu;
            }

        }
        if(!wK) {
            write();
        } else {
            write_wK();
        }
    }
    size_t nsh = primary_->nshell();
    size_t nsh_u = nsh * (nsh + 1) / 2;
    nsh_u = nsh_u * (nsh_u + 1) / 2;
    if(wK) {
        outfile->Printf("  We computed %lu wK shell quartets total.\n",nshqu);
        outfile->Printf("  Whereas there are %lu wK unique shell quartets.\n",nsh_u);
    } else {
        outfile->Printf("  We computed %lu shell quartets total.\n",nshqu);
        outfile->Printf("  Whereas there are %lu unique shell quartets.\n",nsh_u);
    }
    if (nshqu > nsh_u) {
        outfile->Printf("  %7.2f percent of shell quartets recomputed by reordering.\n", (nshqu - nsh_u) / float(nsh_u) * 100);
    }

}

void PKManager::compute_integrals_wK() {
    // Call the same routine with additional branching for each
    // shell. Bit less efficient than re-writing the routine but more maintainable
    compute_integrals(true);
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

void PKManager::integrals_buffering_wK(const double *buffer, unsigned int P, unsigned int Q,
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
            iobuffers_[thread]->fill_values_wK(val, i, j, k, l);
        }
    }

}

// use the vectors of Cleft and Cright to determine which of these densities,
// if any, are symmetric matrices.
void PKManager::form_D_vec(std::vector<SharedMatrix> D, std::vector<SharedMatrix> Cl, std::vector<SharedMatrix> Cr) {

    D_ = D;
    all_sym_ = true;
    symmetric_.clear();
    // Check which density matrices are asymmetric, if any
    for(int N = 0; N < D.size(); ++N) {
        symmetric_.push_back(Cl[N] == Cr[N]);
        all_sym_ = all_sym_ && (Cl[N] == Cr[N]);
    }

    // For debugging
    if(options_.get_bool("PK_ALL_NONSYM")) {
        all_sym_ = false;
        for(int N = 0; N < D.size(); ++N) {
            symmetric_[N] = false;
        }
        outfile->Printf("  All matrices considered asymmetric.\n");
    }

    for (int N = 0; N < D.size(); ++N) {
        // Symmetric density matrix: store triangular form
        if(is_sym(N)) {
            double* Dvec = new double[pk_pairs_];
// All elements are filled, don't need memset
            //::memset((void *)Dvec,'\0',pk_pairs_ * sizeof(double));
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
        // Non-symmetric density matrix: store the whole thing with
        // appropriate factors
        } else {
            double* Dvec = new double[nbf_ * nbf_];
            D_vec_.push_back(Dvec);
            size_t pqval = 0;
            for (int p = 0; p < nbf_; ++p) {
                for(int q = 0; q < nbf_; ++q) {
                    if(p != q) {
                        Dvec[pqval] = D[N]->get(0,p,q);
                    } else {
                        Dvec[pqval] = 0.5 * D[N]->get(0,p,q);
                    }
                    ++pqval;
                }
            }

        }
    }

}

void PKManager::make_J_vec(std::vector<SharedMatrix> J) {
    for(int N = 0; N < J.size(); ++N) {
        // Symmetric density matrix: stores triangular J
        if(is_sym(N)) {
            double* Jvec = new double[pk_pairs_];
            ::memset((void*)Jvec,'\0',pk_pairs_ * sizeof(double));
            JK_vec_.push_back(Jvec);
        // Non-symmetric density matrix: stores NULL pointer
        // as we will use full J.
        // Shared pointer initialized to null
        //TODO Could actually store triangle since J is always symmetric
        } else {
            JK_vec_.push_back(nullptr);
        }
    }
}

void PKManager::get_results(std::vector<SharedMatrix> J,std::string exch) {
    // Now, directly transfer data to resulting matrices
    for(int N = 0; N < J.size(); ++N) {
        // Symmetric density matrix: copy triangle and build full matrix
        if(is_sym(N) && exch != "wK") {
            double *Jp = JK_vec_[N];
            for(int p = 0; p < nbf_; ++p) {
                for(int q = 0; q <= p; ++q) {
                    J[N]->set(0,p,q,*Jp++);
                }
            }
            J[N]->copy_lower_to_upper();
            delete [] JK_vec_[N];
        // Non-symmetric density matrix: result is already in the target
        //TODO: if we store the triangle only, change that.
        } else if (exch == "") {
          double** Jp = J[N]->pointer();
          for(int p = 0; p < nbf_; ++p) {
            Jp[p][p] *= 0.5;
          }
        }
    }
    JK_vec_.clear();

}

void PKManager::form_K(std::vector<SharedMatrix> K) {
    // Right now, this supports both J and K. K asym is
    // formed at the same time than J asym for convenience.
    // The following call only forms K sym.
    form_J(K,"K");
}

void PKManager::form_wK(std::vector<SharedMatrix> wK) {
    // We call form_J with an appropriate flag. This way,
    // we can reuse the code for K. Not optimal but should work.
    form_J(wK,"wK");
}

void PKManager::finalize_D() {
    for(int N = 0; N < D_vec_.size(); ++N) {
        delete [] D_vec_[N];
    }
    D_vec_.clear();
}

PKMgrDisk::PKMgrDisk(std::shared_ptr<PSIO> psio, std::shared_ptr<BasisSet> primary,
          size_t memory, Options &options) : PKManager(primary,memory,options) {
    psio_ = psio;
    AIO_ = std::shared_ptr<AIOHandler>(new AIOHandler(psio_));
    max_batches_ = options.get_int("PK_MAX_BUCKETS");
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

void PKMgrDisk::initialize_wK() {
    prestripe_files_wK();
    print_batches_wK();
    allocate_buffers_wK();
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

    outfile->Printf("  Sizing the integral batches needed.\n");
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
            // The last batch disappeared thus we need to correct batch_for_pq
            for(size_t i = 0; i < batch_for_pq_.size(); ++i) {
              if(batch_for_pq_[i] == lastb) batch_for_pq_[i]--;
            }
        }
    }

    int nbatches = batch_pq_min_.size();
    if (nbatches > max_batches_) {
      outfile->Printf( "  PKJK: maximum number of batches exceeded\n") ;
      outfile->Printf( "  PK computation needs %d batches, max. number: %d\n", nbatches, max_batches_);
      throw PSIEXCEPTION("  PK Failure: max batches exceeded\n");

    }

    // Finally, we want to store the pairs of indices p,q corresponding
    // to min pq and max pq for each batch (useful for non-sym. density matrices)
    // We go beyond the number of pq to make sure we include the last one.
    outfile->Printf("  Building batch lookup table.\n");
    ind_for_pq_[0] = std::make_pair(0,0);
    int nb = 0;
    for(int p = 0; p <= nbf(); ++p) {
        for(int q = 0; q <= p; ++q) {
            pq = INDEX2(p,q);
            if(batch_pq_max_[nb] == pq) {
                ind_for_pq_[pq] = std::make_pair(p,q);
                ++nb;
            }
        }
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

void PKMgrDisk::write_wK() {
    get_buffer()->write_wK(batch_index_min_,batch_index_max_,pk_pairs());
}

void PKMgrDisk::open_PK_file() {
    psio_->open(pk_file_, PSIO_OPEN_OLD);
}

void PKMgrDisk::close_PK_file(bool keep) {
    psio_->close(pk_file_,keep);
}

void PKMgrDisk::prepare_JK(std::vector<SharedMatrix> D, std::vector<SharedMatrix> Cl,
                           std::vector<SharedMatrix> Cr) {
    if (writing()) {
        finalize_PK();
        set_writing(false);
    } else {
        open_PK_file();
    }
    form_D_vec(D,Cl,Cr);
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

void PKMgrDisk::form_J(std::vector<SharedMatrix> J, std::string exch,
                       std::vector<SharedMatrix> K) {
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
        if (exch == "K") {
            label = PKWorker::get_label_K(batch);
        } else if(exch == "wK") {
            label = PKWorker::get_label_wK(batch);
        } else {
            label = PKWorker::get_label_J(batch);
        }
        psio_->read_entry(pk_file_, label, (char *) j_block, batch_size * sizeof(double));

        // Read one entry, use it for all density matrices
        for(int N = 0; N < J.size(); ++N) {
            // Symmetric density matrix, pure triangular
            if(is_sym(N) && exch != "wK") {
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
//DEBUG                        if(!exch && rs == 0) {
//DEBUG                          outfile->Printf("PK int (%lu|%lu) = %20.16f\n",pq,rs,*j_ptr);
//DEBUG                        }
                        J_pq += *j_ptr * (*D_rs);
                        *J_rs += *j_ptr * D_pq;
                        ++D_rs;
                        ++J_rs;
                        ++j_ptr;
                    }
                    J_vec[pq] += J_pq;
                }
            // Non-symmetric density matrix case
            } else if (exch == "" || exch == "wK") {
                double* j_ptr = j_block;
                const int fp = ind_for_pq_[min_pq].first;
                const int fq = ind_for_pq_[min_pq].second;
                const int maxp = ind_for_pq_[max_pq].first;
                if(exch != "wK") {
                    double* D_vec = D_glob_vecs(N);
                    double** J_vec = J[N]->pointer();
                    for(int p = fp; p <= maxp; ++p) {
                      int maxq = (p == maxp) ? ind_for_pq_[max_pq].second : p + 1;
                      int poffs = p * nbf();
                      int q = (p == fp) ? fq : 0;
                      for(;q < maxq; ++q) {
                        int qoffs = q * nbf();
                        for(int r = 0; r <= p; ++r) {
                          int roffs = r * nbf();
                          int maxs = (r == p) ? q : r;
                          for(int s = 0; s <= maxs; ++s) {
//                            if(!exch && INDEX2(r,s) == 0) {
//                              outfile->Printf("PK int (%lu|%lu) = %20.16f\n",INDEX2(p,q),INDEX2(r,s),*j_ptr);
//                              int x = 0;  // For the lolz
//                            }
                              J_vec[p][q] += *j_ptr * (D_vec[roffs + s] + D_vec[s * nbf() + r]);
                              J_vec[q][p] += *j_ptr * (D_vec[roffs + s] + D_vec[s * nbf() + r]);
                              J_vec[r][s] += *j_ptr * (D_vec[poffs + q] + D_vec[qoffs + p]);
                              J_vec[s][r] += *j_ptr * (D_vec[poffs + q] + D_vec[qoffs + p]);
                              ++j_ptr;
                          }

                        }
                      }
                    }
                }
                // Since we just read a batch, might as well compute K
                // Primitive algorithm, just contract integrals with appropriate
                // element on the fly. Might be faster than reading/writing the appropriate
                // PK supermatrix
                if(K.size() || exch == "wK") {
                    double** Dmat = original_D(N)->pointer();
                    double** K_vec;
                    if(exch == "wK") {
                        K_vec = J[N]->pointer();
                    } else {
                        K_vec = K[N]->pointer();
                    }
                    j_ptr = j_block;
                    for(int p = fp; p <= maxp; ++p) {
                      int maxq = (p == maxp) ? ind_for_pq_[max_pq].second : p + 1;
               //       int poffs = p * nbf();
                      int q = (p == fp) ? fq : 0;
                      for(;q < maxq; ++q) {
                 //       int qoffs = q * nbf();
                        for(int r = 0; r <= p; ++r) {
                   //       int roffs = r * nbf();
                          int maxs = (r == p) ? q : r;
                          for(int s = 0; s <= maxs; ++s) {
                          // Need ugly factors for now. A better solution would be great.
                            double fac = 1.0;
                     //       int soffs = s * nbf();
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
                            K_vec[p][r] += (*j_ptr) * fac * Dmat[q][s];
                            K_vec[r][p] += (*j_ptr) * fac * Dmat[s][q];
                            K_vec[q][r] += (*j_ptr) * fac * Dmat[p][s];
                            K_vec[p][s] += (*j_ptr) * fac * Dmat[q][r];
                            K_vec[s][p] += (*j_ptr) * fac * Dmat[r][q];
                            K_vec[r][q] += (*j_ptr) * fac * Dmat[s][p];
                            K_vec[s][q] += (*j_ptr) * fac * Dmat[r][p];
                            K_vec[q][s] += (*j_ptr) * fac * Dmat[p][r];
                            ++j_ptr;
                          }
                        }
                      }
                    }

                }
            }  // end of non-symmetric case
        }  // End of loop over J matrices

        delete [] label;
        delete [] j_block;
    }  // End of batch loop
    get_results(J,exch);
}

void PKMgrDisk::finalize_JK() {
    finalize_D();
    close_PK_file(true);
}

PKMgrReorder::PKMgrReorder(std::shared_ptr<PSIO> psio, std::shared_ptr<BasisSet> primary,
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

// Pre-striping PK file for wK
void PKMgrReorder::prestripe_files_wK() {
    // PK file should still be open at this point for AIO writing of
    // the J/K supermatrices
    for(int batch = 0; batch < batch_ind_min().size(); ++batch) {
        size_t batch_size = batch_ind_max()[batch] - batch_ind_min()[batch];
        // Keep the labels around in a vector
        label_wK_.push_back(PKWorker::get_label_wK(batch));
        AIO()->zero_disk(pk_file(),label_wK_[batch],1,batch_size);
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

void PKMgrReorder::allocate_buffers_wK() {
    // We need to redimension and reset the number of tasks
    size_t mem_per_thread = memory() / nthreads();
    size_t buf_size = mem_per_thread / 2; // For AIO
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
    outfile->Printf("  wK Task number: %lu\n",ntasks());
    outfile->Printf("  wK Buffer size: %lu\n",buf_size);
    outfile->Printf("  wK Buffer per thread: %lu\n",buf_per_thread);

    for(int i = 0; i < nthreads(); ++i) {
        buffer(i)->allocate_wK(buf_size,buf_per_thread);
    }
}

void PKMgrReorder::form_PK() {
    compute_integrals();
    set_writing(true);
}

void PKMgrReorder::form_PK_wK() {
    compute_integrals_wK();
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

PKMgrYoshimine::PKMgrYoshimine(std::shared_ptr<PSIO> psio, std::shared_ptr<BasisSet> primary,
                               size_t memory, Options &options) :
    PKMgrDisk(psio,primary,memory,options) {
    iwl_file_J_ = PSIF_SO_PKSUPER1;
    iwl_file_K_ = PSIF_SO_PKSUPER2;
    iwl_file_wK_ = PSIF_WK_PK;
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
    // We need to check whether iwlsize bytes can actually fit
    // in memory.
    size_t safemem = memory() * 9 / 10;
    size_t nrows = iwlsize / safemem;
    size_t leftover = iwlsize % safemem;
    safemem = std::min(iwlsize, safemem);
    if (nrows > 0 ) {
        AIO()->zero_disk(iwl_file_J_,IWL_KEY_BUF,nrows,safemem);
    }
    AIO()->zero_disk(iwl_file_J_,IWL_KEY_BUF,1,leftover);

    // And we do the same for the file containing K buckets
    // We need at most twice as much IWL buffers since integrals can be
    // pre-sorted in two different buckets for K
    psio()->open(iwl_file_K_, PSIO_OPEN_NEW);
    if (nrows > 0) {
        AIO()->zero_disk(iwl_file_K_, IWL_KEY_BUF, 2*nrows, safemem);
    }
    AIO()->zero_disk(iwl_file_K_, IWL_KEY_BUF, 2, leftover);
}

void PKMgrYoshimine::prestripe_files_wK() {
    psio()->open(iwl_file_wK_, PSIO_OPEN_NEW);
    // Pre-stripe the new file
    // Number of IWL buffers necessary to write all integrals
    size_t num_iwlbuf = pk_size() / ints_per_buf_ + 1;
    // The last buffer of each batch may be partially filled only
    num_iwlbuf += batch_ind_min().size();
    size_t iwlsize_bytes = num_iwlbuf * iwl_int_size_;
    size_t iwlsize = iwlsize_bytes / sizeof(double) + 1;
    AIO()->zero_disk(iwl_file_wK_,IWL_KEY_BUF,1,iwlsize);
}

void PKMgrYoshimine::allocate_buffers() {
    int buf_per_thread = batch_ind_min().size();
    // Factor 2 because we need an address for J and an address for K
    // J and K addresses for the same buckets are stored contiguously, i.e.
    // element [0] is J address for first bucket and element [1] the
    // K address for the first bucket
    // we actually need the boost shared_ptr, the std::shared_ptr does
    // not support array syntax
    std::shared_ptr<std::vector<size_t>> current_pos(new std::vector<size_t>(2 * buf_per_thread));

    (*current_pos)[0] = 0;
    (*current_pos)[1] = 0;
    //Current position of the bucket in the IWL file
    // For K, each integral can potentially go into two buckets, so we
    // give each bucket twice the J bucket's size. It's wasting quite a bit
    // of space, other solution is to explicitly count integrals.
    for(int i = 1; i < buf_per_thread; ++i) {
        size_t batchsize = batch_ind_max()[i - 1] - batch_ind_min()[i - 1];
        size_t iwlperbatch = batchsize / ints_per_buf_ + 1;
        (*current_pos)[2 * i] = iwlperbatch * iwl_int_size_ + (*current_pos)[2 * i - 2];
        (*current_pos)[2 * i + 1] = 2 * (iwlperbatch * iwl_int_size_) + (*current_pos)[2 * i - 1];
    }

    // Ok, now we have the size of a buffer and how many buffers
    // we want for each thread. We can allocate IO buffers.
    for(int i = 0; i < nthreads(); ++i) {
        fill_buffer(SharedPKWrkr(new PKWrkrIWL(primary(),sieve(),AIO(),iwl_file_J_,iwl_file_K_,ints_per_buf_,batch_for_pq(),current_pos)));
    }

}

void PKMgrYoshimine::allocate_buffers_wK() {
    // Need a vector of starting positions for each batch in the wK
    // IWL file
    int bufperthread = batch_ind_min().size();
    std::shared_ptr<std::vector<size_t>> current_pos(new std::vector<size_t>(bufperthread));
    (*current_pos)[0] = 0;
    for(int i = 1; i < bufperthread; ++i) {
        size_t batchsize = batch_ind_max()[i - 1] - batch_ind_min()[i - 1];
        size_t iwlperbatch = batchsize / ints_per_buf_ + 1;
        (*current_pos)[i] = iwlperbatch * iwl_int_size_ + (*current_pos)[i - 1];
    }

    // Now we pass these new positions to the Workers, which will also
    // take care of their setup for wK.
    for(int i = 0; i < nthreads(); ++i) {
        buffer(i)->allocate_wK(current_pos,iwl_file_wK_);
    }

}

void PKMgrYoshimine::form_PK() {
    compute_integrals();

    // Now we need to actually read the integral file and sort
    // the buckets to write them to the PK file
    // We deallocate buffers and synchronize AIO writing in there.

    sort_ints();

}

void PKMgrYoshimine::form_PK_wK() {
    compute_integrals_wK();

    // Now we need to actually read the integral file and sort
    // the buckets to write them to the PK file
    // We deallocate buffers and synchronize AIO writing in there.

    sort_ints_wK();
}

void PKMgrYoshimine::compute_integrals(bool wK) {

    // Get an AO integral factory
    std::shared_ptr<IntegralFactory> intfact(new IntegralFactory(primary()));

    // Get ERI object, one per thread
    std::vector<std::shared_ptr<TwoBodyAOInt> > tb;
    if(!wK) {
        for(int i = 0; i < nthreads(); ++i) {
            tb.push_back(std::shared_ptr<TwoBodyAOInt>(intfact->erd_eri()));
        }
    } else {
        for(int i = 0; i < nthreads(); ++i) {
            tb.push_back(std::shared_ptr<TwoBodyAOInt>(intfact->erf_eri(omega())));
        }
    }

    // Loop over significant shell pairs from ERISieve
    const std::vector< std::pair<int, int> >& sh_pairs = sieve()->shell_pairs();
    size_t npairs = sh_pairs.size();
    // We avoid having one more branch in the loop by moving it outside
    if(!wK) {
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
//DEBUG                    outfile->Printf("Computing shell <%i %i|%i %i>\n", P, Q, R, S);
                    tb[thread]->compute_shell(P, Q, R, S);
                    integrals_buffering(tb[thread]->buffer(),P,Q,R,S);
                }
            }
        } // end of parallelized loop

        // We write all remaining buffers to disk.
        write();
    } else {
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
//DEBUG                    outfile->Printf("Computing shell <%i %i|%i %i>\n", P, Q, R, S);
                    tb[thread]->compute_shell(P, Q, R, S);
                    integrals_buffering_wK(tb[thread]->buffer(),P,Q,R,S);
                }
            }
        } // end of parallelized loop

        // We write all remaining buffers to disk.
        write_wK();

    }

}

void PKMgrYoshimine::compute_integrals_wK() {
    // Simply invoke the compute_integrals routine with a flag
    compute_integrals(true);
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

void PKMgrYoshimine::write_wK() {
    // Each thread buffer has several buckets that may be partially full.
    // Concatenate all internal buckets to the ones of thread buffer 0.
    double val;
    size_t i, j, k, l;
    SharedPKWrkr buf0 = buffer(0);
    SharedPKWrkr buftarget;
    for(int t = 1; t < nthreads(); ++t) {
        buftarget = buffer(t);
        unsigned int nbufs = buftarget->nbuf();
        // Get through all wK buffers
        for(int buf = 0; buf < nbufs; ++buf) {
            while(buftarget->pop_value_wK(buf,val,i,j,k,l)) {
                buf0->insert_value_wK(buf,val,i,j,k,l);
            }
        }

    }
    // By now we should have transferred all remaining integrals to
    // buffer 0, we now flush it.
    buf0->flush_wK();

}

void PKMgrYoshimine::sort_ints(bool wK) {
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
    // wK is done after J/K so we open the file as old
    psio()->open(pk_file(), wK ? PSIO_OPEN_OLD : PSIO_OPEN_NEW);
    //TODO: this function may need renaming
    finalize_PK();
    set_writing(false);
    if(!wK) {
        close_iwl_buckets();
        generate_J_PK(twoel_ints,max_size);
        // Need to reset to zero the two-el integral array
        ::memset((void*)twoel_ints,'\0',max_size * sizeof(double));
        generate_K_PK(twoel_ints,max_size);

    } else {
        close_iwl_buckets_wK();
        generate_wK_PK(twoel_ints,max_size);
    }

    //delete two-el int array
    delete [] twoel_ints;

    psio()->close(pk_file(),1);

}

void PKMgrYoshimine::sort_ints_wK() {
    // We call sort_ints with a flag
    sort_ints(true);
}

void PKMgrYoshimine::close_iwl_buckets() {
    psio()->close(iwl_file_J_, 1);
    psio()->close(iwl_file_K_, 1);
}

void PKMgrYoshimine::close_iwl_buckets_wK() {
    psio()->close(iwl_file_wK_, 1);
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

void PKMgrYoshimine::generate_wK_PK(double *twoel_ints, size_t max_size) {

    IWL inbuf(psio().get(),iwl_file_wK_,0.0, 1, 0);

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
            char* label = PKWorker::get_label_wK(batch);
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

void PKMgrInCore::initialize() {
    print_batches();
    allocate_buffers();
}

void PKMgrInCore::initialize_wK() {
   /// Nothing to do
   print_batches_wK();
}

void PKMgrInCore::print_batches() {
    PKManager::print_batches();
    outfile->Printf("  Performing in-core PK\n");
    int nbufincore = do_wk() ? 3 : 2;
    outfile->Printf("  Using %lu doubles for integral storage.\n", nbufincore * pk_size());
}

void PKMgrInCore::allocate_buffers() {
    // Need to allocate two big arrays
    J_ints_ = std::unique_ptr< double [] >(new double[pk_size()]);
    K_ints_ = std::unique_ptr< double [] >(new double[pk_size()]);
    ::memset((void*) J_ints_.get(), '\0', pk_size() * sizeof(double));
    ::memset((void*) K_ints_.get(), '\0', pk_size() * sizeof(double));
    if(do_wk()) {
        wK_ints_ = std::unique_ptr< double [] >(new double[pk_size()]);
        ::memset((void*) wK_ints_.get(), '\0', pk_size() * sizeof(double));
    }

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
        buf = SharedPKWrkr(new PKWrkrInCore(primary(),sieve(),
              buffer_size,lastbuf,J_ints_.get(),K_ints_.get(),
              wK_ints_.get(), nthreads()));
        fill_buffer(buf);
        set_ntasks(nthreads());
    }
}

void PKMgrInCore::form_PK() {
    compute_integrals();
    if(!do_wk()) {
        finalize_PK();
    }
}

void PKMgrInCore::form_PK_wK() {
    compute_integrals_wK();
    finalize_PK();
}

void PKMgrInCore::write() {
    get_buffer()->finalize_ints(pk_pairs());
}

void PKMgrInCore::write_wK() {
    get_buffer()->finalize_ints_wK(pk_pairs());
}

void PKMgrInCore::finalize_PK() {
    for(int i = 0; i < nthreads(); ++i) {
        buffer(i).reset();
    }
}

void PKMgrInCore::prepare_JK(std::vector<SharedMatrix> D, std::vector<SharedMatrix> Cl,
                             std::vector<SharedMatrix> Cr) {
    form_D_vec(D,Cl,Cr);
}

void PKMgrInCore::form_J(std::vector<SharedMatrix> J, std::string exch, std::vector<SharedMatrix> K) {

    make_J_vec(J);

    for(int N = 0; N < J.size(); ++N) {
        double* j_ptr;
        if(exch == "K") {
            j_ptr = K_ints_.get();
        } else {
            j_ptr = J_ints_.get();
        }
//TODO We should totally parallelize this loop now.
        // Symmetric density matrix case
        if(is_sym(N) && exch != "wK") {
            double* J_vec = JK_glob_vecs(N);
            double* D_vec = D_glob_vecs(N);
            for(size_t pq = 0; pq < pk_pairs(); ++pq) {
                double D_pq = D_vec[pq];
                double *D_rs = D_vec;
                double J_pq = 0.0;
                double *J_rs = J_vec;
                for(size_t rs = 0; rs <= pq; ++rs) {
        //DEBUG               if(!exch && rs == 0) {
        //DEBUG                 outfile->Printf("PK int (%lu|%lu) = %20.16f\n",pq,rs,*j_ptr);
        //DEBUG               }
                    J_pq += *j_ptr * (*D_rs);
                    *J_rs += *j_ptr * D_pq;
                    ++D_rs;
                    ++J_rs;
                    ++j_ptr;
                }
                J_vec[pq] += J_pq;
            }

        //TODO ? Fuse J and K loops ?
        // Non-symmetric density matrix
        } else if (exch == "" || exch == "wK") {
            if(exch == "") {
                double* D_vec = D_glob_vecs(N);
                double** J_vec = J[N]->pointer();
                for(int p = 0;p < nbf(); ++p) {
                  int poffs = p * nbf();
                  for(int q = 0;q <= p; ++q) {
                    int qoffs = q * nbf();
                    for(int r = 0; r <= p; ++r) {
                      int roffs = r * nbf();
                      int maxs = (r == p) ? q : r;
                      for(int s = 0; s <= maxs; ++s) {
                          J_vec[p][q] += *j_ptr * (D_vec[roffs + s] + D_vec[s * nbf() + r]);
                          J_vec[q][p] += *j_ptr * (D_vec[roffs + s] + D_vec[s * nbf() + r]);
                          J_vec[r][s] += *j_ptr * (D_vec[poffs + q] + D_vec[qoffs + p]);
                          J_vec[s][r] += *j_ptr * (D_vec[poffs + q] + D_vec[qoffs + p]);
                          ++j_ptr;
                      }

                    }
                  }
                }
            }
            // Since we just read a batch, might as well compute K
            // Primitive algorithm, just contract integrals with appropriate
            // element on the fly. Might be faster than reading/writing the appropriate
            // PK supermatrix
            if(K.size() || exch =="wK") {
                double** Dmat = original_D(N)->pointer();
                double** K_vec;
                if(exch == "wK") {
                    K_vec = J[N]->pointer();
                    j_ptr = wK_ints_.get();
                } else {
                    K_vec = K[N]->pointer();
                    // We use J supermatrix because it contains every unique integral
                    // K supermatrix has summed some integrals that we need separately
                    j_ptr = J_ints_.get();
                }
                for(int p = 0;p < nbf(); ++p) {
          //        int poffs = p * nbf();
                  for(int q = 0;q <= p; ++q) {
          //          int qoffs = q * nbf();
                    for(int r = 0; r <= p; ++r) {
          //            int roffs = r * nbf();
                      int maxs = (r == p) ? q : r;
                      for(int s = 0; s <= maxs; ++s) {
                      // Need ugly factors for now. A better solution would be great.
                        double fac = 1.0;
          //              int soffs = s * nbf();
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
//DEBUG                       bool prt = (p == 0 && r == 0) || (q == 0 && r == 0) || (p == 0 && s == 0) || (s == 0 && q == 0);
//DEBUG                       if(exch == "wK" && prt) {
//DEBUG                         outfile->Printf("PK int (%d %d|%d %d) = %20.16f\n",p,q,r,s,*j_ptr);
//DEBUG                       }
                        K_vec[p][r] += (*j_ptr) * fac * Dmat[q][s];
                        K_vec[r][p] += (*j_ptr) * fac * Dmat[s][q];
                        K_vec[q][r] += (*j_ptr) * fac * Dmat[p][s];
                        K_vec[p][s] += (*j_ptr) * fac * Dmat[q][r];
                        K_vec[s][p] += (*j_ptr) * fac * Dmat[r][q];
                        K_vec[r][q] += (*j_ptr) * fac * Dmat[s][p];
                        K_vec[s][q] += (*j_ptr) * fac * Dmat[r][p];
                        K_vec[q][s] += (*j_ptr) * fac * Dmat[p][r];
                        ++j_ptr;
                        }
                    }
                  }
                }

            }

        } // End of non-symmetric condition
    } // End of loop over J/K matrices

    get_results(J,exch);

}

void PKMgrInCore::finalize_JK() {
    finalize_D();
}

}  // End namespace pk
}  // End namespace psi
