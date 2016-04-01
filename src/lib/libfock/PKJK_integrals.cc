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

#include <libmints/sointegral_twobody.h>

#include<lib3index/cholesky.h>

#include <sstream>
#include "libparallel/ParallelPrinter.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#include "writers.h"

using namespace std;
using namespace psi;

namespace psi {

void PKJK::integrals_reorder() {

    outfile->Printf(" Computing reordered integrals for PK\n\n");
    int max_buckets = Process::environment.options.get_int("MAX_BUCKETS");

    PKmanager_ = boost::shared_ptr<PK_integrals>(new PK_integrals(primary_, psio_, max_buckets, memory_));

    // Get an AO integral factory
    boost::shared_ptr<IntegralFactory> intfact(new IntegralFactory(primary_));

    // Get ERI object, one per thread
    std::vector<boost::shared_ptr<TwoBodyAOInt> > tb;
    for(int i = 0; i < nthreads_; ++i) {
        tb.push_back(boost::shared_ptr<TwoBodyAOInt>(intfact->erd_eri()));
    }

    // Get all rs index for a given pq to fit in the same
    // batch
    PKmanager_->batch_sizing();

    PKmanager_->print_batches();

    // Create the buffer class, that is storing integrals in buffers in order
    // and taking care of writing to disk.
    PKmanager_->allocate_buffers();

    // Loop over buffer-filling tasks. Initially, we fill a buffer using multiple
    // threads, then write it asynchronously to disk while filling the next buffer.
    //
    // TODO (maybe): Other possibility: each thread has its own buffer, and writes it to disk
    // as the job is complete. We need a signaling mechanism for the end of the task
    // as integrals may be stored at random places in the buffer. We also need
    // to pre-stripe the PK file in this case.

    PKmanager_->open_files(false);
    size_t task_size = 0;
    for (int buf = 0; buf < PKmanager_->buf_ntasks(); ++ buf) {
        timer_on("Actual integral computation");
        // Here we need a vector of tasks to distribute over threads
        task_size = PKmanager_->task_quartets();
        outfile->Printf("The task size is %12zu\n",task_size);
#pragma omp parallel for schedule(dynamic) num_threads(nthreads_)
        for(size_t task = 0; task < task_size; ++task ) {
            int thread = 0;
#ifdef _OPENMP
            thread = omp_get_thread_num();
#endif
            short int P = PKmanager_->P(task);
            short int Q = PKmanager_->Q(task);
            short int R = PKmanager_->R(task);
            short int S = PKmanager_->S(task);
            tb[thread]->compute_shell(P,Q,R,S);
            PKmanager_->integrals_buffering(tb[thread]->buffer(), P, Q, R, S);
        }
        timer_off("Actual integral computation");

        // All shell quartets for the current task are done. We write the ordered integrals
        // to disk
        timer_on("AIO write");
        PKmanager_->write();
        timer_off("AIO write");
    }
    // We want ot deallocate buffers and wait for writing as late as possible
    PKmanager_->set_writing(true);
//    PKmanager_->deallocate_buffers();
//    PKmanager_->close_files();
}


void PKJK::integrals(){
    // We want to get the SO integrals in parallel
    // and write them to disk.
    outfile->Printf(" Computing integrals for PK\n\n");

    // Get an integral factory
    boost::shared_ptr<IntegralFactory> intfact(new IntegralFactory(primary_));

    // Now get SO basis object
    boost::shared_ptr<SOBasisSet> sobasis(new SOBasisSet(primary_,intfact));

    // Get ERI object
    std::vector<boost::shared_ptr<TwoBodyAOInt> > tb;
    for (int i = 0; i < nthreads_; ++i) {
        tb.push_back(boost::shared_ptr<TwoBodyAOInt>(intfact->eri()));
    }
    boost::shared_ptr<TwoBodySOInt> eri(new TwoBodySOInt(tb, intfact));

    // Print out some useful information
    outfile->Printf( "   Calculation information:\n");
    outfile->Printf( "      Number of atoms:                %4d\n", primary_->molecule()->natom());
    outfile->Printf( "      Number of AO shells:            %4d\n", primary_->nshell());
    outfile->Printf( "      Number of SO shells:            %4d\n", sobasis->nshell());
    outfile->Printf( "      Number of primitives:           %4d\n", primary_->nprimitive());
    outfile->Printf( "      Number of atomic orbitals:      %4d\n", primary_->nao());
    outfile->Printf( "      Number of basis functions:      %4d\n\n", primary_->nbf());
    outfile->Printf( "      Number of irreps:               %4d\n", sobasis->nirrep());
    outfile->Printf( "      Integral cutoff                 %4.2e\n", cutoff_);
    outfile->Printf( "      Number of threads:              %4d\n", nthreads_);
    outfile->Printf( "      Number of functions per irrep: [");
    for (int i=0; i<sobasis->nirrep(); ++i) {
        outfile->Printf( "%4d ", sobasis->nfunction_in_irrep(i));
    }
    outfile->Printf( "]\n\n");

    // Open IWL buffer for integral storage
    std::vector<IWLAIOWriter*> writers;
    std::vector<IWLAsync*> buffers;
    boost::shared_ptr<AIOHandler> aio(new AIOHandler(psio_));

    // Small and stupid test code
//
//    psio_address dummy;
//    int nloops = 40;
//    int num[nthreads_][nloops];
//    int num2[nthreads_][nloops];
//    int num3[nthreads_][nloops];
//    int num4[nthreads_][nloops];
//    #pragma omp parallel for schedule(dynamic) num_threads(nthreads_)
//    for (int i = 0; i < nloops; ++i) {
//       int thread = 0;
//       #ifdef _OPENMP
//         thread = omp_get_thread_num();
//       #endif
//       printf("Hi, this is thread %d\n", thread);
//       num[thread][i] = i + thread * 100;
//       num2[thread][i] = num[thread][i] + 1000;
//       num3[thread][i] = num[thread][i] + 2000;
//       num4[thread][i] = num[thread][i] + 3000;
//       #pragma omp critical(POTATO)
//       {
//       printf("Queuing job %d \n", i);
//       aio->write(33, "Dumb key", (char *)  &num[thread][i], sizeof(int), PSIO_ZERO, &dummy);
//       aio->write(33, "Dumb key", (char *) &num2[thread][i], sizeof(int), PSIO_ZERO, &dummy);
//       aio->write(33, "Dumb key", (char *) &num3[thread][i], sizeof(int), PSIO_ZERO, &dummy);
//       aio->write(33, "Dumb key", (char *) &num4[thread][i], sizeof(int), PSIO_ZERO, &dummy);
//       }
//    }
//    aio->synchronize();
//    printf("AIO is now synchronized\n");
//       
//    throw PSIEXCEPTION("Testing asynch. writing\n");
//
//
    // End of small and stupid test code

    for (int nt = 0; nt < nthreads_; ++nt) {
        buffers.push_back(new IWLAsync(psio_, aio, PSIF_SO_TEI));
        writers.push_back(new IWLAIOWriter(*buffers[nt]));
    }

    //std::vector<IWLWriter*> writers;
    //std::vector<IWL*> buffers;
    //boost::shared_ptr<AIOHandler> aio(new AIOHandler(psio_));
    //for (int nt = 0; nt < nthreads_; ++nt) {
    //    buffers.push_back(new IWL(psio_.get(), PSIF_SO_TEI, cutoff_, 0, 0));
    //    writers.push_back(new IWLWriter(*buffers[nt]));
    //}
    // Open the file
    buffers[0]->open_file(0);

    outfile->Printf("      Computing two-electron integrals...");
    //setvbuf (stdout, NULL, _IONBF, BUFSIZ);

    // Instead of using the iterator, we use a big nested loop that we can
    // parallelize with OpenMP

    #pragma omp parallel for schedule(dynamic) num_threads(nthreads_) 
    for(int i = 0; i < sobasis->nshell(); ++i) {
      int num_uniq_pk;
      int p_arr[3];
      int q_arr[3];
      int r_arr[3];
      int s_arr[3];
      int p, q, r, s;
      int thread = 0;
      #ifdef _OPENMP
        thread = omp_get_thread_num();
      #endif
      for(int j = 0; j <= i; ++j) {
        for (int k = 0; k <= j; ++k) {
          for(int l = 0; l <= k; ++l) {
            p_arr[0] = i;
            q_arr[0] = j;
            r_arr[0] = k;
            s_arr[0] = l;
            // Now we take care of all symmetry cases
            if((i == j && i == k) || (j == k && j == l)) {
                num_uniq_pk = 1;
            } else if (i == k || j == l) {
                num_uniq_pk = 2;
                p_arr[1] = i;
                q_arr[1] = k;
                r_arr[1] = j;
                s_arr[1] = l;
            } else if (j == k) {
                num_uniq_pk = 2;
                p_arr[1] = i;
                q_arr[1] = l;
                r_arr[1] = j;
                s_arr[1] = k;
            } else if (i == j || k == l) {
                num_uniq_pk = 2;
                p_arr[1] = i;
                q_arr[1] = k;
                r_arr[1] = j;
                s_arr[1] = l;
            } else {
                num_uniq_pk = 3;
                p_arr[1] = i;
                q_arr[1] = k;
                r_arr[1] = j;
                s_arr[1] = l;

                p_arr[2] = i;
                q_arr[2] = l;
                r_arr[2] = j;
                s_arr[2] = k;
            }

            for(int npk = 0; npk < num_uniq_pk; ++npk) {
              p = p_arr[npk];
              q = q_arr[npk];
              r = r_arr[npk];
              s = s_arr[npk];
              //Sort shells based on AM to save ERI some work doing permutation resorting
              if (sobasis->am(p) < sobasis->am(q)) {
                  std::swap(p,q);
              }
              if (sobasis->am(r) < sobasis->am(s)) {
                  std::swap(r,s);
              }
              if (sobasis->am(p) + sobasis->am(q) >
                      sobasis->am(r) + sobasis->am(s)) {
                  std::swap(p, r);
                  std::swap(q, s);
              }
//              printf("Computing integral <%i %i|%i %i>\n", p, q, r, s);
              eri->compute_shell(p, q, r, s, *writers[thread]);
            }


          }
        }
      }
    }

    unsigned long int total_count = 0;
    // Flush out the buffers now
    // Only one buffer should have the lastbuf_ flag set!
    for(int nt = 1; nt < nthreads_; ++nt) {
        buffers[nt]->flush(0);
        // We want to keep the file
        buffers[nt]->set_keep_flag(true);
        total_count += writers[nt]->count();
    }
    buffers[0]->flush(1);
    buffers[0]->set_keep_flag(true);
    total_count += writers[0]->count();
    aio->synchronize();
    outfile->Printf( "done\n");
    for (int nt = 0; nt < nthreads_; ++ nt) {
        outfile->Printf("Thread %d has computed %lu integrals \n", nt, writers[nt]->count());
    }
    for(int nt = 0; nt < nthreads_; ++nt) {
        delete buffers[nt];
        delete writers[nt];
    }
    psio_->close(PSIF_SO_TEI, 1);

    outfile->Printf( "      Computed %lu non-zero two-electron integrals.\n"
                     "        Stored in file %d.\n\n", total_count, PSIF_SO_TEI);
}

}
