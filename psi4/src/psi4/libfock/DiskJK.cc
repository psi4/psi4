/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
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

#include "psi4/lib3index/3index.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/aiohandler.h"
#include "psi4/libqt/qt.h"
#include "psi4/psi4-dec.h"
#include "psi4/psifiles.h"
#include "psi4/libmints/sieve.h"
#include "psi4/libiwl/iwl.hpp"
#include "jk.h"
#include "jk_independent.h"
#include "link.h"
#include "direct_screening.h"
#include "cubature.h"
#include "points.h"

#include "psi4/lib3index/cholesky.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/sobasis.h"
#include "psi4/libmints/mintshelper.h"
#include <sstream>
#include "psi4/libparallel/ParallelPrinter.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace psi;

namespace psi {


DiskJK::DiskJK(std::shared_ptr<BasisSet> primary, Options& options) :
   JK(primary), options_(options)
{
    common_init();
}
DiskJK::~DiskJK()
{
}
void DiskJK::common_init()
{
}
void DiskJK::print_header() const
{
    if (print_) {
        outfile->Printf( "  ==> DiskJK: Disk-Based J/K Matrices <==\n\n");

        outfile->Printf( "    J tasked:          %11s\n", (do_J_ ? "Yes" : "No"));
        outfile->Printf( "    K tasked:          %11s\n", (do_K_ ? "Yes" : "No"));
        outfile->Printf( "    wK tasked:         %11s\n", (do_wK_ ? "Yes" : "No"));
        outfile->Printf( "    Memory (MB):       %11ld\n", (memory_ *8L) / (1024L * 1024L));
        if (do_wK_)
            outfile->Printf( "    Omega:             %11.3E\n", omega_);
        outfile->Printf( "    Schwarz Cutoff:    %11.0E\n\n", cutoff_);
    }
}
void DiskJK::preiterations()
{
    std::shared_ptr<MintsHelper> mints(new MintsHelper(primary_, options_, 0));
    mints->integrals();
    if(do_wK_)
        mints->integrals_erf(omega_);


    std::shared_ptr<SOBasisSet> bas = mints->sobasisset();

    so2symblk_ = new int[primary_->nbf()];
    so2index_  = new int[primary_->nbf()];
    size_t so_count = 0;
    size_t offset = 0;
    for (int h = 0; h < bas->nirrep(); ++h) {
        for (int i = 0; i < bas->dimension()[h]; ++i) {
            so2symblk_[so_count] = h;
            so2index_[so_count] = so_count-offset;
            ++so_count;
        }
        offset += bas->dimension()[h];
    }
    mints.reset();
}
void DiskJK::compute_JK()
{
    std::shared_ptr<PSIO> psio(new PSIO());
    IWL *iwl = new IWL(psio.get(), PSIF_SO_TEI, cutoff_, 1, 1);
    Label *lblptr = iwl->labels();
    Value *valptr = iwl->values();
    int labelIndex, pabs, qabs, rabs, sabs, prel, qrel, rrel, srel, psym, qsym, rsym, ssym;
    double value;
    bool lastBuffer;
    if(J_.size() == K_.size()){
        do{
            lastBuffer = iwl->last_buffer();
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

                int pqsym = psym ^ qsym;
                int rssym = rsym ^ ssym;
                int qrsym = qsym ^ rsym;
                int pssym = psym ^ ssym;
                int prsym = psym ^ rsym;
                int qssym = qsym ^ ssym;


                for (size_t N = 0; N < J_.size(); N++) {

                    Matrix* J = J_[N].get();
                    Matrix* K = K_[N].get();
                    Matrix* D = D_[N].get();

                    int sym = J->symmetry();

                    /* (pq|rs) */
                    if(pqsym == rssym && pqsym == sym){
                        J->add(rsym, rrel, srel, D->get(psym, prel, qrel) * value);
                    }

                    if(qrsym == pssym && qrsym == sym){
                        K->add(qsym, qrel, rrel, D->get(psym, prel, srel) * value);
                    }

                    if(pabs!=qabs && rabs!=sabs && (pabs!=rabs || qabs!=sabs)){
                        /* (pq|sr) */
                        if(qssym == prsym && qssym == sym){
                            K->add(qsym, qrel, srel, D->get(psym, prel, rrel) * value);
                        }

                        /* (qp|rs) */
                        if(rssym == pqsym && rssym == sym){
                            J->add(rsym, rrel, srel, D->get(qsym, qrel, prel) * value);
                        }

                        if(prsym == qssym && prsym == sym){
                            K->add(psym, prel, rrel, D->get(qsym, qrel, srel) * value);
                        }

                        /* (qp|sr) */
                        if(pssym == qrsym && pssym == sym){
                            K->add(psym, prel, srel, D->get(qsym, qrel, rrel) * value);
                        }

                        /* (rs|pq) */
                        if(pqsym == rssym && pqsym == sym){
                            J->add(psym, prel, qrel, D->get(rsym, rrel, srel) * value);
                        }

                        if(pssym == qrsym && pssym == sym){
                            K->add(ssym, srel, prel, D->get(rsym, rrel, qrel) * value);
                        }

                        /* (sr|pq) */
                        if(pqsym == rssym && pqsym == sym){
                            J->add(psym, prel, qrel, D->get(ssym, srel, rrel) * value);
                        }

                        if(prsym == qssym && prsym == sym){
                            K->add(rsym, rrel, prel, D->get(ssym, srel, qrel) * value);
                        }

                        /* (rs|qp) */
                        if(qssym == prsym && qssym == sym){
                            K->add(ssym, srel, qrel, D->get(rsym, rrel, prel) * value);
                        }

                        /* (sr|qp) */
                        if(qrsym == pssym && qrsym == sym){
                            K->add(rsym, rrel, qrel, D->get(ssym, srel, prel) * value);
                        }
                    }else if(pabs!=qabs && rabs!=sabs && pabs==rabs && qabs==sabs){
                        /* (pq|sr) */
                        if(qssym == prsym && qssym == sym){
                            K->add(qsym, qrel, srel, D->get(psym, prel, rrel) * value);
                        }
                        /* (qp|rs) */
                        if(rssym == pqsym && rssym == sym){
                            J->add(rsym, rrel, srel, D->get(qsym, qrel, prel) * value);
                        }
                        if(prsym == qssym && prsym == sym){
                            K->add(psym, prel, rrel, D->get(qsym, qrel, srel) * value);
                        }

                        /* (qp|sr) */
                        if(pssym == qrsym && pssym == sym){
                            K->add(psym, prel, srel, D->get(qsym, qrel, rrel) * value);
                        }
                    }else if(pabs!=qabs && rabs==sabs){
                        /* (qp|rs) */
                        if(rssym == pqsym && rssym == sym){
                            J->add(rsym, rrel, srel, D->get(qsym, qrel, prel) * value);
                        }

                        if(prsym == qssym && prsym == sym){
                            K->add(psym, prel, rrel, D->get(qsym, qrel, srel) * value);
                        }

                        /* (rs|pq) */
                        if(pqsym == rssym && pqsym == sym){
                            J->add(psym, prel, qrel, D->get(rsym, rrel, srel) * value);
                        }

                        if(pssym == qrsym && pssym == sym){
                            K->add(ssym, srel, prel, D->get(rsym, rrel, qrel) * value);
                        }

                        /* (rs|qp) */
                        if(qssym == prsym && qssym == sym){
                            K->add(ssym, srel, qrel, D->get(rsym, rrel, prel) * value);
                        }
                    }else if(pabs==qabs && rabs!=sabs){
                        /* (pq|sr) */
                        if(qssym == prsym && qssym == sym){
                            K->add(qsym, qrel, srel, D->get(psym, prel, rrel) * value);
                        }

                        /* (rs|pq) */
                        if(pqsym == rssym && pqsym == sym){
                            J->add(psym, prel, qrel, D->get(rsym, rrel, srel) * value);
                        }

                        if(pssym == qrsym && pssym == sym){
                            K->add(ssym, srel, prel, D->get(rsym, rrel, qrel) * value);
                        }

                        /* (sr|pq) */
                        if(pqsym == rssym && pqsym == sym){
                            J->add(psym, prel, qrel, D->get(ssym, srel, rrel) * value);
                        }

                        if(prsym == qssym && prsym == sym){
                            K->add(rsym, rrel, prel, D->get(ssym, srel, qrel) * value);
                        }
                    }else if(pabs==qabs && rabs==sabs && (pabs!=rabs || qabs!=sabs)){
                        /* (rs|pq) */
                        if(pqsym == rssym && pqsym == sym){
                            J->add(psym, prel, qrel, D->get(rsym, rrel, srel) * value);
                        }

                        if(pssym == qrsym && pssym == sym){
                            K->add(ssym, srel, prel, D->get(rsym, rrel, qrel) * value);
                        }
                    }
                }
            } /* end loop through current buffer */
            if(!lastBuffer) iwl->fetch();
        }while(!lastBuffer);

        for (size_t N = 0; N < J_.size(); N++) {
            J_[N]->copy_lower_to_upper();
            if (K_[N]->symmetry()) K_[N]->transpose_this();
        }
    }else{
        // J and K to be handled separately
        do{
            lastBuffer = iwl->last_buffer();
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

                int pqsym = psym ^ qsym;
                int rssym = rsym ^ ssym;
                int qrsym = qsym ^ rsym;
                int pssym = psym ^ ssym;
                int prsym = psym ^ rsym;
                int qssym = qsym ^ ssym;
                // Coulomb terms
                for (size_t N = 0; N < J_.size(); N++) {

                    Matrix* J = J_[N].get();
                    Matrix* D = D_[N].get();

                    int sym = J->symmetry();

                    /* (pq|rs) */
                    if(pqsym == rssym && pqsym == sym){
                        J->add(rsym, rrel, srel, D->get(psym, prel, qrel) * value);
                    }

                    if(pabs!=qabs && rabs!=sabs && (pabs!=rabs || qabs!=sabs)){
                        /* (qp|rs) */
                        if(rssym == pqsym && rssym == sym){
                            J->add(rsym, rrel, srel, D->get(qsym, qrel, prel) * value);
                        }
                        /* (rs|pq) */
                        if(pqsym == rssym && pqsym == sym){
                            J->add(psym, prel, qrel, D->get(rsym, rrel, srel) * value);
                        }
                        /* (sr|pq) */
                        if(pqsym == rssym && pqsym == sym){
                            J->add(psym, prel, qrel, D->get(ssym, srel, rrel) * value);
                        }
                    }else if(pabs!=qabs && rabs!=sabs && pabs==rabs && qabs==sabs){
                        /* (qp|rs) */
                        if(rssym == pqsym && rssym == sym){
                            J->add(rsym, rrel, srel, D->get(qsym, qrel, prel) * value);
                        }
                    }else if(pabs!=qabs && rabs==sabs){
                        /* (qp|rs) */
                        if(rssym == pqsym && rssym == sym){
                            J->add(rsym, rrel, srel, D->get(qsym, qrel, prel) * value);
                        }
                        /* (rs|pq) */
                        if(pqsym == rssym && pqsym == sym){
                            J->add(psym, prel, qrel, D->get(rsym, rrel, srel) * value);
                        }
                    }else if(pabs==qabs && rabs!=sabs){
                        /* (rs|pq) */
                        if(pqsym == rssym && pqsym == sym){
                            J->add(psym, prel, qrel, D->get(rsym, rrel, srel) * value);
                        }
                        /* (sr|pq) */
                        if(pqsym == rssym && pqsym == sym){
                            J->add(psym, prel, qrel, D->get(ssym, srel, rrel) * value);
                        }
                    }else if(pabs==qabs && rabs==sabs && (pabs!=rabs || qabs!=sabs)){
                        /* (rs|pq) */
                        if(pqsym == rssym && pqsym == sym){
                            J->add(psym, prel, qrel, D->get(rsym, rrel, srel) * value);
                        }
                    }
                }

                // Exchange terms
                for (size_t N = 0; N < K_.size(); N++) {

                    Matrix* K = K_[N].get();
                    Matrix* D = D_[N].get();

                    int sym = K->symmetry();

                    /* (pq|rs) */
                    if(qrsym == pssym && qrsym == sym){
                        K->add(qsym, qrel, rrel, D->get(psym, prel, srel) * value);
                    }

                    if(pabs!=qabs && rabs!=sabs && (pabs!=rabs || qabs!=sabs)){
                        /* (pq|sr) */
                        if(qssym == prsym && qssym == sym){
                            K->add(qsym, qrel, srel, D->get(psym, prel, rrel) * value);
                        }
                        /* (qp|rs) */
                        if(prsym == qssym && prsym == sym){
                            K->add(psym, prel, rrel, D->get(qsym, qrel, srel) * value);
                        }
                        /* (qp|sr) */
                        if(pssym == qrsym && pssym == sym){
                            K->add(psym, prel, srel, D->get(qsym, qrel, rrel) * value);
                        }
                        /* (rs|pq) */
                        if(pssym == qrsym && pssym == sym){
                            K->add(ssym, srel, prel, D->get(rsym, rrel, qrel) * value);
                        }
                        /* (sr|pq) */
                        if(prsym == qssym && prsym == sym){
                            K->add(rsym, rrel, prel, D->get(ssym, srel, qrel) * value);
                        }
                        /* (rs|qp) */
                        if(qssym == prsym && qssym == sym){
                            K->add(ssym, srel, qrel, D->get(rsym, rrel, prel) * value);
                        }
                        /* (sr|qp) */
                        if(qrsym == pssym && qrsym == sym){
                            K->add(rsym, rrel, qrel, D->get(ssym, srel, prel) * value);
                        }
                    }else if(pabs!=qabs && rabs!=sabs && pabs==rabs && qabs==sabs){
                        /* (pq|sr) */
                        if(qssym == prsym && qssym == sym){
                            K->add(qsym, qrel, srel, D->get(psym, prel, rrel) * value);
                        }
                        /* (qp|rs) */
                        if(prsym == qssym && prsym == sym){
                            K->add(psym, prel, rrel, D->get(qsym, qrel, srel) * value);
                        }
                        /* (qp|sr) */
                        if(pssym == qrsym && pssym == sym){
                            K->add(psym, prel, srel, D->get(qsym, qrel, rrel) * value);
                        }
                    }else if(pabs!=qabs && rabs==sabs){
                        /* (qp|rs) */
                        if(prsym == qssym && prsym == sym){
                            K->add(psym, prel, rrel, D->get(qsym, qrel, srel) * value);
                        }
                        /* (rs|pq) */
                        if(pssym == qrsym && pssym == sym){
                            K->add(ssym, srel, prel, D->get(rsym, rrel, qrel) * value);
                        }
                        /* (rs|qp) */
                        if(qssym == prsym && qssym == sym){
                            K->add(ssym, srel, qrel, D->get(rsym, rrel, prel) * value);
                        }
                    }else if(pabs==qabs && rabs!=sabs){
                        /* (pq|sr) */
                        if(qssym == prsym && qssym == sym){
                            K->add(qsym, qrel, srel, D->get(psym, prel, rrel) * value);
                        }
                        /* (rs|pq) */
                        if(pssym == qrsym && pssym == sym){
                            K->add(ssym, srel, prel, D->get(rsym, rrel, qrel) * value);
                        }
                        /* (sr|pq) */
                        if(prsym == qssym && prsym == sym){
                            K->add(rsym, rrel, prel, D->get(ssym, srel, qrel) * value);
                        }
                    }else if(pabs==qabs && rabs==sabs && (pabs!=rabs || qabs!=sabs)){
                        /* (rs|pq) */
                        if(pssym == qrsym && pssym == sym){
                            K->add(ssym, srel, prel, D->get(rsym, rrel, qrel) * value);
                        }
                    }
                }
            } /* end loop through current buffer */
            if(!lastBuffer) iwl->fetch();
        }while(!lastBuffer);

        for (size_t N = 0; N < J_.size(); N++) {
            J_[N]->copy_lower_to_upper();
        }
        for (size_t N = 0; N < K_.size(); N++) {
            if (K_[N]->symmetry()) K_[N]->transpose_this();
        }
    }
    iwl->set_keep_flag(1);
    delete iwl;

    if(do_wK_){
        iwl = new IWL(psio.get(), PSIF_SO_ERF_TEI, cutoff_, 1, 1);
        lblptr = iwl->labels();
        valptr = iwl->values();

        do{
            lastBuffer = iwl->last_buffer();
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

                int qrsym = qsym ^ rsym;
                int pssym = psym ^ ssym;
                int prsym = psym ^ rsym;
                int qssym = qsym ^ ssym;

                // Exchange terms
                for (size_t N = 0; N < wK_.size(); N++) {

                    Matrix* wK = wK_[N].get();
                    Matrix* D = D_[N].get();

                    int sym = wK->symmetry();

                    /* (pq|rs) */
                    if(qrsym == pssym && qrsym == sym){
                        wK->add(qsym, qrel, rrel, D->get(psym, prel, srel) * value);
                    }

                    if(pabs!=qabs && rabs!=sabs && (pabs!=rabs || qabs!=sabs)){
                        /* (pq|sr) */
                        if(qssym == prsym && qssym == sym){
                            wK->add(qsym, qrel, srel, D->get(psym, prel, rrel) * value);
                        }
                        /* (qp|rs) */
                        if(prsym == qssym && prsym == sym){
                            wK->add(psym, prel, rrel, D->get(qsym, qrel, srel) * value);
                        }
                        /* (qp|sr) */
                        if(pssym == qrsym && pssym == sym){
                            wK->add(psym, prel, srel, D->get(qsym, qrel, rrel) * value);
                        }
                        /* (rs|pq) */
                        if(pssym == qrsym && pssym == sym){
                            wK->add(ssym, srel, prel, D->get(rsym, rrel, qrel) * value);
                        }
                        /* (sr|pq) */
                        if(prsym == qssym && prsym == sym){
                            wK->add(rsym, rrel, prel, D->get(ssym, srel, qrel) * value);
                        }
                        /* (rs|qp) */
                        if(qssym == prsym && qssym == sym){
                            wK->add(ssym, srel, qrel, D->get(rsym, rrel, prel) * value);
                        }
                        /* (sr|qp) */
                        if(qrsym == pssym && qrsym == sym){
                            wK->add(rsym, rrel, qrel, D->get(ssym, srel, prel) * value);
                        }
                    }else if(pabs!=qabs && rabs!=sabs && pabs==rabs && qabs==sabs){
                        /* (pq|sr) */
                        if(qssym == prsym && qssym == sym){
                            wK->add(qsym, qrel, srel, D->get(psym, prel, rrel) * value);
                        }
                        /* (qp|rs) */
                        if(prsym == qssym && prsym == sym){
                            wK->add(psym, prel, rrel, D->get(qsym, qrel, srel) * value);
                        }
                        /* (qp|sr) */
                        if(pssym == qrsym && pssym == sym){
                            wK->add(psym, prel, srel, D->get(qsym, qrel, rrel) * value);
                        }
                    }else if(pabs!=qabs && rabs==sabs){
                        /* (qp|rs) */
                        if(prsym == qssym && prsym == sym){
                            wK->add(psym, prel, rrel, D->get(qsym, qrel, srel) * value);
                        }
                        /* (rs|pq) */
                        if(pssym == qrsym && pssym == sym){
                            wK->add(ssym, srel, prel, D->get(rsym, rrel, qrel) * value);
                        }
                        /* (rs|qp) */
                        if(qssym == prsym && qssym == sym){
                            wK->add(ssym, srel, qrel, D->get(rsym, rrel, prel) * value);
                        }
                    }else if(pabs==qabs && rabs!=sabs){
                        /* (pq|sr) */
                        if(qssym == prsym && qssym == sym){
                            wK->add(qsym, qrel, srel, D->get(psym, prel, rrel) * value);
                        }
                        /* (rs|pq) */
                        if(pssym == qrsym && pssym == sym){
                            wK->add(ssym, srel, prel, D->get(rsym, rrel, qrel) * value);
                        }
                        /* (sr|pq) */
                        if(prsym == qssym && prsym == sym){
                            wK->add(rsym, rrel, prel, D->get(ssym, srel, qrel) * value);
                        }
                    }else if(pabs==qabs && rabs==sabs && (pabs!=rabs || qabs!=sabs)){
                        /* (rs|pq) */
                        if(pssym == qrsym && pssym == sym){
                            wK->add(ssym, srel, prel, D->get(rsym, rrel, qrel) * value);
                        }
                    }
                }
            } /* end loop through current buffer */
            if(!lastBuffer) iwl->fetch();
        }while(!lastBuffer);

        for (size_t N = 0; N < wK_.size(); N++) {
            if (wK_[N]->symmetry()) wK_[N]->transpose_this();
        }
        iwl->set_keep_flag(1);
        delete iwl;
    }
}
void DiskJK::postiterations()
{
    delete[] so2symblk_;
    delete[] so2index_;
}

}
