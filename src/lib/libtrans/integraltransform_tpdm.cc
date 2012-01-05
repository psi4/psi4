#include "integraltransform.h"
#include <libpsio/psio.hpp>
#include <libciomr/libciomr.h>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <libmints/matrix.h>
#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include "psifiles.h"
#include "ccfiles.h"
#include "mospace.h"
#define EXTERN
#include <libdpd/dpd.gbl>

using namespace psi;

/**
 * Transform the one- and two-particle density matrices from the MO to the SO basis
 */
void
IntegralTransform::backtransform_density()
{
    check_initialized();
    // This limitation can be remedied by accounting for the fact that Pitzer orbital numbering is not
    // dense, so certain quantities must be alloc'd for the full MO space.  It's no limitation, though
    if(frozenOrbitals_ != None)
        throw SanityCheckError("No orbitals can be frozen in density matrix transformations\n",
                               __FILE__, __LINE__);
    // The full MO space must be in the list of spaces used, let's check
    bool allFound = false;
    for(int i = 0; i < spacesUsed_.size(); ++i)
        if(spacesUsed_[i] == MOSPACE_ALL) allFound = true;
    if(!allFound)
        throw PSIEXCEPTION("MOSpace::all must be amongst the spaces passed "
                           "to the integral object's constructor");

    int nActive = nmo_ - nfzv_;
    double *tempSo = new double[nTriSo_];
    double *tempMo = new double[nTriMo_];
    ::memset((void*)tempSo, '\0', nTriSo_*sizeof(double));
    double **tempOPDM = block_matrix(nmo_, nmo_);
    int *order = init_int_array(nmo_);
    // We want to keep Pitzer ordering, so this is just an identity mapping
    for(int n = 0; n < nmo_; ++n) order[n] = n;

    psio_->open(PSIF_AO_OPDM, PSIO_OPEN_OLD);
    psio_->open(PSIF_MO_OPDM, PSIO_OPEN_OLD);
    psio_->open(PSIF_MO_LAG, PSIO_OPEN_OLD);


    if(transformationType_ == Restricted){
        /*
         * Start by transforming the OPDM to the SO basis
         */
        psio_->read_entry(PSIF_MO_OPDM, "MO-basis OPDM", (char *) tempOPDM[0], sizeof(double)*nActive*nActive);
        for(int p = 0; p < nActive; ++p){
          for(int q = 0; q <= p; ++q){
              int P = aCorrToPitzer_[p];
              int Q = aCorrToPitzer_[q];
              size_t PQ = INDEX(P,Q);
              tempMo[PQ] = 0.5 * (tempOPDM[p][q] + tempOPDM[q][p]);
          }
        }
        if(print_>4){
            fprintf(outfile, "The MO basis OPDM");
            print_array(tempMo, nmo_, outfile);
        }

        for(int h = 0, moOffset = 0, soOffset = 0; h < nirreps_; ++h){
            double **pC = Ca_->pointer(h);
            trans_one(mopi_[h], sopi_[h], tempMo, tempSo, pC, moOffset, &(order[soOffset]),true);
            soOffset += sopi_[h];
            moOffset += mopi_[h];
        }
        if(print_>4){
            fprintf(outfile, "The SO basis OPDM");
            print_array(tempSo, nso_, outfile);
        }
        psio_->write_entry(PSIF_AO_OPDM, "SO-basis OPDM", (char *) tempSo, sizeof(double)*nTriSo_);

        /*
         * The Lagrangian
         */
        psio_->read_entry(PSIF_MO_LAG, "MO-basis Lagrangian", (char *) tempOPDM[0], sizeof(double)*nActive*nActive);
        for(int p = 0; p < nActive; ++p){
          for(int q = 0; q <= p; ++q){
              int P = aCorrToPitzer_[p];
              int Q = aCorrToPitzer_[q];
              size_t PQ = INDEX(P,Q);
              tempMo[PQ] = 0.5 * (tempOPDM[p][q] + tempOPDM[q][p]);
          }
        }
        if(print_>4){
            fprintf(outfile, "The MO basis Lagrangian\n");
            print_array(tempMo, nmo_, outfile);
        }

        for(int h = 0, moOffset = 0, soOffset = 0; h < nirreps_; ++h){
            double **pC = Ca_->pointer(h);
            trans_one(mopi_[h], sopi_[h], tempMo, tempSo, pC, moOffset, &(order[soOffset]), true);
            soOffset += sopi_[h];
            moOffset += mopi_[h];
        }
        if(print_>4){
            fprintf(outfile, "The SO basis Lagrangian\n");
            print_array(tempSo, nso_, outfile);
        }
        psio_->write_entry(PSIF_AO_OPDM, "SO-basis Lagrangian", (char *) tempSo, sizeof(double)*nTriSo_);

        /*
         * Now, work on the TPDM
         */
        backtransform_tpdm_restricted();
    }else{
        /*
         * Start by transforming the OPDM to the SO basis
         */
        psio_->read_entry(PSIF_MO_OPDM, "MO-basis Alpha OPDM", (char *) tempOPDM[0], sizeof(double)*nActive*nActive);
        for(int p = 0; p < nActive; ++p){
          for(int q = 0; q <= p; ++q){
              int P = aCorrToPitzer_[p];
              int Q = aCorrToPitzer_[q];
              size_t PQ = INDEX(P,Q);
              tempMo[PQ] = 0.5 * (tempOPDM[p][q] + tempOPDM[q][p]);
          }
        }
        if(print_>4){
            fprintf(outfile, "The MO basis Alpha OPDM");
            print_array(tempMo, nmo_, outfile);
        }

        for(int h = 0, moOffset = 0, soOffset = 0; h < nirreps_; ++h){
            double **pC = Ca_->pointer(h);
            trans_one(mopi_[h], sopi_[h], tempMo, tempSo, pC, moOffset, &(order[soOffset]), true);
            soOffset += sopi_[h];
            moOffset += mopi_[h];
        }
        psio_->read_entry(PSIF_MO_OPDM, "MO-basis Beta OPDM", (char *) tempOPDM[0], sizeof(double)*nActive*nActive);
        for(int p = 0; p < nActive; ++p){
          for(int q = 0; q <= p; ++q){
              int P = bCorrToPitzer_[p];
              int Q = bCorrToPitzer_[q];
              size_t PQ = INDEX(P,Q);
              tempMo[PQ] = 0.5 * (tempOPDM[p][q] + tempOPDM[q][p]);
          }
        }
        if(print_>4){
            fprintf(outfile, "The MO basis Beta OPDM");
            print_array(tempMo, nmo_, outfile);
        }
        for(int h = 0, moOffset = 0, soOffset = 0; h < nirreps_; ++h){
            // Note the final argument here, which tells the code to accumulate the beta contribution into the alpha
            double **pC = Cb_->pointer(h);
            trans_one(mopi_[h], sopi_[h], tempMo, tempSo, pC, moOffset, &(order[soOffset]), true, 1.0);
            soOffset += sopi_[h];
            moOffset += mopi_[h];
        }
        if(print_>4){
            fprintf(outfile, "The SO basis OPDM");
            print_array(tempSo, nso_, outfile);
        }
        psio_->write_entry(PSIF_AO_OPDM, "SO-basis OPDM", (char *) tempSo, sizeof(double)*nTriSo_);

        /*
         * The Lagrangian
         */
        psio_->read_entry(PSIF_MO_LAG, "MO-basis Alpha Lagrangian", (char *) tempOPDM[0], sizeof(double)*nActive*nActive);
        for(int p = 0; p < nActive; ++p){
          for(int q = 0; q <= p; ++q){
              int P = aCorrToPitzer_[p];
              int Q = aCorrToPitzer_[q];
              size_t PQ = INDEX(P,Q);
              tempMo[PQ] = 0.5 * (tempOPDM[p][q] + tempOPDM[q][p]);
          }
        }
        if(print_>4){
            fprintf(outfile, "The MO basis Alpha Lagrangian\n");
            print_array(tempMo, nmo_, outfile);
        }
        for(int h = 0, moOffset = 0, soOffset = 0; h < nirreps_; ++h){
            double **pC = Ca_->pointer(h);
            trans_one(mopi_[h], sopi_[h], tempMo, tempSo, pC, moOffset, &(order[soOffset]), true);
            soOffset += sopi_[h];
            moOffset += mopi_[h];
        }
        psio_->read_entry(PSIF_MO_LAG, "MO-basis Beta Lagrangian", (char *) tempOPDM[0], sizeof(double)*nActive*nActive);
        for(int p = 0; p < nActive; ++p){
          for(int q = 0; q <= p; ++q){
              int P = bCorrToPitzer_[p];
              int Q = bCorrToPitzer_[q];
              size_t PQ = INDEX(P,Q);
              tempMo[PQ] = 0.5 * (tempOPDM[p][q] + tempOPDM[q][p]);
          }
        }
        if(print_>4){
            fprintf(outfile, "The MO basis Beta Lagrangian\n");
            print_array(tempMo, nmo_, outfile);
        }
        for(int h = 0, moOffset = 0, soOffset = 0; h < nirreps_; ++h){
            // Note the final argument here, which tells the code to accumulate the beta contribution into the alpha
            double **pC = Cb_->pointer(h);
            trans_one(mopi_[h], sopi_[h], tempMo, tempSo, pC, moOffset, &(order[soOffset]), true, 1.0);
            soOffset += sopi_[h];
            moOffset += mopi_[h];
        }
        if(print_>4){
            fprintf(outfile, "The SO basis Lagrangian\n");
            print_array(tempSo, nso_, outfile);
        }
        psio_->write_entry(PSIF_AO_OPDM, "SO-basis Lagrangian", (char *) tempSo, sizeof(double)*nTriSo_);

        /*
         * Now, work on the TPDM
         */
        backtransform_tpdm_unrestricted();
    }

    free(order);
    free_block(tempOPDM);
    delete [] tempMo;
    delete [] tempSo;

    psio_->close(PSIF_AO_OPDM, 1);
    psio_->close(PSIF_MO_OPDM, 1);
    psio_->close(PSIF_MO_LAG, 1);
}
