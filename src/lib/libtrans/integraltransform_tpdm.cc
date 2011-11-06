#include "integraltransform.h"
#include <libchkpt/chkpt.hpp>
#include <libpsio/psio.hpp>
#include <libciomr/libciomr.h>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
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
    if(_frozenOrbitals != None)
        throw SanityCheckError("No orbitals can be frozen in density matrix transformations\n",
                               __FILE__, __LINE__);
    // The full MO space must be in the list of spaces used, let's check
    bool allFound = false;
    for(int i = 0; i < _spacesUsed.size(); ++i)
        if(_spacesUsed[i] == MOSPACE_ALL) allFound = true;
    if(!allFound)
        throw PSIEXCEPTION("MOSpace::all must be amongst the spaces passed "
                           "to the integral object's constructor");

    int nActive = _nmo - _nfzv;
    double *tempSo = new double[_nTriSo];
    double *tempMo = new double[_nTriMo];
    double **tempOPDM = block_matrix(_nmo, _nmo);
    int *order = init_int_array(_nmo);
    // We want to keep Pitzer ordering, so this is just an identity mapping
    for(int n = 0; n < _nmo; ++n) order[n] = n;

    _psio->open(PSIF_AO_OPDM, PSIO_OPEN_OLD);
    _psio->open(PSIF_MO_OPDM, PSIO_OPEN_OLD);
    _psio->open(PSIF_MO_LAG, PSIO_OPEN_OLD);


    if(_transformationType == Restricted){
        /*
         * Start by transforming the OPDM to the SO basis
         */
        _psio->read_entry(PSIF_MO_OPDM, "MO-basis OPDM", (char *) tempOPDM[0], sizeof(double)*nActive*nActive);
        for(int p = 0; p < nActive; ++p){
          for(int q = 0; q <= p; ++q){
              int P = _aCorrToPitzer[p];
              int Q = _aCorrToPitzer[q];
              size_t PQ = INDEX(P,Q);
              tempMo[PQ] = 0.5 * (tempOPDM[p][q] + tempOPDM[q][p]);
          }
        }
        if(_print>4){
            fprintf(outfile, "The MO basis OPDM");
            print_array(tempMo, _nmo, outfile);
        }

        for(int h = 0, moOffset = 0, soOffset = 0; h < _nirreps; ++h){
            trans_one(_mopi[h], _sopi[h], tempMo, tempSo, _Ca[h], moOffset, &(order[soOffset]),true);
            soOffset += _sopi[h];
            moOffset += _mopi[h];
        }
        if(_print>4){
            fprintf(outfile, "The SO basis OPDM");
            print_array(tempSo, _nso, outfile);
        }
        _psio->write_entry(PSIF_AO_OPDM, "SO-basis OPDM", (char *) tempSo, sizeof(double)*_nTriSo);

        /*
         * The Lagrangian
         */
        _psio->read_entry(PSIF_MO_LAG, "MO-basis Lagrangian", (char *) tempOPDM[0], sizeof(double)*nActive*nActive);
        for(int p = 0; p < nActive; ++p){
          for(int q = 0; q <= p; ++q){
              int P = _aCorrToPitzer[p];
              int Q = _aCorrToPitzer[q];
              size_t PQ = INDEX(P,Q);
              tempMo[PQ] = 0.5 * (tempOPDM[p][q] + tempOPDM[q][p]);
          }
        }
        if(_print>4){
            fprintf(outfile, "The MO basis Lagrangian\n");
            print_array(tempMo, _nmo, outfile);
        }

        for(int h = 0, moOffset = 0, soOffset = 0; h < _nirreps; ++h){
            trans_one(_mopi[h], _sopi[h], tempMo, tempSo, _Ca[h], moOffset, &(order[soOffset]));
            soOffset += _sopi[h];
            moOffset += _mopi[h];
        }
        if(_print>4){
            fprintf(outfile, "The SO basis Lagrangian\n");
            print_array(tempSo, _nso, outfile);
        }
        _psio->write_entry(PSIF_AO_OPDM, "SO-basis Lagrangian", (char *) tempSo, sizeof(double)*_nTriSo);

        /*
         * Now, work on the TPDM
         */
        backtransform_tpdm_restricted();
    }else{
        /*
         * Start by transforming the OPDM to the SO basis
         */
        _psio->read_entry(PSIF_MO_OPDM, "MO-basis Alpha OPDM", (char *) tempOPDM[0], sizeof(double)*nActive*nActive);
        for(int p = 0; p < nActive; ++p){
          for(int q = 0; q <= p; ++q){
              int P = _aCorrToPitzer[p];
              int Q = _aCorrToPitzer[q];
              size_t PQ = INDEX(P,Q);
              tempMo[PQ] = 0.5 * (tempOPDM[p][q] + tempOPDM[q][p]);
          }
        }
        if(_print>4){
            fprintf(outfile, "The MO basis Alpha OPDM");
            print_array(tempMo, _nmo, outfile);
        }

        for(int h = 0, moOffset = 0, soOffset = 0; h < _nirreps; ++h){
            trans_one(_mopi[h], _sopi[h], tempMo, tempSo, _Ca[h], moOffset, &(order[soOffset]), true);
            soOffset += _sopi[h];
            moOffset += _mopi[h];
        }
        _psio->read_entry(PSIF_MO_OPDM, "MO-basis Beta OPDM", (char *) tempOPDM[0], sizeof(double)*nActive*nActive);
        for(int p = 0; p < nActive; ++p){
          for(int q = 0; q <= p; ++q){
              int P = _bCorrToPitzer[p];
              int Q = _bCorrToPitzer[q];
              size_t PQ = INDEX(P,Q);
              tempMo[PQ] = 0.5 * (tempOPDM[p][q] + tempOPDM[q][p]);
          }
        }
        if(_print>4){
            fprintf(outfile, "The MO basis Beta OPDM");
            print_array(tempMo, _nmo, outfile);
        }
        for(int h = 0, moOffset = 0, soOffset = 0; h < _nirreps; ++h){
            // Note the final argument here, which tells the code to accumulate the beta contribution into the alpha
            trans_one(_mopi[h], _sopi[h], tempMo, tempSo, _Cb[h], moOffset, &(order[soOffset]), true, 1.0);
            soOffset += _sopi[h];
            moOffset += _mopi[h];
        }
        if(_print>4){
            fprintf(outfile, "The SO basis OPDM");
            print_array(tempSo, _nso, outfile);
        }
        _psio->write_entry(PSIF_AO_OPDM, "SO-basis OPDM", (char *) tempSo, sizeof(double)*_nTriSo);

        /*
         * The Lagrangian
         */
        _psio->read_entry(PSIF_MO_LAG, "MO-basis Alpha Lagrangian", (char *) tempOPDM[0], sizeof(double)*nActive*nActive);
        for(int p = 0; p < nActive; ++p){
          for(int q = 0; q <= p; ++q){
              int P = _aCorrToPitzer[p];
              int Q = _aCorrToPitzer[q];
              size_t PQ = INDEX(P,Q);
              tempMo[PQ] = 0.5 * (tempOPDM[p][q] + tempOPDM[q][p]);
          }
        }
        if(_print>4){
            fprintf(outfile, "The MO basis Alpha Lagrangian\n");
            print_array(tempMo, _nmo, outfile);
        }
        for(int h = 0, moOffset = 0, soOffset = 0; h < _nirreps; ++h){
            trans_one(_mopi[h], _sopi[h], tempMo, tempSo, _Ca[h], moOffset, &(order[soOffset]), true);
            soOffset += _sopi[h];
            moOffset += _mopi[h];
        }
        _psio->read_entry(PSIF_MO_LAG, "MO-basis Beta Lagrangian", (char *) tempOPDM[0], sizeof(double)*nActive*nActive);
        for(int p = 0; p < nActive; ++p){
          for(int q = 0; q <= p; ++q){
              int P = _bCorrToPitzer[p];
              int Q = _bCorrToPitzer[q];
              size_t PQ = INDEX(P,Q);
              tempMo[PQ] = 0.5 * (tempOPDM[p][q] + tempOPDM[q][p]);
          }
        }
        if(_print>4){
            fprintf(outfile, "The MO basis Beta Lagrangian\n");
            print_array(tempMo, _nmo, outfile);
        }
        for(int h = 0, moOffset = 0, soOffset = 0; h < _nirreps; ++h){
            // Note the final argument here, which tells the code to accumulate the beta contribution into the alpha
            trans_one(_mopi[h], _sopi[h], tempMo, tempSo, _Cb[h], moOffset, &(order[soOffset]), true, 1.0);
            soOffset += _sopi[h];
            moOffset += _mopi[h];
        }
        if(_print>4){
            fprintf(outfile, "The SO basis Lagrangian\n");
            print_array(tempSo, _nso, outfile);
        }
        _psio->write_entry(PSIF_AO_OPDM, "SO-basis Lagrangian", (char *) tempSo, sizeof(double)*_nTriSo);

        /*
         * Now, work on the TPDM
         */
        backtransform_tpdm_unrestricted();
    }

    free(order);
    free_block(tempOPDM);
    delete [] tempMo;
    delete [] tempSo;

    _psio->close(PSIF_AO_OPDM, 1);
    _psio->close(PSIF_MO_OPDM, 1);
    _psio->close(PSIF_MO_LAG, 1);
}
