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

    if(_transformationType == Restricted){
        /*
         * Start by transforming the OPDM to the SO basis
         */
        double *tempSo = new double[_nTriSo];
        double *tempMo = new double[_nTriMo];
        double **tempOPDM = block_matrix(_nmo, _nmo);

        // Before we start, we need to build a mapping array from correlated (no frozen virtuals)
        // to Pitzer, so that the density is ordered consistently with the MO coefficients
        int nActive = _nmo - _nfzv;
        int *toPitzer = new int[nActive];
        size_t corrCount = 0;
        size_t pitzerOffset = 0;
        // Frozen DOCC
        for(int h = 0; h < _nirreps; ++h){
            for(int n = 0; n < _frzcpi[h]; ++n){
                toPitzer[corrCount++] = pitzerOffset + n;
            }
            pitzerOffset += _mopi[h];
        }
        // Active DOCC
        pitzerOffset = 0;
        for(int h = 0; h < _nirreps; ++h){
            for(int n = _frzcpi[h]; n < _clsdpi[h]; ++n){
                toPitzer[corrCount++] = pitzerOffset + n;
            }
            pitzerOffset += _mopi[h];
        }
        // Active VIRT
        pitzerOffset = 0;
        for(int h = 0; h < _nirreps; ++h){
            for(int n = _clsdpi[h]; n < _mopi[h] - _frzvpi[h]; ++n){
                toPitzer[corrCount++] = pitzerOffset + n;
            }
            pitzerOffset += _mopi[h];
        }

        int *order = init_int_array(_nmo);
        // We want to keep Pitzer ordering, so this is just an identity mapping
        for(int n = 0; n < _nmo; ++n) order[n] = n;

        /*
         * The OPDM
         */
        _psio->open(PSIF_MO_OPDM, PSIO_OPEN_OLD);
        _psio->read_entry(PSIF_MO_OPDM, "MO-basis OPDM", (char *) tempOPDM[0], sizeof(double)*nActive*nActive);
        for(int p = 0; p < nActive; ++p){
          for(int q = 0; q <= p; ++q){
              int P = toPitzer[p];
              int Q = toPitzer[q];
              size_t PQ = INDEX(P,Q);
              tempMo[PQ] = 0.5 * (tempOPDM[p][q] + tempOPDM[q][p]);
          }
        }
        _psio->close(PSIF_MO_OPDM, 1);
        if(_print>4){
            fprintf(outfile, "The MO basis OPDM");
            print_array(tempMo, _nmo, outfile);
        }
        for(int n = 0; n < _nTriSo; ++n) tempSo[n] = 0.0;
        for(int h = 0, moOffset = 0, soOffset = 0; h < _nirreps; ++h){
            trans_one(_mopi[h], _sopi[h], tempMo, tempSo, _Ca[h], moOffset, &(order[soOffset]), true);
            soOffset += _sopi[h];
            moOffset += _mopi[h];
        }
        if(_print>4){
            fprintf(outfile, "The SO basis OPDM");
            print_array(tempSo, _nso, outfile);
        }
        _psio->open(PSIF_AO_OPDM, PSIO_OPEN_OLD);
        _psio->write_entry(PSIF_AO_OPDM, "SO-basis OPDM", (char *) tempSo, sizeof(double)*_nTriSo);
        _psio->close(PSIF_AO_OPDM, 1);

        /*
         * The Lagrangian
         */
        _psio->open(PSIF_MO_LAG, PSIO_OPEN_OLD);
        _psio->read_entry(PSIF_MO_LAG, "MO-basis Lagrangian", (char *) tempOPDM[0], sizeof(double)*nActive*nActive);
        for(int p = 0; p < nActive; ++p){
          for(int q = 0; q <= p; ++q){
              int P = toPitzer[p];
              int Q = toPitzer[q];
              size_t PQ = INDEX(P,Q);
              tempMo[PQ] = 0.5 * (tempOPDM[p][q] + tempOPDM[q][p]);
          }
        }
        _psio->close(PSIF_MO_LAG, 1);
        if(_print>4){
            fprintf(outfile, "The MO basis Lagrangian\n");
            print_array(tempMo, _nmo, outfile);
        }
        for(int n = 0; n < _nTriSo; ++n) tempSo[n] = 0.0;
        for(int h = 0, moOffset = 0, soOffset = 0; h < _nirreps; ++h){
            trans_one(_mopi[h], _sopi[h], tempMo, tempSo, _Ca[h], moOffset, &(order[soOffset]), true);
            soOffset += _sopi[h];
            moOffset += _mopi[h];
        }
        if(_print>4){
            fprintf(outfile, "The SO basis Lagrangian\n");
            print_array(tempSo, _nso, outfile);
        }
        _psio->open(PSIF_AO_OPDM, PSIO_OPEN_OLD);
        _psio->write_entry(PSIF_AO_OPDM, "SO-basis Lagrangian", (char *) tempSo, sizeof(double)*_nTriSo);
        _psio->close(PSIF_AO_OPDM, 1);

        free_block(tempOPDM);
        delete [] tempMo;
        delete [] tempSo;

        /*
         * Now, work on the TPDM
         */
        backtransform_tpdm_restricted();
    }else{
        throw FeatureNotImplemented("libtrans", "Unrestricted density transformations", __FILE__, __LINE__);
    }
}
