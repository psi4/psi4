#include "dcft.h"
#include <libdpd/dpd.h>
#include <libiwl/iwl.hpp>
#include <libtrans/integraltransform.h>
#include "defines.h"

using namespace boost;

namespace psi{ namespace dcft{

/**
 * Updates the MO coefficients, transforms the integrals into both chemists'
 * and physcists' notation.
 */
void
DCFTSolver::transform_integrals()
{
    dpdbuf4 I, Irs, Isr;

    _ints->update_orbitals();
    if(_print > 1){
        fprintf(outfile, "\tTransforming integrals...\n");
        fflush(outfile);
    }
    _ints->set_print(_print - 2 >= 0 ? _print - 2 : 0);

    // Generate the integrals in various spaces in chemists' notation
    _ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::vir);
    _ints->transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::occ, MOSpace::occ);
    _ints->transform_tei(MOSpace::vir, MOSpace::vir, MOSpace::occ, MOSpace::occ);
    _ints->transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::vir, MOSpace::vir);
    if((_options.get_str("AO_BASIS") == "NONE")){
        _ints->transform_tei(MOSpace::vir, MOSpace::vir, MOSpace::vir, MOSpace::vir);
    }

    // The integral object closes this file - we need to re-open it.
    _psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    /*
     * Re-sort the chemists' notation integrals to physisicts' notation
     * (pq|rs) = <pr|qs>
     */
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"),0, "MO Ints (OV|OV)");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[O,O]"), ID("[V,V]"), "MO Ints <OO|VV>");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"),0, "MO Ints <OO|VV>");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, rspq, ID("[V,V]"), ID("[O,O]"), "MO Ints <VV|OO>");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "MO Ints (OV|ov)");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[O,o]"), ID("[V,v]"), "MO Ints <Oo|Vv>");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints (ov|ov)");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[o,o]"), ID("[v,v]"), "MO Ints <oo|vv>");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo|vv>");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, rspq, ID("[v,v]"), ID("[o,o]"), "MO Ints <vv|oo>");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O>=O]+"), ID("[O>=O]+"), 0, "MO Ints (OO|OO)");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[O,O]"), ID("[O,O]"), "MO Ints <OO|OO>");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[o,o]"),
                  ID("[O>=O]+"), ID("[o>=o]+"), 0, "MO Ints (OO|oo)");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[O,o]"), ID("[O,o]"), "MO Ints <Oo|Oo>");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[o,o]"),
                  ID("[O>=O]+"), ID("[o>=o]+"), 0, "MO Ints (OO|oo)");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, rspq, ID("[o,o]"), ID("[O,O]"), "MO Ints (oo|OO)");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o>=o]+"), ID("[o>=o]+"), 0, "MO Ints (oo|oo)");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[o,o]"), ID("[o,o]"), "MO Ints <oo|oo>");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[o,o]"),
                  ID("[V>=V]+"), ID("[o>=o]+"), 0, "MO Ints (VV|oo)");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, sqrp, ID("[o,V]"), ID("[o,V]"), "MO Ints <oV|oV>");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[o,o]"),
                  ID("[V>=V]+"), ID("[o>=o]+"), 0, "MO Ints (VV|oo)");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, rspq, ID("[o,o]"), ID("[V,V]"), "MO Ints (oo|VV)");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[v,v]"),
                  ID("[O>=O]+"), ID("[v>=v]+"), 0, "MO Ints (OO|vv)");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, rspq, ID("[v,v]"), ID("[O,O]"), "MO Ints (vv|OO)");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>=O]+"), ID("[V>=V]+"), 0, "MO Ints (OO|VV)");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[O,V]"), ID("[O,V]"), "MO Ints <OV|OV>");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[v,v]"),
                  ID("[O>=O]+"), ID("[v>=v]+"), 0, "MO Ints (OO|vv)");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[O,v]"), ID("[O,v]"), "MO Ints <Ov|Ov>");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>=o]+"), ID("[v>=v]+"), 0, "MO Ints (oo|vv)");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[o,v]"), ID("[o,v]"), "MO Ints <ov|ov>");
    dpd_buf4_close(&I);

    if(_options.get_str("AO_BASIS") == "NONE"){
        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                      ID("[V>=V]+"), ID("[V>=V]+"), 0, "MO Ints (VV|VV)");
        dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[V,V]"), ID("[V,V]"), "MO Ints <VV|VV>");
        dpd_buf4_close(&I);

        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[v,v]"),
                      ID("[V>=V]+"), ID("[v>=v]+"), 0, "MO Ints (VV|vv)");
        dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[V,v]"), ID("[V,v]"), "MO Ints <Vv|Vv>");
        dpd_buf4_close(&I);

        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[v,v]"),
                      ID("[V>=V]+"), ID("[v>=v]+"), 0, "MO Ints (VV|vv)");
        dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, rspq, ID("[v,v]"), ID("[V,V]"), "MO Ints (vv|VV)");
        dpd_buf4_close(&I);

        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[v,v]"),
                      ID("[v>=v]+"), ID("[v>=v]+"), 0, "MO Ints (vv|vv)");
        dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[v,v]"), ID("[v,v]"), "MO Ints <vv|vv>");
        dpd_buf4_close(&I);
    }

    /*
     * Antisymmetrize the <OV|OV> and <ov|ov> integrals
     */
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV|OV>");
    dpd_buf4_copy(&I, PSIF_LIBTRANS_DPD, "MO Ints <OV|OV> - (OV|OV)");
    dpd_buf4_close(&I);
    dpd_buf4_init(&Irs, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV|OV> - (OV|OV)");
    dpd_buf4_init(&Isr, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
    for(int h = 0; h < _nIrreps; ++h){
        dpd_buf4_mat_irrep_init(&Irs, h);
        dpd_buf4_mat_irrep_init(&Isr, h);
        dpd_buf4_mat_irrep_rd(&Irs, h);
        dpd_buf4_mat_irrep_rd(&Isr, h);
        for(int row = 0; row < Irs.params->rowtot[h]; ++row){
            for(int col = 0; col < Irs.params->coltot[h]; ++col){
                Irs.matrix[h][row][col] -= Isr.matrix[h][row][col];
            }
        }
        dpd_buf4_mat_irrep_wrt(&Irs, h);
        dpd_buf4_mat_irrep_close(&Irs, h);
        dpd_buf4_mat_irrep_close(&Isr, h);
    }

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints <ov|ov>");
    dpd_buf4_copy(&I, PSIF_LIBTRANS_DPD, "MO Ints <ov|ov> - (ov|ov)");
    dpd_buf4_close(&I);
    dpd_buf4_init(&Irs, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints <ov|ov> - (ov|ov)");
    dpd_buf4_init(&Isr, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints (ov|ov)");
    for(int h = 0; h < _nIrreps; ++h){
        dpd_buf4_mat_irrep_init(&Irs, h);
        dpd_buf4_mat_irrep_init(&Isr, h);
        dpd_buf4_mat_irrep_rd(&Irs, h);
        dpd_buf4_mat_irrep_rd(&Isr, h);
        for(int row = 0; row < Irs.params->rowtot[h]; ++row){
            for(int col = 0; col < Irs.params->coltot[h]; ++col){
                Irs.matrix[h][row][col] -= Isr.matrix[h][row][col];
            }
        }
        dpd_buf4_mat_irrep_wrt(&Irs, h);
        dpd_buf4_mat_irrep_close(&Irs, h);
        dpd_buf4_mat_irrep_close(&Isr, h);
    }


    // The integral transformation object also provided us with the Fock matrix
    // in the current basis, from which we can get the new denominator matrices.
    // N.B. These are not neccesarily the eigenvalues, rather they are the diagonal
    // elements of F0 in the current basis; F is diagonal, not F0.
    build_denominators();

    _psio->close(PSIF_LIBTRANS_DPD, 1);
}

/**
 * Consructs the denominators using the diagonal elements of the unperturbed Fock matrix.
 * Also builds the MO coefficient tensors in the alpha/beta, occ/vir spaces.
 */
void
DCFTSolver::build_denominators()
{
    dpdbuf4 D;
    dpdfile2 F;

    double *aF0 = new double[_nTriSo];
    double *bF0 = new double[_nTriSo];

    IWL::read_one(_psio.get(), PSIF_OEI, PSIF_MO_A_FOCK, aF0, _nTriSo, 0, 0, outfile);
    IWL::read_one(_psio.get(), PSIF_OEI, PSIF_MO_B_FOCK, bF0, _nTriSo, 0, 0, outfile);

    double *aOccEvals = new double [_nAOcc];
    double *bOccEvals = new double [_nBOcc];
    double *aVirEvals = new double [_nAVir];
    double *bVirEvals = new double [_nBVir];
    // Pick out the diagonal elements of the Fock matrix, making sure that they are in the order
    // used by the DPD library, i.e. starting from zero for each space and ordering by irrep
    int pitzerOffset = 0;
    int aOccCount = 0, bOccCount = 0, aVirCount = 0, bVirCount = 0, aCount = 0, bCount = 0;
    double **Ca = _chkpt->rd_alpha_scf();
    double **Cb = _chkpt->rd_beta_scf();
    for(int h = 0, soOffset = 0; h < _nIrreps; ++h){
        bCount = aCount = pitzerOffset + _frzcPI[h];
        for(int a = 0; a < _nAOccPI[h]; ++a){
            aOccEvals[aOccCount++] = aF0[INDEX(aCount, aCount)];
            for(int mu = 0; mu < _soPI[h]; ++mu)
                _aOccC[h][mu][a] = Ca[mu+soOffset][aCount];
            ++aCount;
        }
        for(int a = 0; a < _nAVirPI[h]; ++a){
            aVirEvals[aVirCount++] = aF0[INDEX(aCount, aCount)];
            for(int mu = 0; mu < _soPI[h]; ++mu)
                _aVirC[h][mu][a] = Ca[mu+soOffset][aCount];
            ++aCount;
        }
        for(int b = 0; b < _nBOccPI[h]; ++b){
            bOccEvals[bOccCount++] = bF0[INDEX(bCount, bCount)];
            for(int mu = 0; mu < _soPI[h]; ++mu)
                _bOccC[h][mu][b] = Cb[mu+soOffset][bCount];
            ++bCount;
        }
        for(int b = 0; b < _nBVirPI[h]; ++b){
            bVirEvals[bVirCount++] = bF0[INDEX(bCount, bCount)];
            for(int mu = 0; mu < _soPI[h]; ++mu)
                _bVirC[h][mu][b] = Cb[mu+soOffset][bCount];
            ++bCount;
        }
        pitzerOffset += _moPI[h];
        soOffset     += _soPI[h];
    }

//    fprintf(outfile, "All Alpha MOs\n");
//    print_mat(Ca, _nSo, _nMo, outfile);
//    fprintf(outfile, "All Beta MOs\n");
//    print_mat(Cb, _nSo, _nMo, outfile);
//    for(int h = 0; h < _nIrreps; ++h){
//        fprintf(outfile, "Alpha Occ MOs for Irrep %d\n", h);
//        print_mat(_aOccC[h], _soPI[h], _nAOccPI[h], outfile);
//        fprintf(outfile, "Alpha Vir MOs for Irrep %d\n", h);
//        print_mat(_aVirC[h], _soPI[h], _nAVirPI[h], outfile);
//        fprintf(outfile, "Beta Occ MOs for Irrep %d\n", h);
//        print_mat(_bOccC[h], _soPI[h], _nBOccPI[h], outfile);
//        fprintf(outfile, "Beta Vir MOs for Irrep %d\n", h);
//        print_mat(_bVirC[h], _soPI[h], _nBVirPI[h], outfile);
//    }

    free_block(Ca);
    free_block(Cb);

    /*
     * Construct the F0 matrices, which are the standard Fock matrices with the
     * diagonal elements removed
     */
    //Alpha Occupied
    dpd_file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "F0 <O|O>");
    dpd_file2_mat_init(&F);
    int offset = 0;
    for(int h = 0; h < _nIrreps; ++h){
        offset += _frzcPI[h];
        for(int i = 0 ; i < _nAOccPI[h]; ++i){
            for(int j = 0 ; j < _nAOccPI[h]; ++j){
                F.matrix[h][i][j] = (i==j ? 0.0 : aF0[INDEX((i+offset), (j+offset))]);
            }
        }
        offset += _moPI[h];
    }
    dpd_file2_mat_wrt(&F);
    dpd_file2_close(&F);

    //Beta Occupied
    dpd_file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('o'), "F0 <o|o>");
    dpd_file2_mat_init(&F);
    offset = 0;
    for(int h = 0; h < _nIrreps; ++h){
        offset += _frzcPI[h];
        for(int i = 0 ; i < _nBOccPI[h]; ++i){
            for(int j = 0 ; j < _nBOccPI[h]; ++j){
                F.matrix[h][i][j] = (i==j ? 0.0 : bF0[INDEX((i+offset), (j+offset))]);
            }
        }
        offset += _moPI[h];
    }
    dpd_file2_mat_wrt(&F);
    dpd_file2_close(&F);

    //Alpha Virtual
    dpd_file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "F0 <V|V>");
    dpd_file2_mat_init(&F);
    offset = 0;
    for(int h = 0; h < _nIrreps; ++h){
        offset += _nAOccPI[h];
        for(int i = 0 ; i < _nAVirPI[h]; ++i){
            for(int j = 0 ; j < _nAVirPI[h]; ++j){
                F.matrix[h][i][j] = (i==j ? 0.0 : aF0[INDEX((i+offset), (j+offset))]);
            }
        }
        offset += _moPI[h] - _nAOccPI[h];
    }
    dpd_file2_mat_wrt(&F);
    dpd_file2_close(&F);

    //Beta Virtual
    dpd_file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('v'), ID('v'), "F0 <v|v>");
    dpd_file2_mat_init(&F);
    offset = 0;
    for(int h = 0; h < _nIrreps; ++h){
        offset += _nBOccPI[h];
        for(int i = 0 ; i < _nBVirPI[h]; ++i){
            for(int j = 0 ; j < _nBVirPI[h]; ++j){
                F.matrix[h][i][j] = (i==j ? 0.0 : bF0[INDEX((i+offset), (j+offset))]);
            }
        }
        offset += _moPI[h] - _nBOccPI[h];
    }
    dpd_file2_mat_wrt(&F);
    dpd_file2_close(&F);

    delete [] aF0;
    delete [] bF0;

//    fprintf(outfile, "occ a\n\t");
//    for(int n = 0; n < naocc; ++n) fprintf(outfile, "%6.2f ", occ_a[n]);
//    fprintf(outfile, "\nocc b\n\t");
//    for(int n = 0; n < nbocc; ++n) fprintf(outfile, "%6.2f ", occ_b[n]);
//    fprintf(outfile, "\nvir a\n\t");
//    for(int n = 0; n < navir; ++n) fprintf(outfile, "%6.2f ", vir_a[n]);
//    fprintf(outfile, "\nvir b\n\t");
//    for(int n = 0; n < nbvir; ++n) fprintf(outfile, "%6.2f ", vir_b[n]);
//    fprintf(outfile, "\n\n");

    ///////////////////////////////
    // The alpha-alpha spin case //
    ///////////////////////////////
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "D <OO|VV>");
    for(int h = 0; h < _nIrreps; ++h){
        dpd_buf4_mat_irrep_init(&D, h);
        for(int row = 0; row < D.params->rowtot[h]; ++row){
            int i = D.params->roworb[h][row][0];
            int j = D.params->roworb[h][row][1];
            for(int col = 0; col < D.params->coltot[h]; ++col){
                int a = D.params->colorb[h][col][0];
                int b = D.params->colorb[h][col][1];
                D.matrix[h][row][col] = 1.0/
                                    (_regularizer + aOccEvals[i] + aOccEvals[j] - aVirEvals[a] - aVirEvals[b]);
            }
        }
        dpd_buf4_mat_irrep_wrt(&D, h);
        dpd_buf4_mat_irrep_close(&D, h);
    }
    dpd_buf4_close(&D);

    //////////////////////////////
    // The alpha-beta spin case //
    //////////////////////////////
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "D <Oo|Vv>");
    for(int h = 0; h < _nIrreps; ++h){
        dpd_buf4_mat_irrep_init(&D, h);
        for(int row = 0; row < D.params->rowtot[h]; ++row){
            int i = D.params->roworb[h][row][0];
            int j = D.params->roworb[h][row][1];
            for(int col = 0; col < D.params->coltot[h]; ++col){
                int a = D.params->colorb[h][col][0];
                int b = D.params->colorb[h][col][1];
                D.matrix[h][row][col] = 1.0/
                                    (_regularizer + aOccEvals[i] + bOccEvals[j] - aVirEvals[a] - bVirEvals[b]);
            }
        }
        dpd_buf4_mat_irrep_wrt(&D, h);
        dpd_buf4_mat_irrep_close(&D, h);
    }
    dpd_buf4_close(&D);

    /////////////////////////////
    // The beta-beta spin case //
    /////////////////////////////
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "D <oo|vv>");
    for(int h = 0; h < _nIrreps; ++h){
        dpd_buf4_mat_irrep_init(&D, h);
        for(int row = 0; row < D.params->rowtot[h]; ++row){
            int i = D.params->roworb[h][row][0];
            int j = D.params->roworb[h][row][1];
            for(int col = 0; col < D.params->coltot[h]; ++col){
                int a = D.params->colorb[h][col][0];
                int b = D.params->colorb[h][col][1];
                D.matrix[h][row][col] = 1.0/
                                    (_regularizer + bOccEvals[i] + bOccEvals[j] - bVirEvals[a] - bVirEvals[b]);
            }
        }
        dpd_buf4_mat_irrep_wrt(&D, h);
        dpd_buf4_mat_irrep_close(&D, h);
    }
    dpd_buf4_close(&D);

    delete [] aOccEvals;
    delete [] bOccEvals;
    delete [] aVirEvals;
    delete [] bVirEvals;
}

}} // Namespaces




