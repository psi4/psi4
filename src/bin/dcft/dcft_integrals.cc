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
    if(print_ > 1){
        fprintf(outfile, "\tTransforming integrals...\n");
        fflush(outfile);
    }
    _ints->set_print(print_ - 2 >= 0 ? print_ - 2 : 0);

    // Generate the integrals in various spaces in chemists' notation
    _ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::vir);
    _ints->transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::occ, MOSpace::occ);
    _ints->transform_tei(MOSpace::vir, MOSpace::vir, MOSpace::occ, MOSpace::occ);
    _ints->transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::vir, MOSpace::vir);
    if((options_.get_str("AO_BASIS") == "NONE")){
        _ints->transform_tei(MOSpace::vir, MOSpace::vir, MOSpace::vir, MOSpace::vir);
    }

    // The integral object closes this file - we need to re-open it.
    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

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
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, psrq, ID("[O,v]"), ID("[o,V]"), "MO Ints <Ov|oV>");
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

    if(options_.get_str("AO_BASIS") == "NONE"){
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
    for(int h = 0; h < nirrep_; ++h){
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
    for(int h = 0; h < nirrep_; ++h){
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

    psio_->close(PSIF_LIBTRANS_DPD, 1);
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

    double *aOccEvals = new double [nalpha_];
    double *bOccEvals = new double [nbeta_];
    double *aVirEvals = new double [navir_];
    double *bVirEvals = new double [nbvir_];
    // Pick out the diagonal elements of the Fock matrix, making sure that they are in the order
    // used by the DPD library, i.e. starting from zero for each space and ordering by irrep
    int aOccCount = 0, bOccCount = 0, aVirCount = 0, bVirCount = 0;

    //Diagonal elements of the Fock matrix
    //Alpha spin
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0; i < naoccpi_[h]; ++i){
            aOccEvals[aOccCount++] = moFa_->get(h, i, i);
            for(int mu = 0; mu < nsopi_[h]; ++mu)
                aocc_c_->set(h, mu, i, Ca_->get(h, mu, i));
        }

        for(int a = 0; a < navirpi_[h]; ++a){
            aVirEvals[aVirCount++] = moFa_->get(h, naoccpi_[h] + a, naoccpi_[h] + a);
            for(int mu = 0; mu < nsopi_[h]; ++mu)
                avir_c_->set(h, mu, a, Ca_->get(h, mu, naoccpi_[h] + a));
        }        
    }

    //Off-diagonal elements of the Fock matrix
    //Alpha occupied

    dpd_file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "F <O|O>");
    dpd_file2_mat_init(&F);
    int offset = 0;
    for(int h = 0; h < nirrep_; ++h){
        offset += frzcpi_[h];
        for(int i = 0 ; i < naoccpi_[h]; ++i){
            for(int j = 0 ; j < naoccpi_[h]; ++j){
                F.matrix[h][i][j] = moFa_->get(h, i, j);
            }
        }
        offset += nmopi_[h];
    }
    dpd_file2_mat_wrt(&F);
    dpd_file2_close(&F);

    //Alpha Virtual
    dpd_file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "F <V|V>");
    dpd_file2_mat_init(&F);
    offset = 0;
    for(int h = 0; h < nirrep_; ++h){
        offset += naoccpi_[h];
        for(int i = 0 ; i < navirpi_[h]; ++i){
            for(int j = 0 ; j < navirpi_[h]; ++j){
                F.matrix[h][i][j] = moFa_->get(h, i + naoccpi_[h], j + naoccpi_[h]);
            }
        }
        offset += nmopi_[h] - naoccpi_[h];
    }
    dpd_file2_mat_wrt(&F);
    dpd_file2_close(&F);

    //Diagonal elements of the Fock matrix
    //Beta spin
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0; i < nboccpi_[h]; ++i){
            bOccEvals[bOccCount++] = moFb_->get(h, i, i);
            for(int mu = 0; mu < nsopi_[h]; ++mu)
                bocc_c_->set(h, mu, i, Cb_->get(h, mu, i));
        }
        for(int a = 0; a < nbvirpi_[h]; ++a){
            bVirEvals[bVirCount++] = moFb_->get(h, nboccpi_[h] + a, nboccpi_[h] + a);
            for(int mu = 0; mu < nsopi_[h]; ++mu)
                bvir_c_->set(h, mu, a, Cb_->get(h, mu, nboccpi_[h] + a));
        }
    }

    //Off-diagonal elements of the Fock matrix
    //Beta Occupied

    dpd_file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('o'), "F <o|o>");
    dpd_file2_mat_init(&F);
    offset = 0;
    for(int h = 0; h < nirrep_; ++h){
        offset += frzcpi_[h];
        for(int i = 0 ; i < nboccpi_[h]; ++i){
            for(int j = 0 ; j < nboccpi_[h]; ++j){
                F.matrix[h][i][j] = moFb_->get(h, i, j);
            }
        }
        offset += nmopi_[h];
    }
    dpd_file2_mat_wrt(&F);
    dpd_file2_close(&F);

    //Beta Virtual
    dpd_file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('v'), ID('v'), "F <v|v>");
    dpd_file2_mat_init(&F);
    offset = 0;
    for(int h = 0; h < nirrep_; ++h){
        offset += nboccpi_[h];
        for(int i = 0 ; i < nbvirpi_[h]; ++i){
            for(int j = 0 ; j < nbvirpi_[h]; ++j){
                F.matrix[h][i][j] = moFb_->get(h, i + nboccpi_[h], j + nboccpi_[h]);
            }
        }
        offset += nmopi_[h] - nboccpi_[h];
    }
    dpd_file2_mat_wrt(&F);
    dpd_file2_close(&F);

    ///////////////////////////////
    // The alpha-alpha spin case //
    ///////////////////////////////
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "D <OO|VV>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&D, h);
        for(int row = 0; row < D.params->rowtot[h]; ++row){
            int i = D.params->roworb[h][row][0];
            int j = D.params->roworb[h][row][1];
            for(int col = 0; col < D.params->coltot[h]; ++col){
                int a = D.params->colorb[h][col][0];
                int b = D.params->colorb[h][col][1];
                D.matrix[h][row][col] = 1.0/
                                    (regularizer_ + aOccEvals[i] + aOccEvals[j] - aVirEvals[a] - aVirEvals[b]);
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
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&D, h);
        for(int row = 0; row < D.params->rowtot[h]; ++row){
            int i = D.params->roworb[h][row][0];
            int j = D.params->roworb[h][row][1];
            for(int col = 0; col < D.params->coltot[h]; ++col){
                int a = D.params->colorb[h][col][0];
                int b = D.params->colorb[h][col][1];
                D.matrix[h][row][col] = 1.0/
                                    (regularizer_ + aOccEvals[i] + bOccEvals[j] - aVirEvals[a] - bVirEvals[b]);
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
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&D, h);
        for(int row = 0; row < D.params->rowtot[h]; ++row){
            int i = D.params->roworb[h][row][0];
            int j = D.params->roworb[h][row][1];
            for(int col = 0; col < D.params->coltot[h]; ++col){
                int a = D.params->colorb[h][col][0];
                int b = D.params->colorb[h][col][1];
                D.matrix[h][row][col] = 1.0/
                                    (regularizer_ + bOccEvals[i] + bOccEvals[j] - bVirEvals[a] - bVirEvals[b]);
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




