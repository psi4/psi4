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

/*! \file
    \ingroup CCSORT
    \brief Enter brief description of file here
*/

#include "psi4/libdpd/dpd.h"
#include "psi4/psifiles.h"
#include "MOInfo.h"
#include "Params.h"

#include "psi4/liboptions/liboptions.h"
#include "psi4/libqt/qt.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libtrans/mospace.h"
#include "psi4/lib3index/3index.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/petitelist.h"
#include "psi4/libmints/basisset.h"
#include "ccwave.h"

namespace psi { namespace ccenergy {

void CCEnergyWavefunction::form_df_ints(Options &options, int **cachelist, int *cachefiles, dpd_file4_cache_entry *priority)
{
    /*
     * Set up the DF tensor machinery
     */
    std::shared_ptr<BasisSet> dfBasis = get_basisset("DF_BASIS_CC");
    int nocc = doccpi_.sum();
    int nvir = nmo_ - nocc;
    int nbf = basisset_->nbf();
    int nbf2 = nbf*nbf;
    SharedMatrix Ca = Ca_subset("AO");
    DFTensor dfints(basisset_, dfBasis, Ca, nocc, nvir, nocc, nvir, options);
    SharedMatrix Qao = dfints.Qso();

    // Build the symmetrization matrix for the RI basis
    PetiteList petite(dfBasis, integral_, false);
    SharedMatrix dfAOtoSO = petite.aotoso();
    const Dimension &soDim = AO2SO_->colspi();
    const Dimension &dfDim = dfAOtoSO->colspi();
    SharedMatrix symQao(new Matrix(nirrep_, (const int*)dfDim, nbf2));
    double **pQao = Qao->pointer();
    for(int h = 0; h < nirrep_; ++h){
        double **pAOSO = dfAOtoSO->pointer(h);
        double **pSymQao = symQao->pointer(h);
        int nQso = dfAOtoSO->coldim(h);
        int nQao = dfAOtoSO->rowdim(h);
        if(nQso == 0) continue;
        C_DGEMM('t', 'n', nQso, nbf2, nQao, 1.0, pAOSO[0], nQso, pQao[0], nbf2, 0.0, pSymQao[0], nbf2);
    }
    // We're done with the original AO basis integrals now
    Qao.reset();

    double *htints = new double[nbf2];

    /*
     * Set up the DPD machinery
     */
    std::vector<int*> aospaces;
    if(params_.ref == 2){ // UHF
        aospaces.push_back(moinfo_.aoccpi);
        aospaces.push_back(moinfo_.aocc_sym);
        aospaces.push_back(moinfo_.sopi);
        aospaces.push_back(moinfo_.sosym);
        aospaces.push_back(moinfo_.boccpi);
        aospaces.push_back(moinfo_.bocc_sym);
        aospaces.push_back(moinfo_.sopi);
        aospaces.push_back(moinfo_.sosym);
        aospaces.push_back(moinfo_.avirtpi);
        aospaces.push_back(moinfo_.avir_sym);
        aospaces.push_back(moinfo_.bvirtpi);
        aospaces.push_back(moinfo_.bvir_sym);
    }else{ // R(O)HF
        aospaces.push_back(moinfo_.occpi);
        aospaces.push_back(moinfo_.occ_sym);
        aospaces.push_back(moinfo_.sopi);
        aospaces.push_back(moinfo_.sosym);
        aospaces.push_back(moinfo_.virtpi);
        aospaces.push_back(moinfo_.vir_sym);
    }
    int *dforbspi = new int[moinfo_.nirreps];
    int *dummyorbspi = new int[moinfo_.nirreps];
    int count = 0;
    for(int h = 0; h < moinfo_.nirreps; ++h){
        dummyorbspi[h] = 0;
        int norb = dfAOtoSO->coldim(h);
        dforbspi[h] = norb;
        count += norb;
    }
    dummyorbspi[0] = 1;
    int *dforbsym = new int[count];
    int *dummyorbsym = new int[1];
    dummyorbsym[0] = 0;
    count = 0;
    for(int h = 0; h < moinfo_.nirreps; ++h)
        for(int orb = 0; orb < dforbspi[h]; ++orb)
            dforbsym[count++] = h;
    aospaces.push_back(dforbspi);
    aospaces.push_back(dforbsym);
    aospaces.push_back(dummyorbspi);
    aospaces.push_back(dummyorbsym);

    dpd_init(1, moinfo_.nirreps, params_.memory, 0, cachefiles, cachelist, NULL, aospaces.size()/2, aospaces);

    delete [] dforbspi;
    delete [] dforbsym;
    delete [] dummyorbspi;
    delete [] dummyorbsym;

    /*      The IDs of the spaces
     *   R(O)HF
     * Q,D = 43
     * O,V = 27

     *   UHF
     * Q,D 94
     */
    dpdbuf4 I;

    if(params_.ref == 2){
        throw PSIEXCEPTION("UHF Density fitting NYI");
    }else{
        // Transform the AO indices to the SO basis
        global_dpd_->buf4_init(&I, PSIF_CC_OEI, 0, 43, 5, 43, 8, 0, "B(Q|pq)");
        for(int h = 0; h < nirrep_; ++h){
            global_dpd_->buf4_mat_irrep_init(&I, h);

            double **pQao = symQao->pointer(h);
            for(int pq = 0; pq < I.params->rowtot[h]; ++pq){
                for(int Gr=0; Gr < nirrep_; Gr++) {
                    // Transform ( Q | AO AO ) -> ( Q | AO SO )
                    int Gs = h^Gr;
                    int nrows = nbf;
                    int ncols = soDim[Gs];
                    int nlinks = nbf;
                    int rs = I.col_offset[h][Gr];
                    double **pc4a = AO2SO_->pointer(Gs);
                    if(nrows && ncols && nlinks)
                        C_DGEMM('n', 'n', nrows, ncols, nlinks, 1.0, pQao[pq],
                                nlinks, pc4a[0], ncols, 0.0, htints, nbf);
                    // Transform ( Q | AO SO ) -> ( Q | SO SO )
                    nrows = soDim[Gr];
                    ncols = soDim[Gs];
                    nlinks = nbf;
                    rs = I.col_offset[h][Gr];
                    double **pc3a = AO2SO_->pointer(Gr);
                    if(nrows && ncols && nlinks)
                        C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0, pc3a[0], nrows,
                                htints, nbf, 0.0, &I.matrix[h][pq][rs], ncols);
                } /* Gr */
            } /* pq */
            global_dpd_->buf4_mat_irrep_wrt(&I, h);
            global_dpd_->buf4_mat_irrep_close(&I, h);
        }
        // Permute for fast AO basis contractions
        global_dpd_->buf4_sort(&I, PSIF_CC_OEI, rspq, 8, 43, "B(pq|Q)");

#if 1
        // (OV|Q)
        dpdbuf4 OV;
        global_dpd_->buf4_init(&OV, PSIF_CC_OEI, 0, 43, 27, 43, 27, 0, "B(Q|OV)");
        for(int h = 0; h < nirrep_; ++h){
            global_dpd_->buf4_mat_irrep_init(&OV, h);
            global_dpd_->buf4_mat_irrep_init(&I, h);
            global_dpd_->buf4_mat_irrep_rd(&I, h);

            for(int pq = 0; pq < I.params->rowtot[h]; ++pq){
                for(int Gr=0; Gr < nirrep_; Gr++) {
                    // Transform ( Q | SO SO ) -> ( Q | SO V )
                    int Gs = h^Gr;
                    int nrows = moinfo_.sopi[Gr];
                    int ncols = moinfo_.virtpi[Gs];
                    int nlinks = moinfo_.sopi[Gs];
                    int rs = I.col_offset[h][Gr];
                    double **pc4a = moinfo_.Cv[Gs];
                    if(nrows && ncols && nlinks)
                        C_DGEMM('n', 'n', nrows, ncols, nlinks, 1.0,  &I.matrix[h][pq][rs],
                                nlinks, pc4a[0], ncols, 0.0, htints, nbf);
                    // Transform ( Q | SO V ) -> ( Q | O V )
                    nrows = moinfo_.occpi[Gr];
                    ncols = moinfo_.virtpi[Gs];
                    nlinks = moinfo_.sopi[Gr];
                    rs = OV.col_offset[h][Gr];
                    double **pc3a = moinfo_.Co[Gr];
                    if(nrows && ncols && nlinks)
                        C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0, pc3a[0], nrows,
                                htints, nbf, 0.0, &OV.matrix[h][pq][rs], ncols);
                } /* Gr */
            } /* pq */
            global_dpd_->buf4_mat_irrep_wrt(&OV, h);
            global_dpd_->buf4_mat_irrep_close(&OV, h);
            global_dpd_->buf4_mat_irrep_close(&I, h);
        }
        // Permute for fast contractions
        global_dpd_->buf4_sort(&OV, PSIF_CC_OEI, rspq, 27, 43, "B(OV|Q)");
        global_dpd_->buf4_close(&OV);


        // This could be done using the half-transformed intermediate above, but with extra memory cost
        dpdbuf4 VV;
        global_dpd_->buf4_init(&VV, PSIF_CC_OEI, 0, 43, 10, 43, 13, 0, "B(Q|VV)");
        for(int h = 0; h < nirrep_; ++h){
            global_dpd_->buf4_mat_irrep_init(&VV, h);
            global_dpd_->buf4_mat_irrep_init(&I, h);
            global_dpd_->buf4_mat_irrep_rd(&I, h);

            for(int pq = 0; pq < I.params->rowtot[h]; ++pq){
                for(int Gr=0; Gr < nirrep_; Gr++) {
                    // Transform ( Q | SO SO ) -> ( Q | SO V )
                    int Gs = h^Gr;
                    int nrows = moinfo_.sopi[Gr];
                    int ncols = moinfo_.virtpi[Gs];
                    int nlinks = moinfo_.sopi[Gs];
                    int rs = I.col_offset[h][Gr];
                    double **pc4a = moinfo_.Cv[Gs];
                    if(nrows && ncols && nlinks)
                        C_DGEMM('n', 'n', nrows, ncols, nlinks, 1.0,  &I.matrix[h][pq][rs],
                                nlinks, pc4a[0], ncols, 0.0, htints, nbf);
                    // Transform ( Q | SO V ) -> ( Q | V V )
                    nrows = moinfo_.virtpi[Gr];
                    ncols = moinfo_.virtpi[Gs];
                    nlinks = moinfo_.sopi[Gr];
                    rs = VV.col_offset[h][Gr];
                    double **pc3a = moinfo_.Cv[Gr];
                    if(nrows && ncols && nlinks)
                        C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0, pc3a[0], nrows,
                                htints, nbf, 0.0, &VV.matrix[h][pq][rs], ncols);
                } /* Gr */
            } /* pq */
            global_dpd_->buf4_mat_irrep_wrt(&VV, h);
            global_dpd_->buf4_mat_irrep_close(&VV, h);
            global_dpd_->buf4_mat_irrep_close(&I, h);
        }
        // Permute for fast contractions
        global_dpd_->buf4_sort(&VV, PSIF_CC_OEI, rspq, 13, 43, "B(VV|Q)");

        global_dpd_->buf4_close(&VV);
        global_dpd_->buf4_close(&I);

#else
        // (OV|Q)
        dpdbuf4 VV, OV;
        global_dpd_->buf4_init(&VV, PSIF_CC_OEI, 0, 43, 10, 43, 13, 0, "B(Q|VV)");
        global_dpd_->buf4_init(&OV, PSIF_CC_OEI, 0, 43, 27, 43, 27, 0, "B(Q|OV)");
        for(int h = 0; h < nirreps; ++h){
            global_dpd_->buf4_mat_irrep_init(&OV, h);
            global_dpd_->buf4_mat_irrep_init(&VV, h);
            global_dpd_->buf4_mat_irrep_init(&I, h);
            global_dpd_->buf4_mat_irrep_rd(&I, h);

            for(int pq = 0; pq < I.params->rowtot[h]; ++pq){
                for(int Gr=0; Gr < nirreps; Gr++) {
                    // Transform ( Q | SO SO ) -> ( Q | SO V )
                    int Gs = h^Gr;
                    int nrows = moinfo.sopi[Gr];
                    int ncols = moinfo.virtpi[Gs];
                    int nlinks = moinfo.sopi[Gs];
                    int rs = I.col_offset[h][Gr];
                    double **pc4a = moinfo.Cv[Gs];
                    if(nrows && ncols && nlinks)
                        C_DGEMM('n', 'n', nrows, ncols, nlinks, 1.0,  &I.matrix[h][pq][rs],
                                nlinks, pc4a[0], ncols, 0.0, htints, nbf);
                    // Transform ( Q | SO V ) -> ( Q | V V )
                    nrows = moinfo.virtpi[Gr];
                    ncols = moinfo.virtpi[Gs];
                    nlinks = moinfo.sopi[Gr];
                    rs = VV.col_offset[h][Gr];
                    double **pc3a = moinfo.Cv[Gr];
                    if(nrows && ncols && nlinks)
                        C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0, pc3a[0], nrows,
                                htints, nbf, 0.0, &VV.matrix[h][pq][rs], ncols);
                    // Transform ( Q | SO V ) -> ( Q | O V )
                    nrows = moinfo.occpi[Gr];
                    ncols = moinfo.virtpi[Gs];
                    nlinks = moinfo.sopi[Gr];
                    rs = OV.col_offset[h][Gr];
                    pc3a = moinfo.Co[Gr];
                    if(nrows && ncols && nlinks)
                        C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0, pc3a[0], nrows,
                                htints, nbf, 0.0, &OV.matrix[h][pq][rs], ncols);
                } /* Gr */
            } /* pq */
            global_dpd_->buf4_mat_irrep_wrt(&VV, h);
            global_dpd_->buf4_mat_irrep_close(&VV, h);
            global_dpd_->buf4_mat_irrep_wrt(&OV, h);
            global_dpd_->buf4_mat_irrep_close(&OV, h);
            global_dpd_->buf4_mat_irrep_close(&I, h);
        }
        // Permute for fast contractions
        global_dpd_->buf4_sort(&OV, PSIF_CC_OEI, rspq, 27, 43, "B(OV|Q)");
        // Permute for fast contractions
        global_dpd_->buf4_sort(&VV, PSIF_CC_OEI, rspq, 13, 43, "B(VV|Q)");
        global_dpd_->buf4_close(&OV);
        global_dpd_->buf4_close(&VV);
        global_dpd_->buf4_close(&I);
#endif
        delete [] htints;
    }
#if 0
    // Test the resulting integrals
    dpdbuf4 Q1, Q2;
    global_dpd_->buf4_init(&I, PSIF_CC_OEI, 0, uSO, uSO, pSO, pSO, 0 , "DF (pq|rs)");
    global_dpd_->buf4_init(&Q1, PSIF_CC_OEI, 0, QD, uSO, QD, pSO, 0, "B(Q|pq)");
    global_dpd_->buf4_init(&Q2, PSIF_CC_OEI, 0, QD, uSO, QD, pSO, 0, "B(Q|pq)");
    global_dpd_->contract444(&Q1, &Q2, &I, 1, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&Q1);
    global_dpd_->buf4_close(&Q2);
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_CC_OEI, 0, uSO, uSO, pSO, pSO, 0 , "DF (pq|rs)");
    global_dpd_->buf4_print(&I, outfile, 1);
    global_dpd_->buf4_close(&I);
#endif

    dpd_set_default(0);
}


}} // Namespaces
