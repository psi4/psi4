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

/*! \file
    \ingroup CCSORT
    \brief Enter brief description of file here 
*/

#include <libdpd/dpd.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#include "libmints/mints.h"
#include "liboptions/liboptions.h"
#include "libqt/qt.h"
#include "libtrans/integraltransform.h"
#include "libtrans/mospace.h"
#include "lib3index/3index.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

void form_df_ints(Options &options, int **cachelist, int *cachefiles, dpd_file4_cache_entry *priority)
{
    /*
     * Set up the DF tensor machinery
     */
    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
    boost::shared_ptr<Molecule> molecule = Process::environment.molecule();
    boost::shared_ptr<BasisSet> aoBasis = BasisSet::pyconstruct_orbital(molecule,
        "BASIS", options.get_str("BASIS"));
    boost::shared_ptr<BasisSet> dfBasis = BasisSet::pyconstruct_auxiliary(molecule,
        "DF_BASIS_CC", options.get_str("DF_BASIS_CC"), "RIFIT", options.get_str("BASIS"));
    int nmo = wfn->nmo();
    int nocc = wfn->doccpi().sum();
    int nvir = nmo - nocc;
    int nirreps = wfn->nirrep();
    int nbf = aoBasis->nbf();
    int nbf2 = nbf*nbf;
    SharedMatrix Ca = wfn->Ca_subset("AO");
    DFTensor dfints(aoBasis, dfBasis, Ca, nocc, nvir, nocc, nvir, options);
    SharedMatrix Qao = dfints.Qso();

    // Build the symmetrization matrix for the RI basis
    SharedMatrix aoAOtoSO = wfn->aotoso();
    PetiteList petite(dfBasis, wfn->integral(), false);
    SharedMatrix dfAOtoSO = petite.aotoso();
    const Dimension &soDim = aoAOtoSO->colspi();
    const Dimension &dfDim = dfAOtoSO->colspi();
    SharedMatrix symQao(new Matrix(nirreps, (const int*)dfDim, nbf2));
    double **pQao = Qao->pointer();
    for(int h = 0; h < nirreps; ++h){
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
    if(params.ref == 2){ // UHF
        aospaces.push_back(moinfo.aoccpi);
        aospaces.push_back(moinfo.aocc_sym);
        aospaces.push_back(moinfo.sopi);
        aospaces.push_back(moinfo.sosym);
        aospaces.push_back(moinfo.boccpi);
        aospaces.push_back(moinfo.bocc_sym);
        aospaces.push_back(moinfo.sopi);
        aospaces.push_back(moinfo.sosym);
        aospaces.push_back(moinfo.avirtpi);
        aospaces.push_back(moinfo.avir_sym);
        aospaces.push_back(moinfo.bvirtpi);
        aospaces.push_back(moinfo.bvir_sym);
    }else{ // R(O)HF
        aospaces.push_back(moinfo.occpi);
        aospaces.push_back(moinfo.occ_sym);
        aospaces.push_back(moinfo.sopi);
        aospaces.push_back(moinfo.sosym);
        aospaces.push_back(moinfo.virtpi);
        aospaces.push_back(moinfo.vir_sym);
    }
    int *dforbspi = new int[moinfo.nirreps];
    int *dummyorbspi = new int[moinfo.nirreps];
    int count = 0;
    for(int h = 0; h < moinfo.nirreps; ++h){
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
    for(int h = 0; h < moinfo.nirreps; ++h)
        for(int orb = 0; orb < dforbspi[h]; ++orb)
            dforbsym[count++] = h;
    aospaces.push_back(dforbspi);
    aospaces.push_back(dforbsym);
    aospaces.push_back(dummyorbspi);
    aospaces.push_back(dummyorbsym);

    dpd_init(1, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, NULL, aospaces.size()/2, aospaces);

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

    if(params.ref == 2){
        throw PSIEXCEPTION("UHF Density fitting NYI");
    }else{
        // Transform the AO indices to the SO basis
        global_dpd_->buf4_init(&I, PSIF_CC_OEI, 0, 43, 5, 43, 8, 0, "B(Q|pq)");
        for(int h = 0; h < nirreps; ++h){
            global_dpd_->buf4_mat_irrep_init(&I, h);

            double **pQao = symQao->pointer(h);
            for(int pq = 0; pq < I.params->rowtot[h]; ++pq){
                for(int Gr=0; Gr < nirreps; Gr++) {
                    // Transform ( Q | AO AO ) -> ( Q | AO SO )
                    int Gs = h^Gr;
                    int nrows = nbf;
                    int ncols = soDim[Gs];
                    int nlinks = nbf;
                    int rs = I.col_offset[h][Gr];
                    double **pc4a = aoAOtoSO->pointer(Gs);
                    if(nrows && ncols && nlinks)
                        C_DGEMM('n', 'n', nrows, ncols, nlinks, 1.0, pQao[pq],
                                nlinks, pc4a[0], ncols, 0.0, htints, nbf);
                    // Transform ( Q | AO SO ) -> ( Q | SO SO )
                    nrows = soDim[Gr];
                    ncols = soDim[Gs];
                    nlinks = nbf;
                    rs = I.col_offset[h][Gr];
                    double **pc3a = aoAOtoSO->pointer(Gr);
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
        for(int h = 0; h < nirreps; ++h){
            global_dpd_->buf4_mat_irrep_init(&OV, h);
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
                    // Transform ( Q | SO V ) -> ( Q | O V )
                    nrows = moinfo.occpi[Gr];
                    ncols = moinfo.virtpi[Gs];
                    nlinks = moinfo.sopi[Gr];
                    rs = OV.col_offset[h][Gr];
                    double **pc3a = moinfo.Co[Gr];
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
        for(int h = 0; h < nirreps; ++h){
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
