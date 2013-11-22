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
    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    boost::shared_ptr<BasisSet> aoBasis = BasisSet::construct(parser, molecule, "BASIS");
    boost::shared_ptr<BasisSet> dfBasis = BasisSet::construct(parser, molecule, "DF_BASIS_CC");
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
    }else{ // R(O)HF
        aospaces.push_back(moinfo.occpi);
        aospaces.push_back(moinfo.occ_sym);
        aospaces.push_back(moinfo.sopi);
        aospaces.push_back(moinfo.sosym);
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

    // The IDs of the spaces
    int QD = (params.ref == 2 ? 58 : 30);
    int uSO = 5;
    int pSO = 8;

    dpdbuf4 I;

    // Transform the AO indices to the SO basis
    global_dpd_->buf4_init(&I, PSIF_CC_OEI, 0, QD, uSO, QD, pSO, 0, "B(Q|pq)");
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
    global_dpd_->buf4_sort(&I, PSIF_CC_OEI, rspq, pSO, QD, "B(pq|Q)");
    global_dpd_->buf4_close(&I);
    delete [] htints;

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
