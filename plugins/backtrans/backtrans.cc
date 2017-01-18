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

#include "psi4-dec.h"
#include <libmints/mints.h>
#include "libciomr/libciomr.h"
#include <liboptions/liboptions.h>
#include "libdpd/dpd.h"
#include "psifiles.h"
#include "libiwl/iwl.h"
#include "libmints/mints.h"
#include "libtrans/integraltransform.h"
#include "libiwl/iwl.hpp"
#include "libplugin/plugin.h"
#include "psifiles.h"

#define ID(x) ints.DPD_ID(x)
INIT_PLUGIN

using namespace boost;

namespace psi{ namespace backtrans{

extern "C" int
read_options(std::string name, Options &options)
{
    if (name == "BACKTRANS"|| options.read_globals()) {
        /*- The amount of information printed
            to the output file -*/
        options.add_int("PRINT", 1);
    }

    return true;
}

extern "C"
SharedWavefunction backtrans(SharedWavefunction wfn, Options &options)
{
    dpdbuf4 I, G;
    int print = options.get_int("PRINT");

    boost::shared_ptr<PSIO> psio(_default_psio_lib_);


    std::vector<boost::shared_ptr<MOSpace> > spaces;
  
    std::string reference = options.get_str("REFERENCE");

    spaces.push_back(MOSpace::all);
    IntegralTransform ints(wfn,
                           spaces,
                           reference == "RHF" ? IntegralTransform::Restricted : IntegralTransform::Unrestricted,
                           IntegralTransform::DPDOnly,
                           IntegralTransform::QTOrder,
                           IntegralTransform::None,
                           false);
    dpd_set_default(ints.get_dpd_id());
    ints.set_print(print);
    ints.set_keep_dpd_so_ints(1);
    ints.set_keep_iwl_so_ints(1);
    ints.set_write_dpd_so_tpdm(true);
    ints.initialize();
    ints.transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);
    ints.backtransform_density();

    psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    psio->open(PSIF_TPDM_PRESORT, PSIO_OPEN_OLD);
    psio->open(PSIF_SO_PRESORT, PSIO_OPEN_OLD);
    psio->open(PSIF_AO_TPDM, PSIO_OPEN_OLD);


    char **labels = wfn->molecule()->irrep_labels();
    int nso       = wfn->nso();
    int nmo       = wfn->nmo();
    int nIrreps   = wfn->nirrep();
    int *orbspi   = wfn->nmopi();
    int *frzcpi   = wfn->frzcpi();
    int *frzvpi   = wfn->frzvpi();
    int *clsdpi   = wfn->doccpi();
    double eNuc   = wfn->molecule()->nuclear_repulsion_energy();
    double eOne   = 0.0;
    double eTwo   = 0.0;

    outfile->Printf("\n\n\t Irrep  frzcpi doccpi frzvpi  sopi\n");
    for(int h = 0; h < nIrreps; ++h){
        outfile->Printf("\t  %3s    %3d    %3d    %3d    %3d\n",
                labels[h], frzcpi[h], clsdpi[h], frzvpi[h], orbspi[h]);
    }
    int nTriMo = nmo * (nmo + 1) / 2;
    int nTriSo = nso * (nso + 1) / 2;
    double *temp = new double[nTriSo];
 
    /*
     * The SO basis results
     */
    Matrix soT(nIrreps, orbspi, orbspi);
    Matrix soV(nIrreps, orbspi, orbspi);
    IWL::read_one(psio.get(), PSIF_OEI, PSIF_SO_T, temp, nTriSo, 0, 0, "outfile");
    soT.set(temp);
    IWL::read_one(psio.get(), PSIF_OEI, PSIF_SO_V, temp, nTriSo, 0, 0, "outfile");
    soV.set(temp);
    Matrix soOei("SO OEI", nIrreps, orbspi, orbspi);
    soOei.add(soT);
    soOei.add(soV);
    Matrix soOpdm("SO OPDM", nIrreps, orbspi, orbspi);
    psio->open(PSIF_AO_OPDM, PSIO_OPEN_OLD);
    psio->read_entry(PSIF_AO_OPDM, "SO-basis OPDM", (char *) temp, sizeof(double)*nTriSo);
    psio->close(PSIF_AO_OPDM, 1);
    soOpdm.set(temp);
    if(print > 4){
        soOpdm.print();
        soOei.print();
    }
    eOne  = eNuc + soOpdm.vector_dot(soOei);

    global_dpd_->buf4_init(&G, PSIF_AO_TPDM, 0, ID("[n,n]"), ID("[n,n]"),
                  ID("[n>=n]+"), ID("[n>=n]+"), 0, "SO Basis TPDM (nn|nn)");
    global_dpd_->buf4_init(&I, PSIF_SO_PRESORT, 0, ID("[n,n]"), ID("[n,n]"),
                  ID("[n>=n]+"), ID("[n>=n]+"), 0, "SO Ints (nn|nn)");
    eTwo   = global_dpd_->buf4_dot(&I, &G);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&G);

    outfile->Printf("\n\tSO basis results\n");
    outfile->Printf("\tOne electron energy = %16.10f\n", eOne);
    outfile->Printf("\tTwo electron energy = %16.10f\n", eTwo);
    outfile->Printf("\tTotal energy        = %16.10f\n", eOne + eTwo);

    /*
     * MO basis results; these depend on spin case
     */
    const int *aPitzer = ints.alpha_corr_to_pitzer();
    const int *bPitzer = ints.beta_corr_to_pitzer(); // This is the same as the alpha array for RHF references
    if(print > 4)
        for(int n = 0; n < nmo; ++n)
            outfile->Printf("\tAlpha: %3d -> %3d   Beta: %3d -> %3d\n", n, aPitzer[n], n, bPitzer[n]);
    if(reference == "RHF"){
        Matrix moOei("MO OEI", nIrreps, orbspi, orbspi);
        IWL::read_one(psio.get(), PSIF_OEI, PSIF_MO_OEI, temp, nTriMo, 0, 0, "outfile");
        moOei.set(temp);
        Matrix moOpdm("MO OPDM", nIrreps, orbspi, orbspi);
      
      
        psio->open(PSIF_MO_OPDM, PSIO_OPEN_OLD);
        double **tempOPDM = block_matrix(nmo, nmo);
        psio->read_entry(PSIF_MO_OPDM, "MO-basis OPDM", (char *) tempOPDM[0], sizeof(double)*nmo*nmo);
        double *tempMo = new double[nTriMo];
        for(int p = 0; p < nmo; ++p){
          for(int q = 0; q <= p; ++q){
              int P = aPitzer[p];
              int Q = aPitzer[q];
              size_t PQ = INDEX(P,Q);
              tempMo[PQ] = 0.5 * (tempOPDM[p][q] + tempOPDM[q][p]);
          }
        }
        free_block(tempOPDM);
        psio->close(PSIF_MO_OPDM, 1);
        moOpdm.set(tempMo);
        if(print > 4){
            outfile->Printf("The MO basis OPDM, in Pitzer order\n");
            print_array(tempMo, nmo, "outfile");
            moOei.print();
        }
        delete[] tempMo;
      
        eOne = eNuc + moOpdm.vector_dot(moOei);
      
        global_dpd_->buf4_init(&G, PSIF_TPDM_PRESORT, 0, ID("[A,A]"), ID("[A,A]"),
                      ID("[A>=A]+"), ID("[A>=A]+"),0, "MO TPDM (AA|AA)");
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[A,A]"), ID("[A,A]"),
                      ID("[A>=A]+"), ID("[A>=A]+"), 0, "MO Ints (AA|AA)");
        eTwo   = global_dpd_->buf4_dot(&I, &G);
        global_dpd_->buf4_close(&I);
        global_dpd_->buf4_close(&G);
   }else{
        Matrix aMoOei("Alpha MO OEI", nIrreps, orbspi, orbspi);
        IWL::read_one(psio.get(), PSIF_OEI, PSIF_MO_A_OEI, temp, nTriMo, 0, 0, "outfile");
        aMoOei.set(temp);
        Matrix bMoOei("Alpha MO OEI", nIrreps, orbspi, orbspi);
        IWL::read_one(psio.get(), PSIF_OEI, PSIF_MO_B_OEI, temp, nTriMo, 0, 0, "outfile");
        bMoOei.set(temp);
      
        psio->open(PSIF_MO_OPDM, PSIO_OPEN_OLD);
        double **tempOPDM = block_matrix(nmo, nmo);
        Matrix aMoOpdm("Alpha MO OPDM", nIrreps, orbspi, orbspi);
        psio->read_entry(PSIF_MO_OPDM, "MO-basis Alpha OPDM", (char *) tempOPDM[0], sizeof(double)*nmo*nmo);
        double *tempMo = new double[nTriMo];
        for(int p = 0; p < nmo; ++p){
          for(int q = 0; q <= p; ++q){
              int P = aPitzer[p];
              int Q = aPitzer[q];
              size_t PQ = INDEX(P,Q);
              tempMo[PQ] = 0.5 * (tempOPDM[p][q] + tempOPDM[q][p]);
          }
        }
        aMoOpdm.set(tempMo);
        if(print > 4){
            outfile->Printf("The Alpha MO basis OPDM, in Pitzer order\n");
            print_array(tempMo, nmo, "outfile");
            aMoOei.print();
        }
        Matrix bMoOpdm("Beta MO OPDM", nIrreps, orbspi, orbspi);
        psio->read_entry(PSIF_MO_OPDM, "MO-basis Beta OPDM", (char *) tempOPDM[0], sizeof(double)*nmo*nmo);
        for(int p = 0; p < nmo; ++p){
          for(int q = 0; q <= p; ++q){
              int P = bPitzer[p];
              int Q = bPitzer[q];
              size_t PQ = INDEX(P,Q);
              tempMo[PQ] = 0.5 * (tempOPDM[p][q] + tempOPDM[q][p]);
          }
        }
        bMoOpdm.set(tempMo);
        free_block(tempOPDM);
        psio->close(PSIF_MO_OPDM, 1);
        if(print > 4){
            outfile->Printf("The Beta MO basis OPDM, in Pitzer order\n");
            print_array(tempMo, nmo, "outfile");
            bMoOei.print();
        }
        delete[] tempMo;
      
        eOne = eNuc + aMoOpdm.vector_dot(aMoOei) + bMoOpdm.vector_dot(bMoOei);
     
        // Alpha-Alpha 
        global_dpd_->buf4_init(&G, PSIF_TPDM_PRESORT, 0, ID("[A,A]"), ID("[A,A]"),
                      ID("[A>=A]+"), ID("[A>=A]+"),0, "MO TPDM (AA|AA)");
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[A,A]"), ID("[A,A]"),
                      ID("[A>=A]+"), ID("[A>=A]+"), 0, "MO Ints (AA|AA)");
        eTwo = global_dpd_->buf4_dot(&I, &G);
        global_dpd_->buf4_close(&I);
        global_dpd_->buf4_close(&G);
        // Alpha-Beta
        global_dpd_->buf4_init(&G, PSIF_TPDM_PRESORT, 0, ID("[A,A]"), ID("[a,a]"),
                      ID("[A>=A]+"), ID("[a>=a]+"),0, "MO TPDM (AA|aa)");
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[A,A]"), ID("[a,a]"),
                      ID("[A>=A]+"), ID("[a>=a]+"), 0, "MO Ints (AA|aa)");
        eTwo += global_dpd_->buf4_dot(&I, &G);
        global_dpd_->buf4_close(&I);
        global_dpd_->buf4_close(&G);
        // Beta-Beta
        global_dpd_->buf4_init(&G, PSIF_TPDM_PRESORT, 0, ID("[a,a]"), ID("[a,a]"),
                      ID("[a>=a]+"), ID("[a>=a]+"),0, "MO TPDM (aa|aa)");
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[a,a]"), ID("[a,a]"),
                      ID("[a>=a]+"), ID("[a>=a]+"), 0, "MO Ints (aa|aa)");
        eTwo += global_dpd_->buf4_dot(&I, &G);
        global_dpd_->buf4_close(&I);
        global_dpd_->buf4_close(&G);
   }

    outfile->Printf("\n\tMO basis results\n");
    outfile->Printf("\tOne electron energy = %16.10f\n", eOne);
    outfile->Printf("\tTwo electron energy = %16.10f\n", eTwo);
    outfile->Printf("\tTotal energy        = %16.10f\n", eOne + eTwo);

    psio->close(PSIF_TPDM_PRESORT, 1);
    psio->close(PSIF_LIBTRANS_DPD, 1);
    psio->close(PSIF_AO_TPDM, 1);
    psio->close(PSIF_SO_PRESORT, 1);

    delete [] temp;

    return wfn;
}

}} // End Namespace
