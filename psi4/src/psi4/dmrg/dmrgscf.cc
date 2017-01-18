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
#ifdef USING_CheMPS2

//#include <libplugin/plugin.h>
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/writer.h"
#include "psi4/libmints/molecule.h"
#include "psi4/psi4-dec.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/psifiles.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libmints/wavefunction.h"

#include "psi4/libmints/typedefs.h"
//Header above this comment contains typedef std::shared_ptr<psi::Matrix> SharedMatrix;
#include "psi4/libciomr/libciomr.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libfock/jk.h"
#include "psi4/libmints/writer_file_prefix.h"
//Header above allows to obtain "filename.moleculename" with psi::get_writer_file_prefix(std::string name)

#include <stdlib.h>
#include <iostream>
#include <fstream>

#include <chemps2/Irreps.h>
#include <chemps2/Problem.h>
#include <chemps2/CASSCF.h>
#include <chemps2/Initialize.h>
#include <chemps2/EdmistonRuedenberg.h>
#include <chemps2/CASPT2.h>
#include <chemps2/Lapack.h>

using namespace std;

// This allows us to be lazy in getting the spaces in DPD calls
#define ID(x) ints->DPD_ID(x)

//INIT_PLUGIN

namespace psi{ namespace dmrg{

int chemps2_groupnumber(const string SymmLabel){

    int SyGroup = 0;
    bool stopFindGN = false;
    const int magic_number_max_groups_chemps2 = 8;
    do {
        if ( SymmLabel.compare(CheMPS2::Irreps::getGroupName(SyGroup)) == 0 ){ stopFindGN = true; }
        else { SyGroup++; }
    } while (( !stopFindGN ) && ( SyGroup < magic_number_max_groups_chemps2 ));

    (*outfile) << "Psi4 symmetry group was found to be <" << SymmLabel.c_str() << ">." << endl;
    if ( SyGroup >= magic_number_max_groups_chemps2 ){
        (*outfile) << "CheMPS2 did not recognize this symmetry group name. CheMPS2 only knows:" << endl;
        for (int cnt=0; cnt<magic_number_max_groups_chemps2; cnt++){
            (*outfile) << "   <" << (CheMPS2::Irreps::getGroupName(cnt)).c_str() << ">" << endl;
        }
        throw PSIEXCEPTION("CheMPS2 did not recognize the symmetry group name!");
    }
    return SyGroup;

}


void buildJK(SharedMatrix MO_RDM, SharedMatrix MO_JK, SharedMatrix Cmat, std::shared_ptr<JK> myJK, std::shared_ptr<Wavefunction> wfn){

    const int nso    = wfn->nso();
    const int nmo    = wfn->nmo();
    const int nirrep = wfn->nirrep();
    int * nmopi = init_int_array(nirrep);
    int * nsopi = init_int_array(nirrep);
    for ( int h = 0; h < nirrep; ++h ){
        nmopi[h] = wfn->nmopi()[h];
        nsopi[h] = wfn->nsopi()[h];
    }

    // nso can be different from nmo
    SharedMatrix SO_RDM;     SO_RDM = SharedMatrix( new Matrix( "SO RDM",   nirrep, nsopi, nsopi ) );
    SharedMatrix Identity; Identity = SharedMatrix( new Matrix( "Identity", nirrep, nsopi, nsopi ) );
    SharedMatrix SO_JK;       SO_JK = SharedMatrix( new Matrix( "SO JK",    nirrep, nsopi, nsopi ) );
    SharedMatrix work;         work = SharedMatrix( new Matrix( "work",     nirrep, nsopi, nmopi ) );

    work->gemm(false, false, 1.0, Cmat, MO_RDM, 0.0);
    SO_RDM->gemm(false, true, 1.0, work, Cmat, 0.0);

    std::vector<SharedMatrix> & CL = myJK->C_left();
    CL.clear();
    CL.push_back( SO_RDM );

    std::vector<SharedMatrix> & CR = myJK->C_right();
    CR.clear();
    Identity->identity();
    CR.push_back( Identity );

    myJK->set_do_J(true);
    myJK->set_do_K(true);
    myJK->set_do_wK(false);
    myJK->compute();

    SO_JK->copy( myJK->K()[0] );
    SO_JK->scale( -0.5 );
    SO_JK->add( myJK->J()[0] );

    work->gemm(false, false, 1.0, SO_JK, Cmat, 0.0);
    MO_JK->gemm(true, false, 1.0, Cmat, work,  0.0);

}


void copyPSIMXtoCHEMPS2MX( SharedMatrix source, CheMPS2::DMRGSCFindices * iHandler, CheMPS2::DMRGSCFmatrix * target ){

    for (int irrep = 0; irrep < iHandler->getNirreps(); irrep++){
        for (int orb1 = 0; orb1 < iHandler->getNORB(irrep); orb1++){
            for (int orb2 = 0; orb2 < iHandler->getNORB(irrep); orb2++){
                target->set(irrep, orb1, orb2, source->get(irrep, orb1, orb2));
            }
        }
    }

}


/*void copyCHEMPS2MXtoPSIMX( CheMPS2::DMRGSCFmatrix * source, CheMPS2::DMRGSCFindices * iHandler, SharedMatrix target ){
    for (int irrep = 0; irrep < iHandler->getNirreps(); irrep++){
        for (int orb1 = 0; orb1 < iHandler->getNORB(irrep); orb1++){
            for (int orb2 = 0; orb2 < iHandler->getNORB(irrep); orb2++){
                target->set(irrep, orb1, orb2, source->get(irrep, orb1, orb2));
            }
        }
    }

}*/


void buildQmatOCC( CheMPS2::DMRGSCFmatrix * theQmatOCC, CheMPS2::DMRGSCFindices * iHandler, SharedMatrix MO_RDM, SharedMatrix MO_JK, SharedMatrix Cmat, std::shared_ptr<JK> myJK, std::shared_ptr<Wavefunction> wfn ){

    MO_RDM->zero();
    for (int irrep = 0; irrep < iHandler->getNirreps(); irrep++){
        for (int orb = 0; orb < iHandler->getNOCC(irrep); orb++){
            MO_RDM->set(irrep, orb, orb, 2.0);
        }
    }
    buildJK( MO_RDM, MO_JK, Cmat, myJK, wfn );
    copyPSIMXtoCHEMPS2MX( MO_JK, iHandler, theQmatOCC );

}


void buildQmatACT( CheMPS2::DMRGSCFmatrix * theQmatACT, CheMPS2::DMRGSCFindices * iHandler, double * DMRG1DM, SharedMatrix MO_RDM, SharedMatrix MO_JK, SharedMatrix Cmat, std::shared_ptr<JK> myJK, std::shared_ptr<Wavefunction> wfn ){

    MO_RDM->zero();
    const int nOrbDMRG = iHandler->getDMRGcumulative(iHandler->getNirreps());
    for (int irrep = 0; irrep < iHandler->getNirreps(); irrep++){
        const int NOCC = iHandler->getNOCC(irrep);
        const int shift = iHandler->getDMRGcumulative(irrep);
        for (int orb1 = 0; orb1 < iHandler->getNDMRG(irrep); orb1++){
            for (int orb2 = orb1; orb2 < iHandler->getNDMRG(irrep); orb2++){
                const double value = DMRG1DM[ shift + orb1 + nOrbDMRG * ( shift + orb2 ) ];
                MO_RDM->set(irrep, NOCC+orb1, NOCC+orb2, value);
                MO_RDM->set(irrep, NOCC+orb2, NOCC+orb1, value);
            }
        }
    }
    buildJK( MO_RDM, MO_JK, Cmat, myJK, wfn );
    copyPSIMXtoCHEMPS2MX( MO_JK, iHandler, theQmatACT );

}


SharedMatrix print_rdm_ao( CheMPS2::DMRGSCFindices * idx, double * DMRG1DM, SharedMatrix MO_RDM, SharedMatrix Cmat, std::shared_ptr<Wavefunction> wfn ){

    const int num_irreps = idx->getNirreps();
    const int tot_dmrg   = idx->getDMRGcumulative( num_irreps );
    MO_RDM->zero();

    for ( int irrep = 0; irrep < num_irreps; irrep++ ){

        const int NOCC  = idx->getNOCC( irrep );
        const int NACT  = idx->getNDMRG( irrep );
        const int shift = idx->getDMRGcumulative( irrep );

        for ( int occ = 0; occ < NOCC; occ++ ){
            MO_RDM->set( irrep, occ, occ, 2.0 );
        }

        for ( int orb1 = 0; orb1 < NACT; orb1++ ){
            for ( int orb2 = orb1; orb2 < NACT; orb2++ ){
                const double value = 0.5 * ( DMRG1DM[ shift + orb1 + tot_dmrg * ( shift + orb2 ) ] + DMRG1DM[ shift + orb2 + tot_dmrg * ( shift + orb1 ) ] );
                MO_RDM->set( irrep, NOCC + orb1, NOCC + orb2, value );
                MO_RDM->set( irrep, NOCC + orb2, NOCC + orb1, value );
            }
        }
    }

    const int nirrep = wfn->nirrep();
    int * nmopi = init_int_array(nirrep);
    int * nsopi = init_int_array(nirrep);
    for ( int h = 0; h < nirrep; ++h ){
        nmopi[h] = wfn->nmopi()[h];
        nsopi[h] = wfn->nsopi()[h];
    }
    const int nao = wfn->aotoso()->rowspi( 0 );

    SharedMatrix tfo;       tfo = SharedMatrix( new Matrix( num_irreps, nao, nmopi ) );
    SharedMatrix work;     work = SharedMatrix( new Matrix( num_irreps, nao, nmopi ) );
    SharedMatrix AO_RDM; AO_RDM = SharedMatrix( new Matrix( nao, nao ) );

     tfo->gemm( false, false, 1.0, wfn->aotoso(), Cmat,   0.0 );
    work->gemm( false, false, 1.0, tfo,           MO_RDM, 0.0 );

    for ( int ao_row = 0; ao_row < nao; ao_row++ ){
        for ( int ao_col = 0; ao_col < nao; ao_col++ ){
            double value = 0.0;
            for ( int irrep = 0; irrep < num_irreps; irrep++ ){
                for ( int mo = 0; mo < nmopi[ irrep ]; mo++ ){
                    value += work->get( irrep, ao_row, mo ) * tfo->get( irrep, ao_col, mo );
                }
            }
            AO_RDM->set( 0, ao_row, ao_col, value );
        }
    }

    return AO_RDM;

}


void buildHamDMRG( std::shared_ptr<IntegralTransform> ints, std::shared_ptr<MOSpace> Aorbs_ptr, CheMPS2::DMRGSCFmatrix * theTmatrix, CheMPS2::DMRGSCFmatrix * theQmatOCC, CheMPS2::DMRGSCFindices * iHandler, CheMPS2::Hamiltonian * HamDMRG, std::shared_ptr<PSIO> psio, std::shared_ptr<Wavefunction> wfn ){

    ints->update_orbitals();
    // Since we don't regenerate the SO ints, we don't call sort_so_tei, and the OEI are not updated !!!!!
    ints->transform_tei( Aorbs_ptr, Aorbs_ptr, Aorbs_ptr, Aorbs_ptr );
    dpd_set_default(ints->get_dpd_id());
    const int nirrep = wfn->nirrep();

    // Econstant and one-electron integrals
    double Econstant = wfn->molecule()->nuclear_repulsion_energy();
    for (int h = 0; h < iHandler->getNirreps(); h++){
        const int NOCC = iHandler->getNOCC(h);
        for (int froz = 0; froz < NOCC; froz++){
            Econstant += 2 * theTmatrix->get(h, froz, froz) + theQmatOCC->get(h, froz, froz);
        }
        const int shift = iHandler->getDMRGcumulative(h);
        const int NDMRG = iHandler->getNDMRG(h);
        for (int orb1 = 0; orb1 < NDMRG; orb1++){
            for (int orb2 = orb1; orb2 < NDMRG; orb2++){
                HamDMRG->setTmat( shift+orb1, shift+orb2, theTmatrix->get(h, NOCC+orb1, NOCC+orb2) + theQmatOCC->get(h, NOCC+orb1, NOCC+orb2) );
            }
        }
    }
    HamDMRG->setEconst( Econstant );

    // Two-electron integrals
    dpdbuf4 K;
    psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[S,S]"), ID("[S,S]"), ID("[S>=S]+"), ID("[S>=S]+"), 0, "MO Ints (SS|SS)");
    for(int h = 0; h < nirrep; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for(int pq = 0; pq < K.params->rowtot[h]; ++pq){
            const int p = K.params->roworb[h][pq][0];
            const int q = K.params->roworb[h][pq][1];
            for(int rs = 0; rs < K.params->coltot[h]; ++rs){
                const int r = K.params->colorb[h][rs][0];
                const int s = K.params->colorb[h][rs][1];
                HamDMRG->setVmat( p, r, q, s, K.matrix[h][pq][rs] );
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);
    psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

}

void buildTmatrix( CheMPS2::DMRGSCFmatrix * theTmatrix, CheMPS2::DMRGSCFindices * iHandler, std::shared_ptr<PSIO> psio, SharedMatrix Cmat, std::shared_ptr<Wavefunction> wfn ){

    const int nirrep = wfn->nirrep();
    const int nmo    = wfn->nmo();
    const int nTriMo = nmo * (nmo + 1) / 2;
    const int nso    = wfn->nso();
    const int nTriSo = nso * (nso + 1) / 2;
    int * mopi       = init_int_array(nirrep);
    int * sopi       = init_int_array(nirrep);
    for ( int h = 0; h < nirrep; ++h ){
        mopi[h] = wfn->nmopi()[h];
        sopi[h] = wfn->nsopi()[h];
    }
    double * work1   = new double[ nTriSo ];
    double * work2   = new double[ nTriSo ];
    IWL::read_one(psio.get(), PSIF_OEI, PSIF_SO_T, work1, nTriSo, 0, 0, "outfile");
    IWL::read_one(psio.get(), PSIF_OEI, PSIF_SO_V, work2, nTriSo, 0, 0, "outfile");
    for (int n = 0; n < nTriSo; n++){ work1[n] += work2[n]; }
    delete [] work2;

    SharedMatrix soOei; soOei = SharedMatrix( new Matrix("SO OEI", nirrep, sopi, sopi) );
    SharedMatrix half;   half = SharedMatrix( new Matrix(  "Half", nirrep, mopi, sopi) );
    SharedMatrix moOei; moOei = SharedMatrix( new Matrix("MO OEI", nirrep, mopi, mopi) );

    soOei->set( work1 );
    half->gemm(true, false, 1.0, Cmat, soOei, 0.0);
    moOei->gemm(false, false, 1.0, half, Cmat, 0.0);
    delete [] work1;

    copyPSIMXtoCHEMPS2MX( moOei, iHandler, theTmatrix );

}


void fillRotatedTEI_coulomb( std::shared_ptr<IntegralTransform> ints, std::shared_ptr<MOSpace> OAorbs_ptr, CheMPS2::DMRGSCFintegrals * theRotatedTEI, CheMPS2::DMRGSCFindices * iHandler, std::shared_ptr<PSIO> psio, std::shared_ptr<Wavefunction> wfn ){

    ints->update_orbitals();
    // Since we don't regenerate the SO ints, we don't call sort_so_tei, and the OEI are not updated !!!!!
    ints->transform_tei( OAorbs_ptr, OAorbs_ptr, MOSpace::all, MOSpace::all );
    dpd_set_default(ints->get_dpd_id());
    const int nirrep = wfn->nirrep();

    // Two-electron integrals
    dpdbuf4 K;
    psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    // To only process the permutationally unique integrals, change the ID("[A,A]") to ID("[A>=A]+")
    //global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[A,A]"), ID("[A,A]"), ID("[A>=A]+"), ID("[A>=A]+"), 0, "MO Ints (AA|AA)");
    //int buf4_init(dpdbuf4 *Buf, int inputfile, int irrep, int pqnum, int rsnum, int file_pqnum, int file_rsnum, int anti, const char *label);
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[Q,Q]"), ID("[A,A]"), ID("[Q>=Q]+"), ID("[A>=A]+"), 0, "MO Ints (QQ|AA)");
    for(int h = 0; h < nirrep; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for(int pq = 0; pq < K.params->rowtot[h]; ++pq){
            const int p = K.params->roworb[h][pq][0];
            const int q = K.params->roworb[h][pq][1];
            const int psym = K.params->psym[p];
            const int qsym = K.params->qsym[q];
            const int prel = p - K.params->poff[psym];
            const int qrel = q - K.params->qoff[qsym];
            for(int rs = 0; rs < K.params->coltot[h]; ++rs){
                const int r = K.params->colorb[h][rs][0];
                const int s = K.params->colorb[h][rs][1];
                const int rsym = K.params->rsym[r];
                const int ssym = K.params->ssym[s];
                const int rrel = r - K.params->roff[rsym];
                const int srel = s - K.params->soff[ssym];
                theRotatedTEI->set_coulomb( psym, qsym, rsym, ssym, prel, qrel, rrel, srel, K.matrix[h][pq][rs] );
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);
    psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

}


void fillRotatedTEI_exchange( std::shared_ptr<IntegralTransform> ints, std::shared_ptr<MOSpace> OAorbs_ptr, std::shared_ptr<MOSpace> Vorbs_ptr, CheMPS2::DMRGSCFintegrals * theRotatedTEI, CheMPS2::DMRGSCFindices * iHandler, std::shared_ptr<PSIO> psio ){

    ints->update_orbitals();
    ints->transform_tei( Vorbs_ptr, OAorbs_ptr, Vorbs_ptr, OAorbs_ptr );
    dpd_set_default(ints->get_dpd_id());

    // Two-electron integrals
    dpdbuf4 K;
    psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    // To only process the permutationally unique integrals, change the ID("[A,A]") to ID("[A>=A]+")
    //global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[A,A]"), ID("[A,A]"), ID("[A>=A]+"), ID("[A>=A]+"), 0, "MO Ints (AA|AA)");
    //int buf4_init(dpdbuf4 *Buf, int inputfile, int irrep, int pqnum, int rsnum, int file_pqnum, int file_rsnum, int anti, const char *label);
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[T,Q]"), ID("[T,Q]"), ID("[T,Q]"), ID("[T,Q]"), 0, "MO Ints (TQ|TQ)");
    for(int h = 0; h < iHandler->getNirreps(); ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for(int pq = 0; pq < K.params->rowtot[h]; ++pq){
            const int p = K.params->roworb[h][pq][0];
            const int q = K.params->roworb[h][pq][1];
            const int psym = K.params->psym[p];
            const int qsym = K.params->qsym[q];
            const int prel = p - K.params->poff[psym] + iHandler->getNOCC(psym) + iHandler->getNDMRG(psym);
            const int qrel = q - K.params->qoff[qsym];
            for(int rs = 0; rs < K.params->coltot[h]; ++rs){
                const int r = K.params->colorb[h][rs][0];
                const int s = K.params->colorb[h][rs][1];
                const int rsym = K.params->rsym[r];
                const int ssym = K.params->ssym[s];
                const int rrel = r - K.params->roff[rsym] + iHandler->getNOCC(rsym) + iHandler->getNDMRG(rsym);
                const int srel = s - K.params->soff[ssym];
                theRotatedTEI->set_exchange( qsym, ssym, psym, rsym, qrel, srel, prel, rrel, K.matrix[h][pq][rs] );
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);
    psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

}


void copyUNITARYtoPSIMX( CheMPS2::DMRGSCFunitary * unitary, CheMPS2::DMRGSCFindices * iHandler, SharedMatrix target ){

    for (int irrep = 0; irrep < iHandler->getNirreps(); irrep++){
        for (int orb1 = 0; orb1 < iHandler->getNORB(irrep); orb1++){
            for (int orb2 = 0; orb2 < iHandler->getNORB(irrep); orb2++){
                target->set( irrep, orb1, orb2, unitary->getBlock(irrep)[ orb1 + iHandler->getNORB(irrep) * orb2 ] );
            }
        }
    }

}


void update_WFNco( SharedMatrix orig_coeff, CheMPS2::DMRGSCFindices * iHandler, CheMPS2::DMRGSCFunitary * unitary, std::shared_ptr<Wavefunction> wfn, SharedMatrix work1, SharedMatrix work2 ){

    copyUNITARYtoPSIMX( unitary, iHandler, work2 );
    wfn->Ca()->gemm(false, true, 1.0, orig_coeff, work2, 0.0);
    wfn->Cb()->copy(wfn->Ca());

}


SharedWavefunction dmrg(SharedWavefunction wfn, Options& options)
{

    /* This plugin is able to perform a DMRG calculation in a molecular orbital active space. */

    /*******************************
     *   Environment information   *
     *******************************/
    std::shared_ptr<PSIO> psio(_default_psio_lib_); // Grab the global (default) PSIO object, for file I/O
    if (!wfn){ throw PSIEXCEPTION("SCF has not been run yet!"); }

    /*************************
     *   Fetch the options   *
     *************************/

    const int wfn_irrep               = options.get_int("DMRG_IRREP");
    const int wfn_multp               = options.get_int("DMRG_MULTIPLICITY");
    int * dmrg_states                 = options.get_int_array("DMRG_SWEEP_STATES");
    const int ndmrg_states            = options["DMRG_SWEEP_STATES"].size();
    double * dmrg_econv               = options.get_double_array("DMRG_SWEEP_ENERGY_CONV");
    const int ndmrg_econv             = options["DMRG_SWEEP_ENERGY_CONV"].size();
    int * dmrg_maxsweeps              = options.get_int_array("DMRG_SWEEP_MAX_SWEEPS");
    const int ndmrg_maxsweeps         = options["DMRG_SWEEP_MAX_SWEEPS"].size();
    double * dmrg_noiseprefactors     = options.get_double_array("DMRG_SWEEP_NOISE_PREFAC");
    const int ndmrg_noiseprefactors   = options["DMRG_SWEEP_NOISE_PREFAC"].size();
    double * dmrg_dvdson_rtol         = options.get_double_array("DMRG_SWEEP_DVDSON_RTOL");
    const int ndmrg_dvdson_rtol       = options["DMRG_SWEEP_DVDSON_RTOL"].size();
    const bool dmrg_print_corr        = options.get_bool("DMRG_PRINT_CORR");
    const bool mps_chkpt              = options.get_bool("DMRG_MPS_WRITE");
    int * frozen_docc                 = options.get_int_array("RESTRICTED_DOCC");
    int * active                      = options.get_int_array("ACTIVE");
    const double d_convergence        = options.get_double("DMRG_SCF_GRAD_THR");
    const bool dmrg_store_unit        = options.get_bool("DMRG_UNITARY_WRITE");
    const bool dmrg_do_diis           = options.get_bool("DMRG_DIIS");
    const double dmrg_diis_branch     = options.get_double("DMRG_SCF_DIIS_THR");
    const bool dmrg_store_diis        = options.get_bool("DMRG_DIIS_WRITE");
    const int dmrg_max_iter           = options.get_int("DMRG_SCF_MAX_ITER");
    const int dmrg_which_root         = options.get_int("DMRG_EXCITATION");
    const bool dmrg_state_avg         = options.get_bool("DMRG_SCF_STATE_AVG");
    const string dmrg_active_space    = options.get_str("DMRG_SCF_ACTIVE_SPACE");
    const bool dmrg_loc_random        = options.get_bool("DMRG_LOCAL_INIT");
    const bool dmrg_caspt2            = options.get_bool("DMRG_CASPT2_CALC");
    const string dmrg_caspt2_orb      = options.get_str("DMRG_CASPT2_ORBS");
    const bool PSEUDOCANONICAL        = ( dmrg_caspt2_orb.compare("PSEUDOCANONICAL") == 0 ) ? true : false;
    const double dmrg_ipea            = options.get_double("DMRG_CASPT2_IPEA");
    const double dmrg_imag_shift      = options.get_double("DMRG_CASPT2_IMAG");
    const bool dmrg_molden            = options.get_bool("DMRG_MOLDEN_WRITE");
    const bool dmrg_density_ao        = options.get_bool("DMRG_OPDM_AO_PRINT");
    const int dmrg_num_vec_diis       = CheMPS2::DMRGSCF_numDIISvecs;
    const std::string unitaryname     = psi::get_writer_file_prefix( wfn->molecule()->name() ) + ".unitary.h5";
    const std::string diisname        = psi::get_writer_file_prefix( wfn->molecule()->name() ) + ".DIIS.h5";

    /****************************************
     *   Check if the input is consistent   *
     ****************************************/

    const int SyGroup= chemps2_groupnumber( wfn->molecule()->sym_label() );
    const int nmo    = wfn->nmo();
    const int nirrep = wfn->nirrep();
    int * orbspi     = init_int_array(nirrep);
    int * docc       = init_int_array(nirrep);
    int * socc       = init_int_array(nirrep);
    for ( int h = 0; h < nirrep; ++h ){
        orbspi[h] = wfn->nmopi()[h];
        docc[h] = wfn->doccpi()[h];
        socc[h] = wfn->soccpi()[h];
    }
    if ( wfn_irrep<0 )                            { throw PSIEXCEPTION("Option DMRG_IRREP (integer) may not be smaller than zero!"); }
    if ( wfn_multp<1 )                            { throw PSIEXCEPTION("Option DMRG_MULTIPLICITY (integer) should be larger or equal to one: DMRG_MULTIPLICITY = (2S+1) >= 1 !"); }
    if ( ndmrg_states==0 )                        { throw PSIEXCEPTION("Option DMRG_SWEEP_STATES (integer array) should be set!"); }
    if ( ndmrg_econv==0 )                         { throw PSIEXCEPTION("Option DMRG_SWEEP_ENERGY_CONV (double array) should be set!"); }
    if ( ndmrg_maxsweeps==0 )                     { throw PSIEXCEPTION("Option DMRG_SWEEP_MAX_SWEEPS (integer array) should be set!"); }
    if ( ndmrg_noiseprefactors==0 )               { throw PSIEXCEPTION("Option DMRG_SWEEP_NOISE_PREFAC (double array) should be set!"); }
    if ( ndmrg_states!=ndmrg_econv )              { throw PSIEXCEPTION("Options DMRG_SWEEP_STATES (integer array) and DMRG_ECONV (double array) should contain the same number of elements!"); }
    if ( ndmrg_states!=ndmrg_maxsweeps )          { throw PSIEXCEPTION("Options DMRG_SWEEP_STATES (integer array) and DMRG_SWEEP_MAX_SWEEPS (integer array) should contain the same number of elements!"); }
    if ( ndmrg_states!=ndmrg_noiseprefactors )    { throw PSIEXCEPTION("Options DMRG_SWEEP_STATES (integer array) and DMRG_SWEEP_NOISE_PREFAC (double array) should contain the same number of elements!"); }
    if ( ndmrg_states!=ndmrg_dvdson_rtol )        { throw PSIEXCEPTION("Options DMRG_SWEEP_STATES (integer array) and DMRG_SWEEP_DVDSON_RTOL (double array) should contain the same number of elements!"); }
    if ( options["RESTRICTED_DOCC"].size() != nirrep ){ throw PSIEXCEPTION("Option RESTRICTED_DOCC (integer array) should contain as many elements as there are irreps!"); }
    if ( options["ACTIVE"].size()      != nirrep ){ throw PSIEXCEPTION("Option ACTIVE (integer array) should contain as many elements as there are irreps!"); }
    for ( int cnt=0; cnt<ndmrg_states; cnt++ ){
       if ( dmrg_states[cnt] < 2 ){
          throw PSIEXCEPTION("Entries in DMRG_SWEEP_STATES (integer array) should be larger than 1!");
       }
       if ( dmrg_econv[cnt] <= 0.0 ){
          throw PSIEXCEPTION("Entries in DMRG_ECONV (double array) should be positive!");
       }
       if ( dmrg_maxsweeps[cnt] < 1 ){
          throw PSIEXCEPTION("Entries in DMRG_SWEEP_MAX_SWEEPS (integer array) should be positive!");
       }
       if ( dmrg_dvdson_rtol[cnt] <= 0.0 ){
          throw PSIEXCEPTION("Entries in DMRG_SWEEP_DVDSON_RTOL (double array) should be positive!");
       }
    }
    if ( d_convergence<=0.0 )                     { throw PSIEXCEPTION("Option DMRG_SCF_GRAD_THR (double) must be larger than zero!"); }
    if ( dmrg_diis_branch<=0.0 )                  { throw PSIEXCEPTION("Option DMRG_SCF_DIIS_THR (double) must be larger than zero!"); }
    if ( dmrg_max_iter<1 )                        { throw PSIEXCEPTION("Option DMRG_SCF_MAX_ITER (integer) must be larger than zero!"); }
    if ( dmrg_which_root<0 )                      { throw PSIEXCEPTION("Option DMRG_EXCITATION (integer) must be larger than zero!"); }
    if (( dmrg_caspt2 ) && ( dmrg_ipea < 0.0 ))   { throw PSIEXCEPTION("Option DMRG_CASPT2_IPEA (double) must be larger than zero!"); }
    if (( dmrg_molden ) && (( dmrg_caspt2 ) && ( PSEUDOCANONICAL == false ))){
       throw PSIEXCEPTION("Conflicting options: the molden file requires pseudocanonical orbitals, and caspt2 is requested in the active space orbitals.");
    }

    /*******************************************
     *   Create a CheMPS2::ConvergenceScheme   *
     *******************************************/

    CheMPS2::Initialize::Init();
    CheMPS2::ConvergenceScheme * OptScheme = new CheMPS2::ConvergenceScheme( ndmrg_states );
    for (int cnt=0; cnt<ndmrg_states; cnt++){
       OptScheme->set_instruction( cnt, dmrg_states[cnt], dmrg_econv[cnt], dmrg_maxsweeps[cnt], dmrg_noiseprefactors[cnt], dmrg_dvdson_rtol[cnt] );
    }

    /******************************************************************************
     *   Print orbital information; check consistency of frozen_docc and active   *
     ******************************************************************************/

    int * nvirtual = new int[nirrep];
    bool virtualsOK = true;
    for (int cnt=0; cnt<nirrep; cnt++){
       nvirtual[cnt] = orbspi[cnt] - frozen_docc[cnt] - active[cnt];
       if ( nvirtual[cnt] < 0 ){ virtualsOK = false; }
    }
    (*outfile) << "wfn_irrep   = " << wfn_irrep << endl;
    (*outfile) << "wfn_multp   = " << wfn_multp << endl;
    (*outfile) << "numOrbitals = [ " << orbspi[0];
    for (int cnt=1; cnt<nirrep; cnt++){ (*outfile) << " , " << orbspi[cnt];      } (*outfile) << " ]" << endl;
    (*outfile) << "R(O)HF DOCC = [ " << docc[0];
    for (int cnt=1; cnt<nirrep; cnt++){ (*outfile) << " , " << docc[cnt];        } (*outfile) << " ]" << endl;
    (*outfile) << "R(O)HF SOCC = [ " << socc[0];
    for (int cnt=1; cnt<nirrep; cnt++){ (*outfile) << " , " << socc[cnt];        } (*outfile) << " ]" << endl;
    (*outfile) << "frozen_docc = [ " << frozen_docc[0];
    for (int cnt=1; cnt<nirrep; cnt++){ (*outfile) << " , " << frozen_docc[cnt]; } (*outfile) << " ]" << endl;
    (*outfile) << "active      = [ " << active[0];
    for (int cnt=1; cnt<nirrep; cnt++){ (*outfile) << " , " << active[cnt];      } (*outfile) << " ]" << endl;
    (*outfile) << "virtual     = [ " << nvirtual[0];
    for (int cnt=1; cnt<nirrep; cnt++){ (*outfile) << " , " << nvirtual[cnt];    } (*outfile) << " ]" << endl;
    if ( !virtualsOK ){ throw PSIEXCEPTION("For at least one irrep: frozen_docc[ irrep ] + active[ irrep ] > numOrbitals[ irrep ]!"); }

    /*******************************************
     *   Create another bit of DMRG preamble   *
     *******************************************/
    CheMPS2::DMRGSCFindices * iHandler = new CheMPS2::DMRGSCFindices(nmo, SyGroup, frozen_docc, active, nvirtual);
    CheMPS2::DMRGSCFunitary * unitary = new CheMPS2::DMRGSCFunitary(iHandler);
    CheMPS2::DIIS * theDIIS = NULL;
    CheMPS2::DMRGSCFintegrals * theRotatedTEI = new CheMPS2::DMRGSCFintegrals( iHandler );
    const int nOrbDMRG = iHandler->getDMRGcumulative(nirrep);
    double * DMRG1DM = new double[nOrbDMRG * nOrbDMRG];
    double * DMRG2DM = new double[nOrbDMRG * nOrbDMRG * nOrbDMRG * nOrbDMRG];
    CheMPS2::DMRGSCFmatrix * theFmatrix = new CheMPS2::DMRGSCFmatrix( iHandler ); theFmatrix->clear();
    CheMPS2::DMRGSCFmatrix * theQmatOCC = new CheMPS2::DMRGSCFmatrix( iHandler ); theQmatOCC->clear();
    CheMPS2::DMRGSCFmatrix * theQmatACT = new CheMPS2::DMRGSCFmatrix( iHandler ); theQmatACT->clear();
    CheMPS2::DMRGSCFmatrix * theTmatrix = new CheMPS2::DMRGSCFmatrix( iHandler ); theTmatrix->clear();
    CheMPS2::DMRGSCFwtilde * wmattilde  = new CheMPS2::DMRGSCFwtilde( iHandler );
    delete [] nvirtual;

    /***************************************************
     *   Create the active space Hamiltonian storage   *
     ***************************************************/

    int nElectrons = 0;
    for (int cnt=0; cnt<nirrep; cnt++){ nElectrons += 2 * docc[cnt] + socc[cnt]; }
    (*outfile) << "nElectrons  = " << nElectrons << endl;

    // Number of electrons in the active space
    int nDMRGelectrons = nElectrons;
    for (int cnt=0; cnt<nirrep; cnt++){ nDMRGelectrons -= 2 * frozen_docc[cnt]; }
    (*outfile) << "nEl. active = " << nDMRGelectrons << endl;

    // Create the CheMPS2::Hamiltonian --> fill later
    int * orbitalIrreps = new int[ nOrbDMRG ];
    int counterFillOrbitalIrreps = 0;
    for (int h=0; h<nirrep; h++){
       for (int cnt=0; cnt<active[h]; cnt++){ //Only the active space is treated with DMRG-SCF!
          orbitalIrreps[counterFillOrbitalIrreps] = h;
          counterFillOrbitalIrreps++;
       }
    }
    CheMPS2::Hamiltonian * HamDMRG = new CheMPS2::Hamiltonian(nOrbDMRG, SyGroup, orbitalIrreps);
    delete [] orbitalIrreps;

    /* Create the CheMPS2::Problem
       You can fill Ham later, as Problem only keeps a pointer to the Hamiltonian object.
       Since only doubly occupied frozen orbitals are allowed, wfn_multp and wfn_irrep do not change. */
    CheMPS2::Problem * Prob = new CheMPS2::Problem( HamDMRG , wfn_multp-1 , nDMRGelectrons , wfn_irrep );
    if ( !(Prob->checkConsistency()) ){ throw PSIEXCEPTION("CheMPS2::Problem : No Hilbert state vector compatible with all symmetry sectors!"); }
    Prob->SetupReorderD2h(); // Does nothing if group not d2h

    /**************************************
     *   Input is parsed and consistent   *
     *   Start with DMRG                  *
     **************************************/

    SharedMatrix work1; work1 = SharedMatrix( new Matrix("work1", nirrep, orbspi, orbspi) );
    SharedMatrix work2; work2 = SharedMatrix( new Matrix("work2", nirrep, orbspi, orbspi) );
    std::shared_ptr<JK> myJK; myJK = std::shared_ptr<JK>(new DiskJK(wfn->basisset(), options));
    myJK->set_cutoff(0.0);
    myJK->initialize();
    SharedMatrix orig_coeff; orig_coeff = SharedMatrix( new Matrix( wfn->Ca() ) );

    std::vector<int> OAorbs; // Occupied + active
    std::vector<int> Aorbs;  // Only active
    std::vector<int> Vorbs;  // Virtual
    std::vector<int> empty;
    for (int h = 0; h < iHandler->getNirreps(); h++){
       for (int orb = 0; orb < iHandler->getNOCC(h) + iHandler->getNDMRG(h); orb++){
          OAorbs.push_back( iHandler->getOrigNOCCstart(h) + orb );
       }
       for (int orb = 0; orb < iHandler->getNDMRG(h); orb++){
          Aorbs.push_back( iHandler->getOrigNDMRGstart(h) + orb );
       }
       for (int orb = 0; orb < iHandler->getNVIRT(h); orb++){
          Vorbs.push_back( iHandler->getOrigNVIRTstart(h) + orb );
       }
    }
    std::shared_ptr<MOSpace> OAorbs_ptr; OAorbs_ptr = std::shared_ptr<MOSpace>( new MOSpace( 'Q', OAorbs, empty ) );
    std::shared_ptr<MOSpace>  Aorbs_ptr;  Aorbs_ptr = std::shared_ptr<MOSpace>( new MOSpace( 'S',  Aorbs, empty ) );
    std::shared_ptr<MOSpace>  Vorbs_ptr;  Vorbs_ptr = std::shared_ptr<MOSpace>( new MOSpace( 'T',  Vorbs, empty ) );
    std::vector<std::shared_ptr<MOSpace> > spaces;
    spaces.push_back( OAorbs_ptr   );
    spaces.push_back(  Aorbs_ptr   );
    spaces.push_back(  Vorbs_ptr   );
    spaces.push_back( MOSpace::all );
    // CheMPS2 requires RHF or ROHF orbitals.
    std::shared_ptr<IntegralTransform> ints;
    ints = std::shared_ptr<IntegralTransform>( new IntegralTransform( wfn, spaces, IntegralTransform::Restricted ) );
    ints->set_keep_iwl_so_ints( true );
    ints->set_keep_dpd_so_ints( true );
    //ints->set_print(6);

    (*outfile) << "###########################################################" << endl;
    (*outfile) << "###                                                     ###" << endl;
    (*outfile) << "###                       DMRG-SCF                      ###" << endl;
    (*outfile) << "###                                                     ###" << endl;
    (*outfile) << "###            CheMPS2 by Sebastian Wouters             ###" << endl;
    (*outfile) << "###        https://github.com/SebWouters/CheMPS2        ###" << endl;
    (*outfile) << "###   Comput. Phys. Commun. 185 (6), 1501-1514 (2014)   ###" << endl;
    (*outfile) << "###                                                     ###" << endl;
    (*outfile) << "###########################################################" << endl;
    (*outfile) << endl;
    (*outfile) << "Number of variables in the x-matrix = " << unitary->getNumVariablesX() << endl;

    //Convergence variables
    double gradNorm = 1.0;
    double updateNorm = 1.0;
    double * theupdate = new double[ unitary->getNumVariablesX() ];
    for (int cnt=0; cnt<unitary->getNumVariablesX(); cnt++){ theupdate[cnt] = 0.0; }
    double * theDIISparameterVector = NULL;
    double Energy = 1e8;

    int theDIISvectorParamSize = 0;
    int maxlinsize = 0;
    for (int irrep=0; irrep<nirrep; irrep++){
        const int linsize_irrep = iHandler->getNORB(irrep);
        theDIISvectorParamSize += linsize_irrep*(linsize_irrep-1)/2;
        if (linsize_irrep>maxlinsize){ maxlinsize = linsize_irrep; }
    }

    const int nOrbDMRG_pow4    = nOrbDMRG * nOrbDMRG * nOrbDMRG * nOrbDMRG;
    const int unitary_worksize = 4 * maxlinsize * maxlinsize;
    const int sizeWorkMem      = ( nOrbDMRG_pow4 > unitary_worksize ) ? nOrbDMRG_pow4 : unitary_worksize;
    const int tot_dmrg_power6  = nOrbDMRG_pow4 * nOrbDMRG * nOrbDMRG;
    double * mem1 = new double[ sizeWorkMem ];
    double * mem2 = new double[ ( PSEUDOCANONICAL ) ? sizeWorkMem : max( sizeWorkMem, tot_dmrg_power6 ) ];

    CheMPS2::EdmistonRuedenberg * theLocalizer = NULL;
    if ( dmrg_active_space.compare("LOC")==0 ){ theLocalizer = new CheMPS2::EdmistonRuedenberg( HamDMRG->getVmat(), iHandler->getGroupNumber() ); }

    //Load unitary from disk
    if ( dmrg_store_unit ){
        struct stat stFileInfo;
        int intStat = stat( unitaryname.c_str(), &stFileInfo );
        if (intStat==0){ unitary->loadU( unitaryname ); }
    }

    //Load DIIS from disk
    if (( dmrg_do_diis ) && ( dmrg_store_diis )){
        struct stat stFileInfo;
        int intStat = stat( diisname.c_str(), &stFileInfo );
        if (intStat==0){
            if (theDIIS == NULL){
                theDIIS = new CheMPS2::DIIS( theDIISvectorParamSize, unitary->getNumVariablesX(), dmrg_num_vec_diis );
                theDIISparameterVector = new double[ theDIISvectorParamSize ];
            }
            theDIIS->loadDIIS( diisname );
        }
    }

    int nIterations = 0;

    /*****************************
     ***   Actual DMRG loops   ***
     *****************************/
    while ((gradNorm > d_convergence) && (nIterations < dmrg_max_iter)){

        nIterations++;

        //Update the unitary transformation
        if (unitary->getNumVariablesX() > 0){

            std::ofstream capturing;
            std::streambuf * cout_buffer;
            string chemps2filename = outfile_name + ".chemps2";
            (*outfile) << "CheMPS2 output is temporarily written to the file " << chemps2filename << " and will be copied here." << endl;
            capturing.open( chemps2filename.c_str() , ios::trunc ); // truncate
            cout_buffer = cout.rdbuf( capturing.rdbuf() );

            unitary->updateUnitary(mem1, mem2, theupdate, true, true); //multiply = compact = true
            if (( dmrg_do_diis ) && ( updateNorm <= dmrg_diis_branch )){
                if ( dmrg_active_space.compare("NO")==0 ){
                    cout << "DIIS has started. Active space not rotated to NOs anymore!" << endl;
                }
                if ( dmrg_active_space.compare("LOC")==0 ){
                    cout << "DIIS has started. Active space not rotated to localized orbitals anymore!" << endl;
                }
                if (theDIIS == NULL){
                    theDIIS = new CheMPS2::DIIS( theDIISvectorParamSize, unitary->getNumVariablesX(), dmrg_num_vec_diis );
                    theDIISparameterVector = new double[ theDIISvectorParamSize ];
                    unitary->makeSureAllBlocksDetOne(mem1, mem2);
                }
                unitary->getLog(theDIISparameterVector, mem1, mem2);
                theDIIS->appendNew(theupdate, theDIISparameterVector);
                theDIIS->calculateParam(theDIISparameterVector);
                unitary->updateUnitary(mem1, mem2, theDIISparameterVector, false, false); //multiply = compact = false
            }

            cout.rdbuf(cout_buffer);
            capturing.close();
            std::ifstream copying;
            copying.open( chemps2filename , ios::in ); // read only
            if (copying.is_open()){
                string line;
                while( getline( copying, line ) ){ (*outfile) << line << endl; }
                copying.close();
            }
            system(("rm " + chemps2filename).c_str());

        }
        if (( dmrg_store_unit ) && (gradNorm!=1.0)){ unitary->saveU( unitaryname ); }
        if (( dmrg_store_diis ) && (updateNorm!=1.0) && (theDIIS!=NULL)){ theDIIS->saveDIIS( diisname ); }

        //Fill HamDMRG
        update_WFNco( orig_coeff, iHandler, unitary, wfn, work1, work2 );
        buildTmatrix( theTmatrix, iHandler, psio, wfn->Ca(), wfn );
        buildQmatOCC( theQmatOCC, iHandler, work1, work2, wfn->Ca(), myJK, wfn );
        buildHamDMRG( ints, Aorbs_ptr, theTmatrix, theQmatOCC, iHandler, HamDMRG, psio, wfn );

        //Localize the active space and reorder the orbitals within each irrep based on the exchange matrix
        if (( dmrg_active_space.compare("LOC")==0 ) && (theDIIS==NULL)){ //When the DIIS has started: stop

            std::ofstream capturing;
            std::streambuf * cout_buffer;
            string chemps2filename = outfile_name + ".chemps2";
            (*outfile) << "CheMPS2 output is temporarily written to the file " << chemps2filename << " and will be copied here." << endl;
            capturing.open( chemps2filename.c_str() , ios::trunc ); // truncate
            cout_buffer = cout.rdbuf( capturing.rdbuf() );

            theLocalizer->Optimize( mem1, mem2, dmrg_loc_random );
            theLocalizer->FiedlerExchange(maxlinsize, mem1, mem2);
            CheMPS2::CASSCF::fillLocalizedOrbitalRotations(theLocalizer->getUnitary(), iHandler, mem1);
            unitary->rotateActiveSpaceVectors(mem1, mem2);

            cout.rdbuf(cout_buffer);
            capturing.close();
            std::ifstream copying;
            copying.open( chemps2filename , ios::in ); // read only
            if (copying.is_open()){
                string line;
                while( getline( copying, line ) ){ (*outfile) << line << endl; }
                copying.close();
            }
            system(("rm " + chemps2filename).c_str());

            update_WFNco( orig_coeff, iHandler, unitary, wfn, work1, work2 );
            buildTmatrix( theTmatrix, iHandler, psio, wfn->Ca(), wfn );
            buildQmatOCC( theQmatOCC, iHandler, work1, work2, wfn->Ca(), myJK, wfn );
            buildHamDMRG( ints, Aorbs_ptr, theTmatrix, theQmatOCC, iHandler, HamDMRG, psio, wfn );
            (*outfile) << "Rotated the active space to localized orbitals, sorted according to the exchange matrix." << endl;

        }

        //Do the DMRG sweeps, and calculate the 2DM
        {
            std::ofstream capturing;
            std::streambuf * cout_buffer;
            string chemps2filename = outfile_name + ".chemps2";
            (*outfile) << "CheMPS2 output is temporarily written to the file " << chemps2filename << " and will be copied here." << endl;
            capturing.open( chemps2filename.c_str() , ios::trunc ); // truncate
            cout_buffer = cout.rdbuf( capturing.rdbuf() );

            for (int cnt = 0; cnt < nOrbDMRG_pow4; cnt++){ DMRG2DM[ cnt ] = 0.0; } //Clear the 2-RDM (to allow for state-averaged calculations)
            const string psi4TMPpath = PSIOManager::shared_object()->get_default_path();
            CheMPS2::DMRG * theDMRG = new CheMPS2::DMRG(Prob, OptScheme, mps_chkpt, psi4TMPpath);
            for (int state = -1; state < dmrg_which_root; state++){
                if (state > -1){ theDMRG->newExcitation( fabs( Energy ) ); }
                Energy = theDMRG->Solve();
                if ( dmrg_state_avg ){ // When SA-DMRGSCF: 2DM += current 2DM
                    theDMRG->calc2DMandCorrelations();
                    CheMPS2::CASSCF::copy2DMover( theDMRG->get2DM(), nOrbDMRG, DMRG2DM );
                }
                if ((state == -1) && (dmrg_which_root > 0)){ theDMRG->activateExcitations( dmrg_which_root ); }
            }
            if ( !(dmrg_state_avg) ){ // When SS-DMRGSCF: 2DM += last 2DM
                theDMRG->calc2DMandCorrelations();
                CheMPS2::CASSCF::copy2DMover( theDMRG->get2DM(), nOrbDMRG, DMRG2DM );
            }
            if ( dmrg_print_corr ){ theDMRG->getCorrelations()->Print(); }
            if ( CheMPS2::DMRG_storeRenormOptrOnDisk ){ theDMRG->deleteStoredOperators(); }
            delete theDMRG;
            if ((dmrg_state_avg) && (dmrg_which_root > 0)){
                const double averagingfactor = 1.0 / (dmrg_which_root+1);
                for (int cnt = 0; cnt < nOrbDMRG_pow4; cnt++){ DMRG2DM[ cnt ] *= averagingfactor; }
            }
            CheMPS2::CASSCF::setDMRG1DM( nDMRGelectrons, nOrbDMRG, DMRG1DM, DMRG2DM );

            cout.rdbuf(cout_buffer);
            capturing.close();
            std::ifstream copying;
            copying.open( chemps2filename , ios::in ); // read only
            if (copying.is_open()){
                string line;
                while( getline( copying, line ) ){ (*outfile) << line << endl; }
                copying.close();
            }
            system(("rm " + chemps2filename).c_str());
        }

        if (( dmrg_active_space.compare("NO")==0 ) && (theDIIS==NULL)){ //When the DIIS has started: stop
            CheMPS2::CASSCF::copy_active( DMRG1DM, theFmatrix, iHandler, true );
            CheMPS2::CASSCF::block_diagonalize( 'A', theFmatrix, unitary, mem1, mem2, iHandler, true, DMRG2DM, NULL, NULL ); // Unitary is updated and DMRG2DM rotated
            CheMPS2::CASSCF::setDMRG1DM( nDMRGelectrons, nOrbDMRG, DMRG1DM, DMRG2DM );
            update_WFNco( orig_coeff, iHandler, unitary, wfn, work1, work2 );
            buildTmatrix( theTmatrix, iHandler, psio, wfn->Ca(), wfn );
            buildQmatOCC( theQmatOCC, iHandler, work1, work2, wfn->Ca(), myJK, wfn );
            (*outfile) << "Rotated the active space to natural orbitals, sorted according to the NOON." << endl;
        }

        if (dmrg_max_iter == nIterations){
            if ( dmrg_store_unit ){ unitary->saveU( unitaryname ); }
            break;
        }

        buildQmatACT( theQmatACT, iHandler, DMRG1DM, work1, work2, wfn->Ca(), myJK, wfn );
        fillRotatedTEI_coulomb(  ints, OAorbs_ptr, theRotatedTEI, iHandler, psio, wfn );
        fillRotatedTEI_exchange( ints, OAorbs_ptr, Vorbs_ptr,  theRotatedTEI, iHandler, psio );

        {
            std::ofstream capturing;
            std::streambuf * cout_buffer;
            string chemps2filename = outfile_name + ".chemps2";
            (*outfile) << "CheMPS2 output is temporarily written to the file " << chemps2filename << " and will be copied here." << endl;
            capturing.open( chemps2filename.c_str() , ios::trunc ); // truncate
            cout_buffer = cout.rdbuf( capturing.rdbuf() );

            CheMPS2::CASSCF::buildFmat( theFmatrix, theTmatrix, theQmatOCC, theQmatACT, iHandler, theRotatedTEI, DMRG2DM, DMRG1DM);
            CheMPS2::CASSCF::buildWtilde(wmattilde, theTmatrix, theQmatOCC, theQmatACT, iHandler, theRotatedTEI, DMRG2DM, DMRG1DM);
            CheMPS2::CASSCF::augmentedHessianNR(theFmatrix, wmattilde, iHandler, unitary, theupdate, &updateNorm, &gradNorm);

            cout.rdbuf(cout_buffer);
            capturing.close();
            std::ifstream copying;
            copying.open( chemps2filename , ios::in ); // read only
            if (copying.is_open()){
                string line;
                while( getline( copying, line ) ){ (*outfile) << line << endl; }
                copying.close();
            }
            system(("rm " + chemps2filename).c_str());
        }
    }

    outfile->Printf("The DMRG-SCF energy = %3.10f \n", Energy);
    Process::environment.globals["CURRENT ENERGY"] = Energy;
    Process::environment.globals["DMRG-SCF ENERGY"] = Energy;

    if ((( dmrg_molden ) || (( dmrg_caspt2 ) && ( PSEUDOCANONICAL ))) && ( nIterations > 0 )){

        (*outfile) << "################################################" << endl;
        (*outfile) << "###                                          ###" << endl;
        (*outfile) << "###   Rotation to pseudocanonical orbitals   ###" << endl;
        (*outfile) << "###                                          ###" << endl;
        (*outfile) << "################################################" << endl;
        CheMPS2::CASSCF::construct_fock( theFmatrix, theTmatrix, theQmatOCC, theQmatACT, iHandler );
        CheMPS2::CASSCF::block_diagonalize( 'O', theFmatrix, unitary, mem1, mem2, iHandler, false, NULL, NULL, NULL );
        CheMPS2::CASSCF::block_diagonalize( 'A', theFmatrix, unitary, mem1, mem2, iHandler, false, DMRG2DM, NULL, NULL );
        CheMPS2::CASSCF::block_diagonalize( 'V', theFmatrix, unitary, mem1, mem2, iHandler, false, NULL, NULL, NULL );
        CheMPS2::CASSCF::setDMRG1DM( nDMRGelectrons, nOrbDMRG, DMRG1DM, DMRG2DM );
        update_WFNco( orig_coeff, iHandler, unitary, wfn, work1, work2 );
        buildTmatrix( theTmatrix, iHandler, psio, wfn->Ca(), wfn );
        buildQmatOCC( theQmatOCC, iHandler, work1, work2, wfn->Ca(), myJK, wfn );
        buildQmatACT( theQmatACT, iHandler, DMRG1DM, work1, work2, wfn->Ca(), myJK, wfn );
        CheMPS2::CASSCF::construct_fock( theFmatrix, theTmatrix, theQmatOCC, theQmatACT, iHandler );

    }

    if (( dmrg_molden ) && ( nIterations > 0 )){

        SharedVector occupation = SharedVector( new Vector( nirrep, orbspi ) );
        for ( int h = 0; h < nirrep; h++ ){
            const int NOCC = iHandler->getNOCC( h );
            const int NACT = iHandler->getNDMRG( h );
            const int NVIR = iHandler->getNVIRT( h );
            for ( int occ = 0; occ < NOCC; occ++ ){
               occupation->set( h, occ, 1.0 );
            }
            for ( int act = 0; act < NACT; act++ ){
                const double rdmvalue = DMRG1DM[ ( 1 + nOrbDMRG ) * ( iHandler->getDMRGcumulative( h ) + act ) ];
                occupation->set( h, NOCC + act, 0.5 * rdmvalue );
            }
            for ( int vir = 0; vir < NVIR; vir++ ){
                occupation->set( h, NOCC + NACT + vir, 0.0 );
            }
        }

        SharedVector sp_energies = SharedVector( new Vector( nirrep, orbspi ) );
        for ( int h = 0; h < nirrep; h++ ){
            const int NOCC = iHandler->getNOCC( h );
            const int NACT = iHandler->getNDMRG( h );
            const int NVIR = iHandler->getNVIRT( h );
            for ( int orb = 0; orb < orbspi[ h ]; orb++ ){
               sp_energies->set( h, orb, theFmatrix->get( h, orb, orb ) );
            }
        }

        std::shared_ptr<MoldenWriter> molden( new MoldenWriter( wfn ) );
        std::string filename = get_writer_file_prefix( wfn->molecule()->name() ) + ".pseudocanonical.molden";
        outfile->Printf( "Write molden file to %s. \n", filename.c_str() );
        molden->write( filename, wfn->Ca(), wfn->Ca(), sp_energies, sp_energies, occupation, occupation, true );

    }

    if ( dmrg_density_ao ){

        (*outfile) << "############################" << endl;
        (*outfile) << "###                      ###" << endl;
        (*outfile) << "###   DMRG 1-RDM in AO   ###" << endl;
        (*outfile) << "###                      ###" << endl;
        (*outfile) << "############################" << endl;
        (*outfile) << "Please check the molden file for AO basis function information." << endl;
        SharedMatrix AO_RDM = print_rdm_ao( iHandler, DMRG1DM, work1, wfn->Ca(), wfn );
        AO_RDM->print("outfile");

    }

    if (( dmrg_caspt2 ) && ( nIterations > 0 )){

        (*outfile) << "###########################" << endl;
        (*outfile) << "###                     ###" << endl;
        (*outfile) << "###     DMRG-CASPT2     ###" << endl;
        (*outfile) << "###                     ###" << endl;
        (*outfile) << "###########################" << endl;

        buildHamDMRG( ints, Aorbs_ptr, theTmatrix, theQmatOCC, iHandler, HamDMRG, psio, wfn );

        double * contract = new double[ tot_dmrg_power6 ];
        double * three_dm = new double[ tot_dmrg_power6 ];

        {
            std::ofstream capturing;
            std::streambuf * cout_buffer;
            string chemps2filename = outfile_name + ".chemps2";
            (*outfile) << "CheMPS2 output is temporarily written to the file " << chemps2filename << " and will be copied here." << endl;
            capturing.open( chemps2filename.c_str() , ios::trunc ); // truncate
            cout_buffer = cout.rdbuf( capturing.rdbuf() );

            for (int cnt = 0; cnt < nOrbDMRG_pow4; cnt++){ DMRG2DM[ cnt ] = 0.0; } //Clear the 2-RDM (to allow for state-averaged calculations)
            const string psi4TMPpath = PSIOManager::shared_object()->get_default_path();
            CheMPS2::DMRG * theDMRG = new CheMPS2::DMRG(Prob, OptScheme, false, psi4TMPpath); // Rotated orbital space --> do not use checkpoint
            for (int state = -1; state < dmrg_which_root; state++){
                if (state > -1){ theDMRG->newExcitation( fabs( Energy ) ); }
                const double E_CASSCF = theDMRG->Solve();
                if ((state == -1) && (dmrg_which_root > 0)){ theDMRG->activateExcitations( dmrg_which_root ); }
            }
            theDMRG->calc_rdms_and_correlations( true );
            CheMPS2::CASSCF::copy2DMover( theDMRG->get2DM(), nOrbDMRG, DMRG2DM  );        // 2-RDM
            CheMPS2::CASSCF::setDMRG1DM( nDMRGelectrons, nOrbDMRG, DMRG1DM, DMRG2DM );    // 1-RDM
            buildQmatACT( theQmatACT, iHandler, DMRG1DM, work1, work2, wfn->Ca(), myJK, wfn );
            CheMPS2::CASSCF::construct_fock( theFmatrix, theTmatrix, theQmatOCC, theQmatACT, iHandler );
            CheMPS2::CASSCF::copy_active( theFmatrix, mem2, iHandler );                   // Fock
            for ( int cnt = 0; cnt < tot_dmrg_power6; cnt++ ){ contract[ cnt ] = 0.0; }
            for ( int ham_orbz = 0; ham_orbz < nOrbDMRG; ham_orbz++ ){
               theDMRG->Symm4RDM( three_dm, ham_orbz, ham_orbz, false );
               int size = tot_dmrg_power6;
               double f_zz = 0.5 * mem2[ ham_orbz + nOrbDMRG * ham_orbz ];
               int inc1 = 1;
               daxpy_( &size, &f_zz, three_dm, &inc1, contract, &inc1 ); // trace( Fock * 4-RDM )
            }
            if ( PSEUDOCANONICAL == false ){
               for ( int ham_orb1 = 0; ham_orb1 < nOrbDMRG; ham_orb1++ ){
                  for ( int ham_orb2 = ham_orb1 + 1; ham_orb2 < nOrbDMRG; ham_orb2++ ){
                     if ( HamDMRG->getOrbitalIrrep( ham_orb1 ) == HamDMRG->getOrbitalIrrep( ham_orb2 ) ){
                        theDMRG->Symm4RDM( three_dm, ham_orb1, ham_orb2, false );
                        int size = tot_dmrg_power6;
                        double f_12 = 0.5 * ( mem2[ ham_orb1 + nOrbDMRG * ham_orb2 ] + mem2[ ham_orb2 + nOrbDMRG * ham_orb1 ] );
                        int inc1 = 1;
                        daxpy_( &size, &f_12, three_dm, &inc1, contract, &inc1 ); // trace( Fock * 4-RDM )
                     }
                  }
               }
               // CheMPS2::Cumulant::gamma4_fock_contract_ham( Prob, theDMRG->get3DM(), theDMRG->get2DM(), mem2, contract );
            }
            theDMRG->get3DM()->fill_ham_index( 1.0, false, three_dm, 0, nOrbDMRG ); // 3-RDM --> three_dm was used as work space for the constracted 4-RDM
            if (CheMPS2::DMRG_storeMpsOnDisk){        theDMRG->deleteStoredMPS();       }
            if (CheMPS2::DMRG_storeRenormOptrOnDisk){ theDMRG->deleteStoredOperators(); }
            delete theDMRG;

            cout.rdbuf(cout_buffer);
            capturing.close();
            std::ifstream copying;
            copying.open( chemps2filename , ios::in ); // read only
            if (copying.is_open()){
                string line;
                while( getline( copying, line ) ){ (*outfile) << line << endl; }
                copying.close();
            }
            system(("rm " + chemps2filename).c_str());
       }

       if ( PSEUDOCANONICAL == false ){
           (*outfile) << "CASPT2 : Deviation from pseudocanonical = " << CheMPS2::CASSCF::deviation_from_blockdiag( theFmatrix, iHandler ) << endl;
           CheMPS2::CASSCF::block_diagonalize( 'O', theFmatrix, unitary, mem1, mem2, iHandler, false, NULL, NULL, NULL );
           CheMPS2::CASSCF::block_diagonalize( 'A', theFmatrix, unitary, mem1, mem2, iHandler, false, DMRG2DM, three_dm, contract ); // 2-RDM, 3-RDM, and trace( Fock * cu(4)-4-RDM )
           CheMPS2::CASSCF::block_diagonalize( 'V', theFmatrix, unitary, mem1, mem2, iHandler, false, NULL, NULL, NULL );
           CheMPS2::CASSCF::setDMRG1DM( nDMRGelectrons, nOrbDMRG, DMRG1DM, DMRG2DM ); // 1-RDM
           update_WFNco( orig_coeff, iHandler, unitary, wfn, work1, work2 );
           buildTmatrix( theTmatrix, iHandler, psio, wfn->Ca(), wfn );
           buildQmatOCC( theQmatOCC, iHandler, work1, work2, wfn->Ca(), myJK, wfn );
           buildQmatACT( theQmatACT, iHandler, DMRG1DM, work1, work2, wfn->Ca(), myJK, wfn );
           CheMPS2::CASSCF::construct_fock( theFmatrix, theTmatrix, theQmatOCC, theQmatACT, iHandler ); // Fock
       }

       fillRotatedTEI_coulomb(  ints, OAorbs_ptr, theRotatedTEI, iHandler, psio, wfn );
       fillRotatedTEI_exchange( ints, OAorbs_ptr, Vorbs_ptr,  theRotatedTEI, iHandler, psio );

       (*outfile) << "CASPT2 : Norm F - F_pseudocan = " << CheMPS2::CASSCF::deviation_from_blockdiag( theFmatrix, iHandler ) << endl;
       double E_CASPT2 = 0.0;
       {
            std::ofstream capturing;
            std::streambuf * cout_buffer;
            string chemps2filename = outfile_name + ".chemps2";
            (*outfile) << "CheMPS2 output is temporarily written to the file " << chemps2filename << " and will be copied here." << endl;
            capturing.open( chemps2filename.c_str() , ios::trunc ); // truncate
            cout_buffer = cout.rdbuf( capturing.rdbuf() );

            CheMPS2::CASPT2 * myCASPT2 = new CheMPS2::CASPT2( iHandler, theRotatedTEI, theTmatrix, theFmatrix, DMRG1DM, DMRG2DM, three_dm, contract, dmrg_ipea );
            E_CASPT2 = myCASPT2->solve( dmrg_imag_shift );
            delete myCASPT2;

            cout.rdbuf(cout_buffer);
            capturing.close();
            std::ifstream copying;
            copying.open( chemps2filename , ios::in ); // read only
            if (copying.is_open()){
                string line;
                while( getline( copying, line ) ){ (*outfile) << line << endl; }
                copying.close();
            }
            system(("rm " + chemps2filename).c_str());

       }

       delete [] three_dm;
       delete [] contract;

       outfile->Printf("The DMRG-CASPT2 variational correction energy = %3.10f \n", E_CASPT2);
       outfile->Printf("The DMRG-CASPT2 energy = %3.10f \n", Energy + E_CASPT2);
       Process::environment.globals["CURRENT ENERGY"]    = Energy + E_CASPT2;
       Process::environment.globals["DMRG-CASPT2 ENERGY"] = Energy + E_CASPT2;

    }

    delete [] mem1;
    delete [] mem2;
    delete [] theupdate;
    if (theDIISparameterVector!=NULL){ delete [] theDIISparameterVector; }
    if (theLocalizer!=NULL){ delete theLocalizer; }
    if (theDIIS!=NULL){ delete theDIIS; }

    delete wmattilde;
    delete theTmatrix;
    delete theQmatOCC;
    delete theQmatACT;
    delete theFmatrix;
    delete [] DMRG1DM;
    delete [] DMRG2DM;
    delete theRotatedTEI;
    delete unitary;
    delete iHandler;

    delete OptScheme;
    delete Prob;
    delete HamDMRG;

    return wfn;
}

}} // End Namespaces
#endif
