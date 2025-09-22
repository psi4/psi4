/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

// clang-format off
#ifdef USING_CheMPS2

//#include <libplugin/plugin.h>
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/writer.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/psifiles.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libmints/wavefunction.h"
#include "dmrg.h"

#include "psi4/libmints/typedefs.h"
//Header above this comment contains typedef std::shared_ptr<psi::Matrix> SharedMatrix;
#include "psi4/libciomr/libciomr.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libfock/jk.h"
#include "psi4/libmints/writer_file_prefix.h"
//Header above allows to obtain "filename.moleculename" with psi::get_writer_file_prefix(std::string name)

#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <fstream>

#include <chemps2/Irreps.h>
#include <chemps2/Problem.h>
#include <chemps2/CASSCF.h>
#include <chemps2/Initialize.h>
#include <chemps2/EdmistonRuedenberg.h>
#include <chemps2/CASPT2.h>
#include <chemps2/Lapack.h>

// This allows us to be lazy in getting the spaces in DPD calls
#define ID(x) ints.DPD_ID(x)

//INIT_PLUGIN

namespace psi{ namespace dmrg{

int chemps2_groupnumber(const std::string& SymmLabel){

    int SyGroup = 0;
    bool stopFindGN = false;
    const int magic_number_max_groups_chemps2 = 8;
    do {
        if ( SymmLabel.compare(CheMPS2::Irreps::getGroupName(SyGroup)) == 0 ){ stopFindGN = true; }
        else { SyGroup++; }
    } while (( !stopFindGN ) && ( SyGroup < magic_number_max_groups_chemps2 ));

    (*outfile->stream()) << "Psi4 symmetry group was found to be <" << SymmLabel.c_str() << ">." << std::endl;
    if ( SyGroup >= magic_number_max_groups_chemps2 ){
        (*outfile->stream()) << "CheMPS2 did not recognize this symmetry group name. CheMPS2 only knows:" << std::endl;
        for (int cnt=0; cnt<magic_number_max_groups_chemps2; cnt++){
            (*outfile->stream()) << "   <" << (CheMPS2::Irreps::getGroupName(cnt)).c_str() << ">" << std::endl;
        }
        throw PSIEXCEPTION("CheMPS2 did not recognize the symmetry group name!");
    }
    return SyGroup;

}


void buildJK(const Matrix& MO_RDM, Matrix& MO_JK, const Matrix& Cmat, JK& myJK){

    const auto& nsopi = Cmat.rowspi();

    auto Identity = std::make_shared<Matrix>( "Identity", nsopi, nsopi );
    Matrix SO_JK( "SO JK",    nsopi, nsopi );

    auto SO_RDM = MO_RDM.clone();
    SO_RDM->back_transform(Cmat);

    auto & CL = myJK.C_left();
    CL.clear();
    CL.push_back( SO_RDM );

    auto & CR = myJK.C_right();
    CR.clear();
    Identity->identity();
    CR.push_back( Identity );

    myJK.set_do_J(true);
    myJK.set_do_K(true);
    myJK.set_do_wK(false);
    myJK.compute();

    SO_JK.copy( myJK.K()[0] );
    SO_JK.scale( -0.5 );
    SO_JK.add( myJK.J()[0] );

    MO_JK.copy(SO_JK);
    MO_JK.transform(Cmat);

}


void copyPSIMXtoCHEMPS2MX(const Matrix& source, const CheMPS2::DMRGSCFindices& iHandler, CheMPS2::DMRGSCFmatrix& target ){

    for (int irrep = 0; irrep < iHandler.getNirreps(); irrep++){
        for (int orb1 = 0; orb1 < iHandler.getNORB(irrep); orb1++){
            for (int orb2 = 0; orb2 < iHandler.getNORB(irrep); orb2++){
                target.set(irrep, orb1, orb2, source.get(irrep, orb1, orb2));
            }
        }
    }

}


void buildQmatOCC(CheMPS2::DMRGSCFmatrix& theQmatOCC, const CheMPS2::DMRGSCFindices& iHandler, Matrix& MO_RDM, Matrix& MO_JK, const Matrix& Cmat, JK& myJK){

    MO_RDM.zero();
    for (int irrep = 0; irrep < iHandler.getNirreps(); irrep++){
        for (int orb = 0; orb < iHandler.getNOCC(irrep); orb++){
            MO_RDM.set(irrep, orb, orb, 2.0);
        }
    }
    buildJK(MO_RDM, MO_JK, Cmat, myJK);
    copyPSIMXtoCHEMPS2MX(MO_JK, iHandler, theQmatOCC );

}


void buildQmatACT(CheMPS2::DMRGSCFmatrix& theQmatACT, const CheMPS2::DMRGSCFindices& iHandler, const std::vector<double>& DMRG1DM, Matrix& MO_RDM, Matrix& MO_JK, const Matrix& Cmat, JK& myJK){

    MO_RDM.zero();
    const int nOrbDMRG = iHandler.getDMRGcumulative(iHandler.getNirreps());
    for (int irrep = 0; irrep < iHandler.getNirreps(); irrep++){
        const int NOCC = iHandler.getNOCC(irrep);
        const int shift = iHandler.getDMRGcumulative(irrep);
        for (int orb1 = 0; orb1 < iHandler.getNDMRG(irrep); orb1++){
            for (int orb2 = orb1; orb2 < iHandler.getNDMRG(irrep); orb2++){
                const double value = DMRG1DM[ shift + orb1 + nOrbDMRG * ( shift + orb2 ) ];
                MO_RDM.set(irrep, NOCC+orb1, NOCC+orb2, value);
                MO_RDM.set(irrep, NOCC+orb2, NOCC+orb1, value);
            }
        }
    }
    buildJK(MO_RDM, MO_JK, Cmat, myJK);
    copyPSIMXtoCHEMPS2MX(MO_JK, iHandler, theQmatACT);

}


Matrix build_rdm(const CheMPS2::DMRGSCFindices& idx, const std::vector<double>& DMRG1DM){

    const int num_irreps = idx.getNirreps();
    const int tot_dmrg   = idx.getDMRGcumulative( num_irreps );
    std::vector<int> dim(num_irreps);
    for (auto h = 0; h < dim.size(); h++) {
        dim[h] = idx.getNORB(h);
    }
    Dimension nmopi(dim);
    Matrix MO_RDM(nmopi, nmopi);

    for ( int irrep = 0; irrep < num_irreps; irrep++ ){

        const int NOCC  = idx.getNOCC( irrep );
        const int NACT  = idx.getNDMRG( irrep );
        const int shift = idx.getDMRGcumulative( irrep );

        for ( int occ = 0; occ < NOCC; occ++ ){
            MO_RDM.set( irrep, occ, occ, 2.0 );
        }

        for ( int orb1 = 0; orb1 < NACT; orb1++ ){
            for ( int orb2 = orb1; orb2 < NACT; orb2++ ){
                const double value = 0.5 * ( DMRG1DM[ shift + orb1 + tot_dmrg * ( shift + orb2 ) ] + DMRG1DM[ shift + orb2 + tot_dmrg * ( shift + orb1 ) ] );
                MO_RDM.set( irrep, NOCC + orb1, NOCC + orb2, value );
                MO_RDM.set( irrep, NOCC + orb2, NOCC + orb1, value );
            }
        }
    }

    return MO_RDM;
}

Matrix mo_to_ao(Matrix& MO_mat, const Matrix& Cmat, const Matrix& aotoso ){
    const auto& nmopi = Cmat.colspi();
    const int nao = aotoso.rowspi( 0 );

    Matrix AO_mat(nao, nao);

    const auto tfo = linalg::doublet(aotoso, Cmat, false, false);
    const auto work = linalg::doublet(tfo, MO_mat, false, false);

    for ( int ao_row = 0; ao_row < nao; ao_row++ ){
        for ( int ao_col = 0; ao_col < nao; ao_col++ ){
            double value = 0.0;
            for ( int irrep = 0; irrep < MO_mat.nirrep(); irrep++ ){
                for ( int mo = 0; mo < nmopi[ irrep ]; mo++ ){
                    value += work.get( irrep, ao_row, mo ) * tfo.get( irrep, ao_col, mo );
                }
            }
            AO_mat.set( 0, ao_row, ao_col, value );
        }
    }

    return AO_mat;

}


void buildHamDMRG(IntegralTransform& ints, std::shared_ptr<MOSpace> Aorbs_ptr, const CheMPS2::DMRGSCFmatrix& theTmatrix, const CheMPS2::DMRGSCFmatrix& theQmatOCC, const CheMPS2::DMRGSCFindices& iHandler, CheMPS2::Hamiltonian& HamDMRG, std::shared_ptr<PSIO> psio, const Wavefunction& wfn ){

    ints.update_orbitals();
    // Since we don't regenerate the SO ints, we don't call sort_so_tei, and the OEI are not updated !!!!!
    ints.transform_tei( Aorbs_ptr, Aorbs_ptr, Aorbs_ptr, Aorbs_ptr );
    dpd_set_default(ints.get_dpd_id());
    const int nirrep = wfn.nirrep();

    // Econstant and one-electron integrals
    double Econstant = wfn.molecule()->nuclear_repulsion_energy(wfn.get_dipole_field_strength());
    for (int h = 0; h < iHandler.getNirreps(); h++){
        const int NOCC = iHandler.getNOCC(h);
        for (int froz = 0; froz < NOCC; froz++){
            Econstant += 2 * theTmatrix.get(h, froz, froz) + theQmatOCC.get(h, froz, froz);
        }
        const int shift = iHandler.getDMRGcumulative(h);
        const int NDMRG = iHandler.getNDMRG(h);
        for (int orb1 = 0; orb1 < NDMRG; orb1++){
            for (int orb2 = orb1; orb2 < NDMRG; orb2++){
                HamDMRG.setTmat( shift+orb1, shift+orb2, theTmatrix.get(h, NOCC+orb1, NOCC+orb2) + theQmatOCC.get(h, NOCC+orb1, NOCC+orb2) );
            }
        }
    }
    HamDMRG.setEconst( Econstant );

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
                HamDMRG.setVmat( p, r, q, s, K.matrix[h][pq][rs] );
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);
    psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

}

void buildTmatrix(CheMPS2::DMRGSCFmatrix& theTmatrix, const CheMPS2::DMRGSCFindices& iHandler, const Matrix& Cmat, const Wavefunction& wfn){
    auto moOei = *wfn.mintshelper()->so_kinetic()->clone();
    moOei.add(wfn.mintshelper()->so_potential());
    moOei.transform(Cmat);

    copyPSIMXtoCHEMPS2MX( moOei, iHandler, theTmatrix);
}


void fillRotatedTEI_coulomb(IntegralTransform& ints, std::shared_ptr<MOSpace> OAorbs_ptr, CheMPS2::DMRGSCFintegrals& theRotatedTEI, const CheMPS2::DMRGSCFindices& iHandler, std::shared_ptr<PSIO> psio){

    ints.update_orbitals();
    // Since we don't regenerate the SO ints, we don't call sort_so_tei, and the OEI are not updated !!!!!
    ints.transform_tei( OAorbs_ptr, OAorbs_ptr, MOSpace::all, MOSpace::all );
    dpd_set_default(ints.get_dpd_id());

    // Two-electron integrals
    dpdbuf4 K;
    psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    // To only process the permutationally unique integrals, change the ID("[A,A]") to ID("[A>=A]+")
    //global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[A,A]"), ID("[A,A]"), ID("[A>=A]+"), ID("[A>=A]+"), 0, "MO Ints (AA|AA)");
    //int buf4_init(dpdbuf4 *Buf, int inputfile, int irrep, int pqnum, int rsnum, int file_pqnum, int file_rsnum, int anti, const char *label);
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[Q,Q]"), ID("[A,A]"), ID("[Q>=Q]+"), ID("[A>=A]+"), 0, "MO Ints (QQ|AA)");
    for(int h = 0; h < iHandler.getNirreps(); ++h){
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
                theRotatedTEI.set_coulomb( psym, qsym, rsym, ssym, prel, qrel, rrel, srel, K.matrix[h][pq][rs] );
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);
    psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

}


void fillRotatedTEI_exchange(IntegralTransform& ints, std::shared_ptr<MOSpace> OAorbs_ptr, std::shared_ptr<MOSpace> Vorbs_ptr, CheMPS2::DMRGSCFintegrals& theRotatedTEI, const CheMPS2::DMRGSCFindices& iHandler, std::shared_ptr<PSIO> psio ){

    ints.update_orbitals();
    ints.transform_tei( Vorbs_ptr, OAorbs_ptr, Vorbs_ptr, OAorbs_ptr );
    dpd_set_default(ints.get_dpd_id());

    // Two-electron integrals
    dpdbuf4 K;
    psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    // To only process the permutationally unique integrals, change the ID("[A,A]") to ID("[A>=A]+")
    //global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[A,A]"), ID("[A,A]"), ID("[A>=A]+"), ID("[A>=A]+"), 0, "MO Ints (AA|AA)");
    //int buf4_init(dpdbuf4 *Buf, int inputfile, int irrep, int pqnum, int rsnum, int file_pqnum, int file_rsnum, int anti, const char *label);
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[T,Q]"), ID("[T,Q]"), ID("[T,Q]"), ID("[T,Q]"), 0, "MO Ints (TQ|TQ)");
    for(int h = 0; h < iHandler.getNirreps(); ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for(int pq = 0; pq < K.params->rowtot[h]; ++pq){
            const int p = K.params->roworb[h][pq][0];
            const int q = K.params->roworb[h][pq][1];
            const int psym = K.params->psym[p];
            const int qsym = K.params->qsym[q];
            const int prel = p - K.params->poff[psym] + iHandler.getNOCC(psym) + iHandler.getNDMRG(psym);
            const int qrel = q - K.params->qoff[qsym];
            for(int rs = 0; rs < K.params->coltot[h]; ++rs){
                const int r = K.params->colorb[h][rs][0];
                const int s = K.params->colorb[h][rs][1];
                const int rsym = K.params->rsym[r];
                const int ssym = K.params->ssym[s];
                const int rrel = r - K.params->roff[rsym] + iHandler.getNOCC(rsym) + iHandler.getNDMRG(rsym);
                const int srel = s - K.params->soff[ssym];
                theRotatedTEI.set_exchange( qsym, ssym, psym, rsym, qrel, srel, prel, rrel, K.matrix[h][pq][rs] );
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);
    psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

}


// Convert a Unitary object to a Psi matrix.
Matrix UNITARYtoPSIMX(CheMPS2::DMRGSCFunitary&  unitary, const CheMPS2::DMRGSCFindices& iHandler){

    std::vector<int> dim(iHandler.getNirreps());

    for (auto h = 0; h < dim.size(); h++) {
        dim[h] = iHandler.getNORB(h);
    }
    Dimension nmopi(dim);
    Matrix U("Unitary", nmopi, nmopi);

    for (int irrep = 0; irrep < iHandler.getNirreps(); irrep++){
        for (int orb1 = 0; orb1 < iHandler.getNORB(irrep); orb1++){
            for (int orb2 = 0; orb2 < iHandler.getNORB(irrep); orb2++){
                U.set(irrep, orb1, orb2, unitary.getBlock(irrep)[ orb1 + iHandler.getNORB(irrep) * orb2 ] );
            }
        }
    }

    return U;

}


// Update the wavefunction coefficients
void update_WFNco(const Matrix& C_0, const CheMPS2::DMRGSCFindices& iHandler, CheMPS2::DMRGSCFunitary&  unitary, const Wavefunction& wfn){

    auto U = UNITARYtoPSIMX( unitary, iHandler);
    wfn.Ca()->gemm(false, true, 1.0, C_0, U, 0.0);
    wfn.Cb()->copy(wfn.Ca());

}

SharedWavefunction dmrg(SharedWavefunction ref_wfn, Options& options)
{
    auto dmrg = std::make_shared<DMRGSolver>(ref_wfn, options);
    dmrg->compute_energy();
    return dmrg;
}

DMRGSolver::DMRGSolver(SharedWavefunction ref_wfn, Options& options) : Wavefunction(options)
{
    reference_wavefunction_ = ref_wfn;
    shallow_copy(ref_wfn);
    Ca_ = ref_wfn->Ca()->clone();
    Cb_ = ref_wfn->Cb()->clone();
    Da_ = ref_wfn->Da()->clone();
    Db_ = ref_wfn->Db()->clone();
    Fa_ = ref_wfn->Fa()->clone();
    Fb_ = ref_wfn->Fb()->clone();
    same_a_b_orbs_ = true;
}

double DMRGSolver::compute_energy()
{

    /* This plugin is able to perform a DMRG calculation in a molecular orbital active space. */

    /*******************************
     *   Environment information   *
     *******************************/
    std::shared_ptr<PSIO> psio(_default_psio_lib_); // Grab the global (default) PSIO object, for file I/O
    if (!reference_wavefunction_){ throw PSIEXCEPTION("SCF has not been run yet!"); }

    /*************************
     *   Fetch the options   *
     *************************/

    const int wfn_irrep               = options_.get_int("DMRG_IRREP");
    const int wfn_multp               = options_.get_int("DMRG_MULTIPLICITY");
    auto dmrg_states                  = options_.get_int_vector("DMRG_SWEEP_STATES");
    const int ndmrg_states            = options_["DMRG_SWEEP_STATES"].size();
    double * dmrg_econv               = options_.get_double_array("DMRG_SWEEP_ENERGY_CONV");
    const int ndmrg_econv             = options_["DMRG_SWEEP_ENERGY_CONV"].size();
    auto dmrg_maxsweeps               = options_.get_int_vector("DMRG_SWEEP_MAX_SWEEPS");
    const int ndmrg_maxsweeps         = options_["DMRG_SWEEP_MAX_SWEEPS"].size();
    double * dmrg_noiseprefactors     = options_.get_double_array("DMRG_SWEEP_NOISE_PREFAC");
    const int ndmrg_noiseprefactors   = options_["DMRG_SWEEP_NOISE_PREFAC"].size();
    double * dmrg_dvdson_rtol         = options_.get_double_array("DMRG_SWEEP_DVDSON_RTOL");
    const int ndmrg_dvdson_rtol       = options_["DMRG_SWEEP_DVDSON_RTOL"].size();
    const bool dmrg_print_corr        = options_.get_bool("DMRG_PRINT_CORR");
    const bool mps_chkpt              = options_.get_bool("DMRG_MPS_WRITE");
    auto frozen_docc                  = options_.get_int_vector("RESTRICTED_DOCC");
    auto active                       = options_.get_int_vector("ACTIVE");
    const double d_convergence        = options_.get_double("DMRG_SCF_GRAD_THR");
    const bool dmrg_store_unit        = options_.get_bool("DMRG_UNITARY_WRITE");
    const bool dmrg_do_diis           = options_.get_bool("DMRG_DIIS");
    const double dmrg_diis_branch     = options_.get_double("DMRG_SCF_DIIS_THR");
    const bool dmrg_store_diis        = options_.get_bool("DMRG_DIIS_WRITE");
    const int dmrg_max_iter           = options_.get_int("DMRG_SCF_MAX_ITER");
    const int dmrg_which_root         = options_.get_int("DMRG_EXCITATION");
    const bool dmrg_state_avg         = options_.get_bool("DMRG_SCF_STATE_AVG");
    const std::string dmrg_active_space    = options_.get_str("DMRG_SCF_ACTIVE_SPACE");
    const bool dmrg_loc_random        = options_.get_bool("DMRG_LOCAL_INIT");
    const bool dmrg_caspt2            = options_.get_bool("DMRG_CASPT2_CALC");
    const std::string dmrg_caspt2_orb      = options_.get_str("DMRG_CASPT2_ORBS");
    const bool PSEUDOCANONICAL        = ( dmrg_caspt2_orb.compare("PSEUDOCANONICAL") == 0 ) ? true : false;
    const double dmrg_ipea            = options_.get_double("DMRG_CASPT2_IPEA");
    const double dmrg_imag_shift      = options_.get_double("DMRG_CASPT2_IMAG");
    const bool dmrg_molden            = options_.get_bool("DMRG_MOLDEN_WRITE");
    const bool dmrg_density_ao        = options_.get_bool("DMRG_OPDM_AO_PRINT");
    const int dmrg_num_vec_diis       = CheMPS2::DMRGSCF_numDIISvecs;
    const std::string unitaryname     = psi::get_writer_file_prefix( molecule()->name() ) + ".unitary.h5";
    const std::string diisname        = psi::get_writer_file_prefix( molecule()->name() ) + ".DIIS.h5";

    /****************************************
     *   Check if the input is consistent   *
     ****************************************/

    const int SyGroup= chemps2_groupnumber( molecule()->sym_label() );
    const auto& orbspi = nmopi();
    const auto& docc = doccpi();
    const auto& socc = soccpi();

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
    if ( options_["RESTRICTED_DOCC"].size() != nirrep_ ){ throw PSIEXCEPTION("Option RESTRICTED_DOCC (integer array) should contain as many elements as there are irreps!"); }
    if ( options_["ACTIVE"].size()      != nirrep_ ){ throw PSIEXCEPTION("Option ACTIVE (integer array) should contain as many elements as there are irreps!"); }
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
    auto OptScheme = std::make_unique<CheMPS2::ConvergenceScheme>(ndmrg_states);
    for (int cnt=0; cnt<ndmrg_states; cnt++){
       OptScheme->set_instruction( cnt, dmrg_states[cnt], dmrg_econv[cnt], dmrg_maxsweeps[cnt], dmrg_noiseprefactors[cnt], dmrg_dvdson_rtol[cnt] );
    }

    /******************************************************************************
     *   Print orbital information; check consistency of frozen_docc and active   *
     ******************************************************************************/

    std::vector<int> nvirtual(nirrep_);
    bool virtualsOK = true;
    for (int cnt=0; cnt<nirrep_; cnt++){
       nvirtual[cnt] = orbspi[cnt] - frozen_docc[cnt] - active[cnt];
       if ( nvirtual[cnt] < 0 ){ virtualsOK = false; }
    }
    (*outfile->stream()) << "wfn_irrep   = " << wfn_irrep << std::endl;
    (*outfile->stream()) << "wfn_multp   = " << wfn_multp << std::endl;
    (*outfile->stream()) << "numOrbitals = [ " << orbspi[0];
    for (int cnt=1; cnt<nirrep_; cnt++){ (*outfile->stream()) << " , " << orbspi[cnt];      } (*outfile->stream()) << " ]" << std::endl;
    (*outfile->stream()) << "R(O)HF DOCC = [ " << docc[0];
    for (int cnt=1; cnt<nirrep_; cnt++){ (*outfile->stream()) << " , " << docc[cnt];        } (*outfile->stream()) << " ]" << std::endl;
    (*outfile->stream()) << "R(O)HF SOCC = [ " << socc[0];
    for (int cnt=1; cnt<nirrep_; cnt++){ (*outfile->stream()) << " , " << socc[cnt];        } (*outfile->stream()) << " ]" << std::endl;
    (*outfile->stream()) << "frozen_docc = [ " << frozen_docc[0];
    for (int cnt=1; cnt<nirrep_; cnt++){ (*outfile->stream()) << " , " << frozen_docc[cnt]; } (*outfile->stream()) << " ]" << std::endl;
    (*outfile->stream()) << "active      = [ " << active[0];
    for (int cnt=1; cnt<nirrep_; cnt++){ (*outfile->stream()) << " , " << active[cnt];      } (*outfile->stream()) << " ]" << std::endl;
    (*outfile->stream()) << "virtual     = [ " << nvirtual[0];
    for (int cnt=1; cnt<nirrep_; cnt++){ (*outfile->stream()) << " , " << nvirtual[cnt];    } (*outfile->stream()) << " ]" << std::endl;
    if ( !virtualsOK ){ throw PSIEXCEPTION("For at least one irrep: frozen_docc[ irrep ] + active[ irrep ] > numOrbitals[ irrep ]!"); }

    /*******************************************
     *   Create another bit of DMRG preamble   *
     *******************************************/
    CheMPS2::DMRGSCFindices iHandler(nmo_, SyGroup, frozen_docc.data(), active.data(), nvirtual.data());
    CheMPS2::DMRGSCFunitary unitary(&iHandler);
    std::unique_ptr<CheMPS2::DIIS> theDIIS = nullptr;
    CheMPS2::DMRGSCFintegrals theRotatedTEI(&iHandler);
    const int nOrbDMRG = iHandler.getDMRGcumulative(nirrep_);
    std::vector<double> DMRG1DM(nOrbDMRG * nOrbDMRG);
    std::vector<double> DMRG2DM(nOrbDMRG * nOrbDMRG * nOrbDMRG * nOrbDMRG);
    CheMPS2::DMRGSCFmatrix theFmatrix(&iHandler);
    theFmatrix.clear();
    CheMPS2::DMRGSCFmatrix theQmatOCC(&iHandler);
    theQmatOCC.clear();
    CheMPS2::DMRGSCFmatrix theQmatACT(&iHandler);
    theQmatACT.clear();
    CheMPS2::DMRGSCFmatrix theTmatrix(&iHandler);
    theTmatrix.clear();
    CheMPS2::DMRGSCFwtilde wmattilde(&iHandler);

    /***************************************************
     *   Create the active space Hamiltonian storage   *
     ***************************************************/

    int nElectrons = 0;
    for (int cnt=0; cnt<nirrep_; cnt++){ nElectrons += 2 * docc[cnt] + socc[cnt]; }
    (*outfile->stream()) << "nElectrons  = " << nElectrons << std::endl;

    // Number of electrons in the active space
    int nDMRGelectrons = nElectrons;
    for (int cnt=0; cnt<nirrep_; cnt++){ nDMRGelectrons -= 2 * frozen_docc[cnt]; }
    (*outfile->stream()) << "nEl. active = " << nDMRGelectrons << std::endl;

    // Create the CheMPS2::Hamiltonian --> fill later
    std::vector<int> orbitalIrreps(nOrbDMRG);
    int counterFillOrbitalIrreps = 0;
    for (int h=0; h<nirrep_; h++){
       for (int cnt=0; cnt<active[h]; cnt++){ //Only the active space is treated with DMRG-SCF!
          orbitalIrreps[counterFillOrbitalIrreps] = h;
          counterFillOrbitalIrreps++;
       }
    }
    CheMPS2::Hamiltonian HamDMRG(nOrbDMRG, SyGroup, orbitalIrreps.data());

    /* Create the CheMPS2::Problem
       You can fill Ham later, as Problem only keeps a pointer to the Hamiltonian object.
       Since only doubly occupied frozen orbitals are allowed, wfn_multp and wfn_irrep do not change. */
    auto Prob = std::make_unique<CheMPS2::Problem>(&HamDMRG, wfn_multp-1, nDMRGelectrons, wfn_irrep);
    if ( !(Prob->checkConsistency()) ){ throw PSIEXCEPTION("CheMPS2::Problem : No Hilbert state vector compatible with all symmetry sectors!"); }
    Prob->SetupReorderD2h(); // Does nothing if group not d2h

    /**************************************
     *   Input is parsed and consistent   *
     *   Start with DMRG                  *
     **************************************/

    Matrix work1("work1", orbspi, orbspi);
    Matrix work2("work2", orbspi, orbspi);
    std::shared_ptr<JK> myJK = std::make_shared<DiskJK>(basisset(), options_);
    myJK->set_cutoff(0.0);
    myJK->initialize();
    const auto orig_coeff = *Ca()->clone();

    std::vector<int> OAorbs; // Occupied + active
    std::vector<int> Aorbs;  // Only active
    std::vector<int> Vorbs;  // Virtual
    std::vector<int> empty;
    for (int h = 0; h < iHandler.getNirreps(); h++){
       for (int orb = 0; orb < iHandler.getNOCC(h) + iHandler.getNDMRG(h); orb++){
          OAorbs.push_back( iHandler.getOrigNOCCstart(h) + orb );
       }
       for (int orb = 0; orb < iHandler.getNDMRG(h); orb++){
          Aorbs.push_back( iHandler.getOrigNDMRGstart(h) + orb );
       }
       for (int orb = 0; orb < iHandler.getNVIRT(h); orb++){
          Vorbs.push_back( iHandler.getOrigNVIRTstart(h) + orb );
       }
    }
    auto OAorbs_ptr = std::make_shared<MOSpace>('Q', OAorbs, empty);
    auto Aorbs_ptr = std::make_shared<MOSpace>('S', Aorbs, empty);
    auto Vorbs_ptr = std::make_shared<MOSpace>('T', Vorbs, empty);
    std::vector<std::shared_ptr<MOSpace> > spaces;
    spaces.push_back( OAorbs_ptr   );
    spaces.push_back(  Aorbs_ptr   );
    spaces.push_back(  Vorbs_ptr   );
    spaces.push_back( MOSpace::all );
    // CheMPS2 requires RHF or ROHF orbitals.
    IntegralTransform ints(IntegralTransform(shared_from_this(), spaces, IntegralTransform::TransformationType::Restricted));
    ints.set_keep_iwl_so_ints( true );
    ints.set_keep_dpd_so_ints( true );

    (*outfile->stream()) << "###########################################################" << std::endl;
    (*outfile->stream()) << "###                                                     ###" << std::endl;
    (*outfile->stream()) << "###                       DMRG-SCF                      ###" << std::endl;
    (*outfile->stream()) << "###                                                     ###" << std::endl;
    (*outfile->stream()) << "###            CheMPS2 by Sebastian Wouters             ###" << std::endl;
    (*outfile->stream()) << "###        https://github.com/SebWouters/CheMPS2        ###" << std::endl;
    (*outfile->stream()) << "###   Comput. Phys. Commun. 185 (6), 1501-1514 (2014)   ###" << std::endl;
    (*outfile->stream()) << "###                                                     ###" << std::endl;
    (*outfile->stream()) << "###########################################################" << std::endl;
    (*outfile->stream()) << std::endl;
    (*outfile->stream()) << "Number of variables in the x-matrix = " << unitary.getNumVariablesX() << std::endl;

    //Convergence variables
    double gradNorm = 1.0;
    double updateNorm = 1.0;
    std::vector<double> theupdate(unitary.getNumVariablesX(), 0.0);
    std::vector<double> theDIISparameterVector;
    double Energy = 1e8;

    int theDIISvectorParamSize = 0;
    int maxlinsize = 0;
    for (int irrep=0; irrep<nirrep_; irrep++){
        const int linsize_irrep = iHandler.getNORB(irrep);
        theDIISvectorParamSize += linsize_irrep*(linsize_irrep-1)/2;
        if (linsize_irrep>maxlinsize){ maxlinsize = linsize_irrep; }
    }

    const int nOrbDMRG_pow4    = nOrbDMRG * nOrbDMRG * nOrbDMRG * nOrbDMRG;
    const int unitary_worksize = 4 * maxlinsize * maxlinsize;
    const int sizeWorkMem      = ( nOrbDMRG_pow4 > unitary_worksize ) ? nOrbDMRG_pow4 : unitary_worksize;
    const int tot_dmrg_power6  = nOrbDMRG_pow4 * nOrbDMRG * nOrbDMRG;
    std::vector<double> mem1(sizeWorkMem);
    std::vector<double> mem2(( PSEUDOCANONICAL ) ? sizeWorkMem : std::max( sizeWorkMem, tot_dmrg_power6 ));

    std::unique_ptr<CheMPS2::EdmistonRuedenberg> theLocalizer = nullptr;
    if ( dmrg_active_space.compare("LOC")==0 ){ theLocalizer = std::make_unique<CheMPS2::EdmistonRuedenberg>( HamDMRG.getVmat(), iHandler.getGroupNumber() ); }

    //Load unitary from disk
    if ( dmrg_store_unit ){
        struct stat stFileInfo;
        int intStat = stat( unitaryname.c_str(), &stFileInfo );
        if (intStat==0){ unitary.loadU( unitaryname ); }
    }

    //Load DIIS from disk
    if (( dmrg_do_diis ) && ( dmrg_store_diis )){
        struct stat stFileInfo;
        int intStat = stat( diisname.c_str(), &stFileInfo );
        if (intStat==0){
            if (theDIIS == nullptr){
                theDIIS = std::make_unique<CheMPS2::DIIS>( theDIISvectorParamSize, unitary.getNumVariablesX(), dmrg_num_vec_diis );
                theDIISparameterVector = std::vector<double>(theDIISvectorParamSize);
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
        if (unitary.getNumVariablesX() > 0){

            std::ofstream capturing;
            std::streambuf * cout_buffer;
            std::string chemps2filename = outfile_name + ".chemps2";
            (*outfile->stream()) << "CheMPS2 output is temporarily written to the file " << chemps2filename << " and will be copied here." << std::endl;
            capturing.open( chemps2filename.c_str() , std::ios::trunc ); // truncate
            cout_buffer = std::cout.rdbuf( capturing.rdbuf() );

            unitary.updateUnitary(mem1.data(), mem2.data(), theupdate.data(), true, true); //multiply = compact = true
            if (( dmrg_do_diis ) && ( updateNorm <= dmrg_diis_branch )){
                if ( dmrg_active_space.compare("NO")==0 ){
                    std::cout << "DIIS has started. Active space not rotated to NOs anymore!" << std::endl;
                }
                if ( dmrg_active_space.compare("LOC")==0 ){
                    std::cout << "DIIS has started. Active space not rotated to localized orbitals anymore!" << std::endl;
                }
                if (theDIIS == nullptr){
                    theDIIS = std::make_unique<CheMPS2::DIIS>( theDIISvectorParamSize, unitary.getNumVariablesX(), dmrg_num_vec_diis );
                    theDIISparameterVector = std::vector<double>(theDIISvectorParamSize);
                    unitary.makeSureAllBlocksDetOne(mem1.data(), mem2.data());
                }
                unitary.getLog(theDIISparameterVector.data(), mem1.data(), mem2.data());
                theDIIS->appendNew(theupdate.data(), theDIISparameterVector.data());
                theDIIS->calculateParam(theDIISparameterVector.data());
                unitary.updateUnitary(mem1.data(), mem2.data(), theDIISparameterVector.data(), false, false); //multiply = compact = false
            }

            std::cout.rdbuf(cout_buffer);
            capturing.close();
            std::ifstream copying;
            copying.open( chemps2filename , std::ios::in ); // read only
            if (copying.is_open()){
                std::string line;
                while( getline( copying, line ) ){ (*outfile->stream()) << line << std::endl; }
                copying.close();
            }
            system(("rm " + chemps2filename).c_str());

        }
        if (( dmrg_store_unit ) && (gradNorm!=1.0)){ unitary.saveU( unitaryname ); }
        if (( dmrg_store_diis ) && (updateNorm!=1.0) && (theDIIS!=nullptr)){ theDIIS->saveDIIS( diisname ); }

        //Fill HamDMRG
        update_WFNco(orig_coeff, iHandler, unitary, *this);
        buildTmatrix(theTmatrix, iHandler, *Ca_, *this);
        buildQmatOCC(theQmatOCC, iHandler, work1, work2, *Ca_, *myJK);
        buildHamDMRG(ints, Aorbs_ptr, theTmatrix, theQmatOCC, iHandler, HamDMRG, psio, *this);

        //Localize the active space and reorder the orbitals within each irrep based on the exchange matrix
        if (( dmrg_active_space.compare("LOC")==0 ) && (theDIIS==nullptr)){ //When the DIIS has started: stop

            std::ofstream capturing;
            std::streambuf * cout_buffer;
            std::string chemps2filename = outfile_name + ".chemps2";
            (*outfile->stream()) << "CheMPS2 output is temporarily written to the file " << chemps2filename << " and will be copied here." << std::endl;
            capturing.open( chemps2filename.c_str() , std::ios::trunc ); // truncate
            cout_buffer = std::cout.rdbuf( capturing.rdbuf() );

            theLocalizer->Optimize( mem1.data(), mem2.data(), dmrg_loc_random );
            theLocalizer->FiedlerExchange(maxlinsize, mem1.data(), mem2.data());
            CheMPS2::CASSCF::fillLocalizedOrbitalRotations(theLocalizer->getUnitary(), &iHandler, mem1.data());
            unitary.rotateActiveSpaceVectors(mem1.data(), mem2.data());

            std::cout.rdbuf(cout_buffer);
            capturing.close();
            std::ifstream copying;
            copying.open( chemps2filename , std::ios::in ); // read only
            if (copying.is_open()){
                std::string line;
                while( getline( copying, line ) ){ (*outfile->stream()) << line << std::endl; }
                copying.close();
            }
            system(("rm " + chemps2filename).c_str());

            update_WFNco(orig_coeff, iHandler, unitary, *this);
            buildTmatrix(theTmatrix, iHandler, *Ca_, *this);
            buildQmatOCC(theQmatOCC, iHandler, work1, work2, *Ca_, *myJK);
            buildHamDMRG(ints, Aorbs_ptr, theTmatrix, theQmatOCC, iHandler, HamDMRG, psio, *this );
            (*outfile->stream()) << "Rotated the active space to localized orbitals, sorted according to the exchange matrix." << std::endl;

        }

        //Do the DMRG sweeps, and calculate the 2DM
        {
            std::ofstream capturing;
            std::streambuf * cout_buffer;
            std::string chemps2filename = outfile_name + ".chemps2";
            (*outfile->stream()) << "CheMPS2 output is temporarily written to the file " << chemps2filename << " and will be copied here." << std::endl;
            capturing.open( chemps2filename.c_str() , std::ios::trunc ); // truncate
            cout_buffer = std::cout.rdbuf( capturing.rdbuf() );

            for (int cnt = 0; cnt < nOrbDMRG_pow4; cnt++){ DMRG2DM[ cnt ] = 0.0; } //Clear the 2-RDM (to allow for state-averaged calculations)
            const std::string psi4TMPpath = PSIOManager::shared_object()->get_default_path();
            auto theDMRG = std::make_unique<CheMPS2::DMRG>(Prob.get(), OptScheme.get(), mps_chkpt, psi4TMPpath);
            for (int state = -1; state < dmrg_which_root; state++){
                if (state > -1){ theDMRG->newExcitation( std::fabs( Energy ) ); }
                Energy = theDMRG->Solve();
                if ( dmrg_state_avg ){ // When SA-DMRGSCF: 2DM += current 2DM
                    theDMRG->calc2DMandCorrelations();
                    CheMPS2::CASSCF::copy2DMover( theDMRG->get2DM(), nOrbDMRG, DMRG2DM.data() );
                }
                if ((state == -1) && (dmrg_which_root > 0)){ theDMRG->activateExcitations( dmrg_which_root ); }
            }
            if ( !(dmrg_state_avg) ){ // When SS-DMRGSCF: 2DM += last 2DM
                theDMRG->calc2DMandCorrelations();
                CheMPS2::CASSCF::copy2DMover( theDMRG->get2DM(), nOrbDMRG, DMRG2DM.data() );
            }
            if ( dmrg_print_corr ){ theDMRG->getCorrelations()->Print(); }
            if ( CheMPS2::DMRG_storeRenormOptrOnDisk ){ theDMRG->deleteStoredOperators(); }
            theDMRG.reset();
            if ((dmrg_state_avg) && (dmrg_which_root > 0)){
                const double averagingfactor = 1.0 / (dmrg_which_root+1);
                for (int cnt = 0; cnt < nOrbDMRG_pow4; cnt++){ DMRG2DM[ cnt ] *= averagingfactor; }
            }
            CheMPS2::CASSCF::setDMRG1DM( nDMRGelectrons, nOrbDMRG, DMRG1DM.data(), DMRG2DM.data() );

            std::cout.rdbuf(cout_buffer);
            capturing.close();
            std::ifstream copying;
            copying.open( chemps2filename , std::ios::in ); // read only
            if (copying.is_open()){
                std::string line;
                while( getline( copying, line ) ){ (*outfile->stream()) << line << std::endl; }
                copying.close();
            }
            system(("rm " + chemps2filename).c_str());
        }

        if (( dmrg_active_space.compare("NO")==0 ) && (theDIIS==nullptr)){ //When the DIIS has started: stop
            CheMPS2::CASSCF::copy_active( DMRG1DM.data(), &theFmatrix, &iHandler, true );
            CheMPS2::CASSCF::block_diagonalize( 'A', &theFmatrix, &unitary, mem1.data(), mem2.data(), &iHandler, true, DMRG2DM.data(), nullptr, nullptr ); // Unitary is updated and DMRG2DM rotated
            CheMPS2::CASSCF::setDMRG1DM( nDMRGelectrons, nOrbDMRG, DMRG1DM.data(), DMRG2DM.data() );
            update_WFNco(orig_coeff, iHandler, unitary, *this);
            buildTmatrix(theTmatrix, iHandler, *Ca_, *this);
            buildQmatOCC(theQmatOCC, iHandler, work1, work2, *Ca_, *myJK);
            (*outfile->stream()) << "Rotated the active space to natural orbitals, sorted according to the NOON." << std::endl;
        }

        if (dmrg_max_iter == nIterations){
            if ( dmrg_store_unit ){ unitary.saveU( unitaryname ); }
            break;
        }

        buildQmatACT( theQmatACT, iHandler, DMRG1DM, work1, work2, *Ca_, *myJK);
        fillRotatedTEI_coulomb(  ints, OAorbs_ptr, theRotatedTEI, iHandler, psio);
        fillRotatedTEI_exchange( ints, OAorbs_ptr, Vorbs_ptr,  theRotatedTEI, iHandler, psio );

        {
            std::ofstream capturing;
            std::streambuf * cout_buffer;
            std::string chemps2filename = outfile_name + ".chemps2";
            (*outfile->stream()) << "CheMPS2 output is temporarily written to the file " << chemps2filename << " and will be copied here." << std::endl;
            capturing.open( chemps2filename.c_str() , std::ios::trunc ); // truncate
            cout_buffer = std::cout.rdbuf( capturing.rdbuf() );

            CheMPS2::CASSCF::buildFmat(&theFmatrix, &theTmatrix, &theQmatOCC, &theQmatACT, &iHandler, &theRotatedTEI, DMRG2DM.data(), DMRG1DM.data());
            CheMPS2::CASSCF::buildWtilde(&wmattilde, &theTmatrix, &theQmatOCC, &theQmatACT, &iHandler, &theRotatedTEI, DMRG2DM.data(), DMRG1DM.data());
            CheMPS2::CASSCF::augmentedHessianNR(&theFmatrix, &wmattilde, &iHandler, &unitary, theupdate.data(), &updateNorm, &gradNorm);

            std::cout.rdbuf(cout_buffer);
            capturing.close();
            std::ifstream copying;
            copying.open( chemps2filename , std::ios::in ); // read only
            if (copying.is_open()){
                std::string line;
                while( getline( copying, line ) ){ (*outfile->stream()) << line << std::endl; }
                copying.close();
            }
            system(("rm " + chemps2filename).c_str());
        }
    }

    outfile->Printf("The DMRG-SCF energy = %3.10f \n", Energy);
    set_energy(Energy);
    set_scalar_variable("CURRENT ENERGY", Energy);
    set_scalar_variable("DMRG-SCF TOTAL ENERGY", Energy);
    set_module("chemps2");

    if ((( dmrg_molden ) || (( dmrg_caspt2 ) && ( PSEUDOCANONICAL ))) && ( nIterations > 0 )){

        (*outfile->stream()) << "################################################" << std::endl;
        (*outfile->stream()) << "###                                          ###" << std::endl;
        (*outfile->stream()) << "###   Rotation to pseudocanonical orbitals   ###" << std::endl;
        (*outfile->stream()) << "###                                          ###" << std::endl;
        (*outfile->stream()) << "################################################" << std::endl;
        CheMPS2::CASSCF::construct_fock(&theFmatrix, &theTmatrix, &theQmatOCC, &theQmatACT, &iHandler );
        CheMPS2::CASSCF::block_diagonalize( 'O', &theFmatrix, &unitary, mem1.data(), mem2.data(), &iHandler, false, nullptr, nullptr, nullptr );
        CheMPS2::CASSCF::block_diagonalize( 'A', &theFmatrix, &unitary, mem1.data(), mem2.data(), &iHandler, false, DMRG2DM.data(), nullptr, nullptr );
        CheMPS2::CASSCF::block_diagonalize( 'V', &theFmatrix, &unitary, mem1.data(), mem2.data(), &iHandler, false, nullptr, nullptr, nullptr );
        CheMPS2::CASSCF::setDMRG1DM( nDMRGelectrons, nOrbDMRG, DMRG1DM.data(), DMRG2DM.data());
        update_WFNco(orig_coeff, iHandler, unitary, *this);
        buildTmatrix(theTmatrix, iHandler, *Ca(), *this);
        buildQmatOCC(theQmatOCC, iHandler, work1, work2, *Ca(), *myJK);
        buildQmatACT(theQmatACT, iHandler, DMRG1DM, work1, work2, *Ca(), *myJK);
        CheMPS2::CASSCF::construct_fock(&theFmatrix, &theTmatrix, &theQmatOCC, &theQmatACT, &iHandler );

        epsilon_a_->set_name("DMRG Pseudocanonical Orbital Energies");
        for ( int h = 0; h < nirrep_; h++ ){
            const int NOCC = iHandler.getNOCC( h );
            const int NACT = iHandler.getNDMRG( h );
            const int NVIR = iHandler.getNVIRT( h );
            for ( int orb = 0; orb < orbspi[ h ]; orb++ ){
               epsilon_a_->set( h, orb, theFmatrix.get(h, orb, orb ) );
            }
        }
        epsilon_b_ = epsilon_a_;
    }

    if (( dmrg_molden ) && ( nIterations > 0 )){
        // Prevent writing beta orbitals.
        // This may need to be revisited, if we ever have gradients.
        same_a_b_dens_ = true;
    }

    // Reference on what the density is:
    // https://github.com/SebWouters/CheMPS2/issues/83
    density_map_["DMRG MS-AVERAGED (MO)"] = std::make_shared<Matrix>(std::move(build_rdm(iHandler, DMRG1DM)));
    density_map_["DMRG MS-AVERAGED"] = std::make_shared<Matrix>(std::move(mo_to_ao(*density_map_["DMRG MS-AVERAGED (MO)"], *Ca_, *aotoso())));

    if ( dmrg_density_ao ){

        (*outfile->stream()) << "############################" << std::endl;
        (*outfile->stream()) << "###                      ###" << std::endl;
        (*outfile->stream()) << "###   DMRG 1-RDM in AO   ###" << std::endl;
        (*outfile->stream()) << "###                      ###" << std::endl;
        (*outfile->stream()) << "############################" << std::endl;
        (*outfile->stream()) << "Please check the molden file for AO basis function information." << std::endl;
        density_map_["DMRG MS-AVERAGED"]->print("outfile");
    }

    if (( dmrg_caspt2 ) && ( nIterations > 0 )){

        (*outfile->stream()) << "###########################" << std::endl;
        (*outfile->stream()) << "###                     ###" << std::endl;
        (*outfile->stream()) << "###     DMRG-CASPT2     ###" << std::endl;
        (*outfile->stream()) << "###                     ###" << std::endl;
        (*outfile->stream()) << "###########################" << std::endl;

        buildHamDMRG( ints, Aorbs_ptr, theTmatrix, theQmatOCC, iHandler, HamDMRG, psio, *this );

        std::vector<double> contract(tot_dmrg_power6);
        std::vector<double> three_dm(tot_dmrg_power6);

        {
            std::ofstream capturing;
            std::streambuf * cout_buffer;
            std::string chemps2filename = outfile_name + ".chemps2";
            (*outfile->stream()) << "CheMPS2 output is temporarily written to the file " << chemps2filename << " and will be copied here." << std::endl;
            capturing.open( chemps2filename.c_str() , std::ios::trunc ); // truncate
            cout_buffer = std::cout.rdbuf( capturing.rdbuf() );

            for (int cnt = 0; cnt < nOrbDMRG_pow4; cnt++){ DMRG2DM[ cnt ] = 0.0; } //Clear the 2-RDM (to allow for state-averaged calculations)
            const std::string psi4TMPpath = PSIOManager::shared_object()->get_default_path();
            auto theDMRG = std::make_unique<CheMPS2::DMRG>(Prob.get(), OptScheme.get(), false, psi4TMPpath); // Rotated orbital space --> do not use checkpoint
            for (int state = -1; state < dmrg_which_root; state++){
                if (state > -1){ theDMRG->newExcitation( std::fabs( Energy ) ); }
                const double E_CASSCF = theDMRG->Solve();
                if ((state == -1) && (dmrg_which_root > 0)){ theDMRG->activateExcitations( dmrg_which_root ); }
            }
            theDMRG->calc_rdms_and_correlations( true );
            CheMPS2::CASSCF::copy2DMover( theDMRG->get2DM(), nOrbDMRG, DMRG2DM.data());
            CheMPS2::CASSCF::setDMRG1DM( nDMRGelectrons, nOrbDMRG, DMRG1DM.data(), DMRG2DM.data());
            buildQmatACT( theQmatACT, iHandler, DMRG1DM, work1, work2, *Ca_, *myJK);
            CheMPS2::CASSCF::construct_fock(&theFmatrix, &theTmatrix, &theQmatOCC, &theQmatACT, &iHandler );
            CheMPS2::CASSCF::copy_active(&theFmatrix, mem2.data(), &iHandler );                   // Fock
            for ( int cnt = 0; cnt < tot_dmrg_power6; cnt++ ){ contract[ cnt ] = 0.0; }
            for ( int ham_orbz = 0; ham_orbz < nOrbDMRG; ham_orbz++ ){
               theDMRG->Symm4RDM( three_dm.data(), ham_orbz, ham_orbz, false );
               int size = tot_dmrg_power6;
               double f_zz = 0.5 * mem2[ ham_orbz + nOrbDMRG * ham_orbz ];
               int inc1 = 1;
               daxpy_( &size, &f_zz, three_dm.data(), &inc1, contract.data(), &inc1 ); // trace( Fock * 4-RDM )
            }
            if ( PSEUDOCANONICAL == false ){
               for ( int ham_orb1 = 0; ham_orb1 < nOrbDMRG; ham_orb1++ ){
                  for ( int ham_orb2 = ham_orb1 + 1; ham_orb2 < nOrbDMRG; ham_orb2++ ){
                     if ( HamDMRG.getOrbitalIrrep( ham_orb1 ) == HamDMRG.getOrbitalIrrep( ham_orb2 ) ){
                        theDMRG->Symm4RDM( three_dm.data(), ham_orb1, ham_orb2, false );
                        int size = tot_dmrg_power6;
                        double f_12 = 0.5 * ( mem2[ ham_orb1 + nOrbDMRG * ham_orb2 ] + mem2[ ham_orb2 + nOrbDMRG * ham_orb1 ] );
                        int inc1 = 1;
                        daxpy_( &size, &f_12, three_dm.data(), &inc1, contract.data(), &inc1 ); // trace( Fock * 4-RDM )
                     }
                  }
               }
               // CheMPS2::Cumulant::gamma4_fock_contract_ham( Prob, theDMRG->get3DM(), theDMRG->get2DM(), mem2, contract );
            }
            theDMRG->get3DM()->fill_ham_index( 1.0, false, three_dm.data(), 0, nOrbDMRG ); // 3-RDM --> three_dm was used as work space for the constracted 4-RDM
            if (CheMPS2::DMRG_storeMpsOnDisk){        theDMRG->deleteStoredMPS();       }
            if (CheMPS2::DMRG_storeRenormOptrOnDisk){ theDMRG->deleteStoredOperators(); }
            theDMRG.reset();

            std::cout.rdbuf(cout_buffer);
            capturing.close();
            std::ifstream copying;
            copying.open( chemps2filename , std::ios::in ); // read only
            if (copying.is_open()){
                std::string line;
                while( getline( copying, line ) ){ (*outfile->stream()) << line << std::endl; }
                copying.close();
            }
            system(("rm " + chemps2filename).c_str());
       }

       if ( PSEUDOCANONICAL == false ){
           (*outfile->stream()) << "CASPT2 : Deviation from pseudocanonical = " << CheMPS2::CASSCF::deviation_from_blockdiag(&theFmatrix, &iHandler ) << std::endl;
           CheMPS2::CASSCF::block_diagonalize( 'O', &theFmatrix, &unitary, mem1.data(), mem2.data(), &iHandler, false, nullptr, nullptr, nullptr );
           CheMPS2::CASSCF::block_diagonalize( 'A', &theFmatrix, &unitary, mem1.data(), mem2.data(), &iHandler, false, DMRG2DM.data(), three_dm.data(), contract.data() ); // 2-RDM, 3-RDM, and trace( Fock * cu(4)-4-RDM )
           CheMPS2::CASSCF::block_diagonalize( 'V', &theFmatrix, &unitary, mem1.data(), mem2.data(), &iHandler, false, nullptr, nullptr, nullptr );
           CheMPS2::CASSCF::setDMRG1DM( nDMRGelectrons, nOrbDMRG, DMRG1DM.data(), DMRG2DM.data()); // 1-RDM
           update_WFNco(orig_coeff, iHandler, unitary, *this);
           buildTmatrix(theTmatrix, iHandler, *Ca_, *this);
           buildQmatOCC(theQmatOCC, iHandler, work1, work2, *Ca_, *myJK);
           buildQmatACT(theQmatACT, iHandler, DMRG1DM, work1, work2, *Ca_, *myJK);
           CheMPS2::CASSCF::construct_fock(&theFmatrix, &theTmatrix, &theQmatOCC, &theQmatACT, &iHandler ); // Fock
       }

       fillRotatedTEI_coulomb(  ints, OAorbs_ptr, theRotatedTEI, iHandler, psio);
       fillRotatedTEI_exchange( ints, OAorbs_ptr, Vorbs_ptr,  theRotatedTEI, iHandler, psio );

       (*outfile->stream()) << "CASPT2 : Norm F - F_pseudocan = " << CheMPS2::CASSCF::deviation_from_blockdiag(&theFmatrix, &iHandler ) << std::endl;
       double E_CASPT2 = 0.0;
       {
            std::ofstream capturing;
            std::streambuf * cout_buffer;
            std::string chemps2filename = outfile_name + ".chemps2";
            (*outfile->stream()) << "CheMPS2 output is temporarily written to the file " << chemps2filename << " and will be copied here." << std::endl;
            capturing.open( chemps2filename.c_str() , std::ios::trunc ); // truncate
            cout_buffer = std::cout.rdbuf( capturing.rdbuf() );

            auto myCASPT2 = std::make_unique<CheMPS2::CASPT2>(&iHandler, &theRotatedTEI, &theTmatrix, &theFmatrix, DMRG1DM.data(), DMRG2DM.data(), three_dm.data(), contract.data(), dmrg_ipea );
            E_CASPT2 = myCASPT2->solve( dmrg_imag_shift );
            myCASPT2.reset();

            std::cout.rdbuf(cout_buffer);
            capturing.close();
            std::ifstream copying;
            copying.open( chemps2filename , std::ios::in ); // read only
            if (copying.is_open()){
                std::string line;
                while( getline( copying, line ) ){ (*outfile->stream()) << line << std::endl; }
                copying.close();
            }
            system(("rm " + chemps2filename).c_str());

       }

       outfile->Printf("The DMRG-CASPT2 variational correction energy = %3.10f \n", E_CASPT2);
       outfile->Printf("The DMRG-CASPT2 energy = %3.10f \n", Energy + E_CASPT2);
       set_energy(Energy + E_CASPT2);
       set_scalar_variable("CURRENT ENERGY", Energy + E_CASPT2);
       set_scalar_variable("DMRG-CASPT2 TOTAL ENERGY", Energy + E_CASPT2);
       set_module("chemps2");

    }

    return this->scalar_variable("CURRENT ENERGY");
}

SharedVector DMRGSolver::occupation_a() const {
    if (density_map_.find("DMRG MS-AVERAGED (MO)") == density_map_.end()) {
        throw PSIEXCEPTION("Density matrix not yet set.");
    }
    const auto& density = density_map_.at("DMRG MS-AVERAGED (MO)");
    auto occupation = std::make_shared<Vector>(density->rowspi());
    for (size_t h = 0; h < occupation->nirrep(); h++) {
        for (size_t i = 0; i < occupation->dim(h); i++) {
            occupation->set(h, i, 0.5 * density->get(h, i, i));
        }
    }
    return occupation;
}

}} // End Namespaces
#endif
// clang-format on
