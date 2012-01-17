/** Standard library includes */
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sstream>
#include <fstream>
#include <string> 
#include <iomanip>
#include <vector> 

#include <libtrans/integraltransform.h>
#include <libtrans/mospace.h>
#include <libmints/mints.h>
#include <liboptions/liboptions.h>
#include <libdiis/diismanager.h>
#include <libciomr/libciomr.h>

#include "omp2wave.h"
#include "defines.h"

using namespace psi;
using namespace boost;

namespace psi { namespace omp2wave {

OMP2Wave::OMP2Wave(boost::shared_ptr<Wavefunction> reference_wavefunction, Options &options)
    : Wavefunction(options, _default_psio_lib_)
{
    reference_wavefunction_ = reference_wavefunction;
    common_init();
}//

OMP2Wave::~OMP2Wave()
{
}//


void OMP2Wave::common_init()
{
 
	tol_Eod=options_.get_double("E_CONVERGENCE");
	tol_t2=options_.get_double("R_CONVERGENCE");
    tol_grad=options_.get_double("RMS_MOGRAD_CONVERGENCE");
	mograd_max=options_.get_double("MAX_MOGRAD_CONVERGENCE");
        cc_maxiter=options_.get_int("CC_MAXITER");
	mo_maxiter=options_.get_int("MO_MAXITER");
	print_=options_.get_int("PRINT"); 
	cachelev=options_.get_int("CACHELEVEL"); 
	num_vecs=options_.get_int("DIIS_MAX_VECS");
	exp_cutoff=options_.get_int("CUTOFF");
	memory=options_.get_int("MEMORY"); 
	
	step_max=options_.get_double("MO_STEP_MAX");
	lshift_parameter=options_.get_double("LEVEL_SHIFT");
	os_scale=options_.get_double("MP2_OS_SCALE");
	ss_scale=options_.get_double("MP2_SS_SCALE");
	sos_scale=options_.get_double("SOS_SCALE");
	sos_scale2=options_.get_double("SOS_SCALE2");
	
	//wfn=options_.get_str("WFN");
    if (options_["LEVEL_SHIFT"].has_changed()) {
	  level_shift="TRUE";
    }
	//lineq=options_.get_str("LINEQ");
	orth_type=options_.get_str("ORTH_TYPE");
	//stability=options_.get_str("STABILITY");
	opt_method=options_.get_str("OPT_METHOD");
	hess_type=options_.get_str("HESS_TYPE");
	omp2_orb_energy=options_.get_str("OMP2_ORBS_PRINT");
	natorb=options_.get_str("NAT_ORBS");
	reference=options_.get_str("REFERENCE");
	do_scs=options_.get_str("DO_SCS");
	do_sos=options_.get_str("DO_SOS");
	write_mo_coeff=options_.get_str("MO_WRITE");
	read_mo_coeff=options_.get_str("MO_READ");
	
	cutoff = pow(10.0,-exp_cutoff);
	
	
	if (print_ > 0) options_.print();
        title();
	get_moinfo();
	
	// Memory allocation
	HmoA = boost::shared_ptr<Matrix>(new Matrix("MO-basis alpha one-electron ints", nirreps, mopi, mopi));
	HmoB = boost::shared_ptr<Matrix>(new Matrix("MO-basis beta one-electron ints", nirreps, mopi, mopi));
	FockA = boost::shared_ptr<Matrix>(new Matrix("MO-basis alpha Fock matrix", nirreps, mopi, mopi));
	FockB = boost::shared_ptr<Matrix>(new Matrix("MO-basis beta Fock matrix", nirreps, mopi, mopi));
	gamma1corrA = boost::shared_ptr<Matrix>(new Matrix("MO-basis alpha correlation OPDM", nirreps, mopi, mopi));
	gamma1corrB = boost::shared_ptr<Matrix>(new Matrix("MO-basis beta correlation OPDM", nirreps, mopi, mopi));
	g1symmA = boost::shared_ptr<Matrix>(new Matrix("MO-basis alpha OPDM", nirreps, mopi, mopi));
	g1symmB = boost::shared_ptr<Matrix>(new Matrix("MO-basis beta OPDM", nirreps, mopi, mopi));
	GFockA = boost::shared_ptr<Matrix>(new Matrix("MO-basis alpha generalized Fock matrix", nirreps, mopi, mopi));
	GFockB = boost::shared_ptr<Matrix>(new Matrix("MO-basis beta generalized Fock matrix", nirreps, mopi, mopi));
	UorbA = boost::shared_ptr<Matrix>(new Matrix("Alpha MO rotation matrix", nirreps, mopi, mopi));
	UorbB = boost::shared_ptr<Matrix>(new Matrix("Beta MO rotation matrix", nirreps, mopi, mopi));
	KorbA = boost::shared_ptr<Matrix>(new Matrix("K alpha MO rotation", nirreps, mopi, mopi)); 
	KorbB = boost::shared_ptr<Matrix>(new Matrix("K beta MO rotation", nirreps, mopi, mopi)); 
	KsqrA = boost::shared_ptr<Matrix>(new Matrix("K^2 alpha MO rotation", nirreps, mopi, mopi)); 
	KsqrB = boost::shared_ptr<Matrix>(new Matrix("K^2 beta MO rotation", nirreps, mopi, mopi)); 
	HG1A = boost::shared_ptr<Matrix>(new Matrix("Alpha h*g1symm", nirreps, mopi, mopi));
	HG1B = boost::shared_ptr<Matrix>(new Matrix("Beta h*g1symm", nirreps, mopi, mopi));
	WorbA = boost::shared_ptr<Matrix>(new Matrix("Alpha MO gradient matrix", nirreps, mopi, mopi));
	WorbB = boost::shared_ptr<Matrix>(new Matrix("Beta MO gradient matrix", nirreps, mopi, mopi));
	GooA = boost::shared_ptr<Matrix>(new Matrix("Alpha Goo intermediate", nirreps, aoccpiA, aoccpiA));
	GooB = boost::shared_ptr<Matrix>(new Matrix("Beta Goo intermediate", nirreps, aoccpiB, aoccpiB));
	GvvA = boost::shared_ptr<Matrix>(new Matrix("Alpha Gvv intermediate", nirreps, avirtpiA, avirtpiA));
	GvvB = boost::shared_ptr<Matrix>(new Matrix("Beta Gvv intermediate", nirreps, avirtpiB, avirtpiB));

        Molecule& mol = *reference_wavefunction_->molecule().get();
        CharacterTable ct = mol.point_group()->char_table();
        fprintf(outfile,"\tMO spaces per irreps... \n\n"); fflush(outfile);
        fprintf(outfile, "\tIRREP   FC   AOCC  BOCC  AVIR    BVIR  FV \n");
        fprintf(outfile, "\t==========================================\n");                                                 
        for(int h = 0; h < nirrep_; ++h){
         fprintf(outfile, "\t %3s   %3d   %3d   %3d   %3d    %3d   %3d\n",
                             ct.gamma(h).symbol(), frzcpi[h], aoccpiA[h], aoccpiB[h], avirtpiA[h], avirtpiB[h], frzvpi[h]);
        }
        fprintf(outfile,     "\t==========================================\n\n");
	fflush(outfile);

	
    // Alloc ints
    std::vector<boost::shared_ptr<MOSpace> > spaces;
    spaces.push_back(MOSpace::occ);
    spaces.push_back(MOSpace::vir);

    // Perform MO integral transformation on spaces
//    ints = new IntegralTransform(chkpt_, spaces, 
//                           IntegralTransform::Unrestricted,
//                           IntegralTransform::DPDOnly,
//                           IntegralTransform::QTOrder,
//                           IntegralTransform::None,
//                           false);
    ints = new IntegralTransform(reference_wavefunction_, spaces, 
                           IntegralTransform::Unrestricted,
                           IntegralTransform::DPDOnly,
                           IntegralTransform::QTOrder,
                           IntegralTransform::None,
                           false);
                           
                          
    ints->set_print(0);
    ints->set_dpd_id(0);
    ints->set_keep_iwl_so_ints(true);
    ints->set_keep_dpd_so_ints(true);           
    ints->initialize();
    dpd_set_default(ints->get_dpd_id());

}//

void OMP2Wave::title()
{
   fprintf(outfile,"\n");
   fprintf(outfile," ============================================================================== \n");
   fprintf(outfile," ============================================================================== \n");
   fprintf(outfile," ============================================================================== \n");
   fprintf(outfile,"\n");
   fprintf(outfile,"                       OMP2 (OO-MP2)   \n");
   fprintf(outfile,"              Program Written by Ugur Bozkaya\n") ; 
   fprintf(outfile,"              Latest Revision December 30, 2011\n") ;
   fprintf(outfile,"\n");
   fprintf(outfile,"              U. Bozkaya, J. M. Turney, Y. Yamaguchi, H. F. Schaefer,  \n") ;
   fprintf(outfile,"              and C. D. Sherrill, J. Chem. Phys. 135, 104103 (2011). \n") ;
   fprintf(outfile,"\n");
   fprintf(outfile," ============================================================================== \n");
   fprintf(outfile," ============================================================================== \n");
   fprintf(outfile," ============================================================================== \n");
   fprintf(outfile,"\n");
   fflush(outfile);

}//


double OMP2Wave::compute_energy()
{   
        
	// Warnings 
	if (nfrzc != 0 || nfrzv != 0) {
	  fprintf(stderr,  "\tThe OMP2 method has been implemented for only all-electron computations, yet.\n");
	  fprintf(outfile, "\tThe OMP2 method has been implemented for only all-electron computations, yet.\n");
	  fflush(outfile);
          return EXIT_FAILURE;
	}
	
	mo_optimized = 0;
	trans_ints();  
	t2_1st_sc();
	ref_energy();
	mp2_energy();
	Emp2L=Emp2;
        EcorrL=Emp2L-Escf;
	Emp2L_old=Emp2;
	
	fprintf(outfile,"\n \n"); 
	fprintf(outfile,"\tComputing MP2 energy using SCF MOs... \n"); 
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %12.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %12.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %12.14f\n", Eref);
	fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %12.14f\n", Emp2AA);
	fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %12.14f\n", Emp2AB);
	fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %12.14f\n", Emp2BB);
	fprintf(outfile,"\tMP2 Correlation Energy (a.u.)      : %12.14f\n", Ecorr);
	fprintf(outfile,"\tMP2 Total Energy (a.u.)            : %12.14f\n", Emp2);
	fprintf(outfile,"\tScaled_SS Correlation Energy (a.u.): %12.14f\n", Escsmp2AA+Escsmp2BB);
	fprintf(outfile,"\tScaled_OS Correlation Energy (a.u.): %12.14f\n", Escsmp2AB);
	fprintf(outfile,"\tSCS-MP2 Total Energy (a.u.)        : %12.14f\n", Escsmp2);
	fprintf(outfile,"\tSOS-MP2 Total Energy (a.u.)        : %12.14f\n", Esosmp2);
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\n"); 
	fflush(outfile);

	response_pdms();
	GFockmo();
	idp();
	mograd();
        occ_iterations();
	
        if (rms_wog == 0.0 && fabs(DE) >= tol_Eod) {
	  semi_canonic();
	  t2_1st_sc();
        }     

  if (conver == 1) {
        ref_energy();
	mp2_energy();
	
	fprintf(outfile,"\n"); 
	fprintf(outfile,"\tComputing MP2 energy using optimized MOs... \n"); 
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %12.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %12.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %12.14f\n", Eref);
	fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %12.14f\n", Emp2AA);
	fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %12.14f\n", Emp2AB);
	fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %12.14f\n", Emp2BB);
	fprintf(outfile,"\tMP2 Correlation Energy (a.u.)      : %12.14f\n", Ecorr);
	fprintf(outfile,"\tMP2 Total Energy (a.u.)            : %12.14f\n", Emp2);
	fprintf(outfile,"\tScaled_SS Correlation Energy (a.u.): %12.14f\n", Escsmp2AA+Escsmp2BB);
	fprintf(outfile,"\tScaled_OS Correlation Energy (a.u.): %12.14f\n", Escsmp2AB);
	fprintf(outfile,"\tSCS-MP2 Total Energy (a.u.)        : %12.14f\n", Escsmp2);
	fprintf(outfile,"\tSOS-MP2 Total Energy (a.u.)        : %12.14f\n", Esosmp2);
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\n"); 
	fflush(outfile);


	fprintf(outfile,"\n");
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\t================ CONGRATULATIONS OMP2 ITERATIONS ARE CONVERGED =============== \n");
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %12.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %12.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %12.14f\n", Eref);
	fprintf(outfile,"\tOMP2 Correlation Energy (a.u.)     : %12.14f\n", Emp2L-Escf);
	fprintf(outfile,"\tEomp2 - Eref (a.u.)                : %12.14f\n", Emp2L-Eref);
	fprintf(outfile,"\tOMP2 Total Energy (a.u.)           : %12.14f\n", Emp2L);
	fprintf(outfile,"\tSCS-OMP2 Total Energy (a.u.)       : %12.14f\n", Escsmp2);
	fprintf(outfile,"\tSOS-OMP2 Total Energy (a.u.)       : %12.14f\n", Esosmp2);
	fprintf(outfile,"\t============================================================================ \n");
	fprintf(outfile,"\n");
	fflush(outfile);

	chkpt_->wt_etot(Emp2L);
	chkpt_->wt_emp2(Emp2L);
	chkpt_->wt_ecorr(Emp2L-Escf);
	chkpt_->wt_eref(Eref);
	
	if (do_scs == "TRUE") {
	  chkpt_->wt_etot(Escsmp2);
	  chkpt_->wt_ecorr(Escsmp2-Escf);
	}
    
	else if (do_sos == "TRUE") {
	  chkpt_->wt_etot(Esosmp2);
	  chkpt_->wt_ecorr(Esosmp2-Escf);
	}
 
	if (natorb == "TRUE") nbo();
	if (omp2_orb_energy == "TRUE") semi_canonic(); 
	
	// Write MO coefficients to Cmo.psi
        if (write_mo_coeff == "TRUE"){
	  fprintf(outfile,"\n\tWriting MO coefficients in pitzer order to external files CmoA.psi and CmoB.psi...\n");  
	  fflush(outfile);
	  double **C_pitzerA = block_matrix(nso,nmo);
	  double **C_pitzerB = block_matrix(nso,nmo);
	  memset(C_pitzerA[0], 0, sizeof(int)*nso*nmo);
	  memset(C_pitzerB[0], 0, sizeof(int)*nso*nmo);
    
	  //set C_pitzer
	  C_pitzerA = Ca_->to_block_matrix();    
	  C_pitzerB = Cb_->to_block_matrix();    
	
	  // write binary data
	  ofstream OutFile1;
	  OutFile1.open("CmoA.psi", ios::out | ios::binary);
	  OutFile1.write( (char*)C_pitzerA[0], sizeof(double)*nso*nmo);
	  OutFile1.close();
	  
	  // write binary data
	  ofstream OutFile2;
	  OutFile2.open("CmoB.psi", ios::out | ios::binary);
	  OutFile2.write( (char*)C_pitzerB[0], sizeof(double)*nso*nmo);
	  OutFile2.close();  
	  
	  free_block(C_pitzerA);
	  free_block(C_pitzerB);
	}     
	
  }// end if (conver == 1)

        mem_release();
	if (do_scs == "TRUE") return Escsmp2;
	else if (do_sos == "TRUE") return Esosmp2;
        else return Emp2L;

} // end of compute_energy



void OMP2Wave::ref_energy()
{
     double Ehf;     
     Ehf=0.0;
     
     // alpha contribution
     for (int h=0; h<nirreps; h++){
      for (int i=0; i<occpiA[h];i++) {
	Ehf+=HmoA->get(h,i,i) + FockA->get(h,i,i);
      }
    }  
    
    // beta contribution
     for (int h=0; h<nirreps; h++){
      for (int i=0; i<occpiB[h];i++) {
	Ehf+=HmoB->get(h,i,i) + FockB->get(h,i,i);
      }
    }  
    
    Eref=(0.5*Ehf)+Enuc; 
    
} // end of ref_energy


void OMP2Wave::mp2_energy()
{
     dpdbuf4 K, T;
     
     psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
     psio_->open(PSIF_OMP2_DPD, PSIO_OPEN_OLD);
     
     Ecorr=0.0;
     Escsmp2AA=0.0;
     Escsmp2AB=0.0;
     Escsmp2BB=0.0;
     Escsmp2=0.0;
     Esosmp2AB=0.0;
     
     // Compute Energy
     // Alpha-Alpha spin contribution
     dpd_buf4_init(&T, PSIF_OMP2_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO||VV>");
     Ecorr += 0.25 * dpd_buf4_dot(&T, &K);     
     dpd_buf4_close(&T);
     dpd_buf4_close(&K);
     
     Emp2AA=Ecorr;    
     Escsmp2AA = ss_scale * Emp2AA; 
     
     
     // Alpha-Beta spin contribution
     dpd_buf4_init(&T, PSIF_OMP2_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                 ID("[O,o]"), ID("[V,v]"), 0, "T2_1 <Oo|Vv>");
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
     Ecorr += dpd_buf4_dot(&T, &K);     
     dpd_buf4_close(&T);
     dpd_buf4_close(&K);
     
     Emp2AB = Ecorr - Emp2AA;
     Escsmp2AB = os_scale * Emp2AB;  
     if (mo_optimized == 0) Esosmp2AB = sos_scale * Emp2AB; 
     else if (mo_optimized == 1) Esosmp2AB = sos_scale2 * Emp2AB;  
     
     
     // Beta-Beta spin contribution
     dpd_buf4_init(&T, PSIF_OMP2_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1 <oo|vv>");
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo||vv>");
     Ecorr += 0.25 * dpd_buf4_dot(&T, &K);     
     dpd_buf4_close(&T);
     dpd_buf4_close(&K);
     
     Emp2BB = Ecorr - Emp2AA - Emp2AB;  
     Escsmp2BB = ss_scale * Emp2BB;  
     
     Emp2 = Eref + Ecorr;
     Escsmp2 = Eref + Escsmp2AA + Escsmp2AB + Escsmp2BB;
     Esosmp2 = Eref + Esosmp2AB;     
     
     psio_->close(PSIF_LIBTRANS_DPD, 1);
     psio_->close(PSIF_OMP2_DPD, 1);    
     
     Emp2=Eref+Ecorr;
     
} // end of mp2_energy


void OMP2Wave::mp2l_energy()
{      
     //fprintf(outfile,"\n mp2l_energy is starting... \n"); fflush(outfile);
     
     // One-electron contribution
     /*
     double Etpdm_oooo, Etpdm_oovv, Etpdm_ovov, Etpdm_vovo;
     HG1A->zero();
     HG1B->zero();
     HG1A->gemm(false, false, 1.0, HmoA, g1symmA, 0.0);
     HG1B->gemm(false, false, 1.0, HmoB, g1symmB, 0.0);
     Emp2_adm = 0.0;
     Emp2_adm = HG1A->trace() + HG1B->trace() + Enuc;
     Eopdm = HG1A->trace() + HG1B->trace();
     fprintf(outfile,"\n OPDM Contribution               : %12.14f\n", Eopdm); fflush(outfile);
     */
     
    // Two-electron contribution
    dpdbuf4 G, K;
    
    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);  
    psio_->open(PSIF_OMP2_DENSITY, PSIO_OPEN_OLD); 
    
    // OOOO-Block contribution
    // E += G_IJKL <IJ||KL>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "MO Ints <OO||OO>");
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "TPDM <OO|OO>");
    Emp2_adm += dpd_buf4_dot(&G, &K);         
    dpd_buf4_close(&K);
    dpd_buf4_close(&G);   
    
    // E += G_ijkl <ij||kl>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "MO Ints <oo||oo>");
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "TPDM <oo|oo>");
    Emp2_adm += dpd_buf4_dot(&G, &K);  
    dpd_buf4_close(&K);
    dpd_buf4_close(&G);  
    
    // E += 4*G_IjKl <Ij||Kl>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "MO Ints <Oo|Oo>");
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "TPDM <Oo|Oo>");
    Emp2_adm += 4.0 * dpd_buf4_dot(&G, &K); 
    dpd_buf4_close(&K);
    dpd_buf4_close(&G);
    //Etpdm_oooo = Emp2_adm - Enuc - Eopdm;   
    //fprintf(outfile,"\n OOOO-Block Contribution         : %12.14f\n", Etpdm_oooo); fflush(outfile);
    
    
    // OOVV-Block contribution
    // E += 2*G_IJAB <IJ||AB>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO||VV>");
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
    Emp2_adm += 2.0 * dpd_buf4_dot(&G, &K);   
    dpd_buf4_close(&K);
    dpd_buf4_close(&G);    
    
    // E += 2*G_ijab <ij|ab>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo||vv>");
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "TPDM <oo|vv>");
    Emp2_adm += 2.0 * dpd_buf4_dot(&G, &K);   
    dpd_buf4_close(&K);
    dpd_buf4_close(&G); 
    
    // E += 8*G_IjAb <Ij||Ab>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "TPDM <Oo|Vv>");
    Emp2_adm += 8.0 * dpd_buf4_dot(&G, &K);   
    dpd_buf4_close(&K);
    dpd_buf4_close(&G); 
    //Etpdm_oovv = Emp2_adm - Enuc - Eopdm - Etpdm_oooo;
    //fprintf(outfile,"\n OOVV-Block Contribution         : %12.14f\n", Etpdm_oovv); fflush(outfile);

    
    // OVOV-Block contribution
    // E += 4*G_IAJB <IA||JB>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV||OV>");
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
    Emp2_adm += 4.0 * dpd_buf4_dot(&G, &K);     
    dpd_buf4_close(&K);
    dpd_buf4_close(&G);
    
    // E += 4*G_iajb <ia||jb>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints <ov||ov>");
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "TPDM <ov|ov>");
    Emp2_adm += 4.0 * dpd_buf4_dot(&G, &K);     
    dpd_buf4_close(&K);
    dpd_buf4_close(&G);
    
    // E += 4*G_IaJb <Ia||Jb>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "MO Ints <Ov|Ov>");
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "TPDM <Ov|Ov>");
    Emp2_adm += 4.0 * dpd_buf4_dot(&G, &K);     
    dpd_buf4_close(&K);
    dpd_buf4_close(&G);
    
    // VOVO-Block contribution
    // E += 4*G_AiBj <Ai||Bj>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0,ID("[V,o]"), ID("[V,o]"),
                  ID("[V,o]"), ID("[V,o]"), 0, "MO Ints <Vo|Vo>");
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[V,o]"), ID("[V,o]"),
                  ID("[V,o]"), ID("[V,o]"), 0, "TPDM <Vo|Vo>"); 
    Emp2_adm += 4.0 * dpd_buf4_dot(&G, &K);     
    dpd_buf4_close(&K);
    dpd_buf4_close(&G);
    //Etpdm_ovov = Emp2_adm - Enuc - Eopdm - Etpdm_oooo - Etpdm_oovv;
    //fprintf(outfile,"\n OVOV-Block Contribution         : %12.14f\n", Etpdm_ovov); fflush(outfile);
    

    psio_->close(PSIF_LIBTRANS_DPD, 1);  
    psio_->close(PSIF_OMP2_DENSITY, 1);  
  
    EcorrL=Emp2_adm-Escf;        
    Emp2L=Emp2_adm;  
    DE = Emp2L - Emp2L_old;
    Emp2L_old = Emp2L;
    

    //fprintf(outfile,"\n MP2 Total Energy via pdms       : %12.14f\n", Emp2_adm); fflush(outfile);    
    //fprintf(outfile,"\n mp2l_energy is done... \n"); fflush(outfile);
    
} // end of mp2l_energy

void OMP2Wave::nbo()
{
      
fprintf(outfile,"\n  \n");
fprintf(outfile," ============================================================================== \n");
fprintf(outfile," ======================== NBO ANALYSIS ======================================== \n");
fprintf(outfile," ============================================================================== \n");
fprintf(outfile,"\n Diagonalizing one-particle response density matrix... \n");
fprintf(outfile,"\n");
fflush(outfile);
 
      SharedMatrix Udum = boost::shared_ptr<Matrix>(new Matrix("Udum", nirreps, mopi, mopi));
      SharedVector diag = boost::shared_ptr<Vector>(new Vector("Natural orbital occupation numbers", nirreps, mopi));

      // Diagonalizing Alpha-OPDM
      Udum->zero();
      
      //diag->zero();
      for(int h = 0; h < nirreps; h++){
	  for(int i = 0; i < mopi[h]; i++){
	    diag->set(h,i,0.0);
	  }
	}
   
      g1symmA->diagonalize(Udum, diag);
	
      //trace
      //sum=diag->trace();
      sum=0.0;
      for(int h = 0; h < nirreps; h++){
	  for(int i = 0; i < mopi[h]; i++){
	    sum+=diag->get(h,i);
	  }
	}      
      
      fprintf(outfile, "\n Trace of alpha one-particle density matrix: %20.14f \n",  sum);
      fprintf(outfile,"\n");
      fflush(outfile);

      //print
      diag->print();      
      
      
      
      
      // Diagonalizing Beta-OPDM
      Udum->zero();
      
      //diag->zero();
      for(int h = 0; h < nirreps; h++){
	  for(int i = 0; i < mopi[h]; i++){
	    diag->set(h,i,0.0);
	  }
	}
   
      g1symmB->diagonalize(Udum, diag);
	
      //trace
      //sum=diag->trace();
      sum=0.0;
      for(int h = 0; h < nirreps; h++){
	  for(int i = 0; i < mopi[h]; i++){
	    sum+=diag->get(h,i);
	  }
	}      
      
      fprintf(outfile, "\n Trace of beta one-particle density matrix: %20.14f \n",  sum);
      fprintf(outfile,"\n");
      fflush(outfile);

      //print
      diag->print(); 
      
} // end of nbo

double OMP2Wave::dot_product(int dim, double *x, double *y)
{
  double value;
  
  value = 0.0;
  for (int i=0; i < dim; i++) {
    value += x[i] * y[i];
  }
  
  return value; 
    
} // end of dot_product


void OMP2Wave::mem_release()
{   

	delete ints;
	delete [] idprowA;
	delete [] idprowB;
	delete [] idpcolA;
	delete [] idpcolB;
	delete [] idpirrA;
	delete [] idpirrB;
	delete [] wogA;
	delete [] wogB;
	delete [] kappaA;
	delete [] kappaB;
	delete [] kappa_barA;
	delete [] kappa_barB;
	delete [] evalsA;
	delete [] evalsB;
	delete [] evals_c1A;
	delete [] evals_c1B;
	delete [] pitzer2symblk;
	delete [] pitzer2symirrep;
	delete [] c1topitzerA;
	delete [] c1topitzerB;
	delete [] pitzer2c1A;
	delete [] pitzer2c1B;
	delete [] PitzerOffset;
	delete [] sosym;
	delete [] mosym;
	delete [] mosym_c1;
	delete [] occ_offA;
	delete [] occ_offB;
	delete [] vir_offA;
	delete [] vir_offB;
	delete [] occ2symblkA;
	delete [] occ2symblkB;
	delete [] virt2symblkA;
	delete [] virt2symblkB;
	/*
	delete [] c1toqtA;
	delete [] c1toqtB;
	delete [] qt2c1A;
	delete [] qt2c1B;
	*/
	
	if (opt_method == "DIIS") {
	  free_block(errvecsA);
	  free_block(errvecsB);
	  free_block(vecsA);
	  free_block(vecsB);
	}
	
	chkpt_.reset();
	Ca_.reset();
	Cb_.reset();
	Ca_ref.reset();
	Cb_ref.reset();
	Hso.reset();
	Tso.reset();
	Vso.reset();
	HmoA.reset();
	HmoB.reset();
	FockA.reset();
	FockB.reset();
	gamma1corrA.reset();
	gamma1corrB.reset();
	g1symmA.reset();
	g1symmB.reset();
	GFockA.reset();
	GFockB.reset();
	UorbA.reset();
	UorbB.reset();
	KorbA.reset();
	KorbB.reset();
	KsqrA.reset();
	KsqrB.reset();
	HG1A.reset();
	HG1B.reset();
	WorbA.reset();
	WorbB.reset();
	GooA.reset();
	GooB.reset();
	GvvA.reset();
	GvvB.reset();

	//fprintf(outfile,"\n mem_release done. \n"); fflush(outfile);

}//


} }

