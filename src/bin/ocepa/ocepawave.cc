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
#include <libqt/qt.h>

#include "ocepawave.h"
#include "defines.h"

using namespace psi;
using namespace boost;

namespace psi { namespace ocepawave {

OCEPAWave::OCEPAWave(boost::shared_ptr<Wavefunction> reference_wavefunction, Options &options)
    : Wavefunction(options, _default_psio_lib_)
{
    reference_wavefunction_ = reference_wavefunction;
    common_init();
}//

OCEPAWave::~OCEPAWave()
{
}//

/*=======================*/
/*  common_init()        */
/*=======================*/
void OCEPAWave::common_init()
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
	sos_scale=options_.get_double("MP2_SOS_SCALE");
	cepa_os_scale_=options_.get_double("CEPA_OS_SCALE");
	cepa_ss_scale_=options_.get_double("CEPA_SS_SCALE");
	cepa_sos_scale_=options_.get_double("CEPA_SOS_SCALE");
	
	orth_type=options_.get_str("ORTH_TYPE");
	opt_method=options_.get_str("OPT_METHOD");
	hess_type=options_.get_str("HESS_TYPE");
	ocepa_orb_energy=options_.get_str("OCEPA_ORBS_PRINT");
	natorb=options_.get_str("NAT_ORBS");
	reference=options_.get_str("REFERENCE");
	do_scs=options_.get_str("DO_SCS");
	do_sos=options_.get_str("DO_SOS");
	write_mo_coeff=options_.get_str("MO_WRITE");
	read_mo_coeff=options_.get_str("MO_READ");
        lineq=options_.get_str("LINEQ_SOLVER"); 
        compute_cepal=options_.get_str("CEPAL_ENERGY"); 
	level_shift=options_.get_str("DO_LEVEL_SHIFT");
        twopdm_abcd_type=options_.get_str("TPDM_ABCD_TYPE");
	//scs_type_=options_.get_str("SCS_TYPE");
	//sos_type_=options_.get_str("SOS_TYPE");
	wfn_type_=options_.get_str("WFN_TYPE");
	dertype=options_.get_str("DERTYPE");

        if (reference == "RHF" || reference == "RKS") reference_ = "RESTRICTED";
        else if (reference == "UHF" || reference == "UKS") reference_ = "UNRESTRICTED";

	cutoff = pow(10.0,-exp_cutoff);
	if (print_ > 0) options_.print();
        title();
	get_moinfo();
        time4grad = 0;// means i will not compute the gradient 
	

if (reference_ == "RESTRICTED") {
	// Memory allocation
	HmoA = boost::shared_ptr<Matrix>(new Matrix("MO-basis alpha one-electron ints", nirrep_, nmopi_, nmopi_));
	FockA = boost::shared_ptr<Matrix>(new Matrix("MO-basis alpha Fock matrix", nirrep_, nmopi_, nmopi_));
	gamma1corr = boost::shared_ptr<Matrix>(new Matrix("MO-basis correlation OPDM", nirrep_, nmopi_, nmopi_));
	g1symm = boost::shared_ptr<Matrix>(new Matrix("MO-basis OPDM", nirrep_, nmopi_, nmopi_));
	GFock = boost::shared_ptr<Matrix>(new Matrix("MO-basis generalized Fock matrix", nirrep_, nmopi_, nmopi_));
	UorbA = boost::shared_ptr<Matrix>(new Matrix("Alpha MO rotation matrix", nirrep_, nmopi_, nmopi_));
	KorbA = boost::shared_ptr<Matrix>(new Matrix("K alpha MO rotation", nirrep_, nmopi_, nmopi_)); 
	KsqrA = boost::shared_ptr<Matrix>(new Matrix("K^2 alpha MO rotation", nirrep_, nmopi_, nmopi_)); 
	HG1 = boost::shared_ptr<Matrix>(new Matrix("h*g1symm", nirrep_, nmopi_, nmopi_));
	WorbA = boost::shared_ptr<Matrix>(new Matrix("Alpha MO gradient matrix", nirrep_, nmopi_, nmopi_));
	GooA = boost::shared_ptr<Matrix>(new Matrix("Alpha Goo intermediate", nirrep_, aoccpiA, aoccpiA));
	GvvA = boost::shared_ptr<Matrix>(new Matrix("Alpha Gvv intermediate", nirrep_, avirtpiA, avirtpiA));


        Molecule& mol = *reference_wavefunction_->molecule().get();
        CharacterTable ct = mol.point_group()->char_table();
        fprintf(outfile,"\tMO spaces per irreps... \n\n"); fflush(outfile);
        fprintf(outfile, "\tIRREP   FC    OCC   VIR  FV \n");
        fprintf(outfile, "\t==============================\n");                                                 
        for(int h = 0; h < nirrep_; ++h){
         fprintf(outfile, "\t %3s   %3d   %3d   %3d  %3d\n",
                             ct.gamma(h).symbol(), frzcpi_[h], aoccpiA[h], avirtpiA[h], frzvpi_[h]);
        }
        fprintf(outfile,     "\t==============================\n"); 
	fflush(outfile);

    // Alloc ints
    std::vector<boost::shared_ptr<MOSpace> > spaces;
    spaces.push_back(MOSpace::occ);
    spaces.push_back(MOSpace::vir);

    ints = new IntegralTransform(reference_wavefunction_, spaces, 
                           IntegralTransform::Restricted,
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

}  // end if (reference_ == "RESTRICTED") 

    
else if (reference_ == "UNRESTRICTED") {
	// Memory allocation
	HmoA = boost::shared_ptr<Matrix>(new Matrix("MO-basis alpha one-electron ints", nirrep_, nmopi_, nmopi_));
	FockA = boost::shared_ptr<Matrix>(new Matrix("MO-basis alpha Fock matrix", nirrep_, nmopi_, nmopi_));
	gamma1corrA = boost::shared_ptr<Matrix>(new Matrix("MO-basis alpha correlation OPDM", nirrep_, nmopi_, nmopi_));
	g1symmA = boost::shared_ptr<Matrix>(new Matrix("MO-basis alpha OPDM", nirrep_, nmopi_, nmopi_));
	GFockA = boost::shared_ptr<Matrix>(new Matrix("MO-basis alpha generalized Fock matrix", nirrep_, nmopi_, nmopi_));
	UorbA = boost::shared_ptr<Matrix>(new Matrix("Alpha MO rotation matrix", nirrep_, nmopi_, nmopi_));
	KorbA = boost::shared_ptr<Matrix>(new Matrix("K alpha MO rotation", nirrep_, nmopi_, nmopi_)); 
	KsqrA = boost::shared_ptr<Matrix>(new Matrix("K^2 alpha MO rotation", nirrep_, nmopi_, nmopi_)); 
	HG1A = boost::shared_ptr<Matrix>(new Matrix("Alpha h*g1symm", nirrep_, nmopi_, nmopi_));
	WorbA = boost::shared_ptr<Matrix>(new Matrix("Alpha MO gradient matrix", nirrep_, nmopi_, nmopi_));
	GooA = boost::shared_ptr<Matrix>(new Matrix("Alpha Goo intermediate", nirrep_, aoccpiA, aoccpiA));
	GvvA = boost::shared_ptr<Matrix>(new Matrix("Alpha Gvv intermediate", nirrep_, avirtpiA, avirtpiA));

	HmoB = boost::shared_ptr<Matrix>(new Matrix("MO-basis beta one-electron ints", nirrep_, nmopi_, nmopi_));
	FockB = boost::shared_ptr<Matrix>(new Matrix("MO-basis beta Fock matrix", nirrep_, nmopi_, nmopi_));
	gamma1corrB = boost::shared_ptr<Matrix>(new Matrix("MO-basis beta correlation OPDM", nirrep_, nmopi_, nmopi_));
	g1symmB = boost::shared_ptr<Matrix>(new Matrix("MO-basis beta OPDM", nirrep_, nmopi_, nmopi_));
	GFockB = boost::shared_ptr<Matrix>(new Matrix("MO-basis beta generalized Fock matrix", nirrep_, nmopi_, nmopi_));
	UorbB = boost::shared_ptr<Matrix>(new Matrix("Beta MO rotation matrix", nirrep_, nmopi_, nmopi_));
	KorbB = boost::shared_ptr<Matrix>(new Matrix("K beta MO rotation", nirrep_, nmopi_, nmopi_)); 
	KsqrB = boost::shared_ptr<Matrix>(new Matrix("K^2 beta MO rotation", nirrep_, nmopi_, nmopi_)); 
	HG1B = boost::shared_ptr<Matrix>(new Matrix("Beta h*g1symm", nirrep_, nmopi_, nmopi_));
	WorbB = boost::shared_ptr<Matrix>(new Matrix("Beta MO gradient matrix", nirrep_, nmopi_, nmopi_));
	GooB = boost::shared_ptr<Matrix>(new Matrix("Beta Goo intermediate", nirrep_, aoccpiB, aoccpiB));
	GvvB = boost::shared_ptr<Matrix>(new Matrix("Beta Gvv intermediate", nirrep_, avirtpiB, avirtpiB));

        Molecule& mol = *reference_wavefunction_->molecule().get();
        CharacterTable ct = mol.point_group()->char_table();
        fprintf(outfile,"\tMO spaces per irreps... \n\n"); fflush(outfile);
        fprintf(outfile, "\tIRREP   FC   AOCC  BOCC  AVIR    BVIR  FV \n");
        fprintf(outfile, "\t==========================================\n");                                                 
        for(int h = 0; h < nirrep_; ++h){
         fprintf(outfile, "\t %3s   %3d   %3d   %3d   %3d    %3d   %3d\n",
                             ct.gamma(h).symbol(), frzcpi_[h], aoccpiA[h], aoccpiB[h], avirtpiA[h], avirtpiB[h], frzvpi_[h]);
        }
        fprintf(outfile,     "\t==========================================\n\n");
	fflush(outfile);

	
    // Alloc ints
    std::vector<boost::shared_ptr<MOSpace> > spaces;
    spaces.push_back(MOSpace::occ);
    spaces.push_back(MOSpace::vir);

    // Perform MO integral transformation on spaces
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

}// end if (reference_ == "UNRESTRICTED") 
}// end common_init


/*=======================*/
/*  title()              */
/*=======================*/
void OCEPAWave::title()
{
   fprintf(outfile,"\n");
   fprintf(outfile," ============================================================================== \n");
   fprintf(outfile," ============================================================================== \n");
   fprintf(outfile," ============================================================================== \n");
   fprintf(outfile,"\n");
   fprintf(outfile,"                       OCEPA (OO-CEPA)   \n");
   fprintf(outfile,"              Program Written by Ugur Bozkaya,\n") ; 
   fprintf(outfile,"              Latest Revision November 03, 2012.\n") ;
   fprintf(outfile,"\n");
   fprintf(outfile," ============================================================================== \n");
   fprintf(outfile," ============================================================================== \n");
   fprintf(outfile," ============================================================================== \n");
   fprintf(outfile,"\n");
   fflush(outfile);
}//


/*=======================*/
/*  compute_energy()     */
/*=======================*/
double OCEPAWave::compute_energy()
{   
        
	// Warnings 
	if (nfrzc != 0 || nfrzv != 0) {
	  fprintf(stderr,  "\tThe OCEPA method has been implemented for only all-electron computations, yet.\n");
	  fprintf(outfile, "\tThe OCEPA method has been implemented for only all-electron computations, yet.\n");
	  fflush(outfile);
          return EXIT_FAILURE;
	}
	
	mo_optimized = 0;
        timer_on("trans_ints");
	if (reference_ == "RESTRICTED") trans_ints_rhf();  
	else if (reference_ == "UNRESTRICTED") trans_ints_uhf();  
        timer_off("trans_ints");
        timer_on("REF Energy");
	ref_energy();
        timer_off("REF Energy");
        timer_on("T2(1)");
	t2_1st_sc();
        timer_off("T2(1)");
        timer_on("MP2 Energy");
	mp2_energy();
        timer_off("MP2 Energy");
        Ecepa = Emp2;
	EcepaL = Ecepa;
        EcorrL = Ecorr;
	EcepaL_old = Ecepa;

	fprintf(outfile,"\n"); 
	fprintf(outfile,"\tComputing MP2 energy using SCF MOs (Canonical MP2)... \n"); 
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %12.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %12.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %12.14f\n", Eref);
	fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %12.14f\n", Emp2AA);
	fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %12.14f\n", Emp2AB);
	fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %12.14f\n", Emp2BB);
	fprintf(outfile,"\tScaled_SS Correlation Energy (a.u.): %12.14f\n", Escsmp2AA+Escsmp2BB);
	fprintf(outfile,"\tScaled_OS Correlation Energy (a.u.): %12.14f\n", Escsmp2AB);
	fprintf(outfile,"\tSCS-MP2 Total Energy (a.u.)        : %12.14f\n", Escsmp2);
	fprintf(outfile,"\tSOS-MP2 Total Energy (a.u.)        : %12.14f\n", Esosmp2);
	fprintf(outfile,"\tSCSN-MP2 Total Energy (a.u.)       : %12.14f\n", Escsnmp2);
	fprintf(outfile,"\tSCS-MI-MP2 Total Energy (a.u.)     : %12.14f\n", Escsmimp2);
	fprintf(outfile,"\tSCS-MP2-VDW Total Energy (a.u.)    : %12.14f\n", Escsmp2vdw);
	fprintf(outfile,"\tSOS-PI-MP2 Total Energy (a.u.)     : %12.14f\n", Esospimp2);
	fprintf(outfile,"\tMP2 Correlation Energy (a.u.)      : %12.14f\n", Ecorr);
	fprintf(outfile,"\tMP2 Total Energy (a.u.)            : %12.14f\n", Emp2);
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\n"); 
	fflush(outfile);
	Process::environment.globals["MP2 TOTAL ENERGY"] = Emp2;
	Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
	Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
	Process::environment.globals["SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;
	Process::environment.globals["SCS-MI-MP2 TOTAL ENERGY"] = Escsmimp2;
	Process::environment.globals["SCS-MP2-VDW TOTAL ENERGY"] = Escsmp2vdw;
	Process::environment.globals["SOS-PI-MP2 TOTAL ENERGY"] = Esospimp2;

if (wfn_type_ == "OCEPA") {
	response_pdms();
	GFockmo();
	idp();
	mograd();
        occ_iterations();
	
        if (rms_wog == 0.0 && fabs(DE) >= tol_Eod) {
	  semi_canonic();
	  t2_amps();
        }     

  if (conver == 1) {
        ref_energy();
	cepa_energy();
        
        /*
        fprintf(outfile,"\n"); 
	fprintf(outfile,"\tComputing CEPA energy using optimized MOs... \n"); 
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %12.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %12.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %12.14f\n", Eref);
	fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %12.14f\n", EcepaAA);
	fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %12.14f\n", EcepaAB);
	fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %12.14f\n", EcepaBB);
	fprintf(outfile,"\tSCS-CEPA(0) Total Energy (a.u.)    : %12.14f\n", Escscepa);
	fprintf(outfile,"\tSOS-CEPA(0) Total Energy (a.u.)    : %12.14f\n", Esoscepa);
	fprintf(outfile,"\tCEPA(0) Correlation Energy (a.u.)  : %12.14f\n", Ecepa-Escf);
	fprintf(outfile,"\tEcepa - Eref (a.u.)                : %12.14f\n", Ecepa-Eref);
	fprintf(outfile,"\tCEPA(0) Total Energy (a.u.)        : %12.14f\n", Ecepa);
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\n"); 
	fflush(outfile);
        */

	fprintf(outfile,"\n");
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\t================ OCEPA FINAL RESULTS ========================================= \n");
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %12.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %12.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %12.14f\n", Eref);
	fprintf(outfile,"\tSCS-OCEPA(0) Total Energy (a.u.)   : %12.14f\n", Escscepa);
	fprintf(outfile,"\tSOS-OCEPA(0) Total Energy (a.u.)   : %12.14f\n", Esoscepa);
	fprintf(outfile,"\tOCEPA(0) Correlation Energy (a.u.) : %12.14f\n", EcepaL-Escf);
	fprintf(outfile,"\tEocepa - Eref (a.u.)               : %12.14f\n", EcepaL-Eref);
	fprintf(outfile,"\tOCEPA(0) Total Energy (a.u.)       : %12.14f\n", EcepaL);
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\n");
	fflush(outfile);
	
	// Set the global variables with the energies
	Process::environment.globals["OCEPA(0) TOTAL ENERGY"] = EcepaL;
	Process::environment.globals["SCS-OCEPA(0) TOTAL ENERGY"] =  Escscepa;
	Process::environment.globals["SOS-OCEPA(0) TOTAL ENERGY"] =  Esoscepa;
	Process::environment.globals["CURRENT ENERGY"] = EcepaL;
	Process::environment.globals["CURRENT REFERENCE ENERGY"] = Eref;
	Process::environment.globals["CURRENT CORRELATION ENERGY"] = EcepaL-Escf;

        // if scs on	
	if (do_scs == "TRUE") {
	    Process::environment.globals["CURRENT ENERGY"] = Escscepa;

	}
    
        // else if sos on	
	else if (do_sos == "TRUE") {
	         Process::environment.globals["CURRENT ENERGY"] = Esoscepa;
	}

	if (natorb == "TRUE") nbo();
	if (ocepa_orb_energy == "TRUE") semi_canonic(); 
	
	// Write MO coefficients to Cmo.psi
        if (write_mo_coeff == "TRUE"){
	  fprintf(outfile,"\n\tWriting MO coefficients in pitzer order to external file CmoA.psi...\n");  
	  fflush(outfile);
	  double **C_pitzerA = block_matrix(nso_,nmo_);
	  memset(C_pitzerA[0], 0, sizeof(double)*nso_*nmo_);
    
	  //set C_pitzer
	  C_pitzerA = Ca_->to_block_matrix();    
	
	  // write binary data
	  ofstream OutFile1;
	  OutFile1.open("CmoA.psi", ios::out | ios::binary);
	  OutFile1.write( (char*)C_pitzerA[0], sizeof(double)*nso_*nmo_);
	  OutFile1.close();
	  free_block(C_pitzerA);
	 
          if (reference_ == "UNRESTRICTED" ) { 
	      fprintf(outfile,"\n\tWriting MO coefficients in pitzer order to external file CmoB.psi...\n");  
	      fflush(outfile);
	      double **C_pitzerB = block_matrix(nso_,nmo_);
	      memset(C_pitzerB[0], 0, sizeof(double)*nso_*nmo_);

	      //set C_pitzer
	      C_pitzerB = Cb_->to_block_matrix();    

	      // write binary data
	      ofstream OutFile2;
	      OutFile2.open("CmoB.psi", ios::out | ios::binary);
	      OutFile2.write( (char*)C_pitzerB[0], sizeof(double)*nso_*nmo_);
	      OutFile2.close();  
	      free_block(C_pitzerB);
          }
	}     
	
        // Compute Analytic Gradients
        if (dertype == "FIRST") {
            time4grad = 1;
            coord_grad();
        }

  }// end if (conver == 1)
}// end if wfn_type_ == "OCEPA"


else if (wfn_type_ == "CEPA") {
        // Perform CEPA iterations
        cc_iterations();

	fprintf(outfile,"\n");
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\t================ CEPA FINAL RESULTS ========================================== \n");
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %12.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %12.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %12.14f\n", Eref);
	fprintf(outfile,"\tCEPA(0) Correlation Energy (a.u.)  : %12.14f\n", Ecorr);
	fprintf(outfile,"\tCEPA(0) Total Energy (a.u.)        : %12.14f\n", Ecepa);
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\n");
	fflush(outfile);
	
	// Set the global variables with the energies
	Process::environment.globals["CEPA(0) TOTAL ENERGY"] = Ecepa;
	Process::environment.globals["CURRENT ENERGY"] = Ecepa;
	Process::environment.globals["CURRENT REFERENCE ENERGY"] = Eref;
	Process::environment.globals["CURRENT CORRELATION ENERGY"] = Ecorr;
        EcepaL = Ecepa;

}// end else if wfn_type_ == "CEPA"

        mem_release();
        return EcepaL;

} // end of compute_energy


/*=======================*/
/*  ref_energy()         */
/*=======================*/
void OCEPAWave::ref_energy()
{
     double Ehf;     
     Ehf=0.0;

 if (reference_ == "RESTRICTED") {
    for (int h=0; h<nirrep_; h++){
      for (int i=0; i<occpiA[h];i++) {
	Ehf+=HmoA->get(h,i,i) + FockA->get(h,i,i);
      }
    }         
    Eref = Ehf + Enuc;
 }// end rhf
 
 else if (reference_ == "UNRESTRICTED") { 
     // alpha contribution
     for (int h=0; h<nirrep_; h++){
      for (int i=0; i<occpiA[h];i++) {
	Ehf+=HmoA->get(h,i,i) + FockA->get(h,i,i);
      }
    }  
    
    // beta contribution
     for (int h=0; h<nirrep_; h++){
      for (int i=0; i<occpiB[h];i++) {
	Ehf+=HmoB->get(h,i,i) + FockB->get(h,i,i);
      }
    }  
    
    Eref = (0.5 * Ehf) + Enuc; 
 }// end uhf
    
} // end of ref_energy


/*=======================*/
/*  mp2_energy()         */
/*=======================*/
void OCEPAWave::mp2_energy()
{
     dpdbuf4 K, T;
     
     psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
     psio_->open(PSIF_OCEPA_DPD, PSIO_OPEN_OLD);
     
     Ecorr = 0.0;

     Escsmp2AA = 0.0;
     Escsmp2AB = 0.0;
     Escsmp2BB = 0.0;
     Escsmp2 = 0.0;

     Esosmp2AB = 0.0;
     Esosmp2 = 0.0;

     Escsnmp2AA = 0.0;
     Escsnmp2BB = 0.0;
     Escsnmp2 = 0.0;
     
     Escsmimp2AA = 0.0;
     Escsmimp2AB = 0.0;
     Escsmimp2BB = 0.0;
     Escsmimp2 = 0.0;

     Escsmp2vdwAA = 0.0;
     Escsmp2vdwAB = 0.0;
     Escsmp2vdwBB = 0.0;
     Escsmp2vdw = 0.0;

     Esospimp2AB = 0.0;
     Esospimp2 = 0.0;
     
 if (reference_ == "RESTRICTED") {
     // Compute Energy
     // Alpha-Alpha spin contribution
     dpd_buf4_init(&T, PSIF_OCEPA_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2AA <OO|VV>");
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
     Emp2AA = 0.5 * dpd_buf4_dot(&T, &K);     
     dpd_buf4_close(&T);
     Emp2BB = Emp2AA;  
     
     Escsmp2AA = ss_scale * Emp2AA; 
     Escsnmp2AA = 1.76 * Emp2AA; 
     Escsmimp2AA = 1.29 * Emp2AA; 
     Escsmp2vdwAA = 0.50 * Emp2AA; 
     
     Escsmp2BB = ss_scale * Emp2BB;  
     Escsnmp2BB = 1.76 * Emp2BB; 
     Escsmimp2BB = 1.29 * Emp2BB; 
     Escsmp2vdwBB = 0.50 * Emp2BB; 
     
     // Alpha-Beta spin contribution
     dpd_buf4_init(&T, PSIF_OCEPA_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
     Emp2AB = dpd_buf4_dot(&T, &K);     
     dpd_buf4_close(&T);
     dpd_buf4_close(&K);
     
     Escsmp2AB = os_scale * Emp2AB;  
     Esosmp2AB = sos_scale * Emp2AB; 
     Escsmimp2AB = 0.40 * Emp2AB; 
     Escsmp2vdwAB = 1.28 * Emp2AB; 
     Esospimp2AB = 1.40 * Emp2AB; 
     Ecorr = Emp2AA + Emp2BB + Emp2AB;
 
 }// end rhf
 
 else if (reference_ == "UNRESTRICTED") { 
     // Compute Energy
     // Alpha-Alpha spin contribution
     dpd_buf4_init(&T, PSIF_OCEPA_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO||VV>");
     Ecorr += 0.25 * dpd_buf4_dot(&T, &K);     
     dpd_buf4_close(&T);
     dpd_buf4_close(&K);
     
     Emp2AA = Ecorr;    
     Escsmp2AA = ss_scale * Emp2AA; 
     Escsnmp2AA = 1.76 * Emp2AA; 
     Escsmimp2AA = 1.29 * Emp2AA; 
     Escsmp2vdwAA = 0.50 * Emp2AA; 
     
     
     // Alpha-Beta spin contribution
     dpd_buf4_init(&T, PSIF_OCEPA_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                 ID("[O,o]"), ID("[V,v]"), 0, "T2 <Oo|Vv>");
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
     Ecorr += dpd_buf4_dot(&T, &K);     
     dpd_buf4_close(&T);
     dpd_buf4_close(&K);
     
     Emp2AB = Ecorr - Emp2AA;
     Escsmp2AB = os_scale * Emp2AB;  
     Esosmp2AB = sos_scale * Emp2AB; 
     Escsmimp2AB = 0.40 * Emp2AB; 
     Escsmp2vdwAB = 1.28 * Emp2AB; 
     Esospimp2AB = 1.40 * Emp2AB; 
     
     
     // Beta-Beta spin contribution
     dpd_buf4_init(&T, PSIF_OCEPA_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2 <oo|vv>");
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo||vv>");
     Ecorr += 0.25 * dpd_buf4_dot(&T, &K);     
     dpd_buf4_close(&T);
     dpd_buf4_close(&K);
     
     Emp2BB = Ecorr - Emp2AA - Emp2AB;  
     Escsmp2BB = ss_scale * Emp2BB;  
     Escsnmp2BB = 1.76 * Emp2BB; 
     Escsmimp2BB = 1.29 * Emp2BB; 
     Escsmp2vdwBB = 0.50 * Emp2BB; 
     
 }// end uhf

     Emp2 = Eref + Ecorr;
     Escsmp2 = Eref + Escsmp2AA + Escsmp2AB + Escsmp2BB;
     Esosmp2 = Eref + Esosmp2AB;     
     Escsnmp2 = Eref + Escsnmp2AA + Escsnmp2BB;
     Escsmimp2 = Eref + Escsmimp2AA + Escsmimp2AB + Escsmimp2BB;
     Escsmp2vdw = Eref + Escsmp2vdwAA + Escsmp2vdwAB + Escsmp2vdwBB;
     Esospimp2 = Eref + Esospimp2AB;     
     
     psio_->close(PSIF_LIBTRANS_DPD, 1);
     psio_->close(PSIF_OCEPA_DPD, 1);    
     
} // end of mp2_energy


/*=======================*/
/*  cepa_energy()         */
/*=======================*/
void OCEPAWave::cepa_energy()
{
     dpdbuf4 K, T;
     
     psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
     psio_->open(PSIF_OCEPA_DPD, PSIO_OPEN_OLD);
     
     Ecorr = 0.0;
 
 if (reference_ == "RESTRICTED") {
     // Compute Energy
     // Alpha-Alpha spin contribution
     dpd_buf4_init(&T, PSIF_OCEPA_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2AA <OO|VV>");
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
     EcepaAA = 0.5 * dpd_buf4_dot(&T, &K);     
     dpd_buf4_close(&T);
     EcepaBB = EcepaAA;    
     
     // Alpha-Beta spin contribution
     dpd_buf4_init(&T, PSIF_OCEPA_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
     EcepaAB = dpd_buf4_dot(&T, &K);     
     dpd_buf4_close(&T);
     dpd_buf4_close(&K);     

 }// end rhf

 else if (reference_ == "UNRESTRICTED") {
     // Compute Energy
     // Alpha-Alpha spin contribution
     dpd_buf4_init(&T, PSIF_OCEPA_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO||VV>");
     EcepaAA = 0.25 * dpd_buf4_dot(&T, &K);     
     dpd_buf4_close(&T);
     dpd_buf4_close(&K);     
     
     
     // Alpha-Beta spin contribution
     dpd_buf4_init(&T, PSIF_OCEPA_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                 ID("[O,o]"), ID("[V,v]"), 0, "T2 <Oo|Vv>");
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
     EcepaAB = dpd_buf4_dot(&T, &K);     
     dpd_buf4_close(&T);
     dpd_buf4_close(&K);     
 
     // Beta-Beta spin contribution
     dpd_buf4_init(&T, PSIF_OCEPA_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2 <oo|vv>");
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo||vv>");
     EcepaBB = 0.25 * dpd_buf4_dot(&T, &K);     
     dpd_buf4_close(&T);
     dpd_buf4_close(&K);     

 }// end uhf
     
     Ecorr = EcepaAA + EcepaBB + EcepaAB;
     Ecepa = Eref + Ecorr;
     Escscepa = Eref + ((cepa_ss_scale_ * (EcepaAA + EcepaBB)) + (cepa_os_scale_ * EcepaAB));
     Esoscepa = Eref + (cepa_sos_scale_ * EcepaAB);
     
     psio_->close(PSIF_LIBTRANS_DPD, 1);
     psio_->close(PSIF_OCEPA_DPD, 1);    
     
} // end of cepa_energy


/*=======================*/
/*  cepal_energy()        */
/*=======================*/
void OCEPAWave::cepal_energy()
{      
     //fprintf(outfile,"\n cepal_energy is starting... \n"); fflush(outfile);
     // One-electron contribution
     /*
     double Etpdm_oooo, Etpdm_oovv, Etpdm_ovov, Etpdm_vovo;
     HG1A->zero();
     HG1B->zero();
     HG1A->gemm(false, false, 1.0, HmoA, g1symmA, 0.0);
     HG1B->gemm(false, false, 1.0, HmoB, g1symmB, 0.0);
     Ecepa_rdm = 0.0;
     Ecepa_rdm = HG1A->trace() + HG1B->trace() + Enuc;
     Eopdm = HG1A->trace() + HG1B->trace();
     fprintf(outfile,"\n OPDM Contribution               : %12.14f\n", Eopdm); fflush(outfile);
     */
     
    // Two-electron contribution
    dpdbuf4 G, K;
    
    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);  
    psio_->open(PSIF_OCEPA_DENSITY, PSIO_OPEN_OLD); 

 if (reference_ == "RESTRICTED") {
    // OOOO-Block contribution
    // E += 2*G_ijkl <ij|kl>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "MO Ints <OO|OO>");
    dpd_buf4_init(&G, PSIF_OCEPA_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "TPDM <OO|OO>");
    Ecepa_rdm += 2.0 * dpd_buf4_dot(&G, &K);         
    dpd_buf4_close(&K);
    dpd_buf4_close(&G);   
    

    // VVVV-Block contribution
    // E += 2*G_ABCD <AB|CD>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "MO Ints <VV|VV>");
    dpd_buf4_init(&G, PSIF_OCEPA_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "TPDM <VV|VV>");    
    Ecepa_rdm += 2.0 * dpd_buf4_dot(&G, &K);         
    dpd_buf4_close(&K);
    dpd_buf4_close(&G); 
    //double Evvvv = 2.0 * dpd_buf4_dot(&G, &K);         
    //Ecepa_rdm += Evvvv;
    //fprintf(outfile,"\n VVVV Contribution               : %12.14f\n", Evvvv); 
    //fflush(outfile);
    
    
    // OOVV-Block contribution
    // E += 8*G_IJAB <IJ|AB>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
    dpd_buf4_init(&G, PSIF_OCEPA_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
    Ecepa_rdm += 8.0 * dpd_buf4_dot(&G, &K);   
    dpd_buf4_close(&K);
    dpd_buf4_close(&G);    
    
    
    // OVOV-Block contribution
    // E += 4*G_IAJB <IA|JB>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV|OV>");
    dpd_buf4_init(&G, PSIF_OCEPA_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
    Ecepa_rdm += 4.0 * dpd_buf4_dot(&G, &K);     
    dpd_buf4_close(&K);
    dpd_buf4_close(&G);

 }// end rhf

 else if (reference_ == "UNRESTRICTED") {
    // OOOO-Block contribution
    // E += G_IJKL <IJ||KL>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "MO Ints <OO||OO>");
    dpd_buf4_init(&G, PSIF_OCEPA_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "TPDM <OO|OO>");
    Ecepa_rdm += dpd_buf4_dot(&G, &K);         
    dpd_buf4_close(&K);
    dpd_buf4_close(&G);   
    //double Eoooo = dpd_buf4_dot(&G, &K);         
    
    // E += G_ijkl <ij||kl>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "MO Ints <oo||oo>");
    dpd_buf4_init(&G, PSIF_OCEPA_DENSITY, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "TPDM <oo|oo>");
    Ecepa_rdm += dpd_buf4_dot(&G, &K);  
    dpd_buf4_close(&K);
    dpd_buf4_close(&G);  
    //Eoooo += dpd_buf4_dot(&G, &K);         
    
    // E += 4*G_IjKl <Ij||Kl>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "MO Ints <Oo|Oo>");
    dpd_buf4_init(&G, PSIF_OCEPA_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "TPDM <Oo|Oo>");
    Ecepa_rdm += 4.0 * dpd_buf4_dot(&G, &K); 
    dpd_buf4_close(&K);
    dpd_buf4_close(&G);
    //Eoooo += 4.0 * dpd_buf4_dot(&G, &K);         
    //Ecepa_rdm += Eoooo;
    //fprintf(outfile,"\n OOOO Contribution               : %12.14f\n", Eoooo); 
    //fflush(outfile);

     
    // VVVV-Block contribution
    // E += 2*G_ABCD <AB|CD>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "MO Ints <VV|VV>");
    dpd_buf4_init(&G, PSIF_OCEPA_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "TPDM <VV|VV>");    
    Ecepa_rdm += 2.0 * dpd_buf4_dot(&G, &K);         
    dpd_buf4_close(&K);
    dpd_buf4_close(&G); 
    //double Evvvv = 2.0 * dpd_buf4_dot(&G, &K);         

    // E += 2*G_abcd <ab|cd>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v,v]"), ID("[v,v]"), 0, "MO Ints <vv|vv>");
    dpd_buf4_init(&G, PSIF_OCEPA_DENSITY, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v,v]"), ID("[v,v]"), 0, "TPDM <vv|vv>");
    Ecepa_rdm += 2.0 * dpd_buf4_dot(&G, &K);  
    dpd_buf4_close(&K);
    dpd_buf4_close(&G);  
    //Evvvv += 2.0 * dpd_buf4_dot(&G, &K);         
    
    // E += 4*G_AbCd <Ab|Cd>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,v]"), ID("[V,v]"),
                  ID("[V,v]"), ID("[V,v]"), 0, "MO Ints <Vv|Vv>");
    dpd_buf4_init(&G, PSIF_OCEPA_DENSITY, 0, ID("[V,v]"), ID("[V,v]"),
                  ID("[V,v]"), ID("[V,v]"), 0, "TPDM <Vv|Vv>");
    Ecepa_rdm += 4.0 * dpd_buf4_dot(&G, &K); 
    dpd_buf4_close(&K);
    dpd_buf4_close(&G);
    //Ecepa_rdm += Evvvv;
    //Evvvv += 4.0 * dpd_buf4_dot(&G, &K);         
    //fprintf(outfile,"\n VVVV Contribution               : %12.14f\n", Evvvv); 
    //fflush(outfile);
 
    
    // OOVV-Block contribution
    // E += 2*G_IJAB <IJ||AB>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO||VV>");
    dpd_buf4_init(&G, PSIF_OCEPA_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
    Ecepa_rdm += 2.0 * dpd_buf4_dot(&G, &K);   
    dpd_buf4_close(&K);
    dpd_buf4_close(&G);    
    
    // E += 2*G_ijab <ij|ab>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo||vv>");
    dpd_buf4_init(&G, PSIF_OCEPA_DENSITY, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "TPDM <oo|vv>");
    Ecepa_rdm += 2.0 * dpd_buf4_dot(&G, &K);   
    dpd_buf4_close(&K);
    dpd_buf4_close(&G); 
    
    // E += 8*G_IjAb <Ij||Ab>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
    dpd_buf4_init(&G, PSIF_OCEPA_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "TPDM <Oo|Vv>");
    Ecepa_rdm += 8.0 * dpd_buf4_dot(&G, &K);   
    dpd_buf4_close(&K);
    dpd_buf4_close(&G); 

    
    // OVOV-Block contribution
    // E += 4*G_IAJB <IA||JB>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV||OV>");
    dpd_buf4_init(&G, PSIF_OCEPA_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
    Ecepa_rdm += 4.0 * dpd_buf4_dot(&G, &K);     
    dpd_buf4_close(&K);
    dpd_buf4_close(&G);
    
    // E += 4*G_iajb <ia||jb>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints <ov||ov>");
    dpd_buf4_init(&G, PSIF_OCEPA_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "TPDM <ov|ov>");
    Ecepa_rdm += 4.0 * dpd_buf4_dot(&G, &K);     
    dpd_buf4_close(&K);
    dpd_buf4_close(&G);
    
    // E += 4*G_IaJb <Ia||Jb>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "MO Ints <Ov|Ov>");
    dpd_buf4_init(&G, PSIF_OCEPA_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "TPDM <Ov|Ov>");
    Ecepa_rdm += 4.0 * dpd_buf4_dot(&G, &K);     
    dpd_buf4_close(&K);
    dpd_buf4_close(&G);
    
    // VOVO-Block contribution
    // E += 4*G_AiBj <Ai||Bj>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0,ID("[V,o]"), ID("[V,o]"),
                  ID("[V,o]"), ID("[V,o]"), 0, "MO Ints <Vo|Vo>");
    dpd_buf4_init(&G, PSIF_OCEPA_DENSITY, 0, ID("[V,o]"), ID("[V,o]"),
                  ID("[V,o]"), ID("[V,o]"), 0, "TPDM <Vo|Vo>"); 
    Ecepa_rdm += 4.0 * dpd_buf4_dot(&G, &K);     
    dpd_buf4_close(&K);
    dpd_buf4_close(&G);
    
    // OVVO-Block contribution
    // E += 8*G_IaBj <Ia||Bj>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0,ID("[O,v]"), ID("[V,o]"),
                  ID("[O,v]"), ID("[V,o]"), 0, "MO Ints <Ov|Vo>");
    dpd_buf4_init(&G, PSIF_OCEPA_DENSITY, 0, ID("[O,v]"), ID("[V,o]"),
                  ID("[O,v]"), ID("[V,o]"), 0, "TPDM <Ov|Vo>"); 
    Ecepa_rdm += 8.0 * dpd_buf4_dot(&G, &K);     
    dpd_buf4_close(&K);
    dpd_buf4_close(&G);
    
    /*
    Eopdm = HG1A->trace() + HG1B->trace() + Enuc;
    double Etpdm = Ecepa_rdm - Eopdm;
    fprintf(outfile,"\n TPDM Contribution               : %12.14f\n", Etpdm); 
    fflush(outfile);
    */

 }// end uhf

    psio_->close(PSIF_LIBTRANS_DPD, 1);  
    psio_->close(PSIF_OCEPA_DENSITY, 1);  
  
    EcorrL = Ecepa_rdm - Escf;        
    EcepaL = Ecepa_rdm;  
    DE = EcepaL - EcepaL_old;
    EcepaL_old = EcepaL;
    
    //fprintf(outfile,"\n CEPA Total Energy via pdms       : %12.14f\n", Ecepa_rdm); fflush(outfile);    
    //fprintf(outfile,"\n cepal_energy is done... \n"); fflush(outfile);
    
} // end of cepal_energy


/*=======================*/
/*  nbo()                */
/*=======================*/
void OCEPAWave::nbo()
{

fprintf(outfile,"\n  \n");
fprintf(outfile," ============================================================================== \n");
fprintf(outfile," ======================== NBO ANALYSIS ======================================== \n");
fprintf(outfile," ============================================================================== \n");
fprintf(outfile,"\n Diagonalizing one-particle response density matrix... \n");
fprintf(outfile,"\n");
fflush(outfile);

      SharedMatrix Udum = boost::shared_ptr<Matrix>(new Matrix("Udum", nirrep_, nmopi_, nmopi_));
      SharedVector diag = boost::shared_ptr<Vector>(new Vector("Natural orbital occupation numbers", nirrep_, nmopi_));

      // Diagonalizing Alpha-OPDM
      Udum->zero();

      //diag->zero();
      for(int h = 0; h < nirrep_; h++){
	  for(int i = 0; i < nmopi_[h]; i++){
	    diag->set(h,i,0.0);
	  }
	}
 
 if (reference_ == "RESTRICTED") {
      g1symmA->diagonalize(Udum, diag);
	
      //trace
      //sum=diag->trace();
      sum=0.0;
      for(int h = 0; h < nirrep_; h++){
	  for(int i = 0; i < nmopi_[h]; i++){
	    sum+=diag->get(h,i);
	  }
	}
      
      
      fprintf(outfile, "\n Trace of one-particle density matrix: %20.14f \n",  sum);
      fprintf(outfile,"\n");
      fflush(outfile);
 }// end rhf

 else if (reference_ == "UNRESTRICTED") {
      g1symmA->diagonalize(Udum, diag);
	
      //trace
      //sum=diag->trace();
      sum=0.0;
      for(int h = 0; h < nirrep_; h++){
	  for(int i = 0; i < nmopi_[h]; i++){
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
      for(int h = 0; h < nirrep_; h++){
	  for(int i = 0; i < nmopi_[h]; i++){
	    diag->set(h,i,0.0);
	  }
	}
   
      g1symmB->diagonalize(Udum, diag);
	
      //trace
      //sum=diag->trace();
      sum=0.0;
      for(int h = 0; h < nirrep_; h++){
	  for(int i = 0; i < nmopi_[h]; i++){
	    sum+=diag->get(h,i);
	  }
	}      
      
      fprintf(outfile, "\n Trace of beta one-particle density matrix: %20.14f \n",  sum);
      fprintf(outfile,"\n");
      fflush(outfile);

 }// end uhf

      //print
      diag->print();     
      
} // end of nbo


/*=======================*/
/*  mem_release()        */
/*=======================*/
void OCEPAWave::mem_release()
{   
	delete ints;
	delete [] idprowA;
	delete [] idpcolA;
	delete [] idpirrA;
	delete [] pitzer2symblk;
	delete [] pitzer2symirrep;
	delete [] PitzerOffset;
	delete [] sosym;
	delete [] mosym;
	delete [] occ_offA;
	delete [] vir_offA;
	delete [] occ2symblkA;
	delete [] virt2symblkA;
        delete wogA;
	delete kappaA;
	delete kappa_barA;

       if (reference_ == "UNRESTRICTED") {
	delete [] idprowB;
	delete [] idpcolB;
	delete [] idpirrB;
	delete [] occ_offB;
	delete [] vir_offB;
	delete [] occ2symblkB;
	delete [] virt2symblkB;
	delete wogB;
	delete kappaB;
	delete kappa_barB;
      }
	
	if (opt_method == "DIIS" && wfn_type_ == "OCEPA") {
          delete vecsA;
          delete errvecsA;
          if (reference_ == "UNRESTRICTED") delete vecsB;
          if (reference_ == "UNRESTRICTED") delete errvecsB;
	}
	
	chkpt_.reset();

      if (reference_ == "RESTRICTED") {
	Ca_.reset();
	Ca_ref.reset();
	Hso.reset();
	Tso.reset();
	Vso.reset();
	HmoA.reset();
	FockA.reset();
	gamma1corr.reset();
	g1symm.reset();
	GFock.reset();
	UorbA.reset();
	KorbA.reset();
	KsqrA.reset();
	HG1.reset();
	WorbA.reset();
	GooA.reset();
	GvvA.reset();
       }

       else if (reference_ == "UNRESTRICTED") {
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
       }

	//fprintf(outfile,"\n mem_release done. \n"); fflush(outfile);
}//
} }

