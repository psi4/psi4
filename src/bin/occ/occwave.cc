#include <fstream>
#include <math.h>

#include <libtrans/integraltransform.h>
#include <libtrans/mospace.h>

#include "occwave.h"

using namespace psi;
using namespace boost;

namespace psi { namespace occwave{

OCCWave::OCCWave(boost::shared_ptr<Wavefunction> reference_wavefunction, Options &options)
    : Wavefunction(options, _default_psio_lib_)
{
    reference_wavefunction_ = reference_wavefunction;
    common_init();
}//

OCCWave::~OCCWave()
{
}//


void OCCWave::common_init()
{

	tol_Eod=options_.get_double("E_CONVERGENCE");
	tol_t2=options_.get_double("R_CONVERGENCE");
	mograd_max=options_.get_double("MAX_MOGRAD_CONVERGENCE");

        cc_maxiter=options_.get_int("CC_MAXITER");
	mo_maxiter=options_.get_int("MO_MAXITER");
	print_=options_.get_int("PRINT"); 
	cachelev=options_.get_int("CACHELEVEL"); 
	exp_cutoff=options_.get_int("CUTOFF");
        tol_pcg=options_.get_double("PCG_CONVERGENCE");
        pcg_maxiter=options_.get_int("PCG_MAXITER");
	num_vecs=options_.get_int("MO_DIIS_NUM_VECS");
	cc_maxdiis_=options_.get_int("CC_DIIS_MAX_VECS");
	cc_mindiis_=options_.get_int("CC_DIIS_MIN_VECS");
        ep_maxiter=options_.get_int("EP_MAXITER");
	
	step_max=options_.get_double("MO_STEP_MAX");
	lshift_parameter=options_.get_double("LEVEL_SHIFT");
	os_scale=options_.get_double("MP2_OS_SCALE");
	ss_scale=options_.get_double("MP2_SS_SCALE");
	sos_scale=options_.get_double("MP2_SOS_SCALE");
	sos_scale2=options_.get_double("MP2_SOS_SCALE2");
	cepa_os_scale_=options_.get_double("CEPA_OS_SCALE");
	cepa_ss_scale_=options_.get_double("CEPA_SS_SCALE");
	cepa_sos_scale_=options_.get_double("CEPA_SOS_SCALE");
	e3_scale=options_.get_double("E3_SCALE");
	lambda_damping=options_.get_double("MOGRAD_DAMPING");
	
	wfn_type_=options_.get_str("WFN_TYPE");
	orth_type=options_.get_str("ORTH_TYPE");
	opt_method=options_.get_str("OPT_METHOD");
	//hess_type=options_.get_str("HESS_TYPE");
	occ_orb_energy=options_.get_str("OCC_ORBS_PRINT");
	natorb=options_.get_str("NAT_ORBS");
	reference=options_.get_str("REFERENCE");
	do_scs=options_.get_str("DO_SCS");
	do_sos=options_.get_str("DO_SOS");
	write_mo_coeff=options_.get_str("MO_WRITE");
	read_mo_coeff=options_.get_str("MO_READ");
        lineq=options_.get_str("LINEQ_SOLVER"); 
	level_shift=options_.get_str("DO_LEVEL_SHIFT");
	scs_type_=options_.get_str("SCS_TYPE");
	sos_type_=options_.get_str("SOS_TYPE");
	dertype=options_.get_str("DERTYPE");
	pcg_beta_type_=options_.get_str("PCG_BETA_TYPE");
        twopdm_abcd_type=options_.get_str("TPDM_ABCD_TYPE");
	dertype=options_.get_str("DERTYPE");
	pcg_beta_type_=options_.get_str("PCG_BETA_TYPE");
        compute_ccl=options_.get_str("CCL_ENERGY"); 
        orb_resp_solver_=options_.get_str("ORB_RESP_SOLVER"); 
        ip_poles=options_.get_str("IP_POLES"); 
        ea_poles=options_.get_str("EA_POLES"); 
        ep_ip_poles=options_.get_str("EP_IP_POLES"); 
        ep_ea_poles=options_.get_str("EP_EA_POLES"); 
        ekt_ip_=options_.get_str("EKT_IP"); 
        ekt_ea_=options_.get_str("EKT_EA"); 
        orb_opt_=options_.get_str("ORB_OPT"); 

    // Tying orbital convergence to the desired e_conv,
    //   particularly important for sane numerical frequencies by energy
    //   These have been determined by linear fits to a step fn 
    //   based on e_conv on limited numerical tests.
    // The printed value from options_.print() will not be accurate
    //   since newly set orbital conv is not written back to options
    if (options_["RMS_MOGRAD_CONVERGENCE"].has_changed()) {
        tol_grad=options_.get_double("RMS_MOGRAD_CONVERGENCE");
    } 
    else {
        double temp;
        temp = 2.0 - 0.5 * log10(tol_Eod); // I think (U.B) this is the desirable map balancing accuracy and efficiency.
        //temp = 3.0 - 0.5 * log10(tol_Eod); // Lori's old map leads unecessary iterations for the omp2-2 test case.
        //temp = 1.74 - 0.71 * log10(tol_Eod); //OLD map for wfn != OMP2
        if (temp < 5.0) {
            temp = 5.0;
        }
        tol_grad = pow(10.0, -temp);
    }

        if (reference == "RHF" || reference == "RKS") reference_ = "RESTRICTED";
        else if (reference == "UHF" || reference == "UKS" || reference == "ROHF") reference_ = "UNRESTRICTED";

        // Only UHF is allowed for the standard methods, except for MP2
        if (reference == "ROHF" && orb_opt_ == "FALSE" && wfn_type_ != "OMP2") {
           throw PSIEXCEPTION("The ROHF reference is not available for the standard methods (except for MP2)!");
        }

        // Only ROHF-MP2 energy is available, not the gradients
        else if (reference == "ROHF" && orb_opt_ == "FALSE" && dertype != "NONE") {
           throw PSIEXCEPTION("ROHF-MP2 analytic gradients are not available, UHF-MP2 is recommended.");
        }

        if (options_.get_str("DO_DIIS") == "TRUE") do_diis_ = 1;
        else if (options_.get_str("DO_DIIS") == "FALSE") do_diis_ = 0;

	cutoff = pow(10.0,-exp_cutoff);
	if (print_ > 0) options_.print();
        title();
        if (reference == "ROHF") reference_wavefunction_->semicanonicalize();
	get_moinfo();

	
if (reference_ == "RESTRICTED") {
	// Memory allocation
	HmoA = boost::shared_ptr<Matrix>(new Matrix("MO-basis alpha one-electron ints", nirrep_, nmopi_, nmopi_));
	FockA = boost::shared_ptr<Matrix>(new Matrix("MO-basis alpha Fock matrix", nirrep_, nmopi_, nmopi_));
	gamma1corr = boost::shared_ptr<Matrix>(new Matrix("MO-basis alpha correlation OPDM", nirrep_, nmopi_, nmopi_));
	g1symm = boost::shared_ptr<Matrix>(new Matrix("MO-basis alpha OPDM", nirrep_, nmopi_, nmopi_));
	GFock = boost::shared_ptr<Matrix>(new Matrix("MO-basis alpha generalized Fock matrix", nirrep_, nmopi_, nmopi_));
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

        // Compute costs
        //cost_iabc_ = 8 * nooA * nvoA * nvoA * nvoA;
        // compute cost_iabc and cost_abcd
        cost_iabc_ = 0;
        cost_abcd_ = 0;
        for(int h=0; h < nirrep_; h++) {
            cost_iabc_ += (ULI)ov_pairpiAA[h] * (ULI)vv_pairpiAA[h];
            cost_abcd_ += (ULI)vv_pairpiAA[h] * (ULI)vv_pairpiAA[h];
        }
        cost_iabc_ /= (ULI)1024 * (ULI)1024;
        cost_abcd_ /= (ULI)1024 * (ULI)1024;
        cost_iabc_ *= (ULI)sizeof(double);
        cost_abcd_ *= (ULI)sizeof(double);
       
        // print
    if (wfn_type_ == "OMP2") {
        // Print memory
        memory = Process::environment.get_memory();
        memory_mb_ = memory/1000000L;
        fprintf(outfile,"\n\tMemory is %6lu MB \n", memory_mb_); 
        fprintf(outfile,"\tCost of iabc is %6lu MB \n", cost_iabc_); 
        fprintf(outfile,"\tCost of abcd is %6lu MB \n", cost_abcd_); 
	fflush(outfile);
        if (cost_iabc_ < memory_mb_) { 
            incore_iabc_ = 1;
            fprintf(outfile,     "\tSwitching to the incore algoritm for iabc..\n"); 
	    fflush(outfile);
        }
        else { 
            incore_iabc_ = 0;
            fprintf(outfile,     "\tSwitching to the out of core algoritm for iabc..\n"); 
	    fflush(outfile);
        }

        //cost_abcd_ = 8 * nvoA * nvoA * nvoA * nvoA;
        if (cost_abcd_ < memory_mb_) { 
            incore_abcd_ = 1;
            fprintf(outfile,     "\tSwitching to the incore algoritm for abcd..\n"); 
	    fflush(outfile);
        }
        else { 
            incore_abcd_ = 0;
            fprintf(outfile,     "\tSwitching to the out of core algoritm for abcd..\n"); 
	    fflush(outfile);
        }
    }// end if (wfn_type_ == "OMP2") 


    // Alloc ints
    std::vector<boost::shared_ptr<MOSpace> > spaces;
    spaces.push_back(MOSpace::occ);
    spaces.push_back(MOSpace::vir);

if (wfn_type_ == "OMP2" && incore_iabc_ == 0) { 
    ints = new IntegralTransform(reference_wavefunction_, spaces, 
                           IntegralTransform::Restricted,
                           IntegralTransform::IWLAndDPD,
                           IntegralTransform::QTOrder,
                           IntegralTransform::OccOnly,
                           false);
}

else {
    ints = new IntegralTransform(reference_wavefunction_, spaces, 
                           IntegralTransform::Restricted,
                           IntegralTransform::DPDOnly,
                           IntegralTransform::QTOrder,
                           IntegralTransform::OccOnly,
                           false);
}
                           
                          
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
	HmoB = boost::shared_ptr<Matrix>(new Matrix("MO-basis beta one-electron ints", nirrep_, nmopi_, nmopi_));
	FockA = boost::shared_ptr<Matrix>(new Matrix("MO-basis alpha Fock matrix", nirrep_, nmopi_, nmopi_));
	FockB = boost::shared_ptr<Matrix>(new Matrix("MO-basis beta Fock matrix", nirrep_, nmopi_, nmopi_));
	gamma1corrA = boost::shared_ptr<Matrix>(new Matrix("MO-basis alpha correlation OPDM", nirrep_, nmopi_, nmopi_));
	gamma1corrB = boost::shared_ptr<Matrix>(new Matrix("MO-basis beta correlation OPDM", nirrep_, nmopi_, nmopi_));
	g1symmA = boost::shared_ptr<Matrix>(new Matrix("MO-basis alpha OPDM", nirrep_, nmopi_, nmopi_));
	g1symmB = boost::shared_ptr<Matrix>(new Matrix("MO-basis beta OPDM", nirrep_, nmopi_, nmopi_));
	GFockA = boost::shared_ptr<Matrix>(new Matrix("MO-basis alpha generalized Fock matrix", nirrep_, nmopi_, nmopi_));
	GFockB = boost::shared_ptr<Matrix>(new Matrix("MO-basis beta generalized Fock matrix", nirrep_, nmopi_, nmopi_));
	UorbA = boost::shared_ptr<Matrix>(new Matrix("Alpha MO rotation matrix", nirrep_, nmopi_, nmopi_));
	UorbB = boost::shared_ptr<Matrix>(new Matrix("Beta MO rotation matrix", nirrep_, nmopi_, nmopi_));
	KorbA = boost::shared_ptr<Matrix>(new Matrix("K alpha MO rotation", nirrep_, nmopi_, nmopi_)); 
	KorbB = boost::shared_ptr<Matrix>(new Matrix("K beta MO rotation", nirrep_, nmopi_, nmopi_)); 
	KsqrA = boost::shared_ptr<Matrix>(new Matrix("K^2 alpha MO rotation", nirrep_, nmopi_, nmopi_)); 
	KsqrB = boost::shared_ptr<Matrix>(new Matrix("K^2 beta MO rotation", nirrep_, nmopi_, nmopi_)); 
	HG1A = boost::shared_ptr<Matrix>(new Matrix("Alpha h*g1symm", nirrep_, nmopi_, nmopi_));
	HG1B = boost::shared_ptr<Matrix>(new Matrix("Beta h*g1symm", nirrep_, nmopi_, nmopi_));
	WorbA = boost::shared_ptr<Matrix>(new Matrix("Alpha MO gradient matrix", nirrep_, nmopi_, nmopi_));
	WorbB = boost::shared_ptr<Matrix>(new Matrix("Beta MO gradient matrix", nirrep_, nmopi_, nmopi_));
	GooA = boost::shared_ptr<Matrix>(new Matrix("Alpha Goo intermediate", nirrep_, aoccpiA, aoccpiA));
	GooB = boost::shared_ptr<Matrix>(new Matrix("Beta Goo intermediate", nirrep_, aoccpiB, aoccpiB));
	GvvA = boost::shared_ptr<Matrix>(new Matrix("Alpha Gvv intermediate", nirrep_, avirtpiA, avirtpiA));
	GvvB = boost::shared_ptr<Matrix>(new Matrix("Beta Gvv intermediate", nirrep_, avirtpiB, avirtpiB));

        // ROHF-MP2
        if (reference == "ROHF" && orb_opt_ == "FALSE" && wfn_type_ == "OMP2") {
	    t1A = boost::shared_ptr<Matrix>(new Matrix("t_I^A", nirrep_, aoccpiA, avirtpiA)); 
	    t1B = boost::shared_ptr<Matrix>(new Matrix("t_i^a", nirrep_, aoccpiB, avirtpiB)); 
        }

        Molecule& mol = *reference_wavefunction_->molecule().get();
        CharacterTable ct = mol.point_group()->char_table();
        fprintf(outfile,"\tMO spaces per irreps... \n\n"); fflush(outfile);
        fprintf(outfile, "\tIRREP   FC   AOCC  BOCC  AVIR    BVIR  FV \n");
        fprintf(outfile, "\t==========================================\n");                                                 
        for(int h = 0; h < nirrep_; ++h){
         fprintf(outfile, "\t %3s   %3d   %3d   %3d   %3d    %3d   %3d\n",
                             ct.gamma(h).symbol(), frzcpi_[h], aoccpiA[h], aoccpiB[h], avirtpiA[h], avirtpiB[h], frzvpi_[h]);
        }
        fprintf(outfile,     "\t==========================================\n");
	fflush(outfile);

    // Alloc ints
    std::vector<boost::shared_ptr<MOSpace> > spaces;
    spaces.push_back(MOSpace::occ);
    spaces.push_back(MOSpace::vir);

    ints = new IntegralTransform(reference_wavefunction_, spaces, 
                           IntegralTransform::Unrestricted,
                           IntegralTransform::DPDOnly,
                           IntegralTransform::QTOrder,
                           IntegralTransform::OccOnly,
                           false);
                           
                          
    ints->set_print(0);
    ints->set_dpd_id(0);
    ints->set_keep_iwl_so_ints(true);
    ints->set_keep_dpd_so_ints(true);           
    ints->initialize();
    dpd_set_default(ints->get_dpd_id());

}// end if (reference_ == "UNRESTRICTED") 
}// end common_init

void OCCWave::title()
{
   fprintf(outfile,"\n");
   fprintf(outfile," ============================================================================== \n");
   fprintf(outfile," ============================================================================== \n");
   fprintf(outfile," ============================================================================== \n");
   fprintf(outfile,"\n");
   if (wfn_type_ == "OMP2" && orb_opt_ == "TRUE") fprintf(outfile,"                       OMP2 (OO-MP2)   \n");
   else if (wfn_type_ == "OMP2" && orb_opt_ == "FALSE") fprintf(outfile,"                       MP2   \n");
   else if (wfn_type_ == "OMP3" && orb_opt_ == "TRUE") fprintf(outfile,"                       OMP3 (OO-MP3)   \n");
   else if (wfn_type_ == "OMP3" && orb_opt_ == "FALSE") fprintf(outfile,"                       MP3   \n");
   else if (wfn_type_ == "OCEPA" && orb_opt_ == "TRUE") fprintf(outfile,"                       OCEPA (OO-CEPA)   \n");
   else if (wfn_type_ == "OCEPA" && orb_opt_ == "FALSE") fprintf(outfile,"                       CEPA   \n");
   else if (wfn_type_ == "OMP2.5" && orb_opt_ == "TRUE") fprintf(outfile,"                       OMP2.5 (OO-MP2.5)   \n");
   else if (wfn_type_ == "OMP2.5" && orb_opt_ == "FALSE") fprintf(outfile,"                       MP2.5  \n");
   fprintf(outfile,"              Program Written by Ugur Bozkaya,\n") ; 
   fprintf(outfile,"              Latest Revision April 4, 2013.\n") ;
   fprintf(outfile,"\n");
   fprintf(outfile," ============================================================================== \n");
   fprintf(outfile," ============================================================================== \n");
   fprintf(outfile," ============================================================================== \n");
   fprintf(outfile,"\n");
   fflush(outfile);
}//


double OCCWave::compute_energy()
{   
        
	// Warnings 
	if (nfrzc != 0 && orb_opt_ == "TRUE") {
          throw FeatureNotImplemented("Orbital-optimized methods", "Frozen core/virtual", __FILE__, __LINE__);
	}

	else if (nfrzv != 0 && orb_opt_ == "TRUE") {
          throw FeatureNotImplemented("Orbital-optimized methods", "Frozen core/virtual", __FILE__, __LINE__);
	}

	else if (nfrzv != 0 && orb_opt_ == "FALSE") {
          throw FeatureNotImplemented("OCC module standard methods", "Frozen virtual", __FILE__, __LINE__);
	}

	else if (nfrzc != 0 && dertype != "NONE") {
          throw FeatureNotImplemented("OCC module analytic gradients", "Frozen core/virtual", __FILE__, __LINE__);
	}

        // Call the appropriate manager
        if (wfn_type_ == "OMP2" && orb_opt_ == "TRUE") omp2_manager();
        else if (wfn_type_ == "OMP2" && orb_opt_ == "FALSE") mp2_manager();
        else if (wfn_type_ == "OMP3" && orb_opt_ == "TRUE") omp3_manager();
        else if (wfn_type_ == "OMP3" && orb_opt_ == "FALSE") mp3_manager();
        else if (wfn_type_ == "OCEPA" && orb_opt_ == "TRUE") ocepa_manager();
        else if (wfn_type_ == "OCEPA" && orb_opt_ == "FALSE") cepa_manager();
        else if (wfn_type_ == "OMP2.5" && orb_opt_ == "TRUE") omp2_5_manager();
        else if (wfn_type_ == "OMP2.5" && orb_opt_ == "FALSE") mp2_5_manager();
	
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

        // release the memory
        mem_release();

        if (wfn_type_ == "OMP2") return Emp2L;
        else if (wfn_type_ == "OMP3" || wfn_type_ == "OMP2.5") return Emp3L;
        else if (wfn_type_ == "OCEPA") return EcepaL;
        else if (wfn_type_ == "CEPA") return Ecepa;
} // end of compute_energy


void OCCWave::nbo()
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
      g1symm->diagonalize(Udum, diag);
	
      //trace
      //sum=diag->trace();
      sum=0.0;
      for(int h = 0; h < nirrep_; h++){
	  for(int i = 0; i < nmopi_[h]; i++){
	    sum+=diag->get(h,i);
	  }
	}
      
      fprintf(outfile, "\n Trace of one-particle density matrix: %20.14f \n\n",  sum);
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
      
      fprintf(outfile, "\n Trace of alpha one-particle density matrix: %20.14f \n\n",  sum);
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

void OCCWave::mem_release()
{   
	chkpt_.reset();
	delete ints;
	delete [] pitzer2symblk;
	delete [] pitzer2symirrep;
	delete [] PitzerOffset;
	delete [] sosym;
	delete [] mosym;
	delete [] occ_offA;
	delete [] vir_offA;
	delete [] occ2symblkA;
	delete [] virt2symblkA;
        delete [] pitzer2qtA;
        delete [] qt2pitzerA;

      if (reference_ == "RESTRICTED") {
        delete [] oo_pairpiAA;
        delete [] ov_pairpiAA;
        delete [] vv_pairpiAA;
        delete oo_pairidxAA;
        delete vv_pairidxAA;

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
	delete [] occ_offB;
	delete [] vir_offB;
	delete [] occ2symblkB;
	delete [] virt2symblkB;
        delete [] pitzer2qtB;
        delete [] qt2pitzerB;

        if (reference == "ROHF" && orb_opt_ == "FALSE" && wfn_type_ != "OMP2") {
 	    t1A.reset();
	    t1B.reset();
        }

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

