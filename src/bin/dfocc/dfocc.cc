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

#include <fstream>
#include <math.h>
#include <libqt/qt.h>
#include "dfocc.h"
#include "defines.h"

using namespace psi;
using namespace boost;

namespace psi { namespace dfoccwave {

DFOCC::DFOCC(boost::shared_ptr<Wavefunction> reference_wavefunction, Options &options)
    : Wavefunction(options, _default_psio_lib_)
{
    reference_wavefunction_ = reference_wavefunction;
    common_init();
}//

DFOCC::~DFOCC()
{
}//


void DFOCC::common_init()
{

    print_=options_.get_int("PRINT"); 
    if (print_ > 0) options_.print();

    cc_maxiter=options_.get_int("CC_MAXITER");
    mo_maxiter=options_.get_int("MO_MAXITER");
    num_vecs=options_.get_int("MO_DIIS_NUM_VECS");
    cc_maxdiis_=options_.get_int("CC_DIIS_MAX_VECS");
    cc_mindiis_=options_.get_int("CC_DIIS_MIN_VECS");
    exp_cutoff=options_.get_int("CUTOFF");
    pcg_maxiter=options_.get_int("PCG_MAXITER");

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
    tol_Eod=options_.get_double("E_CONVERGENCE");
    tol_t2=options_.get_double("R_CONVERGENCE");
    tol_pcg=options_.get_double("PCG_CONVERGENCE");
    reg_param=options_.get_double("REG_PARAM");

    orth_type=options_.get_str("ORTH_TYPE");
    opt_method=options_.get_str("OPT_METHOD");
    occ_orb_energy=options_.get_str("OCC_ORBS_PRINT");
    natorb=options_.get_str("NAT_ORBS");
    do_scs=options_.get_str("DO_SCS");
    do_sos=options_.get_str("DO_SOS");
    lineq=options_.get_str("LINEQ_SOLVER");
    level_shift=options_.get_str("DO_LEVEL_SHIFT");
    scs_type_=options_.get_str("SCS_TYPE");
    sos_type_=options_.get_str("SOS_TYPE");
    dertype=options_.get_str("DERTYPE");
    orb_resp_solver_=options_.get_str("ORB_RESP_SOLVER");
    ekt_ip_=options_.get_str("EKT_IP");
    ekt_ea_=options_.get_str("EKT_EA");
    reference=options_.get_str("REFERENCE");
    wfn_type_=options_.get_str("WFN_TYPE");
    orb_opt_=options_.get_str("ORB_OPT");
    conv_tei_type=options_.get_str("CONV_TEI_TYPE");
    pcg_beta_type_=options_.get_str("PCG_BETA_TYPE");
    regularization=options_.get_str("REGULARIZATION");
    read_scf_3index=options_.get_str("READ_SCF_3INDEX");
    freeze_core_=options_.get_str("FREEZE_CORE");
    oeprop_=options_.get_str("OEPROP");
    comput_s2_=options_.get_str("COMPUT_S2");

    //title
    title();

    //   Tying orbital convergence to the desired e_conv,
    //   particularly important for the same numerical frequencies by energy
    //   These have been determined by linear fits to a step fn
    //   based on e_conv on limited numerical tests.
    //   The printed value from options_.print() will not be accurate
    //   since newly set orbital conv is not written back to options
    if (options_["RMS_MOGRAD_CONVERGENCE"].has_changed()) {
        tol_grad=options_.get_double("RMS_MOGRAD_CONVERGENCE");
    }
    else {
        double temp;
        temp = (-0.9 * log10(tol_Eod)) - 1.6; 
        if (temp < 4.0) {
            temp = 4.0;
        }
        tol_grad = pow(10.0, -temp);
        //tol_grad = 100.0*tol_Eod; 
        fprintf(outfile,"\tRMS orbital gradient is changed to : %12.2e\n", tol_grad);
        fflush(outfile);
    }

    // Determine the MAXIMUM MOGRAD CONVERGENCE
    if (options_["MAX_MOGRAD_CONVERGENCE"].has_changed()) {
    mograd_max=options_.get_double("MAX_MOGRAD_CONVERGENCE");
    }
    else {
        double temp2;
        temp2 = (-0.8 * log10(tol_grad)) - 0.5;
        if (temp2 < 3.0) {
            temp2 = 3.0;
        }
        mograd_max = pow(10.0, -temp2);
        //mograd_max = 10.0*tol_grad;
        fprintf(outfile,"\tMAX orbital gradient is changed to : %12.2e\n", mograd_max);
        fflush(outfile);
    }

    // Figure out REF
    if (reference == "RHF" || reference == "RKS") reference_ = "RESTRICTED";
    else if (reference == "UHF" || reference == "UKS" || reference == "ROHF") reference_ = "UNRESTRICTED";
    if (reference == "ROHF") reference_wavefunction_->semicanonicalize();

    // Only ROHF-MP2 energy is available, not the gradients
    if (reference == "ROHF" && orb_opt_ == "FALSE" && dertype == "FIRST") {
             throw PSIEXCEPTION("ROHF DF-MP2 analytic gradients are not available, UHF DF-MP2 is recommended.");
    }

    // Frozen Core
    if (freeze_core_ == "TRUE" && orb_opt_ == "FALSE" && dertype == "FIRST") {
             throw PSIEXCEPTION("Frozen core gradients are available only for orbital-optimized methods.");
    }

    // DIIS
    if (options_.get_str("DO_DIIS") == "TRUE") do_diis_ = 1;
    else if (options_.get_str("DO_DIIS") == "FALSE") do_diis_ = 0;

    // Figure out HESSIAN TYPE
    if (options_["HESS_TYPE"].has_changed()) {
        hess_type=options_.get_str("HESS_TYPE");
    }
    else {
         if (reference_ == "RESTRICTED" && freeze_core_ == "FALSE") {
             hess_type = "HF"; 
         }
         else if (reference_ == "RESTRICTED" && freeze_core_ == "TRUE") {
             hess_type = "APPROX_DIAG";
         }
         else if (reference_ == "UNRESTRICTED") {
             hess_type = "APPROX_DIAG";
             fprintf(outfile,"\tMO Hessian type is changed to 'APPROX_DIAG'\n");
             fflush(outfile);
         }
    }

    // Regularization 
    if (regularization == "TRUE") {
        fprintf(outfile,"\n\tNOTE: A regularization procedure will be applied to the method.\n");
        fprintf(outfile,"\tThe regularization parameter is : %12.2f mh\n", reg_param * 1000.0);
        fflush(outfile);
    }

    cutoff = pow(10.0,-exp_cutoff);
    get_moinfo();
    pair_index();


if (reference_ == "RESTRICTED") {
	// Memory allocation
	HmoA = SharedTensor2d(new Tensor2d("MO-basis alpha one-electron ints", nmo_, nmo_));
        FijA = SharedTensor2d(new Tensor2d("Fint <I|J>", naoccA, naoccA));
        FabA = SharedTensor2d(new Tensor2d("Fint <A|B>", navirA, navirA));
        HooA = SharedTensor2d(new Tensor2d("OEI <O|O>", noccA, noccA));
        HovA = SharedTensor2d(new Tensor2d("OEI <O|V>", noccA, nvirA));
        HvoA = SharedTensor2d(new Tensor2d("OEI <V|O>", nvirA, noccA));
        HvvA = SharedTensor2d(new Tensor2d("OEI <V|V>", nvirA, nvirA));

    // if we need PDMs
    if (orb_opt_ == "TRUE" || dertype != "NONE" || oeprop_ == "TRUE") {
        GijA = SharedTensor2d(new Tensor2d("G Intermediate <I|J>", naoccA, naoccA));
        GabA = SharedTensor2d(new Tensor2d("G Intermediate <A|B>", navirA, navirA));
        G1c_oo = SharedTensor2d(new Tensor2d("Correlation OPDM <O|O>", noccA, noccA));
        G1c_vv = SharedTensor2d(new Tensor2d("Correlation OPDM <V|V>", nvirA, nvirA));
        G1 = SharedTensor2d(new Tensor2d("MO-basis OPDM", nmo_, nmo_));
        G1ao = SharedTensor2d(new Tensor2d("AO-basis OPDM", nso_, nso_));
        G1c = SharedTensor2d(new Tensor2d("MO-basis correlation OPDM", nmo_, nmo_));
        GF = SharedTensor2d(new Tensor2d("MO-basis GFM", nmo_, nmo_));
        GFao = SharedTensor2d(new Tensor2d("AO-basis GFM", nso_, nso_));
        GFoo = SharedTensor2d(new Tensor2d("MO-basis GFM <O|O>", noccA, noccA));
        GFvo = SharedTensor2d(new Tensor2d("MO-basis GFM <V|O>", nvirA, noccA));
        GFov = SharedTensor2d(new Tensor2d("MO-basis GFM <O|V>", noccA, nvirA));
        GFvv = SharedTensor2d(new Tensor2d("MO-basis GFM <V|V>", nvirA, nvirA));
        GFtvv = SharedTensor2d(new Tensor2d("MO-basis Complementary GFM <V|V>", nvirA, nvirA));
        WorbA = SharedTensor2d(new Tensor2d("MO-basis alpha MO gradient", nmo_, nmo_));
        UorbA = SharedTensor2d(new Tensor2d("Alpha MO rotation matrix", nmo_, nmo_));
        KorbA = SharedTensor2d(new Tensor2d("Alpha K MO rotation parameters matrix", nmo_, nmo_));
        KsqrA = SharedTensor2d(new Tensor2d("Alpha K^2 MO rotation parameters matrix", nmo_, nmo_));
        AvoA = SharedTensor2d(new Tensor2d("Diagonal MO Hessian <V|O>", nvirA, noccA));
        if (nfrzc > 0) AooA = SharedTensor2d(new Tensor2d("Diagonal MO Hessian <I|FC>", naoccA, nfrzc));
    }

        fprintf(outfile,"\n\tMO spaces... \n\n"); fflush(outfile);
        fprintf(outfile, "\t FC   OCC   VIR   FV \n");
        fprintf(outfile, "\t----------------------\n");                                                 
        fprintf(outfile, "\t%3d  %3d   %3d  %3d\n", nfrzc, naoccA, navirA, nfrzv);
	fflush(outfile);

        // memory requirements
        cost_ampAA = 0;
        cost_ampAA = (ULI)nocc2AA * (ULI)nvir2AA;
        cost_ampAA /= (ULI)1024 * (ULI)1024;
        cost_ampAA *= (ULI)sizeof(double);
        cost_amp = (ULI)3.0 * cost_ampAA;
        memory = Process::environment.get_memory();
        memory_mb = memory/1000000L;
        fprintf(outfile,"\n\tAvailable memory is: %6lu MB \n", memory_mb);
        fprintf(outfile,"\tMinimum required memory for the DFOCC module is: %6lu MB \n", cost_amp);
        fflush(outfile);

}  // end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
	// Memory allocation
	HmoA = SharedTensor2d(new Tensor2d("MO-basis alpha one-electron ints", nmo_, nmo_));
	HmoB = SharedTensor2d(new Tensor2d("MO-basis beta one-electron ints", nmo_, nmo_));
        FijA = SharedTensor2d(new Tensor2d("Fint <I|J>", naoccA, naoccA));
        FijB = SharedTensor2d(new Tensor2d("Fint <i|j>", naoccB, naoccB));
        FabA = SharedTensor2d(new Tensor2d("Fint <A|B>", navirA, navirA));
        FabB = SharedTensor2d(new Tensor2d("Fint <a|b>", navirB, navirB));
        HooA = SharedTensor2d(new Tensor2d("OEI <O|O>", noccA, noccA));
        HooB = SharedTensor2d(new Tensor2d("OEI <o|o>", noccB, noccB));
        HovA = SharedTensor2d(new Tensor2d("OEI <O|V>", noccA, nvirA));
        HovB = SharedTensor2d(new Tensor2d("OEI <o|v>", noccB, nvirB));
        HvoA = SharedTensor2d(new Tensor2d("OEI <V|O>", nvirA, noccA));
        HvoB = SharedTensor2d(new Tensor2d("OEI <v|o>", nvirB, noccB));
        HvvA = SharedTensor2d(new Tensor2d("OEI <V|V>", nvirA, nvirA));
        HvvB = SharedTensor2d(new Tensor2d("OEI <v|v>", nvirB, nvirB));

    // if we need PDMs
    if (orb_opt_ == "TRUE" || dertype != "NONE" || oeprop_ == "TRUE") {
        GijA = SharedTensor2d(new Tensor2d("G Intermediate <I|J>", naoccA, naoccA));
        GijB = SharedTensor2d(new Tensor2d("G Intermediate <i|j>", naoccB, naoccB));
        GabA = SharedTensor2d(new Tensor2d("G Intermediate <A|B>", navirA, navirA));
        GabB = SharedTensor2d(new Tensor2d("G Intermediate <a|b>", navirB, navirB));
        G1c_ooA = SharedTensor2d(new Tensor2d("Correlation OPDM <O|O>", noccA, noccA));
        G1c_ooB = SharedTensor2d(new Tensor2d("Correlation OPDM <o|o>", noccB, noccB));
        G1c_vvA = SharedTensor2d(new Tensor2d("Correlation OPDM <V|V>", nvirA, nvirA));
        G1c_vvB = SharedTensor2d(new Tensor2d("Correlation OPDM <v|v>", nvirB, nvirB));
        G1cA = SharedTensor2d(new Tensor2d("MO-basis alpha correlation OPDM", nmo_, nmo_));
        G1cB = SharedTensor2d(new Tensor2d("MO-basis beta correlation OPDM", nmo_, nmo_));
        G1A = SharedTensor2d(new Tensor2d("MO-basis alpha OPDM", nmo_, nmo_));
        G1B = SharedTensor2d(new Tensor2d("MO-basis beta OPDM", nmo_, nmo_));
        G1ao = SharedTensor2d(new Tensor2d("AO-basis OPDM", nso_, nso_));
        GFA = SharedTensor2d(new Tensor2d("MO-basis alpha GFM", nmo_, nmo_));
        GFB = SharedTensor2d(new Tensor2d("MO-basis beta GFM", nmo_, nmo_));
        GFao = SharedTensor2d(new Tensor2d("AO-basis GFM", nso_, nso_));
        GFooA = SharedTensor2d(new Tensor2d("MO-basis GFM <O|O>", noccA, noccA));
        GFooB = SharedTensor2d(new Tensor2d("MO-basis GFM <o|o>", noccB, noccB));
        GFvoA = SharedTensor2d(new Tensor2d("MO-basis GFM <V|O>", nvirA, noccA));
        GFvoB = SharedTensor2d(new Tensor2d("MO-basis GFM <v|o>", nvirB, noccB));
        GFovA = SharedTensor2d(new Tensor2d("MO-basis GFM <O|V>", noccA, nvirA));
        GFovB = SharedTensor2d(new Tensor2d("MO-basis GFM <o|v>", noccB, nvirB));
        GFvvA = SharedTensor2d(new Tensor2d("MO-basis GFM <V|V>", nvirA, nvirA));
        GFvvB = SharedTensor2d(new Tensor2d("MO-basis GFM <v|v>", nvirB, nvirB));
        GFtvvA = SharedTensor2d(new Tensor2d("MO-basis Complementary GFM <V|V>", nvirA, nvirA));
        GFtvvB = SharedTensor2d(new Tensor2d("MO-basis Complementary GFM <v|v>", nvirB, nvirB));
        WorbA = SharedTensor2d(new Tensor2d("MO-basis alpha MO gradient", nmo_, nmo_));
        WorbB = SharedTensor2d(new Tensor2d("MO-basis beta MO gradient", nmo_, nmo_));
        UorbA = SharedTensor2d(new Tensor2d("Alpha MO rotation matrix", nmo_, nmo_));
        UorbB = SharedTensor2d(new Tensor2d("Beta MO rotation matrix", nmo_, nmo_));
        KorbA = SharedTensor2d(new Tensor2d("Alpha K MO rotation parameters matrix", nmo_, nmo_));
        KorbB = SharedTensor2d(new Tensor2d("Beta K MO rotation parameters matrix", nmo_, nmo_));
        KsqrA = SharedTensor2d(new Tensor2d("Alpha K^2 MO rotation parameters matrix", nmo_, nmo_));
        KsqrB = SharedTensor2d(new Tensor2d("Beta K^2 MO rotation parameters matrix", nmo_, nmo_));
        AvoA = SharedTensor2d(new Tensor2d("Diagonal MO Hessian <V|O>", nvirA, noccA));
        AvoB = SharedTensor2d(new Tensor2d("Diagonal MO Hessian <v|o>", nvirB, noccB));
        if (nfrzc > 0) {
            AooA = SharedTensor2d(new Tensor2d("Diagonal MO Hessian <I|FC>", naoccA, nfrzc));
            AooB = SharedTensor2d(new Tensor2d("Diagonal MO Hessian <i|FC>", naoccA, nfrzc));
        }
    }

        // ROHF-MP2
        if (reference == "ROHF" && wfn_type_ == "DF-OMP2") {
            t1A = SharedTensor2d(new Tensor2d("T1_1 <I|A>", naoccA, navirA));
            t1B = SharedTensor2d(new Tensor2d("T1_1 <i|a>", naoccB, navirB));
            GiaA = SharedTensor2d(new Tensor2d("G Intermediate <I|A>", naoccA, navirA));
            GiaB = SharedTensor2d(new Tensor2d("G Intermediate <i|a>", naoccB, navirB));
            GaiA = SharedTensor2d(new Tensor2d("G Intermediate <A|I>", navirA, naoccA));
            GaiB = SharedTensor2d(new Tensor2d("G Intermediate <a|i>", navirB, naoccB));
            G1c_ovA = SharedTensor2d(new Tensor2d("Correlation OPDM <O|V>", noccA, nvirA));
            G1c_ovB = SharedTensor2d(new Tensor2d("Correlation OPDM <o|v>", noccB, nvirB));
            G1c_voA = SharedTensor2d(new Tensor2d("Correlation OPDM <V|O>", nvirA, noccA));
            G1c_voB = SharedTensor2d(new Tensor2d("Correlation OPDM <v|o>", nvirB, noccB));
        }

        fprintf(outfile,"\n\tMO spaces... \n\n"); fflush(outfile);
        fprintf(outfile, "\t FC   AOCC   BOCC  AVIR   BVIR   FV \n");
        fprintf(outfile, "\t------------------------------------------\n");
        fprintf(outfile, "\t%3d   %3d   %3d   %3d    %3d   %3d\n", nfrzc, naoccA, naoccB, navirA, navirB, nfrzv);
        fflush(outfile);

        // memory requirements
        cost_ampAA = 0;
        cost_ampAA = (ULI)nocc2AA * (ULI)nvir2AA;
        cost_ampAA /= (ULI)1024 * (ULI)1024;
        cost_ampAA *= (ULI)sizeof(double);
        cost_ampBB = (ULI)nocc2BB * (ULI)nvir2BB;
        cost_ampBB /= (ULI)1024 * (ULI)1024;
        cost_ampBB *= (ULI)sizeof(double);
        cost_ampAB = (ULI)nocc2AB * (ULI)nvir2AB;
        cost_ampAB /= (ULI)1024 * (ULI)1024;
        cost_ampAB *= (ULI)sizeof(double);
        cost_amp = MAX0(cost_ampAA, cost_ampBB);
        cost_amp = MAX0(cost_amp, cost_ampAB);
        cost_amp = (ULI)3.0 * cost_amp;
        memory = Process::environment.get_memory();
        memory_mb = memory/1000000L;
        fprintf(outfile,"\n\tAvailable memory is: %6lu MB \n", memory_mb);
        fprintf(outfile,"\tMinimum required memory for the DFOCC module is: %6lu MB \n", cost_amp);
        fflush(outfile);
}// else if (reference_ == "UNRESTRICTED")
	
        //fprintf(outfile,"\tI am here.\n"); fflush(outfile);

}// end common_init

void DFOCC::title()
{
   fprintf(outfile,"\n");
   fprintf(outfile," ============================================================================== \n");
   fprintf(outfile," ============================================================================== \n");
   fprintf(outfile," ============================================================================== \n");
   fprintf(outfile,"\n");
   if (wfn_type_ == "DF-OMP2" && orb_opt_ == "TRUE") fprintf(outfile,"                      DF-OMP2 (DF-OO-MP2)   \n");
   else if (wfn_type_ == "DF-OMP2" && orb_opt_ == "FALSE") fprintf(outfile,"                       DF-MP2   \n");
   else if (wfn_type_ == "DF-OMP3" && orb_opt_ == "TRUE") fprintf(outfile,"                       DF-OMP3 (DF-OO-MP3)   \n");
   else if (wfn_type_ == "DF-OMP3" && orb_opt_ == "FALSE") fprintf(outfile,"                       DF-MP3   \n");
   else if (wfn_type_ == "DF-OCEPA(0)" && orb_opt_ == "TRUE") fprintf(outfile,"                       DF-OCEPA(0) (DF-OO-CEPA)   \n");
   else if (wfn_type_ == "DF-OCEPA(0)" && orb_opt_ == "FALSE") fprintf(outfile,"                       DF-CEPA(0)   \n");
   else if (wfn_type_ == "DF-OMP2.5" && orb_opt_ == "TRUE") fprintf(outfile,"                       DF-OMP2.5 (DF-OO-MP2.5)   \n");
   else if (wfn_type_ == "DF-OMP2.5" && orb_opt_ == "FALSE") fprintf(outfile,"                       DF-MP2.5  \n");
   else if (wfn_type_ == "CD-OMP2" && orb_opt_ == "TRUE") fprintf(outfile,"                      CD-OMP2 (CD-OO-MP2)   \n");
   else if (wfn_type_ == "CD-OMP2" && orb_opt_ == "FALSE") fprintf(outfile,"                       CD-MP2   \n");
   fprintf(outfile,"              Program Written by Ugur Bozkaya\n") ; 
   fprintf(outfile,"              Latest Revision June 17, 2014\n") ;
   fprintf(outfile,"\n");
   fprintf(outfile," ============================================================================== \n");
   fprintf(outfile," ============================================================================== \n");
   fprintf(outfile," ============================================================================== \n");
   fprintf(outfile,"\n");
   fflush(outfile);

}//

void DFOCC::title_grad()
{
   fprintf(outfile,"\n");
   fprintf(outfile," ============================================================================== \n");
   fprintf(outfile," ============================================================================== \n");
   fprintf(outfile," ============================================================================== \n");
   fprintf(outfile,"\n");
   fprintf(outfile,"                         DFGRAD   \n");
   fprintf(outfile,"            A General Analytic Gradients Code   \n");
   fprintf(outfile,"               for Density-Fitted Methods       \n");
   fprintf(outfile,"                   by Ugur Bozkaya\n") ; 
   fprintf(outfile,"              Latest Revision June 17, 2014\n") ;
   fprintf(outfile,"\n");
   fprintf(outfile," ============================================================================== \n");
   fprintf(outfile," ============================================================================== \n");
   fprintf(outfile," ============================================================================== \n");
   fprintf(outfile,"\n");
   fflush(outfile);

}//

double DFOCC::compute_energy()
{   

        // Call the appropriate manager
        do_cd = "FALSE";
        if (wfn_type_ == "DF-OMP2" && orb_opt_ == "TRUE") omp2_manager();
        else if (wfn_type_ == "DF-OMP2" && orb_opt_ == "FALSE") mp2_manager();
        else if (wfn_type_ == "DF-OMP3" && orb_opt_ == "TRUE") omp3_manager();
        else if (wfn_type_ == "DF-OMP3" && orb_opt_ == "FALSE") mp3_manager();
        else if (wfn_type_ == "DF-OCEPA(0)" && orb_opt_ == "TRUE") ocepa_manager();
        else if (wfn_type_ == "DF-OCEPA(0)" && orb_opt_ == "FALSE") cepa_manager();
        else if (wfn_type_ == "DF-OMP2.5" && orb_opt_ == "TRUE") omp2_5_manager();
        else if (wfn_type_ == "DF-OMP2.5" && orb_opt_ == "FALSE") mp2_5_manager();
        else if (wfn_type_ == "CD-OMP2" && orb_opt_ == "TRUE") cd_omp2_manager();
        else if (wfn_type_ == "CD-OMP2" && orb_opt_ == "FALSE") cd_mp2_manager();

        /*
        timer_on("DF-RHF GRAD");
	ref_grad();  
        timer_off("DF-RHF GRAD");
        */

        if (wfn_type_ == "DF-OMP2") Etotal = Emp2L;
        else if (wfn_type_ == "DF-OMP3" || wfn_type_ == "DF-OMP2.5") Etotal = Emp3L;
        else if (wfn_type_ == "DF-OCEPA") Etotal = EcepaL;

        return Etotal;

} // end of compute_energy

}}

