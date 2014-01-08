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

#include <libqt/qt.h>
#include "dfocc.h"

using namespace psi;
using namespace std;


namespace psi{ namespace dfoccwave{
  
//======================================================================
//             CD-OMP2 Manager
//======================================================================             
void DFOCC::cd_omp2_manager()
{
        do_cd = "TRUE";
        time4grad = 0;// means i will not compute the gradient
	mo_optimized = 0;// means MOs are not optimized
	orbs_already_opt = 0;// means orbitals are not optimized yet.
	orbs_already_sc = 0;// menas orbitals are not semicanonical yet.
        timer_on("CD Integrals");
        cd_ints();
        trans_cd();
        timer_off("CD Integrals");

        // memalloc for density intermediates
        Jc = SharedTensor1d(new Tensor1d("DF_BASIS_SCF J_Q", nQ_ref));
        g1Qc = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1_Q", nQ_ref));
        g1Qt = SharedTensor1d(new Tensor1d("DF_BASIS_SCF G1t_Q", nQ_ref));
        g1Q = SharedTensor1d(new Tensor1d("DF_BASIS_CC G1_Q", nQ));
        g1Qt2 = SharedTensor1d(new Tensor1d("DF_BASIS_CC G1t_Q", nQ));

        if (conv_tei_type == "DISK") { 
           tei_oooo_chem_ref();
           tei_ooov_chem_ref();
           tei_oovv_chem_ref();
           tei_ovov_chem_ref(); 
           if (reference_ == "UNRESTRICTED") {
            tei_oooo_phys_ref();
            tei_ooov_phys_ref();
            tei_oovv_phys_ref();
            tei_ovov_phys_ref();
            tei_oooo_anti_symm_ref();
            tei_ooov_anti_symm_ref();
            tei_oovv_anti_symm_ref();
            tei_ovov_anti_symm_ref();
           }

           tei_iajb_chem();
           //tei_ijab_chem();// for Hessian
           if (reference_ == "UNRESTRICTED") {
               tei_ijab_phys();
               tei_ijab_anti_symm();
           }
        }// if (conv_tei_type == "DISK")  
           fock();

        // ROHF REF
        if (reference == "ROHF") t1_1st_sc();
	t2_1st_sc();
	Emp2L=Emp2;
        EcorrL=Emp2L-Escf;
	Emp2L_old=Emp2;
	
	fprintf(outfile,"\n"); 
	if (reference == "ROHF") fprintf(outfile,"\tComputing CD-MP2 energy using SCF MOs (CD-ROHF-MP2)... \n"); 
	else fprintf(outfile,"\tComputing CD-MP2 energy using SCF MOs (Canonical CD-MP2)... \n"); 
	fprintf(outfile,"\t======================================================================= \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	fprintf(outfile,"\tCD-HF Energy (a.u.)                : %20.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	if (reference_ == "UNRESTRICTED") fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	if (reference_ == "UNRESTRICTED") fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	if (reference_ == "UNRESTRICTED") fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	if (reference_ == "UNRESTRICTED") fprintf(outfile,"\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	if (reference_ == "UNRESTRICTED") fprintf(outfile,"\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	if (reference_ == "UNRESTRICTED") fprintf(outfile,"\tCD-SCS-MP2 Total Energy (a.u.)     : %20.14f\n", Escsmp2);
	if (reference_ == "UNRESTRICTED") fprintf(outfile,"\tCD-SOS-MP2 Total Energy (a.u.)     : %20.14f\n", Esosmp2);
	if (reference_ == "UNRESTRICTED") fprintf(outfile,"\tCD-SCSN-MP2 Total Energy (a.u.)    : %20.14f\n", Escsnmp2);
	if (reference == "ROHF") fprintf(outfile,"\tCD-MP2 Singles Energy (a.u.)       : %20.14f\n", Emp2_t1);
	if (reference == "ROHF") fprintf(outfile,"\tCD-MP2 Doubles Energy (a.u.)       : %20.14f\n", Ecorr - Emp2_t1);
	fprintf(outfile,"\tCD-MP2 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
	fprintf(outfile,"\tCD-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
	fprintf(outfile,"\t======================================================================= \n");
	fflush(outfile);
	Process::environment.globals["CURRENT ENERGY"] = Emp2;
	Process::environment.globals["CD-MP2 TOTAL ENERGY"] = Emp2;
	Process::environment.globals["CD-SCS-MP2 TOTAL ENERGY"] = Escsmp2;
	Process::environment.globals["CD-SOS-MP2 TOTAL ENERGY"] = Esosmp2;
	Process::environment.globals["CD-SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;

        Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
        Process::environment.globals["CURRENT CORRELATION ENERGY"] = Emp2 - Escf;
        Process::environment.globals["CD-MP2 CORRELATION ENERGY"] = Emp2 - Escf;
        Process::environment.globals["CD-SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
        Process::environment.globals["CD-SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
        Process::environment.globals["CD-SCSN-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;

        Process::environment.globals["CD-MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        Process::environment.globals["CD-MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA+Emp2BB;

	omp2_opdm();
	omp2_tpdm();
        mp2l_energy();
        separable_tpdm();
	gfock_vo();
	gfock_ov();
        gfock_oo();
        gfock_vv();
	idp();
	mograd();
        occ_iterations();
	
        if (rms_wog <= tol_grad && fabs(DE) >= tol_Eod) {
           orbs_already_opt = 1;
	   if (conver == 1) fprintf(outfile,"\n\tOrbitals are optimized now.\n");
	   else if (conver == 0) { 
                    fprintf(outfile,"\n\tMAX MOGRAD did NOT converged, but RMS MOGRAD converged!!!\n");
	            fprintf(outfile,"\tI will consider the present orbitals as optimized.\n");
           }
	   fprintf(outfile,"\tSwitching to the standard CD-MP2 computation after semicanonicalization of the MOs... \n");
	   fflush(outfile);
	   semi_canonic();
           trans_cd();
        if (conv_tei_type == "DISK") { 
           tei_iajb_chem();
           if (reference_ == "UNRESTRICTED") {
               tei_ijab_phys();
               tei_ijab_anti_symm();
           }
        }// if (conv_tei_type == "DISK") 
        if (conv_tei_type == "DISK") { 
           tei_oooo_chem_ref();
           tei_ooov_chem_ref();
           tei_oovv_chem_ref();
           tei_ovov_chem_ref();
           if (reference_ == "UNRESTRICTED") {
            tei_oooo_phys_ref();
            tei_ooov_phys_ref();
            tei_oovv_phys_ref();
            tei_ovov_phys_ref();
            tei_oooo_anti_symm_ref();
            tei_ooov_anti_symm_ref();
            tei_oovv_anti_symm_ref();
            tei_ovov_anti_symm_ref();
           }
        }// if (conv_tei_type == "DISK") 
           fock();
	   t2_1st_sc();
           conver = 1;
           if (dertype == "FIRST") {
	       omp2_opdm();
	       omp2_tpdm();
               separable_tpdm();
               gfock_vo();
               gfock_ov();
               gfock_oo();
               gfock_vv();
           }
        } 

  if (conver == 1) {
        ref_energy();
	mp2_energy();
        if (orbs_already_opt == 1) Emp2L = Emp2;
	
	fprintf(outfile,"\n"); 
	fprintf(outfile,"\tComputing MP2 energy using optimized MOs... \n"); 
	fprintf(outfile,"\t======================================================================= \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	fprintf(outfile,"\tCD-HF Energy (a.u.)                : %20.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	if (reference_ == "UNRESTRICTED") fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	if (reference_ == "UNRESTRICTED") fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	if (reference_ == "UNRESTRICTED") fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	if (reference_ == "UNRESTRICTED") fprintf(outfile,"\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	if (reference_ == "UNRESTRICTED") fprintf(outfile,"\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	if (reference_ == "UNRESTRICTED") fprintf(outfile,"\tCD-SCS-MP2 Total Energy (a.u.)     : %20.14f\n", Escsmp2);
	if (reference_ == "UNRESTRICTED") fprintf(outfile,"\tCD-SOS-MP2 Total Energy (a.u.)     : %20.14f\n", Esosmp2);
	if (reference_ == "UNRESTRICTED") fprintf(outfile,"\tCD-SCSN-MP2 Total Energy (a.u.)    : %20.14f\n", Escsnmp2);
	fprintf(outfile,"\tCD-MP2 Correlation Energy (a.u.)   : %20.14f\n", Emp2 - Escf);
	fprintf(outfile,"\tCD-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
	fprintf(outfile,"\t======================================================================= \n");
	fflush(outfile);

	fprintf(outfile,"\n");
	fprintf(outfile,"\t======================================================================= \n");
	fprintf(outfile,"\t================ CD-OMP2 FINAL RESULTS ================================ \n");
	fprintf(outfile,"\t======================================================================= \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	fprintf(outfile,"\tCD-HF Energy (a.u.)                : %20.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	if (reference_ == "UNRESTRICTED") fprintf(outfile,"\tCD-SCS-OMP2 Total Energy (a.u.)    : %20.14f\n", Escsmp2);
	if (reference_ == "UNRESTRICTED") fprintf(outfile,"\tCD-SOS-OMP2 Total Energy (a.u.)    : %20.14f\n", Esosmp2);
	if (reference_ == "UNRESTRICTED") fprintf(outfile,"\tCD-SCSN-OMP2 Total Energy (a.u.)   : %20.14f\n", Escsnmp2);
	fprintf(outfile,"\tCD-OMP2 Correlation Energy (a.u.)  : %20.14f\n", Emp2L-Escf);
	fprintf(outfile,"\tEdfomp2 - Eref (a.u.)              : %20.14f\n", Emp2L-Eref);
	fprintf(outfile,"\tCD-OMP2 Total Energy (a.u.)        : %20.14f\n", Emp2L);
	fprintf(outfile,"\t======================================================================= \n");
	fprintf(outfile,"\n");
	fflush(outfile);

	// Set the global variables with the energies
	Process::environment.globals["CD-OMP2 TOTAL ENERGY"] = Emp2L;
	Process::environment.globals["CD-SCS-OMP2 TOTAL ENERGY"] =  Escsmp2;
	Process::environment.globals["CD-SOS-OMP2 TOTAL ENERGY"] =  Esosmp2;
	Process::environment.globals["CD-SCSN-OMP2 TOTAL ENERGY"] = Escsnmp2;
	Process::environment.globals["CURRENT ENERGY"] = Emp2L;
	Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
	Process::environment.globals["CURRENT CORRELATION ENERGY"] = Emp2L-Escf;

        Process::environment.globals["CD-OMP2 CORRELATION ENERGY"] = Emp2L - Escf;
        Process::environment.globals["CD-SCS-OMP2 CORRELATION ENERGY"] =  Escsmp2 - Escf;
        Process::environment.globals["CD-SOS-OMP2 CORRELATION ENERGY"] =  Esosmp2 - Escf;
        Process::environment.globals["CD-SCSN-OMP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;

        // if scs on	
	if (do_scs == "TRUE") {
	    if (scs_type_ == "SCS") {
	       Process::environment.globals["CURRENT ENERGY"] = Escsmp2;
	       Process::environment.globals["CURRENT CORRELATION ENERGY"] = Escsmp2 - Escf;
            }

	    else if (scs_type_ == "SCSN") {
	       Process::environment.globals["CURRENT ENERGY"] = Escsnmp2;
	       Process::environment.globals["CURRENT CORRELATION ENERGY"] = Escsnmp2 - Escf;
            }
	}
    
        // else if sos on	
	else if (do_sos == "TRUE") {
	     if (sos_type_ == "SOS") {
	         Process::environment.globals["CURRENT ENERGY"] = Esosmp2;
  	         Process::environment.globals["CURRENT CORRELATION ENERGY"] = Esosmp2 - Escf;
             }
	}
 
	//if (natorb == "TRUE") nbo();
	//if (occ_orb_energy == "TRUE") semi_canonic(); 

        // Compute Analytic Gradients
        /*
        if (dertype == "FIRST") {
	    fprintf(outfile,"\tAnalytic gradient computation is starting...\n");
	    fflush(outfile);
            coord_grad();
	    fprintf(outfile,"\tNecessary information has been sent to DERIV, which will take care of the rest.\n");
	    fflush(outfile);
        }
        */

  }// end if (conver == 1)
}// end omp2_manager 

//======================================================================
//             CD-MP2 Manager
//======================================================================             
void DFOCC::cd_mp2_manager()
{
        do_cd = "TRUE";
        time4grad = 0;// means i will not compute the gradient
	mo_optimized = 0;// means MOs are not optimized
        timer_on("CD Integrals");
        cd_ints();
        if (dertype == "NONE" && ekt_ip_ == "FALSE" && ekt_ea_ == "FALSE") {
            trans_cd_mp2();
        }
        else trans_cd();
        timer_off("CD Integrals");

        // ROHF REF
        //fprintf(outfile,"\tI am here.\n"); fflush(outfile);
        if (reference == "ROHF") t1_1st_sc();

        if (conv_tei_type == "DISK") { 
            tei_iajb_chem();
            if (reference_ == "UNRESTRICTED") {
                tei_ijab_phys();
                tei_ijab_anti_symm();
            }
        }// if (conv_tei_type == "DISK")  
	//t2_1st_sc();
        //mp2_energy();
        if (dertype == "NONE" && ekt_ip_ == "FALSE" && ekt_ea_ == "FALSE") mp2_direct();
        else {
	     t2_1st_sc();
             mp2_energy();
        }
	Emp2L=Emp2;
        EcorrL=Emp2L-Escf;
	
	fprintf(outfile,"\n"); 
	if (reference == "ROHF") fprintf(outfile,"\tComputing CD-MP2 energy using SCF MOs (CD-ROHF-MP2)... \n"); 
	else fprintf(outfile,"\tComputing CD-MP2 energy using SCF MOs (Canonical CD-MP2)... \n"); 
	fprintf(outfile,"\t======================================================================= \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	fprintf(outfile,"\tCD-HF Energy (a.u.)                : %20.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	if (reference_ == "UNRESTRICTED") fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	if (reference_ == "UNRESTRICTED") fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	if (reference_ == "UNRESTRICTED") fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	if (reference_ == "UNRESTRICTED") fprintf(outfile,"\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	if (reference_ == "UNRESTRICTED") fprintf(outfile,"\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	if (reference_ == "UNRESTRICTED") fprintf(outfile,"\tCD-SCS-MP2 Total Energy (a.u.)     : %20.14f\n", Escsmp2);
	if (reference_ == "UNRESTRICTED") fprintf(outfile,"\tCD-SOS-MP2 Total Energy (a.u.)     : %20.14f\n", Esosmp2);
	if (reference_ == "UNRESTRICTED") fprintf(outfile,"\tCD-SCSN-MP2 Total Energy (a.u.)    : %20.14f\n", Escsnmp2);
	if (reference_ == "ROHF") fprintf(outfile,"\tCD-MP2 Singles Energy (a.u.)       : %20.14f\n", Emp2_t1);
	if (reference_ == "ROHF") fprintf(outfile,"\tCD-MP2 Doubles Energy (a.u.)       : %20.14f\n", Ecorr - Emp2_t1);
	fprintf(outfile,"\tCD-MP2 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
	fprintf(outfile,"\tCD-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
	fprintf(outfile,"\t======================================================================= \n");
	fflush(outfile);
	Process::environment.globals["CURRENT ENERGY"] = Emp2;
	Process::environment.globals["CD-MP2 TOTAL ENERGY"] = Emp2;
	Process::environment.globals["CD-SCS-MP2 TOTAL ENERGY"] = Escsmp2;
	Process::environment.globals["CD-SOS-MP2 TOTAL ENERGY"] = Esosmp2;
	Process::environment.globals["CD-SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;

        Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
        Process::environment.globals["CURRENT CORRELATION ENERGY"] = Emp2 - Escf;
        Process::environment.globals["CD-MP2 CORRELATION ENERGY"] = Emp2 - Escf;
        Process::environment.globals["CD-SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
        Process::environment.globals["CD-SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
        Process::environment.globals["CD-SCSN-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;

        Process::environment.globals["CD-MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        Process::environment.globals["CD-MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA+Emp2BB;

        /*
        // Compute Analytic Gradients
        if (dertype == "FIRST" || ekt_ip_ == "TRUE" || ekt_ea_ == "TRUE") {
	    fprintf(outfile,"\tAnalytic gradient computation is starting...\n");
	    fprintf(outfile,"\tComputing response density matrices...\n");
	    fflush(outfile);
	    omp2_response_pdms();
            fprintf(outfile,"\tComputing off-diagonal blocks of GFM...\n");
            fflush(outfile);
	    gfock();
            fprintf(outfile,"\tForming independent-pairs...\n");
            fflush(outfile);
	    idp2();
            fprintf(outfile,"\tComputing orbital gradient...\n");
            fflush(outfile);
	    mograd();
            coord_grad();

            if (ekt_ip_ == "TRUE" && ekt_ea_ == "TRUE") {
                ekt_ip();
                ekt_ea();
            }

            else if (ekt_ip_ == "TRUE" && ekt_ea_ == "FALSE") {
                ekt_ip();
            }

            else if (ekt_ip_ == "FALSE" && ekt_ea_ == "TRUE") {
                ekt_ea();
            }

            else if (ekt_ip_ == "FALSE" && ekt_ea_ == "FALSE") {
	        fprintf(outfile,"\tNecessary information has been sent to DERIV, which will take care of the rest.\n");
	        fflush(outfile);
            }
        }// if (dertype == "FIRST" || ekt_ip_ == "TRUE" || ekt_ea_ == "TRUE") 
        */

}// end mp2_manager 

}} // End Namespaces


