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
#include "occwave.h"

using namespace boost;
using namespace psi;
using namespace std;


namespace psi{ namespace occwave{
  
//======================================================================
//             OMP2 Manager
//======================================================================             
void OCCWave::omp2_manager()
{
	mo_optimized = 0;
	orbs_already_opt = 0;
	orbs_already_sc = 0;
        timer_on("trans_ints");
	if (reference_ == "RESTRICTED") trans_ints_rhf();  
	else if (reference_ == "UNRESTRICTED") trans_ints_uhf();  
        timer_off("trans_ints");
        timer_on("T2(1)");
	omp2_t2_1st_sc();
        timer_off("T2(1)");
        timer_on("REF Energy");
	ref_energy();
        timer_off("REF Energy");
        timer_on("MP2 Energy");
	omp2_mp2_energy();
        timer_off("MP2 Energy");
	Emp2L=Emp2;
        EcorrL=Emp2L-Escf;
	Emp2L_old=Emp2;
        if (ip_poles == "TRUE") omp2_ip_poles();
        if (ep_ip_poles == "TRUE") ep2_ip();
	
	fprintf(outfile,"\n"); 
	fprintf(outfile,"\tComputing MP2 energy using SCF MOs (Canonical MP2)... \n"); 
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	fprintf(outfile,"\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	fprintf(outfile,"\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	fprintf(outfile,"\tSCS-MP2 Total Energy (a.u.)        : %20.14f\n", Escsmp2);
	fprintf(outfile,"\tSOS-MP2 Total Energy (a.u.)        : %20.14f\n", Esosmp2);
	fprintf(outfile,"\tSCSN-MP2 Total Energy (a.u.)       : %20.14f\n", Escsnmp2);
	fprintf(outfile,"\tSCS-MP2-VDW Total Energy (a.u.)    : %20.14f\n", Escsmp2vdw);
	fprintf(outfile,"\tSOS-PI-MP2 Total Energy (a.u.)     : %20.14f\n", Esospimp2);
	fprintf(outfile,"\tMP2 Correlation Energy (a.u.)      : %20.14f\n", Ecorr);
	fprintf(outfile,"\tMP2 Total Energy (a.u.)            : %20.14f\n", Emp2);
	fprintf(outfile,"\t============================================================================== \n");
	fflush(outfile);
	Process::environment.globals["MP2 TOTAL ENERGY"] = Emp2;
	Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
	Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
	Process::environment.globals["SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;
	Process::environment.globals["SCS-MP2-VDW TOTAL ENERGY"] = Escsmp2vdw;
	Process::environment.globals["SOS-PI-MP2 TOTAL ENERGY"] = Esospimp2;

        Process::environment.globals["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
        Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
        Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
        Process::environment.globals["SCSN-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
        Process::environment.globals["SCS-MP2-VDW CORRELATION ENERGY"] = Escsmp2vdw - Escf;
        Process::environment.globals["SOS-PI-MP2 CORRELATION ENERGY"] = Esospimp2 - Escf;

       Process::environment.globals["MP2 OPPOSITE-SPIN ENERGY"] = Emp2AB;
       Process::environment.globals["MP2 SAME-SPIN ENERGY"] = Emp2AA+Emp2BB;

	omp2_response_pdms();
	gfock();
	idp();
	mograd();
        occ_iterations();
	
        if (rms_wog <= tol_grad && fabs(DE) >= tol_Eod) {
           orbs_already_opt = 1;
	   if (conver == 1) fprintf(outfile,"\n\tOrbitals are optimized now.\n");
	   else if (conver == 0) { 
                    fprintf(outfile,"\n\tMAX MOGRAD did NOEscsmp2 - EscfoT converged, but RMS MOGRAD converged!!!\n");
	            fprintf(outfile,"\tI will consider the present orbitals as optimized.\n");
           }
	   fprintf(outfile,"\tSwitching to the standard MP2 computation after semicanonicalization of the MOs... \n");
	   fflush(outfile);
	   semi_canonic();
	   if (reference_ == "RESTRICTED") trans_ints_rhf();  
	   else if (reference_ == "UNRESTRICTED") trans_ints_uhf();  
	   omp2_t2_1st_sc();
           conver = 1;
           if (dertype == "FIRST") {
               omp2_response_pdms();
	       gfock();
           }
        } 

  if (conver == 1) {
        ref_energy();
	omp2_mp2_energy();
        if (orbs_already_opt == 1) Emp2L = Emp2;

        // Green's function
        if (ip_poles == "TRUE") {
	   if (orbs_already_sc == 0) {
               semi_canonic();
	       if (reference_ == "RESTRICTED") trans_ints_rhf();  
	       else if (reference_ == "UNRESTRICTED") trans_ints_uhf();  
	       omp2_t2_1st_sc();
           }
           omp2_ip_poles();
        }

        if (ep_ip_poles == "TRUE") {
	   if (orbs_already_sc == 0) {
               semi_canonic();
	       if (reference_ == "RESTRICTED") trans_ints_rhf();  
	       else if (reference_ == "UNRESTRICTED") trans_ints_uhf();  
           }
           ep2_ip();
        }

        // EKT
        if (ekt_ip_ == "TRUE") { 
            if (orbs_already_sc == 1) {
	        omp2_response_pdms();
	        gfock();
            }
            gfock_diag();
            ekt_ip();
        }
	
	fprintf(outfile,"\n"); 
	fprintf(outfile,"\tComputing MP2 energy using optimized MOs... \n"); 
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	fprintf(outfile,"\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	fprintf(outfile,"\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	fprintf(outfile,"\tSCS-MP2 Total Energy (a.u.)        : %20.14f\n", Escsmp2);
	fprintf(outfile,"\tSOS-MP2 Total Energy (a.u.)        : %20.14f\n", Esosmp2);
	fprintf(outfile,"\tSCSN-MP2 Total Energy (a.u.)       : %20.14f\n", Escsnmp2);
	fprintf(outfile,"\tSCS-MP2-VDW Total Energy (a.u.)    : %20.14f\n", Escsmp2vdw);
	fprintf(outfile,"\tSOS-PI-MP2 Total Energy (a.u.)     : %20.14f\n", Esospimp2);
	fprintf(outfile,"\tMP2 Correlation Energy (a.u.)      : %20.14f\n", Ecorr);
	fprintf(outfile,"\tMP2 Total Energy (a.u.)            : %20.14f\n", Emp2);
	fprintf(outfile,"\t============================================================================== \n");
	fflush(outfile);


	fprintf(outfile,"\n");
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\t================ OMP2 FINAL RESULTS ========================================== \n");
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	fprintf(outfile,"\tSCS-OMP2 Total Energy (a.u.)       : %20.14f\n", Escsmp2);
	fprintf(outfile,"\tSOS-OMP2 Total Energy (a.u.)       : %20.14f\n", Esosmp2);
	fprintf(outfile,"\tSCSN-OMP2 Total Energy (a.u.)      : %20.14f\n", Escsnmp2);
	fprintf(outfile,"\tSCS-OMP2-VDW Total Energy (a.u.)   : %20.14f\n", Escsmp2vdw);
	fprintf(outfile,"\tSOS-PI-OMP2 Total Energy (a.u.)    : %20.14f\n", Esospimp2);
	fprintf(outfile,"\tOMP2 Correlation Energy (a.u.)     : %20.14f\n", Emp2L-Escf);
	fprintf(outfile,"\tEomp2 - Eref (a.u.)                : %20.14f\n", Emp2L-Eref);
	fprintf(outfile,"\tOMP2 Total Energy (a.u.)           : %20.14f\n", Emp2L);
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\n");
	fflush(outfile);

	// Set the global variables with the energies
	Process::environment.globals["OMP2 TOTAL ENERGY"] = Emp2L;
	Process::environment.globals["SCS-OMP2 TOTAL ENERGY"] =  Escsmp2;
	Process::environment.globals["SOS-OMP2 TOTAL ENERGY"] =  Esosmp2;
	Process::environment.globals["SCSN-OMP2 TOTAL ENERGY"] = Escsnmp2;
	Process::environment.globals["SCS-OMP2-VDW TOTAL ENERGY"] = Escsmp2vdw;
	Process::environment.globals["SOS-PI-OMP2 TOTAL ENERGY"] = Esospimp2;
	Process::environment.globals["CURRENT ENERGY"] = Emp2L;
	Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
	Process::environment.globals["CURRENT CORRELATION ENERGY"] = Emp2L-Escf;

        Process::environment.globals["OMP2 CORRELATION ENERGY"] = Emp2L - Escf;
        Process::environment.globals["SCS-OMP2 CORRELATION ENERGY"] =  Escsmp2 - Escf;
        Process::environment.globals["SOS-OMP2 CORRELATION ENERGY"] =  Esosmp2 - Escf;
        Process::environment.globals["SCSN-OMP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
        Process::environment.globals["SCS-OMP2-VDW CORRELATION ENERGY"] = Escsmp2vdw - Escf;
        Process::environment.globals["SOS-PI-OMP2 CORRELATION ENERGY"] = Esospimp2 - Escf;

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


	    else if (scs_type_ == "SCSVDW") {
	       Process::environment.globals["CURRENT ENERGY"] = Escsmp2vdw;
	       Process::environment.globals["CURRENT CORRELATION ENERGY"] = Escsmp2vdw - Escf;
            }
	}
    
        // else if sos on	
	else if (do_sos == "TRUE") {
	     if (sos_type_ == "SOS") {
	         Process::environment.globals["CURRENT ENERGY"] = Esosmp2;
  	         Process::environment.globals["CURRENT CORRELATION ENERGY"] = Esosmp2 - Escf;
             }

	     else if (sos_type_ == "SOSPI") {
	             Process::environment.globals["CURRENT ENERGY"] = Esospimp2;
	             Process::environment.globals["CURRENT CORRELATION ENERGY"] = Esospimp2 - Escf;
             }
	}
 
	if (natorb == "TRUE") nbo();
	if (occ_orb_energy == "TRUE") semi_canonic(); 

        // Compute Analytic Gradients
        if (dertype == "FIRST") {
	    fprintf(outfile,"\tAnalytic gradient computation is starting...\n");
	    fflush(outfile);
            coord_grad();
	    fprintf(outfile,"\tNecessary information has been sent to DERIV, which will take care of the rest.\n");
	    fflush(outfile);
        }

  }// end if (conver == 1)
}// end omp2_manager 

//======================================================================
//             MP2 Manager
//======================================================================             
void OCCWave::mp2_manager()
{
        time4grad = 0;// means i will not compute the gradient
        timer_on("trans_ints");
        if (reference_ == "RESTRICTED" && dertype == "NONE") trans_ints_rmp2(); 
        else if (reference_ == "RESTRICTED" && dertype == "FIRST") trans_ints_rhf();
	else if (reference_ == "UNRESTRICTED" && dertype == "FIRST") trans_ints_uhf();  
	else if (reference_ == "UNRESTRICTED" && dertype == "NONE" && reference == "ROHF") trans_ints_uhf();  
	else if (reference_ == "UNRESTRICTED" && dertype == "NONE" && reference != "ROHF") trans_ints_ump2();  
        timer_off("trans_ints");
        // ROHF REF
        if (reference == "ROHF") {
        timer_on("T1(1)");
	t1_1st_sc();
        timer_off("T1(1)");
        }// end if (reference == "ROHF") 
        timer_on("T2(1)");
	omp2_t2_1st_sc();
        timer_off("T2(1)");
        Eref = Escf;
        timer_on("MP2 Energy");
	omp2_mp2_energy();
        timer_off("MP2 Energy");
	Emp2L=Emp2;
        EcorrL=Emp2L-Escf;
	Emp2L_old=Emp2;
        if (ip_poles == "TRUE") omp2_ip_poles();
        if (ep_ip_poles == "TRUE") ep2_ip();
	
	fprintf(outfile,"\n"); 
	if (reference == "ROHF") fprintf(outfile,"\tComputing MP2 energy using SCF MOs (ROHF-MP2)... \n"); 
	else fprintf(outfile,"\tComputing MP2 energy using SCF MOs (Canonical MP2)... \n"); 
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	fprintf(outfile,"\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	fprintf(outfile,"\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	fprintf(outfile,"\tSCS-MP2 Total Energy (a.u.)        : %20.14f\n", Escsmp2);
	fprintf(outfile,"\tSOS-MP2 Total Energy (a.u.)        : %20.14f\n", Esosmp2);
	fprintf(outfile,"\tSCSN-MP2 Total Energy (a.u.)       : %20.14f\n", Escsnmp2);
	fprintf(outfile,"\tSCS-MP2-VDW Total Energy (a.u.)    : %20.14f\n", Escsmp2vdw);
	fprintf(outfile,"\tSOS-PI-MP2 Total Energy (a.u.)     : %20.14f\n", Esospimp2);
	if (reference == "ROHF") fprintf(outfile,"\tMP2 Singles Energy (a.u.)          : %20.14f\n", Emp2_t1);
	if (reference == "ROHF") fprintf(outfile,"\tMP2 Doubles Energy (a.u.)          : %20.14f\n", Ecorr - Emp2_t1);
	fprintf(outfile,"\tMP2 Correlation Energy (a.u.)      : %20.14f\n", Ecorr);
	fprintf(outfile,"\tMP2 Total Energy (a.u.)            : %20.14f\n", Emp2);
	fprintf(outfile,"\t============================================================================== \n");
	fflush(outfile);
	Process::environment.globals["CURRENT ENERGY"] = Emp2;
	Process::environment.globals["MP2 TOTAL ENERGY"] = Emp2;
	Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
	Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
	Process::environment.globals["SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;
	Process::environment.globals["SCS-MP2-VDW TOTAL ENERGY"] = Escsmp2vdw;
	Process::environment.globals["SOS-PI-MP2 TOTAL ENERGY"] = Esospimp2;

        Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
        Process::environment.globals["CURRENT CORRELATION ENERGY"] = Emp2 - Escf;
        Process::environment.globals["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
        Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
        Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
        Process::environment.globals["SCSN-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
        Process::environment.globals["SCS-MP2-VDW CORRELATION ENERGY"] = Escsmp2vdw - Escf;
        Process::environment.globals["SOS-PI-MP2 CORRELATION ENERGY"] = Esospimp2 - Escf;

       Process::environment.globals["MP2 OPPOSITE-SPIN ENERGY"] = Emp2AB;
       Process::environment.globals["MP2 SAME-SPIN ENERGY"] = Emp2AA+Emp2BB;

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

	    else if (scs_type_ == "SCSVDW") {
	       Process::environment.globals["CURRENT ENERGY"] = Escsmp2vdw;
	       Process::environment.globals["CURRENT CORRELATION ENERGY"] = Escsmp2vdw - Escf;
            }
	}
    
        // else if sos on	
	else if (do_sos == "TRUE") {
	     if (sos_type_ == "SOS") {
	         Process::environment.globals["CURRENT ENERGY"] = Esosmp2;
  	         Process::environment.globals["CURRENT CORRELATION ENERGY"] = Esosmp2 - Escf;
             }

	     else if (sos_type_ == "SOSPI") {
	             Process::environment.globals["CURRENT ENERGY"] = Esospimp2;
	             Process::environment.globals["CURRENT CORRELATION ENERGY"] = Esospimp2 - Escf;
             }
	}
 
        // Compute Analytic Gradients
        if (dertype == "FIRST" || ekt_ip_ == "TRUE") {
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
            if (ekt_ip_ == "TRUE") ekt_ip();
	    fprintf(outfile,"\tNecessary information has been sent to DERIV, which will take care of the rest.\n");
	    fflush(outfile);
        }

}// end mp2_manager 


//======================================================================
//             OMP3 Manager
//======================================================================             
void OCCWave::omp3_manager()
{
	mo_optimized = 0;
	orbs_already_opt = 0;
	orbs_already_sc = 0;
        time4grad = 0;// means i will not compute the gradient
        timer_on("trans_ints");
	if (reference_ == "RESTRICTED") trans_ints_rhf();  
	else if (reference_ == "UNRESTRICTED") trans_ints_uhf();  
        timer_off("trans_ints");
        timer_on("REF Energy");
	ref_energy();
        timer_off("REF Energy");
        timer_on("T2(1)");
	omp3_t2_1st_sc();
        timer_off("T2(1)");
        timer_on("MP2 Energy");
	omp3_mp2_energy();
        timer_off("MP2 Energy");

	fprintf(outfile,"\n"); 
	fprintf(outfile,"\tComputing MP2 energy using SCF MOs (Canonical MP2)... \n"); 
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	fprintf(outfile,"\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	fprintf(outfile,"\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	fprintf(outfile,"\tSCS-MP2 Total Energy (a.u.)        : %20.14f\n", Escsmp2);
	fprintf(outfile,"\tSOS-MP2 Total Energy (a.u.)        : %20.14f\n", Esosmp2);
	fprintf(outfile,"\tSCSN-MP2 Total Energy (a.u.)       : %20.14f\n", Escsnmp2);
	fprintf(outfile,"\tSCS-MP2-VDW Total Energy (a.u.)    : %20.14f\n", Escsmp2vdw);
	fprintf(outfile,"\tSOS-PI-MP2 Total Energy (a.u.)     : %20.14f\n", Esospimp2);
	fprintf(outfile,"\tMP2 Correlation Energy (a.u.)      : %20.14f\n", Ecorr);
	fprintf(outfile,"\tMP2 Total Energy (a.u.)            : %20.14f\n", Emp2);
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\n"); 
	fflush(outfile);
	Process::environment.globals["MP2 TOTAL ENERGY"] = Emp2;
	Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
	Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
	Process::environment.globals["SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;
	Process::environment.globals["SCS-MP2-VDW TOTAL ENERGY"] = Escsmp2vdw;
	Process::environment.globals["SOS-PI-MP2 TOTAL ENERGY"] = Esospimp2;

        Process::environment.globals["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
        Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
        Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
        Process::environment.globals["SCSN-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
        Process::environment.globals["SCS-MP2-VDW CORRELATION ENERGY"] = Escsmp2vdw - Escf;
        Process::environment.globals["SOS-PI-MP2 CORRELATION ENERGY"] = Esospimp2 - Escf;

        Process::environment.globals["MP2 OPPOSITE-SPIN ENERGY"] = Emp2AB;
        Process::environment.globals["MP2 SAME-SPIN ENERGY"] = Emp2AA+Emp2BB;

        timer_on("T2(2)");
	t2_2nd_sc();
        timer_off("T2(2)");
        timer_on("MP3 Energy");
	mp3_energy();
        timer_off("MP3 Energy");
	Emp3L=Emp3;
        EcorrL=Emp3L-Escf;
	Emp3L_old=Emp3;
        if (ip_poles == "TRUE") omp3_ip_poles();
	
	fprintf(outfile,"\n"); 
	fprintf(outfile,"\tComputing MP3 energy using SCF MOs (Canonical MP3)... \n"); 
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp3AA);
	fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp3AB);
	fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp3BB);
	fprintf(outfile,"\tMP2.5 Correlation Energy (a.u.)    : %20.14f\n", (Emp2 - Escf) + 0.5 * (Emp3-Emp2));
	fprintf(outfile,"\tMP2.5 Total Energy (a.u.)          : %20.14f\n", 0.5 * (Emp3+Emp2));
	fprintf(outfile,"\tSCS-MP3 Total Energy (a.u.)        : %20.14f\n", Escsmp3);
	fprintf(outfile,"\tSOS-MP3 Total Energy (a.u.)        : %20.14f\n", Esosmp3);
	fprintf(outfile,"\tSCSN-MP3 Total Energy (a.u.)       : %20.14f\n", Escsnmp3);
	fprintf(outfile,"\tSCS-MP3-VDW Total Energy (a.u.)    : %20.14f\n", Escsmp3vdw);
	fprintf(outfile,"\tSOS-PI-MP3 Total Energy (a.u.)     : %20.14f\n", Esospimp3);
	fprintf(outfile,"\t3rd Order Energy (a.u.)            : %20.14f\n", Emp3-Emp2);
	fprintf(outfile,"\tMP3 Correlation Energy (a.u.)      : %20.14f\n", Ecorr);
	fprintf(outfile,"\tMP3 Total Energy (a.u.)            : %20.14f\n", Emp3);
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\n"); 
	fflush(outfile);
	Process::environment.globals["MP3 TOTAL ENERGY"] = Emp3;
	Process::environment.globals["SCS-MP3 TOTAL ENERGY"] = Escsmp3;
	Process::environment.globals["SOS-MP3 TOTAL ENERGY"] = Esosmp3;
	Process::environment.globals["SCSN-MP3 TOTAL ENERGY"] = Escsnmp3;
	Process::environment.globals["SCS-MP3-VDW TOTAL ENERGY"] = Escsmp3vdw;
	Process::environment.globals["SOS-PI-MP3 TOTAL ENERGY"] = Esospimp3;

	Process::environment.globals["MP2.5 CORRELATION ENERGY"] = (Emp2 - Escf) + 0.5 * (Emp3-Emp2);
	Process::environment.globals["MP2.5 TOTAL ENERGY"] = 0.5 * (Emp3+Emp2);
	Process::environment.globals["MP3 CORRELATION ENERGY"] = Emp3 - Escf;
	Process::environment.globals["SCS-MP3 CORRELATION ENERGY"] = Escsmp3 - Escf;
	Process::environment.globals["SOS-MP3 CORRELATION ENERGY"] = Esosmp3 - Escf;
	Process::environment.globals["SCSN-MP3 CORRELATION ENERGY"] = Escsnmp3 - Escf;
	Process::environment.globals["SCS-MP3-VDW CORRELATION ENERGY"] = Escsmp3vdw - Escf;
	Process::environment.globals["SOS-PI-MP3 CORRELATION ENERGY"] = Esospimp3 - Escf;

	omp3_response_pdms();
	gfock();
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
	   fprintf(outfile,"\tSwitching to the standard MP3 computation after semicanonicalization of the MOs... \n");
	   fflush(outfile);
	   semi_canonic();
	   if (reference_ == "RESTRICTED") trans_ints_rhf();  
	   else if (reference_ == "UNRESTRICTED") trans_ints_uhf();  
	   omp3_t2_1st_sc();
	   t2_2nd_sc();
           conver = 1;
           if (dertype == "FIRST") {
               omp3_response_pdms();
	       gfock();
           }
        }     

  if (conver == 1) {
        ref_energy();
	omp3_mp2_energy();
	mp3_energy();
        if (orbs_already_opt == 1) Emp3L = Emp3;

        if (ip_poles == "TRUE") {
	   if (orbs_already_sc == 0) {
               semi_canonic();
	       if (reference_ == "RESTRICTED") trans_ints_rhf();  
	       else if (reference_ == "UNRESTRICTED") trans_ints_uhf();  
	       omp3_t2_1st_sc();
	       t2_2nd_sc();
           }
           omp3_ip_poles();
        }

        // EKT
        if (ekt_ip_ == "TRUE") { 
            if (orbs_already_sc == 1) {
	        omp3_response_pdms();
	        gfock();
            }
            gfock_diag();
            ekt_ip();
        }

        fprintf(outfile,"\n"); 
	fprintf(outfile,"\tComputing MP2 energy using optimized MOs... \n"); 
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	fprintf(outfile,"\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	fprintf(outfile,"\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	fprintf(outfile,"\tSCS-MP2 Total Energy (a.u.)        : %20.14f\n", Escsmp2);
	fprintf(outfile,"\tSOS-MP2 Total Energy (a.u.)        : %20.14f\n", Esosmp2);
	fprintf(outfile,"\tSCSN-MP2 Total Energy (a.u.)       : %20.14f\n", Escsnmp2);
	fprintf(outfile,"\tSCS-MP2-VDW Total Energy (a.u.)    : %20.14f\n", Escsmp2vdw);
	fprintf(outfile,"\tSOS-PI-MP2 Total Energy (a.u.)     : %20.14f\n", Esospimp2);
	fprintf(outfile,"\tMP2 Correlation Energy (a.u.)      : %20.14f\n", Ecorr);
	fprintf(outfile,"\tMP2 Total Energy (a.u.)            : %20.14f\n", Emp2);
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\n"); 
	fflush(outfile);

	fprintf(outfile,"\n"); 
	fprintf(outfile,"\tComputing MP3 energy using optimized MOs... \n"); 
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp3AA);
	fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp3AB);
	fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp3BB);
	fprintf(outfile,"\tMP2.5 Correlation Energy (a.u.)    : %20.14f\n", (Emp2 - Escf) + 0.5 * (Emp3-Emp2));
	fprintf(outfile,"\tMP2.5 Total Energy (a.u.)          : %20.14f\n", 0.5 * (Emp3+Emp2));
	fprintf(outfile,"\tSCS-MP3 Total Energy (a.u.)        : %20.14f\n", Escsmp3);
	fprintf(outfile,"\tSOS-MP3 Total Energy (a.u.)        : %20.14f\n", Esosmp3);
	fprintf(outfile,"\tSCSN-MP3 Total Energy (a.u.)       : %20.14f\n", Escsnmp3);
	fprintf(outfile,"\tSCS-MP3-VDW Total Energy (a.u.)    : %20.14f\n", Escsmp3vdw);
	fprintf(outfile,"\tSOS-PI-MP3 Total Energy (a.u.)     : %20.14f\n", Esospimp3);
	fprintf(outfile,"\t3rd Order Energy (a.u.)            : %20.14f\n", Emp3-Emp2);
	fprintf(outfile,"\tMP3 Correlation Energy (a.u.)      : %20.14f\n", Ecorr);
	fprintf(outfile,"\tMP3 Total Energy (a.u.)            : %20.14f\n", Emp3);
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\n"); 
	fflush(outfile);


	fprintf(outfile,"\n");
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\t================ OMP3 FINAL RESULTS ========================================== \n");
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	fprintf(outfile,"\tSCS-OMP3 Total Energy (a.u.)       : %20.14f\n", Escsmp3);
	fprintf(outfile,"\tSOS-OMP3 Total Energy (a.u.)       : %20.14f\n", Esosmp3);
	fprintf(outfile,"\tSCSN-OMP3 Total Energy (a.u.)      : %20.14f\n", Escsnmp3);
	fprintf(outfile,"\tSCS-OMP3-VDW Total Energy (a.u.    : %20.14f\n", Escsmp3vdw);
	fprintf(outfile,"\tSOS-PI-OMP3 Total Energy (a.u.)    : %20.14f\n", Esospimp3);
	fprintf(outfile,"\tOMP3 Correlation Energy (a.u.)     : %20.14f\n", Emp3L-Escf);
	fprintf(outfile,"\tEomp3 - Eref (a.u.)                : %20.14f\n", Emp3L-Eref);
	fprintf(outfile,"\tOMP3 Total Energy (a.u.)           : %20.14f\n", Emp3L);
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\n");
	fflush(outfile);
	
	// Set the global variables with the energies
	Process::environment.globals["OMP3 TOTAL ENERGY"] = Emp3L;
	Process::environment.globals["SCS-OMP3 TOTAL ENERGY"] =  Escsmp3;
	Process::environment.globals["SOS-OMP3 TOTAL ENERGY"] =  Esosmp3;
	Process::environment.globals["SCSN-OMP3 TOTAL ENERGY"] = Escsnmp3;
	Process::environment.globals["SCS-OMP3-VDW TOTAL ENERGY"] = Escsmp3vdw;
	Process::environment.globals["SOS-PI-OMP3 TOTAL ENERGY"] = Esospimp3;
	Process::environment.globals["CURRENT ENERGY"] = Emp3L;
	Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
	Process::environment.globals["CURRENT CORRELATION ENERGY"] = Emp3L-Escf;

	Process::environment.globals["OMP3 CORRELATION ENERGY"] = Emp3L - Escf;
	Process::environment.globals["SCS-OMP3 CORRELATION ENERGY"] =  Escsmp3 - Escf;
	Process::environment.globals["SOS-OMP3 CORRELATION ENERGY"] =  Esosmp3 - Escf;
	Process::environment.globals["SCSN-OMP3 CORRELATION ENERGY"] = Escsnmp3 - Escf;
	Process::environment.globals["SCS-OMP3-VDW CORRELATION ENERGY"] = Escsmp3vdw - Escf;
	Process::environment.globals["SOS-PI-OMP3 CORRELATION ENERGY"] = Esospimp3 - Escf;

        // if scs on	
	if (do_scs == "TRUE") {
	    if (scs_type_ == "SCS") {
	       Process::environment.globals["CURRENT ENERGY"] = Escsmp3;
	       Process::environment.globals["CURRENT CORRELATION ENERGY"] = Escsmp3 - Escf;
            }

	    else if (scs_type_ == "SCSN") {
	       Process::environment.globals["CURRENT ENERGY"] = Escsnmp3;
	       Process::environment.globals["CURRENT CORRELATION ENERGY"] = Escsnmp3 - Escf;
            }

	    else if (scs_type_ == "SCSVDW") {
	       Process::environment.globals["CURRENT ENERGY"] = Escsmp3vdw;
	       Process::environment.globals["CURRENT CORRELATION ENERGY"] = Escsmp3vdw - Escf;
            }
	}
    
        // else if sos on	
	else if (do_sos == "TRUE") {
	     if (sos_type_ == "SOS") {
	         Process::environment.globals["CURRENT ENERGY"] = Esosmp3;
	         Process::environment.globals["CURRENT CORRELATION ENERGY"] = Esosmp3 - Escf;
             }

	     else if (sos_type_ == "SOSPI") {
	             Process::environment.globals["CURRENT ENERGY"] = Esospimp3;
	             Process::environment.globals["CURRENT CORRELATION ENERGY"] = Esospimp3 - Escf;
             }
	}

	if (natorb == "TRUE") nbo();
	if (occ_orb_energy == "TRUE") semi_canonic(); 

        // Compute Analytic Gradients
        if (dertype == "FIRST") {
            time4grad = 1;
	    fprintf(outfile,"\tAnalytic gradient computation is starting...\n");
	    fflush(outfile);
            coord_grad();
	    fprintf(outfile,"\tNecessary information has been sent to DERIV, which will take care of the rest.\n");
	    fflush(outfile);
        }

  }// end if (conver == 1)
}// end omp3_manager 

//======================================================================
//             MP3 Manager
//======================================================================             
void OCCWave::mp3_manager()
{
        time4grad = 0;// means i will not compute the gradient
        timer_on("trans_ints");
	if (reference_ == "RESTRICTED") trans_ints_rhf();  
	else if (reference_ == "UNRESTRICTED") trans_ints_uhf();  
        timer_off("trans_ints");
        Eref = Escf;
        timer_on("T2(1)");
	omp3_t2_1st_sc();
        timer_off("T2(1)");
        timer_on("MP2 Energy");
	omp3_mp2_energy();
        timer_off("MP2 Energy");

	fprintf(outfile,"\n"); 
	fprintf(outfile,"\tComputing MP2 energy using SCF MOs (Canonical MP2)... \n"); 
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	fprintf(outfile,"\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	fprintf(outfile,"\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	fprintf(outfile,"\tSCS-MP2 Total Energy (a.u.)        : %20.14f\n", Escsmp2);
	fprintf(outfile,"\tSOS-MP2 Total Energy (a.u.)        : %20.14f\n", Esosmp2);
	fprintf(outfile,"\tSCSN-MP2 Total Energy (a.u.)       : %20.14f\n", Escsnmp2);
	fprintf(outfile,"\tSCS-MP2-VDW Total Energy (a.u.)    : %20.14f\n", Escsmp2vdw);
	fprintf(outfile,"\tSOS-PI-MP2 Total Energy (a.u.)     : %20.14f\n", Esospimp2);
	fprintf(outfile,"\tMP2 Correlation Energy (a.u.)      : %20.14f\n", Ecorr);
	fprintf(outfile,"\tMP2 Total Energy (a.u.)            : %20.14f\n", Emp2);
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\n"); 
	fflush(outfile);
	Process::environment.globals["MP2 TOTAL ENERGY"] = Emp2;
	Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
	Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
	Process::environment.globals["SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;
	Process::environment.globals["SCS-MP2-VDW TOTAL ENERGY"] = Escsmp2vdw;
	Process::environment.globals["SOS-PI-MP2 TOTAL ENERGY"] = Esospimp2;

        Process::environment.globals["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
        Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
        Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
        Process::environment.globals["SCSN-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
        Process::environment.globals["SCS-MP2-VDW CORRELATION ENERGY"] = Escsmp2vdw - Escf;
        Process::environment.globals["SOS-PI-MP2 CORRELATION ENERGY"] = Esospimp2 - Escf;

        Process::environment.globals["MP2 OPPOSITE-SPIN ENERGY"] = Emp2AB;
        Process::environment.globals["MP2 SAME-SPIN ENERGY"] = Emp2AA+Emp2BB;

        timer_on("T2(2)");
	t2_2nd_sc();
        timer_off("T2(2)");
        timer_on("MP3 Energy");
	mp3_energy();
        timer_off("MP3 Energy");
	Emp3L=Emp3;
        EcorrL=Emp3L-Escf;
	Emp3L_old=Emp3;
        if (ip_poles == "TRUE") omp3_ip_poles();
	
	fprintf(outfile,"\n"); 
	fprintf(outfile,"\tComputing MP3 energy using SCF MOs (Canonical MP3)... \n"); 
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp3AA);
	fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp3AB);
	fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp3BB);
	fprintf(outfile,"\tMP2.5 Correlation Energy (a.u.)    : %20.14f\n", (Emp2 - Escf) + 0.5 * (Emp3-Emp2));
	fprintf(outfile,"\tMP2.5 Total Energy (a.u.)          : %20.14f\n", 0.5 * (Emp3+Emp2));
	fprintf(outfile,"\tSCS-MP3 Total Energy (a.u.)        : %20.14f\n", Escsmp3);
	fprintf(outfile,"\tSOS-MP3 Total Energy (a.u.)        : %20.14f\n", Esosmp3);
	fprintf(outfile,"\tSCSN-MP3 Total Energy (a.u.)       : %20.14f\n", Escsnmp3);
	fprintf(outfile,"\tSCS-MP3-VDW Total Energy (a.u.)    : %20.14f\n", Escsmp3vdw);
	fprintf(outfile,"\tSOS-PI-MP3 Total Energy (a.u.)     : %20.14f\n", Esospimp3);
	fprintf(outfile,"\t3rd Order Energy (a.u.)            : %20.14f\n", Emp3-Emp2);
	fprintf(outfile,"\tMP3 Correlation Energy (a.u.)      : %20.14f\n", Ecorr);
	fprintf(outfile,"\tMP3 Total Energy (a.u.)            : %20.14f\n", Emp3);
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\n"); 
	fflush(outfile);
	Process::environment.globals["CURRENT ENERGY"] = Emp3;
	Process::environment.globals["CURRENT CORRELATION ENERGY"] = Emp3 - Escf;
	Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;

	Process::environment.globals["MP3 TOTAL ENERGY"] = Emp3;
	Process::environment.globals["SCS-MP3 TOTAL ENERGY"] = Escsmp3;
	Process::environment.globals["SOS-MP3 TOTAL ENERGY"] = Esosmp3;
	Process::environment.globals["SCSN-MP3 TOTAL ENERGY"] = Escsnmp3;
	Process::environment.globals["SCS-MP3-VDW TOTAL ENERGY"] = Escsmp3vdw;
	Process::environment.globals["SOS-PI-MP3 TOTAL ENERGY"] = Esospimp3;

	Process::environment.globals["MP2.5 CORRELATION ENERGY"] = (Emp2 - Escf) + 0.5 * (Emp3-Emp2);
	Process::environment.globals["MP2.5 TOTAL ENERGY"] = 0.5 * (Emp3+Emp2);
	Process::environment.globals["MP3 CORRELATION ENERGY"] = Emp3 - Escf;
	Process::environment.globals["SCS-MP3 CORRELATION ENERGY"] = Escsmp3 - Escf;
	Process::environment.globals["SOS-MP3 CORRELATION ENERGY"] = Esosmp3 - Escf;
	Process::environment.globals["SCSN-MP3 CORRELATION ENERGY"] = Escsnmp3 - Escf;
	Process::environment.globals["SCS-MP3-VDW CORRELATION ENERGY"] = Escsmp3vdw - Escf;
	Process::environment.globals["SOS-PI-MP3 CORRELATION ENERGY"] = Esospimp3 - Escf;

        // EKT
        if (ekt_ip_ == "TRUE") { 
	    omp3_response_pdms();
	    gfock();
            gfock_diag();
            ekt_ip();
        }

        // if scs on	
	if (do_scs == "TRUE") {
	    if (scs_type_ == "SCS") {
	       Process::environment.globals["CURRENT ENERGY"] = Escsmp3;
	       Process::environment.globals["CURRENT CORRELATION ENERGY"] = Escsmp3 - Escf;
            }

	    else if (scs_type_ == "SCSN") {
	       Process::environment.globals["CURRENT ENERGY"] = Escsnmp3;
	       Process::environment.globals["CURRENT CORRELATION ENERGY"] = Escsnmp3 - Escf;
            }

	    else if (scs_type_ == "SCSVDW") {
	       Process::environment.globals["CURRENT ENERGY"] = Escsmp3vdw;
	       Process::environment.globals["CURRENT CORRELATION ENERGY"] = Escsmp3vdw - Escf;
            }
	}
    
        // else if sos on	
	else if (do_sos == "TRUE") {
	     if (sos_type_ == "SOS") {
	         Process::environment.globals["CURRENT ENERGY"] = Esosmp3;
	         Process::environment.globals["CURRENT CORRELATION ENERGY"] = Esosmp3 - Escf;
             }

	     else if (sos_type_ == "SOSPI") {
	             Process::environment.globals["CURRENT ENERGY"] = Esospimp3;
	             Process::environment.globals["CURRENT CORRELATION ENERGY"] = Esospimp3 - Escf;
             }
	}

        // Compute Analytic Gradients
        if (dertype == "FIRST" || ekt_ip_ == "TRUE") {
            time4grad = 1;
	    fprintf(outfile,"\tAnalytic gradient computation is starting...\n");
            fprintf(outfile,"\tComputing response density matrices...\n");
            fflush(outfile);
	    omp3_response_pdms();
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
            if (ekt_ip_ == "TRUE") ekt_ip();
	    fprintf(outfile,"\tNecessary information has been sent to DERIV, which will take care of the rest.\n");
	    fflush(outfile);
        }

}// end mp3_manager 


//======================================================================
//             OCEPA Manager
//======================================================================             
void OCCWave::ocepa_manager()
{
	mo_optimized = 0;// means MOs are not optimized yet.
	orbs_already_opt = 0;
	orbs_already_sc = 0;
        time4grad = 0;// means i will not compute the gradient
        timer_on("trans_ints");
	if (reference_ == "RESTRICTED") trans_ints_rhf();  
	else if (reference_ == "UNRESTRICTED") trans_ints_uhf();  
        timer_off("trans_ints");
        timer_on("REF Energy");
	ref_energy();
        timer_off("REF Energy");
        timer_on("T2(1)");
	ocepa_t2_1st_sc();
        timer_off("T2(1)");
        timer_on("MP2 Energy");
	ocepa_mp2_energy();
        timer_off("MP2 Energy");
        Ecepa = Emp2;
	EcepaL = Ecepa;
        EcorrL = Ecorr;
	EcepaL_old = Ecepa;

	fprintf(outfile,"\n"); 
	fprintf(outfile,"\tComputing MP2 energy using SCF MOs (Canonical MP2)... \n"); 
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	fprintf(outfile,"\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	fprintf(outfile,"\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	fprintf(outfile,"\tSCS-MP2 Total Energy (a.u.)        : %20.14f\n", Escsmp2);
	fprintf(outfile,"\tSOS-MP2 Total Energy (a.u.)        : %20.14f\n", Esosmp2);
	fprintf(outfile,"\tSCSN-MP2 Total Energy (a.u.)       : %20.14f\n", Escsnmp2);
	fprintf(outfile,"\tSCS-MP2-VDW Total Energy (a.u.)    : %20.14f\n", Escsmp2vdw);
	fprintf(outfile,"\tSOS-PI-MP2 Total Energy (a.u.)     : %20.14f\n", Esospimp2);
	fprintf(outfile,"\tMP2 Correlation Energy (a.u.)      : %20.14f\n", Ecorr);
	fprintf(outfile,"\tMP2 Total Energy (a.u.)            : %20.14f\n", Emp2);
	fprintf(outfile,"\t============================================================================== \n");
	fflush(outfile);
	Process::environment.globals["MP2 TOTAL ENERGY"] = Emp2;
	Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
	Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
	Process::environment.globals["SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;
	Process::environment.globals["SCS-MP2-VDW TOTAL ENERGY"] = Escsmp2vdw;
	Process::environment.globals["SOS-PI-MP2 TOTAL ENERGY"] = Esospimp2;

	Process::environment.globals["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
	Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
	Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
	Process::environment.globals["SCSN-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
	Process::environment.globals["SCS-MP2-VDW CORRELATION ENERGY"] = Escsmp2vdw - Escf;
	Process::environment.globals["SOS-PI-MP2 CORRELATION ENERGY"] = Esospimp2 - Escf;

        Process::environment.globals["MP2 OPPOSITE-SPIN ENERGY"] = Emp2AB;
        Process::environment.globals["MP2 SAME-SPIN ENERGY"] = Emp2AA+Emp2BB;

	ocepa_response_pdms();
	gfock();
	idp();
	mograd();
        if (rms_wog > tol_grad) occ_iterations();
        else {
           orbs_already_opt = 1;
	   fprintf(outfile,"\n\tOrbitals are already optimized, switching to the canonical CEPA computation... \n");
	   fflush(outfile);
           cepa_iterations();
        }
	
        if (rms_wog <= tol_grad && fabs(DE) >= tol_Eod) {
           orbs_already_opt = 1;
	   if (conver == 1) fprintf(outfile,"\n\tOrbitals are optimized now.\n");
	   else if (conver == 0) { 
                    fprintf(outfile,"\n\tMAX MOGRAD did NOT converged, but RMS MOGRAD converged!!!\n");
	            fprintf(outfile,"\tI will consider the present orbitals as optimized.\n");
           }
	   fprintf(outfile,"\tSwitching to the standard CEPA computation... \n");
	   fflush(outfile);
           ref_energy();
	   cepa_energy();
           Ecepa_old = EcepaL;
           cepa_iterations();
           if (dertype == "FIRST") {
               ocepa_response_pdms();
	       gfock();
           }
        } 

  if (conver == 1) {
        ref_energy();
	cepa_energy();
        if (orbs_already_opt == 1) EcepaL = Ecepa;

        // EKT
        if (ekt_ip_ == "TRUE") { 
            if (orbs_already_opt == 1) {
	        ocepa_response_pdms();
	        gfock();
            }
            gfock_diag();
            ekt_ip();
        }

	fprintf(outfile,"\n");
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\t================ OCEPA FINAL RESULTS ========================================= \n");
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	fprintf(outfile,"\tSCS-OCEPA(0) Total Energy (a.u.)   : %20.14f\n", Escscepa);
	fprintf(outfile,"\tSOS-OCEPA(0) Total Energy (a.u.)   : %20.14f\n", Esoscepa);
	fprintf(outfile,"\tOCEPA(0) Correlation Energy (a.u.) : %20.14f\n", EcepaL-Escf);
	fprintf(outfile,"\tEocepa - Eref (a.u.)               : %20.14f\n", EcepaL-Eref);
	fprintf(outfile,"\tOCEPA(0) Total Energy (a.u.)       : %20.14f\n", EcepaL);
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\n");
	fflush(outfile);
	
	// Set the global variables with the energies
	Process::environment.globals["OCEPA(0) TOTAL ENERGY"] = EcepaL;
	Process::environment.globals["SCS-OCEPA(0) TOTAL ENERGY"] =  Escscepa;
	Process::environment.globals["SOS-OCEPA(0) TOTAL ENERGY"] =  Esoscepa;
	Process::environment.globals["CURRENT ENERGY"] = EcepaL;
	Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
	Process::environment.globals["CURRENT CORRELATION ENERGY"] = EcepaL-Escf;

	Process::environment.globals["OCEPA(0) CORRELATION ENERGY"] = EcepaL - Escf;
	Process::environment.globals["SCS-OCEPA(0) CORRELATION ENERGY"] =  Escscepa - Escf;
	Process::environment.globals["SOS-OCEPA(0) CORRELATION ENERGY"] =  Esoscepa - Escf;

        // if scs on	
	if (do_scs == "TRUE") {
	    Process::environment.globals["CURRENT ENERGY"] = Escscepa;
	    Process::environment.globals["CURRENT CORRELATION ENERGY"] = Escscepa - Escf;

	}
    
        // else if sos on	
	else if (do_sos == "TRUE") {
	         Process::environment.globals["CURRENT ENERGY"] = Esoscepa;
	         Process::environment.globals["CURRENT CORRELATION ENERGY"] = Esoscepa - Escf;
	}

	if (natorb == "TRUE") nbo();
	if (occ_orb_energy == "TRUE") semi_canonic(); 
	
        // Compute Analytic Gradients
        if (dertype == "FIRST") {
            time4grad = 1;
	    fprintf(outfile,"\tAnalytic gradient computation is starting...\n");
	    fflush(outfile);
            coord_grad();
	    fprintf(outfile,"\tNecessary information has been sent to DERIV, which will take care of the rest.\n");
	    fflush(outfile);
        }

  }// end if (conver == 1)
}// end ocepa_manager 


//======================================================================
//             CEPA Manager
//======================================================================             
void OCCWave::cepa_manager()
{
        timer_on("trans_ints");
	if (reference_ == "RESTRICTED") trans_ints_rhf();  
	else if (reference_ == "UNRESTRICTED") trans_ints_uhf();  
        timer_off("trans_ints");
        Eref = Escf;
        timer_on("T2(1)");
	ocepa_t2_1st_sc();
        timer_off("T2(1)");
        timer_on("MP2 Energy");
	ocepa_mp2_energy();
        timer_off("MP2 Energy");
        Ecepa = Emp2;
	Ecepa_old = Emp2;

	fprintf(outfile,"\n"); 
	fprintf(outfile,"\tComputing MP2 energy using SCF MOs (Canonical MP2)... \n"); 
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	fprintf(outfile,"\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	fprintf(outfile,"\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	fprintf(outfile,"\tSCS-MP2 Total Energy (a.u.)        : %20.14f\n", Escsmp2);
	fprintf(outfile,"\tSOS-MP2 Total Energy (a.u.)        : %20.14f\n", Esosmp2);
	fprintf(outfile,"\tSCSN-MP2 Total Energy (a.u.)       : %20.14f\n", Escsnmp2);
	fprintf(outfile,"\tSCS-MP2-VDW Total Energy (a.u.)    : %20.14f\n", Escsmp2vdw);
	fprintf(outfile,"\tSOS-PI-MP2 Total Energy (a.u.)     : %20.14f\n", Esospimp2);
	fprintf(outfile,"\tMP2 Correlation Energy (a.u.)      : %20.14f\n", Ecorr);
	fprintf(outfile,"\tMP2 Total Energy (a.u.)            : %20.14f\n", Emp2);
	fprintf(outfile,"\t============================================================================== \n");
	fflush(outfile);
	Process::environment.globals["MP2 TOTAL ENERGY"] = Emp2;
	Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
	Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
	Process::environment.globals["SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;
	Process::environment.globals["SCS-MP2-VDW TOTAL ENERGY"] = Escsmp2vdw;
	Process::environment.globals["SOS-PI-MP2 TOTAL ENERGY"] = Esospimp2;

	Process::environment.globals["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
	Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
	Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
	Process::environment.globals["SCSN-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
	Process::environment.globals["SCS-MP2-VDW CORRELATION ENERGY"] = Escsmp2vdw - Escf;
	Process::environment.globals["SOS-PI-MP2 CORRELATION ENERGY"] = Esospimp2 - Escf;

        Process::environment.globals["MP2 OPPOSITE-SPIN ENERGY"] = Emp2AB;
        Process::environment.globals["MP2 SAME-SPIN ENERGY"] = Emp2AA+Emp2BB;

        // Perform CEPA iterations
        cepa_iterations();

	fprintf(outfile,"\n");
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\t================ CEPA FINAL RESULTS ========================================== \n");
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	fprintf(outfile,"\tCEPA(0) Correlation Energy (a.u.)  : %20.14f\n", Ecorr);
	fprintf(outfile,"\tCEPA(0) Total Energy (a.u.)        : %20.14f\n", Ecepa);
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\n");
	fflush(outfile);
	
	// Set the global variables with the energies
	Process::environment.globals["CEPA(0) TOTAL ENERGY"] = Ecepa;
	Process::environment.globals["CURRENT ENERGY"] = Ecepa;
	Process::environment.globals["CURRENT REFERENCE ENERGY"] = Eref;
	Process::environment.globals["CURRENT CORRELATION ENERGY"] = Ecorr;
        //EcepaL = Ecepa;
        
        // Compute Analytic Gradients
        if (dertype == "FIRST" || ekt_ip_ == "TRUE") {
            time4grad = 1;
	    fprintf(outfile,"\tAnalytic gradient computation is starting...\n");
            fprintf(outfile,"\tComputing response density matrices...\n");
            fflush(outfile);
	    ocepa_response_pdms();
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
            if (ekt_ip_ == "TRUE") ekt_ip();
	    fprintf(outfile,"\tNecessary information has been sent to DERIV, which will take care of the rest.\n");
	    fflush(outfile);
        }
}// end cepa_manager 


//======================================================================
//             OMP2.5 Manager
//======================================================================             
void OCCWave::omp2_5_manager()
{
	mo_optimized = 0;
	orbs_already_opt = 0;
	orbs_already_sc = 0;
        time4grad = 0;// means i will not compute the gradient
        timer_on("trans_ints");
	if (reference_ == "RESTRICTED") trans_ints_rhf();  
	else if (reference_ == "UNRESTRICTED") trans_ints_uhf();  
        timer_off("trans_ints");
        timer_on("REF Energy");
	ref_energy();
        timer_off("REF Energy");
        timer_on("T2(1)");
	omp3_t2_1st_sc();
        timer_off("T2(1)");
        timer_on("MP2 Energy");
	omp3_mp2_energy();
        timer_off("MP2 Energy");

	fprintf(outfile,"\n"); 
	fprintf(outfile,"\tComputing MP2 energy using SCF MOs (Canonical MP2)... \n"); 
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	fprintf(outfile,"\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	fprintf(outfile,"\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	fprintf(outfile,"\tSCS-MP2 Total Energy (a.u.)        : %20.14f\n", Escsmp2);
	fprintf(outfile,"\tSOS-MP2 Total Energy (a.u.)        : %20.14f\n", Esosmp2);
	fprintf(outfile,"\tSCSN-MP2 Total Energy (a.u.)       : %20.14f\n", Escsnmp2);
	fprintf(outfile,"\tSCS-MP2-VDW Total Energy (a.u.)    : %20.14f\n", Escsmp2vdw);
	fprintf(outfile,"\tSOS-PI-MP2 Total Energy (a.u.)     : %20.14f\n", Esospimp2);
	fprintf(outfile,"\tMP2 Correlation Energy (a.u.)      : %20.14f\n", Ecorr);
	fprintf(outfile,"\tMP2 Total Energy (a.u.)            : %20.14f\n", Emp2);
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\n"); 
	fflush(outfile);
	Process::environment.globals["MP2 TOTAL ENERGY"] = Emp2;
	Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
	Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
	Process::environment.globals["SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;
	Process::environment.globals["SCS-MP2-VDW TOTAL ENERGY"] = Escsmp2vdw;
	Process::environment.globals["SOS-PI-MP2 TOTAL ENERGY"] = Esospimp2;

	Process::environment.globals["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
	Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
	Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
	Process::environment.globals["SCSN-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
	Process::environment.globals["SCS-MP2-VDW CORRELATION ENERGY"] = Escsmp2vdw - Escf;
	Process::environment.globals["SOS-PI-MP2 CORRELATION ENERGY"] = Esospimp2 - Escf;

        Process::environment.globals["MP2 OPPOSITE-SPIN ENERGY"] = Emp2AB;
        Process::environment.globals["MP2 SAME-SPIN ENERGY"] = Emp2AA+Emp2BB;

        timer_on("T2(2)");
	t2_2nd_sc();
        timer_off("T2(2)");
        timer_on("MP3 Energy");
	mp3_energy();
        timer_off("MP3 Energy");
	Emp3L=Emp3;
        EcorrL=Emp3L-Escf;
	Emp3L_old=Emp3;
        if (ip_poles == "TRUE") omp3_ip_poles();
	
	fprintf(outfile,"\n"); 
	fprintf(outfile,"\tComputing MP2.5 energy using SCF MOs (Canonical MP2.5)... \n"); 
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp3AA);
	fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp3AB);
	fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp3BB);
	fprintf(outfile,"\t0.5 Energy Correction (a.u.)       : %20.14f\n", Emp3-Emp2);
	fprintf(outfile,"\tMP2.5 Correlation Energy (a.u.)    : %20.14f\n", Ecorr);
	fprintf(outfile,"\tMP2.5 Total Energy (a.u.)          : %20.14f\n", Emp3);
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\n"); 
	fflush(outfile);
	Process::environment.globals["MP2.5 TOTAL ENERGY"] = Emp3;

	omp3_response_pdms();
	gfock();
        //ccl_energy();
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
	   fprintf(outfile,"\tSwitching to the standard MP2.5 computation after semicanonicalization of the MOs... \n");
	   fflush(outfile);
	   semi_canonic();
	   if (reference_ == "RESTRICTED") trans_ints_rhf();  
	   else if (reference_ == "UNRESTRICTED") trans_ints_uhf();  
	   omp3_t2_1st_sc();
	   t2_2nd_sc();
           conver = 1;
           if (dertype == "FIRST") {
               omp3_response_pdms();
	       gfock();
           }
        }     

  if (conver == 1) {
        ref_energy();
	omp3_mp2_energy();
	mp3_energy();
        if (orbs_already_opt == 1) Emp3L = Emp3;

        if (ip_poles == "TRUE") {
	   if (orbs_already_sc == 0) {
               semi_canonic();
	       if (reference_ == "RESTRICTED") trans_ints_rhf();  
	       else if (reference_ == "UNRESTRICTED") trans_ints_uhf();  
	       omp3_t2_1st_sc();
	       t2_2nd_sc();
           }
           omp3_ip_poles();
        }

        // EKT
        if (ekt_ip_ == "TRUE") { 
            if (orbs_already_sc == 1) {
	        omp3_response_pdms();
	        gfock();
            }
            gfock_diag();
            ekt_ip();
        }

        fprintf(outfile,"\n"); 
	fprintf(outfile,"\tComputing MP2 energy using optimized MOs... \n"); 
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	fprintf(outfile,"\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	fprintf(outfile,"\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	fprintf(outfile,"\tSCS-MP2 Total Energy (a.u.)        : %20.14f\n", Escsmp2);
	fprintf(outfile,"\tSOS-MP2 Total Energy (a.u.)        : %20.14f\n", Esosmp2);
	fprintf(outfile,"\tSCSN-MP2 Total Energy (a.u.)       : %20.14f\n", Escsnmp2);
	fprintf(outfile,"\tSCS-MP2-VDW Total Energy (a.u.)    : %20.14f\n", Escsmp2vdw);
	fprintf(outfile,"\tSOS-PI-MP2 Total Energy (a.u.)     : %20.14f\n", Esospimp2);
	fprintf(outfile,"\tMP2 Correlation Energy (a.u.)      : %20.14f\n", Ecorr);
	fprintf(outfile,"\tMP2 Total Energy (a.u.)            : %20.14f\n", Emp2);
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\n"); 
	fflush(outfile);

	fprintf(outfile,"\n"); 
	fprintf(outfile,"\tComputing MP2.5 energy using optimized MOs... \n"); 
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp3AA);
	fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp3AB);
	fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp3BB);
	fprintf(outfile,"\t0.5 Energy Correction (a.u.)       : %20.14f\n", Emp3-Emp2);
	fprintf(outfile,"\tMP2.5 Correlation Energy (a.u.)    : %20.14f\n", Ecorr);
	fprintf(outfile,"\tMP2.5 Total Energy (a.u.)          : %20.14f\n", Emp3);
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\n"); 
	fflush(outfile);


	fprintf(outfile,"\n");
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\t================ OMP2.5 FINAL RESULTS ======================================== \n");
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	fprintf(outfile,"\tOMP2.5 Correlation Energy (a.u.)   : %20.14f\n", Emp3L-Escf);
	fprintf(outfile,"\tEomp2.5 - Eref (a.u.)              : %20.14f\n", Emp3L-Eref);
	fprintf(outfile,"\tOMP2.5 Total Energy (a.u.)         : %20.14f\n", Emp3L);
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\n");
	fflush(outfile);
	
	// Set the global variables with the energies
	Process::environment.globals["OMP2.5 TOTAL ENERGY"] = Emp3L;
	Process::environment.globals["OMP2.5 CORRELATION ENERGY"] = Emp3L - Escf;
	Process::environment.globals["CURRENT ENERGY"] = Emp3L;
	Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
	Process::environment.globals["CURRENT CORRELATION ENERGY"] = Emp3L-Escf;

	if (natorb == "TRUE") nbo();
	if (occ_orb_energy == "TRUE") semi_canonic(); 

        // Compute Analytic Gradients
        if (dertype == "FIRST") {
            time4grad = 1;
	    fprintf(outfile,"\tAnalytic gradient computation is starting...\n");
	    fflush(outfile);
            coord_grad();
	    fprintf(outfile,"\tNecessary information has been sent to DERIV, which will take care of the rest.\n");
	    fflush(outfile);
        }

  }// end if (conver == 1)
}// end omp2.5_manager 


//======================================================================
//             MP2.5 Manager
//======================================================================             
void OCCWave::mp2_5_manager()
{
        time4grad = 0;// means i will not compute the gradient
        timer_on("trans_ints");
	if (reference_ == "RESTRICTED") trans_ints_rhf();  
	else if (reference_ == "UNRESTRICTED") trans_ints_uhf();  
        timer_off("trans_ints");
        Eref = Escf;
        timer_on("T2(1)");
	omp3_t2_1st_sc();
        timer_off("T2(1)");
        timer_on("MP2 Energy");
	omp3_mp2_energy();
        timer_off("MP2 Energy");

	fprintf(outfile,"\n"); 
	fprintf(outfile,"\tComputing MP2 energy using SCF MOs (Canonical MP2)... \n"); 
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
	fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
	fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
	fprintf(outfile,"\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA+Escsmp2BB);
	fprintf(outfile,"\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
	fprintf(outfile,"\tSCS-MP2 Total Energy (a.u.)        : %20.14f\n", Escsmp2);
	fprintf(outfile,"\tSOS-MP2 Total Energy (a.u.)        : %20.14f\n", Esosmp2);
	fprintf(outfile,"\tSCSN-MP2 Total Energy (a.u.)       : %20.14f\n", Escsnmp2);
	fprintf(outfile,"\tSCS-MP2-VDW Total Energy (a.u.)    : %20.14f\n", Escsmp2vdw);
	fprintf(outfile,"\tSOS-PI-MP2 Total Energy (a.u.)     : %20.14f\n", Esospimp2);
	fprintf(outfile,"\tMP2 Correlation Energy (a.u.)      : %20.14f\n", Ecorr);
	fprintf(outfile,"\tMP2 Total Energy (a.u.)            : %20.14f\n", Emp2);
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\n"); 
	fflush(outfile);
	Process::environment.globals["MP2 TOTAL ENERGY"] = Emp2;
	Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
	Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
	Process::environment.globals["SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;
	Process::environment.globals["SCS-MP2-VDW TOTAL ENERGY"] = Escsmp2vdw;
	Process::environment.globals["SOS-PI-MP2 TOTAL ENERGY"] = Esospimp2;

	Process::environment.globals["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
	Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
	Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
	Process::environment.globals["SCSN-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
	Process::environment.globals["SCS-MP2-VDW CORRELATION ENERGY"] = Escsmp2vdw - Escf;
	Process::environment.globals["SOS-PI-MP2 CORRELATION ENERGY"] = Esospimp2 - Escf;

        Process::environment.globals["MP2 OPPOSITE-SPIN ENERGY"] = Emp2AB;
        Process::environment.globals["MP2 SAME-SPIN ENERGY"] = Emp2AA+Emp2BB;

        timer_on("T2(2)");
	t2_2nd_sc();
        timer_off("T2(2)");
        timer_on("MP3 Energy");
	mp3_energy();
        timer_off("MP3 Energy");
	Emp3L=Emp3;
        EcorrL=Emp3L-Escf;
	Emp3L_old=Emp3;
        if (ip_poles == "TRUE") omp3_ip_poles();
	
	fprintf(outfile,"\n"); 
	fprintf(outfile,"\tComputing MP2.5 energy using SCF MOs (Canonical MP2.5)... \n"); 
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
	fprintf(outfile,"\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
	fprintf(outfile,"\tREF Energy (a.u.)                  : %20.14f\n", Eref);
	fprintf(outfile,"\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp3AA);
	fprintf(outfile,"\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp3AB);
	fprintf(outfile,"\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp3BB);
	fprintf(outfile,"\t0.5 Energy Correction (a.u.)       : %20.14f\n", Emp3-Emp2);
	fprintf(outfile,"\tMP2.5 Correlation Energy (a.u.)    : %20.14f\n", Ecorr);
	fprintf(outfile,"\tMP2.5 Total Energy (a.u.)          : %20.14f\n", Emp3);
	fprintf(outfile,"\t============================================================================== \n");
	fprintf(outfile,"\n"); 
	fflush(outfile);
	Process::environment.globals["MP2.5 TOTAL ENERGY"] = Emp3;
	Process::environment.globals["MP2.5 CORRELATION ENERGY"] = Emp3 - Escf;
	Process::environment.globals["CURRENT ENERGY"] = Emp3L;
	Process::environment.globals["CURRENT REFERENCE ENERGY"] = Eref;
	Process::environment.globals["CURRENT CORRELATION ENERGY"] = Emp3L-Escf;

        // Compute Analytic Gradients
        if (dertype == "FIRST" || ekt_ip_ == "TRUE") {
            time4grad = 1;
	    fprintf(outfile,"\tAnalytic gradient computation is starting...\n");
            fprintf(outfile,"\tComputing response density matrices...\n");
            fflush(outfile);
	    omp3_response_pdms();
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
            if (ekt_ip_ == "TRUE") ekt_ip();
	    fprintf(outfile,"\tNecessary information has been sent to DERIV, which will take care of the rest.\n");
	    fflush(outfile);
        }

}// end omp2.5_manager 


}} // End Namespaces


