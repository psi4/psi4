#include <libmoinfo/libmoinfo.h>
#include <liboptions/liboptions.h>
#include "mrcc.h"
#include "matrix.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{
    extern MOInfo *moinfo;

using namespace std;

void CCMRCC::print_mrccsd_energy(int cycle)
{
  delta_energy = current_energy-old_energy;
  if(cycle==0){
    print_method("\tMultireference Coupled Cluster\n\t\tUsing the DPD Library");
    fprintf(outfile,"\n  ------------------------------------------------------------------------------");
    fprintf(outfile,"\n  @CC Cycle      Energy          Delta E    ||DeltaT1|| ||DeltaT2|| Timing  DIIS");
    fprintf(outfile,"\n  @CC           (Hartree)       (Hartree)                           (Sec)");
    fprintf(outfile,"\n  ------------------------------------------------------------------------------");
  }
  if(cycle>=0){
    fprintf(outfile,"\n  @CC %3d  %18.12f  %11.4e   %8.3e   %8.3e %7.0f",cycle,current_energy,delta_energy,delta_t1_amps,delta_t2_amps,total_time);

    
    if((fabs(delta_energy) < options_.get_double("E_CONVERGENCE")) && (cycle!=0)){
      char star = (options_.get_str("CORR_WFN") == "CCSD") ? '*' : ' ';
      fprintf(outfile,"\n  ------------------------------------------------------------------------------");
      fprintf(outfile,"\n\n%6c%1c Mk-MRCCSD total energy      = %20.12f\n",' ',star,current_energy);
    }
  }else if(cycle==-1){
    char star = ' ';
    if(options_.get_str("CORR_WFN") == "CCSD")
      star = '*';
    fprintf(outfile,"\n\n%6c%1c Mk-MRCCSD total energy      = %20.12f\n",' ',star,current_energy);
    print_eigensystem(moinfo->get_nrefs(),Heff,right_eigenvector);
  }

  fflush(outfile);
}

}} /* End Namespaces */
