/*! \file
    \ingroup INPUT
    \brief Enter brief description of file here 
*/
#define EXTERN
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <psifiles.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <liboptions/liboptions.h>
#include <libchkpt/chkpt.h>
#include "input.h"
#include <physconst.h>
#include "global.h"
#include "defines.h"

namespace psi { namespace input {

extern Options options;

//void start_io(int argc, char *argv[])
void start_io()
{
  int i, errcod;
//int num_extra_args = 0;
//char **extra_args;
//extra_args = (char **) malloc(argc*sizeof(char *));

  keep_chkpt = 0;  
  read_chkpt = 0;
  chkpt_mos = 0;
  chkpt_geom = 0; 
  dont_project_mos = 0; 
  geomdat_geom = 0;
  save_oldcalc = 0;
  overwrite_output = 0;
  no_comshift = 0;
  no_reorient = 0;
			  
//for (i=1; i<argc; i++) {
    /*--- read MOs from checkpoint file and project onto new basis ---*/
    //if (strcmp(argv[i], "--chkptmos") == 0) {
    //read_chkpt = 1;
    //chkpt_mos = 1;
    //}
    read_chkpt = chkpt_mos = options.get_int("CHKPT_MOS");
    /*--- read MOs from checkpoint file and save to a separate file ---*/
    //else if (strcmp(argv[i], "--savemos") == 0) {
    //  read_chkpt = 1;
    //  save_oldcalc = 1;
    //}
    read_chkpt = save_oldcalc = options.get_int("SAVE_MOS");
    /*--- don't project MOs but simply keep them ---*/
    //else if (strcmp(argv[i], "--noproject") == 0) {
    //  dont_project_mos = 1;
    //}
    dont_project_mos = options.get_int("NOPROJECT");
    /*--- read geometry from checkpoint file (in findif calculations) ---*/
    if (options.get_int("CHKPT_GEOM") /* strcmp(argv[i], "--chkptgeom") == 0 */ ) {
      read_chkpt = 0;
      chkpt_geom = 1;
      /* preserve the information about the original reference frame
       * so that properties can be rotated back later by other codes */
      keep_ref_frame = 1;     
      print_lvl = 0;
      cartOn = 1;
      overwrite_output = 0;
    }
	//else if (strcmp(argv[i], "--keepoutput") == 0) {
	//	/* Don't overwrite the output file. This was added for the new psirb driver module. -Jet 30 Jul 07 */
	//	overwrite_output = 0;
	//}
    overwrite_output = !options.get_int("KEEP_OUTPUT");
    /*--- read geometry from geom.dat file (in findif calculations) ---*/
    // JET -> NOT USED IN PSI.DAT
    //if (options.get_int("GEOMDAT") /* strcmp(argv[i], "--geomdat") == 0 */) {
    //  geomdat_geom = 1;
    //  geomdat_entry = atoi(argv[i+1]);  i++;
    //  keep_ref_frame = 1;
    //  cartOn = 1;
    //  overwrite_output = 0;
    //}
    //else if (strcmp(argv[i], "--nocomshift") == 0) {
    //  no_comshift = 1;
    //}
    no_comshift = options.get_int("NO_COM_SHIFT");
    //else if (strcmp(argv[i], "--noreorient") == 0) {
    //  no_reorient = 1;
    //}
    no_reorient = options.get_int("NO_REORIENT");
    //else if (strcmp(argv[i], "--keepchkpt") == 0) {
    //  keep_chkpt = 1;
    //}
    keep_chkpt = options.get_int("KEEP_CHKPT");
    //else {
    //  extra_args[num_extra_args++] = argv[i];
    //}
//}
  
  //errcod = psi_start(&infile,&outfile,&psi_file_prefix,num_extra_args, extra_args, overwrite_output);
  //if (errcod != PSI_RETURN_SUCCESS)
  //  abort();
  //ip_cwk_add(":INPUT");
  tstart();
  //psio_init(); psio_ipv1_config();

//free(extra_args);

  return;
}

void stop_io()
{
  tstop();
 // psio_done();
 // psi_stop(infile,outfile,psi_file_prefix);
}

}} // namespace psi::input

