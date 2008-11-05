/*! \file
    \ingroup INPUT
    \brief Enter brief description of file here 
*/
#define EXTERN
#include <cstdio>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <cstring>
#include "input.h"
#include <physconst.h>
#include "global.h"
#include "defines.h"
#include <liboptions/liboptions.h>

namespace psi { namespace input {

void parsing(Options & options)
{
  int i;
  char tmp_label[80];

  /*--- read MOs from checkpoint file and project onto new basis ---*/
  chkpt_mos   = options.get_bool_option("CHKPT_MOS");

  /*--- read geometry from checkpoint file (in findif calculations) ---*/
  // this keyword may be obseleted in psi4
  chkpt_geom = options.get_bool_option("CHKPT_GEOM");

  /*--- don't project MOs but simply keep them ---*/
  dont_project_mos = options.get_bool_option("NOPROJECT");

  geomdat_geom = options.get_bool_option("GEOMDAT");
  save_oldcalc = 0;
  overwrite_output = 0;
  no_comshift = options.get_bool_option("NO_COM_SHIFT");
  no_reorient = 0;
			  
  // will not survive psi4
  //if (geomdat_geom) 
  //  geomdat_entry = atoi(argv[i+1]);  i++;

  no_reorient = options.get_bool_option("NO_REORIENT");

  read_chkpt = 0;

    /*--- read MOs from checkpoint file and save to a separate file ---*/
   if (options.get_bool_option("SAVE_MOS")) {
      read_chkpt = 1;
      save_oldcalc = 1;
   }

    /*--- read geometry from checkpoint file (in findif calculations) ---*/
   if (chkpt_geom) {
      read_chkpt = 0;
      chkpt_geom = 1;
      /* preserve the information about the original reference frame
 *        * so that properties can be rotated back later by other codes */
      keep_ref_frame = 1;
      print_lvl = 0;
      cartOn = 1;
      overwrite_output = 0;
    }

    if (chkpt_mos) read_chkpt = 1;

	/* Don't overwrite the output file */
	if (options.get_bool_option("KEEP_OUTPUT"))
      overwrite_output = 0;

    /* Keep the chkpt file. */
   keep_chkpt = options.get_bool_option("KEEP_CHKPT");  

   label       = options.get_str_option("LABEL");
   shownorm    = options.get_bool_option("SHOWNORM");
   puream      = options.get_bool_option("PUREAM");
   expert      = options.get_bool_option("EXPERT");
   print_lvl   = options.get_int_option("PRINT");
   subgroup    = options.get_cstr_option("SUBGROUP");
   unique_axis = options.get_cstr_option("UNIQUE_AXIS");
   nfragments  = options.get_int_option("NFRAGMENTS");
   keep_ref_frame = options.get_bool_option("KEEP_REF_FRAME");
   normalize_contractions = options.get_bool_option("NORMALIZE");

   frozen_core = options.get_str_option("FREEZE_CORE");
   nfzv = options.get_int_option("FREEZE_VIRT");

   if (chkpt_mos) read_chkpt = 1;

	// Find where geometry is present in input file
   if (chkpt_geom == 0) {
     /* check to make sure all fragment geometries are present */
	 if (geomdat_geom == 0) {
	   /*No default for these two unless running a findif procedure*/
	   if (ip_exist("ZMAT",0) == 1) {
	     cartOn = 0;
         for (i=1;i<nfragments;++i) {
           sprintf(tmp_label,"ZMAT%d",i+1);
           if (ip_exist(tmp_label,0) == 0)
	         punt("input cannot find all the needed fragment structures!");
         }
       }
	   else if (ip_exist("GEOMETRY",0) == 1) {
	     cartOn = 1;
         for (i=1;i<nfragments;++i) {
           sprintf(tmp_label,"GEOMETRY%d",i+1);
           if (ip_exist(tmp_label,0) == 0)
	         punt("input cannot find all the needed fragment structures!");
         }
       }
	   else
	     punt("Both ZMAT and GEOMETRY are missing!");
	 }
   }
   return;
}

}} // namespace psi::input
