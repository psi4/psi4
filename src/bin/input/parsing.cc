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

void parsing()
{
  int i;
  char tmp_label[80];

  /*--- read MOs from checkpoint file and project onto new basis ---*/
  chkpt_mos   = options["CHKPT_MOS"].to_integer();

  /*--- read geometry from checkpoint file (in findif calculations) ---*/
  // this keyword may be obseleted in psi4
  chkpt_geom = options["CHKPT_GEOM"].to_integer();

  /*--- don't project MOs but simply keep them ---*/
  dont_project_mos = options["NOPROJECT"].to_integer();

  geomdat_geom = options["GEOMDAT"].to_integer();
  save_oldcalc = 0;
  overwrite_output = 0;
  no_comshift = options["NO_COM_SHIFT"].to_integer();
  no_reorient = 0;
			  
  // will not survive psi4
  //if (geomdat_geom) 
  //  geomdat_entry = atoi(argv[i+1]);  i++;

  no_reorient = options["NO_REORIENT"].to_integer();

  read_chkpt = 0;

    /*--- read MOs from checkpoint file and save to a separate file ---*/
   if (options["SAVE_MOS"].to_integer()) {
      read_chkpt = 1;
      save_oldcalc = 1;
   }

    /*--- read geometry from checkpoint file (in findif calculations) ---*/
   if (chkpt_geom) {
      read_chkpt = 0;
      chkpt_geom = 1;
      // preserve the information about the original reference frame
      // so that properties can be rotated back later by other codes
      keep_ref_frame = 1;
      print_lvl = 0;
      cartOn = 1;
      overwrite_output = 0;
    }

    if (chkpt_mos) read_chkpt = 1;

	/* Don't overwrite the output file */
	if (options["KEEP_OUTPUT"].to_integer())
      overwrite_output = 0;

    /* Keep the chkpt file. */
   keep_chkpt = options["KEEP_CHKPT"].to_integer();  

   label       = options["LABEL"].to_string();
   shownorm    = options["SHOWNORM"].to_integer();
   puream      = options["PUREAM"].to_integer();
   expert      = options["EXPERT"].to_integer();
   print_lvl   = options["PRINT"].to_integer();
   subgroup    = options["SUBGROUP"].to_string();
   unique_axis = options["UNIQUE_AXIS"].to_string();
   nfragments  = options["NFRAGMENTS"].to_integer();
   keep_ref_frame = options["KEEP_REF_FRAME"].to_integer();
   normalize_contractions = options["NORMALIZE"].to_integer();

   frozen_core = options["FREEZE_CORE"].to_string();
   nfzv = options["FREEZE_VIRT"].to_integer();

   if (chkpt_mos) read_chkpt = 1;

   // Not sure how to handle these using Options, yet. Options can handle
   // arrays of array. But how is it supposed to know that ZMAT and GEOMETRY
   // can exist as ZMAT%d and GEOMETRY%d ???
   
  // Find where geometry is present in input file
   if (chkpt_geom == 0) {
     /* check to make sure all fragment geometries are present */
     if (geomdat_geom == 0) {
     /*No default for these two unless running a findif procedure*/
       if (ip_exist(const_cast<char*>("ZMAT"),0) == 1) {
         cartOn = 0;
         for (i=1;i<nfragments;++i) {
           sprintf(tmp_label,"ZMAT%d",i+1);
           if (ip_exist(tmp_label,0) == 0)
             punt("input cannot find all the needed fragment structures!");
         }
       }
       else if (ip_exist(const_cast<char*>("GEOMETRY"),0) == 1) {
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
}

}} // namespace psi::input
