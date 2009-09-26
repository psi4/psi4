/*! \file
    \ingroup INPUT
    \brief Enter brief description of file here 
*/
#define EXTERN
#include <cstdio>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <liboptions/liboptions.h>
#include <libchkpt/chkpt.h>
#include <cstring>
#include "input.h"
#include <physconst.h>
#include "global.h"
#include "defines.h"

namespace psi { namespace input {

void parsing(Options & options)
{
   int errcod, i;
   char *guess, tmp_label[80];

     /* Same as --noreorient */
     no_reorient = options.get_bool("NO_REORIENT");

     /* same as --chkptmos */
     chkpt_mos = options.get_bool("CHKPT_MOS");
     if (chkpt_mos) read_chkpt = 1;

     label = options.get_cstr("LABEL");

     /*Flag if the user wants to see that his/her basis set is normalized*/
     shownorm = options.get_bool("SHOWNORM");

     /*Flag if the user wants to have his/her basis set normalized*/
     normalize_contractions = options.get_bool("NORMALIZE");

     puream = options.get_bool("PUREAM");

     expert = options.get_bool("EXPERT");

     /* read print_lvl if chkpt_geom too */
	 print_lvl = options.get_int("PRINT");

     /* allow the user to specify subgroup=C1 for entire findif calc -
        hope this doesn't mess anything up (RAK 9-04) */
	 subgroup = options.get_str("SUBGROUP");
	 
     /*------------------------------------------
       Parse Some boolean information from input
      ------------------------------------------*/

     /* If not reading geometry from chkpt file, we will
	need some options specified in input.dat */
     if (chkpt_geom == 0) {
	   print_lvl = options.get_int("PRINT");

	 /*Get the subgroup label */
	 subgroup = options.get_str("SUBGROUP");
	 
	 /*Get the unique axis*/
	 unique_axis = options.get_str("UNIQUE_AXIS");

     nfragments = options.get_int("NFRAGMENTS");

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

       units = options.get_str("UNITS");

       if ((units == "BOHR") || (units == "AU"))
         conv_factor = 1.0;
       else if (units == "ANGSTROMS" || units == "ANGSTROM")
         conv_factor = 1.0 / _bohr2angstroms;
	   
	   /*Set reference frame to be the frame of the input geometry*/
	   keep_ref_frame = options.get_bool("KEEP_REF_FRAME");
	 }
     }


     /* Check if need to freeze core */
     frozen_core = NULL;
     errcod = ip_string("FREEZE_CORE",&frozen_core,0);
     if (frozen_core == NULL)
       frozen_core = strdup("FALSE");

     /* Check if need to freeze virtuals */
     frozen_virt = NULL;
     errcod = ip_string("FREEZE_VIRT",&frozen_virt,0);
     if (frozen_virt == NULL)
       frozen_virt = strdup("FALSE");

     return;
}

}} // namespace psi::input
