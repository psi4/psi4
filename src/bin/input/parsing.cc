/*! \file 
    \ingroup (INPUT)
    \brief Enter brief description of file here 
*/
#define EXTERN
#include <stdio.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <string.h>
#include "input.h"
#include <physconst.h>
#include "global.h"
#include "defines.h"

namespace psi { namespace input {

void parsing()
{
   int errcod;
   char *guess;

     /* Same as --noreorient */
     errcod = ip_boolean("NO_REORIENT",&no_reorient,0);

     /* same as --chkptmos */
     errcod = ip_boolean("CHKPT_MOS",&chkpt_mos,0);
     if (chkpt_mos) read_chkpt = 1;

     errcod = ip_string("LABEL", &label,0);
     if(errcod != IPE_OK)
       label = "Default PSI label";

     /*Default = false*/
     /*Flag if the user wants to see that his/her basis set is normalized*/
     shownorm = 0;
     errcod = ip_boolean("SHOWNORM",&shownorm,0);

     /*Default = true*/
     /*Flag if the user wants to have his/her basis set normalized*/
     normalize_contractions = 1;
     errcod = ip_boolean("NORMALIZE",&normalize_contractions,0);

     /*Default = false*/
     puream = 0;
     errcod = ip_boolean("PUREAM",&puream,0);

     /* Safe behavior is on by default */
     expert = 0;
     errcod = ip_boolean("EXPERT",&expert,0);

     /* read print_lvl if chkpt_geom too */
	 errcod = ip_data("PRINT","%d",&print_lvl,0);

     /* allow the user to specify subgroup=C1 for entire findif calc -
        hope this doesn't mess anything up (RAK 9-04) */
	 subgroup = NULL;
	 errcod = ip_string("SUBGROUP", &subgroup,0);
	 if (subgroup != NULL && strlen(subgroup) != 0)
	     if ( strcmp(subgroup,"C1") )
		   subgroup = NULL;
	 
     /*------------------------------------------
       Parse Some boolean information from input
      ------------------------------------------*/

     /* If not reading geometry from chkpt file, we will
	need some options specified in input.dat */
     if (chkpt_geom == 0) {
	 /*Default = 1*/
	 print_lvl = 1;
	 errcod = ip_data("PRINT","%d",&print_lvl,0);

	 /*Get the subgroup label */
	 subgroup = NULL;
	 errcod = ip_string("SUBGROUP", &subgroup,0);
	 if (subgroup != NULL && strlen(subgroup) != 0)
	     if (strcmp(subgroup,"C1") && strcmp(subgroup,"C2") && strcmp(subgroup,"CS") && strcmp(subgroup,"CI") &&
		 strcmp(subgroup,"C2V") && strcmp(subgroup,"C2H") && strcmp(subgroup,"D2"))
		 subgroup = NULL;
	 
	 /*Get the unique axis*/
	 unique_axis = NULL;
	 errcod = ip_string("UNIQUE_AXIS", &unique_axis, 0);
	 if (unique_axis != NULL && strlen(unique_axis) != 0)
	     if (strcmp(unique_axis,"X") && strcmp(unique_axis,"Y") && strcmp(unique_axis,"Z"))
		 unique_axis = NULL;

	 if (geomdat_geom == 0) {
	   /*No default for these two unless running a findif procedure*/
	   if (ip_exist("ZMAT",0) == 1)
	     cartOn = 0;
	   else if (ip_exist("GEOMETRY",0) == 1)
	     cartOn = 1;
	   else
	     punt("Both ZMAT and GEOMETRY are missing!");
	   
	   /*Default = ANGSTROMS*/
	   units = strdup("ANGSTROMS");
	   errcod = ip_string("UNITS",&units,0);
	   if (!strcmp(units,"BOHR") || !strcmp(units,"AU"))
	     conv_factor = 1.0;
	   else if (!strcmp(units,"ANGSTROMS") || !strcmp(units,"ANGSTROM"))
	     conv_factor = 1.0 / _bohr2angstroms;
	   else
	     punt("Unrecognized UNITS");
	   
	   /*Set reference frame to be the frame of the input geometry*/
	   keep_ref_frame = 0;
	   errcod = ip_boolean("KEEP_REF_FRAME",&keep_ref_frame,0);
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
