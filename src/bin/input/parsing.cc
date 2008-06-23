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

namespace psi { namespace input {

void parsing()
{
   int errcod, i;
   char tmp_label[80];

   no_reorient = options.get_bool_option("no_reorient");
   chkpt_mos   = options.get_bool_option("chkpt_mos");
   label       = options.get_str("label");
   shownorm    = options.get_bool_option("shownorm");
   puream      = options.get_bool_option("puream");
   expert      = options.get_bool_option("expert");
   print_lvl   = options.get_int_option("print");
   subgroup    = options.get_str_option_with_choices("SUBGROUP");
   unique_axis = options.get_str_option_with_choices("UNIQUE_AXIS");
   nfragments  = options.get_int_option("NFRAGMENTS");
   keep_ref_frame = options.get_bool_option("KEEP_REF_FRAME");
   normalize_contractions = options.get_bool_option("normalize");

   units       = options.get_str_option_with_choices("UNITS");
   if (!strcmp(units,"BOHR") || !strcmp(units,"AU"))
     conv_factor = 1.0;
   else if (!strcmp(units,"ANGSTROMS") || !strcmp(units,"ANGSTROM"))
     conv_factor = 1.0 / _bohr2angstroms;
   else
     punt("Unrecognized UNITS");

   frozen_core = options.get_bool_option("FREEZE_CORE");
   frozen_virt = options.get_bool_option("FREEZE_VIRT");

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
