/*! \file
    \ingroup OEPROP
    \brief Enter brief description of file here 
*/
#define EXTERN
#include <stdexcept>
#include "includes.h"
#include "globals.h"
#include "prototypes.h"

namespace psi { namespace oeprop {

void parsing()
{
  int i,errcod;
  double xmin,xmax,ymin,ymax,zmin,zmax;

  /* Set some defaults for certain wavefunctions */
  update_energy_with_MVD = 0;

//  errcod = ip_string("WFN", &wfn, 0);
  wfn = options.get_str("WFN");

//  if(errcod == IPE_OK) {
//  if (!strcmp(wfn, "CI") || !strcmp(wfn, "DETCI") ||
//      !strcmp(wfn, "CCSD") || !strcmp(wfn, "DETCAS") ||
//      !strcmp(wfn, "CASSCF") || !strcmp(wfn, "RASSCF") ||
//      !strcmp(wfn, "MP2") || !strcmp(wfn, "EOM_CCSD") ||
//      !strcmp(wfn, "CC2") || !strcmp(wfn, "EOM_CC2") ||
//      !strcmp(wfn, "CCSD_MVD"))  {
  if( wfn == "CI" || wfn == "DETCI" || 
      wfn == "CCSD" || wfn == "DETCAS" || 
      wfn == "CASSCF" || wfn == "RASSCF" || 
      wfn == "MP2" || wfn == "EOM_CCSD" || 
      wfn == "CC2" || wfn == "EOM_CC2" || 
      wfn == "CCSD_MVD") {
    read_opdm = 1;
    opdm_file = PSIF_MO_OPDM;
    corr = 0;
//    opdm_basis = (char *) malloc(3*sizeof(char));
//    strcpy(opdm_basis,"MO");
//    opdm_format = (char *) malloc(7*sizeof(char));
//    strcpy(opdm_format,"SQUARE");
  }
//  if (!strcmp(wfn, "SCF_MVD") || !strcmp(wfn, "CCSD_MVD"))
  if( wfn == "SCF_MVD" || wfn == "CCSD_MVD" )
    update_energy_with_MVD = 1;
//  }
//  else {
//    fprintf(outfile,"Incorrect wave function type!\n");
//    abort();
//  }

  transdens = 0;
//  errcod = ip_boolean("TRANSITION_DENSITY", &transdens, 0);
  transdens = options.get_int("TRANSITION_DENSITY");
  if (transdens)
    asymm_opdm = 1;
  else
    asymm_opdm = 0;

//  errcod = ip_string("REFERENCE", &ref, 0);
  ref = options.get_str("REFERENCE");
  
  /* Parsing section */

//  errcod = ip_boolean("READ_OPDM",&read_opdm,0);
  read_opdm = options.get_int("READ_OPDM");
  
  if (read_opdm) {
//    errcod = ip_data("OPDM_FILE","%d",&opdm_file,0);
    opdm_file = options.get_double("OPDM_FILE");
    if ((opdm_file >= PSIO_MAXUNIT) || (opdm_file <= 0))
//      punt("OPDM_FILE out of range");
      throw PsiException("OPDM_FILE out of range", __FILE__, __LINE__);
    if ((opdm_file != 40) && (opdm_file != 79) && 
        (opdm_file != PSIF_MO_OPDM)) {
//      errcod = ip_string("OPDM_BASIS",&opdm_basis,0);
      if(options["OPDM_BASIS"].has_changed())
        opdm_basis = options.get_str("OPDM_BASIS");
//      if ((errcod == IPE_KEY_NOT_FOUND) || (errcod != IPE_OK)) {
      else opdm_basis = "AO";
//        opdm_basis = (char *) malloc(3*sizeof(char));
//        strcpy(opdm_basis,"AO");
//      }
//      errcod = ip_string("OPDM_FORMAT",&opdm_format,0);
      if(options["OPDM_FORMAT"].has_changed())     
        opdm_format = options.get_str("OPDM_FORMAT");
//      if ((errcod == IPE_KEY_NOT_FOUND) || (errcod != IPE_OK)) { 
      else opdm_format = "TRIANG";
//        opdm_format = (char *) malloc(7*sizeof(char));
//        strcpy(opdm_format,"TRIANG");
//      }
    }

//    errcod = ip_boolean("WRTNOS",&wrtnos,0);
//    errcod = ip_boolean("ASYMM_OPDM",&asymm_opdm,0);
    wrtnos = options.get_int("WRTNOS");
    asymm_opdm = options.get_int("ASYMM_OPDM");
  }

//  errcod = ip_boolean("SPIN_PROP",&spin_prop,0);
  spin_prop = options.get_int("SPIN_PROP");
  if (spin_prop && read_opdm)
    spin_prop = 0;
  if (iopen == 0)
    spin_prop = 0;

//  errcod = ip_data("PRINT","%d",&print_lvl,0);
  print_lvl = options.get_double("PRINT");
  if (print_lvl < 0)
    print_lvl = 1;
//  errcod = ip_boolean("PRINT_NOS",&print_nos,0);
  print_nos = options.get_int("PRINT_NOS");

//  errcod = ip_boolean("CORREL_CORR",&corr,0);
  corr = options.get_int("CORREL_CORR");
  /*--- corr should be zero since we are not using Psi2 any longer ---*/
  corr = 0;
  if (corr) {
//    errcod = ip_data("ZVEC_FILE","%d",&zvec_file,0);
    zvec_file = options.get_double("ZVEC_FILE");
    if ((zvec_file >= PSIO_MAXUNIT) || (zvec_file <= 0))
      throw PsiException("ZVEC_FILE out of range", __FILE__, __LINE__);
//      punt("ZVEC_FILE out of range");
//    errcod = ip_boolean("DELETE_ZVEC",&delete_zvec,0);
    delete_zvec = options.get_int("DELETE_ZVEC");
  }

//  errcod = ip_data("MPMAX","%d",&mpmax,0);
  mpmax = options.get_double("MPMAX");
  if (mpmax < 1)
    mpmax = 1;
  else if (mpmax > 3)
    mpmax = 3;


//  if (ip_exist("MP_REF_XYZ",0)) {
  if(options["MP_REF_XYZ"].has_changed()) {
//    ip_count("MP_REF_XYZ",&i,0);
    i = options["MP_REF_XYZ"].size();
    if (i != 3)
      throw PsiException("MP_REF_XYZ must have 3 components", __FILE__, __LINE__);
    mp_ref_xyz = options.get_double_array("MP_REF_XYZ");
//      punt("MP_REF_XYZ must have 3 components");
//    for (i=0;i<3;i++) {
//      errcod = ip_data("MP_REF_XYZ","%lf",&mp_ref_xyz[i],1,i);
//      if (errcod != IPE_OK)
//        punt("Error in the definition of MP_REF_XYZ");
//    }
    mp_ref = -1;         /* mp_ref = -1 means that mp_ref_xyz specified by user */
  }
  else {
//    errcod = ip_data("MP_REF","%d",&mp_ref,0);
    mp_ref = options.get_double("MP_REF");
    if (mp_ref <= 0)             /* Default is COM */
      mp_ref = 1;
  }

//  if (ip_exist("LM_REF_XYZ",0)) {
  if(options["LM_REF_XYZ"].has_changed()) {
    i = options["LM_REF_XYZ"].size();
//    ip_count("LM_REF_XYZ",&i,0);
    if (i != 3)
      throw PsiException("LM_REF_XYZ must have 3 components", __FILE__, __LINE__);
    Lm_ref_xyz = options.get_double_array("LM_REF_XYZ");
//      punt("LM_REF_XYZ must have 3 components");
//    for (i=0;i<3;i++) {
//      errcod = ip_data("LM_REF_XYZ","%lf",&Lm_ref_xyz[i],1,i);
//      if (errcod != IPE_OK)
//        punt("Error in the definition of LM_REF_XYZ");
//    }
  }
  else {
    Lm_ref_xyz[0] = Lm_ref_xyz[1] = Lm_ref_xyz[2] = 0.0;
  }

//  errcod = ip_boolean("NUC_ESP",&nuc_esp,0);
  nuc_esp = options.get_int("NUC_ESP");
  if (spin_prop)
    nuc_esp = 1;

  fine_structure_alpha = 1.0;
//  errcod = ip_data("FINE_STRUCTURE_ALPHA","%lf",&fine_structure_alpha,0);
  fine_structure_alpha = options.get_double("FINE_STRUCTURE_ALPHA");

  QED_darwin = 0;
//  errcod = ip_boolean("QED_DARWIN",&QED_darwin,0);
  QED_darwin = options.get_int("QED_DARWIN");
}

}} // namespace psi::oeprop
