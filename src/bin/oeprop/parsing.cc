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

  errcod = ip_string("WFN", &wfn, 0);

  if (errcod == IPE_OK) {
    if (!strcmp(wfn, "CI") || !strcmp(wfn, "DETCI") ||
        !strcmp(wfn, "CCSD") || !strcmp(wfn, "DETCAS") ||
        !strcmp(wfn, "CASSCF") || !strcmp(wfn, "RASSCF") ||
        !strcmp(wfn, "MP2") || !strcmp(wfn, "EOM_CCSD") ||
        !strcmp(wfn, "CC2") || !strcmp(wfn, "EOM_CC2") ||
        !strcmp(wfn, "CCSD_MVD"))  {
      read_opdm = 1;
      opdm_file = PSIF_MO_OPDM;
      corr = 0;
      opdm_basis = (char *) malloc(3*sizeof(char));
      strcpy(opdm_basis,"MO");
      opdm_format = (char *) malloc(7*sizeof(char));
      strcpy(opdm_format,"SQUARE");
    }
    if (!strcmp(wfn, "SCF_MVD") || !strcmp(wfn, "CCSD_MVD"))
      update_energy_with_MVD = 1;
  }
  else {
    fprintf(outfile,"Incorrect wave function type!\n");
    abort();
  }

  transdens = 0;
  errcod = ip_boolean("TRANSITION_DENSITY", &transdens, 0);
  if (transdens)
    asymm_opdm = 1;
  else
    asymm_opdm = 0;

  errcod = ip_string("REFERENCE", &ref, 0);
  
  if (errcod != IPE_OK) {
    fprintf(outfile,"Incorrect reference!\n");
    abort();
  }
    
  /* Parsing section */

  errcod = ip_boolean("READ_OPDM",&read_opdm,0);
  if (read_opdm) {
    errcod = ip_data("OPDM_FILE","%d",&opdm_file,0);
    if ((opdm_file >= PSIO_MAXUNIT) || (opdm_file <= 0))
      punt("OPDM_FILE out of range");
    if ((opdm_file != 40) && (opdm_file != 79) && 
        (opdm_file != PSIF_MO_OPDM)) {
      errcod = ip_string("OPDM_BASIS",&opdm_basis,0);
      if ((errcod == IPE_KEY_NOT_FOUND) || (errcod != IPE_OK)) {
        opdm_basis = (char *) malloc(3*sizeof(char));
        strcpy(opdm_basis,"AO");
      }
      errcod = ip_string("OPDM_FORMAT",&opdm_format,0);
      if ((errcod == IPE_KEY_NOT_FOUND) || (errcod != IPE_OK)) { 
        opdm_format = (char *) malloc(7*sizeof(char));
        strcpy(opdm_format,"TRIANG");
      }
    }

    errcod = ip_boolean("WRTNOS",&wrtnos,0);
    errcod = ip_boolean("ASYMM_OPDM",&asymm_opdm,0);
  }

  errcod = ip_boolean("SPIN_PROP",&spin_prop,0);
  if (spin_prop && read_opdm)
    spin_prop = 0;
  if (iopen == 0)
    spin_prop = 0;

  errcod = ip_data("PRINT","%d",&print_lvl,0);
  if (print_lvl < 0)
    print_lvl = 1;
  errcod = ip_boolean("PRINT_NOS",&print_nos,0);

  errcod = ip_boolean("CORREL_CORR",&corr,0);
  /*--- corr should be zero since we are not using Psi2 any longer ---*/
  corr = 0;
  if (corr) {
    errcod = ip_data("ZVEC_FILE","%d",&zvec_file,0);
    if ((zvec_file >= PSIO_MAXUNIT) || (zvec_file <= 0))
      punt("ZVEC_FILE out of range");
    errcod = ip_boolean("DELETE_ZVEC",&delete_zvec,0);
  }

  errcod = ip_data("MPMAX","%d",&mpmax,0);
  if (mpmax < 1)
    mpmax = 1;
  else if (mpmax > 3)
    mpmax = 3;


  if (ip_exist("MP_REF_XYZ",0)) {
    ip_count("MP_REF_XYZ",&i,0);
    if (i != 3)
      punt("MP_REF_XYZ must have 3 components");
    for (i=0;i<3;i++) {
      errcod = ip_data("MP_REF_XYZ","%lf",&mp_ref_xyz[i],1,i);
      if (errcod != IPE_OK)
        punt("Error in the definition of MP_REF_XYZ");
    }
    mp_ref = -1;         /* mp_ref = -1 means that mp_ref_xyz specified by user */
  }
  else {
    errcod = ip_data("MP_REF","%d",&mp_ref,0);
    if (mp_ref <= 0)             /* Default is COM */
      mp_ref = 1;
  }

  if (ip_exist("LM_REF_XYZ",0)) {
    ip_count("LM_REF_XYZ",&i,0);
    if (i != 3)
      punt("LM_REF_XYZ must have 3 components");
    for (i=0;i<3;i++) {
      errcod = ip_data("LM_REF_XYZ","%lf",&Lm_ref_xyz[i],1,i);
      if (errcod != IPE_OK)
        punt("Error in the definition of LM_REF_XYZ");
    }
  }
  else {
    Lm_ref_xyz[0] = Lm_ref_xyz[1] = Lm_ref_xyz[2] = 0.0;
  }

  errcod = ip_boolean("NUC_ESP",&nuc_esp,0);
  if (spin_prop)
    nuc_esp = 1;

  fine_structure_alpha = 1.0;
  errcod = ip_data("FINE_STRUCTURE_ALPHA","%lf",&fine_structure_alpha,0);

  QED_darwin = 0;
  errcod = ip_boolean("QED_DARWIN",&QED_darwin,0);
}

}} // namespace psi::oeprop
