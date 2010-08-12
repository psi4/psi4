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

void parsing(Options & options)
{
  int i,errcod;
  double xmin,xmax,ymin,ymax,zmin,zmax, *v3;

  /* Set some defaults for certain wavefunctions */
  update_energy_with_MVD = 0;

  wfn = options.get_str("WFN");

  if (wfn == "CI" || wfn == "DETCI" ||
      wfn == "CCSD" || wfn == "DETCAS" ||
      wfn == "CASSCF" || wfn == "RASSCF" ||
      wfn == "MP2" || wfn == "EOM_CCSD" ||
      wfn == "CC2" || wfn == "EOM_CC2" ||
      wfn == "CCSD_MVD") {
        read_opdm = 1;
        opdm_file = PSIF_MO_OPDM;
        corr = 0;
        opdm_basis = "MO";
        opdm_format = "SQUARE";
  }
  else throw PsiException("Unknown wavefunction type", __FILE__, __LINE__);

  if( wfn == "SCF_MVD" || wfn == "CCSD_MVD" )
    update_energy_with_MVD = 1;

  transdens = options.get_bool("TRANSITION_DENSITY");

  if (transdens)
    asymm_opdm = 1;
  else
    asymm_opdm = 0;

  ref = options.get_str("REFERENCE");
  
  /* Parsing section */

  if (options["READ_OPDM"].has_changed())
    options.get_bool("READ_OPDM");

  if (read_opdm) {
    if ( options["OPDM_FILE"].has_changed() )
      opdm_file = options.get_int("OPDM_FILE");
    if ((opdm_file >= PSIO_MAXUNIT) || (opdm_file <= 0))
      throw PsiException("OPDM_FILE out of range", __FILE__, __LINE__);

    if ((opdm_file != 40) && (opdm_file != 79) && 
        (opdm_file != PSIF_MO_OPDM)) {

      if(options["OPDM_BASIS"].has_changed())
        opdm_basis = options.get_str("OPDM_BASIS");
      else
        opdm_basis = "AO";

      if(options["OPDM_FORMAT"].has_changed())
        opdm_format = options.get_str("OPDM_FORMAT");
      else
        opdm_format = "TRIANG";
    }
    wrtnos = options.get_bool("WRTNOS");
    asymm_opdm = options.get_bool("ASYMM_OPDM");
  }

  spin_prop = options.get_bool("SPIN_PROP");
  if (spin_prop && read_opdm)
    spin_prop = 0;
  if (iopen == 0)
    spin_prop = 0;

  print_lvl = options.get_int("PRINT");

  if (print_lvl < 0) print_lvl = 1;

  print_nos = options.get_bool("PRINT_NOS");

  /*--- corr should be zero since we are not using Psi2 any longer ---*/
  corr = 0;

  mpmax = options.get_int("MPMAX");
  if (mpmax < 1)
    mpmax = 1;
  else if (mpmax > 3)
    mpmax = 3;

  if(options["MP_REF_XYZ"].has_changed()) {
    i = options["MP_REF_XYZ"].size();
    if (i != 3)
      throw PsiException("MP_REF_XYZ must have 3 components", __FILE__, __LINE__);

    v3 = options.get_double_array("MP_REF_XYZ");
    memcpy(mp_ref_xyz, v3, 3*sizeof(double));
    delete [] v3;

    mp_ref = -1;         /* mp_ref = -1 means that mp_ref_xyz specified by user */
  }
  else {
    mp_ref = options.get_int("MP_REF");
    if (mp_ref <= 0)             /* Default is COM */
      mp_ref = 1;
  }

  if(options["LM_REF_XYZ"].has_changed()) {
    i = options["LM_REF_XYZ"].size();
    if (i != 3)
      throw PsiException("LM_REF_XYZ must have 3 components", __FILE__, __LINE__);

    v3 = options.get_double_array("LM_REF_XYZ");
    memcpy(Lm_ref_xyz, v3, 3*sizeof(double));
    delete [] v3;
  }
  else {
    Lm_ref_xyz[0] = Lm_ref_xyz[1] = Lm_ref_xyz[2] = 0.0;
  }

  nuc_esp = options.get_bool("NUC_ESP");
  if (spin_prop)
    nuc_esp = 1;

  fine_structure_alpha = options.get_double("FINE_STRUCTURE_ALPHA");

  QED_darwin = options.get_bool("QED_DARWIN");
}

}} // namespace psi::oeprop
