/*! \file
    \ingroup INTDER
    \brief Enter brief description of file here 
*/
#include "params.h"
#include "displacements.h"
#include "intco.h"

#define EXTERN
#include "globals.h"
#undef EXTERN

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libchkpt/chkpt.h>
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <physconst.h>
#include <libpsio/psio.h>
#include <psifiles.h>
#include <masses.h>

using namespace psi::intder;

namespace psi { namespace intder {
extern Displacements gDisplacements;
extern InternalCoordinates gIntCo;
}}

void Params::parseInputFile()
{
  
  int errcod;
  char *junk;
  int a = 0;
  int i = 0;
  int val = 0;
  double dval = 0.0;
  int num_dummy = 0;
  int dummy_start;  
  
  natom = chkpt_rd_natom();
  natomdummy = chkpt_rd_nallatom();
  ncartesians = 3 * natom;
  atom_dummy = chkpt_rd_atom_dummy();
  num_dummy = natomdummy - natom;

  n_spf = 0;
  
  print_lvl = 3000;  //Note there is a strange print_lev convention details in intder.1
  invtol = 1.0E-10;  //Setting inversion testing criteria for FLIN function
  
  errcod = ip_data("PRINT_LVL","%d",&(print_lvl),0);
  parsePrint();

  standAlone = 0;
  ip_boolean("STANDALONE",&standAlone,0);

  // Form the mapping
  moved_dummy = new int[natomdummy];
  dummy_start = natom;
  for (i=0; i < natomdummy; i++) {
    if (atom_dummy[i]) {
      moved_dummy[i] = dummy_start;
      num_dummy++;
      dummy_start++;
    }
    else {
      moved_dummy[i] = i - num_dummy;
    }
    if(num_dummy)
      fprintf(outfile,"%d\n", i);
  }
  if(num_dummy) {
    fprintf(outfile, "Dummy Atom Relocation:\n");
    for (i=0; i < natomdummy; i++) {
      fprintf(outfile, "\tAtom #%d from input was moved to become Atom #%d\n", i, moved_dummy[i]);
    }
  }
  
  // If intco.dat exists append it to the current input parser.
  simplesPresent = 0;

  ffile_noexit(&fp_intco, "intco.dat", 2);

  if(ip_exist(":INTCO", 0)) {
    ip_cwk_add(":INTCO");
    simplesPresent = 1;
  }
  else if(fp_intco != NULL) {
    ip_append(fp_intco, outfile);
    ip_cwk_add(":INTCO");
    simplesPresent = 1;
  }
  else {
    fprintf(outfile, "\nNo intco's in input.dat or intco.dat!!\n");
    exit(PSI_RETURN_FAILURE); 
  }
  
  fprintf(outfile, "simples_present? %s\n", simplesPresent == 1 ? "true" : "false");
  
  gIntCo.loadInternalCoordinates();  
  nintco = gIntCo.InternalCoordinateSize();
  
  atEquilibrium = 0;
  ip_boolean("EQUIL",&atEquilibrium,0);

  derlvl = 0;
  if(ip_exist("DERLVL",0)) {
    errcod = ip_string("DERLVL", &(junk),0);
    if(errcod != IPE_OK) derlvl = 0;  // Default is to print B matrix
    else if(!strcmp(junk, "FIRST")) derlvl = 1;
    else if(!strcmp(junk, "SECOND")) derlvl = 2;
    else if(!strcmp(junk, "THIRD")) derlvl = 3;
    else if(!strcmp(junk, "FOURTH")) derlvl = 4;
    free(junk);
  }

  highOrder = 0;
  ip_boolean("HIGH_ORDER",&highOrder,0);
  
  transformType = 0;
  if(ip_exist("TRANSFORM",0)) {
    errcod = ip_string("TRANSFORM", &(junk),0);
    if(errcod != IPE_OK) transformType = 0;  // Default is Cartesians to internals
    else if(!strcmp(junk, "C_TO_I")) transformType = 0;
    else if(!strcmp(junk, "I_TO_C")) transformType = 1;
    else if(!strcmp(junk, "I_TO_C_R")) transformType = 2;
    else if(!strcmp(junk, "PROJECTCART")) transformType = 3;
    free(junk);
  }
  
  stop = 0;
  if(ip_exist("STOP",0)) {
    errcod = ip_string("STOP", &(junk),0);
    if(errcod != IPE_OK) numtest = 0;  // Default is to run the program!
    else if(!strcmp(junk, "B_AND_C")) stop = 1;
    else if(!strcmp(junk, "B")) stop = 2;
    else if(!strcmp(junk, "C")) stop = 3;
    free(junk);
  }

  analysisType = 0;
  skipTransform = 0;
  if(ip_exist("FREQ_ANAL",0)) {
    errcod = ip_string("FREQ_ANAL", &(junk),0);
    if(!strcmp(junk, "NONE")) analysisType = 0;
    else if(!strcmp(junk, "INT_CO")) analysisType = 1;
    else if(!strcmp(junk, "CART_CO")) analysisType = 2;
    else if(!strcmp(junk, "BOTH")) analysisType = 3;
    else if(!strcmp(junk, "SQMFC")) analysisType = 4;
    else if(!strcmp(junk, "RXN_PATH")) analysisType = 5;
    else if(!strcmp(junk, "INT_CO_R")) analysisType = 6;
    else {
      printf("Invalid value of input keyword DERTYPE: %s\n", junk);
      exit(PSI_RETURN_FAILURE); 
    }
    ip_boolean("SKIP_T",&skipTransform,0);
    free(junk);
  }
      
  rxnpType = 0;
  if(analysisType == 5) {
    if(!ip_exist("RXNP_TYPE",0) || !ip_exist("RXN_COORD",0)) {
      fprintf(outfile,"RXN_PATH analysis, but either no RXNP_TYPE or RXN_COORD specified!");
      exit(2);
    }
    else {
      errcod = ip_string("RXNP_TYPE", &(junk),0);
      ip_data("RXN_COORD", "%i", &rxnCoord, 0);
    }
    
    if(!strcmp(junk, "CART")) rxnpType = 0;
    else if(!strcmp(junk, "INT")) rxnpType = 1;
    else if(!strcmp(junk, "BOTH")) rxnpType = 2;      
    else {
      printf("Invalid value of input keyword RXNP_TYPE: %s\n", junk);
      exit(PSI_RETURN_FAILURE); 
    }
    free(junk);    
  }

  massTrans = 0;
  ip_boolean("MASS_TRANS",&massTrans,0);
  
  numtest = 0;
  if(ip_exist("NUMTEST",0)) {
    errcod = ip_string("NUMTEST", &(junk),0);
    if(errcod != IPE_OK) numtest = 0;  // Default is no numerical test
    else if(!strcmp(junk, "NONE")) numtest = 0;
    else if(!strcmp(junk, "SECOND")) numtest = 1;
    else if(!strcmp(junk, "THIRD")) numtest = 2;
    else if(!strcmp(junk, "ORTHOG")) {
      numtest = 3;
      stop = 1;
      transformType = 1;
      fprintf(outfile, "Switched to C_to_I-type transformation in order to check orthogonality conditions!\n");
    }
    else if(!strcmp(junk, "INVARIANCE")) numtest = 4;
    else if(!strcmp(junk, "BUBT")) numtest = 5;

    if(ip_exist("INT_INCLUDE",0))
      ip_count("INT_INCLUDE",&a,0);
    if(a > nintco) {
      fprintf(outfile, "More intco's to include from numerical test than there are intco's!\n");
      exit(PSI_RETURN_FAILURE);
    }
    
    intcoInclude = new int[a];

    ip_int_array("INT_INCLUDE", &intcoInclude[i], a);

    free(junk);
  }
       
  ip_boolean("IRINT_R", &IRintensities,0); //Needs to be able to read IR intensities from an input vector INT_DIPDER, 
                                           //or perhaps an external file (18) if true
  
  propertyDimension = 0;
  ip_boolean("PROP_DIM", &propertyDimension,0);
  if(propertyDimension && (atEquilibrium != 1)) {
    printf("Cannot determine vector quantities at non-equilibrium geometry!\nSet EQUIL = TRUE\n");
    exit(PSI_RETURN_FAILURE);
  }
  
  Eckart = 0;
  ip_boolean("ECKART", &Eckart,0);
  
  nmodes = 0;
  ip_boolean("PED", &nmodes,0);

  if(numtest == 5) {
    matrixTest = 2;
    numtest = 0;
  } 
  else if(abs(numtest) == 3) {
    stop = 1;
    (numtest > 0) ? transformType = 1 : transformType = -1;
    //Umm... will it matter that we don't have a default value for matrixTest here?
  }
  else
    matrixTest = 1;  
}      

void Params::checkIsotopes()
{
  int a = 0;
  int i = 0;
  int j = 0;
  double *mass = NULL;
  char line1[80];

  if (!ip_exist("ISOTOPES", 0))
    return;

  ip_count("ISOTOPES", &a, 0);
  if (a != natom) {
    fprintf(outfile, "ISOTOPES array has wrong dimension.\n");
    exit(2);
  }

  mass = new double[natom];
  massTransArray = new double[natom];

  for (i = 0; i < natom; i++) {
    ip_data("ISOTOPES", "%s", line1, 1, i);
    for (j = 0; j < 138; j++) {
      if (strcmp(line1, mass_labels[j]) == 0) {
        mass[i] = atomic_masses[j];
	if(!propertyDimension || transformType > 0)
	  massTransArray[i] = atomic_masses[j];
	else
	  massTransArray[i] = 1.0;
        break;
      }
    }
    if (j == 138) {
      fprintf(outfile, "Isotope label %s is unidentifiable.\n", line1);
      exit(2);
    }
  }

  // Tell all the displacements to use these atomic weights
  gDisplacements.useMasses(mass);
  delete[] mass;
}

void Params::checkMasses()
{
  int a = 0;
  int i = 0;
  double *mass = NULL;
  double val = 0.0;

  if (!ip_exist("MASSES", 0))
    return;

  ip_count("MASSES", &a, 0);
  if (a != natom) {
    fprintf(outfile, "MASSES array has wrong dimension.\n");
    exit(2);
  }

  mass = new double[natom];
  massTransArray = new double[natom];

  for (i = 0; i < natom; i++) {
    ip_data("MASSES", "%lf", &val, 1, i);
    mass[i] = val;
    if(!propertyDimension || transformType > 0)
      massTransArray[i] = val;
    else
      massTransArray[i] = 1.0;
  }
  
  // Tell all the displacements to use these atomic weights
  gDisplacements.useMasses(mass);
  delete[] mass;
}

//For BINVRT and presumably other functions, we need to keep the mass array around
//Should we just integrate this into the checkMasses function?
void Params::writeMasses()
{
  int a = 0;
  int i = 0;
  double val = 0.0;

  if (!ip_exist("MASSES", 0))
    return;

  ip_count("MASSES", &a, 0);
  if (a != natom) {
    fprintf(outfile, "MASSES array has wrong dimension.\n");
    exit(2);
  }

  mass_array = new double[natom];
  for (i = 0; i < natom; i++) {
    if(transformType == 3) 
      mass_array[i] = 1.0;
    else {
      ip_data("MASSES", "%lf", &val, 1, i);
      mass_array[i] = val;
    }
  }
}

//Read internal coordinate values and the internal coordinate gradients. OPTKING does
//not store the gradients in terms of intco's, only in terms of cartesians. Do we do this
//transformation in intder eventually? Could we then hack optking to have optking do this on
//its own?
/*
 void Params::readIntcoGradient()
{
  if(at_equilibrium != 0)
  new 1DintcoGrad double[nintco];
    
    if(at_equilibrium == 0 || derlvl >= 1)
    2DintcoGrad = init_matrix(nintco,nintco);
    
    if(derlvl >= 2)
  C3DMatrix *3DintcoGrad = new C3DMatrix(); //oh shit.
    }  
*/

void Params::parsePrint()
{
  print_array = new int[4];
  print_array[3] = print_lvl / 1000;
  print_array[2] = (print_lvl - (print_array[3] * 1000)) / 100;
  print_array[1] = (print_lvl - (print_array[3] * 1000) - (print_array[2] * 100)) / 10;
  print_array[0] = (print_lvl - (print_array[3] * 1000) - (print_array[2] * 100) - (print_array[1] * 10));
  fprintf(outfile, "\nPrint level is %i %i %i %i\n", print_array[3], print_array[2], print_array[1], print_array[0]);
}

void Params::printParams()
{
  fprintf(outfile, "\n\tPROGRAM OPTIONS: \n\n");
  fprintf(outfile, "\t************************************\n");
  fprintf(outfile, "\tAt Equilibrium: %i\n", atEquilibrium);
  fprintf(outfile, "\tDerivative level : %i\n", derlvl);
  fprintf(outfile, "\tHigher order transformation?: %i\n",highOrder);
  fprintf(outfile, "\tTransform type: %i\n",transformType); 
  fprintf(outfile, "\tStopping point: %i\n",stop);
  fprintf(outfile, "\tFrequency analysis type: %i\n", analysisType);
  fprintf(outfile, "\tSkip Transform?: %i\n", skipTransform);
  fprintf(outfile, "\tReaction path analysis type: %i\n", rxnpType);
  if(rxnpType)
    fprintf(outfile, "\tReaction coordinate: %i\n", rxnCoord);
  fprintf(outfile, "\tMass transformation?: %i\n", massTrans);
  fprintf(outfile, "\tNumerical test type: %i\n", numtest);
  if(numtest) //Needs to be formatted correctly
    fprintf(outfile, "\tINTCO's to numerically test: %i\n", intcoInclude[0]);
  fprintf(outfile, "\tIR intensities?: %i\n", IRintensities);
  fprintf(outfile, "\tTransformation property dimension: %i\n", propertyDimension);
  fprintf(outfile, "\tEckart transformation?: %i\n", Eckart);
  fprintf(outfile, "\tPED?: %i\n", nmodes);
  fprintf(outfile, "\tmatrixTest for Linear dependence: %i", matrixTest);
  fprintf(outfile, "\n\t************************************\n");
}
