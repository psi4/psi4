/*! \file
    \ingroup INTDER
    \brief Enter brief description of file here 
*/
#ifndef _psi_bin_intder_params_h_
#define _psi_bin_intder_params_h_
#include "3dmatrix.h"

namespace psi { namespace intder {

class Params
{
public:
  void checkIsotopes();
  void checkMasses();
  void writeMasses();

  void parseInputFile();
  void parsePrint();
  void printParams();

  void readIntcoGradient();

  int standAlone;  // If this flag exists, intder is running standalone, and all parameters should be read from the
                   // input file (which still needs to be PSI format!)

  int print_lvl;
  int *print_array;

  int natom;
  int natomdummy;
  int ndisp;
  int *atom_dummy;    // Flag for if atom is dummy
  int *moved_dummy;

  double *mass_array;       // Mass array
  double *massTransArray;      // Mass transform array (XT)

  int ncartesians;
  int simplesPresent;
  int nintco;
  int symmCoords;
  int n_spf;
  int n_rcom;

  int massTrans;
  int rcomExist;
  int derlvl;         // Highest order of derivative to be transformed. 0 = evaluation of B. nderv <= 4 
                      // If RCOM is present at_equilibrium+nderv <= 3 
  int highOrder;
  int atEquilibrium; // 0 if at equilibrium, 1 if otherwise
  int transformType;  // 0 = Cart to intern,  1 = intern to cart, 2 = project force constants onto internal space
  int numtest;        // 1 = test B(ij,p) and C(qr,p), 2 = test B(ijk,p) and C(qrs,p), 3 = check orthogonal,
                      // 4 = test Cart. force field for invariance with respect to rotation/translation
                      // 5 = test BuB^T

  int matrixTest;     // ITST variable, to test error of matrix inversions (I think)
  
  int *intcoInclude; // an array of which intco's to include in numerical test
  int analysisType;   // 1 = intern. (gfmat subroutine), 2 = cart (normco subroutine), 3 = both, 4 = SQMFC, 
                      // 5 = rxn path intern, 6 = rxn path cart, 7 = rxn path both, 10+r = reduced dimension space
  int skipTransform;
  int IRintensities;
  int propertyDimension; // 0 = scalar, 1 = vector
  int stop;           // this needs a better name. 1 = stop after B and C mats, 2 = stop after B mats, 3 = stop after cart. projection
  int generateDisp;  // 0 = no. Should we just let optking do this? can optking call intder for bmats? 
                      // 1 = first order, 2 = second order
  int Eckart;         // 1 = read in masses and invoke Eckhart conditions
  int nmodes;         // 0 = TED, 1 = PED  
  int rxnpType;
  int rxnCoord;      // Internal coordinate number that is constrained if analysistype = 5

  double invtol; //Inversion tolerance (hard coded in INTDER2000.f header on page 1)

  double *intcoGrad;
  double **intco2ndDer;
  C3DMatrix *intco3rdDer;
  C3DMatrix *intco4thDerBlock;
  C4DMatrix *intco4thDer;
 
  double **UGF;
  double *FConst1;
  double **FConst2;
  C3DMatrix *FConst3;
  C4DMatrix *FConst4;

};

}} // namespace psi::intder

#endif // header guard

