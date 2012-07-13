/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/

#include <string>

namespace psi { namespace ccdensity {

/* Input parameters for cclambda */
struct Params {
  double tolerance;
  long int memory;
  int cachelev;
  int aobasis;
  int ref;
  int onepdm; /* produce ONLY the onepdm for properties */
  int onepdm_grid_dump; // dump the onepdm on a grid to a dx file
  int relax_opdm;
  int use_zeta;
  int calc_xi;
  int connect_xi;
  int restart;
  int ground;
  int transition; 
  int dertype;
  std::string wfn;
  int nstates;
  int prop_sym;
  int prop_root;
  int prop_all;
  std::string gauge;

  /* these are used by Xi and twopdm code */
  int G_irr; 
  int R_irr; 
  int L_irr;
  double R0;   
  double L0;
  int ael;
  double cceom_energy;
  double overlap1; /* <L1|R1> */
  double overlap2; /* <L2|R2> */
  double RD_overlap; /* Rmnef <mn||ef> */
  double RZ_overlap; /* <R|zeta> */
};

struct RHO_Params {
  int L_irr;  
  int R_irr; 
  int G_irr;
  int L_root; 
  int R_root;
  int L_ground; 
  int R_ground;
  double R0;   
  double L0;
  char L1A_lbl[32];
  char L1B_lbl[32];
  char L2AA_lbl[32];
  char L2BB_lbl[32];
  char L2AB_lbl[32];
  char R1A_lbl[32];
  char R1B_lbl[32];
  char R2AA_lbl[32];
  char R2BB_lbl[32];
  char R2AB_lbl[32];
  double cceom_energy;
  double overlap1; /* <L1|R1> */
  double overlap2; /* <L2|R2> */
  double RD_overlap; /* Rmnef <mn||ef> */
  char DIJ_lbl[10];
  char Dij_lbl[10];
  char DAB_lbl[10];
  char Dab_lbl[10];
  char DIA_lbl[10];
  char Dia_lbl[10];
  char DAI_lbl[10];
  char Dai_lbl[10];
  char opdm_lbl[32];
  char opdm_a_lbl[32];
  char opdm_b_lbl[32];
};

struct TD_Params {
  int irrep;
  int root;
  double R0;
  double cceom_energy;
  char L1A_lbl[32];
  char L1B_lbl[32];
  char L2AA_lbl[32];
  char L2BB_lbl[32];
  char L2AB_lbl[32];
  char R1A_lbl[32];
  char R1B_lbl[32];
  char R2AA_lbl[32];
  char R2BB_lbl[32];
  char R2AB_lbl[32];
  double OS;
  double RS_length;
  double RS_velocity;
};


}} // namespace psi::ccdensity
