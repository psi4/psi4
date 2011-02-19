#ifndef STRUCTS_H
#define STRUCTS_H
#include <cstdio>
namespace psi { namespace sapt { 

struct params {

  int print;

  int diisvec;
  int maxiter;
  double e_conv;
  double d_conv;
  double schwarz;

  int foccA;
  int foccB;

  bool nat_orbs;
  bool nat_orbs_t2;
  double occ_cutoff;

  long int memory;

  };

struct workflow {

  // OEints
  int W_ov;
  int W_oo;
  int W_vv;

  // Storage from CHF Solver
  int save_s;
  int save_chf;

  // MO integrals to be computed and stored on disk
  int g_arar;
  int g_bsbs;

  // Types of amplitudes to be formed and stored on disk
  int t_arar;
  int theta_arar;
  int t_bsbs;
  int theta_bsbs;
  int Y2_ar;
  int Y2_bs;
  int t_ar;
  int t_bs;
  int t2_arar;
  int theta2_arar;
  int t2_bsbs;
  int theta2_bsbs;
  int t_arbs;
  int t_bsar;
  int Y3_ar;
  int Y3_bs;

  // Types of intermediates to be formed and stored on disk
  int mp2_opdm;
  int theta_ar_ar; 
  int theta_bs_bs; 
  int theta2_ar_ar; 
  int theta2_bs_bs; 
  int theta_bs_ar; 
  int theta_ar_bs; 
  int gt_ar_arbs;
  int gt_bs_arbs;

  };

/* This struct holds information from the SCF calculations that is used
 * throughout the SAPT program. */
struct calcinfo {

  int nso; // Number of symmetry orbitals
  int nmo; // Number of molecular orbitals
  int nri; // Number of RI orbitals
  int nrio; // Number of RI orbitals + 3 for dressing
  int nsotri; // Size of the SO lower triangle
  int nmotri; // Size of the MO lower triange
  int ntei; // Number of two electron integrals, SO (probably obsolete now)
  int noccA; // Number of occpied orbitals on monomer A
  int noccB; // Number of occpied orbitals on monomer B
  int nvirA; // Number of virtual orbitals on monomer A
  int nvirB; // Number of virtual orbitals on monomer B
  int NA; // Number of electrons of monomer A
  int NB; // Number of electrons of monomer B

  int *ioff; // Standard ioff array
  int *index2i; // Maps a compound ij to the corresponding i
  int *index2j; // Maps a compound ij to the corresponding j

  double enuc_A; // Nuclear repulsion of monomer A
  double enuc_B; // Nuclear repulsion of monomer B
  double enuc_D; // Nuclear repulsion of the dimer
  double eHF_A; // Hartree Fock energy of monomer A
  double eHF_B; // Hartree Fock energy of monomer B
  double eHF_D; // Hartree Fock energy of the dimer

  double *evalsA; // Hartree Fock orbital energies of monomer A
  double *evalsB; // Hartree Fock orbital energies of monomer B
  double *S; // Matrix of overlap integrals
  double *VA; // Nuclear attraction integrals to nuclei of monomer A
  double *VB; // Nuclear attraction integrals to nuclei of monomer B
  double *diagAA; // Diagonal (AA|P) type fitting integrals (summed)
  double *diagBB; // Diagonal (BB|P) type fitting integrals (summed)

  double **C; // SCF coefficient matrix of dimer
  double **CA; // SCF coefficient matrix of monomer A
  double **CB; // SCF coefficient matrix of monomer B
  double **sA; // Approximate CPHF coefficients
  double **sB; // Approximate CPHF coefficients
  double **CHFA; // CPHF coefficients of monomer A
  double **CHFB; // CPHF coefficients of monomer B

  double **S_AB; // Overlap integrals (first index in monomer A basis, second 
                 // index in monomer B basis)
  double **VABB; // Nuclear attraction integrals to nuclei of monomer A
                 // (first index in monomer B basis, second index in
                 // monomer B basis)
  double **VBAA; // Nuclear attraction integrals to nuclei of monomer B
                 // (first index in monomer A basis, second index in
                 // monomer A basis)
  double **VAAB; // Nuclear attraction integrals to nuclei of monomer A
                 // (first index in monomer A basis, second index in
                 // monomer B basis)
  double **VBAB; // Nuclear attraction integrals to nuclei of monomer B
                 // (first index in monomer A basis, second index in
                 // monomer B basis)
  double **WABS; // Electrostatic potential of monomer A in monomer B basis
  double **WABB;
  double **WASS;
  double **WBAR; // Electrostatic potential of monomer B in monomer A basis
  double **WBAA;
  double **WBRR;

  };

struct noinfo {

  int nvirA;
  int nvirB;

  double disp20;

  double *evalsA;
  double *evalsB;

  double **CA;
  double **CB;

  };

/* This struct stores the results of the SAPT components */
struct results {

  double elst10;
  double exch10_s2;
  double exch10;
  double disp20;
  double ind20;
  double exch_ind20;
  double exch_disp20;

  double elst12;
  double exch11;
  double exch12;
  double ind22;

  double disp21;
  double disp22sdq;
  double disp22t;

  double ind30;
  double exch_ind30;
  double disp30;
  double exch_disp30;
  double ind_disp30;
  double exch_ind_disp30;
  double elst13;

  };

}}
#endif
