#ifndef STRUCTS_H
#define STRUCTS_H
#include <cstdio>
namespace psi { namespace sapt { 

struct params  {
  int t2_restart;
  int df_restart;
  int logfile;
  FILE* logfilename;

  int print;

  int diisvec;
  int maxiter;
  double e_conv;
  double d_conv;
  double schwarz;

  double memory;
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
  double **WBAR; // Electrostatic potential of monomer B in monomer A basis
  };

/* This struct stores the results of the SAPT components */
struct results {

  double hf_int;

  double elst10;

  double exch10;

  double disp20;

  double indr20;
  double indrA_B;
  double indrB_A;

  double exch_indr20;
  double exch_indrA_B;
  double exch_indrB_A;

  double exch_disp20;

  double deltaHF;
  };

}}
#endif
