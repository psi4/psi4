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
  int foccC;

  bool nat_orbs;
  double occ_cutoff;

  long int memory;

  };

struct workflow {

  // DFints
  int save_Jmhalf;

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

  double **Jmhalf;

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
  double exch10;
  double disp20;
  double ind20;
  double exch_ind20;
  double exch_disp20;

  double disp20_os;
  double disp20_ss;
  double exch_disp20_os;
  double exch_disp20_ss;

  double disp20chf;

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

/* This struct holds the information for a three-body SAPT calculation */
struct three_body_info {

  int nso; // Number of symmetry orbitals
  int nmo; // Number of molecular orbitals
  int nrio; // Number of RI orbitals
  int nsotri; // Size of the SO lower triangle
  int nmotri; // Size of the MO lower triange
  int noccA; // Number of occpied orbitals on monomer A
  int noccB; // Number of occpied orbitals on monomer B
  int noccC; // Number of occpied orbitals on monomer C
  int nvirA; // Number of virtual orbitals on monomer A
  int nvirB; // Number of virtual orbitals on monomer B
  int nvirC; // Number of virtual orbitals on monomer C

  int *ioff; // Standard ioff array
  int *index2i; // Maps a compound ij to the corresponding i
  int *index2j; // Maps a compound ij to the corresponding j

  double eHF_A; // Hartree-Fock energy of monomer A
  double eHF_B; // Hartree-Fock energy of monomer B
  double eHF_C; // Hartree-Fock energy of monomer C
  double eHF_AB; // Hartree-Fock energy of dimer AB
  double eHF_AC; // Hartree-Fock energy of dimer AC
  double eHF_BC; // Hartree-Fock energy of dimer BC
  double eHF_ABC; // Hartree-Fock energy of trimer ABC

  double *evalsA; // Hartree Fock orbital energies of monomer A
  double *evalsB; // Hartree Fock orbital energies of monomer B
  double *evalsC; // Hartree Fock orbital energies of monomer C
  double *S; // Matrix of overlap integrals
  double *VA; // Nuclear attraction integrals to nuclei of monomer A
  double *VB; // Nuclear attraction integrals to nuclei of monomer B
  double *VC; // Nuclear attraction integrals to nuclei of monomer C

  double **CA; // SCF coefficient matrix of monomer A
  double **CB; // SCF coefficient matrix of monomer B
  double **CC; // SCF coefficient matrix of monomer C

  double **sA_B; // Approx. CPHF coefficients of monomer A due to monomer B 
  double **sA_C; // Approx. CPHF coefficients of monomer A due to monomer C
  double **sB_A; // Approx. CPHF coefficients of monomer B due to monomer A
  double **sB_C; // Approx. CPHF coefficients of monomer B due to monomer C
  double **sC_A; // Approx. CPHF coefficients of monomer C due to monomer A
  double **sC_B; // Approx. CPHF coefficients of monomer C due to monomer B

  double **CHFA_B; // CPHF coefficients of monomer A due to monomer B 
  double **CHFA_C; // CPHF coefficients of monomer A due to monomer C
  double **CHFB_A; // CPHF coefficients of monomer B due to monomer A
  double **CHFB_C; // CPHF coefficients of monomer B due to monomer C
  double **CHFC_A; // CPHF coefficients of monomer C due to monomer A
  double **CHFC_B; // CPHF coefficients of monomer C due to monomer B

  double **S_AB; // Overlap integrals
  double **S_AC; // Overlap integrals 
  double **S_BC; // Overlap integrals  
  double **S_BA; // Overlap integrals
  double **S_CA; // Overlap integrals
  double **S_CB; // Overlap integrals

  double **VABB; // Nuclear attraction integrals to nuclei of monomer A
  double **VACC; // Nuclear attraction integrals to nuclei of monomer A

  double **VBAA; // Nuclear attraction integrals to nuclei of monomer B
  double **VBCC; // Nuclear attraction integrals to nuclei of monomer B

  double **VCAA; // Nuclear attraction integrals to nuclei of monomer C
  double **VCBB; // Nuclear attraction integrals to nuclei of monomer C

  double **WABS; // Electrostatic potential of monomer A in monomer B basis
  double **WACT; // Electrostatic potential of monomer A in monomer C basis
  double **WBAR; // Electrostatic potential of monomer B in monomer A basis
  double **WBCT; // Electrostatic potential of monomer B in monomer C basis
  double **WCAR; // Electrostatic potential of monomer C in monomer A basis
  double **WCBS; // Electrostatic potential of monomer C in monomer B basis

  double **WABB;
  double **WACC;
  double **WASS;
  double **WATT;

  double **WBAA;
  double **WBCC;
  double **WBRR;
  double **WBTT;

  double **WCAA;
  double **WCBB;
  double **WCRR;
  double **WCSS;
  };

struct three_body_results {

  double e_HF;

  double exch100_s2;
  double exch010_s2;
  double exch001_s2;

  double exch100_s3;
  double exch010_s3;
  double exch001_s3;

  double exch100_s4;
  double exch010_s4;
  double exch001_s4;

  double ind110;
  double ind101;
  double ind011;

  double ind111;
  double ind210;
  double ind021;
  double ind102;
  double ind120;
  double ind012;
  double ind201;

  double ind_disp210;
  double ind_disp021;
  double ind_disp102;
  double ind_disp120;
  double ind_disp012;
  double ind_disp201;

  double disp111;

  double disp3100;
  double disp3010;
  double disp3001;

  double disp211d;
  double disp121d;
  double disp112d;

  double disp211t;
  double disp121t;
  double disp112t;

  double disp220s;
  double disp202s;
  double disp022s;

  double disp220d;
  double disp202d;
  double disp022d;

  double disp220t;
  double disp202t;
  double disp022t;

  double disp220q;
  double disp202q;
  double disp022q;

  double exch_ind200_s2;
  double exch_ind020_s2;
  double exch_ind002_s2;

  double exch_ind110_s2;
  double exch_ind101_s2;
  double exch_ind011_s2;

  double exch_ind200_s3;
  double exch_ind020_s3;
  double exch_ind002_s3;

  double exch_ind110_s3;
  double exch_ind101_s3;
  double exch_ind011_s3;

  double exch_disp200_s2;
  double exch_disp020_s2;
  double exch_disp002_s2;

  double exch_disp110_s2;
  double exch_disp101_s2;
  double exch_disp011_s2;

  double exch_disp200_s3;
  double exch_disp020_s3;
  double exch_disp002_s3;

  double exch_disp110_s3;
  double exch_disp101_s3;
  double exch_disp011_s3;

  double disp20_AB;
  double disp20_AC;
  double disp20_BC;

  double indr20;
  double indrA_B;
  double indrA_C;
  double indrB_A;
  double indrB_C;
  double indrC_A;
  double indrC_B;

  };

}}
#endif
