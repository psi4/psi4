/*! \file
    \ingroup MP2
    \brief Enter brief description of file here 
*/
#ifndef _psi_src_bin_mp2_moinfo_h_
#define _psi_src_bin_mp2_moinfo_h_

struct moinfo {
  int nmo;               /* no. of molecular orbitals */
  int nso;               /* no. of symmetry adapted atomic orbitals */
  int nao;               /* no. of atomic orbitals */
  int nirreps;           /* no. of irreducible representations */
  char **irreplabels;    /* irrep labels */
  int *mopi;             /* all MOs per irrep */
  int ndocc;             /* no. of doubly occupied MOs per irrep */
  int *doccpi;           /* all doubly occupied MOs per irrep */
  int nsocc;             /* no. of singly occupied MOs per irrep */
  int *soccpi;           /* all singly occupied MOs per irrep */
  int nvirt;             /* no. of uoccupied MOs per irrep */
  int *virtpi;           /* all uoccupied MOs per irrep */
  int nfzdocc;           /* no. of frozen occupied MOs */
  int nfzvirt;           /* no. of frozen virtual MOs */
  int nactmo;            /* no. of active MOs */
  int *fzdoccpi;         /* frozen occupied MOs per irrep */
  int *fzvirtpi;         /* frozen virtual MOs per irrep */
  
  int *occpi;            /* occupied MOs per irrep */
  int *virpi;            /* virtual MOs per irrep */
  int *occ_sym;          /* occupied MOs symmetry */
  int *vir_sym;          /* virtual MOs symmetry */
  int *occ_off;          /* occupied orbital offsets within each irrep */ 
  int *vir_off;          /* virtual orbital offsets within each irrep */
  int *qt_occ;           /* CC->QT active occupied reordering array */
  int *qt_vir;           /* CC->QT active virtual reordering array */
  
  int *aoccpi;           /* alpha occupied orbitals per irrep */
  int *boccpi;           /* beta occupied orbitals per irrep */
  int *avirpi;           /* alpha virtual orbitals per irrep */
  int *bvirpi;           /* beta virtual orbitals per irrep */
  int *aocc_sym;         /* alpha occupied MOs symmetry */
  int *bocc_sym;         /* beta occupied MOs symmetry */
  int *avir_sym;         /* alpha virtual MOs symmetry */
  int *bvir_sym;         /* beta virtual MOs symmetry */
  int *aocc_off;         /* alpha occupied orbital offsets within each irrep */
  int *bocc_off;         /* beta occupied orbital offsets within each irrep */
  int *avir_off;         /* alpha virtual orbital offsets within each irrep */
  int *bvir_off;         /* beta virtual orbital offsets within each irrep */
  int *qt_aocc;          /* CC->QT alpha occupied reordering array */
  int *qt_bocc;          /* CC->QT beta occupied reordering array */
  int *qt_avir;          /* CC->QT alpha virtual reordering array*/
  int *qt_bvir;          /* CC->QT beta virtual reordering array */
  
  double Enuc;           /* Nuclear repulsion energy */
  double Escf;           /* SCF energy */
  double Emp2;           /* MP2 energy */
  double emp2_ss;        /* Same-spin correlation energy*/
  double emp2_os;        /* Opposite-spin correlation energy*/
  double escsmp2_ss;     /* Scaled same-spin correlation energy*/
  double escsmp2_os;     /* Scaled opposite-spin correlation energy*/

  double **opdm;         /* One-particle density matrix */
  double **W;            /* Energy-weighted One-particle density matrix */
  double **I;            /* Orbital Lagrangian */
};

#endif /* Header guard */
