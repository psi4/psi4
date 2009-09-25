/*! \file
    \ingroup CIS
    \brief Enter brief description of file here 
*/

namespace psi { namespace cis {

struct MOInfo {
  int nirreps;           /* no. of irreducible representations */
  int nmo;               /* no. of molecular orbitals */
  int nso;               /* no. of symmetry orbitals */
  int *orbspi;           /* no. of MOs per irrep */
  int *orbsym;           /* orbital symmetry (Pitzer/SO) */
  int *clsdpi;           /* no. of closed-shells per irrep excl. frdocc */
  int *openpi;           /* no. of open-shells per irrep */
  int *uoccpi;           /* no. of unoccupied orbitals per irr. ex. fruocc */
  int *frdocc;           /* no. of frozen core orbitals per irrep */
  int *fruocc;           /* no. of frozen unoccupied orbitals per irrep */
  char **labels;         /* irrep labels */
  int *occpi;            /* no. of occupied orbs. (incl. open) per irrep */
  int *aoccpi;           /* no. of alpha occupied orbs. (incl. open) per irrep */
  int *boccpi;           /* no. of beta occupied orbs. (incl. open) per irrep */
  int *virtpi;           /* no. of virtual orbs. (incl. open) per irrep */
  int *avirtpi;          /* no. of alpha virtual orbs. (incl. open) per irrep */
  int *bvirtpi;          /* no. of beta virtual orbs. (incl. open) per irrep */
  int *occ_sym;          /* relative occupied index symmetry */
  int *aocc_sym;         /* relative alpha occupied index symmetry */
  int *bocc_sym;         /* relative beta occupied index symmetry */
  int *vir_sym;          /* relative virtual index symmetry */
  int *avir_sym;         /* relative alpha virtual index symmetry */
  int *bvir_sym;         /* relative beta virtual index symmetry */
  int *occ_off;          /* occupied orbital offsets within each irrep */
  int *aocc_off;         /* occupied alpha orbital offsets within each irrep */
  int *bocc_off;         /* occupied beta orbital offsets within each irrep */
  int *vir_off;          /* virtual orbital offsets within each irrep */
  int *avir_off;         /* virtual alpha orbital offsets within each irrep */
  int *bvir_off;         /* virtual beta orbital offsets within each irrep */
  double enuc;           /* Nuclear repulsion energy */
  double escf;           /* SCF energy (from chkpt) */
  double eref;           /* Reference energy (file100) */
  double **singlet_evals;/* RHF-CIS singlet excitation energies */
  double **singlet_d;    /* RHF-CIS(D) singlet excitation corrections */
  double **singlet_weakp;/* RHF-CIS(D) singlet weak-pair corrections */
  double **triplet_evals;/* RHF-CIS triplet excitation energies */
  double **triplet_d;    /* RHF-CIS triplet excitation corrections */
  double **uhf_evals;    /* UHF-CIS excitation energies */
  double **uhf_d;        /* UHF-CIS excitation corrections */
};

}} // namespace psi::cis
