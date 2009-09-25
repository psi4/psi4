/*! \file
    \ingroup CCSORT
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include "Params.h"
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccsort {

/* 
** fock(): Build the alpha and beta Fock matrices from the
** one-electron integrals/frozen-core operator and active two-electron
** integrals on disk.
**
** TDC, 1996
** Modified to include UHF references, TDC, June 2001.
**
** Notes:
**
** (1) This routine isn't absolutely necessary for RHF and UHF
** references, as one could simply read the Fock eigenvalues from
** PSIF_CHKPT and be happy with that. However, this code is useful as a
** partial check of the integral transformation and sorting routines.
**
** (2) An alternative but currently unused algorithm may be found in 
** fock_build.c.
*/

void fock_rhf(void);
void fock_uhf(void);

void fock(void)
{
  if(params.ref == 2) fock_uhf();
  else fock_rhf();
}

void fock_uhf(void)
{
  int h, nirreps;
  int i, j, I, J, Gi, Gj, IM, JM, MI, MJ;
  int a, b, A, B, Ga, Gb, AM, BM, MA, MB;
  int m, M, Gm;
  int *aoccpi, *boccpi, *aocc_off, *bocc_off;
  int *avirtpi, *bvirtpi, *avir_off, *bvir_off;
  dpdfile2 fIJ, fij, fAB, fab, fIA, fia;
  dpdfile2 hIJ, hij, hAB, hab, hIA, hia;
  dpdbuf4 A_AA, A_BB, A_AB, C_AA, C_BB, C_AB, E_AA, E_BB, E_AB;

  nirreps = moinfo.nirreps;
  aoccpi = moinfo.aoccpi;
  boccpi = moinfo.boccpi;
  avirtpi = moinfo.avirtpi;
  bvirtpi = moinfo.bvirtpi;
  aocc_off = moinfo.aocc_off;
  bocc_off = moinfo.bocc_off;
  avir_off = moinfo.avir_off;
  bvir_off = moinfo.bvir_off;

  dpd_file2_init(&hIJ, CC_OEI, 0, 0, 0, "h(I,J)");
  dpd_file2_init(&hij, CC_OEI, 0, 2, 2, "h(i,j)");
  dpd_file2_init(&hAB, CC_OEI, 0, 1, 1, "h(A,B)");
  dpd_file2_init(&hab, CC_OEI, 0, 3, 3, "h(a,b)");
  dpd_file2_init(&hIA, CC_OEI, 0, 0, 1, "h(I,A)");
  dpd_file2_init(&hia, CC_OEI, 0, 2, 3, "h(i,a)");

  dpd_file2_mat_init(&hIJ);
  dpd_file2_mat_init(&hij);
  dpd_file2_mat_init(&hAB);
  dpd_file2_mat_init(&hab);
  dpd_file2_mat_init(&hIA);
  dpd_file2_mat_init(&hia);

  dpd_file2_mat_rd(&hIJ);
  dpd_file2_mat_rd(&hij);
  dpd_file2_mat_rd(&hAB);
  dpd_file2_mat_rd(&hab);
  dpd_file2_mat_rd(&hIA);
  dpd_file2_mat_rd(&hia);

  dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_init(&fij, CC_OEI, 0, 2, 2, "fij");

  dpd_file2_mat_init(&fIJ);
  dpd_file2_mat_init(&fij);

  for(h=0; h < nirreps; h++) {

    for(i=0; i < aoccpi[h]; i++)
      for(j=0; j < aoccpi[h]; j++)
	fIJ.matrix[h][i][j] = hIJ.matrix[h][i][j];

    for(i=0; i < boccpi[h]; i++)
      for(j=0; j < boccpi[h]; j++)
	fij.matrix[h][i][j] = hij.matrix[h][i][j];
  }

  dpd_buf4_init(&A_AA, CC_AINTS, 0, 0, 0, 0, 0, 1, "A <IJ|KL>");
  dpd_buf4_init(&A_AB, CC_AINTS, 0, 22, 22, 22, 22, 0, "A <Ij|Kl>");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&A_AA, h);
    dpd_buf4_mat_irrep_rd(&A_AA, h);
    dpd_buf4_mat_irrep_init(&A_AB, h);
    dpd_buf4_mat_irrep_rd(&A_AB, h);
    for(Gi=0; Gi < nirreps; Gi++) {
      Gj = Gi; Gm = Gi^h;
      for(i=0; i < aoccpi[Gi]; i++) {
	I = aocc_off[Gi] + i;
	for(j=0; j < aoccpi[Gj]; j++) {
	  J = aocc_off[Gj] + j;
	  for(m=0; m < aoccpi[Gm]; m++) {
	    M = aocc_off[Gm] + m;
	    IM = A_AA.params->rowidx[I][M];
	    JM = A_AA.params->colidx[J][M];
	    fIJ.matrix[Gi][i][j] += A_AA.matrix[h][IM][JM];
	  }
	  for(m=0; m < boccpi[Gm]; m++) {
	    M = bocc_off[Gm] + m;
	    IM = A_AB.params->rowidx[I][M];
	    JM = A_AB.params->colidx[J][M];
	    fIJ.matrix[Gi][i][j] += A_AB.matrix[h][IM][JM];
	  }
	}
      }
    }
    dpd_buf4_mat_irrep_close(&A_AA, h);
    dpd_buf4_mat_irrep_close(&A_AB, h);
  }
  dpd_buf4_close(&A_AA);
  dpd_buf4_close(&A_AB);

  dpd_buf4_init(&A_BB, CC_AINTS, 0, 10, 10, 10, 10, 1, "A <ij|kl>");
  dpd_buf4_init(&A_AB, CC_AINTS, 0, 22, 22, 22, 22, 0, "A <Ij|Kl>");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&A_BB, h);
    dpd_buf4_mat_irrep_rd(&A_BB, h);
    dpd_buf4_mat_irrep_init(&A_AB, h);
    dpd_buf4_mat_irrep_rd(&A_AB, h);
    for(Gi=0; Gi < nirreps; Gi++) {
      Gj = Gi; Gm = Gi^h;
      for(i=0; i < boccpi[Gi]; i++) {
	I = bocc_off[Gi] + i;
	for(j=0; j < boccpi[Gj]; j++) {
	  J = bocc_off[Gj] + j;
	  for(m=0; m < boccpi[Gm]; m++) {
	    M = bocc_off[Gm] + m;
	    IM = A_BB.params->rowidx[I][M];
	    JM = A_BB.params->colidx[J][M];
	    fij.matrix[Gi][i][j] += A_BB.matrix[h][IM][JM];
	  }
	  for(m=0; m < aoccpi[Gm]; m++) {
	    M = aocc_off[Gm] + m;
	    MI = A_AB.params->rowidx[M][I];
	    MJ = A_AB.params->colidx[M][J];
	    fij.matrix[Gi][i][j] += A_AB.matrix[h][MI][MJ];
	  }
	}
      }
    }
    dpd_buf4_mat_irrep_close(&A_BB, h);
    dpd_buf4_mat_irrep_close(&A_AB, h);
  }
  dpd_buf4_close(&A_BB);
  dpd_buf4_close(&A_AB);

  dpd_file2_mat_wrt(&fIJ);
  dpd_file2_mat_wrt(&fij);

  dpd_file2_close(&fIJ);
  dpd_file2_close(&fij);

  dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
  dpd_file2_init(&fab, CC_OEI, 0, 3, 3, "fab");

  dpd_file2_mat_init(&fAB);
  dpd_file2_mat_init(&fab);

  for(h=0; h < nirreps; h++) {

    for(a=0; a < avirtpi[h]; a++)
      for(b=0; b < avirtpi[h]; b++)
	fAB.matrix[h][a][b] = hAB.matrix[h][a][b];

    for(a=0; a < bvirtpi[h]; a++)
      for(b=0; b < bvirtpi[h]; b++)
	fab.matrix[h][a][b] = hab.matrix[h][a][b];
  }

  dpd_buf4_init(&C_AA, CC_CINTS, 0, 20, 20, 20, 20, 0, "C <IA||JB>");
  dpd_buf4_init(&C_AB, CC_CINTS, 0, 26, 26, 26, 26, 0, "C <Ai|Bj>");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&C_AA, h);
    dpd_buf4_mat_irrep_rd(&C_AA, h);
    dpd_buf4_mat_irrep_init(&C_AB, h);
    dpd_buf4_mat_irrep_rd(&C_AB, h);
    for(Ga=0; Ga < nirreps; Ga++) {
      Gb = Ga; Gm = Ga^h;
      for(a=0; a < avirtpi[Ga]; a++) {
	A = avir_off[Ga] + a;
	for(b=0; b < avirtpi[Gb]; b++) {
	  B = avir_off[Gb] + b;
	  for(m=0; m < aoccpi[Gm]; m++) {
	    M = aocc_off[Gm] + m;
	    MA = C_AA.params->rowidx[M][A];
	    MB = C_AA.params->colidx[M][B];
	    fAB.matrix[Ga][a][b] += C_AA.matrix[h][MA][MB];
	  }
	  for(m=0; m < boccpi[Gm]; m++) {
	    M = bocc_off[Gm] + m;
	    AM = C_AB.params->rowidx[A][M];
	    BM = C_AB.params->colidx[B][M];
	    fAB.matrix[Ga][a][b] += C_AB.matrix[h][AM][BM];
	  }
	}
      }
    }
    dpd_buf4_mat_irrep_close(&C_AA, h);
    dpd_buf4_mat_irrep_close(&C_AB, h);
  }
  dpd_buf4_close(&C_AA);
  dpd_buf4_close(&C_AB);

  dpd_buf4_init(&C_BB, CC_CINTS, 0, 30, 30, 30, 30, 0, "C <ia||jb>");
  dpd_buf4_init(&C_AB, CC_CINTS, 0, 24, 24, 24, 24, 0, "C <Ia|Jb>");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&C_BB, h);
    dpd_buf4_mat_irrep_rd(&C_BB, h);
    dpd_buf4_mat_irrep_init(&C_AB, h);
    dpd_buf4_mat_irrep_rd(&C_AB, h);
    for(Ga=0; Ga < nirreps; Ga++) {
      Gb = Ga; Gm = Ga^h;
      for(a=0; a < bvirtpi[Ga]; a++) {
	A = bvir_off[Ga] + a;
	for(b=0; b < bvirtpi[Gb]; b++) {
	  B = bvir_off[Gb] + b;
	  for(m=0; m < boccpi[Gm]; m++) {
	    M = bocc_off[Gm] + m;
	    MA = C_BB.params->rowidx[M][A];
	    MB = C_BB.params->colidx[M][B];
	    fab.matrix[Ga][a][b] += C_BB.matrix[h][MA][MB];
	  }
	  for(m=0; m < aoccpi[Gm]; m++) {
	    M = aocc_off[Gm] + m;
	    MA = C_AB.params->rowidx[M][A];
	    MB = C_AB.params->colidx[M][B];
	    fab.matrix[Ga][a][b] += C_AB.matrix[h][MA][MB];
	  }
	}
      }
    }
    dpd_buf4_mat_irrep_close(&C_BB, h);
    dpd_buf4_mat_irrep_close(&C_AB, h);
  }
  dpd_buf4_close(&C_BB);
  dpd_buf4_close(&C_AB);

  dpd_file2_mat_wrt(&fAB);
  dpd_file2_mat_wrt(&fab);

  dpd_file2_close(&fAB);
  dpd_file2_close(&fab);

  /* Prepare the alpha and beta occ-vir Fock matrix files */
  dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
  dpd_file2_init(&fia, CC_OEI, 0, 2, 3, "fia");
  dpd_file2_mat_init(&fIA);
  dpd_file2_mat_init(&fia);

  /* One-electron (frozen-core) contributions */
  for(h=0; h < nirreps; h++) {

    for(i=0; i < aoccpi[h]; i++) 
      for(a=0; a < avirtpi[h]; a++) 
	fIA.matrix[h][i][a] = hIA.matrix[h][i][a];   

    for(i=0; i < boccpi[h]; i++) 
      for(a=0; a < bvirtpi[h]; a++) 
	fia.matrix[h][i][a] = hia.matrix[h][i][a];   
  }

  /* Two-electron contributions */

  /* Prepare the E integral buffers */
  dpd_buf4_init(&E_AA, CC_EINTS, 0, 21, 0, 21, 0, 1, "E <AI|JK>");
  dpd_buf4_init(&E_AB, CC_EINTS, 0, 26, 22, 26, 22, 0, "E <Ai|Jk>");

  for(h=0; h < nirreps; h++) {

    dpd_buf4_mat_irrep_init(&E_AA, h);
    dpd_buf4_mat_irrep_init(&E_AB, h);
    dpd_buf4_mat_irrep_rd(&E_AA, h);
    dpd_buf4_mat_irrep_rd(&E_AB, h);

    /* Loop over irreps of the target */
    for(Gi=0; Gi < nirreps; Gi++) {
      Ga = Gi; Gm = h^Gi;

      /* Loop over orbitals of the target */
      for(i=0; i < aoccpi[Gi]; i++) {
	I = aocc_off[Gi] + i;
	for(a=0; a < avirtpi[Ga]; a++) {
	  A = avir_off[Ga] + a;

	  for(m=0; m < aoccpi[Gm]; m++) {
	    M = aocc_off[Gm] + m;

	    AM = E_AA.params->rowidx[A][M];
	    IM = E_AA.params->colidx[I][M];

	    fIA.matrix[Gi][i][a] += E_AA.matrix[h][AM][IM];

	  }
	}
      }

      /* Loop over orbitals of the target */
      for(i=0; i < aoccpi[Gi]; i++) {
	I = aocc_off[Gi] + i;
	for(a=0; a < avirtpi[Ga]; a++) {
	  A = avir_off[Ga] + a;

	  for(m=0; m < boccpi[Gm]; m++) {
	    M = bocc_off[Gm] + m;

	    AM = E_AB.params->rowidx[A][M];
	    IM = E_AB.params->colidx[I][M];

	    fIA.matrix[Gi][i][a] += E_AB.matrix[h][AM][IM];

	  }
	}
      }

    }

    dpd_buf4_mat_irrep_close(&E_AA, h);
    dpd_buf4_mat_irrep_close(&E_AB, h);
  }

  dpd_buf4_close(&E_AA);
  dpd_buf4_close(&E_AB);

  dpd_buf4_init(&E_BB, CC_EINTS, 0, 31, 10, 31, 10, 1, "E <ai|jk>");
  dpd_buf4_init(&E_AB, CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");

  for(h=0; h < nirreps; h++) {

    dpd_buf4_mat_irrep_init(&E_BB, h);
    dpd_buf4_mat_irrep_init(&E_AB, h);
    dpd_buf4_mat_irrep_rd(&E_BB, h);
    dpd_buf4_mat_irrep_rd(&E_AB, h);

    /* Loop over irreps of the target */
    for(Gi=0; Gi < nirreps; Gi++) {
      Ga = Gi; Gm = h^Gi;

      /* Loop over orbitals of the target */
      for(i=0; i < boccpi[Gi]; i++) {
	I = bocc_off[Gi] + i;
	for(a=0; a < bvirtpi[Ga]; a++) {
	  A = bvir_off[Ga] + a;

	  for(m=0; m < boccpi[Gm]; m++) {
	    M = bocc_off[Gm] + m;

	    AM = E_BB.params->rowidx[A][M];
	    IM = E_BB.params->colidx[I][M];

	    fia.matrix[Gi][i][a] += E_BB.matrix[h][AM][IM];

	  }
	}
      }
   
      /* Loop over orbitals of the target */
      for(i=0; i < boccpi[Gi]; i++) {
	I = bocc_off[Gi] + i;
	for(a=0; a < bvirtpi[Ga]; a++) {
	  A = bvir_off[Ga] + a;

	  for(m=0; m < aoccpi[Gm]; m++) {
	    M = aocc_off[Gm] + m;

	    MI = E_AB.params->rowidx[M][I];
	    MA = E_AB.params->colidx[M][A];

	    fia.matrix[Gi][i][a] += E_AB.matrix[h][MI][MA];

	  }
	}
      }
    }

    dpd_buf4_mat_irrep_close(&E_BB, h);
    dpd_buf4_mat_irrep_close(&E_AB, h);

  }

  /* Close the E integral buffers */
  dpd_buf4_close(&E_BB);
  dpd_buf4_close(&E_AB);

  /* Close the alpha and beta occ-vir Fock matrix files */
  dpd_file2_mat_wrt(&fIA);
  dpd_file2_mat_wrt(&fia);
  dpd_file2_mat_close(&fIA);
  dpd_file2_mat_close(&fia);
  dpd_file2_close(&fIA);
  dpd_file2_close(&fia);

  dpd_file2_mat_close(&hIJ);
  dpd_file2_mat_close(&hij);
  dpd_file2_mat_close(&hAB);
  dpd_file2_mat_close(&hab);
  dpd_file2_mat_close(&hIA);
  dpd_file2_mat_close(&hia);

  dpd_file2_close(&hIJ);
  dpd_file2_close(&hij);
  dpd_file2_close(&hAB);
  dpd_file2_close(&hab);
  dpd_file2_close(&hIA);
  dpd_file2_close(&hia);
}

void fock_rhf(void)
{
  int h, Gi, Gj, Ga, Gb, Gm;
  int i,j,a,b,m;
  int I, J, A, B, M;
  int IM, JM, MA, MB, AM;
  int nirreps;
  int *occpi, *virtpi;
  int *occ_off, *vir_off;
  int *occ_sym, *vir_sym;
  int *openpi;
  dpdfile2 fIJ, fij, fAB, fab, fIA, fia, Hoo, Hvv, Hov;
  dpdbuf4 AInts_anti, AInts, CInts, CInts_anti, EInts_anti, EInts;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;
  occ_off = moinfo.occ_off; vir_off = moinfo.vir_off;
  occ_sym = moinfo.occ_sym; vir_sym = moinfo.vir_sym;
  openpi = moinfo.openpi;

  dpd_file2_init(&Hoo, CC_OEI, 0, 0, 0, "h(i,j)");
  dpd_file2_init(&Hvv, CC_OEI, 0, 1, 1, "h(a,b)");
  dpd_file2_init(&Hov, CC_OEI, 0, 0, 1, "h(i,a)");
  dpd_file2_mat_init(&Hoo);
  dpd_file2_mat_init(&Hvv);
  dpd_file2_mat_init(&Hov);
  dpd_file2_mat_rd(&Hoo);
  dpd_file2_mat_rd(&Hvv);
  dpd_file2_mat_rd(&Hov);
  
  /* Prepare the alpha and beta occ-occ Fock matrix files */
  dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_init(&fij, CC_OEI, 0, 0, 0, "fij");
  dpd_file2_mat_init(&fIJ);
  dpd_file2_mat_init(&fij);

  /* One-electron (frozen-core) contributions */
  for(h=0; h < nirreps; h++) {

      for(i=0; i < occpi[h]; i++) 
          for(j=0; j < occpi[h]; j++) 
              fIJ.matrix[h][i][j] = Hoo.matrix[h][i][j];   

      for(i=0; i < (occpi[h]-openpi[h]); i++) 
          for(j=0; j < (occpi[h]-openpi[h]); j++) 
              fij.matrix[h][i][j] = Hoo.matrix[h][i][j];   
    }

  /* Two-electron contributions */

  /* Prepare the A integral buffers */
  dpd_buf4_init(&AInts_anti, CC_AINTS, 0, 0, 0, 0, 0, 1, "A <ij|kl>");
  dpd_buf4_init(&AInts, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");

  for(h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(&AInts_anti, h);
      dpd_buf4_mat_irrep_rd(&AInts_anti, h);

      /* Loop over irreps of the target */
      for(Gi=0; Gi < nirreps; Gi++) {
          Gj=Gi; Gm=h^Gi;

          /* Loop over the orbitals of the target */
          for(i=0; i < occpi[Gi]; i++) {
              I = occ_off[Gi] + i;
              for(j=0; j < occpi[Gj]; j++) {
                  J = occ_off[Gj] + j;
                  for(m=0; m < occpi[Gm]; m++) {
                      M = occ_off[Gm] + m;

                      IM = AInts_anti.params->rowidx[I][M];
                      JM = AInts_anti.params->colidx[J][M];

                      fIJ.matrix[Gi][i][j] += AInts_anti.matrix[h][IM][JM];

                    }
                }
            }

          /* Loop over the orbitals of the target */
          for(i=0; i < (occpi[Gi] - openpi[Gi]); i++) {
              I = occ_off[Gi] + i;
              for(j=0; j < (occpi[Gj] - openpi[Gj]); j++) {
                  J = occ_off[Gj] + j;
                  for(m=0; m < (occpi[Gm] - openpi[Gm]); m++) {
                      M = occ_off[Gm] + m;

                      IM = AInts_anti.params->rowidx[I][M];
                      JM = AInts_anti.params->colidx[J][M];

                      fij.matrix[Gi][i][j] += AInts_anti.matrix[h][IM][JM];

                    }
                }
            }

        }
      
      dpd_buf4_mat_irrep_close(&AInts_anti, h);

      dpd_buf4_mat_irrep_init(&AInts, h);
      dpd_buf4_mat_irrep_rd(&AInts, h);

      /* Loop over irreps of the target */
      for(Gi=0; Gi < nirreps; Gi++) {
          Gj=Gi; Gm=h^Gi;

          /* Loop over the orbitals of the target */
          for(i=0; i < occpi[Gi]; i++) {
              I = occ_off[Gi] + i;
              for(j=0; j < occpi[Gj]; j++) {
                  J = occ_off[Gj] + j;
                  for(m=0; m < (occpi[Gm] - openpi[Gm]); m++) {
                      M = occ_off[Gm] + m;

                      IM = AInts.params->rowidx[I][M];
                      JM = AInts.params->colidx[J][M];

                      fIJ.matrix[Gi][i][j] += AInts.matrix[h][IM][JM];

                    }
                }
            }

          /* Loop over the orbitals of the target */
          for(i=0; i < (occpi[Gi] - openpi[Gi]); i++) {
              I = occ_off[Gi] + i;
              for(j=0; j < (occpi[Gj] - openpi[Gj]); j++) {
                  J = occ_off[Gj] + j;
                  for(m=0; m < occpi[Gm]; m++) {
                      M = occ_off[Gm] + m;

                      IM = AInts.params->rowidx[I][M];
                      JM = AInts.params->colidx[J][M];

                      fij.matrix[Gi][i][j] += AInts.matrix[h][IM][JM];

                    }
                }
            }

        }
      
      dpd_buf4_mat_irrep_close(&AInts, h);

    }

  /* Close the A Integral buffers */
  dpd_buf4_close(&AInts_anti);
  dpd_buf4_close(&AInts);
  
  /* Close the alpha and beta occ-occ Fock matrix files */
  dpd_file2_mat_wrt(&fIJ);
  dpd_file2_mat_wrt(&fij);
  dpd_file2_mat_close(&fIJ);
  dpd_file2_mat_close(&fij);
  dpd_file2_close(&fIJ);
  dpd_file2_close(&fij);
  
  /* Prepare the alpha and beta vir-vir Fock matrix files */
  dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
  dpd_file2_init(&fab, CC_OEI, 0, 1, 1, "fab");
  dpd_file2_mat_init(&fAB);
  dpd_file2_mat_init(&fab);

  /* One-electron (frozen-core) contributions */
  for(h=0; h < nirreps; h++) {

      for(a=0; a < (virtpi[h] - openpi[h]); a++) 
          for(b=0; b < (virtpi[h] - openpi[h]); b++) 
              fAB.matrix[h][a][b] = Hvv.matrix[h][a][b];   

      for(a=0; a < virtpi[h]; a++) 
          for(b=0; b < virtpi[h]; b++) 
              fab.matrix[h][a][b] = Hvv.matrix[h][a][b];   
    }

  /* Two-electron contributions */

  /* Prepare the C integral buffers */
  dpd_buf4_init(&CInts_anti, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
  dpd_buf4_init(&CInts, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");

  for(h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(&CInts_anti, h);
      dpd_buf4_mat_irrep_rd(&CInts_anti, h);

      /* Loop over irreps of the target */
      for(Ga=0; Ga < nirreps; Ga++) {
	  Gb = Ga; Gm = h^Ga;

	  /* Loop over orbitals of the target */
	  for(a=0; a < (virtpi[Ga] - openpi[Ga]); a++) {
	      A = vir_off[Ga] + a;
	      for(b=0; b < (virtpi[Gb] - openpi[Gb]); b++) {
		  B = vir_off[Gb] + b;

		  for(m=0; m < occpi[Gm]; m++) {
		      M = occ_off[Gm] + m;

		      MA = CInts_anti.params->rowidx[M][A];
		      MB = CInts_anti.params->colidx[M][B];

		      fAB.matrix[Ga][a][b] += CInts_anti.matrix[h][MA][MB];

		    }
		}
	    }

	  /* Loop over orbitals of the target */
	  for(a=0; a < virtpi[Ga]; a++) {
	      A = vir_off[Ga] + a;
	      for(b=0; b < virtpi[Gb]; b++) {
		  B = vir_off[Gb] + b;

		  for(m=0; m < (occpi[Gm] - openpi[Gm]); m++) {
		      M = occ_off[Gm] + m;

		      MA = CInts_anti.params->rowidx[M][A];
		      MB = CInts_anti.params->colidx[M][B];

		      fab.matrix[Ga][a][b] += CInts_anti.matrix[h][MA][MB];
		    }
		}
	    }
	}

      dpd_buf4_mat_irrep_close(&CInts_anti, h);

      dpd_buf4_mat_irrep_init(&CInts, h);
      dpd_buf4_mat_irrep_rd(&CInts, h);

      /* Loop over irreps of the target */
      for(Ga=0; Ga < nirreps; Ga++) {
	  Gb = Ga; Gm = h^Ga;

	  /* Loop over orbitals of the target */
	  for(a=0; a < (virtpi[Ga] - openpi[Ga]); a++) {
	      A = vir_off[Ga] + a;
	      for(b=0; b < (virtpi[Gb] - openpi[Gb]); b++) {
		  B = vir_off[Gb] + b;

		  for(m=0; m < (occpi[Gm] - openpi[Gm]); m++) {
		      M = occ_off[Gm] + m;

		      MA = CInts.params->rowidx[M][A];
		      MB = CInts.params->colidx[M][B];

		      fAB.matrix[Ga][a][b] += CInts.matrix[h][MA][MB];

		    }
		}
	    }

	  /* Loop over orbitals of the target */
	  for(a=0; a < virtpi[Ga]; a++) {
	      A = vir_off[Ga] + a;
	      for(b=0; b < virtpi[Gb]; b++) {
		  B = vir_off[Gb] + b;

		  for(m=0; m < occpi[Gm]; m++) {
		      M = occ_off[Gm] + m;

		      MA = CInts.params->rowidx[M][A];
		      MB = CInts.params->colidx[M][B];

		      fab.matrix[Ga][a][b] += CInts.matrix[h][MA][MB];
		    }
		}
	    }
	}

      dpd_buf4_mat_irrep_close(&CInts, h);
    }

  /* Close the C integral buffers */
  dpd_buf4_close(&CInts_anti);
  dpd_buf4_close(&CInts);

  /* Close the alpha and beta vir-vir Fock matrix files */
  dpd_file2_mat_wrt(&fAB);
  dpd_file2_mat_wrt(&fab);
  dpd_file2_mat_close(&fAB);
  dpd_file2_mat_close(&fab);
  dpd_file2_close(&fAB);
  dpd_file2_close(&fab);
  
  /* Prepare the alpha and beta occ-vir Fock matrix files */
  dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
  dpd_file2_init(&fia, CC_OEI, 0, 0, 1, "fia");
  dpd_file2_mat_init(&fIA);
  dpd_file2_mat_init(&fia);

  /* One-electron (frozen-core) contributions */
  for(h=0; h < nirreps; h++) {

      for(i=0; i < occpi[h]; i++) 
          for(a=0; a < (virtpi[h] - openpi[h]); a++) 
              fIA.matrix[h][i][a] = Hov.matrix[h][i][a];   

      for(i=0; i < (occpi[h] - openpi[h]); i++) 
          for(a=0; a < virtpi[h]; a++) 
              fia.matrix[h][i][a] = Hov.matrix[h][i][a];   
    }

  /* Two-electron contributions */

  /* Prepare the E integral buffers */
  dpd_buf4_init(&EInts_anti, CC_EINTS, 0, 11, 0, 11, 0, 1, "E <ai|jk>");
  dpd_buf4_init(&EInts, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");

  for(h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(&EInts_anti, h);
      dpd_buf4_mat_irrep_rd(&EInts_anti, h);

      /* Loop over irreps of the target */
      for(Gi=0; Gi < nirreps; Gi++) {
	  Ga = Gi; Gm = h^Gi;

	  /* Loop over orbitals of the target */
	  for(i=0; i < occpi[Gi]; i++) {
	      I = occ_off[Gi] + i;
	      for(a=0; a < (virtpi[Ga] - openpi[Ga]); a++) {
		  A = vir_off[Ga] + a;

		  for(m=0; m < occpi[Gm]; m++) {
		      M = occ_off[Gm] + m;

		      AM = EInts_anti.params->rowidx[A][M];
		      IM = EInts_anti.params->colidx[I][M];

		      fIA.matrix[Gi][i][a] += EInts_anti.matrix[h][AM][IM];

		    }
		}
	    }

	  /* Loop over orbitals of the target */
	  for(i=0; i < (occpi[Gi] - openpi[Gi]); i++) {
	      I = occ_off[Gi] + i;
	      for(a=0; a < virtpi[Ga]; a++) {
		  A = vir_off[Ga] + a;

		  for(m=0; m < (occpi[Gm] - openpi[Gm]); m++) {
		      M = occ_off[Gm] + m;

		      AM = EInts_anti.params->rowidx[A][M];
		      IM = EInts_anti.params->colidx[I][M];

		      fia.matrix[Gi][i][a] += EInts_anti.matrix[h][AM][IM];

		    }
		}
	    }
	}

      dpd_buf4_mat_irrep_close(&EInts_anti, h);

      dpd_buf4_mat_irrep_init(&EInts, h);
      dpd_buf4_mat_irrep_rd(&EInts, h);

      /* Loop over irreps of the target */
      for(Gi=0; Gi < nirreps; Gi++) {
	  Ga = Gi; Gm = h^Gi;

	  /* Loop over orbitals of the target */
	  for(i=0; i < occpi[Gi]; i++) {
	      I = occ_off[Gi] + i;
	      for(a=0; a < (virtpi[Ga] - openpi[Ga]); a++) {
		  A = vir_off[Ga] + a;

		  for(m=0; m < (occpi[Gm] - openpi[Gm]); m++) {
		      M = occ_off[Gm] + m;

		      AM = EInts.params->rowidx[A][M];
		      IM = EInts.params->colidx[I][M];

		      fIA.matrix[Gi][i][a] += EInts.matrix[h][AM][IM];

		    }
		}
	    }

	  /* Loop over orbitals of the target */
	  for(i=0; i < (occpi[Gi] - openpi[Gi]); i++) {
	      I = occ_off[Gi] + i;
	      for(a=0; a < virtpi[Ga]; a++) {
		  A = vir_off[Ga] + a;

		  for(m=0; m < occpi[Gm]; m++) {
		      M = occ_off[Gm] + m;

		      AM = EInts.params->rowidx[A][M];
		      IM = EInts.params->colidx[I][M];

		      fia.matrix[Gi][i][a] += EInts.matrix[h][AM][IM];

		    }
		}
	    }
	}

      dpd_buf4_mat_irrep_close(&EInts, h);

    }

  /* Close the E integral buffers */
  dpd_buf4_close(&EInts_anti);
  dpd_buf4_close(&EInts);

  /* Close the alpha and beta occ-vir Fock matrix files */
  dpd_file2_mat_wrt(&fIA);
  dpd_file2_mat_wrt(&fia);
  dpd_file2_mat_close(&fIA);
  dpd_file2_mat_close(&fia);
  dpd_file2_close(&fIA);
  dpd_file2_close(&fia);

  dpd_file2_mat_close(&Hoo);
  dpd_file2_mat_close(&Hvv);
  dpd_file2_mat_close(&Hov);
  dpd_file2_close(&Hoo);
  dpd_file2_close(&Hvv);
  dpd_file2_close(&Hov);

}

}} // namespace psi::ccsort
