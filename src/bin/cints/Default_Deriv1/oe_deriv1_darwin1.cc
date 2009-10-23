/*! \file
  \ingroup CINTS
  \brief Compute the derivative of the one-electron Darwin relativistic correction.
  \author R.A. King
*/
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <libciomr/libciomr.h>
#include <libint/libint.h>
#include <libderiv/libderiv.h>
#include <libqt/qt.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <psifiles.h>
#include <physconst.h>
#include <psiconfig.h>

#include "defines.h"
#define EXTERN
#include "global.h"

namespace psi {
  namespace cints {

    static const int use_cca_integrals_standard = (PSI_INTEGRALS_STANDARD == 1);

void oe_deriv1_darwin1(void) {
  double **AO_at_nuc, contr, ax, ay, az, ai, energy, *energy_atom;
  double zval, **darwin_deriv, **darwin_deriv_orb, **norm_bf;
  double  **temp_mat, **der_rho, Rab_x, Rab_y, Rab_z, Rab, sign, ***dAO_at_nuc, tval, ***UX;
  int natom, nao, max_am, **xpow_bf, **ypow_bf, **zpow_bf, xpow, ypow, zpow ;
  int i, j, k, l, m, ibf, iatom, iang, nbfi, ig, igmin, igmax, i_1stbf, coord, occ_offset;
  int atom_A, atom_B, atom_C, xyz, nmo, mo, mo_max, pitz_offset, qts_offset, irrep, cnt;
  char *label;

  natom = Molecule.num_atoms;
  nao = BasisSet.num_ao;
  max_am = BasisSet.max_am;
  nmo = MOInfo.num_mo;

  /* determine powers of x,y,z for each AO, and normalization constants */
  xpow_bf = init_int_matrix(max_am+1,(max_am+1)*(max_am+2)/2);
  ypow_bf = init_int_matrix(max_am+1,(max_am+1)*(max_am+2)/2);
  zpow_bf = init_int_matrix(max_am+1,(max_am+1)*(max_am+2)/2);
  norm_bf = block_matrix(max_am+1,(max_am+1)*(max_am+2)/2);
  for(l=0;l<=max_am;l++) {
    ibf = 0;
    for(i = 0; i <= l; i++){
      for(j = 0; j <= i; j++){
        xpow_bf[l][ibf] = l - i;
        ypow_bf[l][ibf] = i - j;
        zpow_bf[l][ibf] = j;
        /* norm_bf is an additional normalization factor for gaussians with l>=d */
        norm_bf[l][ibf] = use_cca_integrals_standard ? 1.0 : sqrt(df[2*l]/(df[2*(l-i)]*df[2*(i-j)]*df[2*j]));
        ibf++;
      }
    }
  }

  /* Compute values of basis functions at each nucleus */
  AO_at_nuc = block_matrix(natom,nao);

  for(i=0;i<BasisSet.num_shells;i++) {      /* loop over shells of bf's */
    iatom = BasisSet.shells[i].center - 1;  /* lookup atom number for shell */
    ax = Molecule.centers[iatom].x;         /* lookup coordinates of atom */
    ay = Molecule.centers[iatom].y;     
    az = Molecule.centers[iatom].z;
    iang  = BasisSet.shells[i].am - 1;      /* iang = l for this shell */
    nbfi = (iang+2)*(iang+1)/2;             /* nbfi = # of bf in shell */
    igmin = BasisSet.shells[i].fprim - 1;   /* locate first primitive in shell */
    i_1stbf = BasisSet.shells[i].fao - 1;   /* locate 1st AO in shell */
    igmax = igmin + BasisSet.shells[i].n_prims - 1; /* last primitive = 1st + # primitives - 1 */

    for(ibf=0;ibf<nbfi;ibf++) {             /* loop over contracted gaussians */
      xpow = xpow_bf[iang][ibf];
      ypow = ypow_bf[iang][ibf];
      zpow = zpow_bf[iang][ibf];
      for(ig=igmin;ig<=igmax;ig++) {        /* loop over primitive functions */
        ai    = BasisSet.cgtos[ig].exp;
        contr = BasisSet.cgtos[ig].ccoeff[iang];
          
        for(k=0;k<natom;k++) {              /* loop over nuclei at which to evaluate functions */
          Rab_x = Molecule.centers[k].x - ax;
          Rab_y = Molecule.centers[k].y - ay;
          Rab_z = Molecule.centers[k].z - az;
          Rab = Rab_x*Rab_x + Rab_y*Rab_y + Rab_z*Rab_z;

          tval = contr * norm_bf[iang][ibf] * exp(-ai * Rab); /* avoid 0^0 */
          for (m=0; m<xpow; ++m)
            tval *= Rab_x;
          for (m=0; m<ypow; ++m)
            tval *= Rab_y;
          for (m=0; m<zpow; ++m)
            tval *= Rab_z;
          AO_at_nuc[k][i_1stbf+ibf] += tval;
        }
      }
    }
  }
  if (UserOptions.print_lvl >= PRINT_DEBUG) {
    fprintf(outfile,"Values of AO's at nuclei\n");
    print_mat(AO_at_nuc,natom,nao,outfile);
  }

  /* Compute Darwin energy */
  energy = 0.0;
  energy_atom = init_array(natom);
  for (iatom=0; iatom<natom; ++iatom) {
    zval = Molecule.centers[iatom].Z_nuc;
    for (i=0; i<nao; ++i)
      for (j=0; j<nao; ++j)
        energy_atom[iatom] += zval * Dens[i][j] * AO_at_nuc[iatom][i] * AO_at_nuc[iatom][j];
  }
  for (iatom=0; iatom<natom; ++iatom) {
    energy_atom[iatom] *= UserOptions.fine_structure_alpha * UserOptions.fine_structure_alpha
      * _pi/(2*_c_au*_c_au);
    energy += energy_atom[iatom];
  }
  if (UserOptions.print_lvl >= PRINT_DEBUG) {
    fprintf(outfile,"One-electron Darwin term per atom\n");
    for (iatom=0; iatom<natom; ++iatom)
      fprintf(outfile,"\tAtom %d: %15.10lf\n",iatom,energy_atom[iatom]);
  }
  free(energy_atom);

  fprintf(outfile,"\n\tOne-electron Darwin energy: %15.10lf\n", energy);
  energy = energy + chkpt_rd_etot();
  chkpt_wt_etot(energy);

  /* Compute AO first-derivative of Darwin integrals */
  /* d[phi_A(R_B)]/dR_Cx = 0, if (A==B) or (A!=B!=C) */
  dAO_at_nuc = init_3d_array(3*natom,natom,nao);
  for(i=0; i<BasisSet.num_shells; i++) {     
    atom_A = BasisSet.shells[i].center - 1;
    ax = Molecule.centers[atom_A].x;       
    ay = Molecule.centers[atom_A].y;
    az = Molecule.centers[atom_A].z;
    iang  = BasisSet.shells[i].am - 1;    
    nbfi = (iang+2)*(iang+1)/2;           
    i_1stbf = BasisSet.shells[i].fao - 1; 
    igmin = BasisSet.shells[i].fprim - 1;   /* locate first primitive in shell */
    igmax = igmin + BasisSet.shells[i].n_prims - 1; /* last primitive = 1st + # primitives - 1 */

    for (atom_B=0; atom_B<natom; ++atom_B) {
      if (atom_B == atom_A) continue;
      Rab_x = Molecule.centers[atom_B].x - ax;
      Rab_y = Molecule.centers[atom_B].y - ay;
      Rab_z = Molecule.centers[atom_B].z - az;
      Rab = Rab_x*Rab_x + Rab_y*Rab_y + Rab_z*Rab_z;

      for(ibf=0;ibf<nbfi;ibf++) {             /* loop over contracted gaussians */
        xpow = xpow_bf[iang][ibf];
        ypow = ypow_bf[iang][ibf];
        zpow = zpow_bf[iang][ibf];
        for(ig=igmin;ig<=igmax;ig++) {        /* loop over primitive functions */
          ai    = BasisSet.cgtos[ig].exp;
          contr = BasisSet.cgtos[ig].ccoeff[iang];

          for (atom_C=0; atom_C<natom; ++atom_C) {  /* Take 3N derivatives wrt R_Cxyz */
            if ((atom_A != atom_B) && (atom_B != atom_C) && (atom_A != atom_C)) continue;

              /* (-1 if B==C) * 2 beta (Rbx-Rax) * phi_A(R_B) */
            if (atom_B==atom_C) sign = -1; else sign = 1;
            tval = 2.0 * sign * ai * contr * norm_bf[iang][ibf] * exp(-ai * Rab); /* avoid 0^0 */
            if (xpow) tval *= pow(Rab_x, xpow);
            if (ypow) tval *= pow(Rab_y, ypow);
            if (zpow) tval *= pow(Rab_z, zpow);
            dAO_at_nuc[3*atom_C+0][atom_B][i_1stbf+ibf] += Rab_x * tval;
            dAO_at_nuc[3*atom_C+1][atom_B][i_1stbf+ibf] += Rab_y * tval;
            dAO_at_nuc[3*atom_C+2][atom_B][i_1stbf+ibf] += Rab_z * tval;

             /* (-1 if A==C) * l_x / (Rbx-Rax) * phi_A(R_B) */
            if (atom_A==atom_C) sign = -1; else sign = 1;
            if (xpow > 0) {
              tval = sign * contr * norm_bf[iang][ibf] * exp(-ai * Rab);
              if (xpow >= 2) tval *= xpow * pow(Rab_x, xpow-1);
              if (ypow) tval *= pow(Rab_y, ypow);
              if (zpow) tval *= pow(Rab_z, zpow);
              dAO_at_nuc[3*atom_C+0][atom_B][i_1stbf+ibf] += tval;
            }
            if (ypow > 0) {
              tval = sign * contr * norm_bf[iang][ibf] * exp(-ai * Rab);
              if (xpow) tval *= pow(Rab_x, xpow);
              if (ypow >= 2) tval *= ypow * pow(Rab_y, ypow-1);
              if (zpow) tval *= pow(Rab_z, zpow);
              dAO_at_nuc[3*atom_C+1][atom_B][i_1stbf+ibf] += tval;
            }
            if (zpow > 0) {
              tval = sign * contr * norm_bf[iang][ibf] * exp(-ai * Rab);
              if (xpow) tval *= pow(Rab_x, xpow);
              if (ypow) tval *= pow(Rab_y, ypow);
              if (zpow >= 2) tval *= zpow * pow(Rab_z, zpow-1);
              dAO_at_nuc[3*atom_C+2][atom_B][i_1stbf+ibf] += tval;
            }
          } /* end loop atom_C */
        } /* end loop ig */
      } /* end loop ibf */ 
    } /* end atom_B */
  } /* end shells */

  if (UserOptions.print_lvl >= PRINT_DEBUG) {
    for (i=0; i<3*natom; ++i) {
      fprintf(outfile,"\n\tOne-electron Darwin integral derivative %d\n",i);
      print_mat(dAO_at_nuc[i],natom,nao,outfile);
    } 
    for (i=0; i<3*natom; ++i) {
      tval = 0.0;
      for (j=0; j<natom; ++j)
        for (k=0; k<nao; ++k)
          tval += dAO_at_nuc[i][j][k];
      fprintf(outfile,"Sum of dAO_at_nuc for derivative %d (analytic): %15.10lf\n", i, tval); 
    }
  }

  /* compute integral derivative part of 1e Darwin gradient 
     = pi/(2c^2) * sum_k Z_k sum_uv Duv [phi_u(k)*dphi_v(k)/dRax+dphi_u(k)/dRax*phi_v(k)]
     = pi/(c^2) * sum_k Z_k sum_uv Duv [phi_u(k)*dphi_v(k)] */
  darwin_deriv = block_matrix(natom,3);
  for (iatom=0; iatom<natom; ++iatom) {
    for (xyz=0; xyz<3; ++xyz) {
      tval = 0.0;
      for (k=0; k<natom; ++k) {
        zval = Molecule.centers[k].Z_nuc;
        for (i=0; i<nao; ++i)
          for (j=0; j<nao; ++j)
            tval += zval * Dens[i][j] * AO_at_nuc[k][i] * dAO_at_nuc[3*iatom+xyz][k][j]; 
      }
      darwin_deriv[iatom][xyz] = UserOptions.fine_structure_alpha * UserOptions.fine_structure_alpha 
        * _pi/(_c_au * _c_au) * tval;
    }
  }

  /* compute orbital derivative part of 1e Darwin gradient */
  UX = (double ***) malloc(3*natom*sizeof(double **));
  for (i=0;i<3*natom;++i)
    UX[i] = block_matrix(nmo,nmo);
  label = (char *) malloc(PSIO_KEYLEN * sizeof(char));
  psio_open(PSIF_CPHF,1);
  for (i=0;i<3*natom;++i) {
    for(j=0; j<PSIO_KEYLEN; j++) label[j] = '\0';
    sprintf(label, "UX(%d)",i);
    psio_read_entry(PSIF_CPHF, label, (char *) &(UX[i][0][0]), nmo*nmo*sizeof(double));
  }
  psio_close(PSIF_CPHF,1);
  free(label);

  temp_mat = block_matrix(nmo,nao);
  der_rho = block_matrix(nao,nao);
  darwin_deriv_orb = block_matrix(natom,3);
  for (iatom=0; iatom<natom; ++iatom) {
    for (xyz=0; xyz<3; ++xyz) {

      mmult(UX[3*iatom+xyz], 1, MOInfo.scf_evec[0], 0,  temp_mat, 0, nmo, nmo, nao, 0);

      zero_mat(der_rho,nao,nao);

      for (i=0; i<nao; ++i) {
        for (j=0; j<nao; ++j) {
          occ_offset = 0;
          cnt = -1;
          for(irrep=0; irrep<Symmetry.nirreps; irrep++) {
            for(m=0; m<MOInfo.clsdpi[irrep]; m++)
              der_rho[i][j] += 2.0 * MOInfo.scf_evec[0][m+occ_offset][i]*temp_mat[m+occ_offset][j];
            occ_offset += MOInfo.orbspi[irrep];
          }
        }
      }

      /* = pi/(2c^2) sum_k Z_k sum_uv [C_u^m d(C_v^m)/dRax + d(C_u^m)/dRax C_v^m]  phi_u(k) phi_v(k) */
      /* = pi/(c^2) sum_k Z_k sum_uv [C_u^m d(C_v^m)/dRax]  phi_u(k) phi_v(k) */
      tval = 0.0;
      for (k=0; k<natom; ++k) {
        zval = Molecule.centers[k].Z_nuc;
        for (i=0; i<nao; ++i)
          for (j=0; j<nao; ++j)
            tval += zval * der_rho[i][j] * AO_at_nuc[k][i] * AO_at_nuc[k][j];
      }
      darwin_deriv_orb[iatom][xyz] = UserOptions.fine_structure_alpha * UserOptions.fine_structure_alpha
        * _pi/(_c_au * _c_au) * tval;
    }
  }

  if (UserOptions.print_lvl >= PRINT_DEBUG) {
    fprintf(outfile,"1e Darwin gradient, integral derivative part\n");
    print_mat(darwin_deriv,natom,3,outfile);
    fprintf(outfile,"1e Darwin gradient, orbital derivative part\n");
    print_mat(darwin_deriv_orb,natom,3,outfile);
  }

  for (iatom=0; iatom<natom; ++iatom)
    for (xyz=0; xyz<3; ++xyz)
      darwin_deriv[iatom][xyz] += darwin_deriv_orb[iatom][xyz];

  if (UserOptions.print_lvl >= PRINT_DEBUG) {
    fprintf(outfile,"\tOne-electron Darwin gradient, total\n");
    print_mat(darwin_deriv,natom,3,outfile);
  }

  for (iatom=0; iatom<natom; ++iatom)
    for (xyz=0; xyz<3; ++xyz)
      Grad[iatom][xyz] += darwin_deriv[iatom][xyz];

  free_block(darwin_deriv_orb);
  free_block(temp_mat);
  free_block(der_rho);
  free_block(darwin_deriv);
  for (i=0;i<3*natom;++i)
    free_block(UX[i]);
  free(UX);

  free_block(norm_bf);
  free_block(AO_at_nuc);
  free_3d_array(dAO_at_nuc,3*natom,natom);

  free_int_matrix(xpow_bf);
  free_int_matrix(ypow_bf);
  free_int_matrix(zpow_bf);

  /* to compute derivative integrals numerically for testing use */
  /* oe_deriv1_darwin1_test(); */
}

}}
