/*! \file
  \ingroup CINTS
  Compute the derivative of the mass-velocity relativistic correction.
  Energy code adapted from oeprop by E.F. Valeev.
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
#include <psifiles.h>
#include <physconst.h>
#include <psiconfig.h>

#include "defines.h"
#define EXTERN
#include "global.h"

#define TEST_OE_DERIV1_MVC_DISP (0.00001)

namespace psi {
  namespace cints {

    static const int use_cca_integrals_standard = (PSI_INTEGRALS_STANDARD == 1);

extern double ***init_box(int a, int b, int c);
extern void free_box(double ***box, int a, int b);
extern void MI_OSrecurs(double pax, double pay, double paz,
                 double pbx, double pby, double pbz, double gamma,
                 int lmaxi, int lmaxj, int maxm);
extern double ***MIX;
extern double ***MIY;
extern double ***MIZ;

void oe_deriv1_mvc_test(void) {
  double  contr, ax, ay, az, bx, by, bz, ai, bj,  energy, **norm_bf;
  double  Rab2, ***UX, tmp, over_pf, ***dmvc_ints, ***pmvc_ints, ***mmvc_ints;
  double x0, y0, z0, tx, ty, tz, t2x, t2y, t2z, gamma, tx1, ty1, tz1;
  double pax, pay, paz, pbx, pby, pbz, px, py, pz, contr_i, contr_j, tval;
  int natom, nao, max_am, **xpow_bf, **ypow_bf, **zpow_bf, xpow, ypow, zpow ;
  int i, j, k, l, m, coord, occ_offset, eq_shell, eq_atoms;
  int lx1, ly1, lz1, lx2, ly2, lz2, atom_C;
  int ibf, iatom, iang, nbfi, ig, igmin, igmax, i_1stbf;
  int jbf, jatom, jang, nbfj, jg, jgmin, jgmax, j_1stbf;
  int xyz, nmo, mo, mo_max, irrep, cnt;
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

  MIX = init_box(max_am+3,max_am+3,0);
  MIY = init_box(max_am+3,max_am+3,0);
  MIZ = init_box(max_am+3,max_am+3,0);
  MIX[0][0][0] = MIY[0][0][0] = MIZ[0][0][0] = 1.0;

  pmvc_ints = init_3d_array(3*natom,nao,nao);
  mmvc_ints = init_3d_array(3*natom,nao,nao);
  dmvc_ints = init_3d_array(3*natom,nao,nao);

  /* compute mass-velocity integrals at +R_Cx*/
  for (atom_C=0; atom_C<natom; ++atom_C) {
    for (xyz=0; xyz<3; ++xyz) {

  for(i=0;i<BasisSet.num_shells;i++) {  /* I-shell loop */ 
    iatom = BasisSet.shells[i].center - 1;
    if (iatom != atom_C) continue;
    ax = Molecule.centers[iatom].x;
    ay = Molecule.centers[iatom].y;
    az = Molecule.centers[iatom].z;

    if (xyz==0) ax += TEST_OE_DERIV1_MVC_DISP;
    else if (xyz==1) ay += TEST_OE_DERIV1_MVC_DISP;
    else if (xyz==2) az += TEST_OE_DERIV1_MVC_DISP;

    iang  = BasisSet.shells[i].am - 1;
    nbfi = (iang+2)*(iang+1)/2;
    igmin = BasisSet.shells[i].fprim - 1;
    i_1stbf = BasisSet.shells[i].fao - 1;
    igmax = igmin + BasisSet.shells[i].n_prims - 1;

    for(j=0;j<BasisSet.num_shells;j++) {  /* J-shell loop */
      jatom = BasisSet.shells[j].center - 1;
      if (jatom == iatom) continue;
      bx = Molecule.centers[jatom].x;
      by = Molecule.centers[jatom].y;
      bz = Molecule.centers[jatom].z;
      jang  = BasisSet.shells[j].am - 1;
      nbfj = (jang+2)*(jang+1)/2;
      jgmin = BasisSet.shells[j].fprim - 1;
      j_1stbf = BasisSet.shells[j].fao - 1;
      jgmax = jgmin + BasisSet.shells[j].n_prims - 1;

      eq_shell = (i == j);
      eq_atoms = (iatom == jatom);

      Rab2 = (bx-ax)*(bx-ax) + (by-ay)*(by-ay) + (bz-az)*(bz-az);

      for(ig=igmin;ig<=igmax;ig++) { /* loop over primitives ig in shell I */
        ai = BasisSet.cgtos[ig].exp;
        contr_i = BasisSet.cgtos[ig].ccoeff[iang];

        for(jg=jgmin;jg<=jgmax;jg++) { /* loop over primitives jg in shell J */
          bj = BasisSet.cgtos[jg].exp;
          contr_j = BasisSet.cgtos[jg].ccoeff[jang];
          gamma = ai + bj;

          tmp = ai * bj * Rab2 / gamma;
          over_pf = contr_i * contr_j * exp(-tmp) * pow(M_PI/gamma,1.5);
          //if (eq_shell && ig != jg)
          //  over_pf *= 2.0;

          /* Computing moment integrals */
          px = (ax*ai + bx*bj)/gamma;
          py = (ay*ai + by*bj)/gamma;
          pz = (az*ai + bz*bj)/gamma;
          if (eq_atoms)
            pax = pay = paz = pbx = pby = pbz = 0.0;
          else {
            pax = px - ax;  pay = py - ay;  paz = pz - az;
            pbx = px - bx;  pby = py - by;  pbz = pz - bz;
          }
          MI_OSrecurs(pax,pay,paz,pbx,pby,pbz,gamma,iang+2,jang+2,0);

          for(ibf=0;ibf<nbfi;ibf++) { /* loop over gaussians ibf */
            lx1 = xpow_bf[iang][ibf];
            ly1 = ypow_bf[iang][ibf];
            lz1 = zpow_bf[iang][ibf];

            for(jbf=0;jbf<nbfj;jbf++) { /* loop over gaussians jbf */
              lx2 = xpow_bf[jang][jbf];
              ly2 = ypow_bf[jang][jbf];
              lz2 = zpow_bf[jang][jbf];

              x0 = MIX[lx1][lx2][0];
              y0 = MIY[ly1][ly2][0];
              z0 = MIZ[lz1][lz2][0];

              /* One-dimensional del-squared integrals */
              tx = bj*(2*lx2+1)*x0 - 2*bj*bj*MIX[lx1][lx2+2][0];
              if (lx2 >= 2)
                tx -= 0.5*lx2*(lx2-1)*MIX[lx1][lx2-2][0];
              ty = bj*(2*ly2+1)*y0 - 2*bj*bj*MIY[ly1][ly2+2][0];
              if (ly2 >= 2)
                ty -= 0.5*ly2*(ly2-1)*MIY[ly1][ly2-2][0];
              tz = bj*(2*lz2+1)*z0 - 2*bj*bj*MIZ[lz1][lz2+2][0];
              if (lz2 >= 2)
                tz -= 0.5*lz2*(lz2-1)*MIZ[lz1][lz2-2][0];

              t2x = 4*(2*lx1+1)*(2*lx2+1)*ai*bj*x0 -
                8*(2*lx1+1)*ai*bj*bj*MIX[lx1][lx2+2][0] -
                8*(2*lx2+1)*ai*ai*bj*MIX[lx1+2][lx2][0] +
                16*ai*ai*bj*bj*MIX[lx1+2][lx2+2][0];
              if (lx1 >= 2)
                t2x += 4*lx1*(lx1-1)*bj*bj*MIX[lx1-2][lx2+2][0] -
                  2*lx1*(lx1-1)*(2*lx2+1)*bj*MIX[lx1-2][lx2][0];
              if (lx2 >= 2)
                t2x += 4*lx2*(lx2-1)*ai*ai*MIX[lx1+2][lx2-2][0] -
                  2*(2*lx1+1)*lx2*(lx2-1)*ai*MIX[lx1][lx2-2][0];
              if ((lx1 >= 2) && (lx2 >= 2))
                t2x += lx1*(lx1-1)*lx2*(lx2-1)*MIX[lx1-2][lx2-2][0];
 
              t2y = 4*(2*ly1+1)*(2*ly2+1)*ai*bj*y0 -
                8*(2*ly1+1)*ai*bj*bj*MIY[ly1][ly2+2][0] -
                8*(2*ly2+1)*ai*ai*bj*MIY[ly1+2][ly2][0] +
                16*ai*ai*bj*bj*MIY[ly1+2][ly2+2][0];
              if (ly1 >= 2)
                t2y += 4*ly1*(ly1-1)*bj*bj*MIY[ly1-2][ly2+2][0] -
                  2*ly1*(ly1-1)*(2*ly2+1)*bj*MIY[ly1-2][ly2][0];
              if (ly2 >= 2)
                t2y += 4*ly2*(ly2-1)*ai*ai*MIY[ly1+2][ly2-2][0] -
                  2*(2*ly1+1)*ly2*(ly2-1)*ai*MIY[ly1][ly2-2][0];
              if ((ly1 >= 2) && (ly2 >= 2))
                t2y += ly1*(ly1-1)*ly2*(ly2-1)*MIY[ly1-2][ly2-2][0];
 
              t2z = 4*(2*lz1+1)*(2*lz2+1)*ai*bj*z0 -
                8*(2*lz1+1)*ai*bj*bj*MIZ[lz1][lz2+2][0] -
                8*(2*lz2+1)*ai*ai*bj*MIZ[lz1+2][lz2][0] +
                16*ai*ai*bj*bj*MIZ[lz1+2][lz2+2][0];
              if (lz1 >= 2)
                t2z += 4*lz1*(lz1-1)*bj*bj*MIZ[lz1-2][lz2+2][0] -
                  2*lz1*(lz1-1)*(2*lz2+1)*bj*MIZ[lz1-2][lz2][0];
              if (lz2 >= 2)
                t2z += 4*lz2*(lz2-1)*ai*ai*MIZ[lz1+2][lz2-2][0] -
                  2*(2*lz1+1)*lz2*(lz2-1)*ai*MIZ[lz1][lz2-2][0];
              if ((lz1 >= 2) && (lz2 >= 2))
                t2z += lz1*(lz1-1)*lz2*(lz2-1)*MIZ[lz1-2][lz2-2][0];

                pmvc_ints[3*atom_C+xyz][i_1stbf+ibf][j_1stbf+jbf] += norm_bf[iang][ibf] * norm_bf[jang][jbf] *
                  over_pf * (t2x*y0*z0 + x0*t2y*z0 + x0*y0*t2z + 8*tx*ty*z0 + 8*tx*y0*tz + 8*x0*ty*tz);
            } /* end loop gaussians jbf */
          } /* end loop gaussians ibf */
        } /* end loop jg primitives */
      } /* end loop ig primitives */
    } /* end loop over j shell */
  } /* end loop over i shell */

  }
  }

  /* compute mass-velocity integrals at -R_Cx*/
  for (atom_C=0; atom_C<natom; ++atom_C) {
    for (xyz=0; xyz<3; ++xyz) {

  for(i=0;i<BasisSet.num_shells;i++) {  /* I-shell loop */ 
    iatom = BasisSet.shells[i].center - 1;
    if (iatom != atom_C) continue;
    ax = Molecule.centers[iatom].x;
    ay = Molecule.centers[iatom].y;
    az = Molecule.centers[iatom].z;

    if (xyz==0) ax -= TEST_OE_DERIV1_MVC_DISP;
    else if (xyz==1) ay -= TEST_OE_DERIV1_MVC_DISP;
    else if (xyz==2) az -= TEST_OE_DERIV1_MVC_DISP;

    iang  = BasisSet.shells[i].am - 1;
    nbfi = (iang+2)*(iang+1)/2;
    igmin = BasisSet.shells[i].fprim - 1;
    i_1stbf = BasisSet.shells[i].fao - 1;
    igmax = igmin + BasisSet.shells[i].n_prims - 1;

    for(j=0;j<BasisSet.num_shells;j++) {  /* J-shell loop */
      jatom = BasisSet.shells[j].center - 1;
      if (jatom == iatom) continue;
      bx = Molecule.centers[jatom].x;
      by = Molecule.centers[jatom].y;
      bz = Molecule.centers[jatom].z;
      jang  = BasisSet.shells[j].am - 1;
      nbfj = (jang+2)*(jang+1)/2;
      jgmin = BasisSet.shells[j].fprim - 1;
      j_1stbf = BasisSet.shells[j].fao - 1;
      jgmax = jgmin + BasisSet.shells[j].n_prims - 1;

      eq_shell = (i == j);
      eq_atoms = (iatom == jatom);

      Rab2 = (bx-ax)*(bx-ax) + (by-ay)*(by-ay) + (bz-az)*(bz-az);

      for(ig=igmin;ig<=igmax;ig++) { /* loop over primitives ig in shell I */
        ai = BasisSet.cgtos[ig].exp;
        contr_i = BasisSet.cgtos[ig].ccoeff[iang];

        for(jg=jgmin;jg<=jgmax;jg++) { /* loop over primitives jg in shell J */
          bj = BasisSet.cgtos[jg].exp;
          contr_j = BasisSet.cgtos[jg].ccoeff[jang];
          gamma = ai + bj;

          tmp = ai * bj * Rab2 / gamma;
          over_pf = contr_i * contr_j * exp(-tmp) * pow(M_PI/gamma,1.5);
          //if (eq_shell && ig != jg)
          //  over_pf *= 2.0;

          /* Computing moment integrals */
          px = (ax*ai + bx*bj)/gamma;
          py = (ay*ai + by*bj)/gamma;
          pz = (az*ai + bz*bj)/gamma;
          if (eq_atoms)
            pax = pay = paz = pbx = pby = pbz = 0.0;
          else {
            pax = px - ax;  pay = py - ay;  paz = pz - az;
            pbx = px - bx;  pby = py - by;  pbz = pz - bz;
          }
          MI_OSrecurs(pax,pay,paz,pbx,pby,pbz,gamma,iang+2,jang+2,0);

          for(ibf=0;ibf<nbfi;ibf++) { /* loop over gaussians ibf */
            lx1 = xpow_bf[iang][ibf];
            ly1 = ypow_bf[iang][ibf];
            lz1 = zpow_bf[iang][ibf];

            for(jbf=0;jbf<nbfj;jbf++) { /* loop over gaussians jbf */
              lx2 = xpow_bf[jang][jbf];
              ly2 = ypow_bf[jang][jbf];
              lz2 = zpow_bf[jang][jbf];

              x0 = MIX[lx1][lx2][0];
              y0 = MIY[ly1][ly2][0];
              z0 = MIZ[lz1][lz2][0];

              /* One-dimensional del-squared integrals */
              tx = bj*(2*lx2+1)*x0 - 2*bj*bj*MIX[lx1][lx2+2][0];
              if (lx2 >= 2)
                tx -= 0.5*lx2*(lx2-1)*MIX[lx1][lx2-2][0];
              ty = bj*(2*ly2+1)*y0 - 2*bj*bj*MIY[ly1][ly2+2][0];
              if (ly2 >= 2)
                ty -= 0.5*ly2*(ly2-1)*MIY[ly1][ly2-2][0];
              tz = bj*(2*lz2+1)*z0 - 2*bj*bj*MIZ[lz1][lz2+2][0];
              if (lz2 >= 2)
                tz -= 0.5*lz2*(lz2-1)*MIZ[lz1][lz2-2][0];

              t2x = 4*(2*lx1+1)*(2*lx2+1)*ai*bj*x0 -
                8*(2*lx1+1)*ai*bj*bj*MIX[lx1][lx2+2][0] -
                8*(2*lx2+1)*ai*ai*bj*MIX[lx1+2][lx2][0] +
                16*ai*ai*bj*bj*MIX[lx1+2][lx2+2][0];
              if (lx1 >= 2)
                t2x += 4*lx1*(lx1-1)*bj*bj*MIX[lx1-2][lx2+2][0] -
                  2*lx1*(lx1-1)*(2*lx2+1)*bj*MIX[lx1-2][lx2][0];
              if (lx2 >= 2)
                t2x += 4*lx2*(lx2-1)*ai*ai*MIX[lx1+2][lx2-2][0] -
                  2*(2*lx1+1)*lx2*(lx2-1)*ai*MIX[lx1][lx2-2][0];
              if ((lx1 >= 2) && (lx2 >= 2))
                t2x += lx1*(lx1-1)*lx2*(lx2-1)*MIX[lx1-2][lx2-2][0];
 
              t2y = 4*(2*ly1+1)*(2*ly2+1)*ai*bj*y0 -
                8*(2*ly1+1)*ai*bj*bj*MIY[ly1][ly2+2][0] -
                8*(2*ly2+1)*ai*ai*bj*MIY[ly1+2][ly2][0] +
                16*ai*ai*bj*bj*MIY[ly1+2][ly2+2][0];
              if (ly1 >= 2)
                t2y += 4*ly1*(ly1-1)*bj*bj*MIY[ly1-2][ly2+2][0] -
                  2*ly1*(ly1-1)*(2*ly2+1)*bj*MIY[ly1-2][ly2][0];
              if (ly2 >= 2)
                t2y += 4*ly2*(ly2-1)*ai*ai*MIY[ly1+2][ly2-2][0] -
                  2*(2*ly1+1)*ly2*(ly2-1)*ai*MIY[ly1][ly2-2][0];
              if ((ly1 >= 2) && (ly2 >= 2))
                t2y += ly1*(ly1-1)*ly2*(ly2-1)*MIY[ly1-2][ly2-2][0];
 
              t2z = 4*(2*lz1+1)*(2*lz2+1)*ai*bj*z0 -
                8*(2*lz1+1)*ai*bj*bj*MIZ[lz1][lz2+2][0] -
                8*(2*lz2+1)*ai*ai*bj*MIZ[lz1+2][lz2][0] +
                16*ai*ai*bj*bj*MIZ[lz1+2][lz2+2][0];
              if (lz1 >= 2)
                t2z += 4*lz1*(lz1-1)*bj*bj*MIZ[lz1-2][lz2+2][0] -
                  2*lz1*(lz1-1)*(2*lz2+1)*bj*MIZ[lz1-2][lz2][0];
              if (lz2 >= 2)
                t2z += 4*lz2*(lz2-1)*ai*ai*MIZ[lz1+2][lz2-2][0] -
                  2*(2*lz1+1)*lz2*(lz2-1)*ai*MIZ[lz1][lz2-2][0];
              if ((lz1 >= 2) && (lz2 >= 2))
                t2z += lz1*(lz1-1)*lz2*(lz2-1)*MIZ[lz1-2][lz2-2][0];

                mmvc_ints[3*atom_C+xyz][i_1stbf+ibf][j_1stbf+jbf] += norm_bf[iang][ibf] * norm_bf[jang][jbf] *
                  over_pf * (t2x*y0*z0 + x0*t2y*z0 + x0*y0*t2z + 8*tx*ty*z0 + 8*tx*y0*tz + 8*x0*ty*tz);
            } /* end loop gaussians jbf */
          } /* end loop gaussians ibf */
        } /* end loop jg primitives */
      } /* end loop ig primitives */
    } /* end loop over j shell */
  } /* end loop over i shell */
  }
  }

  for (i=0;i<3*natom;++i)
    for (j=0;j<nao;++j)
      for (k=0;k<nao;++k)
        dmvc_ints[i][j][k] = (pmvc_ints[i][j][k]-mmvc_ints[i][j][k]) /
         (2.0*TEST_OE_DERIV1_MVC_DISP);

  if (UserOptions.print_lvl >= PRINT_DEBUG) {
    for (i=0; i<3*natom; ++i) {
      fprintf(outfile,"\n\tMass-velocity integral derivatives, numerical (%d)\n",i);
      print_mat(dmvc_ints[i],nao,nao,outfile);
    }
    for (i=0; i<3*natom; ++i) {
      tval = 0.0;
      for (j=0; j<nao; ++j)
        for (k=0; k<nao; ++k)
          tval += dmvc_ints[i][j][k];
      fprintf(outfile,"Sum of lower-triangle MVC ints for derivative %d (numerical): %15.10lf\n", i, tval);
    }
  }

  free_box(MIX, max_am+3,max_am+3);
  free_box(MIY, max_am+3,max_am+3);
  free_box(MIZ, max_am+3,max_am+3);
  free_3d_array(pmvc_ints,3*natom,nao);
  free_3d_array(mmvc_ints,3*natom,nao);
  free_3d_array(dmvc_ints,3*natom,nao);

  free_block(norm_bf);
  free_int_matrix(xpow_bf);
  free_int_matrix(ypow_bf);
  free_int_matrix(zpow_bf);
  return;
}

}}

