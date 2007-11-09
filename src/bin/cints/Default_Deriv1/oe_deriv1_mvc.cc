/*! \file oe_deriv1_mvc.c;
  \ingroup (CINTS)
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
#include <libchkpt/chkpt.h>
#include <psifiles.h>
#include <physconst.h>

#include "defines.h"
#define EXTERN
#include "global.h"

namespace psi {
  namespace CINTS {

double ***MIX, ***MIY, ***MIZ;
void MI_OSrecurs(double pax, double pay, double paz,
                 double pbx, double pby, double pbz, double gamma,
                 int lmaxi, int lmaxj, int maxm);
extern double ***init_box(int a, int b, int c);
extern void ***free_box(double ***box, int a, int b);
extern void oe_deriv1_mvc_test(void);

void oe_deriv1_mvc(void) {
  double **mvc_ints, contr, ax, ay, az, bx, by, bz, ai, bj,  energy, **norm_bf;
  double  Rab2, ***UX, tmp, over_pf, **mvc_deriv_orb, ***dmvc_ints, **mvc_deriv;
  double x0, y0, z0, tx, ty, tz, t2x, t2y, t2z, gamma, **der_rho, **temp_mat, tx1, ty1, tz1;
  double pax, pay, paz, pbx, pby, pbz, px, py, pz, contr_i, contr_j, tval;

  int natom, nao, max_am, **xpow_bf, **ypow_bf, **zpow_bf, xpow, ypow, zpow ;
  int i, j, k, l, m, coord, occ_offset, eq_shell, eq_atoms, datom;
  int lx1, ly1, lz1, lx2, ly2, lz2;
  int ibf, iatom, iang, nbfi, ig, igmin, igmax, i_1stbf;
  int jbf, jatom, jang, nbfj, jg, jgmin, jgmax, j_1stbf, jlim;
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
        norm_bf[l][ibf] = sqrt(df[2*l]/(df[2*(l-i)]*df[2*(i-j)]*df[2*j]));
        ibf++;
      }
    }
  }

  MIX = init_box(max_am+3,max_am+3,0);
  MIY = init_box(max_am+3,max_am+3,0);
  MIZ = init_box(max_am+3,max_am+3,0);
  MIX[0][0][0] = MIY[0][0][0] = MIZ[0][0][0] = 1.0;

  /* compute mass-velocity integrals */
  mvc_ints = block_matrix(nao,nao);
  for(i=0;i<BasisSet.num_shells;i++) {  /* I-shell loop */ 
    iatom = BasisSet.shells[i].center - 1;
    ax = Molecule.centers[iatom].x;
    ay = Molecule.centers[iatom].y;
    az = Molecule.centers[iatom].z;
    iang  = BasisSet.shells[i].am - 1;
    nbfi = (iang+2)*(iang+1)/2;
    igmin = BasisSet.shells[i].fprim - 1;
    i_1stbf = BasisSet.shells[i].fao - 1;
    igmax = igmin + BasisSet.shells[i].n_prims - 1;

    for(j=0;j<=i;j++) {  /* J-shell loop */
      jatom = BasisSet.shells[j].center - 1;
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
        if (eq_shell)
          jgmax = ig;

        for(jg=jgmin;jg<=jgmax;jg++) { /* loop over primitives jg in shell J */
          bj = BasisSet.cgtos[jg].exp;
          contr_j = BasisSet.cgtos[jg].ccoeff[jang];
          gamma = ai + bj;

          tmp = ai * bj * Rab2 / gamma;
          over_pf = contr_i * contr_j * exp(-tmp) * pow(M_PI/gamma,1.5);
          if (eq_shell && ig != jg)
            over_pf *= 2.0;

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
            if (eq_shell) jlim = ibf + 1;
            else jlim = nbfj;

            for(jbf=0;jbf<jlim;jbf++) { /* loop over gaussians jbf */
              lx2 = xpow_bf[jang][jbf];
              ly2 = ypow_bf[jang][jbf];
              lz2 = zpow_bf[jang][jbf];

              x0 = MIX[lx1][lx2][0];
              y0 = MIY[ly1][ly2][0];
              z0 = MIZ[lz1][lz2][0];

              /* One-dimensional KE integrals */
              tx = bj*(2*lx2+1)*x0 - 2*bj*bj*MIX[lx1][lx2+2][0];
              if (lx2 >= 2)
                tx -= 0.5*lx2*(lx2-1)*MIX[lx1][lx2-2][0];
              ty = bj*(2*ly2+1)*y0 - 2*bj*bj*MIY[ly1][ly2+2][0];
              if (ly2 >= 2)
                ty -= 0.5*ly2*(ly2-1)*MIY[ly1][ly2-2][0];
              tz = bj*(2*lz2+1)*z0 - 2*bj*bj*MIZ[lz1][lz2+2][0];
              if (lz2 >= 2)
                tz -= 0.5*lz2*(lz2-1)*MIZ[lz1][lz2-2][0];

              /* Diagonal two-dimensional integrals */
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

                mvc_ints[i_1stbf+ibf][j_1stbf+jbf] += norm_bf[iang][ibf] * norm_bf[jang][jbf] *
                  over_pf * (t2x*y0*z0 + x0*t2y*z0 + x0*y0*t2z + 8*tx*ty*z0 + 8*tx*y0*tz + 8*x0*ty*tz);
            } /* end loop gaussians jbf */
          } /* end loop gaussians ibf */
        } /* end loop jg primitives */
      } /* end loop ig primitives */
    } /* end loop over j shell */
  } /* end loop over i shell */

  if (UserOptions.print_lvl >= PRINT_DEBUG) {
    fprintf(outfile,"Values of mass-velocity integrals in AO basis\n");
    print_mat(mvc_ints,nao,nao,outfile);
  }

  /* compute mass-velocity energy correction from lower triangle of integrals */
  energy = 0.0;
  for (i=0; i<nao; ++i) {
    for (j=0; j<=i; ++j)  {
      if (i == j)
        energy += Dens[i][j] * mvc_ints[i][j];
      else
        energy += 2.0 * Dens[i][j] * mvc_ints[i][j];
    }
  }
  energy *= -1/(8.0 * _c_au * _c_au) *
    UserOptions.fine_structure_alpha * UserOptions.fine_structure_alpha;
  fprintf(outfile,"\n\tMass-Velocity Correction energy: %15.10lf\n\n", energy);
  energy = energy + chkpt_rd_etot();
  chkpt_wt_etot(energy);

  /* compute orbital derivative part of gradient */
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
  mvc_deriv_orb = block_matrix(natom,3);
  for (iatom=0; iatom<natom; ++iatom) {
    for (xyz=0; xyz<3; ++xyz) {
      mmult(UX[3*iatom+xyz], 1, MOInfo.scf_evec[0], 0,  temp_mat, 0, nmo, nmo, nao, 0);
      zero_mat(der_rho,nao,nao);
      for (i=0; i<nao; ++i) {
        for (j=0; j<nao; ++j) {
          occ_offset = 0;
          cnt = -1;
          for(irrep=0; irrep<Symmetry.nirreps; irrep++) {
            for(m=0; m<MOInfo.clsdpi[irrep]; m++) /* 2.0 is for docc orbitals */
              der_rho[i][j] += 2.0 * MOInfo.scf_evec[0][m+occ_offset][i]*temp_mat[m+occ_offset][j];
            occ_offset += MOInfo.orbspi[irrep];
          }
        }
      }
      for (i=0; i<nao; ++i) { /* use lower triangle mvc_ints */
        for (j=0; j<=i; ++j) {
          if (i == j)
            mvc_deriv_orb[iatom][xyz] += 2.0 * der_rho[i][j] * mvc_ints[i][j];
          else
            mvc_deriv_orb[iatom][xyz] += 2.0 * (der_rho[i][j] + der_rho[j][i]) * mvc_ints[i][j];
        }
      }
    }
  }

  /* dE/dRx = -1/8c^2 * [d(Cu)Cv + Cu*d(Cv)] = -1/(8c^2)(2*dCu*Cv) */ 
  for (iatom=0; iatom<natom; ++iatom)
    for (xyz=0; xyz<3; ++xyz)
      mvc_deriv_orb[iatom][xyz] *= -1.0/(8.0*_c_au*_c_au) * 
        UserOptions.fine_structure_alpha * UserOptions.fine_structure_alpha;

  free_block(temp_mat);
  free_block(der_rho);

  /* compute derivative mass-velocity integrals */
  dmvc_ints = init_box(3*natom,nao,nao);
  for (datom=0; datom<natom; ++datom) {
    for (xyz=0; xyz<3; ++xyz) {
      for(i=0;i<BasisSet.num_shells;i++) {  /* I-shell loop */ 
        iatom = BasisSet.shells[i].center - 1;
        if (iatom != datom) continue;
        ax = Molecule.centers[iatom].x;
        ay = Molecule.centers[iatom].y;
        az = Molecule.centers[iatom].z;
        iang  = BasisSet.shells[i].am - 1;
        nbfi = (iang+2)*(iang+1)/2;
        igmin = BasisSet.shells[i].fprim - 1;
        i_1stbf = BasisSet.shells[i].fao - 1;
        igmax = igmin + BasisSet.shells[i].n_prims - 1;
    
        for(j=0;j<BasisSet.num_shells;j++) {  /* J-shell loop */
          jatom = BasisSet.shells[j].center - 1;
          if (iatom == jatom) continue;
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
              if (eq_shell && ig != jg)
                over_pf *= 2.0;
    
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
              MI_OSrecurs(pax,pay,paz,pbx,pby,pbz,gamma,iang+3,jang+2,0);
    
              /* First term: (x-Rx)^l d/dRx[e^(-ai*R^2)] = 2*ai*(x-Rx)^(lx+1)*[e^(-ai*R^2)] */
              for(ibf=0;ibf<nbfi;ibf++) { /* loop over gaussians ibf */
                lx1 = xpow_bf[iang][ibf];
                ly1 = ypow_bf[iang][ibf];
                lz1 = zpow_bf[iang][ibf];
                if (xyz == 0)      lx1 += 1;
                else if (xyz == 1) ly1 += 1;
                else if (xyz == 2) lz1 += 1;

                for(jbf=0;jbf<nbfj;jbf++) { /* loop over gaussians jbf */
                  lx2 = xpow_bf[jang][jbf];
                  ly2 = ypow_bf[jang][jbf];
                  lz2 = zpow_bf[jang][jbf];
    
                  x0 = MIX[lx1][lx2][0];
                  y0 = MIY[ly1][ly2][0];
                  z0 = MIZ[lz1][lz2][0];

                  /* One-dimensional KE integrals */
                  tx = bj*(2*lx2+1)*x0 - 2*bj*bj*MIX[lx1][lx2+2][0];
                  if (lx2 >= 2)
                    tx -= 0.5*lx2*(lx2-1)*MIX[lx1][lx2-2][0];
                  ty = bj*(2*ly2+1)*y0 - 2*bj*bj*MIY[ly1][ly2+2][0];
                  if (ly2 >= 2)
                    ty -= 0.5*ly2*(ly2-1)*MIY[ly1][ly2-2][0];
                  tz = bj*(2*lz2+1)*z0 - 2*bj*bj*MIZ[lz1][lz2+2][0];
                  if (lz2 >= 2)
                    tz -= 0.5*lz2*(lz2-1)*MIZ[lz1][lz2-2][0];

                  /* Diagonal two-dimensional integrals */
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

                  tval = norm_bf[iang][ibf] * norm_bf[jang][jbf] * over_pf *
                    (t2x*y0*z0 + x0*t2y*z0 + x0*y0*t2z + 8*tx*ty*z0 + 8*tx*y0*tz + 8*x0*ty*tz);

                  tval *= 2.0 * ai;

                  dmvc_ints[3*datom+xyz][i_1stbf+ibf][j_1stbf+jbf] += tval;

                } /* end loop gaussians jbf */
              } /* end loop gaussians ibf */
              /* Second term: d/dRx [(x-Rx)^lx] e^(-ai*R^2) = lx*(x-Rx)^(lx-1)*[e^(-ai*R^2)] */
              for(ibf=0;ibf<nbfi;ibf++) { /* loop over gaussians ibf */
                lx1 = xpow_bf[iang][ibf];
                ly1 = ypow_bf[iang][ibf];
                lz1 = zpow_bf[iang][ibf];

                if ( ((xyz==0) && (lx1==0)) ||
                     ((xyz==1) && (ly1==0)) ||
                     ((xyz==2) && (lz1==0)) ) continue;

                if (xyz == 0)      lx1 -= 1;
                else if (xyz == 1) ly1 -= 1;
                else if (xyz == 2) lz1 -= 1;

                for(jbf=0;jbf<nbfj;jbf++) { /* loop over gaussians jbf */
                  lx2 = xpow_bf[jang][jbf];
                  ly2 = ypow_bf[jang][jbf];
                  lz2 = zpow_bf[jang][jbf];
    
                  x0 = MIX[lx1][lx2][0];
                  y0 = MIY[ly1][ly2][0];
                  z0 = MIZ[lz1][lz2][0];

                  /* One-dimensional KE integrals */
                  tx = bj*(2*lx2+1)*x0 - 2*bj*bj*MIX[lx1][lx2+2][0];
                  if (lx2 >= 2)
                    tx -= 0.5*lx2*(lx2-1)*MIX[lx1][lx2-2][0];
                  ty = bj*(2*ly2+1)*y0 - 2*bj*bj*MIY[ly1][ly2+2][0];
                  if (ly2 >= 2)
                    ty -= 0.5*ly2*(ly2-1)*MIY[ly1][ly2-2][0];
                  tz = bj*(2*lz2+1)*z0 - 2*bj*bj*MIZ[lz1][lz2+2][0];
                  if (lz2 >= 2)
                    tz -= 0.5*lz2*(lz2-1)*MIZ[lz1][lz2-2][0];

                  /* Diagonal two-dimensional integrals */
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

                  tval = norm_bf[iang][ibf] * norm_bf[jang][jbf] * over_pf *
                    (t2x*y0*z0 + x0*t2y*z0 + x0*y0*t2z + 8*tx*ty*z0 + 8*tx*y0*tz + 8*x0*ty*tz);

                  if (xyz == 0) tval *= -1.0 * (lx1+1);
                  else if (xyz == 1) tval *= -1.0 * (ly1+1);
                  else if (xyz == 2) tval *= -1.0 * (lz1+1);

                  dmvc_ints[3*datom+xyz][i_1stbf+ibf][j_1stbf+jbf] += tval;
                } /* end loop gaussians jbf */
              } /* end loop gaussians ibf */
            } /* end loop jg primitives */
          } /* end loop ig primitives */
        } /* end loop over j shell */
      } /* end loop over i shell */
    } /* end xyz */
  } /* end iatom */

  if (UserOptions.print_lvl >= PRINT_DEBUG) {
    for (i=0; i<3*natom; ++i) {
      fprintf(outfile,"\n\tMass-velocity integral derivatives, analytic (%d)\n",i);
      print_mat(dmvc_ints[i],nao,nao,outfile);
    }
    for (i=0; i<3*natom; ++i) {
      tval = 0.0;
      for (j=0; j<nao; ++j)
        for (k=0; k<nao; ++k)
          tval += dmvc_ints[i][j][k];
      fprintf(outfile,"Sum of MVC integrals for derivative %d (analytic): %15.10lf\n", i, tval);
    }
  }

  /* compute integral derivative part of gradient  */
  /* -1/(8c^2) Duv [ del^2(d(phi_u)) del^2(phi_v) + del^2(phi_u) del^2(d(phi_v)) */
  /* -1/(4c^2) Duv [ del^2(d(phi_u)) del^2(phi_v) + del^2(phi_u) del^2(d(phi_v)) */
  mvc_deriv = block_matrix(natom,3);
  for (datom=0; datom<natom; ++datom) {
    for (xyz=0; xyz<3; ++xyz) {
      tval = 0.0;
      for (i=0; i<nao; ++i) {
        for (j=0; j<nao; ++j)
          tval += Dens[i][j] * dmvc_ints[3*datom+xyz][i][j];
      }
      mvc_deriv[datom][xyz] = -1.0 / (4.0 * _c_au * _c_au) * tval
        * UserOptions.fine_structure_alpha * UserOptions.fine_structure_alpha;
    }
  }

  if (UserOptions.print_lvl >= PRINT_DEBUG) {
    fprintf(outfile,"Mass-velocity gradient, integral derivative part\n");
    print_mat(mvc_deriv,natom,3,outfile);
    fprintf(outfile,"Mass-velocity gradient, orbital derivative part\n");
    print_mat(mvc_deriv_orb,natom,3,outfile);
  }

  for (iatom=0; iatom<natom; ++iatom)
    for (xyz=0; xyz<3; ++xyz)
      mvc_deriv[iatom][xyz] += mvc_deriv_orb[iatom][xyz];

  if (UserOptions.print_lvl >= PRINT_DEBUG) {
    fprintf(outfile,"\tTotal mass velocity gradient\n");
    print_mat(mvc_deriv,natom,3,outfile);
  }

  for (iatom=0; iatom<natom; ++iatom)
    for (xyz=0; xyz<3; ++xyz)
      Grad[iatom][xyz] += mvc_deriv[iatom][xyz];

  free_box(MIX, max_am+3,max_am+3);
  free_box(MIY, max_am+3,max_am+3);
  free_box(MIZ, max_am+3,max_am+3);

  free_block(mvc_deriv);
  free_block(mvc_deriv_orb);
  free_block(mvc_ints);
  free_box(dmvc_ints, 3*natom, nao);
  for (i=0;i<3*natom;++i)
    free_block(UX[i]);
  free(UX);

  free_block(norm_bf);
  free_int_matrix(xpow_bf);
  free_int_matrix(ypow_bf);
  free_int_matrix(zpow_bf);

  /* to test numerically use the following: */
  /*  oe_deriv1_mvc_test(); */

  return;
}

void MI_OSrecurs(double pax, double pay, double paz,
                 double pbx, double pby, double pbz, double gamma,
                 int lmaxi, int lmaxj, int maxm)
{
  int i,j,k;
  double pp = 1/(2*gamma);

    /* Computing starting integrals for recursive procedure */
                 
  if (maxm > 1)  
    MIX[0][0][2] = MIY[0][0][2] = MIZ[0][0][2] = pp;

    /* Upward recursion in j for i=0 */
  
  for(j=0;j<lmaxj;j++)
    for(k=0;k<=maxm;k++) {
      MIX[0][j+1][k] = pbx*MIX[0][j][k];
      MIY[0][j+1][k] = pby*MIY[0][j][k];
      MIZ[0][j+1][k] = pbz*MIZ[0][j][k];
      if (j>0) {
        MIX[0][j+1][k] += j*pp*MIX[0][j-1][k];
        MIY[0][j+1][k] += j*pp*MIY[0][j-1][k];
        MIZ[0][j+1][k] += j*pp*MIZ[0][j-1][k];
      }
      if (k>0) {
        MIX[0][j+1][k] += k*pp*MIX[0][j][k-1];
        MIY[0][j+1][k] += k*pp*MIY[0][j][k-1];
        MIZ[0][j+1][k] += k*pp*MIZ[0][j][k-1];
      }
    }

    /* Upward recursion in i for all j's */

  for(i=0;i<lmaxi;i++)
    for(j=0;j<=lmaxj;j++)
      for(k=0;k<=maxm;k++) {
        MIX[i+1][j][k] = pax*MIX[i][j][k];
        MIY[i+1][j][k] = pay*MIY[i][j][k];
        MIZ[i+1][j][k] = paz*MIZ[i][j][k];
        if (i>0) {
          MIX[i+1][j][k] += i*pp*MIX[i-1][j][k];
          MIY[i+1][j][k] += i*pp*MIY[i-1][j][k];
          MIZ[i+1][j][k] += i*pp*MIZ[i-1][j][k];
        }
        if (j>0) {
          MIX[i+1][j][k] += j*pp*MIX[i][j-1][k];
          MIY[i+1][j][k] += j*pp*MIY[i][j-1][k];
          MIZ[i+1][j][k] += j*pp*MIZ[i][j-1][k];
        }
        if (k>0) {
          MIX[i+1][j][k] += k*pp*MIX[i][j][k-1];
          MIY[i+1][j][k] += k*pp*MIY[i][j][k-1];
          MIZ[i+1][j][k] += k*pp*MIZ[i][j][k-1];
        }
      }
}

}}

