/*! \file
  \ingroup CINTS
  Computes derivatives of one-electron Darwin integrals numerically.
  \author R.A. King
*/
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>

#include <libciomr/libciomr.h>
#include <libint/libint.h>
#include <libderiv/libderiv.h>
#include <libqt/qt.h>
#include <physconst.h>
#include <psiconfig.h>

#include "defines.h"
#define EXTERN
#include "global.h"

#define TEST_OE_DERIV1_DARWIN1_DISP (0.0005)

namespace psi {
  namespace cints {

    static const int use_cca_integrals_standard = (PSI_INTEGRALS_STANDARD == 1);

void oe_deriv1_darwin1_test(void) {
  double **AO_at_nuc, contr, ax, ay, az, Rab, ai, e_darwin1, zval, **norm_bf;
  double Rab_x, Rab_y, Rab_z, sign, ***pAO_at_nuc, ***mAO_at_nuc, ***dAO_at_nuc, tval;
  int natom, nao, max_am, **xpow_bf, **ypow_bf, **zpow_bf, xpow, ypow, zpow;
  int i, j, k, l, ibf, iatom, iang, nbfi, ig, igmin, igmax, i_1stbf;
  int atom_A, atom_B, atom_C, xyz, ao;

  natom = Molecule.num_atoms;
  nao = BasisSet.num_ao;
  max_am = BasisSet.max_am;

  /* determine powers of x,y,z for each AO, and normalization constants */
  xpow_bf = init_int_matrix(max_am+1,(max_am+1)*(max_am+2)/2);
  ypow_bf = init_int_matrix(max_am+1,(max_am+1)*(max_am+2)/2);
  zpow_bf = init_int_matrix(max_am+1,(max_am+1)*(max_am+2)/2);
  norm_bf = init_matrix(max_am+1,(max_am+1)*(max_am+2)/2);
  for(l=0;l<=max_am;l++) {
    ibf = 0;
    for(i = 0; i <= l; i++){
      for(j = 0; j <= i; j++){
        xpow_bf[l][ibf] = l - i;
        ypow_bf[l][ibf] = i - j;
        zpow_bf[l][ibf] = j;
        norm_bf[l][ibf] = use_cca_integrals_standard ? 1.0 : sqrt(df[2*l]/(df[2*(l-i)]*df[2*(i-j)]*df[2*j]));
        ibf++;
      }
    }
  }

  /* compute AO's at atom positions */
  AO_at_nuc = block_matrix(natom,nao);
  for(i=0;i<BasisSet.num_shells;i++) {
    iatom = BasisSet.shells[i].center - 1;
    ax = Molecule.centers[iatom].x;      
    ay = Molecule.centers[iatom].y;     
    az = Molecule.centers[iatom].z;
    iang  = BasisSet.shells[i].am - 1;  
    nbfi = (iang+2)*(iang+1)/2;         
    igmin = BasisSet.shells[i].fprim - 1;
    i_1stbf = BasisSet.shells[i].fao - 1;
    igmax = igmin + BasisSet.shells[i].n_prims - 1;

    for(ibf=0;ibf<nbfi;ibf++) {         
      xpow = xpow_bf[iang][ibf];
      ypow = ypow_bf[iang][ibf];
      zpow = zpow_bf[iang][ibf];
      for(ig=igmin;ig<=igmax;ig++) {   
        ai    = BasisSet.cgtos[ig].exp;
        contr = BasisSet.cgtos[ig].ccoeff[iang];
        for(k=0;k<natom;k++) {
          Rab_x = Molecule.centers[k].x - ax;
          Rab_y = Molecule.centers[k].y - ay;
          Rab_z = Molecule.centers[k].z - az;
          Rab = Rab_x*Rab_x + Rab_y*Rab_y + Rab_z*Rab_z;

          tval = contr * norm_bf[iang][ibf] * exp(-ai * Rab);
          if (xpow) tval *= pow(Rab_x, xpow);
          if (ypow) tval *= pow(Rab_y, ypow);
          if (zpow) tval *= pow(Rab_z, zpow);
          AO_at_nuc[k][i_1stbf+ibf] += tval;
        }
      }
    }
  }
  //fprintf(outfile,"AO_at_nuc");
  //print_mat(AO_at_nuc,natom,nao,outfile);

  /* compute AO's at +R_Cx */
  pAO_at_nuc = init_3d_array(3*natom,natom,nao);
  for (atom_C=0; atom_C<natom; ++atom_C) {
    for (xyz=0; xyz<3; ++xyz) {

  for(i=0;i<BasisSet.num_shells;i++) {
    iatom = BasisSet.shells[i].center - 1;
    ax = Molecule.centers[iatom].x;      
    ay = Molecule.centers[iatom].y;     
    az = Molecule.centers[iatom].z;
    if (iatom==atom_C) {
      if (xyz==0) ax += TEST_OE_DERIV1_DARWIN1_DISP;
      if (xyz==1) ay += TEST_OE_DERIV1_DARWIN1_DISP;
      if (xyz==2) az += TEST_OE_DERIV1_DARWIN1_DISP;
    }
    iang  = BasisSet.shells[i].am - 1;  
    nbfi = (iang+2)*(iang+1)/2;         
    igmin = BasisSet.shells[i].fprim - 1;
    i_1stbf = BasisSet.shells[i].fao - 1;
    igmax = igmin + BasisSet.shells[i].n_prims - 1;

    for(ibf=0;ibf<nbfi;ibf++) {
      xpow = xpow_bf[iang][ibf];
      ypow = ypow_bf[iang][ibf];
      zpow = zpow_bf[iang][ibf];
      for(ig=igmin;ig<=igmax;ig++) {   
        ai    = BasisSet.cgtos[ig].exp;
        contr = BasisSet.cgtos[ig].ccoeff[iang];
        for(k=0;k<natom;k++) {        
          Rab_x = Molecule.centers[k].x - ax;
          Rab_y = Molecule.centers[k].y - ay;
          Rab_z = Molecule.centers[k].z - az;
          if (k==atom_C) {
            if (xyz==0) Rab_x += TEST_OE_DERIV1_DARWIN1_DISP;
            if (xyz==1) Rab_y += TEST_OE_DERIV1_DARWIN1_DISP;
            if (xyz==2) Rab_z += TEST_OE_DERIV1_DARWIN1_DISP;
          }
          Rab = Rab_x*Rab_x + Rab_y*Rab_y + Rab_z*Rab_z;

          tval = contr * norm_bf[iang][ibf] * exp(-ai * Rab); /* avoid 0^0 */
          if (xpow) tval *= pow(Rab_x, xpow);
          if (ypow) tval *= pow(Rab_y, ypow);
          if (zpow) tval *= pow(Rab_z, zpow);
          pAO_at_nuc[3*atom_C+xyz][k][i_1stbf+ibf] += tval;
        }
      }
    }
  }
  }
  }
//  fprintf(outfile,"pAO_at_nuc[2]");
//  print_mat(pAO_at_nuc[2],natom,nao,outfile);

  /* compute AO's at -R_Cx */
  mAO_at_nuc = init_3d_array(3*natom,natom,nao);
  for (atom_C=0; atom_C<natom; ++atom_C) {
    for (xyz=0; xyz<3; ++xyz) {

  for(i=0;i<BasisSet.num_shells;i++) {
    iatom = BasisSet.shells[i].center - 1;
    ax = Molecule.centers[iatom].x;      
    ay = Molecule.centers[iatom].y;     
    az = Molecule.centers[iatom].z;
    if (iatom==atom_C) {
      if (xyz==0) ax -= TEST_OE_DERIV1_DARWIN1_DISP;
      if (xyz==1) ay -= TEST_OE_DERIV1_DARWIN1_DISP;
      if (xyz==2) az -= TEST_OE_DERIV1_DARWIN1_DISP;
    }
    iang  = BasisSet.shells[i].am - 1;  
    nbfi = (iang+2)*(iang+1)/2;         
    igmin = BasisSet.shells[i].fprim - 1;
    i_1stbf = BasisSet.shells[i].fao - 1;
    igmax = igmin + BasisSet.shells[i].n_prims - 1;

    for(ibf=0;ibf<nbfi;ibf++) {         
      xpow = xpow_bf[iang][ibf];
      ypow = ypow_bf[iang][ibf];
      zpow = zpow_bf[iang][ibf];
      for(ig=igmin;ig<=igmax;ig++) {   
        ai    = BasisSet.cgtos[ig].exp;
        contr = BasisSet.cgtos[ig].ccoeff[iang];
        for(k=0;k<natom;k++) {        
          Rab_x = Molecule.centers[k].x - ax;
          Rab_y = Molecule.centers[k].y - ay;
          Rab_z = Molecule.centers[k].z - az;
          if (k==atom_C) {
            if (xyz==0) Rab_x -= TEST_OE_DERIV1_DARWIN1_DISP;
            if (xyz==1) Rab_y -= TEST_OE_DERIV1_DARWIN1_DISP;
            if (xyz==2) Rab_z -= TEST_OE_DERIV1_DARWIN1_DISP;
          }
          Rab = Rab_x*Rab_x + Rab_y*Rab_y + Rab_z*Rab_z;

          tval = contr * norm_bf[iang][ibf] * exp(-ai * Rab);
          if (xpow) tval *= pow(Rab_x, xpow);
          if (ypow) tval *= pow(Rab_y, ypow);
          if (zpow) tval *= pow(Rab_z, zpow);
          mAO_at_nuc[3*atom_C+xyz][k][i_1stbf+ibf] += tval;
        }
      }
    }
  }
  }
  }
//  fprintf(outfile,"mAO_at_nuc[2]");
//  print_mat(mAO_at_nuc[2],natom,nao,outfile);

  dAO_at_nuc = init_3d_array(3*natom,natom,nao);
  for (i=0;i<3*natom;++i)
    for (k=0;k<natom;++k)
      for (ao=0;ao<nao;++ao)
        dAO_at_nuc[i][k][ao] = (pAO_at_nuc[i][k][ao]-mAO_at_nuc[i][k][ao]) /
         (2.0*TEST_OE_DERIV1_DARWIN1_DISP);

//  for (i=0; i<3*natom; ++i) {
//    fprintf(outfile,"Numerical One-electron Darwin integral derivative %d\n",i);
//    print_mat(dAO_at_nuc[i],natom,nao,outfile);
//  }
  for (i=0; i<3*natom; ++i) {
    tval = 0.0;
    for (j=0; j<natom; ++j)
      for (k=0; k<nao; ++k)
        tval += dAO_at_nuc[i][j][k];
    fprintf(outfile,"Sum of dAO_at_nuc for derivative %d (numerical): %15.10lf\n", i, tval);
  }

  free_block(AO_at_nuc);
  free_3d_array(pAO_at_nuc,3*natom,natom);
  free_3d_array(mAO_at_nuc,3*natom,natom);
  free_3d_array(dAO_at_nuc,3*natom,natom);
  free_int_matrix(xpow_bf);
  free_int_matrix(ypow_bf);
  free_int_matrix(zpow_bf);
}

}}

