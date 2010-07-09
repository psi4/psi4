/*! \file
    \ingroup OPTKING
    \brief generates an empirical guess Hessian from a given set of
      salcs according to Schlegel, Theor. Chim. Acta, 66, 333 (1984) or
      according to Fischer and Almlof, J. Phys. Chem., 96, 9770 (1992).
*/

#define EXTERN
#include "globals.h"
#undef EXTERN
#include "cartesians.h"
#include "simples.h"
#include "salc.h"
#include "opt.h"

#include <cov_radii.h>
#include <libpsio/psio.h>

namespace psi { //namespace optking {

// returns covalent bond length in bohr from atomic numbers
inline double Rcov(int ZA, int ZB) {
  return (cov_radii[ZA] + cov_radii[ZB]) / _bohr2angstroms;
}
inline double Rcov(double ZA, double ZB) {
  return Rcov((int) ZA, (int) ZB);
}

// returns period from atomic number
inline int period(int ZA) {
  if      (ZA ==  1) return 1;
  else if (ZA <= 10) return 2;
  else if (ZA <= 18) return 3;
  else if (ZA <= 36) return 4;
  else               return 5;
}
inline int period(double ZA) {
  return period((int) ZA);
}

void empirical_H(const simples_class & simples, const salc_set &symm, const cartesians &carts) {
  int i, j, k, atomA, atomB, atomC, atomD, simple, count = -1, perA, perB;
  int I, K, L;
  int a, b, natom;
  double rAB, rBC, rABcov, rBCcov, rBDcov, prefactor_i;
  double A, B, C, D, E, r1[3], r2[3], r3[3], tval;
  double *f, *Z, val, *coord, norm_r1, norm_r2, norm_r3;

  if (optinfo.empirical_H == OPTInfo::SCHLEGEL)
    fprintf(outfile,"\nGenerating empirical Hessian (Schlegel)\n");
  else if (optinfo.empirical_H == OPTInfo::FISCHER)
    fprintf(outfile,"\nGenerating empirical Hessian (Fischer & Almlof)\n");

  f = init_array(simples.get_num());      
  coord = carts.get_coord();
  natom = carts.get_natom();
  Z = carts.get_fatomic_num();

  // Form diagonal Hessian in simple internals first
  if (optinfo.empirical_H == OPTInfo::SCHLEGEL) {
    for (i=0;i<simples.get_num(STRE); ++i) {
      atomA = simples.get_atom(STRE, i, 0);
      atomB = simples.get_atom(STRE, i, 1);
      perA = period(carts.get_Z(atomA));
      perB = period(carts.get_Z(atomB));
      A = 1.734;

      if ((perA==1) && (perB==1))
         B = -0.244;
      else if (((perA==1) && (perB==2)) || ((perB==1) && (perA==2)))
         B = 0.352;
      else if ((perA==2) && (perB==2))
         B = 1.085;
      else if (((perA==1) && (perB==3)) || ((perB==1) && (perA==3)))
         B = 0.660;
      else if (((perA==2) && (perB==3)) || ((perB==2) && (perA==3)))
         B = 1.522;
      else
         B = 2.068;

      rAB = carts.R(atomA,atomB);
      // fc in au / bohr^2
      tval = A/((rAB-B)*(rAB-B)*(rAB-B));
      // fc in aJ/Ang^2
      f[++count] = tval * _hartree2aJ / SQR(_bohr2angstroms);
    }
    for (i=0;i<simples.get_num(BEND);++i) {
      atomA = simples.get_atom(BEND, i, 0);
      atomB = simples.get_atom(BEND, i, 1);
      atomC = simples.get_atom(BEND, i, 2);
      if ( ((int) (carts.get_Z(atomA) == 1)) ||
           ((int) (carts.get_Z(atomC) == 1)) )
        val = 0.160;
      else
        val = 0.250;
      f[++count] = val * _hartree2aJ;
    }
    for (i=0;i<simples.get_num(TORS); ++i) {
      atomB = simples.get_atom(TORS, i, 1);
      atomC = simples.get_atom(TORS, i, 2);
      A = 0.0023;
      B = 0.07;
      rBCcov = Rcov(carts.get_Z(atomB), carts.get_Z(atomC));
      rBC = carts.R(atomB,atomC);
      if (rBC > (rBCcov + A/B)) B = 0.0; // keep > 0
      f[++count] = (A - (B*(rBC - rBCcov))) * _hartree2aJ;
    }
    for (i=0;i<simples.get_num(OUT); ++i) {
      atomA = simples.get_atom(OUT, i, 0);
      atomB = simples.get_atom(OUT, i, 1);
      atomC = simples.get_atom(OUT, i, 2);
      atomD = simples.get_atom(OUT, i, 3);
      A = 0.045;
      for (j=0;j<3;++j) {
         r1[j] = coord[3*atomA+j] - coord[3*atomB+j];
         r2[j] = coord[3*atomC+j] - coord[3*atomB+j];
         r3[j] = coord[3*atomD+j] - coord[3*atomB+j];
      }
      norm_r1 = sqrt( SQR(r1[0]) + SQR(r1[1]) + SQR(r1[2]) );
      norm_r2 = sqrt( SQR(r2[0]) + SQR(r2[1]) + SQR(r2[2]) );
      norm_r3 = sqrt( SQR(r3[0]) + SQR(r3[1]) + SQR(r3[2]) );

      tval =  ( r1[0] * (r2[1]*r3[2] - r2[2]*r3[1])
               -r1[1] * (r2[0]*r3[2] - r2[2]*r3[0])
               +r1[2] * (r2[0]*r3[1] - r2[1]*r3[0]) ) / (norm_r1 * norm_r2 * norm_r3);
      // I only want magnitude of this angle for force constant evaluation so make sure not negative
      tval = 1 - fabs(tval);
      f[++count] = A * pow(tval,4) * _hartree2aJ;
    }
    for (i=0;i<simples.get_num(LINB); ++i) {
      f[++count] = 0.10 * _hartree2aJ;
    }
  }
  else if (optinfo.empirical_H == OPTInfo::FISCHER) {
    for (i=0;i<simples.get_num(STRE); ++i) {
      atomA = simples.get_atom(STRE, i, 0);
      atomB = simples.get_atom(STRE, i, 1);

      rAB   = carts.R(atomA,atomB);
      rABcov = Rcov(carts.get_Z(atomA), carts.get_Z(atomB));

      A = 0.3601; B = 1.944;

      tval = A * exp(-B*(rAB - rABcov));
      f[++count] = tval * _hartree2aJ / SQR(_bohr2angstroms);
    }
    for (i=0;i<simples.get_num(BEND); ++i) {
      atomA = simples.get_atom(BEND, i, 0);
      atomB = simples.get_atom(BEND, i, 1);
      atomC = simples.get_atom(BEND, i, 2);

      rABcov = Rcov(carts.get_Z(atomA), carts.get_Z(atomB));
      rBCcov = Rcov(carts.get_Z(atomB), carts.get_Z(atomC));
      rAB = carts.R(atomA, atomB);
      rBC = carts.R(atomB, atomC);

      A = 0.089; B = 0.11; C = 0.44; D = -0.42;

      tval = A + B/(pow(rABcov*rBCcov, D)) * exp(-C*( rAB + rBC - rABcov - rBCcov));
      f[++count] = tval * _hartree2aJ;
    }
    for (i=0;i<simples.get_num(TORS); ++i) {
      atomB = simples.get_atom(TORS, i, 1);
      atomC = simples.get_atom(TORS, i, 2);

      rBCcov = Rcov(carts.get_Z(atomB), carts.get_Z(atomC));
      rBC = carts.R(atomB,atomC);

      A = 0.0015; B = 14.0; C = 2.85; D = 0.57; E = 4.00;

      // count number of additional bonds on central atoms
      double **bonds = simples.bond_connectivity_matrix(natom);
      //print_mat(bonds,natom,natom,outfile);
      L = 0;
      for (j=0; j<natom; ++j) {
        if (j == atomC) continue;
        if (bonds[atomB][j]) ++L;
      }
      for (j=0; j<natom; ++j) {
        if (j == atomB) continue;
        if (bonds[atomC][j]) ++L;
      }
      //fprintf(outfile,"Number of bonds %d\n", L);
      tval = B * pow(L,D) / pow(rBC * rBCcov, E);
      tval = A + tval * exp(-C * (rBC - rBCcov));
      f[++count] = tval * _hartree2aJ;
    }
    for (i=0;i<simples.get_num(OUT);++i) {
      atomA = simples.get_atom(OUT, i, 0);
      atomB = simples.get_atom(OUT, i, 1);
      atomC = simples.get_atom(OUT, i, 2);
      atomD = simples.get_atom(OUT, i, 3);
      val = simples.get_val_A_or_rad(OUT, i); // in rad

      A = 0.0025; B = 0.0061; C = 3.00; D = 4.00; E = 0.80;

      rABcov = Rcov(carts.get_Z(atomA), carts.get_Z(atomB));
      rBCcov = Rcov(carts.get_Z(atomB), carts.get_Z(atomC));
      rBDcov = Rcov(carts.get_Z(atomB), carts.get_Z(atomD));
      carts.R(atomA, atomB);

      tval = B * pow(rBCcov*rBDcov, E) * pow(cos(val),D);
      tval = A + tval * exp(-C*(rAB - rABcov));
      f[++count] = tval * _hartree2aJ;
    }
    for (i=0;i<simples.get_num(LINB);++i) {
      f[++count] = 0.10 * _hartree2aJ;
    }
  }

  // use same guess force constants for interfragment coordinates
  // regardless of chosen empirical_H type 
  for (i=0;i<simples.get_num(FRAG, 1); ++i) { // loop over fragment sets
    if (simples.frag_get_coord_on(i, 0)) {
      int min_a,min_b;
      rAB = 1e6;
      // find minimum distance between two atoms, one in each fragment
      for (a=0; a<simples.get_natom(FRAG, i, FRAG_A); ++a) {
        atomA = simples.get_atom(FRAG, i, a, FRAG_A);
        for (b=0; b<simples.get_natom(FRAG, i, FRAG_B); ++b) {
          atomB = simples.get_atom(FRAG, i, b, FRAG_B);
          tval = carts.R(atomA,atomB); // distance in au
          if (tval < rAB) {
            rAB = tval;
            min_a = atomA;
            min_b = atomB;
          }
        }
      }
      rABcov = Rcov(carts.get_Z(min_a), carts.get_Z(min_b));
      A = 0.3601; B = 1.944;
      tval = A * exp(-B*(rAB - rABcov));
      if (optinfo.frag_dist_rho == 1)
        tval *= pow(rAB,4);
      f[++count] = tval * _hartree2aJ / SQR(_bohr2angstroms);
    }
/* tried to guess fc from nuclear repulsion energy - but couldn't find a reasonable method
    for (I=0; I<6; ++I) {
      if (simples.frag.get_coord_on(i,I)) {
        int xyz, ivect;
        double weight, energy, *disp_coord_plus, *disp_coord_minus, delta = 0.01, first, second;
        energy = nuclear_repulsion(Z,coord);
        disp_coord_plus = new double [3*natom];
        disp_coord_minus = new double [3*natom];
        for (j=0; j<3*natom; ++j)
          disp_coord_plus[j] = disp_coord_minus[j] = coord[j];
        for (a=0; a<simples.frag.get_A_natom(i); ++a) { 
          atomA = simples.frag.get_A_atom(i,a);        
          for (K=0;K<simples.frag.get_A_P(i);++K) {    
            weight = simples.frag.get_A_weight(i,K,a); 
            for (xyz=0;xyz<3;++xyz) {
              tval = weight * simples.frag.get_A_s(i,I,3*K+xyz);
              disp_coord_plus[3*atomA+xyz]  += delta * tval;
              disp_coord_minus[3*atomA+xyz] -= delta * tval;
            }
          }
        }
        for (b=0; b<simples.frag.get_B_natom(i); ++b) {
          atomB = simples.frag.get_B_atom(i,b);       
          for (K=0;K<simples.frag.get_B_P(i);++K) {  
            weight = simples.frag.get_B_weight(i,K,b
            for (xyz=0;xyz<3;++xyz) {
              tval = weight * simples.frag.get_B_s(i,I,3*K+xyz);
              disp_coord_plus[3*atomB+xyz]  += delta * tval;
              disp_coord_minus[3*atomB+xyz] -= delta * tval;
            }
          }
        }
        first = (nuclear_repulsion(Z,disp_coord_plus) - nuclear_repulsion(Z,disp_coord_minus)) / (2*delta);
        first *= _hartree2aJ;
        if (I == 0) first = first / _bohr2angstroms;
        fprintf(outfile,"Nuclear first derivative (aJ/A,rad): %15.10lf\n", first);
        second = nuclear_repulsion(Z,disp_coord_plus) + nuclear_repulsion(Z,disp_coord_minus) - 2*energy;
        second = second / (delta*delta) * _hartree2aJ;
        if (I == 0) second = second / SQR(_bohr2angstroms);
        fprintf(outfile,"Nuclear second derivative (aJ/A^2,rad^2): %15.10lf\n", second);
        tval = delta * first + 0.5 * SQR(delta) * second;
        fprintf(outfile,"Energy change %15.10lf\n", tval);
        delete [] disp_coord_plus;
        delete [] disp_coord_minus;
      f[++count] = 0.001;
      }
    } */
    if (simples.frag_get_coord_on(i, 1))
      f[++count] = 0.001;
    if (simples.frag_get_coord_on(i, 2))
      f[++count] = 0.001;
    if (simples.frag_get_coord_on(i, 3))
      f[++count] = 0.0005;
    if (simples.frag_get_coord_on(i, 4))
      f[++count] = 0.0005;
    if (simples.frag_get_coord_on(i, 5))
      f[++count] = 0.0005;
  }
  free_array(coord);

  //fprintf(outfile,"Diagonal force constants for simple internals\n");
  //print_mat(&f,1,simples.get_num(),outfile);

 // Now transform into salc coordinates U^t H U
  double **intcos;
  double **f_new;

  intcos = block_matrix(symm.get_num(),simples.get_num());
  int id, index;

  for (i=0;i<symm.get_num();++i) {
    prefactor_i = symm.get_prefactor(i);
    for (j=0;j<symm.get_length(i);++j) {
      id = symm.get_simple(i,j);
      index = simples.id_to_index(id);
      intcos[i][index] = prefactor_i * symm.get_coeff(i,j);
    }
  }

  // fprintf(outfile,"Simples to Salc matrix\n");
  // print_mat(intcos,symm.get_num(),simples.get_num(),outfile);

  f_new = block_matrix(symm.get_num(),symm.get_num());
  for (i=0;i<symm.get_num();++i)
    for (j=0;j<symm.get_num();++j)
      for (k=0;k<simples.get_num();++k)
        f_new[i][j] += intcos[i][k] * f[k] * intcos[j][k]; 

 /*** write to PSIF_OPTKING ***/
  open_PSIF();
  psio_write_entry(PSIF_OPTKING, "Symmetric Force Constants",
      (char *) &(f_new[0][0]),symm.get_num()*symm.get_num()*sizeof(double));
  close_PSIF();

  free_array(f);
  free_block(f_new);
  free_block(intcos);
  return;
}

}//}

