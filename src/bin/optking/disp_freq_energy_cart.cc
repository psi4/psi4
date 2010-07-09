/*! \file
    \ingroup OPTKING
    \brief 
 DISP_FREQ_ENERGY_CART makes symmetry-adapted cartesian coordinates and displaces
   along them to setup frequencies by energy points
Formulas for finite differences
3-point - diagonal
  O(1/h^2): [f(1,0) + f(-1,0) - 2f(0,0)]/(h^2)
   [which is same as : [f(2,0) + f(-2,0) - 2f(0,0)]/(4h^2)
3-point - off-diagonal
  O(1/h^2): old-way: [f(1,1) - f(1,-1) - f(-1,1) + f(-1,-1)]/(4h^2)
  O(1/h^2): new-way: [f(1,1)+f(-1,-1)+2f(0,0) -f(1,0) -f(-1,0) -f(0,1) -f(0,-1)]/(2h^2)

5-point formula - diagonal
  O(1/h^4): [-f(-2,0) + 16f(-1,0) + 16f(1,0) - f(2,0) - 30f(0,0)] / (12h^2)
5-point formula - off-diagonal
  O(1/h^4): old-way [ 1f(-2,-2) - 8f(-1,-2) + 8f(+1,-2) - 1f(+2,-2) - 8f(-2,-1)
  + 64f(-1,-1) - 64f(+1,-1) + 8f(+2,-1) + 8f(-2,+1) - 64f(-1,+1) + 64f(+1,+1)
 - 8f(+2,+1) - 1f(-2,+2) + 8f(-1,+2) - 8f(+1,+2) + 1f(+2,+2) / (144h^2)
  O(1/h^4): new-way [ -1f(-1,-2) - 1f(-2,-1) + 9f(-1,-1) - 1f(+1,-1)
    - 1f(-1,1) + 9f(+1,+1) - 1f(+2,+1) - 1f(1,2)
    + 1f(-2,0) - 7f(-1,0)  - 7f(+1,0) + 1f(+2,0)
    + 1f(0,-2) - 7f(0,-1)  - 7f(0,+1) + 1f(0,+2) + 12f(0,0)]/(12h^2)
*/

#define EXTERN
#include "globals.h"
#undef EXTERN
#include "cartesians.h"
#include "simples.h"
#include "salc.h"
#include "opt.h"

#include <libchkpt/chkpt.h>
#include <libqt/qt.h>
#include <libpsio/psio.h>

namespace psi { //namespace optking {

int disp_freq_energy_cart(const cartesians &carts)
{
  int i,j,a,b,I,k,dim, ndisp_all, nsalc_all, natom, atom, xyz, cnt, loner, *ndisp;
  int op, atom2, nirreps, natom_unique, irrep, diag_ind, atom_unique;
  double *coord, energy, ***disp_orig, **displacements, *character, *disp_e;
  double *v, *f, *q, *masses, **constraints, normval, dotval, **constraints_ortho, *Zvals;
  int **ict, print=0,nconstraints, *ua2a;
  int *nsalc_orig, *nsalc, *irrep_per_disp;
  double ***geoms, **moi, *v1_A, *v1_B, *big_zeroes;
  double *com, **X, *evals, tval, tval1, tval2, *disp_grad, **cartrep;
  double ***salc_orig, ***salc;
  psio_address next;
  // print=1;

  coord = carts.get_coord();
  energy = carts.get_energy();
  natom = carts.get_natom();
  masses = carts.get_mass();
  Zvals = carts.get_Zvals();
  nirreps = syminfo.nirreps;

  chkpt_init(PSIO_OPEN_OLD);
  natom_unique = chkpt_rd_num_unique_atom();
  ict = chkpt_rd_ict();
  ua2a = chkpt_rd_ua2a();
  cartrep = chkpt_rd_cartrep();
  chkpt_close();

  if (print) {
    fprintf(outfile,"Character table from syminfo\n");
    for (irrep=0; irrep<nirreps; ++irrep) {
      for (op=0; op<nirreps; ++op)
        fprintf(outfile,"%d",syminfo.ct[irrep][op]);
      fprintf(outfile,"\n");
    }
    fprintf(outfile,"Cartesian rep from chkpt\n");
    print_mat(cartrep,nirreps,9,outfile);
    fprintf(outfile,"ICT table\n");
    for (i=0; i<natom; ++i) {
      for (op=0; op<nirreps; ++op)
        fprintf(outfile," %d ", ict[op][i]);
      fprintf(outfile,"\n");
    }
  }

  /**** Make sure center of mass is at origin ****/
  tval = 0.0;
  com = init_array(3);
  for (i=0; i<natom; ++i) {
    for (xyz=0; xyz<3 ; ++xyz)
      com[xyz] += masses[3*i+xyz] * coord[3*i+xyz];
    tval += masses[3*i];
  }
  for (xyz=0; xyz<3 ; ++xyz)
    com[xyz] = com[xyz] / tval;

  for (i=0; i<natom; ++i)
    for (xyz=0; xyz<3 ; ++xyz)
      coord[3*i+xyz] -= com[xyz];

  if (print) {
    fprintf(outfile,"\n\nCenter of mass\n");
    for (xyz=0; xyz<3 ; ++xyz)
      fprintf(outfile,"%15.10lf", com[xyz]);
    fprintf(outfile,"\n\nNew Geometry wrt COM\n");
    for (i=0; i<natom; ++i)
      fprintf(outfile,"%15.10lf%15.10lf%15.10lf\n",coord[3*i],coord[3*i+1],coord[3*i+2]);
    fflush(outfile);
  }
  free_array(com);

  /**** Determine the rotational coordinates from MOI tensor ****/
  moi = block_matrix(3,3);
  for (i=0; i<natom; ++i) {
    moi[0][0] += masses[3*i] * (SQR(coord[3*i+1]) + SQR(coord[3*i+2]));
    moi[1][1] += masses[3*i] * (SQR(coord[3*i+0]) + SQR(coord[3*i+2]));
    moi[2][2] += masses[3*i] * (SQR(coord[3*i+0]) + SQR(coord[3*i+1]));
    moi[0][1] -= masses[3*i] * coord[3*i+0] * coord[3*i+1];
    moi[0][2] -= masses[3*i] * coord[3*i+0] * coord[3*i+2];
    moi[1][2] -= masses[3*i] * coord[3*i+1] * coord[3*i+2];
  }
  moi[1][0] = moi[0][1];
  moi[2][0] = moi[0][2];
  moi[2][1] = moi[1][2];

  X = block_matrix(3,3);
  evals = init_array(3);
  opt_sq_rsp(3,3,moi,evals,1,X,1.0E-14);
  free_block(moi);
  if (print) {
    fprintf(outfile,"Eigenvectors (X) of MOI tensor\n");
    print_mat(X,3,3,outfile);
  }

  /**** Build matrix with COM and rotational constraints as rows ****/
  constraints = block_matrix(6, 3*natom);
  /* COM constraints */
  for (i=0; i<natom; ++i) {
    constraints[0][3*i]   = 1.0 * sqrt(masses[3*i]);
    constraints[1][3*i+1] = 1.0 * sqrt(masses[3*i+1]);
    constraints[2][3*i+2] = 1.0 * sqrt(masses[3*i+2]);
  }
  /* Rotational constraints */
  for (i=0; i<natom; ++i) {
    tval1 = (coord[3*i+0] * X[0][1]) + (coord[3*i+1] * X[1][1]) + (coord[3*i+2] * X[2][1]);
    tval2 = (coord[3*i+0] * X[0][2]) + (coord[3*i+1] * X[1][2]) + (coord[3*i+2] * X[2][2]);
    for (xyz = 0; xyz < 3; ++xyz)
      constraints[3][3*i+xyz] = (tval1 * X[xyz][2] - tval2 * X[xyz][1]) * sqrt(masses[3*i]);
  }
  for (i=0; i<natom; ++i) {
    tval1 = (coord[3*i+0] * X[0][2]) + (coord[3*i+1] * X[1][2]) + (coord[3*i+2] * X[2][2]);
    tval2 = (coord[3*i+0] * X[0][0]) + (coord[3*i+1] * X[1][0]) + (coord[3*i+2] * X[2][0]);
    for (xyz = 0; xyz < 3; ++xyz)
      constraints[4][3*i+xyz] = (tval1 * X[xyz][0] - tval2 * X[xyz][2]) * sqrt(masses[3*i]);
  }
  for (i=0; i<natom; ++i) {
    tval1 = (coord[3*i+0] * X[0][0]) + (coord[3*i+1] * X[1][0]) + (coord[3*i+2] * X[2][0]);
    tval2 = (coord[3*i+0] * X[0][1]) + (coord[3*i+1] * X[1][1]) + (coord[3*i+2] * X[2][1]);
    for (xyz = 0; xyz < 3; ++xyz)
      constraints[5][3*i+xyz] = (tval1 * X[xyz][1] - tval2 * X[xyz][0]) * sqrt(masses[3*i]);
  }
  if (print) {
    fprintf(outfile,"Raw COM/Rotational Constraints\n");
    print_mat(constraints,6,3*natom,outfile);
  }
  free_block(X);

  /* Remove NULL constraint (if present) and normalize rest of them */
  cnt = -1;
  for (i = 0; i < 6; ++i) {
    dot_array(constraints[i], constraints[i], 3*natom, &normval);
    if (normval > 1.0E-10) {
      ++cnt;
      for (j=0; j<3*natom; ++j)
        constraints[cnt][j] = constraints[i][j] / sqrt(normval) ;
    }
  }
  nconstraints = cnt+1;

  /* Orthogonalize rotations and translations exactly-is this ever necessary?*/
  constraints_ortho = block_matrix(nconstraints,3*natom);
  cnt = 0;
  for (i=0; i<nconstraints; ++i)
    cnt += schmidt_add(constraints_ortho, cnt, 3*natom, constraints[i]);
  for (i=0; i<nconstraints; ++i)
    for (j=0; j<3*natom; ++j)
      constraints[i][j] = constraints_ortho[i][j];
  free_block(constraints_ortho);

  /**** Form symmetry-adapted cartesian vectors ****/
  salc_orig = (double ***) malloc(nirreps*sizeof(double **));
  nsalc_orig = init_int_array(nirreps);

  /* loop over irrep of coordinates */
  for (irrep=0; irrep < nirreps; ++irrep) {
    salc_orig[irrep] = block_matrix(3*natom,3*natom);
    cnt=-1;

    /* loop over unique atoms and xyz coordinates */
    for (atom_unique=0; atom_unique < natom_unique; ++atom_unique) {
      atom = ua2a[atom_unique];
      /* determine if atom is loner, having no symmetry-related partners */
      loner = 1;
      for (op=0; op < nirreps; ++op) { 
        if (atom != ict[op][atom] - 1)
         loner = 0;
      }

      for (xyz=0; xyz<3; ++xyz) {
        ++cnt;
        /* if a loner, don't try to adapt */
        if (loner) {
          if ( get_irrep_xyz( cartrep, xyz ) == irrep) {
            salc_orig[irrep][cnt][3*atom+xyz] = 1.0;
          }
        }
        else {
          diag_ind = xyz*3 + xyz;
          for (op=0; op < nirreps; ++op) { 
            atom2 = ict[op][atom] - 1;
            for (i=0; i<3*natom; ++i)
              salc_orig[irrep][cnt][3*atom2+xyz] += cartrep[op][diag_ind] / syminfo.ct[irrep][op];
          }
        }
      }
    }
    nsalc_orig[irrep] = cnt+1;
  }
  if (print) {
    for (irrep=0; irrep < nirreps; ++irrep) {
      fprintf(outfile,"Cartesian SALCs for Irrep %d\n", irrep);
      print_mat(salc_orig[irrep], nsalc_orig[irrep], 3*natom, outfile);
    }
  }

  /**** Schmidt orthogonalize coordinates to remove COM/ROT constraints ****/
  salc = (double ***) malloc(nirreps*sizeof(double **));
  nsalc = init_int_array(nirreps);
  v = init_array(3*natom);
  nsalc_all = 0;
  for (irrep=0; irrep< nirreps; ++irrep) {
    salc[irrep] = block_matrix(nsalc_orig[irrep],3*natom);
    cnt = 0;

    for (i=0; i<nsalc_orig[irrep]; ++i) {
      zero_array(v,3*natom);

      for (I=0; I<3*natom; I++)
        v[I] = salc_orig[irrep][i][I];
 
      for (j=0; j<nconstraints; j++) {
        dot_array(salc_orig[irrep][i], constraints[j], 3*natom, &dotval) ;
        for (I=0; I<3*natom; I++)
          v[I] -= dotval * constraints[j][I] ;
      }

      dot_array(v, v, 3*natom, &normval);
      if (normval > 1.0E-10) {
        for (j=0; j<3*natom; ++j)
          v[j] = v[j] / sqrt(normval) ;
  
        cnt += schmidt_add(salc[irrep], cnt, 3*natom, v);
      }
    }
    nsalc[irrep] = cnt;
    nsalc_all += cnt;
  }
  for (irrep=0; irrep<nirreps; ++irrep) {
    free_block(salc_orig[irrep]);
  }
  free_block(constraints);
  free_int_array(nsalc_orig);
  free_array(v);

  if (optinfo.freq_irrep != -1) {
    for (irrep=0; irrep < nirreps; ++irrep)
      if (irrep != (optinfo.freq_irrep - 1) ) {
        nsalc_all -= nsalc[irrep];
        nsalc[irrep] = 0;
      }
  }

  if (print) {
    for (irrep=0; irrep < nirreps; ++irrep)
      if (nsalc[irrep]) {
        fprintf(outfile,"Symmetry Adapted Cartesian Coordinates, irrep %d\n", irrep);
        print_mat(salc[irrep], nsalc[irrep], 3*natom, outfile);
      }
  }

  check_coordinates(natom,coord,masses,Zvals,nsalc,salc);

  /* determine number of displacements */
  ndisp_all = 0;
  ndisp = init_int_array(nirreps);
  for (irrep=0; irrep<nirreps; ++irrep) {

    /* diagonal displacements */
    if (irrep == 0) {
      if (optinfo.points == 3)
        ndisp[irrep] = 2 * nsalc[irrep];
      else if (optinfo.points == 5)
        ndisp[irrep] = 4 * nsalc[irrep];
    }
    else {
      if (optinfo.points == 3)
        ndisp[irrep] = nsalc[irrep];
      else if (optinfo.points == 5)
        ndisp[irrep] = 2* nsalc[irrep];
    }

    /* off-diagonal displacements */
    if (optinfo.points == 3)
      ndisp[irrep] += 2 * nsalc[irrep] * (nsalc[irrep] - 1) / 2;
    else if (optinfo.points == 5)
      ndisp[irrep] += 8 * nsalc[irrep] * (nsalc[irrep] - 1) / 2;

    ndisp_all += ndisp[irrep];
  }

  for (irrep=0; irrep<nirreps;++irrep)
    fprintf(outfile,"ndisp[%d]: %d\n",irrep,ndisp[irrep]);
  fflush(outfile);

  /*** construct displaced geometries ***/
  geoms = (double ***) malloc(nirreps*sizeof(double **));
  for (irrep=0; irrep<nirreps; ++irrep) {

    geoms[irrep] = block_matrix(ndisp[irrep],3*natom);
    for (i=0; i<ndisp[irrep]; ++i)
      for (j=0; j<3*natom; ++j)
        geoms[irrep][i][j] = coord[j];

    /* construct diagonal displacements */
    if (irrep == 0) {
      cnt = 0;
      for (i=0; i<nsalc[irrep]; ++i) {
        if (optinfo.points == 3) {
          for (j=0; j < 3*natom; ++j) {
            geoms[irrep][cnt+0][j] -= optinfo.disp_size * salc[irrep][i][j] / sqrt(masses[j]);
            geoms[irrep][cnt+1][j] += optinfo.disp_size * salc[irrep][i][j] / sqrt(masses[j]);
          }
          cnt += 2;
        }
        else if (optinfo.points == 5) {
          for (j=0; j < 3*natom; ++j) {
            geoms[irrep][cnt+0][j] -= 2.0 * optinfo.disp_size * salc[irrep][i][j] / sqrt(masses[j]);
            geoms[irrep][cnt+1][j] -= 1.0 * optinfo.disp_size * salc[irrep][i][j] / sqrt(masses[j]);
            geoms[irrep][cnt+2][j] += 1.0 * optinfo.disp_size * salc[irrep][i][j] / sqrt(masses[j]);
            geoms[irrep][cnt+3][j] += 2.0 * optinfo.disp_size * salc[irrep][i][j] / sqrt(masses[j]);
          }
          cnt += 4;
        }
      }
    }
    else { /* asymmetric diagonal displacements */
      cnt = 0;
      for (i=0; i<nsalc[irrep]; ++i) {
        if (optinfo.points == 3) {
          for (j=0; j < 3*natom; ++j)
            geoms[irrep][cnt][j]   -= optinfo.disp_size * salc[irrep][i][j] / sqrt(masses[j]);
          cnt += 1;
        }
        else if (optinfo.points == 5) {
          for (j=0; j < 3*natom; ++j) {
            geoms[irrep][cnt+0][j] -= 2.0 * optinfo.disp_size * salc[irrep][i][j] / sqrt(masses[j]);
            geoms[irrep][cnt+1][j] -= 1.0 * optinfo.disp_size * salc[irrep][i][j] / sqrt(masses[j]);
          }
          cnt += 2;
        }
      }
    }
    /* add off-diagonal displacements */

    for (i=0; i<nsalc[irrep]; ++i) {
      for (j=0; j<i; ++j) {
        if (optinfo.points == 3) {
          for (k=0; k < 3*natom; ++k) {
            geoms[irrep][cnt+0][k] += ( + salc[irrep][i][k] + salc[irrep][j][k] )
              * optinfo.disp_size / sqrt(masses[k]);
            geoms[irrep][cnt+1][k] += ( - salc[irrep][i][k] - salc[irrep][j][k] )
              * optinfo.disp_size / sqrt(masses[k]);
          }
          cnt += 2;
        }
        else if (optinfo.points == 5) {
          for (k=0; k < 3*natom; ++k) {
            geoms[irrep][cnt+0][k] += ( - 1.0 * salc[irrep][i][k] - 2.0 * salc[irrep][j][k] )
              * optinfo.disp_size / sqrt(masses[k]);
            geoms[irrep][cnt+1][k] += ( - 2.0 * salc[irrep][i][k] - 1.0 * salc[irrep][j][k] )
              * optinfo.disp_size / sqrt(masses[k]);
            geoms[irrep][cnt+2][k] += ( - 1.0 * salc[irrep][i][k] - 1.0 * salc[irrep][j][k] )
              * optinfo.disp_size / sqrt(masses[k]);
            geoms[irrep][cnt+3][k] += ( + 1.0 * salc[irrep][i][k] - 1.0 * salc[irrep][j][k] )
              * optinfo.disp_size / sqrt(masses[k]);
            geoms[irrep][cnt+4][k] += ( - 1.0 * salc[irrep][i][k] + 1.0 * salc[irrep][j][k] )
              * optinfo.disp_size / sqrt(masses[k]);
            geoms[irrep][cnt+5][k] += ( + 1.0 * salc[irrep][i][k] + 1.0 * salc[irrep][j][k] )
              * optinfo.disp_size / sqrt(masses[k]);
            geoms[irrep][cnt+6][k] += ( + 2.0 * salc[irrep][i][k] + 1.0 * salc[irrep][j][k] )
              * optinfo.disp_size / sqrt(masses[k]);
            geoms[irrep][cnt+7][k] += ( + 1.0 * salc[irrep][i][k] + 2.0 * salc[irrep][j][k] )
              * optinfo.disp_size / sqrt(masses[k]);
          }
          cnt += 8;
        }
      }
    }
  }
  free_array(masses);

  if(optinfo.external_energies){
    int counter = 1;
    FILE *fp_dispcart = fopen("dispcart", "w");
    // ACS (11/06) Print the displaced and reference geometry for external programs
    // Use 1 based counting, and remember that the first "displacement" is the reference itself
    fprintf(fp_dispcart,"%d\n",counter++);
    for(int n = 0; n < carts.get_natom(); n++){
      fprintf(fp_dispcart,"%16.10lf %16.10lf %16.10lf\n",coord[3*n],coord[3*n+1],coord[3*n+2]);
    }

    for (irrep=0; irrep<nirreps; ++irrep) {
      for(int disp = 0; disp< ndisp[irrep]; disp++){
        fprintf(fp_dispcart,"%d\n",counter++);
        for(int n = 0; n < carts.get_natom(); n++){
          fprintf(fp_dispcart,"%16.10lf %16.10lf %16.10lf\n",
            geoms[irrep][disp][3*n],geoms[irrep][disp][3*n+1],geoms[irrep][disp][3*n+2]);
        }
      }
    }
    fclose(fp_dispcart);
  }

  /**** write out info to PSIF_OPTKING, geoms and coordinates in irrep order ****/
  open_PSIF();

  big_zeroes = init_array(ndisp_all*3*natom);
  psio_write_entry(PSIF_OPTKING, "OPT: Displaced geometries", (char *) &(big_zeroes[0]),
      ndisp_all*3*natom*sizeof(double));
  psio_write_entry(PSIF_OPTKING, "OPT: Adapted cartesians", (char *) &(big_zeroes[0]),
      ndisp_all*3*natom*sizeof(double));
  psio_write_entry(PSIF_OPTKING, "OPT: Displaced gradients", (char *) &(big_zeroes[0]),
      ndisp_all*3*natom*sizeof(double));
  free_array(big_zeroes);

  next = PSIO_ZERO;
  for (irrep=0; irrep<nirreps; ++irrep) {
    if (ndisp[irrep]) {
      psio_write(PSIF_OPTKING, "OPT: Displaced geometries",
          (char *) &(geoms[irrep][0][0]), ndisp[irrep]*3*natom*sizeof(double), next, &next);
//fprintf(outfile,"Displaced Geometries of irrep %d.\n",irrep);
//print_mat(geoms[irrep],ndisp[irrep],3*natom,outfile);
    }
    free_block(geoms[irrep]);
  }
  free(geoms);

  next = PSIO_ZERO;
  for (irrep=0; irrep<nirreps; ++irrep) {
    if (nsalc[irrep]) {
      psio_write(PSIF_OPTKING, "OPT: Adapted cartesians", (char *) &(salc[irrep][0][0]),
      nsalc[irrep]*3*natom*sizeof(double), next, &next);
    }
    free_block(salc[irrep]);
  }
  free(salc);

  /* make list of irreps for each disp used to make prefixes */
  irrep_per_disp = (int *) malloc(ndisp_all*sizeof(int));
  cnt = -1;
  for (irrep=0; irrep<nirreps; ++irrep) {
    for (i=0; i<ndisp[irrep]; ++i)
      irrep_per_disp[++cnt] = irrep;
  }
  fprintf(outfile,"Irreps per disp: ");
  for (i=0; i<ndisp_all; ++i)
    fprintf(outfile,"%d ",irrep_per_disp[i]);
  fprintf(outfile,"\n");

  /**** Write parameters to PSIO file ****/
  optinfo.disp_num = 0;
  disp_e = init_array(ndisp_all);
  optinfo.micro_iteration = 0;

  psio_write_entry(PSIF_OPTKING, "OPT: Displaced energies",
      (char *) &(disp_e[0]), ndisp_all*sizeof(double));
  psio_write_entry(PSIF_OPTKING, "OPT: Reference geometry",
      (char *) &(coord[0]), 3*natom*sizeof(double));
  psio_write_entry(PSIF_OPTKING, "OPT: Reference energy",
      (char *) &(energy), sizeof(double));
  psio_write_entry(PSIF_OPTKING, "OPT: Current disp_num",
      (char *) &(optinfo.disp_num),sizeof(int));
  psio_write_entry(PSIF_OPTKING, "OPT: Num. of disp.",
      (char *) &(ndisp_all), sizeof(int));
  psio_write_entry(PSIF_OPTKING, "OPT: Disp. per irrep",
      (char *) &(ndisp[0]), nirreps*sizeof(int));
  psio_write_entry(PSIF_OPTKING, "OPT: Irrep per disp",
      (char *) &(irrep_per_disp[0]), ndisp_all*sizeof(int));
  psio_write_entry(PSIF_OPTKING, "OPT: Num. of coord.",
      (char *) &(nsalc_all), sizeof(int));
  psio_write_entry(PSIF_OPTKING, "OPT: Coord. per irrep",
      (char *) &(nsalc[0]), nirreps*sizeof(int));
  psio_write_entry(PSIF_OPTKING, "Micro_iteration",
      (char *) &(optinfo.micro_iteration),sizeof(int));
  close_PSIF();

  for (irrep=0; irrep<nirreps; ++irrep)
    fprintf(outfile,"Number of %s displaced geometries is %d.\n",syminfo.irrep_lbls[irrep],ndisp[irrep]);
  fprintf(outfile,"Total number of displaced geometries is %d.\n",ndisp_all);

  free_int_array(nsalc);
  free_int_array(ndisp);
  free_array(coord);
  free_array(disp_e);
  free_int_array(ua2a);
  return(ndisp_all);
}

}//} /* namespace psi::optking */

