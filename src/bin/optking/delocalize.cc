/*! \file
    \ingroup OPTKING
    \brief This function forms a set of delocalized, symmetry-adapted
     internal coordinates, given a set of simples and their
     s-vectors.  The new coordinates will not mix coordinates of
     different types if mix_types == 1
*/

#define EXTERN
#include "globals.h"
#undef EXTERN
#include "cartesians.h"
#include "simples.h"
#include "salc.h"
#include "opt.h"

namespace psi { namespace optking {

void delocalize(const simples_class & simples, const cartesians & carts) {
  int error,i,j,k,a,b,c,d,id,count,sub_index,row[5],dim[5];
  int rotor_type, degrees_of_freedom, col, natom;
  double **stre_mat, **bend_mat, **tors_mat, **out_mat;
  double **evectst, **evectst_symm, **coord_symm, *fmass;
  double *evals, **evects, **B, **uBt, **BBt, **u, **temp_mat;

  natom = optinfo.natom;

  // Build B matrix for simples
  B = init_matrix(simples.get_num(),natom*3);
  count = -1;
  for (i=0; i<simples.get_num(STRE); ++i) {
    a = simples.get_atom(STRE, i, 0);
    b = simples.get_atom(STRE, i, 1);
    ++count;
    for (k=0;k<3;++k) {
      B[count][3*a+k] += simples.get_s(STRE, i, 0, k);
      B[count][3*b+k] += simples.get_s(STRE, i, 1, k);
    }
  }
  for (i=0; i<simples.get_num(BEND); ++i) {
    a = simples.get_atom(BEND, i, 0);
    b = simples.get_atom(BEND, i, 1);
    c = simples.get_atom(BEND, i, 2);
    ++count;
    for (k=0;k<3;++k) {
      B[count][3*a+k] += simples.get_s(BEND, i, 0, k);
      B[count][3*b+k] += simples.get_s(BEND, i, 1, k);
      B[count][3*c+k] += simples.get_s(BEND, i, 2, k);
    }
  }
  for (i=0; i<simples.get_num(TORS); ++i) {
    a = simples.get_atom(TORS, i, 0);
    b = simples.get_atom(TORS, i, 1);
    c = simples.get_atom(TORS, i, 2);
    d = simples.get_atom(TORS, i, 3);
    ++count;
    for (k=0;k<3;++k) {
      B[count][3*a+k] += simples.get_s(TORS, i, 0, k);
      B[count][3*b+k] += simples.get_s(TORS, i, 1, k);
      B[count][3*c+k] += simples.get_s(TORS, i, 2, k);
      B[count][3*d+k] += simples.get_s(TORS, i, 3, k);
    }
  }
  for (i=0; i<simples.get_num(OUT); ++i) {
    a = simples.get_atom(OUT, i, 0);
    b = simples.get_atom(OUT, i, 1);
    c = simples.get_atom(OUT, i, 2);
    d = simples.get_atom(OUT, i, 3);
    ++count;
    for (k=0;k<3;++k) {
      B[count][3*a+k] += simples.get_s(OUT, i, 0, k);
      B[count][3*b+k] += simples.get_s(OUT, i, 1, k);
      B[count][3*c+k] += simples.get_s(OUT, i, 2, k);
      B[count][3*d+k] += simples.get_s(OUT, i, 3, k);
    }
  }
  for (i=0; i<simples.get_num(LINB); ++i) {
    a = simples.get_atom(LINB, i, 0);
    b = simples.get_atom(LINB, i, 1);
    c = simples.get_atom(LINB, i, 2);
    ++count;
    for (k=0;k<3;++k) {
      B[count][3*a+k] += simples.get_s(LINB, i, 0, k);
      B[count][3*b+k] += simples.get_s(LINB, i, 1, k);
      B[count][3*c+k] += simples.get_s(LINB, i, 2, k);
    }
  }
  //  print_mat2(B,simples.get_num(),natom*3,outfile);

  uBt = init_matrix(3*natom, simples.get_num());
  fmass = carts.get_fmass();
  u = mass_mat(fmass);
  free_array(fmass);

  // Form BBt matrix
  BBt = init_matrix(simples.get_num(),simples.get_num());
  if (optinfo.mix_types) {
    opt_mmult(u,0,B,1,uBt,0,3*natom,3*natom,simples.get_num(),0);

    opt_mmult(B,0,uBt,0,BBt,0,simples.get_num(),3*natom,simples.get_num(),0);

    // opt_mmult(B,0,B,1,BBt,0,simples.get_num(),natom*3,simples.get_num(),0);
  }
  else {
    // Make BBt block diagonal by multiplying only stre*stre, etc.
    dim[0] = simples.get_num(STRE);
    dim[1] = simples.get_num(BEND);
    dim[2] = simples.get_num(TORS);
    dim[3] = simples.get_num(OUT);
    dim[4] = simples.get_num(LINB);
    row[0] = 0;
    row[1] = dim[0];
    row[2] = row[1]+dim[1];
    row[3] = row[2]+dim[2];
    row[4] = row[3]+dim[3];
    double **ptr;
    ptr = (double **) malloc(simples.get_num()*sizeof(double *));
    for (i=0;i<5;++i) {
      for (j=0;j<dim[i];++j) {
        ptr[j] = BBt[row[i]+j] + row[i];
      }
      if (dim[i] != 0)

        opt_mmult(u,0,&(B[row[i]]),1,uBt,0,3*natom,3*natom,dim[i],0);
      opt_mmult(&(B[row[i]]),0,uBt,0,ptr,0,dim[i],3*natom,dim[i],0);

      //  opt_mmult(&(B[row[i]]),0,&(B[row[i]]),1,ptr,0,dim[i],natom*3,dim[i],0);
    }
    free(ptr);
  }
  free_matrix(u);
  free_matrix(B);
  free_matrix(uBt);

  //fprintf(outfile,"The BB^t Matrix:\n");
  //print_mat2(BBt,simples.get_num(),simples.get_num(),outfile);

  // Diagonalize BBt
  evals = init_array(simples.get_num());
  evects = init_matrix(simples.get_num(),simples.get_num());
  opt_sq_rsp(simples.get_num(),simples.get_num(), BBt, evals, 1, evects, 1.0E-14);

  /* check eigenvectors of BBt */
  if (optinfo.print_delocalize) {
    fprintf(outfile,"\n\n+++Delocalized Coordinate Formation+++\n");
    fprintf(outfile,"\n\nChecking eigenvectors of BBt...\n");
    eivout(evects, evals, simples.get_num(), simples.get_num(), outfile);
  }
  temp_mat = init_matrix(simples.get_num(),simples.get_num());
  opt_mmult(BBt,0,evects,0,temp_mat,0,simples.get_num(),simples.get_num(),simples.get_num(),0);
  for (j=0;j<simples.get_num();++j) {
    error = 0;
    for (i=0;i<simples.get_num();++i) {
      if ( fabs(temp_mat[i][j] - evals[j]*evects[i][j]) > 1.0E-13)  error = 1;
      if (error == 1) { fprintf(outfile,"Error in evect %d\n",j); error = 0;}
    }
  }
  free_matrix(temp_mat);
  free_matrix(BBt);

  /* check for proper number of non-redundant coordinates (3n-6) */ 
  num_nonzero = 0;
  for (i=0;i<simples.get_num();++i)
    if( evals[i] > optinfo.ev_tol ) ++num_nonzero;

  chkpt_init(PSIO_OPEN_OLD);
  rotor_type = chkpt_rd_rottype();
  chkpt_close();

  switch (rotor_type) {
    case 3:
      degrees_of_freedom = 3 * natom - 5;
      break;
    case 6:
      degrees_of_freedom = 0;
      break;
    default:
      degrees_of_freedom = 3 * natom - 6;
      break;
  }
  if (num_nonzero == degrees_of_freedom) {
    fprintf(outfile,"\nGood: # of delocalized coordinates = # of degrees");
    fprintf(outfile," of freedom.\n");
  }
  else if (num_nonzero < degrees_of_freedom) {
    fprintf(outfile,"# of delocalized coordinates < # of degrees of ");
    fprintf(outfile,"freedom!\n You may need more simple internals.\n");
  }
  else if (num_nonzero > degrees_of_freedom) {
    fprintf(outfile,"# of delocalized coordinates > # of degrees of ");
    fprintf(outfile,"freedom!\n You might try increasing EV_TOL.\n");
    rm_rotations(simples, carts, num_nonzero, evects);
    fprintf(outfile,"Eigenvectors of BBt after removal of rotations\n");
    eivout(evects, evals, simples.get_num(), simples.get_num(), outfile);
  }


  /* transpose evects matrix, also throw out redundant coordinates
     (eigenvectors corresponding to zero eigenvalues*/

  char aline[MAX_LINELENGTH], **buffer, *err;
  int h;
  irr = init_int_array(simples.get_num());

  evectst = init_matrix(num_nonzero,simples.get_num());

  for (i=0;i<simples.get_num();++i) {
    for (j=simples.get_num()-num_nonzero;j<simples.get_num();++j) {
      evectst[j-simples.get_num()+num_nonzero][i] = evects[i][j];
    }
  }
  if (optinfo.print_delocalize == 1) {
    fprintf(outfile,"\nNon-redundant delocalized internal coordinates");
    fprintf(outfile,"(each row is a coordinate):\n");
    print_mat(evectst,num_nonzero,simples.get_num(),outfile);
  }
  free_array(evals);
  free_matrix(evects);

  /* send evectst to irrep.cc to be properly symmetrized */
  evectst_symm = irrep(simples, evectst);
  free_matrix(evectst);

  if (optinfo.print_delocalize == 1) {
    fprintf(outfile,"\nSymmetrized evects\n");
    print_mat(evectst_symm,num_nonzero,simples.get_num(),outfile);
  }

  // print out coordinates to intco.dat
  opt_ffile(&fp_intco, "intco.dat", 2);
  count = 0;
  for( ; ; ) {
    err = fgets(aline, MAX_LINELENGTH, fp_intco);
    if (err == NULL) break;
    ++count;
  }
  rewind(fp_intco);

  /* read all but the last line of the file into memory... */
  buffer = (char **) malloc((count-1) * sizeof(char *));
  for(i=0; i < count-1; i++) {
    buffer[i] = (char *) malloc(MAX_LINELENGTH * sizeof(char));
    err = fgets(buffer[i], MAX_LINELENGTH, fp_intco);
  }
  rewind(fp_intco);

  /* ...and overwite */
  for(i=0; i < count-1; i++) {
    fprintf(fp_intco, "%s", buffer[i]);
    free(buffer[i]);
  }
  free(buffer);

  // Print out coordinates to intco.dat
  fprintf(fp_intco,"  symm = ( \n");
  for (h=0; h<syminfo.nirreps; ++h) {
    if (h==1) fprintf(fp_intco, "  )\n  asymm = (\n");
    for (i=0;i<num_nonzero;++i) {
      if (irr[i] == h) {

        // print irrep label
        fprintf(fp_intco,"    (");
        fprintf(fp_intco,"\"%s\"",syminfo.clean_irrep_lbls[irr[i]]);

        // print simples ids
        fprintf(fp_intco," (");
        for (col=0, j=0;j<simples.get_num();++j) {
          if ( fabs(evectst_symm[i][j]) > 1.0E-10 ) {
            fprintf(fp_intco," %d",simples.index_to_id(j));
            if (col == 20) { fprintf(fp_intco,"\n    "); col = -1; }
            ++col;
          }
        }
        fprintf(fp_intco,")\n");

        // print coefficients of salc
        fprintf(fp_intco,"   (");
        for (col=0, j=0;j<simples.get_num();++j) {
          if ( fabs(evectst_symm[i][j]) > 1.0E-10 ) {
           fprintf(fp_intco,"%15.7lf",evectst_symm[i][j]);
          //fprintf(fp_intco,"%15.10lf",evectst_symm[i][j]);
          if (col == 7) { fprintf(fp_intco,"\n    "); col = -1; }
          ++col;
         }
        }
        fprintf(fp_intco,"))\n");
      }
    }
  }

  fprintf(fp_intco,"  )\n)\n"); 
  fflush(fp_intco);
  fclose(fp_intco);

  free_matrix(evectst_symm);
  return;
}

// removes rotations - returns number of rotations removed
void rm_rotations(const simples_class & simples, const cartesians & carts,
    int &num_nonzero, double **evects) {

  int i, j, k, a, b, c, d, ivect, cnt;
  double *disp_coord, scale=0.001, *fatomic_num, energy, disp_energy;
  double *coord, rot_tol = 1.0E-10;

  fatomic_num = carts.get_fatomic_num();
  coord = carts.get_coord();

  energy = nuclear_repulsion(fatomic_num,coord);

  disp_coord = new double [3*optinfo.natom];
  for (ivect=0; ivect<num_nonzero; ++ivect) {

    for (i=0;i<(3*optinfo.natom);++i)
      disp_coord[i] = coord[i];

    cnt = -1;
    for (i=0; i<simples.get_num(STRE); ++i) {
      a = simples.get_atom(STRE, i, 0);
      b = simples.get_atom(STRE, i, 1);
      ++cnt;
      for (k=0;k<3;++k) {
        disp_coord[3*a+k] += scale * evects[cnt][ivect] * simples.get_s(STRE, i, 0, k);
        disp_coord[3*b+k] += scale * evects[cnt][ivect] * simples.get_s(STRE, i, 1, k);
      }
    }
    for (i=0; i<simples.get_num(BEND); ++i) {
      a = simples.get_atom(BEND, i, 0);
      b = simples.get_atom(BEND, i, 1);
      c = simples.get_atom(BEND, i, 2);
      ++cnt;
      for (k=0;k<3;++k) {
        disp_coord[3*a+k] += scale * evects[cnt][ivect] * simples.get_s(BEND, i, 0, k);
        disp_coord[3*b+k] += scale * evects[cnt][ivect] * simples.get_s(BEND, i, 1, k);
        disp_coord[3*c+k] += scale * evects[cnt][ivect] * simples.get_s(BEND, i, 2, k);
      }
    }
    for (i=0; i<simples.get_num(TORS); ++i) {
      a = simples.get_atom(TORS, i, 0);
      b = simples.get_atom(TORS, i, 1);
      c = simples.get_atom(TORS, i, 2);
      d = simples.get_atom(TORS, i, 3);
      ++cnt;
      for (k=0;k<3;++k) {
        disp_coord[3*a+k] += scale * evects[cnt][ivect] * simples.get_s(TORS, i, 0, k);
        disp_coord[3*b+k] += scale * evects[cnt][ivect] * simples.get_s(TORS, i, 1, k);
        disp_coord[3*c+k] += scale * evects[cnt][ivect] * simples.get_s(TORS, i, 2, k);
        disp_coord[3*d+k] += scale * evects[cnt][ivect] * simples.get_s(TORS, i, 3, k);
      }
    }
    for (i=0; i<simples.get_num(OUT); ++i) {
      a = simples.get_atom(OUT, i, 0);
      b = simples.get_atom(OUT, i, 1);
      c = simples.get_atom(OUT, i, 2);
      d = simples.get_atom(OUT, i, 3);
      ++cnt;
      for (k=0;k<3;++k) {
        disp_coord[3*a+k] += scale * evects[cnt][ivect] * simples.get_s(OUT, i, 0, k); 
        disp_coord[3*b+k] += scale * evects[cnt][ivect] * simples.get_s(OUT, i, 1, k);
        disp_coord[3*c+k] += scale * evects[cnt][ivect] * simples.get_s(OUT, i, 2, k);
        disp_coord[3*d+k] += scale * evects[cnt][ivect] * simples.get_s(OUT, i, 3, k);
      }
    }
    disp_energy = nuclear_repulsion(fatomic_num, disp_coord);
    // fprintf(outfile,"dispenergy - energy: %16.12lf\n",disp_energy - energy);
    if (fabs(disp_energy - energy) < rot_tol) {
      fprintf(outfile,"rotational coordinate eliminated");
      for (i=ivect+1; i<num_nonzero; ++i) 
        for (j=0; j<simples.get_num(); ++j) 
          evects[i-1][j] = evects[i][j];

      --ivect;
      --num_nonzero;
    }
  }
  free_array(coord);
  delete [] disp_coord;
  delete [] fatomic_num;
  return;
}

}} /* namespace psi::optking */

