/*###########################################################################*/
/*! 
  \file
  \ingroup EXTREMA
  \brief zmat transformation functions

  z-matrix derived class member functions for transformations between 
  z-matrix and cartesian coordinates. */
/*						Joseph P. Kenny 11/29/01
  ##########################################################################*/

#define EXTERN
#include "extrema.h"
#include "inline.h"

using namespace psi::extrema;

/*---------------------------------------------------------------------------*/
/*! \fn zmat::compute_B()
  \brief Computes the B matrix.

  Interface for <b>internals::B_row...()</b> functions. */
/*---------------------------------------------------------------------------*/

void zmat::compute_B() {  

    int i, j, pos=0;
    double *B_row0, *B_row1, *B_row2;

  for(i=1;i<num_entries;++i) {
      if(i==1) {
	  B_row0 = B_row_bond(carts, i, simples[pos].get_bond()-1);
	  B[pos] = B_row0;
          ++pos;
	}
      if(i==2) {
	  B_row0 = B_row_bond(carts, i, simples[pos].get_bond()-1);
	  B_row1 = B_row_angle(carts, i, simples[pos+1].get_bond()-1, 
			       simples[pos+1].get_angle()-1);
	  B[pos] = B_row0;
	  B[pos+1] = B_row1;
	  pos += 2;
	}
      if(i>2) {
	  B_row0 = B_row_bond(carts, i, simples[pos].get_bond()-1);
	  B_row1 = B_row_angle(carts, i, simples[pos+1].get_bond()-1, 
			       simples[pos+1].get_angle()-1);
	  B_row2 = B_row_tors(carts, i, simples[pos+2].get_bond()-1, 
			      simples[pos+2].get_angle()-1, 
			      simples[pos+2].get_tors()-1);
	  B[pos] = B_row0;
	  B[pos+1] = B_row1;
	  B[pos+2] = B_row2;
	  pos += 3;
	}
    }

  /*form u*/
  int k=0;
  for(j=0;j<num_entries;++j) {
      if(strncmp(felement[j],"X\0",2)) {
	 u[3*j][3*j] = 1.0 / masses[k]; 
	 u[3*j+1][3*j+1] = 1.0 /masses[k];
	 u[3*j+2][3*j+2] = 1.0 / masses[k];
	 ++k;
       }
       else if (!strncmp(felement[j],"X\0",2)) {
	  u[3*j][3*j] = u[3*j+1][3*j+1] = u[3*j+2][3*j+2]= 1.0;
	  }
  }

  /*form B_red, the reduced dimension B matrix*/

  pos=0;
  for(i=0;i<fnum_coords;++i) 
      if(first_unique[i] && simples[i].get_opt() ) {
	  for(j=0;j<(3*num_entries);++j)
	      B_red[pos][j] = B[i][j];
	  ++pos;
      }  
      
  return;
}



/*---------------------------------------------------------------------------*/
/*! \fn zmat::grad_trans()
  \brief Transforms gradients from cartesians to internals.

  Interface for <b>internals::grad_trans()</b>.  Gradients are calculated
  for redundant coordinates and their equivalency is checked. */
/*---------------------------------------------------------------------------*/

void zmat::grad_trans() {
    
    internals::grad_trans();

    int i, j, div, p=0;

    if(print_lvl>NORMAL_PRINT) {
	fprintf(outfile,"\n  Internal coordinate gradients (a.u):");
	for(i=0;i<fnum_coords;++i)
	    fprintf(outfile,"\n  %i %8s: %15.10lf",
		    i+1,simples[i].get_label(),fgrads[i]);
    }
    fprintf(outfile,"\n");

    double sum;
    for(i=0;i<fnum_coords;++i) 
	if(first_unique[i] && simples[i].get_opt()) {
	    sum = fgrads[i];
	    div=1;
	    for(j=(i+1);j<fnum_coords;++j)
		if(!strcmp(simples[i].get_label(),simples[j].get_label()) &&
		    strcmp(simples[i].get_label(),"")) {
		    sum += fgrads[j];
		    ++div;
		    if( fabs(fgrads[i]-fgrads[j]) > EQUIV_GRAD ) {
			fprintf(outfile,"\n  WARNING: gradients %d and %d",
				i+1,j+1);
			fprintf(outfile," should be equal, differ by %lf",
				fabs(fgrads[i]-fgrads[j]));
		    }
		}
	    grads[p] = sum/((double) div);
	    ++p;
	}
    
    return;
}



/*---------------------------------------------------------------------------*/
/*! \fn zmat::back_transform()
  \brief Computes cartesians corresponding to the current z-matrix.
 
  Interface for <b>internals::back_transform()</b>. */
/*---------------------------------------------------------------------------*/

void zmat::back_transform() {
    
    /*create full coordinate vectors*/
    double *fcoord_new;
    fcoord_new = init_array(fnum_coords);
    
    int i;
    for(i=0;i<fnum_coords;++i) { 
	fcoord_new[i] = simples[i].get_val();
    }

    internals::back_transform(fcoord_new, fcoord_old);

    free(fcoord_old);
    free(fcoord_new);
    return;
}



/*---------------------------------------------------------------------------*/
/*! \fn  zmat::cart_to_internal(double* z_array)
  \brief  Computes z-matrix from cartesian coordinates. 
  \param z_array values of z-matrix coordinates  */
/*--------------------------------------------------------------------------*/

void zmat :: cart_to_internal(double** z_array) {

    int i, j,pos=0;
    double *temp1, *temp2, temp_num, div, n1, n2;
    double tnum;

    temp1 = init_array(3);
    temp2 = init_array(3);



    for(i=1;i<num_entries;++i) {


	if(i==1) {
	    (*z_array)[pos] = compute_bond(carts, i,simples[pos].get_bond()-1);
	    ++pos ;
	}


	if(i==2){

		(*z_array)[pos] = 
		    compute_bond( carts, i, simples[pos].get_bond()-1);
		(*z_array)[pos+1] = 
		    compute_angle( carts, i, simples[pos+1].get_bond()-1,
				   simples[pos+1].get_angle()-1);
	    pos += 2;
	}


	if(i>2) {
	    
	   (*z_array)[pos] = 
	       compute_bond( carts, i, simples[pos].get_bond()-1);
	   (*z_array)[pos+1] = 
	       compute_angle( carts, i, simples[pos+1].get_bond()-1,
			      simples[pos+1].get_angle()-1);
	   (*z_array)[pos+2] = 
	       compute_torsion( carts, i, simples[pos+2].get_bond()-1,
				simples[pos+2].get_angle()-1,
				simples[pos+2].get_tors()-1);
	   pos += 3;
	}
		
    }

    for(i=0;i<fnum_coords;++i) {
	if( (simples[i].get_type() == 2) && (simples[i].get_val() < 0.0) ) 
	    (*z_array)[i] *= -1.0;
    }

    return;
}								     
