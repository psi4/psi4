/*! \file 
    \ingroup (INPUT)
    \brief Enter brief description of file here 
*/
#define EXTERN
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include "input.h"
#include "global.h"
#include "defines.h"

namespace psi { namespace input {

/* read_zmat(): parse the zmat = () keyword in the input file.
**
** Written by EFV, original date unknown
**
** Modified 8/03 by TDC to allow a simpler z-matrix format without all the 
** "extra" parentheses.  This should allow users to cut-and-paste z-matrices from
** other programs for use in PSI.
**
** The following H2O examples should both be valid:
**
** %original format
** zmat = (
**   (o)
**   (h 1 r)
**   (h 1 r 2 a)
** )
**
** %new format (8/03)
** zmat = (
**   o
**   h 1 r
**   h 1 r 2 a
** )
*/

/*this is used to hold ZVARS info*/
struct definition{
      char* variable;
      double value;
};

/*zmat parser, code is found at bottom of this file*/
void parse_zmat(int i, int position, double *value,struct definition *array, 
		int num_vals, int zval_exist);
void parse_zmat_simple(int i, int position, double *value,struct definition *array, 
		       int num_vals, int zval_exist);

void read_zmat()
{
  int i, j, k, m, a, b, c, errcod, value_set, zvar_exist;
  char *buffer;
  int num_vals, entry_length=0, atomcount, fatomcount;
  int A, B, C, D;
  int linearOn = 1;	/* Flag indicating the code is working on the linear fragment in the begginning of the Z-matrix */
  double rAB, rBC, rCD, thetaABC, thetaBCD, phiABCD, val, norm1, norm2;
  double cosABC, sinABC, cosBCD, sinBCD, cosABCD, sinABCD;
  double eAB[3], eBC[3], ex[3], ey[3], ez[3];
  double Z = 0.0;

  char *name;

  int zmat_len, simple_zmat, simple_zvars;

  /*these are passed to parsing function to tell it what position in the row to look at*/
  int bond=2, angle=4, tors=6;

  /* array to store variable definitions */
  struct definition *def_arr;
  
  /*need this to avoid seg fault if character string used but no zvar given*/
  zvar_exist = ip_exist("ZVARS",0);    

  /* must first determine "type" of z-matrix input: nested or simple */
  simple_zmat = 1;
  ip_count("ZMAT", &zmat_len,0);
  for(i=0; i < zmat_len; i++) {  
    ip_count("ZMAT", &entry_length,1,i);
    if(entry_length > 1) {
      simple_zmat = 0;
    }
  }
  
  /* Read number of lines and count atoms in ZMAT */

  num_atoms = num_allatoms = 0;
  if(simple_zmat) {
    if(zmat_len >= 1) {
      num_allatoms = 1;
      errcod = ip_string("ZMAT",&buffer,1,0);
      if(errcod != IPE_OK) punt("Problem with the Z-matrix");
      if(strcmp(buffer,"X")) { free(buffer); num_atoms++; }
    }
    if(zmat_len >= 4) {
      num_allatoms = 2;
      errcod = ip_string("ZMAT",&buffer,1,1);
      if(errcod != IPE_OK) punt("Problem with the Z-matrix");
      if(strcmp(buffer,"X")) { free(buffer); num_atoms++; }
    }
    if(zmat_len >= 9) {
      num_allatoms = 3;
      errcod = ip_string("ZMAT",&buffer,1,4);
      if(errcod != IPE_OK) punt("Problem with the Z-matrix");
      if(strcmp(buffer,"X")) { free(buffer); num_atoms++; }
    }
    if(zmat_len > 9) { 
      num_allatoms = 3 + (zmat_len-9)/7;
      if((zmat_len-9)%7 != 0) punt("Error in z-matrix input!");

      for(i=0; i < num_allatoms-3; i++) {
	errcod = ip_string("ZMAT",&buffer,1,9+(i*7));
	if(errcod != IPE_OK) punt("Problem with the Z-matrix");
	if(strcmp(buffer,"X")) { free(buffer); num_atoms++; }
      }
    }
  }
  else {
    ip_count("ZMAT",&num_allatoms,0);
    if (num_allatoms == 0)
      punt("Z-matrix is empty!");
    for(i=0;i<num_allatoms;i++){
      errcod = ip_string("ZMAT",&buffer,2,i,0);
      if (errcod != IPE_OK)
	punt("Problem with the Z-matrix");
      if (strcmp(buffer,"X")) {
	free(buffer);
	num_atoms++;
      }
    }
  }
  if (num_atoms == 0) punt("Z-matrix contains no atoms!");
  if(num_atoms > MAXATOM) punt("Z-matrix contains more atoms than allowed!");

  /*-----------------------
    Allocate global arrays
    -----------------------*/
  full_geom = block_matrix(num_allatoms,3);
  geometry = (double **) malloc(num_atoms*sizeof(double *));
  atom_dummy = (int *) malloc(sizeof(int)*num_allatoms);
  /* see chkpt.h for info about z_entry structure */
  z_geom = (struct z_entry *) malloc(sizeof(struct z_entry)*num_allatoms); 
  element = (char **) malloc(sizeof(char *)*num_atoms);
  full_element = (char **) malloc(sizeof(char *)*num_allatoms);
  elemsymb_charges = init_array(num_atoms);

  atomcount = 0;
  fatomcount = 0;

  /* read in zvars */
  errcod = 0;
  if( ip_exist("ZVARS",0) ) {
    errcod += ip_count("ZVARS",&num_vals,0);

    simple_zvars = 1;
    for(i=0; i < num_vals; i++) {
      ip_count("ZVARS", &entry_length,1,i);
      if(entry_length > 1) simple_zvars = 0;
    }

    if(simple_zvars) {
      /* there must be an even number of entries */
      if(num_vals%2) punt("Problem with number of simple ZVARS entries.");
      num_vals /= 2;
      def_arr = (struct definition *) malloc( num_vals * sizeof(struct definition) );

      for(i=0; i < num_vals; i++) {
	errcod += ip_string("ZVARS", &def_arr[i].variable, 1, 2*i);
	errcod += ip_data("ZVARS", "%lf", &def_arr[i].value, 1, 2*i+1);
      }
    }
    else {
      def_arr = (struct definition *) malloc( num_vals * sizeof(struct definition) );
      for(i=0;i<num_vals;++i) {
        errcod += ip_string("ZVARS",&def_arr[i].variable,2,i,0);
        errcod += ip_data("ZVARS","%lf",&def_arr[i].value,2,i,1);
      }
    }

    if(errcod > 0)
      punt("Problem parsing ZVARS");
  }

/*   if(simple_zvars) exit(1); */

  for (i=0;i<num_allatoms;++i) {

    /* Process a line of ZMAT and write to z_geom array */

    if(!simple_zmat) {
      ip_count("ZMAT",&entry_length,1,i);
      if ( ((i == 0) && (entry_length != 1)) ||
	   ((i == 1) && (entry_length != 3)) ||
	   ((i == 2) && (entry_length != 5)) ||
	   ((i  > 2) && (entry_length != 7)) ) {
	fprintf(outfile,"  Line %d of ZMAT has a wrong number of entries.\n",i+1);
	punt("Invalid ZMAT");
      }
    }

    if (i == 0) {					/*	1st atom */
      full_geom[i][0] = 0.0;
      full_geom[i][1] = 0.0;
      full_geom[i][2] = 0.0;

      z_geom[0].bond_atom = z_geom[0].angle_atom = z_geom[0].tors_atom = -1;
      z_geom[0].bond_opt  = z_geom[0].angle_opt  = z_geom[0].tors_opt  = -1;
      z_geom[0].bond_val  = z_geom[0].angle_val  = z_geom[0].tors_val  = -999.9;
      z_geom[0].bond_label[0] = z_geom[0].angle_label[0] = z_geom[0].tors_label[0] = ' ';
    }

    else if (i == 1) {					/*	2nd atom */

      if(!simple_zmat) {
	ip_data("ZMAT","%d",&a,2,i,1);
	parse_zmat(i,bond,&rAB,def_arr,num_vals,zvar_exist);
      }
      else {
	ip_data("ZMAT", "%d", &a, 1, 2);
	parse_zmat_simple(i,bond,&rAB,def_arr,num_vals,zvar_exist);
      }
	
      z_geom[1].bond_atom = a;
      z_geom[1].angle_atom = z_geom[1].tors_atom = -1;
      z_geom[1].angle_opt = z_geom[1].tors_opt = -1; 
      z_geom[1].bond_val = rAB * conv_factor;
      z_geom[1].angle_val = z_geom[1].tors_val = -999.9;

      if (rAB < ZERO_BOND_DISTANCE)
	punt("Invalid bond length in atom 2.");
      full_geom[i][0] = 0.0;
      full_geom[i][1] = 0.0;
      full_geom[i][2] = rAB;
    }

    else if (i == 2) {					/*	3rd atom */
      if(!simple_zmat) {
	ip_data("ZMAT","%d",&a,2,i,1);
	ip_data("ZMAT","%d",&b,2,i,3);
      }
      else {
	ip_data("ZMAT", "%d", &a, 1, 5);
	ip_data("ZMAT", "%d", &b, 1, 7);
      }
      if ( ((a == 2) && (b == 1)) ||
	   ((a == 1) && (b == 2)) ) {

	if(!simple_zmat) parse_zmat(i,bond,&rBC,def_arr,num_vals,zvar_exist);
	else parse_zmat_simple(i,bond,&rBC,def_arr,num_vals,zvar_exist);
	    
	if (rBC <= ZERO_BOND_DISTANCE) {
	  fprintf(outfile,"  Invalid bond length in atom 3.\n");
	  punt("Invalid ZMAT");
	}

	if(!simple_zmat) parse_zmat(i,angle,&thetaABC,def_arr,num_vals,zvar_exist);
	else parse_zmat_simple(i,angle,&thetaABC,def_arr,num_vals,zvar_exist);
	
	if (thetaABC <= ZERO_BOND_ANGLE) {
	  fprintf(outfile,"  Invalid bond angle in atom 3.\n");
	  punt("Invalid ZMAT");
	}
	z_geom[2].bond_atom = a;
	z_geom[2].bond_val = rBC * conv_factor;
	z_geom[2].angle_atom = b;
	z_geom[2].angle_val = thetaABC;
	z_geom[2].tors_atom = -1;
	z_geom[2].tors_val = -999.9;
	z_geom[2].tors_opt = -1;

	thetaABC = thetaABC*M_PI/180.0;
           
	if (a == 2) {				/*	ABC case */
	  full_geom[i][0] = full_geom[a-1][0] + rBC*sin(thetaABC);
	  full_geom[i][1] = full_geom[a-1][1];
	  full_geom[i][2] = full_geom[a-1][2] - rBC*cos(thetaABC);
	}
	else {					/*	BAC case */
	  full_geom[i][0] = full_geom[a-1][0] + rBC*sin(thetaABC);
	  full_geom[i][1] = full_geom[a-1][1];
	  full_geom[i][2] = full_geom[a-1][2] + rBC*cos(thetaABC);
	}
      }
      else {
	fprintf(outfile,"  Problem in atom 3 in zmat.\n");
	punt("Invalid ZMAT");;
      }
    }

    else { 
      if(!simple_zmat) {
	ip_data("ZMAT","%d",&c,2,i,1);
	ip_data("ZMAT","%d",&b,2,i,3);
	ip_data("ZMAT","%d",&a,2,i,5);
      }
      else {
	ip_data("ZMAT","%d",&c,1,9+((i-3)*7)+1);
	ip_data("ZMAT","%d",&b,1,9+((i-3)*7)+3);
	ip_data("ZMAT","%d",&a,1,9+((i-3)*7)+5);
      }
      a -= 1; b -= 1; c -= 1;

      if ( (a == b) || (b == c) || (a == c) ||
	   (a >= i) || (b >= i) || (c >= i) ) {
	fprintf(outfile,"  Problem in atom %d of zmat.\n",i);
	punt("Invalid ZMAT");
      }

      if(!simple_zmat) parse_zmat(i,bond,&rCD,def_arr,num_vals,zvar_exist);
      else parse_zmat_simple(i,bond,&rCD,def_arr,num_vals,zvar_exist);
	
      if (rCD <= ZERO_BOND_DISTANCE) {
	fprintf(outfile,"  Invalid bond length in atom %d.\n",i+1);
	punt("Invalid ZMAT");
      }

      if(!simple_zmat) parse_zmat(i,angle,&thetaBCD,def_arr,num_vals,zvar_exist);
      else parse_zmat_simple(i,angle,&thetaBCD,def_arr,num_vals,zvar_exist);
	
      if (thetaBCD <= ZERO_BOND_ANGLE) {
	fprintf(outfile,"  Invalid bond angle in atom %d.\n",i+1);
	punt("Invalid ZMAT");
      }

      if(!simple_zmat) parse_zmat(i,tors,&phiABCD,def_arr,num_vals,zvar_exist);
      else parse_zmat_simple(i,tors,&phiABCD,def_arr,num_vals,zvar_exist);
	   
      z_geom[i].bond_atom = c+1;
      z_geom[i].angle_atom = b+1;
      z_geom[i].tors_atom = a+1;
      z_geom[i].bond_val = rCD * conv_factor;
      z_geom[i].angle_val = thetaBCD;
      z_geom[i].tors_val = phiABCD;
      
      thetaBCD = thetaBCD * M_PI/180.0; 
      phiABCD = phiABCD * M_PI/180.0;


      /* If you want to have linear fragment defined in 
	 the beginning of the Z-matrix - fine, but you still have to
	 supply the third "dummy" atom and a value for the dihedral angle*/

      if ( ((full_geom[i-1][0] - LINEAR_CUTOFF) < 0.0) && linearOn ) {
	if ((full_geom[c][2] - full_geom[b][2]) > 0.0) {
	  full_geom[i][0] = rCD*sin(thetaBCD);
	  full_geom[i][1] = 0.0;
	  full_geom[i][2] = full_geom[c][2] - rCD*cos(thetaBCD);
	}
	else {
	  full_geom[i][0] = rCD*sin(thetaBCD);
	  full_geom[i][1] = 0.0;
	  full_geom[i][2] = full_geom[c][2] + rCD*cos(thetaBCD);
	}
      }

      /* Here starts a piece of code for nonlinear ABC fragment */

      else {
	/* Set "linearity" Flag to false so that the current "if" test 
	   is "false" ever since */
	linearOn = 0;
	   
	/* Note : x,y,z - unit vectors defining a coordinate 
	   system associated with atom C; BC defines z, ABC defines 
	   the xz-plane */
	unit_vec(full_geom[b],full_geom[a],eAB);	/* U1 in the old zmat */
	unit_vec(full_geom[c],full_geom[b],ez);	/* U2 in the old zmat */
	cosABC = -dot_prod(ez,eAB);	/* ez is essentially eBC */
	sinABC = sqrt(1 - (cosABC * cosABC) );
	   
	if ( (sinABC - LINEAR_CUTOFF) < 0.0 ) {
	  fprintf(outfile,"  Dihedral angle in line %d is defined with respect to a linear fragment.\n",i+1);
	  punt("Invalid ZMAT");
	}
	   
	cross_prod(eAB,ez,ey);		/* ey is U3 in the old zmat */
	for(j=0;j<3;j++)			/* normalization of ey */
	  ey[j] /= sinABC;
	cross_prod(ey,ez,ex);		/* ex is U4 in the old zmat */

	/* Intermediates - to avoid calling 
	   trig. functions multiple times */
	cosBCD = cos(thetaBCD);
	sinBCD = sin(thetaBCD);
	cosABCD = cos(phiABCD);
	sinABCD = sin(phiABCD);
	   
	for (j=0;j<3;j++)
	  full_geom[i][j] = full_geom[c][j] + rCD * 
	    ( - ez[j] * cosBCD + 
	      ex[j] * sinBCD * cosABCD + 
	      ey[j] * sinBCD * sinABCD );
      }
    }

    if(!simple_zmat) errcod = ip_string("ZMAT",&buffer,2,i,0);
    else {
      if(i==0) errcod = ip_string("ZMAT",&buffer,1,0);
      else if(i==1) errcod = ip_string("ZMAT",&buffer,1,1);
      else if(i==2) errcod = ip_string("ZMAT",&buffer,1,4);
      else if(i>2) errcod = ip_string("ZMAT",&buffer,1,9+((i-3)*7));
    }

    if (strcmp(buffer,"X")) {
      atom_num(buffer, &Z);
      free(buffer);
      elemsymb_charges[atomcount] = Z;
      element[atomcount] = elem_name[(int)Z];
      full_element[fatomcount] = elem_name[(int)Z];
      geometry[atomcount] = full_geom[i];
      atom_dummy[fatomcount] = 0;
      ++atomcount;
    }
    else {
      full_element[fatomcount] = "X";
      free(buffer);
      atom_dummy[fatomcount] = 1;
    }
    ++fatomcount;
  }

  for(i=0;i<num_allatoms;++i) {
    full_geom[i][0] *=conv_factor;
    full_geom[i][1] *=conv_factor;
    full_geom[i][2] *=conv_factor;
  }

  read_charges();

  return;
}



/*******************************************************************************************************************

  parse_zmat

  function used to parse z-matrix input

  position should be 2 for bond, 4 for angle, 6 for torsion (corresponds to postion of value/variable in zmat row)
  since I didn't want to make more stuff global this is pretty ugly, but it cleans things up in read_zmat();
  ******************************************************************************************************************/

void parse_zmat(int i, int position, double *value, struct definition
*array, int num_vals, int zvar_exist) {

   int j, a, value_set;
   double r;
   char *temp_string, *dollar;

   value_set = 0;
   ip_string("ZMAT",&temp_string,2,i,position);
   if( !isalnum( temp_string[0] ) &&
       temp_string[0] != '-' &&
       temp_string[0] != '+' &&
       temp_string[0] != '.')
     punt("z-matrix entry doesn't begin with letter or number");
   dollar = strchr(temp_string,'$');
   if( dollar != NULL ) {
     if(position == 2) 
       z_geom[i].bond_opt = 1;
     else if(position == 4)
       z_geom[i].angle_opt = 1;
     else if(position == 6)
       z_geom[i].tors_opt = 1;
     *dollar = '\0';
   }
   else {
     if(position == 2) 
       z_geom[i].bond_opt = 0;
     else if(position == 4)
       z_geom[i].angle_opt = 0;
     else if(position == 6)
       z_geom[i].tors_opt = 0;
   }
   
   if( isdigit( temp_string[0] ) ||
       ( (temp_string[0] == '-') && isdigit(temp_string[1]) )  ||
       ( (temp_string[0] == '+') && isdigit(temp_string[1]) )  ||
       ( (temp_string[0] == '.') && isdigit(temp_string[1]) ) ) { 
     *value = atof( temp_string );
     ++value_set;
   }
   if( isalpha( temp_string[0] ) && zvar_exist ) {
     for(j=0;j<num_vals;++j) 
       if( !strcmp(temp_string,array[j].variable) ) {
	 *value = array[j].value;
	 ++value_set;
       }
   }
   if(value_set != 1 ) {
     fprintf(outfile,"  Problem with variable definition from line %i\n",i+1);
     punt("Invalid ZMAT");
   }
   
   if( isalpha( temp_string[0] ) ) {
     if(position == 2) 
       strcpy( z_geom[i].bond_label, temp_string );
     if(position == 4)
       strcpy( z_geom[i].angle_label, temp_string );
     if(position == 6)
       strcpy( z_geom[i].tors_label, temp_string );
   }
   else {
     if(position == 2)
       z_geom[i].bond_label[0] = '\0';
     if(position == 4)
       z_geom[i].angle_label[0] = '\0';
     if(position == 6)
       z_geom[i].tors_label[0] = '\0';
   }
   
   free(temp_string);	
   
   
   return;
}

void parse_zmat_simple(int i, int position, double *value, struct definition *array,
		       int num_vals, int zvar_exist)
{
  int j, value_set;
  char *temp_string, *dollar;

  value_set=0;

  if(i==1) ip_string("ZMAT", &temp_string, 1, 3);
  else if(i==2) ip_string("ZMAT", &temp_string, 1, 4+position);
  else ip_string("ZMAT", &temp_string, 1, 9+((i-3)*7)+position);

  if(!isalnum(temp_string[0]) && temp_string[0] != '-' && 
     temp_string[0] != '+' && temp_string[0] != '.')
    punt("Z-matrix entry doesn't begin with a letter or number.");

  dollar = strchr(temp_string, '$');
  if(dollar != NULL) {
    if(position == 2) z_geom[i].bond_opt = 1;
    else if(position == 4) z_geom[i].angle_opt = 1;
    else if(position == 6) z_geom[i].tors_opt = 1;
    *dollar = '\0';
  }
  else {
    if(position == 2) z_geom[i].bond_opt = 0;
    else if(position == 4) z_geom[i].angle_opt = 0;
    else if(position == 6) z_geom[i].tors_opt = 0;
   }

  if ( isdigit(temp_string[0]) || 
       ( (temp_string[0] == '-') && isdigit(temp_string[1]) )  ||
       ( (temp_string[0] == '+') && isdigit(temp_string[1]) )  ||
       ( (temp_string[0] == '.') && isdigit(temp_string[1]) ) ) { 
    *value = atof(temp_string);
    ++value_set;
  }

  if(isalpha(temp_string[0]) && zvar_exist) {
    for(j=0; j < num_vals; j++)
      if(!strcmp(temp_string, array[j].variable)) {
	*value = array[j].value;
	++value_set;
      }
  }
  if(!value_set) {
    fprintf(outfile, "  Problem with variable definition for atom %d.\n", i+1);
    punt("Invalid Z-matrix.");
  }

  if(isalpha(temp_string[0])) {
    if(position == 2) strcpy(z_geom[i].bond_label, temp_string);
    else if(position == 4) strcpy(z_geom[i].angle_label, temp_string);
    else if(position == 6) strcpy(z_geom[i].tors_label, temp_string);
  }
  else {
    if(position == 2) z_geom[i].bond_label[0] = '\0';
    else if(position == 2) z_geom[i].angle_label[0] = '\0';
    else if(position == 6) z_geom[i].tors_label[0] = '\0';
  }

  free(temp_string);

  return;
}

}} // namespace psi::input
