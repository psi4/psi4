/*! \file
    \ingroup INPUT
    \brief Enter brief description of file here 
*/
#define EXTERN
#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <cmath>
#include <cstring>
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <libqt/qt.h>
#include <libchkpt/chkpt.h>
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
		int num_vals, int zval_exist, char *zmat_lbl);
void parse_zmat_simple(int i, int position, double *value,struct definition *array, 
		       int num_vals, int zval_exist, char *zmat_lbl);

void read_zmat()
{
  int i, j, k, m, a, b, c, errcod, value_set, zvar_exist;
  char *buffer;
  int num_vals, entry_length, atomcount, all_atomcount;
  int A, B, C, D, f;
  int linearOn = 1;	/* Flag indicating a linear fragment in the beginning of the Z-matrix */
  double rAB, rBC, rCD, thetaABC, thetaBCD, phiABCD, val, norm1, norm2;
  double cosABC, sinABC, cosBCD, sinBCD, cosABCD, sinABCD;
  double eAB[3], eBC[3], ex[3], ey[3], ez[3];
  double Z = 0.0;
  char *name, **zvar_lbl, **zmat_lbl;
  int zmat_len, simple_zmat, simple_zvars;

  /*these are passed to parsing function to tell it what position in the row to look at*/
  int bond=2, angle=4, tors=6;
  struct definition *def_arr;

  frag_num_atoms = (int *) malloc(nfragments*sizeof(int));
  frag_num_allatoms = (int *) malloc(nfragments*sizeof(int));
  frag_atom = (int *) malloc(nfragments*sizeof(int));
  frag_allatom = (int *) malloc(nfragments*sizeof(int));
  nref_per_fragment = (int *) malloc(nfragments*sizeof(int));
  ref_pts_lc = (double ***) malloc(nfragments*sizeof(double **));
  num_atoms = 0;
  num_allatoms = 0;
      
  zvar_lbl = (char **) malloc(nfragments*sizeof(char *));
  zmat_lbl = (char **) malloc(nfragments*sizeof(char *));

  for (f=0; f<nfragments; ++f) { /* loop over fragments */
    zvar_lbl[f] = (char *) malloc(10*sizeof(char)); 
    zmat_lbl[f] = (char *) malloc(10*sizeof(char)); 
    if (f == 0) {
      sprintf(zvar_lbl[f],"ZVARS");
      sprintf(zmat_lbl[f],"ZMAT");
    }
    else { 
      sprintf(zvar_lbl[f],"ZVARS%d",f+1);
      sprintf(zmat_lbl[f],"ZMAT%d",f+1);
    }

    /* must first determine "type" of z-matrix input: nested or simple */
    simple_zmat = 1;
    ip_count(zmat_lbl[f], &zmat_len,0);
    entry_length = 0;
    for(i=0; i < zmat_len; i++) {  
      ip_count(zmat_lbl[f], &entry_length,1,i);
      if(entry_length > 1) {
        simple_zmat = 0;
      }
    }

    /* Read number of lines and count atoms in ZMAT */
    frag_num_atoms[f] = frag_num_allatoms[f] = 0;
    if(simple_zmat) {
      if(zmat_len >= 1) {
        frag_num_allatoms[f] = 1;
        errcod = ip_string(zmat_lbl[f],&buffer,1,0);
        if(errcod != IPE_OK) punt("Problem with the Z-matrix");
        if(strcmp(buffer,"X")) { free(buffer); frag_num_atoms[f]++; }
      }
      if(zmat_len >= 4) {
        frag_num_allatoms[f] = 2;
        errcod = ip_string(zmat_lbl[f],&buffer,1,1);
        if(errcod != IPE_OK) punt("Problem with the Z-matrix");
        if(strcmp(buffer,"X")) { free(buffer); frag_num_atoms[f]++; }
      }
      if(zmat_len >= 9) {
        frag_num_allatoms[f] = 3;
        errcod = ip_string(zmat_lbl[f],&buffer,1,4);
        if(errcod != IPE_OK) punt("Problem with the Z-matrix");
        if(strcmp(buffer,"X")) { free(buffer); frag_num_atoms[f]++; }
      }
      if(zmat_len > 9) { 
        frag_num_allatoms[f] = 3 + (zmat_len-9)/7;
        if((zmat_len-9)%7 != 0) punt("Error in z-matrix input!");
  
        for(i=0; i < frag_num_allatoms[f]-3; i++) {
          errcod = ip_string(zmat_lbl[f],&buffer,1,9+(i*7));
          if(errcod != IPE_OK) punt("Problem with the Z-matrix");
          if(strcmp(buffer,"X")) { free(buffer); frag_num_atoms[f]++; }
        }
      }
    }
    else { /* nested ZMAT */
      ip_count(zmat_lbl[f],&frag_num_allatoms[f],0);
      if (frag_num_allatoms[f] == 0)
        punt("A ZMAT is empty!");
      for(i=0;i<frag_num_allatoms[f];i++){
        errcod = ip_string(zmat_lbl[f],&buffer,2,i,0);
        if (errcod != IPE_OK)
          punt("Problem with the Z-matrix");
        if (strcmp(buffer,"X"))
          frag_num_atoms[f]++;
        free(buffer);
      }
    }
    if (frag_num_atoms[f] == 0) punt("A Z-matrix has no non-dummy atoms!");
    num_atoms += frag_num_atoms[f];
    num_allatoms += frag_num_allatoms[f];
  }

  if(num_atoms > MAXATOM) punt("There are more atoms than allowed!");

  frag_allatom[0] = frag_atom[0] = 0;
  for (f=1; f<nfragments; ++f) {
    frag_atom[f] = frag_atom[f-1] + frag_num_atoms[f-1];
    frag_allatom[f] = frag_allatom[f-1] + frag_num_allatoms[f-1];
  }

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
  all_atomcount = 0;
  for (f=0; f<nfragments; ++f) {

    /*need this to avoid seg fault if character string used but no zvar given*/
    zvar_exist = ip_exist(zvar_lbl[f],0);    
    errcod = 0;
    if (zvar_exist) { /* read in ZVARS */
      errcod += ip_count(zvar_lbl[f],&num_vals,0);
      simple_zvars = 1;
      entry_length = 0;
      for (i=0; i<num_vals; i++) {
        ip_count(zvar_lbl[f], &entry_length,1,i);
        if(entry_length>1) simple_zvars = 0;
      }
      if(simple_zvars) {
        if(num_vals%2) punt("Problem with number of simple ZVARS entries.");
        num_vals /= 2;
        def_arr = (struct definition *) malloc( num_vals * sizeof(struct definition) );
        for(i=0; i < num_vals; i++) {
          errcod += ip_string(zvar_lbl[f], &def_arr[i].variable, 1, 2*i);
          errcod += ip_data(zvar_lbl[f], "%lf", &def_arr[i].value, 1, 2*i+1);
        }
      }
      else { /* nested ZVARS */
        def_arr = (struct definition *) malloc( num_vals * sizeof(struct definition) );
        for(i=0;i<num_vals;++i) {
          errcod += ip_string(zvar_lbl[f],&def_arr[i].variable,2,i,0);
          errcod += ip_data(zvar_lbl[f],"%lf",&def_arr[i].value,2,i,1);
        }
      }
      if(errcod > 0)
        punt("Problem parsing ZVARS");
    }

    for (i=0; i<frag_num_allatoms[f]; ++i) {
      /* Process a line of ZMAT and write to z_geom array */
      if (!simple_zmat) {
        entry_length = 0;
        ip_count(zmat_lbl[f],&entry_length,1,i);
        if ( ((i == 0) && (entry_length != 1)) ||
             ((i == 1) && (entry_length != 3)) ||
             ((i == 2) && (entry_length != 5)) ||
             ((i  > 2) && (entry_length != 7)) ) {
          fprintf(outfile,"  Line %d of ZMAT has a wrong number of entries.\n",i+1);
          punt("Invalid ZMAT");
        }
      }

      if (i == 0) {					/*	1st atom */
        full_geom[all_atomcount][0] = 0.0;
        full_geom[all_atomcount][1] = 0.0;
        full_geom[all_atomcount][2] = 0.0;

z_geom[all_atomcount].bond_atom= z_geom[all_atomcount].angle_atom= z_geom[all_atomcount].tors_atom= -1;
z_geom[all_atomcount].bond_opt = z_geom[all_atomcount].angle_opt = z_geom[all_atomcount].tors_opt = -1;
z_geom[all_atomcount].bond_val = z_geom[all_atomcount].angle_val = z_geom[all_atomcount].tors_val = -999.9;
z_geom[all_atomcount].bond_label[0]= z_geom[all_atomcount].angle_label[0]= z_geom[all_atomcount].tors_label[0]=' ';
      }

      else if (i == 1) {					/*	2nd atom */
        if(!simple_zmat) {
          ip_data(zmat_lbl[f],"%d",&a,2,i,1);
          parse_zmat(i,bond,&rAB,def_arr,num_vals,zvar_exist,zmat_lbl[f]);
        }
        else {
          ip_data(zmat_lbl[f], "%d", &a, 1, 2);
          parse_zmat_simple(i,bond,&rAB,def_arr,num_vals,zvar_exist,zmat_lbl[f]);
        }
	
        z_geom[all_atomcount].bond_atom = a + frag_allatom[f];
        z_geom[all_atomcount].angle_atom= z_geom[all_atomcount].tors_atom = -1;
        z_geom[all_atomcount].angle_opt = z_geom[all_atomcount].tors_opt  = -1; 
        z_geom[all_atomcount].bond_val  = rAB * conv_factor;
        z_geom[all_atomcount].angle_val = z_geom[all_atomcount].tors_val  = -999.9;

        if (rAB < ZERO_BOND_DISTANCE)
          punt("Invalid bond length in atom 2.");
        full_geom[all_atomcount][0] = 0.0;
        full_geom[all_atomcount][1] = 0.0;
        full_geom[all_atomcount][2] = rAB;
      }

      else if (i == 2) {					/*	3rd atom */
        if(!simple_zmat) {
          ip_data(zmat_lbl[f],"%d",&a,2,i,1);
          ip_data(zmat_lbl[f],"%d",&b,2,i,3);
        }
        else {
          ip_data(zmat_lbl[f], "%d", &a, 1, 5);
          ip_data(zmat_lbl[f], "%d", &b, 1, 7);
        }
        if ( ((a == 2) && (b == 1)) ||
             ((a == 1) && (b == 2)) ) {

          if(!simple_zmat) parse_zmat(i,bond,&rBC,def_arr,num_vals,zvar_exist,zmat_lbl[f]);
          else parse_zmat_simple(i,bond,&rBC,def_arr,num_vals,zvar_exist, zmat_lbl[f]);
	    
          if (rBC <= ZERO_BOND_DISTANCE) {
            fprintf(outfile,"  Invalid bond length in atom 3.\n");
            punt("Invalid ZMAT");
          }

          if(!simple_zmat) parse_zmat(i,angle,&thetaABC,def_arr,num_vals,zvar_exist,zmat_lbl[f]);
          else parse_zmat_simple(i,angle,&thetaABC,def_arr,num_vals,zvar_exist,zmat_lbl[f]);
	
          //if (thetaABC <= ZERO_BOND_ANGLE) { // does not permit 0 bond angle with atoms 1-2-3
          if (thetaABC < 0.0) {
            fprintf(outfile,"  Invalid bond angle in atom 3.\n");
            punt("Invalid ZMAT");
          }
          z_geom[all_atomcount].bond_atom = a + frag_allatom[f];
          z_geom[all_atomcount].bond_val  = rBC * conv_factor;
          z_geom[all_atomcount].angle_atom= b + frag_allatom[f];
          z_geom[all_atomcount].angle_val = thetaABC;
          z_geom[all_atomcount].tors_atom = -1;
          z_geom[all_atomcount].tors_val  = -999.9;
          z_geom[all_atomcount].tors_opt  = -1;

          thetaABC = thetaABC*M_PI/180.0;
           
	      if (a == 2) {				/*	ABC case */
            full_geom[all_atomcount][0] = full_geom[a+frag_allatom[f]-1][0] + rBC*sin(thetaABC);
            full_geom[all_atomcount][1] = full_geom[a+frag_allatom[f]-1][1];
            full_geom[all_atomcount][2] = full_geom[a+frag_allatom[f]-1][2] - rBC*cos(thetaABC);
	      }
          else {					/*	BAC case */
            full_geom[all_atomcount][0] = full_geom[a+frag_allatom[f]-1][0] + rBC*sin(thetaABC);
            full_geom[all_atomcount][1] = full_geom[a+frag_allatom[f]-1][1];
            full_geom[all_atomcount][2] = full_geom[a+frag_allatom[f]-1][2] + rBC*cos(thetaABC);
          }
        }
        else {
	      fprintf(outfile,"  Problem in atom 3 in zmat.\n");
	      punt("Invalid ZMAT");;
        }
      }

      else { 
        if(!simple_zmat) {
          ip_data(zmat_lbl[f],"%d",&c,2,i,1);
          ip_data(zmat_lbl[f],"%d",&b,2,i,3);
          ip_data(zmat_lbl[f],"%d",&a,2,i,5);
        }
        else {
          ip_data(zmat_lbl[f],"%d",&c,1,9+((i-3)*7)+1);
          ip_data(zmat_lbl[f],"%d",&b,1,9+((i-3)*7)+3);
          ip_data(zmat_lbl[f],"%d",&a,1,9+((i-3)*7)+5);
        }
        a -= 1; b -= 1; c -= 1;
  
        if ( (a == b) || (b == c) || (a == c) ||
             (a >= i) || (b >= i) || (c >= i) ) {
          fprintf(outfile,"  Problem in atom %d of zmat.\n",i);
          punt("Invalid ZMAT");
        }
  
        if(!simple_zmat) parse_zmat(i,bond,&rCD,def_arr,num_vals,zvar_exist,zmat_lbl[f]);
        else parse_zmat_simple(i,bond,&rCD,def_arr,num_vals,zvar_exist,zmat_lbl[f]);
	  
        if (rCD <= ZERO_BOND_DISTANCE) {
          fprintf(outfile,"  Invalid bond length in atom %d.\n",i+1);
          punt("Invalid ZMAT");
        }
  
        if(!simple_zmat) parse_zmat(i,angle,&thetaBCD,def_arr,num_vals,zvar_exist,zmat_lbl[f]);
        else parse_zmat_simple(i,angle,&thetaBCD,def_arr,num_vals,zvar_exist,zmat_lbl[f]);
	  
        if (thetaBCD <= ZERO_BOND_ANGLE) {
          fprintf(outfile,"  Invalid bond angle in atom %d.\n",i+1);
          punt("Invalid ZMAT");
        }
  
        if(!simple_zmat) parse_zmat(i,tors,&phiABCD,def_arr,num_vals,zvar_exist,zmat_lbl[f]);
        else parse_zmat_simple(i,tors,&phiABCD,def_arr,num_vals,zvar_exist,zmat_lbl[f]);
	     
        z_geom[all_atomcount].bond_atom  = c+1+frag_allatom[f];
        z_geom[all_atomcount].angle_atom = b+1+frag_allatom[f];
        z_geom[all_atomcount].tors_atom  = a+1+frag_allatom[f];
        z_geom[all_atomcount].bond_val   = rCD * conv_factor;
        z_geom[all_atomcount].angle_val  = thetaBCD;
        z_geom[all_atomcount].tors_val   = phiABCD;
        
        thetaBCD = thetaBCD * M_PI/180.0; 
        phiABCD = phiABCD * M_PI/180.0;

        /* If you want to have linear fragment defined in 
	    the beginning of the Z-matrix - fine, but you still have to
	    supply the third "dummy" atom and a value for the dihedral angle*/
        a += frag_allatom[f]; /* switch to full index */
        b += frag_allatom[f];
        c += frag_allatom[f];

        if ( ((full_geom[all_atomcount-1][0] - LINEAR_CUTOFF) < 0.0) && linearOn ) {
          if ((full_geom[c][2] - full_geom[b][2]) > 0.0) {
            full_geom[all_atomcount][0] = rCD*sin(thetaBCD);
            full_geom[all_atomcount][1] = 0.0;
            full_geom[all_atomcount][2] = full_geom[c][2] - rCD*cos(thetaBCD);
          }
          else {
            full_geom[all_atomcount][0] = rCD*sin(thetaBCD);
            full_geom[all_atomcount][1] = 0.0;
            full_geom[all_atomcount][2] = full_geom[c][2] + rCD*cos(thetaBCD);
          }
        }

        else { /* code for nonlinear ABC fragment */
          /* Set "linearity" Flag to false */
          linearOn = 0;
	   
	      /* Note : x,y,z - unit vectors defining a coordinate 
	      system associated with atom C; BC defines z, ABC defines the xz-plane */

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
	        full_geom[all_atomcount][j] = full_geom[c][j] + rCD * 
	          ( - ez[j] * cosBCD + ex[j] * sinBCD * cosABCD + ey[j] * sinBCD * sinABCD );
        }
      }

      if(!simple_zmat) errcod = ip_string(zmat_lbl[f],&buffer,2,i,0);
      else {
        if(i==0) errcod = ip_string(zmat_lbl[f],&buffer,1,0);
        else if(i==1) errcod = ip_string(zmat_lbl[f],&buffer,1,1);
        else if(i==2) errcod = ip_string(zmat_lbl[f],&buffer,1,4);
        else if(i>2) errcod = ip_string(zmat_lbl[f],&buffer,1,9+((i-3)*7));
      }

      if (strcmp(buffer,"X")) {
        atom_num(buffer, &Z);
        free(buffer);
        elemsymb_charges[atomcount] = Z;
        element[atomcount] = elem_name[(int)Z];
        full_element[all_atomcount] = elem_name[(int)Z];
        geometry[atomcount] = full_geom[all_atomcount];
        atom_dummy[all_atomcount] = 0;
        ++atomcount;
      }
      else {
        static const char* dummy_elem_name = "X";
        full_element[all_atomcount] = const_cast<char*>(dummy_elem_name);
        free(buffer);
        atom_dummy[all_atomcount] = 1;
      }
      ++all_atomcount;
    }
  }

  for(i=0;i<num_allatoms;++i) {
    full_geom[i][0] *=conv_factor;
    full_geom[i][1] *=conv_factor;
    full_geom[i][2] *=conv_factor;
  }

  read_charges();

  fprintf(outfile,"Coordinates after reading z-matrices\n");
  print_mat(geometry, num_atoms, 3, outfile);

  return;
}



/*******************************************************************************************************************

  parse_zmat

  function used to parse z-matrix input

  position should be 2 for bond, 4 for angle, 6 for torsion (corresponds to postion of value/variable in zmat row)
  since I didn't want to make more stuff global this is pretty ugly, but it cleans things up in read_zmat();
  ******************************************************************************************************************/

void parse_zmat(int i, int position, double *value, struct definition
*array, int num_vals, int zvar_exist, char *zmat_lbl) {

   int j, a, value_set;
   double r;
   char *temp_string, *dollar;

   value_set = 0;
   ip_string(zmat_lbl,&temp_string,2,i,position);
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
   
   {
     const char* tstring = isalpha(temp_string[0]) ? temp_string : "";
     if(position == 2) strncpy( z_geom[i].bond_label, tstring, CHKPT_ZMAT_LABEL_LEN);
     if(position == 4) strncpy( z_geom[i].angle_label, tstring, CHKPT_ZMAT_LABEL_LEN);
     if(position == 6) strncpy( z_geom[i].tors_label, tstring, CHKPT_ZMAT_LABEL_LEN);
   }
   
   free(temp_string);	
   
   return;
}

void parse_zmat_simple(int i, int position, double *value, struct definition *array,
		       int num_vals, int zvar_exist, char *zmat_lbl)
{
  int j, value_set;
  char *temp_string, *dollar;

  value_set=0;

  if(i==1) ip_string(zmat_lbl, &temp_string, 1, 3);
  else if(i==2) ip_string(zmat_lbl, &temp_string, 1, 4+position);
  else ip_string(zmat_lbl, &temp_string, 1, 9+((i-3)*7)+position);

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

  {
    const char* tstring = isalpha(temp_string[0]) ? temp_string : "";
    if(position == 2) strncpy( z_geom[i].bond_label, tstring, CHKPT_ZMAT_LABEL_LEN);
    if(position == 4) strncpy( z_geom[i].angle_label, tstring, CHKPT_ZMAT_LABEL_LEN);
    if(position == 6) strncpy( z_geom[i].tors_label, tstring, CHKPT_ZMAT_LABEL_LEN);
  }

  free(temp_string);

  return;
}

}} // namespace psi::input
