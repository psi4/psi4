/*###########################################################################*/
/*! 
  \file
  \ingroup EXTREMA
  \brief Member function definitions. */
/*						Joseph P. Kenny 12/07/01
  ###########################################################################*/

#define EXTERN
#include "extrema.h"
#include "masses.h"


using namespace psi::extrema;

/*---------------------------------------------------------------------------*/
/*! \fn coord_base_carts::coord_base_carts() 
  \brief Reads cartesian info.
/*---------------------------------------------------------------------------*/

coord_base_carts :: coord_base_carts() {

    int i, num, p;
    char **temp_names;

    /* if dummy atoms are present derived classes will reset num_entries */
    num_atoms = num_entries = chkpt_rd_natom();
 
    carts = init_array(3*num_atoms);
    c_grads = init_array(3*num_atoms);
    masses = init_array(num_atoms);
    e_names = (char**) malloc(num_atoms * sizeof(char*) );

    temp_names = (char **) chkpt_rd_felement();
    num = chkpt_rd_nallatom();
    p=-1;
    for(i=0;i<num;++i) {
	if(strncmp(temp_names[i],"X\0",2))
	   e_names[++p] = temp_names[i];
	else
	   free( temp_names[i] );
    }

    double** t_carts;
    t_carts = chkpt_rd_geom();
    p=-1;
    for(i=0;i<(3*num_atoms);++i) {
	++p;
	if(fabs(t_carts[0][p])>1.0e-15)
	    carts[i] = t_carts[0][p];
    }
    free(t_carts);

    coord_base_carts::read_file11();

    return;
}



/*---------------------------------------------------------------------------*/
/*! \fn coord_base_carts::read_file11()
  \brief Reads cartesian gradients from file11. */
/*---------------------------------------------------------------------------*/

void coord_base_carts :: read_file11() {

    int i, natom, count = 1, continue_flag = 1;
    char label[133]; char line1[133]; char *tmp_ptr;
    FILE *fp_11;
    double an,x,y,z,energy;
    
    ffile(&fp_11,"file11.dat",2);
    if (fp_11 == NULL) {
	punt("Could not open file11.dat");
    }
    
    tmp_ptr = fgets(label, MAX_LINELENGTH, fp_11);
    
    if (tmp_ptr == NULL) {
	punt("Touble reading first line of file11.dat");
    }
    
    fgets(line1, MAX_LINELENGTH, fp_11);
    if (sscanf(line1, "%d %lf", &natom, &energy) != 2) {
	punt("Trouble reading natoms and energy from file11.dat");
    }
    
    if(natom!=num_atoms)
	punt("Numbers of atoms differ in file11 and chkpt");
    
    rewind(fp_11);
    
    while ( fgets(label, MAX_LINELENGTH, fp_11) != NULL ) {
	
	fgets(line1, MAX_LINELENGTH, fp_11);
	sscanf(line1, "%d %lf", &natom, &energy);
	
	/*read in one chunk at a time*/ 
	for (i=0; i<num_atoms ; i++) {
	    if(fscanf(fp_11, "%lf %lf %lf %lf", &an, &x, &y, &z) != 4) {
		punt("Trouble reading cartesian coordinates from file11.dat");
	    }
	    masses[i]=an2masses[(int) an];
	}
	
	for (i=0; i<num_atoms ; i++) {
	    if(fscanf(fp_11, "%lf %lf %lf", &x, &y, &z) != 3) {
		punt("Trouble reading gradients from file11.dat");
	    }
	    c_grads[3*i] = x;
	    c_grads[3*i+1] = y;
	    c_grads[3*i+2] = z;
	}
	
	fgets(line1, MAX_LINELENGTH, fp_11);
	
    }

    fclose(fp_11);
    return;
}



/*---------------------------------------------------------------------------*/
/*! \fn coord_base_carts::write_chkpt()
  \brief Writes cartesians to file 30. */
/*---------------------------------------------------------------------------*/

void coord_base_carts::write_chkpt() {
    
    int i,j, pos=-1;
    double** cart_matrix;
    cart_matrix = block_matrix(num_atoms,3);    

    for(i=0;i<num_atoms;++i) 
	for(j=0;j<3;++j) 
	    cart_matrix[i][j] = carts[++pos];

    chkpt_wt_geom(cart_matrix);

    return;
}



/*---------------------------------------------------------------------------*/
/*! \fn coord_base_carts::print_carts(double)
  \brief Prints cartesians.

  \param conv conversion factor;
  either <b>1.0</b> for bohr or <b>_bohr2angstroms</b> for angstroms */
/*---------------------------------------------------------------------------*/

void coord_base_carts :: print_carts(double conv) {

    int i, j;
    double **temp;

    temp = init_matrix(num_entries,3);
    for(i=0;i<num_entries;++i) 
	for(j=0;j<3;++j) 
	    temp[i][j] = carts[3*i+j]*conv;
    if(conv==1.0)
	fprintf(outfile,"\n  Cartesian Coordinates (bohr):\n");
    else
	fprintf(outfile,"\n  Cartesian Coordinates (angstroms):\n");
     fprintf(outfile,"                       x              y         ");
    fprintf(outfile,"       z\n");
    fprintf(outfile,"                --------------- --------------- ");
    fprintf(outfile,"---------------\n");
    for(i=0;i<num_entries;++i)
	fprintf(outfile,"  %12s  %15.10lf %15.10lf %15.10lf\n",
		e_names[i], temp[i][0], temp[i][1], temp[i][2]);
    free_matrix(temp,num_entries);

    return;
}



/*--------------------------------------------------------------------------*/
/*! \fn coord_base_carts::print_c_grads()
  \brief Prints cartesian gradients (Hartree/Bohr). */
/*--------------------------------------------------------------------------*/

void coord_base_carts :: print_c_grads() {
	
    int i, j;

    double **cgtemp;
    cgtemp = init_matrix(num_entries,3);
	
    for(i=0;i<num_entries;++i) 
	for(j=0;j<3;++j) 
	    cgtemp[i][j] = c_grads[3*i+j];
    fprintf(outfile,"\n  Cartesian Gradients (a.u):\n");
    fprintf(outfile,"                       x              y         ");
    fprintf(outfile,"       z\n");
    fprintf(outfile,"                --------------- --------------- ");
    fprintf(outfile,"---------------\n");
    for(i=0;i<num_entries;++i)
	fprintf(outfile,"  %12s  %15.10lf %15.10lf %15.10lf\n",
		e_names[i], cgtemp[i][0], cgtemp[i][1], cgtemp[i][2]);
    free_matrix(cgtemp,num_entries);
    
    return;
}







