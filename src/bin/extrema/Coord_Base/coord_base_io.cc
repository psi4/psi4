/*###########################################################################*/
/*! 
  \file
  \ingroup EXTREMA
  \brief File read and write functions common to all coordinate types. */ 
/*						Joseph P. Kenny 11/28/01
  ###########################################################################*/

#define EXTERN
#include"extrema.h"

using namespace psi::extrema;

/*--------------------------------------------------------------------------*/ 
/*! \fn coord_base::parse_input()
  \brief Parses input for info common to all coordinate types. */
/*---------------------------------------------------------------------------*/

void coord_base :: parse_input() {

    char *buffer;

    print_lvl = NORMAL_PRINT+1;
    errcod = ip_data("PRINT","%d",&print_lvl,0); 
    errcod = ip_data("EXTREMA_PRINT","%d",&print_lvl,0);
    fprintf(outfile,"\n  PRINT:         %d",print_lvl);
   
    grad_max = 6;
    errcod = ip_data("GRAD_MAX","%d",&grad_max,0);
    fprintf(outfile,"\n  GRAD_MAX:      %d",grad_max);

    update = "BFGS";
    errcod = ip_string("UPDATE",(char **) &update,0);
    fprintf(outfile,"\n  UPDATE:        %s", update);

    angle_abort = 1;
    errcod = ip_boolean("ANGLE_ABORT",&angle_abort,0);
    fprintf(outfile,"\n  ANGLE_ABORT:   %d", angle_abort);

    // if( ip_exist("SYMMETRY",0) ) {
    //	errcod = ip_string("SYMMETRY",&symmetry,0);
    //	fprintf(outfile,"\n  SYMMETRY:        %s", symmetry);
    //   }
    //else
    //	symmetry = chkpt_rd_sym_label();

    /* do_deriv=1 here only means derivative is possibility
       must read opt.dat to see when it is time to do it for real */
    do_deriv = 0;
    if(ip_exist("DERTYPE",0)) {
	errcod = ip_string("DERTYPE",&buffer,0);
	if( !strncmp(buffer,"FIRST",5) && !strncmp(buffer,"SECOND",6) )
	    do_deriv = 1;
    }
    
    /* see last comment */
    if(ip_exist("OPT",0)) {
	ip_boolean("OPT",&do_opt,0);
    }
         
    return;
}



/*---------------------------------------------------------------------------*/
/*! \fn coord_base::read_opt()
  \brief Reads from opt.dat()
  
  Reads previous coordinates, gradients, and hessian inverse from opt.dat. */
/*---------------------------------------------------------------------------*/

void coord_base :: read_opt() {

	FILE *opt_ptr;

	int i, j, error;

	ffile_noexit(&opt_ptr,"opt.dat",2);
	if( opt_ptr != NULL ) {      
	    
	    ip_done();
	    ip_set_uppercase(1);
	    ip_initialize(opt_ptr,outfile);
	    ip_cwk_add(":OPT_INFO");

	    ip_data("ITERATION","%d",&iteration,0);
	    ++iteration;
      
	    /*read old coordinate vector*/
	    error = 0;
	    error += !ip_exist("COORD",0);  
	    ip_count("COORD",&num_coords,0);
	  
	    for (i=0;i<num_coords;++i) 
	     	error += ip_data("COORD","%lf",&coords_old[i],1,i);
		
	    if(error != 0)
		punt("Problem reading old coordinate values from opt.dat");
  
	    /*read old gradient vector*/
	    error += !ip_exist("GRAD",0);

	    for (i=0;i<num_coords;++i) {
	        error += ip_data("GRAD","%lf",&grads_old[i],1,i);
	    } 

	    if(error != 0)
		punt("Problem reading old gradient from opt.dat");
	    
	    /*read hmat, the inverse of the hessian*/
	    error += (!ip_exist("HINV",0));
	    
	    for (i=0;i<num_coords;++i) {
		for (j=0;j<num_coords;++j) {
		    error += ip_data("HINV","%lf", &Hi_old[i][j], 2, i, j);
		}
	    }
	    
	    fclose(opt_ptr);
      
	    if(error != 0)
		punt("Problem reading old hessian from opt.dat");
        }
        else iteration=1;

	fprintf(outfile,"\n\n  Beginning iteration: %d\n",iteration);

	return;
}



/*---------------------------------------------------------------------------*/
/*! \fn coord_base::write_opt()
  \brief Writes to opt.dat()

  Writes coordinates, gradients, and hessian inverse to opt.dat. */
/*---------------------------------------------------------------------------*/

void coord_base :: write_opt() {

    int place, i, r;
    
    FILE *opt_ptr;

    ip_done();
    ffile(&opt_ptr,"opt.dat",0);
    ip_set_uppercase(1);
    ip_initialize(opt_ptr,outfile);
    
    fprintf(opt_ptr,"opt_info: (\n\n");

    /* write number of coords */
    fprintf(opt_ptr,"  num_coords = %d\n\n",num_coords);    
    
    /*write iteration number*/
    fprintf(opt_ptr,"  iteration = %d\n\n",iteration);
    
    /*write coordinate vector*/
    place = 0;
    fprintf(opt_ptr,"  coord = ( ");
    for (i=0;i<num_coords;++i) {
	if( place==8 ) {
	    place = 0;
	    fprintf(opt_ptr,"\n            ");
	}
	fprintf(opt_ptr,"%lf  ",coord_write[i]);
	++place;
    }
    fprintf(opt_ptr,")\n\n");

    /*write gradient vector*/
    place = 0;
    fprintf(opt_ptr,"  grad = ( ");
    for (i=0;i<num_coords;++i) {
	if( place==8 ) {
	    place = 0;
	    fprintf(opt_ptr,"\n            ");
	}
	fprintf(opt_ptr,"%.20lf  ",grads[i]);
	++place;
    }
    fprintf(opt_ptr,")\n\n");
    
    /*write Hessian one row at a time*/
    place=0;
    fprintf(opt_ptr,"  hinv = ( ");
    for(r=0;r<num_coords;++r) {
	if( place==0 )
	    fprintf(opt_ptr,"\n         ( ");
	for (i=0;i<num_coords;++i) {
	    if( place==8 ) {
		place = 0;
		fprintf(opt_ptr,"\n           ");
	    }
	    fprintf(opt_ptr,"%lf  ",Hi[r][i]);
	    ++place;
	}
	fprintf(opt_ptr,")\n         ");
	place = 0;
    }
    fprintf(opt_ptr,")\n\n");
    
    fprintf(opt_ptr,"          )\n");
    
    fclose(opt_ptr);
    return;
}





