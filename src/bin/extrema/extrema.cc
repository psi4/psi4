/*###########################################################################*/
/*! \file
    \ingroup EXTREMA
  \brief Main function and related small functions for extrema.

  Provides main function and just enough input parsing to know what
  method drivers main should call. */

  /*! \fn void main()
    \brief Main function for extrmema.
    Initializes top-level class objects and calls method driver functions. */
/*						Joseph P. Kenny 11/29/01
  ##########################################################################*/

#include "extrema.h"

namespace psi { namespace extrema {
int get_coord_type();
void print_intro();
void start_io(int argc, char *argv[]);
void stop_io();
}}


int main(int argc, char *argv[]) {
    using namespace psi::extrema;
    
    start_io(argc,argv);
    print_intro();
    coord_type = get_coord_type();

    /*------------
      CARTESIANS
      ----------*/
    if(coord_type==1) {
	//carts c_obj;
    }
  

    /*---------
      ZMATRIX
      -------*/
    else if(coord_type==2) {
	zmat z_obj;
	z_obj.optimize();
    }


    /*-------------
      DELOCALIZED
      -----------*/
    else if(coord_type==3) {
	deloc d_obj;
	d_obj.optimize();
    }
    
    stop_io();
    return 0;
}


namespace psi { namespace extrema {

/*-----------------------------------------------------------------------------
  get_coord_type

  print intro and determine coordinate type
  ---------------------------------------------------------------------------*/

int get_coord_type() {

  char *buffer;
  
  if(ip_exist("COORDINATES",0)) {
      errcod = ip_string("COORDINATES", &buffer,0);
      if( !strcmp(buffer,"CARTESIANS") ) 
	  punt("Cartesians not available");
      else if( !strcmp(buffer,"ZMATRIX") ) {
	  if( ip_exist("ZMAT",0) ) {
	      fprintf(outfile,"\n  Using z-matrix coordinates\n");
	      coord_type = ZMAT_TYPE; }
	  else
	      punt("Can't find z-matrix");
      }
      else if( !strcmp(buffer,"DELOCALIZED") ) {
	  fprintf(outfile,"\n Using delocalized internal coordinates\n");
	  coord_type = DELOC_TYPE; }
      else 
	  punt("Problem determining coordinate type");
      free(buffer);
  }    
  else {
      if( ip_exist("ZMAT",0) ) {
        fprintf(outfile,"\n  Defaulting to z-matrix coordinates\n");
        coord_type = ZMAT_TYPE; }
      else if( ip_exist("GEOMETRY",0) ) {
        fprintf(outfile,"\n  Defaulting to delocalized internal coordinates\n");
        coord_type = DELOC_TYPE; }
      else 
	punt("Problem determining coordinate type");
  }

  return coord_type;
}



void print_intro() {

     tstart(outfile);
     fprintf(outfile,"                  --------------------------------------------\n");
     fprintf(outfile,"                                   EXTREMA \n");
     fprintf(outfile,"                        Joseph P. Kenny and Rollin King \n");
     fprintf(outfile,"                  --------------------------------------------\n\n");

  return;
}



void start_io(int argc, char *argv[]) {
   
    psi_start(&infile,&outfile,&psi_file_prefix,argc-1,argv+1,0); 
    ip_cwk_add(":INPUT");
    ip_cwk_add(":EXTREMA");
    psio_init(); psio_ipv1_config();
    chkpt_init(PSIO_OPEN_OLD);

    return;
}



void stop_io() {

    chkpt_close();
    tstop(outfile);
    psi_stop(infile,outfile,psi_file_prefix);
 
    return;
}

}} // namespace psi::extrema
