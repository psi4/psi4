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
#include <psi4-dec.h>
#include <string>

namespace psi { namespace extrema {
int get_coord_type();
void print_intro();
void start_io(int argc, char *argv[]);
void stop_io();


int extrema(int argc, char *argv[]) {
    
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
    return Success;
}



/*-----------------------------------------------------------------------------
  get_coord_type

  print intro and determine coordinate type
  ---------------------------------------------------------------------------*/

int get_coord_type() {

  std::string buffer;

  if(options["COORDINATES"].has_changed()) {
    buffer = options.get_str("COORDINATES");
    if( (buffer == "CARTESIANS") ) 
      punt("Cartesians not available");
    else if( (buffer == "ZMATRIX") ) {
      if( options["ZMAT"].has_changed() ) {
        fprintf(outfile,"\n  Using z-matrix coordinates\n");
        coord_type = ZMAT_TYPE; }
      else
        punt("Can't find z-matrix");
    }
    else if( (buffer == "DELOCALIZED") ) {
      fprintf(outfile,"\n Using delocalized internal coordinates\n");
      coord_type = DELOC_TYPE; }
    else 
      punt("Problem determining coordinate type");
  }    
  else {
    if( options["ZMAT"].has_changed() ) {
      fprintf(outfile,"\n  Defaulting to z-matrix coordinates\n");
      coord_type = ZMAT_TYPE; }
    else if( options["GEOMETRY"].has_changed() ) {
      fprintf(outfile,"\n  Defaulting to delocalized internal coordinates\n");
      coord_type = DELOC_TYPE; }
    else 
      punt("Problem determining coordinate type");
  }

  return coord_type;
}



void print_intro() {

     tstart();
     fprintf(outfile,"                  --------------------------------------------\n");
     fprintf(outfile,"                                   EXTREMA \n");
     fprintf(outfile,"                        Joseph P. Kenny and Rollin King \n");
     fprintf(outfile,"                  --------------------------------------------\n\n");

  return;
}



void start_io(int argc, char *argv[]) {
   
    chkpt_init(PSIO_OPEN_OLD);

    return;
}



void stop_io() {

    chkpt_close();
    tstop();
 
    return;
}

}} // namespace psi::extrema
