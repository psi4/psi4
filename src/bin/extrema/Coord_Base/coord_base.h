#ifndef _psi_bin_extrema_coordbase_h_
#define _psi_bin_extrema_coordbase_h_

namespace psi { namespace extrema {

/*###########################################################################*/
/*! 
  \file
  \ingroup EXTREMA
  \brief Coordinate base class declarations 
  (coord_base and coord_base_carts).*/

/*! \class coord_base_carts
  \brief Cartesian information needed for all coordinate representations.

  All coordinate representations are formed from the basic cartesian
  coordinates and gradients read from chkpt and file11.  This basic
  information is read and stored in this class. */
/*  						Joseph P. Kenny 12/08/01
  ###########################################################################*/

class coord_base_carts {
  
  protected:

    /* these should belong to coord_base, but are needed 
       to read and print cartesians */
    int num_atoms, /*!< number of atoms */ 
	num_entries; /*!< number of atoms + dummy atoms */


    double *carts, /*!< cartesian coordinate array */
	*c_grads, /*!< cartesian gradient array */
	*masses; /*!< masses array, should be member of coord_base class,
		   but it's read in with cartesian info */

    char *symmetry; /*!< symmetry from input.dat or chkpt */
    char **e_names;  /*!< element names (no dummy atoms) */
	
    virtual void print_carts(double conv);
    virtual void print_c_grads();
    void read_file11();
    virtual void write_chkpt();   

  public:

    coord_base_carts();
    ~coord_base_carts(){
	int i;
	free(carts);
	free(c_grads);
	free(masses);
	for(i=0;i<num_atoms;++i)
	    free(e_names[i]);
	free(e_names);
	return;
    }
};



/*###########################################################################*/
/*! \class coord_base
  \brief First level of abstract coordinate classes.

  The <b>coord_base</b> class contains data and functions common to 
  all coordinate types.  All coordinate types derive from this class.
  Member data includes generic coordinates and gradients, cartesians and 
  cartesian gradients, a generic Hessian, information form previous iterations
  ("old" variables), and basic user suppied parameters.  Generic functions
  for coordinate data manipulations are members of this class.

  "generic" variables hold values and no information regarding
  coordinates to which they correspond.  Classes deriving from this 
  class determine the actual coordinate type and are responsible for 
  proper handling of these variables.  */
/*  						Joseph P. Kenny 11/29/01
  ###########################################################################*/

class coord_base : protected coord_base_carts,  protected math_tools {
  
  protected:

    int iteration, /*!< current iteration */ 
	num_coords, /*!< number of coordinates which are actually optimized */
	grad_max, /*!< max allowable gradient is 10^-(grad_max) */
	print_lvl, /*!< print level */
    	do_deriv, /*!< are we doing derivatives */
	do_opt,   /*!< are we doing optimization */
	angle_abort; /*!< die if bad angle */

    double *coords, /*!< generic coordinate array */ 
	*grads, /*!< generic gradient array */
	*atomic_nums, /*!< atomic numbers */
	**Hi, /*!< generic inverse hessian matrix */
	*coords_old, /*!< generic coordinates from previous iteration */
	*grads_old, /*!< generic gradients from previous iteration */
	**Hi_old, /*!< generic hessian inverse */
	*coord_write; /*!< holds coordinate values prior to optimization step
			until opt.dat is written */

    /*! \note "generic" variables hold values and no information regarding
       coordinates to which they correspond.  Classes deriving from this 
       class determine the actual coordinate type and are responsible for 
       proper handling of these variables */ 

    const char *update; /*!< the hessian inverse update method */
    char **felement; /*!< the full list of element names 
		      (including dummy atoms) */

    coord_base();
    ~coord_base(){

	free(coords);
	free(grads);
	free_matrix(Hi, num_coords);
	free(coords_old);
	free(grads_old);
	free_matrix(Hi_old,num_coords);
	free(coord_write);
	free(atomic_nums);
	int i;
        for(i=0;i<num_entries;++i)
	    free(felement[i]);

	return;
    }
    void mem_alloc();
    void parse_input();
    void print_Hi();
    void print_H();
    virtual void read_opt();
    virtual void write_opt();
    void update_Hi();
    void grad_test();
    virtual void initial_Hi()=0;
    void H_test();
};

/*! \fn coord_base::~coord_base()
  \brief Destructor.  Memory is freed and io is stopped */

/*! \fn coord_base::initial_Hi()
  \brief Derived classes must provide a inverse Hessian guess function. */

/*! \fn coord_base::write_chkpt()
  \brief Derived classes may override this function and write necessary 
  information to chkpt */

}} // namespace psi::extrema

#endif // header guard
