#ifndef _psi_bin_extrema_internals_h_
#define _psi_bin_extrema_internals_h_

namespace psi { namespace extrema {

/*###########################################################################*/
/*! 
  \file
  \ingroup EXTREMA
  \brief Derived class <b>internals</b> declaration. */

/*! \class internals
  \brief Abstract internal coordinate derived class.

  Member functions and data common to internal coordinate types.  Primarily,
  this class contains the B matrix and related matrices, functions for
  computing these matrices, and functions which apply these matrices in
  gradient and coordinate transformations 

  "full" variables refer to the complete set of coordinates, including
  any equivalent coordinates.  In the iterative back transformation to 
  cartesians, the complete set of coordinates is utilized.
*/

/*						Joseph P. Kenny 11/29/01
  ###########################################################################*/

class internals : protected coord_base {

  protected:

    int fnum_coords, /*!< full number of coordinates */
	B_dim; /*!< fnum_coords for zmat, num_coords for deloc */
    double **full_geom, /*!< full geometry (including dummy atoms) */
	*fcoords, /*!< full set of internal coordinate values */
	*fcoords_old, /*!< full set of previous internal coordinate values */
	*fgrads, /*!< full set of internal coordinate gradients */
	**B,/*!< the B matrix */
	**B_red, /*!< reduced dimension B matrix (no equivalent coords) */
	**G, /*!< the G=B.B^t matrix */
	**A, /*!< the A=u.B^t.G matrix */
	**u; /*!< the u matrix (triplets of inverse atomic masses) */

    internals() { return; };
    void mem_alloc();
    ~internals() {
	free_matrix(full_geom,num_entries);
	free_matrix(B_red,num_coords);
        free_matrix(B,fnum_coords);
	free_matrix(G,fnum_coords);
	free_matrix(A,3*num_entries);
	free_matrix(u,3*num_entries);
	return;
    }
    virtual void compute_B() = 0;
    virtual void cart_to_internal(double**) = 0;
    void back_transform(double*, double*);
    double *B_row_bond(double*, int, int);
    double *B_row_angle(double*, int, int, int);
    double *B_row_tors(double*, int, int, int, int);
    void print_B();
    void compute_G();
    void compute_A();
    virtual void grad_trans();
};

/*! \fn internals::internals()
  \brief Dummy constructor.  
  
  This is a dummy constructor which exists to make compilers happy.  It does
  nothing.  The actual constructor may depend on data which is not know
  when a derived class is initialized and may be called later in that
  class's constructor. */

/*! \fn internals::~internals()
  \brief Destructor.  Frees memory. */

/*! \fn internals::compute_B()
  \brief Derived classes must provide a function to compute the B matrix 
  though calls to <b>B_row...()</b> functions. */

/*! \fn internals::cart_to_internal(double *cg)
  \brief Derived classes must provide a function to transform cartesian
  coordinates into internal coordinates.  This is used by  
  <b>back_transform()</b>. 
  \param double *cg the cartesian geometry for which internals are desired. */

}} // namespace psi::extrema

#endif // header guard
