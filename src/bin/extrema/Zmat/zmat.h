#ifndef _psi_bin_extrema_zmat_h_
#define _psi_bin_extrema_zmat_h_

namespace psi { namespace extrema {

/*###########################################################################*/
/*! 
  \file
  \ingroup EXTREMA
  \brief Header file for the zmat class. */

/*! \class zmat
   \brief Derived class for z-matrix coordinates.

   A top level class which allows z-matrix
   manipulations by driving member functions of the abstract classes.
   Currently optimization is the only available method. */

/*###########################################################################*/


class zmat : protected internals {

  protected:

    struct z_entry *z_geom;
    /*!< an array of z_entry structs, defined in 
	chkpt.h, a z-matrix is read/written to chkpt as an array of z_entry 
	structs */
    
    simple *simples;
    /*!< an array of simple classes, the simple concrete
       class provides interfaces to simple internal coordinate information */ 

    int *first_unique;
    /*!< 1 if a simple intenal from the full set 
      of simples is the first unique coordinate, 0 otherwise, 
      only unique coordinates need to be optimized */

    int *pos_neg_pairs;
    /*!< positive negative torsion angle pairs are 
      indicated by matching positive integers, for keeping torsion
      angle pairs matching */

    double *fcoord_old;
    /*!< full coordinates prior to optimization step, needed for
      back transformation to cartesians */

    double bond_lim; /*!< maximum bond displacement in bohr */
    double angle_lim; /*!< maximum angle displacement in radians */


  void compute_B(void);
  void cart_to_internal(double**);
  void print_internals(void);
  void initial_Hi(void);
  void grad_trans();
  void parse_input();
  void write_chkpt();
  void newton_step();
  void back_transform();
  void print_carts(double);
  void print_c_grads();

  public:

    zmat();
    ~zmat() {
        int i; 
	free(z_geom); 
	free(simples);
	free(first_unique);
	free(pos_neg_pairs);
	return; }
    void optimize();

};

/*! \fn zmat::~zmat() 
  \brief zmat destructor 

  Frees allocated memory.*/

}} // namespace psi::extrema

#endif // header guard
