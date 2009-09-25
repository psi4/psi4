#ifndef _psi_bin_extrema_deloc_h_
#define _psi_bin_extrema_deloc_h_

namespace psi { namespace extrema {

/*###########################################################################*/
/*! 
  \file
  \ingroup EXTREMA
  \brief Header file for the deloc (delocalized internals) class. */

/*! \class deloc
   \brief Derived class for delocalized internal coordinates.

   A top level class which allows delocalized internal coordinate
   manipulations by driving member functions of the abstract classes. */
/*###########################################################################*/


class deloc : public internals {

  private:

    int	degrees_of_freedom, /*!<molecular degrees of freedom*/
	num_simples, /*!< number of generated simple internals */
	num_nonzero, /*< number of nonzero eigenvalues of G*/
	*deloc_irrep, /*< irrep numbers for coordinates */
	count_array[4], /*!< numbers of each simple internal type*/
	**bonds; /*!< matrix indicating which atoms are bonded */

    double **deloc_define, /*!< delocalized internal coordinate definitions*/ 
	ev_tol; /*!< tolerance for nonzero eigenvalue of G matrix */

    void parse_input();
    void init_simples();
    void set_simples();
    void ir_project();
    int get_bond_index(int at, int bo);
    int get_angle_index(int at, int bo, int an);
    int get_torsion_index(int at, int bo, int an, int to);
    void compute_B_simple();
    void compute_B();
    void compute_block_G();
    void cart_to_internal(double**);
    void initial_Hi();
    void read_opt();
    void newton_step();
    void back_transform();
    void write_opt();
    void read_bonds();
    void print_simples();

    char* point_group;

    simple *simples;

  public:

    deloc();
    ~deloc() {
	free(point_group);
	free(simples);
	free_int_matrix(bonds);
	return; 
    }
    void optimize();

};

/*! \fn deloc::~deloc() 
  \brief Frees memory. */

}} // namespace psi::extrema

#endif // header guard
