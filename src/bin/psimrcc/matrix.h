#ifndef _psi_src_bin_psimrcc_matrix_h_
#define _psi_src_bin_psimrcc_matrix_h_
/***************************************************************************
 *  PSIMRCC : Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/

#include <libutil/memory_manager.h>
#include <vector>
#include <string>

namespace psi{ namespace psimrcc{

class CCIndex;

/**
  @author Francesco Evangelista <frank@ccc.uga.edu>
*/
class CCMatrix{
  typedef std::vector<std::pair<int,int> > intpairvec;
  typedef std::vector<double>              DoubleVec;
  typedef std::vector<size_t>              Size_tVec;
  typedef std::vector<bool>                BoolVec;
public:
  ///////////////////////////////////////////////////////////////////////////////
  // Class Constructor and Destructor
  ///////////////////////////////////////////////////////////////////////////////
  CCMatrix(std::string& str,CCIndex* left_index,CCIndex* right_index);
  ~CCMatrix();

  ///////////////////////////////////////////////////////////////////////////////
  // Class Interface
  ///////////////////////////////////////////////////////////////////////////////
  // Functions for scalars
  void         add_scalar(double val);
  void         set_scalar(double val);
  double       get_scalar();

  bool         is_out_of_core();
  bool         is_out_of_core(int h)             const {return(out_of_core[h]);}

  // Functions to get the properties of a matrix
  std::string&  get_label()                             {return(label);}
  std::string&  get_index_label()                       {return(index_label);}
  size_t       get_memory2()                      const {return(memory2);}
  size_t       get_memorypi2(int h)               const {return(memorypi2[h]);}
  int          get_reference()                   const {return(reference);}
  bool         is_integral()                     const {return(integral);}
  bool         is_antisymmetric()                const {return(antisymmetric);}
  bool         is_chemist()                      const {return(chemist_notation);}
  bool         is_fock()                         const {return(fock);}
  int          get_symmetry()                    const {return(symmetry);}

  // Functions to access the indexing and the matrix elements
  CCIndex*     get_left()                        const {return(left);}
  CCIndex*     get_right()                       const {return(right);}
  size_t       get_left_pairpi(int h)            const {return(left_pairpi[h]);}
  size_t       get_right_pairpi(int h)           const {return(right_pairpi[h]);}
  size_t       get_block_sizepi(int h)           const {return(block_sizepi[h]);}
  double**     operator[](int h)                 const {return(matrix[h]);}
  double***    get_matrix()                            {naccess++;return(matrix);}

  // Access the matrix elements
  double       get_two_address_element(short p, short q);
  void         set_two_address_element(short p, short q,double value);
  void         add_two_address_element(short p, short q,double value);
  double       get_four_address_element(short p, short q, short r, short s);
  void         set_four_address_element(short p, short q, short r, short s,double value);
  void         add_four_address_element(short p, short q, short r, short s,double value);

  void         add_six_address_element(short i, short j, short k, short a, short b, short c,double value);
  void         add_six_address_element_abc(short i, short j, short k, size_t abc,double value);
  void         add_six_address_element_ijk(size_t ijk, short a, short b, short c,double value);
  double       get_six_address_element(short i, short j, short k, short a, short b, short c);
  void         add_six_address_element_Pij(short i, short j, short k, short a, short b, short c,double value);
  void         add_six_address_element_Pij_abc(short i, short j, short k, size_t abc,double value);
  void         add_six_address_element_Pik(short i, short j, short k, short a, short b, short c,double value);
  void         add_six_address_element_Pjk(short i, short j, short k, short a, short b, short c,double value);
  void         add_six_address_element_Pjk_abc(short i, short j, short k, size_t abc,double value);
  void         add_six_address_element_Pab(short i, short j, short k, short a, short b, short c,double value);
  void         add_six_address_element_Pab_ijk(size_t ijk, short a, short b, short c,double value);
  void         add_six_address_element_Pbc(short i, short j, short k, short a, short b, short c,double value);
  void         add_six_address_element_Pbc_ijk(size_t ijk, short a, short b, short c,double value);
  void         add_six_address_element_Pij_k(short i, short j, short k, size_t abc, double value);
  void         add_six_address_element_Pijk(short i, short j, short k, short a, short b, short c,double value);
  void         add_six_address_element_Pab_c(size_t ijk, short a, short b, short c, double value);
  void         add_six_address_element_Pij_Pab(short i, short j, short k, short a, short b, short c,double value);
  void         add_six_address_element_Pjk_Pbc(short i, short j, short k, short a, short b, short c,double value);
  void         add_six_address_element_Pij_k_Pa_bc(short i, short j, short k, short a, short b, short c,double value);
  void         add_six_address_element_Pi_jk_Pab_c(short i, short j, short k, short a, short b, short c,double value);
  void         add_six_address_element_Pi_jk_Pa_bc(short i, short j, short k, short a, short b, short c,double value);
  // Access the MO indices of a matrix element
  void         get_two_indices(short*& pq,int irrep, int i, int j);
  void         get_two_indices_pitzer(short*& pq,int irrep, int i, int j);
  void         get_four_indices(short*& pqrs,int irrep, int i, int j);
  void         get_four_indices_pitzer(short*& pqrs,int irrep, int i, int j);

  // Matrix operations
  void         add_numerical_factor(double factor);
  void         add_numerical_factor(double factor, int h);
  void         scale(double factor);
  void         scale(double factor, int h);
  void         zero_matrix();
  void         zero_matrix_block(int h);
  void         zero_two_diagonal();
  void         zero_right_four_diagonal();
  void         zero_left_four_diagonal();
  void         zero_non_doubly_occupied();
  void         zero_non_external();
  void         element_by_element_product(double factor,CCMatrix* B_Matrix,CCMatrix* C_Matrix,int h);
  void         element_by_element_division(double factor,CCMatrix* B_Matrix,CCMatrix* C_Matrix,int h);
  void         element_by_element_addition(double factor,CCMatrix* B_Matrix,int h);
  void         tensor_product(std::string& reindexing,double factor,CCMatrix* B_Matrix,CCMatrix* C_Matrix);
  static double dot_product(CCMatrix* B_Matrix, CCMatrix* C_Matrix, int h);

  // Very Special (VS) Matrix operations



  // Printing
  void         print();
  void         print_dpdmatrix(int n, FILE *out);

  // Memory
  bool         is_allocated();
  bool         is_block_allocated(int h);
  void         allocate_memory();
  void         allocate_block(int h);
  void         free_memory();
  void         free_block(int h);
  int          get_naccess()    {return(naccess);}

  // IO
  void         load();
  void         load_irrep(int h);
  void         dump_to_disk();
  void         dump_to_disk(int first_irrep,int last_irrep);
  void         dump_block_to_disk(int h);
  void         write_block_to_disk(int h);
  void         read_from_disk();
  void         read_from_disk(int first_irrep,int last_irrep);
  void         read_block_from_disk(int h);
  size_t       read_strip_from_disk(int h, int strip, double* buffer);
private:
  ///////////////////////////////////////////////////////////////////////////////
  // Class private functions
  ///////////////////////////////////////////////////////////////////////////////
  std::string  compute_index_label();
  ///////////////////////////////////////////////////////////////////////////////
  // Class data
  ///////////////////////////////////////////////////////////////////////////////
  std::string  label;              // The matrix label
  std::string  index_label;        // The index label
  int          nirreps;            // The number of irreps
  int          reference;          // The reference zeroth-order wavefunction
  double***    matrix;             // Pointer to the allocated memory
                                   // matrix[irrep][left_pair][right_pair]
  CCIndex*     left;               // Pointer to the left indexing scheme
  CCIndex*     right;              // Pointer to the right indexing scheme
  int          symmetry;           // Symmetry of the indices
  size_t*      block_sizepi;       // Size of a subblock of matrix per irrep
  size_t*      left_pairpi;        // Left indexing tuples per irrep
  size_t*      right_pairpi;       // Right indexing tuple per irrep
  bool         integral;           // Is this a two electron integral?
  bool         chemist_notation;   // Is this a two electron integral in chemist notation?
  bool         antisymmetric;      // Is this an antisymmetric two electron integral?
  bool         fock;               // Is this a fock matrix?
  size_t       memory2;             // Memory required for storage in bytes
  Size_tVec    memorypi2;           // Memory required for storage in bytes
  BoolVec      out_of_core;        // Is this irrep stored on disk?
  int          naccess;            // How many times you have called get_matrix();
public:
  static double fraction_of_memory_for_buffer;
};

}} /* End Namespaces */

#endif // _psi_src_bin_psimrcc_matrix_h_
