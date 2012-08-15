#ifndef _psi_src_bin_psimrcc_ccindex_h
#define _psi_src_bin_psimrcc_ccindex_h

/*! \file
    \ingroup (PSIMRCC)
    \brief   This class is used to store n-tuples of MOs indices (p,q,r,..)
*/
/*********************************************************
  CCIndex Class
  1) Purpose
    This class is used to store n-tuples of MOs indices (p,q,r,..)
  2) Use
  3) Details
    The pairs are stored according to the irreps, for example

          A1           A2          B1       B2
    ----------------------------------------------
    |(0,0),(0,1)| (4,5),(7,9)  | .......|........|
    ----------------------------------------------
  4) Uses
    MOInfo class
    STL <vector>

*********************************************************/

#include <libutil/libutil.h>
#include <vector>

namespace psi{
    extern FILE *outfile;
    namespace psimrcc{

typedef std::vector<std::vector<int> > vecvecint;

/**
	@author Francesco Evangelista <frank@ccc.uga.edu>
*/

class CCIndex{
  typedef std::vector<int*>         pIntVec;
  typedef std::vector<size_t>  Size_tVec;
  typedef std::vector<size_t*> pSize_tVec;
public:
  ///////////////////////////////////////////////////////////////////////////////
  // Class Constructor and Destructor
  ///////////////////////////////////////////////////////////////////////////////
  CCIndex(std::string str);
  ~CCIndex();

  ///////////////////////////////////////////////////////////////////////////////
  // Class Interface
  ///////////////////////////////////////////////////////////////////////////////
  // Print routines
  void        print();
  // Get the numeber of tuples and the numeber of indices per tuple
  int         get_ntuples()                               {return(ntuples);}
  int         get_nelements()                             {return(nelements);}
  std::string get_label()                                 {return(label);}

  // Get the tuples
  short**     get_tuples()                                {return(tuples);}
  short*      get_tuple(int i)                            {return(tuples[i]);}

  // Get the tuples irrep structure
  size_t      get_first(int i)                            {return(first[i]);}
  size_t      get_last(int i)                             {return(last[i]);}
  size_t      get_pairpi(int i)                           {return(tuplespi[i]);}
  size_t      get_tuplespi(int i)                         {return(tuplespi[i]);}
  Size_tVec&  get_pairpi()                                {return(tuplespi);}
  Size_tVec&  get_first()                                 {return(first);}
  Size_tVec&  get_last()                                  {return(last);}
  Size_tVec&  get_tuplespi()                              {return(tuplespi);}


  // Given a set of number retrieve the corresponding tuple index relative to the tuple's irrep
  size_t      get_tuple_abs_index(short p)                    {return(first[one_index_to_irrep[p]] + one_index_to_tuple_rel_index[p]);}
  size_t      get_tuple_abs_index(short p, short q)           {return(first[two_index_to_irrep[p][q]] + two_index_to_tuple_rel_index[p][q]);}
  size_t      get_tuple_abs_index(short p, short q, short r)  {return(first[three_index_to_irrep[p][q][r]] + three_index_to_tuple_rel_index[p][q][r]);}

  // Given a set of number retrieve the corresponding tuple index relative to the tuple's irrep
  size_t      get_tuple_rel_index(short p)                    {return(one_index_to_tuple_rel_index[p]);}
  size_t      get_tuple_rel_index(short p, short q)           {return(two_index_to_tuple_rel_index[p][q]);}
  size_t      get_tuple_rel_index(short p, short q, short r)  {return(three_index_to_tuple_rel_index[p][q][r]);}

  // Given a set of number retrieve the corresponding tuple irrep
  int         get_tuple_irrep(short p)                    {return(one_index_to_irrep[p]);}
  int         get_tuple_irrep(short p, short q)           {return(two_index_to_irrep[p][q]);}
  int         get_tuple_irrep(short p, short q, short r)  {return(three_index_to_irrep[p][q][r]);}

  vecvecint&   get_indices_to_pitzer()                    {return(indices_to_pitzer);}

  size_t*     get_one_index_to_tuple_rel_index()          {return(one_index_to_tuple_rel_index);};
  size_t**    get_two_index_to_tuple_rel_index()          {return(two_index_to_tuple_rel_index);};
  size_t***   get_three_index_to_tuple_rel_index()        {return(three_index_to_tuple_rel_index);};
  int*        get_one_index_to_irrep()                    {return(one_index_to_irrep);};
  int**       get_two_index_to_irrep()                    {return(two_index_to_irrep);};
  int***      get_three_index_to_irrep()                  {return(three_index_to_irrep);};

  int**       get_element_irrep()                         {return(element_irrep);}

private:
  ///////////////////////////////////////////////////////////////////////////////
  // Class private functions
  ///////////////////////////////////////////////////////////////////////////////
  void        init();
  void        cleanup();
  void        make_zero_index();
  void        make_one_index();
  void        make_two_index();
  void        make_three_index();
  ///////////////////////////////////////////////////////////////////////////////
  // Class data
  ///////////////////////////////////////////////////////////////////////////////
  // Type                           // Name
  std::string                       label;                  // The label of this object
  int                               nirreps;                // The number of irreps
  int                               nelements;              // Number of elements
  Size_tVec                         dimension;              // Size of the elements space
  vecvecint                         mospi;                  // Size of each irrep of the elements space
  vecvecint                         first_mos;              // Last mo of each irrep of the elements space
  vecvecint                         indices_to_pitzer;      // Mapping from the mos of this space to Pitzer
  bool                              greater_than_or_equal;  // >= tuples
  bool                              greater_than;           // >  tuples
  size_t                            ntuples;                // Number of tuples
  short**                           tuples;                 // The tuples ordered as matrix : tuples[number][element]
  Size_tVec                         first;                  // First pair of irrep
  Size_tVec                         last;                   // Last  pair of irrep
  Size_tVec                         tuplespi;               // Number of tuples for irrep
  size_t*                           one_index_to_tuple_rel_index;     // 1->tuple mapping array
  size_t**                          two_index_to_tuple_rel_index;     // 2->tuple mapping array
  size_t***                         three_index_to_tuple_rel_index;   // 3->tuple mapping array
  int*                              one_index_to_irrep;     // 1->irrep mapping array
  int**                             two_index_to_irrep;     // 2->irrep mapping array
  int***                            three_index_to_irrep;   // 3->irrep mapping array
  int**                             element_irrep;          // Irrep of each element
};

}} /* End Namespaces */

#endif // _psi_src_bin_psimrcc_ccindex_h
