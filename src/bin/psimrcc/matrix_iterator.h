// #ifndef _psi_src_bin_psimrcc_ccmatrix_iterator_h
// #define _psi_src_bin_psimrcc_ccmatrix_iterator_h
// 
// /*! \file    index_iterator.h
//     \ingroup (PSIMRCC)
//     \brief   This class is used to iterate over n-tuples of MOs indices (p,q,r,..)
// */
// 
// #include <string>
// 
// #include <index.h>
// 
// namespace psi{ namespace psimrcc{
// 
// class CCMatrixIterator{
// public:
//   ///////////////////////////////////////////////////////////////////////////////
//   // Class Constructor and Destructor
//   ///////////////////////////////////////////////////////////////////////////////
//   CCMatrixIterator(string left_str,string right_str,int irrep);
//   ~CCMatrixIterator();
//   ///////////////////////////////////////////////////////////////////////////////
//   // Class Methods
//   ///////////////////////////////////////////////////////////////////////////////
// 
//   bool operator++()
//   {
//     if(rel < max_index_in_irrep_ - 1){
//       // Increment index
//       ++rel;
//       ++abs;
//       if(nelements_>0){
//         ind[0] = tuples[abs][0];
//         ind_sym[0] = element_irrep[0][ind[0]];
//       }
//       if(nelements_>1){
//         ind[1] = tuples[abs][1];
//         ind_sym[1] = element_irrep[1][ind[1]];
//       }
//       if(nelements_>2){
//         ind[2] = tuples[abs][2];
//         ind_sym[2] = element_irrep[2][ind[2]];
//       }
//       return(true);
//     }else{
//       // Check if it reached the last non empty irrep
//       sym = next_non_empty_irrep(sym);
//       if(sym == -1)
//         return(false);
//       else{
//         max_index_in_irrep_ = ccindex_->get_pairpi(sym);
//         rel = 0;
//         ++abs;
//         if(nelements_>0){
//           ind[0] = tuples[abs][0];
//           ind_sym[0] = element_irrep[0][ind[0]];
//         }
//         if(nelements_>1){
//           ind[1] = tuples[abs][1];
//           ind_sym[1] = element_irrep[1][ind[1]];
//         }
//         if(nelements_>2){
//           ind[2] = tuples[abs][2];
//           ind_sym[2] = element_irrep[2][ind[2]];
//         }
//         return(true);
//       }
//     }
//   }
//   void        reset();
// 
//   ///////////////////////////////////////////////////////////////////////////////
//   // Class Public Members
//   ///////////////////////////////////////////////////////////////////////////////
//   int                               sym;
//   size_t                            rel;
//   size_t                            abs;
//   short*                            ind;
//   int*                              ind_sym;
// private:
//   ///////////////////////////////////////////////////////////////////////////////
//   // Class private functions
//   ///////////////////////////////////////////////////////////////////////////////
//   void        init();
//   void        cleanup();
//   int         next_non_empty_irrep(int n);
//   ///////////////////////////////////////////////////////////////////////////////
//   // Class data
//   ///////////////////////////////////////////////////////////////////////////////
//   // Type                           // Name
//   CCIndex*                          ccindex_left;
//   CCIndex*                          ccindex_right;
//   int                               nelements_;
//   int**                             element_irrep;
//   short**                           tuples;
//   size_t                            max_index_in_irrep_;
// 
// protected:
//   static int                        nirreps_;
// };
// 
// }} /* End Namespaces */
// 
// #endif // _psi_src_bin_psimrcc_ccindex_iterator_h
