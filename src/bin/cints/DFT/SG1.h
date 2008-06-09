#ifndef _psi_src_bin_cints_DFT_SG1_h
#define _psi_src_bin_cints_DFT_SG1_h

/*! \file SG1.h
    \ingroup CINTS
    \brief Enter brief description of file here 

 These are the first Partitioning parameters for the SG-1 Grid
 It is in this format so that I don't have to use if statements */

namespace psi { namespace CINTS {
  double SG1alpha1[] = 
    { 0.25, 0.1667, 0.10};
  
  double SG1alpha2[] = 
    { 0.50, 0.50, 0.40};
  
  double SG1alpha3[] = 
    { 1.0, 0.90, 0.80};
  
  double SG1alpha4[] = 
    { 4.50, 3.50, 2.50};
  
  int SG1angular[]=
    { 6, 38, 86, 194, 86 };
  
  int SG1a2param[]=
    { 0, 0, 0,
      1, 1, 1, 1, 1, 1, 1, 1,
      2, 2, 2, 2, 2, 2, 2, 2,};
};





};
#endif
