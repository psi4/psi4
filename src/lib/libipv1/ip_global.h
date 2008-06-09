/*! \file
    \ingroup IPV1
    \brief Enter brief description of file here 
*/

#ifndef _psi_src_lib_libipv1_ipglobal_h_
#define _psi_src_lib_libipv1_ipglobal_h_

#ifdef __cplusplus
extern "C" {
#endif

 /* Global data for the ip routines. */

#ifdef _IP_ALLOCATE_GLOBAL_
#define EXTERN
#define INITIALIZE(x,y) x=y
#else
#define EXTERN extern
#define INITIALIZE(x,y) x
#endif

  /*  The input file. */
EXTERN FILE INITIALIZE(*ip_in,NULL);

  /*  The output file. */
EXTERN FILE INITIALIZE(*ip_out,NULL);

  /* The parsed input. */
EXTERN ip_keyword_tree_t INITIALIZE(*ip_tree,NULL);
EXTERN ip_keyword_tree_t INITIALIZE(*sub_tree,NULL);

  /* Pointers within ip_tree to the current keywords
   * (up and down are unused). */
EXTERN ip_keyword_tree_list_t INITIALIZE(*ip_cwk,NULL);

 /* This is 1 if everything is to be converted into uppercase. */
EXTERN int INITIALIZE(ip_uppercase,0);

 /* This is 1 if keywords are to be shown as they are searched. */
EXTERN int ip_keyword;

#ifdef __cplusplus
}
#endif

#endif /* header guard */
