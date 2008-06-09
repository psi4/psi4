/*! \file
    \ingroup IPV1
    \brief Enter brief description of file here 
*/

#ifndef _psi_src_lib_libipv1_iptypes_h_
#define _psi_src_lib_libipv1_iptypes_h_

#ifdef __cplusplus
extern "C" {
#endif

#define IP_UNDEFINED  0
#define IP_ARRAY      1
#define IP_SCALAR     2

/* This is a tree with all of the keywords and values. */
struct ip_keyword_tree_struct {
  char *keyword;
  struct ip_keyword_tree_struct *across; /* Circular list. */
  struct ip_keyword_tree_struct *up;    /* Terminated by NULL. */
  struct ip_keyword_tree_struct *down;  /* Terminated by NULL. */
  struct ip_value_struct *value;
  };

struct ip_keyword_tree_list_struct {
  struct ip_keyword_tree_struct *kt;
  struct ip_keyword_tree_list_struct *p;
  };

struct ip_value_struct {
  int type;
  union {
    struct ip_array_struct *array;
    char *scalar;
    } v;
  };

struct ip_array_struct {
  int n;
  struct ip_value_struct **values;
  };

typedef struct ip_value_struct ip_value_t;
typedef struct ip_array_struct ip_array_t;
typedef struct ip_keyword_tree_struct ip_keyword_tree_t;
typedef struct ip_keyword_tree_list_struct ip_keyword_tree_list_t;

#ifdef __cplusplus
}
#endif

#endif /* header guard */
