/*! \file 
    \ingroup (IPV1)
    \brief Enter brief description of file here 
*/

/*! \defgroup IPV1 libipv1: The Input Parsing Library */

#include <stdio.h>
#include <stdlib.h>
#include "tmpl.h"
#include "ip_types.h"
#include "ip_global.h"

#include "ip_alloc.gbl"
#include "ip_alloc.lcl"

#include "ip_error.gbl"

extern "C" {

ip_keyword_tree_t *ip_alloc_keyword_tree(void)
{
  ip_keyword_tree_t *result;

  result = (ip_keyword_tree_t *) malloc(sizeof(ip_keyword_tree_t));
  if (!result) {
    perror("ip_alloc_keyword_tree: malloc failed");
    ip_error(NULL);
    }

  result->up = NULL;
  result->down = NULL;
  result->across = NULL;
  result->keyword = NULL;
  result->value = NULL;

  return result;
  }

void ip_free_keyword_tree(ip_keyword_tree_t *tree)
{
  ip_keyword_tree_t *I,*nextI;

  if (!tree) return;

  /* Free all the keyword_trees in the across list. */
  I=tree;
  do {
    /* Free the sub trees first. */
    ip_free_keyword_tree(I->down);
    free(I->keyword);
    /* free any values */
    
    nextI = I->across;
    ip_free_value(I->value);

    /* Zero out I (so if I accidently use it again I'll get SEGV). */
    I->down = NULL;
    I->keyword = NULL;
    I->across = NULL;
    I->up = NULL;
    I->value = NULL;

    free(I);
    } while ((I = nextI) != tree);

  }

ip_value_t *ip_alloc_value(void)
{
  ip_value_t *result;

  result = (ip_value_t *) malloc(sizeof(ip_value_t));
  if (!result) {
    perror("ip_alloc_value: malloc failed");
    ip_error(NULL);
    }
  result->type = IP_UNDEFINED;
  result->v.scalar = NULL;
  return result;
  }

void ip_free_value(ip_value_t *value)
{
  if (!value) return;

  if (value->type == IP_ARRAY) {
    ip_free_array(value->v.array);
    value->v.array = NULL;
    }
  else if (value->type == IP_SCALAR) {
    free(value->v.scalar);
    value->v.scalar = NULL;
    }

  free(value);
  }

void ip_free_array(ip_array_t *array)
{
  int i;

  for (i=0; i<array->n; i++) {
    ip_free_value(array->values[i]);
    }

  free(array->values);
  free(array);
  }

}; /* extern "C" */
