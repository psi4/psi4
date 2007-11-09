/*! \file 
    \ingroup (IPV1)
    \brief Enter brief description of file here 
*/
/* This file provides the routines to do the initial parse of the input
 * file. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tmpl.h"
#include "ip_types.h"
#define _IP_ALLOCATE_GLOBAL_
#include "ip_global.h"
#include "y.tab.h"

#include "ip_read.gbl"
#include "ip_read.lcl"

#include "ip_error.gbl"
#include "ip_print.gbl"
#include "ip_alloc.gbl"
#include "ip_cwk.gbl"
#include "ip_data.gbl"

extern "C" {

/* The input for yacc. */
extern FILE *yyin;
extern void yyparse(void);

/* Set up static variables. */
/*static ip_keyword_tree_t *sub_tree = NULL;*/

/* Set the ip_uppercase global. */
void ip_set_uppercase(int uc)
{
  if (uc) ip_uppercase = 1;
  else ip_uppercase = 0;
  }

/* Initialize the ip routines.  This involves parsing the entire file and
 * converting it into an internal representation. */
/* in = the input file. */
/* out = the output file. */
void ip_initialize(FILE *in, FILE *out)
{

  if (in)  ip_in = yyin = in;
  else     ip_in = yyin = stdin;

  if (out) ip_out = out;
  else     ip_out = stdout;

  /* If ip_tree is not NULL, then ip_initialize has been called twice,
   * with a ip_done inbetween. Call ip_done now. */
  if (ip_tree) {
    ip_warn("ip_initialize has been called twice without an ip_done");
    ip_done();
    }

  yyparse();

  /* The initial cwk list is nothing. */
  ip_cwk_clear();

  ip_internal_values();

}

/* Continue adding to the ip_tree with, presumably, another input file. 
 * This should be called after ip_initialize has been called with different
 * input file.  Multiple calls to ip_append, with different input files,
 * are allowed. */
/* in = the input file. */
/* out = the output file. */
void ip_append(FILE *in, FILE *out)
{

  if (in)  ip_in = yyin = in;
  else     ip_in = yyin = stdin;

  if (out) ip_out = out;
  else     ip_out = stdout;

  /* If ip_tree is NULL, then ip_initialize has probably not been called,
   * but there may have been no data in the first file so we just have
   * to hope the user of these routines knows what is going on. */

/*  08/06/1998 EdV - comment: ip_push_keyword() shifts sub_tree a level 
    down before checking whether the parsed keyword already exists. 
    Therefore sub_tree->down MUST refer to the top of 
    the keyword tree, or ip_tree itself.
    
    if (sub_tree != ip_tree) {
    ip_error("ip_append: sub_tree != ip_tree - impossible");
    }*/

  if (ip_tree) { /* some input has been parsed */
    sub_tree = ip_alloc_keyword_tree();
    sub_tree->down = ip_tree;
  }

  yyparse();

  ip_internal_values();
  }

/* Set up internal ip variables based on the input that has been read in. */
void ip_internal_values(void)
{
  int errcod;

  errcod = ip_boolean(":ip:keyword",&ip_keyword,0);
  if (errcod) ip_keyword = 0;

  }

/* Free all of the data. */
void ip_done(void)
{
  ip_free_keyword_tree(ip_tree);
  ip_tree = NULL;
  sub_tree = NULL;
  ip_in = NULL;
  ip_out = NULL;
  }

void ip_push_keyword(char *keyword)
{
  ip_keyword_tree_t *I, *new_keyword;

  /* If this is the first keyword, then create the tree. */
  if (!ip_tree) {
    sub_tree = ip_tree = ip_alloc_keyword_tree();
    sub_tree->across = sub_tree;
    sub_tree->keyword = keyword;
    return;
    }

  /* This is not the first keyword, so descend the tree. */

  /* If sub_tree is at the top (NULL), then move to ip_tree. */
  if (!sub_tree) {
    sub_tree = ip_tree;
    }
  /* If there is not already a sub_tree->down, then create it. */
  else if (!sub_tree->down) {
    sub_tree->down = ip_alloc_keyword_tree();
    sub_tree->down->across = sub_tree->down;
    sub_tree->down->up = sub_tree;

    sub_tree = sub_tree->down;
    sub_tree->keyword = keyword;
    return;
    }
  /* Descend the tree, but keep track of where we were. */
  else {
    sub_tree = sub_tree->down;
    }

  /* Does the keyword exist in the current sub tree? */
  I=sub_tree;
  do {

    if (!strcmp(I->keyword,keyword)) {
      /* We found it. */
      sub_tree = I;
      free(keyword);
      return;
      }

    } while ((I = I->across) != sub_tree);

  /* We could not find it -- create a new entry. */

  new_keyword = ip_alloc_keyword_tree();
  new_keyword->across = sub_tree->across;
  new_keyword->keyword = keyword;
  sub_tree->across = new_keyword;

  new_keyword->up = sub_tree->up;

  /* Move us down to the new keyword. */
  sub_tree = new_keyword;
  }

void ip_pop_keyword(void)
{
  /* Make sure we aren't already on top. */
  if (!sub_tree) {
    ip_error("ip_pop_keyword: tried to pop above top");
    }
  sub_tree = sub_tree->up;
  }

void ip_assign_value(ip_value_t *value)
{

  /* If sub_tree is still NULL then we have a problem. */
  if (!sub_tree) {
    ip_error("tried to put a keyword at the top level - impossible");
    }

  /* Check for duplicate definitions. */
  if (sub_tree->value) {
#   ifdef DUP_WARN
    /* Warn the user about duplicate definitions. */
    ip_warn("duplicate definition of the following keyword:");
    ip_print_keyword(ip_out,sub_tree);
    fprintf(ip_out,"\n");
    ip_warn("the new value will be ignored");
#   endif /* DUP_WARN */
    ip_free_value(value);
    }
  else sub_tree->value = value;
  }

ip_value_t *ip_scalar(char *scalar)
{
  ip_value_t *result;

  result = ip_alloc_value();
  result->type = IP_SCALAR;
  result->v.scalar = scalar;

  return result;
  }

/* This adds an element to the end of an array, arrayval. */
ip_value_t *ip_array(ip_value_t *arrayval, ip_value_t *newval)
{
  int i;
  ip_value_t *result, **oldvalues;

  /* If arrayval is NULL, then this is to be the first entry in the
   * array. */
  if (!arrayval) {
    result = ip_alloc_value();
    result->type = IP_ARRAY;
    result->v.array = (ip_array_t *) malloc(sizeof(ip_array_t));
    result->v.array->n = 1;
    result->v.array->values = (ip_value_t **) malloc(sizeof(ip_value_t *));
    result->v.array->values[0] = newval;
    return result;
    }

  if (arrayval->type != IP_ARRAY) {
    ip_error("ip_array: tried to extend a nonarray datum");
    }

  result = arrayval;
  oldvalues = result->v.array->values;
  result->v.array->values =
          (ip_value_t **) malloc((result->v.array->n+1)*sizeof(ip_value_t *));
  for (i=0; i<result->v.array->n; i++) {
    result->v.array->values[i] = oldvalues[i];
    }
  free(oldvalues);
  result->v.array->values[result->v.array->n] = newval;
  result->v.array->n++;

  return result;
  }

}; /* extern "C" */
