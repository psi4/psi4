/*! \file
    \ingroup R12
    \brief Enter brief description of file here 
*/
/*-------------------------------------------------------------------------
  This is not real DERIV_LVL, it is just a constant that reflects the fact
  that t-integrals require incremented angular momentum (effectively, 
  derivatives of ERIs).
 -------------------------------------------------------------------------*/
#define MAX_AM 16
#define DERIV_LVL 1
#define DEFAULT_NEW_AM 6
#define DEFAULT_OPT_AM 6
#define DEFAULT_MAX_CLASS_SIZE 785

typedef struct {

  /* Twice the maximum AM for which manager routines need to be generated */
  int new_am;

  /* Twice the maximum AM for which workers need to be generated */
  int opt_am;

  /* Twice the AM for which manager routines are already
     generated and stored in an existing library */
  int old_am;

  /* Max number of integrals to be processed in one worker. If a class
   is larger than this then split the worker into several. */ 
  int max_class_size;

} Libr12Params_t;
