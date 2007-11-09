/*! \file 
    \ingroup (DERIV)
    \brief Enter brief description of file here 
*/
#define MAX_AM 16
#define DERIV_LVL 2
#define DEFAULT_NEW_AM1 6  /*--- first derivatives of up to f-functions ---*/
#define DEFAULT_NEW_AM2 4  /*--- first derivatives of up to d-functions ---*/
#define DEFAULT_NEW_AM12 4  /*--- second derivatives of up to d-functions ---*/
#define EMIT_DERIV2_MANAGERS 0  /*--- whether to produce manager functions that compute
				   second derivatives only ---*/


typedef struct {

  /* Twice the maximum AM for which deriv1/deriv2/deriv12 manager routines need to be generated */
  int new_am;

  /* Twice the maximum AM for which deriv1 manager routines need to be generated */
  int new_am1;

#if EMIT_DERIV2_MANAGERS
  /* Twice the maximum AM for which deriv2 manager routines need to be generated */
  int new_am2;
#endif

  /* Twice the maximum AM for which deriv12 manager routines need to be generated */
  int new_am12;

  /* The maximum AM + 1 for which machine-generated VRR workers are present in libint.a */
  int opt_am;

  /* Twice the AM for which manager routines are already
     generated and stored in an existing library */
  int old_am;

  /* The maximum AM for which VRR workers will be made inline functions
     Determined by libint.a */
  int max_am_to_inline_vrr_worker;

  /* The maximum AM for which DERIV workers will be made inline functions */
  int max_am_to_inline_deriv_worker;

  /* The maximum AM of the VRR managers which will inline VRR workers
     Determined by libint.a */
  int max_am_manager_to_inline_vrr_worker;

  /* The maximum AM of the VRR managers which will inline DERIV workers */
  int max_am_manager_to_inline_deriv_worker;

  /* The maximum AM for which HRR workers will be made inline functions
     Determined by libint.a */
  int max_am_to_inline_hrr_worker;

  /* The maximum AM for which D1HRR workers will be made inline functions */
  int max_am_to_inline_d1hrr_worker;

  /* The maximum AM of the HRR/VRR managers which will inline HRR workers */
  int max_am_manager_to_inline_hrr_worker;

  /* The maximum AM of the HRR managers which will inline D1HRR workers */
  int max_am_manager_to_inline_d1hrr_worker;

  /* The maximum AM for which VRR managers will be inlined into HRR managers */
  int max_am_to_inline_vrr_manager;

} LibderivParams_t;
