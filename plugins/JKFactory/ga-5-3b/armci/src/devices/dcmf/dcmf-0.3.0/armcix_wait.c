#if HAVE_CONFIG_H
#   include "config.h"
#endif
/* begin_generated_IBM_copyright_prolog                             */
/*                                                                  */
/* ---------------------------------------------------------------- */
/* (C)Copyright IBM Corp.  2007, 2008                               */
/* IBM BSD License                                                  */
/* ---------------------------------------------------------------- */
/*                                                                  */
/* end_generated_IBM_copyright_prolog                               */
/**
 * \file armci/src/armcix/dcmf/armcix_wait.c
 * \brief DCMF ARMCI Extension wait functions.
 */

#include "armcix_impl.h"

/**
 * \brief DCMF ARMCI Extension blocking wait operation for a specifc request
 *
 * The armcix_opaque_t structure is an opaque object contains a
 * armcix_dcmf_opaque_t structure which is used to maintain DCMF
 * ARMCIX state information for an operation in progress.
 *
 * This function invokes DCMF_Messager_advance() until the operation
 * completes and its associated callback is invoked and decrements the
 * active count.
 *
 * \param[in] cmpl_info Pointer to the ARMCIX opaque object
 *
 * \todo define return values
 * \return 0
 *
 * \see armcix_dcmf_opaque_t
 */
int ARMCIX_Wait (armcix_opaque_t * cmpl_info)
{
  DCMF_CriticalSection_enter (0);
  armcix_dcmf_opaque_t * dcmf = (armcix_dcmf_opaque_t *) cmpl_info;
  while (dcmf->active) DCMF_Messager_advance();
  DCMF_CriticalSection_exit  (0);

  return 0;
}

/**
 * \brief DCMF ARMCI Extension blocking wait operation for all requests to a specific process
 *
 * This function invokes DCMF_Messager_advance() until all operations to the
 * specified process complete and the associated callbacks are invoked and
 * decrements the active count.
 *
 * \param[in] proc Remote process rank
 *
 * \todo define return values
 * \return 0
 *
 * \see ARMCIX_DCMF_Connection_t
 * \see __connection
 */
int ARMCIX_WaitProc (int proc)
{
  DCMF_CriticalSection_enter (0);
  while (__connection[proc].active) DCMF_Messager_advance();
  DCMF_CriticalSection_exit  (0);

  return 0;
}

/**
 * \brief DCMF ARMCI Extension blocking wait operation for all requests to all processes
 *
 * This function invokes DCMF_Messager_advance() until all operations to all
 * processes complete and the associated callbacks are invoked to
 * decrement the global active count.
 *
 * \todo define return values
 * \return 0
 *
 * \see ARMCIX_DCMF_Connection_t
 * \see __global_connection
 */
int ARMCIX_WaitAll ()
{
  DCMF_CriticalSection_enter (0);
  while (__global_connection.active) DCMF_Messager_advance();
  DCMF_CriticalSection_exit  (0);

  return 0;
}
