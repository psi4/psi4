/* begin_generated_IBM_copyright_prolog                             */
/*                                                                  */
/* ---------------------------------------------------------------- */
/* (C)Copyright IBM Corp.  2007, 2008                               */
/* IBM BSD License                                                  */
/* ---------------------------------------------------------------- */
/*                                                                  */
/* end_generated_IBM_copyright_prolog                               */
/**
 * \file armci/src/armcix/dcmf/armcix_impl.h
 * \brief DCMF ARMCI Extension implementation interface.
 */
#ifndef __armci_src_x_armcix_impl_h
#define __armci_src_x_armcix_impl_h

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "dcmf.h"
#include "armcix.h"

/* verify that the version of the installed dcmf library is compatible */
#if (DCMF_VERSION_RELEASE == 0)
  #if (DCMF_VERSION_MAJOR == 3)
    #if (DCMF_VERSION_MINOR < 0)
      #error Incompatible dcmf minor version
    #endif
  #else
    #error Incompatible dcmf major version
  #endif
#else
  #error Incompatible dcmf release version
#endif


/*
#define BLOCKING_OPERATIONS_REQUIRE_FENCE
#warning remove the previous #define if blocking put/acc operations do not require a fence
*/
typedef struct ARMCIX_DCMF_Connection_t
{
  DCMF_Request_t request;      /**< \todo lazy allocate the request object?  */
  unsigned active;             /**< Number of active messages to this peer   */
  unsigned peer;               /**< Maximum system size = 2^32               */
  struct
  {
    unsigned origin;
    unsigned target;
  } sequence;
  struct
  {
    unsigned watermark;
    unsigned origin     :1;
    unsigned target     :1;
    unsigned unused     :30;
  } fence;
  unsigned unused0;
  unsigned unused1;
  DCMF_Memregion_t local_mem_region;
  DCMF_Memregion_t remote_mem_region;
}
ARMCIX_DCMF_Connection_t __attribute__ ((__aligned__ (16)));

typedef struct ARMCIX_DCMF_Request_t
{
  DCMF_Request_t                 request;
  DCQuad                         quad[7];
} ARMCIX_DCMF_Request_t __attribute__ ((__aligned__ (16)));


typedef struct armcix_dcmf_opaque_t
{
  ARMCIX_DCMF_Connection_t * connection;
  unsigned                   active;
}
armcix_dcmf_opaque_t;


static inline void armcix_dcmf_compile_time_assert ()
{
  COMPILE_TIME_ASSERT(sizeof(armcix_dcmf_opaque_t)<=sizeof(armcix_opaque_t));
}

extern ARMCIX_DCMF_Connection_t __global_connection;
extern ARMCIX_DCMF_Connection_t * __connection;

static inline size_t armcix_dcmf_va_to_offset (DCMF_Memregion_t * mr, void * va)
{
  size_t bytes;
  void * base;
  DCMF_Memregion_query (mr, &bytes, &base);
  return ((size_t)va) - ((size_t)base);
}



/**
 * \brief Generic decrement callback
 *
 * \param[in] clientdata The variable to decrement
 */
void ARMCIX_DCMF_cb_decrement (void * clientdata, DCMF_Error_t *err);

/**
 * \brief Callback function for non-blocking operations
 *
 * \param[in] clientdata The non-blocking handle to complete
 */
void ARMCIX_DCMF_NbOp_cb_done (void * clientdata, DCMF_Error_t *err);

/**
 * \brief Allocate a request from the free request pool
 *
 * Attempt to increase the size of the request pool if a free request is not
 * available. Otherwise, if the maximum request pool size has been reached,
 * block until a request completes and becomes available.
 *
 * \param[in] cb_free Callback to invoke with the request is free'd
 *
 * \return A free request
 *
 * \see ARMCIX_DCMF_request_free
 */
ARMCIX_DCMF_Request_t * ARMCIX_DCMF_request_allocate (DCMF_Callback_t cb_free);

/**
 * \brief Release a request into the free request pool
 *
 * The callback associated with the request is invoked before the
 * request is released.
 *
 * \see ARMCIX_DCMF_request_allocate
 */
void ARMCIX_DCMF_request_free (ARMCIX_DCMF_Request_t * request);


/**
 * \brief Register the DCMF ARMCI Extention get operation.
 *
 * \see DCMF_Get_register
 */
void ARMCIX_DCMF_Get_register ();

/**
 * \brief Register the DCMF ARMCI Extention put operation.
 *
 * \param[in]  connection_array Connection array
 *
 * \see DCMF_Send_register
 */
void ARMCIX_DCMF_Put_register (ARMCIX_DCMF_Connection_t * connection_array);

/**
 * \brief Register the DCMF ARMCI Extention accumulate operation.
 *
 * \param[in]  connection_array Connection array
 *
 * \see DCMF_Send_register
 */
void ARMCIX_DCMF_Acc_register (ARMCIX_DCMF_Connection_t * connection_array);

/**
 * \brief Register the DCMF ARMCI Extention fence operation.
 *
 * \param[in]  connection_array Connection array
 *
 * \see DCMF_Control_register
 */
void ARMCIX_DCMF_Fence_register (ARMCIX_DCMF_Connection_t * connection_array);


/**
 * \brief Register the DCMF ARMCI Extention rmw operation.
 *
 * \see DCMF_Control_register
 * \see DCMF_Send_register
 */
void ARMCIX_DCMF_Rmw_register ();







#endif
