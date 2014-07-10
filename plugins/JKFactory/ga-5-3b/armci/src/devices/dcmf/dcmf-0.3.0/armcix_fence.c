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
 * \file armci/src/armcix/dcmf/armcix_fence.c
 * \brief DCMF ARMCI Extension for fence operations.
 */

#include "armcix_impl.h"
#include <stdio.h>

DCMF_Protocol_t __fence_rts_protocol;
DCMF_Protocol_t __fence_ack_protocol;

/**
 * \brief DCMF ARMCI Extention receive short fence request callback
 *
 * \see ARMCIX_Fence
 * \see ARMCIX_AllFence
 * \see ARMCIX_DCMF_ReceiveFenceAck
 * \see DCMF_RecvSendShort
 */
void ARMCIX_DCMF_ReceiveFenceRequest (void           * clientdata,
                                      const DCQuad   * msginfo,
                                      unsigned         count,
                                      unsigned         peer,
                                      const char     * src,
                                      unsigned         bytes)
{
  DCMF_Callback_t * cb = (DCMF_Callback_t *) msginfo;
  DCMF_Control_t info;
  DCMF_Callback_t * ack = (DCMF_Callback_t *) &info;
  *ack = *cb;

  DCMF_Control (&__fence_ack_protocol,
                DCMF_SEQUENTIAL_CONSISTENCY,
                peer,
                &info);
}

/**
 * \brief Receive a fence control message.
 *
 * The fence message type is either a fence \e request or a fence \e ack.
 *
 * \param[in] clientdata Registered clientdata, the armci connection array
 * \param[in] info       Fence control information
 * \param[in] peer       Rank of the node that sent this control message
 *
 * \see ARMCIX_Fence
 * \see ARMCIX_AllFence
 * \see ARMCIX_DCMF_ReceiveFenceRequest
 * \see DCMF_RecvControl
 */
void ARMCIX_DCMF_ReceiveFenceAck (void                 * clientdata,
                                  const DCMF_Control_t * info,
                                  unsigned               peer)
{
  DCMF_Callback_t * cb = (DCMF_Callback_t *) info;
  if (cb->function)
    cb->function(cb->clientdata, NULL);
}


/**
 * \brief Register the DCMF ARMCI Extention fence operation.
 *
 * \param[in]  connection_array Connection array
 *
 * \see DCMF_Control_register
 */
void ARMCIX_DCMF_Fence_register (ARMCIX_DCMF_Connection_t * connection_array)
{
  DCMF_CriticalSection_enter (0);

  DCMF_Send_Configuration_t send_configuration = {
    DCMF_DEFAULT_SEND_PROTOCOL,
    DCMF_DEFAULT_NETWORK,
    ARMCIX_DCMF_ReceiveFenceRequest,
    connection_array,
    NULL,
    NULL
  };
  DCMF_Send_register (&__fence_rts_protocol, &send_configuration);

  DCMF_Control_Configuration_t configuration = {
    DCMF_DEFAULT_CONTROL_PROTOCOL,
    DCMF_DEFAULT_NETWORK,
    ARMCIX_DCMF_ReceiveFenceAck,
    connection_array
  };
  DCMF_Control_register (&__fence_ack_protocol, &configuration);

  DCMF_CriticalSection_exit (0);
}


/**
 * \brief Point-to-point fence operation.
 *
 * Blocks until all active messages between the local node and the remote
 * node have completed and acknowledged by the remote node.
 *
 * \param[in] proc       Rank of the remote node to fence
 *
 * \see ARMCIX_AllFence
 * \see ARMCIX_DCMF_ReceiveFenceRequest
 * \see ARMCIX_DCMF_ReceiveFenceAck
 */
void ARMCIX_Fence (int proc)
{
  DCMF_CriticalSection_enter (0);

  DCMF_Request_t request;
  volatile unsigned active = 1;
  DCQuad quad;
  DCMF_Callback_t * cb = (DCMF_Callback_t *) &quad;
  cb->function   = ARMCIX_DCMF_cb_decrement;
  cb->clientdata = (void *) &active;
  DCMF_Send ( &__fence_rts_protocol,
              &request,
              (DCMF_Callback_t) { NULL, NULL },
              DCMF_SEQUENTIAL_CONSISTENCY,
              proc,
              0,
              NULL,
              (DCQuad *) &quad,
              1);

  while (active) DCMF_Messager_advance ();

  DCMF_CriticalSection_exit (0);
}

/**
 * \brief Global fence operation.
 *
 * Blocks until all active messages between the local node and all remote
 * nodes have completed and acknowledged by the remote node.
 *
 * \see ARMCIX_Fence
 * \see ARMCIX_DCMF_ReceiveFenceRequest
 * \see ARMCIX_DCMF_ReceiveFenceAck
 */
void ARMCIX_AllFence ()
{
  DCMF_CriticalSection_enter (0);

  unsigned size = DCMF_Messager_size ();
  unsigned peer;

  volatile unsigned active = 0;
  DCQuad quad;
  DCMF_Callback_t * cb = (DCMF_Callback_t *) &quad;
  cb->function   = ARMCIX_DCMF_cb_decrement;
  cb->clientdata = (void *) &active;

  DCMF_Callback_t cb_null = { NULL, NULL };
  DCMF_Callback_t cb_done = { (void (*)(void *, DCMF_Error_t *))ARMCIX_DCMF_request_free, NULL };
  for (peer = 0; peer < size; peer++)
  {
    ARMCIX_DCMF_Request_t * new_request = ARMCIX_DCMF_request_allocate (cb_null);
    cb_done.clientdata = new_request;

    active++;
    DCMF_Send ( &__fence_rts_protocol,
                &(new_request->request),
                cb_done,
                DCMF_SEQUENTIAL_CONSISTENCY,
                peer,
                0,
                NULL,
                (DCQuad *) &quad,
                1);

    while (active) DCMF_Messager_advance ();
  }

  DCMF_CriticalSection_exit (0);
}

void ARMCIX_Barrier ()
{
#warning implement ARMCIX_Barrier ?
assert (0);
}
