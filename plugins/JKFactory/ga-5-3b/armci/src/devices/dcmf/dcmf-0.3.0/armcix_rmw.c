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
 * \file armci/src/armcix/dcmf/armcix_rmw.c
 * \brief DCMF ARMCI Extension for rmw operations.
 */

#include "armcix_impl.h"
#include <stdio.h>

typedef union ARMCIX_DCMF_RMWRequest_t
{
  DCQuad quad[2];
  struct
  {
    int        op;
    void     * ploc;
    void     * prem;
    int        extra;
    unsigned * active;
  };
}
ARMCIX_DCMF_RMWRequest_t __attribute__ ((__aligned__ (16)));

typedef union ARMCIX_DCMF_RMWResponse_t
{
  DCQuad quad;
  struct
  {
    int        op;
    void     * ploc;
    unsigned * active;
    union
    {
      int      ival;
      long     lval;
    };
  };
}
ARMCIX_DCMF_RMWResponse_t __attribute__ ((__aligned__ (16)));


DCMF_Protocol_t __rmw_request_protocol;
DCMF_Protocol_t __rmw_response_protocol;



void ARMCIX_DCMF_RecvRMWRequest (void           * clientdata,
                                 const DCQuad   * msginfo,
                                 unsigned         count,
                                 unsigned         peer,
                                 const char     * src,
                                 unsigned         bytes)
{
  ARMCIX_DCMF_RMWRequest_t * info = (ARMCIX_DCMF_RMWRequest_t *) msginfo;

  /* Initialize the RMW response data                                        */
  ARMCIX_DCMF_RMWResponse_t response;
  response.op = info->op;
  response.ploc = info->ploc;
  response.active = info->active;

  //fprintf (stderr, "ARMCIX_DCMF_RecvRMWRequest() - info->op == %d, info->ploc == %p, info->prem == %p, info->extra == %d, peer == %d\n", info->op, info->ploc, info->prem, info->extra, peer);

  switch (info->op) 
  {
    case ARMCI_FETCH_AND_ADD:
    case ARMCI_SWAP:
      response.ival = *((int *) info->prem);
      break;
    case ARMCI_FETCH_AND_ADD_LONG:
    case ARMCI_SWAP_LONG:
      response.lval = *((long *) info->prem);
      break;
    default: 
      armci_die("rmw: operation not supported",info->op);
      break;
  }

  /* Send the rmw response.                                                  */
  DCMF_Control (&__rmw_response_protocol,
                DCMF_SEQUENTIAL_CONSISTENCY,
                peer,
                (DCMF_Control_t *) &response);

  /* Perform the swap or add                                                 */
  switch (info->op) 
  {
    case ARMCI_FETCH_AND_ADD:
      *((int *) info->prem) += info->extra;
      break;
    case ARMCI_FETCH_AND_ADD_LONG:
      *((long *) info->prem) += info->extra;
      break;
    case ARMCI_SWAP:
      *((int *) info->prem) = info->extra;
      break;
    case ARMCI_SWAP_LONG:
      *((long *) info->prem) = info->extra;
      break;
    default: 
      armci_die("rmw: operation not supported",info->op);
      break;
  }
}

/**
 * \brief Receive a read-modify-write response message.
 *
 * \param[in] clientdata Registered clientdata
 * \param[in] info       Read-Modify-Write data sent from the target node
 * \param[in] peer       Rank of the node that sent this control message
 *
 * \see ARMCIX_DCMF_Connection_t
 * \see ARMCIX_DCMF_Control_t
 * \see DCMF_RecvControl
 */
void ARMCIX_DCMF_ReceiveRMWResponse (void                 * clientdata,
                                     const DCMF_Control_t * info,
                                     unsigned               peer)
{
  ARMCIX_DCMF_RMWResponse_t * response   = (ARMCIX_DCMF_RMWResponse_t *) info;

  //fprintf (stderr, "ARMCIX_DCMF_ReceiveRMWResponse() - response->op == %d, response->ploc == %p, response->ival == %d, response->lval == %d, response->active == %p\n", response->op, response->ploc, response->ival, response->lval, response->active);

  switch (response->op)
  {
    case ARMCI_FETCH_AND_ADD:
    case ARMCI_SWAP:
      *((int *) response->ploc) = response->ival;
      break;
    case ARMCI_FETCH_AND_ADD_LONG:
    case ARMCI_SWAP_LONG:
      *((long *) response->ploc) = response->lval;
      break;
    default: 
      armci_die("rmw: operation not supported",response->op);
      break;
  }

  // Complete the rmw operation
  *(response->active) = 0;
}


/**
 * \brief Register the DCMF ARMCI Extention rmw operation.
 *
 * \see DCMF_Control_register
 * \see DCMF_Send_register
 */
void ARMCIX_DCMF_Rmw_register ()
{
  DCMF_CriticalSection_enter (0);

  DCMF_Send_Configuration_t request_configuration = {
    DCMF_DEFAULT_SEND_PROTOCOL,
    DCMF_DEFAULT_NETWORK,
    ARMCIX_DCMF_RecvRMWRequest,
    NULL,
    NULL,
    NULL
  };
  DCMF_Send_register (&__rmw_request_protocol, &request_configuration);

  DCMF_Control_Configuration_t response_configuration = {
    DCMF_DEFAULT_CONTROL_PROTOCOL,
    DCMF_DEFAULT_NETWORK,
    ARMCIX_DCMF_ReceiveRMWResponse,
    NULL
  };
  DCMF_Control_register (&__rmw_response_protocol, &response_configuration);

  DCMF_CriticalSection_exit (0);
}


/**
 * \brief ARMCI Extension blocking read-modify-write operation.
 *
 * \param[in] op
 * \param[in] ploc
 * \param[in] prem
 * \param[in] extra
 * \param[in] proc
 *
 * \retval ???
 */
int ARMCIX_Rmw (int op, int * ploc, int * prem, int extra, int proc)
{
  DCMF_CriticalSection_enter (0);

  volatile unsigned active = 1;
  
  //fprintf (stderr, "ARMCIX_Rmw() - op == %d, ploc == %p, prem == %p, extra == %d, proc == %d\n", op, ploc, prem, extra, proc);

  /* Initialize the RMW request data                                        */
  ARMCIX_DCMF_RMWRequest_t info;
  info.op = op;
  info.ploc = ploc;
  info.prem = prem;
  switch (op)
  {
    case ARMCI_FETCH_AND_ADD:
    case ARMCI_FETCH_AND_ADD_LONG:
      info.extra = extra;
      break;
    case ARMCI_SWAP:
    case ARMCI_SWAP_LONG:
      info.extra = *ploc;
      break;
    default: 
      armci_die("rmw: operation not supported",op);
      break;
  }

  info.active = (unsigned *)&active;

  DCMF_Request_t request;
  DCMF_Callback_t cb_wait = { NULL, NULL };

  DCMF_Send ( &__rmw_request_protocol,
              &request,
              cb_wait,
              DCMF_SEQUENTIAL_CONSISTENCY,
              proc,
              0,
              NULL,
              (DCQuad *)&info,
              2);

  //fprintf (stderr, "ARMCIX_Rmw() > active == %d (&active == %p)\n", active, &active);
  while (active) DCMF_Messager_advance ();
  //fprintf (stderr, "ARMCIX_Rmw() < active == %d (&active == %p)\n", active, &active);

  DCMF_CriticalSection_exit  (0);

  return 0;
}
