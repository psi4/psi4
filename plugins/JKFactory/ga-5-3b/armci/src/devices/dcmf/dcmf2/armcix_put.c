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
 * \file armci/src/x/dcmf/armcix_put.c
 * \brief DCMF ARMCI Extension for put operations.
 */

#include "armcix_impl.h"

DCMF_Protocol_t __put_protocol;

/**
 * \brief Register the DCMF ARMCI Extention put operation.
 *
 * \param[in]  connection_array Connection array
 *
 * \see DCMF_Send_register
 */
void ARMCIX_DCMF_Put_register (ARMCIX_DCMF_Connection_t * connection_array)
{
  DCMF_CriticalSection_enter (0);

  DCMF_Put_Configuration_t put_configuration = { DCMF_DEFAULT_PUT_PROTOCOL };
  DCMF_Put_register (&__put_protocol, &put_configuration);

  DCMF_CriticalSection_exit (0);
}

/**
 * \brief ARMCI Extension blocking put operation.
 *
 * \param[in] src       Source buffer on the local node
 * \param[in] dst       Destination buffer on the remote node
 * \param[in] bytes     Number of bytes to transfer
 * \param[in] proc      Remote node rank
 *
 * \return ???
 */
int ARMCIX_Put( void * src, void * dst, int bytes, int proc)
{
  DCMF_CriticalSection_enter (0);

  volatile unsigned active = 1;
  DCMF_Callback_t cb_wait = { ARMCIX_DCMF_cb_decrement, (void *)&active };
  DCMF_Request_t request;

  DCMF_Memregion_t * src_memregion = &__connection[proc].local_mem_region;
  DCMF_Memregion_t * dst_memregion = &__connection[proc].remote_mem_region;

  DCMF_Result result =
    DCMF_Put (&__put_protocol,
              &request,
              cb_wait,
              DCMF_SEQUENTIAL_CONSISTENCY,
              proc,
              bytes,
              src_memregion,
              dst_memregion,
              armcix_dcmf_va_to_offset (src_memregion, src),
              armcix_dcmf_va_to_offset (dst_memregion, dst));

#ifdef BLOCKING_OPERATIONS_REQUIRE_FENCE
  ARMCIX_Fence (proc);
#else
  while (active) DCMF_Messager_advance ();
#endif

  DCMF_CriticalSection_exit  (0);

  return (result != DCMF_SUCCESS);
}


/**
 * \brief ARMCI Extension non-blocking put operation.
 *
 * \param[in] src       Source buffer on the local node
 * \param[in] dst       Destination buffer on the remote node
 * \param[in] bytes     Number of bytes to transfer
 * \param[in] proc      Remote node rank
 * \param[in] nb_handle ARMCI non-blocking handle
 *
 * \return ???
 */
int ARMCIX_NbPut (void * src, void * dst, int bytes, int proc, armci_ihdl_t nb_handle)
{
  DCMF_CriticalSection_enter (0);

  armcix_dcmf_opaque_t * dcmf = (armcix_dcmf_opaque_t *) &nb_handle->cmpl_info;
  dcmf->active = 1;
  dcmf->connection = &__connection[proc];

  __connection[proc].active++;
  __global_connection.active++;

  DCMF_Callback_t cb_free = { ARMCIX_DCMF_NbOp_cb_done, nb_handle };
  ARMCIX_DCMF_Request_t * new_request = ARMCIX_DCMF_request_allocate (cb_free);
  DCMF_Callback_t cb_done = { (void(*)(void *)) ARMCIX_DCMF_request_free, new_request };

  DCMF_Memregion_t * src_memregion = &__connection[proc].local_mem_region;
  DCMF_Memregion_t * dst_memregion = &__connection[proc].remote_mem_region;

  DCMF_Result result =
    DCMF_Put (&__put_protocol,
              &(new_request->request),
              cb_done,
              DCMF_SEQUENTIAL_CONSISTENCY,
              proc,
              bytes,
              src_memregion,
              dst_memregion,
              armcix_dcmf_va_to_offset (src_memregion, src),
              armcix_dcmf_va_to_offset (dst_memregion, dst));

  DCMF_CriticalSection_exit  (0);

  return (result != DCMF_SUCCESS);
}



/**
 * \brief ARMCI Extension blocking vector put operation.
 *
 * \param[in] darr      Descriptor array
 * \param[in] len       Length of descriptor array
 * \param[in] proc      Remote process(or) ID
 *
 * \return ???
 */
int ARMCIX_PutV (armci_giov_t * darr, int len, int proc)
{
  armci_ireq_t nb_request;
  armci_ihdl_t nb_handle = (armci_ihdl_t) &nb_request;
  ARMCIX_NbPutV (darr, len, proc, nb_handle);

#ifdef BLOCKING_OPERATIONS_REQUIRE_FENCE
  ARMCIX_Fence (proc);
#else
  ARMCIX_Wait (&nb_handle->cmpl_info);
#endif

  return 0;
}


/**
 * \brief ARMCI Extension non-blocking vector put operation.
 *
 * \param[in] darr      Descriptor array
 * \param[in] len       Length of descriptor array
 * \param[in] proc      Remote process(or) ID
 * \param[in] nb_handle ARMCI non-blocking handle
 *
 * \return ???
 */
int ARMCIX_NbPutV (armci_giov_t * darr, int len, int proc, armci_ihdl_t nb_handle)
{
  DCMF_CriticalSection_enter (0);

  //fprintf (stderr, "ARMCIX_NbPutV() >> len=%d, proc=%d\n", len, proc);

  // Calculate the number of requests
  unsigned n = 0;
  unsigned i, j;
  for (i = 0; i < len; i++)
    for (j = 0; j < darr[i].ptr_array_len; j++)
      n++;

  armcix_dcmf_opaque_t * dcmf = (armcix_dcmf_opaque_t *) &nb_handle->cmpl_info;
  dcmf->connection = &__connection[proc];
  dcmf->active = n;

  __connection[proc].active += n;
  __global_connection.active += n;

  DCMF_Memregion_t * src_memregion = &__connection[proc].local_mem_region;
  DCMF_Memregion_t * dst_memregion = &__connection[proc].remote_mem_region;

  DCMF_Callback_t cb_free = { ARMCIX_DCMF_NbOp_cb_done, nb_handle };
  DCMF_Callback_t cb_done = { (void(*)(void *)) ARMCIX_DCMF_request_free, NULL };
  for (i = 0; i < len; i++)
  {
    for (j = 0; j < darr[i].ptr_array_len; j++)
    {
      //fprintf (stderr, "ARMCIX_NbPutV() -- src=%p, dst=%p, bytes=%d\n", darr[i].src_ptr_array[j], darr[i].dst_ptr_array[j], darr[i].bytes);
      ARMCIX_DCMF_Request_t * new_request = ARMCIX_DCMF_request_allocate (cb_free);
      cb_done.clientdata = new_request;

      DCMF_Put (&__put_protocol,
                &(new_request->request),
                cb_done,
                DCMF_SEQUENTIAL_CONSISTENCY,
                proc,
                darr[i].bytes,
                src_memregion,
                dst_memregion,
                armcix_dcmf_va_to_offset (src_memregion, darr[i].src_ptr_array[j]),
                armcix_dcmf_va_to_offset (dst_memregion, darr[i].dst_ptr_array[j]));
    }
  }
  //fprintf (stderr, "ARMCIX_NbPutV() <<\n");

  DCMF_CriticalSection_exit  (0);

  return 0;
}



unsigned ARMCIX_DCMF_PutS_recurse (void * src_ptr, int * src_stride_arr, 
                                   void * dst_ptr, int * dst_stride_arr, 
                                   int * seg_count, int stride_levels, int proc,
                                   armci_ihdl_t nb_handle)
{
  unsigned num_requests = 0;

  //fprintf (stderr, "ARMCIX_DCMF_PutS_recurse() >> \n");

  if (stride_levels == 0)
  {
    //fprintf (stderr, "ARMCIX_DCMF_PutS_recurse() dst=%p, src=%p, bytes=%d, nb_handle=%p\n", dst_ptr, src_ptr, seg_count[0], nb_handle);

    DCMF_Callback_t cb_free = { ARMCIX_DCMF_NbOp_cb_done, nb_handle };
    ARMCIX_DCMF_Request_t * new_request = ARMCIX_DCMF_request_allocate (cb_free);
    DCMF_Callback_t cb_done = { (void(*)(void *)) ARMCIX_DCMF_request_free, new_request };

    DCMF_Memregion_t * src_memregion = &__connection[proc].local_mem_region;
    DCMF_Memregion_t * dst_memregion = &__connection[proc].remote_mem_region;

    DCMF_Put (&__put_protocol,
              &(new_request->request),
              cb_done,
              DCMF_SEQUENTIAL_CONSISTENCY,
              proc,
              seg_count[0],
              src_memregion,
              dst_memregion,
              armcix_dcmf_va_to_offset (src_memregion, src_ptr),
              armcix_dcmf_va_to_offset (dst_memregion, dst_ptr));

    num_requests++;
  }
  else
  {
    char * src_tmp = (char *) src_ptr;
    char * dst_tmp = (char *) dst_ptr;
    unsigned i;
    for (i = 0; i < seg_count[stride_levels]; i++)
    {
      num_requests += ARMCIX_DCMF_PutS_recurse (src_tmp, src_stride_arr, 
                                                dst_tmp, dst_stride_arr, 
                                                seg_count, (stride_levels-1), proc,
                                                nb_handle);

      src_tmp += src_stride_arr[(stride_levels-1)];
      dst_tmp += dst_stride_arr[(stride_levels-1)];
    }
  }

  //fprintf (stderr, "ARMCIX_DCMF_PutS_recurse() << num_requests = %d\n", num_requests);

  return num_requests;
}


/**
 * \brief ARMCI Extension blocking strided put operation.
 *
 * \param[in] src_ptr        pointer to 1st segment at source
 * \param[in] src_stride_arr array of strides at source
 * \param[in] dst_ptr        pointer to 1st segment at destination
 * \param[in] dst_stride_arr array of strides at destination
 * \param[in] seg_count      number of segments at each stride levels: count[0]=bytes
 * \param[in] stride_levels  number of stride levels
 * \param[in] proc           remote process(or) ID
 *
 * \return ???
 */
int ARMCIX_PutS (void * src_ptr, int * src_stride_arr, 
                 void * dst_ptr, int * dst_stride_arr, 
                 int * seg_count, int stride_levels, int proc)
{
  armci_ireq_t nb_request;
  armci_ihdl_t nb_handle = (armci_ihdl_t) &nb_request;
  ARMCIX_NbPutS (src_ptr, src_stride_arr, dst_ptr, dst_stride_arr,
                 seg_count, stride_levels, proc, nb_handle);

#ifdef BLOCKING_OPERATIONS_REQUIRE_FENCE
  ARMCIX_Fence (proc);
#else
  ARMCIX_Wait (&nb_handle->cmpl_info);
#endif

  return 0;
}

/**
 * \brief ARMCI Extension non-blocking strided put operation.
 *
 * \param[in] src_ptr        pointer to 1st segment at source
 * \param[in] src_stride_arr array of strides at source
 * \param[in] dst_ptr        pointer to 1st segment at destination
 * \param[in] dst_stride_arr array of strides at destination
 * \param[in] seg_count      number of segments at each stride levels: count[0]=bytes
 * \param[in] stride_levels  number of stride levels
 * \param[in] proc           remote process(or) ID
 * \param[in] nb_handle      ARMCI non-blocking handle
 *
 * \return ???
 */
int ARMCIX_NbPutS (void * src_ptr, int * src_stride_arr, 
                   void * dst_ptr, int * dst_stride_arr, 
                   int * seg_count, int stride_levels, int proc,
                   armci_ihdl_t nb_handle)
{
  DCMF_CriticalSection_enter (0);

  // Calculate the number of requests
  unsigned i;
  unsigned n = 1;
  for (i = 0; i < stride_levels; i++) n = n * seg_count[i+1];

  armcix_dcmf_opaque_t * dcmf = (armcix_dcmf_opaque_t *) &nb_handle->cmpl_info;
  dcmf->connection = &__connection[proc];
  dcmf->active = n;

  __connection[proc].active += n;
  __global_connection.active += n;

  unsigned count;
  count = ARMCIX_DCMF_PutS_recurse (src_ptr, src_stride_arr, 
                                    dst_ptr, dst_stride_arr, 
                                    seg_count, stride_levels, proc,
                                    nb_handle);

  //fprintf (stderr, "ARMCIX_NbPutS() -- n=%d == count=%d\n", n, count);
  assert (n == count);

  DCMF_CriticalSection_exit  (0);

  return 0;
}
