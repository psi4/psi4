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
 * \file armci/src/x/dcmf/armcix_get.c
 * \brief DCMF ARMCI Extension for get operations.
 */

#include "armcix_impl.h"

DCMF_Protocol_t __get_protocol;

/**
 * \brief Register the DCMF ARMCI Extention get operation.
 *
 * \see DCMF_Get_register
 */
void ARMCIX_DCMF_Get_register ()
{
  DCMF_Get_Configuration_t configuration = {
    DCMF_DEFAULT_GET_PROTOCOL
  };
  DCMF_Get_register (&__get_protocol, &configuration);
}


/**
 * \brief ARMCI Extension blocking get operation.
 *
 * \param[in] src       Source buffer on the remote node
 * \param[in] dst       Destination buffer on the local node
 * \param[in] bytes     Number of bytes to transfer
 * \param[in] proc      Remote node rank
 *
 * \retval 0 Success
 */
int ARMCIX_Get(void * src, void * dst, int bytes, int proc)
{
  DCMF_CriticalSection_enter (0);

  volatile unsigned active = 1;
  DCMF_Callback_t cb_wait = { ARMCIX_DCMF_cb_decrement, (void *)&active };
  DCMF_Request_t request;

  DCMF_Memregion_t * src_memregion = &__connection[proc].remote_mem_region;
  DCMF_Memregion_t * dst_memregion = &__connection[proc].local_mem_region;

  DCMF_Result result =
    DCMF_Get (&__get_protocol,
              &request,
              cb_wait,
              DCMF_SEQUENTIAL_CONSISTENCY,
              proc,
              bytes,
              src_memregion,
              dst_memregion,
              armcix_dcmf_va_to_offset (src_memregion, src),
              armcix_dcmf_va_to_offset (dst_memregion, dst));

  while (active) DCMF_Messager_advance ();

  DCMF_CriticalSection_exit  (0);

  return (result != DCMF_SUCCESS);
}


/**
 * \brief ARMCI Extension non-blocking get operation.
 *
 * \param[in] src       Source buffer on the remote node
 * \param[in] dst       Destination buffer on the local node
 * \param[in] bytes     Number of bytes to transfer
 * \param[in] proc      Remote node rank
 * \param[in] nb_handle ARMCI non-blocking handle
 *
 * \return ???
 */
int ARMCIX_NbGet (void * src, void * dst, int bytes, int proc, armci_ihdl_t nb_handle)
{
  DCMF_CriticalSection_enter (0);

  armcix_dcmf_opaque_t * dcmf = (armcix_dcmf_opaque_t *) &nb_handle->cmpl_info;
  dcmf->active = 1;
  dcmf->connection = &__connection[proc];

  //fprintf (stderr, "ARMCIX_NbGet() dst=%p, src=%p, bytes=%d, request=%p\n", dst, src, bytes, &(dcmf->request));

  __connection[proc].active++;
  __global_connection.active++;

  DCMF_Callback_t cb_free = { ARMCIX_DCMF_NbOp_cb_done, nb_handle };
  ARMCIX_DCMF_Request_t * new_request = ARMCIX_DCMF_request_allocate (cb_free);
  DCMF_Callback_t cb_done = { (void(*)(void *)) ARMCIX_DCMF_request_free, new_request };

  DCMF_Memregion_t * src_memregion = &__connection[proc].remote_mem_region;
  DCMF_Memregion_t * dst_memregion = &__connection[proc].local_mem_region;

  DCMF_Result result =
    DCMF_Get (&__get_protocol,
              &(new_request->request),
              cb_done,
              DCMF_SEQUENTIAL_CONSISTENCY,
              proc,
              bytes,
              src_memregion,
              dst_memregion,
              armcix_dcmf_va_to_offset (src_memregion, src),
              armcix_dcmf_va_to_offset (dst_memregion, dst));

  DCMF_CriticalSection_exit (0);

  return (result != DCMF_SUCCESS);
}



/**
 * \brief ARMCI Extension blocking vector get operation.
 *
 * \param[in] darr      Descriptor array
 * \param[in] len       Length of descriptor array
 * \param[in] proc      Remote process(or) ID
 *
 * \return ???
 */
int ARMCIX_GetV (armci_giov_t * darr, int len, int proc)
{
  armci_ireq_t nb_request;
  armci_ihdl_t nb_handle = (armci_ihdl_t) &nb_request;
  ARMCIX_NbGetV (darr, len, proc, nb_handle);
  ARMCIX_Wait (&nb_handle->cmpl_info);

  return 0;
}


/**
 * \brief ARMCI Extension non-blocking vector get operation.
 *
 * \param[in] darr      Descriptor array
 * \param[in] len       Length of descriptor array
 * \param[in] proc      Remote process(or) ID
 * \param[in] nb_handle ARMCI non-blocking handle
 *
 * \return ???
 */
int ARMCIX_NbGetV (armci_giov_t * darr, int len, int proc, armci_ihdl_t nb_handle)
{
  DCMF_Result result = DCMF_ERROR;

  DCMF_CriticalSection_enter (0);

  //fprintf (stderr, "ARMCIX_NbGetV() >> len=%d, proc=%d\n", len, proc);

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

  //fprintf (stderr, "ARMCIX_NbGetV() -- n=%d, dcmf->active=%d, __connection[%d].active=%d, __global_connection.active=%d\n", n, dcmf->active, proc, __connection[proc].active, __global_connection.active);

  DCMF_Memregion_t * src_memregion = &__connection[proc].remote_mem_region;
  DCMF_Memregion_t * dst_memregion = &__connection[proc].local_mem_region;

  DCMF_Callback_t cb_free = { ARMCIX_DCMF_NbOp_cb_done, nb_handle };
  DCMF_Callback_t cb_done = { (void(*)(void *)) ARMCIX_DCMF_request_free, NULL };
  for (i = 0; i < len; i++)
  {
    for (j = 0; j < darr[i].ptr_array_len; j++)
    {
      //fprintf (stderr, "ARMCIX_NbGetV() -- src=%p, dst=%p, bytes=%d\n", darr[i].src_ptr_array[j], darr[i].dst_ptr_array[j], darr[i].bytes);
      ARMCIX_DCMF_Request_t * new_request = ARMCIX_DCMF_request_allocate (cb_free);
      cb_done.clientdata = new_request;

      result =
        DCMF_Get (&__get_protocol,
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

  //fprintf (stderr, "ARMCIX_NbGetV() << result=%d\n", result);
  DCMF_CriticalSection_exit  (0);

  return (result != DCMF_SUCCESS);
}


unsigned ARMCIX_DCMF_GetS_recurse (void * src_ptr, int * src_stride_arr, 
                                   void * dst_ptr, int * dst_stride_arr, 
                                   int * seg_count, int stride_levels, int proc,
                                   armci_ihdl_t nb_handle)
{
  unsigned num_requests = 0;

  //fprintf (stderr, "ARMCIX_DCMF_GetS_recurse() >> \n");

  if (stride_levels == 0)
  {
    //fprintf (stderr, "ARMCIX_DCMF_GetS_recurse() dst=%p, src=%p, bytes=%d, request=%p\n", dst_ptr, src_ptr, seg_count[0], request);

    DCMF_Callback_t cb_free = { ARMCIX_DCMF_NbOp_cb_done, nb_handle };
    ARMCIX_DCMF_Request_t * new_request = ARMCIX_DCMF_request_allocate (cb_free);
    DCMF_Callback_t cb_done = { (void(*)(void *)) ARMCIX_DCMF_request_free, new_request };

    DCMF_Memregion_t * src_memregion = &__connection[proc].remote_mem_region;
    DCMF_Memregion_t * dst_memregion = &__connection[proc].local_mem_region;

    DCMF_Get (&__get_protocol,
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
      num_requests += ARMCIX_DCMF_GetS_recurse (src_tmp, src_stride_arr, 
                                                dst_tmp, dst_stride_arr, 
                                                seg_count, (stride_levels-1), proc,
                                                nb_handle);

      src_tmp += src_stride_arr[(stride_levels-1)];
      dst_tmp += dst_stride_arr[(stride_levels-1)];
    }
  }

  //fprintf (stderr, "ARMCIX_DCMF_GetS_recurse() << num_requests = %d\n", num_requests);

  return num_requests;
}


/**
 * \brief ARMCI Extension blocking strided get operation.
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
int ARMCIX_GetS (void * src_ptr, int * src_stride_arr, 
                 void * dst_ptr, int * dst_stride_arr, 
                 int * seg_count, int stride_levels, int proc)
{
  armci_ireq_t nb_request;
  armci_ihdl_t nb_handle = (armci_ihdl_t) &nb_request;
  ARMCIX_NbGetS (src_ptr, src_stride_arr, dst_ptr, dst_stride_arr,
                 seg_count, stride_levels, proc, nb_handle);
  ARMCIX_Wait (&nb_handle->cmpl_info);

  return 0;
}

/**
 * \brief ARMCI Extension non-blocking strided get operation.
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
int ARMCIX_NbGetS (void * src_ptr, int * src_stride_arr, 
                   void * dst_ptr, int * dst_stride_arr, 
                   int * seg_count, int stride_levels, int proc,
                   armci_ihdl_t nb_handle)
{
  DCMF_CriticalSection_enter (0);

  assert (nb_handle != NULL);
  if (stride_levels > 0)
  {
    assert (src_stride_arr);
    assert (dst_stride_arr);
  }

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
  count = ARMCIX_DCMF_GetS_recurse (src_ptr, src_stride_arr, 
                                    dst_ptr, dst_stride_arr, 
                                    seg_count, stride_levels, proc,
                                    nb_handle);

  assert (n == count);

  DCMF_CriticalSection_exit  (0);

  return 0;
}
