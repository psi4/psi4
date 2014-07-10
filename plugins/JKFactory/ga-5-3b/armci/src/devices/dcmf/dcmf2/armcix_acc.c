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
 * \file armci/src/x/dcmf/armcix_acc.c
 * \brief DCMF ARMCI Extension for accumulate operations.
 */

#include "armcix_impl.h"

typedef struct {
    float real;
    float imag;
} complex_t;

typedef struct {
    double real;
    double imag;
} dcomplex_t;

typedef union ARMCIX_DCMF_AccInfo_t
{
  DCQuad         raw[2];
  struct
  {
    void       * dst;
    int          datatype;
    unsigned     bytes;
    union
    {
      int        ival;
      long       lval;
      float      fval;
      double     dval;
      complex_t  cplxval;
      dcomplex_t dcplxval;
    };
  };
}
ARMCIX_DCMF_AccInfo_t __attribute__ ((__aligned__ (16)));

typedef struct ARMCIX_DCMF_AccNbInfo_t
{
  ARMCIX_DCMF_AccInfo_t      info;
  ARMCIX_DCMF_Connection_t * connection;
  void                     * buffer;
}
ARMCIX_DCMF_AccNbInfo_t;

#define ACCUMULATE( DTYPE, scale, elems, src, dst) {\
    int j;\
    DTYPE *a =(DTYPE *)(dst);\
    DTYPE *b =(DTYPE *)(src);\
    DTYPE alpha = *(DTYPE *)(scale);\
    for(j=0;j<(elems);j++)a[j] += alpha*b[j];\
}
        
#define CPL_ACCUMULATE( DTYPE, scale, elems, src, dst) {\
    int j;\
    DTYPE *a =(DTYPE *)(dst);\
    DTYPE *b =(DTYPE *)(src);\
    DTYPE alpha = *(DTYPE *)(scale);\
    for(j=0;j<(elems);j++){\
        a[j].real += alpha.real*b[j].real - alpha.imag*b[j].imag;\
        a[j].imag += alpha.imag*b[j].real + alpha.real*b[j].imag;\
    }\
}

DCMF_Protocol_t __acc_protocol;

/**
 * \brief DCMF ARMCI Extention receive short accumulate operation callback
 *
 * \see DCMF_RecvSendShort
 */
void ARMCIX_DCMF_RecvAcc1 (void           * clientdata,
                           const DCQuad   * msginfo,
                           unsigned         count,
                           unsigned         peer,
                           const char     * src,
                           unsigned         bytes)
{
  //ARMCIX_DCMF_Connection_t * connection = (ARMCIX_DCMF_Connection_t *) clientdata;
  ARMCIX_DCMF_AccInfo_t * info = (ARMCIX_DCMF_AccInfo_t *) msginfo;

  switch (info->datatype)
  {
    case ARMCI_ACC_INT:
      ACCUMULATE( int, &info->ival, bytes/sizeof(int), src, info->dst);
      break;
    case ARMCI_ACC_DBL:
      ACCUMULATE( double, &info->dval, bytes/sizeof(double), src, info->dst);
      break;
    case ARMCI_ACC_FLT:
      ACCUMULATE( float, &info->fval, bytes/sizeof(float), src, info->dst);
      break;
    case ARMCI_ACC_CPL:
      CPL_ACCUMULATE( complex_t, &info->cplxval, bytes/sizeof(complex_t), src, info->dst);
      break;
    case ARMCI_ACC_DCP:
      CPL_ACCUMULATE( dcomplex_t, &info->dcplxval, bytes/sizeof(dcomplex_t), src, info->dst);
      break;
    case ARMCI_ACC_LNG:
      ACCUMULATE( long, &info->lval, bytes/sizeof(long), src, info->dst);
      break;
    default:
      assert (0);
      break;
  }
}


void ARMCIX_DCMF_AccReceiveComplete (ARMCIX_DCMF_AccNbInfo_t * nbinfo)
{
  ARMCIX_DCMF_AccInfo_t * info = (ARMCIX_DCMF_AccInfo_t *) &nbinfo->info;

  switch (info->datatype)
  {
    case ARMCI_ACC_INT:
      ACCUMULATE( int, &info->ival, info->bytes/sizeof(int), nbinfo->buffer, info->dst);
      break;
    case ARMCI_ACC_DBL:
      ACCUMULATE( double, &info->dval, info->bytes/sizeof(double), nbinfo->buffer, info->dst);
      break;
    case ARMCI_ACC_FLT:
      ACCUMULATE( float, &info->fval, info->bytes/sizeof(float), nbinfo->buffer, info->dst);
      break;
    case ARMCI_ACC_CPL:
      CPL_ACCUMULATE( complex_t, &info->cplxval, info->bytes/sizeof(complex_t), nbinfo->buffer, info->dst);
      break;
    case ARMCI_ACC_DCP:
      CPL_ACCUMULATE( dcomplex_t, &info->dcplxval, info->bytes/sizeof(dcomplex_t), nbinfo->buffer, info->dst);
      break;
    case ARMCI_ACC_LNG:
      ACCUMULATE( long, &info->lval, info->bytes/sizeof(long), nbinfo->buffer, info->dst);
      break;
    default:
      assert (0);
      break;
  }

  free (nbinfo);
}


/**
 * \brief DCMF ARMCI Extention receive accumulate operation callback
 *
 * \see DCMF_RecvSend
 */
DCMF_Request_t * ARMCIX_DCMF_RecvAcc2 (void             * clientdata,
                                       const DCQuad     * msginfo,
                                       unsigned           count,
                                       unsigned           peer,
                                       unsigned           sndlen,
                                       unsigned         * rcvlen,
                                       char            ** rcvbuf,
                                       DCMF_Callback_t  * cb_done)
{
  ARMCIX_DCMF_Connection_t * connection = (ARMCIX_DCMF_Connection_t *) clientdata;

  ARMCIX_DCMF_AccNbInfo_t * nbinfo = (ARMCIX_DCMF_AccNbInfo_t *) malloc (sizeof(ARMCIX_DCMF_AccNbInfo_t) + sndlen);

  memcpy (&nbinfo->info, msginfo, sizeof(ARMCIX_DCMF_AccInfo_t));
  nbinfo->connection = &connection[peer];
  nbinfo->buffer = (void *)(((char *) nbinfo) + sizeof(ARMCIX_DCMF_AccNbInfo_t));

  *rcvlen = sndlen;
  *rcvbuf = nbinfo->buffer;

  cb_done->function   = (void *) ARMCIX_DCMF_AccReceiveComplete;
  cb_done->clientdata = (void *) nbinfo;

  return &connection[peer].request;
}


/**
 * \brief Register the DCMF ARMCI Extention accumulate operation.
 *
 * \param[in]  connection_array Connection array
 *
 * \see DCMF_Send_register
 */
void ARMCIX_DCMF_Acc_register (ARMCIX_DCMF_Connection_t * connection_array)
{
  DCMF_CriticalSection_enter (0);

  DCMF_Send_Configuration_t configuration = {
    DCMF_DEFAULT_SEND_PROTOCOL,
    ARMCIX_DCMF_RecvAcc1,
    connection_array,
    ARMCIX_DCMF_RecvAcc2,
    connection_array
  };
  DCMF_Send_register (&__acc_protocol, &configuration);

  DCMF_CriticalSection_exit (0);
}

/**
 * \brief ARMCI Extension blocking accumulate operation.
 *
 * \param[in] datatype  accumulate datatype (operation code)
 * \param[in] scale     opaque pointer to the scaling factor for accumulate
 * \param[in] src       Source buffer on the local node
 * \param[in] dst       Destination buffer on the remote node
 * \param[in] bytes     Number of bytes to transfer
 * \param[in] proc      Remote node rank
 *
 * \return ???
 */
int ARMCIX_Acc (int datatype, void * scale, void * src, void * dst, int bytes, int proc)
{
  DCMF_CriticalSection_enter (0);

  volatile unsigned active = 1;
  DCMF_Callback_t cb_wait = { ARMCIX_DCMF_cb_decrement, (void *)&active };
  DCMF_Request_t request;

  ARMCIX_DCMF_AccInfo_t info;
  info.dst = dst;
  info.bytes = bytes;
  info.datatype = datatype;
  switch (datatype)
  {
    case ARMCI_ACC_INT:
      info.ival = *((int *)scale);
      break;
    case ARMCI_ACC_DBL:
      info.dval = *((double *)scale);
      break;
    case ARMCI_ACC_FLT:
      info.fval = *((float *)scale);
      break;
    case ARMCI_ACC_CPL:
      info.cplxval.real = ((complex_t *)scale)->real;
      info.cplxval.imag = ((complex_t *)scale)->imag;
      break;
    case ARMCI_ACC_DCP:
      info.dcplxval.real = ((dcomplex_t *)scale)->real;
      info.dcplxval.imag = ((dcomplex_t *)scale)->imag;
      break;
    case ARMCI_ACC_LNG:
      info.lval = *((long *)scale);
      break;
    default:
      assert (0);
      break;
  }

  DCMF_Send ( &__acc_protocol,
              &request,
              cb_wait,
              DCMF_SEQUENTIAL_CONSISTENCY,
              proc,
              bytes,
              (char *) src,
              (DCQuad *) &info,
              2);

#ifdef BLOCKING_OPERATIONS_REQUIRE_FENCE
  ARMCIX_Fence (proc);
#else
  while (active) DCMF_Messager_advance ();
#endif

  DCMF_CriticalSection_exit  (0);

  return 0;
}


/**
 * \brief ARMCI Extension non-blocking accumulate operation.
 *
 * \param[in] datatype  accumulate datatype (operation code)
 * \param[in] scale     opaque pointer to the scaling factor for accumulate
 * \param[in] src       Source buffer on the local node
 * \param[in] dst       Destination buffer on the remote node
 * \param[in] bytes     Number of bytes to transfer
 * \param[in] proc      Remote node rank
 * \param[in] nb_handle ARMCI non-blocking handle
 *
 * \return ???
 */
int ARMCIX_NbAcc (int datatype, void * scale, void * src, void * dst, int bytes, int proc, armci_ihdl_t nb_handle)
{
  DCMF_CriticalSection_enter (0);

  armcix_dcmf_opaque_t * dcmf = (armcix_dcmf_opaque_t *) &nb_handle->cmpl_info;
  dcmf->active++;
  dcmf->connection = &__connection[proc];

  __connection[proc].active++;
  __global_connection.active++;

  DCMF_Callback_t cb_free = { ARMCIX_DCMF_NbOp_cb_done, nb_handle };
  ARMCIX_DCMF_Request_t * new_request = ARMCIX_DCMF_request_allocate (cb_free);
  DCMF_Callback_t cb_done = { (void(*)(void *)) ARMCIX_DCMF_request_free, new_request };

  ARMCIX_DCMF_AccInfo_t * info = (ARMCIX_DCMF_AccInfo_t *) &(new_request->quad[0]);
  info->dst = dst;
  info->bytes = bytes;
  info->datatype = datatype;
  switch (datatype)
  {
    case ARMCI_ACC_INT:
      info->ival = *((int *)scale);
      break;
    case ARMCI_ACC_DBL:
      info->dval = *((double *)scale);
      break;
    case ARMCI_ACC_FLT:
      info->fval = *((float *)scale);
      break;
    case ARMCI_ACC_CPL:
      info->cplxval.real = ((complex_t *)scale)->real;
      info->cplxval.imag = ((complex_t *)scale)->imag;
      break;
    case ARMCI_ACC_DCP:
      info->dcplxval.real = ((dcomplex_t *)scale)->real;
      info->dcplxval.imag = ((dcomplex_t *)scale)->imag;
      break;
    case ARMCI_ACC_LNG:
      info->lval = *((long *)scale);
      break;
    default:
      assert (0);
      break;
  }

  DCMF_Send ( &__acc_protocol,
              &(new_request->request),
              cb_done,
              DCMF_SEQUENTIAL_CONSISTENCY,
              proc,
              bytes,
              (char *) src,
              (DCQuad *) info,
              2);

  DCMF_CriticalSection_exit  (0);

  return 0;
}


/**
 * \brief ARMCI Extension blocking vector accumulate operation.
 *
 * \todo something goofy with AccV .. should be able to replace with combination of
 *       ARMCIX_AccV and ARMCIX_Wait, but that causes test-ibm.x to hang. Maybe
 *       related to interrupts?
 *
 * \param[in] datatype accumulate datatype (operation code)
 * \param[in] scale opaque pointer to the scaling factor for accumulate
 * \param[in] darr descriptor array
 * \param[in] len length of the descriptor array
 * \param[in] proc process(or) ID
 */
int ARMCIX_AccV (int datatype, void * scale, armci_giov_t * darr, int len, int proc)
{
#if 0
#error causes test-ibm.x to hang!
  armci_ireq_t nb_request;
  armci_ihdl_t nb_handle = (armci_ihdl_t) &nb_request;
  ARMCIX_NbAccV (datatype, scale, darr, len, proc, nb_handle);
#warning remove this ARMCIX_Fence() and implement some sort of ack scheme.
  ARMCIX_Fence (proc);
  ARMCIX_Wait (&nb_handle->cmpl_info);
#else
  DCMF_CriticalSection_enter (0);

  // Calculate the number of requests
  unsigned n = 0;
  unsigned i, j;
  for (i = 0; i < len; i++)
    for (j = 0; j < darr[i].ptr_array_len; j++)
      n++;

  volatile unsigned active = n;
  DCMF_Callback_t cb_wait = { ARMCIX_DCMF_cb_decrement, (void *)&active };
  DCMF_Request_t request[n];

  ARMCIX_DCMF_AccInfo_t info;
  info.datatype = datatype;
  switch (datatype)
  {
    case ARMCI_ACC_INT:
      info.ival = *((int *)scale);
      break;
    case ARMCI_ACC_DBL:
      info.dval = *((double *)scale);
      break;
    case ARMCI_ACC_FLT:
      info.fval = *((float *)scale);
      break;
    case ARMCI_ACC_CPL:
      info.cplxval.real = ((complex_t *)scale)->real;
      info.cplxval.imag = ((complex_t *)scale)->imag;
      break;
    case ARMCI_ACC_DCP:
      info.dcplxval.real = ((dcomplex_t *)scale)->real;
      info.dcplxval.imag = ((dcomplex_t *)scale)->imag;
      break;
    case ARMCI_ACC_LNG:
      info.lval = *((long *)scale);
      break;
    default:
      assert (0);
      break;
  }

  for (i = 0; i < len; i++)
  {
    info.bytes = darr[i].bytes;
    for (j = 0; j < darr[i].ptr_array_len; j++)
    {
      info.dst = darr[i].dst_ptr_array[j];
      DCMF_Send ( &__acc_protocol,
                  &request[--n],
                  cb_wait,
                  DCMF_SEQUENTIAL_CONSISTENCY,
                  proc,
                  info.bytes,
                  (char *) darr[i].src_ptr_array[j],
                  (DCQuad *) &info,
                  2);
    }
  }

#ifdef BLOCKING_OPERATIONS_REQUIRE_FENCE
  ARMCIX_Fence (proc);
#else
  // Poll until all accumulate messages have been sent.
  while (active) DCMF_Messager_advance ();
#endif

  DCMF_CriticalSection_exit  (0);
#endif
  return 0;
}


/**
 * \brief ARMCI Extension non-blocking vector accumulate operation.
 *
 * \param[in] datatype  accumulate datatype (operation code)
 * \param[in] scale     opaque pointer to the scaling factor for accumulate
 * \param[in] darr      Descriptor array
 * \param[in] len       Length of descriptor array
 * \param[in] proc      Remote process(or) ID
 * \param[in] nb_handle ARMCI non-blocking handle
 *
 * \return ???
 */
int ARMCIX_NbAccV (int datatype, void * scale, armci_giov_t * darr, int len, int proc, armci_ihdl_t nb_handle)
{
  DCMF_CriticalSection_enter (0);

  // Calculate the number of requests
  unsigned n = 0;
  unsigned i, j;
  for (i = 0; i < len; i++)
    for (j = 0; j < darr[i].ptr_array_len; j++)
      n++;

  armcix_dcmf_opaque_t * dcmf = (armcix_dcmf_opaque_t *) &nb_handle->cmpl_info;
  dcmf->connection = &__connection[proc];
  dcmf->active += n;

  __connection[proc].active += n;
  __global_connection.active += n;
#if 0
  ARMCIX_DCMF_AccInfo_t info;
  info.datatype = datatype;
  switch (datatype)
  {
    case ARMCI_ACC_INT:
      info.ival = *((int *)scale);
      break;
    case ARMCI_ACC_DBL:
      info.dval = *((double *)scale);
      break;
    case ARMCI_ACC_FLT:
      info.fval = *((float *)scale);
      break;
    case ARMCI_ACC_CPL:
      info.cplxval.real = ((complex_t *)scale)->real;
      info.cplxval.imag = ((complex_t *)scale)->imag;
      break;
    case ARMCI_ACC_DCP:
      info.dcplxval.real = ((dcomplex_t *)scale)->real;
      info.dcplxval.imag = ((dcomplex_t *)scale)->imag;
      break;
    case ARMCI_ACC_LNG:
      info.lval = *((long *)scale);
      break;
    default:
      assert (0);
      break;
  }
#endif

  for (i = 0; i < len; i++)
  {
    //info.bytes = darr[i].bytes;
    for (j = 0; j < darr[i].ptr_array_len; j++)
    {
      DCMF_Callback_t cb_free = { ARMCIX_DCMF_NbOp_cb_done, nb_handle };
      ARMCIX_DCMF_Request_t * new_request = ARMCIX_DCMF_request_allocate (cb_free);
      DCMF_Callback_t cb_done = { (void(*)(void *)) ARMCIX_DCMF_request_free, new_request };

      ARMCIX_DCMF_AccInfo_t * info = (ARMCIX_DCMF_AccInfo_t *) &(new_request->quad[0]);
      //info->dst = dst;
      //info->bytes = bytes;
      info->datatype = datatype;
      info->bytes = darr[i].bytes;
      switch (datatype)
      {
        case ARMCI_ACC_INT:
          info->ival = *((int *)scale);
          break;
        case ARMCI_ACC_DBL:
          info->dval = *((double *)scale);
          break;
        case ARMCI_ACC_FLT:
          info->fval = *((float *)scale);
          break;
        case ARMCI_ACC_CPL:
          info->cplxval.real = ((complex_t *)scale)->real;
          info->cplxval.imag = ((complex_t *)scale)->imag;
          break;
        case ARMCI_ACC_DCP:
          info->dcplxval.real = ((dcomplex_t *)scale)->real;
          info->dcplxval.imag = ((dcomplex_t *)scale)->imag;
          break;
        case ARMCI_ACC_LNG:
          info->lval = *((long *)scale);
          break;
        default:
          assert (0);
          break;
      }

      info->dst = darr[i].dst_ptr_array[j];
      DCMF_Send ( &__acc_protocol,
                  &(new_request->request),
                  cb_done,
                  DCMF_SEQUENTIAL_CONSISTENCY,
                  proc,
                  info->bytes,
                  (char *) darr[i].src_ptr_array[j],
                  (DCQuad *) info,
                  2);
    }
  }

  DCMF_CriticalSection_exit  (0);

  return 0;
}




unsigned ARMCIX_DCMF_AccS_recurse (int datatype, void * scale,
                                   void * src_ptr, int * src_stride_arr,
                                   void * dst_ptr, int * dst_stride_arr,
                                   int * seg_count, int stride_levels, int proc,
                                   armci_ihdl_t nb_handle)
{
  unsigned num_requests = 0;

  //fprintf (stderr, "ARMCIX_DCMF_AccS_recurse() >> \n");

  if (stride_levels == 0)
  {
    DCMF_Callback_t cb_free = { ARMCIX_DCMF_NbOp_cb_done, nb_handle };
    ARMCIX_DCMF_Request_t * new_request = ARMCIX_DCMF_request_allocate (cb_free);
    DCMF_Callback_t cb_done = { (void(*)(void *)) ARMCIX_DCMF_request_free, new_request };

    ARMCIX_DCMF_AccInfo_t * info = (ARMCIX_DCMF_AccInfo_t *) &(new_request->quad[0]);
    //info->dst = dst;
    //info->bytes = bytes;
    info->datatype = datatype;
    switch (datatype)
    {
      case ARMCI_ACC_INT:
        info->ival = *((int *)scale);
        break;
      case ARMCI_ACC_DBL:
        info->dval = *((double *)scale);
        break;
      case ARMCI_ACC_FLT:
        info->fval = *((float *)scale);
        break;
      case ARMCI_ACC_CPL:
        info->cplxval.real = ((complex_t *)scale)->real;
        info->cplxval.imag = ((complex_t *)scale)->imag;
        break;
      case ARMCI_ACC_DCP:
        info->dcplxval.real = ((dcomplex_t *)scale)->real;
        info->dcplxval.imag = ((dcomplex_t *)scale)->imag;
        break;
      case ARMCI_ACC_LNG:
        info->lval = *((long *)scale);
        break;
      default:
        assert (0);
        break;
      }
#if 0
    ARMCIX_DCMF_AccInfo_t info;
    info.datatype = datatype;
    switch (datatype)
    {
      case ARMCI_ACC_INT:
        info.ival = *((int *)scale);
        break;
      case ARMCI_ACC_DBL:
        info.dval = *((double *)scale);
        break;
      case ARMCI_ACC_FLT:
        info.fval = *((float *)scale);
        break;
      case ARMCI_ACC_CPL:
        info.cplxval.real = ((complex_t *)scale)->real;
        info.cplxval.imag = ((complex_t *)scale)->imag;
        break;
      case ARMCI_ACC_DCP:
        info.dcplxval.real = ((dcomplex_t *)scale)->real;
        info.dcplxval.imag = ((dcomplex_t *)scale)->imag;
        break;
      case ARMCI_ACC_LNG:
        info.lval = *((long *)scale);
        break;
      default:
        assert (0);
        break;
    }

    DCMF_Callback_t cb_free = { ARMCIX_DCMF_NbOp_cb_done, nb_handle };
    DCMF_Request_t * new_request = ARMCIX_DCMF_request_allocate (cb_free);
    DCMF_Callback_t cb_done = { (void(*)(void *)) ARMCIX_DCMF_request_free, new_request };
#endif
    info->bytes = seg_count[0];
    info->dst = dst_ptr;

    DCMF_Send ( &__acc_protocol,
                &(new_request->request),
                cb_done,
                DCMF_SEQUENTIAL_CONSISTENCY,
                proc,
                info->bytes,
                (char *) src_ptr,
                (DCQuad *) info,
                2);

    num_requests++;
  }
  else
  {
    char * src_tmp = (char *) src_ptr;
    char * dst_tmp = (char *) dst_ptr;
    unsigned i;
    for (i = 0; i < seg_count[stride_levels]; i++)
    {
      num_requests += ARMCIX_DCMF_AccS_recurse (datatype, scale,
                                                src_tmp, src_stride_arr,
                                                dst_tmp, dst_stride_arr,
                                                seg_count, (stride_levels-1), proc,
                                                nb_handle);

      src_tmp += src_stride_arr[(stride_levels-1)];
      dst_tmp += dst_stride_arr[(stride_levels-1)];
    }
  }

  //fprintf (stderr, "ARMCIX_DCMF_AccS_recurse() << num_requests = %d\n", num_requests);

  return num_requests;
}


/**
 * \brief ARMCI Extension blocking strided accumulate operation.
 *
 * \param[in] datatype       accumulate datatype (operation code)
 * \param[in] scale          opaque pointer to the scaling factor for accumulate
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
int ARMCIX_AccS (int datatype, void * scale,
                 void * src_ptr, int * src_stride_arr, 
                 void * dst_ptr, int * dst_stride_arr, 
                 int * seg_count, int stride_levels, int proc)
{
#if 0
#error causes test-ibm.x to hang!
  armci_ireq_t nb_request;
  armci_ihdl_t nb_handle = (armci_ihdl_t) &nb_request;
  ARMCIX_NbAccS (datatype, scale,
                 src_ptr, src_stride_arr,
                 dst_ptr, dst_stride_arr,
                 seg_count, stride_levels, proc,
                 nb_handle);
#warning remove this ARMCIX_Fence() and implement some sort of ack scheme.
  ARMCIX_Fence (proc);
  ARMCIX_Wait (&nb_handle->cmpl_info);
#else
  DCMF_CriticalSection_enter (0);

  //fprintf (stderr, "ARMCIX_AccS() >> \n");
  //fprintf (stderr, "ARMCIX_AccS() -- __connection[%d].sequence.origin=%d, __connection[%d].active=%d, __global_connection.active=%d\n", proc, __connection[proc].sequence.origin, proc, __connection[proc].active, __global_connection.active);

  // Calculate the number of requests
  unsigned i;
  unsigned n = 1;
  for (i = 0; i < stride_levels; i++) n = n * seg_count[i+1];

  armci_ireq_t nb_handle;
  armcix_dcmf_opaque_t * dcmf = (armcix_dcmf_opaque_t *) &nb_handle.cmpl_info;
  dcmf->connection = &__connection[proc];
  dcmf->active = n;

  __connection[proc].active += n;
  __global_connection.active += n;

  unsigned count;
  count = ARMCIX_DCMF_AccS_recurse (datatype, scale,
                                    src_ptr, src_stride_arr,
                                    dst_ptr, dst_stride_arr,
                                    seg_count, stride_levels, proc,
                                    (armci_ihdl_t) &nb_handle);

#ifdef BLOCKING_OPERATIONS_REQUIRE_FENCE
  ARMCIX_Fence (proc);
#else
  assert (n == count);
  while (dcmf->active) DCMF_Messager_advance ();
#endif

  //fprintf (stderr, "ARMCIX_AccS() << \n");

  DCMF_CriticalSection_exit  (0);
#endif
  return 0;
}

/**
 * \brief ARMCI Extension non-blocking strided accumulate operation.
 *
 * \param[in] datatype       accumulate datatype (operation code)
 * \param[in] scale          opaque pointer to the scaling factor for accumulate
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
int ARMCIX_NbAccS (int datatype, void * scale,
                   void * src_ptr, int * src_stride_arr, 
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
  dcmf->active += n;

  __connection[proc].active += n;
  __global_connection.active += n;

  unsigned count;
  count = ARMCIX_DCMF_AccS_recurse (datatype, scale,
                                    src_ptr, src_stride_arr, 
                                    dst_ptr, dst_stride_arr, 
                                    seg_count, stride_levels, proc,
                                    nb_handle);

  assert (n == count);

  DCMF_CriticalSection_exit  (0);

  return 0;
}
