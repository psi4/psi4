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
 * \file armci/src/armcix/dcmf/armcix_impl.c
 * \brief DCMF ARMCI Extension implementation.
 */

#include "armcix_impl.h"
#if HAVE_STRINGS_H
#   include <strings.h>
#endif

ARMCIX_DCMF_Connection_t __global_connection;
ARMCIX_DCMF_Connection_t * __connection;
DCMF_Memregion_t __local_mem_region;

volatile unsigned __memregions_to_receive;

typedef struct ARMCIX_DCMF_RequestInfo_t
{
  ARMCIX_DCMF_Request_t              info;
  DCMF_Callback_t                    cb_free;
  struct ARMCIX_DCMF_RequestInfo_t * next;
  unsigned                           unused;
} ARMCIX_DCMF_RequestInfo_t __attribute__ ((__aligned__ (16)));

typedef struct ARMCIX_DCMF_RequestPool_t
{
  ARMCIX_DCMF_RequestInfo_t * head;
  unsigned max;
  unsigned current;
  unsigned increment;
} ARMCIX_DCMF_RequestPool_t;


ARMCIX_DCMF_RequestPool_t __armcix_dcmf_requestpool;

void ARMCIX_DCMF_request_print (char * label)
{
  char str[1024];
  if (label == NULL) str[0] = 0;
  else snprintf (str, 1024, "[%s] ", label);

  fprintf (stderr, "%s__armcix_dcmf_requestpool { head = %p, max = %d, current = %d, increment = %d }\n", str, __armcix_dcmf_requestpool.head, __armcix_dcmf_requestpool.max, __armcix_dcmf_requestpool.current, __armcix_dcmf_requestpool.increment);

  ARMCIX_DCMF_RequestInfo_t * p = __armcix_dcmf_requestpool.head;
  while (p != NULL)
  {
    fprintf (stderr, "    (%p)->next = %p\n", p, p->next);
    p = p->next;
  }
}

void ARMCIX_DCMF_request_initialize (unsigned max, unsigned increment)
{
  unsigned count = max;
  if (increment > 0 && increment < max) count = increment;

  __armcix_dcmf_requestpool.head = (ARMCIX_DCMF_RequestInfo_t *) malloc (sizeof(ARMCIX_DCMF_RequestInfo_t) * count);
  assert (__armcix_dcmf_requestpool.head!=NULL);

  __armcix_dcmf_requestpool.max = max;
  __armcix_dcmf_requestpool.current = count;
  __armcix_dcmf_requestpool.increment = increment;

  unsigned i;
  for (i=1; i<count; i++) __armcix_dcmf_requestpool.head[i-1].next = & __armcix_dcmf_requestpool.head[i];
  __armcix_dcmf_requestpool.head[count-1].next = NULL;

  //ARMCIX_DCMF_request_print ("init");
}

ARMCIX_DCMF_Request_t * ARMCIX_DCMF_request_allocate (DCMF_Callback_t cb_free)
{
  //ARMCIX_DCMF_request_print ("allocate");

  if (__armcix_dcmf_requestpool.head == NULL)
  {
    if (__armcix_dcmf_requestpool.current < __armcix_dcmf_requestpool.max)
    {
      // Allocate a new block of request objects and add them to the request pool.
      __armcix_dcmf_requestpool.head = 
        (ARMCIX_DCMF_RequestInfo_t *) malloc (sizeof(ARMCIX_DCMF_RequestInfo_t) * __armcix_dcmf_requestpool.increment);
      assert (__armcix_dcmf_requestpool.head!=NULL);

      __armcix_dcmf_requestpool.current += __armcix_dcmf_requestpool.increment;
      unsigned i;
      for (i=1; i<__armcix_dcmf_requestpool.increment; i++)
        __armcix_dcmf_requestpool.head[i-1].next = & __armcix_dcmf_requestpool.head[i];
      __armcix_dcmf_requestpool.head[__armcix_dcmf_requestpool.increment-1].next = NULL;
      //fprintf (stderr, "ARMCIX_DCMF_request_allocate() .. allocate a new block of requests (current = %d -> %d)\n", previous, __armcix_dcmf_requestpool.current);
    }
    else
    {
      // The request pool has already reached its maximum size, advance until a request is freed.
      do
      {
        DCMF_Messager_advance ();
      } while (__armcix_dcmf_requestpool.head == NULL);
    }
  }

  // Get the next free request object from the request pool, and set the
  // request pool pointer to the next available request object.
  ARMCIX_DCMF_RequestInfo_t * _request = (ARMCIX_DCMF_RequestInfo_t *) __armcix_dcmf_requestpool.head;
  __armcix_dcmf_requestpool.head = _request->next;

  // Initialize the new request object before return
  _request->cb_free = cb_free;
  _request->next = NULL;

  return (ARMCIX_DCMF_Request_t *) _request;
}

void ARMCIX_DCMF_request_free (ARMCIX_DCMF_Request_t * request)
{
  ARMCIX_DCMF_RequestInfo_t * _request = (ARMCIX_DCMF_RequestInfo_t *) request;

  // Invoke the "free" callback if it is specified.
  if (_request->cb_free.function != NULL)
    _request->cb_free.function (_request->cb_free.clientdata, NULL);

  // Return the request to the free request pool.
  _request->next = __armcix_dcmf_requestpool.head;
  __armcix_dcmf_requestpool.head = _request;
}

/**
 * \brief Generic decrement callback
 *
 * \param[in] clientdata Address of the variable to decrement
 */
void ARMCIX_DCMF_cb_decrement (void * clientdata, DCMF_Error_t *err)
{
  unsigned * value = (unsigned *) clientdata;
  (*value)--;
}

/**
 * \brief Callback function for non-blocking operations
 *
 * \param[in] clientdata The non-blocking handle to complete
 */
void ARMCIX_DCMF_NbOp_cb_done (void * clientdata, DCMF_Error_t *err)
{
  armci_ihdl_t nb_handle = (armci_ihdl_t) clientdata;
  armcix_dcmf_opaque_t * dcmf = (armcix_dcmf_opaque_t *) &nb_handle->cmpl_info;
  dcmf->active--;
  dcmf->connection->active--;
  __global_connection.active--;
}

/**
 * \brief DCMF ARMCI Extention receive memory region short callback
 *
 * \see DCMF_RecvSend
 */
void ARMCIX_DCMF_RecvMemregion1 (void           * clientdata,
                                 const DCQuad   * msginfo,
                                 unsigned         count,
                                 unsigned         peer,
                                 const char     * src,
                                 unsigned         bytes)
{
  ARMCIX_DCMF_Connection_t * connection = (ARMCIX_DCMF_Connection_t *) clientdata;
  memcpy (&connection[peer].remote_mem_region, src, bytes);
  __memregions_to_receive--;
}


/**
 * \brief DCMF ARMCI Extention receive memory region long callback
 *
 * \see DCMF_RecvSend
 */
DCMF_Request_t * ARMCIX_DCMF_RecvMemregion2 (void             * clientdata,
                                             const DCQuad     * msginfo,
                                             unsigned           count,
                                             unsigned           peer,
                                             unsigned           sndlen,
                                             unsigned         * rcvlen,
                                             char            ** rcvbuf,
                                             DCMF_Callback_t  * cb_done)
{
  assert(0);
  ARMCIX_DCMF_Connection_t * connection = (ARMCIX_DCMF_Connection_t *) clientdata;

  *rcvlen = sndlen;
  *rcvbuf = (char *) &connection[peer].remote_mem_region;

  cb_done->function   = (void (*)(void *, DCMF_Error_t *))free; // still works, for now.
  cb_done->clientdata = (void *) malloc (sizeof (DCMF_Request_t));

  return cb_done->clientdata;
}

void ARMCIX_DCMF_Connection_initialize ()
{
  DCMF_CriticalSection_enter(0);

  __global_connection.peer = (unsigned) -1;

  unsigned rank = DCMF_Messager_rank ();
  unsigned size = DCMF_Messager_size ();
  posix_memalign ((void **)&__connection, 16, sizeof(ARMCIX_DCMF_Connection_t) * size);
  bzero ((void *)__connection, sizeof(ARMCIX_DCMF_Connection_t) * size);

  void * base  = NULL;
  size_t bytes = (size_t) -1;

  unsigned i;
  for (i = 0; i < size; i++)
  {
    __connection[i].peer = i;
#warning fix memregion setup to handle non-global address space pinning.
    //DCMF_Result result =
      DCMF_Memregion_create (&__connection[i].local_mem_region,
                             &bytes, (size_t) -1, NULL, 0);
  }

  // Register a send protocol to exchange memory regions
  DCMF_Protocol_t send_protocol;
  DCMF_Send_Configuration_t send_configuration = {
    DCMF_DEFAULT_SEND_PROTOCOL,
    DCMF_DefaultNetwork,
    ARMCIX_DCMF_RecvMemregion1,
    __connection,
    ARMCIX_DCMF_RecvMemregion2,
    __connection
  };
  DCMF_Send_register (&send_protocol, &send_configuration);

  DCMF_Request_t request;
  volatile unsigned active;
  DCMF_Callback_t cb_done = { ARMCIX_DCMF_cb_decrement, (void *) &active };

  // Exchange the memory regions
  __memregions_to_receive = size;
  for (i = 0; i < size; i++)
  {
    unsigned peer = (rank+i)%size;
    active = 1;
    DCMF_Send (&send_protocol,
               &request,
               cb_done,
               DCMF_SEQUENTIAL_CONSISTENCY,
               peer,
               sizeof(DCMF_Memregion_t),
               (char *) &__connection[peer].local_mem_region,
               (DCQuad *) NULL,
               0);
    while (active) DCMF_Messager_advance();
  }
  while (__memregions_to_receive) DCMF_Messager_advance();

  DCMF_CriticalSection_exit(0);
}


static inline int
ENV_Bool(char * env, int * dval)
{
  int result = *dval;
  if(env != NULL)
    {
      if (strcmp(env, "0") == 0)
        result = 0;
      else if  (strcmp(env, "0") == 1)
        result = 1;
    }
  return *dval = result;
}

static inline int
ENV_Int(char * env, int * dval)
{
  int result = *dval;
  if(env != NULL)
    {
      result = (int) strtol((const char *)env, NULL, 10);
    }
  return *dval = result;
}

/**
 * \brief Initialize the DCMF ARMCI resources
 */
int ARMCIX_Init ()
{
  DCMF_CriticalSection_enter(0);

  DCMF_Messager_initialize ();

  ARMCIX_DCMF_Connection_initialize ();

  /* Determine request pool defaults */
  int ARMCIX_DCMF_REQUESTPOOL_MAX = 1000;
  ENV_Int (getenv ("ARMCIX_DCMF_REQUESTPOOL_MAX"), &ARMCIX_DCMF_REQUESTPOOL_MAX);
  int ARMCIX_DCMF_REQUESTPOOL_INC = 0;
  ENV_Int (getenv ("ARMCIX_DCMF_REQUESTPOOL_INC"), &ARMCIX_DCMF_REQUESTPOOL_INC);
  ARMCIX_DCMF_request_initialize (ARMCIX_DCMF_REQUESTPOOL_MAX, ARMCIX_DCMF_REQUESTPOOL_INC);



  ARMCIX_DCMF_Get_register ();

  ARMCIX_DCMF_Put_register (__connection);

  ARMCIX_DCMF_Acc_register (__connection);

  ARMCIX_DCMF_Fence_register (__connection);

  ARMCIX_DCMF_Rmw_register ();

  /* Determine interrupt mode */
  int interrupts = 1;
  ENV_Bool (getenv ("DCMF_INTERRUPT"),  &interrupts);
  ENV_Bool (getenv ("DCMF_INTERRUPTS"), &interrupts);

  DCMF_Configure_t config;
  memset (&config, 0x00, sizeof(DCMF_Configure_t));
  config.interrupts = (interrupts==0)?DCMF_INTERRUPTS_OFF:DCMF_INTERRUPTS_ON;
  DCMF_Messager_configure (&config, &config);

  DCMF_Messager_configure (NULL, &config);

  //ARMCIX_DCMF_request_print ("after armcix_init");

  DCMF_CriticalSection_exit(0);

  return 0;
}
