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
 * \file armci/src/armcix/dcmf/armcix_lock.c
 * \brief DCMF ARMCI Extension for lock operations.
 */

#include "armcix_impl.h"

typedef struct ARMCIX_DCMF_Lockwaiter_t
{
  unsigned   peer;
  void     * start;
  void     * end;
  unsigned * update;
}
ARMCIX_DCMF_Lockwaiter_t;

typedef struct ARMCIX_DCMF_Lockwaitq_t
{
  ARMCIX_DCMF_Lockwaiter_t * queue;
  unsigned                   size;
  unsigned                   head;
  unsigned                   tail;
}
ARMCIX_DCMF_Lockwaitq_t;

typedef enum
{
  ARMCIX_DCMF_Lock_Request,
  ARMCIX_DCMF_Lock_Ack,
  ARMCIX_DCMF_Lock_Release
}
ARMCIX_DCMF_Lockinfo_operation_t;

typedef struct ARMCIX_DCMF_Lockinfo_t
{
  ARMCIX_DCMF_Lockinfo_operation_t op;
  unsigned     * update;
  union
  {
    struct
    {
      void     * start;
      void     * end;
    } request;
    struct
    {
      unsigned   slot;
      unsigned   unused;
    } key;
  };
}
ARMCIX_DCMF_Lockinfo_t __attribute__ ((__aligned__ (16)));

/** Lock operation DCMF protocol global object. */
DCMF_Protocol_t         __lock_req_protocol;
DCMF_Protocol_t         __lock_ack_protocol;

/**
 * \brief Origin node active lock slots.
 */
#if (MAX_SLOTS > 256)
#error MAX_SLOTS must be less than or equal to 256
#endif
unsigned char              * __lock_slot;

/** Target node pending lock queue. */
ARMCIX_DCMF_Lockwaitq_t __lock_pending;

/**
 * \brief Acquire a lock on a memory address range.
 *
 * Inspects the memory lock table and, if there is an open slot and the
 * requested memory address range does not conflict with an existing lock,
 * updates the memory locak table with the lock information and returns
 * the slot of the acquired lock.
 *
 * \todo Can this be moved up to a more platform-neutral place?
 *
 * \param[in] local_memlock_table The memory lock table
 * \param[in] pstart              The start of the memory address range
 * \param[in] pend                The end of the memory address range
 *
 * \return    The memory lock table slot of the successfully acquired lock, 
 *            otherwise \c -1
 *
 * \see ARMCIX_release_lock
 */
int ARMCIX_aquire_lock (memlock_t * local_memlock_table, void * pstart, void * pend)
{
  int acquired_lock = -1;

  /* inspect the table */
  unsigned conflict = 0; 
  unsigned slot = 0;
  for(; slot < MAX_SLOTS; slot++)
    {
      /* nonzero starting address means the slot is occupied */ 
      if(local_memlock_table[slot].start == NULL)
        {
          /* remember a free slot to store address range */
          acquired_lock = slot;
        }
      else
        {
          /* check for conflict: overlap between stored and current range */
          if ( (pstart >= local_memlock_table[slot].start &&
                pstart <= local_memlock_table[slot].end) ||
               (pend >= local_memlock_table[slot].start &&
                pend <= local_memlock_table[slot].end) )
            {
              conflict = 1;
              break;
            }
        }
    }

  if (acquired_lock != -1 && !conflict)
    {
      /* acquired the memory lock: enter address into the table */
      local_memlock_table[acquired_lock].start = pstart;
      local_memlock_table[acquired_lock].end = pend;
    }
  else
    {
      acquired_lock = -1;
    }

  return acquired_lock;
}

/**
 * \brief Release a lock of the memory address range in the specified lock slot.
 *
 * Clears the memory address range in specified lock slot of the memory lock
 * table.
 *
 * \todo Can this be moved up to a more platform-neutral place?
 *
 * \param[in] local_memlock_table The memory lock table
 * \param[in] slot                The memory lock slot to release
 *
 * \see ARMCIX_aquire_lock
 */
void ARMCIX_release_lock (memlock_t * local_memlock_table, unsigned slot)
{
  local_memlock_table[slot].start = NULL;
  local_memlock_table[slot].end   = NULL;
}

/**
 * \brief Receive a lock control message.
 *
 * The lock message type is either a lock \e request, \e acknowledgement, 
 * or \e release.
 *
 * For a lock request message the local memlock table is inspected for an 
 * available lock slot.  If a slot is found a lock acknowledgement message
 * is immediately sent to the peer node. Otherwise, the lock request
 * information is saved in a pending queue.
 *
 * The lock acknowledgment message is sent from the target node when the
 * origin node has acquired a lock of a resource on the target node. The
 * update variable is unset which will end the polling advance operation
 * inside the blocking ARMCIX_Lockmem() function.
 *
 * The lock release message frees the specified slot and attempts to acquire
 * a lock for the first lock waiter in the pending queue. If the lock is 
 * aqcuired then a lock acknowledgement is sent to the waiter.
 *
 * \param[in] clientdata Registered clientdata, the armci connection array
 * \param[in] info       Lock control information
 * \param[in] peer       Rank of the node that sent this control message
 *
 * \see DCMF_RecvControl
 */
void ARMCIX_DCMF_RecvLockMessage (void                 * clientdata,
                                  const DCMF_Control_t * info,
                                  unsigned               peer)
{
  /* Get the lock information sent with this lock message                    */
  ARMCIX_DCMF_Lockinfo_t * lockinfo = (ARMCIX_DCMF_Lockinfo_t *) info;

  switch (lockinfo->op)
  {
    case ARMCIX_DCMF_Lock_Request:
      {
        unsigned queue_request = 1;
        if (__lock_pending.head == __lock_pending.tail)
        {
          /* There are no pending lock requests. Attempt to acquire the      */
          /* request lock.                                                   */
          int lock = ARMCIX_aquire_lock ((memlock_t *) clientdata, lockinfo->request.start, lockinfo->request.end);
          if (lock != -1)
          {
            /* The lock was acquired. Send a lock acknowledgement message to */
            /* the origin node with the slot of the aquired lock and the     */
            /* address of the update variable on the origin node.            */
            ARMCIX_DCMF_Lockinfo_t ack;
            ack.op             = ARMCIX_DCMF_Lock_Ack;
            ack.update         = lockinfo->update;
            ack.key.slot       = lock;
            ack.key.unused     = 0;

            DCMF_Control (&__lock_ack_protocol,
                          DCMF_SEQUENTIAL_CONSISTENCY,
                          peer,
                         (DCMF_Control_t *) &ack);

            queue_request = 0;
          }
        }

        if (queue_request)
        {
          /* There were pending lock requests, or the acquire lock attempt   */
          /* failed. Queue the lock request until another node releases a    */
          /* lock.                                                           */
          ARMCIX_DCMF_Lockwaiter_t * tail = &__lock_pending.queue[__lock_pending.tail];
          tail->peer   = peer;
          tail->start  = lockinfo->request.start;
          tail->end    = lockinfo->request.end;
          tail->update = lockinfo->update;

          /* Advance the tail index                                            */
          __lock_pending.tail = (__lock_pending.tail+1)%__lock_pending.size;
        }
      }
      break;

    case ARMCIX_DCMF_Lock_Ack:

      /* Save the lock slot until the lock is released.                      */
      __lock_slot[peer] = lockinfo->key.slot;

      /* Set the update variable to zero. This will break the polling loop   */
      /* in ARMCIX_Lockmem() and allow that function to return.              */
      *(lockinfo->update) = 0;
      break;

    case ARMCIX_DCMF_Lock_Release:
      {
        /* Release the lock in this slot.                                    */
        ARMCIX_release_lock ((memlock_t *) clientdata, lockinfo->key.slot);

        if (__lock_pending.head != __lock_pending.tail)
        {
          /* There is a pending lock request. Attempt to acquire the lock as */
          /* specified in the pending queue now that the previous lock has   */
          /* been released.                                                  */
          ARMCIX_DCMF_Lockwaiter_t * head = &__lock_pending.queue[__lock_pending.head];
          int lock = ARMCIX_aquire_lock ((memlock_t *) clientdata, head->start, head->end);
          if (lock != -1)
          {
            /* The aquire lock attempt was successful. Send a lock           */
            /* acknowledgement message to the origin node of the pending     */
            /* lock request with the slot of the aquired lock and the        */
            /* address of the update variable on the origin node.            */
            ARMCIX_DCMF_Lockinfo_t ack;
            ack.op             = ARMCIX_DCMF_Lock_Ack;
            ack.update         = head->update;
            ack.key.slot       = lock;
            ack.key.unused     = 0;

            DCMF_Control (&__lock_ack_protocol,
                          DCMF_SEQUENTIAL_CONSISTENCY,
                          head->peer,
                          (DCMF_Control_t *) &ack);

           /* Advance the head index (dequeue the pending lock request)      */
            __lock_pending.head = (__lock_pending.head+1)%__lock_pending.size;
          }
        }
      }
      break;
    default:
      /* Bad! Bad! */
      assert (0);
      break;
  }
}


/**
 * \brief DCMF ARMCI Extention receive short lock request callback
 *
 * \see ARMCIX_DCMF_RecvLockAck
 * \see DCMF_RecvSendShort
 */
void ARMCIX_DCMF_RecvLockRequest (void           * clientdata,
                                  const DCQuad   * msginfo,
                                  unsigned         count,
                                  unsigned         peer,
                                  const char     * src,
                                  unsigned         bytes)
{
  ARMCIX_DCMF_RecvLockMessage (clientdata, (const DCMF_Control_t *) msginfo, peer);
}


/**
 * \brief Initialize the ARMCI Extention lock resources.
 *
 * Register the DCMF Control protocol used to pass lock messages between nodes.
 * Allocate a lock pending queue for use when a lock request is received by
 * the receive callback function ARMCIX_DCMF_RecvLockMessage() and the resource
 * is currently allocated to another node.
 *
 * \param[in]  local_memlock_table memlock table
 *
 * \see ARMCIX_DCMF_Lockwaitq_t
 * \see ARMCIX_DCMF_Lockwaiter_t
 * \see ARMCIX_DCMF_RecvLockMessage
 * \see DCMF_Control_register
 */
void ARMCIX_init_memlock (memlock_t * local_memlock_table)
{
  DCMF_CriticalSection_enter (0);

  DCMF_Send_Configuration_t send_configuration = {
    DCMF_DEFAULT_SEND_PROTOCOL,
    DCMF_DEFAULT_NETWORK,
    ARMCIX_DCMF_RecvLockRequest,
    local_memlock_table,
    NULL,
    NULL
  };
  DCMF_Send_register (&__lock_req_protocol, &send_configuration);

  DCMF_Control_Configuration_t ctrl_configuration = {
    DCMF_DEFAULT_CONTROL_PROTOCOL,
    DCMF_DEFAULT_NETWORK,
    ARMCIX_DCMF_RecvLockMessage,
    local_memlock_table
  };
  DCMF_Control_register (&__lock_ack_protocol, &ctrl_configuration);

  unsigned msize = DCMF_Messager_size ();
  unsigned qsize = sizeof(ARMCIX_DCMF_Lockwaiter_t) * msize+1;
  __lock_pending.size = msize;
  __lock_pending.queue = (ARMCIX_DCMF_Lockwaiter_t *) malloc (qsize);
  __lock_pending.head = 0;
  __lock_pending.tail = 0;

  __lock_slot = (unsigned char *) malloc (sizeof(unsigned char) * msize);
  memset(__lock_slot, 0x00, sizeof(unsigned char) * msize);

  DCMF_CriticalSection_exit (0);
}

/**
 * \brief ARMCI Extension blocking memory lock operation.
 *
 * Send a lock request to the remote node and block until the lock has been
 * acquired on the remote node.
 *
 * \param[in] pstart The start virtual address of the range of memory to lock.
 * \param[in] pend   The end virtual address of the range of memory to lock.
 * \param[in] proc   Remote process(or) ID
 *
 * \see ARMCIX_DCMF_Lockinfo_t
 * \see ARMCIX_DCMF_RecvLockMessage
 * \see DCMF_Control
 */
void ARMCIX_Lockmem (void * pstart, void * pend, int proc)
{
  DCMF_CriticalSection_enter (0);

  volatile unsigned active = 1;

  ARMCIX_DCMF_Lockinfo_t info;
  info.op             = ARMCIX_DCMF_Lock_Request;
  info.update         = (unsigned *)&active;
  info.request.start  = pstart;
  info.request.end    = pend;

  DCMF_Request_t request;
  DCMF_Send ( &__lock_req_protocol,
              &request,
              (DCMF_Callback_t) { NULL, NULL },
              DCMF_SEQUENTIAL_CONSISTENCY,
              proc,
              0,
              NULL,
              (DCQuad *) &info,
              sizeof(ARMCIX_DCMF_Lockinfo_t)/sizeof(DCQuad));

  while (active) DCMF_Messager_advance ();

  DCMF_CriticalSection_exit  (0);
}


/**
 * \brief ARMCI Extension release memory lock operation.
 *
 * Send a lock release message to the remote node. This is a \e fire-and-forget
 * operation because the node does not block for an acknowledgement that the
 * lock release was successful.
 *
 * \param[in] proc   Remote process(or) ID
 *
 * \see ARMCIX_DCMF_Lockinfo_t
 * \see ARMCIX_DCMF_RecvLockMessage
 * \see DCMF_Control
 */
void ARMCIX_Unlockmem (int proc)
{
  DCMF_CriticalSection_enter (0);

  ARMCIX_DCMF_Lockinfo_t info;
  info.op              = ARMCIX_DCMF_Lock_Release;
  info.key.slot        = (unsigned) __lock_slot[proc];
  info.key.unused      = 0;

  DCMF_Request_t request;
  volatile unsigned active = 1;
  DCMF_Send ( &__lock_req_protocol,
              &request,
              (DCMF_Callback_t) { ARMCIX_DCMF_cb_decrement, (void *)&active },
              DCMF_SEQUENTIAL_CONSISTENCY,
              proc,
              0,
              NULL,
              (DCQuad *) &info,
              sizeof(ARMCIX_DCMF_Lockinfo_t)/sizeof(DCQuad));

  while (active) DCMF_Messager_advance ();

  DCMF_CriticalSection_exit  (0);
}
