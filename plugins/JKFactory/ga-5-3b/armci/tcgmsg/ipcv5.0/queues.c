#if HAVE_CONFIG_H
#   include "config.h"
#endif

extern void free(void *ptr);

#include "tcgmsgP.h"

static const long false = 0;
static const long true  = 1;

extern void Busy(int);

extern long async_send(SendQEntry *);


/**
 * Given a nodeid return a unqiue integer constructed by
 * combining it with the value of a counter
 */
static long NextMsgID(long node)
{
    static long id = 0;
    static long mask = (1<<20)-1;

    id = (id + 1) & mask;
    if (id == 0) id = 1;

    return (node << 20) + id;
}


/**
 * Given an id from NextMsgID extract the node
 */
static long NodeFromMsgID(long msgid)
{
    long node = msgid >> 20;

    if (node < 0 || node > NNODES_())
        Error("NodeFromMsgID: invalid msgid", msgid);

    return node;
}


/**
 * Flush as many messages as possible without blocking from
 * the send q to the specified node.
 */
static void flush_send_q_node(long node)
{
    while (TCGMSG_proc_info[node].sendq) {

        if (!async_send(TCGMSG_proc_info[node].sendq)) {
            /* Send is incomplete ... stop processing this q*/
            break;
        }
        else {
            SendQEntry *tmp = TCGMSG_proc_info[node].sendq;

            TCGMSG_proc_info[node].sendq = (SendQEntry *) TCGMSG_proc_info[node].sendq->next;
            if (tmp->free_buf_on_completion)
                (void) free(tmp->buf);
            tmp->active = false;    /* Matches NewSendQEntry() */
        }
    }
}


/**
 * Flush as many messages as possible without blocking
 * from all of the send q's.
 */
void flush_send_q()
{
    long node;
    long nproc = NNODES_();

    for (node=0; node<nproc; node++)
        if (TCGMSG_proc_info[node].sendq)
            flush_send_q_node(node);
}    


/**
 * Return 0 if the message operation is incomplete.
 * Return 1 if the message operation is complete.
 */
long msg_status(long msgid)
{
    long node = NodeFromMsgID(msgid);
    SendQEntry *entry;
    long status = 1;

    flush_send_q();

    /* Attempt to find the msgid in the message q.  If it is not
       there then the send is complete */

    for (entry=TCGMSG_proc_info[node].sendq; entry; entry=(SendQEntry *) entry->next) {
        if (entry->msgid == msgid) {
            status = 0;
            break;
        }
    }

    return status;
}


/**
 * Wait for the operation referred to by msgid to complete.
 */
void msg_wait(long msgid)
{
    long nspin = 0;
    long spinlim = 1000000;

    while (!msg_status(msgid)) {
        nspin++;
        if (nspin < spinlim)
            Busy(100);
        else 
            usleep(1);
    }
}


static SendQEntry *NewSendQEntry(void)
{
    SendQEntry *new = TCGMSG_sendq_ring;

    if (new->active)
        Error("NewSendQEntry: too many outstanding sends\n", 0L);

    TCGMSG_sendq_ring = (SendQEntry *) TCGMSG_sendq_ring->next_in_ring;

    new->active = true;

    return new;
}


long msg_async_snd(long type, char *buf, long lenbuf, long node)
{
    long msgid;
    SendQEntry *entry;

    if (node<0 || node>=TCGMSG_nnodes)
        Error("msg_async_send: node is out of range", node);

    if (node == TCGMSG_nodeid)
        Error("msg_async_send: cannot send to self", node);

    msgid = NextMsgID(node);
    entry = NewSendQEntry();

    /* Insert a new entry into the q */

    entry->tag   = TCGMSG_proc_info[node].n_snd++; /* Increment tag */
    entry->msgid = msgid;
    entry->type  = type;
    entry->buf   = buf;
    entry->free_buf_on_completion = 0;
    entry->lenbuf= lenbuf;
    entry->node  = node;
    entry->next  = (SendQEntry *) 0;
    entry->written = 0;
    entry->buffer_number = 0;

    /* Attach to the send q */

    if (!TCGMSG_proc_info[node].sendq)
        TCGMSG_proc_info[node].sendq = entry;
    else {
        SendQEntry *cur = TCGMSG_proc_info[node].sendq;

        while (cur->next)
            cur = cur->next;
        cur->next = entry;
    }

    /* Attempt to flush the send q */

    flush_send_q();

    return msgid;
}


/**
 * synchronous send of message to a process
 *
 * long *type     = user defined integer message type (input)
 * char *buf      = data buffer (input)
 * long *lenbuf   = length of buffer in bytes (input)
 * long *node     = node to send to (input)
 *
 * for zero length messages only the header is sent
 */
void msg_snd(long type, char *buf, long lenbuf, long node)
{
    msg_wait(msg_async_snd(type, buf, lenbuf, node));
}
