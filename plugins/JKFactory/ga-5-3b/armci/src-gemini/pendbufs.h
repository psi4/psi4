/** @file Split buffer implementation. 
 * @author Sriram Krishnamoorthy
 * 
 * Supports multiple short/immediate buffers posted per client and a
 * client-independent number of buffers to handle large messages. 
 */
#ifndef _PENDBUFS_H_
#define _PENDBUFS_H_

#if defined(PEND_BUFS)

#include "armcip.h"
#include "request.h"


/**The buf should be the first field in immbuf_t and pendbuf_t. For
   example, look at openib.c:armci_rcv_req and maybe other places*/
typedef struct immbuf_t {
  char *buf; /*immediate buffer[IMMBUF_LEN]*/
/*   IMMBUF_NW_T fields; */
  IMMBUF_NW_T
  struct immbuf_t *immbuf_list_next;
} immbuf_t;

typedef struct pendbuf_t {
  char *buf; /*pending buffer[PENDBUF_LEN]*/
/*   PENDBUF_NW_T fields; */
  PENDBUF_NW_T
  int status; /*<Finite State Machine status*/
  int avail;  /*<Buffer is available*/
  immbuf_t *vbuf; /*<Immediate buffer corresponding to this pending buffer*/
  struct pendbuf_t *order_next,*order_prev; /*Order among pending buffers processing requests from same client*/
  int commit_me; /*This buffer needs to be committed (by calling
		   progresson it)*/
} pendbuf_t;

int get_immbuf_len(); 
int get_immbuf_num();
int get_pendbuf_len(); 
int get_pendbuf_num();
#define IMM_BUF_LEN get_immbuf_len()
#define PENDING_BUF_LEN get_pendbuf_len()
#define IMM_BUF_NUM get_immbuf_num()
#define PENDING_BUF_NUM get_pendbuf_num()

/*Initialize number and length of pending and initialized
 *  buffers. Called when initializing the network connections
 */
void armci_pbuf_init_buffer_env();

/*----------------Stuff to enforce server-side ordering----------------*/

enum PendBufOrderingRule {
  ONE_PBUF_MESG,         /**<Only one pending mesg progressed at any time*/
  ONE_PBUF_MESG_PER_PROC,/**<One pending mesg per client proc*/
  ACC_NO_ORDER,          /**<Consecutive ACCs from same client proc
			    can overlap*/
  PUTACC_SPLIT_ORDER     /**<ACC_NO_ORDER+concurrent get of data from
			    client for non-immediate PUTs*/
};


extern pendbuf_t *serv_pendbuf_arr;   /**<Array of pending buffers*/

void armci_pendbuf_service_req(immbuf_t *immbuf);
void armci_pendbuf_init();
void armci_pendbuf_done_put(int pbufid);
void armci_pendbuf_done_get(int pbufid);

#endif /*PEND_BUFS*/
#endif /*_PENDBUFS_H_*/

