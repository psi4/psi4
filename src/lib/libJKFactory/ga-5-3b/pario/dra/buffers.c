#if HAVE_CONFIG_H
#   include "config.h"
#endif

/** @file
 * Buffer Manager for managing buffers in any application
 */

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#if HAVE_STRINGS_H
#   include <strings.h>
#endif

#include "buffers.h"

/*#define STATBUF 1 */ /* set if static buffers will be used */
/*#define DEBUG 1 */

#define MAX_CTXT 5

int ctxt_count = 0;


/** buffer management initialization routine */
void buffer_init(
        buf_context_t *ctxt, int nbuf, int buf_size, void (*fptr)(char*))
{
    int i;
    long diff;
    if (nbuf < 1 || nbuf > MAXBUF) {
        printf("Too many (or too few) buffers requested, using default number (%d) of buffers", DEFBUF);
        nbuf = DEFBUF;
    }

    /* create a context */
    ctxt->ctxt_id = ctxt_count++;
    if (ctxt->ctxt_id > MAX_CTXT) {
        printf("Max number of contexts reached!\n");
        return;
    }

    ctxt->nbuf = nbuf;
    ctxt->size = buf_size;
    ctxt->buf = (_buffer_t *)malloc(sizeof(_buffer_t) * nbuf);
    ctxt->fptr = fptr;

#ifdef STATBUF
    double buffers[MAXBUF][DBL_BUF_SIZE];
    for (i = 0; i < nbuf; i++) { 
        ctxt->buf[i].buffer = (char *) buffers[i];
        bzero(buffers[i], sizeof(buffers[i]));
    }
#else /* STATBUF */

    /* get buffer memory */
    for (i = 0; i < nbuf; i++) {
        ctxt->buf[i].buffer = (char*) malloc((buf_size + ALIGN-1) *sizeof(double));

        if (ctxt->buf[i].buffer == NULL) {
            printf("Could not allocate memory for buffers!\n");
            return;
        }
        bzero(ctxt->buf[i].buffer, sizeof(ctxt->buf[i].buffer));

        /* align buffer address */
        diff = ((long)(ctxt->buf[i].buffer)) % (sizeof(double)*ALIGN);
        if(diff) 
        {
            ctxt->buf[i].align_off = (int) (sizeof(double)*ALIGN - diff);
        } else 
        {
            ctxt->buf[i].align_off = (int) 0;
        }
        ctxt->buf[i].buffer += ctxt->buf[i].align_off;

    }
#endif /* STATBUF */
    for (i = 0; i < nbuf; i++) {
        ctxt->buf[i].active = 0;
        ctxt->buf[i].group_id = 0;
    }
#ifdef DEBUG
    printf("Created a context\n\n");
#endif
}


/** internal function to return empty buffer handle */
int get_buf_hdl(buf_context_t *ctxt)
{
    int i;
    for (i = 0; i < ctxt->nbuf; i++) {
        if (ctxt->buf[i].active == 0) {
            ctxt->buf[i].active = 1;
            return i;
        }
    }
    return -1;
}


char* get_buf(buf_context_t *ctxt, int call_id)
{
    int hdl;
    char *buf;

    hdl = get_buf_hdl(ctxt);
    if (hdl == -1) {
        int cur_buf;
        /* no buffer is available, wait for the oldest buffer to become free */
        cur_buf = (ctxt->last_buf + 1) % ctxt->nbuf;

        /* free this buffer by calling the callback function provided by the application */
        ctxt->fptr(ctxt->buf[cur_buf].buffer); 
        hdl = cur_buf;
    }
    buf = ctxt->buf[hdl].buffer;
    ctxt->buf[hdl].buf_hdl = hdl;
    ctxt->buf[hdl].call_id = call_id;
    ctxt->buf[hdl].active = 1;
    ctxt->last_buf = hdl;
#ifdef DEBUG
    printf("Giving a buffer with internal handle: %d\n", hdl);
#endif
    return (buf);
}


/** function to free a buffer */
void free_buf(buf_context_t *ctxt, char *buf)
{
    int i;
    for (i = 0; i < ctxt->nbuf; i++) {
        if (ctxt->buf[i].buffer == buf) {
            ctxt->buf[i].active = 0;
            break;
        }
    }
}


/** function to complete an entire call */
void buf_complete_call(buf_context_t *ctxt, int call_id)
{
    int i;
#ifdef DEBUG
    printf("Completing call with call id: %d\n", call_id);
#endif
    for (i = 0; i < ctxt->nbuf; i++) {
        if (ctxt->buf[i].call_id == call_id && ctxt->buf[i].active == 1) {
            /* force completion by calling user function */
            ctxt->fptr(ctxt->buf[i].buffer);
            ctxt->buf[i].active = 0;
        }
    }
}


/** function to return the call_id associated with a particular buffer */
int buf_get_call_id(buf_context_t *ctxt, char *buf)
{
    int i;
    for (i = 0; i < ctxt->nbuf; i++)
        if (ctxt->buf[i].buffer == buf)
            return(ctxt->buf[i].call_id);
    printf("Buf_man error: Cannot find call_id for this buffer\n");
    return -1;
}


/**
 * Function to return an array of buffers associated with a call_id
 * the last two parameters are output
 */
int get_bufs_of_call_id(
        buf_context_t *ctxt, int call_id, int *n_buf, char *bufs[])
{
    int i, count = 0;

    for (i = 0; i < ctxt->nbuf; i++)
        if (ctxt->buf[i].call_id == call_id) {
            bufs[count++] = ctxt->buf[i].buffer;
        }
    *n_buf = count;
    if (*n_buf == 0) {
#ifdef DEBUG 
        printf("Buf_man: No active buffer found for call_id %d\n", call_id);
#endif
        return -1; /* no buffer found */
    }

    return 0; /* success */
}


/** terminates an application context */
void buf_terminate(buf_context_t *ctxt)
{
#ifndef STATBUF
    int i;
    for (i = 0; i < ctxt->nbuf; i++) {
        ctxt->buf[i].buffer -= ctxt->buf[i].align_off;
        free(ctxt->buf[i].buffer);
    }

#endif

    free(ctxt->buf);
    ctxt_count--; /* this context can be reallocated */
}
