#ifndef BUFFERS_H_ /* so that this file is not included twice */
#define BUFFERS_H_

/** buffer size --- adjust to be a multiplicity of the
    striping factor in a parallel filesystem */
#define DBL_BUF_SIZE 131072

#define BUF_SIZE     (DBL_BUF_SIZE*sizeof(double))
#define INT_BUF_SIZE (BUF_SIZE/sizeof(int))

/* alignment factor for the internal buffer */
#define ALIGN 16

#define MAXBUF 16 /* max # of buffers that can be used */
#define DEFBUF 4 /* default # of buffers */

/** internal buffer structure */
typedef struct {
    char *buffer;
    int align_off; /**< caching alignment offset */
    int buf_hdl;   /**< buffer handle; index of the buffer in the ctxt */
    int group_id;  /**< identify callback function to use to release buffer */
    int call_id;   /**< id to be used to complete an entire  call */
    int active;    /**< if the buffer active or not */
} _buffer_t;

/** structure to create application context */
typedef struct {
    int ctxt_id;
    _buffer_t *buf;      /**< will be allocated nbuf buffers*/
    int nbuf;
    int size;            /**< in bytes */
    void (*fptr)(char*); /**< array of pointers to functions provided
                              by the application */
    int last_buf;        /* utility variable;
                            contains the last buf assigned in this ctxt */
} buf_context_t;

void buffer_init(buf_context_t *ctxt, int nbuf, int buf_size, void (*fptr)(char*));
char *get_buf(buf_context_t *ctxt, int call_id);
void buf_terminate(buf_context_t *ctxt);
void buf_complete_call(buf_context_t *ctxt, int call_id);
int buf_get_call_id(buf_context_t *ctxt, char *buf);
int get_bufs_of_call_id(buf_context_t *ctxt, int call_id, int *n_buf, char *bufs[]);
void free_buf(buf_context_t *ctxt, char *buf);

#endif /* BUFFERS_H_ */
