#ifndef SFF2C_H_
#define SFF2C_H_

#define sf_close_          F77_FUNC_(sf_close,SF_CLOSE)
#define sf_create_         F77_FUNC_(sf_create,SF_CREATE)
#define sf_create_suffix_  F77_FUNC_(sf_create_suffix,SF_CREATE_SUFFIX)
#define sf_destroy_        F77_FUNC_(sf_destroy,SF_DESTROY)
#define sf_errmsg_         F77_FUNC_(sf_errmsg,SF_ERRMSG)
#define sf_open_           F77_FUNC_(sf_open,SF_OPEN)
#define sf_read_           F77_FUNC_(sf_read,SF_READ)
#define sf_rwtor_          F77_FUNC_(sf_rwtor,SF_RWTOR)
#define sf_waitall_        F77_FUNC_(sf_waitall,SF_WAITALL)
#define sf_wait_           F77_FUNC_(sf_wait,SF_WAIT)
#define sf_write_          F77_FUNC_(sf_write,SF_WRITE)
#define sf_fsync_          F77_FUNC_(sf_fsync,SF_FSYNC)

extern Integer sfi_create(char *fname, SFsize_t *size_hard_limit, SFsize_t *size_soft_limit, SFsize_t *req_size, Integer *handle);
extern Integer sfi_create_suffix(char *fname, SFsize_t *size_hard_limit, SFsize_t *size_soft_limit, SFsize_t *req_size, Integer *handle, Integer *suffix);
extern void    sfi_errmsg(int code, char *msg);

extern Integer FATR sf_close_(Integer *s_a);
extern Integer FATR sf_destroy_(Integer *s_a);
extern Integer FATR sf_open_(Integer *s_a);
extern Integer FATR sf_read_(Integer *s_a, SFsize_t *offset, SFsize_t *bytes, char *buffer, Integer *req_id);
extern Integer FATR sf_rwtor_(Integer *s_a);
extern Integer FATR sf_waitall_(Integer *list, Integer *num);
extern Integer FATR sf_wait_(Integer *req_id);
extern Integer FATR sf_write_(Integer *s_a, SFsize_t *offset, SFsize_t *bytes, char *buffer, Integer *req_id);
extern Integer FATR sf_fsync_(Integer *s_a);

#ifdef F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
extern Integer sf_create(char *fname, SFsize_t *size_hard_limit, SFsize_t *size_soft_limit, SFsize_t *req_size, Integer *handle, int len);
extern Integer sf_create_suffix(char *fname, SFsize_t *size_hard_limit, SFsize_t *size_soft_limit, SFsize_t *req_size, Integer *handle, Integer *suffix, int len);
extern void    sf_errmsg(Integer *code, char *msg, int len);
#else /* F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS */
extern Integer sf_create(char *fname, int len, SFsize_t *size_hard_limit, SFsize_t *size_soft_limit, SFsize_t *req_size, Integer *handle);
extern Integer sf_create_suffix(char *fname, int len, SFsize_t *size_hard_limit, SFsize_t *size_soft_limit, SFsize_t *req_size, Integer *handle, Integer *suffix);
extern void    sf_errmsg(Integer *code, char *msg, int len);
#endif /* F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS */

#endif /* SFF2C_H_ */
