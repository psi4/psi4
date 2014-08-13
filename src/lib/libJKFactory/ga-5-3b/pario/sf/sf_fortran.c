#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STRING_H
#   include <string.h>
#endif

#include "sf.h"
#include "sff2c.h"
#include "farg.h"

#define MAX_NAME 256
static char cname[MAX_NAME+1];

Integer FATR sf_create_(
#ifdef F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
        char *fname,
        SFsize_t *size_hard_limit,
        SFsize_t *size_soft_limit,
        SFsize_t *req_size,
        Integer *handle,
        int len
#else
        char *fname,
        int len,
        SFsize_t *size_hard_limit,
        SFsize_t *size_soft_limit,
        SFsize_t *req_size,
        Integer *handle
#endif
        )
{
    ga_f2cstring(fname, len, cname, MAX_NAME);
    return sfi_create(cname, size_hard_limit, size_soft_limit, req_size, handle);
}


Integer FATR sf_create_suffix_(
#ifdef F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
        char *fname,
        SFsize_t *size_hard_limit,
        SFsize_t *size_soft_limit,
        SFsize_t *req_size,
        Integer *handle,
        Integer *suffix,
        int len
#else
        char *fname,
        int len,
        SFsize_t *size_hard_limit,
        SFsize_t *size_soft_limit,
        SFsize_t *req_size,
        Integer *handle,
        Integer *suffix
#endif
        )
{
    ga_f2cstring(fname, len, cname, MAX_NAME);
    return sfi_create_suffix(cname, size_hard_limit, size_soft_limit,
            req_size, handle, suffix);
}

/*****************************************************************************/

void FATR sf_errmsg_(Integer *code, char *msg, int msglen)
{
    char buf[80];

    sfi_errmsg((int) *code, buf);

    (void) ga_c2fstring(buf, msg, msglen);
}
