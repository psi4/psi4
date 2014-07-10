#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>
#include <string.h>

#include "comex.h"
#include "comex_impl.h"
#include "openib.h"

// Device specific implementation
//
void* COMEXD_register_memory(void *buf, int len)
{
    return openib_register_memory(buf, len);
}

int COMEXD_deregister_memory(void *buf)
{
    return openib_deregister_memory(buf);
}

int COMEXD_put_nbi(void *src, void *dst, int bytes, int proc)
{
    if (proc == l_state.rank) {
        (void) memcpy(dst, src, bytes);
        return 0;
    }

    return openib_put_nbi(src, dst, bytes, proc);
}

int COMEXD_get_nbi(void *src, void *dst, int bytes, int proc)
{
    if (proc == l_state.rank) {
        (void) memcpy(dst, src, bytes);
        return 0;
    }

    return openib_get_nbi(src, dst, bytes, proc);
}

void COMEXD_network_lock(int proc)
{
    openib_network_lock(proc);
}

void COMEXD_network_unlock(int proc)
{
    openib_network_unlock(proc);
}

int COMEXD_waitproc(int proc)
{
    return openib_waitproc(proc);    
}

int COMEXD_waitall()
{
    openib_waitall();    
    return 0;
}

int COMEXD_initialize()
{
    return openib_initialize();
}

int COMEXD_finalize()
{
    return openib_finalize();
}
