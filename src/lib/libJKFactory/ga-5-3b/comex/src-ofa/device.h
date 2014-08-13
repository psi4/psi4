#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>
#include <string.h>

// Device specific functions
void* COMEXD_register_memory(void *buf, int len);

int COMEXD_deregister_memory(void *buf);

int COMEXD_put_nbi(void *src, void *dst, int bytes, int proc);

int COMEXD_get_nbi(void *src, void *dst, int bytes, int proc);

void COMEXD_network_lock(int proc);

void COMEXD_network_unlock(int proc);

int COMEXD_waitproc(int proc);

int COMEXD_waitall();

int COMEXD_initialize();

int COMEXD_finalize();

