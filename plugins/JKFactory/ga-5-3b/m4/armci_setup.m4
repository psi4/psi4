# ARMCI_SETUP
# -----------
# ARMCI is a sensitive piece of code! It mostly depends on the network
# selection, but on occasion the type of system is also important (aka
# TARGET). A lot of this code was ported from the old GNUmakefile, for
# better or worse.
AC_DEFUN([ARMCI_SETUP],
[AC_REQUIRE([GA_ARMCI_NETWORK])
AS_CASE([$ga_armci_network],
[BGML], [
    AC_DEFINE([NO_SHM], [1], [Disables use of shared memory])
    AC_DEFINE([USE_MALLOC], [1], [TODO])
    AC_DEFINE([BGML], [1], [TODO])
    ],
[CRAY_SHMEM], [
    AC_DEFINE([CLUSTER], [1], [TODO])
    AC_DEFINE([CRAY_XT], [1], [TODO])
    AC_DEFINE([CRAY_SHMEM], [1], [TODO])
    ],
[DCMF], [
    AC_DEFINE([NO_SHM], [1], [Disables use of shared memory])
    AC_DEFINE([ARMCIX], [1], [Defined when using DCMF])
    AC_DEFINE([DCMF], [1], [Defined when using DCMF])
    ],
[ELAN3], [
    AC_DEFINE([PTHREADS], [1], [TODO])
    AC_DEFINE([SERVER_THREAD], [1], [TODO])
    AC_DEFINE([_REENTRANT], [1], [TODO])
    AC_DEFINE([ACC_SMP], [1], [TODO])
    AC_DEFINE([QUADRICS], [1], [TODO])
    AS_IF([test x$LIBELAN_NATTACH != x],
        [AC_DEFINE([MULTI_CTX], [1], [TODO])],
        [AC_DEFINE([ALLOC_MUNMAP], [1], [TODO])])
    ],
[ELAN4], [
    AC_DEFINE([PTHREADS], [1], [TODO])
    AC_DEFINE([SERVER_THREAD], [1], [TODO])
    AC_DEFINE([_REENTRANT], [1], [TODO])
    AC_DEFINE([ACC_SMP], [1], [TODO])
    AC_DEFINE([QUADRICS], [1], [TODO])
    AC_DEFINE([DOELAN4], [1], [TODO])
    ],
[GEMINI], [
    AC_DEFINE([LIBONESIDED], [1], [for Gemini])
    AC_DEFINE([CRAY_UGNI], [1], [for Gemini])
    AC_DEFINE([GEMINI], [1], [for Gemini])
    ],
[GM], [
    AC_DEFINE([PTHREADS], [1], [TODO])
    AC_DEFINE([DATA_SERVER], [1], [TODO])
    AC_DEFINE([ALLOW_PIN], [1], [TODO])
    AC_DEFINE([SERVER_THREAD], [1], [TODO])
    AC_DEFINE([_REENTRANT], [1], [TODO])
    AC_DEFINE([GM], [1], [TODO])
    ],
[LAPI], [
    AC_DEFINE([LAPI], [1], [TODO])
    AS_IF([test x$LAPI_RDMA != x],
        [AC_DEFINE([ALLOW_PIN], [1], [TODO])
         AC_DEFINE([LAPI_RDMA], [1], [TODO])])
    ],
[MELLANOX], [
    AC_DEFINE([PTHREADS], [1], [TODO])
    AC_DEFINE([DATA_SERVER], [1], [TODO])
    AC_DEFINE([SERVER_THREAD], [1], [TODO])
    AC_DEFINE([_REENTRANT], [1], [TODO])
    AC_DEFINE([VAPI], [1], [TODO])
    AC_DEFINE([ALLOW_PIN], [1], [TODO])
    AC_DEFINE([MELLANOX], [1], [TODO])
    ],
[MPI_MT], [
    AC_DEFINE([PTHREADS], [1], [TODO])
    AC_DEFINE([DATA_SERVER], [1], [TODO])
    AC_DEFINE([SERVER_THREAD], [1], [TODO])
    AC_DEFINE([MPI_MT], [1], [TODO])
    ],
[MPI_SPAWN], [
    AC_DEFINE([DATA_SERVER], [1], [TODO])
    AC_DEFINE([MPI_SPAWN], [1], [TODO])
    AS_IF([test x$ARMCI_SPAWN_CRAY_XT != x],
        [AC_DEFINE([SPAWN_CRAY_XT], [1], [TODO])])
    ],
[OPENIB], [
    AC_DEFINE([PTHREADS], [1], [TODO])
    AC_DEFINE([DATA_SERVER], [1], [TODO])
    AC_DEFINE([SERVER_THREAD], [1], [TODO])
    AC_DEFINE([_REENTRANT], [1], [TODO])
    AC_DEFINE([VAPI], [1], [TODO])
    AC_DEFINE([ALLOW_PIN], [1], [TODO])
    AC_DEFINE([PEND_BUFS], [1], [TODO])
    AC_DEFINE([OPENIB], [1], [TODO])
    ],
[PORTALS], [
    AC_DEFINE([DATA_SERVER], [1], [TODO])
    AC_DEFINE([_REENTRANT], [1], [TODO])
    AC_DEFINE([SERVER_THREAD], [1], [TODO])
    AC_DEFINE([CRAY_XT], [1], [TODO])
    AC_DEFINE([PORTALS], [1], [TODO])
    ],
[SOCKETS], [
    AC_DEFINE([DATA_SERVER], [1], [TODO])
    AC_DEFINE([SOCKETS], [1], [TODO])
    ],
[VIA], [
    AC_DEFINE([PTHREADS], [1], [TODO])
    AC_DEFINE([DATA_SERVER], [1], [TODO])
    AC_DEFINE([SERVER_THREAD], [1], [TODO])
    AC_DEFINE([_REENTRANT], [1], [TODO])
    AC_DEFINE([VIA], [1], [TODO])
    ]
)
AS_IF([test x$REPORT_SHMMAX != x],
    [AC_DEFINE([REPORT_SHMMAX], [1], [TODO])])
AS_IF([test x$thread_safety = xyes],
    [AC_DEFINE([POSIX_THREADS], [1], [TODO])
     AC_DEFINE([_REENTRANT], [1], [TODO])])
])dnl
