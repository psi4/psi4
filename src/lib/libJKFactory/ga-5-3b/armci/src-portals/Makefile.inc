libarmci_la_SOURCES += src-portals/acc.h
libarmci_la_SOURCES += src-portals/aggregate.c
libarmci_la_SOURCES += src-portals/armci.c
libarmci_la_SOURCES += src-portals/armcip.h
libarmci_la_SOURCES += src-portals/armci_portals.c
libarmci_la_SOURCES += src-portals/armci_portals.h
libarmci_la_SOURCES += src-portals/atomics-i386.h
libarmci_la_SOURCES += src-portals/buffers.c
libarmci_la_SOURCES += src-portals/caccumulate.c
libarmci_la_SOURCES += src-portals/ccopy.c
libarmci_la_SOURCES += src-portals/clusterinfo.c
libarmci_la_SOURCES += src-portals/copy.h
libarmci_la_SOURCES += src-portals/ds-shared.c
libarmci_la_SOURCES += src-portals/fence.c
libarmci_la_SOURCES += src-portals/kr_malloc.c
libarmci_la_SOURCES += src-portals/kr_malloc.h
libarmci_la_SOURCES += src-portals/locks.c
libarmci_la_SOURCES += src-portals/locks.h
libarmci_la_SOURCES += src-portals/memlock.c
libarmci_la_SOURCES += src-portals/memlock.h
libarmci_la_SOURCES += src-portals/memory.c
libarmci_la_SOURCES += src-portals/message.c
libarmci_la_SOURCES += src-portals/mutex.c
libarmci_la_SOURCES += src-portals/pack.c
libarmci_la_SOURCES += src-portals/pendbufs.h
libarmci_la_SOURCES += src-portals/portals_cp.c
libarmci_la_SOURCES += src-portals/portals_ds.c
libarmci_la_SOURCES += src-portals/portals.c
libarmci_la_SOURCES += src-portals/request.c
libarmci_la_SOURCES += src-portals/request.h
libarmci_la_SOURCES += src-portals/rmw.c
libarmci_la_SOURCES += src-portals/rtinfo.c
libarmci_la_SOURCES += src-portals/semaphores.c
libarmci_la_SOURCES += src-portals/semaphores.h
libarmci_la_SOURCES += src-portals/shmalloc.h
libarmci_la_SOURCES += src-portals/shmem.c
libarmci_la_SOURCES += src-portals/shmem.h
libarmci_la_SOURCES += src-portals/shmlimit.c
libarmci_la_SOURCES += src-portals/shmlimit.h
libarmci_la_SOURCES += src-portals/signaltrap.c
libarmci_la_SOURCES += src-portals/signaltrap.h
libarmci_la_SOURCES += src-portals/sockets.h
libarmci_la_SOURCES += src-portals/spawn.c
libarmci_la_SOURCES += src-portals/spinlock.h
libarmci_la_SOURCES += src-portals/strided.c
libarmci_la_SOURCES += src-portals/utils.h
libarmci_la_SOURCES += src-portals/vector.c
if PORTALS_ENABLE_NEW_MALLOC
libarmci_la_SOURCES += src-portals/new_memory.c
endif
if MSG_COMMS_MPI
libarmci_la_SOURCES += src-portals/groups.c
endif
if ARMCI_ENABLE_GPC_CALLS
libarmci_la_SOURCES += src-portals/gpc.c
endif
if THREAD_SAFE
libarmci_la_SOURCES += src-portals/threads.c
libarmci_la_SOURCES += src-portals/utils.c
endif

include_HEADERS += src-portals/armci.h
include_HEADERS += src-portals/gpc.h
include_HEADERS += src-portals/message.h

AM_CPPFLAGS += -I$(top_srcdir)/src-portals
AM_CPPFLAGS += -I$(top_srcdir)/src/include
