#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_SETJMP_H
#   include <setjmp.h>
#endif
#if HAVE_SYS_TYPES_H
#   include <sys/types.h>
#endif
#if HAVE_SYS_SYSCALL_H
#   include <sys/syscall.h>
#endif
#if HAVE_SYS_MMAN_H
#   include <sys/mman.h>
#endif
#if HAVE_SYS_PARAM_H
#   include <sys/param.h>
#endif
#if HAVE_SYS_WAIT_H
#   include <sys/wait.h>
#endif
#if HAVE_SIGNAL_H
#   include <signal.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#if HAVE_STRING_H
#   include <string.h>
#endif
#if HAVE_MATH_H
#   include <math.h>
#endif
#if HAVE_ASSERT_H
#   include <assert.h>
#endif
#if HAVE_UNISTD_H
#   include <unistd.h>
#endif
#if HAVE_FCNTL_H
#   include <fcntl.h>
#endif
#if HAVE_ERRNO_H
#   include <errno.h>
#endif
#if HAVE_DIRENT_H
#   include <dirent.h>
#endif
#if HAVE_STDARG_H
#   include <stdarg.h>
#endif
#define FILE_DS FILE *

FILE_DS armci_storage_fopen(char *filename);
FILE_DS armci_storage_fopen_RONLY(char *filename);
void armci_storage_fclose(FILE_DS filed);

int armci_storage_read_ptr(FILE_DS file_d,void *ptr,size_t size,off_t ofs);
int armci_storage_read_pages(FILE_DS file_d, unsigned long first_page,
                unsigned long *page_arr, unsigned long page_arr_sz,int pagesize,
                off_t ofs);
int armci_storage_write_ptr(FILE_DS file_d,void *ptr,size_t size,off_t ofs);
int armci_storage_write_pages(FILE_DS file_d, unsigned long first_page,
                unsigned long *page_arr, unsigned long page_arr_sz,int pagesize,
                off_t ofs);
