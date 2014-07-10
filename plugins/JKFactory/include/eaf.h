#ifndef EAF_H 
#define EAH_H

#if 0
/* This section used by both C and Fortran */
#endif

#define   EAF_RW -1
#define   EAF_W  -2
#define   EAF_R  -3

#ifndef EAF_FORTRAN

/* This section used by only C */

/* This to ensure size_t is defined */
#include <stdio.h>
#include <sys/types.h>

typedef double eaf_off_t;

int  EAF_Aread(int fd, eaf_off_t offset, void *buf, size_t bytes, int *req_id);
int  EAF_Awrite(int fd, eaf_off_t offset, const void *buf, size_t bytes, int *req_id);
int  EAF_Close(int fd);
int  EAF_Delete(const char *fname);
int  EAF_Eof(int code);
void EAF_Errmsg(int code, char *msg);
int  EAF_Length(int fd, eaf_off_t *length);
int  EAF_Length(int fd, eaf_off_t *length);
int  EAF_Open(const char *fname, int type, int *fd);
void EAF_Print_stats(int fd);
int  EAF_Probe(int id, int *status);
int  EAF_Read(int fd, eaf_off_t offset, void *buf, size_t bytes);
int  EAF_Stat(const char *path, int *avail_kb, char *fstype, int fslen);
int  EAF_Truncate(int fd, eaf_off_t length);
int  EAF_Wait(int fd, int id);
int  EAF_Write(int fd, eaf_off_t offset, const void *buf, size_t bytes);

#endif
#endif
