#ifndef SF_H_
#define SF_H_

typedef double SFsize_t;

#include "typesf2c.h"

extern int SF_Create(char* fname, SFsize_t size_hard_limit,
        SFsize_t size_soft_limit, SFsize_t req_size, int *handle);
 
extern int SF_Create_suffix(char* fname, SFsize_t size_hard_limit,
        SFsize_t size_soft_limit, SFsize_t req_size, int *handle, int suffix);

extern int SF_Destroy(int handle);

extern int SF_Rwtor(int handle);

extern int SF_Open(int handle);

extern int SF_Close(int handle);

extern int SF_Fsync(int handle);

extern int SF_Write(int handle, SFsize_t offset, SFsize_t bytes,
        char *buffer, int *req_id);

extern int SF_Read(int handle, SFsize_t offset, SFsize_t bytes,
        char *buffer, int *req_id);

extern int SF_Wait(int *req_id);

extern int SF_Waitall(int *list, int num);

extern void SF_Errmsg(int code, char *msg);

#endif /* SF_H_ */
