#ifndef _psi_src_lib_libciomr_types_h_
#define _psi_src_lib_libciomr_types_h_

#include "iomrparam.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef
struct {
  int junk;
  } r_async_t;

typedef
struct {
  int junk;
  } s_async_t;

typedef
struct {
  char *path;
#if BUFF
  FILE *stream;
#else
  int stream;
#endif
  } sequential_volume_t;

typedef
struct {
  int n;
  int blocksize;
  int last_ioop;
  int verbose;
  PSI_FPTR next;
  PSI_FPTR previous_size; /* added for error checking with unsigned's */ 
  int unit; /* The unit number. */
  PSI_FPTR incount;  /* The number of bytes read. */
  PSI_FPTR outcount;  /* The number of bytes written. */
  sequential_volume_t v[MAX_VOLUME];
  } sequential_t;

typedef
struct {
  int junk;
  } ram_t;
  

typedef
struct {
  int method;
  union {
    r_async_t *r_async;
    s_async_t *s_async;
    sequential_t *sequential;
    ram_t *ram;
    } ptr;
  } ioFILE_t;

  sequential_t *
  sequential_ioopen(char* baseparam,int unit);
  void
  sequential_ioclos(sequential_t* ud, int status);
  void
  sequential_iordr(sequential_t* ud, char* buffer,PSI_FPTR first,int length);
  void
  sequential_iowrr(sequential_t *ud,char* buffer,PSI_FPTR first,int length);
  void
  sequential_iordwrr(char* caller,int ioop,sequential_t* ud,char* buffer,PSI_FPTR first,int length);
  PSI_FPTR sequential_iosize(sequential_t* ud);
  r_async_t * r_async_ioopen(char *param, int unit);
  void r_async_ioclos(r_async_t *ud, int status);
  void r_async_iordr(r_async_t *ud, char *buffer, PSI_FPTR first, int length);
  void r_async_iowrr(r_async_t *ud, char *buffer, PSI_FPTR first, int length);
  PSI_FPTR r_async_iosize(r_async_t *ud);
  s_async_t *s_async_ioopen(char *param, int unit);
  void s_async_ioclos(s_async_t *ud, int status);
  void s_async_iordr(s_async_t *ud, char *buffer, PSI_FPTR first, int length);
  void s_async_iowrr(s_async_t *ud, char *buffer, PSI_FPTR first, int length);
  PSI_FPTR s_async_iosize(s_async_t *ud);
  ram_t * ram_ioopen(char *param, int unit);
  void ram_ioclos(ram_t *ud, int status);
  void ram_iordr(ram_t *ud, char *buffer, PSI_FPTR first, int length);
  void ram_iowrr(ram_t *ud, char *buffer, PSI_FPTR first, int length);
  PSI_FPTR ram_iosize(ram_t *ud);



#ifdef __cplusplus
}
#endif

#endif /* header guard */
