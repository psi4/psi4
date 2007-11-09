/*! \file 
    \ingroup (IPV1)
    \brief Enter brief description of file here 
*/
#ifndef _psi_src_lib_libipv1_tmpl_h_
#define _psi_src_lib_libipv1_tmpl_h_

#ifdef __cplusplus
extern "C" {
#endif

#ifdef DEC
#define NO_CONST
#endif

#ifdef NCUBE_V2
#define NO_CONST
#endif

#ifdef NO_CONST
#define CONST
#else
#define CONST const
#endif

#define VOID void
#define VOID_PTR void *

#ifdef PROTO_LIBRARY

VOID_PTR malloc(size_t);
VOID_PTR realloc(VOID_PTR,size_t);
VOID free(VOID_PTR);

FILE *fopen(CONST char *, CONST char *);
int   fclose(FILE *);
int fprintf(FILE *, CONST char *, ... );
int printf(CONST char *, ... );
int sprintf(char *, CONST char *, ... );

VOID perror(CONST char *);

VOID exit(int);

int strcmp(CONST char *, CONST char *);
int strlen(CONST char *);
char *strcpy(char *, CONST char *);
char *strcat(char *, CONST char *);

/* These are needed for getc and putc (to avoid sun3 gcc warnings). */
int _filbuf(void), _flsbuf(void);

int yylex(void);

int yyparse(void);

/* This is for gcc. */
void abort( void );

#endif /* PROTO_LIBRARY */

#ifdef __cplusplus
}
#endif

#endif /* header guard */

