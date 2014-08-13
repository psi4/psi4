#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Header: /tmp/hpctools/ga/tcgmsg/ipcv5.0/pfilecopy.c,v 1.5 2004-04-01 02:23:05 manoj Exp $ */

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STRINGS_H
#   include <strings.h>
#endif
#if HAVE_STRING_H
#   include <string.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif

/* extern void free(void *ptr); */

#include "msgtypesc.h"
#include "sndrcv.h"
#include "tcgmsgP.h"

/**
 * Process node0 has a file (assumed unopened) named fname.
 * This file will be copied to all other processes which must
 * simultaneously invoke pfilecopy. Since the processes may be
 * using the same directory one probably ought to make sure
 * that each process uses a different name in the call.
 *
 *    e.g.
 *
 *    on node 0    pfilecopy(99, 0, 'argosin')
 *    on node 1    pfilecopy(99, 0, 'argosin_001')
 *    on node 2    pfilecopy(99, 0, 'argosin_002')
 */
void tcgi_pfilecopy(long *type, long *node0, char *filename)
{
    char *buffer;
    FILE *file;
    long length, nread=32768, len_nread=sizeof(long);
    long typenr = (*type & 32767) | MSGINT;   /* Force user type integer */
    long typebuf =(*type & 32767) | MSGCHR;

    if (!(buffer = malloc((unsigned) nread)))
        Error("pfilecopy: failed to allocate the I/O buffer",nread);

    if (*node0 == NODEID_()) {

        /* I have the original file ... open and check its size */

        if ((file = fopen(filename,"r")) == (FILE *) NULL) {
            (void) fprintf(stderr,"me=%ld, filename = %s.\n",
                           (long)NODEID_(),filename);
            Error("pfilecopy: node0 failed to open original file", *node0);
        }

        /* Quick sanity check on the length */

        (void) fseek(file, 0L, (int) 2);   /* Seek to end of file */
        length = ftell(file);              /* Find the length of file */
        (void) fseek(file, 0L, (int) 0);   /* Seek to beginning of file */
        if ( (length<0) || (length>1e12) )
            Error("pfilecopy: the file length is -ve or very big", length);

        /* Send the file in chunks of nread bytes */

        while (nread) {
            nread = fread(buffer, 1, (int) nread, file);
            BRDCST_(&typenr, (char *) &nread, &len_nread, node0);
            typenr++;
            if (nread) {
                BRDCST_(&typebuf, buffer, &nread, node0);
                typebuf++;
            }
        }
    }
    else {

        /* Open the file for the duplicate */

        if ((file = fopen(filename,"w+")) == (FILE *) NULL) {
            (void) fprintf(stderr,"me=%ld, filename = %s.\n",
                           (long)NODEID_(),filename);
            Error("pfilecopy: failed to open duplicate file", *node0);
        }

        /* Receive data and write to file */

        while (nread) {
            BRDCST_(&typenr, (char *) &nread, &len_nread, node0);
            typenr++;
            if (nread) {
                BRDCST_(&typebuf, buffer, &nread, node0);
                typebuf++;
                if (nread != fwrite(buffer, 1, (int) nread, file))
                    Error("pfilecopy: error data to duplicate file", nread);
            }
        }
    }

    /* Tidy up the stuff we have been using */

    (void) fflush(file);
    (void) fclose(file);
    (void) free(buffer);
}

/** The original C interface to PFCOPY_. */
void PFILECOPY_(long *type, long *node0, char *filename)
{
    tcgi_pfilecopy(type, node0, filename);
}

void PFCOPY_(long *type, long *node0, char *fname, int len)
{
    /* Fortran wrapper around pfilecopy */

    char *filename;

#ifdef DEBUG 
    (void) printf("me=%d, type=%d, node0=%d, fname=%x, fname=%.8s, len=%d\n",
                  NODEID_(), *type, *node0, fname, fname, len);
#endif 

    /* Strip trailing blanks off the file name */

    while ((len > 0) && (fname[len-1] == ' '))
        len--;
    if (len <= 0)
        Error("pfcopy_: file name length is toast", (long) len);

    /* Generate a NULL terminated string */

    filename = malloc( (unsigned) (len+1) );
    if (filename) {
        (void) bcopy(fname, filename, len);
        filename[len] = '\0';
    }
    else
        Error("PFCOPY_: failed to malloc space for filename", (long) len);

    /* Now call the C routine to do the work */

    tcgi_pfilecopy(type, node0, filename);

    (void) free(filename);
}
