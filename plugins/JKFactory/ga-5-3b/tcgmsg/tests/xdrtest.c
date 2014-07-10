#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDIO_H
#include <stdio.h>
#endif

#if HAVE_RPC_TYPES_H
#include <rpc/types.h>
#endif
#if HAVE_RPC_XDR_H
#include <rpc/xdr.h>
#endif
#if HAVE_STRING_H
#include <string.h>
#endif

static char *xdrbuf;

static XDR xdrs;

int main(int argc, char **argv)
{
    long data[4];
    u_int len;
    long *temp=data;

    if (argc != 2) {
        return 1;
    }

    xdrbuf = malloc(4096);

    if (strcmp(argv[1], "encode") == 0) {
        xdrmem_create(&xdrs, xdrbuf, 4096, XDR_ENCODE);
        (void) fprintf(stderr," encode xdr_setpos=%d\n",
                       xdr_setpos(&xdrs, (u_int) 0));
        (void) scanf("%ld %ld %ld %ld", data, data+1, data+2, data+3);
        (void) fprintf(stderr,"encode Input longs %ld, %ld, %ld, %ld\n",
                       data[0], data[1], data[2], data[3]);
        len = 4;
        (void) fprintf(stderr,"encode xdr_array=%d\n",
                       xdr_array(&xdrs, (char **) &temp, &len, (u_int) 4096,
                           (u_int) sizeof(long), (xdrproc_t)xdr_long));
        len = 4*4 + 4;
        (void) fprintf(stderr,"encode len=%lu\n", (long unsigned)len);
        (void) fwrite(&len, 4, 1, stdout);
        (void) fwrite(xdrbuf, 1, len, stdout);
        (void) fprintf(stderr,"encode data written\n");
        return 0;
    }
    else {
        xdrmem_create(&xdrs, xdrbuf, 4096, XDR_DECODE);
        (void) fprintf(stderr," decode xdr_setpos=%d\n",
                       xdr_setpos(&xdrs, (u_int) 0));
        (void) fread(&len, 4, 1, stdin);
        (void) fprintf(stderr,"decode len=%lu\n", (long unsigned)len);
        (void) fread(xdrbuf, 1, len, stdin);
        (void) fprintf(stderr,"decode data read\n");
        (void) fprintf(stderr,"decode xdr_array=%d\n",
                       xdr_array(&xdrs, (char **) &temp, &len, (u_int) 4096,
                           (u_int) sizeof(long), (xdrproc_t)xdr_long));
        (void) fprintf(stderr,"decode Input longs %ld, %ld, %ld, %ld\n",
                       data[0], data[1], data[2], data[3]);
        return 0;
    }
}
