#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <stdlib.h>

#include "macommon.h"
#include "global.h"
#include "dra.h"
#include "drap.h"

#define LEN 10

int main(int argc, char **argv)
{
    int from, to, type;
    Integer idata[LEN];
#if 0
    int fd;
#endif
    Integer i, ii, imax, offset, status;
    DoublePrecision ddata[LEN];

    if(argc < 2){
        printf("program prints data from a binary file to screen\n");
        printf("usage: dbg_read.x <filename> <type> <from> <to>\n");
        printf("type: 1 - integer 2 - double \n <from> <to> -range of elements (0, ...)\n");
        return(1);
    }

    type = atoi(argv[2]);
    from = atoi(argv[3]);
    to   = atoi(argv[4]);

    if(from < 0 || to < from) {printf("range error\n"); return 1;}
#if 0
    if(!(fd = dra_el_open(argv[1],DRA_R))){printf("not found\n"); return 1;} 
#else
    /* TODO This must be an old test program using an old API...
     * consider removing this program. */
    return 1;
#endif

    switch (type){

        case 1: 
            for(i=from; i<= to; i+= LEN){
                imax = PARIO_MIN(i+LEN-1,to);
                offset = sizeof(Integer)*i;
#if 0
                status=dra_el_read(idata, sizeof(Integer), imax -i+1, fd, offset);
#else
                status = 1;
#endif
                if(!status)printf("error read failed\n");
                for(ii=0;ii<imax-i+1;ii++) printf("%4ld ",idata[ii]);
                printf("\n");
            }
            break;
        case 2: 
            for(i=from; i<= to; i+= LEN){ 
                imax = PARIO_MIN(i+LEN-1,to);
                offset = sizeof(DoublePrecision)*i;
#if 0
                status=dra_el_read(ddata, sizeof(DoublePrecision), imax -i+1, fd, offset);
#else
                status = 1;
#endif
                if(!status)printf("error read failed\n");
                for(ii=0;ii<imax-i+1;ii++) printf("%f ",ddata[ii]);
                printf("\n");
            }
            break;
        default: printf("type error\n"); return 1;
    }
}
