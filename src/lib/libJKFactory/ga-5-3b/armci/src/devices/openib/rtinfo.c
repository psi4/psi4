#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STRING_H
#   include <string.h>
#endif

/* determine number of CPUs on the current SMP node- Linux version for now */

int armci_getnumcpus(void)
{
    int numproc = 0;
    FILE* fp;
    char line[80];
    fp = fopen("/proc/cpuinfo","r");

    if(fp == NULL) 
        return -1;
    
    while(!feof(fp)){
        fgets(line,80,fp);
        if(strncmp(line,"processor",9)==0) 
            numproc++;
    }
    fclose(fp);
    return numproc;
}

