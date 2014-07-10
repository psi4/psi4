#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_SYS_TYPES_H
#   include <sys/types.h>
#endif
#if HAVE_SYS_SYSCALL_H
#   include <sys/syscall.h>
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
#if HAVE_STDARG_H
#   include <stdarg.h>
#endif
#include "armci_storage.h"
extern int armci_me;
extern void armci_die(char *msg, int code);
#define DEBUG 0
FILE_DS armci_storage_fopen(char *filename)
{
    FILE_DS file_d;
    remove ( filename );
    file_d = fopen(filename,"w+");
    if(file_d==NULL) armci_die("armci_storage_fopen(): cannot open file",0);
    if(DEBUG); printf("\n%d:filed=%p %s\n",armci_me,file_d,filename);
    return(file_d);
}
FILE_DS armci_storage_fopen_RONLY(char *filename)
{
    FILE_DS file_d;
    file_d = fopen(filename,"r");
    if(file_d==NULL)
       armci_die("armci_storage_fopen_RONLY():cannot open file",0);
    if(DEBUG); printf("\n%d:filed=%p %s\n",armci_me,file_d,filename);
    return(file_d);
}
void armci_storage_fclose(FILE_DS filed)
{
    int rc = fclose(filed);
    if(rc!=0)armci_die("armci_storage_fclose(): cannot close file",0);
}

int armci_storage_read_ptr(FILE_DS file_d,void *ptr,size_t size,off_t ofs)
{
    int rc=0,orc=0,isize=size;

    if(isize<0) armci_die("armci_storage_read_ptr(): Invalid size(<0)", isize);
    
    rc = fseek(file_d,ofs,SEEK_SET);
    if(rc)armci_die("fseek failed in armci_storage_read_ptr",rc);

    while(orc!=isize){
       printf("ptr=%p ofs=%ld size=%ld filed=%p\n", ptr,ofs,size,file_d);
       rc = fread(ptr,1,size,file_d);
       if(!rc) armci_die("armci_storage_read_ptr(): fread failed",0);
       
       orc+=rc;
       if(DEBUG); printf("\n%d:read %d so far of %d\n",armci_me,orc,isize);
       if(orc!=isize){
         ptr+=rc;
         size-=rc;
       }
    }
    rc = fseek(file_d,0,SEEK_SET);
    if(rc)armci_die("fseek failed in armci_storage_read_ptr",rc);
    return 0;
}

int armci_storage_read_pages(FILE_DS file_d, unsigned long first_page,
                unsigned long *page_arr, unsigned long page_arr_sz,int pagesize,
                off_t ofs)
{
    int i,rc=0;
    /*this can be heavily optimized*/ 
    for(i=0;i<page_arr_sz;i++){
       off_t ofs1 = ofs+(page_arr[i]-first_page)*pagesize;
       void *ptr = (void *)(page_arr[i]*pagesize);
       rc=armci_storage_read_ptr(file_d,ptr,pagesize,ofs1);
    }
    return rc;
}

int armci_storage_write_ptr(FILE_DS file_d,void *ptr,size_t size,off_t ofs)
{
    int rc=0,orc=0,isize=size;
    if(DEBUG)printf("in storage_write_ptr %ld is ofs\n",ofs); 
    rc = fseek(file_d,ofs,SEEK_SET);
    if(rc)armci_die("fseek failed in armci_storage_write_ptr",rc);

    while(orc!=isize){
       rc = fwrite(ptr,1,size,file_d);
       orc+=rc;
       if(DEBUG)printf("\n%d:wrote %d so far of %d\n",armci_me,orc,isize);
       if(orc!=isize){
         ptr+=rc;
         size-=rc;
       }
    }
    rc = fseek(file_d,0,SEEK_SET); /* for file sync. data flushed to disk  */
    if(rc)armci_die("fseek failed in armci_storage_write_ptr",rc);
    return 0;
}

int armci_storage_write_pages(FILE_DS file_d, unsigned long first_page,
                unsigned long *page_arr, unsigned long page_arr_sz,int pagesize,
                off_t ofs)
{
    int i,rc=0;
    /*this can be heavily optimized*/ 
    for(i=0;i<page_arr_sz;i++){
       off_t ofs1 = ofs+(page_arr[i]-first_page)*pagesize;
       void *ptr = (void *)(page_arr[i]*pagesize);
       rc = armci_storage_write_ptr(file_d,ptr,pagesize,ofs1);
    }
    return 0;
}
