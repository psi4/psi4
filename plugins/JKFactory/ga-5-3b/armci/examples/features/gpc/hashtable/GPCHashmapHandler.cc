/* $Id:  */
#include <stdio.h>

#define DEBUG 0

#ifdef WIN32
#  include <windows.h>
#  define sleep(x) Sleep(1000*(x))
#else
#  include <unistd.h>
#endif

#include "armci.h"
#define ARMCI_ENABLE_GPC_CALLS
#include "gpc.h"

#include "Hash_common.h"
#include "GPCHashmap.h"

static GPCHashmap *gGPCHashmap=NULL;

#if DEBUG
char* g_GPCops[6] = {"HASHMAP_CREATE", "HASHMAP_DESTROY", "HASHMAP_INSERT",
                      "HASHMAP_PRINT", "HASHMAP_GET", "HASHMAP_REHASH"};
#endif

int gpc_disthashmap_handler(int to, int from, void *hdr,   int hlen,
                            void *data,  int dlen,
                            void *rhdr,  int rhlen, int *rhsize,
                            void *rdata, int rdlen, int *rdsize,
                            int rtype)
{
hash_hdr_t *lhdr;
int rank; MP_MYID(&rank);
 
     lhdr = (hash_hdr_t*)ARMCI_Gpc_translate(hdr,to,from);

     const char *buf = (const char*)data;
     size_t bufsize  = dlen;
     int hash_op     = lhdr->hash_op;

     // for correctness ?
     *rhsize = sizeof(int);
     *rdsize = sizeof(int);
     
#if DEBUG
     int me; MP_MYID(&me);
     printf("%d: In GPC handler: %s from=%d to=%d\n",
            me, g_GPCops[hash_op-HASHMAP_CREATE], from, to); fflush(stdout);
#endif
     
     switch(hash_op) {
       case HASHMAP_CREATE:
          gGPCHashmap = new GPCHashmap();
          gGPCHashmap->create();
         break;

       case HASHMAP_DESTROY:
         gGPCHashmap->destroy();
         if(gGPCHashmap != NULL) delete gGPCHashmap;
         break;

       case HASHMAP_INSERT:
         gGPCHashmap->insert(buf, bufsize);
         // piggy-back the globalIds
         gGPCHashmap->getGlobalIds(buf, bufsize, (int*)rdata);
         break;

       case HASHMAP_PRINT:
         gGPCHashmap->print();
         break;

       case HASHMAP_REHASH:
         *rdsize = sizeof(int);
         gGPCHashmap->rehash((int*)rdata);
         break;

       default:
         ARMCI_Error("gpc_disthashmap_handler(): Invalid hashmap operation",0);
     }

     return GPC_DONE;
}

void gpc_disthashmap_exec_nb(int hash_op, char *buf, size_t bufsize,
                             int proc, int gpc_handle, gpc_hdl_t *nbh)
{
int hlen;
hash_hdr_t header;
int rheader;

#if DEBUG
     int me; MP_MYID(&me);
     printf("%d: Executing GPC: (%s) buf=%p bufsize=%ld proc=%d\n",
            me, g_GPCops[hash_op-HASHMAP_CREATE], buf, bufsize, proc);
     fflush(stdout);
#endif 
 
     header.hash_op = hash_op;
     header.buf = buf;
     header.bufsize = bufsize;
     hlen=sizeof(header);

     if(nbh!= NULL) ARMCI_Gpc_init_handle(nbh);

#if 0
     if(ARMCI_Gpc_exec(gpc_handle, proc, &header, hlen, buf, bufsize,
                       &rheader, sizeof(int), &rdata, sizeof(int), nbh))
        ARMCI_Error("gpc_disthashmap_exec_nb(): ARMCI_Gpc_exec failed", 0);
#endif

     if(ARMCI_Gpc_exec(gpc_handle, proc, &header, hlen, buf, bufsize,
                       &rheader, sizeof(int), buf, bufsize, nbh))
       ARMCI_Error("gpc_disthashmap_exec_nb(): ARMCI_Gpc_exec failed", 0);
}

void gpc_disthashmap_exec_wait(gpc_hdl_t *nbh)
{
    if(nbh!= NULL) ARMCI_Gpc_wait(nbh);
}

void gpc_disthashmap_exec(int hash_op, char *buf, size_t bufsize,
                          int proc, int gpc_handle)
{

#if NON_BLOCKING // non blocking GPC is flaky
gpc_hdl_t nbh;
     gpc_disthashmap_exec_nb(hash_op, buf, bufsize, proc, gpc_handle, &nbh);
     gpc_disthashmap_exec_wait(&nbh);
#else
     gpc_disthashmap_exec_nb(hash_op, buf, bufsize, proc, gpc_handle, NULL);
     gpc_disthashmap_exec_wait(NULL);
#endif
     /*
int hlen;
hash_hdr_t header;
int rheader, rdata;
     
     header.hash_op = hash_op;
     header.buf = buf;
     header.bufsize = bufsize;
     hlen=sizeof(header);

     ARMCI_Gpc_init_handle(&nbh);

     if(ARMCI_Gpc_exec(gpc_handle, proc, &header, hlen, buf, bufsize,
                       &rheader, sizeof(int), &rdata, sizeof(int), &nbh))
       ARMCI_Error("gpc_disthashmap_exec(): ARMCI_Gpc_exec failed", 0);
     
     ARMCI_Gpc_wait(&nbh);
     */     
}

