#if HAVE_CONFIG_H
#   include "config.h"
#endif

/*interfaces for saving and restoring GA's */
#include "globalp.h"
#include "message.h"
#include "base.h"
#include "macdecls.h"
#include "ga_ckpt.h"
#include "ga-papi.h"
#include "ga-wapi.h"

#define DEBUG 0

void ga_set_spare_procs(int *spare)
{
    extern int ga_spare_procs;
    ga_spare_procs=*spare;
    armci_set_spare_procs(*spare);
}

int ga_icheckpoint_init(Integer *gas, int num)
{
    int rid,i,hdl;
    armci_ckpt_ds_t ckptds;
    extern ARMCI_Group* ga_get_armci_group_(int);
    /*code needs to be written to check if all the gas have same pgroup*/
    (void)ARMCI_Ckpt_create_ds(&ckptds,2*num);
    for(i=0;i<num*2;i=i+2){
    hdl = gas[i/2]+GA_OFFSET;
       printf("\n%d:i=%d hdl=%d gas=%d %d %d",pnga_nodeid(),i,hdl,gas[i/2],i/2,GA[hdl].p_handle);fflush(stdout);
       ckptds.ptr_arr[i]=&GA[hdl];
       ckptds.sz[i]=sizeof(global_array_t);
       ckptds.saveonce[i]=1;
       ckptds.ptr_arr[i+1]=GA[hdl].ptr[pnga_nodeid()];
       ckptds.sz[i+1]=GA[hdl].size;
    }
    hdl = gas[0]+GA_OFFSET;
    if(GA[hdl].p_handle >=1)
       rid = ARMCI_Ckpt_init(NULL,ga_get_armci_group_(GA[hdl].p_handle),0,0,&ckptds);
    else
       rid = ARMCI_Ckpt_init(NULL,ARMCI_Get_world_group(),0,0,&ckptds);
    for(i=0;i<num;i=i+1){
       int hdl = gas[i]+GA_OFFSET;
       GA[hdl].record_id = rid;
    }
}


/*get the list of changed pages from touched_page_array and rewrite the 
 * changed pages*/
int ga_icheckpoint(Integer *gas, int num)
{
    int i,rc,rid;
    int hdl = gas[0]+GA_OFFSET;
    /*code needs to be written to make sure all gas have same rid*/
    rid = GA[hdl].record_id;
    ARMCI_Ckpt(rid);
}


int ga_irecover(int rid)
{
    int rc;
    /*restore state*/
    /*if longjmp things are hosed */
    if(rid == 0){
       printf("\n%d:in recover\n",pnga_nodeid());
       armci_irecover(rid,1);
    }
    else
      armci_irecover(rid,0);
    /*set the default GA group again which means*/
    set_ga_group_is_for_ft(1);
    /*create the new list just replace the bad guy*/
    return(1);
}


int ga_icheckpoint_finalize(int g_a)
{
    /*free the record id, make it available for next use*/
    /*close the files used for checkpointing*/
}
