#if HAVE_CONFIG_H
#   include "config.h"
#endif
#define  EXTERN
#include "armcip.h"

int PARMCI_Wait(armci_hdl_t* usr_hdl)
{
    armci_ihdl_t nb_handle = (armci_ihdl_t)usr_hdl;
    int success=0;
    int direct = SAMECLUSNODE(nb_handle->proc);

#ifdef BGML
    assert(nb_handle->cmpl_info);
    BGML_Wait(&(nb_handle->count));
    return(success);
#else

    if(direct) {
        return(success);
    }

    if(nb_handle) {
        if(nb_handle->agg_flag) {
            armci_agg_complete(nb_handle, UNSET);
            return (success);
        }
    }

    if(nb_handle){

#ifdef ARMCI_NB_WAIT

        if(nb_handle->tag==0){
            ARMCI_NB_WAIT(nb_handle->cmpl_info);
            return(success);
        }
#if defined(LAPI) || defined(ALLOW_PIN) || defined(ARMCIX)
        if(nb_handle->tag!=0 && nb_handle->bufid==NB_NONE){
            ARMCI_NB_WAIT(nb_handle->cmpl_info);
            return(success);
        }
#endif

#endif
#ifdef COMPLETE_HANDLE
        COMPLETE_HANDLE(nb_handle->bufid,nb_handle->tag,(&success));
#endif
    }
#endif

    return(success);
}

/** 
 *  * implicit handle 
 *   */
static armci_hdl_t armci_nb_handle[ARMCI_MAX_IMPLICIT];/*implicit non-blocking handle*/
static char hdl_flag[ARMCI_MAX_IMPLICIT];
static int impcount=0;


armci_hdl_t *armci_set_implicit_handle (int op, int proc) 
{
    armci_ihdl_t nbh;
    int i=impcount%ARMCI_MAX_IMPLICIT;
    if(hdl_flag[i]=='1')
        PARMCI_Wait(&armci_nb_handle[i]);

    nbh = (armci_ihdl_t)&armci_nb_handle[i];
#ifdef BGML
    nbh->count=0;
#endif
    nbh->tag   = GET_NEXT_NBTAG();
    nbh->op    = op;
    nbh->proc  = proc;
    nbh->bufid = NB_NONE;
    nbh->agg_flag = 0;
    hdl_flag[i]='1';
    ++impcount;
    return &armci_nb_handle[i];
}

/* wait for all non-blocking operations to finish */
int PARMCI_WaitAll (void) {
#ifdef BGML
    BGML_WaitAll();
#elif ARMCIX
    ARMCIX_WaitAll ();
#else
    int i;
    if(impcount) {
        for(i=0; i<ARMCI_MAX_IMPLICIT; i++) {
            if(hdl_flag[i] == '1') {
                PARMCI_Wait(&armci_nb_handle[i]);
                hdl_flag[i]='0';
            }
        }
    }
    impcount=0;
#endif
    return 0;
}  

/* wait for all non-blocking operations to a particular process to finish */
int PARMCI_WaitProc (int proc) {
#ifdef BGML
    BGML_WaitProc(proc);
#elif ARMCIX
    ARMCIX_WaitProc (proc);
#else
    int i;
    if(impcount) {
        for(i=0; i<ARMCI_MAX_IMPLICIT; i++) {
            if(hdl_flag[i]=='1' &&
                    ((armci_ihdl_t)&armci_nb_handle[i])->proc==proc) {
                PARMCI_Wait(&armci_nb_handle[i]);
                hdl_flag[i]='0';
            }
        }
    }
#endif
    return 0;
}


int PARMCI_Test(armci_hdl_t *usr_hdl)
{
    armci_ihdl_t nb_handle = (armci_ihdl_t)usr_hdl;
    int success=0;
#ifdef BGML
    success=(int)nb_handle->count;
#else
    int direct=SAMECLUSNODE(nb_handle->proc);
    if(direct)return(success);
    if(nb_handle) {
        if(nb_handle->agg_flag) {
            armci_die("test for aggregate handle not yet implemented\n",0);
        }
    }
    if(nb_handle){
#ifdef ARMCI_NB_TEST
        if(nb_handle->tag==0){
            ARMCI_NB_TEST(nb_handle->cmpl_info,&success);
            return(success);
        }
#ifdef LAPI
        if(nb_handle->tag!=0 && nb_handle->bufid==NB_NONE){
            ARMCI_NB_TEST(nb_handle->cmpl_info,&success);
            return(success);
        }
#endif
#endif
#ifdef TEST_HANDLE
        TEST_HANDLE(nb_handle->bufid,nb_handle->tag,(&success));
#endif
    }
#endif
    return(success);
}

