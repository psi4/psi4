#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include "globalp.h"
#include "base.h"
#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#define DEBUG 0

/*The structure of gai_nbhdl_t is (this is our internal handle)*/
#define NUM_HDLS 20
typedef struct {
    unsigned int ihdl_index:8;
    unsigned int ga_nbtag:24;
}gai_nbhdl_t;


/*Each element in the armci handle linked list is of type ga_armcihdl_t*/
typedef struct struct_armcihdl_t{
    armci_hdl_t* handle;
    struct struct_armcihdl_t *next;
    struct struct_armcihdl_t *previous;
    int index;
    int ga_hdlarr_index;
}ga_armcihdl_t;


/*We create an array of type ga_nbhdl_array_t. Each of the elements in this
  array is the head of the armcihandle linked list that is associated with
  each GA call.
*/
typedef struct{
    ga_armcihdl_t *ahandle;
    int count;
    int ga_nbtag;
} ga_nbhdl_array_t;


/*fills up the armci_hdl_t entries in ga_armcihdl_t */
static armci_hdl_t hdl_array[NUM_HDLS];


/*index of the following array goes into ihdl_index. while waiting for a
 *non-bloking ga call, we first check if
 *(list_element_array[inbhandle->ihdl_index].ga_nbtag == inbhandle->ga_nbtag)
 *if it is, then we complete all the armci handles in the linked list this
 *points to.
*/
static ga_nbhdl_array_t ga_ihdl_array[NUM_HDLS];


/*this is the array of linked list elements. */
static ga_armcihdl_t list_element_array[NUM_HDLS] = {
{&(hdl_array[0]), NULL,NULL,0, -1 },{&(hdl_array[1]), NULL,NULL, 1,-1 }, 
{&(hdl_array[2]), NULL,NULL,2, -1 },{&(hdl_array[3]), NULL,NULL, 3,-1 },
{&(hdl_array[4]), NULL,NULL,4, -1 },{&(hdl_array[5]), NULL,NULL, 5,-1 }, 
{&(hdl_array[6]), NULL,NULL,6, -1 },{&(hdl_array[7]), NULL,NULL, 7,-1 }, 
{&(hdl_array[8]), NULL,NULL,8, -1 },{&(hdl_array[9]), NULL,NULL, 9,-1 },
{&(hdl_array[10]),NULL,NULL,10,-1 },{&(hdl_array[11]),NULL,NULL,11,-1 },
{&(hdl_array[12]),NULL,NULL,12,-1 },{&(hdl_array[13]),NULL,NULL,13,-1 },
{&(hdl_array[14]),NULL,NULL,14,-1 },{&(hdl_array[15]),NULL,NULL,15,-1 },
{&(hdl_array[16]),NULL,NULL,16,-1 },{&(hdl_array[17]),NULL,NULL,17,-1 },
{&(hdl_array[18]),NULL,NULL,18,-1 },{&(hdl_array[19]),NULL,NULL,19,-1 }};






static int nextIHAelement=-1; /*oldest ga_ihdl_array element*/
static int nextLEAelement=-1; /*oldest list_element_array element*/
static int ihdl_array_avail[NUM_HDLS]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
static int list_ele_avail[NUM_HDLS]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

/*\ a unique tag everytime
\*/
static unsigned int ga_nb_tag;
unsigned int get_next_tag(){
    return((++ga_nb_tag));
}

/*\ the only way to complete a list element! 
 *  does the basic list operation: remove element, update previous and next
 *  links of the previous and next elements in the linked list
 *  prev==null => this was the element pointed by the head(ie, first element).
\*/
static void clear_list_element(int index){
ga_armcihdl_t *listele,*prev,*next;
    if(DEBUG){
       printf("\n%ld:clearing handle %d\n",(long)GAme,index);fflush(stdout);
    }
    listele = &(list_element_array[index]);

    /*first wait for the armci handle */
    ARMCI_Wait(listele->handle);

    /*set prev and next links of my prev element and my next element*/
    prev=listele->previous;
    next = listele->next;
    if(prev)
       prev->next = next;
    else
       ga_ihdl_array[listele->ga_hdlarr_index].ahandle=next;
    if(next)
       next->previous = prev; 

    /*since one element from the linked list is completed, update the count*/
    ga_ihdl_array[listele->ga_hdlarr_index].count--;

    /*reset the prev and next pointers and initialize the handle*/
    listele->next=NULL;
    listele->previous=NULL;
    ARMCI_INIT_HANDLE(listele->handle);
    list_ele_avail[index]=1;
}


/*\ Get the next available list element from the list element array, if 
 *  nothing is available, free element with index nextLEAelement
\*/
ga_armcihdl_t* get_armcihdl(){
int i;
ga_armcihdl_t *ret_handle;

    /*first see if an element from the list_ele_arr is already available */
    for(i=0;i<NUM_HDLS;i++)
       if(list_ele_avail[i]){
         list_ele_avail[i]=0;
         ARMCI_INIT_HANDLE(list_element_array[i].handle);
         if(DEBUG){
           printf("\n%ld:found a free handle %d\n",(long)GAme,i);fflush(stdout);
         }
         return(&(list_element_array[i]));
       }

    /*nothing is available so best element to clear is nextLEAelement(LRU)*/
    if(nextLEAelement==-1)
       nextLEAelement=0;
    if(DEBUG){
       printf("\n%ld:have to clear handle %d\n",(long)GAme,nextLEAelement);
       fflush(stdout);
    }
    clear_list_element(nextLEAelement);
    list_ele_avail[nextLEAelement]=0;

    ret_handle=&(list_element_array[nextLEAelement]);

    /*update the LRU element index */
    nextLEAelement = (nextLEAelement+1)%NUM_HDLS;
    return(ret_handle);
}

/*\ Input is the index to the ga_ihdl_array that has the head of the list.
 *  This function waits for all the elements in the list.
\*/
static void free_armci_handle_list(int elementtofree){
ga_armcihdl_t *first = ga_ihdl_array[elementtofree].ahandle,*next;
    /*call clear_list_element for every element in the list*/
    while(first!=NULL){
       next=first->next;
       clear_list_element(first->index);
       first=next;
    }

    /*reset the head of the list for reuse*/
    ga_ihdl_array[elementtofree].count=0;
    ga_ihdl_array[elementtofree].ga_nbtag=0;
    ga_ihdl_array[elementtofree].ahandle=NULL;
    ihdl_array_avail[elementtofree]=1;
}
      

/*\ Add the armci handle list element to the end of the list.
\*/
static void add_armcihdl_to_list(ga_armcihdl_t *listelement, int headindex){
ga_armcihdl_t *first=ga_ihdl_array[headindex].ahandle;

    ga_ihdl_array[headindex].count++;
    listelement->ga_hdlarr_index = headindex;
    if(ga_ihdl_array[headindex].ahandle==NULL){
       ga_ihdl_array[headindex].ahandle=listelement;
       listelement->previous= NULL;
       return;
    }
    while(first->next!=NULL){
       first=first->next;
    }
    first->next=listelement;
    listelement->previous=first;
}


/*\ Complete the list of armci handles associated with a particular GA request.
 *  specific=-1 means free the next available one. other values complete the
 *  armci handle list pointed to by head at that "specific" element.
\*/
static int get_GAnbhdl_element(int specific){
int elementtofree,i;
    if(specific!=-1)elementtofree=specific;
    else {
       for(i=0;i<NUM_HDLS;i++)
         if(ihdl_array_avail[i]){
           ihdl_array_avail[i]=0;
           return(i);
         }
       if(nextIHAelement==-1)  
         nextIHAelement=0;       
       elementtofree=nextIHAelement;
       nextIHAelement = (elementtofree+1)%NUM_HDLS;
    }
    free_armci_handle_list(elementtofree);
    return(elementtofree);
}


/*\ called from ga_put/get before a call to every non-blocking armci request. 
\*/
armci_hdl_t* get_armci_nbhandle(Integer *nbhandle){
gai_nbhdl_t *inbhandle = (gai_nbhdl_t *)nbhandle;
ga_armcihdl_t *ret_handle;
    if(inbhandle->ihdl_index == (NUM_HDLS+1)){
       inbhandle->ihdl_index = get_GAnbhdl_element(-1);
       inbhandle->ga_nbtag = get_next_tag();
       ga_ihdl_array[(inbhandle->ihdl_index)].ga_nbtag=inbhandle->ga_nbtag; 
    }
    ret_handle = get_armcihdl(); 
    add_armcihdl_to_list(ret_handle,inbhandle->ihdl_index);
    return(ret_handle->handle);
}

/*\ the wait routine which is called inside nga_nbwait and ga_nbwait
\*/ 
int nga_wait_internal(Integer *nbhandle){
gai_nbhdl_t *inbhandle = (gai_nbhdl_t *)nbhandle;
int retval = 0;
    if(inbhandle->ihdl_index==(NUM_HDLS+1))retval=0;
    else if(inbhandle->ga_nbtag !=ga_ihdl_array[inbhandle->ihdl_index].ga_nbtag)
       retval=0;
    else
       free_armci_handle_list(inbhandle->ihdl_index);
    
    return(retval);
}


static int test_list_element(int index){
ga_armcihdl_t *listele;
    if(DEBUG){
       printf("\n%ld:clearing handle %d\n",(long)GAme,index);fflush(stdout);
    }
    listele = &(list_element_array[index]);

    return (ARMCI_Test(listele->handle));
}

static int test_armci_handle_list(int elementtofree){
ga_armcihdl_t *first = ga_ihdl_array[elementtofree].ahandle,*next;
 int done = 1; 
    /*call clear_list_element for every element in the list*/
    while(first!=NULL){
       next=first->next;
       if (test_list_element(first->index) == 0) {
	 done = 0;
	 break;
       }
       first=next;
    }
    return (done);
}

/*\ the test routine which is called inside nga_nbtest
\*/ 
int nga_test_internal(Integer *nbhandle){
gai_nbhdl_t *inbhandle = (gai_nbhdl_t *)nbhandle;
int retval = 0;
    if(inbhandle->ihdl_index==(NUM_HDLS+1))retval=0;
    else if(inbhandle->ga_nbtag !=ga_ihdl_array[inbhandle->ihdl_index].ga_nbtag)
       retval=0;
    else
       return (test_armci_handle_list(inbhandle->ihdl_index));
    
    return(retval);
}

/*\ unlike in ARMCI, user doesnt have to initialize handle in GA.
 *  it is done by the get/put call instead.
\*/
void ga_init_nbhandle(Integer *nbhandle)
{
gai_nbhdl_t *inbhandle = (gai_nbhdl_t *)nbhandle;
    inbhandle->ihdl_index=(NUM_HDLS+1);
}
