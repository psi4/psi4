/*! \file
    \ingroup R12
    \brief Enter brief description of file here 
*/
/*------------------------------------------------------------------------------------------------------
  Builds a library of functions-calls for applying HRR and VRR
 ------------------------------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libint/libint.h>
#include "mem_man.h"
#include "build_libr12.h"
#include <libint/constants.h>
#define MAXNODE 3000
#define NONODE -1000000
#define SSR12SSNODE -1111 /* This is a special node - (ss||ss) */
#define NUMPARENTS 21
#define NUMCHILDREN 6
#define NUMGRTTYPES 4     /* There are 4 types classes to be evaluated - ERIs (g), r12 (r), [r12,t1] (t1), [r12,t2] (t2) */
#define MOFFSET     11    /* This is to take care of a paossibility that class 0 will be someone's child. Otherwise
			     it would be confused with (ss|ss)^0 */

static int last_hrr_node = 0;      /* Global pointer to the last node on the HRR stack */
static int last_vrr_node = 0;      /* Global pointer to the last node on the VRR stack */

extern FILE *outfile, *hrr_header, *init_code;
extern int libr12_stack_size[MAX_AM/2+1];
extern Libr12Params_t Params;

typedef struct node{
  int A, B, C, D;         /* Angular momenta on centers A and C */
  int m;
  int grt_type;           /* Type of the class - 0 = ERI (g), 1 = r12 (r), 2 = [r12,t1] (t1), 3 = [r12,t2] (t2) */
  int size;               /* Class size in double words */
  int pointer;
  int children[NUMCHILDREN];        /* Up to 8 children of the class */
  int parents_counter;
  int num_parents;        /* Number of parents */
  int parents[NUMPARENTS];         /* Pointers to parents */
  int llink;              /* Pointer to a class computed right before computing this one */
  int rlink;              /* Pointer to a class to be computed after this one is */
  int marked;             /* Flag indicating that this node has been computed */
  int target;             /* Flag indicating that this node is among targets */
  } class;

static int first_hrr_to_compute = 0; /* Number of the first class to be computed
				    (pointer to the beginning of the linked list) */
static int first_vrr_to_compute = 0; /* Number of the first class to be computed
				    (pointer to the beginning of the linked list) */

/*--- This is the maximum ang. momentum allowed for any (intermediate) classes ---*/
#define LMAX_AM LIBINT_MAX_AM
static int hrr_hash_table[NUMGRTTYPES][2*LMAX_AM][2*LMAX_AM][2*LMAX_AM][2*LMAX_AM];
static int vrr_hash_table[NUMGRTTYPES][2*LMAX_AM][2*LMAX_AM][4*LMAX_AM];

void mark_hrr_parents(int n, class *allnodes, int rent);
void mark_vrr_parents(int n, class *allnodes, int rent);
int mk_hrr_node(class node, class *allnodes, int new);
int mk_vrr_node(class node, class *allnodes, int new);
int alloc_mem_vrr(class *nodes);
int alloc_mem_hrr(class *nodes);

void emit_grt_order()
{
  int old_am = Params.old_am;
  int new_am = Params.new_am;
  int opt_am = Params.opt_am;

  int i, j, k, l;
  int la, lc, lc_min, ld, ld_max, ld_min;
  int lb, lb_min, lb_max;
  int current_highest_am;
  int base_mem, hrr_mem, vrr_mem;
  int child0, child1, child;
  int num_children;
  int offset;
  class nodes[MAXNODE];    /* Stack of nodes */
  class *hrr_nodes = &(nodes[0]);
  class *vrr_nodes;
  int target_data;
  int done;
  int max_stack_size = 0;
  int target_hrr_nodes[NUMGRTTYPES];         /* Array of unique targets on the HRR graph */ 
  int num_hrr_targets;
  int target_vrr_nodes[1000];
  int num_vrr_targets;
  char hrr_code_name[] = "hrr_grt_order_0000.cc";
  char hrr_function_name[] = "hrr_grt_order_0000";
  char vrr_code_name[] = "vrr_grt_order_0000.cc";
  char vrr_function_name[] = "vrr_grt_order_0000";
  FILE *hrr_code, *vrr_code;


  for(la=0;la<=new_am;la++) {
    lb_max = la/2;
    lb_min = (la <= new_am/2) ? 0 : la - new_am/2;
    lc_min = la;
    for(lb=lb_min;lb<=lb_max;lb++) {
    for(lc=lc_min;lc<=new_am;lc++) {
      ld_max = lc/2;
      ld_min = (lc <= new_am/2) ? 0 : lc - new_am/2;
      for(ld=ld_min;ld<=ld_max;ld++) {

      current_highest_am = (la-lb > lb) ? la-lb : lb;
      current_highest_am = (current_highest_am > lc-ld) ? current_highest_am : lc-ld;
      current_highest_am = (current_highest_am > ld) ? current_highest_am : ld;

      /*---------------------------------------------------------------
	Form code and function names for HRR and VRR ordering routines
       ---------------------------------------------------------------*/
      hrr_function_name[14] = am_letter[la-lb];
      hrr_function_name[15] = am_letter[lb];
      hrr_function_name[16] = am_letter[lc-ld];
      hrr_function_name[17] = am_letter[ld];
      vrr_function_name[14] = am_letter[la-lb];
      vrr_function_name[15] = am_letter[lb];
      vrr_function_name[16] = am_letter[lc-ld];
      vrr_function_name[17] = am_letter[ld];
      hrr_code_name[14] = am_letter[la-lb];
      hrr_code_name[15] = am_letter[lb];
      hrr_code_name[16] = am_letter[lc-ld];
      hrr_code_name[17] = am_letter[ld];
      hrr_code = fopen(hrr_code_name,"w");
      vrr_code_name[14] = am_letter[la-lb];
      vrr_code_name[15] = am_letter[lb];
      vrr_code_name[16] = am_letter[lc-ld];
      vrr_code_name[17] = am_letter[ld];
      vrr_code = fopen(vrr_code_name,"w");

      /*-----------------------------------
	Write the overhead to the HRR code
       -----------------------------------*/
      fprintf(hrr_code,"#include <stdio.h>\n");
      fprintf(hrr_code,"#include <string.h>\n");
      fprintf(hrr_code,"#include <libint/libint.h>\n");
      fprintf(hrr_code,"#include \"libr12.h\"\n");
      fprintf(hrr_code,"#include <libint/hrr_header.h>\n\n");
      fprintf(hrr_code,"#include \"r12_hrr_header.h\"\n\n");
      fprintf(hrr_code,"extern void %s(Libr12_t *, prim_data *);\n\n",vrr_function_name);
      fprintf(hrr_code,"  /* Computes (%c%c|%c%c) integrals for linear R12-methods */\n\n",
	      am_letter[la-lb],am_letter[lb],am_letter[lc-ld],am_letter[ld]);
      fprintf(hrr_code,"void %s(Libr12_t *Libr12, int num_prim_comb)\n{\n",hrr_function_name);
      fprintf(hrr_code," prim_data *Data = Libr12->PrimQuartet;\n");
      fprintf(hrr_code," REALTYPE *int_stack = Libr12->int_stack;\n");
      fprintf(hrr_code," int i,j;\n REALTYPE tmp, *target;\n\n");

      /*--------------------------------------------------------------
	Include the function into the hrr_header.h and init_libint.cc
       --------------------------------------------------------------*/
      fprintf(hrr_header,"void %s(Libr12_t *, int);\n",hrr_function_name);
      fprintf(init_code,"  build_r12_grt[%d][%d][%d][%d] = %s;\n",la-lb,lb,lc-ld,ld,hrr_function_name);

      /*-----------------------------------
	Write the overhead to the VRR code
       -----------------------------------*/
      fprintf(vrr_code,"#include <stdio.h>\n");
      fprintf(vrr_code,"#include <libint/libint.h>\n");
      fprintf(vrr_code,"#include \"libr12.h\"\n");
      fprintf(vrr_code,"#include <libint/vrr_header.h>\n");
      fprintf(vrr_code,"#include \"r12_vrr_header.h\"\n\n");
      fprintf(vrr_code,"void %s(Libr12_t *Libr12, prim_data *Data)\n{\n",vrr_function_name);

      
      /*--------------------------------------------------
	Starting at the target node(s) set up an HRR graph.
       --------------------------------------------------*/

      last_hrr_node = 0;
      num_hrr_targets=0;
      for(i=0;i<NUMGRTTYPES;i++) {
	target_hrr_nodes[num_hrr_targets] = last_hrr_node;
	hrr_nodes[last_hrr_node].A = la-lb;
	hrr_nodes[last_hrr_node].B = lb;
	hrr_nodes[last_hrr_node].C = lc-ld;
	hrr_nodes[last_hrr_node].D = ld;
	hrr_nodes[last_hrr_node].m = 0;
	hrr_nodes[last_hrr_node].grt_type = i;
	first_hrr_to_compute = last_hrr_node;
	k = mk_hrr_node(hrr_nodes[last_hrr_node], hrr_nodes, 1);
	if (k == first_hrr_to_compute) { /* If the node hasn't been added to the tree before */
	  if (num_hrr_targets) {
	    hrr_nodes[target_hrr_nodes[num_hrr_targets-1]].llink = first_hrr_to_compute;
	    hrr_nodes[first_hrr_to_compute].rlink = target_hrr_nodes[num_hrr_targets-1];
	    hrr_nodes[first_hrr_to_compute].llink = -1;
	  }
	  else {
	    hrr_nodes[first_hrr_to_compute].rlink = -1;
	    hrr_nodes[first_hrr_to_compute].llink = -1;
	  }
	  num_hrr_targets++;
	}
	/* Mark the class as a target unless it's an (ss|ss) type class */
	if (k >= 0)
	  hrr_nodes[k].target = 1;
      }

	/*-------------------------------------------
	  Traverse the graph starting at each target
	 -------------------------------------------*/
      for(i=0;i<num_hrr_targets;i++) {
	j = target_hrr_nodes[i];
	for(k=0; k<NUMCHILDREN; k++)
	  if(hrr_nodes[j].children[k] >= 0){
	    mark_hrr_parents(hrr_nodes[j].children[k], hrr_nodes, j);
	  }
      }

      
      init_mem(1);

      /*---------------------------------------------------------------
	Allocate and zero out space for classes to be generated by VRR
	Those do not include (ss|ss) and (ss||ss)
       ---------------------------------------------------------------*/
      for(i=last_hrr_node-1;i>=0;i--)
	if (hrr_nodes[i].B == 0 && hrr_nodes[i].D == 0) {
	  hrr_nodes[i].marked = 1;
	  /*--- do not allocate space for (ss|ss) and (ss||ss) ---*/
	  if (hrr_nodes[i].A == 0 && hrr_nodes[i].C == 0 &&
	       (hrr_nodes[i].grt_type == 0 || hrr_nodes[i].grt_type == 1))
	    continue;
	  hrr_nodes[i].pointer = get_mem(hrr_nodes[i].size);
	  fprintf(hrr_code," Libr12->");
	  switch(hrr_nodes[i].grt_type) {
	  case 0:
	      fprintf(hrr_code,"gvrr_classes");
	      break;
	  case 1:
	      fprintf(hrr_code,"rvrr_classes");
	      break;
	  case 2:
	      fprintf(hrr_code,"t1vrr_classes");
	      break;
	  case 3:
	      fprintf(hrr_code,"t2vrr_classes");
	      break;
	  }
	  fprintf(hrr_code,"[%d][%d] = int_stack + %d;\n",
		  hrr_nodes[i].A,hrr_nodes[i].C,hrr_nodes[i].pointer);
	}
      base_mem = get_total_memory();
      fprintf(hrr_code," memset(int_stack,0,%d*sizeof(REALTYPE));\n\n",base_mem);
      fprintf(hrr_code," Libr12->r12vrr_stack = int_stack + %d;\n",base_mem);
      
      
      
      /*----------------------------
	Build the HRR call sequence
       ----------------------------*/
      if (lb != 0 || ld != 0) {
	target_data = alloc_mem_hrr(hrr_nodes);
      }


      hrr_mem = get_total_memory();
      if (max_stack_size < hrr_mem)
	max_stack_size = hrr_mem;
      fprintf(hrr_code," for(i=0;i<num_prim_comb;i++) {\n");
      fprintf(hrr_code,"   %s(Libr12, Data);\n",vrr_function_name);
      fprintf(hrr_code,"   Data++;\n }\n\n");

      
      /*------------------------------------------------
	Evaluate the HRR tree for each class
       ------------------------------------------------*/
      /*--- If we have non-(ss|ss) class - perform the standard procedure ---*/
      j = first_hrr_to_compute;
      if (last_hrr_node > 0)
      do {
	fprintf(hrr_code, " /*--- compute (%c%c|",
		am_letter[hrr_nodes[j].A],am_letter[hrr_nodes[j].B]);
	switch (hrr_nodes[j].grt_type) {
	case 0:
	    break;
	case 1:
	    fprintf(hrr_code,"|");
	    break;
	case 2:
	    fprintf(hrr_code,"[r12,T1]|");
	    break;
	case 3:
	    fprintf(hrr_code,"[r12,T2]|");
	    break;
	}
	fprintf(hrr_code,"%c%c) ---*/\n",
		am_letter[hrr_nodes[j].C],am_letter[hrr_nodes[j].D]);


	if (hrr_nodes[j].B > 0 || hrr_nodes[j].D > 0) {
	  /*--- compute the number of children ---*/
	  num_children = 0;
	  for(i=0;i<NUMCHILDREN;i++)
	    if (hrr_nodes[j].children[i] >= 0)
	      num_children++;
	    
	  if (hrr_nodes[j].B == 0 && hrr_nodes[j].D != 0) {
	    if (hrr_nodes[j].grt_type == 3)
	      fprintf(hrr_code, "   t2hrr3_build_%c%c(Libr12->ShellQuartet.CD,Libr12->ShellQuartet.AC,int_stack+%d,",
		      am_letter[hrr_nodes[j].C], am_letter[hrr_nodes[j].D], hrr_nodes[j].pointer);
	    else
	      fprintf(hrr_code, "   hrr3_build_%c%c(Libr12->ShellQuartet.CD,int_stack+%d,",
		      am_letter[hrr_nodes[j].C], am_letter[hrr_nodes[j].D], hrr_nodes[j].pointer);
	  }
	  else if (hrr_nodes[j].B != 0) {
	    if (hrr_nodes[j].grt_type == 2)
	      fprintf(hrr_code, "   t1hrr1_build_%c%c(Libr12->ShellQuartet.AB,Libr12->ShellQuartet.AC,int_stack+%d,",
		      am_letter[hrr_nodes[j].A], am_letter[hrr_nodes[j].B], hrr_nodes[j].pointer);
	    else
	      fprintf(hrr_code, "   hrr1_build_%c%c(Libr12->ShellQuartet.AB,int_stack+%d,",
		      am_letter[hrr_nodes[j].A], am_letter[hrr_nodes[j].B], hrr_nodes[j].pointer);
	  }
	
	  for(i=0;i<NUMCHILDREN;i++)
	    if (hrr_nodes[j].children[i] >= 0) {
	      child = hrr_nodes[j].children[i];
	      fprintf(hrr_code, " int_stack+%d,", hrr_nodes[child].pointer);
	    }
	    
	  if (hrr_nodes[j].B == 0 && hrr_nodes[j].D != 0) {
	    if (hrr_nodes[j].grt_type == 3)
	      fprintf(hrr_code, " %d, %d);\n", hrr_nodes[j].A,hrr_nodes[j].B);
	    else
	      fprintf(hrr_code, " %d);\n", io[1+hrr_nodes[j].A]*io[1+hrr_nodes[j].B]);
	  }
	  else if (hrr_nodes[j].B != 0) {
	    if (hrr_nodes[j].grt_type == 2)
	      fprintf(hrr_code, " %d, %d);\n", hrr_nodes[j].C,hrr_nodes[j].D);
	    else
	      fprintf(hrr_code, " %d);\n", io[1+hrr_nodes[j].C]*io[1+hrr_nodes[j].D]);
	  }
	}
	
	/* Pass the "target" quartets to CINTS */
        if (hrr_nodes[j].target) {
	  fprintf(hrr_code,"     Libr12->te_ptr[%d] = int_stack + %d;\n",hrr_nodes[j].grt_type,hrr_nodes[j].pointer);
	}
	j = hrr_nodes[j].rlink;
      } while (j != -1);
      
      fprintf(hrr_code,"\n}\n");
      fclose(hrr_code);
      printf("Done with %s\n",hrr_code_name);
      for(i=0;i<last_hrr_node;i++) {
	hrr_nodes[i].llink = -1;
        hrr_nodes[i].rlink = -1;
      }


      /*----------------------------
	Zero out the hashing tables
       ----------------------------*/
      for(i=0;i<NUMGRTTYPES;i++)
	for(j=0;j<2*LMAX_AM;j++)
	  for(k=0;k<2*LMAX_AM;k++)
	    memset(vrr_hash_table[i][j][k],0,(4*LMAX_AM)*sizeof(int));
      for(i=0;i<NUMGRTTYPES;i++)
	for(j=0;j<2*LMAX_AM;j++)
	  for(k=0;k<2*LMAX_AM;k++)
	    for(l=0;l<2*LMAX_AM;l++)
	      memset(hrr_hash_table[i][j][k][l],0,(2*LMAX_AM)*sizeof(int));

      
      /*------------------------------------------------------------------
	Now generate the VRR graph using the (e0|f0) type classes present
	in the HRR graph as "potential" targets
       ------------------------------------------------------------------*/

      vrr_nodes = &(hrr_nodes[last_hrr_node]);
      last_vrr_node = 0;
      for(i=0;i<last_hrr_node;i++)
	if (hrr_nodes[i].B == 0 && hrr_nodes[i].D == 0) {
	  vrr_nodes[last_vrr_node].A = hrr_nodes[i].A;
	  vrr_nodes[last_vrr_node].B = hrr_nodes[i].B;
	  vrr_nodes[last_vrr_node].C = hrr_nodes[i].C;
	  vrr_nodes[last_vrr_node].D = hrr_nodes[i].D;
	  vrr_nodes[last_vrr_node].m = hrr_nodes[i].m;
	  vrr_nodes[last_vrr_node].grt_type = hrr_nodes[i].grt_type;
	  k = mk_vrr_node(vrr_nodes[last_vrr_node], vrr_nodes, 1);
	  if (k >= 0)
	    vrr_nodes[k].target = 1;
	}
      /*--- Now find the true targets (nodes with 0 parents) ---*/
      num_vrr_targets = 0;
      for(i=0;i<last_vrr_node;i++)
	if (vrr_nodes[i].num_parents == 0) {
	  target_vrr_nodes[num_vrr_targets] = i;
	  if (num_vrr_targets > 0) {
	    vrr_nodes[target_vrr_nodes[num_vrr_targets-1]].llink = i;
	    vrr_nodes[i].rlink = target_vrr_nodes[num_vrr_targets-1];
	    vrr_nodes[i].llink = -1;
	  }
	  else {
	    vrr_nodes[i].rlink = -1;
	    vrr_nodes[i].llink = -1;
	  }
	  num_vrr_targets++;
	}
      first_vrr_to_compute = target_vrr_nodes[num_vrr_targets-1];

      /* Traverse the graph starting at each target */
      for(i=0;i<num_vrr_targets;i++) {
	j = target_vrr_nodes[i];
	for(k=0; k<NUMCHILDREN; k++){
	  if(vrr_nodes[j].children[k] >= 0){
	    mark_vrr_parents(vrr_nodes[j].children[k], vrr_nodes, j);
	  }
	}
      }

      init_mem(1);

      /* Build the call sequence */
      target_data = alloc_mem_vrr(vrr_nodes);
      vrr_mem = base_mem + get_total_memory();
      if (max_stack_size < vrr_mem)
	max_stack_size = vrr_mem;
      fprintf(vrr_code," REALTYPE *r12vrr_stack = Libr12->r12vrr_stack;\n");
      fprintf(vrr_code," REALTYPE *tmp, *target_ptr;\n int i, am[2];\n\n");
      
      j = first_vrr_to_compute;
      do {
	fprintf(vrr_code, " /*--- compute (%c%c|",
		am_letter[vrr_nodes[j].A],am_letter[vrr_nodes[j].B]);
	switch (vrr_nodes[j].grt_type) {
	case 0:
	    break;
	case 1:
	    fprintf(vrr_code,"|");
	    break;
	case 2:
	    fprintf(vrr_code,"[r12,T1]|");
	    break;
	case 3:
	    fprintf(vrr_code,"[r12,T2]|");
	    break;
	}
	fprintf(vrr_code,"%c%c)^%d ---*/\n",
		am_letter[vrr_nodes[j].C],am_letter[vrr_nodes[j].D],vrr_nodes[j].m);

	/*---------------------------------------------------------
	  Decide which routine to use to compute the current class
	 ---------------------------------------------------------*/
	switch (vrr_nodes[j].grt_type) {
	case 0:
	    if (vrr_nodes[j].A <= LIBINT_OPT_AM && vrr_nodes[j].C <= LIBINT_OPT_AM)
	      fprintf(vrr_code, " _BUILD_%c0%c0(Data,", am_letter[vrr_nodes[j].A], am_letter[vrr_nodes[j].C]);
	    else {
	      fprintf(vrr_code, " am[0] = %d;  am[1] = %d;\n", vrr_nodes[j].A, vrr_nodes[j].C);
	      fprintf(vrr_code, " vrr_build_xxxx(am,Data,");
	    }
	    break;
	case 1:
	    if (vrr_nodes[j].A <= opt_am && vrr_nodes[j].C <= opt_am)
	      fprintf(vrr_code, " _R_BUILD_%c0%c0(Data,", am_letter[vrr_nodes[j].A], am_letter[vrr_nodes[j].C]);
	    else {
	      fprintf(vrr_code, " am[0] = %d;  am[1] = %d;\n", vrr_nodes[j].A, vrr_nodes[j].C);
	      fprintf(vrr_code, " r_vrr_build_xxxx(am,Data,");
	    }
	    break;
	case 2:
	    if (vrr_nodes[j].A <= opt_am && vrr_nodes[j].C <= opt_am)
	      fprintf(vrr_code, " _T1_BUILD_%c0%c0(Data,", am_letter[vrr_nodes[j].A], am_letter[vrr_nodes[j].C]);
	    else {
	      fprintf(vrr_code, " am[0] = %d;  am[1] = %d;\n", vrr_nodes[j].A, vrr_nodes[j].C);
	      fprintf(vrr_code, " t1_vrr_build_xxxx(am,Data,");
	    }
	    break;
	case 3:
	    if (vrr_nodes[j].A <= opt_am && vrr_nodes[j].C <= opt_am)
	      fprintf(vrr_code, " _T2_BUILD_%c0%c0(Data,", am_letter[vrr_nodes[j].A], am_letter[vrr_nodes[j].C]);
	    else {
	      fprintf(vrr_code, " am[0] = %d;  am[1] = %d;\n", vrr_nodes[j].A, vrr_nodes[j].C);
	      fprintf(vrr_code, " t2_vrr_build_xxxx(am,Data,");
	    }
	    break;
	}
	switch (vrr_nodes[j].grt_type) {
	case 2:
	    fprintf(vrr_code,"&(Libr12->ShellQuartet),");
	    break;
	case 3:
	    fprintf(vrr_code,"&(Libr12->ShellQuartet),");
	    break;
	}
	fprintf(vrr_code, "r12vrr_stack+%d", vrr_nodes[j].pointer);
	num_children = (vrr_nodes[j].grt_type == 1) ? 6 : 5;
	for(k=0; k<num_children; k++){
	    if(vrr_nodes[j].children[k] >= 0) /*--- this child is a "real" class ---*/
	      fprintf(vrr_code, ", r12vrr_stack+%d", vrr_nodes[vrr_nodes[j].children[k]].pointer);
	    else if (vrr_nodes[j].children[k] == NONODE) /*--- this child is a no-class ---*/
	      fprintf(vrr_code, ", NULL");
	    else if (vrr_nodes[j].children[k] == SSR12SSNODE) /*--- this is a (ss||ss) ---*/
	      fprintf(vrr_code, ", &(Data->ss_r12_ss)");
	    else /*--- this is a (ss|ss)^m ---*/
	      fprintf(vrr_code, ", Data->F+%d", (-1)*(vrr_nodes[j].children[k] + MOFFSET));
	}
	fprintf(vrr_code, ");\n");

	/*-----------------------------------------------
	  If this derivative class is one of the targets
	  copy it to a location pointed by deriv_classes
	  to be used by the calling hrr_order routine
	 -----------------------------------------------*/
	if (vrr_nodes[j].target == 1) {
	  fprintf(vrr_code, " tmp = r12vrr_stack + %d;\n", vrr_nodes[j].pointer);
	  fprintf(vrr_code, " target_ptr = Libr12->");
	  switch (vrr_nodes[j].grt_type) {
	  case 0:
	      fprintf(vrr_code,"g");
	      break;
	  case 1:
	      fprintf(vrr_code,"r");
	      break;
	  case 2:
	      fprintf(vrr_code,"t1");
	      break;
	  case 3:
	      fprintf(vrr_code,"t2");
	      break;
	  }
	  fprintf(vrr_code,"vrr_classes[%d][%d];\n",
		  vrr_nodes[j].A,vrr_nodes[j].C);
	  fprintf(vrr_code, " for(i=0;i<%d;i++)\n",vrr_nodes[j].size);
	  fprintf(vrr_code, "   target_ptr[i] += tmp[i];\n\n");
	}
	else
	  fprintf(vrr_code, "\n");

	j = vrr_nodes[j].rlink;
      } while (j != -1);
      fprintf(vrr_code, "\n}\n\n");
      fclose(vrr_code);
      printf("Done with %s\n",vrr_code_name);
      for(i=0;i<last_vrr_node;i++) {
	vrr_nodes[i].llink = -1;
        vrr_nodes[i].rlink = -1;
      }

      /* compare this max_stack_size to the libint_stack_size for this angular momentum */
      if (libr12_stack_size[current_highest_am] < max_stack_size)
	libr12_stack_size[current_highest_am] = max_stack_size;

      max_stack_size = 0;
      
      }
    }
    }
  }
}


/* Recursive function that build a hybrid HRR subgraph given the parent */

int mk_hrr_node(class node, class *allnodes, int new)
{

  int i, j, k, l;
  class O[NUMCHILDREN];
  int subnodes = 0;
  int thisnode;
  int rlink, llink;
  int made = 0;

/*  if (node.A == 0 && node.B == 0 && node.C == 0 && node.D == 0)
    return -1;*/

  /* Search for the parent node on stack
     If it's not there - we'll add it to the end of the stack */
  thisnode = last_hrr_node;
  /* it's already placed on the stack allnodes - make sure children don't get created again (made = 1) */
  if (hrr_hash_table[node.grt_type][node.A][node.B][node.C][node.D]) {
    i = hrr_hash_table[node.grt_type][node.A][node.B][node.C][node.D] - 1;
    thisnode = i;
    made = 1;
  }

  /* it's not computed, add it, and make it the first to compute! */
  if(!made){
    allnodes[thisnode].A = node.A;
    allnodes[thisnode].B = node.B;
    allnodes[thisnode].C = node.C;
    allnodes[thisnode].D = node.D;
    allnodes[thisnode].m = node.m;
    allnodes[thisnode].grt_type = node.grt_type;
    hrr_hash_table[node.grt_type][node.A][node.B][node.C][node.D] = thisnode + 1;
    allnodes[thisnode].num_parents = 0;
    allnodes[thisnode].parents_counter = 0;
    allnodes[thisnode].marked = 0;
    allnodes[thisnode].pointer = 0;
    memset(allnodes[thisnode].parents,0,NUMPARENTS*sizeof(int));
    for(i=0;i<NUMCHILDREN;i++)
      allnodes[thisnode].children[i] = NONODE;
    allnodes[thisnode].size = io[1+node.A]*io[1+node.B]*io[1+node.C]*io[1+node.D];
    allnodes[thisnode].target = 0;
    /* We just added a node ..*/
    last_hrr_node++;
    /* If stack is overfull - exit */
    if(last_hrr_node >= MAXNODE) {
      printf(" Maximum stack size is reached. Change MAXNODE and recompile.\n\n");
      exit(1);
    }
  }

  /* If the parent class wasn't on stack already (!new) - increase the parent counter */
  if(!new){
    allnodes[thisnode].num_parents++;
    allnodes[thisnode].parents_counter++;
    if (allnodes[thisnode].num_parents > NUMPARENTS) {
      printf("Number of parents exceeds the limit\n");
      exit(1);
    }
  }


  /* now make all child nodes */
  if (!made) {
    if(node.B){
      O[0].A = node.A+1;
      O[0].B = node.B-1;
      O[0].C = node.C;
      O[0].D = node.D;
      O[0].m = node.m;
      O[0].grt_type = node.grt_type;
      allnodes[thisnode].children[0] = mk_hrr_node(O[0], allnodes, made);
      O[1].A = node.A;
      O[1].B = node.B-1;
      O[1].C = node.C;
      O[1].D = node.D;
      O[1].m = node.m;
      O[1].grt_type = node.grt_type;
      allnodes[thisnode].children[1] = mk_hrr_node(O[1], allnodes, made);
      /* Special case - [r12.T1] */
      if (node.grt_type == 2) {
	O[2].A = node.A+1;
	O[2].B = node.B-1;
	O[2].C = node.C;
	O[2].D = node.D;
	O[2].m = node.m;
	O[2].grt_type = 0;
	allnodes[thisnode].children[2] = mk_hrr_node(O[2], allnodes, made);
        O[3].A = node.A;
	O[3].B = node.B-1;
	O[3].C = node.C+1;
	O[3].D = node.D;
	O[3].m = node.m;
	O[3].grt_type = 0;
	allnodes[thisnode].children[3] = mk_hrr_node(O[3], allnodes, made);
	O[4].A = node.A;
	O[4].B = node.B-1;
	O[4].C = node.C;
	O[4].D = node.D;
	O[4].m = node.m;
	O[4].grt_type = 0;
	allnodes[thisnode].children[4] = mk_hrr_node(O[4], allnodes, made);
      }
    }
    else if(node.D){
      O[0].A = node.A;
      O[0].B = node.B;
      O[0].C = node.C+1;
      O[0].D = node.D-1;
      O[0].m = node.m;
      O[0].grt_type = node.grt_type;
      allnodes[thisnode].children[0] = mk_hrr_node(O[0], allnodes, made);
      O[1].A = node.A;
      O[1].B = node.B;
      O[1].C = node.C;
      O[1].D = node.D-1;
      O[1].m = node.m;
      O[1].grt_type = node.grt_type;
      allnodes[thisnode].children[1] = mk_hrr_node(O[1], allnodes, made);
      /* Special case - [r12.T2] */
      if (node.grt_type == 3) {
	O[2].A = node.A;
	O[2].B = node.B;
	O[2].C = node.C+1;
	O[2].D = node.D-1;
	O[2].m = node.m;
	O[2].grt_type = 0;
	allnodes[thisnode].children[2] = mk_hrr_node(O[2], allnodes, made);
        O[3].A = node.A+1;
	O[3].B = node.B;
	O[3].C = node.C;
	O[3].D = node.D-1;
	O[3].m = node.m;
	O[3].grt_type = 0;
	allnodes[thisnode].children[3] = mk_hrr_node(O[3], allnodes, made);
	O[4].A = node.A;
	O[4].B = node.B;
	O[4].C = node.C;
	O[4].D = node.D-1;
	O[4].m = node.m;
	O[4].grt_type = 0;
	allnodes[thisnode].children[4] = mk_hrr_node(O[4], allnodes, made);
      }
    }
  }

  return thisnode;

}



/* Recursive function that builds a hybrid VRR subgraph given the parent */

int mk_vrr_node(class node, class *allnodes, int new)
{

  int i, j, k, l;
  class O[NUMCHILDREN];
  int subnodes = 0;
  int thisnode;
  int made = 0;

  /* If it's not a derivative class - do some checks to see if need to proceed */
  if (node.grt_type == 0 && node.A + node.B + node.C + node.D == 0)
    return (-1)*node.m - MOFFSET;
  else if (node.grt_type == 1 && node.A + node.B + node.C + node.D == 0)
    return SSR12SSNODE;

  /* Search for the parent node on stack
     If it's not there - we'll add it to the end of the stack */
  thisnode = last_vrr_node;
  /* it's already placed on the stack allnodes - make sure children don't get created again (made = 1) */
  if (vrr_hash_table[node.grt_type][node.A][node.C][node.m]) {
    i = vrr_hash_table[node.grt_type][node.A][node.C][node.m] - 1;
    thisnode = i;
    made = 1;
  }

  /* it's not computed, add it, and make it the first to compute! */
  if(!made){
    allnodes[thisnode].A = node.A;
    allnodes[thisnode].B = node.B;
    allnodes[thisnode].C = node.C;
    allnodes[thisnode].D = node.D;
    allnodes[thisnode].m = node.m;
    allnodes[thisnode].grt_type = node.grt_type;
    vrr_hash_table[node.grt_type][node.A][node.C][node.m] = thisnode + 1;
    allnodes[thisnode].num_parents = 0;
    allnodes[thisnode].parents_counter = 0;
    allnodes[thisnode].marked = 0;
    allnodes[thisnode].target = 0;
    allnodes[thisnode].pointer = 0;
    memset(allnodes[thisnode].parents,0,NUMPARENTS*sizeof(int));
    for(i=0;i<NUMCHILDREN;i++)
      allnodes[thisnode].children[i] = NONODE;
    allnodes[thisnode].size = io[1+node.A]*io[1+node.B]*io[1+node.C]*io[1+node.D];
    allnodes[thisnode].llink = -1;
    allnodes[thisnode].rlink = -1;
    /* We just added a node ..*/
    last_vrr_node++;
    /* If stack is overfull - exit */
    if(last_vrr_node+last_hrr_node >= MAXNODE) {
      printf(" Maximum stack size is reached. Change MAXNODE and recompile.\n\n");
      exit(1);
    }
  }

  /* If the parent class wasn't on stack already (!new) - increase the parent counter */
  if(!new){
    allnodes[thisnode].num_parents++;
    allnodes[thisnode].parents_counter++;
    if (allnodes[thisnode].num_parents > NUMPARENTS) {
      printf("Number of parents exceeds the limit\n");
      exit(1);
    }
  }


  /* now make all child nodes */
  if (!made) {
    /* regular ERI */
    if (node.grt_type == 0) {
      if(node.A){              /*--- (a0|c0 ---*/
    	O[0].A = node.A-1;
	O[0].B = 0;
	O[0].C = node.C;
	O[0].D = 0;
	O[0].m = node.m;
	O[0].grt_type = 0;
	allnodes[thisnode].children[0] = mk_vrr_node(O[0], allnodes, made);
	O[1].A = node.A-1;
	O[1].B = 0;
	O[1].C = node.C;
	O[1].D = 0;
	O[1].m = node.m+1;
	O[1].grt_type = 0;
	allnodes[thisnode].children[1] = mk_vrr_node(O[1], allnodes, made);
	if(node.A>1){
	  O[2].A = node.A-2;
	  O[2].B = 0;
	  O[2].C = node.C;
	  O[2].D = 0;
	  O[2].m = node.m;
	  O[2].grt_type = 0;
	  allnodes[thisnode].children[2] = mk_vrr_node(O[2], allnodes, made);
	  O[3].A = node.A-2;
	  O[3].B = 0;
	  O[3].C = node.C;
	  O[3].D = 0;
	  O[3].m = node.m+1;
	  O[3].grt_type = 0;
	  allnodes[thisnode].children[3] = mk_vrr_node(O[3], allnodes, made);
	}
	if(node.C){
	  O[4].A = node.A-1;
	  O[4].B = 0;
	  O[4].C = node.C-1;
	  O[4].D = 0;
	  O[4].m = node.m+1;
	  O[4].grt_type = 0;
	  allnodes[thisnode].children[4] = mk_vrr_node(O[4], allnodes, made);
	}
      }
      else if(node.C){     /*--- (00|c0) ---*/
	O[0].A = node.A;
	O[0].B = 0;
	O[0].C = node.C-1;
	O[0].D = 0;
	O[0].m = node.m;
	O[0].grt_type = 0;
	allnodes[thisnode].children[0] = mk_vrr_node(O[0], allnodes, made);
	O[1].A = node.A;
	O[1].B = 0;
	O[1].C = node.C-1;
	O[1].D = 0;
	O[1].m = node.m+1;
	O[1].grt_type = 0;
	allnodes[thisnode].children[1] = mk_vrr_node(O[1], allnodes, made);
	if(node.C>1){
	  O[2].A = node.A;
	  O[2].B = 0;
	  O[2].C = node.C-2;
	  O[2].D = 0;
	  O[2].m = node.m;
	  O[2].grt_type = 0;
	  allnodes[thisnode].children[2] = mk_vrr_node(O[2], allnodes, made);
	  O[3].A = node.A;
	  O[3].B = 0;
	  O[3].C = node.C-2;
	  O[3].D = 0;
	  O[3].m = node.m+1;
	  O[3].grt_type = 0;
	  allnodes[thisnode].children[3] = mk_vrr_node(O[3], allnodes, made);
	}
      }
    }
    /* Integral of r12 */
    else if (node.grt_type == 1) {
      if(node.A){        /*--- (a0||c0) ---*/
	O[0].A = node.A-1;
	O[0].B = 0;
	O[0].C = node.C;
	O[0].D = 0;
	O[0].m = 0;
	O[0].grt_type = 1;
	allnodes[thisnode].children[0] = mk_vrr_node(O[0], allnodes, made);
	if (node.A > 1) {
	  O[1].A = node.A-2;
	  O[1].B = 0;
	  O[1].C = node.C;
	  O[1].D = 0;
	  O[1].m = 0;
	  O[1].grt_type = 1;
	  allnodes[thisnode].children[1] = mk_vrr_node(O[1], allnodes, made);
	}
	O[2].A = node.A;
	O[2].B = 0;
	O[2].C = node.C;
	O[2].D = 0;
	O[2].m = 0;
	O[2].grt_type = 0;
	allnodes[thisnode].children[2] = mk_vrr_node(O[2], allnodes, made);
	O[3].A = node.A-1;
	O[3].B = 0;
	O[3].C = node.C;
	O[3].D = 0;
	O[3].m = 0;
	O[3].grt_type = 0;
	allnodes[thisnode].children[3] = mk_vrr_node(O[3], allnodes, made);
	if (node.A > 1) {
	  O[4].A = node.A-2;
	  O[4].B = 0;
	  O[4].C = node.C;
	  O[4].D = 0;
	  O[4].m = 0;
	  O[4].grt_type = 0;
	  allnodes[thisnode].children[4] = mk_vrr_node(O[4], allnodes, made);
	}
	if(node.C){
	  O[5].A = node.A-1;
	  O[5].B = 0;
	  O[5].C = node.C-1;
	  O[5].D = 0;
	  O[5].m = 0;
	  O[5].grt_type = 0;
	  allnodes[thisnode].children[5] = mk_vrr_node(O[5], allnodes, made);
	}
      }
      else if (node.C){         /*--- (00||c0) ---*/
	O[0].A = node.A;
	O[0].B = 0;
	O[0].C = node.C-1;
	O[0].D = 0;
	O[0].m = 0;
	O[0].grt_type = 1;
	allnodes[thisnode].children[0] = mk_vrr_node(O[0], allnodes, made);
	if (node.C > 1) {
	  O[1].A = node.A;
	  O[1].B = 0;
	  O[1].C = node.C-2;
	  O[1].D = 0;
	  O[1].m = 0;
	  O[1].grt_type = 1;
	  allnodes[thisnode].children[1] = mk_vrr_node(O[1], allnodes, made);
	}
	O[2].A = node.A;
	O[2].B = 0;
	O[2].C = node.C;
	O[2].D = 0;
	O[2].m = 0;
	O[2].grt_type = 0;
	allnodes[thisnode].children[2] = mk_vrr_node(O[2], allnodes, made);
	O[3].A = node.A;
	O[3].B = 0;
	O[3].C = node.C-1;
	O[3].D = 0;
	O[3].m = 0;
	O[3].grt_type = 0;
	allnodes[thisnode].children[3] = mk_vrr_node(O[3], allnodes, made);
	if(node.C > 1){
	  O[4].A = node.A;
	  O[4].B = 0;
	  O[4].C = node.C-2;
	  O[4].D = 0;
	  O[4].m = 0;
	  O[4].grt_type = 0;
	  allnodes[thisnode].children[4] = mk_vrr_node(O[4], allnodes, made);
	}
      }
    }
    /*--- Integral of [r12,T1] ---*/
    else if (node.grt_type == 2) {
      O[0].A = node.A;
      O[0].B = 0;
      O[0].C = node.C;
      O[0].D = 0;
      O[0].m = 0;
      O[0].grt_type = 0;
      allnodes[thisnode].children[0] = mk_vrr_node(O[0], allnodes, made);
      O[1].A = node.A+1;
      O[1].B = 0;
      O[1].C = node.C;
      O[1].D = 0;
      O[1].m = 0;
      O[1].grt_type = 0;
      allnodes[thisnode].children[1] = mk_vrr_node(O[1], allnodes, made);
      O[2].A = node.A;
      O[2].B = 0;
      O[2].C = node.C+1;
      O[2].D = 0;
      O[2].m = 0;
      O[2].grt_type = 0;
      allnodes[thisnode].children[2] = mk_vrr_node(O[2], allnodes, made);
      if (node.A) {
	O[3].A = node.A-1;
	O[3].B = 0;
	O[3].C = node.C+1;
	O[3].D = 0;
	O[3].m = 0;
	O[3].grt_type = 0;
	allnodes[thisnode].children[3] = mk_vrr_node(O[3], allnodes, made);
	O[4].A = node.A-1;
	O[4].B = 0;
	O[4].C = node.C;
	O[4].D = 0;
	O[4].m = 0;
	O[4].grt_type = 0;
	allnodes[thisnode].children[4] = mk_vrr_node(O[4], allnodes, made);
      }
    }
    /*--- Integrals of [r12,T2] ---*/
    else if (node.grt_type == 3) {
      O[0].A = node.A;
      O[0].B = 0;
      O[0].C = node.C;
      O[0].D = 0;
      O[0].m = 0;
      O[0].grt_type = 0;
      allnodes[thisnode].children[0] = mk_vrr_node(O[0], allnodes, made);
      O[1].A = node.A;
      O[1].B = 0;
      O[1].C = node.C+1;
      O[1].D = 0;
      O[1].m = 0;
      O[1].grt_type = 0;
      allnodes[thisnode].children[1] = mk_vrr_node(O[1], allnodes, made);
      O[2].A = node.A+1;
      O[2].B = 0;
      O[2].C = node.C;
      O[2].D = 0;
      O[2].m = 0;
      O[2].grt_type = 0;
      allnodes[thisnode].children[2] = mk_vrr_node(O[2], allnodes, made);
      if (node.C > 0) {
	O[3].A = node.A+1;
	O[3].B = 0;
	O[3].C = node.C-1;
	O[3].D = 0;
	O[3].m = 0;
	O[3].grt_type = 0;
	allnodes[thisnode].children[3] = mk_vrr_node(O[3], allnodes, made);
	O[4].A = node.A;
	O[4].B = 0;
	O[4].C = node.C-1;
	O[4].D = 0;
	O[4].m = 0;
	O[4].grt_type = 0;
	allnodes[thisnode].children[4] = mk_vrr_node(O[4], allnodes, made);
      }
    }
  }

  return thisnode;

}




/* Make hrr_nodes[rent] a parent of hrr_nodes[n] and proceed recursively */

void mark_hrr_parents(int n, class *allnodes, int rent)
{
  int i;
  int *tmp;

  /* handle case where it's in the parent list already */
  for(i=allnodes[n].num_parents-1; i>=allnodes[n].parents_counter; i--)
    if(rent==allnodes[n].parents[i]) return;

  /* if the parent rent is not in the list - add it to the list! */
  i = --allnodes[n].parents_counter;
  allnodes[n].parents[i] = rent;
  /* hits from all of the parents has been received - schedule it for computation and mark all of its children */
  if (i == 0 && (allnodes[n].B != 0 || allnodes[n].D != 0)) {
    /*--- take it out of the list if it's in there already ---*/
    if (allnodes[n].llink != -1) {
      allnodes[allnodes[n].llink].rlink = allnodes[n].rlink;
    }
    if (allnodes[n].rlink != -1) {
      allnodes[allnodes[n].rlink].llink = allnodes[n].llink;
    }
    /*--- put it in the beginning ---*/
    allnodes[n].llink = -1;
    allnodes[n].rlink = first_hrr_to_compute;
    allnodes[first_hrr_to_compute].llink = n;
    first_hrr_to_compute = n;

    for(i=0; i<NUMCHILDREN; i++)
      if(allnodes[n].children[i] >= 0)
	mark_hrr_parents(allnodes[n].children[i], allnodes, n);
  }
  return;
}


/* Make vrr_nodes[rent] a parent of vrr_nodes[n] and proceed recursively */

void mark_vrr_parents(int n, class *allnodes, int rent)
{
  int i;
  int *tmp;

  /* handle case where it's in there already */
  for(i=allnodes[n].num_parents-1; i>=allnodes[n].parents_counter; i--)
    if(rent==allnodes[n].parents[i]) return;


  /* if the parent rent is not in the list - add it to the list! */
  i = --allnodes[n].parents_counter;
  allnodes[n].parents[i] = rent;
  /* hits from all of the parents has been received - schedule it for computation and mark all of its children */
  if (i == 0) {
    /*--- take it out of the list if it's in there already ---*/
    if (allnodes[n].llink != -1) {
      allnodes[allnodes[n].llink].rlink = allnodes[n].rlink;
    }
    if (allnodes[n].rlink != -1) {
      allnodes[allnodes[n].rlink].llink = allnodes[n].llink;
    }
    /*--- put it in the beginning ---*/
    allnodes[n].llink = -1;
    allnodes[n].rlink = first_vrr_to_compute;
    if (first_vrr_to_compute >= 0)
      allnodes[first_vrr_to_compute].llink = n;
    first_vrr_to_compute = n;

    for(i=0; i<NUMCHILDREN; i++)
      if(allnodes[n].children[i] >= 0)
	mark_vrr_parents(allnodes[n].children[i], allnodes, n);

  }
  return;
}



/* This functions controls memory placement of computed classes on the CINTS stack */

int alloc_mem_hrr(class *nodes)
{
  int i, j, k, l;
  int size;
  int child;
  int free_it;

  j = first_hrr_to_compute;
  do{
    /* Node to compute */
    if (nodes[j].marked == 0) {
      nodes[j].marked = 1; /* sign that it has been passed */
      nodes[j].pointer = get_mem(nodes[j].size); /* Allocate memory for it on a CINTS-provided stack */
    }
      
    /* Figure out which children can be freed,
       i.e. which children are not targets and have all parents marked */
    for(k=0; k<NUMCHILDREN; k++){
      child = nodes[j].children[k];
      if(child >= 0)
	if (nodes[child].target == 0) {
	  free_it = 1;
	  for(l=0; l<nodes[child].num_parents; l++)
	    if(!nodes[nodes[child].parents[l]].marked)
	      free_it = 0;
	  if(free_it)
	    free_mem(nodes[child].pointer, nodes[child].size);
	}
    }
    j = nodes[j].rlink;
  } while (j != -1);
  
  return nodes[0].pointer;
}


/* This functions controls memory placement of computed classes on the CINTS stack */

int alloc_mem_vrr(class *nodes)
{
  int i, j, k, l;
  int size;
  int child;
  int free_it;

  /* Mark all nodes as not computed */
  for(i=0; i<last_vrr_node; i++)
    nodes[i].marked = 0;

  j = first_vrr_to_compute;
  do{
    /* Node to compute */
    nodes[j].marked = 1; /* sign that it has been passed */
    nodes[j].pointer = get_mem(nodes[j].size); /* Allocate memory for it on a CINTS-provided stack */

    /* Figure out which children can be freed,
       i.e. which children have all parents marked */
    for(k=0; k<NUMCHILDREN; k++){
      child = nodes[j].children[k];
      if(child >= 0){
        free_it = 1;
        for(l=0; l<nodes[child].num_parents; l++){
          if(!nodes[nodes[child].parents[l]].marked) {
            free_it = 0;
	  }
	}
        if(free_it){
          free_mem(nodes[child].pointer, nodes[child].size);
	}
      }
    }
    j = nodes[j].rlink;
  } while (j != -1);
  
  return nodes[0].pointer;
}

