/*! \file
    \ingroup INT
    \brief Enter brief description of file here 
*/
/*------------------------------------------------------------------------------------------------------
  Builds a library of functions-calls for applying HRR and VRR
 ------------------------------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mem_man.h"
#include "build_libint.h"
#include <libint/constants.h>
#define MAXNODE 20000
#define NONODE -1000000
#define SSSSNODE -1111

static int last_hrr_node = 0;      /* Global pointer to the last node on the HRR stack */
static int last_vrr_node = 0;      /* Global pointer to the last node on the VRR stack */

extern FILE *outfile, *hrr_header, *init_code;
extern int libint_stack_size[MAX_AM/2+1];
extern LibintParams_t Params;

typedef struct hrr_node{
  int A, B, C, D;               /* Angular momenta on centers A and C */
  int size;               /* Class size in double words */
  int pointer;
  int children[2];        /* Up to 5 children of the class */
  int parents_counter;
  int num_parents;        /* Number of parents */
  int parents[5];         /* Pointers to parents */
  int llink;              /* Pointer to a class computed right before computing this one */
  int rlink;              /* Pointer to a class to be computed after this one is */
  int marked;             /* Flag indicating that this node has been computed */
  int target;             /* Flag indicating that this node is among targets */
  } hrr_class;

typedef struct vrr_node{
  int A, C;               /* Angular momenta on centers A and C */
  int m;                  /* Auxiliary index */
  int size;               /* Class size in double words */
  int pointer;
  int children[5];        /* Up to 5 children of the class */
  int parents_counter;
  int num_parents;        /* Number of parents */
  int parents[9];         /* Pointers to parents */
  int llink;              /* Pointer to a class computed right before computing this one */
  int rlink;              /* Pointer to a class to be computed after this one is */
  int marked;
  int target;
  } vrr_class;


static int first_hrr_to_compute = 0; /* Number of the first class to be computed
				    (pointer to the beginning of the linked list) */
static int first_vrr_to_compute = 0; /* Number of the first class to be computed
				    (pointer to the beginning of the linked list) */

static int hrr_hash_table[MAX_AM+1][MAX_AM+1][MAX_AM+1][MAX_AM+1];
static int vrr_hash_table[MAX_AM+1][MAX_AM+1][2*MAX_AM+1];

void emit_order()
{
  int old_am = Params.old_am;
  int new_am = Params.new_am;
  int opt_am = Params.opt_am;
  int am_to_inline_into_hrr = Params.max_am_to_inline_vrr_manager;
  int am_to_inline_vrr = Params.max_am_manager_to_inline_vrr_worker;
  int am_to_inline_hrr = Params.max_am_manager_to_inline_hrr_worker;
  int to_inline_into_hrr, to_inline_vrr, to_inline_hrr;

  int i, j, k, l;
  int la, lc, lc_min, ld, ld_max, ld_min;
  int lb, lb_min, lb_max;
  int current_highest_am, max_node_am;
  int base_mem, hrr_mem, vrr_mem;
  hrr_class hrr_nodes[MAXNODE];    /* Stack of HRR nodes */
  vrr_class vrr_nodes[MAXNODE];    /* Stack of VRR nodes */
  int target_data;
  int done;
  int max_stack_size = 0;
  int num_vrr_targets;
  int target_vrr_nodes[MAX_AM*MAX_AM];
  char hrr_code_name[] = "hrr_order_0000.cc";
  char hrr_function_name[] = "hrr_order_0000";
  char vrr_code_name[] = "vrr_order_0000.cc";
  char vrr_function_name[] = "vrr_order_0000";
  char inline_vrr_list_name[] = "inline_vrr_order_0000.h";
  char inline_hrr_list_name[] = "inline_hrr_order_0000.h";
  FILE *hrr_code, *vrr_code, *inline_vrr_list, *inline_hrr_list;

  for(la=0;la<=new_am;la++) {
    lb_max = la/2;
    lb_min = (la <= new_am/2) ? 0 : la - new_am/2;
    lc_min = (la == 0) ? 1 : la;
    for(lb=lb_min;lb<=lb_max;lb++) {
    for(lc=lc_min;lc<=new_am;lc++) {
      ld_max = lc/2;
      ld_min = (lc <= new_am/2) ? 0 : lc - new_am/2;
      for(ld=ld_min;ld<=ld_max;ld++) {

      current_highest_am = (la-lb > lb) ? la-lb : lb;
      current_highest_am = (current_highest_am > lc-ld) ? current_highest_am : lc-ld;
      current_highest_am = (current_highest_am > ld) ? current_highest_am : ld;
      to_inline_into_hrr = (current_highest_am <= am_to_inline_into_hrr) ? 1 : 0;
      to_inline_vrr = (current_highest_am <= am_to_inline_vrr) ? 1 : 0;
      to_inline_hrr = (current_highest_am <= am_to_inline_hrr) ? 1 : 0;

      /*---------------------------------------------------------------
	Form code and function names for HRR and VRR ordering routines
       ---------------------------------------------------------------*/
      sprintf(hrr_function_name,"hrr_order_%c%c%c%c",
	      am_letter[la-lb],am_letter[lb],
	      am_letter[lc-ld],am_letter[ld]);
      sprintf(vrr_function_name,"vrr_order_%c%c%c%c",
	      am_letter[la-lb],am_letter[lb],
	      am_letter[lc-ld],am_letter[ld]);
      sprintf(hrr_code_name,"%s.cc",hrr_function_name);
      if (to_inline_into_hrr)
	sprintf(vrr_code_name,"%s.h",vrr_function_name);
      else
	sprintf(vrr_code_name,"%s.cc",vrr_function_name);
      sprintf(inline_vrr_list_name,"inline_vrr_order_%c%c%c%c.h",
	      am_letter[la-lb],am_letter[lb],
	      am_letter[lc-ld],am_letter[ld]);
      sprintf(inline_hrr_list_name,"inline_hrr_order_%c%c%c%c.h",
	      am_letter[la-lb],am_letter[lb],
	      am_letter[lc-ld],am_letter[ld]);
      hrr_code = fopen(hrr_code_name,"w");
      vrr_code = fopen(vrr_code_name,"w");
      inline_vrr_list = fopen(inline_vrr_list_name,"w");
      inline_hrr_list = fopen(inline_hrr_list_name,"w");

      /*-----------------------------------
	Write the overhead to the HRR code
       -----------------------------------*/
      fprintf(hrr_code,"#include <stdio.h>\n");
      fprintf(hrr_code,"#include <string.h>\n");
      fprintf(hrr_code,"#include <psi4-dec.h>\n");
      fprintf(hrr_code,"#include \"libint.h\"\n");
      if (to_inline_hrr) {
	fprintf(hrr_code,"#define INLINE_HRR_WORKER\n");
	fprintf(hrr_code,"#include \"%s\"\n",inline_hrr_list_name);
      }
      fprintf(hrr_code,"#include \"hrr_header.h\"\n\n");
      if (to_inline_into_hrr)
	fprintf(hrr_code,"#include \"%s\"\n",vrr_code_name);
      else
	fprintf(hrr_code,"extern void vrr_order_%c%c%c%c(Libint_t*, prim_data*);\n\n",
		am_letter[la-lb],am_letter[lb],am_letter[lc-ld],am_letter[ld]);
      fprintf(hrr_code,"  /* Computes quartets of (%c%c|%c%c) integrals */\n\n",
	      am_letter[la-lb],am_letter[lb],am_letter[lc-ld],am_letter[ld]);
      fprintf(hrr_code,"REALTYPE *%s(Libint_t *Libint, int num_prim_comb)\n{\n",hrr_function_name);
      fprintf(hrr_code," prim_data *Data = Libint->PrimQuartet;\n");
      fprintf(hrr_code," REALTYPE *int_stack = Libint->int_stack;\n");
      fprintf(hrr_code," int i;\n\n");
      
      /*-------------------------------------------------------------
	Include the function into the hrr_header.h and init_libint.c
       -------------------------------------------------------------*/
      fprintf(hrr_header,"REALTYPE *%s(Libint_t *, int);\n",hrr_function_name);
      fprintf(init_code,"  build_eri[%d][%d][%d][%d] = %s;\n",la-lb,lb,lc-ld,ld,hrr_function_name);

      /*-----------------------------------
	Write the overhead to the VRR code
       -----------------------------------*/
      fprintf(vrr_code,"#include <stdio.h>\n");
      fprintf(vrr_code,"#include <psi4-dec.h>\n");
      fprintf(vrr_code,"#include \"libint.h\"\n");
      if (to_inline_vrr) {
	fprintf(vrr_code,"#define INLINE_VRR_WORKER\n");
	fprintf(vrr_code,"#include \"%s\"\n",inline_vrr_list_name);
      }
      fprintf(vrr_code,"#include \"vrr_header.h\"\n\n");
      fprintf(vrr_code,"  /* Computes quartets necessary to compute (%c%c|%c%c) integrals */\n\n",
	      am_letter[la-lb],am_letter[lb],am_letter[lc-ld],am_letter[ld]);
      if (to_inline_into_hrr)
	fprintf(vrr_code,"inline ");
      fprintf(vrr_code,"void %s(Libint_t * Libint, prim_data *Data)\n{\n",vrr_function_name);

      /*----------------------------
	Zero out the hashing tables
       ----------------------------*/
      for(i=0;i<=MAX_AM;i++)
	for(j=0;j<=MAX_AM;j++) {
	  memset(vrr_hash_table[i][j],0,(2*MAX_AM+1)*sizeof(int));
	  for(k=0;k<=MAX_AM;k++)
	    memset(hrr_hash_table[i][j][k],0,(MAX_AM+1)*sizeof(int));
	}

      
      
      /* (a b|c d) */
      hrr_nodes[0].A = la-lb;
      hrr_nodes[0].B = lb;
      hrr_nodes[0].C = lc-ld;
      hrr_nodes[0].D = ld;
      hrr_nodes[0].llink = -1;
      hrr_nodes[0].rlink = -1;
      first_hrr_to_compute = 0;
      last_hrr_node = 0;
      mk_hrr_node(hrr_nodes[0], hrr_nodes, 0);
      hrr_nodes[0].target = 1;

      /*-------------------------------------------
	Traverse the graph starting at each target
       -------------------------------------------*/
      for(k=0; k<2; k++)
	if(hrr_nodes[0].children[k]>0)
	  mark_hrr_parents(hrr_nodes[0].children[k], hrr_nodes, 0);

      /*-----------------------------------------------------------------
	Empirical rule as with what heap size to start memory allocation
       -----------------------------------------------------------------*/
      init_mem(1);
      
      /*---------------------------------------------------------------
	Allocate and zero out space for classes to be generated by VRR
       ---------------------------------------------------------------*/
      for(i=last_hrr_node-1;i>=0;i--) {
	if (hrr_nodes[i].B == 0 && hrr_nodes[i].D == 0) {
	  hrr_nodes[i].marked = 1;
	  hrr_nodes[i].pointer = get_mem(hrr_nodes[i].size);
	  fprintf(hrr_code," Libint->vrr_classes[%d][%d] = int_stack + %d;\n",
		  hrr_nodes[i].A,hrr_nodes[i].C,hrr_nodes[i].pointer);
	}
      }
      base_mem = get_total_memory();
      fprintf(hrr_code," memset(int_stack,0,%d*sizeof(REALTYPE));\n\n",base_mem);
      fprintf(hrr_code," Libint->vrr_stack = int_stack + %d;\n",base_mem);
      
      /*----------------------------
	Build the HRR call sequence
       ----------------------------*/
      target_data = alloc_mem_hrr(hrr_nodes);
      hrr_mem = get_total_memory();
      if (max_stack_size < hrr_mem)
	max_stack_size = hrr_mem;
      fprintf(hrr_code," for(i=0;i<num_prim_comb;i++) {\n");
      fprintf(hrr_code,"   vrr_order_%c%c%c%c(Libint, Data);\n",
	      am_letter[la-lb],am_letter[lb],am_letter[lc-ld],am_letter[ld]);
      fprintf(hrr_code,"   Data++;\n }\n");
      
      j = first_hrr_to_compute;
      do {
	fprintf(hrr_code, " /*--- compute (%c%c|%c%c) ---*/\n",
		am_letter[hrr_nodes[j].A],am_letter[hrr_nodes[j].B],
		am_letter[hrr_nodes[j].C],am_letter[hrr_nodes[j].D]);
	if (hrr_nodes[j].B == 0 && hrr_nodes[j].D != 0) {
	  fprintf(hrr_code, " hrr3_build_%c%c(Libint->CD,int_stack+%d,int_stack+%d,",
		  am_letter[hrr_nodes[j].C], am_letter[hrr_nodes[j].D], hrr_nodes[j].pointer,
		  hrr_nodes[hrr_nodes[j].children[0]].pointer);
	  fprintf(hrr_code, "int_stack+%d,%d);\n", hrr_nodes[hrr_nodes[j].children[1]].pointer,
		  io[1+hrr_nodes[j].A]*io[1+hrr_nodes[j].B]);
	  /* Add this function to the list of inlined functions if necessary */
	  max_node_am = (hrr_nodes[j].C > hrr_nodes[j].D) ? hrr_nodes[j].C : hrr_nodes[j].D;
	  if (to_inline_hrr && max_node_am <= Params.max_am_to_inline_hrr_worker)
	    fprintf(inline_hrr_list,"#include \"hrr3_build_%c%c.h\"\n", am_letter[hrr_nodes[j].C], am_letter[hrr_nodes[j].D]);
	}
	else if (hrr_nodes[j].B != 0) {
	  fprintf(hrr_code, " hrr1_build_%c%c(Libint->AB,int_stack+%d,int_stack+%d,",
		  am_letter[hrr_nodes[j].A], am_letter[hrr_nodes[j].B], hrr_nodes[j].pointer,
		  hrr_nodes[hrr_nodes[j].children[0]].pointer);
	  fprintf(hrr_code, "int_stack+%d,%d);\n", hrr_nodes[hrr_nodes[j].children[1]].pointer,
		  io[1+hrr_nodes[j].C]*io[1+hrr_nodes[j].D]);
	  /* Add this function to the list of inlined functions if necessary */
	  max_node_am = (hrr_nodes[j].A > hrr_nodes[j].B) ? hrr_nodes[j].A : hrr_nodes[j].B;
	  if (to_inline_hrr && max_node_am <= Params.max_am_to_inline_hrr_worker)
	    fprintf(inline_hrr_list,"#include \"hrr1_build_%c%c.h\"\n", am_letter[hrr_nodes[j].A], am_letter[hrr_nodes[j].B]);
	}
	j = hrr_nodes[j].rlink;
      } while (j != -1);
      fprintf(hrr_code," return int_stack+%d;}\n",target_data);
      fclose(hrr_code);
      fclose(inline_hrr_list);
      printf("Done with %s\n",hrr_code_name);
      for(i=0;i<last_hrr_node;i++) {
	hrr_nodes[i].llink = 0;
        hrr_nodes[i].rlink = 0;
      }


      /*------------------------------------------------------------------
	Now generate the VRR graph using the (e0|f0) type classes present
	in the HRR graph as "potential" targets
       ------------------------------------------------------------------*/
      
      last_vrr_node = 0;
      num_vrr_targets = 0;
      for(i=0;i<last_hrr_node;i++)
	if (hrr_nodes[i].B == 0 && hrr_nodes[i].D == 0) {
	  if ((j = vrr_hash_table[hrr_nodes[i].A][hrr_nodes[i].C][0]) == 0) {
	    target_vrr_nodes[num_vrr_targets] = last_vrr_node;
	    vrr_nodes[last_vrr_node].A = hrr_nodes[i].A;
	    vrr_nodes[last_vrr_node].C = hrr_nodes[i].C;
	    vrr_nodes[last_vrr_node].m = 0;
	    vrr_nodes[last_vrr_node].llink = -1;
	    vrr_nodes[last_vrr_node].rlink = -1;
	    if (num_vrr_targets) {
	      vrr_nodes[target_vrr_nodes[num_vrr_targets-1]].llink = last_vrr_node;
	      vrr_nodes[last_vrr_node].rlink = target_vrr_nodes[num_vrr_targets-1];
	      vrr_nodes[last_vrr_node].llink = -1;
	    }
	    else {
	      vrr_nodes[last_vrr_node].rlink = -1;
	      vrr_nodes[last_vrr_node].llink = -1;
	    }
	    first_vrr_to_compute = last_vrr_node;
	    mk_vrr_node(vrr_nodes[last_vrr_node], vrr_nodes, 0);
	    vrr_nodes[first_vrr_to_compute].target = 1;
	    num_vrr_targets++;
	  }
	  else
	    vrr_nodes[j-1].target = 1;
	  if (first_vrr_to_compute == last_vrr_node && i == last_hrr_node-1)
	    punt("Edward, you fucked up\n");
	}

      /* Traverse the graph starting at each target */
      for(i=0;i<num_vrr_targets;i++) {
	j = target_vrr_nodes[i];
	for(k=0; k<5; k++){
	  if(vrr_nodes[j].children[k]>0){
	    mark_vrr_parents(vrr_nodes[j].children[k], vrr_nodes, j);
	  }
	}
      }

      /*-----------------------------------------------------------------
	Empirical rule as with what size heap to start memory allocation
       -----------------------------------------------------------------*/
      init_mem(1);


      /* Build the call sequence */
      target_data = alloc_mem_vrr(vrr_nodes);
      vrr_mem = base_mem + get_total_memory();
      if (max_stack_size < vrr_mem)
	max_stack_size = vrr_mem;
      fprintf(vrr_code," REALTYPE *vrr_stack = Libint->vrr_stack;\n");
      fprintf(vrr_code," REALTYPE *tmp, *target_ptr;\n int i, am[2];\n");
      
      j = first_vrr_to_compute;
      do {
	if (vrr_nodes[j].A <= opt_am && vrr_nodes[j].C <= opt_am) {
	  fprintf(vrr_code, " _BUILD_%c0%c0(Data,", am_letter[vrr_nodes[j].A], am_letter[vrr_nodes[j].C]);
	  /* Add this function to the list of inlined functions if necessary */
	  max_node_am = (vrr_nodes[j].A > vrr_nodes[j].C) ? vrr_nodes[j].A : vrr_nodes[j].C;
	  if (to_inline_vrr && max_node_am <= Params.max_am_to_inline_vrr_worker)
	    fprintf(inline_vrr_list,"#include \"build_%c0%c0.h\"\n", am_letter[vrr_nodes[j].A], am_letter[vrr_nodes[j].C]);
	}
	else {
	  fprintf(vrr_code, " am[0] = %d;  am[1] = %d;\n", vrr_nodes[j].A, vrr_nodes[j].C);
	  fprintf(vrr_code, " vrr_build_xxxx(am,Data,");
	}
	fprintf(vrr_code, "vrr_stack+%d", vrr_nodes[j].pointer);
	for(k=0; k<5; k++){
	  if(vrr_nodes[j].children[k] > 0)
	    fprintf(vrr_code, ", vrr_stack+%d", vrr_nodes[vrr_nodes[j].children[k]].pointer);
	  else if (vrr_nodes[j].children[k] == NONODE)
	    fprintf(vrr_code, ", NULL");
	  else
	    fprintf(vrr_code, ", Data->F+%d", (-1)*vrr_nodes[j].children[k]);
	}
	fprintf(vrr_code, ");\n");
	if (vrr_nodes[j].target == 1) {
	  fprintf(vrr_code, "   tmp = vrr_stack + %d;\n", vrr_nodes[j].pointer);
	  fprintf(vrr_code, "   target_ptr = Libint->vrr_classes[%d][%d];\n",vrr_nodes[j].A,vrr_nodes[j].C);
	  fprintf(vrr_code, "   for(i=0;i<%d;i++)\n",(vrr_nodes[j].A+1)*(vrr_nodes[j].A+2)*(vrr_nodes[j].C+1)*(vrr_nodes[j].C+2)/4);
	  fprintf(vrr_code, "     target_ptr[i] += tmp[i];\n");
	}
	j = vrr_nodes[j].rlink;
      } while (j != -1);
      fprintf(vrr_code, "\n}\n\n");
      fclose(vrr_code);
      fclose(inline_vrr_list);
      printf("Done with %s\n",vrr_code_name);
      for(i=0;i<last_vrr_node;i++) {
	vrr_nodes[i].llink = 0;
        vrr_nodes[i].rlink = 0;
      }
      
      
      /* compare this max_stack_size to the libint_stack_size for
	 this angular momentum  and reset it */
      if (libint_stack_size[current_highest_am] < max_stack_size)
	libint_stack_size[current_highest_am] = max_stack_size;

      max_stack_size = 0;

      }
    }
    }
  }
  return;
}


/* Recursive function that build the HRR subgraph given the parent */

int mk_hrr_node(hrr_class node, hrr_class *allnodes, int new)
{

  int i, j, k, l;
  hrr_class O[2];
  int subnodes = 0;
  int thisnode;
  int rlink, llink;
  int made = 0;

  /* Search for the parent node on stack
     If it's not there - we'll add it to the end of the stack */
  thisnode = last_hrr_node;
  /* it's already placed on the stack allnodes - make sure children don't get created again (made = 1) */
  if (hrr_hash_table[node.A][node.B][node.C][node.D]) {
    i = hrr_hash_table[node.A][node.B][node.C][node.D] - 1;
    thisnode = i;
    made = 1;
  }

  /* it's not computed, add it, and make it the first to compute! */
  if(!made){
    allnodes[thisnode].A = node.A;
    allnodes[thisnode].B = node.B;
    allnodes[thisnode].C = node.C;
    allnodes[thisnode].D = node.D;
    hrr_hash_table[node.A][node.B][node.C][node.D] = thisnode + 1;
    allnodes[thisnode].num_parents = 0;
    allnodes[thisnode].parents_counter = 0;
    allnodes[thisnode].marked = 0;
    allnodes[thisnode].pointer = 0;
    memset(allnodes[thisnode].parents,0,5*sizeof(int));
    allnodes[thisnode].children[0] = NONODE;
    allnodes[thisnode].children[1] = NONODE;
    allnodes[thisnode].size = io[1+node.A]*io[1+node.B]*io[1+node.C]*io[1+node.D];
    allnodes[thisnode].target = 0;
    /* We just added a node ..*/
    last_hrr_node++;
    /* If stack is overfull - exit */
    if(last_hrr_node==MAXNODE) {
      printf(" Maximum stack size is reached. Change MAXNODE and recompile.\n\n");
      exit(1);
    }
  }

  /* If the parent class wasn't on stack already (!new) - increase the parent counter */
  if(!new){
    allnodes[thisnode].num_parents++;
    allnodes[thisnode].parents_counter++;
    }


  /* now make all child nodes */
  if (!made) {
  if(node.B){
    O[0].A = node.A+1;
    O[0].B = node.B-1;
    O[0].C = node.C;
    O[0].D = node.D;
    allnodes[thisnode].children[0] = 
            mk_hrr_node(O[0], allnodes, made);
    O[1].A = node.A;
    O[1].B = node.B-1;
    O[1].C = node.C;
    O[1].D = node.D;
    allnodes[thisnode].children[1] =
	    mk_hrr_node(O[1], allnodes, made);
  }
  else if(node.D){
    O[0].A = node.A;
    O[0].B = node.B;
    O[0].C = node.C+1;
    O[0].D = node.D-1;
    allnodes[thisnode].children[0] = 
            mk_hrr_node(O[0], allnodes, made);
    O[1].A = node.A;
    O[1].B = node.B;
    O[1].C = node.C;
    O[1].D = node.D-1;
    allnodes[thisnode].children[1] = 
            mk_hrr_node(O[1], allnodes, made);
  }
  }

  return thisnode;

}


/* Recursive function that build the VRR subgraph given the parent */

int mk_vrr_node(vrr_class node, vrr_class *allnodes, int new)
{

  int i, j, k, l;
  vrr_class O[5];
  int subnodes = 0;
  int thisnode;
  int rlink, llink;
  int made = 0;

  /* Are there any children? */
  if(node.A+node.C == 0) return (-1)*node.m;

  /* Search for the parent node on stack
     If it's not there - we'll add it to the end of the stack */
  thisnode = last_vrr_node;
  /* it's already placed on the stack allnodes - make sure children don't get created again (made = 1) */
  if (vrr_hash_table[node.A][node.C][node.m]) {
    i = vrr_hash_table[node.A][node.C][node.m] - 1;
    thisnode = i;
    made = 1;
  }

  /* it's not computed, add it, and make it the first to compute! */
  if(!made){
    allnodes[thisnode].A = node.A;
    allnodes[thisnode].C = node.C;
    allnodes[thisnode].m = node.m;
    vrr_hash_table[node.A][node.C][node.m] = thisnode + 1;
    allnodes[thisnode].num_parents = 0;
    allnodes[thisnode].parents_counter = 0;
    allnodes[thisnode].marked = 0;
    allnodes[thisnode].target = 0;
    allnodes[thisnode].pointer = 0;
    memset(allnodes[thisnode].parents,0,10*sizeof(int));
    allnodes[thisnode].children[0] = NONODE;
    allnodes[thisnode].children[1] = NONODE;
    allnodes[thisnode].children[2] = NONODE;
    allnodes[thisnode].children[3] = NONODE;
    allnodes[thisnode].children[4] = NONODE;
    allnodes[thisnode].size = io[1+node.A]*io[1+node.C];
    /* We just added a node ..*/
    last_vrr_node++;
    /* If stack is overfull - exit */
    if(last_vrr_node==MAXNODE) {
      printf(" Maximum stack size is reached. Change MAXNODE and recompile.\n\n");
      exit(1);
    }
  }

  /* If the parent class wasn't on stack already (!new) - increase the parent counter */
  if(!new){
    allnodes[thisnode].num_parents++;
    allnodes[thisnode].parents_counter++;
    }


  /* now make all child nodes */
  if (!made) {
  if(node.A){
    O[0].A = node.A-1;
    O[0].C = node.C;
    O[0].m = node.m;
    allnodes[thisnode].children[0] = 
            mk_vrr_node(O[0], allnodes, made);
    O[1].A = node.A-1;
    O[1].C = node.C;
    O[1].m = node.m+1;
    allnodes[thisnode].children[1] = 
            mk_vrr_node(O[1], allnodes, made);
    if(node.A>1){
      O[2].A = node.A-2;
      O[2].C = node.C;
      O[2].m = node.m;
      allnodes[thisnode].children[2] = 
            mk_vrr_node(O[2], allnodes, made);
      O[3].A = node.A-2;
      O[3].C = node.C;
      O[3].m = node.m+1;
      allnodes[thisnode].children[3] = 
            mk_vrr_node(O[3], allnodes, made);
      }
    if(node.C){
      O[4].A = node.A-1;
      O[4].C = node.C-1;
      O[4].m = node.m+1;
      allnodes[thisnode].children[4] = 
            mk_vrr_node(O[4], allnodes, made);
      }
    }
  else if(node.C){
    O[0].A = node.A;
    O[0].C = node.C-1;
    O[0].m = node.m;
    allnodes[thisnode].children[0] = 
            mk_vrr_node(O[0], allnodes, made);
    O[1].A = node.A;
    O[1].C = node.C-1;
    O[1].m = node.m+1;
    allnodes[thisnode].children[1] = 
            mk_vrr_node(O[1], allnodes, made);
    if(node.C>1){
      O[2].A = node.A;
      O[2].C = node.C-2;
      O[2].m = node.m;
      allnodes[thisnode].children[2] = 
            mk_vrr_node(O[2], allnodes, made);
      O[3].A = node.A;
      O[3].C = node.C-2;
      O[3].m = node.m+1;
      allnodes[thisnode].children[3] = 
            mk_vrr_node(O[3], allnodes, made);
      }
    }
  }

  return thisnode;

}


/* Make hrr_nodes[rent] a parent of hrr_nodes[n] and proceed recursively */

int mark_hrr_parents(int n, hrr_class *allnodes, int rent)
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
    allnodes[n].llink = -1;
    allnodes[n].rlink = first_hrr_to_compute;
    allnodes[first_hrr_to_compute].llink = n;
    first_hrr_to_compute = n;

    for(i=0; i<2; i++)
      if(allnodes[n].children[i]>0)
	mark_hrr_parents(allnodes[n].children[i], allnodes, n);
  }
  return;
}


/* Make vrr_nodes[rent] a parent of vrr_nodes[n] and proceed recursively */

int mark_vrr_parents(int n, vrr_class *allnodes, int rent)
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
    allnodes[n].llink = -1;
    allnodes[n].rlink = first_vrr_to_compute;
    allnodes[first_vrr_to_compute].llink = n;
    first_vrr_to_compute = n;

    for(i=0; i<5; i++)
      if(allnodes[n].children[i]>0)
	mark_vrr_parents(allnodes[n].children[i], allnodes, n);

  }
  return;
}



/* This functions controls memory placement of computed classes on the CINTS stack */

int alloc_mem_hrr(hrr_class *nodes)
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
    for(k=0; k<2; k++){
      child = nodes[j].children[k];
      if(child>0)
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

int alloc_mem_vrr(vrr_class *nodes)
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
    for(k=0; k<5; k++){
      child = nodes[j].children[k];
      if(child>0){
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
