/*! \file
    \ingroup DERIV
    \brief Enter brief description of file here 
*/
/*------------------------------------------------------------------------------------------------------
  Builds a library of functions-calls for applying HRR and VRR
 ------------------------------------------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <libint/libint.h>
#include "mem_man.h"
#include "build_libderiv.h"
#define MAXNODE 10000
#define MAX_NUM_TARGET_VRR_NODES 2000
#define NONODE -1000000
#define NUMPARENTS 40

#define INDEX(a,b) (((a) > (b)) ? io[(a)] + (b) : io[(b)] + (a))

static int last_hrr_node = 0;      /* Global pointer to the last node on the HRR stack */
static int last_vrr_node = 0;      /* Global pointer to the last node on the VRR stack */

extern FILE *outfile, *d1hrr_header, *init_code;
extern int libderiv12_stack_size[MAX_AM/2+1];
extern LibderivParams_t Params;

typedef struct node{
  int A, B, C, D;         /* Angular momenta on centers A and C */
  int m;
  int deriv_lvl;          /* Derivative level */
  int deriv_ind[12];      /* Derivative indices along each of nuclear coordinates */
  int size;               /* Class size in double words */
  int pointer;
  int children[8];        /* Up to 8 children of the class */
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
static int hrr_hash_table[2*LMAX_AM][2*LMAX_AM][2*LMAX_AM][2*LMAX_AM];
static int vrr_hash_table[2*LMAX_AM][2*LMAX_AM][4*LMAX_AM];

static int mk_dhrr_node(class node, class *allnodes, int new);
static int mk_deriv_node(class node, class *allnodes, int new);
static void mark_dhrr_parents(int n, class *allnodes, int rent);
static void mark_vrr_parents(int n, class *allnodes, int rent);
static void mark_parents(int n, class *allnodes, int rent);
static int alloc_mem_dhrr(class *nodes);
static int alloc_mem_vrr(class *nodes);

static void get_deriv_indices(class *node,int * di, int *dj);

void emit_deriv12_managers()
{
  int new_am = Params.new_am12;
  int old_am = Params.old_am;
  int opt_am = Params.opt_am;
  int am_to_inline_into_hrr = Params.max_am_to_inline_vrr_manager;
  int am_to_inline_vrr = Params.max_am_manager_to_inline_vrr_worker;
  int am_to_inline_hrr = Params.max_am_manager_to_inline_hrr_worker;
  int am_to_inline_deriv = Params.max_am_manager_to_inline_deriv_worker;
  int am_to_inline_d1hrr = Params.max_am_manager_to_inline_d1hrr_worker;
  int to_inline_into_hrr, to_inline_vrr, to_inline_hrr, to_inline_deriv, to_inline_d1hrr;

  int i, j, k, l;
  int di, dj, dk, dl;
  int la, lc, lc_min, ld, ld_max, ld_min;
  int lb, lb_min, lb_max;
  int current_highest_am;
  int max_node_am;
  int last_mem;
  int child0, child1, child;
  int num_children;
  int offset;

  class nodes[MAXNODE];    /* Stack of nodes */
  class *hrr_nodes = &(nodes[0]);
  class *vrr_nodes;
  int target_data;
  int done;
  int max_stack_size = 0;

  int target_hrr_nodes[144+12];         /* List of unique targets on the top of HRR tree (144 second derivatives + 12 first derivatives + 1 ERI) */
  int target_d2hrr_nodes[12][12];         /* This is used to refer to second derivative targets only */ 
  int target_d1hrr_nodes[12];             /* This is used to refer to first derivative targets only */ 
  int num_hrr_targets;
  int first_empty_slot;

  int target_vrr_nodes[MAX_NUM_TARGET_VRR_NODES];
  int num_vrr_targets;
  const char am_letter[] = "0pdfghiklmnoqrtuvwxyz";
  const char cart_comp[] = "XYZ";
  char hrr_code_name[80];
  char hrr_function_name[80];
  char vrr_code_name[80];
  char vrr_function_name[80];
  char inline_vrr_list_name[80];
  char inline_hrr_list_name[80];
  FILE *hrr_code, *vrr_code, *inline_vrr_list, *inline_hrr_list;
  static int io[] = {1,3,6,10,15,21,28,36,45,55,66,78,91,105,120,136,153,171,190,210};


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
      to_inline_into_hrr = (current_highest_am <= am_to_inline_into_hrr) ? 1 : 0;
      to_inline_vrr = (current_highest_am <= am_to_inline_vrr) ? 1 : 0;
      to_inline_hrr = (current_highest_am <= am_to_inline_hrr) ? 1 : 0;
      to_inline_deriv = (current_highest_am <= am_to_inline_deriv) ? 1 : 0;
      to_inline_d1hrr = (current_highest_am <= am_to_inline_d1hrr) ? 1 : 0;

      /*---------------------------------------------------------------
	Form code and function names for HRR and VRR ordering routines
       ---------------------------------------------------------------*/
      sprintf(hrr_function_name,"d12hrr_order_%c%c%c%c",
	      am_letter[la-lb],
	      am_letter[lb],
	      am_letter[lc-ld],
	      am_letter[ld]);
      sprintf(vrr_function_name,"d12vrr_order_%c%c%c%c",
	      am_letter[la-lb],
	      am_letter[lb],
	      am_letter[lc-ld],
	      am_letter[ld]);
      sprintf(hrr_code_name,"%s.cc",hrr_function_name);
      if (to_inline_into_hrr)
	sprintf(vrr_code_name,"%s.h",vrr_function_name);
      else
	sprintf(vrr_code_name,"%s.cc",vrr_function_name);
      sprintf(inline_vrr_list_name,"inline_%s.h",vrr_function_name);
      sprintf(inline_hrr_list_name,"inline_%s.h",hrr_function_name);
      hrr_code = fopen(hrr_code_name,"w");
      vrr_code = fopen(vrr_code_name,"w");
      inline_vrr_list = fopen(inline_vrr_list_name,"w");
      inline_hrr_list = fopen(inline_hrr_list_name,"w");

      /*-----------------------------------
	Write the overhead to the HRR code
       -----------------------------------*/
      fprintf(hrr_code,"#include <stdio.h>\n");
      fprintf(hrr_code,"#include <string.h>\n");
      fprintf(hrr_code,"#include <libint/libint.h>\n");
      fprintf(hrr_code,"#include \"libderiv.h\"\n");
      if (to_inline_hrr)
	fprintf(hrr_code,"#define INLINE_HRR_WORKER\n");
      if (to_inline_d1hrr)
	fprintf(hrr_code,"#define INLINE_D1HRR_WORKER\n");
      if (to_inline_hrr || to_inline_d1hrr)
	fprintf(hrr_code,"#include \"%s\"\n",inline_hrr_list_name);
      fprintf(hrr_code,"#include <libint/hrr_header.h>\n\n");
      fprintf(hrr_code,"#include \"d1hrr_header.h\"\n\n");
      if (to_inline_into_hrr)
	fprintf(hrr_code,"#include \"%s\"\n",vrr_code_name);
      else
	fprintf(hrr_code,"extern void %s(Libderiv_t *, prim_data *);\n\n",vrr_function_name);
      fprintf(hrr_code,"  /* Computes derivatives of (%c%c|%c%c) integrals */\n\n",
	      am_letter[la-lb],am_letter[lb],am_letter[lc-ld],am_letter[ld]);
      fprintf(hrr_code,"void %s(Libderiv_t *Libderiv, int num_prim_comb)\n{\n",hrr_function_name);
      fprintf(hrr_code," prim_data *Data = Libderiv->PrimQuartet;\n");
      fprintf(hrr_code," double *int_stack = Libderiv->int_stack;\n");
      fprintf(hrr_code," double *zero_stack = Libderiv->zero_stack;\n");
      fprintf(hrr_code," int i,j;\n double tmp, *target;\n\n");

      /*-----------------------------------------------------------------
	Include the function into the dhrr_header.h and init_libderiv.cc
       -----------------------------------------------------------------*/
      fprintf(d1hrr_header,"void %s(Libderiv_t *, int);\n",hrr_function_name);
      fprintf(init_code,"  build_deriv12_eri[%d][%d][%d][%d] = %s;\n",la-lb,lb,lc-ld,ld,hrr_function_name);

      /*-----------------------------------
	Write the overhead to the VRR code
       -----------------------------------*/
      fprintf(vrr_code,"#include <stdio.h>\n");
      fprintf(vrr_code,"#include <libint/libint.h>\n");
      fprintf(vrr_code,"#include \"libderiv.h\"\n");
      if (to_inline_vrr)
	fprintf(vrr_code,"#define INLINE_VRR_WORKER\n");
      if (to_inline_deriv)
	fprintf(vrr_code,"#define INLINE_DERIV_WORKER\n");
      if (to_inline_hrr)
	fprintf(vrr_code,"#define INLINE_HRR_WORKER\n");
      if (to_inline_vrr || to_inline_deriv || to_inline_hrr)
	fprintf(vrr_code,"#include \"%s\"\n",inline_vrr_list_name);
      fprintf(vrr_code,"#include <libint/vrr_header.h>\n");
      fprintf(vrr_code,"#include <libint/hrr_header.h>\n");
      fprintf(vrr_code,"#include \"deriv_header.h\"\n\n");
      fprintf(vrr_code,"  /* Computes quartets necessary to compute derivatives of (%c%c|%c%c) integrals */\n\n",
	      am_letter[la-lb],am_letter[lb],am_letter[lc-ld],am_letter[ld]);
      if (to_inline_into_hrr)
	fprintf(vrr_code,"inline ");
      fprintf(vrr_code,"void %s(Libderiv_t *Libderiv, prim_data *Data)\n{\n",vrr_function_name);

      
      /*--------------------------------------------------
	Starting at the target node(s) set up an HRR graph.
       --------------------------------------------------*/
      last_hrr_node = 0;
      num_hrr_targets=0;

      /*--- First add second derivative ERIs ---*/
      for(di=0;di<12;di++) { /* Now USING translational invariance here */
	if (di<3 || di>5)
	  for(dj=di;dj<12;dj++)
	    if (dj<3 || dj>5)
	      {
	      target_hrr_nodes[num_hrr_targets] = last_hrr_node;
	      hrr_nodes[last_hrr_node].A = la-lb;
	      hrr_nodes[last_hrr_node].B = lb;
	      hrr_nodes[last_hrr_node].C = lc-ld;
	      hrr_nodes[last_hrr_node].D = ld;
	      hrr_nodes[last_hrr_node].m = 0;
	      hrr_nodes[last_hrr_node].deriv_lvl = 2;
	      for(dk=0;dk<12;dk++)
		hrr_nodes[last_hrr_node].deriv_ind[dk] = 0;
	      hrr_nodes[last_hrr_node].deriv_ind[di] += 1;
	      hrr_nodes[last_hrr_node].deriv_ind[dj] += 1;
	      first_hrr_to_compute = last_hrr_node;
	      k = mk_dhrr_node(hrr_nodes[last_hrr_node], hrr_nodes, 1);
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
	      hrr_nodes[k].target = 1;
	    }
      }

#if 1
      /*--- Now add first derivative ERIs ---*/
      for(di=0;di<12;di++) 
        if (di<3 || di>5) {
	target_hrr_nodes[num_hrr_targets] = last_hrr_node;
	hrr_nodes[last_hrr_node].A = la-lb;
	hrr_nodes[last_hrr_node].B = lb;
	hrr_nodes[last_hrr_node].C = lc-ld;
	hrr_nodes[last_hrr_node].D = ld;
	hrr_nodes[last_hrr_node].m = 0;
	hrr_nodes[last_hrr_node].deriv_lvl = 1;
	for(dk=0;dk<12;dk++)
	  hrr_nodes[last_hrr_node].deriv_ind[dk] = 0;
	hrr_nodes[last_hrr_node].deriv_ind[di] += 1;
	first_hrr_to_compute = last_hrr_node;
	k = mk_dhrr_node(hrr_nodes[last_hrr_node], hrr_nodes, 1);
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
	hrr_nodes[k].target = 1;
	}
#endif

	/*-------------------------------------------
	  Traverse the graph starting at each target
	 -------------------------------------------*/
      for(i=0;i<num_hrr_targets;i++) {
	j = target_hrr_nodes[i];
	for(k=0; k<8; k++)
	  if(hrr_nodes[j].children[k]>0){
	    mark_dhrr_parents(hrr_nodes[j].children[k], hrr_nodes, j);
	  }
      }

      /* Empirically determined that it's better to start with a small stack */
      init_mem(1);

      
      /*---------------------------------------------------------------
	Allocate and zero out space for classes to be generated by VRR
       ---------------------------------------------------------------*/
      for(i=last_hrr_node-1;i>=0;i--) {
	if (hrr_nodes[i].B == 0 && hrr_nodes[i].D == 0) {
	  hrr_nodes[i].marked = 1;
	  hrr_nodes[i].pointer = get_mem(hrr_nodes[i].size);
	  get_deriv_indices(&hrr_nodes[i],&di,&dj);
	  if (hrr_nodes[i].deriv_lvl == 0)
	    fprintf(hrr_code," Libderiv->dvrr_classes[%d][%d] = int_stack + %d;\n",
		    hrr_nodes[i].A, hrr_nodes[i].C, hrr_nodes[i].pointer);
	  else if (hrr_nodes[i].deriv_lvl == 1) {
	    fprintf(hrr_code," Libderiv->deriv_classes[%d][%d][%d] = int_stack + %d;\n",
		    hrr_nodes[i].A, hrr_nodes[i].C, di, hrr_nodes[i].pointer);
	  }
	  else if (hrr_nodes[i].deriv_lvl == 2) {
	    fprintf(hrr_code," Libderiv->deriv2_classes[%d][%d][%d] = int_stack + %d;\n",
		    hrr_nodes[i].A, hrr_nodes[i].C, di*12+dj, hrr_nodes[i].pointer);
	  }
	  else {
	    punt("deriv_lvl = 3\n");
	  }
	}
      }
      fprintf(hrr_code," memset(int_stack,0,%lu);\n\n",get_total_memory()*sizeof(double));

      
      /*----------------------------
	Build the HRR call sequence
       ----------------------------*/
      if (lb != 0 || ld != 0) {
	target_data = alloc_mem_dhrr(hrr_nodes);
      }


      last_mem = get_total_memory();
      fprintf(hrr_code," Libderiv->dvrr_stack = int_stack + %d;\n",last_mem);
      fprintf(hrr_code," for(i=0;i<num_prim_comb;i++) {\n");
      fprintf(hrr_code,"   %s(Libderiv, Data);\n",vrr_function_name);
      fprintf(hrr_code,"   Data++;\n }\n\n");

      
      /*------------------------------------------------
	Evaluate the HRR tree for each derivative class
       ------------------------------------------------*/
      /*--- If we have non-(ss|ss) class - perform the standard procedure ---*/
      j = first_hrr_to_compute;
      do {
	fprintf(hrr_code, " /*--- compute (%c%c|%c%c) ---*/\n",
		am_letter[hrr_nodes[j].A],am_letter[hrr_nodes[j].B],
		am_letter[hrr_nodes[j].C],am_letter[hrr_nodes[j].D]);

	if (hrr_nodes[j].B > 0 || hrr_nodes[j].D > 0) {
	  /*--- compute the number of children ---*/
	  num_children = 0;
	  for(i=0;i<8;i++)
	    if (hrr_nodes[j].children[i] > 0)
	      num_children++;
	    
	  if (hrr_nodes[j].B == 0 && hrr_nodes[j].D != 0) {
	    offset = 6;
	    if (num_children > 2) {
	      fprintf(hrr_code, "   d1hrr3_build_%c%c(Libderiv->CD,int_stack+%d,",
		      am_letter[hrr_nodes[j].C], am_letter[hrr_nodes[j].D], hrr_nodes[j].pointer);
	      /* Add this function to the list of inlined functions if necessary */
	      max_node_am = (hrr_nodes[j].C > hrr_nodes[j].D) ? hrr_nodes[j].C : hrr_nodes[j].D;
	      if (to_inline_d1hrr && max_node_am <= Params.max_am_to_inline_d1hrr_worker)
		fprintf(inline_hrr_list,"#include \"d1hrr3_build_%c%c.h\"\n", am_letter[hrr_nodes[j].C], am_letter[hrr_nodes[j].D]);
	    }
	    else {
	      fprintf(hrr_code, "   hrr3_build_%c%c(Libderiv->CD,int_stack+%d,",
		      am_letter[hrr_nodes[j].C], am_letter[hrr_nodes[j].D], hrr_nodes[j].pointer);
	      /* Add this function to the list of inlined functions if necessary */
	      max_node_am = (hrr_nodes[j].C > hrr_nodes[j].D) ? hrr_nodes[j].C : hrr_nodes[j].D;
	      if (to_inline_hrr && max_node_am <= Params.max_am_to_inline_hrr_worker)
		fprintf(inline_hrr_list,"#include <libint/hrr3_build_%c%c.h>\n",
			am_letter[hrr_nodes[j].C], am_letter[hrr_nodes[j].D]);
	    }
	  }
	  else if (hrr_nodes[j].B != 0) {
	    offset = 0;
	    if (num_children > 2) {
	      fprintf(hrr_code, "   d1hrr1_build_%c%c(Libderiv->AB,int_stack+%d,",
		      am_letter[hrr_nodes[j].A], am_letter[hrr_nodes[j].B], hrr_nodes[j].pointer);
	      /* Add this function to the list of inlined functions if necessary */
	      max_node_am = (hrr_nodes[j].A > hrr_nodes[j].B) ? hrr_nodes[j].A : hrr_nodes[j].B;
	      if (to_inline_d1hrr && max_node_am <= Params.max_am_to_inline_d1hrr_worker)
		fprintf(inline_hrr_list,"#include \"d1hrr1_build_%c%c.h\"\n", am_letter[hrr_nodes[j].A], am_letter[hrr_nodes[j].B]);
	    }
	    else {
	      fprintf(hrr_code, "   hrr1_build_%c%c(Libderiv->AB,int_stack+%d,",
		      am_letter[hrr_nodes[j].A], am_letter[hrr_nodes[j].B], hrr_nodes[j].pointer);
	      /* Add this function to the list of inlined functions if necessary */
	      max_node_am = (hrr_nodes[j].A > hrr_nodes[j].B) ? hrr_nodes[j].A : hrr_nodes[j].B;
	      if (to_inline_hrr && max_node_am <= Params.max_am_to_inline_hrr_worker)
		fprintf(inline_hrr_list,"#include <libint/hrr1_build_%c%c.h>\n", am_letter[hrr_nodes[j].A], am_letter[hrr_nodes[j].B]);
	    }
	  }
	
	  /*--- If the first child is one of VRR derivative classes - need to compute its location ---*/
	  child0 = hrr_nodes[j].children[0];
	  fprintf(hrr_code, "int_stack+%d,",hrr_nodes[child0].pointer);

	  /*--- If the second child is one of VRR derivative classes - need to compute its location ---*/
	  child1 = hrr_nodes[j].children[1];
	  fprintf(hrr_code, "int_stack+%d,",hrr_nodes[child1].pointer);

	  /*--- Now go through the rest of the children ---*/
	  if (num_children > 2)
	    for(i=0;i<6;i++) {
	      if (hrr_nodes[j].children[i+2] > 0) {
		child = hrr_nodes[j].children[i+2];
		fprintf(hrr_code, " %.1lf, int_stack+%d, ",
			(double) hrr_nodes[j].deriv_ind[offset+i],
			hrr_nodes[child].pointer);
	      }
	      else
		fprintf(hrr_code, " 0.0, zero_stack,");
	    }
	  
	    
	  if (hrr_nodes[j].B == 0 && hrr_nodes[j].D != 0)
	    fprintf(hrr_code, "%d);\n", io[hrr_nodes[j].A]*io[hrr_nodes[j].B]);
	  else if (hrr_nodes[j].B != 0)
	    fprintf(hrr_code, "%d);\n", io[hrr_nodes[j].C]*io[hrr_nodes[j].D]);

	}
	
	/* Pass the "target" quartets to CINTS */
        if (hrr_nodes[j].target) {
	  get_deriv_indices(&hrr_nodes[j],&di,&dj);
	  if (hrr_nodes[j].deriv_lvl == 1)
	    fprintf(hrr_code,"     Libderiv->ABCD[%d] = int_stack + %d;\n",di,hrr_nodes[j].pointer);
	  else if (hrr_nodes[j].deriv_lvl == 2)
	    fprintf(hrr_code,"     Libderiv->ABCD[%d] = int_stack + %d;\n",12+di*12+dj
		    /*INDEX(di,dj)*/
,hrr_nodes[j].pointer);
	  else
	    punt("deriv_lvl = 3\n");
	}
	j = hrr_nodes[j].rlink;
      } while (j != -1);
      
      fprintf(hrr_code,"\n}\n",target_data);
      fclose(hrr_code);
      fclose(inline_hrr_list);
      printf("Done with %s\n",hrr_code_name);
      for(i=0;i<last_hrr_node;i++) {
	hrr_nodes[i].llink = 0;
        hrr_nodes[i].rlink = 0;
      }

      /*----------------------------
	Zero out the hashing tables
       ----------------------------*/
      for(i=0;i<2*LMAX_AM;i++)
	for(j=0;j<2*LMAX_AM;j++)
	  memset(vrr_hash_table[i][j],0,(4*LMAX_AM)*sizeof(int));
      for(i=0;i<2*LMAX_AM;i++)
	for(j=0;j<2*LMAX_AM;j++)
	  for(k=0;k<2*LMAX_AM;k++)
	    memset(hrr_hash_table[i][j][k],0,(2*LMAX_AM)*sizeof(int));

      
      /*------------------------------------------------------------------
	Now generate the VRR graph using the (e0|f0) type classes present
	in the HRR graph as "potential" targets
       ------------------------------------------------------------------*/

      vrr_nodes = &(hrr_nodes[last_hrr_node]);
      last_vrr_node = 0;
      num_vrr_targets = 0;
      for(i=0;i<last_hrr_node;i++)
	if (hrr_nodes[i].B == 0 && hrr_nodes[i].D == 0) {
	  target_vrr_nodes[num_vrr_targets] = last_vrr_node;
	  vrr_nodes[last_vrr_node].A = hrr_nodes[i].A;
	  vrr_nodes[last_vrr_node].B = hrr_nodes[i].B;
	  vrr_nodes[last_vrr_node].C = hrr_nodes[i].C;
	  vrr_nodes[last_vrr_node].D = hrr_nodes[i].D;
	  vrr_nodes[last_vrr_node].m = hrr_nodes[i].m;
	  vrr_nodes[last_vrr_node].deriv_lvl = hrr_nodes[i].deriv_lvl;
	  memcpy(vrr_nodes[last_vrr_node].deriv_ind,hrr_nodes[i].deriv_ind,12*sizeof(int));
	  first_empty_slot = last_vrr_node;
	  k = mk_deriv_node(vrr_nodes[last_vrr_node], vrr_nodes, 1);
	  if (k == first_empty_slot) { /* If the node hasn't been added to the tree before */
	    if (num_vrr_targets) {
	      vrr_nodes[target_vrr_nodes[num_vrr_targets-1]].llink = k;
	      vrr_nodes[k].rlink = target_vrr_nodes[num_vrr_targets-1];
	      vrr_nodes[k].llink = -1;
	    }
	    else {
	      vrr_nodes[k].rlink = -1;
	      vrr_nodes[k].llink = -1;
	    }
	    num_vrr_targets++;
	    first_vrr_to_compute = k;
	  }
	  vrr_nodes[k].target = 1;
	  if (first_vrr_to_compute == last_vrr_node && i == last_hrr_node-1)
	    punt("Edward, you fucked up\n");
	}

      /* Traverse the graph starting at each target */
      if (num_vrr_targets >= MAX_NUM_TARGET_VRR_NODES)
	punt("Number of VRR targets has exceeded the hardwired maximum MAX_NUM_TARGET_VRR_NODES\n");
      for(i=0;i<num_vrr_targets;i++) {
	j = target_vrr_nodes[i];
	for(k=0; k<5; k++){
	  if(vrr_nodes[j].children[k]>0){
	    mark_parents(vrr_nodes[j].children[k], vrr_nodes, j);
	  }
	}
      }

      /* Empirically determined that it's better to start with a small stack */
      init_mem(1);

      /* Build the call sequence */
      target_data = alloc_mem_vrr(vrr_nodes);
      last_mem += get_total_memory();
      if (max_stack_size < last_mem)
	max_stack_size = last_mem;
      fprintf(vrr_code," double *dvrr_stack = Libderiv->dvrr_stack;\n double *tmp, *target_ptr;\n");
      fprintf(vrr_code," int i, am[2];\n");
      
      j = first_vrr_to_compute;
      do {
	fprintf(vrr_code, " /* compute (%d %d | %d %d) m=%d deriv level %d */\n",
		vrr_nodes[j].A,vrr_nodes[j].B,
		vrr_nodes[j].C,vrr_nodes[j].D,
		vrr_nodes[j].m, vrr_nodes[j].deriv_lvl);
	fprintf(vrr_code, " /* deriv_ind: %d %d %d  %d %d %d  %d %d %d  %d %d %d */\n",
		vrr_nodes[j].deriv_ind[0], vrr_nodes[j].deriv_ind[1], vrr_nodes[j].deriv_ind[2],
		vrr_nodes[j].deriv_ind[3], vrr_nodes[j].deriv_ind[4], vrr_nodes[j].deriv_ind[5],
		vrr_nodes[j].deriv_ind[6], vrr_nodes[j].deriv_ind[7], vrr_nodes[j].deriv_ind[8],
		vrr_nodes[j].deriv_ind[9], vrr_nodes[j].deriv_ind[10], vrr_nodes[j].deriv_ind[11]);

	/*---------------------------------------------------------
	  Decide which routine to use to compute the current class
	 ---------------------------------------------------------*/
	if (vrr_nodes[j].deriv_lvl) { /*--- use build_deriv ---*/
	  fprintf(vrr_code, " deriv_build_",
		  am_letter[vrr_nodes[j].A],am_letter[vrr_nodes[j].B],
		  am_letter[vrr_nodes[j].C],am_letter[vrr_nodes[j].D]);
	  get_deriv_indices(&vrr_nodes[j],&di,&dj);
	  if (vrr_nodes[j].deriv_lvl == 2 && vrr_nodes[j].target == 0)
	    punt("There's a non-target second derivative class\n");

	  switch (di) {
	  case 0: case 1: case 2:
	      fprintf(vrr_code, "A%c_%c(Data,%d,",cart_comp[di],am_letter[vrr_nodes[j].A],
		      io[vrr_nodes[j].B]*io[vrr_nodes[j].C]*io[vrr_nodes[j].D]);
	      /* Add this function to the list of inlined functions if necessary */
	      max_node_am = vrr_nodes[j].A;
	      if (to_inline_deriv && max_node_am <= Params.max_am_to_inline_deriv_worker)
		fprintf(inline_vrr_list,"#include \"deriv_build_A%c_%c.h\"\n",
		        cart_comp[di], am_letter[vrr_nodes[j].A]);
	      break;

	  case 3: case 4: case 5:
	      fprintf(vrr_code, "B%c_%c(Data,%d,%d,",cart_comp[di-3],am_letter[vrr_nodes[j].B],
		      io[vrr_nodes[j].A],io[vrr_nodes[j].C]*io[vrr_nodes[j].D]);
	      /* Add this function to the list of inlined functions if necessary */
	      max_node_am = vrr_nodes[j].B;
	      if (to_inline_deriv && max_node_am <= Params.max_am_to_inline_deriv_worker)
		fprintf(inline_vrr_list,"#include \"deriv_build_B%c_%c.h\"\n",
		        cart_comp[di-3], am_letter[vrr_nodes[j].B]);
	      break;

	  case 6: case 7: case 8:
	      fprintf(vrr_code, "C%c_%c(Data,%d,%d,",cart_comp[di-6],am_letter[vrr_nodes[j].C],
		      io[vrr_nodes[j].A]*io[vrr_nodes[j].B],io[vrr_nodes[j].D]);
	      /* Add this function to the list of inlined functions if necessary */
	      max_node_am = vrr_nodes[j].C;
	      if (to_inline_deriv && max_node_am <= Params.max_am_to_inline_deriv_worker)
		fprintf(inline_vrr_list,"#include \"deriv_build_C%c_%c.h\"\n",
		        cart_comp[di-6], am_letter[vrr_nodes[j].C]);
	      break;

	  case 9: case 10: case 11:
	      fprintf(vrr_code, "D%c_%c(Data,%d,",cart_comp[di-9],am_letter[vrr_nodes[j].D],
		      io[vrr_nodes[j].A]*io[vrr_nodes[j].B]*io[vrr_nodes[j].C]);
	      /* Add this function to the list of inlined functions if necessary */
	      max_node_am = vrr_nodes[j].D;
	      if (to_inline_deriv && max_node_am <= Params.max_am_to_inline_deriv_worker)
		fprintf(inline_vrr_list,"#include \"deriv_build_D%c_%c.h\"\n",
		        cart_comp[di-9], am_letter[vrr_nodes[j].D]);
	      break;
	  }
	  
	  fprintf(vrr_code, "dvrr_stack+%d", vrr_nodes[j].pointer);
	  for(k=0; k<2; k++){
	    if(vrr_nodes[j].children[k] > 0)
	      fprintf(vrr_code, ", dvrr_stack+%d", vrr_nodes[vrr_nodes[j].children[k]].pointer);
	    else if (vrr_nodes[j].children[k] == NONODE)
	      fprintf(vrr_code, ", NULL");
	    else
	      fprintf(vrr_code, ", Data->F+%d", (-1)*vrr_nodes[j].children[k]);
	  }
	  fprintf(vrr_code, ");\n");
	}
	else if (vrr_nodes[j].B + vrr_nodes[j].D > 0) { /*--- build_hrr ---*/
	  if (vrr_nodes[j].B == 0 && vrr_nodes[j].D != 0) {
	    fprintf(vrr_code, " hrr3_build_%c%c(Libderiv->CD,dvrr_stack+%d,dvrr_stack+%d,",
		    am_letter[vrr_nodes[j].C], am_letter[vrr_nodes[j].D], vrr_nodes[j].pointer,
		    vrr_nodes[vrr_nodes[j].children[0]].pointer);
	    if (vrr_nodes[j].children[1] > 0)
	      fprintf(vrr_code, "dvrr_stack+%d,%d);\n\n", vrr_nodes[vrr_nodes[j].children[1]].pointer,
		      io[vrr_nodes[j].A]*io[vrr_nodes[j].B]);
	    else
	      fprintf(vrr_code, "Data->F,%d);\n\n", io[vrr_nodes[j].A]*io[vrr_nodes[j].B]);
	    /* Add this function to the list of inlined functions if necessary */
	    max_node_am = (vrr_nodes[j].C > vrr_nodes[j].D) ? vrr_nodes[j].C : vrr_nodes[j].D;
	    if (to_inline_hrr && max_node_am <= Params.max_am_to_inline_hrr_worker)
	      fprintf(inline_vrr_list,"#include <libint/hrr3_build_%c%c.h>\n", am_letter[vrr_nodes[j].C], am_letter[vrr_nodes[j].D]);
	  }
	  else if (vrr_nodes[j].B != 0) {
	    fprintf(vrr_code, " hrr1_build_%c%c(Libderiv->AB,dvrr_stack+%d,dvrr_stack+%d,",
		    am_letter[vrr_nodes[j].A], am_letter[vrr_nodes[j].B], vrr_nodes[j].pointer,
		    vrr_nodes[vrr_nodes[j].children[0]].pointer);
	    if (vrr_nodes[j].children[1] > 0)
	      fprintf(vrr_code, "dvrr_stack+%d,%d);\n\n", vrr_nodes[vrr_nodes[j].children[1]].pointer,
		      io[vrr_nodes[j].C]*io[vrr_nodes[j].D]);
	    else
	      fprintf(vrr_code, "Data->F,%d);\n", io[vrr_nodes[j].C]*io[vrr_nodes[j].D]);
	    /* Add this function to the list of inlined functions if necessary */
	    max_node_am = (vrr_nodes[j].A > vrr_nodes[j].B) ? vrr_nodes[j].A : vrr_nodes[j].B;
	    if (to_inline_hrr && max_node_am <= Params.max_am_to_inline_hrr_worker)
	      fprintf(inline_vrr_list,"#include <libint/hrr1_build_%c%c.h>\n", am_letter[vrr_nodes[j].A], am_letter[vrr_nodes[j].B]);
	  }
	}
	else { /*--- build_vrr ---*/
	  if (vrr_nodes[j].A <= LIBINT_OPT_AM && vrr_nodes[j].C <= LIBINT_OPT_AM) {
	    fprintf(vrr_code, " _BUILD_%c0%c0(Data,", am_letter[vrr_nodes[j].A], am_letter[vrr_nodes[j].C]);
	    /* Add this function to the list of inlined functions if necessary */
	    max_node_am = (vrr_nodes[j].A > vrr_nodes[j].C) ? vrr_nodes[j].A : vrr_nodes[j].C;
	    if (to_inline_vrr && max_node_am <= Params.max_am_to_inline_vrr_worker)
	      fprintf(inline_vrr_list,"#include <libint/build_%c0%c0.h>\n", am_letter[vrr_nodes[j].A], am_letter[vrr_nodes[j].C]);
	  }
	  else {
	    fprintf(vrr_code, " am[0] = %d;  am[1] = %d;\n", vrr_nodes[j].A, vrr_nodes[j].C);
	    fprintf(vrr_code, " vrr_build_xxxx(am,Data,");
	  }
	  fprintf(vrr_code, "dvrr_stack+%d", vrr_nodes[j].pointer);
	  for(k=0; k<5; k++){
	    if(vrr_nodes[j].children[k] > 0)
	      fprintf(vrr_code, ", dvrr_stack+%d", vrr_nodes[vrr_nodes[j].children[k]].pointer);
	    else if (vrr_nodes[j].children[k] == NONODE)
	      fprintf(vrr_code, ", NULL");
	    else
	      fprintf(vrr_code, ", Data->F+%d", (-1)*vrr_nodes[j].children[k]);
	  }
	  fprintf(vrr_code, ");\n");
	}

	/*-----------------------------------------------
	  If this derivative class is one of the targets
	  copy it to a location pointed by deriv_classes
	  to be used by the calling hrr_order routine
	 -----------------------------------------------*/
	if (vrr_nodes[j].target == 1) {
	  fprintf(vrr_code, " tmp = dvrr_stack + %d;\n", vrr_nodes[j].pointer);
	  get_deriv_indices(&vrr_nodes[j],&di,&dj);
	  if (vrr_nodes[j].deriv_lvl == 0)
	    fprintf(vrr_code, " target_ptr = Libderiv->dvrr_classes[%d][%d];\n",
		    vrr_nodes[j].A, vrr_nodes[j].C);
	  else if (vrr_nodes[j].deriv_lvl == 1)
	    fprintf(vrr_code, " target_ptr = Libderiv->deriv_classes[%d][%d][%d];\n",
		    vrr_nodes[j].A, vrr_nodes[j].C, di);
	  else if (vrr_nodes[j].deriv_lvl == 2)
	    fprintf(vrr_code, " target_ptr = Libderiv->deriv2_classes[%d][%d][%d];\n",
		    vrr_nodes[j].A, vrr_nodes[j].C, di*12+dj);
	  fprintf(vrr_code, " for(i=0;i<%d;i++)\n",vrr_nodes[j].size);
	  fprintf(vrr_code, "   target_ptr[i] += tmp[i];\n\n");
	}
	else
	  fprintf(vrr_code, "\n");

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
      
      /* compare this max_stack_size to the libint_stack_size for this angular momentum */
      if (libderiv12_stack_size[current_highest_am] < max_stack_size)
	libderiv12_stack_size[current_highest_am] = max_stack_size;

      max_stack_size = 0;
      }

    }
    }
  }
  return;
}


/* Recursive function that build the derivative HRR subgraph given the parent */

static int mk_dhrr_node(class node, class *allnodes, int new)
{

  int i, j, k, l;
  class O[8];
  int subnodes = 0;
  int thisnode;
  int rlink, llink;
  int made = 0;
  static int io[] = {1,3,6,10,15,21,28,36,45,55,66,78,91,105,120,136,153,171,190,210};

  if (node.A == 0 && node.B == 0 && node.C == 0 && node.D == 0 && node.deriv_lvl == 0)
    return -1;

  /* Search for the parent node on stack
     If it's not there - we'll add it to the end of the stack */
  thisnode = last_hrr_node;
  /* it's already placed on the stack allnodes - make sure children don't get created again (made = 1) */
/*  if (node.deriv_lvl != 0) {*/
    for(i=0;i<last_hrr_node;i++)
      if (allnodes[i].deriv_lvl == node.deriv_lvl &&
	  allnodes[i].A == node.A &&
	  allnodes[i].B == node.B &&
	  allnodes[i].C == node.C &&
	  allnodes[i].D == node.D &&
	  allnodes[i].deriv_ind[0] == node.deriv_ind[0] &&
	  allnodes[i].deriv_ind[1] == node.deriv_ind[1] &&
	  allnodes[i].deriv_ind[2] == node.deriv_ind[2] &&
	  allnodes[i].deriv_ind[3] == node.deriv_ind[3] &&
	  allnodes[i].deriv_ind[4] == node.deriv_ind[4] &&
	  allnodes[i].deriv_ind[5] == node.deriv_ind[5] &&
	  allnodes[i].deriv_ind[6] == node.deriv_ind[6] &&
	  allnodes[i].deriv_ind[7] == node.deriv_ind[7] &&
	  allnodes[i].deriv_ind[8] == node.deriv_ind[8] &&
	  allnodes[i].deriv_ind[9] == node.deriv_ind[9] &&
	  allnodes[i].deriv_ind[10] == node.deriv_ind[10] &&
	  allnodes[i].deriv_ind[11] == node.deriv_ind[11]) {
	thisnode = i;
	made = 1;
	break;
      }
/*  }*/

  /* it's not computed, add it, and make it the first to compute! */
  if(!made){
    allnodes[thisnode].A = node.A;
    allnodes[thisnode].B = node.B;
    allnodes[thisnode].C = node.C;
    allnodes[thisnode].D = node.D;
    allnodes[thisnode].m = node.m;
    allnodes[thisnode].deriv_lvl = node.deriv_lvl;
    allnodes[thisnode].deriv_ind[0] = node.deriv_ind[0];
    allnodes[thisnode].deriv_ind[1] = node.deriv_ind[1];
    allnodes[thisnode].deriv_ind[2] = node.deriv_ind[2];
    allnodes[thisnode].deriv_ind[3] = node.deriv_ind[3];
    allnodes[thisnode].deriv_ind[4] = node.deriv_ind[4];
    allnodes[thisnode].deriv_ind[5] = node.deriv_ind[5];
    allnodes[thisnode].deriv_ind[6] = node.deriv_ind[6];
    allnodes[thisnode].deriv_ind[7] = node.deriv_ind[7];
    allnodes[thisnode].deriv_ind[8] = node.deriv_ind[8];
    allnodes[thisnode].deriv_ind[9] = node.deriv_ind[9];
    allnodes[thisnode].deriv_ind[10] = node.deriv_ind[10];
    allnodes[thisnode].deriv_ind[11] = node.deriv_ind[11];
    allnodes[thisnode].num_parents = 0;
    allnodes[thisnode].parents_counter = 0;
    allnodes[thisnode].marked = 0;
    allnodes[thisnode].pointer = 0;
    memset(allnodes[thisnode].parents,0,NUMPARENTS*sizeof(int));
    allnodes[thisnode].children[0] = NONODE;
    allnodes[thisnode].children[1] = NONODE;
    allnodes[thisnode].children[2] = NONODE;
    allnodes[thisnode].children[3] = NONODE;
    allnodes[thisnode].children[4] = NONODE;
    allnodes[thisnode].children[5] = NONODE;
    allnodes[thisnode].children[6] = NONODE;
    allnodes[thisnode].children[7] = NONODE;
    allnodes[thisnode].size = io[node.A]*io[node.B]*io[node.C]*io[node.D];
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
    O[0].deriv_lvl = node.deriv_lvl;
    memcpy(O[0].deriv_ind,node.deriv_ind,12*sizeof(int));
    allnodes[thisnode].children[0] = 
            mk_dhrr_node(O[0], allnodes, made);
    O[1].A = node.A;
    O[1].B = node.B-1;
    O[1].C = node.C;
    O[1].D = node.D;
    O[1].m = node.m;
    O[1].deriv_lvl = node.deriv_lvl;
    memcpy(O[1].deriv_ind,node.deriv_ind,12*sizeof(int));
    allnodes[thisnode].children[1] =
	    mk_dhrr_node(O[1], allnodes, made);
    for(i=0;i<6;i++)
      if (node.deriv_ind[i]) {
	O[2+i].A = node.A;
	O[2+i].B = node.B-1;
	O[2+i].C = node.C;
	O[2+i].D = node.D;
	O[2+i].m = node.m;
	O[2+i].deriv_lvl = node.deriv_lvl-1;
	memcpy(O[2+i].deriv_ind,node.deriv_ind,12*sizeof(int));
	O[2+i].deriv_ind[i]--;
	allnodes[thisnode].children[2+i] =
	        mk_dhrr_node(O[2+i], allnodes, made);
      }
  }
  else if(node.D){
    O[0].A = node.A;
    O[0].B = node.B;
    O[0].C = node.C+1;
    O[0].D = node.D-1;
    O[0].m = node.m;
    O[0].deriv_lvl = node.deriv_lvl;
    memcpy(O[0].deriv_ind,node.deriv_ind,12*sizeof(int));
    allnodes[thisnode].children[0] = 
            mk_dhrr_node(O[0], allnodes, made);
    O[1].A = node.A;
    O[1].B = node.B;
    O[1].C = node.C;
    O[1].D = node.D-1;
    O[1].m = node.m;
    O[1].deriv_lvl = node.deriv_lvl;
    memcpy(O[1].deriv_ind,node.deriv_ind,12*sizeof(int));
    allnodes[thisnode].children[1] = 
            mk_dhrr_node(O[1], allnodes, made);
    for(i=0;i<6;i++)
      if (node.deriv_ind[i+6]) {
	O[2+i].A = node.A;
	O[2+i].B = node.B;
	O[2+i].C = node.C;
	O[2+i].D = node.D-1;
	O[2+i].m = node.m;
	O[2+i].deriv_lvl = node.deriv_lvl-1;
	memcpy(O[2+i].deriv_ind,node.deriv_ind,12*sizeof(int));
	O[2+i].deriv_ind[i+6]--;
	allnodes[thisnode].children[2+i] =
	        mk_dhrr_node(O[2+i], allnodes, made);
      }
  }
  }

  return thisnode;

}


/* Recursive function that builds a hybrid HRR/VRR/deriv subgraph given the parent */

static int mk_deriv_node(class node, class *allnodes, int new)
{

  int i, j, k, l;
  int di;
  class O[5];
  int subnodes = 0;
  int thisnode;
  int rlink, llink;
  int made = 0;
  static int io[] = {1,3,6,10,15,21,28,36,45,55,66,78,91,105,120,136,153,171,190,210};

  /* If it's not a derivative class - do some checks to see if need to proceed */
  if (node.deriv_lvl == 0 && node.A + node.B + node.C + node.D == 0)
    return (-1)*node.m;

  /* Search for the parent node on stack
     If it's not there - we'll add it to the end of the stack */
  thisnode = last_vrr_node;
  /* it's already placed on the stack allnodes - make sure children don't get created again (made = 1) */
  if (node.deriv_lvl != 0) {
    for(i=0;i<last_vrr_node;i++)
      if (allnodes[i].deriv_lvl == node.deriv_lvl &&
	  allnodes[i].A == node.A &&
	  allnodes[i].B == node.B &&
	  allnodes[i].C == node.C &&
	  allnodes[i].D == node.D &&
	  allnodes[i].deriv_ind[0] == node.deriv_ind[0] &&
	  allnodes[i].deriv_ind[1] == node.deriv_ind[1] &&
	  allnodes[i].deriv_ind[2] == node.deriv_ind[2] &&
	  allnodes[i].deriv_ind[3] == node.deriv_ind[3] &&
	  allnodes[i].deriv_ind[4] == node.deriv_ind[4] &&
	  allnodes[i].deriv_ind[5] == node.deriv_ind[5] &&
	  allnodes[i].deriv_ind[6] == node.deriv_ind[6] &&
	  allnodes[i].deriv_ind[7] == node.deriv_ind[7] &&
	  allnodes[i].deriv_ind[8] == node.deriv_ind[8] &&
	  allnodes[i].deriv_ind[9] == node.deriv_ind[9] &&
	  allnodes[i].deriv_ind[10] == node.deriv_ind[10] &&
	  allnodes[i].deriv_ind[11] == node.deriv_ind[11]) {
	thisnode = i;
	made = 1;
	break;
      }
  }
  else if (node.B + node.D != 0) {
    if (hrr_hash_table[node.A][node.B][node.C][node.D]) {
      i = hrr_hash_table[node.A][node.B][node.C][node.D] - 1;
      thisnode = i;
      made = 1;
    }
  }
  else if (vrr_hash_table[node.A][node.C][node.m]) {
    i = vrr_hash_table[node.A][node.C][node.m] - 1;
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
    allnodes[thisnode].deriv_lvl = node.deriv_lvl;
    allnodes[thisnode].deriv_ind[0] = node.deriv_ind[0];
    allnodes[thisnode].deriv_ind[1] = node.deriv_ind[1];
    allnodes[thisnode].deriv_ind[2] = node.deriv_ind[2];
    allnodes[thisnode].deriv_ind[3] = node.deriv_ind[3];
    allnodes[thisnode].deriv_ind[4] = node.deriv_ind[4];
    allnodes[thisnode].deriv_ind[5] = node.deriv_ind[5];
    allnodes[thisnode].deriv_ind[6] = node.deriv_ind[6];
    allnodes[thisnode].deriv_ind[7] = node.deriv_ind[7];
    allnodes[thisnode].deriv_ind[8] = node.deriv_ind[8];
    allnodes[thisnode].deriv_ind[9] = node.deriv_ind[9];
    allnodes[thisnode].deriv_ind[10] = node.deriv_ind[10];
    allnodes[thisnode].deriv_ind[11] = node.deriv_ind[11];
    if (node.deriv_lvl == 0) {
      if (node.B + node.D == 0)
	vrr_hash_table[node.A][node.C][node.m] = thisnode + 1;
      else
	hrr_hash_table[node.A][node.B][node.C][node.D] = thisnode + 1;
    }
    allnodes[thisnode].num_parents = 0;
    allnodes[thisnode].parents_counter = 0;
    allnodes[thisnode].marked = 0;
    allnodes[thisnode].target = 0;
    allnodes[thisnode].pointer = 0;
    memset(allnodes[thisnode].parents,0,NUMPARENTS*sizeof(int));
    allnodes[thisnode].children[0] = NONODE;
    allnodes[thisnode].children[1] = NONODE;
    allnodes[thisnode].children[2] = NONODE;
    allnodes[thisnode].children[3] = NONODE;
    allnodes[thisnode].children[4] = NONODE;
    allnodes[thisnode].size = io[node.A]*io[node.B]*io[node.C]*io[node.D];
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
    /* derivative ERI */
    if (node.deriv_lvl) {
      for(i=0;i<12;i++) {
	if (node.deriv_ind[i] != 0) {
	  switch (i) {
	  case 0: case 1: case 2:
	      O[0].A = node.A+1;
	      O[0].B = node.B;
	      O[0].C = node.C;
	      O[0].D = node.D;
	      O[0].m = 0;
	      O[0].deriv_lvl = node.deriv_lvl-1;
	      for(di=0;di<12;di++)
		O[0].deriv_ind[di] = node.deriv_ind[di];
	      O[0].deriv_ind[i] = node.deriv_ind[i]-1;
	      allnodes[thisnode].children[0] = mk_deriv_node(O[0], allnodes, made);
	      if (node.A > 0) {
		O[1].A = node.A-1;
		O[1].B = node.B;
		O[1].C = node.C;
		O[1].D = node.D;
		O[1].m = 0;
		O[1].deriv_lvl = node.deriv_lvl-1;
		for(di=0;di<12;di++)
		  O[1].deriv_ind[di] = node.deriv_ind[di];
		O[1].deriv_ind[i] = node.deriv_ind[i]-1;
		allnodes[thisnode].children[1] = mk_deriv_node(O[1], allnodes, made);
	      }
	      break;

	  case 3: case 4: case 5:
	      O[0].A = node.A;
	      O[0].B = node.B+1;
	      O[0].C = node.C;
	      O[0].D = node.D;
	      O[0].m = 0;
	      O[0].deriv_lvl = node.deriv_lvl-1;
	      for(di=0;di<12;di++)
		O[0].deriv_ind[di] = node.deriv_ind[di];
	      O[0].deriv_ind[i] = node.deriv_ind[i]-1;
	      allnodes[thisnode].children[0] = mk_deriv_node(O[0], allnodes, made);
	      if (node.B > 0) {
		O[1].A = node.A;
		O[1].B = node.B-1;
		O[1].C = node.C;
		O[1].D = node.D;
		O[1].m = 0;
		O[1].deriv_lvl = node.deriv_lvl-1;
		for(di=0;di<12;di++)
		  O[1].deriv_ind[di] = node.deriv_ind[di];
		O[1].deriv_ind[i] = node.deriv_ind[i]-1;
		allnodes[thisnode].children[1] = mk_deriv_node(O[1], allnodes, made);
	      }
	      break;

	  case 6: case 7: case 8:
	      O[0].A = node.A;
	      O[0].B = node.B;
	      O[0].C = node.C+1;
	      O[0].D = node.D;
	      O[0].m = 0;
	      O[0].deriv_lvl = node.deriv_lvl-1;
	      for(di=0;di<12;di++)
		O[0].deriv_ind[di] = node.deriv_ind[di];
	      O[0].deriv_ind[i] = node.deriv_ind[i]-1;
	      allnodes[thisnode].children[0] = mk_deriv_node(O[0], allnodes, made);
	      if (node.C > 0) {
		O[1].A = node.A;
		O[1].B = node.B;
		O[1].C = node.C-1;
		O[1].D = node.D;
		O[1].m = 0;
		O[1].deriv_lvl = node.deriv_lvl-1;
		for(di=0;di<12;di++)
		  O[1].deriv_ind[di] = node.deriv_ind[di];
		O[1].deriv_ind[i] = node.deriv_ind[i]-1;
		allnodes[thisnode].children[1] = mk_deriv_node(O[1], allnodes, made);
	      }
	      break;

	  case 9: case 10: case 11:
	      O[0].A = node.A;
	      O[0].B = node.B;
	      O[0].C = node.C;
	      O[0].D = node.D+1;
	      O[0].m = 0;
	      O[0].deriv_lvl = node.deriv_lvl-1;
	      for(di=0;di<12;di++)
		O[0].deriv_ind[di] = node.deriv_ind[di];
	      O[0].deriv_ind[i] = node.deriv_ind[i]-1;
	      allnodes[thisnode].children[0] = mk_deriv_node(O[0], allnodes, made);
	      if (node.D > 0) {
		O[1].A = node.A;
		O[1].B = node.B;
		O[1].C = node.C;
		O[1].D = node.D-1;
		O[1].m = 0;
		O[1].deriv_lvl = node.deriv_lvl-1;
		for(di=0;di<12;di++)
		  O[1].deriv_ind[di] = node.deriv_ind[di];
		O[1].deriv_ind[i] = node.deriv_ind[i]-1;
		allnodes[thisnode].children[1] = mk_deriv_node(O[1], allnodes, made);
	      }
	      break;
	  }
	  break;
	}
      }
    }
    /* HRR case */
    else if (node.B + node.D != 0) {
      if(node.B){
	O[0].A = node.A+1;
	O[0].B = node.B-1;
	O[0].C = node.C;
	O[0].D = node.D;
	O[0].m = node.m;
	O[0].deriv_lvl = 0;
	memset(O[0].deriv_ind,0,12*sizeof(int));
	allnodes[thisnode].children[0] = mk_deriv_node(O[0], allnodes, made);
	O[1].A = node.A;
	O[1].B = node.B-1;
	O[1].C = node.C;
	O[1].D = node.D;
	O[1].m = node.m;
	O[1].deriv_lvl = 0;
	memset(O[1].deriv_ind,0,12*sizeof(int));
	allnodes[thisnode].children[1] = mk_deriv_node(O[1], allnodes, made);
      }
      else if(node.D){
	O[0].A = node.A;
	O[0].B = node.B;
	O[0].C = node.C+1;
	O[0].D = node.D-1;
	O[0].m = node.m;
	O[0].deriv_lvl = 0;
	memset(O[0].deriv_ind,0,12*sizeof(int));
	allnodes[thisnode].children[0] = mk_deriv_node(O[0], allnodes, made);
	O[1].A = node.A;
	O[1].B = node.B;
	O[1].C = node.C;
	O[1].D = node.D-1;
	O[1].m = node.m;
	O[1].deriv_lvl = 0;
	memset(O[1].deriv_ind,0,12*sizeof(int));
	allnodes[thisnode].children[1] = mk_deriv_node(O[1], allnodes, made);
      }
    }
    /* VRR case */
    else {
      if(node.A){
	O[0].A = node.A-1;
	O[0].B = 0;
	O[0].C = node.C;
	O[0].D = 0;
	O[0].m = node.m;
	O[0].deriv_lvl = 0;
	memset(O[0].deriv_ind,0,12*sizeof(int));
	allnodes[thisnode].children[0] = mk_deriv_node(O[0], allnodes, made);
	O[1].A = node.A-1;
	O[1].B = 0;
	O[1].C = node.C;
	O[1].D = 0;
	O[1].m = node.m+1;
	O[1].deriv_lvl = 0;
	memset(O[1].deriv_ind,0,12*sizeof(int));
	allnodes[thisnode].children[1] = mk_deriv_node(O[1], allnodes, made);
	if(node.A>1){
	  O[2].A = node.A-2;
	  O[2].B = 0;
	  O[2].C = node.C;
	  O[2].D = 0;
	  O[2].m = node.m;
	  O[2].deriv_lvl = 0;
	  memset(O[2].deriv_ind,0,12*sizeof(int));
	  allnodes[thisnode].children[2] = mk_deriv_node(O[2], allnodes, made);
	  O[3].A = node.A-2;
	  O[3].B = 0;
	  O[3].C = node.C;
	  O[3].D = 0;
	  O[3].m = node.m+1;
	  O[3].deriv_lvl = 0;
	  memset(O[3].deriv_ind,0,12*sizeof(int));
	  allnodes[thisnode].children[3] = mk_deriv_node(O[3], allnodes, made);
	}
	if(node.C){
	  O[4].A = node.A-1;
	  O[4].B = 0;
	  O[4].C = node.C-1;
	  O[4].D = 0;
	  O[4].m = node.m+1;
	  O[4].deriv_lvl = 0;
	  memset(O[4].deriv_ind,0,12*sizeof(int));
	  allnodes[thisnode].children[4] = mk_deriv_node(O[4], allnodes, made);
	}
      }
      else if(node.C){
	O[0].A = node.A;
	O[0].B = 0;
	O[0].C = node.C-1;
	O[0].D = 0;
	O[0].m = node.m;
	O[0].deriv_lvl = 0;
	memset(O[0].deriv_ind,0,12*sizeof(int));
	allnodes[thisnode].children[0] = mk_deriv_node(O[0], allnodes, made);
	O[1].A = node.A;
	O[1].B = 0;
	O[1].C = node.C-1;
	O[1].D = 0;
	O[1].m = node.m+1;
	O[1].deriv_lvl = 0;
	memset(O[1].deriv_ind,0,12*sizeof(int));
	allnodes[thisnode].children[1] = mk_deriv_node(O[1], allnodes, made);
	if(node.C>1){
	  O[2].A = node.A;
	  O[2].B = 0;
	  O[2].C = node.C-2;
	  O[2].D = 0;
	  O[2].m = node.m;
	  O[2].deriv_lvl = 0;
	  memset(O[2].deriv_ind,0,12*sizeof(int));
	  allnodes[thisnode].children[2] = mk_deriv_node(O[2], allnodes, made);
	  O[3].A = node.A;
	  O[3].B = 0;
	  O[3].C = node.C-2;
	  O[3].D = 0;
	  O[3].m = node.m+1;
	  O[3].deriv_lvl = 0;
	  memset(O[3].deriv_ind,0,12*sizeof(int));
	  allnodes[thisnode].children[3] = mk_deriv_node(O[3], allnodes, made);
	}
      }
    }
  }

  return thisnode;

}



/* Make hrr_nodes[rent] a parent of hrr_nodes[n] and proceed recursively */

static void mark_dhrr_parents(int n, class *allnodes, int rent)
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

    for(i=0; i<8; i++)
      if(allnodes[n].children[i]>0)
	mark_dhrr_parents(allnodes[n].children[i], allnodes, n);
  }
  return;
}


/* Make vrr_nodes[rent] a parent of vrr_nodes[n] and proceed recursively */

static void mark_vrr_parents(int n, class *allnodes, int rent)
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


static void mark_parents(int n, class *allnodes, int rent)
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

static int alloc_mem_dhrr(class *nodes)
{
  int i, j, k, l;
  int size;
  int child;
  int free_it;
  static int io[] = {1,3,6,10,15,21,28,36,45,55,66,78,91,105,120,136,153,171,190,210};

  j = first_hrr_to_compute;
  do{
    /* Node to compute */
    if (nodes[j].marked == 0) {
      nodes[j].marked = 1; /* sign that it has been passed */
      nodes[j].pointer = get_mem(nodes[j].size); /* Allocate memory for it on a CINTS-provided stack */
    }
      
    /* Figure out which children can be freed,
       i.e. which children are not targets and have all parents marked */
    for(k=0; k<8; k++){
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

static int alloc_mem_vrr(class *nodes)
{
  int i, j, k, l;
  int size;
  int child;
  int free_it;
  static int io[] = {1,3,6,10,15,21,28,36,45,55,66,78,91,105,120,136,153,171,190,210};

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
    if (j >= last_vrr_node)
      punt("Linked list goes out of boundary\n");
  } while (j != -1);
  
  return nodes[0].pointer;
}


static void get_deriv_indices(class *node,int *di, int *dj)
{
  int i, j;
  i = j = 0;
  
  switch(node->deriv_lvl) {

  case 0:
    break;

  case 1:
    for(i=0;i<12;i++)
      if (node->deriv_ind[i] != 0) {
	break;
      }
    i = i;
    j = 0;
    break;

  case 2:
    for(i=0;i<12;i++)
      if (node->deriv_ind[i] != 0)
	if (node->deriv_ind[i] == 2) {
	  j = i;
	  break;
	}
	else {
	  for(j=i+1;j<12;j++)
	    if (node->deriv_ind[j] != 0)
	      break;
	  break;
	}
    break;
  }
  
  *di = i; *dj = j;
  return;
}
