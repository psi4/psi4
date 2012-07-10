/*! \file
**  \ingroup DETCI
**  \brief Code to form the CI space as strings with all single replacements.
**
** C. David Sherrill
** Center for Computational Quantum Chemistry
** University of Georgia
** June 1995
**
*/

#include <cstdlib>
#include <cstdio>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "structs.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace detci {

/* DEFINES */
#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))
#define MAXIJ 100000
extern unsigned char ***Occs;


/* GLOBALS for this module */
struct level *sbgr_tr_head;
int sbgr_tr_orbs;
int **sbgr_tr_out;
int sbgr_tr_walks;
int *sbgr_tr_alist;

int *O, *U, *T ;
int **Tij, **Toij, *Tcnt;
unsigned int **Tidx;
signed char **Tsgn;


extern int subgr_lex_addr(struct level *head, int *occs, int nel, int norb);
extern int og_lex_addr(struct olsen_graph *Graph, int *occs, int nel,
      int *listnum);

/* FUNCTION PROTOTYPES this module */
void subgr_trav_init(struct level *head, int ci_orbs, int **outarr, int
   walks);
void subgr_traverse(int i, int j);
void form_stringwr(struct stringwr *strlist, int *occs, int N,
      int num_ci_orbs, struct stringgraph *subgraph, struct olsen_graph
      *Graph, int first_orb_active);
void og_form_repinfo(struct stringwr *string, int num_ci_orbs, 
      struct olsen_graph *Graph, int first_orb_active);
void init_stringwr_temps(int nel, int num_ci_orbs, int nsym);
void free_stringwr_temps(int nsym);




/*
** stringlist():  This function forms the list of strings with their
**    single replacements using the Olsen Graph structures.
**
*/
void stringlist(struct olsen_graph *Graph, struct stringwr **slist)
{

   int i;
   int nirreps, irrep, code, ncodes, walk, listnum, addr;
   int nel_expl;
   struct stringgraph *subgraph;
   int **outarr, *occs;

   nel_expl = Graph->num_el_expl;
   ncodes = Graph->subgr_per_irrep;
   nirreps = Graph->nirreps;

   outarr = init_int_matrix(nel_expl, Graph->max_str_per_irrep);
   occs = init_int_array(nel_expl);

   if (!Parameters.repl_otf) {
      init_stringwr_temps(Graph->num_el_expl, Graph->num_orb, nirreps * ncodes);
      }

   Occs = (unsigned char ***) malloc (nirreps * ncodes * sizeof(unsigned
      char **));

   for (irrep=0,listnum=0; irrep < nirreps; irrep++) {
      for (code = 0; code < ncodes; code++,listnum++) { 

         Occs[listnum] = NULL;

         subgraph = Graph->sg[irrep] + code;
         if (!subgraph->num_strings) continue;

         Occs[listnum] = (unsigned char **) malloc (subgraph->num_strings * 
            sizeof(unsigned char *));
         for (i=0; i<subgraph->num_strings; i++) 
            Occs[listnum][i] = (unsigned char *) malloc (nel_expl * 
                                sizeof(unsigned char));

         slist[listnum] = (struct stringwr *)
            malloc (subgraph->num_strings * sizeof(struct stringwr));
   
         subgr_trav_init(subgraph->lvl, Graph->num_orb, outarr, 0);
         subgr_traverse(0, 0);

         for (walk=0; walk<subgraph->num_strings; walk++) {
            for (i=0; i<nel_expl; i++) occs[i]= outarr[i][walk];
            addr = subgr_lex_addr(subgraph->lvl,occs,nel_expl,Graph->num_orb); 

            if (addr < 0) {
               printf("(stringlist): Impossible string addr\n");
               }

            for (i=0; i<nel_expl; i++) 
               Occs[listnum][addr][i] = (unsigned char) occs[i];

            form_stringwr(slist[irrep * ncodes + code], occs, 
               nel_expl, Graph->num_orb, subgraph, Graph, 
               Graph->num_cor_orbs);
            }
         } /* end loop over subgraph codes */
      } /* end loop over irreps */

   /* free the stringwr scratch space */
   if (!Parameters.repl_otf) {
      free_stringwr_temps(nirreps * ncodes);
      }

   free_int_matrix(outarr);
   free(occs);
   free(sbgr_tr_alist);
}



void subgr_trav_init(struct level *head, int ci_orbs, int **outarr, int
   walks)
{
   sbgr_tr_head = head;
   sbgr_tr_orbs = ci_orbs;
   sbgr_tr_out = outarr;
   sbgr_tr_walks = walks;
   sbgr_tr_alist = init_int_array(ci_orbs+1);
}


void subgr_traverse(int i, int j)
{
   int m,n;
   int k0, k1;

   sbgr_tr_alist[i] = sbgr_tr_head[i].a[j];

   if (i == sbgr_tr_orbs) {
      for (m=1,n=0; m<=sbgr_tr_orbs; m++) {
         if (sbgr_tr_alist[m] != sbgr_tr_alist[m-1]) 
            sbgr_tr_out[n++][sbgr_tr_walks] = m-1;
         }
      sbgr_tr_walks++;
      return;
      }

   k0 = sbgr_tr_head[i].k[0][j];
   k1 = sbgr_tr_head[i].k[1][j];
   if (k0) subgr_traverse(i+1, k0-1);
   if (k1) subgr_traverse(i+1, k1-1);

}


/*
** form_stringwr(): Make the string with replacements list.
**    This version uses the Olsen Graph structures. 
**
** Parameters:
**   strlist          = list of strings (may be different for alpha and beta)
**   occs             = list of occupied orbitals (N long)
**   N                = number of electrons explicitly included
**   num_ci_orbs      = number of active CI orbitals
**   subgraph         = the subgraph including the given walk
**   Graph            = the olsengraph for the alpha or beta strings
**   first_orb_active = first alp/bet orb active
*/
void form_stringwr(struct stringwr *strlist, int *occs, int N,
      int num_ci_orbs, struct stringgraph *subgraph, struct olsen_graph
      *Graph, int first_orb_active)
{
   unsigned char *occlist;
   unsigned int addr;
   int i;
   struct stringwr *node;
   
   occlist = (unsigned char *) malloc (N * sizeof(unsigned char));
   if (occlist == NULL) {
      throw PsiException("(form_stringwr): Malloc error",__FILE__,__LINE__);
      }

   for (i=0; i<N; i++) {
      occlist[i] = occs[i] ;
      }

   addr = subgr_lex_addr(subgraph->lvl, occs, N, num_ci_orbs);
 
   node = strlist + addr;
   node->occs = occlist;

   if (!Parameters.repl_otf) {
      og_form_repinfo(node, num_ci_orbs, Graph, first_orb_active);
      }
}


/*
** og_form_repinfo(): This function is a slightly-modified version of the
**    form_repinfo() function, which forms the single replacement info 
**    for each string.  This version uses the Olsen Graph structures.
**
*/
void og_form_repinfo(struct stringwr *string, int num_ci_orbs, 
      struct olsen_graph *Graph, int first_orb_active)
{

   int nel, p, q, i, j, k, l, ij, oij;
   int nlists, listnum, strlistnum, nsym, jused, nfzc;
   int diagcnt=0; 
   static int *diagij = NULL;
   static int *diagoij = NULL; 
   unsigned int cnt, stringridx;
   int ridx;
   signed char sgn;
  
   nel = Graph->num_el_expl;
   nfzc = Graph->num_fzc_orbs;
   nsym = Graph->nirreps;
   nlists = Graph->subgr_per_irrep * nsym;

   /* Zero out the counting array and O, U, and T */
   zero_int_array(Tcnt, nlists);
   zero_int_array(O, nel+1);
   zero_int_array(U, num_ci_orbs-nel);
   zero_int_array(T, nel);

   for (i=0; i<nel; i++) O[i] = (string->occs)[i] ;
   stringridx = og_lex_addr(Graph, O, nel, &strlistnum);

   /* set up the ij indices for the 'diagonal' entries E_ii           */
   /* this assumes that the values in array O are strictly increasing */
   if (diagij == NULL) diagij = init_int_array(nel);
   if (diagoij == NULL) diagoij = init_int_array(nel);
   for (i=0; i<nel; i++) {
      j = O[i];
      diagij[i] = ioff[j] + j;
      diagoij[i] = j * num_ci_orbs + j;
      }

   /* do the inactive E_ii's first */
   for (i=0; i<first_orb_active; i++) {
      p = (string->occs)[i];
      q = (string->occs)[i];
      ij = ioff[p] + q;
      oij = p * num_ci_orbs + q;
      sgn = (signed char) 1;

      cnt = Tcnt[strlistnum];
      Tij[strlistnum][cnt] = ij;
      Toij[strlistnum][cnt] = oij;
      Tidx[strlistnum][cnt] = stringridx;
      Tsgn[strlistnum][cnt] = sgn;
      Tcnt[strlistnum] += 1; 
      diagcnt++;
      }

   /* now do the single replacements.  first form arrays O and U */

   for (i=0,j=0,k=0; i<num_ci_orbs; i++) {
      if (O[j] != i) U[k++] = i ;
      else j++ ;
      }

   for (i=0; i<nel; i++) {
      q = O[i];
      if (q < first_orb_active) continue ;

      for (j=0; j<num_ci_orbs - nel; j++) {
         p = U[j];
         ij = ioff[MAX0(p,q)] + MIN0(p,q);
         oij = p * num_ci_orbs + q;


         /* arrange occupied list in order */
         for (k=0, jused=0, l=0; k<nel; k++) {
           if (!jused && U[j] < O[k]) { jused=1; T[l++] = U[j]; }
           if (k != i) T[l++] = O[k] ;  /* was else if 8/27/94 */
           }
         if (!jused) T[l++] = U[j] ;
         ridx = og_lex_addr(Graph, T, nel, &listnum);
         
         if (ridx >= 0) {
            while (diagcnt < nel && ij > diagij[diagcnt]) {
               cnt = Tcnt[strlistnum];
               Tij[strlistnum][cnt] = diagij[diagcnt];
               Toij[strlistnum][cnt] = diagoij[diagcnt];
               Tidx[strlistnum][cnt] = stringridx;
               Tsgn[strlistnum][cnt] = (signed char) 1;
               Tcnt[strlistnum] += 1; 
               diagcnt++;
               }

            /* get the sign: 
             * there are i occupied orbitals before the one to be annihilated
             * contributing (-1)^i to the sign.  Need number of occupied orbs
             * (not counting any annihilated) before the one to be created
             * contributing (-1)^l to the sign.  Sign = (-1)^(i+l)
             */
            for (k=0,l=0; k<nel; k++) {
              if (O[k] < U[j] && k != i) l++;
              if (O[k] >= U[j]) break;
              }
            sgn = ((i + l) % 2) ? -1 : 1 ;

            cnt = Tcnt[listnum];
            Tij[listnum][cnt] = ij;
            Toij[listnum][cnt] = oij;
            Tidx[listnum][cnt] = ridx;
            Tsgn[listnum][cnt] = sgn;
            Tcnt[listnum] += 1;
            }
         }

      } /* end loop over electrons */

   while (diagcnt < nel) {
      sgn = (signed char) 1;
      listnum = strlistnum;
      cnt = Tcnt[strlistnum];
      Tij[strlistnum][cnt] = diagij[diagcnt];
      Toij[strlistnum][cnt] = diagoij[diagcnt];
      Tidx[strlistnum][cnt] = stringridx;
      Tsgn[strlistnum][cnt] = sgn;
      Tcnt[strlistnum] += 1; 
      diagcnt++;
      }

   /* now write the info in the T matrices */
   string->cnt = init_int_array(nlists);
   string->ij = (int **) malloc(sizeof(int *) * nlists);
   string->oij = (int **) malloc(sizeof(int *) * nlists);
   string->ridx = (unsigned int **) malloc(sizeof(unsigned int *) * nlists);
   string->sgn = (signed char **) malloc(sizeof(signed char *) * nlists);


   for (i=0; i<nlists; i++) {
      string->cnt[i] = cnt = Tcnt[i];
      string->ij[i] = NULL;
      string->oij[i] = NULL;
      string->ridx[i] = NULL;
      string->sgn[i] = NULL;
      if (cnt) {
         string->ij[i] = init_int_array(cnt);
         string->oij[i] = init_int_array(cnt);
         string->ridx[i] = (unsigned int *) malloc(cnt * 
            sizeof(unsigned int));
         string->sgn[i] = (signed char *) malloc(cnt * 
            sizeof(signed char));

         for (k=0; k<cnt; k++) {
            for (l=0,p=0,q=MAXIJ; l<cnt; l++) { 
               ij = Tij[i][l];
               if (ij <= q) { q = ij; p = l; }
               } 
            string->ij[i][k] = Tij[i][p];
            string->oij[i][k] = Toij[i][p];
            string->ridx[i][k] = Tidx[i][p];
            string->sgn[i][k] = Tsgn[i][p];
            Tij[i][p] = MAXIJ;
            }
         }
      } /* end loop over i */


}      


void init_stringwr_temps(int nel, int num_ci_orbs, int nsym)
{
   int maxcnt, i, j;

   O = init_int_array(nel+1) ;
 /*   U = init_int_array(num_ci_orbs - nel) ; */
   /* MLL and CDS +1 in case of no virtual alpha electrons */
   U = init_int_array(num_ci_orbs - nel + 1) ;
   T = init_int_array(nel) ;
   Tcnt = init_int_array(nsym);
   maxcnt = nel * num_ci_orbs; /* num single replacements inc. self-repl */
   Tij = (int **) malloc(sizeof(int *) * nsym);
   Toij = (int **) malloc(sizeof(int *) * nsym);
   Tidx = (unsigned int **) malloc(sizeof(unsigned int *) * nsym);
   Tsgn = (signed char **) malloc(sizeof(signed char *) * nsym);
   
   for (i=0; i<nsym; i++) {
      Tij[i] = init_int_array(maxcnt);
      Toij[i] = init_int_array(maxcnt);
      Tidx[i] = (unsigned int *) malloc(sizeof(unsigned int) * maxcnt);
      Tsgn[i] = (signed char *) malloc(sizeof(signed char) * maxcnt); 
      }
}


void free_stringwr_temps(int nsym)
{
   int i,j;

   free(O);
   free(U);
   free(T);
   free(Tcnt);
   for (i=0; i<nsym; i++) {
      free(Tij[i]);
      free(Toij[i]);
      free(Tidx[i]);
      free(Tsgn[i]);
      }
   free(Tij);
   free(Toij);
   free(Tidx);
   free(Tsgn);
}


}} // namespace psi::detci

