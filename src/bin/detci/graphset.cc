/*! \file
**  \ingroup DETCI
**  \brief Routines needed to maintain the GraphSet Object.  
**
** The graph for alpha/beta strings has a different subgraph for each value 
** of irrep and RAS I hole/RAS III particle/RAS IV particle combination.  
**
** C. David Sherrill, May 1996
** Based on previous code by David Sherrill, 1994-5
**
*/

#define EXTERN

#include <cstdlib>
#include <cstdio>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "structs.h"
#include "globals.h"
#include <iostream>
#include "odometer.h"

namespace psi { namespace detci {

extern struct stringwr **alplist;  
extern struct stringwr **betlist; 

extern void gs_stringlist(struct graph_set *AG, struct stringwr **slist);
extern void gs_set_ciblks(struct graph_set *AG, struct graph_set *BG);
extern void print_ci_space(struct stringwr *strlist, int num_strings,
   int nirreps, int strtypes, int nel, FILE *outfile);

#define UOC_BIT 1
#define OCC_BIT 2
#define N_RAS_SPACES 4

//#define DEBUG 

/* GLOBALS THIS MODULE */
int **Tij, **Toij;
unsigned int **Tridx;
signed char **Tsgn;


/* FUNCTION PROTOTYPES for this module */
void graphset(struct graph_set *GraphSet, int ci_orbs, int num_el, 
      int nirreps, int *orbsym, int ras1_lvl, int ras1_min, int ras1_max, 
      int ras3_lvl, int ras3_max, int num_fzc_orbs, int num_cor_orbs,
      int ras4_lvl, int ras4_max, int ras34_max);
void gs_add_walk(int ras1_idx, int ras3_num, int ras4_num, int *occs, 
      int nel_expl, int norb, int nirreps, int num_fzc_orbs,
      struct graph_set *GraphSet);
void gs_fill(int num_el, int norb, int nirreps, int num_fzc_orbs,
      struct graph_set *GraphSet);
int gs_glex_addr(struct fastgraph *Graph, int *occs, int nel);
void gs_print(struct graph_set *GraphSet, FILE *outfile);
void gs_stringlist(struct graph_set *GraphSet, struct stringwr **slist);
void gs_init_repinfo_temps(int nel, int norbs);
void gs_free_repinfo_temps(void);



/*
** formstrings(): Form the alpha and beta string graphs.  Use a different
**    graph for every irrep * (num el in RAS I) * (num el in RAS III)
**
*/
void formstrings(void)
{
   int i, nlists, nirreps, ncodes;
   int irrep, code, listnum;
   int *occs;

   AlphaGraph = new graph_set[1];

   // Make the graph
   graphset(AlphaGraph, CalcInfo.num_ci_orbs, CalcInfo.num_alp,
      CalcInfo.nirreps, CalcInfo.orbsym,
      Parameters.a_ras1_lvl, Parameters.a_ras1_min, Parameters.a_ras1_max,
      Parameters.ras3_lvl, Parameters.a_ras3_max,
      CalcInfo.num_fzc_orbs, CalcInfo.num_cor_orbs,
      Parameters.ras4_lvl, Parameters.a_ras4_max, Parameters.a_ras34_max);

   if (Parameters.print_lvl > 3)
      gs_print(AlphaGraph, outfile) ;

   ncodes = AlphaGraph->num_codes;
   nirreps = AlphaGraph->nirreps;
   nlists = AlphaGraph->num_graphs;
   alplist = (struct stringwr **) malloc(nlists * sizeof(struct stringwr *));
   for (i=0; i<nlists; i++) alplist[i] = NULL;

   gs_stringlist(AlphaGraph, alplist);

   if (Parameters.print_lvl>=4) {
      for (listnum=0; listnum<AlphaGraph->num_graphs; listnum++) {
         fprintf(outfile, "Strings for list %2d, irrep=%d\n",
            listnum, AlphaGraph->graph_irrep[listnum]);
         print_ci_space(alplist[listnum],
            AlphaGraph->Graph[listnum].num_strings,
            nirreps, nlists, AlphaGraph->num_el_expl, outfile) ;
         }
      }

   /* for beta string graph if necessary */
   if (CalcInfo.iopen) {

      BetaGraph = new graph_set[1];

      graphset(BetaGraph, CalcInfo.num_ci_orbs, CalcInfo.num_bet,
          CalcInfo.nirreps, CalcInfo.orbsym,
          Parameters.b_ras1_lvl,  Parameters.b_ras1_min, Parameters.b_ras1_max,
          Parameters.ras3_lvl, Parameters.b_ras3_max,
          CalcInfo.num_fzc_orbs, CalcInfo.num_cor_orbs,
          Parameters.ras4_lvl, Parameters.b_ras4_max, Parameters.b_ras3_max);

      if (Parameters.print_lvl > 1)
         og_print(BetaGraph, outfile) ;

      ncodes = BetaGraph->num_codes;
      nirreps = BetaGraph->nirreps;
      nlists = BetaGraph->num_graphs;
      betlist = (struct stringwr **) malloc(nlists * sizeof(struct stringwr *));
      for (i=0; i<nlists; i++) betlist[i] = NULL;

      gs_stringlist(BetaGraph, betlist);

      if (Parameters.print_lvl>=4) {
         for (listnum=0; listnum<AlphaGraph->num_graphs; listnum++) {
            fprintf(outfile, "Strings for list %2d, irrep=%d\n",
               listnum, AlphaGraph->graph_irrep[listnum]);
            print_ci_space(alplist[listnum],
               AlphaGraph->Graph[listnum].num_strings,
               nirreps, nlists, AlphaGraph->num_el_expl, outfile) ;
            }
         }
      } /* end if(iopen) */

   else {
      betlist = alplist;
      BetaGraph = AlphaGraph;
      }

   /* get number of alpha/beta strings, ref symmetry, etc */
   gs_set_ciblks(AlphaG, BetaG) ;

}



/*
** graphset(): Form an Olsen/Roos graph for the alpha/beta strings.
** 
** Unfortunately, I found it necessary to split the DRT graph into many
**    subgraphs, as described in Olsen & Roos 1988.  Each graph must 
**    contain walks of only one irrep, and with only one value for the 
**    number of electrons at the RAS I level.  All these graphs are 
**    stored in a collective container, struct olsen_graph.  For open
**    shell systems, there will be two of these, one for alpha and one
**    for beta.  For closed shell systems only the alpha graph is stored.
**
** Even more unfortunately, each subgraph must also be labeled according
**    to the number of electrons in RAS III, or so I believe I figured
**    out at one time.  That's how I'll do it now.  [CDS 5/95]
**
** Even more more unfortunately, I'm going to split it according to 
**    number of electrons in RAS IV!!!  That's just crazy.  [CDS 8/95]
**
** For a full CI, make all strings with a given irrep belong to the same
**    graph.  [CDS 6/95]
**
** New, faster, simpler version based on previous OlsenGraph version.
**    [CDS 5/96]
**
** Parameters:
**    GraphSet     =  struct graph_set to hold all subgraphs
**    ci_orbs      =  number of orbitals explicitly in CI 
**    num_el       =  number of electrons for the relevant string
**    nirreps      =  number of irreducible representations
**    orbsym       =  orbital symmetry array
**    ras1_lvl     =  last level in RAS I
**    ras1_min     =  min number of electrons at RAS I level for the string 
**                    (nb this definition includes core electrons!)
**    ras1_max     =  maximum number of electrons at RAS I level 
**    ras3_lvl     =  first level in RAS III
**    ras3_max     =  max number of electrons in RAS III _for the string_
**    num_fzc_orbs = number of frozen core orbitals
**    num_cor_orbs = number of restricted core orbitals
**    ras4_lvl     =  first level of the new RAS IV 
**    ras4_max     =  max number of electrons in RAS IV for the string
**    ras34_max    =  max number of electrons in RAS III and IV
**
** Returns: none
*/
void graphset(struct graph_set *GraphSet, int ci_orbs, int num_el, 
      int nirreps, int *orbsym, int ras1_lvl, int ras1_min, int ras1_max, 
      int ras3_lvl, int ras3_max, int num_fzc_orbs, int num_cor_orbs,
      int ras4_lvl, int ras4_max, int ras34_max)
{
   Odometer Ras1, Ras2, Ras3, Ras4;
   int n1, n2, n3, n4;
   int n1max, n1min;
   int max_el_ras1;
   int *occs, *array1, *array2, *array3, *array4, **encode_tmp;
   int i, j, ij, k, l;
   int orbs_frozen, fzc_sym=0, code=0, num_el_expl;


   #ifdef DEBUG
   fprintf(outfile, "ras1_lvl = %d   ras1_min = %d  ras1_max = %d\n",
      ras1_lvl, ras1_min, ras1_max) ;
   fprintf(outfile, "ras3_lvl = %d   ras3_max = %d\n", ras3_lvl, ras3_max) ;
   #endif

   // Go ahead and set the occupations of the frozen orbs 
   occs = init_int_array(num_el) ;
   for (i=0; i<num_cor_orbs; i++) occs[i] = i ;
 
   orbs_frozen = num_fzc_orbs + num_cor_orbs;

   for (i=0; i<num_fzc_orbs; i++) fzc_sym ^= orbsym[i];

   // go ahead and make room for the occupations of RAS I, II, and III
   array1 = init_int_array(num_el) ;
   array2 = init_int_array(num_el) ;
   array3 = init_int_array(num_el) ;
   array4 = init_int_array(num_el) ;

   // Initialize the Graph data structure
   GraphSet->num_el = num_el;
   num_el_expl = num_el - num_fzc_orbs;
   GraphSet->num_el_expl = num_el_expl;
   GraphSet->num_orb = ci_orbs ;
   GraphSet->num_fzc_orbs = num_fzc_orbs;
   GraphSet->num_cor_orbs = num_cor_orbs;
   GraphSet->fzc_sym = fzc_sym;
   GraphSet->orbsym = init_int_array(ci_orbs);
   for (i=0; i<ci_orbs; i++) {
      GraphSet->orbsym[i] = orbsym[i+num_fzc_orbs];
      }
   GraphSet->ras1_lvl = ras1_lvl;
   GraphSet->ras1_min = ras1_min;
   GraphSet->ras1_max = ras1_max;
   GraphSet->ras3_lvl = ras3_lvl; 
   GraphSet->ras3_max = ras3_max;
   GraphSet->ras4_lvl = ras4_lvl;
   GraphSet->ras4_max = ras4_max;
   GraphSet->ras34_max = ras34_max;
   GraphSet->nirreps = nirreps;
   GraphSet->str_per_irrep = init_int_array(nirreps);
   n1max = ras1_max - orbs_frozen; 
   n1min = ras1_min - orbs_frozen;

   GraphSet->decode = (int ***) malloc ((ras1_max - ras1_min + 1) *
      sizeof(int **));
   for (i=0; i<(ras1_max - ras1_min + 1); i++) {
      GraphSet->decode[i] = init_int_matrix(ras3_max + 1, ras4_max + 1);
      }

   encode_tmp = init_int_matrix(3, (ras1_max - ras1_min + 1) * (ras3_max + 1)
                   * (ras4_max + 1));

   /* need to know how many possible RAS I/RAS III/RAS IV combinations */
   if (!Parameters.fci_strings) {
      for (i=ras1_max; i>=ras1_min; i--) {
         for (j=0; j<=ras3_max; j++) {
            for (k=0; k<=ras4_max; k++) {
               if ((i+j+k<=num_el) && (num_el-i-j-k<=ras3_lvl-ras1_lvl-1)
                   && (j+k <= ras34_max) && (!(Parameters.r4s && k>=2 &&
                   ras1_max - i > Parameters.ex_lvl))) {
                  GraphSet->decode[i-ras1_min][j][k] = code;
                  encode_tmp[0][code] = i - num_fzc_orbs; 
                  encode_tmp[1][code] = j;
                  encode_tmp[2][code] = k;
                  code++;
                  }
               else GraphSet->decode[i-ras1_min][j][k] = -1;
               }
            }
         }
      }
   else { /* all strings w/ given irrep belong to same graph for FCI */
      for (i=ras1_max; i>=ras1_min; i--) {
         for (j=0; j<=ras3_max; j++) {
            for (k=0; k<=ras4_max; k++) {
               if ((i+j+k<=num_el) && (num_el-i-j-k<=ras3_lvl-ras1_lvl-1)
                   && (j+k <= ras34_max)) {
                  GraphSet->decode[i-ras1_min][j][k] = 0;
                  }
               else GraphSet->decode[i-ras1_min][j][k] = -1;
               }
            }
         }
      code = 1;
      }

   GraphSet->encode = init_int_matrix(3,code);
   for (i=0; i<code; i++) {
      GraphSet->encode[0][i] = encode_tmp[0][i];
      GraphSet->encode[1][i] = encode_tmp[1][i];
      GraphSet->encode[2][i] = encode_tmp[2][i];
      } 
   free_int_matrix(encode_tmp);
 
 
   GraphSet->num_codes = code;
   GraphSet->AllGraph = (struct fastgraph **) malloc (nirreps * code * 
      sizeof(struct fastgraph *));
   for (i=0; i<nirreps * code; i++) {
      GraphSet->AllGraph[i] = (struct fastgraph *) malloc (sizeof(struct
         fastgraph);
      GraphSet->AllGraph[i]->data = NULL;
      GraphSet->AllGraph[i]->num_strings = 0;
      }

   // loop over the possible number of e- in RAS I (n1) and III (n3).
   // and now IV (n4)
   // the number of electrons in RAS II is defined via 
   //    n2 = num_el_this_spin - n1 - n3 - n4
   //
   // Employ the very useful Generalized Odometer
   //

   for (n1 = n1max; n1 >= n1min; n1--) {
      Ras1.resize(n1) ;
      Ras1.set_min_lex(num_cor_orbs) ;
      Ras1.set_max_lex(ras1_lvl) ;

      for (n3 = 0; n3 <= ras3_max; n3++) {

         Ras3.resize(n3) ;
         Ras3.set_min_lex(ras3_lvl) ;
         /* Ras3.set_max_lex(ci_orbs-1) ; */
         Ras3.set_max_lex(ras4_lvl-1);

         for (n4 = 0; n4 <= ras4_max && n4 <= ras34_max - n3; n4++) {

            n2 = num_el - orbs_frozen - n1 - n3 - n4;
            if (n2 < 0 || n2 > ras3_lvl - ras1_lvl - 1) continue ; 

            /* CDS 8/24/95 */
            if (Parameters.r4s && n4 >= 2 && n1max - n1 > Parameters.ex_lvl) 
               continue;

            #ifdef DEBUG
            fprintf(outfile, "n1 = %d, n2 = %d, n3 = %d, n4 = %d\n", 
               n1, n2, n3, n4) ;
            if (n2 < 0) printf("Error: n2 < 0 in form_strings()\n") ;
            #endif

            Ras2.resize(n2) ;
            Ras4.resize(n4) ; 
            Ras2.set_min_lex(ras1_lvl+1) ;
            Ras2.set_max_lex(ras3_lvl-1) ;
            Ras4.set_min_lex(ras4_lvl);
            Ras4.set_max_lex(ci_orbs-1);

            Ras1.reset() ; Ras2.reset() ; Ras3.reset() ;  Ras4.reset() ;

            do {
               Ras1.get_value(array1) ;
               do {
                  Ras2.get_value(array2) ;
                  do {
                     Ras3.get_value(array3) ;
                     do {
                        Ras4.get_value(array4) ;
                        for (i=n1-1, j = num_cor_orbs; i>=0; i--)
                           occs[j++] = array1[i] ;
                        for (i=n2-1; i>=0; i--)
                           occs[j++] = array2[i] ;
                        for (i=n3-1; i>=0; i--)
                           occs[j++] = array3[i] ;
                        for (i=n4-1; i>=0; i--)
                           occs[j++] = array4[i] ;
                        
                        // print out occupations for debugging
                        #ifdef DEBUG
                        for (i=0; i<num_el - num_fzc_orbs; i++) 
                           fprintf(outfile, "%2d ", occs[i]) ;
                        fprintf(outfile, "\n") ;
                        #endif
                  
                        // add this walk to the graph
                        og_add_walk(n1-n1min, n3, n4, occs, num_el_expl,
                           ci_orbs, nirreps, num_fzc_orbs, GraphSet) ;

                        Ras4.increment_lex() ;
                     } while (!Ras4.at_min()) ;
                     Ras3.increment_lex() ;
                     } while (!Ras3.at_min()) ;
                  Ras2.increment_lex() ;
                  } while (!Ras2.at_min()) ;
               Ras1.increment_lex() ;
               } while (!Ras1.at_min()) ;

            } /* end loop over n4 */

         }  /* end loop over n3 */

      } /* end loop over n1 */

   free(array1);
   free(array2);
   free(array3);
   free(array4);
   free(occs);

   /* fill up the graphs from the ki's */
   gs_fill(num_el_expl, ci_orbs, nirreps, num_fzc_orbs, GraphSet);

   /* compact the graphs (go from AllGraph to Graph array) */
   GraphSet->num_str = 0;
   for (i=0,ij=0,k=0; i<nirreps; i++) {
      for (j=0; j<code; j++,ij++) {
         if (l = GraphSet->AllGraph[ij]->num_strings) {
            GraphSet->num_str += l;
            GraphSet->Graph[k] = GraphSet->AllGraph[ij];
            GraphSet->graph_irrep[k] = i;
            GraphSet->graph_code[k] = j;
            k++;
            }
         else free(GraphSet->AllGraph[ij]);
         }
      }
   GraphSet->num_graphs = k;
   /* at this point, AllGraph pointers are no longer valid: use Graph */

   return; 
}



/*
** gs_add_walk(): Add a walk to a subgraph within the Olsen/Roos scheme
**    of subgraphs.  Uses struct graph_set.
**
** Parameters:
**    ras1_idx     = Number of e- in RAS I - minimum # of e- in RAS I
**                   for the given _string_
**    ras3_num     = Number of electrons in RAS III
**    ras4_num     = Number of electrons in RAS IV
**    occs         = array listing orbital each electron occupies 
**    nel_expl     = number of electrons explicitly treated (i.e. minus
**                   all implicitly treated frozen core electrons)
**    norb         = number of orbitals _explicitly_ included 
**    nirreps      = number of irreps
**    num_fzc_orbs = number of frozen core orbitals
**    Graph        = Olsen Graph structure containing all subgraphs for a 
**                   given electron spin (alpha or beta)
**
** Returns: none
**
** Note: The newidx code excludes frozen core electrons from the num_el
**       factor
*/
void gs_add_walk(int ras1_idx, int ras3_num, int ras4_num, int *occs, 
      int nel_expl, int norb, int nirreps, int num_fzc_orbs,
      struct olsen_graph *GraphSet)
{
   int i;
   int irrep, ncodes;
   int cur_el=0, orb=0;
   struct fastgraph *graph;
   int code;
   int *orbsym;

   orbsym = GraphSet->orbsym;
   irrep = GraphSet->fzc_sym;

   /* figure out the irrep for this walk */
   for (i=0; i<nel_expl; i++) {
      orb = occs[i] ;
      irrep ^= orbsym[orb] ;
      }

   /* get a pointer to the right subgraph */
   code = GraphSet->decode[ras1_idx][ras3_num][ras4_num];
   ncodes = GraphSet->num_codes;
   graph = GraphSet->AllGraph[irrep * ncodes + code];

   if (graph == NULL) {
      printf("Error (gs_add_walk): NULL subgraph pointer\n");
      return;
      }
   if (code < 0) {
      printf("Error (gs_add_walk): negative RAS code\n");
      return;
      }

   /* make sure that graph's data has been malloc'd */
   if (graph->data == NULL) {
      graph->data = init_int_matrix(nel_expl+1, norb+1);
      }

   /* loop over all (explicitly included) orbitals */
   for (i=0; i<norb; i++) {

      if (occs[cur_el] == i) { /* if the current orbital is occupied */
         graph->data[cur_el][i] |= OCC_BIT;
         cur_el++ ;
         }
      
      else {
         graph->data[cur_el][i] |= UOC_BIT;
         }

      } /* end loop over i */

   graph->num_strings++;

}



/*
** gs_fill(): Fill out the Olsen-Roos subgraph DRT's.  So far, all we have
**    is a list of preliminary links.
**
** Parameters:
**    nel      = total of explicit electrons for string
**    norb     = number of CI orbitals
**    nirreps  = number of irreps
**    Graph    = GraphSet structure containing all subgraphs for a given
**               electron spin (alpha or beta)
**
** Returns: none
**
** David Sherrill, May 1996
*/
void gs_fill(int nel, int norb, int nirreps, struct graph_set *GraphSet)
{
   int gnum, orb, el;
   int **xmat;

   /* allocate scratch matrices */
   xmat = init_int_matrix(nel+1, norb+1);

   /* loop over graphs */
   for (gnum=0; gnum<GraphSet->num_codes * nirreps; gnum++) {

      graph = GraphSet->AllGraph[gnum];
      if (graph==NULL) {
         throw PsiException("(gs_fill): Error, get NULL graph pointer!",__FILE__,__LINE__);
         }
      
      /* calculate vertex weights x first */
      zero_mat(xmat, nel+1, norb+1);
      xmat[0][0] = 1;

      for (orb=1; orb<=norb; orb++) {
         for (el=0; el<=nel; el++) {
            if (graph->data[el][orb-1] & UOC_BIT) 
               xmat[el][orb] += xmat[el][orb-1];
            if (el > 0 && graph->data[el-1][orb-1] & OCC_BIT)
               xmat[el][orb] += xmat[el-1][orb-1];
            }
         }

      /* now we have x, so calculate y */
      zero_mat(graph->data, nel+1, norb+1);
      for (orb=0; orb<norb; orb++) {
         for (el=0; el<nel; el++) {
            graph->data[el][orb] = xmat[el+1][orb];
            }
         }

      /* check the value of num_strings */
      if (graph->num_strings != xmat[nel][norb]) {
         std::string str = "(gs_fill): num_strings != x[nel][norb] for graph ";
         str += str += static_cast<std::ostringstream*>( &(std::ostringstream() << gnum) )->str();
         throw PsiException(str,__FILE__,__LINE__);
         }

      } /* end loop over graphs */

   /* free scratch matrices */
   free_matrix(xmat, nel+1);

}


void gs_print(struct graphset *GraphSet, FILE *outfile)
{

   int ras1_min, ras1_max, ras3_max, ras4_max, code;
   int i,j,k;
   struct fastgraph *graph;

   ras1_min = GraphSet->ras1_min;
   ras1_max = GraphSet->ras1_max;
   ras3_max = GraphSet->ras3_max;
   ras4_max = GraphSet->ras4_max;

   fprintf(outfile,"\nGraphSet:\n");
   fprintf(outfile,"%3c%2d Electrons\n",' ',GraphSet->num_el);
   fprintf(outfile,"%3c%2d Frozen core orbitals\n",' ',GraphSet->num_fzc_orbs);
   fprintf(outfile,"%3c%2d Restricted core orbs\n",' ',GraphSet->num_cor_orbs);
   fprintf(outfile,"%3c%2d Explicit electrons\n",' ',GraphSet->num_el_expl);
   fprintf(outfile,"%3c%2d Explicit Orbitals\n",' ',GraphSet->num_orb);
   fprintf(outfile,"%3c%2d RAS I level\n",' ',GraphSet->ras1_lvl);
   fprintf(outfile,"%3c%2d RAS I minimum\n",' ',ras1_min);
   fprintf(outfile,"%3c%2d RAS I maximum\n",' ',ras1_max);
   fprintf(outfile,"%3c%2d RAS III level\n",' ',GraphSet->ras3_lvl);
   fprintf(outfile,"%3c%2d RAS III maximum\n",' ',ras3_max);
   fprintf(outfile,"%3c%2d RAS IV maximum\n",' ',ras4_max);
   fprintf(outfile,"%3c%2d Number of irreps\n",' ',GraphSet->nirreps);
   fprintf(outfile,"%3c%2d Number of codes\n",' ',
      Graph->num_codes);
   fprintf(outfile,"%3c%2d Max strings in irrep\n", ' ', 
      Graph->max_str_per_irrep);
   fprintf(outfile,"%3c%2d Strings in total\n\n", ' ', GraphSet->num_str);

   fprintf(outfile, "\n");
   for (i=ras1_min; i<=ras1_max; i++) {
      for (j=0; j<=ras3_max; j++) {
         for (k=0; k<=ras4_max; k++) {
            if ((code = GraphSet->decode[i-ras1_min][j][k]) >= 0) {
               fprintf(outfile, "%5cDecode (%2d,%2d,%2d) = %3d\n",' ',
                  i,j,k,code);
               } 
            }
         }
      }

   fprintf(outfile, "\n%4cString Digraphs\n", ' ');
   for (i=0; i<GraphSet->num_graphs; i++) {
      graph = GraphSet->Graph[i];
      fprintf(outfile, "%6cGraph %3d (Code=%2d,Irrep=%1d): %4d strings, 
         offset = %4d\n", ' ', i, GraphSet->graph_code[i], 
         GraphSet->graph_irrep[i], graph->num_strings, 
         GraphSet->graph_offset[i]);

      print_int_mat(outfile, graph->data, GraphSet->num_el_expl+1, 
         GraphSet->num_orb+1);

      fprintf(outfile, "\n");
      }       

   fprintf(outfile, "\n");
   fflush(outfile);
}


void gs_stringlist(struct graph_set *GraphSet, struct stringwr **slist)
{
   Odometer Ras1, Ras2, Ras3, Ras4;
   int n1, n2, n3, n4;
   int n1max, n1min, orbs_frozen, ci_orbs;
   int *occs, *array1, *array2, *array3, *array4;
   int ras1_lvl, ras3_lvl, ras4_lvl, ras3_max, ras4_max, ras34_max;
   int i, num_el_expl, irrep, code, ncodes, gnum, snum;

   /* Regenerate all strings, store occs, and get replacement info */
   num_el_expl = GraphSet->num_el_expl;
   GraphSet->Occs = (unsigned char ***) malloc (GraphSet->num_graphs *
      sizeof(unsigned char **));
   for (i=0; i<GraphSet->num_graphs; i++) {
      GraphSet->Occs[i] = (unsigned char **) malloc (sizeof(unsigned char *)
         * GraphSet->Graph[i]->num_strings);
      for (j=0; j<GraphSet->Graph[i]->num_strings; j++) {
         GraphSet->Occs[i][j] = (unsigned char *) malloc (num_el_expl *
            sizeof(unsigned char));
         }
      }

   ncodes = GraphSet->num_codes;
   occs = init_int_array(num_el_expl) ;
   for (i=0; i<num_cor_orbs; i++) occs[i] = i;

   orbs_frozen = num_fzc_orbs + num_cor_orbs;
   ci_orbs = GraphSet->num_orb;

   array1 = init_int_array(num_el);
   array2 = init_int_array(num_el);
   array3 = init_int_array(num_el);
   array4 = init_int_array(num_el);

   n1max = GraphSet->ras1_max - orbs_frozen; 
   n1min = GraphSet->ras1_min - orbs_frozen;
   ras1_lvl = GraphSet->ras1_lvl;
   ras3_lvl = GraphSet->ras3_lvl;
   ras4_lvl = GraphSet->ras4_lvl;
   ras3_max = GraphSet->ras3_max;
   ras4_max = GraphSet->ras4_max;
   ras34_max = GraphSet->ras34_max;

   // loop over the possible number of e- in RAS I (n1) and III (n3).
   // and now IV (n4)
   // the number of electrons in RAS II is defined via 
   //    n2 = num_el_this_spin - n1 - n3 - n4
   //
   // Employ the very useful Generalized Odometer
   //

   for (n1 = n1max; n1 >= n1min; n1--) {
      Ras1.resize(n1) ;
      Ras1.set_min_lex(num_cor_orbs) ;
      Ras1.set_max_lex(ras1_lvl) ;

      for (n3 = 0; n3 <= ras3_max; n3++) {

         Ras3.resize(n3) ;
         Ras3.set_min_lex(ras3_lvl) ;
         /* Ras3.set_max_lex(ci_orbs-1) ; */
         Ras3.set_max_lex(ras4_lvl-1);

         for (n4 = 0; n4 <= ras4_max && n4 <= ras34_max - n3; n4++) {

            n2 = num_el - orbs_frozen - n1 - n3 - n4;
            if (n2 < 0 || n2 > ras3_lvl - ras1_lvl - 1) continue ; 

            /* CDS 8/24/95 */
            if (Parameters.r4s && n4 >= 2 && n1max - n1 > Parameters.ex_lvl) 
               continue;

            #ifdef DEBUG
            fprintf(outfile, "n1 = %d, n2 = %d, n3 = %d, n4 = %d\n", 
               n1, n2, n3, n4) ;
            if (n2 < 0) printf("Error: n2 < 0 in form_strings()\n") ;
            #endif

            Ras2.resize(n2) ;
            Ras4.resize(n4) ; 
            Ras2.set_min_lex(ras1_lvl+1) ;
            Ras2.set_max_lex(ras3_lvl-1) ;
            Ras4.set_min_lex(ras4_lvl);
            Ras4.set_max_lex(ci_orbs-1);

            Ras1.reset() ; Ras2.reset() ; Ras3.reset() ;  Ras4.reset() ;

            do {
               Ras1.get_value(array1) ;
               do {
                  Ras2.get_value(array2) ;
                  do {
                     Ras3.get_value(array3) ;
                     do {
                        Ras4.get_value(array4) ;
                        for (i=n1-1, j = num_cor_orbs; i>=0; i--)
                           occs[j++] = array1[i] ;
                        for (i=n2-1; i>=0; i--)
                           occs[j++] = array2[i] ;
                        for (i=n3-1; i>=0; i--)
                           occs[j++] = array3[i] ;
                        for (i=n4-1; i>=0; i--)
                           occs[j++] = array4[i] ;
                        
                        // print out occupations for debugging
                        #ifdef DEBUG
                        for (i=0; i<num_el - num_fzc_orbs; i++) 
                           fprintf(outfile, "%2d ", occs[i]) ;
                        fprintf(outfile, "\n") ;
                        #endif
                  
                        // add this walk to the Occs array

                        irrep = GraphSet->fzc_sym;
                        for (i=0; i<num_el_expl; i++) irrep ^= orbsym[i]; 
                        code = GraphSet->decode[n1-n1min][n3][n4];
                        gnum = GraphSet->AllGraph2Graph[irrep * ncodes + code];
                        snum = gs_glex_addr(GraphSet->Graph[gnum], occs, 
                           num_el_expl);

                        for (i=0; i<num_el_expl; i++) 
                           GraphSet.Occs[gnum][snum][i] = (unsigned char) 
                              occs[i];

                        gs_form_stringwr(slist[gnum], occs, num_el_expl,
                           GraphSet->num_orb, GraphSet->Graph[gnum],
                           GraphSet, GraphSet->num_cor_orbs);

                        Ras4.increment_lex() ;
                     } while (!Ras4.at_min()) ;
                     Ras3.increment_lex() ;
                     } while (!Ras3.at_min()) ;
                  Ras2.increment_lex() ;
                  } while (!Ras2.at_min()) ;
               Ras1.increment_lex() ;
               } while (!Ras1.at_min()) ;

            } /* end loop over n4 */

         }  /* end loop over n3 */

      } /* end loop over n1 */

   free(array1);
   free(array2);
   free(array3);
   free(array4);
   free(occs);


}


int gs_glex_addr(struct fastgraph *graph, int *occs, int nel)
{
   int addr = 0;
   int *yptr;

   yptr = graph->data;

   for (i=0; i<nel; i++, yptr++) {
      j = occs[i];
      addr += yptr[j];
      }
}


/*
** gs_form_stringwr(): Make the string with replacements list.
**    This version uses the graph_set structures. 
**
** Parameters:
**   strlist          = list of strings (may be different for alpha and beta)
**   occs             = list of occupied orbitals (N long)
**   N                = number of electrons explicitly included
**   num_ci_orbs      = number of active CI orbitals
**   graph            = the graph including the given walk
**   GraphSet         = the graph_set for the alpha or beta strings
**   first_orb_active = first alp/bet orb active
*/
void form_stringwr(struct stringwr *strlist, int *occs, int N,
      int num_ci_orbs, int gnum, struct fastgraph *graph, struct graph_set
      *GraphSet, int first_orb_active)
{
   unsigned char *occlist;
   unsigned int addr;
   int i;
   struct stringwr *node;
   
   occlist = (unsigned char *) malloc (N * sizeof(unsigned char));
   if (occlist == NULL) {
      throw PsiException("(gs_form_stringwr): Malloc error",__FILE__,__LINE__);
      }

   for (i=0; i<N; i++) {
      occlist[i] = occs[i] ;
      }

   addr = gs_glex_addr(graph, occs, N);
 
   node = strlist + addr;
   node->occs = occlist;

   if (!Parameters.repl_otf) {
      gs_init_repinfo_temps(GraphSet->num_el_expl, GraphSet->num_orbs);
      gs_form_repinfo(node, GraphSet, occs, gnum, first_orb_active);
      gs_free_repinfo_temps();
      }
}


void gs_form_repinfo(struct stringwr *string, struct graph_set *GraphSet, 
      int *occs, int sgraph, int first_orb_active)
{
   int i, j, ngraphs, nel, nras, cnt, tgraph;
   int *I_n[N_RAS_SPACES];
   int *J_n[N_RAS_SPACES];
   int *D_n[N_RAS_SPACES];
   int Isym, Jsym;

   ngraphs = GraphSet->num_graphs;
   nel = GraphSet->num_el_expl;
   nras = GraphSet->num_ras_spaces;
   Isym = GraphSet->graph_irrep[sgraph];
   Icode = GraphSet->graph_code[sgraph];

   for (i=0; i<nras; i++) {
      I_n[i] = GraphSet->encode2[i][Icode];
      if (I_n[i] < 0) {
         printf("(gs_form_repinfo): Got less than 0 e- in a partition!\n");
         return;
         }
      }

   string->cnt = init_int_array(ngraphs);
   string->ij = (int **) malloc (sizeof(int *) * ngraphs);
   string->oij = (int **) malloc(sizeof(int *) * ngraphs);
   string->ridx = (unsigned int **) malloc(sizeof(unsigned int *) * ngraphs);
   string->sgn = (signed char **) malloc(sizeof(signed char *) * ngraphs);

   for (tgraph=0; tgraph<ngraphs; tgraph++) {
      Jcode = GraphSet->graph_code[tgraph];

      for (i=0; i<nras; i++) {
         J_n[i] = GraphSet->encode2[i][Jcode];
         if (J_n[i] < 0) {
            printf("(gs_form_repinfo): Got less than 0 e- in a partition!\n");
            return;
            }
         D_n[i] = J_n[i] - I_n[i];
         }

      /* are these ok? */
      if ((i = abs(D_n1) + abs(D_n2) + abs(D_n3) + abs(D_n4)) > 2) {
         string->cnt[tgraph] = 0;
         continue;
         }

      Jsym = GraphSet->graph_irrep[tgraph];
      ijsym = Isym ^ Jsym;

      if (i==0) {
         cnt = s2bgen1(GraphSet->Graph[tgraph], occs, I_n, ijsym, nel, 
            GraphSet->num_orb, nras, GraphSet->raslevels, 
            GraphSet->orbsym);
         }
      else {
         cnt = s2bgen2(occs, ijsym, nel, nras, GraphSet->raslevels, tgraph);
         }

      string->cnt[tgraph] = cnt;
      if (cnt) {
         string->ij[i] = init_int_array(cnt);
         string->oij[i] = init_int_array(cnt);
         string->ridx[i] = (unsigned int *) malloc(cnt * sizeof(unsigned int));
         string->sgn[i] = (signed char *) malloc(cnt * sizeof(signed char));
         /* could sort by ij, only takes a few ll (see stringlist.c l.369) */
         for (i=0; i<cnt; i++) {
            string->ij[i] = Jij[i];
            string->oij[i] = Joij[i];
            string->ridx[i] = Jridx[i];
            string->sgn[i] = Jsgn[i];
            }
         }
      else {
         string->ij[tgraph] = NULL; 
         string->oij[tgraph] = NULL; 
         string->ridx[tgraph] = NULL; 
         string->sgn[tgraph] = NULL; 
         }

      } /* end loop over tgraph */

}


int s2bgen1(struct fastgraph *graph, int *occs, int *I_n, int ijsym, int nel, 
      int norb, int nras, int *raslevels, int *orbsym)
{
   int i, j, ij, oij, k, orb, ridx;
   int tloc, tloc2, iused, hole, abshole, part;
   int T[MAX_EL], ecnt[N_RAS_SPACES], ras_occs[N_RAS_SPACES][MAX_EL];
   int cnt = 0;  /* how many singlerepls found for this string */

   for (i=0; i<nel; i++) {
      orb = occs[i];
      for (j=0; j<nras; j++) {
         if (orb < raslelvels[j]) {
            k = ecnt[j];
            ecnt[j] += 1;
            ras_occs[j][k] = orb;
            break; 
            }
         }
      }

   #ifdef DEBUG
   for (i=0; i<nras; i++) {
      if (I_n[i] != ecnt[i]) {
         throw PsiException("(s2bgen1): I_n != ecnt",__FILE__,__LINE__);
         }
      }
   #endif DEBUG
   
   /* do diagonals first */
   if (ijsym == 0) {
      ridx = gs_glex_addr(graph, occs, nel);
      for (k=0; k<nel; k++) {
         i = occs[k];
         ij = ioff[i] + i;
         oij = i * norb + i;
         Tij[cnt] = ij;
         Toij[cnt] = oij;
         Tsgn[cnt] = 1;
         Tridx[cnt] = ridx; 
         cnt++;
         }
      }

       
   /* Done with diagonals.  Now loop over RAS subspaces */
   /* loop over excited electrons/orbitals */

   for (ras=0,abshole=0; ras<nras; ras++) {
   
      /* arrange the occupied list in order */
      for (i=0,tloc=0; i<ras; i++) { 
         for (j=0; j<ecnt[i]; j++) {
            T[tloc++] = ras_occs[i][j];
            }
         }
      k = tloc + ecnt[ras];
      for (i=ras+1; i<nras; i++) {
         for (j=0; j<ecnt[i]; j++) {
            T[k++] = ras_occs[i][j];
            }
         }
      
      for (hole=0; hole<ecnt[ras]; hole++,abshole++) { 
 
         j = ras_occs[ras][hole];
         if (j < CalcInfo.num_cor_orbs) continue;
         jsym = orbsym[j];
         isym = ijsym ^ jsym;

         for (part=0,tloc2=tloc; part<ras_opi[ras][isym]; part++) {
            i = ras_orbs[ras][isym][part]; 
            
            /* make sure i is not occupied already */
            for (k=0,l=0; k<ecnt[ras]; k++) {  
               if (i == ras_occs[ras][k]) {l=1; break;}
               }
            if (l==1) continue;

            /* continue to arrange the occupied list in order */
            for (k=0,iused=0; k<ecnt[ras]; k++) {
               if (!iused && i < ras_occs[ras][k]) {
                  iused=1; 
                  hops = tloc2; 
                  T[tloc2++] = i;
                  }
               if (k != hole) {
                  T[tloc2++] = ras_occs[ras][k];
                  }
               }
            if (!iused) {
               hops = tloc2;
               T[tloc2++] = i;
               }


            ridx = gs_glex_addr(graph, T, nel);

            /* get the sign */
            sgn = ((abshole + hops) % 2) ? -1 : 1;
            ij = INDEX(i,j);
            oij = i * norb + j;
            /* store these results and increment the counter */
            Tij[cnt] = ij;
            Toij[cnt] = oij;
            Tsgn[cnt] = sgn;
            Tridx[cnt] = ridx;
            cnt++;

            } /* end loop over particles */
         } /* end loop over holes */
      } /* end loop over RAS subspaces */

   return(cnt);
}


int s2bgen2(struct fastgraph *graph, int *occs, int *I_n, int ijsym, int nel, 
      int norb, int nras, int *raslevels, int *orbsym)
{
   int i, j, ij, oij, k, orb, ridx;
   int hole, abshole, part;
   int T[MAX_EL], ecnt[N_RAS_SPACES], ras_occs[N_RAS_SPACES][MAX_EL];
   int cnt = 0;  /* how many singlerepls found for this string */
   int **ras_orbs, *ras_opi, *ras_occs_excite, *ras_occs_virt;

   for (i=0; i<nel; i++) {
      orb = occs[i];
      for (j=0; j<nras; j++) {
         if (orb < raslelvels[j]) {
            k = ecnt[j];
            ecnt[j] += 1;
            ras_occs[j][k] = orb;
            break; 
            }
         }
      }

   #ifdef DEBUG
   for (i=0; i<nras; i++) {
      if (I_n[i] != ecnt[i]) {
         throw PsiException("(s2bgen1): I_n != ecnt",__FILE__,__LINE__);
         }
      }
   #endif DEBUG

   ras_occs_excite = ras_occs[down]; 
   ras_occs_virt = ras_occs[up];
   ras_opi = CalcInfo.ras_opi[up];
   ras_orbs = CalcInfo.ras_orbs[up];

   abshole = 0;
   for (i=0; i<down; i++) abshole += ecnt[k];

   /* arrange the part of T we already know about (if any) */
   for (i=0,tloc=0; i<nras; i++) {
      if (i == down) { tlocd = tloc;  tloc += ecnt[i] - 1; }
      else if (i == up) { tlocu = tloc;  tloc += ecnt[i] + 1; }
      else {
         for (j=0; j<ecnt[i]; j++) {
            T[tloc++] = ras_occs[i][j];
            } 
         }
      }


   for (hole=0; hole<ecnt[down]; hole++,abshole++) {

      j = ras_occs_excite[hole];
      if (j < CalcInfo.num_cor_orbs) continue;
      jsym = orbsym[j];
      isym = ijsym ^ jsym;
      for (part=0; part<ras_opi[isym]; part++) {
         i = ras_orbs[isym][part];

         /* make sure i is not occupied already */
         for (k=0,l=0; k<ecnt[up]; k++) {
            if (i == ras_occs_virt[k]) {l=1; break;}
            }
         if (l==1) continue;

         /* arrange the occupied list in order */
         for (k=0,l=tlocd; k<ecnt[down]; k++) {
            if (k != hole) T[l++] = ras_occs_excite[k];
            }

         for (k=0,iused=0,l=tlocu; k<ecnt[up]; k++) {
            if (!iused && i < ras_occs_virt[k]) { iused=1; hops=l; T[l++]=i; }
            T[l++] = ras_occs_virt[k];
            }
         if (!iused) { hops = l; T[l++] = i; }

         ridx = gs_glex_addr(graph, T, nel);

         /* get the sign */
         sgn = ((abshole + hops) % 2) ? -1 : 1;

         ij = INDEX(i,j);
         oij = i * norb + j;

         /* store these results and increment the counter */
         Tij[cnt] = ij;
         Toij[cnt] = oij;
         Tsgn[cnt] = sgn;
         Tridx[cnt] = ridx;
         cnt++;

         } /* end loop over particles */
      } /* end loop over holes */

   return(cnt);

}


void gs_init_repinfo_temps(int nel, int norbs)
{
   int maxcnt;

   maxcnt = nel * norbs;

   Tij = init_int_array(maxcnt);
   Toij = init_int_array(maxcnt);
   Tridx = (unsigned int *) malloc(sizeof(unsigned int) * maxcnt);
   Tsgn = (signed char *) malloc(sizeof(signed char) * maxcnt); 
}


void gs_free_repinfo_temps(void)
{
   free(Tij);
   free(Toij);
   free(Tridx);
   free(Tsgn);
}

}} // namespace psi::detci

