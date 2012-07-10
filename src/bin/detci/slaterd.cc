/*! \file
    \ingroup DETCI
    \brief Enter brief description of file here 
*/

#include <cstdio>
#include <cstdlib> /* was libc.h */
/* gcc 2.7.0 doesn't like #include <cstring> */
#include "slaterd.h"

#include <psi4-dec.h>

namespace psi { namespace detci {

extern double get_twoel(int i, int j, int k, int l);
extern double get_onel(int i, int j);
extern int calc_orb_diff(int cnt, unsigned char *I, unsigned char *J, 
   int *I_alpha_diff, int *J_alpha_diff, int *sign, int *same, 
   int extended);
extern void common_orbs(int *same_alpha, int *same_beta, int cnt_alpha,
   int cnt_beta, int *common_docc, int *common_alpha_socc, 
   int *common_beta_socc, int *cnt_docc, int *cnt_alpha_socc, 
   int *cnt_beta_socc);


void SlaterDeterminant::set(unsigned int na, unsigned char *alpoccs, 
      unsigned int nb, unsigned char *betoccs)
{
   int i;

   if (nalp != na) {
      if (Occs[0] != NULL) free(Occs[0]);
      Occs[0] = (unsigned char *) malloc (sizeof(unsigned char) * na); 
      nalp = na;
      }
   if (nbet != nb) {
      if (Occs[1] != NULL) free(Occs[1]);
      Occs[1] = (unsigned char *) malloc (sizeof(unsigned char) * nb);
      nbet = nb;
      }
   
   for (i=0; i<nalp; i++) {
      Occs[0][i] = alpoccs[i];
      }
   for (i=0; i<nbet; i++) {
      Occs[1][i] = betoccs[i];
      }
}



void SlaterDeterminant::print(void)
{
   print(stdout);
}



void SlaterDeterminant::print(FILE *outfile)
{
   int i;

   fprintf(outfile, "Alpha string: ");
   for (i=0; i<nalp; i++) {
      fprintf(outfile, "%3d ", Occs[0][i]);
      }
   fprintf(outfile, "\n");

   fprintf(outfile, "Beta string : ");
   for (i=0; i<nbet; i++) {
      fprintf(outfile, "%3d ", Occs[1][i]);
      }
   fprintf(outfile, "\n");
}


void SlaterDeterminant::print_config(FILE *outfile)
{
   int i=0, j=0;

   while ((i < nalp) && (j < nbet)) {
      if (Occs[0][i] == Occs[1][j]) {
         fprintf(outfile, "%dX ", Occs[0][i]+1);
         i++; j++;
         }
      else if (Occs[0][i] < Occs[1][j]) {
         fprintf(outfile, "%dA ", Occs[0][i]+1);
         i++;
         } 
      else if (Occs[0][i] > Occs[1][j]) {
         fprintf(outfile, "%dB ", Occs[1][j]+1);
         j++;
         }
      }

   if (i < j) {
      while (i < nalp) {
         fprintf(outfile, "%dA ", Occs[0][i]+1);
         i++;
         }
      }
   else if (i > j) {
      while (j < nbet) {
         fprintf(outfile, "%dB ", Occs[1][j]+1);
         j++;
         }
      }

   fprintf(outfile, "\n") ;

}


SlaterDeterminant& SlaterDeterminant::operator=(const SlaterDeterminant& s)
{
   if (nalp != s.nalp) {
      if (Occs[0] != NULL) free(Occs[0]);
      Occs[0] = (unsigned char *) malloc (sizeof(unsigned char) * s.nalp);   
      }   
   if (nbet != s.nbet) {
      if (Occs[1] != NULL) free(Occs[1]);
      Occs[1] = (unsigned char *) malloc (sizeof(unsigned char) * s.nbet);
      }
   set(s.nalp, s.Occs[0], s.nbet, s.Occs[1]);
   return(*this);
}


int operator ==(SlaterDeterminant& s1, SlaterDeterminant& s2) 
{
   int i;

   if (s1.nalp != s2.nalp || s1.nbet != s2.nbet) return(0);

   for (i=0; i<s1.nalp; i++) {
      if (s1.Occs[0][i] != s2.Occs[0][i]) return(0);
      }
   for (i=0; i<s1.nbet; i++) {
      if (s1.Occs[1][i] != s2.Occs[1][i]) return(0);
      }

   return(1);
}
   

double matrix_element(SlaterDeterminant* I, SlaterDeterminant* J)
{
   static int first_call = 1;
   static int *same_alpha;
   static int *same_beta;
   static int *common_docc;
   static int *common_alp_socc;
   static int *common_bet_socc;
   static int init_nalp, init_nbet;
   static int **I_diff,**J_diff;
   int *common_socc1, *common_socc2;
   int sign=0;
   int nalp, nbet, nalp_same, nbet_same;
   int cnt_docc=0, cnt_alp_socc=0, cnt_bet_socc=0;
   int total_diff=-4, alpha_diff=-2, beta_diff=-2;
   int diffspin, samespin, common1, common2;
   int i, j, k, l, a, b, m;
   double val = 0.0;
  
   if (I->nalp != J->nalp || I->nbet != J->nbet) {
      fprintf(stderr,"(matrix_element): unequal length alp/bet strings!\n");
      return(0.0);
      }

   nalp = I->nalp; 
   nbet = I->nbet;

   if (first_call) {
      same_alpha = (int *) malloc(sizeof(int) * nalp);
      same_beta = (int *) malloc(sizeof(int) * nbet);
      common_alp_socc = (int *) malloc(sizeof(int) * nalp);
      common_bet_socc = (int *) malloc(sizeof(int) * nbet);
      common_docc = (int *) malloc(sizeof(int) * nalp); 
      I_diff = (int **) malloc(sizeof(int *) * 2);
      I_diff[0] = (int *) malloc(sizeof(int) * nalp);
      I_diff[1] = (int *) malloc(sizeof(int) * nalp);
      J_diff = (int **) malloc(sizeof(int *) * 2);
      J_diff[0] = (int *) malloc(sizeof(int) * nalp);
      J_diff[1] = (int *) malloc(sizeof(int) * nalp);
      init_nalp = nalp;
      init_nbet = nbet;
      first_call = 0;
      }

   // make sure that subsequent calls don't change the number of alpha or
   // beta electrons without re-initializing the static arrays
   else {
      if ((nalp != init_nalp) || (nbet != init_nbet)) {
         throw PsiException("(matrix_element): nalp/nbet != init_nalp/nbet",__FILE__,__LINE__);
         }
      
      }

   alpha_diff = calc_orb_diff(nalp, I->Occs[0], J->Occs[0], I_diff[0],
      J_diff[0], &sign, same_alpha, 0);
   beta_diff = calc_orb_diff(nbet, I->Occs[1], J->Occs[1], I_diff[1],
      J_diff[1], &sign, same_beta, 0);

   total_diff = alpha_diff + beta_diff ;
   nalp_same = nalp - alpha_diff ;
   nbet_same = nbet - beta_diff ;

   #ifdef DEBUG
   printf(" alpha_diff = %d\n", alpha_diff) ;
   printf(" beta_diff = %d\n", beta_diff) ;
   printf(" total_diff = %d\n", total_diff) ;
   printf(" sign = %d\n", sign) ;

   for (i=0; i<alpha_diff; i++)
      printf(" I_diff[0][%d] = %d ", i, I_diff[0][i]+1) ;
      printf("\n") ;
   for (i=0; i<alpha_diff; i++)
      printf(" J_diff[0][%d] = %d ", i, J_diff[0][i]+1) ;
      printf("\n\n") ;

   for (i=0; i<beta_diff; i++)
      printf(" I_diff[1][%d] = %d ", i, I_diff[1][i]+1) ;
      printf("\n") ;
   for (i=0; i<beta_diff; i++)
      printf(" J_diff[1][%d] = %d ", i, J_diff[1][i]+1) ;
      printf("\n\n") ;
   for (i=0; i<nalp_same; i++)
      printf(" same_alpha[%d] = %d", i, same_alpha[i]+1) ;
      printf("\n\n") ;
   for (i=0; i<nbet_same; i++)
      printf(" same_beta[%d] = %d", i, same_beta[i]+1) ;
      printf("\n\n") ;
   #endif

   if ((alpha_diff == -2) || (beta_diff == -2)) {
      fprintf(stderr, "(matrix_element): Problem with calc_orb_diff.\n");
      fprintf(stderr, "  Returns -2 value. \n");
      } 

   else if ((alpha_diff == -1) || (beta_diff == -1) || total_diff > 2) {
      return(0.0);
      }

   else if (total_diff == 2) {
      
      if (alpha_diff == 1) {     /* Case 1: 1 in alpha and 1 in beta */

         /* assign i,j,k,l for <ij||kl> */
         i = I_diff[0][0];
         j = I_diff[1][0];
         k = J_diff[0][0];
         l = J_diff[1][0];
         
         #ifdef PRINT_INTS
         if (sign % 2) printf("-"); 
         #endif

         val = get_twoel(i,k,j,l);
         if (sign % 2) val = -val;

         return(val); 
         } /* end Case 1 */ 

      else if ((alpha_diff == 2) || (beta_diff == 2)) { /* Case 2: 2 in alpha */

         if (alpha_diff == 2)
            diffspin = 0; 
         else diffspin = 1;

         /* assign <ij||kl> */
         i = I_diff[diffspin][0];
         j = I_diff[diffspin][1];
         k = J_diff[diffspin][0];
         l = J_diff[diffspin][1]; 

         #ifdef PRINT_INTS
         if (sign % 2) printf("- [ "); 
         #endif

         val = get_twoel(i,k,j,l);
         val -= get_twoel(i,l,j,k);

         #ifdef PRINT_INTS
         if (sign % 2) printf(" ] "); 
         #endif
         
         if (sign % 2) val = -val;
 
         return(val);
         } /* end else if for differ by 2 in alpha or beta */

      else {
         throw PsiException("Error (matrix_element): total_diff != alpha_diff + beta_diff",__FILE__,__LINE__);
         }
 
      } /* end else if for differ by 2 spin orbitals */


   /* Differ by 1 spin orbital */
   else if (total_diff == 1) {
   
      common_orbs(same_alpha, same_beta, nalp_same, nbet_same, common_docc,
         common_alp_socc, common_bet_socc, &cnt_docc, &cnt_alp_socc, 
         &cnt_bet_socc); 
 
      #ifdef DEBUG
      printf("cnt_docc = %d\n", cnt_docc); 
      printf("cnt_alp_socc = %d\n", cnt_alp_socc); 
      printf("cnt_bet_socc = %d\n", cnt_bet_socc); 

      for (i=0; i<cnt_docc; i++)
         printf("common_docc[%d] = %d\n", i, common_docc[i]+1) ; 
         printf("\n") ; 
      for (i=0; i<cnt_alp_socc; i++)
         printf("common_alp_socc[%d] = %d\n", i, common_alp_socc[i]+1) ; 
         printf("\n") ; 
      for (i=0; i<cnt_bet_socc; i++)
         printf("common_bet_socc[%d] = %d\n", i, common_bet_socc[i]+1) ;
         printf("\n") ; 
      #endif

      if (alpha_diff == 1) { /* Case 1: 1 in alpha */   
         diffspin = 0; 
         samespin = 1; 
         common1 = cnt_alp_socc; 
         common2 = cnt_bet_socc; 
         common_socc1 = common_alp_socc; 
         common_socc2 = common_bet_socc; 
         } /* end Case 1: 1 in alpha */

      else if (beta_diff == 1) { /* Case 1: 1 in beta */
         diffspin = 1;
         samespin = 0;
         common1 = cnt_bet_socc; 
         common2 = cnt_alp_socc; 
         common_socc1 = common_bet_socc; 
         common_socc2 = common_alp_socc; 
         } /* end Case 1: 1 in beta */

      else {
         throw PsiException("(matrix_element): impossible case abdiff",__FILE__,__LINE__);
         }

      i = I_diff[diffspin][0];
      j = J_diff[diffspin][0];

      val = get_onel(i,j);

      #ifdef PRINT_INTS
      if (sign % 2) 
         printf(" - [ \n");
      #endif

      for (b=0; b<cnt_docc; b++) { /* looping over all common docc electrons. */
         m = common_docc[b];

         #ifdef PRINT_INTS
         printf("2 * ");
         #endif

         val += 2.0 * get_twoel(i,j,m,m);
        
         #ifdef PRINT_INTS
         printf("-");
         #endif 

         val -= get_twoel(i,m,m,j);
         } /* end loop over all common docc electrons */

      for (b=0; b<common2; b++) { /* looping over all beta or alpha elec.;
                                    beta if differs in one by alpha */
         m = common_socc2[b];
         val += get_twoel(i,j,m,m);
         }          

      for (b=0; b<common1; b++) { /* looping over all alpha or beta elec.;
                                    alpha if differs in one by alpha */
         m = common_socc1[b];

         val += get_twoel(i,j,m,m);
        
         #ifdef PRINT_INTS
         printf("-");
         #endif

         val -= get_twoel(i,m,m,j);
         }

      if (sign % 2) val = -val;

      #ifdef PRINT_INTS
      if (sign % 2) 
         printf(" ] \n");
      #endif

      return(val);

      } /* end if (total_diff == 1) for Case 1: Differ by 1 spin orbital */


   else if (total_diff == 0) { 

      common_orbs(same_alpha, same_beta, nalp_same, nbet_same, common_docc,
         common_alp_socc, common_bet_socc, &cnt_docc, &cnt_alp_socc,
         &cnt_bet_socc); 
     
      val = 0.0;

      #ifdef PRINT_INTS
      if (sign % 2) printf("- [ \n");
      printf(" 2 * [ ");
      #endif

      /* get one-electron integrals */

      for (b=0; b<cnt_docc; b++) {    
         m = common_docc[b];
         val += 2.0 * get_onel(m,m);
         }

      #ifdef PRINT_INTS
      printf(" ] ");
      #endif

      for (b=0; b<cnt_alp_socc; b++) { 
         m = common_alp_socc[b];
         val += get_onel(m,m);
         }

      for (b=0; b<cnt_bet_socc; b++) { 
         m = common_bet_socc[b];
         val += get_onel(m,m);
         }

      #ifdef PRINT_INTS
      printf("\n");
      #endif

      /* get two-electron integrals */

      for (a=0; a<cnt_docc; a++) { 
         i = common_docc[a];

         for (b=0; b<cnt_alp_socc; b++) { 
            j = common_alp_socc[b];
            
            #ifdef PRINT_INTS
            printf("2 * ");
            #endif

            val += 2.0 * get_twoel(i,i,j,j);      

            #ifdef PRINT_INTS
            printf("- ");
            #endif

            val -= get_twoel(i,j,i,j); 
            }

         for (b=0; b<cnt_bet_socc; b++) {
            j = common_bet_socc[b];

            #ifdef PRINT_INTS
            printf("2 * ");
            #endif

            val += 2.0 * get_twoel(i,i,j,j);      

            #ifdef PRINT_INTS
            printf("- ");
            #endif

            val -= get_twoel(i,j,i,j);
            }

         val += get_twoel(i,i,i,i);

         for (b=a+1; b<cnt_docc; b++) {

            j = common_docc[b];

            #ifdef PRINT_INTS
            printf("4 * ");
            #endif

            val += 4.0 * get_twoel(i,i,j,j);

            #ifdef PRINT_INTS
            printf("- 2 * ");
            #endif

            val -= 2.0 * get_twoel(i,j,j,i);

            } /* end loop over pairs of doccs */

         } /* end loop a over doccs */

      /* loop over unique pairs of alpha electrons and alpha with beta */
      for (a=0; a<cnt_alp_socc; a++) {

         i = common_alp_socc[a];

         for (b=a+1; b<cnt_alp_socc; b++) {

            j = common_alp_socc[b];
            
            val += get_twoel(i,i,j,j);

            #ifdef PRINT_INTS
            printf("- ");
            #endif

            val -= get_twoel(i,j,j,i);
            }

         for (b=0; b<cnt_bet_socc; b++) {
            j = common_bet_socc[b];
            val += get_twoel(i,i,j,j);
            }
         }


      /* loop over unique pairs of beta electrons */
      for (a=0; a<cnt_bet_socc; a++) {
         for (b=a+1; b<cnt_bet_socc; b++) {

            i = common_bet_socc[a];
            j = common_bet_socc[b];
            
            val += get_twoel(i,i,j,j);

            #ifdef PRINT_INTS
            printf("- ");
            #endif

            val -= get_twoel(i,j,j,i);
            }
         }

      if (sign % 2) val = -val;

      #ifdef PRINT_INTS
      if (sign % 2) 
         printf(" ] \n");
      #endif

      return(val);

      } /* end total_diff == 0 case */

   else {
      throw PsiException("(matrix_element): Impossible case for total_diff!",__FILE__,__LINE__);
      }

   return(0.0);
}


#ifdef STANDALONE
main() {
   unsigned char string1[2], string2[2], string3[2];

   SlaterDeterminant A;
   SlaterDeterminant B;
   SlaterDeterminant C;
   double value = 0.0 ;

   string1[0] = (unsigned char) 0;
   string1[1] = (unsigned char) 1;
   string2[0] = (unsigned char) 0;
   string2[1] = (unsigned char) 2;
   string3[0] = (unsigned char) 0;
   string3[1] = (unsigned char) 3;

   A.set(2, string1, 2, string1);
   B.set(2, string1, 2, string3);
 
   printf("Slater determinant A \n");
   A.print() ;
   printf("Slater determinant B \n");
   B.print() ;
   printf("Matrix element = %lf\n", matrix_element(&A,&B));
}

#ifdef PRINT_INTS
double get_twoel(int i, int j, int k, int l)
{
   printf("(%d %d | %d %d) ", i, j, k, l);
   return(0.0);
}

double get_onel(int i, int j)
{
   printf("(%d|h|%d) ", i, j);
   return(0.0) ;
}
#endif

#endif

}} // namespace psi::detci

