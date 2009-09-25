/*! \file
    \ingroup DETCI
    \brief Enter brief description of file here 
*/

/* Slater Determinant Class
**
** Based on a previous implementation by David Sherrill using bitsrings 
** from the g++ library and on a symbolic Slater's rules program written
** in C by Matt Leininger
**
** February 7, 1996
**
** C. David Sherrill and Matthew L. Leininger
** Center for Computational Quantum Chemistry
** University of Georgia
** Athens, GA 30606
** 
** Need to #include <cstdio> before this file
**
** Assume number of alpha electrons is greater than or equal to the number
** of beta electrons
**
** Currently matrix_element() uses static temp arrays which are never free'd.
** The arrays should always be small, so this shouldn't be a problem.  The
** malloc'ing is only done once, no matter how many times matrix_element is
** called.
*/


#ifndef _psi_src_bin_detci_slaterd_h
#define _psi_src_bin_detci_slaterd_h

namespace psi { namespace detci {

class SlaterDeterminant {

   protected:
      unsigned nalp;
      unsigned nbet;
      unsigned char *Occs[2];

   public:
      SlaterDeterminant() { nalp=0; nbet=0; Occs[0]=NULL; Occs[1]=NULL; }
      ~SlaterDeterminant() { 
         if (Occs[0] != NULL) free(Occs[0]);
         if (Occs[1] != NULL) free(Occs[1]);
         }
      void set(unsigned int nalp, unsigned char *alpoccs, 
         unsigned int nbet, unsigned char *betoccs);
      void print(void);
      void print(FILE *outfile);
      void print_config(FILE *outfile) ;
      SlaterDeterminant& operator=(const SlaterDeterminant& s) ;
      friend int operator==(SlaterDeterminant& s1, SlaterDeterminant& s2) ;
      friend double matrix_element(SlaterDeterminant* I, SlaterDeterminant* J);
};

}} // namespace psi::detci

#endif // header guard

