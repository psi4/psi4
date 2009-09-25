/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/

namespace psi { namespace ccdensity {

void V_build(void);
void Gijkl(void);
void Gabcd(void);
void Gibja(void);
void Gijka(void);
void Gciab(void);
void Gijab(void);

/* twopdm(): Computes all contributions to the two-particle density
** matrix for CC-like wave functions.  
**
** Note that the contractions evaluated in the functions below
** actually build the bra-ket symmetrized two-particle density:
**
** Gamma'(pq,rs) = 1/2 [Gamma(pq,rs) + Gamma(rs,pq)],
**
** where Gamma(pq,rs) is the original, non-bra-ket-symmetric
** expression.  This is done to satisfy the 
**
** TDC, July 2002
*/

void twopdm(void)
{
/*  V_build(); */
  Gijkl();
  Gabcd();
  Gijka();
  Gciab();
  Gibja();
  Gijab();
}


}} // namespace psi::ccdensity
