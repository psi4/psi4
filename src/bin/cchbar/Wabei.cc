/*! \file
    \ingroup CCHBAR
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstring>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cchbar {

/* Wabei(): Computes all contributions to the Wabei HBAR matrix
** elements, whose spin-orbital definition is:
**
** Wabei = <ab||ei> - Fme t(mi,ab) + t(i,f) <ab||ef> 
**            (I)         (II)             (IIIa)
**   - P(ab) t(i,f) t(m,b) <am||ef> + 1/2 t(i,f) t(mn,ab) <mn||ef>
**               (IIIb)                            (IIIc)
**   + 1/2 P(ab) t(i,f) t(m,a) t(n,b) <mn||ef>
**                   (IIId)
**   + 1/2 t(mn,ab) <mn||ei> + 1/2 P(ab) t(m,a) t(n,b) <mn||ei>
**             (IVa)                         (IVb)
**   + P(ab) t(mi,fb) <am||ef> - P(ab) t(m,a) <mb||ei>
**               (V)                    (VI)
**   - P(ab) t(m,a) t(ni,fb) <mn||ef>
**                 (VII)
**
**  [cf. Gauss and Stanton, JCP 103, 3561-3577 (1995)]
*/

void Wabei_RHF(void);
void Wabei_ROHF(void);
void WABEI_UHF(void);
void Wabei_UHF(void);
void WAbEi_UHF(void);
void WaBeI_UHF(void);

void Wabei_build(void)
{
  if(params.ref == 0) Wabei_RHF();
  else if(params.ref == 1) Wabei_ROHF();
  else if(params.ref == 2) {
    WABEI_UHF();
    Wabei_UHF();
    WAbEi_UHF();
    WaBeI_UHF();
  }
}

}} // namespace psi::cchbar
