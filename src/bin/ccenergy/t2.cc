/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libdpd/dpd.h>
#include "Params.h"
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

void DT2(void), FaetT2(void), FmitT2(void), WmnijT2(void), WmbejT2(void);
void BT2(void), ZT2(void), FT2(void), ET2(void), CT2(void), dijabT2(void);
void BT2_AO(void);
void status(const char *, FILE *);
void FT2_CC2(void);

void t2_build(void)
{
dpdbuf4 tIjAb;
double dotval;

  DT2();
  if(params.print & 2) status("<ij||ab> -> T2", outfile);

  if(params.wfn != "CC2" || params.wfn != "EOM_CC2") { /* skip all this is wfn=CC2 */

    FaetT2();
    FmitT2();
    if(params.print & 2) status("F -> T2", outfile);

    WmnijT2();
    if(params.print & 2) status("Wmnij -> T2", outfile);

#ifdef TIME_CCENERGY
    timer_on("BT2", outfile);
#endif
    if(params.aobasis == "DISK" || params.aobasis == "DIRECT")
      BT2_AO();
    else BT2();
    if(params.print & 2) status("<ab||cd> -> T2", outfile);
#ifdef TIME_CCENERGY
    timer_off("BT2", outfile);
#endif

    ZT2();
    if(params.print & 2) status("Z -> T2", outfile);

#ifdef TIME_CCENERGY
    timer_on("FT2", outfile);
#endif
    FT2();
    if(params.print & 2) status("<ia||bc> -> T2", outfile);
#ifdef TIME_CCENERGY
    timer_off("FT2", outfile);
#endif

    ET2();
    if(params.print & 2) status("<ij||ka> -> T2", outfile);

#ifdef TIME_CCENERGY
    timer_on("WmbejT2", outfile);
#endif
    WmbejT2();
    if(params.print & 2) status("Wmbej -> T2", outfile);
#ifdef TIME_CCENERGY
    timer_off("WmbejT2", outfile);
#endif

#ifdef TIME_CCENERGY
    timer_on("CT2", outfile);
#endif
    CT2();
    if(params.print & 2) status("<ia||jb> -> T2", outfile);
#ifdef TIME_CCENERGY
    timer_off("CT2", outfile);
#endif
  }
  else { /* For CC2, just include (FT2)c->T2 */
    FT2_CC2();
  }

}
}} // namespace psi::ccenergy
