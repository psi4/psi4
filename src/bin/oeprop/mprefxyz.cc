/*! \file
    \ingroup OEPROP
    \brief Enter brief description of file here 
*/
#define EXTERN
#include "includes.h"
#include "prototypes.h"
#include "globals.h"
#include <masses.h>

namespace psi { namespace oeprop {

void compute_mp_ref_xyz()
{
  int i;
  double sum,sumx,sumy,sumz;
  double qnty;
  
  sum = sumx = sumy = sumz = 0.0;
  switch (mp_ref) {
  
    case 1: 	/* The center of mass */
      for(i=0;i<natom;i++) {
        qnty = an2masses[(int)zvals[i]];
        sum += qnty;
        sumx += qnty * geom[i][0];
        sumy += qnty * geom[i][1];
        sumz += qnty * geom[i][2];
      }
      break;

    case 2:	/* The origin */
      sum = 1.0;
      break;
    
    case 3:	/* The center of Mulliken's electronic charge */
      for(i=0;i<natom;i++) {
        qnty = zvals[i] - qnet[i];
        sum += qnty;
        sumx += qnty * geom[i][0];               
        sumy += qnty * geom[i][1];               
        sumz += qnty * geom[i][2];
      }
      break;
    
    case 4:	/* The center of nuclear charge */
      for(i=0;i<natom;i++) {
        qnty = zvals[i];
        sum += qnty;
        sumx += qnty * geom[i][0];
        sumy += qnty * geom[i][1];
        sumz += qnty * geom[i][2];
      }
      break;
    
    case 5:	/* The center of net charge */
      for(i=0;i<natom;i++) {
        qnty = qnet[i];
        sum += qnty;
        sumx += qnty * geom[i][0];
        sumy += qnty * geom[i][1];
        sumz += qnty * geom[i][2];
      }
      break;

  }
  
  mp_ref_xyz[0] = sumx/sum;
  mp_ref_xyz[1] = sumy/sum;
  mp_ref_xyz[2] = sumz/sum;

}

}} // namespace psi::oeprop
