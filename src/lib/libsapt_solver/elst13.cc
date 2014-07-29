/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#include "sapt2p3.h"

namespace psi { namespace sapt {

void SAPT2p3::elst13()
{
  double e_elst130 = elst130(wBAA_,wBRR_,CHFA_,PSIF_SAPT_AMPS,
    "qAA Density Matrix","qRR Density Matrix","Y3 AR Amplitudes",
    foccA_,noccA_,nvirA_);

  if (debug_) {
    psi::fprintf(outfile,"    Elst130,r           = %18.12lf H\n",e_elst130);
  }

  double e_elst103 = elst130(wABB_,wASS_,CHFB_,PSIF_SAPT_AMPS,
    "qBB Density Matrix","qSS Density Matrix","Y3 BS Amplitudes",
    foccB_,noccB_,nvirB_);

  if (debug_) {
    psi::fprintf(outfile,"    Elst103,r           = %18.12lf H\n\n",e_elst103);
  }

  e_elst13_ = e_elst130 + e_elst103;

  if (print_) {
    psi::fprintf(outfile,"    Elst13,r            = %18.12lf H\n",e_elst13_);
    fflush(outfile);
  }
}

double SAPT2p3::elst130(double **wBAA, double **wBRR, double **CHFA, 
  int ampfile, const char *pAAlabel, const char *pRRlabel, 
  const char *Ylabel, int foccA, int noccA, int nvirA)
{
  int aoccA = noccA - foccA;
  
  double **pAA = block_matrix(aoccA,aoccA);
  psio_->read_entry(ampfile,pAAlabel,(char *) pAA[0],
    sizeof(double)*aoccA*aoccA);

  double **pRR = block_matrix(nvirA,nvirA);
  psio_->read_entry(ampfile,pRRlabel,(char *) pRR[0],
    sizeof(double)*nvirA*nvirA);

  double **yAR = block_matrix(aoccA,nvirA); 
  psio_->read_entry(ampfile,Ylabel,(char *) yAR[0],
    sizeof(double)*aoccA*nvirA);

  double e1 = 0.0, e2 = 0.0, e3 = 0.0;

  for (int a=0; a<aoccA; a++) {
    e1 -= 4.0*C_DDOT(aoccA,pAA[a],1,&(wBAA[a+foccA][foccA]),1);
  }

  e2 += 4.0*C_DDOT(nvirA*nvirA,pRR[0],1,wBRR[0],1);
  e3 += 4.0*C_DDOT(aoccA*nvirA,yAR[0],1,CHFA[foccA],1);

  free_block(pAA);
  free_block(pRR);
  free_block(yAR);

  if (debug_) {
    psi::fprintf(outfile,"\n    Elst13_1            = %18.12lf H\n",e1);
    psi::fprintf(outfile,"    Elst13_2            = %18.12lf H\n",e2);
    psi::fprintf(outfile,"    Elst13_3            = %18.12lf H\n",e3);
  }

  return(e1+e2+e3);
}

}}
