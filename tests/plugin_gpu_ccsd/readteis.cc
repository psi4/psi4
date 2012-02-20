#include<psifiles.h>
#include<libiwl/iwl.h>
#include"gpu_ccsd.h"

using namespace psi;
typedef unsigned long int ULI;
namespace psi{
  void ReadTEIs(double*tei,Options&options);
  void ReadInts(struct iwlbuf *Buf,Options&options,double*tei);
}

namespace psi{

void ReadTEIs(double*tei,Options&options){
  struct iwlbuf Buf; 
  // initialize buffer
  iwl_buf_init(&Buf,PSIF_MO_TEI,0.0,1,1);

  fprintf(outfile,"\n  Begin CC integral sort\n\n");

  //sort
  ReadInts(&Buf,options,tei);

  iwl_buf_close(&Buf,1);
}
/**
  * read two-electron integrals
  */
void ReadInts(struct iwlbuf *Buf,Options&options,double*tei){

  double val;
  ULI lastbuf;
  Label *lblptr;
  Value *valptr;
  ULI nocc,idx, p, q, r, s, pq, rs, pqrs;

  lblptr  = Buf->labels;
  valptr  = Buf->values;
  lastbuf = Buf->lastbuf;

  /**
    * first buffer (read in when Buf was initialized)
    */
  for (idx=4*Buf->idx; Buf->idx<Buf->inbuf; Buf->idx++) {
      p = (ULI) lblptr[idx++];
      q = (ULI) lblptr[idx++];
      r = (ULI) lblptr[idx++];
      s = (ULI) lblptr[idx++];

      pq   = Position(p,q);
      rs   = Position(r,s);
      pqrs = Position(pq,rs);
      tei[pqrs] = (double)valptr[Buf->idx];
  }

  /**
    * now do the same for the rest of the buffers
    */
  while(!lastbuf){
      iwl_buf_fetch(Buf);
      lastbuf = Buf->lastbuf;
      for (idx=4*Buf->idx; Buf->idx<Buf->inbuf; Buf->idx++) {

          p = (ULI) lblptr[idx++];
          q = (ULI) lblptr[idx++];
          r = (ULI) lblptr[idx++];
          s = (ULI) lblptr[idx++];

          pq   = Position(p,q);
          rs   = Position(r,s);
          pqrs = Position(pq,rs);
          tei[pqrs] = (double)valptr[Buf->idx];
      }
  }
}

} // end of namespace
  
