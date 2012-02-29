#ifndef SORT_H
#define SORT_H

typedef unsigned long int ULI;
namespace boost{
template<class T>class shared_ptr;
}
namespace psi{
  struct integral{
    ULI ind;
    double val;
  };
  void Sort(struct iwlbuf *Buf,int nfzc,int nfzv,int norbs,int ndoccact,int nvirt);
  void ijkl_terms(double val,ULI pq,ULI rs,ULI p,ULI q,ULI r,ULI s,ULI o,ULI&nijkl,struct integral*ijkl);
  void ijak_terms(double val,ULI p,ULI q,ULI r,ULI s,ULI o,ULI v,ULI&nijak,struct integral*ijak);
  void ijak2_terms(double val,ULI p,ULI q,ULI r,ULI s,ULI o,ULI v,ULI&nijak2,struct integral*ijak2);
  void klcd_terms(double val,ULI pq,ULI rs,ULI p,ULI q,ULI r,ULI s,ULI o,ULI v,ULI&nklcd,struct integral*klcd);
  void akjc_terms(double val,ULI p,ULI q,ULI r,ULI s,ULI o,ULI v,ULI&nklcd,struct integral*klcd);
  void abci3_terms(double val,ULI p,ULI q,ULI r,ULI s,ULI o,ULI v,ULI&nabci3,struct integral*abci3);
  void abci5_terms(double val,ULI p,ULI q,ULI r,ULI s,ULI o,ULI v,ULI&nabci5,struct integral*abci5);
  void abcd1_terms(double val,ULI pq,ULI rs,ULI p,ULI q,ULI r,ULI s,ULI o,ULI v,ULI&nabcd1,struct integral*abcd1);
  void abcd2_terms(double val,ULI pq,ULI rs,ULI p,ULI q,ULI r,ULI s,ULI o,ULI v,ULI&nabcd2,struct integral*abcd2);

  void SortBlock(ULI nelem,ULI blockdim,struct integral*buffer,double*tmp,ULI PSIFILE,char*string,ULI maxdim);
};


#endif
