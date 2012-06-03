#include"psi4-dec.h"
#include<psifiles.h>
#include<libiwl/iwl.h>
#include <libpsio/psio.hpp>
#include<libtrans/integraltransform.h>
#include<libtrans/mospace.h>
#include<libmints/wavefunction.h>
#include<libmints/matrix.h>

#ifdef _OPENMP
    #include<omp.h>
#else
    #define omp_get_wtime() 0.0
#endif

#include"coupledpair.h"
#include"sort.h"
#include"blas.h"

using namespace psi;
namespace psi{ namespace cepa{
  long int Position(long int i,long int j);
}}

namespace psi{ namespace cepa{

void SortNoVVVV(struct iwlbuf *Buf,int nfzc,int nfzv,int norbs,int ndoccact,int nvirt,bool islocal,Options&options);

void TransformIntegrals(int nfzc,int nfzv,int norbs,int ndoccact,int nvirt,bool islocal,Options&options){

  struct iwlbuf Buf; 
  // initialize buffer

  fflush(outfile);
  fprintf(outfile,"\n");
  fprintf(outfile, "        *******************************************************\n");
  fprintf(outfile, "        *                                                     *\n");
  fprintf(outfile, "        *        CEPA Integral Transformation and Sort        *\n");
  fprintf(outfile, "        *                                                     *\n");
  fprintf(outfile, "        *******************************************************\n");
  fprintf(outfile,"\n\n");
  fflush(outfile);

  boost::shared_ptr<Wavefunction> wfn =
       Process::environment.reference_wavefunction();

  boost::shared_ptr<PSIO> psio(_default_psio_lib_);
  std::vector<boost::shared_ptr<MOSpace> > spaces;

  spaces.push_back(MOSpace::occ);
  spaces.push_back(MOSpace::all);
  IntegralTransform ints(wfn, spaces, IntegralTransform::Restricted,
           IntegralTransform::IWLOnly, IntegralTransform::QTOrder,
           IntegralTransform::None, false);
  ints.set_keep_dpd_so_ints(1);
  ints.set_keep_iwl_so_ints(1);
  ints.initialize();

  fprintf(outfile,"        ==> Transform (ox|xx) integrals <==\n");
  fprintf(outfile,"\n");
  ints.transform_tei(MOSpace::occ, MOSpace::all, MOSpace::all, MOSpace::all);
  fprintf(outfile,"\n");

  fprintf(outfile,"        ==> Sort (ox|xx) integrals <==\n");
  fprintf(outfile,"\n");
  iwl_buf_init(&Buf,PSIF_MO_TEI,0.0,1,1);
  SortNoVVVV(&Buf,nfzc,nfzv,norbs,ndoccact,nvirt,islocal,options);
  iwl_buf_close(&Buf,1);
  fprintf(outfile,"\n");

  // (vv|vv)
  fprintf(outfile,"        ==> (AB|CD) <==\n");
  fprintf(outfile,"\n");
  fprintf(outfile,"        (AB|CD) contraction will be integral direct.\n");
  fprintf(outfile,"\n");

}

/**
  * Sort MO integrals, but not the VVVV ones since they weren't generated.
  */
void SortNoVVVV(struct iwlbuf *Buf,int nfzc,int nfzv,int norbs,int ndoccact,int nvirt,bool islocal,Options&options){

  double val;
  ULI o = ndoccact;
  ULI v = nvirt;
  ULI fstact = nfzc;
  ULI lstact = norbs-nfzv;

  ULI lastbuf;
  Label *lblptr;
  Value *valptr;
  ULI nocc,idx, p, q, r, s, pq, rs, pqrs;

  lblptr = Buf->labels;
  valptr = Buf->values;

  lastbuf = Buf->lastbuf;

  // buckets for integrals:
  struct integral *ijkl,*klcd,*akjc,*abci3;
  struct integral *abci5,*abcd1,*abcd2,*ijak,*ijak2;

  // available memory:
  ULI memory = Process::environment.get_memory();

  // 8 bytes for tmp, 16 for integral struct
  ULI maxelem = memory / (sizeof(double) + sizeof(struct integral));
  if (maxelem > o*v*v*v) maxelem = o*v*v*v;

  fprintf(outfile,"        CC integral sort will use                   %7.2lf mb\n",
         maxelem*(sizeof(double) + sizeof(struct integral))/1024./1024.);
  if (maxelem <o*v*v*v){
     fprintf(outfile,"       (for most efficient sort, increase memory by %7.2lf mb)\n",
         (o*v*v*v-maxelem)*(sizeof(double) + sizeof(struct integral))/1024./1024.);
  }
  fprintf(outfile,"\n");


  struct integral*integralbuffer;
  ULI nelem = maxelem/7 - 20;
  if ((nelem+20)*7>maxelem)
     integralbuffer= new integral[(nelem+7)*10];
  else
     integralbuffer= new integral[maxelem];

  ijkl  = integralbuffer;
  ijak  = integralbuffer+(nelem+20);
  klcd  = integralbuffer+(nelem+20)*2;
  akjc  = integralbuffer+(nelem+20)*3;
  abci3 = integralbuffer+(nelem+20)*4;
  abci5 = integralbuffer+(nelem+20)*5;
  ijak2 = integralbuffer+(nelem+20)*6;

  boost::shared_ptr<PSIO> psio(new PSIO());

  psio_address ijkl_addr  = PSIO_ZERO;
  psio_address klcd_addr  = PSIO_ZERO;
  psio_address akjc_addr  = PSIO_ZERO;
  psio_address ijak_addr  = PSIO_ZERO;
  psio_address ijak2_addr = PSIO_ZERO;
  psio_address abci3_addr = PSIO_ZERO;
  psio_address abci5_addr = PSIO_ZERO;

  psio->open(PSIF_DCC_IJKL,PSIO_OPEN_NEW);
  psio->close(PSIF_DCC_IJKL,1);
  psio->open(PSIF_DCC_IAJB,PSIO_OPEN_NEW);
  psio->close(PSIF_DCC_IAJB,1);
  psio->open(PSIF_DCC_IJAK,PSIO_OPEN_NEW);
  psio->close(PSIF_DCC_IJAK,1);
  psio->open(PSIF_DCC_IJAK2,PSIO_OPEN_NEW);
  psio->close(PSIF_DCC_IJAK2,1);
  psio->open(PSIF_DCC_ABCI3,PSIO_OPEN_NEW);
  psio->close(PSIF_DCC_ABCI3,1);
  psio->open(PSIF_DCC_ABCI5,PSIO_OPEN_NEW);
  psio->close(PSIF_DCC_ABCI5,1);
  psio->open(PSIF_DCC_IJAB,PSIO_OPEN_NEW);
  psio->close(PSIF_DCC_IJAB,1);
  ULI nijkl=0;
  ULI totalnijkl=0;
  ULI nijak2=0;
  ULI totalnijak2=0;
  ULI nijak=0;
  ULI totalnijak=0;
  ULI nklcd=0;
  ULI totalnklcd=0;
  ULI nakjc=0;
  ULI totalnakjc=0;
  ULI nabci3=0;
  ULI totalnabci3=0;
  ULI nabci5=0;
  ULI totalnabci5=0;

  fprintf(outfile,"        Initial sort........");fflush(outfile);
  /**
    * first buffer (read in when Buf was initialized)
    */
  for (idx=4*Buf->idx; Buf->idx<Buf->inbuf; Buf->idx++) {
      p = (ULI) lblptr[idx++];
      q = (ULI) lblptr[idx++];
      r = (ULI) lblptr[idx++];
      s = (ULI) lblptr[idx++];

      if (p < fstact || q < fstact || r < fstact || s < fstact) continue;
      if (p > lstact || q > lstact || r > lstact || s > lstact) continue;
      p -= fstact;
      q -= fstact;
      r -= fstact;
      s -= fstact;

      pq   = Position(p,q);
      rs   = Position(r,s);
      pqrs = Position(pq,rs);

      nocc = 0;
      if (p<o) nocc++;
      if (q<o) nocc++;
      if (r<o) nocc++;
      if (s<o) nocc++;

      // which type of integral?

      if (nocc==4){
         val = (double)valptr[Buf->idx];
         ijkl_terms(val,pq,rs,p,q,r,s,o,nijkl,ijkl);

         if (nijkl>=nelem){
            psio->open(PSIF_DCC_IJKL,PSIO_OPEN_OLD);
            psio->write(PSIF_DCC_IJKL,"E2ijkl",(char*)&ijkl[0],nijkl*sizeof(struct integral),ijkl_addr,&ijkl_addr);
            psio->close(PSIF_DCC_IJKL,1);
            totalnijkl+=nijkl;
            nijkl=0;
         }
      }
      else if (nocc==3){
         val = (double)valptr[Buf->idx];
         ijak_terms(val,p,q,r,s,o,v,nijak,ijak);
         if (nijak>=nelem){
            psio->open(PSIF_DCC_IJAK,PSIO_OPEN_OLD);
            psio->write(PSIF_DCC_IJAK,"E2ijak",(char*)&ijak[0],nijak*sizeof(struct integral),ijak_addr,&ijak_addr);
            psio->close(PSIF_DCC_IJAK,1);
            totalnijak+=nijak;
            nijak=0;
         }
         ijak2_terms(val,p,q,r,s,o,v,nijak2,ijak2);
         if (nijak2>=nelem){
            psio->open(PSIF_DCC_IJAK2,PSIO_OPEN_OLD);
            psio->write(PSIF_DCC_IJAK2,"E2ijak2",(char*)&ijak2[0],nijak2*sizeof(struct integral),ijak2_addr,&ijak2_addr);
            psio->close(PSIF_DCC_IJAK2,1);
            totalnijak2+=nijak2;
            nijak2=0;
         }
      } 
      else if (nocc==2){
         val = (double)valptr[Buf->idx];

         if (p<o && q>=o || p>=o && q<o){
            klcd_terms(val,pq,rs,p,q,r,s,o,v,nklcd,klcd);

            if (nklcd>=nelem){
               psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
               psio->write(PSIF_DCC_IAJB,"E2iajb",(char*)&klcd[0],nklcd*sizeof(struct integral),klcd_addr,&klcd_addr);
               psio->close(PSIF_DCC_IAJB,1);
               totalnklcd+=nklcd;
               nklcd=0;
            }
         }
         else{
            akjc_terms(val,p,q,r,s,o,v,nakjc,akjc);
            if (nakjc>=nelem){
               psio->open(PSIF_DCC_IJAB,PSIO_OPEN_OLD);
               psio->write(PSIF_DCC_IJAB,"E2ijab",(char*)&akjc[0],nakjc*sizeof(struct integral),akjc_addr,&akjc_addr);
               psio->close(PSIF_DCC_IJAB,1);
               totalnakjc+=nakjc;
               nakjc=0;
             }
         }
      }
      else if (nocc==1){
         val = (double)valptr[Buf->idx];
         abci3_terms(val,p,q,r,s,o,v,nabci3,abci3);
         if (nabci3>=nelem){
            psio->open(PSIF_DCC_ABCI3,PSIO_OPEN_OLD);
            psio->write(PSIF_DCC_ABCI3,"E2abci3",(char*)&abci3[0],nabci3*sizeof(struct integral),abci3_addr,&abci3_addr);
            psio->close(PSIF_DCC_ABCI3,1);
            totalnabci3+=nabci3;
            nabci3=0;
         }
         abci5_terms(val,p,q,r,s,o,v,nabci5,abci5);
         if (nabci5>=nelem){
            psio->open(PSIF_DCC_ABCI5,PSIO_OPEN_OLD);
            psio->write(PSIF_DCC_ABCI5,"E2abci5",(char*)&abci5[0],nabci5*sizeof(struct integral),abci5_addr,&abci5_addr);
            psio->close(PSIF_DCC_ABCI5,1);
            totalnabci5+=nabci5;
            nabci5=0;
         }
      }
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

          if (p < fstact || q < fstact || r < fstact || s < fstact) continue;
          if (p > lstact || q > lstact || r > lstact || s > lstact) continue;
          p -= fstact;
          q -= fstact;
          r -= fstact;
          s -= fstact;

          pq   = Position(p,q);
          rs   = Position(r,s);
          pqrs = Position(pq,rs);

          nocc = 0;
          if (p<o) nocc++;
          if (q<o) nocc++;
          if (r<o) nocc++;
          if (s<o) nocc++;

          // which type of integral?

          if (nocc==4){
             val = (double)valptr[Buf->idx];
             ijkl_terms(val,pq,rs,p,q,r,s,o,nijkl,ijkl);

             if (nijkl>=nelem){
                psio->open(PSIF_DCC_IJKL,PSIO_OPEN_OLD);
                psio->write(PSIF_DCC_IJKL,"E2ijkl",(char*)&ijkl[0],nijkl*sizeof(struct integral),ijkl_addr,&ijkl_addr);
                psio->close(PSIF_DCC_IJKL,1);
                totalnijkl+=nijkl;
                nijkl=0;
             }
          }
          else if (nocc==3){
             val = (double)valptr[Buf->idx];
             ijak_terms(val,p,q,r,s,o,v,nijak,ijak);
             if (nijak>=nelem){
                psio->open(PSIF_DCC_IJAK,PSIO_OPEN_OLD);
                psio->write(PSIF_DCC_IJAK,"E2ijak",(char*)&ijak[0],nijak*sizeof(struct integral),ijak_addr,&ijak_addr);
                psio->close(PSIF_DCC_IJAK,1);
                totalnijak+=nijak;
                nijak=0;
             }
             ijak2_terms(val,p,q,r,s,o,v,nijak2,ijak2);
             if (nijak2>=nelem){
                psio->open(PSIF_DCC_IJAK2,PSIO_OPEN_OLD);
                psio->write(PSIF_DCC_IJAK2,"E2ijak2",(char*)&ijak2[0],nijak2*sizeof(struct integral),ijak2_addr,&ijak2_addr);
                psio->close(PSIF_DCC_IJAK2,1);
                totalnijak2+=nijak2;
                nijak2=0;
             }
          } 
          else if (nocc==2){
             val = (double)valptr[Buf->idx];

             if (p<o && q>=o || p>=o && q<o){
                klcd_terms(val,pq,rs,p,q,r,s,o,v,nklcd,klcd);

                if (nklcd>=nelem){
                   psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
                   psio->write(PSIF_DCC_IAJB,"E2iajb",(char*)&klcd[0],nklcd*sizeof(struct integral),klcd_addr,&klcd_addr);
                   psio->close(PSIF_DCC_IAJB,1);
                   totalnklcd+=nklcd;
                   nklcd=0;
                }
             }
             else{
                akjc_terms(val,p,q,r,s,o,v,nakjc,akjc);

                if (nakjc>=nelem){
                   psio->open(PSIF_DCC_IJAB,PSIO_OPEN_OLD);
                   psio->write(PSIF_DCC_IJAB,"E2ijab",(char*)&akjc[0],nakjc*sizeof(struct integral),akjc_addr,&akjc_addr);
                   psio->close(PSIF_DCC_IJAB,1);
                   totalnakjc+=nakjc;
                   nakjc=0;
                }
             }
          }
          else if (nocc==1){
             val = (double)valptr[Buf->idx];
             abci3_terms(val,p,q,r,s,o,v,nabci3,abci3);
             if (nabci3>=nelem){
                psio->open(PSIF_DCC_ABCI3,PSIO_OPEN_OLD);
                psio->write(PSIF_DCC_ABCI3,"E2abci3",(char*)&abci3[0],nabci3*sizeof(struct integral),abci3_addr,&abci3_addr);
                psio->close(PSIF_DCC_ABCI3,1);
                totalnabci3+=nabci3;
                nabci3=0;
             }
             abci5_terms(val,p,q,r,s,o,v,nabci5,abci5);
             if (nabci5>=nelem){
                psio->open(PSIF_DCC_ABCI5,PSIO_OPEN_OLD);
                psio->write(PSIF_DCC_ABCI5,"E2abci5",(char*)&abci5[0],nabci5*sizeof(struct integral),abci5_addr,&abci5_addr);
                psio->close(PSIF_DCC_ABCI5,1);
                totalnabci5+=nabci5;
                nabci5=0;
             }
          }
      }
  }
  fprintf(outfile,"done.\n\n");fflush(outfile);
  /**
    * write any leftover bits that might not have been dumped to disk
    */
  if (nabci5!=0){
     psio->open(PSIF_DCC_ABCI5,PSIO_OPEN_OLD);
     psio->write(PSIF_DCC_ABCI5,"E2abci5",(char*)&abci5[0],nabci5*sizeof(struct integral),abci5_addr,&abci5_addr);
     psio->close(PSIF_DCC_ABCI5,1);
     totalnabci5+=nabci5;
     nabci5=0;
  }
  if (nabci3!=0){
     psio->open(PSIF_DCC_ABCI3,PSIO_OPEN_OLD);
     psio->write(PSIF_DCC_ABCI3,"E2abci3",(char*)&abci3[0],nabci3*sizeof(struct integral),abci3_addr,&abci3_addr);
     psio->close(PSIF_DCC_ABCI3,1);
     totalnabci3+=nabci3;
     nabci3=0;
  }
  if (nakjc!=0){
     psio->open(PSIF_DCC_IJAB,PSIO_OPEN_OLD);
     psio->write(PSIF_DCC_IJAB,"E2ijab",(char*)&akjc[0],nakjc*sizeof(struct integral),akjc_addr,&akjc_addr);
     psio->close(PSIF_DCC_IJAB,1);
     totalnakjc+=nakjc;
     nakjc=0;
  }
  if (nklcd!=0){
     psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
     psio->write(PSIF_DCC_IAJB,"E2iajb",(char*)&klcd[0],nklcd*sizeof(struct integral),klcd_addr,&klcd_addr);
     psio->close(PSIF_DCC_IAJB,1);
     totalnklcd+=nklcd;
     nklcd=0;
  }
  if (nijkl!=0){
     psio->open(PSIF_DCC_IJKL,PSIO_OPEN_OLD);
     psio->write(PSIF_DCC_IJKL,"E2ijkl",(char*)&ijkl[0],nijkl*sizeof(struct integral),ijkl_addr,&ijkl_addr);
     psio->close(PSIF_DCC_IJKL,1);
     totalnijkl+=nijkl;
     nijkl=0;
  }
  if (nijak!=0){
     psio->open(PSIF_DCC_IJAK,PSIO_OPEN_OLD);
     psio->write(PSIF_DCC_IJAK,"E2ijak",(char*)&ijak[0],nijak*sizeof(struct integral),ijak_addr,&ijak_addr);
     psio->close(PSIF_DCC_IJAK,1);
     totalnijak+=nijak;
     nijak=0;
  }
  if (nijak2!=0){
     psio->open(PSIF_DCC_IJAK2,PSIO_OPEN_OLD);
     psio->write(PSIF_DCC_IJAK2,"E2ijak2",(char*)&ijak2[0],nijak2*sizeof(struct integral),ijak2_addr,&ijak2_addr);
     psio->close(PSIF_DCC_IJAK2,1);
     totalnijak2+=nijak2;
     nijak2=0;
  }

  /**
    * sort values in each of the files
    */
  double *tmp;
  tmp = new double[maxelem];

  fprintf(outfile,"        Sort (IJ|KL)........");fflush(outfile);
  SortBlock(totalnijkl,o*o*o*o,integralbuffer,tmp,PSIF_DCC_IJKL,"E2ijkl",maxelem);
  fprintf(outfile,"done.\n");fflush(outfile);
  fprintf(outfile,"        Sort (IJ|KA) 1/2....");fflush(outfile);
  SortBlock(totalnijak,o*o*o*v,integralbuffer,tmp,PSIF_DCC_IJAK,"E2ijak",maxelem);
  fprintf(outfile,"done.\n");fflush(outfile);
  fprintf(outfile,"        Sort (IJ|KA) 2/2....");fflush(outfile);
  SortBlock(totalnijak2,o*o*o*v,integralbuffer,tmp,PSIF_DCC_IJAK2,"E2ijak2",maxelem);
  fprintf(outfile,"done.\n");fflush(outfile);
  fprintf(outfile,"        Sort (IA|JB)........");fflush(outfile);
  SortBlock(totalnklcd,o*o*v*v,integralbuffer,tmp,PSIF_DCC_IAJB,"E2iajb",maxelem);
  fprintf(outfile,"done.\n");fflush(outfile);
  fprintf(outfile,"        Sort (IJ|AB)........");fflush(outfile);
  SortBlock(totalnakjc,o*o*v*v,integralbuffer,tmp,PSIF_DCC_IJAB,"E2ijab",maxelem);
  fprintf(outfile,"done.\n");fflush(outfile);
  fprintf(outfile,"        Sort (IA|BC) 1/2....");fflush(outfile);
  SortBlock(totalnabci3,o*v*v*v,integralbuffer,tmp,PSIF_DCC_ABCI3,"E2abci3",maxelem);
  fprintf(outfile,"done.\n");fflush(outfile);
  fprintf(outfile,"        Sort (IA|BC) 1/2....");fflush(outfile);
  SortBlock(totalnabci5,o*v*v*v,integralbuffer,tmp,PSIF_DCC_ABCI5,"E2abci5",maxelem);
  fprintf(outfile,"done.\n");fflush(outfile);

  delete integralbuffer;
  delete tmp;
}

}} // end of namespace
