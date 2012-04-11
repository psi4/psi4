#include"psi4-dec.h"
#include<psifiles.h>
#include<libiwl/iwl.h>
#include <libpsio/psio.hpp>

#include"ccsd.h"
#include"sort.h"
#include"blas.h"

using namespace psi;

namespace psi{

void OutOfCoreSort(int nfzc,int nfzv,int norbs,int ndoccact,int nvirt,bool islocal,Options&options){
  struct iwlbuf Buf; 
  // initialize buffer
  iwl_buf_init(&Buf,PSIF_MO_TEI,0.0,1,1);

  fprintf(outfile,"\n        Begin CC integral sort\n\n");

  //sort
  Sort(&Buf,nfzc,nfzv,norbs,ndoccact,nvirt,islocal,options);

  iwl_buf_close(&Buf,1);
}
/**
  * out-of-core integral sort.  requires 2o^2v^2 doubles +o^2v^2 ULIs storage
  */
void Sort(struct iwlbuf *Buf,int nfzc,int nfzv,int norbs,int ndoccact,int nvirt,bool islocal,Options&options){

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
  struct integral *ijkl,*klcd,*akjc,*abci1,*abci3;
  struct integral *abci5,*abcd1,*abcd2,*ijak,*ijak2;

  // available memory:
  ULI memory = Process::environment.get_memory();

  // 8 bytes for tmp, 16 for integral struct
  ULI maxelem = memory / (sizeof(double) + sizeof(struct integral));
  if (maxelem > v*(v+1)/2*v*(v+1)/2) maxelem = v*(v+1)/2*v*(v+1)/2;

  fprintf(outfile,"        CC integral sort will use                   %7.2lf mb\n",
         maxelem*(sizeof(double) + sizeof(struct integral))/1024./1024.);
  if (maxelem <v*(v+1)/2*v*(v+1)/2){
     fprintf(outfile,"       (for most efficient sort, increase memory by %7.2lf mb)\n",
         (v*(v+1)/2*v*(v+1)/2-maxelem)*(sizeof(double) + sizeof(struct integral))/1024./1024.);
  }
  fprintf(outfile,"\n");


  struct integral*integralbuffer;
  ULI nelem = maxelem/10 - 20;
  if ((nelem+20)*10>maxelem)
     integralbuffer= new integral[(nelem+20)*10];
  else
     integralbuffer= new integral[maxelem];

  ijkl  = integralbuffer;
  ijak  = integralbuffer+(nelem+20);
  klcd  = integralbuffer+(nelem+20)*2;
  akjc  = integralbuffer+(nelem+20)*3;
  abci1 = integralbuffer+(nelem+20)*4;
  abci3 = integralbuffer+(nelem+20)*5;
  abci5 = integralbuffer+(nelem+20)*6;
  abcd1 = integralbuffer+(nelem+20)*7;
  abcd2 = integralbuffer+(nelem+20)*8;
  ijak2 = integralbuffer+(nelem+20)*9;

  boost::shared_ptr<PSIO> psio(new PSIO());

  psio_address ijkl_addr  = PSIO_ZERO;
  psio_address klcd_addr  = PSIO_ZERO;
  psio_address akjc_addr  = PSIO_ZERO;
  psio_address ijak_addr  = PSIO_ZERO;
  psio_address ijak2_addr = PSIO_ZERO;
  psio_address abci1_addr = PSIO_ZERO;
  psio_address abci2_addr = PSIO_ZERO;
  psio_address abci3_addr = PSIO_ZERO;
  psio_address abci5_addr = PSIO_ZERO;
  psio_address abcd1_addr = PSIO_ZERO;
  psio_address abcd2_addr = PSIO_ZERO;

  psio->open(PSIF_IJKL,PSIO_OPEN_NEW);
  psio->close(PSIF_IJKL,1);
  psio->open(PSIF_KLCD,PSIO_OPEN_NEW);
  psio->close(PSIF_KLCD,1);
  psio->open(PSIF_IJAK,PSIO_OPEN_NEW);
  psio->close(PSIF_IJAK,1);
  psio->open(PSIF_IJAK2,PSIO_OPEN_NEW);
  psio->close(PSIF_IJAK2,1);
  psio->open(PSIF_ABCI,PSIO_OPEN_NEW);
  psio->close(PSIF_ABCI,1);
  psio->open(PSIF_ABCI2,PSIO_OPEN_NEW);
  psio->close(PSIF_ABCI2,1);
  psio->open(PSIF_ABCI3,PSIO_OPEN_NEW);
  psio->close(PSIF_ABCI3,1);
  psio->open(PSIF_ABCD1,PSIO_OPEN_NEW);
  psio->close(PSIF_ABCD1,1);
  psio->open(PSIF_ABCD2,PSIO_OPEN_NEW);
  psio->close(PSIF_ABCD2,1);
  psio->open(PSIF_AKJC2,PSIO_OPEN_NEW);
  psio->close(PSIF_AKJC2,1);
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
  ULI nabci1=0;
  ULI totalnabci1=0;
  ULI nabci3=0;
  ULI totalnabci3=0;
  ULI nabci5=0;
  ULI totalnabci5=0;
  ULI nabcd1=0;
  ULI totalnabcd1=0;
  ULI nabcd2=0;
  ULI totalnabcd2=0;

  fprintf(outfile,"        Initial sort.....");fflush(outfile);
  /**
    * first buffer (read in when Buf was initialized)
    */
  for (idx=4*Buf->idx; Buf->idx<Buf->inbuf; Buf->idx++) {
      p = (ULI) lblptr[idx++];
      q = (ULI) lblptr[idx++];
      r = (ULI) lblptr[idx++];
      s = (ULI) lblptr[idx++];

      if (islocal){
         if (p < fstact || q < fstact || r < fstact || s < fstact) continue;
         if (p > lstact || q > lstact || r > lstact || s > lstact) continue;
         p -= fstact;
         q -= fstact;
         r -= fstact;
         s -= fstact;
      }

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
            psio->open(PSIF_IJKL,PSIO_OPEN_OLD);
            psio->write(PSIF_IJKL,"E2ijkl",(char*)&ijkl[0],nijkl*sizeof(struct integral),ijkl_addr,&ijkl_addr);
            psio->close(PSIF_IJKL,1);
            totalnijkl+=nijkl;
            nijkl=0;
         }
      }
      else if (nocc==3){
         val = (double)valptr[Buf->idx];
         ijak_terms(val,p,q,r,s,o,v,nijak,ijak);
         if (nijak>=nelem){
            psio->open(PSIF_IJAK,PSIO_OPEN_OLD);
            psio->write(PSIF_IJAK,"E2ijak",(char*)&ijak[0],nijak*sizeof(struct integral),ijak_addr,&ijak_addr);
            psio->close(PSIF_IJAK,1);
            totalnijak+=nijak;
            nijak=0;
         }
         ijak2_terms(val,p,q,r,s,o,v,nijak2,ijak2);
         if (nijak2>=nelem){
            psio->open(PSIF_IJAK2,PSIO_OPEN_OLD);
            psio->write(PSIF_IJAK2,"E2ijak2",(char*)&ijak2[0],nijak2*sizeof(struct integral),ijak2_addr,&ijak2_addr);
            psio->close(PSIF_IJAK2,1);
            totalnijak2+=nijak2;
            nijak2=0;
         }
      } 
      else if (nocc==2){
         val = (double)valptr[Buf->idx];

         if (p<o && q>=o || p>=o && q<o){
            klcd_terms(val,pq,rs,p,q,r,s,o,v,nklcd,klcd);

            if (nklcd>=nelem){
               psio->open(PSIF_KLCD,PSIO_OPEN_OLD);
               psio->write(PSIF_KLCD,"E2klcd",(char*)&klcd[0],nklcd*sizeof(struct integral),klcd_addr,&klcd_addr);
               psio->close(PSIF_KLCD,1);
               totalnklcd+=nklcd;
               nklcd=0;
            }
         }
         else{
            akjc_terms(val,p,q,r,s,o,v,nakjc,akjc);
            if (nakjc>=nelem){
               psio->open(PSIF_AKJC2,PSIO_OPEN_OLD);
               psio->write(PSIF_AKJC2,"E2akjc2",(char*)&akjc[0],nakjc*sizeof(struct integral),akjc_addr,&akjc_addr);
               psio->close(PSIF_AKJC2,1);
               totalnakjc+=nakjc;
               nakjc=0;
             }
         }
      }
      else if (nocc==1){
         val = (double)valptr[Buf->idx];
         abci1_terms(val,p,q,r,s,o,v,nabci1,abci1);
         if (nabci1>=nelem){
            psio->open(PSIF_ABCI,PSIO_OPEN_OLD);
            psio->write(PSIF_ABCI,"E2abci",(char*)&abci1[0],nabci1*sizeof(struct integral),abci1_addr,&abci1_addr);
            psio->close(PSIF_ABCI,1);
            totalnabci1+=nabci1;
            nabci1=0;
         }
         abci3_terms(val,p,q,r,s,o,v,nabci3,abci3);
         if (nabci3>=nelem){
            psio->open(PSIF_ABCI3,PSIO_OPEN_OLD);
            psio->write(PSIF_ABCI3,"E2abci3",(char*)&abci3[0],nabci3*sizeof(struct integral),abci3_addr,&abci3_addr);
            psio->close(PSIF_ABCI3,1);
            totalnabci3+=nabci3;
            nabci3=0;
         }
         abci5_terms(val,p,q,r,s,o,v,nabci5,abci5);
         if (nabci5>=nelem){
            psio->open(PSIF_ABCI2,PSIO_OPEN_OLD);
            psio->write(PSIF_ABCI2,"E2abci2",(char*)&abci5[0],nabci5*sizeof(struct integral),abci5_addr,&abci5_addr);
            psio->close(PSIF_ABCI2,1);
            totalnabci5+=nabci5;
            nabci5=0;
         }
      }
      else if (nocc==0){
         val = (double)valptr[Buf->idx];
         if (options.get_bool("VABCD_PACKED")){
            abcd1_terms(val,pq,rs,p,q,r,s,o,v,nabcd1,abcd1);
            if (nabcd1>=nelem){
               psio->open(PSIF_ABCD1,PSIO_OPEN_OLD);
               psio->write(PSIF_ABCD1,"E2abcd1",(char*)&abcd1[0],nabcd1*sizeof(struct integral),abcd1_addr,&abcd1_addr);
               psio->close(PSIF_ABCD1,1);
               totalnabcd1+=nabcd1;
               nabcd1=0;
            }
            abcd2_terms(val,pq,rs,p,q,r,s,o,v,nabcd2,abcd2);
            if (nabcd2>=nelem){
               psio->open(PSIF_ABCD2,PSIO_OPEN_OLD);
               psio->write(PSIF_ABCD2,"E2abcd2",(char*)&abcd2[0],nabcd2*sizeof(struct integral),abcd2_addr,&abcd2_addr);
               psio->close(PSIF_ABCD2,1);
               totalnabcd2+=nabcd2;
               nabcd2=0;
            }
         }else{
            abcd3_terms(val,pq,rs,p,q,r,s,o,v,nabcd1,abcd1);
            if (nabcd1>=nelem){
               psio->open(PSIF_ABCD1,PSIO_OPEN_OLD);
               psio->write(PSIF_ABCD1,"E2abcd1",(char*)&abcd1[0],nabcd1*sizeof(struct integral),abcd1_addr,&abcd1_addr);
               psio->close(PSIF_ABCD1,1);
               totalnabcd1+=nabcd1;
               nabcd1=0;
            }
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

          if (islocal){
             if (p < fstact || q < fstact || r < fstact || s < fstact) continue;
             if (p > lstact || q > lstact || r > lstact || s > lstact) continue;
             p -= fstact;
             q -= fstact;
             r -= fstact;
             s -= fstact;
          }

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
                psio->open(PSIF_IJKL,PSIO_OPEN_OLD);
                psio->write(PSIF_IJKL,"E2ijkl",(char*)&ijkl[0],nijkl*sizeof(struct integral),ijkl_addr,&ijkl_addr);
                psio->close(PSIF_IJKL,1);
                totalnijkl+=nijkl;
                nijkl=0;
             }
          }
          else if (nocc==3){
             val = (double)valptr[Buf->idx];
             ijak_terms(val,p,q,r,s,o,v,nijak,ijak);
             if (nijak>=nelem){
                psio->open(PSIF_IJAK,PSIO_OPEN_OLD);
                psio->write(PSIF_IJAK,"E2ijak",(char*)&ijak[0],nijak*sizeof(struct integral),ijak_addr,&ijak_addr);
                psio->close(PSIF_IJAK,1);
                totalnijak+=nijak;
                nijak=0;
             }
             ijak2_terms(val,p,q,r,s,o,v,nijak2,ijak2);
             if (nijak2>=nelem){
                psio->open(PSIF_IJAK2,PSIO_OPEN_OLD);
                psio->write(PSIF_IJAK2,"E2ijak2",(char*)&ijak2[0],nijak2*sizeof(struct integral),ijak2_addr,&ijak2_addr);
                psio->close(PSIF_IJAK2,1);
                totalnijak2+=nijak2;
                nijak2=0;
             }
          } 
          else if (nocc==2){
             val = (double)valptr[Buf->idx];

             if (p<o && q>=o || p>=o && q<o){
                klcd_terms(val,pq,rs,p,q,r,s,o,v,nklcd,klcd);

                if (nklcd>=nelem){
                   psio->open(PSIF_KLCD,PSIO_OPEN_OLD);
                   psio->write(PSIF_KLCD,"E2klcd",(char*)&klcd[0],nklcd*sizeof(struct integral),klcd_addr,&klcd_addr);
                   psio->close(PSIF_KLCD,1);
                   totalnklcd+=nklcd;
                   nklcd=0;
                }
             }
             else{
                akjc_terms(val,p,q,r,s,o,v,nakjc,akjc);

                if (nakjc>=nelem){
                   psio->open(PSIF_AKJC2,PSIO_OPEN_OLD);
                   psio->write(PSIF_AKJC2,"E2akjc2",(char*)&akjc[0],nakjc*sizeof(struct integral),akjc_addr,&akjc_addr);
                   psio->close(PSIF_AKJC2,1);
                   totalnakjc+=nakjc;
                   nakjc=0;
                }
             }
          }
          else if (nocc==1){
             val = (double)valptr[Buf->idx];
             abci1_terms(val,p,q,r,s,o,v,nabci1,abci1);
             if (nabci1>=nelem){
                psio->open(PSIF_ABCI,PSIO_OPEN_OLD);
                psio->write(PSIF_ABCI,"E2abci",(char*)&abci1[0],nabci1*sizeof(struct integral),abci1_addr,&abci1_addr);
                psio->close(PSIF_ABCI,1);
                totalnabci1+=nabci1;
                nabci1=0;
             }
             abci3_terms(val,p,q,r,s,o,v,nabci3,abci3);
             if (nabci3>=nelem){
                psio->open(PSIF_ABCI3,PSIO_OPEN_OLD);
                psio->write(PSIF_ABCI3,"E2abci3",(char*)&abci3[0],nabci3*sizeof(struct integral),abci3_addr,&abci3_addr);
                psio->close(PSIF_ABCI3,1);
                totalnabci3+=nabci3;
                nabci3=0;
             }
             abci5_terms(val,p,q,r,s,o,v,nabci5,abci5);
             if (nabci5>=nelem){
                psio->open(PSIF_ABCI2,PSIO_OPEN_OLD);
                psio->write(PSIF_ABCI2,"E2abci2",(char*)&abci5[0],nabci5*sizeof(struct integral),abci5_addr,&abci5_addr);
                psio->close(PSIF_ABCI2,1);
                totalnabci5+=nabci5;
                nabci5=0;
             }
          }
          else if (nocc==0){
             val = (double)valptr[Buf->idx];
             if (options.get_bool("VABCD_PACKED")){
                abcd1_terms(val,pq,rs,p,q,r,s,o,v,nabcd1,abcd1);
                if (nabcd1>=nelem){
                   psio->open(PSIF_ABCD1,PSIO_OPEN_OLD);
                   psio->write(PSIF_ABCD1,"E2abcd1",(char*)&abcd1[0],nabcd1*sizeof(struct integral),abcd1_addr,&abcd1_addr);
                   psio->close(PSIF_ABCD1,1);
                   totalnabcd1+=nabcd1;
                   nabcd1=0;
                }
                abcd2_terms(val,pq,rs,p,q,r,s,o,v,nabcd2,abcd2);
                if (nabcd2>=nelem){
                   psio->open(PSIF_ABCD2,PSIO_OPEN_OLD);
                   psio->write(PSIF_ABCD2,"E2abcd2",(char*)&abcd2[0],nabcd2*sizeof(struct integral),abcd2_addr,&abcd2_addr);
                   psio->close(PSIF_ABCD2,1);
                   totalnabcd2+=nabcd2;
                   nabcd2=0;
                }
             }else{
                abcd3_terms(val,pq,rs,p,q,r,s,o,v,nabcd1,abcd1);
                if (nabcd1>=nelem){
                   psio->open(PSIF_ABCD1,PSIO_OPEN_OLD);
                   psio->write(PSIF_ABCD1,"E2abcd1",(char*)&abcd1[0],nabcd1*sizeof(struct integral),abcd1_addr,&abcd1_addr);
                   psio->close(PSIF_ABCD1,1);
                   totalnabcd1+=nabcd1;
                   nabcd1=0;
                }
             }
          }

      }

  }
  fprintf(outfile,"done.\n\n");fflush(outfile);
  /**
    * write any leftover bits that might not have been dumped to disk
    */
  if (nabcd2!=0){
     psio->open(PSIF_ABCD2,PSIO_OPEN_OLD);
     psio->write(PSIF_ABCD2,"E2abcd2",(char*)&abcd2[0],nabcd2*sizeof(struct integral),abcd2_addr,&abcd2_addr);
     psio->close(PSIF_ABCD2,1);
     totalnabcd2+=nabcd2;
     nabcd2=0;
  }
  if (nabcd1!=0){
     psio->open(PSIF_ABCD1,PSIO_OPEN_OLD);
     psio->write(PSIF_ABCD1,"E2abcd1",(char*)&abcd1[0],nabcd1*sizeof(struct integral),abcd1_addr,&abcd1_addr);
     psio->close(PSIF_ABCD1,1);
     totalnabcd1+=nabcd1;
     nabcd1=0;
  }
  if (nabci5!=0){
     psio->open(PSIF_ABCI2,PSIO_OPEN_OLD);
     psio->write(PSIF_ABCI2,"E2abci2",(char*)&abci5[0],nabci5*sizeof(struct integral),abci5_addr,&abci5_addr);
     psio->close(PSIF_ABCI2,1);
     totalnabci5+=nabci5;
     nabci5=0;
  }
  if (nabci3!=0){
     psio->open(PSIF_ABCI3,PSIO_OPEN_OLD);
     psio->write(PSIF_ABCI3,"E2abci3",(char*)&abci3[0],nabci3*sizeof(struct integral),abci3_addr,&abci3_addr);
     psio->close(PSIF_ABCI3,1);
     totalnabci3+=nabci3;
     nabci3=0;
  }
  if (nabci1!=0){
     psio->open(PSIF_ABCI,PSIO_OPEN_OLD);
     psio->write(PSIF_ABCI,"E2abci",(char*)&abci1[0],nabci1*sizeof(struct integral),abci1_addr,&abci1_addr);
     psio->close(PSIF_ABCI,1);
     totalnabci1+=nabci1;
     nabci1=0;
  }
  if (nakjc!=0){
     psio->open(PSIF_AKJC2,PSIO_OPEN_OLD);
     psio->write(PSIF_AKJC2,"E2akjc2",(char*)&akjc[0],nakjc*sizeof(struct integral),akjc_addr,&akjc_addr);
     psio->close(PSIF_AKJC2,1);
     totalnakjc+=nakjc;
     nakjc=0;
  }
  if (nklcd!=0){
     psio->open(PSIF_KLCD,PSIO_OPEN_OLD);
     psio->write(PSIF_KLCD,"E2klcd",(char*)&klcd[0],nklcd*sizeof(struct integral),klcd_addr,&klcd_addr);
     psio->close(PSIF_KLCD,1);
     totalnklcd+=nklcd;
     nklcd=0;
  }
  if (nijkl!=0){
     psio->open(PSIF_IJKL,PSIO_OPEN_OLD);
     psio->write(PSIF_IJKL,"E2ijkl",(char*)&ijkl[0],nijkl*sizeof(struct integral),ijkl_addr,&ijkl_addr);
     psio->close(PSIF_IJKL,1);
     totalnijkl+=nijkl;
     nijkl=0;
  }
  if (nijak!=0){
     psio->open(PSIF_IJAK,PSIO_OPEN_OLD);
     psio->write(PSIF_IJAK,"E2ijak",(char*)&ijak[0],nijak*sizeof(struct integral),ijak_addr,&ijak_addr);
     psio->close(PSIF_IJAK,1);
     totalnijak+=nijak;
     nijak=0;
  }
  if (nijak2!=0){
     psio->open(PSIF_IJAK2,PSIO_OPEN_OLD);
     psio->write(PSIF_IJAK2,"E2ijak2",(char*)&ijak2[0],nijak2*sizeof(struct integral),ijak2_addr,&ijak2_addr);
     psio->close(PSIF_IJAK2,1);
     totalnijak2+=nijak2;
     nijak2=0;
  }

  /**
    * sort values in each of the files
    */
  double *tmp;
  tmp = new double[maxelem];

  fprintf(outfile,"        IJKL block.......");fflush(outfile);
  SortBlock(totalnijkl,o*o*o*o,integralbuffer,tmp,PSIF_IJKL,"E2ijkl",maxelem);
  fprintf(outfile,"done.\n");fflush(outfile);
  fprintf(outfile,"        IJAK block 1/2...");fflush(outfile);
  SortBlock(totalnijak,o*o*o*v,integralbuffer,tmp,PSIF_IJAK,"E2ijak",maxelem);
  fprintf(outfile,"done.\n");fflush(outfile);
  fprintf(outfile,"        IJAK block 2/2...");fflush(outfile);
  SortBlock(totalnijak2,o*o*o*v,integralbuffer,tmp,PSIF_IJAK2,"E2ijak2",maxelem);
  fprintf(outfile,"done.\n");fflush(outfile);
  fprintf(outfile,"        KCLD block.......");fflush(outfile);
  SortBlock(totalnklcd,o*o*v*v,integralbuffer,tmp,PSIF_KLCD,"E2klcd",maxelem);
  fprintf(outfile,"done.\n");fflush(outfile);
  fprintf(outfile,"        KLCD block.......");fflush(outfile);
  SortBlock(totalnakjc,o*o*v*v,integralbuffer,tmp,PSIF_AKJC2,"E2akjc2",maxelem);
  fprintf(outfile,"done.\n");fflush(outfile);
  fprintf(outfile,"        ABCI block 1/3...");fflush(outfile);
  SortBlock(totalnabci1,o*v*v*v,integralbuffer,tmp,PSIF_ABCI,"E2abci",maxelem);
  fprintf(outfile,"done.\n");fflush(outfile);
  fprintf(outfile,"        ABCI block 2/3...");fflush(outfile);
  SortBlock(totalnabci3,o*v*v*v,integralbuffer,tmp,PSIF_ABCI3,"E2abci3",maxelem);
  fprintf(outfile,"done.\n");fflush(outfile);
  fprintf(outfile,"        ABCI block 3/3...");fflush(outfile);
  SortBlock(totalnabci5,o*v*v*v,integralbuffer,tmp,PSIF_ABCI2,"E2abci2",maxelem);
  fprintf(outfile,"done.\n");fflush(outfile);


  if (options.get_bool("VABCD_PACKED")){
     fprintf(outfile,"        ABCD block 1/2...");fflush(outfile);
     SortBlock(totalnabcd1,v*(v+1)/2*v*(v+1)/2,integralbuffer,tmp,PSIF_ABCD1,"E2abcd1",maxelem);
     fprintf(outfile,"done.\n");fflush(outfile);
     fprintf(outfile,"        ABCD block 2/2...");fflush(outfile);
     SortBlock(totalnabcd2,v*(v+1)/2*v*(v+1)/2,integralbuffer,tmp,PSIF_ABCD2,"E2abcd2",maxelem);
     fprintf(outfile,"done.\n");fflush(outfile);
  }else{
     fprintf(outfile,"        ABCD block    ...");fflush(outfile);
     SortBlock(totalnabcd1,v*v*v*v,integralbuffer,tmp,PSIF_ABCD1,"E2abcd1",maxelem);
     fprintf(outfile,"done.\n");fflush(outfile);
  }
  fprintf(outfile,"\n");

  delete integralbuffer;

  double *tmp2;
  tmp2 = new double[maxelem];
  /**
    *  Sort ABCI2 integrals (actually, just ABCI2-2*ABCI3)
    */

  ULI nbins,binsize,lastbin;
  for (ULI i=1; i<=o*v*v*v; i++){
      if (maxelem>=(double)o*v*v*v/i){
         binsize = o*v*v*v/i;
         if (i*binsize < o*v*v*v) binsize++;
         nbins = i;
         break;
      }
  }
  lastbin = o*v*v*v - (nbins-1)*binsize;
  psio->open(PSIF_ABCI3,PSIO_OPEN_OLD);
  psio->open(PSIF_ABCI2,PSIO_OPEN_OLD);
  if (islocal || options.get_bool("TRIPLES_LOW_MEMORY")){
     psio->open(PSIF_ABCI4,PSIO_OPEN_NEW);
  }

  abci2_addr = PSIO_ZERO;
  abci3_addr = PSIO_ZERO;
  abci5_addr = PSIO_ZERO;
  psio_address abci4_addr = PSIO_ZERO;
  for (ULI i=0; i<nbins-1; i++){
      psio->read(PSIF_ABCI3,"E2abci3",(char*)&tmp[0],binsize*sizeof(double),abci3_addr,&abci3_addr);
      psio->read(PSIF_ABCI2,"E2abci2",(char*)&tmp2[0],binsize*sizeof(double),abci5_addr,&abci5_addr);
      // this is for the local triples
      if (islocal || options.get_bool("TRIPLES_LOW_MEMORY")){
         psio->write(PSIF_ABCI4,"E2abci4",(char*)&tmp2[0],binsize*sizeof(double),abci4_addr,&abci4_addr);
      }
      F_DAXPY(binsize,-2.0,tmp,1,tmp2,1);
      psio->write(PSIF_ABCI2,"E2abci2",(char*)&tmp2[0],binsize*sizeof(double),abci2_addr,&abci2_addr);
  }
  psio->read(PSIF_ABCI3,"E2abci3",(char*)&tmp[0],lastbin*sizeof(double),abci3_addr,&abci3_addr);
  psio->read(PSIF_ABCI2,"E2abci2",(char*)&tmp2[0],lastbin*sizeof(double),abci5_addr,&abci5_addr);
  if (islocal || options.get_bool("TRIPLES_LOW_MEMORY")){
     psio->write(PSIF_ABCI4,"E2abci4",(char*)&tmp2[0],lastbin*sizeof(double),abci4_addr,&abci4_addr);
  }
  F_DAXPY(lastbin,-2.0,tmp,1,tmp2,1);
  psio->write(PSIF_ABCI2,"E2abci2",(char*)&tmp2[0],lastbin*sizeof(double),abci2_addr,&abci2_addr);
  psio->close(PSIF_ABCI2,1);
  psio->close(PSIF_ABCI3,1);
  if (islocal || options.get_bool("TRIPLES_LOW_MEMORY")){
     psio->close(PSIF_ABCI4,1);
  }

  /**
    *  Combine ABCD1 and ABCD2 integrals if SJS packing
    */
  if (options.get_bool("VABCD_PACKED")){
     for (ULI i=1; i<=v*(v+1)/2*v*(v+1)/2; i++){
         if (maxelem>=(double)v*(v+1)/2*v*(v+1)/2/i){
            binsize = v*(v+1)/2*v*(v+1)/2/i;
            if (i*binsize < v*(v+1)/2*v*(v+1)/2) binsize++;
            nbins = i;
            break;
         }
     }
     lastbin = v*(v+1)/2*v*(v+1)/2 - (nbins-1)*binsize;
     psio->open(PSIF_ABCD1,PSIO_OPEN_OLD);
     psio->open(PSIF_ABCD2,PSIO_OPEN_OLD);
     psio_address abcd1_again = PSIO_ZERO;
     psio_address abcd1_new = PSIO_ZERO;
     psio_address abcd2_new = PSIO_ZERO;
     abcd1_addr = abcd2_addr = PSIO_ZERO;
     for (ULI i=0; i<nbins-1; i++){
         psio->read(PSIF_ABCD1,"E2abcd1",(char*)&tmp[0],binsize*sizeof(double),abcd1_addr,&abcd1_addr);
         psio->read(PSIF_ABCD2,"E2abcd2",(char*)&tmp2[0],binsize*sizeof(double),abcd2_addr,&abcd2_addr);
         F_DAXPY(binsize,-1.0,tmp2,1,tmp,1);
         psio->write(PSIF_ABCD2,"E2abcd2",(char*)&tmp[0],binsize*sizeof(double),abcd2_new,&abcd2_new);
         psio->read(PSIF_ABCD1,"E2abcd1",(char*)&tmp[0],binsize*sizeof(double),abcd1_again,&abcd1_again);
         F_DAXPY(binsize,1.0,tmp2,1,tmp,1);
         psio->write(PSIF_ABCD1,"E2abcd1",(char*)&tmp[0],binsize*sizeof(double),abcd1_new,&abcd1_new);
     }
     psio->read(PSIF_ABCD1,"E2abcd1",(char*)&tmp[0],lastbin*sizeof(double),abcd1_addr,&abcd1_addr);
     psio->read(PSIF_ABCD2,"E2abcd2",(char*)&tmp2[0],lastbin*sizeof(double),abcd2_addr,&abcd2_addr);
     F_DAXPY(lastbin,-1.0,tmp2,1,tmp,1);
     psio->write(PSIF_ABCD2,"E2abcd2",(char*)&tmp[0],lastbin*sizeof(double),abcd2_new,&abcd2_new);
     psio->read(PSIF_ABCD1,"E2abcd1",(char*)&tmp[0],lastbin*sizeof(double),abcd1_again,&abcd1_again);
     F_DAXPY(lastbin,1.0,tmp2,1,tmp,1);
     psio->write(PSIF_ABCD1,"E2abcd1",(char*)&tmp[0],lastbin*sizeof(double),abcd1_new,&abcd1_new);
     psio->close(PSIF_ABCD1,1);
     psio->close(PSIF_ABCD2,1);
  }

  delete tmp;
  delete tmp2;
}

void klcd_terms(double val,ULI pq,ULI rs,ULI p,ULI q,ULI r,ULI s,ULI o,ULI v,ULI&nklcd,struct integral*klcd){
  ULI k,l,c,d;
  if (p<o){
     k=p;
     c=q-o;
     if (r<o){
        l=r;
        d=s-o;
     }
     else{
        d=r-o;
        l=s;
     }
  }
  else{
     c=p-o;
     k=q;
     if (r<o){
        l=r;
        d=s-o;
     }
     else{
        d=r-o;
        l=s;
     }
  }
  klcd[nklcd].ind   = k*o*v*v+c*o*v+l*v+d;
  klcd[nklcd++].val = val;
  if (pq!=rs){
     klcd[nklcd].ind   = l*o*v*v+d*o*v+k*v+c;
     klcd[nklcd++].val = val;
  }
}
void ijkl_terms(double val,ULI pq,ULI rs,ULI p,ULI q,ULI r,ULI s,ULI o,ULI&nijkl,struct integral*ijkl){
  if (p==q){
     if (r==s){
        ijkl[nijkl].ind   = p*o*o*o+r*o*o+q*o+s;
        ijkl[nijkl++].val = val;
        if (pq!=rs){
           ijkl[nijkl].ind   = r*o*o*o+p*o*o+s*o+q;
           ijkl[nijkl++].val = val;
        }
     }else{
        ijkl[nijkl].ind   = p*o*o*o+r*o*o+q*o+s;
        ijkl[nijkl++].val = val;
        ijkl[nijkl].ind   = p*o*o*o+s*o*o+q*o+r;
        ijkl[nijkl++].val = val;
        if (pq!=rs){
           ijkl[nijkl].ind   = r*o*o*o+p*o*o+s*o+q;
           ijkl[nijkl++].val = val;
           ijkl[nijkl].ind   = s*o*o*o+p*o*o+r*o+q;
           ijkl[nijkl++].val = val;
        }
     }
  }else{
     if (r==s){
        ijkl[nijkl].ind   = p*o*o*o+r*o*o+q*o+s;
        ijkl[nijkl++].val = val;
        ijkl[nijkl].ind   = q*o*o*o+r*o*o+p*o+s;
        ijkl[nijkl++].val = val;
        if (pq!=rs){
           ijkl[nijkl].ind   = r*o*o*o+p*o*o+s*o+q;
           ijkl[nijkl++].val = val;
           ijkl[nijkl].ind   = r*o*o*o+q*o*o+s*o+p;
           ijkl[nijkl++].val = val;
        }
     }else{
        ijkl[nijkl].ind   = p*o*o*o+r*o*o+q*o+s;
        ijkl[nijkl++].val = val;
        ijkl[nijkl].ind   = q*o*o*o+r*o*o+p*o+s;
        ijkl[nijkl++].val = val;
        ijkl[nijkl].ind   = p*o*o*o+s*o*o+q*o+r;
        ijkl[nijkl++].val = val;
        ijkl[nijkl].ind   = q*o*o*o+s*o*o+p*o+r;
        ijkl[nijkl++].val = val;
        if (pq!=rs){
           ijkl[nijkl].ind   = r*o*o*o+p*o*o+s*o+q;
           ijkl[nijkl++].val = val;
           ijkl[nijkl].ind   = r*o*o*o+q*o*o+s*o+p;
           ijkl[nijkl++].val = val;
           ijkl[nijkl].ind   = s*o*o*o+p*o*o+r*o+q;
           ijkl[nijkl++].val = val;
           ijkl[nijkl].ind   = s*o*o*o+q*o*o+r*o+p;
           ijkl[nijkl++].val = val;
        }
     }
  }
}


void akjc_terms(double val,ULI p,ULI q,ULI r,ULI s,ULI o,ULI v,ULI&nakjc,struct integral*akjc){
  ULI a,k,j,c;
  if (p>=o){
     a=p-o;
     c=q-o;
     j=r;
     k=s;
  }else{
     j=p;
     k=q;
     a=r-o;
     c=s-o;
  }
  if (j==k){
     if (a==c){
        akjc[nakjc].ind   = k*o*v*v + c*o*v + j*v + a;
        akjc[nakjc++].val = val;
     }
     else{
        akjc[nakjc].ind   = k*o*v*v + c*o*v + j*v + a;
        akjc[nakjc++].val = val;
        akjc[nakjc].ind   = k*o*v*v + a*o*v + j*v + c;
        akjc[nakjc++].val = val;
     }
  }
  else{
     if (a==c){
        akjc[nakjc].ind   = k*o*v*v + c*o*v + j*v + a;
        akjc[nakjc++].val = val;
        akjc[nakjc].ind   = j*o*v*v + c*o*v + k*v + a;
        akjc[nakjc++].val = val;
     }
     else{
        akjc[nakjc].ind   = k*o*v*v + c*o*v + j*v + a;
        akjc[nakjc++].val = val;
        akjc[nakjc].ind   = j*o*v*v + c*o*v + k*v + a;
        akjc[nakjc++].val = val;
        akjc[nakjc].ind   = k*o*v*v + a*o*v + j*v + c;
        akjc[nakjc++].val = val;
        akjc[nakjc].ind   = j*o*v*v + a*o*v + k*v + c;
        akjc[nakjc++].val = val;
     }
  }
}
void ijak2_terms(double val,ULI p,ULI q,ULI r,ULI s,ULI o,ULI v,ULI&nijak2,struct integral*ijak2){
  ULI i,j,a,k;
  if (p>=o){
     a=p-o;
     i=q;
     j=r;
     k=s;
  }else if (q>=o){
     i=p;
     a=q-o;
     j=r;
     k=s;
  }else if (r>=o){
     a=r-o;
     i=s;
     j=p;
     k=q;
  }else if (s>=o){
     i=r;
     a=s-o;
     j=p;
     k=q;
  }
  ijak2[nijak2].ind   = j*o*o*v + a*o*o + k*o + i;
  ijak2[nijak2++].val = val;
  if (k!=j){
     ijak2[nijak2].ind   = k*o*o*v + a*o*o + j*o + i;
     ijak2[nijak2++].val = val;
  }
}
void ijak_terms(double val,ULI p,ULI q,ULI r,ULI s,ULI o,ULI v,ULI&nijak,struct integral*ijak){
  ULI i,j,a,k;
  if (p>=o){
     a=p-o;
     i=q;
     j=r;
     k=s;
  }else if (q>=o){
     i=p;
     a=q-o;
     j=r;
     k=s;
  }else if (r>=o){
     a=r-o;
     i=s;
     j=p;
     k=q;
  }else if (s>=o){
     i=r;
     a=s-o;
     j=p;
     k=q;
  }
  ijak[nijak].ind   = j*o*o*v + i*o*v + k*v + a;
  ijak[nijak++].val = val;
  if (k!=j){
     ijak[nijak].ind   = k*o*o*v + i*o*v + j*v + a;
     ijak[nijak++].val = val;
  }
}
void abci5_terms(double val,ULI p,ULI q,ULI r,ULI s,ULI o,ULI v,ULI&nabci5,struct integral*abci5){
  ULI i,a,b,c;
  if (p<o){
     i=p;
     b=q-o;
     a=r-o;
     c=s-o;
  }else if (q<o){
     i=q;
     b=p-o;
     a=r-o;
     c=s-o;
  }else if (r<o){
     i=r;
     b=s-o;
     a=p-o;
     c=q-o;
  }else if (s<o){
     i=s;
     b=r-o;
     a=p-o;
     c=q-o;
  }
  abci5[nabci5].ind   = a*v*v*o + b*v*o + i*v + c;
  abci5[nabci5++].val = val;
  if (a!=c){
     abci5[nabci5].ind   = c*v*v*o + b*v*o + i*v + a;
     abci5[nabci5++].val = val;
  }
}
void abci3_terms(double val,ULI p,ULI q,ULI r,ULI s,ULI o,ULI v,ULI&nabci3,struct integral*abci3){
  ULI a,f,m,e;
  if (p<o){
     m=p;
     e=q-o;
     a=r-o;
     f=s-o;
  }else if (q<o){
     m=q;
     e=p-o;
     a=r-o;
     f=s-o;
  }else if (r<o){
     m=r;
     e=s-o;
     a=p-o;
     f=q-o;
  }else if (s<o){
     m=s;
     e=r-o;
     a=p-o;
     f=q-o;
  }
  abci3[nabci3].ind   = a*v*v*o + f*v*o + m*v + e;
  abci3[nabci3++].val = val;
  if (a!=f){
     abci3[nabci3].ind   = f*v*v*o + a*v*o + m*v + e;
     abci3[nabci3++].val = val;
  }
}
void abci1_terms(double val,ULI p,ULI q,ULI r,ULI s,ULI o,ULI v,ULI&nabci1,struct integral*abci1){
  ULI i,a,b,c;
  if (p<o){
     i=p;
     b=q-o;
     a=r-o;
     c=s-o;
  }else if (q<o){
     i=q;
     b=p-o;
     a=r-o;
     c=s-o;
  }else if (r<o){
     i=r;
     b=s-o;
     a=p-o;
     c=q-o;
  }else if (s<o){
     i=s;
     b=r-o;
     a=p-o;
     c=q-o;
  }
  abci1[nabci1].ind   = i*v*v*v + a*v*v + b*v + c;
  abci1[nabci1++].val = val;
  if (a!=c){
     abci1[nabci1].ind   = i*v*v*v + c*v*v + b*v + a;
     abci1[nabci1++].val = val;
  }
}
/**
  * ABCD-type integrals, because of weird SJS packing, are really 
  * confusing to sort.  I couldn't think of an analytic way to do
  * this, so I resorted to brute force.
  */
void abcd2_terms(double val,ULI pq,ULI rs,ULI p,ULI q,ULI r,ULI s,ULI o,ULI v,ULI&nabcd2,struct integral*abcd2){
  ULI ind3,a,b,c,d,ind1,ind2,index,flag;
  ULI nvals,vals[16];
  nvals = 0;
  
  a = p-o;
  d = q-o;
  b = r-o;
  c = s-o;
  ind1 = Position(a,b);
  ind2 = Position(c,d);
  if ((b<=a && d<=c) || (a<=b && c<=d)){
     ind3 = ind1*v*(v+1)/2+ind2;
     flag=1;
     for (index=0; index<nvals; index++){
         if (vals[index]==ind3){
            flag=0;
            break;
         }
     }
     if (flag){
        abcd2[nabcd2].ind   = ind3;
        abcd2[nabcd2++].val = val;
        vals[nvals++] = ind3;
     }
     if (ind1!=ind2){
        ind3 = ind2*v*(v+1)/2+ind1;
        flag=1;
        for (index=0; index<nvals; index++){
            if (vals[index]==ind3){
               flag=0;
               break;
            }
        }
        if (flag){
           abcd2[nabcd2].ind   = ind3;
           abcd2[nabcd2++].val = val;
           vals[nvals++] = ind3;
        }
     }
  }
  d = p-o;
  a = q-o;
  b = r-o;
  c = s-o;
  ind1 = Position(a,b);
  ind2 = Position(c,d);
  if ((b<=a && d<=c) || (a<=b && c<=d)){
     ind3 = ind1*v*(v+1)/2+ind2;
     flag=1;
     for (index=0; index<nvals; index++){
         if (vals[index]==ind3){
            flag=0;
            break;
         }
     }
     if (flag){
        abcd2[nabcd2].ind   = ind3;
        abcd2[nabcd2++].val = val;
        vals[nvals++] = ind3;
     }
     if (ind1!=ind2){
        ind3 = ind2*v*(v+1)/2+ind1;
        flag=1;
        for (index=0; index<nvals; index++){
            if (vals[index]==ind3){
               flag=0;
               break;
            }
        }
        if (flag){
           abcd2[nabcd2].ind   = ind3;
           abcd2[nabcd2++].val = val;
           vals[nvals++] = ind3;
        }
     }
  }
  a = p-o;
  d = q-o;
  c = r-o;
  b = s-o;
  ind1 = Position(a,b);
  ind2 = Position(c,d);
  if ((b<=a && d<=c) || (a<=b && c<=d)){
     ind3 = ind1*v*(v+1)/2+ind2;
     flag=1;
     for (index=0; index<nvals; index++){
         if (vals[index]==ind3){
            flag=0;
            break;
         }
     }
     if (flag){
        abcd2[nabcd2].ind   = ind3;
        abcd2[nabcd2++].val = val;
        vals[nvals++] = ind3;
     }
     if (ind1!=ind2){
        ind3 = ind2*v*(v+1)/2+ind1;
        flag=1;
        for (index=0; index<nvals; index++){
            if (vals[index]==ind3){
               flag=0;
               break;
            }
        }
        if (flag){
           abcd2[nabcd2].ind   = ind3;
           abcd2[nabcd2++].val = val;
           vals[nvals++] = ind3;
        }
     }
  }
  d = p-o;
  a = q-o;
  c = r-o;
  b = s-o;
  ind1 = Position(a,b);
  ind2 = Position(c,d);
  if ((b<=a && d<=c) || (a<=b && c<=d)){
     ind3 = ind1*v*(v+1)/2+ind2;
     flag=1;
     for (index=0; index<nvals; index++){
         if (vals[index]==ind3){
            flag=0;
            break;
         }
     }
     if (flag){
        abcd2[nabcd2].ind   = ind3;
        abcd2[nabcd2++].val = val;
        vals[nvals++] = ind3;
     }
     if (ind1!=ind2){
        ind3 = ind2*v*(v+1)/2+ind1;
        flag=1;
        for (index=0; index<nvals; index++){
            if (vals[index]==ind3){
               flag=0;
               break;
            }
        }
        if (flag){
           abcd2[nabcd2].ind   = ind3;
           abcd2[nabcd2++].val = val;
           vals[nvals++] = ind3;
        }
     }
  }
}
void abcd3_terms(double val,ULI pq,ULI rs,ULI p,ULI q,ULI r,ULI s,ULI o,ULI v,ULI&nabcd3,struct integral*abcd3){
  p -= o;
  q -= o;
  r -= o;
  s -= o;
  if (p==q){
     if (r==s){
        abcd3[nabcd3].ind   = p*v*v*v+r*v*v+q*v+s;
        abcd3[nabcd3++].val = val;
        if (pq!=rs){
           abcd3[nabcd3].ind   = r*v*v*v+p*v*v+s*v+q;
           abcd3[nabcd3++].val = val;
        }
     }else{
        abcd3[nabcd3].ind   = p*v*v*v+r*v*v+q*v+s;
        abcd3[nabcd3++].val = val;
        abcd3[nabcd3].ind   = p*v*v*v+s*v*v+q*v+r;
        abcd3[nabcd3++].val = val;
        if (pq!=rs){
           abcd3[nabcd3].ind   = r*v*v*v+p*v*v+s*v+q;
           abcd3[nabcd3++].val = val;
           abcd3[nabcd3].ind   = s*v*v*v+p*v*v+r*v+q;
           abcd3[nabcd3++].val = val;
        }
     }
  }else{
     if (r==s){
        abcd3[nabcd3].ind   = p*v*v*v+r*v*v+q*v+s;
        abcd3[nabcd3++].val = val;
        abcd3[nabcd3].ind   = q*v*v*v+r*v*v+p*v+s;
        abcd3[nabcd3++].val = val;
        if (pq!=rs){
           abcd3[nabcd3].ind   = r*v*v*v+p*v*v+s*v+q;
           abcd3[nabcd3++].val = val;
           abcd3[nabcd3].ind   = r*v*v*v+q*v*v+s*v+p;
           abcd3[nabcd3++].val = val;
        }
     }else{
        abcd3[nabcd3].ind   = p*v*v*v+r*v*v+q*v+s;
        abcd3[nabcd3++].val = val;
        abcd3[nabcd3].ind   = q*v*v*v+r*v*v+p*v+s;
        abcd3[nabcd3++].val = val;
        abcd3[nabcd3].ind   = p*v*v*v+s*v*v+q*v+r;
        abcd3[nabcd3++].val = val;
        abcd3[nabcd3].ind   = q*v*v*v+s*v*v+p*v+r;
        abcd3[nabcd3++].val = val;
        if (pq!=rs){
           abcd3[nabcd3].ind   = r*v*v*v+p*v*v+s*v+q;
           abcd3[nabcd3++].val = val;
           abcd3[nabcd3].ind   = r*v*v*v+q*v*v+s*v+p;
           abcd3[nabcd3++].val = val;
           abcd3[nabcd3].ind   = s*v*v*v+p*v*v+r*v+q;
           abcd3[nabcd3++].val = val;
           abcd3[nabcd3].ind   = s*v*v*v+q*v*v+r*v+p;
           abcd3[nabcd3++].val = val;
        }
     }
  }
}
void abcd1_terms(double val,ULI pq,ULI rs,ULI p,ULI q,ULI r,ULI s,ULI o,ULI v,ULI&nabcd1,struct integral*abcd1){
  ULI ind3,a,b,c,d,ind1,ind2,index,flag;
  ULI nvals,vals[16];
  nvals = 0;
  
  a = p-o;
  c = q-o;
  b = r-o;
  d = s-o;
  ind1 = Position(a,b);
  ind2 = Position(c,d);
  if ((b<=a && d<=c) || (a<=b && c<=d)){
     ind3 = ind1*v*(v+1)/2+ind2;
     flag=1;
     for (index=0; index<nvals; index++){
         if (vals[index]==ind3){
            flag=0;
            break;
         }
     }
     if (flag){
        abcd1[nabcd1].ind   = ind3;
        abcd1[nabcd1++].val = val;
        vals[nvals++] = ind3;
     }
     if (ind1!=ind2){
        ind3 = ind2*v*(v+1)/2+ind1;
        flag=1;
        for (index=0; index<nvals; index++){
            if (vals[index]==ind3){
               flag=0;
               break;
            }
        }
        if (flag){
           abcd1[nabcd1].ind   = ind3;
           abcd1[nabcd1++].val = val;
           vals[nvals++] = ind3;
        }
     }
  }
  c = p-o;
  a = q-o;
  b = r-o;
  d = s-o;
  ind1 = Position(a,b);
  ind2 = Position(c,d);
  if ((b<=a && d<=c) || (a<=b && c<=d)){
     ind3 = ind1*v*(v+1)/2+ind2;
     flag=1;
     for (index=0; index<nvals; index++){
         if (vals[index]==ind3){
            flag=0;
            break;
         }
     }
     if (flag){
        abcd1[nabcd1].ind   = ind3;
        abcd1[nabcd1++].val = val;
        vals[nvals++] = ind3;
     }
     if (ind1!=ind2){
        ind3 = ind2*v*(v+1)/2+ind1;
        flag=1;
        for (index=0; index<nvals; index++){
            if (vals[index]==ind3){
               flag=0;
               break;
            }
        }
        if (flag){
           abcd1[nabcd1].ind   = ind3;
           abcd1[nabcd1++].val = val;
           vals[nvals++] = ind3;
        }
     }
  }
  a = p-o;
  c = q-o;
  d = r-o;
  b = s-o;
  ind1 = Position(a,b);
  ind2 = Position(c,d);
  if ((b<=a && d<=c) || (a<=b && c<=d)){
     ind3 = ind1*v*(v+1)/2+ind2;
     flag=1;
     for (index=0; index<nvals; index++){
         if (vals[index]==ind3){
            flag=0;
            break;
         }
     }
     if (flag){
        abcd1[nabcd1].ind   = ind3;
        abcd1[nabcd1++].val = val;
        vals[nvals++] = ind3;
     }
     if (ind1!=ind2){
        ind3 = ind2*v*(v+1)/2+ind1;
        flag=1;
        for (index=0; index<nvals; index++){
            if (vals[index]==ind3){
               flag=0;
               break;
            }
        }
        if (flag){
           abcd1[nabcd1].ind   = ind3;
           abcd1[nabcd1++].val = val;
           vals[nvals++] = ind3;
        }
     }
  }
  c = p-o;
  a = q-o;
  d = r-o;
  b = s-o;
  ind1 = Position(a,b);
  ind2 = Position(c,d);
  if ((b<=a && d<=c) || (a<=b && c<=d)){
     ind3 = ind1*v*(v+1)/2+ind2;
     flag=1;
     for (index=0; index<nvals; index++){
         if (vals[index]==ind3){
            flag=0;
            break;
         }
     }
     if (flag){
        abcd1[nabcd1].ind   = ind3;
        abcd1[nabcd1++].val = val;
        vals[nvals++] = ind3;
     }
     if (ind1!=ind2){
        ind3 = ind2*v*(v+1)/2+ind1;
        flag=1;
        for (index=0; index<nvals; index++){
            if (vals[index]==ind3){
               flag=0;
               break;
            }
        }
        if (flag){
           abcd1[nabcd1].ind   = ind3;
           abcd1[nabcd1++].val = val;
           vals[nvals++] = ind3;
        }
     }
  }

}

void SortBlock(ULI nelem,ULI blockdim,struct integral*buffer,double*tmp,ULI PSIFILE,char*string,ULI maxdim){
  boost::shared_ptr<PSIO> psio(new PSIO());
  // does the block fit in core?
  if (nelem<=maxdim && blockdim<=maxdim){
     psio->open(PSIFILE,PSIO_OPEN_OLD);
     psio->read_entry(PSIFILE,string,(char*)&buffer[0],nelem*sizeof(struct integral));
     psio->close(PSIFILE,0);

     memset((void*)tmp,'\0',blockdim*sizeof(double));
     for (ULI j=0; j<nelem; j++){
         tmp[buffer[j].ind] = buffer[j].val;
     }

     psio->open(PSIFILE,PSIO_OPEN_NEW);
     psio->write_entry(PSIFILE,string,(char*)&tmp[0],blockdim*sizeof(double));
     psio->close(PSIFILE,1);
  }else{
     ULI nbins,binsize,lastbin;
     for (ULI i=1; i<=blockdim; i++){
         if (maxdim>=(double)blockdim/i){
            binsize = blockdim/i;
            if (i*binsize < blockdim) binsize++;
            nbins = i;
            break;
         }
     }
     lastbin = blockdim - (nbins-1)*binsize;

     ULI initialnbins,initialbinsize,initiallastbin;
     for (ULI i=1; i<=nelem; i++){
         if (maxdim>=(double)nelem/i){
            initialbinsize = nelem/i;
            if (i*initialbinsize < nelem) initialbinsize++;
            initialnbins = i;
            break;
         }
     }
     initiallastbin = nelem - (initialnbins-1)*initialbinsize;

     psio_address*addr,addr1,addr2;
     addr = new psio_address[nbins];
     addr1 = PSIO_ZERO;

     psio->open(PSIFILE,PSIO_OPEN_OLD);

     psio->open(PSIF_TEMP,PSIO_OPEN_NEW);
     addr2=PSIO_ZERO;
     for (ULI k=0; k<nbins; k++){
         addr1=PSIO_ZERO;
         memset((void*)tmp,'\0',binsize*sizeof(double));
         for (ULI i=0; i<initialnbins-1; i++){
             psio->read(PSIFILE,string,(char*)&buffer[0],initialbinsize*sizeof(struct integral),addr1,&addr1);
             for (ULI j=0; j<initialbinsize; j++){
                 if (buffer[j].ind < (k+1)*binsize && buffer[j].ind>=k*binsize){
                     tmp[buffer[j].ind-k*binsize] = buffer[j].val;
                 }
             }
         }
         psio->read(PSIFILE,string,(char*)&buffer[0],initiallastbin*sizeof(struct integral),addr1,&addr1);
         for (ULI j=0; j<initiallastbin; j++){
             if (buffer[j].ind < (k+1)*binsize && buffer[j].ind>=k*binsize){
                 tmp[buffer[j].ind-k*binsize] = buffer[j].val;
             }
         }
         psio->write(PSIF_TEMP,string,(char*)&tmp[0],binsize*sizeof(double),addr2,&addr2);
     }
     psio->close(PSIFILE,1);
     psio->close(PSIF_TEMP,1);

     psio->rename_file(PSIF_TEMP,PSIFILE);

     delete addr;
  }
}

} // end of namespace
  
