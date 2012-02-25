#include"psi4-dec.h"
#include<psifiles.h>
#include"blas.h"
#include"cim.h"
#include<libmints/wavefunction.h>
#include<libmints/matrix.h>

using namespace psi;
using namespace boost;

namespace psi{

/*
 * build central domains based on Fock matrix elements
 */
void CIM::OccupiedDomains(){

  fprintf(outfile,"\n");
  fprintf(outfile,"  ==> Define Clusters of Interacting Occupied Orbitals <==\n");
  fprintf(outfile,"\n");

  // threshold for central domains
  thresh1 = options_.get_double("THRESH1");

  int**pdomains = (int**)malloc(ndoccact*sizeof(int*));
  int*pdomainsizes = (int*)malloc(ndoccact*sizeof(int));
  for (int i=0; i<ndoccact; i++){
      pdomains[i] = (int*)malloc(ndoccact*sizeof(int));
      pdomainsizes[i] = 0;
      for (int j=0; j<ndoccact; j++){
          pdomains[i][j] = -999;
      }
  }
  for (int i=0; i<ndoccact; i++){
      for (int j=0; j<ndoccact; j++){
          if (fabs(Fock[i+nfzc][j+nfzc]) >= thresh1){
             pdomains[i][j] = j;
             pdomainsizes[i]++;
          }
      }
  }

  // collapse redundant domains:
 
  central = (int**)malloc(ndoccact*sizeof(int*));
  ncentral = (int*)malloc(ndoccact*sizeof(int));
  for (int i=0; i<ndoccact; i++){
      central[i] = (int*)malloc(ndoccact*sizeof(int));
      ncentral[i] = 1;
      central[i][0] = i;
  }

  int ndomains = ndoccact;
  int *skip = (int*)malloc(ndoccact*sizeof(int));
  for (int i=0; i<ndoccact; i++) skip[i] = 0;
  for (int i=0; i<ndoccact; i++){
      if (skip[i]) continue;
      // check if domain i fits in domain j
      for (int j=0; j<ndoccact; j++){
          if (i==j)    continue;
          if (skip[j]) continue;
          int count = 0;
          for (int k=0; k<ndoccact; k++){
              int dum  = pdomains[i][k] - pdomains[j][k];
              if (dum==0 && pdomains[i][k]!=-999) count++;
          }
          if (count==pdomainsizes[i]){
             skip[i] = 1;
             central[j][ncentral[j]++] = i;
             ndomains--;
          }
      }
  }

  // remove central orbitals from central mo domain:
  for (int i=0; i<ndoccact; i++){
      if (skip[i]) continue;
      for (int j=0; j<ncentral[i]; j++){
          pdomains[i][central[i][j]] = -999;
      }
  }

  // threshold for environmental domains
  env = (int**)malloc(ndoccact*sizeof(int*));
  nenv = (int*)malloc(ndoccact*sizeof(int));
  for (int i=0; i<ndoccact; i++){
      env[i] = (int*)malloc(ndoccact*sizeof(int));
      nenv[i] = 0;
      for (int j=0; j<ndoccact; j++){
          env[i][j] = -999;
      }
  }
  thresh2 = options_.get_double("THRESH2");
  // loop over unique clusters
  for (int i=0; i<ndoccact; i++){
      if (skip[i]) continue;
      // loop over all occupied orbitals
      for (int j=0; j<ndoccact; j++){
          int isenv = 0;
          int redundant = 0;
          // check interaction of orbital j with central orbitals in cluster i
          for (int k=0; k<ncentral[i]; k++){
              if (central[i][k]==j){
                 redundant=1;
                 break;
              }
              if (fabs(Fock[central[i][k]+nfzc][j+nfzc]) >= thresh2){
                 isenv=1;
              }
          }
          if (redundant) continue;
          // check interaction of orbital j with central mo domain in cluster i
          for (int k=0; k<ndoccact; k++){
              if (pdomains[i][k]==-999) continue;
              if (pdomains[i][k]==j){
                 redundant=1;
                 break;
              }
              if (fabs(Fock[pdomains[i][k]+nfzc][j+nfzc]) >= thresh2){
                 isenv=1;
              }
          }
          if (redundant) continue;
          // add orbital to environmental domain
          if (isenv)
             env[i][j] = j;
      }
  }

  // print out clusters
  for (int i=0; i<ndoccact; i++){
      if (skip[i]) continue;
      fprintf(outfile,"  Cluster %3i:\n",i);
      int count=0;
      int first=1;
      for (int j=0; j<ncentral[i]; j++){
          count++;
          if (count==1){
             if (first){
                first=0;
                fprintf(outfile,"     Central orbitals:      ");
             }else{
                fprintf(outfile,"                            ");
             }
          }
          fprintf(outfile," %3i",central[i][j]);
          if (count==8){
             fprintf(outfile,"\n");
             count=0;
          }
      }
      fprintf(outfile,"\n");
      count=0;
      first=1;
      for (int j=0; j<ndoccact; j++){
          if (pdomains[i][j]!=-999){
             count++;
             if (count==1){
                if (first){
                   first=0;
                   fprintf(outfile,"     Central domain:        ");
                }else{
                   fprintf(outfile,"                            ");
                }
             }
             fprintf(outfile," %3i",pdomains[i][j]);
          }
          if (count==8){
             fprintf(outfile,"\n");
             count=0;
          }
      }
      fprintf(outfile,"\n");
      count=0;
      first=1;
      for (int j=0; j<ndoccact; j++){
          if (env[i][j]!=-999){
             count++;
             if (count==1){
                if (first){
                   first=0;
                   fprintf(outfile,"     Environmental domain:  ");
                }else{
                   fprintf(outfile,"                            ");
                }
             }
             fprintf(outfile," %3i",env[i][j]);
          }
          if (count==8){
             fprintf(outfile,"\n");
             count=0;
          }
      }
      fprintf(outfile,"\n");
      fprintf(outfile,"\n");
  }
  fflush(outfile);


  // TODO
  /*fprintf(outfile,"\n");
  fprintf(outfile,"  Collapse Redundant Domains:\n");
  fprintf(outfile,"\n");

  int**tempdomains=(int**)malloc(ndoccact*sizeof(double*));
  int*tempdomainsizes=(int*)malloc(ndoccact*sizeof(double));
  for (int i=0; i<ndoccact; i++){
      tempdomains[i] = (int*)malloc(ndoccact*sizeof(int));
      tempdomainsizes[i] = 0;
      for (int j=0; j<ndoccact; j++){
          tempdomains[i][j] = -999;
      }
  }
  for (int i=0; i<ndoccact; i++){
      if (skip[i]) continue;
      for (int j=0; j<ncentral[i]; j++){
          tempdomains[i][central[i][j]] = central[i][j];
          tempdomainsizes[i]++;
      }
      for (int j=0; j<ndoccact; j++){
          if (pdomains[i][j]==-999) continue;
          tempdomains[i][j] = pdomains[i][j];
          tempdomainsizes[i]++;
      }
      for (int j=0; j<ndoccact; j++){
          if (env[i][j]==-999) continue;
          tempdomains[i][j] = env[i][j];
          tempdomainsizes[i]++;
      }
  }
  
  ndomains = 0;
  for (int i=0; i<ndoccact; i++){
      if (!skip[i]) ndomains++;
  }
  printf("ndomains %5i\n",ndomains);

  for (int i=0; i<ndoccact; i++){
      if (skip[i]) continue;
      // check if domain i fits in domain j
      for (int j=0; j<ndoccact; j++){
          if (i==j)    continue;
          if (skip[j]) continue;
          int count = 0;
          for (int k=0; k<ndoccact; k++){
              int dum  = tempdomains[i][k] - tempdomains[j][k];
              if (dum==0 && tempdomains[i][k]!=-999) count++;
          }
          if (count==tempdomainsizes[i]){
             skip[i] = 1;
             for (int k=0; k<ncentral[i]; k++){
                 // is the orbital already a central orbital?
                 int redundant=0;
                 for (int l=0; l<ncentral[j]; l++){
                     if (central[i][k]==central[j][l]){
                        redundant=1;
                        break;
                     }
                 }
                 if (redundant) continue;
                 // add to central orbitals of cluster j
                 central[j][ncentral[j]]    = central[i][k];
                 pdomains[j][central[i][k]] = -999;
                 env[j][central[i][k]]      = -999;
                 ncentral[j]++;
             }
          }
      }
  }
  ndomains = 0;
  for (int i=0; i<ndoccact; i++){
      if (!skip[i]) ndomains++;
  }
  printf("ndomains %5i\n",ndomains);


  fprintf(outfile,"\n");
  fprintf(outfile,"  Full Occupied Domains:\n");
  fprintf(outfile,"\n");

  for (int i=0; i<ndoccact; i++){
      if (skip[i]) continue;
      fprintf(outfile,"  Cluster %3i:\n",i);
      int count=0;
      int first=1;
      for (int j=0; j<ncentral[i]; j++){
          count++;
          if (count==1){
             if (first){
                first=0;
                fprintf(outfile,"     Central orbitals:      ");
             }else{
                fprintf(outfile,"                            ");
             }
          }
          fprintf(outfile," %3i",central[i][j]);
          if (count==8){
             fprintf(outfile,"\n");
             count=0;
          }
      }
      fprintf(outfile,"\n");
      count=0;
      first=1;
      for (int j=0; j<ndoccact; j++){
          if (pdomains[i][j]!=-999){
             count++;
             if (count==1){
                if (first){
                   first=0;
                   fprintf(outfile,"     Central domain:        ");
                }else{
                   fprintf(outfile,"                            ");
                }
             }
             fprintf(outfile," %3i",pdomains[i][j]);
          }
          if (count==8){
             fprintf(outfile,"\n");
             count=0;
          }
      }
      fprintf(outfile,"\n");
      count=0;
      first=1;
      for (int j=0; j<ndoccact; j++){
          if (env[i][j]!=-999){
             count++;
             if (count==1){
                if (first){
                   first=0;
                   fprintf(outfile,"     Environmental domain:  ");
                }else{
                   fprintf(outfile,"                            ");
                }
             }
             fprintf(outfile," %3i",env[i][j]);
          }
          if (count==8){
             fprintf(outfile,"\n");
             count=0;
          }
      }
      fprintf(outfile,"\n");
      fprintf(outfile,"\n");
  }
  fflush(outfile);*/

}

}

