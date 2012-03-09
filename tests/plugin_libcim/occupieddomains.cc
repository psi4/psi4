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
 * this algorithm is called dual-environment cim in
 * J. Chem. Phys. 135, 104111 (2011)
 */
void CIM::OccupiedDomains(){

  // all orbitals:         domain,domainsize
  // central orbitals:     central, ncentral
  // central mo domain:    modomain, nmodomain
  // environmental domain: env, nenv

  fprintf(outfile,"  ==> Define Clusters of Interacting Occupied Orbitals <==\n");
  fprintf(outfile,"\n");

  // threshold for central domains
  thresh1 = options_.get_double("THRESH1");

  domain = (int**)malloc(ndoccact*sizeof(int*));
  domainsize = (int*)malloc(ndoccact*sizeof(int));
  for (int i=0; i<ndoccact; i++){
      domain[i] = (int*)malloc(ndoccact*sizeof(int));
      domainsize[i] = 0;
      for (int j=0; j<ndoccact; j++){
          domain[i][j] = isempty;
      }
  }
  for (int i=0; i<ndoccact; i++){
      for (int j=0; j<ndoccact; j++){
          if (fabs(Fock[i+nfzc][j+nfzc]) >= thresh1){
             domain[i][j] = j;
             domainsize[i]++;
          }
      }
  }

  // list of central orbitals 
  central = (int**)malloc(ndoccact*sizeof(int*));
  ncentral = (int*)malloc(ndoccact*sizeof(int));
  for (int i=0; i<ndoccact; i++){
      central[i] = (int*)malloc(ndoccact*sizeof(int));
      for (int j=0; j<ndoccact; j++){
          central[i][j] = isempty;
      }
      ncentral[i] = 1;
      central[i][i] = i;
  }
  ndomains = ndoccact;

  // list of central mo domains
  modomain = (int**)malloc(ndoccact*sizeof(int*));
  nmodomain = (int*)malloc(ndoccact*sizeof(int));
  for (int i=0; i<ndoccact; i++){
      modomain[i] = (int*)malloc(ndoccact*sizeof(int));
      nmodomain[i] = 0;
  }
  for (int i=0; i<ndoccact; i++){
      for (int j=0; j<ndoccact; j++){
          modomain[i][j] = domain[i][j];
      }
  }
  // remove central orbitals from central mo domain:
  for (int i=0; i<ndoccact; i++){
      modomain[i][i] = isempty;
      nmodomain[i] = 0;
      for (int j=0; j<ndoccact; j++){
          if (modomain[i][j]!=isempty) nmodomain[i]++;
      }
      domainsize[i]=0;
      for (int j=0; j<ndoccact; j++){
          if (domain[i][j]!=isempty) domainsize[i]++;
      }
  }


  // collapse redundant domains:
  skip = (int*)malloc(ndoccact*sizeof(int));
  for (int i=0; i<ndoccact; i++) skip[i] = 0;
  for (int i=0; i<ndoccact; i++){
      if (skip[i]) continue;
      // check if domain i fits in domain j
      for (int j=0; j<ndoccact; j++){
          if (i==j)    continue;
          if (skip[j]) continue;
          int count = 0;
          for (int k=0; k<ndoccact; k++){
              if (domain[i][k]==isempty) continue;
              if (domain[i][k]==domain[j][k]) count++;
          }
          // add central orbital from {i} to {j}
          // also remove that orbital from the central mo domain of {j}
          if (count==domainsize[i]){
             skip[i] = 1;
             central[j][i] = i;
             ncentral[j]++;
             modomain[j][i] = isempty;
             ndomains--;
          }
      }
  }

  // threshold for environmental domains
  env = (int**)malloc(ndoccact*sizeof(int*));
  nenv = (int*)malloc(ndoccact*sizeof(int));
  for (int i=0; i<ndoccact; i++){
      env[i] = (int*)malloc(ndoccact*sizeof(int));
      nenv[i] = 0;
      for (int j=0; j<ndoccact; j++){
          env[i][j] = isempty;
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
          // check interaction of orbital j with all orbitals in cluster {i}
          for (int k=0; k<ndoccact; k++){
              if (domain[i][k]==isempty) continue;
              if (k==j){
                 redundant=1;
                 break;
              }
              if (fabs(Fock[k+nfzc][j+nfzc]) >= thresh2){
                 isenv=1;
              }
          }
          if (redundant) continue;
          // add orbital to environmental domain
          if (isenv){
             env[i][j] = j;
          }
      }
  }
  for (int i=0; i<ndoccact; i++){
      if (skip[i]) continue;
      domainsize[i]=0;
      for (int j=0; j<ndoccact; j++){
          if (env[i][j]!=isempty) domain[i][j] = j;
          if (domain[i][j]!=isempty) domainsize[i]++;
      }
  }

  // collapse again
  for (int i=0; i<ndoccact; i++){     // {I}
      if (skip[i]) continue;
      int tossi=0;
      for (int j=0; j<ndoccact; j++){ // {J}
          if (i==j) continue;
          if (skip[j]) continue;
          int count = 0;
          for (int k=0; k<ndoccact; k++){
              if (domain[i][k]==isempty) continue;
              if (domain[i][k]==domain[j][k]) count++;
          }
          // if {I} fits in {J}, add centrals of {I} to {J} and remove
          // those orbitals from the central mo and environmental domains
          if (count==domainsize[i]){
             for (int k=0; k<ndoccact; k++){
                 if (central[i][k]==isempty) continue;
                 central[j][k]  = k;
                 modomain[j][k] = isempty;
                 env[j][k]      = isempty;
             }
             tossi=1;
          }
          
      }
      if (tossi) skip[i] = 1;
  }

  // count orbitals in each class:
  double cccost  = 0.0;
  double tcost   = 0.0;
  ndomains = 0;
  for (int i=0; i<ndoccact; i++){
      if (skip[i]) continue;
      ncentral[i] = 0;
      nmodomain[i] = 0;
      nenv[i] = 0;
      ndomains++;
      for (int j=0; j<ndoccact; j++){
          if (central[i][j]!=isempty)  ncentral[i]++;
          if (modomain[i][j]!=isempty) nmodomain[i]++;
          if (env[i][j]!=isempty)      nenv[i]++;
      }
      double dum = (double)domainsize[i]/ndoccact;
      cccost  += pow(dum,6);
      tcost   += pow(dum,7);
  }

  fprintf(outfile,"  Cluster |");
  fprintf(outfile," Central Orbitals     |");
  fprintf(outfile," Central MO Domain    |");
  fprintf(outfile," Environmental Domain\n");
  fprintf(outfile,"  ----------------------------------------------------------------------------\n");
  // print clusters:
  int clusternum=-1;
  for (int i=0; i<ndoccact; i++){
      if (skip[i]) continue;
      clusternum++;

      int nrowscentral  = ncentral[i]  / 5;
      int nrowsmodomain = nmodomain[i] / 5;
      int nrowsenv      = nenv[i]      / 5;
      if (ncentral[i]%5>0)  nrowscentral++;
      if (nmodomain[i]%5>0) nrowsmodomain++;
      if (nenv[i]%5>0)      nrowsenv++;
      int nrows = nrowscentral;
      if (nrowsmodomain>nrows) nrows = nrowscentral;
      if (nrowsmodomain>nrows) nrows = nrowsmodomain;
      if (nrowsenv>nrows)      nrows = nrowsenv;
      int nspacescentral  = nrowscentral  * 5 - ncentral[i];
      int nspacesmodomain = nrowsmodomain * 5 - nmodomain[i];
      int nspacesenv      = nrowsenv      * 5 - nenv[i];

      int lastcentral=-1;
      int lastmodomain=-1;
      int lastenv=-1;

      
      for (int row = 0; row<nrows; row++){


          // cluster number
          if (row==0){
             fprintf(outfile,"      %3i |",clusternum);
          }
          else{
             fprintf(outfile,"          |");
          }

          // central orbitals
          if (nrowscentral<row+1){
             fprintf(outfile,"                      |");
          }
          else{
             if (row<nrowscentral-1){
                int count=0;
                for (int j=lastcentral+1; j<ndoccact; j++){
                    if (central[i][j]!=isempty){
                       fprintf(outfile," %3i",j);
                       count++;
                    }
                    if (count==5){
                       lastcentral = j;
                       break;
                    }
                }
             }
             else{
                for (int j=lastcentral+1; j<ndoccact; j++){
                    if (central[i][j]!=isempty){
                       fprintf(outfile," %3i",j);
                    }
                }
                for (int j=0; j<nspacescentral; j++){
                    fprintf(outfile,"    ");
                }
             }
             fprintf(outfile,"  |");
          }
          // central mo domain
          if (nrowsmodomain<row+1){
             fprintf(outfile,"                      |");
          }
          else{
             if (row<nrowsmodomain-1){
                int count=0;
                for (int j=lastmodomain+1; j<ndoccact; j++){
                    if (modomain[i][j]!=isempty){
                       fprintf(outfile," %3i",j);
                       count++;
                    }
                    if (count==5){
                       lastmodomain = j;
                       break;
                    }
                }
             }
             else{
                for (int j=lastmodomain+1; j<ndoccact; j++){
                    if (modomain[i][j]!=isempty){
                       fprintf(outfile," %3i",j);
                    }
                }
                for (int j=0; j<nspacesmodomain; j++){
                    fprintf(outfile,"    ");
                }
             }
             fprintf(outfile,"  |");
          }
          // environmental doain
          if (nrowsenv<row+1){
             fprintf(outfile,"\n");
          }
          else{
             if (row<nrowsenv-1){
                int count=0;
                for (int j=lastenv+1; j<ndoccact; j++){
                    if (env[i][j]!=isempty){
                       fprintf(outfile," %3i",j);
                       count++;
                    }
                    if (count==5){
                       lastenv = j;
                       break;
                    }
                }
             }
             else{
                for (int j=lastenv+1; j<ndoccact; j++){
                    if (env[i][j]!=isempty){
                       fprintf(outfile," %3i",j);
                    }
                }
             }
             fprintf(outfile,"\n");
          }
         
          
      }
      fprintf(outfile,"  ----------------------------------------------------------------------------\n");
  }
  fflush(outfile);
  fprintf(outfile,"\n");
  fprintf(outfile,"  ==> Estimated Cost (Relative to Canonical Calculations) <==\n");
  fprintf(outfile,"\n");
  fprintf(outfile,"      CCSD iterations:  %7.4lf\n",cccost);
  fprintf(outfile,"      (T) contribution: %7.4lf\n",tcost);
  fprintf(outfile,"\n");
}

}

