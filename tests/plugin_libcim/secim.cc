#include"psi4-dec.h"
#include<psifiles.h>
#include<physconst.h>
#include"blas.h"
#include"cim.h"
#include<libmints/wavefunction.h>
#include<libmints/mints.h>
#include<libmints/matrix.h>
#include<libmints/molecule.h>

using namespace psi;
using namespace boost;

namespace psi{

/*
 * build central domains based on Fock matrix elements
 * this algorithm is called single-environment cim in
 * J. Chem. Phys. 135, 104111 (2011)
 */
void CIM::SECIM(){

  // all orbitals:         domain,domainsize
  // central orbitals:     central, ncentral
  // environmental domain: env, nenv
  // here (as opposed to in decim), the
  // central mo domain doesn't exist

  fprintf(outfile,"  ==> Define Clusters of Interacting Occupied Orbitals <==\n");
  fprintf(outfile,"\n");

  // VDW radii int angstroms 
  // from supplemental material from J. Phys. Chem A 111, 2193 (2007)
  double * VDWRadius = (double*)malloc(100*sizeof(double));
  VDWRadius[1]  = 0.30; // H
  VDWRadius[2]  = 1.16; // He
  VDWRadius[3]  = 1.23; // Li
  VDWRadius[4]  = 0.89; // Be
  VDWRadius[5]  = 0.88; // B
  VDWRadius[6]  = 0.77; // C
  VDWRadius[7]  = 0.70; // N
  VDWRadius[8]  = 0.66; // O
  VDWRadius[9]  = 0.58; // F
  VDWRadius[10] = 0.55; // Ne
  VDWRadius[11] = 1.40; // Na
  VDWRadius[12] = 1.36; // Mg
  VDWRadius[13] = 1.25; // Al
  VDWRadius[14] = 1.17; // Si
  VDWRadius[15] = 1.10; // P
  VDWRadius[16] = 1.11; // S
  VDWRadius[17] = 0.99; // Cl
  VDWRadius[18] = 1.58; // Ar
 
  // get list of heavy atoms
  boost::shared_ptr<Molecule>molecule_ = Process::environment.molecule();
  int natom = molecule_->natom();


  // maximum possible number of domains:
  maxndomains = natom;

  int *heavyatom = (int*)malloc(natom*sizeof(int));
  for (int i=0; i<natom; i++){
      if (molecule_->fZ(i)>1.0) heavyatom[i] = 1;
      else                      heavyatom[i] = 0;
  }
  int ** atoms = (int**)malloc(natom*sizeof(int*));
  for (int i=0; i<natom; i++){
      atoms[i] = (int*)malloc(natom*sizeof(int));
      for (int j=0; j<natom; j++){
          atoms[i][j] = isempty;
      }
  }
  // list of heavy atoms and nearby hydrogens
  skip = (int*)malloc(maxndomains*sizeof(int));
  for (int i=0; i<maxndomains; i++) skip[i] = 0;
  for (int i=natom; i<maxndomains; i++) skip[i] = 1;
  for (int i=0; i<natom; i++){
      if (!heavyatom[i]){
         skip[i]=1;
         continue;
      }
      atoms[i][i] = i;
      double Xx = molecule_->x(i);
      double Xy = molecule_->y(i);
      double Xz = molecule_->z(i);
      // which hydrogens are nearby (within R(H)+R(X)+0.168 A ) ?
      for (int j=0; j<natom; j++){
          if (heavyatom[j]) continue;
          double Hx = molecule_->x(j);
          double Hy = molecule_->y(j);
          double Hz = molecule_->z(j);
          double rHX = (Xx-Hx)*(Xx-Hx) + (Xy-Hy)*(Xy-Hy) + (Xz-Hz)*(Xz-Hz);
          rHX  = sqrt(rHX);
          rHX *= _bohr2angstroms;
          if (rHX <= VDWRadius[(int)(molecule_->fZ(i))]+VDWRadius[(int)(molecule_->fZ(j))]+0.168){
             atoms[i][j] = j;
          }
      }
  }
  free(VDWRadius);

  // list of central orbitals 
  central = (int**)malloc(maxndomains*sizeof(int*));
  ncentral = (int*)malloc(maxndomains*sizeof(int));
  for (int i=0; i<maxndomains; i++){
      central[i] = (int*)malloc(ndoccact*sizeof(int));
      for (int j=0; j<ndoccact; j++){
          central[i][j] = isempty;
      }
      ncentral[i] = 0;
  }
  domain = (int**)malloc(maxndomains*sizeof(int*));
  domainsize = (int*)malloc(maxndomains*sizeof(int));
  for (int i=0; i<maxndomains; i++){
      domain[i] = (int*)malloc(ndoccact*sizeof(int));
      domainsize[i] = 0;
      for (int j=0; j<ndoccact; j++){
          domain[i][j] = isempty;
      }
  }

  // mulliken charges for each atom - use OEProp
  // need to overwrite Da in wfn ...
  boost::shared_ptr<psi::Wavefunction> 
      ref = Process::environment.reference_wavefunction();
  SharedMatrix SaveDa(ref->Da());

  double**C = boys->Clmo->pointer();
  for (int i=0; i<ndoccact; i++){
      SharedMatrix Da_so(SaveDa->clone());
      double**Da_so_pointer = Da_so->pointer();
      // Da is the density associated with LMO i
      for (int mu=0; mu<nso; mu++){
          for (int nu=0; nu<nso; nu++){
              Da_so_pointer[mu][nu] = C[mu][i+nfzc]*C[nu][i+nfzc];
          }
      }
      // compute Mulliken charges:
      boost::shared_ptr<OEProp> oe(new OEProp());
      boost::shared_ptr<Vector> Qa = 
          oe->compute_mulliken_charges_custom_Da(Da_so);
      double * Qa_pointer = Qa->pointer();

      // look at mulliken charges due to LMO i
      for (int j=0; j<natom; j++){
          if (!heavyatom[j]) continue;
          // check if orbital i is central in cluster j
          for (int k=0; k<natom; k++){
              if (atoms[j][k]==isempty) continue;
              if (2.0*Qa_pointer[k]>0.3){
                 central[j][i] = i;
                 domain[j][i] = i;
                 ncentral[j]++;
                 domainsize[j]++;
              }
          }
      }
  }

  // environmental orbitals
  env = (int**)malloc(maxndomains*sizeof(int*));
  nenv = (int*)malloc(maxndomains*sizeof(int));
  for (int i=0; i<maxndomains; i++){
      env[i] = (int*)malloc(ndoccact*sizeof(int));
      nenv[i] = 0;
      for (int j=0; j<ndoccact; j++){
          env[i][j] = isempty;
      }
  }
  // threshold for environmental orbitals
  thresh1 = options_.get_double("THRESH1");
  //thresh1 = sqrt(thresh1);
  for (int i=0; i<natom; i++){
      if (skip[i]) continue;
      for (int j=0; j<ndoccact; j++){
          if (domain[i][j]==isempty) continue;
          for (int k=0; k<ndoccact; k++){
              if (j==k) continue;
              if (fabs(Fock[j+nfzc][k+nfzc]) >= thresh1){
                 if (domain[i][k]!=isempty) continue;
                 env[i][k] = k;
                 nenv[i]++;
              }
          }
      }
  }
  for (int i=0; i<natom; i++){
      if (skip[i]) continue;
      for (int j=0; j<ndoccact; j++){
          if (env[i][j]==isempty) continue;
          domain[i][j] = env[i][j];
      }
  }


  ndomains = natom;

  // list of central mo domains (not used in secim)
  modomain = (int**)malloc(maxndomains*sizeof(int*));
  nmodomain = (int*)malloc(maxndomains*sizeof(int));
  for (int i=0; i<maxndomains; i++){
      modomain[i] = (int*)malloc(ndoccact*sizeof(int));
      nmodomain[i] = 0;
  }
  for (int i=0; i<maxndomains; i++){
      for (int j=0; j<ndoccact; j++){
          modomain[i][j] = isempty;
      }
  }


  // check domainsizes
  for (int i=0; i<natom; i++){
      domainsize[i]=0;
      for (int j=0; j<ndoccact; j++){
          if (domain[i][j]==isempty) continue;
          domainsize[i]++;
      }
  }

  // collapse redundant domains
  for (int i=0; i<natom; i++){     // {I}
      if (skip[i]) continue;
      int tossi=0;
      for (int j=0; j<natom; j++){ // {J}
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
  for (int i=0; i<maxndomains; i++){
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
  for (int i=0; i<maxndomains; i++){
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

