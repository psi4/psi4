#include<libmints/mints.h>
#include"blas.h"
#include"cepa.h"

using namespace psi;

namespace psi{
void OPDM(boost::shared_ptr<psi::CoupledPair>cepa,Options&options);
double Normalize(long int o,long int v,double*t1,double*t2,int cepa_level);
void BuildD1(long int o,long int v,double*t1,double*ta,double*tb,double c0,double*D1);

void OPDM(boost::shared_ptr<psi::CoupledPair>cepa,Options&options){

  long int o = cepa->ndoccact;
  long int v = cepa->nvirt;

  // if t2 was stored on disk, grab it.
  if (cepa->t2_on_disk){
     boost::shared_ptr<PSIO> psio(new PSIO());
     psio->open(PSIF_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_T2,"t2",(char*)&cepa->tempv[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_T2,1);
     cepa->tb = cepa->tempv;
  }

  // Normalize wave function and return leading coefficient 
  double c0 = Normalize(o,v,cepa->t1,cepa->tb,cepa->cepa_level);

  // Build 1-RDM
  double*D1 = (double*)malloc((o+v)*(o+v)*sizeof(double));
  BuildD1(o,v,cepa->t1,cepa->integrals,cepa->tb,c0,D1);

  // Call oeprop
  boost::shared_ptr<OEProp> oe(new OEProp());
  boost::shared_ptr<Wavefunction> wfn =
    Process::environment.reference_wavefunction();
  boost::shared_ptr<Matrix> Ca = wfn->Ca();

  std::stringstream ss;
  ss << cepa->cepa_type;
  std::stringstream ss_a;
  ss_a << ss.str() << " alpha";

  // pass opdm to oeprop
  SharedMatrix opdm_a(new Matrix(ss_a.str(), Ca->colspi(), Ca->colspi()));
  double** opdmap = opdm_a->pointer();
  for (long int i=0; i<cepa->nfzc; i++){
      for (long int j=0; j<cepa->nfzc; j++){
          opdmap[i][j] = 0.0;
      }
      opdmap[i][i] = 1.0;
  }
  for (long int i=0; i<o+v; i++){
      for (long int j=0; j<o+v; j++){
          opdmap[i+cepa->nfzc][j+cepa->nfzc] = D1[i*(o+v)+j];
      }
  }
  oe->set_Da_mo(opdm_a);

  std::stringstream oeprop_label;
  oeprop_label << cepa->cepa_type;
  oe->set_title(oeprop_label.str());
  oe->add("DIPOLE");
  oe->add("MULLIKEN_CHARGES");
  oe->add("NO_OCCUPATIONS");
  if (options.get_int("PRINT") > 1) {
      oe->add("QUADRUPOLE");
  }
   
  fprintf(outfile, "\n");
  fprintf(outfile, "  ==> Properties %s <==\n", ss.str().c_str());
  oe->compute();

  char*line = (char*)malloc(100*sizeof(char));

  std::stringstream ss2;
  if (cepa->cepa_level == -1)      sprintf(line,"CISD DIPOLE X");
  else if (cepa->cepa_level == -2) sprintf(line,"ACPF DIPOLE X");
  else if (cepa->cepa_level == -3) sprintf(line,"AQCC DIPOLE X");
  ss2 << oeprop_label.str() << " DIPOLE X";
  Process::environment.globals[line] =
    Process::environment.globals[ss2.str()];

  if (cepa->cepa_level == -1)      sprintf(line,"CISD DIPOLE Y");
  else if (cepa->cepa_level == -2) sprintf(line,"ACPF DIPOLE Y");
  else if (cepa->cepa_level == -3) sprintf(line,"AQCC DIPOLE Y");
  ss2.str(std::string());
  ss2 << oeprop_label.str() << " DIPOLE Y";
  Process::environment.globals[line] =
    Process::environment.globals[ss2.str()];

  if (cepa->cepa_level == -1)      sprintf(line,"CISD DIPOLE Z");
  else if (cepa->cepa_level == -2) sprintf(line,"ACPF DIPOLE Z");
  else if (cepa->cepa_level == -3) sprintf(line,"AQCC DIPOLE Z");
  ss2.str(std::string());
  ss2 << oeprop_label.str() << " DIPOLE Z";
  Process::environment.globals[line] =
    Process::environment.globals[ss2.str()];

  if (options.get_int("PRINT")>1){

     if (cepa->cepa_level == -1)      sprintf(line,"CISD QUADRUPOLE XX");
     else if (cepa->cepa_level == -2) sprintf(line,"ACPF QUADRUPOLE XX");
     else if (cepa->cepa_level == -3) sprintf(line,"AQCC QUADRUPOLE XX");
     ss2.str(std::string());
     ss2 << oeprop_label.str() << " QUADRUPOLE XX";
     Process::environment.globals["line"] =
       Process::environment.globals[ss2.str()];

     if (cepa->cepa_level == -1)      sprintf(line,"CISD QUADRUPOLE YY");
     else if (cepa->cepa_level == -2) sprintf(line,"ACPF QUADRUPOLE YY");
     else if (cepa->cepa_level == -3) sprintf(line,"AQCC QUADRUPOLE YY");
     ss2.str(std::string());
     ss2 << oeprop_label.str() << " QUADRUPOLE YY";
     Process::environment.globals[line] =
       Process::environment.globals[ss2.str()];

     if (cepa->cepa_level == -1)      sprintf(line,"CISD QUADRUPOLE ZZ");
     else if (cepa->cepa_level == -2) sprintf(line,"ACPF QUADRUPOLE ZZ");
     else if (cepa->cepa_level == -3) sprintf(line,"AQCC QUADRUPOLE ZZ");
     ss2.str(std::string());
     ss2 << oeprop_label.str() << " QUADRUPOLE ZZ";
     Process::environment.globals[line] =
       Process::environment.globals[ss2.str()];

     if (cepa->cepa_level == -1)      sprintf(line,"CISD QUADRUPOLE XY");
     else if (cepa->cepa_level == -2) sprintf(line,"ACPF QUADRUPOLE XY");
     else if (cepa->cepa_level == -3) sprintf(line,"AQCC QUADRUPOLE XY");
     ss2.str(std::string());
     ss2 << oeprop_label.str() << " QUADRUPOLE XY";
     Process::environment.globals[line] =
       Process::environment.globals[ss2.str()];

     if (cepa->cepa_level == -1)      sprintf(line,"CISD QUADRUPOLE XZ");
     else if (cepa->cepa_level == -2) sprintf(line,"ACPF QUADRUPOLE XZ");
     else if (cepa->cepa_level == -3) sprintf(line,"AQCC QUADRUPOLE XZ");
     ss2.str(std::string());
     ss2 << oeprop_label.str() << " QUADRUPOLE XZ";
     Process::environment.globals[line] =
       Process::environment.globals[ss2.str()];

     if (cepa->cepa_level == -1)      sprintf(line,"CISD QUADRUPOLE YZ");
     else if (cepa->cepa_level == -2) sprintf(line,"ACPF QUADRUPOLE YZ");
     else if (cepa->cepa_level == -3) sprintf(line,"AQCC QUADRUPOLE YZ");
     ss2.str(std::string());
     ss2 << oeprop_label.str() << " QUADRUPOLE YZ";
     Process::environment.globals[line] =
       Process::environment.globals[ss2.str()];
  }

  free(D1);
}

// build the 1-electron density
void BuildD1(long int o,long int v,double*t1,double*ta,double*tb,double c0,double*D1){
  long int id,i,j,k,a,b,c,sg,p,count,sg2,nmo=o+v;
  double sum;
  memset((void*)D1,'\0',(o+v)*(o+v)*sizeof(double));
  double*tempd = (double*)malloc(v*v*sizeof(double));

  F_DCOPY(o*o*v*v,tb,1,ta,1);
  for (a=0,id=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              for (j=0; j<o; j++){
                  ta[id++] -= tb[b*o*o*v+a*o*o+i*o+j];
              }
          }
      }
  }


  // D(a,b)
  F_DGEMM('t','n',v,v,o*o*v,1.0,tb,o*o*v,tb,o*o*v,0.0,tempd,v);
  F_DGEMM('t','n',v,v,o*o*v,0.5,ta,o*o*v,ta,o*o*v,1.0,tempd,v);
  F_DGEMM('t','n',v,v,o,1,t1,o,t1,o,1.0,tempd,v);
  for (a=0; a<v; a++){
      for (b=0; b<v; b++){
          D1[(a+o)*nmo+(b+o)] = tempd[a*v+b];
      }
  }
 
  // D(i,j)
  F_DGEMM('n','t',o,o,o*v*v,-1.0,tb,o,tb,o,0.0,tempd,o);
  F_DGEMM('n','t',o,o,o*v*v,-0.5,ta,o,ta,o,1.0,tempd,o);
  F_DGEMM('n','t',o,o,v,-1.0,t1,o,t1,o,1.0,tempd,o);
  for (i=0; i<o; i++){
      for (j=0; j<o; j++){
          D1[i*nmo + j] = tempd[i*o+j];
      }
      D1[i*nmo+i] += 1.0;
  }

  // hmmm ... i could do this as a dgemv, but i'd have sort ta and tb
  for (i=0; i<o; i++){
      for (a=0; a<v; a++){
          sum = t1[a*o+i]*c0;
          for (j=0; j<o; j++){
              for (b=0; b<v; b++){
                  sum += t1[b*o+j]*tb[a*o*o*v+b*o*o+i*o+j];
                  sum += t1[b*o+j]*ta[a*o*o*v+b*o*o+i*o+j];
              }
          }
          D1[i*nmo+a+o] = D1[(a+o)*nmo+i] = sum;
      }
  }

  free(tempd);

}

// normalize the wave function and return the leading coefficient
double Normalize(long int o,long int v,double*t1,double*t2,int cepa_level){

  if (cepa_level == 0) return 1.0;

  double nrm = 1.0;
  double fac = 1.0;
  if (cepa_level == -2)      fac = 1.0/o;
  else if (cepa_level == -3) fac = 1.0-(2.0*o-2.0)*(2.0*o-3.0) / (2.0*o*(2.0*o-1.0));

  long int i,j,a,b,id;
  double dum,sum=0.;

  for (a=0,id=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              for (j=0; j<o; j++){
                  dum  = t2[id++];
                  sum -= dum*dum;
                  dum -= t2[b*o*o*v+a*o*o+i*o+j];
                  sum -= .5*dum*dum;
              }
          }
      }
  }
  for (a=0,id=0; a<v; a++){
     for (i=0; i<o; i++){
          dum  = t1[id++];
          sum -= 2.*dum*dum;
      }
  }
  nrm = sqrt(1.0-fac*sum);

  for (i=0; i<o*o*v*v; i++) t2[i] /= nrm;
  for (i=0; i<o*v; i++)     t1[i] /= nrm;

  return 1.0/nrm;
}


}// end of namespace psi

