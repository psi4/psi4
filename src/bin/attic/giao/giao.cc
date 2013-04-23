/*! \defgroup GIAO giao: Gauge-Including Atomic Orbitals */

/*! 
** \file
** \ingroup GIAO
** \brief Gauge-Including Atomic Orbitals
*/


#include <stdexcept>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libiwl/iwl.h>
#include <psifiles.h>

extern "C" {
  FILE* infile;
  FILE* outfile;
  char* psi_file_prefix;
}

namespace psi { namespace giao {

static const int IOFFMAX = 1024;
static int ioff[IOFFMAX];
static const double CUTOFF = 1.0E-10;
#define INDEX42(ab,cd,X) (ab)*(X)*(X)+(cd)
#define INDEX44(a,b,c,d,X) INDEX42((a)*(X)+(b),(c)*(X)+(d),(X))

void init_misc();
int iwl_buf_rd_all_noperm(struct iwlbuf *Buf, double *ints, int nbf, int unsymm, 
                          int printflg);
void transform_1e_2(double** SS, int ns, double** U, int nu, double** UU);
double* transform_2e_4(double* SSSS, int ns, double** U, int nu, double* scr);
void tildeize_1e(double** H0, double** H1, double** S1, int nbf, double** Ht1);
void tildeize_2e(double* g0, double* g1, double** S1, int nbf, double* gt1);
void form_hamiltonian_deriv(double** Ht1, double* dgt, int nmo, int ndocc, double** dFt);
void write_2e_iwl(double* ints, int nbf, int filenum);
void print_1e(const char* label, double** ints, int nbf);
void print_2e(const char* label, double* ints, int nbf, double cutoff);
void print_2e_canon(const char* label, double* ints, int nbf, double cutoff);

}} // namespace psi::giao

int main(int argc, char* argv[])
{
  using namespace psi::giao;
  
  int errcod = psi_start(&infile,&outfile,&psi_file_prefix,argc-1,argv+1,0);
  if (errcod != PSI_RETURN_SUCCESS)
    throw std::runtime_error("main -- psi_start failed");
  tstart(outfile);
  psio_init(); psio_ipv1_config();
  chkpt_init(PSIO_OPEN_OLD);
  init_misc();
  
  int nmo = chkpt_rd_nmo();
  int nso = chkpt_rd_nso();
  int nbf = nso;
  int nao = chkpt_rd_nao();
  bool puream = (nbf != nao);
  int* clsdpi = chkpt_rd_clsdpi();
  int ndocc = clsdpi[0];
  double** Cso = chkpt_rd_scf();
  double** Uso2ao = chkpt_rd_usotao();
  double** Cao = block_matrix(nao,nmo);
  mmult(Uso2ao,1,Cso,0,Cao,0,nao,nso,nmo,0);
  free_block(Uso2ao);

  struct iwlbuf g_Buf, dgdBx_Buf, dgdBy_Buf, dgdBz_Buf;
//  iwl_buf_init(&g_Buf, PSIF_SO_TEI, 0.0, 1, 1);
  iwl_buf_init(&g_Buf, 47, 0.0, 1, 1);
  iwl_buf_init(&dgdBx_Buf, PSIF_AO_DGDBX, 0.0, 1, 1);
  iwl_buf_init(&dgdBy_Buf, PSIF_AO_DGDBY, 0.0, 1, 1);
  iwl_buf_init(&dgdBz_Buf, PSIF_AO_DGDBZ, 0.0, 1, 1);

  int nijkl = nao*nao*nao*nao;
  double* g_ao = new double[nijkl];
  double* dgdBx_ao = new double[nijkl];
  double* dgdBy_ao = new double[nijkl];
  double* dgdBz_ao = new double[nijkl];
//  iwl_buf_rd_all_noperm(&g_Buf, g_ao, nso, 1, 0);
  iwl_buf_rd_all_noperm(&g_Buf, g_ao, nao, 0, 0);
  iwl_buf_rd_all_noperm(&dgdBx_Buf, dgdBx_ao, nao, 0, 0);
  iwl_buf_rd_all_noperm(&dgdBy_Buf, dgdBy_ao, nao, 0, 0);
  iwl_buf_rd_all_noperm(&dgdBz_Buf, dgdBz_ao, nao, 0, 0);
  
  iwl_buf_close(&g_Buf, 1);
  iwl_buf_close(&dgdBx_Buf, 1);
  iwl_buf_close(&dgdBy_Buf, 1);
  iwl_buf_close(&dgdBz_Buf, 1);
  
  // Make MO basis 2-e integrals
  double* scratch_buf = new double[nijkl];
  //print_2e("Canonical list of ERIs", g_ao, nao, CUTOFF);
//  double* g = transform_2e_4(g_ao, nao, Cao, nmo, scratch_buf);
  double* g = transform_2e_4(g_ao, nso, Cso, nmo, scratch_buf);
  //print_2e_canon("Canonical list of ERIs", g, nmo, CUTOFF);
  double* dgdBx = transform_2e_4(dgdBx_ao, nao, Cao, nmo, scratch_buf);
  //print_2e("dg/dBx (GIAO MO basis)", dgdBx, nmo, CUTOFF);
  double* dgdBy = transform_2e_4(dgdBy_ao, nao, Cao, nmo, scratch_buf);
  //print_2e("dg/dBy (GIAO MO basis)", dgdBy, nmo, CUTOFF);
  double* dgdBz = transform_2e_4(dgdBz_ao, nao, Cao, nmo, scratch_buf);
  //print_2e("dg/dBz (GIAO MO basis)", dgdBz, nmo, CUTOFF);
  
  // Make MO basis dS/dBi
  double** dSdBx_ao = block_matrix(nao,nao);
  double** dSdBy_ao = block_matrix(nao,nao);
  double** dSdBz_ao = block_matrix(nao,nao);
  int nbra = nao * nao;
  iwl_rdone(PSIF_OEI, PSIF_AO_DSDB_X, dSdBx_ao[0], nbra, 0, 0, outfile);
  iwl_rdone(PSIF_OEI, PSIF_AO_DSDB_Y, dSdBy_ao[0], nbra, 0, 0, outfile);
  iwl_rdone(PSIF_OEI, PSIF_AO_DSDB_Z, dSdBz_ao[0], nbra, 0, 0, outfile);
  double** dSdBx = block_matrix(nmo,nmo);
  double** dSdBy = block_matrix(nmo,nmo);
  double** dSdBz = block_matrix(nmo,nmo);
  transform_1e_2(dSdBx_ao, nao, Cao, nmo, dSdBx);
  transform_1e_2(dSdBy_ao, nao, Cao, nmo, dSdBy);
  transform_1e_2(dSdBz_ao, nao, Cao, nmo, dSdBz);
  free_block(dSdBx_ao);
  free_block(dSdBy_ao);
  free_block(dSdBz_ao);
  print_1e("dS/dBx (GIAO MO basis)", dSdBx, nmo);
  print_1e("dS/dBy (GIAO MO basis)", dSdBy, nmo);
  print_1e("dS/dBz (GIAO MO basis)", dSdBz, nmo);

  // Make MO basis dh/dBi
  double** dHdBx_ao = block_matrix(nao,nao);
  double** dHdBy_ao = block_matrix(nao,nao);
  double** dHdBz_ao = block_matrix(nao,nao);
  iwl_rdone(PSIF_OEI, PSIF_AO_DHDB_X, dHdBx_ao[0], nbra, 0, 0, outfile);
  iwl_rdone(PSIF_OEI, PSIF_AO_DHDB_Y, dHdBy_ao[0], nbra, 0, 0, outfile);
  iwl_rdone(PSIF_OEI, PSIF_AO_DHDB_Z, dHdBz_ao[0], nbra, 0, 0, outfile);
  double** dHdBx = block_matrix(nmo,nmo);
  double** dHdBy = block_matrix(nmo,nmo);
  double** dHdBz = block_matrix(nmo,nmo);
  transform_1e_2(dHdBx_ao, nao, Cao, nmo, dHdBx);
  transform_1e_2(dHdBy_ao, nao, Cao, nmo, dHdBy);
  transform_1e_2(dHdBz_ao, nao, Cao, nmo, dHdBz);
  free_block(dHdBx_ao);
  free_block(dHdBy_ao);
  free_block(dHdBz_ao);
  print_1e("dH/dBx (GIAO MO basis)", dHdBx, nmo);
  print_1e("dH/dBy (GIAO MO basis)", dHdBy, nmo);
  print_1e("dH/dBz (GIAO MO basis)", dHdBz, nmo);

  // Make MO basis h
  int ntri = ioff[nso];
  double* V_so = new double[ntri];
  double* T_so = new double[ntri];
  iwl_rdone(PSIF_OEI, PSIF_SO_T, T_so, ntri, 0, 0, outfile);
  iwl_rdone(PSIF_OEI, PSIF_SO_V, V_so, ntri, 0, 0, outfile);
  double** H_so = block_matrix(nso,nso);
  int ij=0;
  for(int i=0; i<nso; i++)
    for(int j=0; j<=i; j++, ij++) {
      double value = T_so[ij] + V_so[ij];
      H_so[i][j] = H_so[j][i] = value;
    }
  delete[] T_so; delete[] V_so;
  double** H = block_matrix(nmo,nmo);
  transform_1e_2(H_so, nso, Cso, nmo, H);
  free_block(H_so);
  print_1e("H (MO basis)", H, nmo);

  // For comparison make L in MO basis
  double** Lx_ao = block_matrix(nao,nao);
  double** Ly_ao = block_matrix(nao,nao);
  double** Lz_ao = block_matrix(nao,nao);
  nbra = nao * nao;
  iwl_rdone(PSIF_OEI, PSIF_AO_LX, Lx_ao[0], nbra, 0, 0, outfile);
  iwl_rdone(PSIF_OEI, PSIF_AO_LY, Ly_ao[0], nbra, 0, 0, outfile);
  iwl_rdone(PSIF_OEI, PSIF_AO_LZ, Lz_ao[0], nbra, 0, 0, outfile);
  double** Lx = block_matrix(nmo,nmo);
  double** Ly = block_matrix(nmo,nmo);
  double** Lz = block_matrix(nmo,nmo);
  transform_1e_2(Lx_ao, nao, Cao, nmo, Lx);
  transform_1e_2(Ly_ao, nao, Cao, nmo, Ly);
  transform_1e_2(Lz_ao, nao, Cao, nmo, Lz);
  free_block(Lx_ao);
  free_block(Ly_ao);
  free_block(Lz_ao);
  print_1e("Lx (MO basis)", Lx, nmo);
  print_1e("Ly (MO basis)", Ly, nmo);
  print_1e("Lz (MO basis)", Lz, nmo);
  delete[] Lx;
  delete[] Ly;
  delete[] Lz;

  // Form F and dF/dBi (just to check against ACES2 and other programs) in MO basis from dH/dBi and dg/dBi
  double** F = block_matrix(nmo,nmo);
  double** dFdBx = block_matrix(nmo,nmo);
  double** dFdBy = block_matrix(nmo,nmo);
  double** dFdBz = block_matrix(nmo,nmo);
  form_hamiltonian_deriv(H,g,nmo,ndocc,F);
  form_hamiltonian_deriv(dHdBx,dgdBx,nmo,ndocc,dFdBx);
  form_hamiltonian_deriv(dHdBy,dgdBy,nmo,ndocc,dFdBy);
  form_hamiltonian_deriv(dHdBz,dgdBz,nmo,ndocc,dFdBz);
  print_1e("F (MO basis)", F, nmo);
  print_1e("dF/dBx (GIAO MO basis)", dFdBx, nmo);
  print_1e("dF/dBy (GIAO MO basis)", dFdBy, nmo);
  print_1e("dF/dBz (GIAO MO basis)", dFdBz, nmo);
  free_block(dFdBx);
  free_block(dFdBy);
  free_block(dFdBz);

  // Now form dH~/dBi in MO basis from H, dS/dBi and dH/dBi
  double** dHtdBx = block_matrix(nmo,nmo);
  double** dHtdBy = block_matrix(nmo,nmo);
  double** dHtdBz = block_matrix(nmo,nmo);
  tildeize_1e(H,dHdBx,dSdBx,nmo,dHtdBx);
  tildeize_1e(H,dHdBy,dSdBy,nmo,dHtdBy);
  tildeize_1e(H,dHdBz,dSdBz,nmo,dHtdBz);
  print_1e("dHt/dBx (GIAO MO basis)", dHtdBx, nmo);
  print_1e("dHt/dBy (GIAO MO basis)", dHtdBy, nmo);
  print_1e("dHt/dBz (GIAO MO basis)", dHtdBz, nmo);
  
  // Write out dH~/dBi in MO basis to PSIF_OEI
  iwl_wrtone(PSIF_OEI, "dh~/dBx (MO Basis)", nmo*nmo, dHtdBx[0]);
  iwl_wrtone(PSIF_OEI, "dh~/dBy (MO Basis)", nmo*nmo, dHtdBy[0]);
  iwl_wrtone(PSIF_OEI, "dh~/dBz (MO Basis)", nmo*nmo, dHtdBz[0]);
  
  // intermediate cleanup
  free_block(dHdBx);
  free_block(dHdBy);
  free_block(dHdBz);

  // Also form dg~/dBi in MO basis from g, dS/dBi and dg/dBi
  nijkl = nmo*nmo*nmo*nmo;
  double* dgtdBx = new double[nijkl];
  double* dgtdBy = new double[nijkl];
  double* dgtdBz = new double[nijkl];
  tildeize_2e(g,dgdBx,dSdBx,nmo,dgtdBx);
  tildeize_2e(g,dgdBy,dSdBy,nmo,dgtdBy);
  tildeize_2e(g,dgdBz,dSdBz,nmo,dgtdBz);
  //print_2e("dgt/dBx (GIAO MO basis)", dgtdBx, nmo, CUTOFF);
  //print_2e("dgt/dBy (GIAO MO basis)", dgtdBy, nmo, CUTOFF);
  //print_2e("dgt/dBz (GIAO MO basis)", dgtdBz, nmo, CUTOFF);  
  delete[] dgdBx;
  delete[] dgdBy;
  delete[] dgdBz;
  
  // Write out dg~/dBi in MO basis to IWL files
  write_2e_iwl(dgtdBx, nmo, 48);
  write_2e_iwl(dgtdBy, nmo, 49);
  write_2e_iwl(dgtdBz, nmo, 50);

  // Finally form dF~/dBi in MO basis from dH~/dBi and dg~/dBi
  double** dFtdBx = block_matrix(nmo,nmo);
  double** dFtdBy = block_matrix(nmo,nmo);
  double** dFtdBz = block_matrix(nmo,nmo);
  form_hamiltonian_deriv(dHtdBx,dgtdBx,nmo,ndocc,dFtdBx);
  form_hamiltonian_deriv(dHtdBy,dgtdBx,nmo,ndocc,dFtdBy);
  form_hamiltonian_deriv(dHtdBz,dgtdBx,nmo,ndocc,dFtdBz);
  print_1e("dFt/dBx (GIAO MO basis)", dFtdBx, nmo);
  print_1e("dFt/dBy (GIAO MO basis)", dFtdBy, nmo);
  print_1e("dFt/dBz (GIAO MO basis)", dFtdBz, nmo);
  iwl_wrtone(PSIF_OEI, PSIF_MO_DFDB_X, nmo*nmo, dFtdBx[0]);
  iwl_wrtone(PSIF_OEI, PSIF_MO_DFDB_Y, nmo*nmo, dFtdBy[0]);
  iwl_wrtone(PSIF_OEI, PSIF_MO_DFDB_Z, nmo*nmo, dFtdBz[0]);

  free_block(F);
  free_block(H);
  free_block(dFtdBx);
  free_block(dFtdBy);
  free_block(dFtdBz);
  free_block(dHtdBx);
  free_block(dHtdBy);
  free_block(dHtdBz);
  free_block(dSdBx);
  free_block(dSdBy);
  free_block(dSdBz);
  delete[] dgtdBx;
  delete[] dgtdBy;
  delete[] dgtdBz;
  delete[] g;

  free_block(Cso);
  free_block(Cao);
  
  chkpt_close();
  psio_done();
  tstop(outfile);
  psi_stop(infile,outfile,psi_file_prefix);
  exit(0);
}

extern "C" const char* gprgid()
{
  const char* prgid = "GIAO";
  return prgid;
}

namespace psi { namespace giao {

void init_misc()
{
  ioff[0] = 0;
  for (int i=1; i<IOFFMAX; i++){
    ioff[i] = ioff[i-1] + i;
  }
}


int iwl_buf_rd_all_noperm(struct iwlbuf *Buf, double *ints, int nbf, int unsymm, 
                          int printflg)
{  
  Label* lblptr = Buf->labels;
  Value* valptr = Buf->values;
  
  int lastbuf = Buf->lastbuf;
  
  for (int idx=4*Buf->idx; Buf->idx<Buf->inbuf; Buf->idx++) {
    long int p = (long int) fabs((long int) lblptr[idx++]);
    long int q = (long int) lblptr[idx++];
    long int r = (long int) lblptr[idx++];
    long int s = (long int) lblptr[idx++];

    long int pqrs = INDEX44(p,q,r,s,nbf);
    
    double value = (double) valptr[Buf->idx];
    ints[pqrs] = value;
    if (unsymm) {
      ints[INDEX44(q,p,r,s,nbf)] = value;
      ints[INDEX44(p,q,s,r,nbf)] = value;
      ints[INDEX44(q,p,s,r,nbf)] = value;
      ints[INDEX44(r,s,p,q,nbf)] = value;
      ints[INDEX44(s,r,p,q,nbf)] = value;
      ints[INDEX44(r,s,q,p,nbf)] = value;
      ints[INDEX44(s,r,q,p,nbf)] = value;
    }
    
    if (printflg) 
      fprintf(outfile, "<%2d %2d %2d %2d [[%3ld]] = %20.10lf\n",
	      p, q, r, s, pqrs, ints[pqrs]) ;
    
  } /*! end loop through current buffer */
  
   /*! read new PSI buffers */
  while (!lastbuf) {
    iwl_buf_fetch(Buf);
    lastbuf = Buf->lastbuf;
    
    for (int idx=4*Buf->idx; Buf->idx<Buf->inbuf; Buf->idx++) {
      long int p = (long int) fabs((long int) lblptr[idx++]);
      long int q = (long int) lblptr[idx++];
      long int r = (long int) lblptr[idx++];
      long int s = (long int) lblptr[idx++];

      long int pqrs = INDEX44(p,q,r,s,nbf);

      double value = (double) valptr[Buf->idx];
      ints[pqrs] = value;
      if (unsymm) {
        ints[INDEX44(q,p,r,s,nbf)] = value;
        ints[INDEX44(p,q,s,r,nbf)] = value;
        ints[INDEX44(q,p,s,r,nbf)] = value;
        ints[INDEX44(r,s,p,q,nbf)] = value;
        ints[INDEX44(s,r,p,q,nbf)] = value;
        ints[INDEX44(r,s,q,p,nbf)] = value;
        ints[INDEX44(s,r,q,p,nbf)] = value;
      }

      if (printflg) 
	fprintf(outfile, "<%d %d %d %d [[%ld]] = %20.10lf\n",
		p, q, r, s, pqrs, ints[pqrs]) ;
      
    } /*! end loop through current buffer */
    
  } /*! end loop over reading buffers */
  
  return(0); /*! we must have reached the last buffer at this point */
}


void transform_1e_2(double** SS, int ns, double** U, int nu, double** UU)
{
  double** SU = block_matrix(ns, nu);
  
  mmult(SS,0,U,0,SU,0,ns,ns,nu,0);
  mmult(U,1,SU,0,UU,0,nu,ns,nu,0);
  
  free_block(SU);
}

double* transform_2e_4(double* SSSS, int ns, double** U, int nu, double* scr)
{
  // First step
  double* USSS = scr;
  int nsss = ns*ns*ns;
  memset(USSS,0,sizeof(double)*nsss*nu);
  int ssss=0;
  for(int i=0; i<ns; i++) {
    for(int jkl=0; jkl<nsss; jkl++, ssss++) {
      double value = SSSS[ssss];
      int usss = jkl;
      for(int p=0; p<nu; p++, usss+=nsss) {
        USSS[usss] += value * U[i][p];
      }
    }
  }

  // Second step
  double* UUSS = SSSS;
  int nss = ns*ns;
  memset(UUSS,0,sizeof(double)*nss*nu*nu);
  int usss=0;
  for(int p=0; p<nu; p++) {
    for(int j=0; j<ns; j++) {
      for(int kl=0; kl<nss; kl++, usss++) {
        double value = USSS[usss];
        int uuss = p*nu*nss + kl;
        for(int q=0; q<nu; q++, uuss+=nss) {
          UUSS[uuss] += value * U[j][q];
        }
      }
    }
  }

  // Third step
  double* UUUS = scr;
  int nuu = nu*nu;
  memset(UUUS,0,sizeof(double)*ns*nuu*nu);
  int uuss=0;
  for(int pq=0; pq<nuu; pq++) {
    for(int k=0; k<ns; k++) {
      for(int l=0; l<ns; l++, uuss++) {
        double value = UUSS[uuss];
        int uuus = pq*nu*ns + l;
        for(int r=0; r<nu; r++, uuus+=ns) {
          UUUS[uuus] += value * U[k][r];
        }
      }
    }
  }

  // Fourth step
  double* UUUU = SSSS;
  int nuuu = nuu*nu;
  memset(UUUU,0,sizeof(double)*nuu*nuu);
  int uuus=0;
  for(int pqr=0; pqr<nuuu; pqr++) {
    for(int l=0; l<ns; l++, uuus++) {
      double value = UUUS[uuus];
      int uuuu = pqr*nu;
      for(int q=0; q<nu; q++, uuuu++) {
        UUUU[uuuu] += value * U[l][q];
      }
    }
  }

  return UUUU;
}

void tildeize_1e(double** H0, double** H1, double** S1, int nbf, double** Ht1)
{
/*  double** S1xH0 = block_matrix(nbf,nbf);
  mmult(S1,0,H0,0,S1xH0,0,nbf,nbf,nbf,0);
  double** H0xS1t = block_matrix(nbf,nbf);
  mmult(H0,0,S1,1,H0xS1t,0,nbf,nbf,nbf,0);
  for(int i=0; i<nbf; i++)
    for(int j=0; j<nbf; j++) {
      Ht1[i][j] = H1[i][j] - 0.5 * (S1xH0[i][j] - H0xS1t[i][j]);
    }
  free_block(S1xH0);
  free_block(H0xS1t);*/
  
  for(int i=0; i<nbf; i++) {
    for(int j=0; j<nbf; j++) {
      
      double value = 0.0;
      for(int p=0; p<nbf; p++)
        value += S1[i][p] * H0[p][j] - S1[j][p] * H0[i][p];
      
      Ht1[i][j] = H1[i][j] - 0.5 * value;
    }
  }
  
  return;
}

void tildeize_2e(double* g0, double* g1, double** S1, int nbf, double* gt1)
{
  long int ijkl = 0;
  for(long int i=0; i<nbf; i++) {
    for(long int j=0; j<nbf; j++) {
      for(long int k=0; k<nbf; k++) {
         for(long int l=0; l<nbf; l++, ijkl++) {
         
           double value = 0.0;
        
           for(long int p=0; p<nbf; p++) {
             long int pjkl = INDEX44(p,j,k,l,nbf);
             long int ipkl = INDEX44(i,p,k,l,nbf);
             long int ijpl = INDEX44(i,j,p,l,nbf);
             long int ijkp = INDEX44(i,j,k,p,nbf);
             value += S1[i][p] * g0[pjkl] - S1[j][p] * g0[ipkl]
                    + S1[k][p] * g0[ijpl] - S1[l][p] * g0[ijkp];
           }
           
           gt1[ijkl] = g1[ijkl] - 0.5 * value;
        }
      }
    }
  }

  return;
}

void form_hamiltonian_deriv(double** Ht1, double* gt1, int nmo, int ndocc, double** Ft1)
{
  for(long int p=0; p<nmo; p++) {
    for(long int q=0; q<nmo; q++) {
      
      double g = 0.0;
      for(long int o=0; o<ndocc; o++) {
        long int pqoo = INDEX44(p,q,o,o,nmo);
        long int pooq = INDEX44(p,o,o,q,nmo);
        g += 2.0*gt1[pqoo] - gt1[pooq];
      }
      Ft1[p][q] = Ht1[p][q] + g;

    }
  }

  return;
}

void write_2e_iwl(double* ints, int nbf, int filenum)
{
  struct iwlbuf IWLBuf;
  iwl_buf_init(&IWLBuf, filenum, 0.0, 0, 0);
  
  long int ijkl = 0;
  for(long int i=0; i<nbf; i++) {
    for(long int j=0; j<nbf; j++) {
      for(long int k=0; k<nbf; k++) {
         for(long int l=0; l<nbf; l++, ijkl++) {
           double value = ints[ijkl];
           iwl_buf_wrt_val(&IWLBuf, i, j, k, l, value, 0, NULL, 0);
        }
      }
    }
  }
  
  iwl_buf_flush(&IWLBuf, 1);
  iwl_buf_close(&IWLBuf, 1);

  return;
}

void print_1e(const char* label, double** ints, int nbf)
{
  fprintf(outfile, "  -%s:\n", label);
  print_mat(ints,nbf,nbf,outfile);
  return;
}

void print_2e(const char* label, double* ints, int nbf, double cutoff)
{
  fprintf(outfile, "  -%s:\n", label);

  for(int i=0; i<nbf; i++) {
    for(int j=0; j<nbf; j++) {
      for(int k=0; k<nbf; k++) {
        for(int l=0; l<nbf; l++) {
          int ijkl = INDEX44(i,j,k,l,nbf);
          double value = ints[ijkl];

          if (fabs(value) > cutoff)
            fprintf(stdout, ">%d %d %d %d [%d] = %20.10lf\n",
                    i, j, k, l, ijkl, value);          
        }
      }
    }
  }

  return;
}

void print_2e_canon(const char* label, double* ints, int nbf, double cutoff)
{
  fprintf(outfile, "  -%s:\n", label);
  
  for(int i=0; i<nbf; i++) {
    for(int j=0; j<=i; j++) {
      int ij = ioff[i]+j;
      for(int k=0; k<=i; k++) {
        int lmax = (i == k) ? j : k;
        for(int l=0; l<=lmax; l++) {
          int kl = ioff[k] + l;
          int ijkl = INDEX44(i,j,k,l,nbf);
          double value = ints[ijkl];

          if (fabs(value) > cutoff)
            fprintf(outfile, ">%d %d %d %d [%d] [%d] = %20.10lf\n",
                    i, j, k, l, ij, kl, value);          
        }
      }
    }
  }

  return;
}

}} // namespace psi:giao
