/*! \file
    \ingroup DETCAS
    \brief Enter brief description of file here 
*/
#include <cstdlib>
#include <cstdio>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "globaldefs.h"
#include "globals.h"

namespace psi { namespace detcas {

void calc_F_act(double *F_act, int nmo, int firstact, int lastact, 
                double **onepdm, double *tei);
void check_F_act(double *F_act, int nmo, int firstact, int lastact, 
                double **onepdm, double *tei);
void test_lag(int nbf, int ncore, int npop, double *onei,
           double *tei, double **opdm, double *tpdm);
void test_fzc(int nbf, int ncore, double *onei, double *tei);
void test_lag2(int nbf, int ncore, int npop, double *onei,
           double *tei, double **opdm, double *tpdm);

void form_F_act(void)
{
  int ncore;

  /* Form the intermediates we need */
  CalcInfo.F_act = init_array(CalcInfo.nbstri);
  ncore = CalcInfo.num_fzc_orbs + CalcInfo.num_cor_orbs;
  calc_F_act(CalcInfo.F_act, CalcInfo.nmo, ncore,  CalcInfo.npop, 
             CalcInfo.opdm, CalcInfo.twoel_ints);
  /*
  check_F_act(CalcInfo.F_act, CalcInfo.nmo, ncore,  CalcInfo.npop, 
             CalcInfo.opdm, CalcInfo.twoel_ints);
  test_lag(CalcInfo.nmo, ncore, CalcInfo.npop, CalcInfo.onel_ints,
           CalcInfo.twoel_ints, CalcInfo.opdm, CalcInfo.tpdm);
  test_fzc(CalcInfo.nmo, ncore, CalcInfo.onel_ints, CalcInfo.twoel_ints);

  test_lag2(CalcInfo.nmo, ncore, CalcInfo.npop, CalcInfo.onel_ints,
           CalcInfo.twoel_ints, CalcInfo.opdm, CalcInfo.tpdm);
  */

  if (Params.print_lvl > 3) {
    fprintf(outfile, "\nActive Fock matrix:\n");
    print_array(CalcInfo.F_act, CalcInfo.nmo, outfile);
  }

}


/*
** calc_F_act
**
** This forms the "active" Fock matrix as defined by eq. 13b of 
** Siegbahn, Heiberg, Roos, and Levy, Physica Scripta 21, 323 (1980).
**
** F_act_{pq} = \sum_{uv}^{active} \gamma_{uv} [ (pq|uv) - 1/2 (pu|qv) ] 
**
** This thing looks symmetric so I will use that in the calculation
** I am assuming that the pairs (p,q) are always given such that p>=q
**
*/
void calc_F_act(double *F_act, int nmo, int firstact, int lastact, 
                double **onepdm, double *tei)
{

  int p,q,pq,u,v,pu,uv,qv,pquv,puqv;
  double val, I1, I2;

  for (p=0; p<nmo; p++) {
    for (q=0; q<=p; q++) {
      pq = ioff[p] + q;

      val = 0.0;

      for (u=firstact; u<lastact; u++) {
        pu = INDEX(p,u);
        for (v=firstact; v<lastact; v++) {
          uv = INDEX(u,v);
          qv = INDEX(q,v);
          pquv = INDEX(pq,uv);
          puqv = INDEX(pu,qv);  
          I1 = tei[pquv];
          I2 = tei[puqv];
          val += onepdm[u][v] * (I1 - 0.5 * I2);
         }
      }

      F_act[pq] = val;
    }
  }
}

void check_F_act(double *F_act, int nmo, int firstact, int lastact, 
                double **onepdm, double *tei)
{
  int p, q, pq, u, v, uv, pquv;
  int pu, qv, puqv;

  double gamma, I1, I2, val, sum = 0.0;

  /* this is a stupid test to check my F_act routine */

  /* hardwire the test ... check F^{act}_{6,0} */
  p = 6; q = 0;
  pq = ioff[p] + q;

  for (u=1; u<6; u++) { 
    for (v=1; v<=u; v++) {
      uv = ioff[u] + v;
      pquv = ioff[pq] + uv;
      gamma = onepdm[u][v];
      I1 = tei[pquv];
      pu = ioff[p] + u;
      qv = ioff[v] + q;
      puqv = ioff[pu] + qv;
      I2 = tei[puqv];

      val = 4.0 * gamma * I1;
      sum += val;

      if (val != 0.0) {
        fprintf(outfile, "gamma[%d][%d] = %12.6lf\n", u, v, gamma);
        fprintf(outfile, "TEI[%2d %2d %2d %2d] = %12.6lf\n", p,q,u,v,I1);
        fprintf(outfile, "contrib = %12.6lf, sum = %12.6lf\n\n", val, sum);
      }

      val = -4.0 * gamma * 0.5 * I2;
      sum += val;

      if (val != 0.0) {
        fprintf(outfile, "TEI[%2d %2d %2d %2d] = %12.6lf\n", p,u,q,v,I2);
        fprintf(outfile, "contrib = %12.6lf, sum = %12.6lf\n\n", val, sum);
      }
    }
  }

  fprintf(outfile, "Final F_act{0,6} = %12.6lf\n", sum / 4.0);
 
}

void test_lag(int nbf, int ncore, int npop, double *oei,
           double *tei, double **opdm, double *tpdm)
{
  int p, q, r, s, t, pr, pq, qr, st, prst, qrst;
  double *lag;

  lag = init_array(nbf * (nbf + 1) / 2);
 
  for (p=0; p<nbf; p++) {
    for (q=0; (q<=p && q<npop); q++) {
      pq = ioff[p] + q;

      for (r=0; r<npop; r++) {
        pr = INDEX(p,r);
        lag[pq] += opdm[q][r] * oei[pr];

        for (s=0; s<npop; s++) {
          for (t=0; t<npop; t++) {
            qr = INDEX(q,r);
            st = INDEX(s,t);
            prst = INDEX(pr,st);
            qrst = INDEX(qr,st); 
            lag[pq] += tei[prst] * tpdm[qrst]; 
          }
        }
      }
    }
  }

  fprintf(outfile, "\nTest lag:\n");
  print_array(lag, nbf, outfile);
  free(lag);
}


void test_fzc(int nbf, int ncore, double *oei, double *tei)
{
  int p, q, pq, k, kk, pqkk, pk, qk, pkqk;
  double *fzc_op;
 
  fzc_op = init_array(nbf * (nbf+1) / 2);
  for (p=0,pq=0; p<nbf; p++) {
    for(q=0; q<=p; q++,pq++) {
      fzc_op[pq] = oei[pq];
      for (k=0; k<ncore; k++) {
        kk = ioff[k] + k;
        pqkk = INDEX(pq,kk);
        pk = INDEX(p,k);
        qk = INDEX(q,k);
        pkqk = INDEX(pk,qk);
        fzc_op[pq] += (2.0 * tei[pqkk] - tei[pkqk]);
      } 
    }
  }

  print_array(fzc_op, nbf, outfile);
  free(fzc_op);

}


void test_lag2(int nbf, int ncore, int npop, double *oei,
           double *tei, double **opdm, double *tpdm)
{
  int p, q, r, s, t, pr, pq, qr, st, prst, qrst;
  double val, sum=0.0;

  /* hardwire the test ... check L_{6,0} */
  p = 6; q = 0;

  pq = ioff[p] + q;

  for (r=0; r<npop; r++) {
    pr = INDEX(p,r);
    val = 2.0 * opdm[q][r] * oei[pr];
    fprintf(outfile,
       "2.0 * opdm[%d][%d] (%12.6lf) x oei[%d][%d] (%12.6lf) = %12.6lf\n",
       q, r, opdm[q][r], p, r, oei[pr], val);
    sum += val;

    for (s=0; s<npop; s++) {
      for (t=0; t<npop; t++) {
        qr = INDEX(q,r);
        st = INDEX(s,t);
        prst = INDEX(pr,st);
        qrst = INDEX(qr,st); 
        val = 2.0 * tei[prst] * tpdm[qrst]; 
        if (val != 0.0) {
          fprintf(outfile, "tei[%d %d %d %d]  = %12.6lf\n", p, r, s, t,
                  tei[prst]);
          fprintf(outfile, "tpdm[%d %d %d %d] = %12.6lf\n", q, r, s, t,
                  tpdm[qrst]);
          fprintf(outfile, "contrib = %12.6lf, sum = %12.6lf\n", val, sum);
        }
        sum += val;
      }
    }
  }

}


}} // end namespace psi::detcas

