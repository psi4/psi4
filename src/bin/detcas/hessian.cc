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

/*
** form_appx_diag_mo_hess
**
** Calculates the approximate diagonal MO Hessian according to
** G. Chaban, M. W. Schmidt, and M. S. Gordon, Theor. Chem. Acc.,
** 97, 88-95 (1997)
** 
** Also uses active-active formulae from G. Chaban, personal communication
** (derived from H. J. Aa. Jensen and H. Agren, Chem. Phys. 104, 229 (1986))
**
** I am assuming that the pairs (p,q) are always given such that p>=q
**
** C. David Sherrill
** April 1998
** 
** Updated with active-active parts, March 2004
*/
void form_appx_diag_mo_hess(int npairs, int *ppair, int *qpair, 
  double *F_core, double *tei, double **opdm, double *tpdm, double *F_act,
  int firstact, int lastact, double *hess)
{

  int pair, p, q, pq, pp, qq;
  int ia;
  int i,ii,a,aa,t,tt,u,tu,v,w,vw,uv,tuvw;
  int qv,puvw,quvw,pv,pu,qu,pupv,quqv,ppuv,qquv,pvqu,puqv,pquv;
  double value;

  fprintf(outfile, "Forming approximate diagonal orbital Hessian\n");

  /* loop over the independent pairs */
  for (pair=0; pair<npairs; pair++) {
    p = ppair[pair];
    q = qpair[pair];
    pq = ioff[p] + q;
    pp = ioff[p] + p;
    qq = ioff[q] + q;
  
    /* H_{ai,ai}, i.e., inactive virt/inactive occ */
    if (p >= lastact && q < firstact) {
      hess[pair] = 4.0 * (F_core[pp] + F_act[pp] - F_core[qq] - F_act[qq]);
    }

    /* H_{at,at}, i.e., inactive virt with active orb */
    else if (p >= lastact && q >= firstact) {
      a = p;  t = q;
      aa = ioff[a] + a;

      hess[pair] = 2.0 * opdm[t][t] * (F_core[aa] + F_act[aa]);

      value = 0.0;
      for (u=firstact; u<lastact; u++) {
        tu = INDEX(t,u);
        value += opdm[t][u] * F_core[tu];
      }
      hess[pair] -= 2.0 * value;

      value = 0.0;
      for (u=firstact; u<lastact; u++) {
        tu = INDEX(t,u);
 
        /* loop over active v,w ... later restrict to symmetry allowed */
        for (v=firstact; v<lastact; v++) {
          for (w=firstact; w<lastact; w++) {
            vw = INDEX(v,w);
            tuvw = INDEX(tu,vw); 
            value += tpdm[tuvw] * tei[tuvw];
          }
        }
         
      }
      hess[pair] -= 2.0 * value; 
       
    } 

    /* H_{ti,ti}, i.e., active orb with inactive occ */
    else if (p >= firstact && q < firstact) {
      t = p;  i = q;
      tt = ioff[t] + t;
      ii = ioff[i] + i; 
      
      hess[pair] = 2.0 * opdm[t][t] * (F_core[ii] + F_act[ii]);
      hess[pair] += 4.0 * (F_core[tt] + F_act[tt] - F_core[ii] - F_act[ii]);

      value = 0.0;
      for (u=firstact; u<lastact; u++) {
        tu = INDEX(t,u);
        value += opdm[t][u] * F_core[tu];
      }
      hess[pair] -= 2.0 * value;

      value = 0.0;
      for (u=firstact; u<lastact; u++) {
        tu = INDEX(t,u);
 
        /* loop over active v,w ... later restrict to symmetry allowed */
        for (v=firstact; v<lastact; v++) {
          for (w=firstact; w<lastact; w++) {
            vw = INDEX(v,w);
            tuvw = INDEX(tu,vw); 
            value += tpdm[tuvw] * tei[tuvw];
          }
        }
         
      }
      hess[pair] -= 2.0 * value; 
    }

    /* 
       H_{t1t2,t1t2}, i.e., active-active (happens for MCSCF, not CASSCF)
       expression given by personal communication from Galina Chaban,
       derived from H. J. Aa. Jensen, H. Agren, Chem. Phys. 104, 229 (1982),
       appendix A.
    */
    else if (p >= firstact && q < lastact) {
      hess[pair] = 2.0 * (  opdm[q][q] * F_core[pp]   
                          + opdm[p][p] * F_core[qq]
                          - 2.0 * opdm[p][q] * F_core[pq] );
      value = 0.0;
      for (u=firstact; u<lastact; u++) {
        pu = INDEX(p,u);
        qu = INDEX(q,u);
        value += opdm[p][u] * F_core[pu] + opdm[q][u] * F_core[qu];
        /* should be able to reduce work below by factor of 2 */
        for (w=firstact; w<lastact; w++) {
          for (v=firstact; v<lastact; v++) {
          vw = INDEX(v,w);
          puvw = INDEX(pu,vw);
          quvw = INDEX(qu,vw);
          value += tpdm[puvw] * tei[puvw];
          value += tpdm[quvw] * tei[quvw]; 
          }
        }
      } 
      hess[pair] -= value * 2.0;

      value = 0.0;
      for (u=firstact; u<lastact; u++) {
        pu = INDEX(p,u);
        qu = INDEX(q,u);
        for (v=firstact; v<lastact; v++) {
          pv = INDEX(p,v);
          qv = INDEX(q,v);
          uv = INDEX(u,v);
          pupv = INDEX(pu,pv);
          quqv = INDEX(qu,qv);
          value += 2.0 * tpdm[pupv] * tei[quqv];
          value += 2.0 * tpdm[quqv] * tei[pupv]; 
          ppuv = INDEX(pp,uv);
          qquv = INDEX(qq,uv);
          value += tpdm[ppuv] * tei[qquv] + tpdm[qquv] * tei[ppuv];
          pvqu = INDEX(pv,qu);
          puqv = INDEX(pu,qv);
          value -= 4.0 * tpdm[pvqu] * tei[puqv];
          pquv = INDEX(pq,uv);
          value -= 2.0 * tpdm[pquv] * tei[pquv];
        }
      }
      hess[pair] += 2.0 * value;
      if (Params.scale_act_act != 1.0) hess[pair] *= Params.scale_act_act;
    } /* end case H_{t1t2,t1t2} */


    else {
      fprintf(outfile, 
             "(form_diag_mo_hess): Error, unrecognized class of indep pair\n");
    }

  } /* end loop over pairs */

}



/*
** form_diag_mo_hess
**
** Calculates the exact diagonal MO Hessian assuming CASSCF
**
** See Galina Chaban et al., Theor Chem Acc (1997) 97:88-95
** and references therein for mathematics of Newton-Raphson
** approach to MCSCF/CASSCF
**
** Supplemented with active-active parts from G. Chaban, personal
** communication
**
** I am assuming that the pairs (p,q) are always given such that p>=q
**
** G. O. Hyde
** January 2002
**
** Active-active parts added by C. D. Sherrill, March 2004
*/
void form_diag_mo_hess(int npairs, int *ppair, int *qpair, 
  double *F_core, double *tei, double **opdm, double *tpdm, double *F_act,
  int firstact, int lastact, double *hess)
{

  int pair, p, q, pq, pp, qq, pqpq, ppqq;
  int i,ii,a,aa,t,tt,u,tu,v,w,vw,tuvw;
  int au, uv, av, tv, ttuv, aauv, tvtu, avau;
  int ui, vi, ti, uvii, uivi, uiti, tuii;
  int qv,puvw,quvw,pv,pu,qu,pupv,quqv,ppuv,qquv,pvqu,puqv,pquv;
  int delta;
  double value;

  fprintf(outfile, "Forming diagonal orbital Hessian\n");

  /* loop over the independent pairs */
  for (pair=0; pair<npairs; pair++) {
    p = ppair[pair];
    q = qpair[pair];
    pq = ioff[p] + q;
    pp = ioff[p] + p;
    qq = ioff[q] + q;
    pqpq = ioff[pq] + pq;
    ppqq = ioff[pp] + qq; 
 
    /* H_{ai,ai}, i.e., inactive virt/inactive occ */
    if (p >= lastact && q < firstact) {
      hess[pair] = 4.0 * (F_core[pp] + F_act[pp] - F_core[qq] - F_act[qq]
                          + 3.0 * tei[pqpq] - tei[ppqq]);
    }

    /* H_{at,at}, i.e., inactive virt with active orb */
    else if (p >= lastact && q >= firstact) {
      a = p;  t = q;
      aa = ioff[a] + a;
      tt = ioff[t] + t;

      hess[pair] = 2.0 * opdm[t][t] * (F_core[aa]);

      value = 0.0;
      for (u=firstact; u<lastact; u++) {
        tu = INDEX(t,u);
        value += opdm[t][u] * F_core[tu];
      }
      hess[pair] -= 2.0 * value;

      value = 0.0;
      for (u=firstact; u<lastact; u++) {
        tu = INDEX(t,u);
 
        /* loop over active v,w ... later restrict to symmetry allowed */
        for (v=firstact; v<lastact; v++) {
          for (w=firstact; w<lastact; w++) {
            vw = INDEX(v,w);
            tuvw = INDEX(tu,vw); 
            value += tpdm[tuvw] * tei[tuvw];
          }
        }
         
      }
      hess[pair] -= 2.0 * value; 

      value = 0.0;
      for (u=firstact; u<lastact; u++) {
        tu = INDEX(t,u);
        au = INDEX(a,u);
        for (v=firstact; v<lastact; v++) {
          uv = INDEX(u,v);
          av = INDEX(a,v);
          tv = INDEX(t,v);
          ttuv = INDEX(tt,uv);
          aauv = INDEX(aa,uv);
          tvtu = INDEX(tv,tu);
          avau = INDEX(av,au);
          value += ( ( tpdm[ttuv] * tei[aauv] ) + 
                     ( 2.0 * tpdm[tvtu] * tei[avau] ) );
        }
      }
      hess[pair] += 2.0 * value;
 
    } 

    /* H_{ti,ti}, i.e., active orb with inactive occ */
    else if (p >= firstact && q < firstact) {
      t = p;  i = q;
      tt = ioff[t] + t;
      ii = ioff[i] + i; 
      ti = ioff[t] + i;
      
      hess[pair] = 2.0 * opdm[t][t] * (F_core[ii]);
      hess[pair] += 4.0 * (F_core[tt] + F_act[tt] - F_core[ii] - F_act[ii]);

      value = 0.0;
      for (u=firstact; u<lastact; u++) {
        tu = INDEX(t,u);
        value += opdm[t][u] * F_core[tu];
      }
      hess[pair] -= 2.0 * value;

      value = 0.0;
      for (u=firstact; u<lastact; u++) {
        tu = INDEX(t,u);
 
        /* loop over active v,w ... later restrict to symmetry allowed */
        for (v=firstact; v<lastact; v++) {
          for (w=firstact; w<lastact; w++) {
            vw = INDEX(v,w);
            tuvw = INDEX(tu,vw); 
            value += tpdm[tuvw] * tei[tuvw];
          }
        }
         
      }
      hess[pair] -= 2.0 * value; 

      value = 0.0;
      for (u=firstact; u<lastact; u++) {
        tu = INDEX(t,u);
        ui = INDEX(u,i);
        for (v=firstact; v<lastact; v++) {
          uv = INDEX(u,v);
          tv = INDEX(t,v);
          vi = INDEX(v,i);
          ttuv = INDEX(tt,uv);
          uvii = INDEX(uv,ii);
          tvtu = INDEX(tv,tu);
          uivi = INDEX(ui,vi);
          value += ( ( tpdm[ttuv] * tei[uvii] ) + 
                     ( 2.0 * tpdm[tvtu] * tei[uivi] ) );
        }
      }
      hess[pair] += 2.0 * value;

      value = 0.0;
      for (u=firstact; u<lastact; u++) {
        tu = INDEX(t,u);
        ui = INDEX(u,i);
        uiti = INDEX(ui,ti);
        tuii = INDEX(tu,ii);
          /* Create delta(t,u) */
          delta = 0.0;
          if( t == u )
            delta = 1.0;
          else
            delta = 0.0;
        value += ( ( delta - opdm[t][u] ) *
                   ( ( 3.0 * tei[uiti] ) - tei[tuii] ) );
      }
      hess[pair] += 4.0 * value;
    }

    /* 
       H_{t1t2,t1t2}, i.e., active-active (happens for MCSCF, not CASSCF)
       expression given by personal communication from Galina Chaban,
       derived from H. J. Aa. Jensen, H. Agren, Chem. Phys. 104, 229 (1982),
       appendix A.
    */
    else if (p >= firstact && q < lastact) {
      hess[pair] = 2.0 * (  opdm[q][q] * F_core[pp]   
                          + opdm[p][p] * F_core[qq]
                          - 2.0 * opdm[p][q] * F_core[pq] );
      value = 0.0;
      for (u=firstact; u<lastact; u++) {
        pu = INDEX(p,u);
        qu = INDEX(q,u);
        value += opdm[p][u] * F_core[pu] + opdm[q][u] * F_core[qu];
        /* should be able to reduce work below by factor of 2 */
        for (w=firstact; w<lastact; w++) {
          for (v=firstact; v<lastact; v++) {
          vw = INDEX(v,w);
          puvw = INDEX(pu,vw);
          quvw = INDEX(qu,vw);
          value += tpdm[puvw] * tei[puvw];
          value += tpdm[quvw] * tei[quvw]; 
          }
        }
      } 
      hess[pair] -= value * 2.0;

      value = 0.0;
      for (u=firstact; u<lastact; u++) {
        pu = INDEX(p,u);
        qu = INDEX(q,u);
        for (v=firstact; v<lastact; v++) {
          pv = INDEX(p,v);
          qv = INDEX(q,v);
          uv = INDEX(u,v);
          pupv = INDEX(pu,pv);
          quqv = INDEX(qu,qv);
          value += 2.0 * tpdm[pupv] * tei[quqv];
          value += 2.0 * tpdm[quqv] * tei[pupv]; 
          ppuv = INDEX(pp,uv);
          qquv = INDEX(qq,uv);
          value += tpdm[ppuv] * tei[qquv] + tpdm[qquv] * tei[ppuv];
          pvqu = INDEX(pv,qu);
          puqv = INDEX(pu,qv);
          value -= 4.0 * tpdm[pvqu] * tei[puqv];
          pquv = INDEX(pq,uv);
          value -= 2.0 * tpdm[pquv] * tei[pquv];
        }
      }
      hess[pair] += 2.0 * value;
      if (Params.scale_act_act != 1.0) hess[pair] *= Params.scale_act_act;
    } /* end case H_{t1t2,t1t2} */

    else {
      fprintf(outfile, 
             "(form_diag_mo_hess): Error, unrecognized class of indep pair\n");
    }

  } /* end loop over pairs */

}

/*
** form_full_mo_hess
**
** Calculates the full MO Hessian d^2E / (Theta_{pq} Theta_{rs})
** for independent pairs (p,q) and (r,s).  Assume independent pairs
** are coming in such that p>=q.
**
** d^2 E / (Theta_{pq} Theta_{rs}) = 
**   y_{pqrs} - y_{qprs} - y_{pqsr} + y_{qpsr}
**   + 1/2 delta_{ps} ( x_{qr} + x_{rq} )
**   - 1/2 delta_{qs} ( x_{pr} + x_{rp} )
**   - 1/2 delta_{pr} ( x_{qs} + x_{sq} )
**   + 1/2 delta_{qr} ( x_{ps} + x_{sp} )
**
** x is the Lagrangian
** y_{pqrs} = gamma_{qs} 
**          + \sum_{mn} ( 2 Gamma_{qsmn} (pr|mn) + 4 Gamma_{qmsn} (pm|rn) )
**
** note: indices q and s must be populated for y to be nonzero
**
** Based on notes by Yukio Yamaguchi (MCSCF eq. 22)
** We need twice his value to be consistent with the "actual" gradient
** (his eq. 21, not his eq. 27) or to match G. Chaban, M. W. Schmidt,
** and M. S. Gordon, Theor. Chem. Acc. 97, 88 (1997).
**
** C. David Sherrill
** September 2003
*/
void form_full_mo_hess(int npairs, int *ppair, int *qpair, double *oei,
  double *tei, double **opdm, double *tpdm, double **lag,
  double **hess)
{
  int npop, nmo;
  int pair1, pair2, p, q, r, s, m, n;
  int qs, ps, pr, qr, mn, qm, sn, pm, rn;
  int qsmn, prmn, qmsn, pmrn, psmn, pmsn, qrmn, qmrn;  
  double sum, ysum;

  fprintf(outfile, "Forming full MCSCF orbital Hessian\n");
 
  nmo = CalcInfo.nmo;
  npop = CalcInfo.npop;

  /* loop over the pairs of independent pairs */
  for (pair1=0; pair1<npairs; pair1++) {
    p = ppair[pair1];
    q = qpair[pair1];

    for (pair2=0; pair2<=pair1; pair2++) {
      r = ppair[pair2];
      s = qpair[pair2];

      /* compute the hessian contribution for p,q,r,s */
      sum = 0.0;
      if (p == s) sum += 0.5 * (lag[q][r] + lag[r][q]);
      if (q == s) sum -= 0.5 * (lag[p][r] + lag[r][p]);
      if (p == r) sum -= 0.5 * (lag[q][s] + lag[s][q]);
      if (q == r) sum += 0.5 * (lag[p][s] + lag[s][p]);

      pr = INDEX(p,r);
      qs = INDEX(q,s);
      ps = INDEX(p,s);
      qr = INDEX(q,r);

      /* compute y_{pqrs} */
      ysum = 0.0;
      if (q < npop && s < npop) {
        ysum += oei[pr] * opdm[q][s];
        for (m=0; m<npop; m++) {
          qm = INDEX(q,m);
          pm = INDEX(p,m);
          for (n=0; n<npop; n++) {
            mn = INDEX(m,n);
            sn = INDEX(s,n);
            rn = INDEX(r,n);
            qsmn = INDEX(qs,mn);
            qmsn = INDEX(qm,sn);
            prmn = INDEX(pr,mn);
            pmrn = INDEX(pm,rn); 
            ysum += tpdm[qsmn] * tei[prmn];
            ysum += 2.0 * tpdm[qmsn] * tei[pmrn];
          }
        }
      }
      sum += ysum;

      /* compute y_{qprs} */
      ysum = 0.0;
      if (p < npop && s < npop) {
        ysum += oei[qr] * opdm[p][s];
        for (m=0; m<npop; m++) {
          pm = INDEX(p,m);
          qm = INDEX(q,m);
          for (n=0; n<npop; n++) {
            mn = INDEX(m,n);
            sn = INDEX(s,n);
            rn = INDEX(r,n);
            psmn = INDEX(ps,mn);
            pmsn = INDEX(pm,sn);
            qrmn = INDEX(qr,mn);
            qmrn = INDEX(qm,rn); 
            ysum += tpdm[psmn] * tei[qrmn];
            ysum += 2.0 * tpdm[pmsn] * tei[qmrn];
          }
        }
      }
      sum -= ysum;

      /* compute y_{pqsr} */
      ysum = 0.0;
      if (q < npop && r < npop) {
        ysum += oei[ps] * opdm[q][r];
        for (m=0; m<npop; m++) {
          qm = INDEX(q,m);
          pm = INDEX(p,m);
          for (n=0; n<npop; n++) {
            mn = INDEX(m,n);
            rn = INDEX(r,n);
            sn = INDEX(s,n);
            qrmn = INDEX(qr,mn);
            qmrn = INDEX(qm,rn);
            psmn = INDEX(ps,mn);
            pmsn = INDEX(pm,sn); 
            ysum += tpdm[qrmn] * tei[psmn];
            ysum += 2.0 * tpdm[qmrn] * tei[pmsn];
          }
        }
      }
      sum -= ysum;

      /* compute y_{qpsr} */
      ysum = 0.0;
      if (p < npop && r < npop) {
        ysum += oei[qs] * opdm[p][r];
        for (m=0; m<npop; m++) {
          pm = INDEX(p,m);
          qm = INDEX(q,m);
          for (n=0; n<npop; n++) {
            mn = INDEX(m,n);
            rn = INDEX(r,n);
            sn = INDEX(s,n);
            prmn = INDEX(pr,mn);
            pmrn = INDEX(pm,rn);
            qsmn = INDEX(qs,mn);
            qmsn = INDEX(qm,sn); 
            ysum += tpdm[prmn] * tei[qsmn];
            ysum += 2.0 * tpdm[pmrn] * tei[qmsn];
          }
        }
      }
      sum += ysum;

      hess[pair1][pair2] = 2.0 * sum;
      hess[pair2][pair1] = 2.0 * sum;
    } /* end loop over pair2 */
  } /* end loop over pair1 */

}

/*
** form_diag_mo_hess_yy
**
** Calculates the diagonal MO Hessian d^2E / Theta_{pq}^2
** for independent pair (p,q).  Assume independent pairs
** are coming in such that p>=q.
**
** d^2 E / Theta_{pq}^2
**   y_{pqpq} - y_{qppq} - y_{pqqp} + y_{qpqp}
**   - 1/2 delta_{qq} ( x_{pp} + x_{pp} )
**   - 1/2 delta_{pp} ( x_{qq} + x_{qq} )
**
** x is the Lagrangian
** y_{pqrs} = gamma_{qs} 
**          + \sum_{mn} ( 2 Gamma_{qsmn} (pr|mn) + 4 Gamma_{qmsn} (pm|rn) )
**
** note: indices q and p must be populated for y to be nonzero
**
** Based on notes by Yukio Yamaguchi (MCSCF eq. 22)
**
** C. David Sherrill
** September 2003
*/
void form_diag_mo_hess_yy(int npairs, int *ppair, int *qpair, double *oei,
  double *tei, double **opdm, double *tpdm, double **lag,
  double *hess)
{
  int npop, nmo;
  int pair1, pair2, p, q, m, n;
  int pq, pp, qq, mn, pn, qn, pm, qm;
  int qqmn, ppmn, qmqn, pmpn, pqmn, qpmn, pmqn, qmpn;
  double sum, ysum;
 
  fprintf(outfile, "Forming diagonal MCSCF orbital Hessian (YY)\n");

  nmo = CalcInfo.nmo;
  npop = CalcInfo.npop;

  /* loop over the pairs of independent pairs */
  for (pair1=0; pair1<npairs; pair1++) {
    p = ppair[pair1];
    q = qpair[pair1];

    pq = INDEX(p,q);
    pp = ioff[p] + p;
    qq = ioff[q] + q;

    /* compute the hessian contribution for p,q,p,q */
    sum = - lag[p][p] - lag[q][q];

    /* compute y_{pqpq} */
    ysum = 0.0;
    if (q < npop) {
      ysum += oei[pp] * opdm[q][q];
       
      for (m=0; m<npop; m++) {
        pm = INDEX(p,m);
        qm = INDEX(q,m);
        for (n=0; n<npop; n++) {
          mn = INDEX(m,n);
          pn = INDEX(p,n);
          qn = INDEX(q,n);
          qqmn = INDEX(qq,mn);
          ppmn = INDEX(pp,mn);
          qmqn = INDEX(qm,qn);
          pmpn = INDEX(pm,pn);
          ysum += tpdm[qqmn] * tei[ppmn];
          ysum += 2.0 * tpdm[qmqn] * tei[pmpn];
        }
      }
    }
    sum += ysum;

    /* compute y_{qppq} = y_{pqqp} */
    ysum = 0.0;
    if (p < npop && q < npop) {
      ysum += oei[pq] * opdm[p][q];
      for (m=0; m<npop; m++) {
        pm = INDEX(p,m);
        qm = INDEX(q,m);
        for (n=0; n<npop; n++) {
          mn = INDEX(m,n);
          qn = INDEX(q,n);
          pn = INDEX(p,n);
          pqmn = INDEX(pq,mn);
          pmqn = INDEX(pm,qn);
          qmpn = INDEX(qm,pn); 

          ysum += tpdm[pqmn] * tei[pqmn];
          ysum += 2.0 * tpdm[pmqn] * tei[qmpn];
        }
      }
    }
    sum -= 2.0 * ysum;

    /* compute y_{qpqp} */
    ysum = 0.0;
    if (p < npop) {
      ysum += oei[qq] * opdm[p][p];
      for (m=0; m<npop; m++) {
        pm = INDEX(p,m);
        qm = INDEX(q,m);
        for (n=0; n<npop; n++) {
          mn = INDEX(m,n);
          pn = INDEX(p,n);
          qn = INDEX(q,n);
          ppmn = INDEX(pp,mn);
          qqmn = INDEX(qq,mn);
          pmpn = INDEX(pm,pn);
          qmqn = INDEX(qm,qn);

          ysum += tpdm[ppmn] * tei[qqmn];
          ysum += 2.0 * tpdm[pmpn] * tei[qmqn];
        }
      }
    }
    sum += ysum;

    hess[pair1] = 2.0 * sum;
  } /* end loop over pair1 */

}

}} // end namespace psi::detcas

