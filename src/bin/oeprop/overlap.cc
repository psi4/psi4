/*! \file
    \ingroup OEPROP
    \brief Enter brief description of file here 
*/
#define EXTERN
#include "includes.h"
#include "globals.h"
#include "prototypes.h"

namespace psi { namespace oeprop {

void compute_overlap()
{
  int i,j,k,l,ig,jg,ibf,jbf,ib,jb,jlim,kk,ll;
  int iatom, jatom, iang, jang, i_1stbf, j_1stbf, nbfi, nbfj;
  int igmin, igmax, jgmin, jgmax;
  int ixm, iym, izm, jxm, jym, jzm, iind, jind;
  int eq_shell, eq_atoms;
  double ax, ay, az, bx, by, bz, xab, yab, zab, rab2;
  double ai, bj, gamma, over_pf, tmp;
  double px, py, pz, pax, pay, paz, pbx, pby, pbz, pcx, pcy, pcz;
  double x0, y0, z0;
  

		/* Initialize some intermediates */
  
  S = init_array(natri);
  MIX = init_box(lmax+1,lmax+1,1);
  MIY = init_box(lmax+1,lmax+1,1);
  MIZ = init_box(lmax+1,lmax+1,1);
  MIX[0][0][0] = MIY[0][0][0] = MIZ[0][0][0] = 1.0;

	/* I-shell loop */

  for(i=0;i<nshell;i++) {
    iatom = snuc[i] - 1;
    ax = geom[iatom][0];
    ay = geom[iatom][1];
    az = geom[iatom][2];
    iang = stype[i]-1;
    izm = 1;
    iym = iang+1;
    ixm = iym*iym;
    i_1stbf = sloc[i];
    nbfi = (iang+2)*(iang+1)/2;
    igmin = sprim[i] - 1;
    igmax = igmin + snumg[i] - 1;
    
    
    	/* J-shell loop */
    
    for(j=0;j<=i;j++) {
      jatom = snuc[j] - 1;
      bx = geom[jatom][0];
      by = geom[jatom][1];
      bz = geom[jatom][2];
      jang = stype[j]-1;
      jzm = 1;
      jym = jang+1;
      jxm = jym*jym;
      j_1stbf = sloc[j];
      nbfj = (jang+2)*(jang+1)/2;
      jgmin = sprim[j] - 1;
      jgmax = jgmin + snumg[j] - 1;
      
      eq_shell = (i == j);
      eq_atoms = (iatom == jatom);
      
      xab = ax - bx;
      yab = ay - by;
      zab = az - bz;
      rab2 = xab*xab + yab*yab + zab*zab;
      
      for(ig=igmin;ig<=igmax;ig++) {
        ai = exps[ig];
        if (eq_shell)
          jgmax = ig;
        for(jg=jgmin;jg<=jgmax;jg++) {
          bj = exps[jg];
          gamma = ai + bj;
          tmp = ai*bj*rab2/gamma;
          
          over_pf = contr[ig]*contr[jg]*exp(-tmp)*pow(M_PI/gamma,1.5);
          if (eq_shell && ig != jg)
            over_pf *= 2.0;
          
          
          /* Compute P and intermediates */
          
          px = (ax*ai + bx*bj)/gamma;
          py = (ay*ai + by*bj)/gamma;
          pz = (az*ai + bz*bj)/gamma;
          if (eq_atoms)
            pax = pay = paz = pbx = pby = pbz = 0.0;
          else {
            pax = px - ax;  pay = py - ay;  paz = pz - az;
            pbx = px - bx;  pby = py - by;  pbz = pz - bz;
          }

		/* Computing moment integrals */

          MI_OSrecurs(pax,pay,paz,pbx,pby,pbz,gamma,iang,jang,0);

		/* Computing elec. potential integrals -
		   see Obara and Saika paper */
          for(ibf=0;ibf<nbfi;ibf++) {
            if (eq_shell)
              jlim = ibf + 1;
            else
              jlim = nbfj;
            for(jbf=0;jbf<jlim;jbf++) {
              ib = i_1stbf + ibf - 1; jb = j_1stbf + jbf - 1;
              
              x0 = MIX[xpow_bf[iang][ibf]][xpow_bf[jang][jbf]][0];
              y0 = MIY[ypow_bf[iang][ibf]][ypow_bf[jang][jbf]][0];
              z0 = MIZ[zpow_bf[iang][ibf]][zpow_bf[jang][jbf]][0];
              
              
              S[ioff[ib]+jb] += over_pf*x0*y0*z0*norm_bf[iang][ibf]*norm_bf[jang][jbf];
              
            }
          }
        }
      }


    }
  }

  if (print_lvl >= PRINTOVERLAPLEVEL) {
    fprintf(outfile,"\tOverlap matrix in AO basis :\n");
    print_array(S,nbfao,outfile);
    fprintf(outfile,"\n\n");
  }

	/* Cleaning up */
  
  free_box(MIX,lmax+1,lmax+1);
  free_box(MIY,lmax+1,lmax+1);
  free_box(MIZ,lmax+1,lmax+1);
}

}} // namespace psi::oeprop
