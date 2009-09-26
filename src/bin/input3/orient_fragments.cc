/*! \file 
    \ingroup (ORIENT_FRAGMENTS)
    \brief arrange multiple fragments based on interfragment coordinates 
*/
#define EXTERN
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libipv1/ip_lib.h>
#include <libchkpt/chkpt.h>
#include <cmath>
#include "input.h"
#include "global.h"
#include "defines.h"
#include "physconst.h"

namespace psi { namespace input {

void orient_fragments()
{
  int i, j, errcod, f, ival, pts, xyz, frag_dist_rho=0;
  int *P; // number of reference points in each fragment
  char ref_pts_lbl[80], error_message[80], frag_lbl[80];
  double R_AB, theta_A, theta_B, phi_A, phi_B, tau, tval, norm;
  double **geom_A, **geom_B;

  ip_boolean("FRAGMENT_DISTANCE_INVERSE", &(frag_dist_rho),0);

  /* a reference point for a fragment is a linear combination of atom positions - each
    fragment needs 3 reference points, unless the fragment has 1 atoms, 2 atoms, or is linear */
  P = (int *) malloc(nfragments*sizeof(int));
  for (f=0; f<nfragments; ++f) {
    if (frag_num_atoms[f]==1)      P[f] = 1;
    else if (frag_num_atoms[f]==2) P[f] = 2;
    else P[f] = 3;
    ref_pts_lc[f] = block_matrix(P[f],frag_num_atoms[f]);
  }

  for (f=0; f<nfragments; ++f) {
    if (f == 0)
      sprintf(ref_pts_lbl,"GEOMETRY_REF_POINTS");
    else
      sprintf(ref_pts_lbl,"GEOMETRY%d_REF_POINTS",f+1);

    /* Read in reference points linear combinations - normalize them if necessary */
    errcod = ip_count(ref_pts_lbl, &ival, 0);
    if (ival < P[f]) {
      fprintf(outfile,"Caution: Too few reference pts for full geometry specification, unless fragment %d is non-linear\n",f+1);
      P[f] = ival;
    }
    else if (ival > P[f]) {
      sprintf(error_message,"%s requires %d reference points!", ref_pts_lbl, P[f]);
      punt(error_message);
    }
    for (i=0; i<P[f]; ++i) {
      errcod = ip_count(ref_pts_lbl, &ival, 1, i);
      if (ival != frag_num_atoms[f]) {
        sprintf(error_message,"Fragment: %d, Reference point: %d requires %d entries!", f+1, i+1, frag_num_atoms[f]);
        punt(error_message);
      }
      for (j=0; j<frag_num_atoms[f]; ++j)
        ip_data(ref_pts_lbl, "%lf", &ref_pts_lc[f][i][j], 2, i, j);
    }
    for (i=0; i<P[f]; ++i) {
      tval = 0.0;
      for (j=0; j<frag_num_atoms[f]; ++j)
        tval += fabs(ref_pts_lc[f][i][j]);
      tval = 1.0/tval;
      for (j=0; j<frag_num_atoms[f]; ++j)
        ref_pts_lc[f][i][j] *= tval;
    }
    fprintf(outfile,"Linear combinations which specify reference points for fragment %d\n", f);
    print_mat(ref_pts_lc[f],P[f],frag_num_atoms[f],outfile);
  }
 
  for (f=1; f<nfragments; ++f) { /* leave first fragment alone */

    /* read in interfragment coordinates A->B; B will be moved to fit these */
    sprintf(frag_lbl, "FRAGMENTS%d%d", f, f+1);
    if (!ip_exist(frag_lbl,0))
      punt("You need to specify interfragment coordinates.");
    ip_count(frag_lbl, &ival, 0);
    if (ival != 6)
      punt("Interfragment coordinates should be given in sets of 6 (add 0's if necessary)");
    ip_data(frag_lbl, "%lf", &R_AB   , 1, 0);
    if (frag_dist_rho)
      R_AB = 1.0/R_AB * conv_factor; /* convert to au - allow 1/R coordinate */
    else 
      R_AB = R_AB * conv_factor;
    ip_data(frag_lbl, "%lf", &theta_A, 1, 1);
    ip_data(frag_lbl, "%lf", &theta_B, 1, 2);
    ip_data(frag_lbl, "%lf", &tau    , 1, 3);
    ip_data(frag_lbl, "%lf", &phi_A  , 1, 4);
    ip_data(frag_lbl, "%lf", &phi_B  , 1, 5);
  fprintf(outfile,"\nGiven %d-%d interfragment coordinates:\n", f, f+1);
  fprintf(outfile,"\t(1/)R_AB:%10.5f, theta_A:%10.5f, theta_B:%10.5f\n", R_AB, theta_A, theta_B);
  fprintf(outfile,"\t     tau:%10.5f,   phi_A:%10.5f,   phi_B:%10.5f\n", tau, phi_A, phi_B);

    /* put geometries of fragments in clean, contiguous blocks */
    geom_A = block_matrix(frag_num_atoms[f-1],3);
    geom_B = block_matrix(frag_num_atoms[f],3);
    for (xyz=0; xyz<3; ++xyz) {
      for (i=0; i<frag_num_atoms[f-1];++i)
        geom_A[i][xyz] = geometry[frag_atom[f-1]+i][xyz];
      for (i=0; i<frag_num_atoms[f];++i)
        geom_B[i][xyz] = geometry[frag_atom[f]+i][xyz];
    }

    /* function places geom_B in system of A with interfragment coordinates */
    orient_fragment(frag_num_atoms[f-1], frag_num_atoms[f], P[f-1], P[f], geom_A, geom_B,
    ref_pts_lc[f-1], ref_pts_lc[f], R_AB, theta_A, theta_B, tau, phi_A, phi_B, outfile);

    for (i=0;i<frag_num_atoms[f];++i)
      for (xyz=0; xyz<3; ++xyz)
        geometry[frag_atom[f]+i][xyz] = geom_B[i][xyz];

    free_block(geom_A);
    free_block(geom_B);
  }

  nref_per_fragment = P;
  return;
}


}} // namespace psi::input
