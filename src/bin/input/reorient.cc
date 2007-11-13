/*! \file 
    \ingroup (INPUT)
    \brief Enter brief description of file here 
*/
#define EXTERN
#include <stdio.h>
#include <libciomr/libciomr.h>
#include <stdlib.h>
#include <math.h>
#include "input.h"
#include "global.h"
#include "defines.h"
#include <physconst.h>
#include <masses.h>

namespace {
  void if_to_invert_axis(double* v1, int* must_invert, int* should_invert, double* maxproj);
}

namespace psi { namespace input {

void reorient()
{
  int i,j;
  int degen;
  int deg_IM1, deg_IM2;
  int nspher_set;
  int prox_i, prox_j;
  int axis, axis1, axis2, xyz;
  int nshould, nmust, must_invert[3], should_invert[3];
  double Xcm = 0.0;
  double Ycm = 0.0;
  double Zcm = 0.0;
  double mass = 0.0;
  double tmp, abs, rel;
  double Zmax, r, min_ij;
  double **IT, *IM, **ITAxes;
  double *v1, *v2, *v3, **R;
  double vabs, maxproj[3];
  double **sset_geom, *sset_dist;
  double origin[] = {0.0, 0.0, 0.0};

  double cos_theta, sin_theta, theta;
  double cos_phix, cos_phiy, phix, sin_phix;
  double v2norm;
  double **Rz, **Rx, **Rzt;

  const double im2rotconst = 0.25/(M_PI*_c_au*_bohr2cm);

    v1 = init_array(3);
    v2 = init_array(3);
    v3 = init_array(3);
    IT = block_matrix(3,3);
    IM = init_array(3);
    ITAxes = block_matrix(3,3);

    for(i=0;i<num_atoms;i++) {
      tmp = an2masses[(int)nuclear_charges[i]];
      Xcm += tmp*geometry[i][0];
      Ycm += tmp*geometry[i][1];
      Zcm += tmp*geometry[i][2];
      mass += tmp;
    }
    Xcm /= mass; Ycm /= mass; Zcm /= mass;

    if (!no_comshift) {
      /*full geom center-of-mass shift*/
      for(i=0;i<num_allatoms;i++) {
	full_geom[i][0] -= Xcm;
	full_geom[i][1] -= Ycm;
	full_geom[i][2] -= Zcm;
      }
    }


    /*Computing inertia tensor, moments of inertia, and principal axes*/
    if (num_atoms > 1) {
      for(i=0;i<num_atoms;i++) {
        tmp = an2masses[(int)nuclear_charges[i]]/_au2amu;
        IT[0][0] += tmp*(geometry[i][1]*geometry[i][1] + geometry[i][2]*geometry[i][2]);
        IT[1][1] += tmp*(geometry[i][0]*geometry[i][0] + geometry[i][2]*geometry[i][2]);
        IT[2][2] += tmp*(geometry[i][0]*geometry[i][0] + geometry[i][1]*geometry[i][1]);
	IT[0][1] -= tmp*geometry[i][0]*geometry[i][1];
	IT[0][2] -= tmp*geometry[i][0]*geometry[i][2];
	IT[1][2] -= tmp*geometry[i][1]*geometry[i][2];
      }
      IT[1][0] = IT[0][1];
      IT[2][0] = IT[0][2];
      IT[2][1] = IT[1][2];
      sq_rsp(3,3,IT,IM,1,ITAxes,1.0E-14);
      IM[0] = fabs(IM[0]); /*Fixing numerical errors in the linear case*/
      fprintf(outfile,"\n  -Rotational constants (cm-1) :\n");
      if (IM[0] < ZERO_MOMENT_INERTIA) /* Linear molecule */
	fprintf(outfile,"    A = **********  ");
      else   /* Regular molecule */
	fprintf(outfile,"    A = %10.5lf  ",im2rotconst/IM[0]);
      if (IM[1] < ZERO_MOMENT_INERTIA)  /* Atom */
	fprintf(outfile,"B = **********  C = **********\n");
      else /* molecule */
	fprintf(outfile,"B = %10.5lf  C = %10.5lf\n",im2rotconst/IM[1],im2rotconst/IM[2]);

    
      /*Computing degeneracy*/
      degen = 0;
      for(i=0;i<2;i++)
	for(j=i+1;j<3 && degen<2;j++) {
	  abs = fabs(IM[i] - IM[j]);
	  tmp = (IM[i] > IM[j]) ? IM[i] : IM[j];
	  if (abs > 1.0E-14)
            rel = abs/tmp;
          else
	    rel = 0.0;
	  if (rel < ZERO_MOMENT_INERTIA) {
	    degen++;
	    deg_IM1 = i;
	    deg_IM2 = j;
	  }
	}

      /*Ensuring the righthandedness of the reference coordinate system*/
      v1[0] = ITAxes[0][1];
      v1[1] = ITAxes[1][1];
      v1[2] = ITAxes[2][1];
      v2[0] = ITAxes[0][2];
      v2[1] = ITAxes[1][2];
      v2[2] = ITAxes[2][2];
      cross_prod(v1,v2,v3);
      ITAxes[0][0] = v3[0];
      ITAxes[1][0] = v3[1];
      ITAxes[2][0] = v3[2];

      /*--------------------------------------------------------
	Check if the original frame is "sufficiently close"
	to the principal frame. If true, then try to minimize
	the necessary reorientation using heuristic rules:
	1. If an axis V has 2 zero cartesian projections
	   and the sign of the other projection is negative
	   the axis needs to be inverted
	2. In an axis V has 1 zero cartesian projection
	   and the largest (in absolute magnitude) of non-zero
	   projections is negative then the axis may need to
	   be inverted
	3. If an axis V has no zero cartesian projections
	   and the largest (in absolute magnitude) of non-zero
	   projections is negative then the axis may beed to
	   be inverted
	This procedure is needed to ensure that in
	finite difference calculation there's not too much
	reorientation (this is a hack and may need to be
	rethought later).
	-EFV
       -------------------------------------------------------*/
      nmust = 0;
      nshould = 0;
      for(axis=0; axis<3; axis++) {

        v1[0] = ITAxes[0][axis];
        v1[1] = ITAxes[1][axis];
        v1[2] = ITAxes[2][axis];

        if_to_invert_axis(v1,&(must_invert[axis]),&(should_invert[axis]),&(maxproj[axis]));
        nmust += must_invert[axis];
        nshould += should_invert[axis];
	
      }

      R = block_matrix(3,3);
      if (nmust == 2) {
	for(axis=0; axis<3; axis++) {
	  if (must_invert[axis])
	    R[axis][axis] = -1.0;
	  else
	    R[axis][axis] = 1.0;
	}
      }
      else if (nmust == 1 && nshould > 0) {
	if (nshould == 2) {
	  for(axis=0; axis<3; axis++)
	    if (should_invert[axis]) {
	      axis1 = axis;
	      axis++;
	      break;
	    }
	  for(; axis<3; axis++)
	    if (should_invert[axis]) {
	      axis2 = axis;
	      break;
	    }
	  if (fabs(maxproj[axis1]) > fabs(maxproj[axis2])) {
	    nshould = 1;
	    should_invert[axis2] = 0;
	  }
	  else {
	    nshould = 1;
	    should_invert[axis1] = 0;
	  }
	}

	for(axis=0; axis<3; axis++) {
	  if (must_invert[axis])
	    R[axis][axis] = -1.0;
	  else if (should_invert[axis])
	    R[axis][axis] = -1.0;
	  else
	    R[axis][axis] = 1.0;
	}
      }
      else if (nmust == 3) {
	R[0][0] = -1.0;
	R[1][1] = -1.0;
	R[2][2] = 1.0;
      }
      else if (nmust == 0 && nshould > 1) {
	if (nshould == 3) {
	  tmp = fabs(maxproj[0]);
	  i = 0;
	  for(axis=1; axis<3; axis++) {
	    if (fabs(maxproj[axis]) < fabs(tmp)) {
	      i = axis;
	      tmp = fabs(maxproj[axis]);
	    }
	  }
	  should_invert[i] = 0;
	  nshould = 2;
	}
	for(axis=0; axis<3; axis++) {
	  if (should_invert[axis])
	    R[axis][axis] = -1.0;
	  else
	    R[axis][axis] = 1.0;
	}
      }
      else {
	R[0][0] = 1.0;
	R[1][1] = 1.0;
	R[2][2] = 1.0;
      }

      /*Reorient if degen == 0 (asymmetric top case).
	Otherwise hope user knows what he/she's doing and leave it as it is.*/
      if (degen == 0 && !no_reorient) {
	rotate_full_geom(ITAxes);
       	rotate_full_geom(R);
      }
      free_block(R);
      
      /*-------------------------------------------------------------
	For linear or symmetric tops need to rotate so that
	the nondegenerate moment of inertia is along z:
	1. v1 - nondegenerate
	2. cos(Theta) = v1.dot.Z
	3. v2 = v1.cross.Z
	4. cos(Phi_x) = v2.dot.X
	5. cos(Phi_y) = v2.dot.Y
	6. From 4 and 5 determine quadrant
	7. Rotate around Z by Phi_x or -Phi_x so that v2 is along X
	8. Rotate around X by Theta
	9. Rotate around Z by -Phi_x or Phi_x

	Added 4/26/2003 by EFV.
       ------------------------------------------------------------*/
      if (degen == 1 && !no_reorient) {

        int must_invert, should_invert, unique_axis;
        double maxproj, invert_pfac;
        
	if (deg_IM1 + deg_IM2 == 3)
          unique_axis = 0;
        else
          unique_axis = 2;
          
        v1[0] = ITAxes[0][unique_axis];
        v1[1] = ITAxes[1][unique_axis];
        v1[2] = ITAxes[2][unique_axis];
        
        /* Figure out if the unique axis needs to be inverted
           If yes - just invert that axis directly, since other axes
           are not even used */
        if_to_invert_axis(v1, &must_invert, &should_invert, &maxproj);
        if (must_invert || should_invert)
          invert_pfac = -1.0;
        else
          invert_pfac = 1.0;

        v1[0] *= invert_pfac;
        v1[1] *= invert_pfac;
        v1[2] *= invert_pfac;
        
	cos_theta = v1[2];
	if ( (1.0 - fabs(cos_theta)) > ZERO_MOMENT_INERTIA) {
	  theta = acos(cos_theta);
	  sin_theta = sin(theta);
	  
	  v3[0] = 0.0;
	  v3[1] = 0.0;
	  v3[2] = 1.0;
	  cross_prod(v1,v3,v2);
	  v2norm = sqrt(v2[0] * v2[0] +
			v2[1] * v2[1] +
			v2[2] * v2[2]);
	  v2[0] /= v2norm;
	  v2[1] /= v2norm;
	  v2[2] /= v2norm;
	  
	  cos_phix = v2[0];
	  cos_phiy = v2[1];
	  phix = acos(cos_phix);
	  
	  if (cos_phiy > 0.0) {
	    phix *= -1.0; 
	  }
	  sin_phix = sin(phix);

	  Rz = block_matrix(3,3);
	  Rz[2][2] = 1.0;
	  Rz[0][0] = cos_phix;
	  Rz[1][1] = cos_phix;
	  Rz[0][1] = sin_phix;
	  Rz[1][0] = -sin_phix;
	  rotate_full_geom(Rz);
	  free_block(Rz);

	  Rx = block_matrix(3,3);
	  Rx[0][0] = 1.0;
	  Rx[1][1] = cos_theta;
	  Rx[2][2] = cos_theta;
	  Rx[1][2] = sin_theta;
	  Rx[2][1] = -sin_theta;
	  rotate_full_geom(Rx);
	  free_block(Rx);

	  Rzt = block_matrix(3,3);
	  Rzt[2][2] = 1.0;
	  Rzt[0][0] = cos_phix;
	  Rzt[1][1] = cos_phix;
	  Rzt[0][1] = -sin_phix;
	  Rzt[1][0] = sin_phix;
	  rotate_full_geom(Rzt);
	  free_block(Rzt);
	}
      }

      /*If degen=0 (asymmetric top) - do nothing
	   degen=1 (linear molecule or symmetric top) - do nothing
           denen=2 (atom or spherical top) - do lots of stuff - see below*/
      if (degen == 0) {
	fprintf(outfile,"    It is an asymmetric top.\n");
	rotor = asymmtop;
      }
      else if (degen == 1)
	switch (deg_IM1 + deg_IM2) {
	  case 3: /*B and C are degenerate - linear or prolate symm. top.*/
		  if (IM[0] < ZERO_MOMENT_INERTIA) {
		    fprintf(outfile,"    It is a linear molecule.\n");
		    rotor = linear;
		  }
		  else {
		    fprintf(outfile,"    It is a prolate symmetric top.\n");
		    rotor = symmtop;
		  }
		  break;
	  case 1: /*A and B are degenerate - oblate top.
		    C is the unique axis. Do nothing.*/
	          fprintf(outfile,"    It is an oblate symmetric top.\n");
		  rotor = symmtop;
	          break;
	}
      else if (degen == 2) { /*Man, this piece of code is nasty!!!*/
	fprintf(outfile,"    It is a spherical top.\n");
	rotor = sphtop;
	Zmax = 0;
	for(i=0;i<num_atoms;i++) /*Finding the heaviest type of atoms positioned NOT in the origin*/
	  if (Zmax < (int)nuclear_charges[i] && sqrt(dot_prod(geometry[i],geometry[i])) > ZERO)
	    Zmax = (int)nuclear_charges[i];
	/*Find subset of heaviest atoms at the same distance from the origin - so called spherical set
	  and store their geometries in sset_geom[][] */
	r = 0.0;
	for(i=0;i<num_atoms && r < ZERO;i++)
	  r = (nuclear_charges[i] == Zmax) ? sqrt(dot_prod(geometry[i],geometry[i])) : 0.0;
	sset_geom = init_matrix(num_atoms,3);
	sset_geom[0][0] = geometry[i-1][0];
	sset_geom[0][1] = geometry[i-1][1];
	sset_geom[0][2] = geometry[i-1][2];
	nspher_set = 1;
	for(j=i;j<num_atoms;j++)
	  if (nuclear_charges[j] == Zmax) {
	    tmp = sqrt(dot_prod(geometry[i],geometry[i]));
	    if ( fabs(tmp - r) < ZERO) {
	      sset_geom[nspher_set][0] = geometry[j][0];
	      sset_geom[nspher_set][1] = geometry[j][1];
	      sset_geom[nspher_set][2] = geometry[j][2];
	      nspher_set++;
	    }
	  }
	/*Find the unique orientation. If nspher_set = 4, 6, or 8 - we are dealing with a tetrahedron,
	  octahedron or cube respectively - and treat these special cases rather easily.
	  In general case, the procedure is as follows:
	  1. Pick an atom (1)
	  2. Pick an atom (2) not related to (1) by inversion
	  3. Loop over the spherical set and find all atoms {(3)}
	     so that dist(1)-(2) == dist(1)-(3)
	  4. switch(1+Number of atoms in {(3)})
	       1,2 - pick different (2)
	       3 - atom(1) is on a C3 axis and (2) and {(3)} are related by C3
	       4 - either atom(1) is on a C2 axis or there're two C3's next to each other
	       5 - find (4) in {(3)} that dist(1)-(4) == dist(1)-(2)
	           1, 2, 4 are related and perpendicular to C3
	  5. Do that again trying to find either 2 C3 axes or C3 and C2 or two C2's
	  6. Let's go ... */
        if (!no_reorient)
	switch(nspher_set) {
	  case 4: /*Tetrahedron*/
	          median_vec(sset_geom[0], sset_geom[1], v1);
		  median_vec(sset_geom[1], sset_geom[2], v2);
		  cross_prod(v1, v2, v3);
		  vectors_to_matrix(v1, v2, v3, ITAxes);
		  rotate_full_geom(ITAxes);
		  break;

	  case 6: /*Octahedron*/
  	          if (!inv_related(sset_geom[0], sset_geom[1]))
		    median_vec(sset_geom[0], sset_geom[1], v1);
		  else 
		    median_vec(sset_geom[0], sset_geom[2], v1);
		  unit_vec(sset_geom[0],origin,v2);
		  cross_prod(v1, v2, v3);
		  unit_vec(v3,origin,v3);
		  cross_prod(v3, v1, v2);
		  vectors_to_matrix(v1, v2, v3, ITAxes);
		  rotate_full_geom(ITAxes);
	          break;

	  case 8: /*Cube*/
	          sset_dist = init_array(nspher_set*(nspher_set+1)/2);
		  calc_distance(sset_geom,sset_dist,nspher_set);
	          min_ij = sset_dist[ioff[1]+0];
		  prox_i = 1;
		  for(i=2;i<nspher_set-2;i++)
		    if (min_ij > (tmp = sset_dist[ioff[i]+0]) && !inv_related(sset_geom[0],sset_geom[i])) {
		      min_ij = tmp;
		      prox_i = i;
		      break;
		    }
		  for(j=prox_i+1;j<nspher_set;j++)
		    if (fabs(min_ij - sset_dist[ioff[j]+0]) < ZERO) {
		      prox_j = j;
		      break;
		    }
		  unit_vec(sset_geom[0], sset_geom[prox_i], v1);
		  unit_vec(sset_geom[0], sset_geom[prox_j], v2);
		  cross_prod(v1, v2, v3);
		  vectors_to_matrix(v1, v2, v3, ITAxes);
		  rotate_full_geom(ITAxes);
		  free(sset_dist);
	          break;

	 default: /*General case*/
/*	          sset_dist = init_array(nspher_set*(nspher_set+1)/2);
		  calc_distance(sset_geom,sset_dist,nspher_set);
		  a1 = 0;
		  a2 = (inv_related(sset_geom[a1],sset_geom[a1+1])) ? a1+2 : a1+1;
		  r = sset_dist[ioff[a2]+i1];
		  num_a3 = 0;
		  for(i=a2+1;i<nspher_set;i++)
		    if (fabs(r - sset_dist[ioff[i]+a1]) < ZERO) {
		      a3[num_a3] = i;
		      num_a3++;
		    }
		  switch (num_a3 + 1) {
		    case 1:
		    case 2: for(i=a2+1;i<nspher_set;i++)
			      if () {
			        a3[num_a3] = i;
				num_a3++;
			      }
			      */
	          break;
	}
	free_matrix(sset_geom,num_atoms);
      }
    }
    else if (num_atoms == 1) { /* Atom */
      fprintf(outfile,"    It is a spherical top.\n");
      rotor = atom;
    }
    else if (num_atoms <= 0) { /* ??? */
      punt("Fewer than 1 atom");
    }

    free(v1);
    free(v2);
    free(v3);
    free(IM);
    free_block(IT);
    free_block(ITAxes);
}

}} // namespace psi::input

namespace {
void if_to_invert_axis(double* v1, int* must_invert, int* should_invert, double* maxproj)
{
  using namespace psi::input;
  int xyz, nzero;
  double vabs;

  *maxproj = 0.0;
  *must_invert = 0;
  *should_invert = 0;
  
  nzero = 0;
  
  for(xyz=0; xyz<3; xyz++) {

    vabs = fabs(v1[xyz]);

    if (vabs < ZERO)
      nzero++;

    if (vabs > fabs(*maxproj)) {
      *maxproj = v1[xyz];
    }

  }

  if (nzero == 2) {
    if (*maxproj < 0.0)
      *must_invert = 1;
    else
      *must_invert = 0;
  }
  else if (nzero < 2) {
    if (*maxproj < 0.0)
      *should_invert = 1;
    else
      *should_invert = 0;
  }
}
} // namespace
