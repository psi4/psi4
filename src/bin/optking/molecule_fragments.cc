/*! \file molecule_fragments.cc
    \ingroup optking
    \brief functions to handle fragments
*/

#include "molecule.h"

#include <cmath>
#include <iostream>
#include <sstream>

#include "linear_algebra.h"
#include "v3d.h"
#include "print.h"
#include "atom_data.h"
#include "physconst.h"

#define EXTERN
#include "globals.h"

#if defined(OPTKING_PACKAGE_PSI)
#include <libmints/molecule.h>
#endif

namespace opt {

using namespace v3d;

// This function manipulates a molecule containing a single fragment.
//
// If fragment_mode == SINGLE, then this function adds connectivity between the
// closest atom members of disconnected groups of atoms to form one superfragment.
// It does not itself add additional internal coordinates; these may be added by
// "add_simples_by_connectivity()" .
//
// If fragment_mode == MULTI, then this function splits one fragment into more
// than one, according to the connectivity previously established for the fragment.
void MOLECULE::fragmentize(void) {
  int i, j, xyz;

  if (fragments.size() != 1) return;

  int natom = fragments[0]->g_natom();
  const bool * const * const connectivity = fragments[0]->g_connectivity_pointer();

  // each row is (potentially) a fragment
  bool **frag_atoms = init_bool_matrix(natom,natom);
  int nfrag=0, ifrag=0;
  bool completed_fragment;
  int start_search = 0;    // fragment includes this atom
  bool more_fragments = true;

  for (ifrag=0; more_fragments==true; ++ifrag) {
    ++nfrag;
    frag_atoms[ifrag][start_search] = true; //connected to self

    do {
      completed_fragment = true;
      for (i=0; i<natom; ++i) {
        if (frag_atoms[ifrag][i]) {
          for (j=0; j<natom; ++j) {
            if (!frag_atoms[ifrag][j] && connectivity[i][j]) { // a new connection
              frag_atoms[ifrag][j] = true;
              completed_fragment = false;
            }
          }
        }
      }
    }
    while (!completed_fragment); // while still finding connections

    more_fragments = false;
    // are there any atoms not in a fragment yet
    for (j=0; j<natom; ++j) {
      bool atom_in_fragment = false;
      for (i=0; i<nfrag; ++i)
        if (frag_atoms[i][j]) atom_in_fragment = true;
      if (!atom_in_fragment) {
        start_search = j; // first atom of new fragment
        more_fragments = true;
        break;
      }
    }
  }

  for (i=0; i<nfrag; ++i) {
    fprintf(outfile, "\tDetected frag with atoms: ");
    for (j=0; j<natom; ++j)
      if (frag_atoms[i][j]) fprintf(outfile," %d", j+1);
    fprintf(outfile,"\n"); 
  }

  // Do nothing.  Atoms are all happily connected.
  if (nfrag == 1) {
    ;
  }
  else { // either connect fragment or split it up

    // make array of number of atoms in each fragment
    int *frag_natom = init_int_array(nfrag);
    for (i=0; i<nfrag; ++i)
      for (j=0; j<natom; ++j)
        if (frag_atoms[i][j]) ++frag_natom[i];

    // make lookup array of first atom of each fragment
    int *frag_offset = init_int_array(nfrag);
    for (i=1; i<nfrag; ++i)
      frag_offset[i] = frag_offset[i-1] + frag_natom[i-1];

    // add connectivity between the closest atoms of disconnected fragments
    // to make one complete superfragment
    if (Opt_params.fragment_mode == OPT_PARAMS::SINGLE) {
      GeomType geom = fragments[0]->g_geom_const_pointer();
      double tval, min;

      for (int f2=0; f2<nfrag; ++f2) {
        for (int f1=0; f1<f2; ++f1) {
          min = 1.0e12;

          for (int f1_atom=0; f1_atom<frag_natom[f1]; ++f1_atom) {
            for (int f2_atom=0; f2_atom<frag_natom[f2]; ++f2_atom) {
              tval = v3d_dist(geom[frag_offset[f1]+f1_atom], geom[frag_offset[f2]+f2_atom]);
              if (tval < min) {
                min = tval;
                i = frag_offset[f1]+f1_atom;
                j = frag_offset[f2]+f2_atom;
              }
            }
          }
          fragments[0]->connect(i,j); // connect closest atoms
          // Now check for possibly symmetry related atoms which are just as close
          // we need them all to avoid symmetry breaking.
          for (int f1_atom=0; f1_atom<frag_natom[f1]; ++f1_atom) {
            for (int f2_atom=0; f2_atom<frag_natom[f2]; ++f2_atom) {
              tval = v3d_dist(geom[frag_offset[f1]+f1_atom], geom[frag_offset[f2]+f2_atom]);
              if (fabs(tval - min) < 1.0e-14) {
                i = frag_offset[f1]+f1_atom;
                j = frag_offset[f2]+f2_atom;
                fragments[0]->connect(i,j);
              }
            }
          }
        }
      }
    }
    // create fragment objects for distinct atom groups
    else if (Opt_params.fragment_mode == OPT_PARAMS::MULTI) {
      double **geom = fragments[0]->g_geom();
      double **grad = fragments[0]->g_grad();
      double *Z     = fragments[0]->g_Z();
  
      for (ifrag=0; ifrag<nfrag; ++ifrag) {
    
        double *Z_frag = init_array(frag_natom[ifrag]);
        double **geom_frag = init_matrix(frag_natom[ifrag], 3);
        double **grad_frag = init_matrix(frag_natom[ifrag], 3);
  
        for (i=0; i<frag_natom[ifrag]; ++i) {
          Z_frag[i] = Z[frag_offset[ifrag] + i];
          for (xyz=0; xyz<3; ++xyz) {
            geom_frag[i][xyz] = geom[frag_offset[ifrag] + i][xyz];
            grad_frag[i][xyz] = grad[frag_offset[ifrag] + i][xyz];
          }
        }
  
        FRAG * one_frag = new FRAG(frag_natom[ifrag], Z_frag, geom_frag);
        one_frag->set_grad(grad_frag);
        free_matrix(grad_frag);

        // update connectivity and add to list
        one_frag->update_connectivity_by_distances();
        fragments.push_back(one_frag);
  
      }
  
      delete fragments[0];
      fragments.erase(fragments.begin()); // remove original, first 
    
      free_array(Z);
      free_matrix(geom);
      free_matrix(grad);
    }
    free_int_array(frag_natom);
    free_int_array(frag_offset);
  }
  free_bool_matrix(frag_atoms);
  return;
}

// add interfragment coordinates
// for now, coordinates for fragments in order 1-2-3-
void MOLECULE::add_interfragment(void) {
  int nA, nB;               // fragment natom
  const double * const * A; // fragment geometries
  const double * const * B;
  const bool * const * cA;  // fragment connectivities
  const bool * const * cB;
  int A1, A2, A3, B1, B2, B3;
  double tval, min;
  int ndA, ndB; // num of reference atoms on each fragment
  char error_msg[100];
  double **weight_A=NULL, **weight_B=NULL;
  FRAG *Afrag, *Bfrag;

  if (fragments.size() == 1) return;

  if (Opt_params.interfragment_mode == OPT_PARAMS::FIXED)
    fprintf(outfile,"\tInterfragment coordinate reference points to be selected from closest atoms and neighbors.\n");
  else if (Opt_params.interfragment_mode == OPT_PARAMS::PRINCIPAL_AXES)
    fprintf(outfile,"\tInterfragment coordinate reference points to be determined by principal axes.\n");

  for (int frag_i=0; frag_i<(fragments.size()-1); ++frag_i) {

    Afrag = fragments[frag_i];
    Bfrag = fragments[frag_i+1];

    A  = Afrag->g_geom_const_pointer();
    nA = Afrag->g_natom();
    cA = Afrag->g_connectivity_pointer();

    B  = Bfrag->g_geom_const_pointer();
    nB = Bfrag->g_natom();
    cB = Bfrag->g_connectivity_pointer();

    if (Opt_params.interfragment_mode == OPT_PARAMS::FIXED) {

      // A1 and B1 will be closest atoms between fragments
      min = 1e9;
      for (int iA=0; iA < nA; ++iA) {
        for (int iB=0; iB < nB; ++iB) {
          tval = v3d_dist(A[iA],B[iB]);
          if (tval < min) {
            min = tval; 
            A1 = iA;
            B1 = iB;
          }
        }
      }
      ndA = ndB = 1;

      fprintf(outfile,"\tNearest atoms on two fragments are %d and %d.\n",
        g_atom_offset(frag_i)+A1+1, g_atom_offset(frag_i+1)+B1+1);

      // A2 is bonded to A1, but A2-A1-B1 must not be collinear
      for (int iA=0; iA < nA; ++iA) {
        if (cA[iA][A1]) {
          if (v3d_angle(B[B1],A[A1],A[iA], tval)) {
            if (tval > Opt_params.interfrag_collinear_tol*_pi && tval < (1-Opt_params.interfrag_collinear_tol)*_pi) {
              A2 = iA;
              ++ndA;
              break;
            }
          }
        }
      }
      if (ndA == 1 && nA > 1) {
        fprintf(outfile, "Fragment A has >1 atoms but no non-collinear atom found bonded to %d", A1+1);
        sprintf(error_msg, "Fragment A has >1 atoms but no non-collinear atom found bonded to %d", A1+1);
        INTCO_EXCEPT(error_msg, true);
      }
  
      // B2 is bonded to B1, but A1-B1-B2 must not be collinear
      for (int iB=0; iB < nB; ++iB) {
        if (cB[iB][B1]) {
          if (v3d_angle(A[A1],B[B1],B[iB],tval)) {
            if (tval > Opt_params.interfrag_collinear_tol*_pi && tval < (1-Opt_params.interfrag_collinear_tol)*_pi) {
              B2 = iB;
              ++ndB;
              break;
            }
          }
        }
      }
      if (ndB == 1 && nB > 1) {
        fprintf(outfile, "Fragment B has >1 atoms but no non-collinear atom found bonded to %d", B1+1);
        sprintf(error_msg, "Fragment B has >1 atoms but no non-collinear atom found bonded to %d", B1+1);
        INTCO_EXCEPT(error_msg,true);
      }
  
      if (ndA == 2) { // we were able to locate a suitable A2
        // A3 is bonded to A2, but A3-A2-A1 must not be collinear
        for (int iA=0; iA < nA; ++iA) {
          if (iA != A1 && cA[iA][A2]) {
            if (v3d_angle(A[A1],A[A2],A[iA], tval)) {
              if (tval > Opt_params.interfrag_collinear_tol*_pi && tval < (1-Opt_params.interfrag_collinear_tol)*_pi) {
                A3 = iA;
                ++ndA;
                break;
              }
            }
          }
        }
        // if we couldn't find a 3rd atom bonded to A2, then look for a 3rd atom bonded to A1
        if (ndA != 3) {
          for (int iA=0; iA < nA; ++iA) {
            if (iA != A2 && cA[iA][A1]) {
              if (v3d_angle(A[A1],A[A2],A[iA], tval)) {
                if (tval > Opt_params.interfrag_collinear_tol*_pi && tval < (1-Opt_params.interfrag_collinear_tol)*_pi) {
                  A3 = iA;
                  ++ndA;
                  break;
                }
              }
            }
          }
        }
      }
  
      if (ndB == 2) { // we were able to locate a suitable B2
      // B3 is bonded to B2, but B3-B2-B1 must not be collinear
        for (int iB=0; iB < nB; ++iB) {
          if (iB != B1 && cB[iB][B2]) {
            if (v3d_angle(B[B1],B[B2],B[iB],tval)) {
              if (tval > Opt_params.interfrag_collinear_tol*_pi && tval < (1-Opt_params.interfrag_collinear_tol)*_pi) {
                B3 = iB;
                ++ndB;          
                break;        
              }
            }
          }
        }
        // if we couldn't find a 3rd atom bonded to B2, then look for a 3rd atom bonded to B1
        if (ndB != 3) { 
          for (int iB=0; iB < nB; ++iB) {
            if (iB != B2 && cB[iB][B1]) {
              if (v3d_angle(B[B1],B[B2],B[iB], tval)) {
                if (tval > Opt_params.interfrag_collinear_tol*_pi && tval < (1-Opt_params.interfrag_collinear_tol)*_pi) {
                  B3 = iB;
                  ++ndB;
                  break;
                }
              }
            }
          }
        }
      }
      // default weights are simply 1 to produce the reference points A1, A2, etc.
      weight_A = init_matrix(3, nA);
      weight_A[0][A1] = 1.0;
      weight_A[1][A2] = 1.0;
      weight_A[2][A3] = 1.0;
  
      weight_B = init_matrix(3, nB);
      weight_B[0][B1] = 1.0;
      weight_B[1][B2] = 1.0;
      weight_B[2][B3] = 1.0;
  
      if (Opt_params.print_lvl >= 3) {
        fprintf(outfile, "\tReference points are linear combination on fragment A\n");
        print_matrix(outfile, weight_A, 3, nA);
        fprintf(outfile, "\tReference points are linear combination on fragment B\n");
        print_matrix(outfile, weight_B, 3, nB);
      }

      INTERFRAG * one_IF = new INTERFRAG(Afrag, Bfrag, frag_i, frag_i+1, weight_A, weight_B, ndA, ndB);
      interfragments.push_back(one_IF);

    } // fixed interfragment coordinates
    else if (Opt_params.interfragment_mode == OPT_PARAMS::PRINCIPAL_AXES) {

      // ref point A[0] and B[0] will be the centers of mass
      // ref points A[1/2] and B[1/2] will on on principal axes
      // nothing to compute now
      if (nA == 1)
        ndA = 1;
      else if (nA == 2) // TODO check linearity
        ndA = 2;
      else 
        ndA = 3;

      if (nB == 1)
        ndB = 1;
      else if (nB == 2)
        ndB = 2;
      else
         ndA = 3;

      weight_A = weight_B = NULL;

      INTERFRAG * one_IF = new INTERFRAG(Afrag, Bfrag, frag_i, frag_i+1, NULL, NULL, ndA, ndB, true);
      interfragments.push_back(one_IF);

    }
  }

  fflush(outfile);
}

// Check to see if displacement along any of the interfragment modes breakes
// the symmetry of the molecule.  If so, freeze it.  This is a hack for now.
// will it work?  RAK 3-2012
void MOLECULE::freeze_interfragment_asymm(void) {
  double **coord_orig = g_geom_2D();
  double disp_size = 0.1;

  fprintf(outfile,"\tChecking interfragment coordinates for ones that break symmetry.\n");
  fflush(outfile);

  for (int I=0; I<interfragments.size(); ++I) {
    double **B = interfragments[I]->compute_B(); // ->g_nintco() X (3*atom A)+3(natom_B)

    int iA = interfragments[I]->g_A_index();
    int iB = interfragments[I]->g_B_index();
    int nA = interfragments[I]->g_natom_A();
    int nB = interfragments[I]->g_natom_B();

    for (int i=0; i<interfragments[I]->g_nintco(); ++i) {
      bool symmetric_intco = true;

      double **coord = matrix_return_copy(coord_orig, g_natom(), 3);

      for (int atom_a=0; atom_a<nA; ++atom_a)
        for (int xyz=0; xyz<3; ++xyz)
          coord[g_atom_offset(iA)+atom_a][xyz] += disp_size * B[i][3*atom_a+xyz];

      for (int atom_b=0; atom_b<nB; ++atom_b)
        for (int xyz=0; xyz<3; ++xyz)
          coord[g_atom_offset(iB)+atom_b][xyz] += disp_size * B[i][3*atom_b+xyz];


#if defined(OPTKING_PACKAGE_PSI)
      psi::Process::environment.molecule()->set_geometry(coord);
      symmetric_intco = psi::Process::environment.molecule()->valid_atom_map();
#elif defined(OPTKING_PACKAGE_QCHEM)
  // not implemented yet
#endif
      if (symmetric_intco)
        fprintf(outfile,"\tInterfragment coordinate %d, %d is symmetric.\n", I+1, i+1);
      else {
        fprintf(outfile,"\tInterfragment coordinate %d, %d breaks symmetry - freezing.\n", I+1, i+1);
        fflush(outfile);
        interfragments[I]->freeze(i);
      }
      free(coord);
    }
    free_matrix(B);
  }

#if defined(OPTKING_PACKAGE_PSI)
      psi::Process::environment.molecule()->set_geometry(coord_orig);
#endif

  return;
}

} // namespace opt

