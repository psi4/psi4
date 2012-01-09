/*! \file
    \ingroup OPTKING
    \brief fd_freq_1(): compute frequencies from gradients
*/

#include "findif.h"

#include <boost/python.hpp>
#include <boost/python/list.hpp>
using namespace boost;
using namespace boost::python;

#include <physconst.h>

namespace psi { namespace findif {

//int iE1(std::vector<int> & Ndisp_pi, std::vector< std::vector<int> > & salcs_pi, int pts,
// int irrep, int ii, int jj, int disp_i, int disp_j);

PsiReturnType fd_freq_1(Options &options, const boost::python::list& grad_list, int freq_irrep_only) {
  int pts = options.get_int("POINTS");
  double disp_size = options.get_double("DISP_SIZE");

  const boost::shared_ptr<Molecule> mol = psi::Process::environment.molecule();
  int Natom = mol->natom();
  boost::shared_ptr<MatrixFactory> fact;
  CdSalcList salc_list(mol, fact);
  int Nirrep = salc_list.nirrep();
  int Nsalc_all = salc_list.ncd();
  double tval;

  // *** Build vectors that list indices of salcs for each irrep
  std::vector< std::vector<int> > salcs_pi;
  for (int h=0; h<Nirrep; ++h)
    salcs_pi.push_back( std::vector<int>() );
  for (int i=0; i<Nsalc_all; ++i)
    salcs_pi[salc_list[i].irrep()].push_back(i);

  // Now remove irreps that are not requested
  if (freq_irrep_only != -1) {
    for (int h=0; h<Nirrep; ++h)
      if (h != freq_irrep_only)
        salcs_pi[h].clear();
  }

  // Count displacements.
  std::vector<int> Ndisp_pi (Nirrep);
  if (pts == 3)
    Ndisp_pi[0] = 2 * salcs_pi[0].size();
  else if (pts == 5)
    Ndisp_pi[0] = 4 * salcs_pi[0].size();

  for (int h=1; h<Nirrep; ++h) {
    if (pts == 3)
      Ndisp_pi[h] = salcs_pi[h].size();
    else if (pts == 5)
      Ndisp_pi[h] = 2* salcs_pi[h].size();
  }
  int Ndisp_all = 0;
  for (int h=0; h<Nirrep; ++h)
    Ndisp_all += Ndisp_pi[h];

  fprintf(outfile,"\n-------------------------------------------------------------\n\n");

  fprintf(outfile, "  Computing second-derivative from gradients using projected, \n");
  fprintf(outfile, "  symmetry-adapted, cartesian coordinates (fd_freq_1).\n");

  fprintf(outfile,"\tNumber of displacements per irrep:\n");
  for (int h=0; h<Nirrep; ++h)
    fprintf(outfile,"\t  Irrep %d: %d\n", h+1, Ndisp_pi[h]);

  fprintf(outfile, "\t%d gradients passed in, including the reference geometry.\n", (int) len(grad_list));
  // We are passing in the reference geometry at the moment; though we are not using
  // its gradient.  Could be removed later.
  if ((int) len(grad_list) != Ndisp_all+1) { // last gradient is the reference, non-displaced one
    fprintf(outfile,"gradients.size() is %d\n", (int) len(grad_list));
    fprintf(outfile,"Ndisp_all is %d\n", Ndisp_all);
    throw PsiException("FINDIF: Incorrect number of gradients passed in!",__FILE__,__LINE__);  }

  fprintf(outfile,"\tGenerating complete list of displacements from unique ones.\n");
  fflush(outfile);

  // *** Generate complete list of displacements from unique ones
  boost::shared_ptr<PointGroup> pg = mol->point_group();
  CharacterTable ct = mol->point_group()->char_table();
  int order = ct.order();

  int **atom_map = compute_atom_map(mol);  // how atoms are mapped to other atoms by operations
  fprintf(outfile,"atom map:\n");
  for (int i=0; i<Natom; ++i) {
    for (int j=0; j<order; ++j)
      fprintf(outfile,"%4d", atom_map[i][j]);
    fprintf(outfile,"\n");
  }

  // Extract the symmetric gradients.
  std::vector< SharedMatrix > gradients;
  for (int i=0; i<Ndisp_pi[0]; ++i)
    gradients.push_back( (SharedMatrix) extract< SharedMatrix >(grad_list[i]) );

  if (options.get_int("PRINT") > 2) {
    fprintf(outfile,"\tSymmetric gradients\n");
    for (int i=0; i<gradients.size(); ++i)
      gradients[i]->print();
  }

  int disp_cnt = Ndisp_pi[0]; // step through current list, for non-symmetric gradients

//  SharedMatrix zero_grad(new Matrix(Natom, 3));

  for (int h=1; h<Nirrep; ++h) { // loop over asymmetric irreps
    int op_disp; // operation that takes + displacement to - displacement
    if (Ndisp_pi[h] == 0) continue;

    IrreducibleRepresentation gamma = ct.gamma(h);

    fprintf(outfile,"Characters for irrep %d\n", h);
    for (int i=0; i<order; ++i)
      fprintf(outfile," %10.5lf", gamma.character(i));
    fprintf(outfile,"\n");

    for (op_disp=0; op_disp<order; ++op_disp) //Find operation that takes + to - displacement.
      if (gamma.character(op_disp) == -1)
        break;
    fprintf(outfile,"\tOperation %d takes plus displacements of irrep %s to minus ones.\n",
      op_disp+1, gamma.symbol());

    SymmetryOperation so = ct.symm_operation(op_disp); // get 3x3 matrix representation of operation

fprintf(outfile,"So\n");
for (int xyz=0; xyz<3; ++xyz)
  fprintf(outfile,"\t %10.5lf %10.5lf %10.5lf\n", so[xyz][0], so[xyz][1], so[xyz][2]);

    // Loop over coordinates of that irrep.
    for (int coord=0; coord<Ndisp_pi[h]; ++coord) {

      gradients.push_back( (SharedMatrix) extract< SharedMatrix >(grad_list[disp_cnt]) );

      fprintf(outfile,"Original displaced gradient:\n");
      gradients[disp_cnt]->print();

      for (int atom=0; atom<Natom; ++atom) {

        int atom2 = atom_map[atom][op_disp]; // how this atom transforms under this op.

        SharedMatrix new_grad(new Matrix(Natom, 3));

        tval = 0.0;
        for (int xyz=0; xyz<3; ++xyz) {
          for (int xyz2=0; xyz2<3; ++xyz2)
            tval += so(xyz, xyz2) * gradients[disp_cnt]->get(atom,xyz2);
          new_grad->set(atom2, xyz, tval);
        }

        gradients.push_back(new_grad);
      }
      if (pts == 5) { // There are 2 displacements.  Do the second one too.
        ;
      }
      ++disp_cnt; //step through original gradient list
      fprintf(outfile,"Transformed displaced gradient:\n");
      gradients.back()->print();
    }

    // Compute forces in internal coordinates.
    // g_q = - (BuBt)^-1 B u f_x.
    // Gradients should be mass-weighted.
    // In mass-weighted coordinates, B vectors are orthonormal.
    // Instead of mass-weighting gradients, I'll just divide the masses into B before multiplying.
    // g_q = c * f_x / sqrt(masses)

    // Build Bu^1/2 matrix for this irrep (or mass-weight gradients)
/*
    double **B_irr = block_matrix(dim_q, 3*Natom);
  
    for (int i=0; i<salcs_pi[h].size(); ++i) {
      int salc_i = salcs_pi[h][i];
      for (int c=0; c<salc_list[salc_i].ncomponent(); ++c) {
        int a          = salc_list[salc_i].component(c).atom;
        int xyz        = salc_list[salc_i].component(c).xyz;
        double coef    = salc_list[salc_i].component(c).coef;
        B_irr[i][3*a+xyz] = coef / sqrt(mol->mass(a));
      }
    }

    disp_start = 0;
    for (int i=0; i<h; ++i)
      disp_start += Ndisp_pi[i];

    // put gradients for this irrep in block_matrix
    double **g_x = block_matrix(Ndisp_pi[h], 3*Natom);

    for (int disp=disp_start; disp<Ndisp_pi[h]; ++disp) {
      for (int a=0; a<Natom; ++a)
        for (int xyz=0; xyz<3; ++xyz)
          g_x[disp][3*a+xyz] = gradients[disp]->get(h, a, xyz); // irrep necessary ?
    }

    // compute g_q ; gradients in internal coordinates B g_x^t = g_q^t -> g_x B^t = g_q
    double **g_q = block_matrix(Ndisp_pi[h], dim_q);

    C_DGEMM('n', 't', Ndisp_pi[h], dim_q, 3*Natom, 1.0, g_x[0], 3*Natom,
      B_irr[0], 3*Natom, 0, g_q[0], dim_q);

    free_block(g_x);
    free_block(B_irr);

    // Compute Hessian via finite differences
*/

//assume pts = 3 for the moment
/*
    if (options.get_int("PRINT") > 3) {
      fprintf(outfile, "\tForce Constants for irrep %d in mass-weighted, ", h+1);
      fprintf(outfile, "symmetry-adapted cartesian coordinates.\n");
      mat_print(H_irr, salcs_pi[h].size(), salcs_pi[h].size(), outfile);
    }

    // diagonalize force constant matrix
    int dim = salcs_pi[h].size(); 
    double *evals= init_array(dim);
    double **evects = block_matrix(dim, dim);

    sq_rsp(dim, dim, H_irr, evals, 3, evects, 1e-14);

    for (int i=0; i<salcs_pi[h].size(); ++i) {
      // copy eigenvector into row
      double *v = init_array(3*Natom);
      for (int x=0; x<3*Natom; ++x)
        v[x] = normal_irr[x][i];
      VIBRATION *vib = new VIBRATION(h, evals[i], v);
      modes.push_back(vib);
    }

    free(evals);
    free_block(evects);
    free_block(H_irr);
    free_block(normal_irr);

    free_block(g_q);

    //disp_start += Ndisp_pi[h];
*/
  }

  delete_atom_map(atom_map, mol);

  if (options.get_int("PRINT") > 2) {
    fprintf(outfile, "Non-symmetric gradients all.\n");
    if (options.get_int("PRINT") > 2)
      for (int i=Ndisp_pi[0]; i<gradients.size(); ++i)
        gradients[i]->print();
  }

/*
  print_vibrations(modes);

  for (int i=0; i<modes.size(); ++i)
    delete modes[i];
  modes.clear();
*/

  return Success;
}

/* iE1() returns index for the energy of a displacement, according to the order
generated in fd_geoms_freq_1()
ii and jj are coordinates, displaced by quantized steps disp_i and disp_j
disp_i,disp_j are {-1,0,+1} for a three-point formula
disp_j,disp_j are {-2,-1,0,+1,+2} for a five-point formula
It is assumed that ii >= jj .
For diagonal displacements disp_j=0 and jj is arbitrary/meaningless.
*/

int iE1(std::vector<int> & Ndisp_pi, std::vector< std::vector<int> > & salcs_pi, int pts,
  int irrep, int ii, int jj, int disp_i, int disp_j) {

  int ndiag_this_irrep;
  int rval=-1;

  // compute starting location of displacements for this irrep
  int start_irr = 0;
  for (int h=0; h<irrep; ++h)
    start_irr += Ndisp_pi[h];

  if (pts == 3) {
    if (disp_j == 0) {  // diagonal; all diagonals at beginning of irrep 
      if (irrep == 0) {
        if (disp_i == -1)
          rval = 2*ii;   // f(-1, 0)
        else if (disp_i == +1)
          rval = 2*ii+1; // f(+1, 0)
      }
      else if (disp_i == -1 || disp_i == +1)
        rval = start_irr + ii;     // f(+1,0) = f(-1, 0)
    }
    else {    // off_diagonal
      if (irrep == 0)
        ndiag_this_irrep = 2 * salcs_pi[0].size();
      else
        ndiag_this_irrep = salcs_pi[irrep].size();

      int ij_pair = 2 * ((ii*(ii-1))/2 + jj);

      if      (disp_i == +1 && disp_j == +1) rval = start_irr + ndiag_this_irrep + ij_pair;
      else if (disp_i == -1 && disp_j == -1) rval = start_irr + ndiag_this_irrep + ij_pair + 1;
    }
  }
  else if (pts == 5) {
    if (disp_j == 0) {   // diagonal 
      if (irrep == 0) {
        if (disp_i == -2)      rval = start_irr + 4*ii;     // f(-2, 0)
        else if (disp_i == -1) rval = start_irr + 4*ii+1;   // f(-1, 0)
        else if (disp_i ==  1) rval = start_irr + 4*ii+2;   // f(+1, 0)
        else if (disp_i ==  2) rval = start_irr + 4*ii+3;   // f(+2, 0)
      }
      else { // irrep != 0
        if      (disp_i == -2 || disp_i == 2) rval = start_irr + 2*ii;     // f(-2, 0)
        else if (disp_i == -1 || disp_i == 1) rval = start_irr + 2*ii+1;   // f(-1, 0)
      }
    }
    else {   //off-diagonal
      if (irrep == 0)
        ndiag_this_irrep = 4 * salcs_pi[0].size() ;
      else
        ndiag_this_irrep = 2 * salcs_pi[irrep].size() ;

      int ij_pair = 8 * ((ii*(ii-1))/2 + jj);

      if      (disp_i == -1 && disp_j == -2) rval = start_irr + ndiag_this_irrep + ij_pair;
      else if (disp_i == -2 && disp_j == -1) rval = start_irr + ndiag_this_irrep + ij_pair+1;
      else if (disp_i == -1 && disp_j == -1) rval = start_irr + ndiag_this_irrep + ij_pair+2;
      else if (disp_i == +1 && disp_j == -1) rval = start_irr + ndiag_this_irrep + ij_pair+3;
      else if (disp_i == -1 && disp_j == +1) rval = start_irr + ndiag_this_irrep + ij_pair+4;
      else if (disp_i == +1 && disp_j == +1) rval = start_irr + ndiag_this_irrep + ij_pair+5;
      else if (disp_i == +2 && disp_j == +1) rval = start_irr + ndiag_this_irrep + ij_pair+6;
      else if (disp_i == +1 && disp_j == +2) rval = start_irr + ndiag_this_irrep + ij_pair+7;
    }
  }

  if (rval < 0) {
    fprintf(outfile,"Problem finding displaced energy.\n");
    throw PsiException("FINDIF: Problem finding displaced energy.",__FILE__,__LINE__);
  }
  //fprintf(outfile,"irrep: %d, ii: %d, jj: %d, disp_i: %d, disp_j: %d, rval: %d \n",
  //  irrep, ii, jj, disp_i, disp_j, rval);
  return rval;
}

}}

