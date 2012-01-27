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

PsiReturnType fd_freq_1(Options &options, const boost::python::list& grad_list, int freq_irrep_only) {
  int pts = options.get_int("POINTS");
  double disp_size = options.get_double("DISP_SIZE");
  int print_lvl = options.get_int("PRINT");

  const boost::shared_ptr<Molecule> mol = psi::Process::environment.molecule();
  int Natom = mol->natom();
  boost::shared_ptr<MatrixFactory> fact;
  CdSalcList salc_list(mol, fact, 0xF, true, true);
  int Nirrep = salc_list.nirrep();

  // *** Build vectors that list indices of salcs for each irrep
  std::vector< std::vector<int> > salcs_pi;
  for (int h=0; h<Nirrep; ++h)
    salcs_pi.push_back( std::vector<int>() );
  for (int i=0; i<salc_list.ncd(); ++i)
    salcs_pi[salc_list[i].irrep()].push_back(i);

  // Now remove irreps that are not requested
  if (freq_irrep_only != -1) {
    for (int h=0; h<Nirrep; ++h)
      if (h != freq_irrep_only)
        salcs_pi[h].clear();
  }

  // Determine total num of salcs and where each irrep starts
  int Nsalc_all = salcs_pi[0].size();
  int salc_irr_start[8];
  salc_irr_start[0] = 0;
  for (int h=1; h<Nirrep; ++h) {
    Nsalc_all += salcs_pi[h].size();
    salc_irr_start[h] = salc_irr_start[h-1] + salcs_pi[h-1].size();
  }

  // ** Count displacements
  // symmetric displacements:
  std::vector<int> Ndisp_pi (Nirrep);
  if (pts == 3)
    Ndisp_pi[0] = 2 * salcs_pi[0].size();
  else if (pts == 5)
    Ndisp_pi[0] = 4 * salcs_pi[0].size();

  // asymmetric displacements:
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
  fprintf(outfile, "  symmetry-adapted, cartesian coordinates (fd_freq_1).\n\n");

  fprintf(outfile, "  %d gradients passed in, including the reference geometry.\n", (int) len(grad_list));

  // We are passing in the reference geometry at the moment, though we are not using
  // its gradient.  Could be removed later.

  if ((int) len(grad_list) != Ndisp_all+1) { // last gradient is the reference, non-displaced one
    fprintf(outfile,"gradients.size() is %d\n", (int) len(grad_list));
    fprintf(outfile,"Ndisp_all is %d\n", Ndisp_all);
    throw PsiException("FINDIF: Incorrect number of gradients passed in!",__FILE__,__LINE__);
  }

  // *** Generate complete list of gradients from unique ones.
  fprintf(outfile,"\tGenerating complete list of displacements from unique ones.\n\n");

  boost::shared_ptr<PointGroup> pg = mol->point_group();
  CharacterTable ct = mol->point_group()->char_table();
  int order = ct.order();

  // atom_map, how atoms are mapped to other atoms by operations
  int **atom_map = compute_atom_map(mol);
  if (print_lvl >= 3) {
    fprintf(outfile,"The atom map:\n");
    for (int i=0; i<Natom; ++i) {
      fprintf(outfile,"\t %d : ", i);
      for (int j=0; j<order; ++j)
        fprintf(outfile,"%4d", atom_map[i][j]);
      fprintf(outfile,"\n");
    }
    fprintf(outfile,"\n");
  }

  // Extract the symmetric gradients.
  std::vector< SharedMatrix > gradients;
  for (int i=0; i<Ndisp_pi[0]; ++i)
    gradients.push_back( (SharedMatrix) extract< SharedMatrix >(grad_list[i]) );

  if (print_lvl >= 3) {
    fprintf(outfile,"\tSymmetric gradients\n");
    for (int i=0; i<gradients.size(); ++i)
      gradients[i]->print();
  }

  // Extract the asymmetric gradients, one at a time and determine the gradient of the
  // non-computed displacements.

  //double tval;
  int disp_cnt = Ndisp_pi[0]; // step through original list of gradients for non-symmetric ones

//  SharedMatrix zero_grad(new Matrix(Natom, 3));

  for (int h=1; h<Nirrep; ++h) { // loop over asymmetric irreps

    if (Ndisp_pi[h] == 0) continue;

    IrreducibleRepresentation gamma = ct.gamma(h);

    if (print_lvl >= 3) {
      fprintf(outfile,"Characters for irrep %d\n", h);
      for (int i=0; i<order; ++i)
        fprintf(outfile," %5.1lf", gamma.character(i));
      fprintf(outfile,"\n");
    }

    // Find operation that takes + to - displacement.
    int op_disp;
    for (op_disp=0; op_disp<order; ++op_disp)
      if (gamma.character(op_disp) == -1)
        break;
    fprintf(outfile,"\tOperation %d takes plus displacements of irrep %s to minus ones.\n",
      op_disp+1, gamma.symbol());

    // Get 3x3 matrix representation of operation.
    SymmetryOperation so = ct.symm_operation(op_disp);

fprintf(outfile,"SO\n");
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

        for (int xyz=0; xyz<3; ++xyz) {
          double tval = 0.0;
          for (int xyz2=0; xyz2<3; ++xyz2)
            tval += so(xyz, xyz2) * gradients[disp_cnt]->get(atom,xyz2);
          new_grad->set(atom2, xyz, tval);
        }

        gradients.push_back(new_grad);
      }
      ++disp_cnt; //step through original gradient list
fprintf(outfile,"Transformed displaced gradient:\n");
gradients.back()->print();
    }

  } // end loop over irreps

  delete_atom_map(atom_map, mol);

  // Fix number of displacements for full list.
  for (int h=0; h<Nirrep; ++h) {
    if (pts == 3)
      Ndisp_pi[0] = 2 * salcs_pi[0].size();
    else if (pts == 5)
      Ndisp_pi[0] = 4 * salcs_pi[0].size();
  }
  Ndisp_all = 0;
  for (int h=0; h<Nirrep; ++h)
    Ndisp_all += Ndisp_pi[h];

  // Mass-weight all the gradients g_xm = 1/sqrt(m) g_x
/*  divide by masses in B matrix instead
  for (int i=0; i<Ndisp_all; ++i) {
    double **disp = gradients[i]->pointer();
    for (int a=0; a<Natom; ++a)
      for (int xyz=0; xyz<3; ++xyz)
        disp[a][xyz] /= sqrt(mol->mass(a));
  }
*/

  char **irrep_lbls = mol->irrep_labels();
  double **H_irr[8];

  std::vector<VIBRATION *> modes;

  for (int h=0; h<Nirrep; ++h) {

    // To store gradients in SALC displacement coordinates.
    double **grads_adapted = block_matrix(Ndisp_pi[h], 3*Natom);
    
    // Build B matrix / sqrt(masses).
    SharedMatrix B_irr_shared = salc_list.matrix_irrep(h);
    double **B_irr = B_irr_shared->pointer();

    for (int salc=0; salc<salcs_pi[h].size(); ++salc)
      for (int a=0; a<Natom; ++a)
        for (int xyz=0; xyz<3; ++xyz)
          B_irr[salc][3*a+xyz] /= sqrt(mol->mass(a));

    // Compute forces in internal coordinates, g_q = G_inv B u g_x
    // In this case, B = c * masses^(1/2).  =>  G=I.
    // Thus, g_q = c * g_x / sqrt(masses) or B g_x = g_q.

    for (int disp=0; disp<Ndisp_pi[h]; ++h)
      for (int salc=0; salc<salcs_pi[h].size(); ++salc)
        for (int a=0; a<Natom; ++a)
          for (int xyz=0; xyz<3; ++xyz)
            grads_adapted[disp][salc] += B_irr[salc][3*a+xyz] * gradients[disp]->get(a,xyz);

    //C_DGEMM('n', 't', Ndisp_pi[h], dim_q, 3*Natom, 1.0, g_x[0], 3*Natom,
    //  B_irr[0], 3*Natom, 0, g_q[0], dim_q);

    //** Construct force constant matrix from finite differences of forces
    H_irr[h] = init_matrix(salcs_pi[h].size(),salcs_pi[h].size());

    if (pts == 3) { // Hij = fj(i+1) - fj(i-1) / (2h)

      for (int i=0; i<salcs_pi[h].size(); ++i)
        for (int j=0; j<salcs_pi[h].size(); ++j)
          H_irr[h][i][j] = (grads_adapted[2*i+1][j] - grads_adapted[2*i][j]) / (2.0 * disp_size);

    }
    else if (pts == 5) { // fj(i-2) - fj(i+2) - 8fj(i-1) + 8fj(i+1)  / (12h)

      for (int i=0; i<salcs_pi[h].size(); ++i)
        for (int j=0; j<salcs_pi[h].size(); ++j)
          H_irr[h][i][j] = (  1.0 * grads_adapted[4*i][j]   - 1.0 * grads_adapted[4*i+1][j]
                            - 8.0 * grads_adapted[4*i+2][j] + 8.0 * grads_adapted[4*i+3][j] )
                            / (12.0 * disp_size);

    }

    if (print_lvl >= 3) {
      fprintf(outfile, "\n\tForce Constants for irrep %s in mass-weighted, ", irrep_lbls[h]);
      fprintf(outfile, "symmetry-adapted cartesian coordinates.\n");
      mat_print(H_irr[h], salcs_pi[h].size(), salcs_pi[h].size(), outfile);
    }

    // diagonalize force constant matrix
    int dim = salcs_pi[h].size();
    double *evals= init_array(dim);
    double **evects = block_matrix(dim, dim);

    sq_rsp(dim, dim, H_irr[h], evals, 3, evects, 1e-14);

    // Bu^1/2 * evects -> normal mode
    double **normal_irr = block_matrix(3*Natom, dim);
    C_DGEMM('t', 'n', 3*Natom, dim, dim, 1.0, B_irr[0], 3*Natom, evects[0],
      dim, 0, normal_irr[0], dim);

    if (print_lvl >= 2) {
      fprintf(outfile,"\n\tNormal coordinates (mass-weighted) for irrep %s:\n", irrep_lbls[h]);
      eivout(normal_irr, evals, 3*Natom, dim, outfile);
    }

    for (int i=0; i<salcs_pi[h].size(); ++i) {
      double *v = init_array(3*Natom);
      for (int x=0; x<3*Natom; ++x)
        v[x] = normal_irr[x][i];
      VIBRATION *vib = new VIBRATION(h, evals[i], v);
      modes.push_back(vib);
    }

    free(evals);
    free_block(evects);
    free_block(normal_irr);
  }

  // This print function also saves frequencies in wavefunction.
  print_vibrations(modes);

  for (int i=0; i<modes.size(); ++i)
    delete modes[i];
  modes.clear();

  // Build complete hessian for transformation to cartesians
  double **H = block_matrix(Nsalc_all, Nsalc_all);

  for (int h=0; h<Nirrep; ++h)
    for (int i=0; i<salcs_pi[h].size(); ++i) {
      int start = salc_irr_start[h];
      for (int j=0; j<=i; ++j)
        H[start + i][start + j] = H[start + j][start + i] = H_irr[h][i][j];
    }

  for (int h=0; h<Nirrep; ++h)
    if (salcs_pi[h].size()) free_block(H_irr[h]);

  // Transform Hessian into cartesian coordinates
  if (print_lvl >= 3) {
    fprintf(outfile, "\n\tFull force constant matrix in mass-weighted SALCS.\n");
    mat_print(H, Nsalc_all, Nsalc_all, outfile);
  }

  // Build Bu^-1/2 matrix for the whole Hessian
  SharedMatrix B_shared = salc_list.matrix();
  double **B = B_shared->pointer();

  // un mass-weighted below
  //for (int i=0; i<Nsalc_all; ++i) 
    //for (int a=0; a<Natom; ++a)
      //for (int xyz=0; xyz<3; ++xyz)
        //B[i][3*a+xyz] *= sqrt(mol->mass(a));

  double **Hx = block_matrix(3*Natom, 3*Natom);

  // Hx = Bt H B
  for (int i=0; i<Nsalc_all; ++i)
    for (int j=0; j<Nsalc_all; ++j)
      for (int x1=0; x1<3*Natom; ++x1)
        for (int x2=0; x2 <= x1; ++x2)
          Hx[x1][x2] += B[i][x1] * H[i][j] * B[j][x2];

  for (int x1=0; x1<3*Natom; ++x1)
    for (int x2=0; x2 < x1; ++x2)
      Hx[x2][x1] = Hx[x1][x2];

  free_block(H);

  if (print_lvl >= 3) {
    fprintf(outfile, "\n\tForce Constants in mass-weighted cartesian coordinates.\n");
    mat_print(Hx, 3*Natom, 3*Natom, outfile);
  }

  // Un-mass-weight Hessian
  for (int x1=0; x1<3*Natom; ++x1)
    for (int x2=0; x2<3*Natom; ++x2)
      Hx[x1][x2] *= sqrt(mol->mass(x1/3)) * sqrt(mol->mass(x2/3));

  if (print_lvl >= 3) {
    fprintf(outfile, "\n\tForce Constants in cartesian coordinates.\n");
    mat_print(Hx, 3*Natom, 3*Natom, outfile);
  }

  FILE *of_Hx = fopen("psi.file15.dat","w");
  fprintf(of_Hx,"%5d", Natom);
  fprintf(of_Hx,"%5d\n", 6*Natom);

  int cnt = -1;
  for (int i=0; i<3*Natom; ++i) {
    for (int j=0; j<3*Natom; ++j) {
      fprintf(of_Hx, "%20.10lf", Hx[i][j]);
      if (++cnt == 2) {
        fprintf(of_Hx,"\n");
        cnt = -1;
      }
    }
  }

  fclose(of_Hx);
  free_block(Hx);

  fprintf(outfile,"\n-------------------------------------------------------------\n");

  return Success;
}

}}

