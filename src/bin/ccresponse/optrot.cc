/*! \file
  \ingroup ccresponse
  \brief Compute optical rotations using the CC linear response formalism.

  Optical rotation is determined from negative of the the mixed electric-
  dipole/magnetic-dipole polarizability, - << mu ; m >>, where mu is the
  electric-dipole vector operator and m is the magnetic-dipole vector
  operator.  We may choose between two representations of mu: length where mu
  = -r, and velocity, where mu = -p.  The velocity representation tensors are
  origin independent, but length representation tensors are not.  The
  "modified" velocity gauge involves subtraction of the zero-frequency <<
  mu(p); m>> tensor from the non-zer-frequency tensor.  Furthermore, if one
  chooses both representations (gauge = both), then the code also computes the
  origin-depenedence vector of the length-representation optical rotation,
  which requires the - << mu(r) ; m >> tensor.

  NB: The optical rotation expression requires the negative of all these
  tensors (much like dipole polarizabilities).

  The CC linear response tensors are computed as:

  << A ; B >>_w = 1/2 C(+/-w) P[A(-w), B(w)] *
      [<0|(1+L) [ABAR,X(B,w)]|0> + 1/2 <0|(1+L) [[HBAR,X(A,-w)],X(B,w)]|0>]

 = 1/2 [<0|(1+L) [ABAR,X(B,w)]|0> + 1/2 <0|(1+L) [[HBAR,X(A,-w)],X(B,w)]|0>]
 + 1/2 [<0|(1+L) [B*BAR,X(A*,w)]|0> + 1/2 <0|(1+L) [[HBAR,X(B*,-w)],X(A*,w)]|0>]

  If w=0, this becomes:

  << A ; B >>_0 = 1/2 C *
      [<0|(1+L) [ABAR,X(B,0)]|0> + 1/2 <0|(1+L) [[HBAR,X(A,0)],X(B,0)]|0>]

 = 1/2[<0|(1+L) [ABAR,X(B,0)]|0> + 1/2 <0|(1+L) [[HBAR,X(A,0)],X(B,0)]|0>]
 + 1/2[<0|(1+L) [B*BAR,X(A*,0)]|0> + 1/2 <0|(1+L) [[HBAR,X(B*,0)],X(A*,0)]|0>]

  Special notes: mu(p) and m are pure-imaginary operators, which means that a
  factor of i is included implicitly in the corresponding integrals. In the
  computation of the << mu(p); m >> tensor, where both operators carry i, an
  extra factor of -1 must appear in the computation of the tensors.  (This
  happens in both tersm of << A ; B >>_w above because complex conjugation of
  both opertors simultaneously doesn't change the overall sign of the tensor.
  Since we require the negative of the response tensor, the linresp() calls
  below for <<P;L>> have +1/2 prefactors rather than the usual -1/2.

  -TDC, 4/09
*/
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include <physconst.h>
#include <libmints/molecule.h>
#include <psi4-dec.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccresponse {

void pertbar(const char *pert, int irrep, int anti);
void compute_X(const char *pert, int irrep, double omega);
void linresp(double *tensor, double A, double B,
             const char *pert_x, int x_irrep, double omega_x,
             const char *pert_y, int y_irrep, double omega_y);

void optrot(void)
{
  double ***tensor_rl, ***tensor_pl, ***tensor_rp, **tensor0;
  double **tensor_rl0, **tensor_rl1, **tensor_pl0, **tensor_pl1;
  char **cartcomp, pert[32], pert_x[32], pert_y[32];
  int alpha, beta, i, j, k;
  double TrG_rl, TrG_pl, M, nu, bohr2a4, m2a, hbar, prefactor;
  double *rotation_rl, *rotation_pl, *rotation_rp, *rotation_mod, **delta;
  char lbl1[32], lbl2[32], lbl3[32];
  int compute_rl=0, compute_pl=0;

  /* Booleans for convenience */
  if(params.gauge == "LENGTH" || params.gauge == "BOTH") compute_rl=1;
  if(params.gauge == "VELOCITY" || params.gauge == "BOTH") compute_pl=1;

  cartcomp = (char **) malloc(3 * sizeof(char *));
  cartcomp[0] = strdup("X");
  cartcomp[1] = strdup("Y");
  cartcomp[2] = strdup("Z");

  tensor_rl = (double ***) malloc(params.nomega * sizeof(double **));
  tensor_pl = (double ***) malloc(params.nomega * sizeof(double **));
  tensor_rp = (double ***) malloc(params.nomega * sizeof(double **));
  for(i=0; i < params.nomega; i++) {
    tensor_rl[i] = block_matrix(3,3);
    tensor_pl[i] = block_matrix(3,3);
    tensor_rp[i] = block_matrix(3,3);
  }
  tensor0 = block_matrix(3,3);
  tensor_rl0 = block_matrix(3,3);
  tensor_rl1 = block_matrix(3,3);
  tensor_pl0 = block_matrix(3,3);
  tensor_pl1 = block_matrix(3,3);
  rotation_rl = init_array(params.nomega);
  rotation_pl = init_array(params.nomega);
  rotation_rp = init_array(params.nomega);
  rotation_mod = init_array(params.nomega);
  delta = block_matrix(params.nomega,3);

  if(compute_pl) {

    /* compute the zero-frequency Rosenfeld tensor for Koch's modified
       velocity optical rotation */
    /* NOTE!!!  The complex conjugation required for the response
       function is handled *implicitly* in the code below because the
       velocity-gauge version of the optical rotation tensor contains
       *two* pure-imaginary operators, p and m.  Thus, even though the
       corresponding perturbed wave functions change sign when we take
       complex conjugates of the p and m operators, the signs cancel,
       giving an identical contribution to the total *zero-frequency*
       Rosenfeld tensor as the original non-complex-conjugate term. */

    sprintf(lbl1, "<<P;L>>_(%5.3f)", 0.0);
    if(!params.restart || !psio_tocscan(CC_INFO, lbl1)) {

      for(alpha=0; alpha < 3; alpha++) {
        sprintf(pert, "P_%1s", cartcomp[alpha]);
        pertbar(pert, moinfo.mu_irreps[alpha], 1);
        compute_X(pert, moinfo.mu_irreps[alpha], 0);

        sprintf(pert, "L_%1s", cartcomp[alpha]);
        pertbar(pert, moinfo.l_irreps[alpha], 1);
        compute_X(pert, moinfo.l_irreps[alpha], 0);
      }

      fprintf(outfile, "\n\tComputing %s tensor.\n", lbl1); fflush(outfile);
      for(alpha=0; alpha < 3; alpha ++) {
        for(beta=0; beta < 3; beta++) {
          sprintf(pert_x, "P_%1s", cartcomp[alpha]);
          sprintf(pert_y, "L_%1s", cartcomp[beta]);
          linresp(&tensor0[alpha][beta], +1.0, 0.0,
                  pert_x, moinfo.mu_irreps[alpha], 0.0,
                  pert_y, moinfo.l_irreps[beta], 0.0);
        }
      }
      psio_write_entry(CC_INFO, lbl1, (char *) tensor0[0], 9*sizeof(double));

      /* Clean up disk space */
      psio_close(CC_LR, 0);
      psio_open(CC_LR, 0);

      for(j=CC_TMP; j <= CC_TMP11; j++) {
        psio_close(j,0);
        psio_open(j,0);
      }
    }
    else {
      fprintf(outfile, "Using %s tensor found on disk.\n", lbl1);
      psio_read_entry(CC_INFO, lbl1, (char *) tensor0[0], 9*sizeof(double));
    }

    if (params.wfn == "CC2")
      fprintf(outfile, "\n     CC2 Optical Rotation Tensor (Velocity Gauge): %s\n", lbl1);
    else if(params.wfn == "CCSD")
      fprintf(outfile, "\n    CCSD Optical Rotation Tensor (Velocity Gauge): %s\n", lbl1);

    fprintf(outfile, "  -------------------------------------------------------------------------\n");
    fprintf(outfile,   "   Evaluated at omega = 0.00 E_h (Inf nm, 0.0 eV, 0.0 cm-1)\n");
    fprintf(outfile, "  -------------------------------------------------------------------------\n");
    mat_print(tensor0, 3, 3, outfile);

  }

  for(i=0; i < params.nomega; i++) {
    zero_mat(tensor_rl0, 3, 3);
    zero_mat(tensor_rl1, 3, 3);
    zero_mat(tensor_pl0, 3, 3);
    zero_mat(tensor_pl1, 3, 3);

    sprintf(lbl1, "1/2 <<Mu;L>>_(%5.3f)", params.omega[i]);
    sprintf(lbl2, "1/2 <<P;L>>_(%5.3f)", params.omega[i]);
    if(!params.restart || ((compute_rl && !psio_tocscan(CC_INFO,lbl1)) || (compute_pl && !psio_tocscan(CC_INFO,lbl2)))) {

      /* prepare the dipole-length and/or dipole-velocity integrals */
      if(compute_rl) {
        for(alpha=0; alpha < 3; alpha++) {
          sprintf(pert, "Mu_%1s", cartcomp[alpha]);
          pertbar(pert, moinfo.mu_irreps[alpha], 0);
        }
      }
      if(compute_pl) {
        for(alpha=0; alpha < 3; alpha++) {
          sprintf(pert, "P_%1s", cartcomp[alpha]);
          pertbar(pert, moinfo.mu_irreps[alpha], 1);
        }
      }

      /* prepare the magnetic-dipole integrals */
      for(alpha=0; alpha < 3; alpha++) {
        sprintf(pert, "L_%1s", cartcomp[alpha]);
        pertbar(pert, moinfo.l_irreps[alpha], 1);
      }

      /* Compute the +omega magnetic-dipole and -omega electric-dipole CC wave functions */
      for(alpha=0; alpha < 3; alpha++) {
        if(compute_rl) {
          sprintf(pert, "Mu_%1s", cartcomp[alpha]);
          compute_X(pert, moinfo.mu_irreps[alpha], -params.omega[i]);
        }

        if(compute_pl) {
          sprintf(pert, "P_%1s", cartcomp[alpha]);
          compute_X(pert, moinfo.mu_irreps[alpha], -params.omega[i]);
        }

        sprintf(pert, "L_%1s", cartcomp[alpha]);
        compute_X(pert, moinfo.l_irreps[alpha], params.omega[i]);
      }

      fprintf(outfile, "\n");
      if(compute_rl) {
        fprintf(outfile, "\tComputing %s tensor.\n", lbl1); fflush(outfile);
        for(alpha=0; alpha < 3; alpha++) {
          for(beta=0; beta < 3; beta++) {
            sprintf(pert_x, "Mu_%1s", cartcomp[alpha]);
            sprintf(pert_y, "L_%1s", cartcomp[beta]);
            linresp(&tensor_rl0[alpha][beta], -0.5, 0.0,
                    pert_x, moinfo.mu_irreps[alpha], -params.omega[i],
                    pert_y, moinfo.l_irreps[beta], params.omega[i]);
          }
        }
        psio_write_entry(CC_INFO, lbl1, (char *) tensor_rl0[0], 9*sizeof(double));
      }
      if(compute_pl) {
        fprintf(outfile, "\tComputing %s tensor.\n", lbl2); fflush(outfile);
        for(alpha=0; alpha < 3; alpha++) {
          for(beta=0; beta < 3; beta++) {
            sprintf(pert_x, "P_%1s", cartcomp[alpha]);
            sprintf(pert_y, "L_%1s", cartcomp[beta]);
            linresp(&tensor_pl0[alpha][beta], +0.5, 0.0,
                    pert_x, moinfo.mu_irreps[alpha], -params.omega[i],
                    pert_y, moinfo.l_irreps[beta], params.omega[i]);
          }
        }
        psio_write_entry(CC_INFO, lbl2, (char *) tensor_pl0[0], 9*sizeof(double));
      }

      /* Clean up disk space */
      if(params.gauge != "BOTH") { /* don't clean up if we want both gauges */
        psio_close(CC_LR, 0);
        psio_open(CC_LR, 0);
      }

      for(j=CC_TMP; j <= CC_TMP11; j++) {
        psio_close(j,0);
        psio_open(j,0);
      }
    }
    else {
      fprintf(outfile, "\n");
      if(compute_rl) {
        fprintf(outfile, "\tUsing %s tensor found on disk.\n", lbl1); fflush(outfile);
        psio_read_entry(CC_INFO, lbl1, (char *) tensor_rl0[0], 9*sizeof(double));
      }
      if(compute_pl) {
        fprintf(outfile, "\tUsing %s tensor found on disk.\n", lbl2); fflush(outfile);
        psio_read_entry(CC_INFO, lbl2, (char *) tensor_pl0[0], 9*sizeof(double));
      }

    }

    sprintf(lbl1, "1/2 <<Mu;L*>>_(%5.3f)", params.omega[i]);
    sprintf(lbl2, "1/2 <<P*;L*>>_(%5.3f)", params.omega[i]);
    sprintf(lbl3, "<<P;Mu>>_(%5.3f)", params.omega[i]);
    if(!params.restart || ((compute_rl && !psio_tocscan(CC_INFO,lbl1)) || (compute_pl && !psio_tocscan(CC_INFO,lbl2)))) {

      /* prepare the dipole-length or dipole-velocity integrals */
      if(compute_rl) {
        for(alpha=0; alpha < 3; alpha++) {
          sprintf(pert, "Mu_%1s", cartcomp[alpha]);
          pertbar(pert, moinfo.mu_irreps[alpha], 0);
        }
      }
      if(compute_pl) {
        for(alpha=0; alpha < 3; alpha++) {
          sprintf(pert, "P*_%1s", cartcomp[alpha]);
          pertbar(pert, moinfo.mu_irreps[alpha], 1);
        }
      }

      /* prepare the complex-conjugate of the magnetic-dipole integrals */
      for(alpha=0; alpha < 3; alpha++) {
        sprintf(pert, "L*_%1s", cartcomp[alpha]);
        pertbar(pert, moinfo.l_irreps[alpha], 1);
      }

      /* Compute the -omega magnetic-dipole and +omega electric-dipole CC wave functions */
      for(alpha=0; alpha < 3; alpha++) {
        if(compute_rl) {
          sprintf(pert, "Mu_%1s", cartcomp[alpha]);
          compute_X(pert, moinfo.mu_irreps[alpha], params.omega[i]);
        }
        if(compute_pl) {
          sprintf(pert, "P*_%1s", cartcomp[alpha]);
          compute_X(pert, moinfo.mu_irreps[alpha], params.omega[i]);
        }

        sprintf(pert, "L*_%1s", cartcomp[alpha]);
        compute_X(pert, moinfo.l_irreps[alpha], -params.omega[i]);
      }

      fprintf(outfile, "\n");
      if(compute_rl) {
        fprintf(outfile, "\tComputing %s tensor.\n", lbl1); fflush(outfile);
        for(alpha=0; alpha < 3; alpha++) {
          for(beta=0; beta < 3; beta++) {
            sprintf(pert_x, "Mu_%1s", cartcomp[alpha]);
            sprintf(pert_y, "L*_%1s", cartcomp[beta]);
            linresp(&tensor_rl1[alpha][beta], -0.5, 0.0,
                    pert_x, moinfo.mu_irreps[alpha], params.omega[i],
                    pert_y, moinfo.l_irreps[beta], -params.omega[i]);
          }
        }
        psio_write_entry(CC_INFO, lbl1, (char *) tensor_rl1[0], 9*sizeof(double));
      }
      if(compute_pl) {
        fprintf(outfile, "\tComputing %s tensor.\n", lbl2); fflush(outfile);
        for(alpha=0; alpha < 3; alpha++) {
          for(beta=0; beta < 3; beta++) {
            sprintf(pert_x, "P*_%1s", cartcomp[alpha]);
            sprintf(pert_y, "L*_%1s", cartcomp[beta]);
            linresp(&tensor_pl1[alpha][beta], +0.5, 0.0,
                    pert_x, moinfo.mu_irreps[alpha], params.omega[i],
                    pert_y, moinfo.l_irreps[beta], -params.omega[i]);
          }
        }
        psio_write_entry(CC_INFO, lbl2, (char *) tensor_pl1[0], 9*sizeof(double));
      }

      if(params.gauge == "BOTH") {
        fprintf(outfile, "\tComputing %s tensor.\n", lbl3); fflush(outfile);
        for(alpha=0; alpha < 3; alpha++) {
          for(beta=0; beta < 3; beta++) {
            sprintf(pert_x, "P_%1s", cartcomp[alpha]);
            sprintf(pert_y, "Mu_%1s", cartcomp[beta]);
            linresp(&tensor_rp[i][alpha][beta],-0.5, 0.0,
                    pert_x, moinfo.mu_irreps[alpha], -params.omega[i],
                    pert_y, moinfo.mu_irreps[beta], params.omega[i]);
          }
        }
        for(alpha=0; alpha < 3; alpha++) {
          for(beta=0; beta < 3; beta++) {
            sprintf(pert_x, "P*_%1s", cartcomp[alpha]);
            sprintf(pert_y, "Mu_%1s", cartcomp[beta]);
            linresp(&tensor_rp[i][alpha][beta],-0.5, 1.0,
                    pert_x, moinfo.mu_irreps[alpha], params.omega[i],
                    pert_y, moinfo.mu_irreps[beta], -params.omega[i]);
          }
        }
        psio_write_entry(CC_INFO, lbl3, (char *) tensor_rp[i][0], 9*sizeof(double));
      }

      /* Clean up disk space */
      psio_close(CC_LR, 0);
      psio_open(CC_LR, 0);

      for(j=CC_TMP; j <= CC_TMP11; j++) {
        psio_close(j,0);
        psio_open(j,0);
      }

    }
    else {
      fprintf(outfile, "\n");
      if(compute_rl) {
        fprintf(outfile, "\tUsing %s tensor found on disk.\n", lbl1); fflush(outfile);
        psio_read_entry(CC_INFO, lbl1, (char *) tensor_rl1[0], 9*sizeof(double));
      }
      if(compute_pl) {
        fprintf(outfile, "\tUsing %s tensor found on disk.\n", lbl2); fflush(outfile);
        psio_read_entry(CC_INFO, lbl2, (char *) tensor_pl1[0], 9*sizeof(double));
      }
      if(params.gauge == "BOTH") {
        fprintf(outfile, "\tUsing %s tensor found on disk.\n", lbl3); fflush(outfile);
        psio_read_entry(CC_INFO, lbl3, (char *) tensor_rp[i][0], 9*sizeof(double));
      }
    }

    /* sum the two 1/2 tensors */
    for(j=0; j < 3; j++)
      for(k=0; k < 3; k++) {
        if(compute_rl) tensor_rl[i][j][k] = tensor_rl0[j][k] + tensor_rl1[j][k];
        if(compute_pl) tensor_pl[i][j][k] = tensor_pl0[j][k] + tensor_pl1[j][k];
      }

    /* compute the specific rotation */
    for(j=0,M=0.0; j < moinfo.natom ;j++) M += Process::environment.molecule()->mass(j);  /* amu */
    nu = params.omega[i]; /* hartree */
    bohr2a4 = _bohr2angstroms * _bohr2angstroms * _bohr2angstroms * _bohr2angstroms;
    m2a = _bohr2angstroms * 1.0e-10;
    hbar = _h/(2.0 * _pi);
    prefactor = 1.0e-2 * hbar/(_c * 2.0 * _pi * _me * m2a * m2a);
    prefactor *= prefactor;
    prefactor *= 288.0e-30 * _pi * _pi * _na * bohr2a4;

    if(compute_rl) {
      if (params.wfn == "CC2")
        fprintf(outfile, "\n            CC2 Optical Rotation Tensor (Length Gauge):\n");
      else if(params.wfn == "CCSD")
        fprintf(outfile, "\n           CCSD Optical Rotation Tensor (Length Gauge):\n");

      fprintf(outfile, "  -------------------------------------------------------------------------\n");
      fprintf(outfile,   "   Evaluated at omega = %8.6f E_h (%6.2f nm, %5.3f eV, %8.2f cm-1)\n", params.omega[i],
              (_c*_h*1e9)/(_hartree2J*params.omega[i]), _hartree2ev*params.omega[i],
              _hartree2wavenumbers*params.omega[i]);
      fprintf(outfile, "  -------------------------------------------------------------------------\n");
      mat_print(tensor_rl[i], 3, 3, outfile);

      TrG_rl = (tensor_rl[i][0][0] + tensor_rl[i][1][1] + tensor_rl[i][2][2])/(3.0 * params.omega[i]);

      rotation_rl[i] = prefactor * TrG_rl * nu * nu / M;
      fprintf(outfile, "\n   Specific rotation using length-gauge electric-dipole Rosenfeld tensor.\n");
      fprintf(outfile, "\t[alpha]_(%5.3f) = %10.5f deg/[dm (g/cm^3)]\n", params.omega[i], rotation_rl[i]);
    }

    if(compute_pl) {

      if (params.wfn == "CC2")
        fprintf(outfile, "\n          CC2 Optical Rotation Tensor (Velocity Gauge):\n");
      else if(params.wfn == "CCSD")
        fprintf(outfile, "\n         CCSD Optical Rotation Tensor (Velocity Gauge):\n");

      fprintf(outfile, "  -------------------------------------------------------------------------\n");
      fprintf(outfile,   "   Evaluated at omega = %8.6f E_h (%6.2f nm, %5.3f eV, %8.2f cm-1)\n", params.omega[i],
              (_c*_h*1e9)/(_hartree2J*params.omega[i]), _hartree2ev*params.omega[i],
              _hartree2wavenumbers*params.omega[i]);
      fprintf(outfile, "  -------------------------------------------------------------------------\n");
      mat_print(tensor_pl[i], 3, 3, outfile);

      TrG_pl = (tensor_pl[i][0][0] + tensor_pl[i][1][1] + tensor_pl[i][2][2])/(3.0 * params.omega[i]);
      TrG_pl /= params.omega[i];

      rotation_pl[i] = prefactor * TrG_pl * nu * nu / M;
      fprintf(outfile, "\n   Specific rotation using velocity-gauge electric-dipole Rosenfeld tensor.\n");
      fprintf(outfile, "\t[alpha]_(%5.3f) = %10.5f deg/[dm (g/cm^3)]\n", params.omega[i], rotation_pl[i]);

      /* subtract the zero-frequency beta tensor */
      for(j=0; j < 3; j++)
        for(k=0; k < 3; k++)
          tensor_pl[i][j][k] -= tensor0[j][k];

      if (params.wfn == "CC2")
        fprintf(outfile, "\n        CC2 Optical Rotation Tensor (Modified Velocity Gauge):\n");
      else if(params.wfn == "CCSD")
        fprintf(outfile, "\n        CCSD Optical Rotation Tensor (Modified Velocity Gauge):\n");

      fprintf(outfile, "  -------------------------------------------------------------------------\n");
      fprintf(outfile,   "   Evaluated at omega = %8.6f E_h (%6.2f nm, %5.3f eV, %8.2f cm-1)\n", params.omega[i],
              (_c*_h*1e9)/(_hartree2J*params.omega[i]), _hartree2ev*params.omega[i],
              _hartree2wavenumbers*params.omega[i]);
      fprintf(outfile, "  -------------------------------------------------------------------------\n");
      mat_print(tensor_pl[i], 3, 3, outfile);

      /* compute the specific rotation */
      TrG_pl = (tensor_pl[i][0][0] + tensor_pl[i][1][1] + tensor_pl[i][2][2])/(3.0 * params.omega[i]);
      TrG_pl /= params.omega[i];

      rotation_mod[i] = prefactor * TrG_pl * nu * nu / M;
      fprintf(outfile, "\n   Specific rotation using modified velocity-gauge Rosenfeld tensor.\n");
      fprintf(outfile, "\t[alpha]_(%5.3f) = %10.5f deg/[dm (g/cm^3)]\n", params.omega[i], rotation_mod[i]);
    }

    if(params.gauge == "BOTH") {
      delta[i][0] = prefactor * (tensor_rp[i][1][2] - tensor_rp[i][2][1]) * nu * nu / M;
      delta[i][1] = prefactor * (tensor_rp[i][2][0] - tensor_rp[i][0][2]) * nu * nu / M;
      delta[i][2] = prefactor * (tensor_rp[i][0][1] - tensor_rp[i][1][0]) * nu * nu / M;
      delta[i][0] /= 6.0 * params.omega[i];
      delta[i][1] /= 6.0 * params.omega[i];
      delta[i][2] /= 6.0 * params.omega[i];
      fprintf(outfile, "\n   Origin-dependence vector for length-gauge rotation deg/[dm (g/cm^3)]/bohr.\n");
      fprintf(outfile, "     Delta_x = %6.2f   Delta_y = %6.2f   Delta_z = %6.2f\n",
              delta[i][0], delta[i][1], delta[i][2]);
    }
  } /* loop i over nomega */

  if(params.nomega > 1) {  /* print a summary table for multi-wavelength calcs */

    if(compute_rl) {
      fprintf(outfile, "\n   ------------------------------------------\n");
      if (params.wfn == "CC2")
        fprintf(outfile,   "       CC2 Length-Gauge Optical Rotation\n");
      else
        fprintf(outfile,   "       CCSD Length-Gauge Optical Rotation\n");
      fprintf(outfile,   "   ------------------------------------------\n");

      if(params.gauge == "BOTH") {
        fprintf(outfile,   "       Omega           alpha                        Delta\n");
        fprintf(outfile,   "    E_h      nm   deg/[dm (g/cm^3)]        deg/[dm (g/cm^3)]/bohr\n");
        fprintf(outfile,   "   -----   ------ ------------------  ----------------------------------\n");
        fprintf(outfile,   "                                          x           y           z      \n");
        for(i=0; i < params.nomega; i++)
          fprintf(outfile, "   %5.3f   %6.2f      %10.5f    %10.5f  %10.5f  %10.5f\n", params.omega[i], (_c*_h*1e9)/(_hartree2J*params.omega[i]),
                  rotation_rl[i], delta[i][0], delta[i][1], delta[i][2]);
      }
      else {
        fprintf(outfile,   "       Omega           alpha\n");
        fprintf(outfile,   "    E_h      nm   deg/[dm (g/cm^3)]\n");
        fprintf(outfile,   "   -----   ------ ------------------\n");
        for(i=0; i < params.nomega; i++)
          fprintf(outfile, "   %5.3f   %6.2f      %10.5f\n", params.omega[i], (_c*_h*1e9)/(_hartree2J*params.omega[i]),
                  rotation_rl[i]);
      }
    }

    if(compute_pl) {
      fprintf(outfile, "\n   ------------------------------------------------------\n");
      if (params.wfn == "CC2")
        fprintf(outfile,   "            CC2 Velocity-Gauge Optical Rotation\n");
      else
        fprintf(outfile,   "            CCSD Velocity-Gauge Optical Rotation\n");
      fprintf(outfile,   "   ------------------------------------------------------\n");
      fprintf(outfile,   "       Omega           alpha (deg/[dm (g/cm^3)]\n");
      fprintf(outfile, "\n    E_h      nm   Velocity-Gauge  Modified Velocity-Gauge\n");
      fprintf(outfile,   "   -----   ------ --------------  -----------------------\n");
      for(i=0; i < params.nomega; i++)
        fprintf(outfile, "   %5.3f   %6.2f   %10.5f          %10.5f\n", params.omega[i], (_c*_h*1e9)/(_hartree2J*params.omega[i]), rotation_pl[i], rotation_mod[i]);
    }
  }

  free(rotation_rl);
  free(rotation_pl);
  free(rotation_rp);
  free(rotation_mod);
  free_block(delta);

  for(i=0; i < params.nomega; i++) {
    free_block(tensor_rl[i]);
    free_block(tensor_pl[i]);
    free_block(tensor_rp[i]);
  }
  free(tensor_rl);
  free(tensor_pl);
  free(tensor_rp);
  free_block(tensor0);
  free_block(tensor_rl0);
  free_block(tensor_rl1);
  free_block(tensor_pl0);
  free_block(tensor_pl1);

  free(cartcomp[0]);
  free(cartcomp[1]);
  free(cartcomp[2]);
  free(cartcomp);
}

}} // namespace psi::ccresponse
