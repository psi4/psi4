/*!
  \ingroup optking
  \file opt_data.h : header for structure that holds optimization data
*/

#ifndef _opt_opt_data_h_
#define _opt_opt_data_h_

#include <fstream>
#include <vector>

#include "package.h"

#include "linear_algebra.h"
#include "molecule.h"
#include "print.h"

#include "io.h"

namespace opt {

using namespace std;

// data for one optimization step
class STEP_DATA {

  double *f_q;        // Internal coordinate forces
  double *geom;          // cartesian coordinate values
  double energy;      // total energy
  double DE_predicted; // energy drop predicted for next step
  double *unit_step;  // unit vector in direction of step in the basis of internal coordinates
  double dq_norm;     // norm of step in internal coordinates
  double dq_gradient; // gradient along step 
  double dq_hessian;  // hessian along step
  double *dq;         // step in internal coordinates

  public:
    //STEP_DATA(ifstream & fin, int Nintco, int Ncart);  // read in date for one step

    STEP_DATA(int Nintco, int Ncart);   // allocate memory only

    ~STEP_DATA(); // free memory

    // save geometry and energy
    void save_geom_energy(double *geom_in, double energy_in, int Ncart);

    // save rest of stuff
    void save_step_info(double DE_predicted_in, double *unit_step_in, double dq_norm_in,
     double dq_gradient_in, double dq_hessian_in, int Nintco);

    // functions to read and write a step to the binary file
    void write(int istep, int Nintco, int Ncart);
    // read step from binary file
    void read(int istep, int Nintco, int Ncart);

    // functions to retrieve data
    double *g_forces_pointer(void) const { return f_q; }
    double *g_geom_const_pointer(void) const { return geom; }
    double *g_dq_pointer(void) const { return dq; }
    double g_energy(void) const { return energy; }
    double g_DE_predicted(void) const { return DE_predicted; }
    double g_dq_norm(void) const { return dq_norm; }
    double g_dq_gradient(void) const { return dq_gradient; }
    double g_dq_hessian(void) const { return dq_hessian; }
};

// data for an optimization
class OPT_DATA {
  int Nintco;        // num. of internal coordinates
  int Ncart;         // num. of cartesian coordinates
  double **H;        // Hessian matrix
  int iteration;     // num. of current iteration, 0, 1, 2, ...
                     // # of previous steps of data stored should be == iteration
  int consecutive_backsteps; // # of consecutive steps backwards, if any
  double *rfo_eigenvector;  // for RFO root-following
  std::vector<STEP_DATA *> steps; 

  public:

    // allocates memory for this step; reads in old steps from binary file
    OPT_DATA(int Nintco_in, int Ncart_in);

    // free memory
    ~OPT_DATA();

    // write data to binary file
    void write(void);

    // save geometry and energy to current (last) step
    void save_geom_energy(double *geom_in, double energy_in) {
      steps[steps.size()-1]->save_geom_energy(geom_in, energy_in, Ncart);
    }

    // save rest of stuff to current (last) step
    void save_step_info(double DE_predicted_in, double *unit_step_in, double dq_norm_in,
        double dq_gradient_in, double dq_hessian_in) {
      steps[steps.size()-1]->save_step_info(DE_predicted_in, unit_step_in, dq_norm_in,
        dq_gradient_in, dq_hessian_in, Nintco);
    }

    // return (pointers) to current-step data
    int g_iteration(void) const { return iteration; }
    double **g_H_pointer(void) { return H; }
    double g_energy(void) const { return steps[steps.size()-1]->g_energy(); }
    double *g_rfo_eigenvector_pointer(void) const { return rfo_eigenvector; }
    // return dimension of Hessian matrix
    int g_nintco(void) const { return Nintco; }

    void set_rfo_eigenvector(double *evect_in) {
      for (int i=0; i<Nintco+1; ++i)
        rfo_eigenvector[i] = evect_in[i];
    }

    // step data
    double *g_forces_pointer(void) const {
      return steps[steps.size()-1]->g_forces_pointer();
    }
    double *g_dq_pointer(void) const {
      return steps[steps.size()-1]->g_dq_pointer();
    }
    // return energy from the previous step (last entry - 1)
    double g_last_energy(void) const {
      if (steps.size() > 1)
        return steps[steps.size()-2]->g_energy();
      else return 0.0;
    }

    // return predicted energy change at the previous step (last entry - 1)
    double g_last_DE_predicted(void) const {
      if (steps.size() > 1)
        return steps[steps.size()-2]->g_DE_predicted();
      else return 0.0;
    }

    // return pointers to arbitrary-step data (pass in index starting at 0 ...)
    double g_energy(int i) const {
      return steps[i]->g_energy();
    }
    double *g_forces_pointer(int i) const {
      return steps.at(i)->g_forces_pointer();
    }
    double *g_geom_const_pointer(int i) const {
      return steps.at(i)->g_geom_const_pointer();
    }
    double *g_dq_pointer(int i) const {
      return steps.at(i)->g_dq_pointer();
    }
    double g_dq_norm(int i) const {
      return steps.at(i)->g_dq_norm();
    }
    double g_dq_gradient(int i) const {
      return steps.at(i)->g_dq_gradient();
    }
    double g_dq_hessian(int i) const {
      return steps.at(i)->g_dq_hessian();
    }

    bool previous_step_report(void) const;

    // check convergence of current step
    bool conv_check(opt::MOLECULE &) const;

    // summarize optimization up til now
    void summary(void) const;

    // perform Hessian update
    void H_update(opt::MOLECULE & mol);

    // read in cartesian Hessian
    double ** read_cartesian_H(void) const;

    // return number of steps present
    int nsteps(void) const { return steps.size(); }

    void decrement_iteration(void) { --iteration; }
    void increment_consecutive_backsteps(void) { ++consecutive_backsteps; }
    void reset_consecutive_backsteps(void) { consecutive_backsteps = 0; }
    int g_consecutive_backsteps(void) { return consecutive_backsteps; }

    void erase_last_step(void) { // free last step
      delete steps.back();
      steps.erase(steps.end()-1);
    }
    void erase_step(int i) {
      delete steps[i];
      steps.erase(steps.begin() + i);
    }
    void reset_iteration_to_size(void) {
      iteration = steps.size() + 1;
    }
};

}

#endif

