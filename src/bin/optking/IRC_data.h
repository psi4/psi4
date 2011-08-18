/*!
  \ingroup optking
  \file IRC_data.h : header for structure that holds IRC data
*/

#ifndef _opt_irc_data_h_
#define _opt_irc_data_h_

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
class IRC_POINT {

  int coord_step;
  double *q_pivot;//pivot point for step
  double *q;           // internal coordinate values
  double *x;           // cartesian coordinate values
  double *f_q;         // internal coordinate forces
  double *f_x;         // cartesian coordinate forces
  double energy;       // total energy

  public:

    IRC_POINT(int coord_in, double *q_p_in, double *q_in, double *x_in, double *f_q_in, double *f_x_in, double E_in)
    {
      coord_step = coord_in;
      q_pivot = q_p_in;
      q = q_in;
      x = x_in;
      f_q = f_q_in;
      f_x = f_x_in;
      energy = E_in;
    }

    ~IRC_POINT() {free_array(q_pivot); free_array(q); free_array(x); free_array(f_q); free_array(f_x);} // free memory

    int g_coord_step(void) const { return coord_step; }
    double *g_q_pivot(void) const { return q_pivot; }
    double *g_q(void) const { return q; }
    double *g_x(void) const { return x; }
    double *g_f_q(void) const { return f_q; }
    double *g_f_x(void) const { return f_x; }
    double g_energy(void) const { return energy; }

};

// IRC reaction path data
class IRC_DATA {
  int coord_step;
  int sphere_step;

  std::vector<IRC_POINT *> steps; 

  public:

    IRC_DATA() {
      sphere_step = 0;
    };

    // free memory
    ~IRC_DATA() {
      for (int i=0; i<steps.size(); ++i)
        delete steps[i];
      steps.clear();
    }

    int size(void) { return steps.size(); }

    int g_next_coord_step(void) {
      if(size() == 0)
        return 0;
      else
        return steps.back()->g_coord_step() + 1;
    }

    double *g_q_pivot(void) const
    {
      return steps[steps.size()-1]->g_q_pivot();
    }
    double *g_q(void) const
    {
      return steps[steps.size()-1]->g_q();
    }
    double *g_x(void) const
    {
      return steps[steps.size()-1]->g_x();
    }
    double *g_f_q(void) const
    {
      return steps[steps.size()-1]->g_f_q();
    }
    double *g_f_x(void) const
    {
      return steps[steps.size()-1]->g_f_x();
    }

    const IRC_POINT &g_step(int index) const { return *(steps[index]); }

    //use like irc_step->g_step(2)->g_x();

    // check convergence of current step
//    bool conv_check(void) const;

    // summarize optimization up til now
//    void summary(void) const;

    void add_irc_point(int coord_in, double *q_p_in, double *q_in, double *x_in, double *f_q_in, double *f_x_in, double E_in)
    {
      IRC_POINT *onepoint = new IRC_POINT(coord_in, q_p_in, q_in, x_in, f_q_in, f_x_in, E_in);
      steps.push_back(onepoint);
    }

    friend void MOLECULE::irc_step(void);

    void point_converged(opt::MOLECULE &mol);
};

}

#endif

