/*
  This file is part of MADNESS.
  
  Copyright (C) 2007,2010 Oak Ridge National Laboratory
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
  
  For more information please contact:
  
  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367
  
  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680
  
  $Id: electronicstructureparams.h 2173 2011-02-23 21:40:46Z justus.c79@gmail.com $
*/
#ifndef ELECTRONICSTRUCTUREPARAMS_H_
#define ELECTRONICSTRUCTUREPARAMS_H_

#include <mra/mra.h>
#include <mra/vmra.h>
#include <misc/ran.h>
#include <misc/misc.h>
#include "mentity.h"
#include "molecularbasis.h"

namespace madness {

struct ElectronicStructureParams
{
  // Size of the cubic box (this needs to change)
  double L;
  // Number of electrons
  int nelec;
  // 1 - LDA; 2 - Hartree-Fock
  int functional;
  // Low value in the BSH / Coulomb fit
  double lo;
  // Spin-polarized
  bool spinpol;
  // Periodic system
  bool periodic;
  // Maximum number of interations
  int maxits;
  // Is source function a nuclear potential or a nuclear charge density?
  bool ispotential;
  // Thresh
  double thresh;
  // Order of wavelets
  int waveorder;
  // Max thresh
  double maxthresh;
  // Max order of wavelets
  int maxwaveorder;
  // Number of empty states
  int nempty;
  // Smearing parameter
  double smear;
  // Total number of bands
  int nbands;
  // Size of k-mesh (hardcoded for 3-d)
  int ngridk0, ngridk1, ngridk2;
  // Maximum occupation
  double maxocc;
  // Read k-points?
  bool kpoints;
  // Fractional coordinates?
  bool fractional;
  // Maximum size of subspace
  unsigned int maxsub;
  // maxrotn
  double maxrotn;
  // Solve for canonical orbitals?
  bool canon;
  // Don't use solver = 0; full KAIN = 1; k-point KAIN = 2 
  int solver;
  // k-mesh offset
  double koffset0, koffset1, koffset2;
  // initial basis set
  std::string basis;
  // number of IO nodes
  int nio;
  // restart calculation; no restart = 0; restart fully = 1;
  // restart using only density = 2;
  int restart;
  // total amount of electronic charge
  double ncharge;
  // width for smearing
  double swidth;
  // print matrices
  bool print_matrices;
  // plot orbitals
  bool plotorbs;
  // convergence criterion for residual
  double rcriterion;
  
  template <typename Archive>
  void serialize(Archive& ar) {
      ar & L & nelec & functional & lo & spinpol &
        periodic & maxits & ispotential & thresh &
        waveorder & maxthresh & maxwaveorder & nempty &
        smear & nbands & ngridk0 & ngridk1 & ngridk2 &
        maxocc & kpoints & fractional & maxsub & 
        maxrotn & canon & solver & koffset0 & koffset1 & 
        koffset2 & basis & nio & restart & ncharge & 
        swidth & print_matrices & plotorbs & rcriterion;
  }

  ElectronicStructureParams()
  {
    L = 10.0;
    nelec = 1;
    functional = 1;
    lo = 1e-4;
    smear = 0.001;
    spinpol = false;
    periodic = false;
    ispotential = false;
    maxits = 100;
    thresh = 1e-6;
    waveorder = 8;
    maxthresh = 1e-6;
    maxwaveorder = 8;
    nempty = 2;
    ngridk0 = 1; ngridk1 = 1; ngridk2 = 1;
    maxocc = 2.0;
    nbands = nelec/maxocc + nempty;
    kpoints = false;
    fractional = false;
    maxsub = 1;
    maxrotn = 0.5;
    canon = true;
    solver = 1;
    koffset0 = 0.0;
    koffset1 = 0.0;
    koffset2 = 0.0;
    basis = "sto-3g";
    nio = 1;
    restart = 0;
    ncharge = 0;
    swidth = 0.001;
    print_matrices = true;
    plotorbs = false;
    rcriterion = 1e-4;
  }

  void read_file(const std::string& filename)
  {
    std::ifstream f(filename.c_str());
    position_stream(f, "dft");
    std::string s;
    bool bnelec = false;
    while (f >> s)
    {
      if (s == "end")
      {
        break;
      }
      else if (s == "nelec")
      {
        f >> nelec;
        bnelec = true;
      }
      else if (s == "solver")
      {
        f >> solver;
      }
      else if (s == "L")
      {
        f >> L;
      }
      else if (s == "functional")
      {
        f >> functional;
      }
      else if (s == "basis")
      {
        f >> basis;
      }
      else if (s == "lo")
      {
        f >> lo;
      }
      else if (s == "swidth")
      {
        f >> swidth;
      }
      else if (s == "nio")
      {
        f >> nio;
      }
      else if (s == "restart")
      {
        f >> restart;
      }
      else if (s == "spinpol")
      {
        std::string tempstr;
        f >> tempstr;
        if (tempstr == "true")
        {
          spinpol = true;
        }
        else if (tempstr == "false")
        {
          spinpol = false;
        }
        else
        {
          MADNESS_EXCEPTION("input error -- spinpol", 0);
        }
      }
      else if (s == "canon")
      {
        std::string tempstr;
        f >> tempstr;
        if (tempstr == "true")
        {
          canon = true;
        }
        else if (tempstr == "false")
        {
          canon = false;
        }
        else
        {
          MADNESS_EXCEPTION("input error -- canon", 0);
        }
      }
      else if (s == "periodic")
      {
        std::string tempstr;
        f >> tempstr;
        if (tempstr == "true")
        {
          periodic = true;
        }
        else if (tempstr == "false")
        {
          periodic = false;
        }
        else
        {
          MADNESS_EXCEPTION("input error -- periodic", 0);
        }
      }
      else if (s == "ispotential")
      {
        std::string tempstr;
        f >> tempstr;
        if (tempstr == "true")
        {
          ispotential = true;
        }
        else if (tempstr == "false")
        {
          ispotential = false;
        }
        else
        {
          MADNESS_EXCEPTION("input error -- ispotential", 0);
        }
      }
      else if (s == "maxits")
      {
        f >> maxits;
      }
      else if (s == "maxsub")
      {
        f >> maxsub;
      }
      else if (s == "maxrotn")
      {
        f >> maxrotn;
      }
      else if (s == "thresh")
      {
        f >> maxthresh;
        thresh = maxthresh;
      }
      else if (s == "waveorder")
      {
        f >> maxwaveorder;
        waveorder = maxwaveorder;
      }
      else if (s == "nempty")
      {
        f >> nempty;
      }
      else if (s == "kpoints")
      {
        std::string tempstr;
        f >> tempstr;
        if (tempstr == "true")
        {
          kpoints = true;
        }
        else if (tempstr == "false")
        {
          kpoints = false;
        }
        else
        {
          MADNESS_EXCEPTION("input error -- kpoints", 0);
        }
      }
      else if (s == "fractional")
      {
        std::string tempstr;
        f >> tempstr;
        if (tempstr == "true")
        {
          fractional = true;
        }
        else if (tempstr == "false")
        {
          fractional = false;
        }
        else
        {
          MADNESS_EXCEPTION("input error -- fractional", 0);
        }
      }
      else if (s == "ngridk")
      {
        f >> ngridk0; f >> ngridk1; f >> ngridk2;
      }
      else if (s == "koffset")
      {
        f >> koffset0; f >> koffset1; f >> koffset2;
      }
      else if (s == "print_matrices") {
          print_matrices = true;
      }
      else if (s == "noprint_matrices") {
          print_matrices = false;
      }
      else if (s == "plotorbs")
      {
        std::string tempstr;
        f >> tempstr;
        if (tempstr == "true")
        {
          plotorbs = true;
        }
        else if (tempstr == "false")
        {
          plotorbs = false;
        }
        else
        {
          MADNESS_EXCEPTION("input error -- plotorbs", 0);
        }
      }
      else if (s == "rcriterion")
      {
        f >> rcriterion;
      }
      else
      {
        std::cout << "esolver: unrecognized input keyword " << s << std::endl;
        MADNESS_EXCEPTION("input error", 0);
      }
    }
    // No spin polarization
    //if (spinpol = true) MADNESS_EXCEPTION("spinpol not implemented", 0);
    // nelec is required
    //if (!bnelec) MADNESS_EXCEPTION("nelec required", 0);
//    // maximum occupation
//    maxocc = (spinpol) ? 1.0 : 2.0;
//    // compute total number of bands
//    nbands = nelec/maxocc + nempty;
    // kpoints only for periodic
    if (kpoints && !periodic)
      MADNESS_EXCEPTION("input error -- k-points only valid with periodic calculation", 0);
  }

  void set_molecular_info(const MolecularEntity& mentity, const AtomicBasisSet& aobasis) {
      lo = mentity.smallest_length_scale();
  }
};

}
#endif /* ELECTRONICSTRUCTUREPARAMS_H_ */
