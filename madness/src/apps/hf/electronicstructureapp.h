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
  
  $Id: electronicstructureapp.h 2173 2011-02-23 21:40:46Z justus.c79@gmail.com $
*/
/*
 * electronicstructureapp.h
 *
 *  Created on: Nov 5, 2008
 *      Author: eh7
 */

#ifndef ELECTRONICSTRUCTUREAPP_H_
#define ELECTRONICSTRUCTUREAPP_H_

#include <mra/mra.h>
#include <misc/ran.h>
#include "electronicstructureparams.h"
//#include "poperator.h"
#include "libxc.h"
#include "complexfun.h"
#include "esolver.h"

class LevelPmap : public WorldDCPmapInterface< Key<3> > {
private:
    const int nproc;
public:
    LevelPmap() : nproc(0) {};

    LevelPmap(World& world) : nproc(world.nproc()) {}

    /// Find the owner of a given key
    ProcessID owner(const Key<3>& key) const {
        Level n = key.level();
        if (n == 0) return 0;
        hashT hash;
        if (n <= 3 || (n&0x1)) hash = key.hash();
        else hash = key.parent().hash();
        //hashT hash = key.hash();
        return hash%nproc;
    }
};

class MolecularPotentialFunctor : public FunctionFunctorInterface<double,3> {
private:
    const MolecularEntity& _mentity;
public:
    MolecularPotentialFunctor(const MolecularEntity& mentity)
        : _mentity(mentity)
    {}

    double operator()(const coordT& x) const
    {
      return _mentity.nuclear_attraction_potential(x[0], x[1], x[2]);
    }
};

class MolecularNuclearChargeDensityFunctor : public FunctionFunctorInterface<double,3> {
private:
    const MolecularEntity& _mentity;
    const double R;
    const bool periodic;
    const std::vector<coordT> _specialpts;
public:
    MolecularNuclearChargeDensityFunctor(const MolecularEntity& mentity, const double& R,
        const bool& periodic, const std::vector<coordT>& specialpts)
      : _mentity(mentity), R(R), periodic(periodic), _specialpts(specialpts) {
    }

    virtual std::vector<coordT> special_points() const
    {
      return _specialpts;
    }

    virtual Level special_level()
    {
      return 10;
    }

    double operator()(const coordT& x) const
    {
        double big = 0.5*R + 6.0*_mentity.smallest_length_scale();
        // Only one contribution at any point due to the short
        // range of the nuclear charge density
        if (periodic)
        {
            for (int xr = -1; xr <= 1; xr += 1)
            {
                double xx = x[0] + xr*R;
                if (xx < big && xx > -big)
                {
                    for (int yr = -1; yr <= 1; yr += 1)
                    {
                        double yy = x[1] + yr*R;
                        if (yy < big && yy > -big)
                        {
                            for (int zr = -1; zr <= 1; zr += 1)
                            {
                                double zz = x[2] + zr*R;
                                if (zz < big && zz > -big)
                                {
                                    return _mentity.nuclear_charge_density(xx, yy, zz);
                                }
                            }
                        }
                    }
                }
            }
        }
        else
        {
            return _mentity.nuclear_charge_density(x[0], x[1], x[2]);
        }
        return 0.0;
    }
};

#define NTRANS 8

class AtomicBasisFunctor : public FunctionFunctorInterface<std::complex<double>,3> {
private:
  const AtomicBasisFunction aofunc;
  const double R;
  const double rangesq;
  const bool periodic;
  const KPoint kpt;
  coordT r;
  std::vector<coordT> _specialpts;
//  const int NTRANS = 2;
  Vector<std::complex<double>,2*NTRANS+1> tx;
  Vector<std::complex<double>,2*NTRANS+1> ty;
  Vector<std::complex<double>,2*NTRANS+1> tz;
public:
  AtomicBasisFunctor(const AtomicBasisFunction& aofunc, double R, 
                     bool periodic, const KPoint kpt)
      : aofunc(aofunc), R(R), rangesq(aofunc.rangesq()), periodic(periodic), kpt(kpt)
  {
    double x, y, z;
    aofunc.get_coords(x,y,z);
    
    r[0]=x; r[1]=y; r[2]=z;
    _specialpts=std::vector<coordT>(1,r);

    for (int ir = -NTRANS; ir <= NTRANS; ir += 1)
    {
      tx[ir+NTRANS] = exp(std::complex<double>(0.0, kpt.k[0]*ir * R));
      ty[ir+NTRANS] = exp(std::complex<double>(0.0, kpt.k[1]*ir * R));
      tz[ir+NTRANS] = exp(std::complex<double>(0.0, kpt.k[2]*ir * R));
    }
}

  virtual std::vector<coordT> special_points() const
  {
    return _specialpts;
  }

  std::complex<double> operator()(const coordT& x) const
  {
    std::complex<double> value = 0.0;
    if (periodic) {
        for (int xx=-NTRANS; xx<=NTRANS; xx++)  {
            const double xxR = xx*R + x[0] -r[0];
            const double xxRsq = xxR*xxR;
            if (xxRsq < rangesq) {
                for (int yy=-NTRANS; yy<=NTRANS; yy++) {
                    const double yyR = yy*R + x[1] - r[1];
                    const double yyRsq = xxRsq + yyR*yyR;
                    if (yyRsq < rangesq) {
                        for (int zz=-NTRANS; zz<=NTRANS; zz++)  {
                            double ao = aofunc(xx*R+x[0], yy*R+x[1], zz*R+x[2]);
                            if (fabs(ao) > 1e-8) {
                                std::complex<double> t1 = tx[xx+NTRANS]*ty[yy+NTRANS]*tz[zz+NTRANS];
            double kx0 = kpt.k[0] * x[0];
            double kx1 = kpt.k[1] * x[1];
            double kx2 = kpt.k[2] * x[2];
            std::complex<double> t2 = exp(std::complex<double>(0.0, -kx0 - kx1 - kx2));
                                value += t1 * t2 * ao;
          }
        }
      }
    }
            }
        }
    }
    else  {
      value = aofunc(x[0], x[1], x[2]);
    }
    return value;
  }
};

double rsquared(const coordT& r) {
    return r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
}


#endif /* ELECTRONICSTRUCTUREAPP_H_ */
