//
// rep.cc
//
// Additional modifications made by Justin Turney <jturney@ccqc.uga.edu>
// for use in PSI4.
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Edward Seidl <seidl@janed.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#include <cmath>

#include <libmints/pointgrp.h>

using namespace std;
using namespace psi;

/////////////////////////////////////////////////////////////////////////

SymRep::SymRep(int i) :
  n(i)
{
  zero();
}

SymRep::SymRep(const SymmetryOperation& so) :
  n(3)
{
  memset(d,0,sizeof(double)*25);
  for (int i=0; i < 3; i++)
    for (int j=0; j < 3; j++)
      d[i][j] = so[i][j];
}

SymRep::~SymRep()
{
  n=0;
}

SymRep::operator SymmetryOperation() const
{
  if (n != 3) {
    // ExEnv::err0() << indent << "SymRep::operator SymmetryOperation(): "
    //      << "trying to cast to symop when n == " << n << endl;
    // abort();
      throw PSIEXCEPTION("SymRep::operator SymmetryOperation(): trying to cast to symop when n != 3");
  }

  SymmetryOperation so;

  for (int i=0; i < 3; i++)
    for (int j=0; j < 3; j++)
      so[i][j] = d[i][j];

  return so;
}

SymRep
SymRep::operate(const SymRep& r) const
{
  if (r.n != n) {
    // ExEnv::err0() << indent << "SymRep::operate(): dimensions don't match: "
    //      << r.n << " != " << n << endl;
    // abort();
      throw PSIEXCEPTION("SymRep::operate(): dimensions don't match");
  }

  SymRep ret(n);

  for (int i=0; i < n; i++) {
    for (int j=0; j < n; j++) {
      double t=0;
      for (int k=0; k < n; k++)
        t += r[i][k] * d[k][j];
      ret[i][j] = t;
    }
  }

  return ret;
}

SymRep
SymRep::transform(const SymRep& r) const
{
  int i,j,k;

  if (r.n != n) {
    // ExEnv::err0() << indent
    //      << "SymRep::symm_transform(): dimensions don't match: "
    //      << r.n << " != " << n << endl;
    // abort();
    throw PSIEXCEPTION("SymRep::operate(): dimensions don't match");
  }

  SymRep ret(n), foo(n);

  // foo = r * d
  for (i=0; i < n; i++) {
    for (j=0; j < n; j++) {
      double t=0;
      for (k=0; k < n; k++)
        t += r[i][k] * d[k][j];
      foo[i][j] = t;
    }
  }

  // ret = (r*d)*r~ = foo*r~
  for (i=0; i < n; i++) {
    for (j=0; j < n; j++) {
      double t=0;
      for (k=0; k < n; k++)
        t += foo[i][k]*r[j][k];
      ret[i][j]=t;
    }
  }

  return ret;
}

void
SymRep::sigma_h()
{
  unit();

  if (n==3) {
    d[2][2] = -1.0;
  } else if (n==5) {
    d[3][3] = d[4][4] = -1.0;
  }
}

void
SymRep::sigma_xz()
{
  unit();

  if (n==2 || n==3 || n==4) {
    d[1][1] = -1.0;
    if (n==4)
      d[2][2] = -1.0;
  } else if (n==5) {
    d[2][2] = d[4][4] = -1.0;
  }
}

void
SymRep::sigma_yz()
{
  unit();

  if (n==2 || n==3 || n==4) {
    d[0][0] = -1.0;
    if (n==4)
      d[3][3] = -1.0;
  } else if (n==5) {
    d[2][2] = d[3][3] = -1.0;
  }
}

void
SymRep::rotation(int nt)
{
  double theta = (nt) ? 2.0*M_PI/nt : 2.0*M_PI;
  rotation(theta);
}

void
SymRep::rotation(double theta)
{
  zero();

  double ctheta = cos(theta);
  double stheta = sin(theta);
  double c2theta = cos(2*theta);
  double s2theta = sin(2*theta);

  switch (n) {
  case 1:
    d[0][0] = 1.0;
    break;

  case 3:
    d[0][0] = ctheta;
    d[0][1] = stheta;
    d[1][0] = -stheta;
    d[1][1] = ctheta;
    d[2][2] = 1.0;
    break;

  case 4:
  case 2:
    d[0][0] = ctheta;
    d[0][1] = stheta;
    d[1][0] = -stheta;
    d[1][1] = ctheta;

    // this is ok since d is hardwired
    d[2][2] = c2theta;
    d[2][3] = -s2theta;
    d[3][2] = s2theta;
    d[3][3] = c2theta;
    break;

  case 5:
    d[0][0] = 1.0;

    d[1][1] = c2theta;
    d[1][2] = s2theta;
    d[2][1] = -s2theta;
    d[2][2] = c2theta;

    d[3][3] = ctheta;
    d[3][4] = -stheta;
    d[4][3] = stheta;
    d[4][4] = ctheta;
    break;

  default:
    // ExEnv::err0() << indent << "SymRep::rotation(): n > 5 (" << n << ")\n";
    // abort();
    throw PSIEXCEPTION("SymRep::rotation(): n > 5");
  }

}

void
SymRep::c2_x()
{
  i();

  if (n==2 || n==3 || n==4) {
    d[0][0] = 1.0;
    if (n==4)
      d[3][3] = 1.0;
  } else if (n==5) {
    d[0][0] = d[1][1] = d[4][4] = 1.0;
  }
}

void
SymRep::c2_y()
{
  i();

  if (n==2 || n==3 || n==4) {
    d[1][1] = 1.0;
    if (n==4)
      d[2][2] = 1.0;
  } else if (n==5) {
    d[0][0] = d[1][1] = d[3][3] = 1.0;
  }
}

void
SymRep::c2_z()
{
  i();

  if (n==2 || n==3 || n==4) {
    d[1][1] = 1.0;
    if (n==4)
      d[2][2] = 1.0;
  } else if (n==5) {
    d[0][0] = d[1][1] = d[3][3] = 1.0;
  }
}

// void
// SymRep::print(ostream& os) const
// {
//   int i;
//
//   os << indent;
//   for (i=0; i < n; i++) os << scprintf("%11d",i+1);
//   os << endl;
//
//   for (i=0; i < n; i++) {
//     os << indent << scprintf("%3d ",i+1);
//     for (int j=0; j < n; j++)
//       os << scprintf(" %10.7f",d[i][j]);
//     os << endl;
//   }
//   os << endl;
// }

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
