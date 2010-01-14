//
// tetra.cc
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
#include <cstring>

#include <libmints/pointgrp.h>

using namespace psi;

// these are the operations which make up T
static void
t_ops(SymmetryOperation *symop)
{
  // identity
  symop[0].E();

  // C2(x)
  symop[9].c2_x();

  // C2(y)
  symop[10].c2_y();

  // C2(z)
  symop[11].rotation((double)M_PI);

  // a = ( 1, 1, 1)
  // b = (-1,-1, 1)
  // c = ( 1,-1,-1)
  // d = (-1, 1,-1)
  // C3 (a)
  symop[1][0][2] =  1.0;
  symop[1][1][0] =  1.0;
  symop[1][2][1] =  1.0;

  // C3 (b)
  symop[2] = symop[1].transform(symop[11]);

  // C3 (c)
  symop[3] = symop[1].transform(symop[9]);

  // C3 (d)
  symop[4] = symop[1].transform(symop[10]);

  // C3^2 (a)
  symop[5][0][1] =  1.0;
  symop[5][1][2] =  1.0;
  symop[5][2][0] =  1.0;

  // C3^2 (b)
  symop[6] = symop[5].transform(symop[11]);

  // C3^2 (c)
  symop[7] = symop[5].transform(symop[9]);

  // C3^2 (d)
  symop[8] = symop[5].transform(symop[10]);
}

// this gives us the operations in Td which aren't in T.
static void
td_ops(SymmetryOperation *symop)
{
  // S4 (x)
  symop[0][0][0] = -1.0;
  symop[0][1][2] = -1.0;
  symop[0][2][1] =  1.0;

  // S4^3 (x)
  symop[1][0][0] = -1.0;
  symop[1][1][2] =  1.0;
  symop[1][2][1] = -1.0;

  // S4 (y)
  symop[2][0][2] =  1.0;
  symop[2][1][1] = -1.0;
  symop[2][2][0] = -1.0;

  // S4^3 (y)
  symop[3][0][2] = -1.0;
  symop[3][1][1] = -1.0;
  symop[3][2][0] =  1.0;

  // S4 (z)
  symop[4][0][1] = -1.0;
  symop[4][1][0] =  1.0;
  symop[4][2][2] = -1.0;

  // S4^3 (z)
  symop[5][0][1] =  1.0;
  symop[5][1][0] = -1.0;
  symop[5][2][2] = -1.0;

  // a = ( 1, 1, 1)
  // b = (-1,-1, 1)
  // c = ( 1,-1,-1)
  // d = (-1, 1,-1)
  // sigma (ac)
  symop[6][0][0] =  1.0;
  symop[6][1][2] =  1.0;
  symop[6][2][1] =  1.0;

  // sigma (bd)
  symop[7][0][0] =  1.0;
  symop[7][1][2] = -1.0;
  symop[7][2][1] = -1.0;

  // sigma (ad)
  symop[8][0][2] =  1.0;
  symop[8][1][1] =  1.0;
  symop[8][2][0] =  1.0;

  // sigma (bc)
  symop[9][0][2] = -1.0;
  symop[9][1][1] =  1.0;
  symop[9][2][0] = -1.0;

  // sigma (ab)
  symop[10][0][1] =  1.0;
  symop[10][1][0] =  1.0;
  symop[10][2][2] =  1.0;

  // sigma (dc)
  symop[11][0][1] = -1.0;
  symop[11][1][0] = -1.0;
  symop[11][2][2] =  1.0;
}

////////////////////////////////////////////////////////////////////////////

void
CharacterTable::t()
{
  // t_ops gives us all the symmetry operations we need
  t_ops(symop);

  int i;

  gamma_[0].init(g,1,"A");
  for (i=0; i < g; i++)
    gamma_[0].rep[i][0][0] = 1.0;

  IrreducibleRepresentation& ire = gamma_[1];
  ire.init(g,2,"E");
  ire.complex_=1;

  IrreducibleRepresentation& irt = gamma_[2];
  irt.init(g,3,"T");
  irt.nrot_ = 1;
  irt.ntrans_ = 1;

  // the symmetry operation matrices give us a basis for irrep T
  for (i=0; i < g; i++)
    irt.rep[i] = symop[i];

  // identity
  ire.rep[0].E();

  // 4 C3's
  ire.rep[1].rotation(2.0*(double)M_PI/3.0);
  ire.rep[2] = ire.rep[1];
  ire.rep[3] = ire.rep[1];
  ire.rep[4] = ire.rep[1];

  ire.rep[5] = ire.rep[1].operate(ire.rep[1]);
  ire.rep[6] = ire.rep[5];
  ire.rep[7] = ire.rep[5];
  ire.rep[8] = ire.rep[5];
  
  // 3 C2's
  ire.rep[9].unit();
  ire.rep[10].unit();
  ire.rep[11].unit();

}

void
CharacterTable::th()
{
  int i,j;

  SymmetryOperation so;
  so.i();
  
  t_ops(symop);
  for (i=0; i < 12; i++)
    symop[i+12] = symop[i].operate(so);
  
  gamma_[0].init(g,1,"Ag");
  gamma_[1].init(g,1,"Au");

  for (i=0; i < 12; i++) {
    gamma_[0].rep[i][0][0] = 1.0;
    gamma_[1].rep[i][0][0] = 1.0;

    gamma_[0].rep[i+12][0][0] =  1.0;
    gamma_[1].rep[i+12][0][0] = -1.0;
  }

  IrreducibleRepresentation& ireg = gamma_[2];
  IrreducibleRepresentation& ireu = gamma_[3];

  IrreducibleRepresentation& irtg = gamma_[4];
  IrreducibleRepresentation& irtu = gamma_[5];

  ireg.init(g,2,"Eg");
  ireu.init(g,2,"Eu");
  ireg.complex_=1;
  ireu.complex_=1;

  irtg.init(g,3,"Tg");
  irtu.init(g,3,"Tu");
  irtg.nrot_=1;
  irtu.ntrans_=1;

  // the symmetry operation matrices form a basis for Tu.  Tg(g)=Tu(g) for
  // the proper rotations, and = -Tu(g) for the improper ones
  for (i=0; i < 12; i++) {
    irtg.rep[i] = symop[i];
    irtu.rep[i] = symop[i];

    irtg.rep[i+12] = symop[i];
    irtu.rep[i+12] = symop[i+12];
  }
    
  // identity
  ireg.rep[0].E();
  
  // 4 C3's
  ireg.rep[1].rotation(2.0*(double)M_PI/3.0);
  ireg.rep[2] = ireg.rep[1];
  ireg.rep[3] = ireg.rep[1];
  ireg.rep[4] = ireg.rep[1];

  // 4 C3^2's
  ireg.rep[5] = ireg.rep[1].operate(ireg.rep[1]);
  ireg.rep[6] = ireg.rep[5];
  ireg.rep[7] = ireg.rep[5];
  ireg.rep[8] = ireg.rep[5];

  // 3 C2's
  ireg.rep[9].unit();
  ireg.rep[10].unit();
  ireg.rep[11].unit();

  SymRep sr(2);
  sr.i();
  
  for (j=0; j < 12; j++) {
    ireu.rep[j] = ireg.rep[j];
    ireg.rep[j+12] = ireg.rep[j];
    ireu.rep[j+12] = ireg.rep[j].operate(sr);
  }
}

void
CharacterTable::td()
{
  // first get the T operations, then the Td operations
  t_ops(symop);
  td_ops(&symop[12]);
  
  int i;
  
  gamma_[0].init(g,1,"A1");
  gamma_[1].init(g,1,"A2");

  for (i=0; i < 12; i++) {
    gamma_[0].rep[i][0][0] = 1.0;
    gamma_[1].rep[i][0][0] = 1.0;

    gamma_[0].rep[i+12][0][0] =  1.0;
    gamma_[1].rep[i+12][0][0] = -1.0;
  }

  IrreducibleRepresentation& ire = gamma_[2];
  ire.init(g,2,"E");

  IrreducibleRepresentation& irt1 = gamma_[3];
  IrreducibleRepresentation& irt2 = gamma_[4];

  irt1.init(g,3,"T1");
  irt2.init(g,3,"T2");
  irt1.nrot_ = 1;
  irt2.ntrans_ = 1;

  // the symmetry operation matrices form a basis for T2.  T1(g)=T2(g) for
  // the proper rotations, and = -T2(g) for the improper ones
  SymmetryOperation so;
  so.i();
  
  for (i=0; i < 12; i++) {
    irt1.rep[i] = symop[i];
    irt2.rep[i] = symop[i];
    irt1.rep[i+12] = symop[i+12].operate(so);
    irt2.rep[i+12] = symop[i+12];
  }
  
  // identity
  ire.rep[0].E();

  // 4 C3's
  ire.rep[1].rotation(2.0*(double)M_PI/3.0);
  ire.rep[2] = ire.rep[1];
  ire.rep[3] = ire.rep[1];
  ire.rep[4] = ire.rep[1];

  // 4 C3^2's
  ire.rep[5] = ire.rep[1].operate(ire.rep[1]);
  ire.rep[6] = ire.rep[5];
  ire.rep[7] = ire.rep[5];
  ire.rep[8] = ire.rep[5];

  // 3 C2's
  ire.rep[9].unit();
  ire.rep[10].unit();
  ire.rep[11].unit();

  // 6 S4's
  ire.rep[12].c2_x();
  ire.rep[13].c2_x();

  ire.rep[14] = ire.rep[12].operate(ire.rep[1]);
  ire.rep[15] = ire.rep[14];
  
  ire.rep[16] = ire.rep[14].operate(ire.rep[1]);
  ire.rep[17] = ire.rep[16];

  for (i=18; i < 24; i++)
    ire.rep[i] = ire.rep[i-6];
}

void
CharacterTable::o()
{
  int i;
  
  // first get the T operations, then the O operations
  t_ops(symop);
  td_ops(&symop[12]);
  
  SymmetryOperation so;
  so.i();

  for (i=12; i < 24; i++)
    symop[i] = symop[i].operate(so);
  
  gamma_[0].init(g,1,"A1");
  gamma_[1].init(g,1,"A2");

  for (i=0; i < 12; i++) {
    gamma_[0].rep[i][0][0] = 1.0;
    gamma_[1].rep[i][0][0] = 1.0;

    gamma_[0].rep[i+12][0][0] =  1.0;
    gamma_[1].rep[i+12][0][0] = -1.0;
  }

  IrreducibleRepresentation& ire = gamma_[2];
  ire.init(g,2,"E");

  IrreducibleRepresentation& irt1 = gamma_[3];
  IrreducibleRepresentation& irt2 = gamma_[4];

  irt1.init(g,3,"T1");
  irt2.init(g,3,"T2");
  irt1.nrot_ = 1;
  irt1.ntrans_ = 1;

  // the symmetry operation matrices form a basis for T1.  T2(g)=T1(g) for
  // the proper rotations, and = -T1(g) for the improper ones
  
  for (i=0; i < 12; i++) {
    irt1.rep[i] = symop[i];
    irt2.rep[i] = symop[i];
    irt1.rep[i+12] = symop[i+12];
    irt2.rep[i+12] = symop[i+12].operate(so);
  }
  
  // identity
  ire.rep[0].E();

  // 4 C3's
  ire.rep[1].rotation(2.0*(double)M_PI/3.0);
  ire.rep[2] = ire.rep[1];
  ire.rep[3] = ire.rep[1];
  ire.rep[4] = ire.rep[1];

  // 4 C3^2's
  ire.rep[5] = ire.rep[1].operate(ire.rep[1]);
  ire.rep[6] = ire.rep[5];
  ire.rep[7] = ire.rep[5];
  ire.rep[8] = ire.rep[5];

  // 3 C2's
  ire.rep[9].unit();
  ire.rep[10].unit();
  ire.rep[11].unit();

  // 6 C4's
  ire.rep[12].c2_x();
  ire.rep[13].c2_x();

  ire.rep[14] = ire.rep[12].operate(ire.rep[1]);
  ire.rep[15] = ire.rep[14];
  
  ire.rep[16] = ire.rep[14].operate(ire.rep[1]);
  ire.rep[17] = ire.rep[16];

  // 6 C2's
  for (i=18; i < 24; i++)
    ire.rep[i] = ire.rep[i-6];
}

void CharacterTable::oh()
{
  int i,j;
  
  SymmetryOperation so;
  so.i();
  
  // first get the T operations, then the O operations, then the Th
  // operations, then the Td operations
  t_ops(symop);
  td_ops(&symop[36]);

  for (i=0; i < 12; i++) {
    symop[i+24] = symop[i].operate(so);
    symop[i+12] = symop[i+36].operate(so);
  }
  
  gamma_[0].init(g,1,"A1g");
  gamma_[1].init(g,1,"A2g");
  gamma_[5].init(g,1,"A1u");
  gamma_[6].init(g,1,"A2u");

  for (i=0; i < 12; i++) {
    gamma_[0].rep[i][0][0] = 1.0;
    gamma_[1].rep[i][0][0] = 1.0;
    gamma_[5].rep[i][0][0] = 1.0;
    gamma_[6].rep[i][0][0] = 1.0;

    gamma_[0].rep[i+12][0][0] =  1.0;
    gamma_[1].rep[i+12][0][0] = -1.0;
    gamma_[5].rep[i+12][0][0] =  1.0;
    gamma_[6].rep[i+12][0][0] = -1.0;

    gamma_[0].rep[i+24][0][0] =  1.0;
    gamma_[1].rep[i+24][0][0] =  1.0;
    gamma_[5].rep[i+24][0][0] = -1.0;
    gamma_[6].rep[i+24][0][0] = -1.0;

    gamma_[0].rep[i+36][0][0] =  1.0;
    gamma_[1].rep[i+36][0][0] = -1.0;
    gamma_[5].rep[i+36][0][0] = -1.0;
    gamma_[6].rep[i+36][0][0] =  1.0;
  }

  // the symmetry operation matrices form a basis for T1u.  T2u(g)=T1u(g) for
  // the proper rotations, and = -T1(g) for the improper ones.
  // T1g(g)=T1u(g) for the O part, and = -T1u(g) for the ixO part.
  // T2g(g)=T1g(g) for proper rotations and =-T1g(g) for improper

  gamma_[3].init(g,3,"T1g");
  gamma_[4].init(g,3,"T2g");
  gamma_[8].init(g,3,"T1u");
  gamma_[9].init(g,3,"T2u");

  gamma_[3].nrot_=1;
  gamma_[8].ntrans_=1;
  
  for (i=0; i < 12; i++) {
    gamma_[3].rep[i] = symop[i];
    gamma_[4].rep[i] = symop[i];
    gamma_[8].rep[i] = symop[i];
    gamma_[9].rep[i] = symop[i];
    
    gamma_[3].rep[i+12] = symop[i+12];
    gamma_[4].rep[i+12] = symop[i+12].operate(so);
    gamma_[8].rep[i+12] = symop[i+12];
    gamma_[9].rep[i+12] = symop[i+12].operate(so);
    
    gamma_[3].rep[i+24] = symop[i+24].operate(so);
    gamma_[4].rep[i+24] = symop[i+24].operate(so);
    gamma_[8].rep[i+24] = symop[i+24];
    gamma_[9].rep[i+24] = symop[i+24];
    
    gamma_[3].rep[i+36] = symop[i+36].operate(so);
    gamma_[4].rep[i+36] = symop[i+36];
    gamma_[8].rep[i+36] = symop[i+36];
    gamma_[9].rep[i+36] = symop[i+36].operate(so);
  }

  IrreducibleRepresentation& ireg = gamma_[2];
  IrreducibleRepresentation& ireu = gamma_[7];

  ireg.init(g,2,"Eg");
  ireu.init(g,2,"Eu");
    
  // identity
  ireg.rep[0].E();
  
  // 4 C3's
  ireg.rep[1].rotation(2.0*(double)M_PI/3.0);
  ireg.rep[2] = ireg.rep[1];
  ireg.rep[3] = ireg.rep[1];
  ireg.rep[4] = ireg.rep[1];

  // 4 C3^2's
  ireg.rep[5] = ireg.rep[1].operate(ireg.rep[1]);
  ireg.rep[6] = ireg.rep[5];
  ireg.rep[7] = ireg.rep[5];
  ireg.rep[8] = ireg.rep[5];

  // 3 C2's
  ireg.rep[9].unit();
  ireg.rep[10].unit();
  ireg.rep[11].unit();

  // 6 C4's
  ireg.rep[12].c2_x();
  ireg.rep[13].c2_x();

  ireg.rep[14] = ireg.rep[12].operate(ireg.rep[1]);
  ireg.rep[15] = ireg.rep[14];
  
  ireg.rep[16] = ireg.rep[14].operate(ireg.rep[1]);
  ireg.rep[17] = ireg.rep[16];

  // 6 C2's
  for (i=18; i < 24; i++)
    ireg.rep[i] = ireg.rep[i-6];

  SymRep sr(2);
  sr.i();
  
  for (j=0; j < 24; j++) {
    ireu.rep[j] = ireg.rep[j];
    ireg.rep[j+24] = ireg.rep[j];
    ireu.rep[j+24] = ireg.rep[j].operate(sr);
  }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
