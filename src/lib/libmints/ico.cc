//
// ico.cc --- implementation of icosahedral operations
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

// #include <util/misc/math.h>
#include <cmath>
#include <string.h>

#include <libmints/pointgrp.h>

using namespace psi;

// these are the operations which make up T
static void
i_ops(SymRep *t1rep, SymRep *t2rep, SymRep *grep, SymRep *hrep)
{
  int i;
  
  // identity
  t1rep[0].E();
  t2rep[0].E();
  grep[0].E();
  hrep[0].E();
    
  //
  // 12 C5's
  //
  // first the 2 C5's about the z axis
  t1rep[1].rotation(2.0*(double)M_PI/5.0);
  t1rep[2].rotation(8.0*(double)M_PI/5.0);
  
  t2rep[1] = t1rep[1].operate(t1rep[1]);
  t2rep[2] = t1rep[2].operate(t1rep[2]);

  grep[1].rotation(2.0*(double)M_PI/5.0);
  grep[2].rotation(8.0*(double)M_PI/5.0);
  
  hrep[1].rotation(2.0*(double)M_PI/5.0);
  hrep[2].rotation(8.0*(double)M_PI/5.0);
   
  // form rotation matrices for the C3 axis about the zx axis (these were
  // taken from turbomole version 2, which claims they were sort of inherited
  // from hondo
  SymRep t1so(3);
  SymRep gso(4);
  SymRep hso(5);

  double c2p5 = cos(2.0*(double)M_PI/5.0);
  double s2p5 = sin(2.0*(double)M_PI/5.0);
  double cosd = s2p5/((1.0-c2p5)*sqrt(3.0));
  double cosd2 = cosd*cosd;
  double sind2 = 1.0 - cosd2;
  double sind = sqrt(sind2);

  t1so[0][0] =  1.0 - 1.5*cosd2;
  t1so[1][0] =  0.5*sqrt(3.0)*cosd;
  t1so[2][0] =  1.5*cosd*sind;
  t1so[0][1] = -0.5*sqrt(3.0)*cosd;
  t1so[1][1] = -0.5;
  t1so[2][1] =  0.5*sqrt(3.0)*sind;
  t1so[0][2] =  1.5*cosd*sind;
  t1so[1][2] = -0.5*sqrt(3.0)*sind;
  t1so[2][2] =  1.0 - 1.5*sind2;

  gso[0][0] = (3.0*sqrt(5.0)+5.0)/20.0;
  gso[0][1] = cosd*sqrt(3.0)*(sqrt(5.0)-1.0)/4.0;
  gso[0][2] = 3.0*sqrt(5.0)/10.0;
  gso[0][3] = -sqrt(5.0-2.0*sqrt(5.0))*sqrt(5.0)/10.0;
  gso[1][0] = -gso[0][1];
  gso[1][1] = (1-sqrt(5.0))/4.0;
  gso[1][2] = cosd*sqrt(3.0)/2.0;
  gso[1][3] = cosd*sqrt(5-2*sqrt(5.0))*sqrt(3.0)/2.0;
  gso[2][0] = gso[0][2];
  gso[2][1] = -gso[1][2];
  gso[2][2] = (5-3*sqrt(5.0))/20.0;
  gso[2][3] = sqrt(5.0-2*sqrt(5.0))*(sqrt(5.0)+5)/20;
  gso[3][0] = -gso[0][3];
  gso[3][1] = gso[1][3];
  gso[3][2] = -gso[2][3];
  gso[3][3] = (sqrt(5.0)+1)/4.0;

  hso[0][0] = -1.0/5.0;
  hso[0][4] = sqrt(3.0)*(sqrt(5.0)+1)/10.0;
  hso[0][3] = 3.0*cosd*(3.0*sqrt(5.0)-5.0)/10.0;
  hso[0][2] = 3.0*cosd*(5.0-sqrt(5.0))/10.0;
  hso[0][1] = sqrt(3.0)*(sqrt(5.0)-1.0)/10.0;
  hso[4][0] = hso[0][4];
  hso[4][4] = (2.0*sqrt(5.0)+1.0)/10.0;
  hso[4][3] = sqrt(3.0)*cosd*(5.0-2.0*sqrt(5.0))/10.0;
  hso[4][2] = sqrt(3.0)*cosd*(5.0-3.0*sqrt(5.0))/5.0;
  hso[4][1] = 2.0/5.0;
  hso[3][0] = -hso[0][3];
  hso[3][4] = -hso[4][3];
  hso[3][3] = -1.0/2.0;
  hso[3][2] = 0.0;
  hso[3][1] = sqrt(3.0)*cosd*(5.0-sqrt(5.0))/5.0;
  hso[2][0] = -hso[0][2];
  hso[2][4] = -hso[4][2];
  hso[2][3] = 0.0;
  hso[2][2] = -1.0/2.0;
  hso[2][1] = -sqrt(3.0)*sqrt(5.0)*cosd/10.0;
  hso[1][0] = hso[0][1];
  hso[1][4] = hso[4][1];
  hso[1][3] = -hso[3][1];
  hso[1][2] = -hso[2][1];
  hso[1][1] = (1.0-2.0*sqrt(5.0))/10.0;
  
  // now rotate the first C5's by 2pi/3 degrees about the zx axis (sort of)
  t1rep[3] = t1rep[1].transform(t1so);
  t1rep[4] = t1rep[2].transform(t1so);

  grep[3] = grep[1].transform(gso);
  grep[4] = grep[2].transform(gso);

  hrep[3] = hrep[1].transform(hso);
  hrep[4] = hrep[2].transform(hso);

  // rotate twice to get the first one aligned along the x axis
  t1rep[3] = t1rep[3].transform(t1rep[1]).transform(t1rep[1]);
  t1rep[4] = t1rep[4].transform(t1rep[1]).transform(t1rep[1]);

  grep[3] = grep[3].transform(grep[1]).transform(grep[1]);
  grep[4] = grep[4].transform(grep[1]).transform(grep[1]);

  hrep[3] = hrep[3].transform(hrep[1]).transform(hrep[1]);
  hrep[4] = hrep[4].transform(hrep[1]).transform(hrep[1]);

  t2rep[3] = t1rep[4].operate(t1rep[4]);
  t2rep[4] = t1rep[3].operate(t1rep[3]);

  t2rep[13] = t1rep[2];
  t2rep[14] = t1rep[1];

  t2rep[15] = t1rep[3];
  t2rep[16] = t1rep[4];
  
  // and then rotate those by 2pi/5 about the z axis 4 times
  for (i=5; i < 13; i++) {
    t1rep[i] = t1rep[i-2].transform(t1rep[1]);
    grep[i] = grep[i-2].transform(grep[1]);
    hrep[i] = hrep[i-2].transform(hrep[1]);

    t2rep[i] = t2rep[i-2].transform(t2rep[1]);
    t2rep[i+12] = t2rep[i+10].transform(t2rep[1]);
  }

  //
  // 12 C5^2's
  //
  // get these from operating on each of the C5's with itself
  for (i=13; i < 25; i++) {
    t1rep[i] = t1rep[i-12].operate(t1rep[i-12]);
    grep[i] = grep[i-12].operate(grep[i-12]);
    hrep[i] = hrep[i-12].operate(hrep[i-12]);
  }

  //
  // 20 C3's
  //
  // first we have 2 C3's about the zx axis
  t1rep[25] = t1so;
  t1rep[26] = t1so.operate(t1so);
  
  grep[25] = gso;
  grep[26] = gso.operate(gso);
  
  hrep[25] = hso;
  hrep[26] = hso.operate(hso);
  
  // and then rotate those by 2pi/5 about the z axis 4 times
  for (i=27; i < 35; i++) {
    t1rep[i] = t1rep[i-2].transform(t1rep[1]);
    grep[i] = grep[i-2].transform(grep[1]);
    hrep[i] = hrep[i-2].transform(hrep[1]);
  }

  // now rotate one of the above C3's by 2pi/3 about the zx axis
  t1rep[35] = t1rep[27].transform(t1so);
  t1rep[36] = t1rep[28].transform(t1so);

  grep[35] = grep[27].transform(gso);
  grep[36] = grep[28].transform(gso);
                      
  hrep[35] = hrep[27].transform(hso);
  hrep[36] = hrep[28].transform(hso);

  // and then rotate those by 2pi/5 about the z axis 4 times
  for (i=37; i < 45; i++) {
    t1rep[i] = t1rep[i-2].transform(t1rep[1]);
    grep[i] = grep[i-2].transform(grep[1]);
    hrep[i] = hrep[i-2].transform(hrep[1]);
  }

  t2rep[25] = t1rep[35];
  t2rep[26] = t1rep[36];
  
  for (i=27; i < 35; i++)
    t2rep[i] = t2rep[i-2].transform(t2rep[1]);
  
  t2rep[35] = t1rep[26];
  t2rep[36] = t1rep[25];
  
  for (i=37; i < 45; i++)
    t2rep[i] = t2rep[i-2].transform(t2rep[1]);

  //
  // 15 C2's
  //
  // first we have a C2 about the y axis
  t1rep[45][0][0] = -1.0;
  t1rep[45][1][1] =  1.0;
  t1rep[45][2][2] = -1.0;

  t2rep[45] = t1rep[45];
  
  grep[45][0][0] = -1.0;
  grep[45][1][1] =  1.0;
  grep[45][2][2] = -1.0;
  grep[45][3][3] =  1.0;
  
  hrep[45][0][0] =  1.0;
  hrep[45][1][1] =  1.0;
  hrep[45][2][2] = -1.0;
  hrep[45][3][3] = -1.0;
  hrep[45][4][4] =  1.0;
  
  // and rotate that by 2pi/5 about the z axis 4 times
  for (i=46; i < 50; i++) {
    t1rep[i] = t1rep[i-1].transform(t1rep[1]);
    t2rep[i] = t2rep[i-1].transform(t2rep[1]);
    grep[i] = grep[i-1].transform(grep[1]);
    hrep[i] = hrep[i-1].transform(hrep[1]);
  }

  // now take the C2 about the y axis and rotate it by 2pi/3 about the zx axis
  t1rep[50] = t1rep[45].transform(t1so);
  grep[50] = grep[45].transform(gso);
  hrep[50] = hrep[45].transform(hso);

  // align this c2 along the x axis
  t1rep[50] = t1rep[50].transform(t1rep[2]).transform(t1rep[2]);
  grep[50] = grep[50].transform(grep[2]).transform(grep[2]);
  hrep[50] = hrep[50].transform(hrep[2]).transform(hrep[2]);

  // and rotate that by 2pi/5 about the z axis 4 times
  for (i=51; i < 55; i++) {
    t1rep[i] = t1rep[i-1].transform(t1rep[1]);
    grep[i] = grep[i-1].transform(grep[1]);
    hrep[i] = hrep[i-1].transform(hrep[1]);
  }

  // finally, take a C2 about the y axis, and rotate it by 2pi/3 about the
  // xz axis, and align it along the x axis
  t1rep[55] = t1rep[45].transform(t1rep[35]).transform(t1rep[1]);
  grep[55] = grep[45].transform(grep[35]).transform(grep[1]);
  hrep[55] = hrep[45].transform(hrep[35]).transform(hrep[1]);

  // and then rotate that by 2pi/5 about the z axis 4 times
  for (i=56; i < 60; i++) {
    t1rep[i] = t1rep[i-1].transform(t1rep[1]);
    grep[i] = grep[i-1].transform(grep[1]);
    hrep[i] = hrep[i-1].transform(hrep[1]);
  }

  t2rep[50] = t1rep[55];
  t2rep[55] = t1rep[50];
  
  for (i=51; i < 55; i++) {
    t2rep[i] = t2rep[i-1].transform(t2rep[1]);
    t2rep[i+5] = t2rep[i+4].transform(t2rep[1]);
  }
}

void
CharacterTable::i()
{
  int i;

  IrreducibleRepresentation& ira = gamma_[0];
  IrreducibleRepresentation& ir1 = gamma_[1];
  IrreducibleRepresentation& ir2 = gamma_[2];
  IrreducibleRepresentation& irg = gamma_[3];
  IrreducibleRepresentation& irh = gamma_[4];

  ira.init(g,1,"A");
  ir1.init(g,3,"T1");
  ir2.init(g,3,"T2");
  irg.init(g,4,"G");
  irh.init(g,5,"H");

  // i_ops gives us all the symmetry operations we need
  i_ops(ir1.rep, ir2.rep, irg.rep, irh.rep);
    
  ir1.nrot_ = 1;
  ir1.ntrans_ = 1;

  for (i=0; i < g; i++) {
    ira.rep[i][0][0] = 1.0;
    symop[i] = ir1.rep[i];
  }
}


void CharacterTable::ih()
{
  int i;

  IrreducibleRepresentation& irag = gamma_[0];
  IrreducibleRepresentation& ir1g = gamma_[1];
  IrreducibleRepresentation& ir2g = gamma_[2];
  IrreducibleRepresentation& irgg = gamma_[3];
  IrreducibleRepresentation& irhg = gamma_[4];

  IrreducibleRepresentation& irau = gamma_[5];
  IrreducibleRepresentation& ir1u = gamma_[6];
  IrreducibleRepresentation& ir2u = gamma_[7];
  IrreducibleRepresentation& irgu = gamma_[8];
  IrreducibleRepresentation& irhu = gamma_[9];

  irag.init(g,1,"Ag");
  ir1g.init(g,3,"T1g");
  ir2g.init(g,3,"T2g");
  irgg.init(g,4,"Gg");
  irhg.init(g,5,"Hg");

  irau.init(g,1,"Au");
  ir1u.init(g,3,"T1u");
  ir2u.init(g,3,"T2u");
  irgu.init(g,4,"Gu");
  irhu.init(g,5,"Hu");

  // i_ops gives us all the symmetry operations we need
  i_ops(ir1g.rep, ir2g.rep, irgg.rep, irhg.rep);
    
  ir1g.nrot_ = 1;
  ir1u.ntrans_ = 1;

  SymRep ti(3), gi(4), hi(5);
  ti.i();
  gi.i();
  hi.i();
  
  for (i=0; i < g/2; i++) {
    irag.rep[i][0][0] = 1.0;
    irau.rep[i][0][0] = 1.0;

    irag.rep[i+60][0][0] =  1.0;
    irau.rep[i+60][0][0] = -1.0;

    ir1g.rep[i+60] = ir1g.rep[i];
    ir2g.rep[i+60] = ir2g.rep[i];
    irgg.rep[i+60] = irgg.rep[i];
    irhg.rep[i+60] = irhg.rep[i];
    
    ir1u.rep[i] = ir1g.rep[i];
    ir2u.rep[i] = ir2g.rep[i];
    irgu.rep[i] = irgg.rep[i];
    irhu.rep[i] = irhg.rep[i];
    
    ir1u.rep[i+60] = ir1g.rep[i].operate(ti);
    ir2u.rep[i+60] = ir2g.rep[i].operate(ti);
    irgu.rep[i+60] = irgg.rep[i].operate(gi);
    irhu.rep[i+60] = irhg.rep[i].operate(hi);
    
    symop[i] = ir1u.rep[i];
    symop[i+60] = ir1u.rep[i+60];
  }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
