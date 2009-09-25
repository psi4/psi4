/*! \file
    \ingroup INTDER
    \brief Enter brief description of file here 
*/
#include "intco.h"
#include "params.h"
#include "displacements.h"
#include "atom.h"
#include "bmat.h"
#define EXTERN
#include "globals.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>
#include <physconst.h>
#include <psifiles.h>

using namespace psi::intder;

namespace psi { namespace intder {
extern InternalCoordinates gIntCo;
extern Params gParams;
extern Displacements gDisplacements;
extern BMat gBMat;
}}

void InternalCoordinate::printInfo()
{
  fprintf(outfile, "\tEmpty internal coordinate.\n");
}

void Stretch::printInfo()
{
  fprintf(outfile, "\tID#%d STRE (%d %d)\n", id, atomA+1, atomB+1);
}

// Conversion of VECT1, page 13
void Stretch::Vect(int disp, int k1, int k2, double *s1, double *dist, int pflag)
{
  int i = 0;
  double c = 0.0;
  Atom* a1 = NULL;
  Atom* a2 = NULL;
  *dist = 0;

  // Get the atoms involved
  a1 = gDisplacements.displacement(disp)->atom(k1);
  a2 = gDisplacements.displacement(disp)->atom(k2);

  s1[0] = a1->getX() - a2->getX();
  s1[1] = a1->getY() - a2->getY();
  s1[2] = a1->getZ() - a2->getZ();

  *dist = s1[0]*s1[0];
  *dist += s1[1]*s1[1];
  *dist += s1[2]*s1[2];
  *dist = sqrt(*dist);
  c = 1.0 / *dist;

  s1[0] *= c;
  s1[1] *= c;
  s1[2] *= c;
  
  if(pflag)
    fprintf(outfile, "Stretch (%d, %d) S: [%10.7f %10.7f %10.7f]\tDistance (Ang): %lf\n", k1, k2, s1[0], s1[1], s1[2], *dist);
}

void Bend::printInfo()
{
  fprintf(outfile, "\tID#%d BEND (%d %d %d)\n", id, atomA+1, atomB+1, atomC+1);
}

// Conversion of VECT2, A-B-C bending page 13
void Bend::Vect(int disp, int k1, int k2, int k3, double *s1, double *s2, double *s3, double *theta, int pflag)
{
  int i = 0;
  double c = 0.0;
  double r12, r32;
  double cos_prod;
  double sin_prod;
  double *e23;
  double *e21;

  *theta = 0;

  e21 = new double[3];
  e23 = new double[3];

  Stretch::Vect(disp, k1, k2, e21, &r12);
  Stretch::Vect(disp, k3, k2, e23, &r32);
    
  cos_prod = dot_prod(e21,e23);
  sin_prod = sqrt(1.0 - pow(cos_prod,2));
  *theta = acos(cos_prod);
     
  for(i = 0; i < 3; i++) {
    s1[i] = (cos_prod * e21[i] - e23[i]) / (r12 * sin_prod);
    s3[i] = (cos_prod * e23[i] - e21[i]) / (r32 * sin_prod);
    s2[i] = -s1[i] - s3[i];
  }

  if(pflag)
    fprintf(outfile, "Bend (%d, %d, %d) \t\tAngle (deg): %lf\n", k1, k2, k3, *theta*180/_pi);
  //fprintf(outfile, "\t\t S1: [%10.7f %10.7f %10.7f]\n", s1[0], s1[1], s1[2]);
  //fprintf(outfile, "\t\t S2: [%10.7f %10.7f %10.7f]\n", s2[0], s2[1], s2[2]);
  //fprintf(outfile, "\t\t S3: [%10.7f %10.7f %10.7f]\n", s3[0], s3[1], s3[2]);
  

  delete[] e21;
  delete[] e23;
}

void Linear1::printInfo()
{
  fprintf(outfile, "\tID#%d LIN1 (%d %d %d %d)\n", id, atomA+1, atomB+1, atomC+1, atomD+1);
}

void Linear1::Vect(int disp, int k1, int k2, int k3, int k4, double *s1, double *s2, double *s3, double *theta, int pflag)
{
  int i;
  double *e21;
  double *e23;
  double *e2m;
  double *ea;
  double *s4;
  double *s5;
  double t21, t23, w, w1, w2, stheta, ctheta, ttheta;

  e21 = new double[3];
  e23 = new double[3];
  e2m = new double[3];
  ea  = new double[3];

  Stretch::Vect(disp, k1, k2, e21, &t21);
  Stretch::Vect(disp, k3, k2, e23, &t23);
  Stretch::Vect(disp, k4, k2, e2m, &w);

  w1 = dot_prod(e21, e2m);
  w2 = sqrt(1.0 - pow(w1, 2.0));

  ea[0] = (w1 * e21[0] - e2m[0]) / w2;
  ea[1] = (w1 * e21[1] - e2m[1]) / w2;
  ea[2] = (w1 * e21[2] - e2m[2]) / w2;

  delete[] e2m;
  e2m = vect_prod(e23, e21);
  stheta = dot_prod(ea, e2m);
  *theta = asin(stheta);
  ctheta = cos(*theta);
  ttheta = stheta/ctheta;

  s4 = vect_prod(ea, e23);
  s5 = vect_prod(ea, e21);

  for (i=0; i<3; i++) {
    s1[i] =  (s4[i]/ctheta - e21[i]*ttheta) / t21;
    s3[i] = -(s5[i]/ctheta + e23[i]*ttheta) / t23;
    s2[i] = -(s1[i] + s3[i]);
  }

  if(pflag)
    fprintf(outfile, "Linear1bend (%i, %i, %i) dummy = %i\t\tAngle (deg): %lf\n", k1, k2, k3, k4, *theta*180/_pi);

  free(e2m);
  delete[] e21;
  delete[] e23;
  delete[] e2m;
  delete[] ea;
  delete[] s4;
  delete[] s5;
}

void OutOfPlane::printInfo()
{
  fprintf(outfile, "\tID#%d OUT  (%d %d %d %d)\n", id, atomA+1, atomB+1, atomC+1, atomD+1);
}

//VECT5
void OutOfPlane::Vect(int disp, int k1, int k2, int k3, int k4, double *s1, double *s2, double *s3, double *s4, double *theta, int pflag)
{
  int i = 0;
  double c = 0.0;
  double r12, r32, r42;
  double cprod1, cprod2;
  double *vprod1, *vprod2, *vprod3;
  double *e24, *e23, *e21;
  double phi;
  double *k3sbend;
  double *v2;
  double *k4sbend;

  e21 = new double[3];
  e23 = new double[3];
  e24 = new double[3];
  k3sbend = new double[3];
  k4sbend = new double[3];
  v2 = new double[3];

  Stretch::Vect(disp, k1, k2, e21, &r12);
  Stretch::Vect(disp, k3, k2, e23, &r32);
  Stretch::Vect(disp, k4, k2, e24, &r42);

  vprod1 = vect_prod(e23, e24);

  cprod1 = dot_prod(e21, e23);
  cprod2 = dot_prod(e21, e24);

  vprod2 = vect_prod(e24, e21);
  vprod3 = vect_prod(e21, e23);

  Bend::Vect(disp, k3, k2, k4, k3sbend, v2, k4sbend, &phi);
  *theta = asin(dot_prod(e21, vprod1)/sin(phi));
  if ((cprod1 + cprod2) > 0.0)  {
    if (*theta >= 0) 
      *theta = (_pi - *theta);
    else 
      *theta = (-_pi - *theta);
  }

  for (i=0; i<3; i++) {
    s1[i] = (vprod1[i] / (cos(*theta) * sin(phi)) - (e21[i] * tan(*theta))) / r12;
    s3[i] = (1 / (r32 * cos(*theta) * sin(phi))) * (vprod2[i] + k4sbend[i] * r42 * sin(*theta));
    s4[i] = (1 / (r42 * cos(*theta) * sin(phi))) * (vprod3[i] + k3sbend[i] * r32 * sin(*theta));
    s2[i] = -(s1[i] + s3[i] + s4[i]); 
  }

  if(pflag)
    fprintf(outfile, "Out (%d, %d, %d, %d) S: [%10.7f %10.7f %10.7f]\tAngle (deg): %lf\n", k1, k2, k3, k4, s1[0], s1[1], s1[2], *theta*180/_pi);

  free(vprod2);
  free(vprod3);
  delete[] e21;
  delete[] e23;
  delete[] e24;
  delete[] k3sbend;
  delete[] k4sbend;
  delete[] v2;
}

void Torsion::printInfo()
{
  fprintf(outfile, "\tID#%d TORS (%d %d %d %d)\n", id, atomA+1, atomB+1, atomC+1, atomD+1);
}

void Torsion::Vect(int disp, int k1, int k2, int k3, int k4, double *s1, double *s2, double *s3, double *s4, double *theta, int pflag)
{
  double *e21, *e32, *e43;
  double *s5, *s6;
  double t21, t32, t43;
  double w1, w2, w3, w4, w5, w6, cp2, cp3, sp2, sp3;
  int i;

  e21 = new double[3];
  e32 = new double[3];
  e43 = new double[3];

  Stretch::Vect(disp, k1, k2, e21, &t21);
  Stretch::Vect(disp, k2, k3, e32, &t32);
  Stretch::Vect(disp, k3, k4, e43, &t43);

  s5 = vect_prod(e21, e32);
  s6 = vect_prod(e43, e32);

  w2 = dot_prod(e21, e32);
  w3 = dot_prod(e32, e43);

  cp2 = -w2;
  cp3 = -w3;
  sp2 = sqrt(1.0 - pow(cp2, 2.0));
  sp3 = sqrt(1.0 - pow(cp3, 2.0));

  w2 = dot_prod(e21, s6);
  w3 = dot_prod(s5, s6);

  w3 = -w3;
  *theta = asin(w2 / (sp2*sp3));
  if (w3 < 0.0)
    *theta = _pi - (*theta);

  w1 = 1.0 / (t21*sp2*sp2);
  w2 = 1.0 / (t43*sp3*sp3);

  for (i=0; i<3; i++) {
    s1[i] = -w1 * s5[i];
    s4[i] = -w2 * s6[i];
  }

  w3 = (t32-t21*cp2)*w1/t32;
  w4 = cp3/(t32*sp3*sp3);
  w5 = (t32-t43*cp3)*w2/t32;
  w6 = cp2/(t32*sp2*sp2);

  for (i=0; i<3; i++) {
    s2[i] = w3*s5[i] + w4*s6[i];
    s3[i] = w5*s6[i] + w6*s5[i];
  }

  if(pflag)
    fprintf(outfile, "Torsional (%d %d %d %d)\t\tAngle (deg): %lf\n", k1, k2, k3, k4, *theta*180/_pi);

  delete[] e21;
  delete[] e32;
  delete[] e43;
  free(s5);
  free(s6);
}

void Spf::printInfo()
{
  fprintf(outfile, "\tID#%d SPF  (%d %d)\n", id, atomA+1, atomB+1);
}

void LinearX::printInfo()
{
  fprintf(outfile, "\tID#%d LINX (%d %d %d %d)\n", id, atomA+1, atomB+1, atomC+1, atomD+1);
}

void LinearX::Vect(int disp, int k1, int k2, int k3, int k4, double *s1, double *s2, double *s3, double *s4, double *theta, int pflag)
{
  double *e1, *e2, *e3, *e4, *e32, *e34;
  double **h11, **h21, **h31, **h22, **h32, **h33;
  double **h3;
  double r23, r43, phi123;
  double tout, cosy;
  double dprod34_3;
  int i, j;

  e1 = new double[3];
  e2 = new double[3];
  e3 = new double[3];
  e4 = new double[3];
  e32 = new double[3];
  e34 = new double[3];
  h3 = init_matrix(3, 3);
  h11 = init_matrix(3, 3);
  h21 = init_matrix(3, 3);
  h31 = init_matrix(3, 3);
  h22 = init_matrix(3, 3);
  h32 = init_matrix(3, 3);
  h33 = init_matrix(3, 3);

  Stretch::Vect(disp, k2, k3, e32, &r23);
  Stretch::Vect(disp, k4, k3, e34, &r43);

  Bend::Vect(disp, k1, k2, k3, e1, e2, e3, &phi123);
  
  dprod34_3 = dot_prod(e34, e3);
  *theta = -(r23 / dprod34_3);

  Stretch::Hijs(disp, k3, k4, h3);
  Bend::Hijs(disp, k1, k2, k3, h11, h21, h31, h22, h32, h33);

  for (i=0; i<3; i++) {
    s1[i] = 0.0;
    s4[i] = 0.0;
    s2[i] = *theta*e32[i]/r23;

    for (j=0; j<3; j++) {
      s1[i] = s1[i] - r23*e34[j]*h31[j][i];
      s2[i] = s2[i] - r23*e34[j]*h32[j][i];
      s4[i] = s4[i] - r23*e3[j]*h3[j][i];
    }
    s3[i] = -s1[i] - s2[i] - s4[i];
  }

  if(pflag)
    fprintf(outfile, "LinearX (%d %d %d %d)\t\tAngle (deg): %lf\n", k1, k2, k3, k4, *theta*180/_pi);

  delete[] e1;
  delete[] e2;
  delete[] e3;
  delete[] e4;
  delete[] e32;
  delete[] e34;
  free(h3);
  free(h11);
  free(h21);
  free(h31);
  free(h22);
  free(h32);
  free(h33);
}

void LinearY::printInfo()
{
  fprintf(outfile, "\tID#%d LINY (%d %d %d %d)\n", id, atomA+1, atomB+1, atomC+1, atomD+1);
}

void LinearY::Vect(int disp, int k1, int k2, int k3, int k4, double *s1, double *s2, double *s3, double *s4, double *theta, int pflag)
{
  double *e1, *e2, *e3, *e4;
  double tout, cosy;
  int i;

  e1 = new double[3];
  e2 = new double[3];
  e3 = new double[3];
  e4 = new double[3];

  OutOfPlane::Vect(disp, k1, k2, k3, k4, e4, e3, e2, e1, &tout);
  *theta = -sin(tout);
  cosy = cos(tout);

  for (i=0; i<3; i++) {
    s2[i] = -cosy * e2[i];
    s3[i] = -cosy * e3[i];
    s4[i] = -cosy * e4[i];
    s1[i] = -cosy * e1[i];
  }

  if(pflag)
    fprintf(outfile, "LinearY (%d %d %d %d)\t\tAngle (deg): %lf\n", k1, k2, k3, k4, *theta*180/_pi);

  delete[] e1;
  delete[] e2;
  delete[] e3;
  delete[] e4;
}

void Rcom::printInfo()
{
  fprintf(outfile, "\tID#%d RCOM (%d %d %d %d)\n", id, atomA+1, atomB+1);
}

InternalCoordinates::~InternalCoordinates()
{
  int index;

  for (index = 0; index < vectorInternals.size(); index++)
    delete vectorInternals[index];
}

void InternalCoordinates::loadInternalCoordinates()
{
  int i = 0, j = 0;
  int cnt = 0;
  int a, b, c, d, e;  // temp storage variables
  
  // Should we be reading the intco section?
  if (!gParams.simplesPresent)
    return;

  if (ip_exist("STRE", 0)) {
    ip_count("STRE", &cnt, 0);
    fprintf(outfile, "%d stretches found.\n", cnt);

    for (i=0; i < cnt; i++) {
      ip_count("STRE", &j, 1, i);

      if (j != 3) {
        fprintf(outfile, "STRE %d is of wrong dimension.\n", i+1);
        fprintf(outfile, "Found %d needed 3\n", j);
        exit(2);
      }

      ip_data("STRE", "%d", &a, 2, i, 0);
      ip_data("STRE", "%d", &b, 2, i, 1);
      ip_data("STRE", "%d", &c, 2, i, 2);

      // Add to the coordinate list here, atoms need to be zero based
      // Reorders the internal coordinates so that dummy atoms are last.
      b--;
      c--;
      addInternalCoordinate(new Stretch(a, gParams.moved_dummy[b], gParams.moved_dummy[c]));
    }
  }

  if (ip_exist("BEND", 0)) {
    ip_count("BEND", &cnt, 0);
    fprintf(outfile, "%d bends found.\n", cnt);

    for (i=0; i<cnt; i++) {
      ip_count("BEND", &j, 1, i);

      if (j != 4) {
        fprintf(outfile, "BEND %d is of wrong dimension.\n", i+1);
        fprintf(outfile, "Found %d needed 4\n", j);
        exit(2);
      }

      ip_data("BEND", "%d", &a, 2, i, 0);
      ip_data("BEND", "%d", &b, 2, i, 1);
      ip_data("BEND", "%d", &c, 2, i, 2);
      ip_data("BEND", "%d", &d, 2, i, 3);

      // Add to the coordinate list here, atoms need to be zero based
      b--;
      c--;
      d--;
      addInternalCoordinate(new Bend(a, gParams.moved_dummy[b], gParams.moved_dummy[c], gParams.moved_dummy[d]));
    }
  }

  if (ip_exist("LIN1", 0)) {
    ip_count("LIN1", &cnt, 0);
    fprintf(outfile, "%d linear1s found.\n", cnt);

    for (i=0; i<cnt; i++) {
      ip_count("LIN1", &j, 1, i);

      if (j != 5) {
        fprintf(outfile, "LIN1 %d is of wrong dimension.\n", i+1);
        fprintf(outfile, "Found %d needed 5\n", j);
        exit(2);
      }

      ip_data("LIN1", "%d", &a, 2, i, 0);
      ip_data("LIN1", "%d", &b, 2, i, 1);
      ip_data("LIN1", "%d", &c, 2, i, 2);
      ip_data("LIN1", "%d", &d, 2, i, 3);
      ip_data("LIN1", "%d", &e, 2, i, 4);

      // Add to the coordinate list here, atoms need to be zero based
      b--;
      c--;
      d--;
      e--;
      addInternalCoordinate(new Linear1(a, gParams.moved_dummy[b], gParams.moved_dummy[c], gParams.moved_dummy[d], gParams.moved_dummy[e]));
    }
  }

  if (ip_exist("OUT", 0)) {
    ip_count("OUT", &cnt, 0);
    fprintf(outfile, "%d out-of-planes found.\n", cnt);

    for (i=0; i<cnt; i++) {
      ip_count("OUT", &j, 1, i);

      if (j != 5) {
        fprintf(outfile, "OUT %d is of wrong dimension.\n", i+1);
        fprintf(outfile, "Found %d needed 5\n", j);
        exit(2);
      }

      ip_data("OUT", "%d", &a, 2, i, 0);
      ip_data("OUT", "%d", &b, 2, i, 1);
      ip_data("OUT", "%d", &c, 2, i, 2);
      ip_data("OUT", "%d", &d, 2, i, 3);
      ip_data("OUT", "%d", &e, 2, i, 4);

      // Add to the coordinate list here, atoms need to be zero based
      b--;
      c--;
      d--;
      e--;
      addInternalCoordinate(new OutOfPlane(a, gParams.moved_dummy[b], gParams.moved_dummy[c], gParams.moved_dummy[d], gParams.moved_dummy[e]));
    }
  }

  if (ip_exist("TORS", 0)) {
    ip_count("TORS", &cnt, 0);
    fprintf(outfile, "%d torsionals found.\n", cnt);

    for (i=0; i<cnt; i++) {
      ip_count("TORS", &j, 1, i);

      if (j != 5) {
        fprintf(outfile, "TORS %d is of wrong dimension.\n", i+1);
        fprintf(outfile, "Found %d needed 5\n", j);
        exit(2);
      }

      ip_data("TORS", "%d", &a, 2, i, 0);
      ip_data("TORS", "%d", &b, 2, i, 1);
      ip_data("TORS", "%d", &c, 2, i, 2);
      ip_data("TORS", "%d", &d, 2, i, 3);
      ip_data("TORS", "%d", &e, 2, i, 4);

      // Add to the coordinate list here, atoms need to be zero based
      b--;
      c--;
      d--;
      e--;
      addInternalCoordinate(new Torsion(a, gParams.moved_dummy[b], gParams.moved_dummy[c], gParams.moved_dummy[d], gParams.moved_dummy[e]));
    }
  }

  if (ip_exist("SPF", 0)) {
    ip_count("SPF", &cnt, 0);
    fprintf(outfile, "%d SPFs found.\n", cnt);

    for (i=0; i<cnt; i++) {
      ip_count("SPF", &j, 1, i);

      if (j != 3) {
        fprintf(outfile, "SPF %d is of wrong dimension.\n", i+1);
        fprintf(outfile, "Found %d needed 3\n", j);
        exit(2);
      }

      ip_data("SPF", "%d", &a, 2, i, 0);
      ip_data("SPF", "%d", &b, 2, i, 1);
      ip_data("SPF", "%d", &c, 2, i, 2);

      // Add to the coordinate list here, atoms need to be zero based
      b--;
      c--;
      addInternalCoordinate(new Spf(a, gParams.moved_dummy[b], gParams.moved_dummy[c]));
    }
  }

  if (ip_exist("LINX", 0)) {
    ip_count("LINX", &cnt, 0);
    fprintf(outfile, "%d linear X's found.\n", cnt);

    for (i=0; i<cnt; i++) {
      ip_count("LINX", &j, 1, i);

      if (j != 5) {
        fprintf(outfile, "LINX %d is of wrong dimension.\n", i+1);
        fprintf(outfile, "Found %d needed 5\n", j);
        exit(2);
      }

      ip_data("LINX", "%d", &a, 2, i, 0);
      ip_data("LINX", "%d", &b, 2, i, 1);
      ip_data("LINX", "%d", &c, 2, i, 2);
      ip_data("LINX", "%d", &d, 2, i, 3);
      ip_data("LINX", "%d", &e, 2, i, 4);

      // Add to the coordinate list here, atoms need to be zero based
      b--;
      c--;
      d--;
      e--;
      addInternalCoordinate(new LinearX(a, gParams.moved_dummy[b], gParams.moved_dummy[c], gParams.moved_dummy[d], gParams.moved_dummy[e]));
    }
  }

  if (ip_exist("LINY", 0)) {
    ip_count("LINY", &cnt, 0);
    fprintf(outfile, "%d linear Y's found.\n", cnt);

    for (i=0; i<cnt; i++) {
      ip_count("LINY", &j, 1, i);

      if (j != 5) {
        fprintf(outfile, "LINY %d is of wrong dimension.\n", i+1);
        fprintf(outfile, "Found %d needed 5\n", j);
        exit(2);
      }

      ip_data("LINY", "%d", &a, 2, i, 0);
      ip_data("LINY", "%d", &b, 2, i, 1);
      ip_data("LINY", "%d", &c, 2, i, 2);
      ip_data("LINY", "%d", &d, 2, i, 3);
      ip_data("LINY", "%d", &e, 2, i, 4);

      // Add to the coordinate list here, atoms need to be zero based
      b--;
      c--;
      d--;
      e--;
      addInternalCoordinate(new LinearY(a, gParams.moved_dummy[b], gParams.moved_dummy[c], gParams.moved_dummy[d], gParams.moved_dummy[e]));
    }
  }

  if (ip_exist("RCOM", 0)) {
    ip_count("RCOM", &cnt, 0);
    fprintf(outfile, "%d RCOMs found.\n", cnt);

    for (i=0; i < cnt; i++) {
      ip_count("RCOM", &j, 1, i);

      if (j != 3) {
        fprintf(outfile, "RCOM %d is of wrong dimension.\n", i+1);
        fprintf(outfile, "Found %d needed 3\n", j);
        exit(2);
      }

      ip_data("RCOM", "%d", &a, 2, i, 0);
      ip_data("RCOM", "%d", &b, 2, i, 1);
      ip_data("RCOM", "%d", &c, 2, i, 2);

      // Add to the coordinate list here, atoms need to be zero based
      b--;
      c--;
      addInternalCoordinate(new Rcom(a, gParams.moved_dummy[b], gParams.moved_dummy[c]));
    }
  }
  
  printInternalCoordinates();
}

int InternalCoordinates::intcoSwitch(double disp, int r, int *Aptr, int *Bptr, int *Cptr, int *Dptr, 
				      double *s1, double *s2, double *s3, double *s4, double *vptr)
{
  double value;
  int i, j;
  int atomA, atomB, atomC, atomD;
  int nIntcoAtoms;

  InternalCoordinate *pPointer;
  Stretch *pStre;
  Bend *pBend;
  Linear1 *pLin1;
  OutOfPlane *pOut;
  Torsion *pTors;
  Spf *pSpf;
  LinearX *pLinX;
  LinearY *pLinY;
  Rcom *pRcom;
  
  pPointer = gIntCo.vectorInternals[r];
  
  switch(pPointer->getType())
    {
    case STRE:
      pStre = (Stretch *)pPointer;
      atomA = (int)pStre->atomA;
      atomB = (int)pStre->atomB;
      Stretch::Vect((int)disp, atomA, atomB, (double*)s1, &value,1);
      for(i = 0; i < 3; i++)
	s2[i] = -(s1[i]);
      nIntcoAtoms = 2;
      break;
    case BEND:
      pBend = (Bend *)pPointer;
      atomA = (int)pBend->atomA;
      atomB = (int)pBend->atomB;
      atomC = (int)pBend->atomC;
      Bend::Vect((int)disp, atomA, atomB, atomC, (double*)s1, (double*)s2, (double*)s3, &value,1);
      nIntcoAtoms = 3;
      break;
    case LIN1:
      pLin1 = (Linear1 *)pPointer;
      atomA = (int)pLin1->atomA;
      atomB = (int)pLin1->atomB;
      atomC = (int)pLin1->atomC;
      atomD = (int)pLin1->atomD;
      Linear1::Vect((int)disp, atomA, atomB, atomC, atomD, (double*)s1, (double*)s2, (double*)s3, &value,1);
      nIntcoAtoms = 3;
      break;
    case OUT:
      pOut = (OutOfPlane *)pPointer;
      atomA = (int)pOut->atomA;
      atomB = (int)pOut->atomB;
      atomC = (int)pOut->atomC;
      atomD = (int)pOut->atomD;
      OutOfPlane::Vect((int)disp, atomA, atomB, atomC, atomD, (double*)s1, (double*)s2, (double*)s3, (double*)s4, &value,1);
      nIntcoAtoms = 4;
      break;
    case TORS:
      pTors = (Torsion *)pPointer;
      atomA = (int)pTors->atomA;
      atomB = (int)pTors->atomB;
      atomC = (int)pTors->atomC;
      atomD = (int)pTors->atomD;
      Torsion::Vect((int)disp, atomA, atomB, atomC, atomD, (double*)s1, (double*)s2, (double*)s3, (double*)s4, &value,1);
      nIntcoAtoms = 4;
      break;
    case SPF:
      break;
    case LINX:
      pLinX = (LinearX *)pPointer;
      atomA = (int)pLinX->atomA;
      atomB = (int)pLinX->atomB;
      atomC = (int)pLinX->atomC;
      atomD = (int)pLinX->atomD;
      LinearX::Vect((int)disp, atomA, atomB, atomC, atomD, (double*)s1, (double*)s2, (double*)s3, (double*)s4, &value,1);
      nIntcoAtoms = 4;
      break;
    case LINY:
      pLinY = (LinearY *)pPointer;
      atomA = (int)pLinY->atomA;
      atomB = (int)pLinY->atomB;
      atomC = (int)pLinY->atomC;
      atomD = (int)pLinY->atomD;
      LinearY::Vect((int)disp, atomA, atomB, atomC, atomD, (double*)s1, (double*)s2, (double*)s3, (double*)s4, &value,1);
      break;
      nIntcoAtoms = 4;
    }
  
  *Aptr = atomA;
  *Bptr = atomB;
  *Cptr = atomC;
  *Dptr = atomD;

  *vptr = value;

  return nIntcoAtoms;
}


void InternalCoordinates::addInternalCoordinate(InternalCoordinate* ico)
{
  vectorInternals.push_back(ico);
}

void InternalCoordinates::printInternalCoordinates()
{
  int index;

  fprintf(outfile, "List of internal coordinates found:\n");

  for (index=0; index < vectorInternals.size(); index++)
    vectorInternals[index]->printInfo();
}

void InternalCoordinates::printSingleIntCo(int i)
{
  vectorInternals[i]->printInfo();
   if(vectorInternals[i]->getType() == STRE) 
    fprintf(outfile, "\tS = %lf\n", gBMat.SVectArray[i]); 
   else
     fprintf(outfile, "\tS = %lf\n", gBMat.SVectArray[i] * 180 / _pi); 
   //Would be neat if I could figure out how to save the SVect value to intco class... need to think about this
}

int InternalCoordinates::InternalCoordinateSize()
{
  return vectorInternals.size();
}
