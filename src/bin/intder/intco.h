/*! \file
    \ingroup INTDER
    \brief Enter brief description of file here 
*/
#ifndef _psi_bin_intder_intco_h_
#define _psi_bin_intder_intco_h_

#include <vector>
#include "3dmatrix.h"

namespace psi { namespace intder {

enum IntCoType { STRE, BEND, LIN1, OUT, TORS, SPF, LINX, LINY, RCOM, SYMM };

class InternalCoordinate 
{
 protected:
  IntCoType type;
  int id;
  double value;

 public:
  IntCoType getType()
    { return type; }
  
  void setID(int i)
    { id = i; }
  int getID()
    { return id; }
  void setValue(double v) 
    { value = v;}
  double getValue()
   { return value; }
  
  virtual void printInfo(); 
};

class Stretch : public InternalCoordinate
{
public:
  int atomA, atomB;

  Stretch(int myid, int aa, int ab)
    { id = myid; atomA = aa; atomB = ab; type = STRE; }
  Stretch()
    { id = -1; atomA = -1; atomB = -1; type = STRE; }

  void printInfo();

  // Conversion of Dr. Allen's VECT1 subroutine from Fortran
  static void Vect(int disp, int k1, int k2, double *v1, double *dist, int pflag = 0);
	//HIJS1 from intder2000.f
  static void Hijs(int disp, int k1, int k2, double **h11);
	//HIJKS1 from intder2000.f
	static void Hijks(int disp, int k1, int k2, C3DMatrix *h111);
  //H4TH1 from intder2000.f
  static void H4th(int disp, int k1, int k2, C4DMatrix *h1111);
	//H5TH1 from intder2000.f
	static void H5th(int disp, int k1, int k2, C5DMatrix *h11111);

};

class Bend : public InternalCoordinate
{
public:
  int atomA, atomB, atomC;

  Bend(int myid, int aa, int ab, int ac)
    { id = myid; atomA = aa; atomB = ab; atomC = ac; type = BEND; }
  Bend()
    { id = -1; atomA = -1; atomB = -1; atomC = -1; type = BEND; }

  void printInfo();

  static void Vect(int disp, int k1, int k2, int k3, double *s1, double *s2, double *s3, double *theta, int pflag = 0);
  static void Hijs(int disp, int k1, int k2, int k3, 
                   double **h11, double **h21, double **h31, 
                   double **h22, double **h32, double **h33);
	//HIJKS2
	static void Hijks(int disp, int k1, int k2, int k3, C3DMatrix *h111, C3DMatrix *h112,
                    C3DMatrix *h113, C3DMatrix *h123, C3DMatrix *h221, C3DMatrix *h222,
                    C3DMatrix *h223, C3DMatrix *h331, C3DMatrix *h332, C3DMatrix *h333);
  static void H4th(int disp, int k1, int k2, int k3, C4DMatrix *h1111, C4DMatrix *h1112, C4DMatrix *h1113,
                    C4DMatrix *h1122, C4DMatrix *h1123, C4DMatrix *h1133, C4DMatrix *h1222,
                    C4DMatrix *h1223, C4DMatrix *h1233, C4DMatrix *h1333, C4DMatrix *h2222,
                    C4DMatrix *h2223, C4DMatrix *h2233, C4DMatrix *h2333, C4DMatrix *h3333);
};

class Linear1 : public InternalCoordinate
{
public:
  int atomA, atomB, atomC, atomD;

  Linear1(int myid, int aa, int ab, int ac, int ad)
    { id = myid; atomA = aa; atomB = ab; atomC = ac; atomD = ad; type = LIN1; }
  Linear1()
    { id = -1; atomA = -1; atomB = -1; atomC = -1; atomD = -1; type = LIN1; }

  void printInfo();

  // Conversion of Dr. Allen's VECT4 subroutine.
  // NOTE: The LIN2 of INTDER2004 is now considered LIN1
  //       The LIN1 of INTDER2004 is no longer supported.
  static void Vect(int disp, int k1, int k2, int k3, int k4, double *s1, double *s2, double *s3, 
		   double *theta, int pflag = 0);
  static void Hijs(int disp, int k1, int k2, int k3, int k4, 
                   double **h11, double **h21, double **h31, 
                   double **h22, double **h32, double **h33);
	//HIJKS from intder2000.f
  static void Hijks(int disp, int k1, int k2, int k3, int k4,
                    C3DMatrix *h111, C3DMatrix *h112, C3DMatrix *h113, C3DMatrix *h123,
                    C3DMatrix *h221, C3DMatrix *h222, C3DMatrix *h223, C3DMatrix *h331,
                    C3DMatrix *h332, C3DMatrix *h333);

};

class OutOfPlane : public InternalCoordinate
{
public:
  int atomA, atomB, atomC, atomD;

  OutOfPlane(int myid, int aa, int ab, int ac, int ad)
    { id = myid; atomA = aa; atomB = ab; atomC = ac; atomD = ad; type = OUT; }
  OutOfPlane()
    { id = -1; atomA = -1; atomB = -1; atomC = -1; atomD = -1; type = OUT; }

  void printInfo();
  static void Vect(int disp, int k1, int k2, int k3, int k4, double *s1, double *s2, double *s3, double *s4,
		   double *theta,int pflag = 0);
  static void Hijs(int disp, int k1, int k2, int k3, int k4,
                   double **h11, double **h21, double **h31, double **h41,
                   double **h22, double **h32, double **h42, double **h33,
                   double **h43, double **h44);
  static void Hijks(int disp, int k1, int k2, int k3, int k4, C3DMatrix *h111, C3DMatrix *h112,
                    C3DMatrix *h221, C3DMatrix *h222, C3DMatrix *h113, C3DMatrix *h123,
                    C3DMatrix *h223, C3DMatrix *h331, C3DMatrix *h332, C3DMatrix *h333,
                    C3DMatrix *h411, C3DMatrix *h421, C3DMatrix *h422, C3DMatrix *h431,
                    C3DMatrix *h432, C3DMatrix *h433, C3DMatrix *h441, C3DMatrix *h442,
                    C3DMatrix *h443, C3DMatrix *h444);
};

class Torsion : public InternalCoordinate
{
public:
  int atomA, atomB, atomC, atomD;
  Torsion(int myid, int aa, int ab, int ac, int ad)
    { id = myid; atomA = aa; atomB = ab; atomC = ac; atomD = ad; type = TORS; }
  Torsion()
    { id = -1; atomA = -1; atomB = -1; atomC = -1; atomD = -1; type = TORS; }

  void printInfo();

  static void Vect(int disp, int k1, int k2, int k3, int k4, double *s1, double *s2, double *s3, 
		   double *s4, double *theta, int pflag = 0);
	
	static void Hijs(int disp, int k1, int k2, int k3, int k4, double **h11, 
                   double **h21, double **h31, double **h41, double **h22, 
                   double **h32, double **h42, double **h33, double **h43, 
                   double **h44);
	static void Hijks(int disp, int k1, int k2, int k3, int k4, C3DMatrix *h111, C3DMatrix *h112,
                    C3DMatrix *h221, C3DMatrix *h222, C3DMatrix *h113, C3DMatrix *h123,
                    C3DMatrix *h223, C3DMatrix *h331, C3DMatrix *h332, C3DMatrix *h333,
                    C3DMatrix *h411, C3DMatrix *h421, C3DMatrix *h422, C3DMatrix *h431,
                    C3DMatrix *h432, C3DMatrix *h433, C3DMatrix *h441, C3DMatrix *h442,
                    C3DMatrix *h443, C3DMatrix *h444);
};

class Spf : public InternalCoordinate
{
public:
  int atomA, atomB;

  Spf(int myid, int aa, int ab)
    { id = myid; atomA = aa; atomB = ab; type = SPF; }
  Spf()
    { id = -1; atomA = -1; atomB = -1; type = SPF; }

  void printInfo();
  
};

class LinearX : public InternalCoordinate
{
public:
  int atomA, atomB, atomC, atomD;

  LinearX(int myid, int aa, int ab, int ac, int ad)
    { id = myid; atomA = aa; atomB = ab; atomC = ac; atomD = ad; type = LINX; }
  LinearX()
    { id = -1; atomA = -1; atomB = -1; atomC = -1; atomD = -1; type = LINX; }

  void printInfo();
	//VECT8
  static void Vect(int disp, int k1, int k2, int k3, int k4, double *s1, double *s2, double *s3, 
		   double *s4, double *theta,int pflag = 0);
//HIJS8 from intder2000.f
	static void Hijs(int disp, int k1, int k2, int k3, int k4, double **h11, double **h21,
                 double **h31, double **h41, double **h22, double **h32, double **h42,
                 double **h33, double **h43, double **h44);
//HIJKS8 from intder2000.f
        static void Hijks(int disp, int k1, int k2, int k3, int k4, C3DMatrix *h111, C3DMatrix *h112,
                    C3DMatrix *h221, C3DMatrix *h222, C3DMatrix *h113, C3DMatrix *h123,
                    C3DMatrix *h223, C3DMatrix *h331, C3DMatrix *h332, C3DMatrix *h333,
                    C3DMatrix *h411, C3DMatrix *h421, C3DMatrix *h422, C3DMatrix *h431,
                    C3DMatrix *h432, C3DMatrix *h433, C3DMatrix *h441, C3DMatrix *h442,
                    C3DMatrix *h443, C3DMatrix *h444);
};

class LinearY : public InternalCoordinate
{
public:
  int atomA, atomB, atomC, atomD;

  LinearY(int myid, int aa, int ab, int ac, int ad)
    { id = myid; atomA = aa; atomB = ab; atomC = ac; atomD = ad; type = LINY; }
  LinearY()
    { id = -1; atomA = -1; atomB = -1; atomC = -1; atomD = -1; type = LINY; }

  void printInfo();

//VECT9
  static void Vect(int disp, int k1, int k2, int k3, int k4, double *s1, double *s2, double *s3, 
		   double *s4, double *theta,int pflag = 0);
//HIJS9 from intder2000.f
	static void Hijs(int disp, int k1, int k2, int k3, int k4, double **h11, double **h21,
                   double **h31, double **h41, double **h22, double **h32, double **h42,
                   double **h33, double **h43, double **h44);
        static void Hijks(int disp, int k1, int k2, int k3, int k4, C3DMatrix *h111, C3DMatrix *h112,
                    C3DMatrix *h221, C3DMatrix *h222, C3DMatrix *h113, C3DMatrix *h123,
                    C3DMatrix *h223, C3DMatrix *h331, C3DMatrix *h332, C3DMatrix *h333,
                    C3DMatrix *h411, C3DMatrix *h421, C3DMatrix *h422, C3DMatrix *h431,
                    C3DMatrix *h432, C3DMatrix *h433, C3DMatrix *h441, C3DMatrix *h442,
                    C3DMatrix *h443, C3DMatrix *h444);
};

class Rcom : public InternalCoordinate
{
public:
  int atomA, atomB;

  Rcom(int myid, int aa, int ab)
    { id = myid; atomA = aa; atomB = ab; type = RCOM; }
  Rcom()
    { id = -1; atomA = -1; atomB = -1; type = RCOM; }

  void printInfo();
};

class InternalCoordinates
{
public:

  std::vector<InternalCoordinate*> vectorInternals;
  InternalCoordinates()
    { };
  ~InternalCoordinates();

  void addInternalCoordinate(InternalCoordinate* ico);
  void printInternalCoordinates();

  void loadInternalCoordinates();
  int InternalCoordinateSize();
  void printSingleIntCo(int);
  int intcoSwitch(double, int, int*, int*, int*, int*, double*, double*, double*, double*, double*);
};

}} // namespace psi::intder

#endif // header guard

