/*! \file
    \ingroup INTDER
    \brief Enter brief description of file here 
*/
#ifndef _psi_bin_intder_transform_h_
#define _psi_bin_intder_transform_h_

#include <vector>
#include "3dmatrix.h"

#ifdef EXTERN
#undef EXTERN
#define EXTERN extern
#else
#define EXTERN
#endif

namespace psi { namespace intder {

class Transform
{

public:

  void der_transform(int, int, double **);

/* Numerical testing stuff
  double *BMatrow;
  C3DMatrix *tempYRow;

  double **SR2Matrix;  //Numerical 2nd derivative test matrix
  C3DMatrix *SR3Matrix;
  C4DMatrix *SR4Matrix;

  double **H11;
  double **H21;
  double **H22;
  double **H31;
  double **H32;
  double **H33;
  double **H41;
  double **H42;
  double **H43;
  double **H44;
  
  C3DMatrix *H111;
  C3DMatrix *H112;
  C3DMatrix *H221;
  C3DMatrix *H222;
  C3DMatrix *H113;
  C3DMatrix *H123;
  C3DMatrix *H223;
  C3DMatrix *H331;
  C3DMatrix *H332;
  C3DMatrix *H333;
  C3DMatrix *H411;
  C3DMatrix *H421;
  C3DMatrix *H422;
  C3DMatrix *H431;
  C3DMatrix *H432;
  C3DMatrix *H433;
  C3DMatrix *H441;
  C3DMatrix *H442;
  C3DMatrix *H443;
  C3DMatrix *H444;
  */

  Transform();
  ~Transform();

  /* Numerical testing stuff
  void ThirdDerinit();
  void FourthDerinit();

  void BRow(double *, double, int, double *);
  void XRow(double **, double, int);
  void YRow(C3DMatrix *, double, int);
  void StoreElement(double*, int, int, double*, double*);
  void StoreElement(double*, int, int, int, double*, double*, double*);
  void StoreElement(double*, int, int, int, int, double*, double*, double*);
  void StoreElement(double*, int, int, int, int, double*, double*, double*, double*);
  */
  void open_PSIF();
  void close_PSIF();
};

/*
class SecondDerivative : public Transform
{
public:

  SecondDerivative();
  ~SecondDerivative();
    
  void i_to_c();

  void EvaluateH2(int, int, double **);
  void EvaluateH3(int, int, int, double **, double **, double **, double **, double **, double **);
  void EvaluateH4(int, int, int, int, double **, double **, double **, double **, double **, double **, double **, double **, double **, double **);

  void NumTest();
  void SRTest1();

  void secondDerivativeOut(int, double **);
  void secondDerivativeIn(int, double **);

private:

  static const double graddisp = 1.0E-4;

};


class ThirdDerivative : public Transform
{
protected:


public:

  ThirdDerivative();
  ~ThirdDerivative();

  void i_to_c();

  void EvaluateH2(int, int, C3DMatrix*, C3DMatrix*);
  void EvaluateH3(int, int, int, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*);
  void EvaluateH4(int, int, int, int, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*);
  
  void FillSR2(int, int, C3DMatrix*);
  void FillSR3(int, int, int, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*);
  void FillSR4(int, int, int, int, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*, C3DMatrix*);
  
  void NumTest();
  void SRTest2();
  
  void thirdDerivativeOut(int, C3DMatrix*);
  void thirdDerivativeIn(int, C3DMatrix*);

 private:

  static const double graddisp = 1.0E-4;

};

class FourthDerivative : public Transform
{
protected:


public:
  FourthDerivative();
  ~FourthDerivative();

  void i_to_c();

  void NumTest();

  void thirdDerivativeBlockOut(int, int, C3DMatrix*);
  void thirdDerivativeBlockIn(int, int, C3DMatrix*);

 private:

  static const double graddisp = 1.0E-4;

};
*/

}} // namespace psi::intder

#endif // header guard
