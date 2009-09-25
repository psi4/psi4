/*! \file
    \ingroup INTDER
    \brief Enter brief description of file here 
*/
#ifndef _psi_bin_intder_bmat_h_
#define _psi_bin_intder_bmat_h_

namespace psi { namespace intder {

class BMat
{
 
public:
  double **BMatrix;
  double **BMatSave;
  double *SVectArray;
  double **AMatrix;
  int disp;

  BMat();
  ~BMat();

  void init();

  void make_BMat();
  void invert_BMat();
  void BRow(double*, double, int, double*);
  void StoreElement(double*, int, int, double*, double*);
  void StoreElement(double*, int, int, int, double*, double*, double*);
  void StoreElement(double*, int, int, int, int, double*, double*, double*);
  void StoreElement(double*, int, int, int, int, double*, double*, double*, double*);
  void print_BMat();
};

}} // namespace psi::intder

#endif // header guard
