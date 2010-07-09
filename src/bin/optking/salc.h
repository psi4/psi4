/*! \file
    \ingroup OPTKING
    \brief Class for salcs
*/

#ifndef _psi3_bin_optking_salc_h_
#define _psi3_bin_optking_salc_h_

#define EXTERN
#include "globals.h"
#undef EXTERN
#include "cartesians.h"

#include <libipv1/ip_lib.h>

namespace psi { //namespace optking {

class salc_class {
    int length;
    double prefactor;
    int *simple;
    double *coeff;
    char *label;
  public:
    salc_class() { 
      simple = new int[MAX_SALC_LENGTH];
      coeff = new double[MAX_SALC_LENGTH];
      label = new char[MAX_LINELENGTH];
    }
    ~salc_class() {
      // fprintf(stdout, "destructing salc class\n");
      delete [] simple;
      delete [] coeff;
      delete [] label;
    }
    void print() const {
      int i, col = 0;
        fprintf(outfile,"    (");
        fprintf(outfile,"\"%s\"",label);
        for (col=0, i=0;i<length;++i, ++col) {
          if (i == 0)
             fprintf(outfile," (%d",simple[i]);
          else
             fprintf(outfile," %d", simple[i]);
          if (col == 15) { fprintf(outfile,"\n     "); col=0; }
        }
        fprintf(outfile,")\n    ");
        for (col=0, i=0;i<length;++i, ++col) {
          if (i == 0)
             fprintf(outfile," (%5.2f", coeff[i]);
          else
             fprintf(outfile," %5.2f", coeff[i]);
          if (col == 9) { fprintf(outfile,"\n     "); col=0; }
        }
        fprintf(outfile,"))\n");
    }
    void set_length(int new_length) {
      if (new_length < MAX_SALC_LENGTH)
        length = new_length;
      else {
        fprintf(outfile,"SALC length exceeds MAX_SALC_LENGTH\n");
        exit(2);
      }
    }
    int get_length() const { return length; }
    void set_simple(int pos, int i) {
      if (pos < length)
         simple[pos] = i;
      else {
         fprintf(outfile,"Position in SALC array exceed SALC length\n");
         exit(2);
      }
    }
    int  get_simple(int pos) const {
      if (pos >= length) {
         fprintf(outfile,"Position in SALC array exceeds SALC length\n");
         exit(2);
      }
      return simple[pos];
    }
    void set_coeff(int pos, double new_coeff) {
      if (pos >= length) {
         fprintf(outfile,"Position in SALC array exceeds SALC length\n");
         exit(2);
      }
      coeff[pos] = new_coeff;
    }
    double get_coeff(int pos) const {
      if (pos >= length) {
         fprintf(outfile,"Position in SALC array exceeds SALC length\n");
         exit(2);
      }
      return coeff[pos];
    }
    char *get_label() const { return label; }
    void set_prefactor(double new_prefactor) { prefactor = new_prefactor; }
    double get_prefactor() const { return prefactor; }
    void set_label(char *new_label);
};

class salc_set {
   int num;
   salc_class *salc_array;
   char *name;
  public:
    salc_set(const char *keyword);
    salc_set();
    ~salc_set() {
      // fprintf(stdout,"destructing salc_set\n");
      delete [] salc_array;
      delete [] name ;
    }
    void print() const {
      int i;
      if (num > 0) {
        fprintf(outfile,"\n  %s = (\n",name);
        for (i=0; i < num; ++i) {
    // fprintf(outfile,"prefactor: %lf\n", salc_array[i].get_prefactor());
           salc_array[i].print();
        }
        fprintf(outfile,"  )\n");
      }
      return;
    }
    void set_num(int i) {num = i;}
    int  get_num(void) const { return num;}
    void set_coeff(int index, int pos, double new_coeff) {
      salc_array[index].set_coeff(pos, new_coeff);
    }
    double get_coeff(int index, int pos) const {
      return salc_array[index].get_coeff(pos);
    }
    int get_length(int index) const {
      return salc_array[index].get_length();
    }
    int get_simple(int index, int pos) const {
      return salc_array[index].get_simple(pos);
    }
    double get_prefactor(int index) const {
      return salc_array[index].get_prefactor();
    }
    char *get_label(int index) const {
      return salc_array[index].get_label();
    }
};

}//} /* namespace psi::optking */

#endif
