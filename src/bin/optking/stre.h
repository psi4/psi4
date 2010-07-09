/*! \file
    \ingroup OPTKING
    \brief Class for stretches
*/

#ifndef _psi3_bin_optking_stre_h_
#define _psi3_bin_optking_stre_h_

namespace psi { //namespace optking {

class stre_class {

  private:

    int id;
    int A;
    int B;
    double val; /* length of bond */
    double s_A[3];  /* The s vector for atom A (Xa-Xb) */
    double s_B[3];  /* The s vector for atom B (Xb-Xa) */

  public:

    friend class simples_class;

    stre_class(void){ }

    ~stre_class() { }

    stre_class(int id_in, int A_in, int B_in, double val_in = 0.0) {
      id = id_in;
      if (A_in < B_in) { // canonical is A < B
        A = A_in;
        B = B_in;
      }
      else {
        A = B_in;
        B = A_in;
      }
      val = val_in;
    }

    void print(FILE *fp_out, bool print_vals) const {
      if (print_vals)
        fprintf(fp_out,"    (%d %d %d) (%.8lf)\n", id, A+1, B+1, val);
      else
        fprintf(fp_out,"    (%d %d %d)\n", id, A+1, B+1);
    }

    void print_s(FILE *fp_out) const {
      fprintf(fp_out,"S vector for stretch %d %d: atom A\n", A+1, B+1);
      fprintf(fp_out,"(%16.10f,%16.10f,%16.10f)\n", s_A[0], s_A[1], s_A[2]);
      fprintf(fp_out,"S vector for stretch %d %d: atom B\n", A+1, B+1);
      fprintf(fp_out,"(%16.10f,%16.10f,%16.10f)\n", s_B[0], s_B[1], s_B[2]);
    }

    void set_id(int i) { id = i;}
    void set_A(int i)  { A = i;}
    void set_B(int i)  { B = i;}
    void set_val(double length) { val = length;}
    void set_s_A(double s_A0, double s_A1, double s_A2) {
      s_A[0] = s_A0; s_A[1] = s_A1; s_A[2] = s_A2;
    }
    void set_s_B(double s_B0, double s_B1, double s_B2) {
      s_B[0] = s_B0; s_B[1] = s_B1; s_B[2] = s_B2;
    }

    int  get_id(void) const { return id;}
    int  get_A(void) const  { return A;}
    int  get_B(void) const  { return B;}

    double get_val(void) const { return val;}
    double get_val_A_or_rad(void) const  { return val; }

    double get_s_A(int i) const { return s_A[i]; }
    double get_s_B(int i) const { return s_B[i]; }

    int  get_atom(int a) const  {
      if (a==0) return A;
      else if (a==1) return B;
      else throw("stre_class::get_atom() : atom index must be 0 or 1.\n");
    }

    double get_s(int atom, int xyz) const  {
      if ( xyz < 0 || xyz > 2) throw ("stre_class::get_s() : xyz must be 0, 1 or 2");
      if (atom==0) return s_A[xyz];
      else if (atom==1) return s_B[xyz];
      else throw("stre_class::get_s() : atom index must be 0 or 1.\n");
    }

    // takes geometry in au; stores bond length in Angstroms
    void compute(double *geom) {
      double tmp;
      tmp = SQR( geom[3*A+0] - geom[3*B+0] ) +
            SQR( geom[3*A+1] - geom[3*B+1] ) +
            SQR( geom[3*A+2] - geom[3*B+2] );
      val = sqrt(tmp)*_bohr2angstroms;
    }

    // takes geometry in au; stores s vectors which are same as unit e vectors
    void compute_s(double *geom) {
      int j;
      double eBA[3], norm;

      for (j=0;j<3;++j)
        eBA[j] = geom[3*A+j] - geom[3*B+j];

      norm = sqrt( SQR(eBA[0]) + SQR(eBA[1]) + SQR(eBA[2]) );
      scalar_div(norm,eBA);

      s_A[0] = eBA[0];
      s_A[1] = eBA[1];
      s_A[2] = eBA[2];
      s_B[0] = -1.0*eBA[0];
      s_B[1] = -1.0*eBA[1];
      s_B[2] = -1.0*eBA[2];
    }

    bool operator==(const stre_class & s2) const {
      if ( this->A == s2.A && this->B == s2.B)
        return true;
      else if ( this->A == s2.B && this->B == s2.A)
        return true;
      else
        return false;
    }

};

}//} /* namespace psi::optking */

#endif
