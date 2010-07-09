/*! \file
    \ingroup OPTKING
    \brief Class for angle bends
*/

#ifndef _psi3_bin_optking_bend_h_
#define _psi3_bin_optking_bend_h_

namespace psi { //namespace optking {

class bend_class {

  private:

    int id;
    int A;
    int B;
    int C;
    double val;
    double s_A[3]; /* The s vector for atom A */
    double s_B[3];
    double s_C[3];

  public:

    friend class simples_class;

    bend_class(void) { }

    ~bend_class() { }

    bend_class(int id_in, int A_in, int B_in, int C_in, double val_in = 0.0) {
      id = id_in;
      B = B_in;
      if (A_in < C_in) { //canonical order is A < C
        A = A_in;
        C = C_in;
      }
      else {
        A = C_in;
        C = A_in;
      }
      val = val_in;
    }

    void print(FILE *fp_out, bool print_vals) const {
      if (print_vals)
        fprintf(fp_out,"    (%d %d %d %d) (%.8lf)\n", id,A+1,B+1,C+1,val);
      else 
        fprintf(fp_out,"    (%d %d %d %d)\n", id,A+1,B+1,C+1);
    }

    void print_s(FILE *fp_out) const {
      fprintf(fp_out,"S vector for bend %d %d %d: atom A\n", A+1, B+1, C+1);
      fprintf(fp_out,"(%16.10f,%16.10f,%16.10f)\n", s_A[0], s_A[1], s_A[2]);
      fprintf(fp_out,"S vector for bend %d %d %d: atom B\n", A+1, B+1, C+1);
      fprintf(fp_out,"(%16.10f,%16.10f,%16.10f)\n", s_B[0], s_B[1], s_B[2]);
      fprintf(fp_out,"S vector for bend %d %d %d: atom C\n", A+1, B+1, C+1);
      fprintf(fp_out,"(%16.10f,%16.10f,%16.10f)\n", s_C[0], s_C[1], s_C[2]);
    }

    void set_id(int i){ id = i;}
    void set_A(int i) { A = i;}
    void set_B(int i) { B = i;}
    void set_C(int i) { C = i;}
    void set_val(double new_val) { val = new_val;}
    void set_s_A(double s_A0, double s_A1, double s_A2) {
         s_A[0] = s_A0; s_A[1] = s_A1; s_A[2] = s_A2;
    }
    void set_s_B(double s_B0, double s_B1, double s_B2) {
         s_B[0] = s_B0; s_B[1] = s_B1; s_B[2] = s_B2;
    }
    void set_s_C(double s_C0, double s_C1, double s_C2) {
         s_C[0] = s_C0; s_C[1] = s_C1; s_C[2] = s_C2;
    }

    int  get_id(void) const { return id;}
    int  get_A(void) const  { return A;}
    int  get_B(void) const  { return B;}
    int  get_C(void) const  { return C;}

    double get_val(void) const  { return val;}
    double get_val_A_or_rad(void) const  {
      return (val/180*_pi);
    }

    double get_s_A(int i) const  { return s_A[i]; }
    double get_s_B(int i) const { return s_B[i]; }
    double get_s_C(int i) const { return s_C[i]; }
    int  get_atom(int a) const  {
      if (a==0) return A;
      else if (a==1) return B;
      else if (a==2) return C;
      else throw("bend_class::get_atom : atom index must be 0, 1 or 2.\n");
    }

    double get_s(int atom, int xyz) const  {
      if ( xyz < 0 || xyz > 2) throw ("bend_class::get_s() : xyz must be 0, 1 or 2");
      if (atom==0) return s_A[xyz];
      else if (atom==1) return s_B[xyz];
      else if (atom==2) return s_C[xyz];
      else throw("bend_class::get_s() : atom index must be 0, 1, or 2");
    }


    // takes geometry in au; stores value in degrees
    void compute(double *geom) {
      int i,j;
      double rBA,rBC,eBA[3],eBC[3],tmp[3],dotprod;
    
      for (j=0;j<3;++j) {
        eBA[j] = geom[3*A+j] - geom[3*B+j];
        eBC[j] = geom[3*C+j] - geom[3*B+j];
      }
    
      rBA = sqrt( SQR(eBA[0])+SQR(eBA[1])+SQR(eBA[2]) );
      rBC = sqrt( SQR(eBC[0])+SQR(eBC[1])+SQR(eBC[2]) );
    
      scalar_div(rBA,eBA);
      scalar_div(rBC,eBC);
    
      dot_array(eBA,eBC,3,&dotprod);
    
      if (dotprod > 1.0) val = 0.0;
      else if (dotprod < -1.0) val = 180.0;
      else val = acos(dotprod) * 180.0 / _pi;
    
      return;
    }

    // takes geometry in au; stores s vectors in 1/Angstroms
    // length of s_A is equal to 1/R(AB) [See Wilson, page 56]
    // length of s_C is equal to 1/R(BC)
    // length of s_B = -s_A - s_C
    void compute_s(double *geom) {
      int j;
      double val_rad, rBA, rBC;
      double eBA[3], eBC[3], tmp[3];
    
      val_rad = val * _pi / 180.0;
    
      for (j=0;j<3;++j) {
        eBA[j] = geom[3*A+j] - geom[3*B+j];
        eBC[j] = geom[3*C+j] - geom[3*B+j];
      }
    
      rBA = sqrt( SQR(eBA[0]) + SQR(eBA[1]) + SQR(eBA[2]) );
      rBC = sqrt( SQR(eBC[0]) + SQR(eBC[1]) + SQR(eBC[2]) );

      scalar_div(rBA, eBA);
      scalar_div(rBC, eBC);

      for (j=0; j<3; ++j) {
        s_A[j] = (eBA[j]*cos(val_rad) - eBC[j]) / (rBA*sin(val_rad));

        s_C[j] = (eBC[j]*cos(val_rad) - eBA[j]) / (rBC*sin(val_rad));
      }

      for (j=0; j<3; ++j)
        s_B[j] = - s_A[j] - s_C[j];
      //  s_B[j] = ((rBA - rBC*cos(val_rad))*eBA[j] + (rBC-rBA*cos(val_rad))*eBC[j])
      //            / (rBA * rBC * sin(val_rad));

      for (j=0; j<3; ++j) {
        s_A[j] /= _bohr2angstroms;
        s_B[j] /= _bohr2angstroms;
        s_C[j] /= _bohr2angstroms;
      }

      return;
    }

    bool operator==(const bend_class & s2) const {
      if ( this->A == s2.A && this->B == s2.B && this->C == s2.C)
        return true;
      else if ( this->A == s2.C && this->B == s2.B && this->C == s2.A)
        return true;
      else
        return false;
    };

};

}//} /* namespace psi::optking */

#endif
