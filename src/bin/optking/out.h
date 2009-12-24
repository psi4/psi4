/*! \file
    \ingroup OPTKING
    \brief Class for out-of-plane angles
*/

#ifndef _psi3_bin_optking_out_h_
#define _psi3_bin_optking_out_h_

namespace psi { namespace optking {

class out_class {

  private:

    int id;
    int A;
    int B;
    int C;
    int D;
    double val;
    double s_A[3]; // The s vector for atom A
    double s_B[3];
    double s_C[3];
    double s_D[3];

  public:

    friend class simples_class;

    out_class(void) { }

    ~out_class() { }
 
    out_class(int id_in, int A_in, int B_in, int C_in, int D_in, double val_in = 0.0) {
      id = id_in;
      A = A_in;
      B = B_in;
      if (C_in < D_in) { // canonical is C < D
        C = C_in;
        D = D_in;
        val = val_in;
      } 
      else {
        C = D_in;
        D = C_in;
        val = -1.0 * val_in;
      }
    }

    void print(FILE *fp_out, bool print_vals) const {
      if (print_vals) 
        fprintf(fp_out,"    (%d %d %d %d %d) (%.8lf)\n", id,A+1,B+1,C+1,D+1,val);
      else
        fprintf(fp_out,"    (%d %d %d %d %d)\n", id,A+1,B+1,C+1,D+1);
    }

    void print_s(FILE *fp_out) const {
      fprintf(fp_out,"S vector for out %d %d %d %d:atom A\n", A+1, B+1, C+1, D+1);
      fprintf(fp_out,"(%16.10f,%16.10f,%16.10f)\n", s_A[0], s_A[1], s_A[2]);
      fprintf(fp_out,"S vector for out %d %d %d %d:atom B\n", A+1, B+1, C+1, D+1);
      fprintf(fp_out,"(%16.10f,%16.10f,%16.10f)\n", s_B[0], s_B[1], s_B[2]);
      fprintf(fp_out,"S vector for out %d %d %d %d:atom C\n", A+1, B+1, C+1, D+1);
      fprintf(fp_out,"(%16.10f,%16.10f,%16.10f)\n", s_C[0], s_C[1], s_C[2]);
      fprintf(fp_out,"S vector for out %d %d %d %d:atom D\n", A+1, B+1, C+1, D+1);
      fprintf(fp_out,"(%16.10f,%16.10f,%16.10f)\n", s_D[0], s_D[1], s_D[2]);
    }

    void set_id(int i){ id = i;}
    void set_A(int i) { A = i;}
    void set_B(int i) { B = i;}
    void set_C(int i) { C = i;}
    void set_D(int i) { D = i;}
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
    void set_s_D(double s_D0, double s_D1, double s_D2) {
         s_D[0] = s_D0; s_D[1] = s_D1; s_D[2] = s_D2;
    }

    int  get_id(void) const{ return id;}
    int  get_A(void) const { return A;}
    int  get_B(void) const { return B;}
    int  get_C(void) const { return C;}
    int  get_D(void) const { return D;}
    double get_val(void) const { return val;}
    double get_val_A_or_rad(void) const  { 
      return (val/180*_pi);
    }

    double get_s_A(int i) const { return s_A[i]; }
    double get_s_B(int i) const { return s_B[i]; }
    double get_s_C(int i) const { return s_C[i]; }
    double get_s_D(int i) const { return s_D[i]; }

    int  get_atom(int a) const  {
      if (a==0) return A;
      else if (a==1) return B;
      else if (a==2) return C;
      else if (a==3) return D;
      else throw("out_class::get_atom : atom index must be 0, 1, 2 or 3.\n");
    }

    double get_s(int atom, int xyz) const  {
      if ( xyz < 0 || xyz > 2) throw ("out_class::get_s() : xyz must be 0, 1 or 2");
      if (atom==0) return s_A[xyz];
      else if (atom==1) return s_B[xyz];
      else if (atom==2) return s_C[xyz];
      else if (atom==3) return s_D[xyz];
      else throw("out_class::get_s() : atom index must be 0, 1, 2, or 3");
    }

    void compute(double *geom) {
      int j;
      double rBA, rBC, rBD, phi_CBD = 0.0, dotprod = 0.0, angle;
      double eBA[3], eBC[3], eBD[3], tmp[3];
    
      for (j=0;j<3;++j) {
        eBA[j] = geom[3*A+j] - geom[3*B+j];
        eBC[j] = geom[3*C+j] - geom[3*B+j];
        eBD[j] = geom[3*D+j] - geom[3*B+j];
      }
    
      rBA = sqrt( SQR(eBA[0]) + SQR(eBA[1]) + SQR(eBA[2]) );
      rBC = sqrt( SQR(eBC[0]) + SQR(eBC[1]) + SQR(eBC[2]) );
      rBD = sqrt( SQR(eBD[0]) + SQR(eBD[1]) + SQR(eBD[2]) );
    
      scalar_div(rBA,eBA);
      scalar_div(rBC,eBC);
      scalar_div(rBD,eBD);
    
      dot_array(eBC,eBD,3,&phi_CBD);
    
      if (phi_CBD > 1.0) phi_CBD = 0.0;
      else if (phi_CBD < -1.0) phi_CBD = _pi ;
      else phi_CBD = acos(phi_CBD) ;
    
      cross_product(eBC,eBD,tmp);
    
      dot_array(tmp,eBA,3,&dotprod);
    
      if (sin(phi_CBD) > optinfo.sin_phi_denominator_tol) dotprod = dotprod / sin(phi_CBD) ;
      else dotprod = 0.0 ;
    
      if (dotprod > 1.0) angle = _pi / 2.0;
      else if (dotprod < -1.0) angle = -1.0 * _pi / 2.0000;
      else angle = asin(dotprod) ;
    
      val = angle * 180.0 / _pi;
      return;
    }

    // takes geometry in au; stores s vectors in 1/Angstroms
    void compute_s(double *geom) {
      int j;
      double rBA, rBC, rBD, phi_CBD = 0.0, val_rad;
      double eBA[3], eBC[3], eBD[3], tmp[3], tmp2[3], tmp3[3], temp;
    
      val_rad = val * _pi / 180.0;
    
      for (j=0;j<3;++j) {
        eBA[j] = geom[3*A+j] - geom[3*B+j];
        eBC[j] = geom[3*C+j] - geom[3*B+j];
        eBD[j] = geom[3*D+j] - geom[3*B+j];
      }
    
      rBA = sqrt( SQR(eBA[0]) + SQR(eBA[1]) + SQR(eBA[2]) );
      rBC = sqrt( SQR(eBC[0]) + SQR(eBC[1]) + SQR(eBC[2]) );
      rBD = sqrt( SQR(eBD[0]) + SQR(eBD[1]) + SQR(eBD[2]) );
    
      scalar_div(rBA,eBA);
      scalar_div(rBC,eBC);
      scalar_div(rBD,eBD);
    
      dot_array(eBC,eBD,3,&phi_CBD);
    
      if (phi_CBD > 1.0) phi_CBD = 0.0;
      else if (phi_CBD < -1.0) phi_CBD = _pi;
      else phi_CBD = acos(phi_CBD);
    
      cross_product(eBC,eBD,tmp);
      scalar_div(cos(val_rad)*sin(phi_CBD),tmp);
      for (j=0;j<3;++j) 
         tmp2[j] = tan(val_rad) * eBA[j];
      for (j=0;j<3;++j) 
         tmp3[j] = (tmp[j] - tmp2[j])/rBA;
      s_A[0] = tmp3[0];
      s_A[1] = tmp3[1];
      s_A[2] = tmp3[2];
    
      cross_product(eBD,eBA,tmp);
      scalar_div(cos(val_rad)*sin(phi_CBD),tmp);
      for (j=0;j<3;++j)
        tmp2[j] = cos(phi_CBD) * eBD[j];
      for (j=0;j<3;++j)
        tmp3[j] = eBC[j] - tmp2[j];
      scalar_mult(tan(val_rad)/SQR(sin(phi_CBD)),tmp3,3);
      for (j=0;j<3;++j)
         tmp2[j] = (tmp[j] - tmp3[j])/rBC;
      s_C[0] = tmp2[0];
      s_C[1] = tmp2[1];
      s_C[2] = tmp2[2];
    
      cross_product(eBA,eBC,tmp);
      scalar_div(cos(val_rad)*sin(phi_CBD),tmp);
      for (j=0;j<3;++j)
        tmp2[j] = cos(phi_CBD) * eBC[j];
      for (j=0;j<3;++j)
        tmp3[j] = eBD[j] - tmp2[j];
      scalar_mult(tan(val_rad)/SQR(sin(phi_CBD)),tmp3,3);

      for (j=0;j<3;++j)
        s_D[j] = (tmp[j] - tmp3[j])/rBD;

      for (j=0;j<3;++j)
        s_B[j]  = (-1.0) * s_A[j] - s_C[j] - s_D[j];

      for (j=0;j<3;++j) {
        s_A[j] /= _bohr2angstroms;
        s_B[j] /= _bohr2angstroms;
        s_C[j] /= _bohr2angstroms;
        s_D[j] /= _bohr2angstroms;
      }

      return;
    }

    bool operator==(const out_class & s2) const {
      if ( this->A == s2.A && this->B == s2.B && this->C == s2.C && this->D == s2.D)
        return true;
      else
        return false;
    };

};

}} /* namespace psi::optking */

#endif
