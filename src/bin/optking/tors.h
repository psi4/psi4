/*! \file
    \ingroup OPTKING
    \brief Class for torsions
*/

#ifndef _psi3_bin_optking_tors_h_
#define _psi3_bin_optking_tors_h_

namespace psi { //namespace optking {

class tors_class {

  private:

    int id;
    int A;
    int B;
    int C;
    int D;
    double val;
    double s_A[3];
    double s_B[3];
    double s_C[3];
    double s_D[3];
    int near_180;
        // +1 if approaching 180
        //  0 if OK
        // -1 if approaching -180

  public:

    friend class simples_class;

    tors_class(void){
      near_180 = 0;
    }

    ~tors_class(){ }

    tors_class(int id_in, int A_in, int B_in, int C_in, int D_in, double val_in = 0.0) {
      id = id_in;
      if (A_in < D_in) { // canonical order is A < D
        A = A_in;
        B = B_in;
        C = C_in;
        D = D_in;
        val = val_in;
      }
      else {
        A = D_in;
        B = C_in;
        C = B_in;
        D = A_in;
        val = -1.0 * val_in;
      }
      near_180 = 0;
    }

    void print(FILE *fp_out, bool print_vals) const {
      if (print_vals)
        fprintf(fp_out,"    (%d %d %d %d %d) (%.8lf)\n", id, A+1, B+1, C+1, D+1, val);
      else
        fprintf(fp_out,"    (%d %d %d %d %d)\n", id, A+1, B+1, C+1, D+1);
    }
 
    void print_s(FILE *fp_out) const {
      fprintf(fp_out,"S vector for tors %d %d %d %d:atom A\n", A+1, B+1, C+1, D+1);
      fprintf(fp_out,"(%16.10f,%16.10f,%16.10f)\n", s_A[0], s_A[1], s_A[2]);
      fprintf(fp_out,"S vector for tors %d %d %d %d:atom B\n", A+1, B+1, C+1, D+1);
      fprintf(fp_out,"(%16.10f,%16.10f,%16.10f)\n", s_B[0], s_B[1], s_B[2]);
      fprintf(fp_out,"S vector for tors %d %d %d %d:atom C\n", A+1, B+1, C+1, D+1);
      fprintf(fp_out,"(%16.10f,%16.10f,%16.10f)\n", s_C[0], s_C[1], s_C[2]);
      fprintf(fp_out,"S vector for tors %d %d %d %d:atom D\n", A+1, B+1, C+1, D+1);
      fprintf(fp_out,"(%16.10f,%16.10f,%16.10f)\n", s_D[0], s_D[1], s_D[2]);
    }

    void set_id(int i){ id = i;}
    void set_A(int i) { A = i;}
    void set_B(int i) { B = i;}
    void set_C(int i) { C = i;}
    void set_D(int i) { D = i;}
    void set_val(double new_val) { val = new_val;}
    void set_near_180(int new_near_180) { near_180 = new_near_180;}
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

    int get_near_180(void) const { return near_180;}
    double get_s_A(int i) const { return s_A[i]; }
    double get_s_B(int i) const { return s_B[i]; }
    double get_s_C(int i) const { return s_C[i]; }
    double get_s_D(int i) const { return s_D[i]; }

    int  get_atom(int a) const  {
      if (a==0) return A;
      else if (a==1) return B;
      else if (a==2) return C;
      else if (a==3) return D;
      else throw("tors_class::get_atom : atom index must be 0, 1, 2 or 3.\n");
    }

    double get_s(int atom, int xyz) const  {
      if ( xyz < 0 || xyz > 2) throw ("tors_class::get_s() : xyz must be 0, 1 or 2");
      if (atom==0) return s_A[xyz];
      else if (atom==1) return s_B[xyz];
      else if (atom==2) return s_C[xyz];
      else if (atom==3) return s_D[xyz];
      else throw("tors_class::get_s() : atom index must be 0, 1, 2, or 3");
    }


    // take geometry in au; store angle in degrees
    void compute(double *geom) {
      int j,k;
      double rAB, rBC, rCD, phi_123, phi_234, dotprod;
      double eAB[3], eBC[3], eCD[3], tmp[3], tmp2[3], tmp3[3];
    
      for (j=0;j<3;++j) {
        eAB[j] = geom[3*B+j] - geom[3*A+j];
        eBC[j] = geom[3*C+j] - geom[3*B+j];
        eCD[j] = geom[3*D+j] - geom[3*C+j];
      }
    
      rAB = sqrt( SQR(eAB[0]) + SQR(eAB[1]) + SQR(eAB[2]) );
      rBC = sqrt( SQR(eBC[0]) + SQR(eBC[1]) + SQR(eBC[2]) );
      rCD = sqrt( SQR(eCD[0]) + SQR(eCD[1]) + SQR(eCD[2]) );
    
      scalar_div(rAB,eAB);
      scalar_div(rBC,eBC);
      scalar_div(rCD,eCD);
    
      phi_123 = phi_234 = 0.0;
      for (j=0;j<3;++j) {
         phi_123 += (-1.0 * eAB[j]) * eBC[j];
         phi_234 += (-1.0 * eBC[j]) * eCD[j];
      }
    
      if (phi_123 > 1.0)
        phi_123 = 0.0;
      else if (phi_123 < -1.0)
        phi_123 = _pi;
      else phi_123 = acos(phi_123);
    
      if (phi_234 > 1.0)
        phi_234 = 0.0;
      else if (phi_234 < -1.0)
        phi_234 = _pi;
      else phi_234 = acos(phi_234);
    
      cross_product(eAB,eBC,tmp);
      cross_product(eBC,eCD,tmp2);
      dot_array(tmp,tmp2,3,&dotprod);

      if ((sin(phi_123) > optinfo.sin_phi_denominator_tol) &&
          (sin(phi_234) > optinfo.sin_phi_denominator_tol)) {
         dotprod /= sin(phi_123);
         dotprod /= sin(phi_234);
      }
      else dotprod = 2.0 ;
    
      if (dotprod > optinfo.cos_tors_near_1_tol) val = 0.0 ;
      else if (dotprod < optinfo.cos_tors_near_neg1_tol) val = 180.0 ;
      else {
        val = acos(dotprod) ;
        // determine sign of torsions
        cross_product(eBC,eCD,tmp);
        dot_array(eAB,tmp,3,&dotprod);
        if (dotprod < 0) val *= -1;
        val *= 180.0 / _pi;
      }

      // extend domain of torsions so delta(vals) can be calculated
      if ((near_180 == -1) && (val > FIX_NEAR_180)) {
         //fprintf(outfile,"near_180=%d angle %15.10lf angle %15.10lf\n",
         //    near_180, val, -180.0 - (180.0 - val) );
        val = -180.0 - (180.0 - val);
      }
      else if ((near_180 == +1) && (val < -1*FIX_NEAR_180)) {
        //fprintf(outfile,"get_near_180(%d)=%d angle %15.10lf angle %15.10lf\n",
        //    i, get_near_180(i), angle, +180.0 + (180.0 + angle) );
        val = +180.0 + (180.0 + val);
      }
      return;
    }

    void fix_near_180(void) {
      if ( val > FIX_NEAR_180)
        near_180 = +1;
      else if ( val < -1*FIX_NEAR_180)
        near_180 = -1;
      else
        near_180 = 0;
      return;
    }

    // take geometry in au; store s vectors in 1/Angstroms
    void compute_s(double *geom) {
      int i,j;
      double rAB,rBC,rCD;
      double eAB[3], eBC[3], eCD[3], tmp[3], tmp2[3];
      double phiABC, phiBCD;
    
      for (j=0;j<3;++j) {
        eAB[j] = geom[3*B+j] - geom[3*A+j];
        eBC[j] = geom[3*C+j] - geom[3*B+j];
        eCD[j] = geom[3*D+j] - geom[3*C+j];
      }
    
      rAB = sqrt( SQR(eAB[0]) + SQR(eAB[1]) + SQR(eAB[2]) );
      rBC = sqrt( SQR(eBC[0]) + SQR(eBC[1]) + SQR(eBC[2]) );
      rCD = sqrt( SQR(eCD[0]) + SQR(eCD[1]) + SQR(eCD[2]) );
    
      scalar_div(rAB,eAB);
      scalar_div(rBC,eBC);
      scalar_div(rCD,eCD);
    
      phiABC = phiBCD = 0.0;
      for (j=0;j<3;++j) {
         phiABC += (-1.0 * eAB[j]) * eBC[j];
         phiBCD += (-1.0 * eBC[j]) * eCD[j];
      }
    
      phiABC = acos(phiABC);
      phiBCD = acos(phiBCD);
    
      cross_product(eAB,eBC,tmp);
      scalar_div(-1.0 * rAB * SQR(sin(phiABC)),tmp);
      s_A[0] = tmp[0];
      s_A[1] = tmp[1];
      s_A[2] = tmp[2];
    
      cross_product(eAB,eBC,tmp);
      scalar_mult((rBC-rAB*cos(phiABC))/(rBC*rAB*SQR(sin(phiABC))),tmp,3);
      cross_product(eCD,eBC,tmp2);
      scalar_mult(cos(phiBCD)/(rBC*SQR(sin(phiBCD))),tmp2,3);
      s_B[0] = tmp[0]+tmp2[0];
      s_B[1] = tmp[1]+tmp2[1];
      s_B[2] = tmp[2]+tmp2[2];
    
      cross_product(eCD,eBC,tmp);
      scalar_mult((rBC-rCD*cos(phiBCD))/(rBC*rCD*SQR(sin(phiBCD))),tmp,3);
      cross_product(eAB,eBC,tmp2);
      scalar_mult(cos(phiABC)/(rBC*SQR(sin(phiABC))),tmp2,3);
      s_C[0] = tmp[0]+tmp2[0];
      s_C[1] = tmp[1]+tmp2[1];
      s_C[2] = tmp[2]+tmp2[2];
    
      cross_product(eCD,eBC,tmp);
      scalar_div(-1.0*rCD*SQR(sin(phiBCD)),tmp);
      s_D[0] = tmp[0];
      s_D[1] = tmp[1];
      s_D[2] = tmp[2];

      for (j=0; j<3; ++j) {
        s_A[j] /= _bohr2angstroms;
        s_B[j] /= _bohr2angstroms;
        s_C[j] /= _bohr2angstroms;
        s_D[j] /= _bohr2angstroms;
      }
      return;
    }

    bool operator==(const tors_class & s2) const {
      if ( this->A == s2.A && this->B == s2.B && this->C == s2.C && this->D == s2.D)
        return true;
      /* allow reverse torsions ?
      else if ( this->A == s2.D && this->B == s2.C && this->C == s2.B && this->D == s2.A)
        return true; */
      else
        return false;
    };

};

}//} /* namespace psi::optking */

#endif
