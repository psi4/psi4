/*! \file
    \ingroup OPTKING
    \brief Class for pairs of orthogonal linear bend coordinates;
 lin with linval1 is toward an axis perpendicular with A-B-C with maximum
 overlap with the y-axis;
 lin with linval2 is orthogonal to this axis
*/

#ifndef _psi3_bin_optking_linb_h_
#define _psi3_bin_optking_linb_h_

#include <libchkpt/chkpt.h>

namespace psi { //namespace optking {

class linb_class {

 private:

    int id;
    int A;
    int B;
    int C;
    int linval; // 1 or 2 for lin1 or lin2
    double val;
    double s_A[3]; /* The s vector for atom A */
    double s_B[3];
    double s_C[3];
    double dummy[3]; // x,y,z of implicit dummy atom to orient linb

  public:

    friend class simples_class;

    linb_class(void) { }

    ~linb_class() { }

    linb_class(int id_in, int A_in, int B_in, int C_in, int linval_in, double val_in = 0.0) {
      id = id_in;
      B = B_in;
      // cononical is A < C 
      if (A_in < C_in) {
        A = A_in;
        C = C_in;
      }
      else {
        A = C_in;
        C = A_in;
      }
      if ((linval_in != 1) && (linval_in != 2))
        throw("linb_class() : unallowed val of linval");
      linval = linval_in;
      val = val_in;
    }

    void print(FILE *fp_out, int print_val) const {
      if (print_val)
        fprintf(fp_out,"    (%d %d %d %d %d) (%.8lf)\n", id,A+1,B+1,C+1,linval,val);
      else 
        fprintf(fp_out,"    (%d %d %d %d)\n", id,A+1,B+1,C+1);
    }

    void print_s(FILE *fp_out) const {
      fprintf(fp_out,"S vector for linb %d %d %d %d: atom A\n", A+1, B+1, C+1, linval);
      fprintf(fp_out,"(%16.10f,%16.10f,%16.10f)\n", s_A[0], s_A[1], s_A[2]);
      fprintf(fp_out,"S vector for linb %d %d %d %d: atom B\n", A+1, B+1, C+1, linval);
      fprintf(fp_out,"(%16.10f,%16.10f,%16.10f)\n", s_B[0], s_B[1], s_B[2]);
      fprintf(fp_out,"S vector for linb %d %d %d %d: atom C\n", A+1, B+1, C+1, linval);
      fprintf(fp_out,"(%16.10f,%16.10f,%16.10f)\n", s_C[0], s_C[1], s_C[2]);
    }

    void set_id(int i){ id = i;}
    void set_A(int i) { A = i;}
    void set_B(int i) { B = i;}
    void set_C(int i) { C = i;}
    void set_linval(int i){ linval = i;}
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
    void set_dummy(double x, double y, double z) {
      dummy[0] = x;  dummy[1] = y; dummy[2] = z;
    }


    int  get_id(void) const { return id;}
    int  get_A(void) const { return A;}
    int  get_B(void) const { return B;}
    int  get_C(void) const { return C;}
    int  get_linval(void) const { return linval;}
    double get_val(void) const { return val;}
    double get_val_A_or_rad(void) const  { 
      return (val/180*_pi);
    }

    double get_s_A(int i) const { return s_A[i]; }
    double get_s_B(int i) const { return s_B[i]; }
    double get_s_C(int i) const { return s_C[i]; }
    double get_dummy(int i) const { return dummy[i]; }

    int  get_atom(int a) const  {
      if (a==0) return A;
      else if (a==1) return B;
      else if (a==2) return C;
      else throw("linb_class::get_atom : atom index must be 0, 1 or 2.\n");
    }

    double get_s(int atom, int xyz) const  {
      if ( xyz < 0 || xyz > 2) throw ("linb_class::get_s() : xyz must be 0, 1 or 2");
      if (atom==0) return s_A[xyz];
      else if (atom==1) return s_B[xyz];
      else if (atom==2) return s_C[xyz];
      else throw("linb_class::get_s() : atom index must be 0, 1, or 2");
    }

    void compute(double *geom) {
      int j;
      double rBA,rBC,rBD,eBA[3],eBC[3],eBD[3],tmp[3],dotprod;
      double angle_ABD, angle_CBD, disp_size, F[3];
      double rBF, eBF[3];
 
      dummy[0] = dummy[1] = dummy[2] = 0;

      // first try placing dummy atoms according to default axes
      if (optinfo.dummy_axis_1 == 0)      dummy[0] = 1;
      else if (optinfo.dummy_axis_1 == 1) dummy[1] = 1;
      else if (optinfo.dummy_axis_1 == 2) dummy[2] = 1;

      // determine direction B->A
      for (j=0;j<3;++j)
        eBA[j] = geom[3*A+j] - geom[3*B+j];
      rBA = sqrt( SQR(eBA[0]) + SQR(eBA[1]) + SQR(eBA[2]) );
      scalar_div(rBA, eBA);

      // see if chosen dummy axis is on fragment line
      dot_array(eBA,dummy,3,&dotprod);
      if ((1 - fabs(dotprod)) < 1.e-5) {
        fprintf(outfile,"linear bend fragment is pointed toward dummy atom\n");
        fprintf(outfile,"try dummy_axis_1 = {1,2,3}\n");
        throw("linb_class::compute atoms pointed at dummy atom\n");
      }

      // 2nd dummy atom F (for linval == 2) is at (B + eBD x eBA)
      for (j=0;j<3;++j) {
        eBD[j] = dummy[j] - geom[3*B+j];
        eBA[j] = geom[3*A+j] - geom[3*B+j];
      }
      rBD = sqrt( SQR(eBD[0])+SQR(eBD[1])+SQR(eBD[2]) );
      rBA = sqrt( SQR(eBA[0])+SQR(eBA[1])+SQR(eBA[2]) );
      scalar_div(rBD,eBD);
      scalar_div(rBA,eBA);
      cross_product(eBD,eBA,eBF);

      // move dummy atoms far,far away
      disp_size = 1.0E9;

      if (linval == 2) { // dummy atom = F
        dummy[0] = geom[ 3*B+0 ] + disp_size * eBF[0];
        dummy[1] = geom[ 3*B+1 ] + disp_size * eBF[1];
        dummy[2] = geom[ 3*B+2 ] + disp_size * eBF[2];
      }
      else { // 1st dummy atom (for linval == 1) is at (B + eBA x eBF) 
        cross_product(eBA,eBF,eBD);
        dummy[0] = geom[ 3*B+0 ] + disp_size * eBD[0];
        dummy[1] = geom[ 3*B+1 ] + disp_size * eBD[1];
        dummy[2] = geom[ 3*B+2 ] + disp_size * eBD[2];
      }

      // angle = <A-B-D/F + <C-B-D/F
      // compute val of A-B-D/F
      for (j=0;j<3;++j) {
        eBA[j] = geom[3*A+j] - geom[3*B+j];
        eBC[j] = dummy[j] - geom[3*B+j];
      }
      rBA = sqrt( SQR(eBA[0])+SQR(eBA[1])+SQR(eBA[2]) );
      rBC = sqrt( SQR(eBC[0])+SQR(eBC[1])+SQR(eBC[2]) );
      scalar_div(rBA,eBA);
      scalar_div(rBC,eBC);
      dot_array(eBA,eBC,3,&dotprod);
      if (dotprod > 1.0) angle_ABD = 0.0;
      else if (dotprod < -1.0) angle_ABD = _pi;
      else angle_ABD = acos(dotprod)*180.0/_pi;
     //fprintf(outfile,"angle(A-B-D/F): %20.15lf\n",angle_ABD);

      // compute val of C-B-D/F
      for (j=0;j<3;++j) {
        eBA[j] = geom[3*C+j] - geom[3*B+j];
        eBC[j] = dummy[j] - geom[3*B+j];
      }
      rBA = sqrt( SQR(eBA[0])+SQR(eBA[1])+SQR(eBA[2]) );
      rBC = sqrt( SQR(eBC[0])+SQR(eBC[1])+SQR(eBC[2]) );
      scalar_div(rBA,eBA);
      scalar_div(rBC,eBC);
      dot_array(eBA,eBC,3,&dotprod);
      if (dotprod > (1.0-MIN_LIN_COS)) angle_CBD = 0.0;
      else if (dotprod < (-1.0+MIN_LIN_COS)) angle_CBD = _pi;
      else angle_CBD = acos(dotprod)*180.0/_pi;
     //fprintf(outfile,"angle(<C-B-D/F): %20.15lf\n",angle_CBD);

      val = angle_ABD+angle_CBD;
      return;
    }


    // s vectors point in direction of increasing internal coordinate val
    // so A and C retreat from D and B advances
    void compute_s(double *geom) {
      int j;
      double rAD, rBD, rCD, eAD[3], eBD[3], eCD[3];

      for (j=0;j<3;++j) {
        eAD[j] = dummy[j] - geom[3*A+j];
        eBD[j] = dummy[j] - geom[3*B+j];
        eCD[j] = dummy[j] - geom[3*C+j];
      }
      rAD = sqrt( SQR(eAD[0])+SQR(eAD[1])+SQR(eAD[2]) );
      rBD = sqrt( SQR(eBD[0])+SQR(eBD[1])+SQR(eBD[2]) );
      rCD = sqrt( SQR(eCD[0])+SQR(eCD[1])+SQR(eCD[2]) );
      scalar_div(rAD, eAD);
      scalar_div(rBD, eBD);
      scalar_div(rCD, eCD);

      // A and C go away from D
      set_s_A(-eAD[0], -eAD[1], -eAD[2]);
      set_s_C(-eCD[0], -eCD[1], -eCD[2]);
      // B goes toward D
      set_s_B( eBD[0],  eBD[1],  eBD[2]);

      return;
    }

    bool operator==(const linb_class & s2) const {
      if ( this->A == s2.A && this->B == s2.B && this->C == s2.C && this->linval == s2.linval)
        return true;
      else
        return false;
    };

};

}//} /* namespace psi::optking */

#endif
