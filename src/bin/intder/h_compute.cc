/*! \file
    \ingroup INTDER
    \brief Enter brief description of file here 
*/
#include "intco.h"
#include "params.h"
#include "displacements.h"
#include "atom.h"
#define EXTERN
#include "globals.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <physconst.h>

using namespace psi::intder;

namespace psi { namespace intder {
extern InternalCoordinates gIntCo;
extern Params gParams;
extern Displacements gDisplacements;
}}

void Stretch::Hijs(int disp, int k1, int k2, double **h11)
{
  int i, j;
  double t21;
  double *s1;

  s1 = new double[3];
  Stretch::Vect(disp, k1, k2, s1, &t21);

  t21 = 1 / t21;

  // Took the 7 "do" loops of Dr. Allen's and fit them into 2 "for" loops
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++) {
       h11[i][j] = -s1[i] * s1[j];

       if (i == j)
         h11[i][i] += 1;

       h11[i][j] * t21;
    }
}

//HIJKS1 from intder2000.f
void Stretch::Hijks(int disp, int k1, int k2, C3DMatrix *h111)
{
  int i, j, k;
  double **h11;
  double *s1; //v1 from intder2000.f
  double t21;
                                                                                           
  s1 = new double[3];
  h11 = init_matrix(3,3);
                                                                                           
  Stretch::Vect(disp, k1, k2, s1, &t21); //VECT1
	Stretch::Hijs(disp, k1, k2, h11);			//HIJS1

	for(k=0; k<3; k++)
		for(j=k; j<3; j++)
			for(i=j; i<3; i++)
				h111->Set(i, j, k, -(s1[i]*h11[k][j]+s1[j]*h11[k][i]+s1[k]*h11[j][i])/t21);

	fill3b(3, 3, h111);

  delete[] s1;
  free_matrix(h11,3);                                                                                     
}

void Stretch::H4th(int disp, int k1, int k2, C4DMatrix *h1111)
{
  double **h11 = init_matrix(3,3);
  C3DMatrix *h111 = new C3DMatrix(3,3,3);
  double *s1 = new double[3];
  double t21;
  double F;
  int l, k, j, i;

  Stretch::Vect(disp, k1, k2, s1, &t21);
  Stretch::Hijs(disp, k1, k2, h11);
  Stretch::Hijks(disp, k1, k2, h111);

  for (l=0; l<3; l++) {
    for (k=0; k<=l; k++) {
      for (j=0; j<=k; j++) {
        for (i=0; i<=j; i++) {
          F = h11[i][l]*h11[k][j] + h11[j][l]*h11[k][i] + h11[i][j]*h11[k][l];
          F+= s1[i]*h111->Get(j,k,l) + s1[j]*h111->Get(i,k,l) + s1[k]*h111->Get(i,j,l);
          F+= s1[l]*h111->Get(i,j,k);
          h1111->Set(i,j,k,l, -F/t21);
        }
       }
     }
   }

   fill4a(3, 3, h1111);

   delete h111;
   delete[] s1;
   free_matrix(h11, 3);
}

//H5TH1 from intder2000.f
void Stretch::H5th(int disp, int k1, int k2, C5DMatrix *h11111)
{
	int i, j, k, l, m;
	double *v1, **h11, t21;
	C3DMatrix h111(3,3,3);
	C4DMatrix h1111(3,3,3,3);
	
	v1 = new double[3];
	h11 = init_matrix(3, 3);

	Stretch::Vect(disp, k1, k2, v1, &t21); //VECT1
	Stretch::Hijs(disp, k1, k2, h11);			 //HIJS1
	Stretch::H4th(disp, k1, k2, &h1111);	 //H4TH1

	for(m=0; m<3; m++) {
		for(l=0; l<3; l++) {
			for(k=0; k<3; k++) {
				for(j=0; j<3; j++) {
					for(i=0; i<3; i++) {
					  h11111->Set(i, j, k, l, m, h11[i][l]*h111.Get(k, j, m) + h11[j][l]*h111.Get(k, i, m) 
                                    + h11[i][j]*h111.Get(k, l, m) + h11[k][j]*h111.Get(i, l, m)
                                    + h11[k][i]*h111.Get(j, l, m) + h11[k][l]*h111.Get(i, j, m)
                                    + h11[i][m]*h111.Get(k, j, l) + h11[j][m]*h111.Get(k, i, l)
                                    + h11[k][k]*h111.Get(i, j, l) + h11[l][k]*h111.Get(k, i, j)
                                    + v1[i]*h1111.Get(j, k, l, m) + v1[j]*h1111.Get(i, k, l, m)
                                    + v1[k]*h1111.Get(i, j, l, m) + v1[l]*h1111.Get(i, j, k, m)
                                    + v1[m]*h1111.Get(i, j, k, l));
						h11111->Set(i, j, k, l, m, -h11111->Get(i, j, k, l, m)/t21);
					}
				}
			}
		}
	}
//sum over m=1,3 is handled in fill4a(int, int, C5DMatrix)
	fill4a(3, 3, h11111);
	delete[] v1;
	free_matrix(h11, 3);
}


void Bend::Hijs(int disp, int k1, int k2, int k3, double **h11, double **h21, double **h31, double **h22, double **h32, double **h33)
{
  int i, j;
  double *s1, *s2, *s3, *e21, *e23;
  double **h11a, **h33a;
  double phi, t21, t23;
  double sphi, ctphi;
  double w1, w2, w3, w4, w5;

  s1 = new double[3];
  s2 = new double[3];
  s3 = new double[3];
  e21 = new double[3];
  e23 = new double[3];
  h11a = init_matrix(3, 3);
  h33a = init_matrix(3, 3);

  Bend::Vect(disp, k1, k2, k3, s1, s2, s3, &phi);
  Stretch::Vect(disp, k1, k2, e21, &t21);
  Stretch::Vect(disp, k3, k2, e23, &t23);

  Stretch::Hijs(disp, k1, k2, h11a);
  Stretch::Hijs(disp, k3, k2, h33a);

  sphi = sin(phi);
  ctphi = cos(phi) / sphi;
  w1 = ctphi;
  w2 = 1 / t21;
  w3 = w1 * w2;
  w4 = 1 / t23;
  w5 = w1 * w4;

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      h11[i][j] = h11a[i][j]*w3 - s1[i]*s1[j]*w1 - (e21[i]*s1[j] + s1[i]*e21[j])*w2;
      h33[i][j] = h33a[i][j]*w5 - s3[i]*s3[j]*w1 - (e23[i]*s3[j] + s3[i]*e23[j])*w4;
    }
  }

  w3 = 1 / (t21 * sphi);
  for (j = 0; j < 3; j++) {
    w4 = w2*e21[j] + w1*s1[j];
    for (i = 0; i < 3; i++) {
      h31[i][j] = -h33a[i][j]*w3 - s3[i]*w4;
      h21[i][j] = -(h11[i][j] + h31[i][j]);
      h32[i][j] = -(h31[i][j] + h33[i][j]);
    }
  }

  for (j = 0; j < 3; j++)
    for (i = 0; i < 3; i++)
      h22[i][j] = -(h21[j][i] + h32[i][j]);

  delete[] s1;
  delete[] s2;
  delete[] s3;
  delete[] e21;
  delete[] e23;
  free(h11a);
  free(h33a);
}

//HIJKS2 from intder2000.f
void Bend::Hijks(int disp, int k1, int k2, int k3, C3DMatrix *h111, C3DMatrix *h112,
                 C3DMatrix *h113, C3DMatrix *h123, C3DMatrix *h221, C3DMatrix *h222,
                 C3DMatrix *h223, C3DMatrix *h331, C3DMatrix *h332, C3DMatrix *h333)
{
	int i, j, k;
	double *s1, *s2, *s3, *e21, *e23;
	double **h11, **h21, **h31, **h22, **h32, **h33;
	C3DMatrix *h111a = new C3DMatrix(3,3,3);
	C3DMatrix *h333a = new C3DMatrix(3,3,3);
	double **h11a, **h33a;
	double phi, t21, t23, sphi, ctphi;
	double w1, w2, w3, w4, w5, w6, w7, w8, w9, w10;

	s1 = new double(3);
	s2 = new double(3);
	s3 = new double(3);
	e21 = new double(3);
	e23 = new double(3);
	h11 = init_matrix(3, 3);
	h21 = init_matrix(3, 3);
	h31 = init_matrix(3, 3);
	h22 = init_matrix(3, 3);
	h32 = init_matrix(3, 3);
	h33 = init_matrix(3, 3);
	h11a = init_matrix(3,3);
	h33a = init_matrix(3,3);

	Bend::Vect(disp, k1, k2, k3, s1, s2, s3, &phi);	//VECT1
	Stretch::Vect(disp, k1, k2, e21, &t21);					//VECT2
	Stretch::Vect(disp, k3, k2, e23, &t23);					//VECT3

	Stretch::Hijs(disp, k1, k2, h11a);	//HIJS1
	Stretch::Hijs(disp, k3, k2, h33a);	//HIJS1
	Bend::Hijs(disp, k1, k2, k3, h11, h21, h31, h22, h32, h33);	//HIJS2
	Stretch::Hijks(disp, k1, k2, h111a);//HIJKS1
	Stretch::Hijks(disp, k3, k2, h333a);//HIJKS1

	sphi = sin(phi);
	ctphi = cos(phi)/sphi;
	w1 = 1/t21;
	w2 = 1/t23;
	w3 = ctphi*w1;
	w4 = ctphi*w2;

	for(k=0; k<3; k++) {
		w5 = s1[k]*ctphi + e21[k]*w1;
		w6 = e21[k]*w3;
		w7 = s1[k]*w1;
		w8 = s3[k]*ctphi + e23[k]*w2;
		w9 = e23[k]*w4;
		w10 = s3[k]*w2;
		for(j=0; j<3; j++) {
			for(i=0; i<3; i++) {
				h221->Set(i, j, k, w5*h11[i][j] + s1[i]*s1[j]*w6 + h11a[i][j]*w7);
				h223->Set(i, j, k, w8*h33[i][j] + s3[i]*s3[j]*w9 + h33a[i][j]*w10);
			}
		}
	}

	for(k=0; k<3; k++) {
		for(j=k; j<3; j++) {
			for(i=j; j<3; i++) {
				h111->Set(i, j, k, -(h221->Get(i, j, k) + h221->Get(j, k, i) + h221->Get(i, k, j)
                            + s1[i]*s1[j]*s1[k] + h111a->Get(i, j, k)*w3));
				h333->Set(i, j, k, -(h223->Get(i, j, k) + h223->Get(j, k, i) + h223->Get(i, k, j)
                            + s3[i]*s3[j]*s3[k] + h333a->Get(i, j, k)*w4));
			}
		}
	}

	fill3b(3, 3, h111);
	fill3b(3, 3, h333);

	for(i=0; i<3; i++) {
		w3 = s1[i]*ctphi + e21[i]*w1;
		w4 = s3[i]*ctphi + e23[i]*w2;
		for(j=0; j<3; j++) {
			for(k=0; k<3; k++) {
				h221->Set(i, j, k, w3*h31[k][j]);
				h223->Set(i, j, k, w4*h31[j][k]);
			}
		}
	}
	
	w3 = 1/(sphi*sphi);

	for(k=0; k<3; k++) {
		for(j=0; j<3; j++) {
			for(i=0; i<3; i++) {
				h113->Set(i, j, k, s3[k]*(s1[i]*s1[j] - h11a[i][j]*w1)*w3 - h221->Get(i, j, k)
                           - h221->Get(j, i, k));
				h331->Set(i, j, k, s1[k]*(s3[i]*s3[j] - h33a[i][j]*w2)*w3 - h223->Get(i, j, k)
                           - h223->Get(j, i, k));
			}
		}
	}

	for(k=0; k<3; k++) {
		for(j=0; j<3; j++) {
			for(i=0; i<3; i++) {
				h123->Set(i, j, k, -(h331->Get(j, k, i) + h113->Get(i, j, k)));
				h112->Set(i, j, k, -(h111->Get(i, j, k) + h113->Get(i, j, k)));
				h332->Set(i, j, k, -(h333->Get(i, j, k) + h331->Get(i, j, k)));
			}
		}
	}

	for(k=0; k<3; k++) {
		for(j=0; j<3; j++) {
			for(i=0; i<3; i++) {
				h221->Set(j, k, i, -(h123->Get(i, j, k) + h112->Get(i, k, j)));
				h223->Set(j, k, i, -(h331->Get(i, j, k) + h123->Get(j, k, i)));
			}
		}
	}

	for(k=0; k<3; k++)
		for(j=0; j<3; j++)
			for(i=0; i<3; i++)
				h222->Set(i, j, k, -(h223->Get(j, k, i) + h221->Get(j, k, i)));

	delete[] s1;
	delete[] s2;
	delete[] s3;
	delete[] e21;
	delete[] e23;
	delete h111a;
	delete h333a;
	free_matrix(h11, 3);
	free_matrix(h21, 3);
	free_matrix(h31, 3);
	free_matrix(h22, 3);
	free_matrix(h32, 3);
	free_matrix(h33, 3);
	free_matrix(h11a, 3);
	free_matrix(h33a, 3);
}

void Bend::H4th(int disp, int k1, int k2, int k3, C4DMatrix *h1111, C4DMatrix *h1112, C4DMatrix *h1113,
                C4DMatrix *h1122, C4DMatrix *h1123, C4DMatrix *h1133, C4DMatrix *h1222,
                C4DMatrix *h1223, C4DMatrix *h1233, C4DMatrix *h1333, C4DMatrix *h2222,
                C4DMatrix *h2223, C4DMatrix *h2233, C4DMatrix *h2333, C4DMatrix *h3333)
{
  double *s1 = new double[3];
  double *s2 = new double[3];
  double *s3 = new double[3];
  double **h11 = init_matrix(3,3);
  double **h21 = init_matrix(3,3);
  double **h31 = init_matrix(3,3);
  double **h22 = init_matrix(3,3);
  double **h32 = init_matrix(3,3);
  double **h33 = init_matrix(3,3);
  C3DMatrix h111(3,3,3), h112(3,3,3), h113(3,3,3), h123(3,3,3), h221(3,3,3), h222(3,3,3);
  C3DMatrix h223(3,3,3), h331(3,3,3), h332(3,3,3), h333(3,3,3);
  C4DMatrix q1111(3,3,3,3), q3333(3,3,3,3);
  C5DMatrix q11111(3,3,3,3,3), q33333(3,3,3,3,3);

  double phi, cscp, cotp, r1, r3;
  int m, l, k, j, i;

  Bend::Vect(disp, k1, k2, k3, s1, s2, s3, &phi);
  Bend::Hijs(disp, k1, k2, k3, h11, h21, h31, h22, h32, h33);
  Bend::Hijks(disp, k1, k2, k3, &h111, &h112, &h113, &h123, &h221, &h222, &h223, &h331, &h332, &h333);

  cscp = 1.0 / sin(phi);
  cotp = cos(phi)*cscp;

  for (l=0; l<3; l++) {
    for (k=0; k<3; k++) {
      for (j=0; j<3; j++) {
        for (i=0; i<3; i++) {
          h1111->Set(i,j,k,l, (s1[i]*s1[j]*s1[k]*s1[l] - h111.Get(j,k,l)*s1[i] -
                               h111.Get(i,j,k)*s1[l]   - h111.Get(i,j,l)*s1[k] - h111.Get(i,k,l)*s1[j] -
                               h11[i][j]*h11[k][l]     - h11[i][l]*h11[k][j]   - h11[i][k]*h11[j][l]) * cotp +
                               h11[i][j]*s1[k]*s1[l]   + h11[i][k]*s1[j]*s1[l] +
                               h11[i][l]*s1[j]*s1[k]   + h11[j][k]*s1[i]*s1[l] +
                               h11[j][l]*s1[k]*s1[i]   + h11[k][l]*s1[j]*s1[i]);
          h1113->Set(i,j,k,l, (s1[i]*s1[j]*s1[k]*s3[l] - h113.Get(j,k,l)*s1[i] -
                               h111.Get(i,j,k)*s3[l]   - h113.Get(i,j,l)*s1[k] - h113.Get(i,k,l)*s1[j] -
                               h11[i][j]*h31[l][k]     - h31[l][i]*h11[k][j]   - h11[i][k]*h31[l][j]) * cotp +
                               h11[i][j]*s1[k]*s3[l]   + h11[i][k]*s1[j]*s3[l] +
                               h31[l][i]*s1[j]*s1[k]   + h11[j][k]*s1[i]*s3[l] +
                               h31[l][j]*s1[k]*s1[i]   + h31[l][k]*s1[j]*s1[i]);
          h1133->Set(i,j,k,l, (s1[i]*s1[j]*s3[k]*s3[l] - h331.Get(l,k,j)*s1[i] -
                               h113.Get(i,j,k)*s3[l]   - h113.Get(i,j,l)*s3[k] - h331.Get(l,k,j)*s1[j] -
                               h11[i][j]*h33[l][k]     - h31[l][i]*h31[k][j]   - h31[k][i]*h31[l][j]) * cotp +
                               h11[i][j]*s3[k]*s3[l]   + h31[k][i]*s1[j]*s3[l] +
                               h31[l][i]*s1[j]*s3[k]   + h31[k][j]*s1[i]*s3[l] +
                               h31[l][j]*s3[k]*s1[i]   + h33[l][k]*s1[j]*s1[i]);
          h1333->Set(i,j,k,l, (s1[i]*s3[j]*s3[k]*s3[l] - h333.Get(l,k,j)*s1[i] -
                               h331.Get(k,j,i)*s3[l]   - h331.Get(l,j,i)*s3[k] - h331.Get(l,k,j)*s3[j] -
                               h31[j][i]*h33[l][k]     - h31[l][i]*h33[k][j]   - h31[k][i]*h33[l][j]) * cotp +
                               h31[j][i]*s3[k]*s3[l]   + h31[k][i]*s3[j]*s3[l] +
                               h31[l][i]*s3[j]*s3[k]   + h33[k][j]*s1[i]*s3[l] +
                               h33[l][j]*s3[k]*s1[i]   + h33[l][k]*s3[j]*s1[i]);
          h3333->Set(i,j,k,l, (s3[i]*s3[j]*s3[k]*s3[l] - h333.Get(l,k,j)*s3[i] -
                               h333.Get(k,j,i)*s3[l]   - h333.Get(l,j,i)*s3[k] - h333.Get(l,k,j)*s3[j] -
                               h33[j][i]*h33[l][k]     - h33[l][i]*h33[k][j]   - h33[k][i]*h33[l][j]) * cotp +
                               h33[j][i]*s3[k]*s3[l]   + h33[k][i]*s3[j]*s3[l] +
                               h33[l][i]*s3[j]*s3[k]   + h33[k][j]*s3[i]*s3[l] +
                               h33[l][j]*s3[k]*s3[i]   + h33[l][k]*s3[j]*s3[i]);
        }
      }
    }
  }

  Stretch::Vect(disp, k1, k2, s1, &r1);
  Stretch::Vect(disp, k3, k2, s3, &r3);
  Stretch::Hijs(disp, k1, k2, h11);
  Stretch::Hijs(disp, k3, k2, h33);
  Stretch::Hijks(disp, k1, k2, &h111);
  Stretch::Hijks(disp, k3, k2, &h333);
  Stretch::H4th(disp, k1, k2, &q1111);
  Stretch::H4th(disp, k3, k2, &q3333);
  Stretch::H5th(disp, k1, k2, &q11111);
  Stretch::H5th(disp, k3, k2, &q33333);

  for (m=0; m<3; m++) {
    for (l=0; l<3; l++) {
      for (k=0; k<3; k++) {
        for (j=0; j<3; j++) {
          for (i=0; i<3; i++) {
            h1111->Set(i,j,k,l, h1111->Get(i,j,k,l) - cscp*q11111.Get(m,i,j,k,l)*s3[m]);
            h1113->Set(i,j,k,l, h1113->Get(i,j,k,l) - cscp*q1111.Get(m,i,j,k)*h33[l][m]);
            h1133->Set(i,j,k,l, h1133->Get(i,j,k,l) - cscp*h111.Get(m,i,j)*h333.Get(k,l,m));
            h1333->Set(i,j,k,l, h1333->Get(i,j,k,l) - cscp*h11[m][i]*q3333.Get(j,k,l,m));
            h3333->Set(i,j,k,l, h3333->Get(i,j,k,l) - cscp*s1[m]*q33333.Get(i,j,k,l,m));
          }
        }
      }
    }
  }

  for (l=0; l<3; l++) {
    for (k=0; k<3; k++) {
      for (j=0; j<3; j++) {
        for (i=0; i<3; i++) {
          h1112->Set(i,j,k,l, -h1111->Get(i,j,k,l) - h1113->Get(i,j,k,l));
          h1123->Set(i,j,k,l, -h1113->Get(i,j,k,l) - h1133->Get(i,j,k,l));
          h1233->Set(i,j,k,l, -h1133->Get(i,j,k,l) - h1333->Get(i,j,k,l));
          h2333->Set(i,j,k,l, -h1333->Get(i,j,k,l) - h3333->Get(i,j,k,l));
        }
      }
    }
  }

  for (l=0; l<3; l++) {
    for (k=0; k<3; k++) {
      for (j=0; j<3; j++) {
        for (i=0; i<3; i++) {
          h1122->Set(i,j,k,l, -h1112->Get(i,j,k,l) - h1123->Get(i,j,k,l));
          h1223->Set(i,j,k,l, -h1123->Get(i,j,k,l) - h1233->Get(i,j,k,l));
          h2233->Set(i,j,k,l, -h1233->Get(i,j,k,l) - h2333->Get(i,j,k,l));
        }
      }
    }
  }

  for (l=0; l<3; l++) {
    for (k=0; k<3; k++) {
      for (j=0; j<3; j++) {
        for (i=0; i<3; i++) {
          h1222->Set(i,j,k,l, -h1122->Get(i,j,k,l) - h1223->Get(i,j,k,l));
          h2223->Set(i,j,k,l, -h1223->Get(i,j,k,l) - h2233->Get(i,j,k,l));
        }
      }
    }
  }

  for (l=0; l<3; l++) {
    for (k=0; k<3; k++) {
      for (j=0; j<3; j++) {
        for (i=0; i<3; i++) {
          h2222->Set(i,j,k,l, -h1222->Get(i,j,k,l) - h2223->Get(i,j,k,l));
        }
      }
    }
  }

  free_matrix(h11, 3);
  free_matrix(h21, 3);
  free_matrix(h31, 3);
  free_matrix(h22, 3);
  free_matrix(h32, 3);
  free_matrix(h33, 3);
  delete[] s1;
  delete[] s2;
  delete[] s3;
}


void Linear1::Hijs(int disp, int k1, int k2, int k3, int k4, double **h11, double **h21, double **h31, double **h22, double **h32, double **h33)
{
  double *s1, *s2, *s3, *e21, *e23, *ea;
  double **h11a, **h33a, **em;
  double th, t21, t23, d, tanth, costh;
  double it21, it23;
  int j, i, k;
  Atom *a4 = NULL;

  s1  = new double[3];
  s2  = new double[3];
  s3  = new double[3];
  e21 = new double[3];
  e23 = new double[3];
  ea  = new double[3];

  h11a = init_matrix(3, 3);
  h33a = init_matrix(3, 3);
  em   = init_matrix(3, 3);

  Linear1::Vect(disp, k1, k2, k3, k4, s1, s2, s3, &th);
  Stretch::Vect(disp, k1, k2, e21, &t21);
  Stretch::Vect(disp, k3, k2, e23, &t23);

  Stretch::Hijs(disp, k1, k2, h11a);
  Stretch::Hijs(disp, k3, k2, h33a);

  a4 = gDisplacements.displacement(disp)->atom(k4);
  ea[0] = a4->getX();
  ea[1] = a4->getY();
  ea[2] = a4->getZ();

  d = 1/(sqrt(dot_prod(ea, ea)));
  ea[0] *= d;
  ea[1] *= d;
  ea[2] *= d;

  tanth = tan(th);
  costh = cos(th);

  // off diagonal
  em[1][0] = ea[2];
  em[2][0] = -ea[1];
  em[2][1] = ea[0];
  em[0][1] = -em[1][0];
  em[0][2] = -em[2][0];
  em[1][2] = -em[2][1];

  // diagonal
  em[0][0] = em[1][1] = em[2][2] = 0.0;

  for (j=0; j<3; j++) {
    for (i=0; i<3; i++) {
      h22[i][j] = 0.0;
      for (k=0; k<3; k++) {
        h22[i][j] += em[i][k]*h33a[k][j];
      }
    }
  }

  it21 = 1/t21;
  it23 = 1/t23;

  for (j=0; j<3; j++) {
    for (i=0; i<3; i++) {
      h11[i][j] = (-h11a[i][j]*it21 + s1[i]*s1[j])*tanth - (e21[i]*s1[j] + s1[i]*e21[j])*it21;
      h31[i][j] = (h22[j][i]/costh - s3[i]*e21[j])*it21 + s3[i]*s1[j]*tanth;
      h33[i][j] = (-h33a[i][j]*it23 + s3[i]*s3[j])*tanth - (e23[i]*s3[j] + s3[i]*e23[j])*it23;
      h21[i][j] = -(h11[i][j] + h31[i][j]);
      h32[i][j] = -(h31[i][j] + h33[i][j]);
    }
  }

  for (j=0; j<3; j++) {
    for (i=0; i<3; i++) {
      h22[i][j] = -(h21[j][i] + h32[i][j]);
    }
  }

  delete s1;
  delete s2;
  delete s3;
  delete e21;
  delete e23;
  delete ea;

  free_matrix(h11a, 3);
  free_matrix(h33a, 3);
  free_matrix(em, 3);
}

void Linear1::Hijks(int disp, int k1, int k2, int k3, int k4, 
                    C3DMatrix *h111, C3DMatrix *h112, C3DMatrix *h113, C3DMatrix *h123,
                    C3DMatrix *h221, C3DMatrix *h222, C3DMatrix *h223, C3DMatrix *h331,
                    C3DMatrix *h332, C3DMatrix *h333)
{
	int i, j, k;
	double *s1, *s2, *s3, *e21, *e23;
	double **h11, **h21, **h31, **h22, **h32, **h33;
	C3DMatrix *h111a = new C3DMatrix(3,3,3);
	C3DMatrix *h333a = new C3DMatrix(3,3,3);
	double **h11a, **h33a;
	double theta, t21, t23;
	Atom *a4 = NULL;
	double d, *ea, tanth, costh, w1, w2, w3, w4, w5, w6;

	s1 = new double(3);
	s2 = new double(3);
	s3 = new double(3);
	e21 = new double(3);
	e23 = new double(3);
	ea = new double(3);

	h11 = init_matrix(3,3);
	h21 = init_matrix(3,3);
	h31 = init_matrix(3,3);
	h22 = init_matrix(3,3);
	h32 = init_matrix(3,3);
	h33 = init_matrix(3,3);
	h11a = init_matrix(3,3);
	h33a = init_matrix(3,3);

	Linear1::Vect(disp, k1, k2, k3, k4, s1, s2, s3, &theta); //VECT3
	Stretch::Vect(disp, k1, k2, e21, &t21);	//VECT1
	Stretch::Vect(disp, k3, k2, e23, &t23);	//VECT1

	Stretch::Hijs(disp, k1, k2, h11a);
	Stretch::Hijs(disp, k2, k2, h33a);

  a4 = gDisplacements.displacement(disp)->atom(k4);
  ea[0] = a4->getX();
  ea[1] = a4->getY();
  ea[2] = a4->getZ();
                                                                                    
  d = 1/(sqrt(dot_prod(ea, ea)));
  ea[0] *= d;
  ea[1] *= d;
  ea[2] *= d;

	Linear1::Hijs(disp, k1, k2, k3, k4, h11, h21, h31, h22, h32, h33);  //HIJS3
	Stretch::Hijks(disp, k1, k2, h111a); //HIJKS1
	Stretch::Hijks(disp, k3, k2, h333a); //HIJKS1

	tanth = tan(theta);
	costh = cos(theta);
	w1 = 1/t21;
	w2 = 1/t23;
	w3 = tanth*w1;
	w4 = tanth*w2;

	for(k=0; k<3; k++) {
		for(j=0; j<3; j++) {
			for(i=0; i<3; i++) {
				h221->Set(i, j, k, h11[i][j]*(s1[k]*tanth - e21[k]/t21) 
                           + s1[k]*s1[j]*e21[i]*tanth/t21
                           - (h11a[i][j]*s1[k])/t21);
				h223->Set(i, j, k, h33[i][j]*(s3[k]*tanth - e23[k]/t23)
                           + s3[k]*s3[j]*e23[i]*tanth/t23
                           - (h33a[i][j]*s3[k])/t23);
			}
		}
	}

  for(k=0; k<3; k++) {
    for(j=0; j<3; j++) {
      for(i=0; i<3; i++) {
				h111->Set(i, j, k, h221->Get(i, j, k) + h221->Get(j, k, j) + h221->Get(k, i, j)
                           + s1[i]*s1[j]*s1[k] - h111a->Get(i, j, k)*w3);
				h333->Set(i, j, k, h223->Get(i, j, k) + h223->Get(j, k, i) + h223->Get(k, i, j)
                           + s3[i]*s3[j]*s3[k] - h333a->Get(i, j, k)*w4);
			}
		}
	}

	fill3b(3, 3, h111);
	fill3b(3, 3, h333);

	for(i=0; i<3; i++) {
		w5 = s1[i]*tanth - e21[i]*w1;
		w6 = s3[i]*tanth - e23[i]*w2;
		for(j=0; j<3; j++) {
			for(k=0; k<3; k++) {
				h221->Set(i, j, k, w5*h31[k][j]);
				h223->Set(i, j, k, w6*h31[j][k]);
			}
		}
	}
	
	w5 = 1/(costh*costh);

	for(k=0; k<3; k++) {
		for(j=0; j<3; j++) {
			for(i=0; i<3; i++) {
				h113->Set(i, j, k, s3[k]*(s1[i]*s1[j] - h11a[i][j]*w1)*w5 
                           + h221->Get(i, j, k) + h221->Get(j, i, k));
				h331->Set(i, j, k, s1[k]*(s3[i]*s3[j] - h33a[i][j]*w2)*w5
                           + h223->Get(i, j, k) + h223->Get(j, i, k));
			}
		}
	}

  for(k=0; k<3; k++) {
    for(j=0; j<3; j++) {
      for(i=0; i<3; i++) {
				h123->Set(i, j, k, -(h331->Get(j, k, i) + h113->Get(i, j, k)));
				h112->Set(i, j, k, -(h111->Get(i, j, k) + h113->Get(i, j, k)));
				h332->Set(i, j, k, -(h333->Get(i, j, k) + h331->Get(i, j, k)));
			}
		}
	}

  for(k=0; k<3; k++) {
    for(j=0; j<3; j++) {
      for(i=0; i<3; i++) {
				h221->Set(j, k, i, -(h123->Get(i, j, k) + h112->Get(i, k, j)));
				h223->Set(j, k, i, -(h332->Get(i, j, k) + h123->Get(j, k, i)));
			}
		}
	}

  for(k=0; k<3; k++)
    for(j=0; j<3; j++)
      for(i=0; i<3; i++)
				h222->Set(i, j, k, -(h223->Get(j, k, i) + h221->Get(j, k, i)));

	delete[] s1;
	delete[] s2;
	delete[] s3;
	delete[] e21;
	delete[] e23;
	delete[] ea;

	free_matrix(h11, 3);
	free_matrix(h21, 3);
  free_matrix(h31, 3);
  free_matrix(h22, 3);
  free_matrix(h32, 3);
  free_matrix(h33, 3);
  free_matrix(h11a, 3);
  free_matrix(h33a, 3);
  
}


//THIS IS HIJS7 in intder2000.f
void OutOfPlane::Hijs(int disp, int k1, int k2, int k3, int k4,
                      double **h11, double **h21, double **h31, double **h41,
                      double **h22, double **h32, double **h42, double **h33,
                      double **h43, double **h44)
{
  int i, j;
  double *s1  = new double[3];
  double *s2  = new double[3];
  double *s3  = new double[3];
  double *s4  = new double[3];
  double *s5;
  double *s6;
  double *bp3 = new double[3];
  double *bp4 = new double[3];
  double *e21 = new double[3];
  double *e23 = new double[3];
  double *e24 = new double[3];
  double t21, t23, t24, phi, gamma;
  double **hp43, **hp44, **cp31, **cp41, **cp43;
  double sp, cp, tp, sg, cg, tg;
  double c21, c23, c24, ctp, c11;
  double c312, c311, c313, c411;
  double c3, c4;
  double c331, c332;
  double c441, c442, c431, c432, c434, c435, c436;
  double xj;

  hp43 = init_matrix(3, 3);
  hp44 = init_matrix(3, 3);
  cp31 = init_matrix(3, 3);
  cp41 = init_matrix(3, 3);
  cp43 = init_matrix(3, 3);

  Stretch::Vect(disp, k1, k2, e21, &t21);
  Stretch::Vect(disp, k3, k2, e23, &t23);
  Stretch::Vect(disp, k4, k2, e24, &t24);

  Bend::Vect(disp, k3, k2, k4, bp3, s2, bp4, &phi);

  OutOfPlane::Vect(disp, k1, k2, k3, k4, s1, s2, s3, s4, &gamma);

  Bend::Hijs(disp, k3, k2, k4, h11, h21, h43, h22, h32, hp44);

  s5 = vect_prod(e23, e24);
  s6 = vect_prod(e24, e21);

  mat1(cp31, e24);
  mat1(cp41, e23);
  mat1(cp43, e21); 

  sp = sin(phi);
  cp = cos(phi);
  tp = sp/cp;
  sg = sin(gamma);
  cg = cos(gamma);
  tg = sp/cp;

  c21  = 1/t21;
  c23  = 1/t23;
  c24  = 1/t24;
  ctp  = 1/tp;
  c11  = tg*c21*c21;
  c312 = c21/(cg*sp);
  c311 = c312*c23;
  c313 = c312*ctp;
  c411 = c312*c24;
  c3   = c23/sp;
  c4   = c24/sp;
  c331 = t24*c3;
  c332 = c331*tg;
  c441 = t23*c4;
  c442 = c441*tg;
  c431 = c3*c24/cg;
  c432 = tg;
  c434 = tg*c3;
  c435 = t24*c3;
  c436 = c435*tg;

  for (j=0; j<3; j++) {
    for (i=j; i<3; i++) {
      h11[i][j] = s1[j]*(tg*s1[i] - e21[i]*c21) - s1[i]*e21[j]*c21 + e21[i]*e21[j]*c11;
      if (i == j)
        h11[i][j] -= c11;
      h11[j][i] = h11[i][j];

      h33[i][j] = s3[i]*bp4[j]*c331 + hp43[j][i]*c332 + s3[j]*(tg*s3[i]-e23[i]*c23-bp3[i]*ctp);
      h33[j][i] = h33[i][j];

      h44[i][j] = s4[i]*bp3[j]*c441 + hp43[i][j]*c442 + s4[j]*(tg*s4[i] - e24[i]*c24 - bp4[i]*ctp);
      h44[j][i] = h44[i][j];
    }
  }

  for (j=0; j<3; j++) {
    xj = tg*s1[j] - e21[j]*c21;
    for (i=0; i<3; i++) {
      h31[i][j] = s3[i]*xj - cp31[i][j]*c311 - e23[i]*s5[j]*c311 - bp3[i]*s5[j]*c313;
      h41[i][j] = s4[i]*xj + cp41[i][j]*c411 - e24[i]*s5[j]*c411 - bp4[i]*s5[j]*c313;
      h21[i][j] = -(h11[i][j] + h31[i][j] + h41[i][j]);
      h43[i][j] = (cp43[j][i] - e24[i]*s6[j])*c431 + s3[j]*s4[i]*c432 - s3[j]*bp4[i]*ctp + e24[i]*bp4[j]*c434 + s4[i]*bp4[j]*c435 + hp44[i][j]*c436;
    }
  }

  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      h32[i][j] = -(h31[i][j] + h33[i][j] + h43[j][i]);
      h42[i][j] = -(h41[i][j] + h43[i][j] + h44[i][j]);
      h22[i][j] = -(h21[j][i] + h32[i][j] + h42[i][j]);
    }
  }

  delete[] s1;
  delete[] s2;
  delete[] s3;
  delete[] s4;
  delete[] bp3;
  delete[] bp4;
  delete[] e21;
  delete[] e23;
  delete[] e24;
                                                                                                                                                             
  free_matrix(hp43, 3);
  free_matrix(hp44, 3);
  free_matrix(cp31, 3);
  free_matrix(cp41, 3);
  free_matrix(cp43, 3);
}

void OutOfPlane::Hijks(int disp, int k1, int k2, int k3, int k4, C3DMatrix *h111, C3DMatrix *h112, 
		       C3DMatrix *h221, C3DMatrix *h222, C3DMatrix *h113, C3DMatrix *h123, 
		       C3DMatrix *h223, C3DMatrix *h331, C3DMatrix *h332, C3DMatrix *h333,
		       C3DMatrix *h411, C3DMatrix *h421, C3DMatrix *h422, C3DMatrix *h431,
		       C3DMatrix *h432, C3DMatrix *h433, C3DMatrix *h441, C3DMatrix *h442,
		       C3DMatrix *h443, C3DMatrix *h444)
{
  int i,j,k;
  double *s1, *s2, *s3, *s4, *e21, *e23, *e24;
  double *bp3, *bp4;
  double **hp33, **hp43, **hp44;
  double **ht11;
  double **cp21, **cp24, *cp2124;
  double **h11, **h21, **h31, **h41, **h22, **h32, **h42, **h33, **h43, **h44;

  double t21, t23, t24, phi, gamma, cg, sg, tg, cp, sp, tp, ctp, c1;
  double s2g, c3, c4, c2p, c1111, c1112, c3331, c3332, c3333, c3335, c3334;
  double c4411, c4412, c431, c4311, c4312, c4313, c4431, c4442, c4441;
  
  C3DMatrix *ht111 = new C3DMatrix(3,3,3);
  C3DMatrix *hp334 = new C3DMatrix(3,3,3);
  C3DMatrix *hp443 = new C3DMatrix(3,3,3);
  C3DMatrix *prod = new C3DMatrix(3,3,3);
  
  s1 = new double[3];
  s2 = new double[3];
  s3 = new double[3];
  s4 = new double[3];
  e21 = new double[3];
  e23 = new double[3];
  e24 = new double[3];
  bp3 = new double[3];
  bp4 = new double[3];
  
  h11 = init_matrix(3, 3);
  h21 = init_matrix(3, 3);
  h31 = init_matrix(3, 3);
  h41 = init_matrix(3, 3);
  h22 = init_matrix(3, 3);
  h32 = init_matrix(3, 3);
  h42 = init_matrix(3, 3);
  h33 = init_matrix(3, 3);
  h43 = init_matrix(3, 3);
  h44 = init_matrix(3, 3);
	
  hp33 = init_matrix(3,3);
  hp43 = init_matrix(3,3);
  hp44 = init_matrix(3,3);
  ht11 = init_matrix(3,3);
  cp21 = init_matrix(3,3);
  cp24 = init_matrix(3,3);
  cp24 = init_matrix(3,3);
  cp2124 = new double[3];

  Stretch::Vect(disp, k1, k2, e21, &t21);
  Stretch::Vect(disp, k3, k2, e23, &t23);
  Stretch::Vect(disp, k4, k2, e24, &t24);
  Bend::Vect(disp, k3, k2, k4, bp3, s2, bp4, &phi); 
  OutOfPlane::Vect(disp, k1, k2, k3, k4, s1, s2, s3, s4, &gamma);
  
  Stretch::Hijs(disp, k1, k2, ht11);
  Bend::Hijs(disp, k3, k2, k4, hp33, h32, hp43, h22, h42, hp44);
  OutOfPlane::Hijs(disp, k1, k2, k3, k4, h11, h21, h31, h41, h22, h32, h42, h33, h43, h44);
  Stretch::Hijks(disp, k1, k2, ht111);
  Bend::Hijks(disp, k3, k2, k4, h111, h112, hp334, h123, h221, h222, h223, hp443, h332, h333);

  mat1(cp21, e21);
  mat1(cp24, e24);
    
  cp2124 = vect_prod(e21, e24);
  prod = tripro();

  cg = cos(gamma);
  sg = sin(gamma);
  tg = tan(gamma);
  cp = cos(phi);
  sp = sin(phi);
  tp = tan(phi);
  ctp = cp / sp;
  c1 = 1 / t21;
  s2g = 1 / (pow(cg,2));
  c3 = 1 / t23;
  c4 = 1 / t24;
  c2p = 1 / (pow(sg,2));
  c1111 = tg * c1;
  c1112 = s2g * c1;
  c3331 = t24 * c3 / sp;
  c3332 = c3331 * tg;
  c3333 = c3331 * s2g;
  c3335 = pow(c3,2);
  c3334 = 2.0 * c3335;
  c4411 = t23 * c4 / sp;
  c4412 = c4411 * s2g;
  c431 = c3 * c4 / (cg * sp);
  c4311 = c431 * c1;
  c4312 = c431 * tg;
  c4313 = c3333 * c4;
  c4431 = c4411 * tg;
  c4442 = pow(c4,2);
  c4441 = 2.0 * c4442;
  
  for(i = 0; i < 3; i++)
    for(j = 0; j < i; j++)
      for(k = 0; k < j; k++) {
	h111->Set(i,j,k, s2g * s1[i] * s1[j] * s1[k] + tg * h11[i][j] * s1[k] + tg * h11[i][k] * s1[j] + c1111 * 
		  (e21[i] * s1[j] * s1[k] - ht111->Get(i,j,k)) - c1112 * ht11[j][k] * s1[i] - c1 * 
		  (e21[i] * h11[j][k] + e21[j] * h11[i][k] + e21[k] * h11[i][j] + ht11[i][j] * s1[k] + ht11[i][k] * s1[j]));
	h333->Set(i,j,k, c3331 * 
		  (h33[i][j] * bp4[k] + hp43[k][i] * s3[j] - (s3[j] * bp4[k] + tg * hp43[k][j]) *
		   (c3 * e23[i] + ctp * bp3[i])) + c3332 * hp334->Get(i,j,k) + c3333 * hp43[k][j] * s3[i] + h33[i][k] * 
		  (tg * s3[j] - ctp * bp3[j] - c3 * e23[j]) + s3[k] * 
		  (s3[i] * s3[j] * s2g + h33[i][j] * tg + bp3[i] * bp3[j] * c2p - hp33[i][j] * ctp + e23[i] * e23[j] * c3334));
	h444->Set(i,j,k, c4411 * (h44[i][j] * bp3[k] + hp43[i][k] * s4[j] - 
				  (s4[j] * bp3[k] + hp43[j][k] * tg) * (e24[i] * c4 + bp4[i] * ctp)) +
		  c4431 * hp443->Get(i,j,k) + c4412 * hp43[j][k] * s4[i] + h44[i][k] * 
		  (tg * s4[j] - ctp * bp4[j] - c4 * e24[j]) + s4[k] * 
		  (s4[i] * s4[j] * s2g + h44[i][j] * tg * bp4[i] * bp4[j] * c2p - hp44[i][j] * ctp + e24[i] * e24[j] * c4441));
	if(i == j) {
	  h333->PlusEq(i,j,k, -(s3[k] * c3335));
	  h444->PlusEq(i,j,k, -(s4[k] * c4442));
	}
      }

  fill3b(3,3, h111);
  fill3b(3,3, h333);
  fill3b(3,3, h444);

  for(k = 0; k < 3; k++)
    for(j = 0; j < k; j++)
      for(i = 0; i < 3; i++) {
	h113->Set(j,k,i, s2g * s3[i] * s1[j] * s1[k] + tg * h31[i][j] * s1[k] + tg * h31[i][k] * s1[j] - c1 *
		  (e21[j] * h31[i][k] + e21[k] * h31[i][j]) - c1112 * ht11[j][k] * s3[i]);
	h113->Set(k,j,i, h113->Get(j,k,i));
	h411->Set(i,j,k, s2g * s4[i] * s1[j] * s1[k] + tg * h41[i][j] * s1[k] + tg * h41[i][k] * s1[j] - c1 *
		  (e21[j] * h41[i][k] + e21[k] * h41[i][j]) - c1112 * ht11[j][k] * s4[i]);
	h411->Set(i,k,j, h411->Get(i,j,k));
	h433->Set(i,j,k, c3331 * (h43[i][j] * bp4[k] + hp44[i][k] * s3[j] + (c4 * e24[i] - ctp*bp4[i]) *
				  (s3[j] * bp4[k] + tg * hp43[k][j]))+c3332 * hp443->Get(i,k,j)+c3333 * hp43[k][j] * s4[i] + h43[i][k] *
		  (tg * s3[j] - ctp * bp3[j] - c3 * e23[j]) + s3[k] * (tg * h43[i][j] + s2g * s4[i] * s3[j] - ctp * hp43[i][j] +
								       c2p * bp4[i] * bp3[j]));
	h433->Set(i,k,j,h433->Get(i,j,k));
      }

  for(i = 0; i < 3; i++)
    for(j = 0; j < i; j++)
      for(k = 0; k < 3; k++) {
	h331->Set(i,j,k, c3331 * h31[i][k] * bp4[j] + c3333*hp43[j][i] * s1[k] + (tg * s3[i] - ctp * bp3[i] - c3 * e23[i]) *
		  h31[j][k] + tg * h31[i][k] * s3[j] + s2g * s3[i] * s3[j] * s1[k]);
	h331->Set(j,i,k, h331->Get(i,j,k));
	h441->Set(i,j,k, c4411 * h41[i][k] * bp3[j] + c4412 * hp43[i][j] * s1[k] + (tg * s4[i] - ctp * bp4[i] - c4 * e24[i]) *
		  h41[j][k] + tg * h41[i][k] * s4[j] + s2g * s4[i] * s4[j] * s1[k]);
	h441->Set(j,i,k, h441->Get(i,j,k));
	h443->Set(i,j,k, c4411 * (h43[j][k] * bp3[i] + hp33[i][k] * s4[j] + (c3 * e23[k] - ctp * bp3[k]) *
				  (bp3[i] * s4[j] + tg * hp43[j][i])) + c4431 * hp334->Get(i,k,j)+c4412 * hp43[j][i] * s3[k] + h43[i][k]*
		  (tg  * s4[j] - ctp * bp4[j] - c4 * e24[j]) + s4[i] * (tg * h43[j][k] + s2g * s4[j] * s3[k] - ctp * hp43[j][k] + 
									c2p * bp4[j] * bp3[k]));
	h443->Set(j,i,k, h443->Get(i,j,k));
      }

  for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++)
      for(k = 0; k < 3; k++) {
	h431->Set(i,j,k, c4311 * (prod->Get(k,j,i) + e21[k] * cp21[i][j] + e24[i] * cp24[j][k] - e24[i] * e21[k] * cp2124[j]) +
		  c4312 * s1[k] * (e24[i] * cp2124[j] - cp21[i][j]) + tg * h41[i][k] * s3[j] + s2g * s4[i] * s3[j] * s1[k] + h31[j][k] * 
		  (tg * s4[i] - ctp * bp4[i]) + c4313 * e24[i] * bp4[j] * s1[k] +c3331* h41[i][k] * bp4[j] +c3333 * hp44[i][j] * s1[k]);
      }

  for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++)
      for(k = 0; k < 3; k++) {
	h112->Set(i,j,k, -(h111->Get(i,j,k) + h113->Get(i,j,k) + h441->Get(k,i,j)));
	h421->Set(i,j,k, -(h411->Get(i,j,k) + h431->Get(i,j,k) + h441->Get(i,j,k)));
	h123->Set(i,j,k, -(h113->Get(i,j,k) + h331->Get(j,k,i) + h431->Get(j,k,i)));
	h332->Set(i,j,k, -(h331->Get(i,j,k) + h333->Get(i,j,k) + h433->Get(k,i,j)));
	h432->Set(i,j,k, -(h431->Get(i,j,k) + h433->Get(i,j,k) + h443->Get(i,k,j)));
	h442->Set(i,j,k, -(h441->Get(i,j,k) + h443->Get(i,j,k) + h444->Get(i,j,k)));
      }
  
  for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++)
      for(k = 0; k < 3; k++) {
	h221->Set(i,j,k, -(h112->Get(i,k,j) + h123->Get(k,j,i) + h421->Get(i,j,k)));
	h223->Set(i,j,k, -(h123->Get(i,j,k) + h332->Get(i,k,j) + h432->Get(i,k,j)));
	h422->Set(i,j,k, -(h421->Get(i,j,k) + h432->Get(i,k,j) + h442->Get(i,k,j)));	
      }

  for(i = 0; i < 3; i++)
    for(j = 0; j < i; j++)
      for(k = 0; k < j; k++) 
	h222->Set(i,j,k, -(h221->Get(i,j,k) + h223->Get(i,j,k) + h442->Get(k,i,j)));

  fill3b(3,3,h222);

  delete[] s1;
  delete[] s2;
  delete[] s3;
  delete[] s4;
  delete[] e21;
  delete[] e23;
  delete[] e24;
  delete[] bp3;
  delete[] bp4;
  
  free_matrix(h11, 3);
  free_matrix(h21, 3);
  free_matrix(h31, 3);
  free_matrix(h41, 3);
  free_matrix(h22, 3);
  free_matrix(h32, 3);
  free_matrix(h42, 3);
  free_matrix(h33, 3);
  free_matrix(h43, 3);
  free_matrix(h44, 3);
  free_matrix(hp33,3);
  free_matrix(hp43,3);
  free_matrix(hp44,3);
  free_matrix(ht11,3);
  free_matrix(cp21,3);
  free_matrix(cp24,3);
  free_matrix(cp24,3);
  delete[] cp2124;
}

void Torsion::Hijs(int disp, int k1, int k2, int k3, int k4, double **h11, 
                   double **h21, double **h31, double **h41, double **h22, 
                   double **h32, double **h42, double **h33, double **h43, 
                   double **h44)
{
	int i, j;
	double *s1, *s2, *s3, *s4, *e21, *e23, *e34;
	double *bp21, *bp22, *bp23, *bp32, *bp34;
	double xx, xy, w1, w2, w3, w4, c1, c2, c3, c4, c5, c6, c7;
	double theta; //W from intder2000.f
	double t21, t23, t34, p2, p3, x1, x2, x3, y1, y2, y3;

	 s1 = new double[3];
   s2 = new double[3];
   s3 = new double[3];
	 s4 = new double[3];
   e21 = new double[3];
   e23 = new double[3];
	 e34 = new double[3];
	 bp21 = new double[3];
   bp22 = new double[3];
   bp23 = new double[3];
   bp32 = new double[3];
	 bp34 = new double[3];

	Torsion::Vect(disp, k1, k2, k3, k4, s1, s2, s3, s4, &theta);
	Stretch::Vect(disp, k1, k2, e21, &t21);
	Stretch::Vect(disp, k3, k2, e23, &t23);
	Stretch::Vect(disp, k4, k3, e34, &t34);
	Bend::Vect(disp, k1, k2, k3, bp21, bp22, bp23, &p2);
	Bend::Vect(disp, k2, k3, k4, bp32, s2, bp34, &p3);

  xx = sin(p2);
  xy = sin(p3);
  xx = t21*xx*xx;
  xy = t34*xy*xy;
  w1 = 1/(t21*xx);
  w2 = 1/(t23*xx);
  w3 = 1/(t34*xy);
  w4 = 1/(t23*xy);
	for(j=0; j<3; j++) {
		for(i=0; i<3; i++) {
     	h11[i][j] = -h11[i][j]*w1;
     	h31[i][j] = h31[i][j]*w2;
    	h44[i][j] = h44[i][j]*w3;
     	h42[i][j] = -h42[i][j]*w4;
  	}
	}
	xx = cos(p2)/sin(p2);
  xy = cos(p3)/sin(p3);
  for(i=0; i<3; i++) {
		w1 = 2*(e21[i]/t21 + bp21[i]*xx);
		w2 = (e23[i]/t23 + 2*bp23[i]*xx);
		w3 = 2*(e34[i]/t34 + bp34[i]*xy);
		w4 = (e23[i]/t23 - 2*bp32[i]*xy);
		for(j=0; j<3; j++) {
			h11[i][j] = h11[i][j] - w1*s1[j];
			h31[i][j]  = h31[i][j] - w2*s1[j];
			h44[i][j] = h44[i][j] - w3*s4[j];
			h42[j][i] = h42[j][i]  + w4*s4[j];
		}
	}
	for(j=0; j<3; j++) {
		for(i=0; i<3; i++) {
			h41[i][j] = 0;
			h21[i][j] = -(h11[i][j] + h31[i][j]);
			h43[i][j] = -(h44[i][j] +  h42[i][j]);
		}
	}

	x1 = t21/t23;
	y1 = t34/t23;
	x2 = cos(p2);
	y2 = sin(p2);
	x3 = cos(p3);
	y3 = sin(p3);
	c1 = x1*x2 - 1;
	c2 = -x3*y1;
	c3 = -x2/t23;
	c4 = -x1*y2;
	c5 = x1*x2/t23;
	c6 = y1*y3;
	c7 = -y1*x3/t23;

	for(i=0; i<3; i++) {
		w1 = c3*e21[i] + c4*bp22[i] + c5*e23[i];
		w2 = c6*bp32[i] + c7*e23[i];
		for(j=0; j<3; j++)
			h22[i][j] = c1*h21[i][j] + c2*h42[j][i] + w1*s1[j] + w2*s4[j];
	}
	
	for(j=0; j<3; j++)
		for(i=0; i<3; i++)
			h32[i][j] = -(h21[j][i] + h22[i][j] + h42[i][j]);
	for(j=0; j<3; j++)
		for(i=0; i<3; i++)
			h33[i][j] = -(h31[i][j] + h32[i][j] + h43[j][i]);

	delete[] s1;
	delete[] s2;
	delete[] s3;
	delete[] s4;
	delete[] e21;
	delete[] e23;
	delete[] e34;
	delete[] bp21;
	delete[] bp22;
	delete[] bp23;
	delete[] bp32;
	delete[] bp34;
}

//Exact translation of HIJKS6 from inter2000.f
void Torsion::Hijks(int disp, int k1, int k2, int k3, int k4, C3DMatrix *h111, C3DMatrix *h112, 
                    C3DMatrix *h221, C3DMatrix *h222, C3DMatrix *h113, C3DMatrix *h123,
                    C3DMatrix *h223, C3DMatrix *h331, C3DMatrix *h332, C3DMatrix *h333,
                    C3DMatrix *h411, C3DMatrix *h421, C3DMatrix *h422, C3DMatrix *h431,
                    C3DMatrix *h432, C3DMatrix *h433, C3DMatrix *h441, C3DMatrix *h442,
                    C3DMatrix *h443, C3DMatrix *h444)
{
	int i, j, k;
	double *v1, *v2, *v3, *v4;
	double *e21, *e23, *e34, *bp21, *bp22, *bp23, *bp32, *bp33, *bp34;
	double **h11, **h21, **h31, **h41, **h22, **h32, **h42, **h33, **h43, **h44;
	double t21, t23, t34, p2, p3, tau, w1, w2, w3, w4, w5, w6;
	double c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16;

	v1 = new double[3];
	v2 = new double[3];
	v3 = new double[3];
	v4 = new double[3];
	e21 = new double[3];
	e23 = new double[3];
	e34 = new double[3];
	bp21 = new double[3];
	bp22 = new double[3];
	bp23 = new double[3];
	bp33 = new double[3];
	bp34 = new double[3];

	h11 = init_matrix(3, 3);
	h21 = init_matrix(3, 3);
	h31 = init_matrix(3, 3);
	h41 = init_matrix(3, 3);
	h22 = init_matrix(3, 3);
	h32 = init_matrix(3, 3);
	h42 = init_matrix(3, 3);
	h33 = init_matrix(3, 3);
	h43 = init_matrix(3, 3);
	h44 = init_matrix(3, 3);

	Torsion::Vect(disp, k1, k2, k3, k4, v1, v2, v3, v4, &tau); //VECT6
	Stretch::Vect(disp, k1, k2, e21, &t21); //VECT1
	Stretch::Vect(disp, k3, k2, e23, &t23);	//VECT1
	Stretch::Vect(disp, k4, k3, e34, &t34);	//VECT1
	Bend::Vect(disp, k1, k2, k3, bp21, bp22, bp23, &p2);	//VECT2
	Bend::Vect(disp, k2, k3, k4, bp32, bp33, bp34, &p3);	//VECT2

	mat1(h32, e23);
	mat1(h21, e21);
	mat1(h43, e34);

	c1 = 1/t21;
	c2 = 1/t34;
	c3 = 1/t23;
	c4 = sin(p2);
	c5 = cos(p2);
	c6 = c5/c4;
	c7 = sin(p3);
	c8 = cos(p3);
	c9 = c8/c7;
	c10 = 1/(c4*c4);
	c11 = 1/(c7*c7);
	c12 = c1*c1;
	c13 = c2*c2;
	c14 = c3*c3;
	c15 = t21*c3;
	c16 = t34*c3;

	w1 = 2*c6;
	w2 = 2*c9;
	w3 = 2*c1;
	w4 = 2*c2;
	w5 = c5*c3;
	w6 = c8*c3;

	for(k=0; k<3; k++){
		h411->Set(1, 1, k, e21[k]*c1 + bp21[k]*c6);
		h411->Set(1, 2, k, e34[k]*c2 + bp34[k]*c9);
		h411->Set(1, 3, k, e23[k]*c3 + bp23[k]*w1);
		h411->Set(2, 1, k, -e23[k]*c3 + bp32[k]*w2);
		h411->Set(2, 2, k, e21[k]*w3 + e23[k]*c3 - bp22[k]*w1);
		h411->Set(2, 3, k, e34[k]*w4 - e23[k]*c3 - bp33[k]*w2);
		h411->Set(3, 1, k, e23[k]*w5 + bp23[k]*c4);
		h411->Set(3, 2, k, -e23[k]*w6 + bp32[k]*c7);
	} 

	for(k=0; k<3; k++) {
		for(i=0; i<3; i++)
			v2[i] = 0;
		v2[k] = 1;
		mat1(h22, v2);
		for(j=0; j<3; j++)
			for(i=0; i<3; i++)
				h421->Set(i, j, k, h22[i][j]);
	}

	w1 = 2*c10*c12;
	w2 = 2*c11*c13;

	for(k=0; k<3; k++) {
		w3 = w1*h411->Get(1, 1, k);
		w4 = w2*h411->Get(1, 2, k);
		for(j=0; j<k; j++)	{
			for(i=0; i<j; i++)	{
				h111->Set(i, j, k, w3*h32[i][j]);
				h444->Set(i, j, k, -w4*h32[i][j]);
			}
		}	
	}

	w1 = c10*c12;
	w2 = c11*c13;
	w3 = w1*c3;
	w4 = w2*c3;

	for(k=0; k<3; k++)	{
		w5 = w1*h411->Get(1, 3, k);
		w6 = w2*h411->Get(2, 3, k);
		for(j=0; j<3; j++)	{
			for(i=0; i<3; i++)	{
				h113->Set(i, j, k, w5*h32[i][j] - w3*h421->Get(i, j, k));
				h442->Set(i, j, k, -w6*h32[i][j] - w4*h421->Get(i, j, k));
			}
		}
	}

	w1 = c1*c3*c10;
	w2 = c2*c3*c11;

	for(k=0; k<3; k++)	{
		w5 = w1*h411->Get(2, 2, k);
		w6 = w2*h411->Get(2, 3, k);
		for(j=0; j<3; j++)	{
			for(i=0; i<3; i++)	{
				h123->Set(i, k, j, -w5*h21[i][j] + w3*h421->Get(i, j, k));
				h432->Set(i, k, j, -w6*h43[i][j] + w4*h421->Get(i, j, k));
			}
		}
	}

	Bend::Hijs(disp, k1, k2, k3, h11, h21, h31, h22, h32, h33);	//HIJS2
	Stretch::Hijs(disp, k1, k2, h44);	//HIJS1
	Stretch::Hijs(disp, k2, k3, h42);	//HIJS1

	w1 = 2*c1;
	w2 = 2*c12;
	
	for(k=0; k<3; k++)	{
		for(i=0; i<k; i++)
			h43[i][k] = 2*(w1*h44[i][k] + c6*h11[i][k] - c10*bp21[i]*bp21[k]);
		h43[k][k] -= w2;
	}

	for(k=0; k<3; k++)
		for(j=0; j<k; j++)
			for(i=0; i<j; i++)
				h111->PlusEq(i, j, k, - v1[j]*h43[i][k]);

	w1 = 2*c6;
	w2 = 2*c10;
	w3 = 2*c3;

	for(k=0; k<3; k++)
		for(i=0; i<3; i++)
			h43[i][k] = h31[k][i]*w1 - bp21[i]*bp23[k]*w2;

	for(k=0; k<3; k++)
		for(j=0; j<3; j++)
			for(i=0; i<3; i++)
				h113->PlusEq(i, j, k, -v1[j]*h43[i][k]);

	for(k=0; k<3; k++)	{
		for(j=0; j<3; j++)
			h43[j][k] = w3*h42[j][k] - w1*h32[k][j] + w2*bp22[j]*bp23[k];
		h43[k][k] -= c14;
	}

	for(k=0; k<3; k++)
		for(j=0; j<3; j++)
			for(i=0; i<3; i++)
				h123->PlusEq(i, j, k, v1[i]*h43[j][k]);

	w1 = c4*c3;
	w2 = c4*c15;
	w3 = c5*c15;
	w4 = w3*c3;
	w5 = c3*c15;

	for(k=0; k<3; k++)
		for(i=0; i<3; i++)
			h43[i][k] = -e21[i]*bp23[k]*w1 + h32[k][i]*w2 + bp22[i]*bp23[k]*w3 - h42[i][k]*w4;

	for(k=0; k<3; k++)
		for(i=0; i<3; i++)
			for(j=0; j<3; j++)
				h223->Set(i, j, k, -v1[j]*(h43[i][k] + e23[i]*w5*h411->Get(3, 1, k)));

	w1 = c3*c4*c15;
	w2 = c5*c14;

	for(k=0; k<3; k++)
		for(i=0; i<3; i++)
			h43[i][k] += -e23[k]*bp22[i]*w1 + e23[k]*w2*(c15*e23[i] - e21[i]);

	for(k=0; k<3; k++)
		for(j=0; j<3; j++)
			for(i=0; i<3; i++)
				h332->Set(i, j, k, v1[j]*h43[k][i]);

	Bend::Hijs(disp, k2, k3, k4, h22, h32, h42, h33, h43, h44); //HIJS2
	Stretch::Hijs(disp, k4, k3, h11); //HIJS1
	Stretch::Hijs(disp, k3, k2, h31); //HIJS1
		
	w1 = 2*c2;
	w2 = 2*c13;

	for(k=0; k<3; k++)	{
		for(i=0; i<k; i++)
			h21[i][k] = 2*(w1*h11[i][k] + c9*h44[i][k] - c11*bp34[i]*bp34[k]);
		h21[k][k] -= w2;
	}

	for(k=0; k<3; k++)
		for(j=0; j<k; j++)
			for(i=0; i<j; i++)
				h444->PlusEq(i, j, k, -v4[j]*h21[i][k]);

	w1 = 2*c9;
	w2 = 2*c11;
	w3 = 2*c3;

	for(k=0; k<3; k++)
		for(i=0; i<3; i++)
			h21[i][k] = w1*h42[i][k] - w2*bp34[i]*bp32[k];

	for(k=0; k<3; k++)
		for(j=0; j<3; j++)
			for(i=0; i<3; i++)
				h442->PlusEq(i, j, k, -v4[j]*h21[i][k]);

	for(k=0; k<3; k++)	{
		for(j=0; j<3; j++)	
			h21[j][k] = w3*h31[j][k] - w1*h32[j][k] + w2*bp33[j]*bp32[k];
		h21[k][k] -= c14;
	}

	for(k=0; k<3; k++)
		for(j=0; j<3; j++)
			for(i=0; i<3; i++)
				h432->PlusEq(i, j, k, v4[i]*h21[j][k]);

	w1 = c7*c3;
	w2 = c7*c16;
	w3 = c8*c16;
	w4 = w3*c3;
	w5 = t34*c14;

	for(k=0; k<3; k++)
		for(i=0; i<3; i++)
			h21[i][k] = -e34[i]*bp32[k]*w1 + h32[i][k]*w2 + bp33[i]*bp32[k]*w3 - h31[i][k]*w4;

	for(k=0; k<3; k++)
		for(i=0; i<3; i++)
			for(j=0; j<3; j++)
				h332->PlusEq(i, j, k, -v4[j]*(h21[i][k] - e23[i]*w5*h411->Get(3, 2, k)));

	w1 = c3*c7*c16;
	w2 = c8*c14;

	for(k=0; k<3; k++)
		for(i=0; i<3; i++)
			h21[i][k] += e23[k]*bp33[i]*w1 + e23[k]*w2*(e34[i] + c16*e23[i]);

	
	for(k=0; k<3; k++)
		for(j=0; j<3; j++)
			for(i=0; i<3; i++)
				h223->PlusEq(i, j, k, v4[j]*h21[k][i]);

	Torsion::Hijs(disp, k1, k2, k3, k4, h11, h21, h31, h41, h22, h32, h42, h33, h43, h44);//HIJS6
	
	for(k=0; k<3; k++)	{
		for(j=0; j<k; j++)	{
			for(i=0; i<j; i++)	{
				h111->PlusEq(i, j, k, -2*h11[j][k]*h411->Get(1, 1, i));
				h444->PlusEq(i, j, k, -2*h44[j][k]*h411->Get(1, 2, i));
			}
		}
	}

	fill3a(3, 3, h111);
	fill3a(3, 3, h444);

	for(i=0; i<3; i++)	{
		w1 = 2*h411->Get(1, 1, i);
		w2 = 2*h411->Get(1, 2, i);
		for(j=0; j<3; j++)	{
			for(k=0; k<3; k++)	{
				h113->PlusEq(i, j, k, -w1*h31[k][j]);
				h442->PlusEq(i, j, k, -w2*h42[j][k]);
				h123->PlusEq(i, j, k, -h21[j][i]*h411->Get(1, 3, k));
				h432->PlusEq(i, j, k, -h43[i][j]*h411->Get(2, 1, k));
			}
		}
	}

	w4 = c5*c15;
	w1 = w4 - 1;
	w2 = c8*c16;
	w3 = w2 - 1;

	for(k=0; k<3; k++)	{
		for(j=0; j<3; j++)	{
			for(i=0; i<3; i++)	{
				h223->PlusEq(i, j, k, w1*h123->Get(j, i, k) - w2*h432->Get(j, k, i)
                              - c15*h21[i][j]*h411->Get(3, 1, k));
				h332->PlusEq(i, j, k, w3*h432->Get(j, i, k) - w4*h123->Get(j, k, i)
                              - c16*h43[j][i]*h411->Get(3, 2, k));
			}
		}
	}
	
	for(k=0; k<3; k++)	{
		for(j=0; j<3; j++)	{
			w1 = c16*(h43[j][k] - c3*v4[j]*e23[k]);
			w2 = c15*(h21[k][j] + c3*v1[j]*e23[k]);
			for(i=0; i<3; i++)	{
				h223->PlusEq(i, j, k, w1*h411->Get(3, 2, i));
				h332->PlusEq(i, j, k, w2*h411->Get(3, 1, i));
			}
		}
	}

	w1 = c5*c3;
	w2 = c4*c15;
	w3 = c5*t21*c14;
	w4 = c8*c3;
	w5 = c7*c16;
	w6 = c8*t34*c14;

	for(k=0; k<3; k++)	{
		h411->Set(1, 1, k, w5*bp33[k] + w6*e23[k] + w4*e34[k]);
		h411->Set(1, 2, k, w2*bp22[k] - w3*e23[k] + w1*e21[k]);
		h411->Set(1, 3, k, -w1*e21[k] - w2*bp22[k] + w3*e23[k]);
		h411->Set(2, 1, k, -w4*e34[k] - w5*bp33[k] - w6*e23[k]);
	}

	for(k=0; k<3; k++)	{
		for(j=0; j<3; j++)	{
			for(i=0; i<3; i++)	{
				h223->PlusEq(i, j, k, h42[j][i]*h411->Get(1, 1, k));
				h332->PlusEq(i, j, k, h31[i][j]*h411->Get(1, 2, k));
			}
		}
	}

	for(k=0; k<3; k++)	{
		for(j=0; j<3; j++)	{
			w1 = h31[k][j] - c3*v1[j]*e23[k];
			w2 = h42[j][k] + c3*v4[j]*e23[k];
			for(i=0; i<3; i++)	{
				h223->PlusEq(i, j, k, w1*h411->Get(1, 3, i));
				h332->PlusEq(i, j, k, w2*h411->Get(2, 1, i));
			}	
		}
	}

	for(k=0; k<3; k++)	{
		for(j=0; j<3; j++)	{
			for(i=0; i<3; i++)	{
				h411->Set(i, j, k, 0);
				h421->Set(i, j, k, 0);
				h431->Set(i, j, k, 0);
				h441->Set(i, j, k, 0);
				h443->Set(j, k, i, -h444->Get(i, j, k) - h442->Get(j, k, i));
				h112->Set(i, j, k, -h111->Get(i, j, k) - h113->Get(i, j, k));
			}
		}	
	}

	for(k=0; k<3; k++)	{
		for(j=0; j<3; j++)	{
			for(i=0; i<3; i++)	{
				h433->Set(k, i, j, -h443->Get(i, k, j) - h432->Get(k, j, i));
				h221->Set(i, j, k, -h112->Get(i, k, j) - h123->Get(k, j, i));
			}
		}
	}
	
	for(k=0; k<3; k++)	{
		for(j=0; j<3; j++)	{
			for(i=0; i<3; i++)	{
				h422->Set(k, i, j, -h432->Get(k, i, j) - h442->Get(i, k, j));
				h331->Set(i, j, k, -h123->Get(k, i, j) - h113->Get(i, k, j));
			}
		}
	}

	for(k=0; k<3; k++)	{
		for(j=0; j<k; j++)	{
			for(i=0; i<j; i++)	{
				h222->Set(i, j, k, -h221->Get(i, j, k) - h223->Get(i, j, k) - h422->Get(k, i, j));
				h333->Set(i, j, k, -h331->Get(i, j, k) - h332->Get(i, j, k) - h433->Get(k, i, j));
			}
		}
	}

	fill3a(3, 3, h222);
	fill3a(3, 3, h333);

	delete[] v1;
	delete[] v2;
	delete[] v3;
	delete[] v4;
	delete[] e21;
  delete[] e23;
  delete[] e34;
  delete[] bp21;
  delete[] bp22;
  delete[] bp23;
  delete[] bp33;
  delete[] bp34;

	free_matrix(h11, 3);
  free_matrix(h21, 3);
  free_matrix(h31, 3);
  free_matrix(h41, 3);
  free_matrix(h22, 3);
  free_matrix(h32, 3);
  free_matrix(h42, 3);
  free_matrix(h33, 3);
  free_matrix(h43, 3);
  free_matrix(h44, 3);
}


//HIJS8 from intder2000.f
void LinearX::Hijs(int disp, int k1, int k2, int k3, int k4, double **h11, double **h21,
                   double **h31, double **h41, double **h22, double **h32, double **h42,
                   double **h33, double **h43, double **h44)
{
	int i, j, k;
	double *e2, *e4, *q1, *q2, *q3;
	double **e22, **e44;
	double **q11, **q21, **q31, **q33, **q22, **q32;
	double w, t32, t34, t123, t;
	C3DMatrix *q111 = new C3DMatrix(3,3,3);
  C3DMatrix *q222 = new C3DMatrix(3,3,3);
  C3DMatrix *q333 = new C3DMatrix(3,3,3);
  C3DMatrix *q444 = new C3DMatrix(3,3,3);
  C3DMatrix *q112 = new C3DMatrix(3,3,3);
  C3DMatrix *q223 = new C3DMatrix(3,3,3);
  C3DMatrix *q331 = new C3DMatrix(3,3,3);
  C3DMatrix *q221 = new C3DMatrix(3,3,3);
  C3DMatrix *q332 = new C3DMatrix(3,3,3);
  C3DMatrix *q113 = new C3DMatrix(3,3,3);
  C3DMatrix *q123 = new C3DMatrix(3,3,3);
	
	e2 = new double[3];
	e4 = new double[3];
	q1 = new double[3];
	q2 = new double[3];
	q3 = new double[3];
	e22 = init_matrix(3, 3);
	e44 = init_matrix(3, 3);
	q11 = init_matrix(3, 3);
	q21 = init_matrix(3, 3);
	q31 = init_matrix(3, 3);
	q33 = init_matrix(3, 3);
	q22 = init_matrix(3, 3);
	q32 = init_matrix(3, 3);

	Stretch::Vect(disp, k2, k3, e2, &t32); //VECT1
	Stretch::Vect(disp, k4, k3, e4, &t34); //VECT1
	Bend::Vect(disp, k1, k2, k3, q1, q2, q3, &t123); //VECT2

	t =	dot_prod(e4, q3);
	w = -t32*t;

	Stretch::Hijs(disp, k4, k3, e44); //HIJS1
	Bend::Hijs(disp, k1, k2, k3, q11, q21, q31, q22, q32, q33); //HIJS2
	Stretch::Hijs(disp, k2, k3, e22); //HIJS1
	Stretch::Hijks(disp, k4, k3, q444); //HIJKS1

	for(j=0; j<3; j++) {
		for(k=0; k<3; k++) {
			h44[j][k] = 0;
			for(i=0; i<3; i++)
				h44[j][k] -= t32*q444->Get(i, j, k)*q3[i];
		}
	}
	
	Bend::Hijks(disp, k1, k2, k3, q111, q112, q113, q123, q221, q222, q223, q331, q332, q333); //HIJKS2
	
	for(k=0; k<3; k++) {
		for(j=0; j<3; j++) {
			h41[j][k] = h42[j][k] = h11[j][k] = h21[j][k] = 0;
			h22[j][k] = w*e22[j][k]/t32;
			for(i=0; i<3; i++) {
				h11[j][k] += -t32*e4[i]*q113->Get(j, k, i);
				h21[j][k] += -e4[i]*(e2[j]*q31[i][k] + t32*q123->Get(k, j, i));
				h22[j][k] += -e4[i]*(e2[j]*q32[i][k] + e2[k]*q32[i][j] + t32*q223->Get(j, k, i));
				h41[j][k] += -t32*e44[i][j]*q31[i][k];
				h42[j][k] += -e44[i][j]*(t32*q32[i][k] + e2[k]*q3[i]);
			}
		}
	}

	for(j=0; j<3; j++) {
		for(k=0; k<3; k++) {
			h31[j][k] = -h11[j][k] - h21[j][k] - h41[j][k];
			h32[j][k] = -h21[k][j] - h22[j][k] - h42[j][k];
			h43[j][k] = -h41[j][k] - h42[j][k] - h44[j][k];
		}
	}

	for(j=0; j<3; j++)
		for(k=0; k<3; k++)
			h33[j][k] = -h31[j][k] - h32[j][k] - h43[k][j];//should this really be h43[k][j]?

	delete q111;
	delete q222;
	delete q333;
	delete q444;
	delete q112;
	delete q223;
	delete q331;
	delete q221;
	delete q332;
	delete q113;
	delete q123;
	delete[] e2;
	delete[] e4;
	delete[] q1;
	delete[] q2;
	delete[] q3;
	free_matrix(e22, 3);
  free_matrix(e44, 3);
  free_matrix(q11, 3);
  free_matrix(q21, 3);
  free_matrix(q31, 3);
  free_matrix(q33, 3);
  free_matrix(q22, 3);
  free_matrix(q32, 3);
}

//HIJKS8 from intder2000.f
void LinearX::Hijks(int disp, int k1, int k2, int k3, int k4, C3DMatrix *h111, C3DMatrix *h112, 
                    C3DMatrix *h221, C3DMatrix *h222, C3DMatrix *h113, C3DMatrix *h123, 
                    C3DMatrix *h223, C3DMatrix *h331, C3DMatrix *h332, C3DMatrix *h333,
                    C3DMatrix *h411, C3DMatrix *h421, C3DMatrix *h422, C3DMatrix *h431,
                    C3DMatrix *h432, C3DMatrix *h433, C3DMatrix *h441, C3DMatrix *h442,
                    C3DMatrix *h443, C3DMatrix *h444)
{
  int i, j, k, l;
  double *q1, *q2, *q3, *q4, *qb, *qc;
  double **q11, **q21, **q31, **q41, **q22, **q32, **q42, **q33, **q43, **q44;
  double r23, r34, phi, w;
  
  C3DMatrix *q111 = new C3DMatrix(3, 3, 3);
  C3DMatrix *q112 = new C3DMatrix(3, 3, 3);
  C3DMatrix *q113 = new C3DMatrix(3, 3, 3);
  C3DMatrix *q123 = new C3DMatrix(3, 3, 3);
  C3DMatrix *q221 = new C3DMatrix(3, 3, 3);
  C3DMatrix *q222 = new C3DMatrix(3, 3, 3);
  C3DMatrix *q223 = new C3DMatrix(3, 3, 3);
  C3DMatrix *q331 = new C3DMatrix(3, 3, 3);
  C3DMatrix *q332 = new C3DMatrix(3, 3, 3);
  C3DMatrix *q444 = new C3DMatrix(3, 3, 3);
	
	C4DMatrix *q1111 = new C4DMatrix(3, 3, 3, 3);
  C4DMatrix *q1112 = new C4DMatrix(3, 3, 3, 3);
  C4DMatrix *q1113 = new C4DMatrix(3, 3, 3, 3);
  C4DMatrix *q1122 = new C4DMatrix(3, 3, 3, 3);
  C4DMatrix *q1123 = new C4DMatrix(3, 3, 3, 3);
  C4DMatrix *q1133 = new C4DMatrix(3, 3, 3, 3);
  C4DMatrix *q1222 = new C4DMatrix(3, 3, 3, 3);
  C4DMatrix *q1223 = new C4DMatrix(3, 3, 3, 3);
  C4DMatrix *q1233 = new C4DMatrix(3, 3, 3, 3);
  C4DMatrix *q1333 = new C4DMatrix(3, 3, 3, 3);
  C4DMatrix *q2222 = new C4DMatrix(3, 3, 3, 3);
  C4DMatrix *q2223 = new C4DMatrix(3, 3, 3, 3);
  C4DMatrix *q2233 = new C4DMatrix(3, 3, 3, 3);
  C4DMatrix *q2333 = new C4DMatrix(3, 3, 3, 3);
  C4DMatrix *q4444 = new C4DMatrix(3, 3, 3, 3);

	q1 = new double[3];
	q2 = new double[3];
	q3 = new double[3];
	q4 = new double[3];
	qb = new double[3];
	qc = new double[3];

	q11 = init_matrix(3, 3);
	q21 = init_matrix(3, 3);
	q31 = init_matrix(3, 3);
	q41 = init_matrix(3, 3);
	q22 = init_matrix(3, 3);
	q32 = init_matrix(3, 3);
	q42 = init_matrix(3, 3);
	q33 = init_matrix(3, 3);
	q43 = init_matrix(3, 3);
	q44 = init_matrix(3, 3);

	Stretch::Vect(disp, k2, k3, qb, &r23); //VECT1
	Stretch::Vect(disp, k4, k3, qc, &r34); //VECT1
	Bend::Vect(disp, k1, k2, k3, q1, q2, q3, &phi); //VECT2
	Bend::Hijs(disp, k1, k2, k3, q11, q21, q31, q22, q32, q33);  //HIJS2
	Stretch::Hijks(disp, k4, k3, q444); //HIJKS1
	Stretch::H4th(disp, k4, k3, q4444); //H4TH1 (don't forget to write this)

	for(i=0; i<3; i++) {
		for(j=0; j<3; j++) {
			for(k=0; k<3; k++) {
				h444->Set(i, j, k, 0);
				h441->Set(i, j, k, 0);
				h442->Set(i, j, k, 0);
				for(l=0; l<3; l++) {
					h444->PlusEq(i, j, k, -r23*q4444->Get(l, i, j, k)*q3[l]);
					h441->PlusEq(i, j, k, -r23*q444->Get(l, i, j)*q31[l][k]);
					h442->PlusEq(i, j, k, -r23*q444->Get(l, i, j)*q32[l][k]);
				}
			}
		}
	}

	Stretch::Hijs(disp, k4, k3, q44); //HIJS2
	Bend::Hijks(disp, k1, k2, k3, q111, q112, q113, q123, q221, q222, q223, q331, q332, q444);//HIJKS2
	Bend::H4th(disp, k1, k2, k3, q1111, q1112, q1113, q1122, q1123, q1133, q1222, q1223, q1233, q1333, q2222, q2223, q2233, q2333, q4444); //H4TH2

	for(i=0; i<3; i++) {
		for(j=0; j<3; j++) {
			for(k=0; k<3; k++) {
				h111->Set(i, j, k, 0);
				h112->Set(i, j, k, 0);
				h221->Set(i, j, k, 0);
				h222->Set(i, j, k, 0);
				h421->Set(i, j, k, 0);
				h422->Set(i, j, k, 0);
				h411->Set(i, j, k, 0);
				for(l=0; l<3; l++) {
					h111->PlusEq(i, j, k, -r23*q1113->Get(i, j, k, l)*qc[l]);
					h112->PlusEq(i, j, k, -r23*q1123->Get(i, j, k, l)*qc[l]);
					h221->PlusEq(i, j, k, -r23*q1223->Get(k, j, i, l)*qc[l]);
					h222->PlusEq(i, j, k, -r23*q2223->Get(i, j, k, l)*qc[l]);
					h421->PlusEq(i, j, k, -r23*q123->Get(k, j, l)*q44[l][i]);
					h422->PlusEq(i, j, k, -r23*q223->Get(j, k, l)*q44[l][i]);
					h411->PlusEq(i, j, k, -r23*q113->Get(j, k, l)*q44[l][i]);
				}
			}	
		}
	}
	
	LinearX::Hijs(disp, k1, k2, k3, k4, q11, q21, q31, q41, q22, q32, q42, q33, q43, q44);//HIJS8
	
	for(i=0; i<3; i++) {
		for(j=0; j<3; j++) {
			for(k=0; k<3; k++) {
				h442->PlusEq(i, j, k, q44[i][j]*qb[k]/r23);
				h421->PlusEq(i, j, k, q41[i][k]*qb[j]/r23);
				h112->PlusEq(i, j, k, q11[i][j]*qb[k]/r23);
				h222->PlusEq(i, j, k, (q22[j][k]*qb[i] + (q22[i][k]*qb[j] + q22[i][j]*qb[k]))/r23);
				h221->PlusEq(i, j, k, (q21[i][k]*qb[j] + q21[j][k]*qb[i])/r23);
				h442->PlusEq(i, j, k, (q42[i][k]*qb[j] + q42[i][j]*qb[k])/r23);
			}
		}
	}
		
	LinearX::Vect(disp, k1, k2, k3, k4, q1, q2, q3, q4, &w); //VECT8
	Stretch::Hijs(disp, k2, k3, q22); //HIJS1
	Stretch::Hijks(disp, k2, k3, q222); //HIJKS1

	for(i=0; i<3; i++) {
		for(j=0; j<3; j++) {
			for(k=0; k<3; k++) {
				h422->PlusEq(i, j, k, (q22[j][k]*q4[i] - 2*qb[j]*qb[k]*q4[i]/r23)/r23);
				h221->PlusEq(i, j, k, (q22[i][j]*q1[k] - 2*qb[i]*qb[j]*q1[k]/r23)/r23);
				h222->PlusEq(i, j, k, (q22[i][j]*q2[k] + q22[j][k]*q2[i] + q22[i][k]*q2[j])/r23
                     + 6*w*qb[i]*qb[j]*qb[k]/(r23*r23*r23)
                     - (qb[i]*qb[j]*q2[k] + qb[j]*qb[k]*q2[i]* + qb[i]*qb[k]*q2[j])*2/(r23*r23)
                     + w*q222->Get(i, j, k)/r23 
                     - 2*w*(q22[i][j]*qb[k] + q22[i][k]*qb[j] + q22[j][k]*qb[i])/(r23*r23));
			}
		}
	}

	for(i=0; i<3; i++) {
		for(j=0; j<3; j++) {
			for(k=0; k<3; k++) {
				h223->Set(i, j, k, -h222->Get(i, j, k) - h221->Get(i, j, k) - h422->Get(k, i, j));
				h113->Set(i, j, k, -h112->Get(i, j, k) - h111->Get(i, j, k) - h411->Get(k, i, j));
				h123->Set(i, j, k, -h112->Get(i, k, j) - h221->Get(k, j, i) - h421->Get(k, j, i));
				h443->Set(i, j, k, -h442->Get(i, j, k) - h441->Get(i, j, k) - h444->Get(i, j, k));
				h431->Set(i, j, k, -h421->Get(i, j, k) - h411->Get(i, j, k) - h441->Get(i, j, k));
				h432->Set(i, j, k, -h422->Get(i, j, k) - h421->Get(i, k, j) - h442->Get(i, j, k));
			}
		}
	}

	for(i=0; i<3; i++) {
		for(j=0; j<3; j++) {
			for(k=0; k<3; k++) {
				h331->Set(i, j, k, -h431->Get(i, j, k) - h123->Get(k, i, j) - h113->Get(i, k, j));
				h332->Set(i, j, k, -h432->Get(i, j, k) - h223->Get(i, k, j) - h123->Get(i, k, j));
				h433->Set(i, j, k, -h431->Get(i, j, k) - h432->Get(i, j, k) - h443->Get(k, i, j));
			}
		}
	}

	for(i=0; i<3; i++)
		for(j=0; j<3; j++)
			for(k=0; k<3; k++)
				h333->Set(i, j, k, -h433->Get(i, j, k) - h331->Get(j, k, i) - h332->Get(j, k, i));

  delete q111;
  delete q112;
  delete q113;
  delete q123;
  delete q221;
  delete q222;
  delete q223;
  delete q331;
  delete q332;
  delete q444;

	delete q1111;
  delete q1112;
  delete q1113;
  delete q1122;
  delete q1123;
	delete q1133;
  delete q1222;
  delete q1223;
  delete q1233;
  delete q1333;
  delete q2222;
  delete q2223;
  delete q2233;
  delete q2333;
  delete q4444;

	delete[] q1;
	delete[] q2;
	delete[] q3;
	delete[] q4;
	delete[] qb;
	delete[] qc;

	free_matrix(q11, 3);
  free_matrix(q21, 3);
  free_matrix(q31, 3);
  free_matrix(q41, 3);
  free_matrix(q22, 3);
  free_matrix(q32, 3);
  free_matrix(q42, 3);
  free_matrix(q33, 3);
  free_matrix(q43, 3);
  free_matrix(q44, 3);
} 


//HIJS9 from intder2000.f
void LinearY::Hijs(int disp, int k1, int k2, int k3, int k4, double **h11, double **h21, 
                   double **h31, double **h41, double **h22, double **h32, double **h42, 
                   double **h33, double **h43, double **h44)
{
	int i, j, k;
	double *e1, *e2, *e3, *e4;
	double **q11, **q12, **q13, **q14, **q22, **q23, **q24, **q33, **q34, **q44;
	double cosy, w, tout;

	e1 = new double[3];
	e2 = new double[3];
	e3 = new double[3];
	e4 = new double[3];
	q11 = init_matrix(3, 3);
	q12 = init_matrix(3, 3);
  q13 = init_matrix(3, 3);
  q14 = init_matrix(3, 3);
  q22 = init_matrix(3, 3);
  q23 = init_matrix(3, 3);
  q24 = init_matrix(3, 3);
  q33 = init_matrix(3, 3);
  q34 = init_matrix(3, 3);
  q44 = init_matrix(3, 3);

	OutOfPlane::Vect(disp, k4, k3, k2, k1, e4, e3, e2, e1, &tout); //VECT5
	w = sin(tout);
	cosy = cos(tout);
	
	OutOfPlane::Hijs(disp, k4, k3, k2, k1, q44, q34, q24, q14, q33, q23, q13, q22, q12, q11);//HIJS7
	for(k=0; k<3; k++) {
		for(j=0; j<3; j++) {
			h22[j][k] = -w*e2[j]*e2[k] - cosy*q22[k][j];
			h32[j][k] = -w*e3[j]*e2[k] - cosy*q23[k][j];
			h42[j][k] = -w*e4[j]*e2[k] - cosy*q24[k][j];
			h33[j][k] = -w*e3[j]*e3[k] - cosy*q33[k][j];
			h43[j][k] = -w*e4[j]*e3[k] - cosy*q34[k][j];
			h44[j][k] = -w*e4[j]*e4[k] - cosy*q44[k][j];
			h41[j][k] = -w*e4[j]*e1[k] - cosy*q14[k][j];
			h31[j][k] = -w*e3[j]*e1[k] - cosy*q13[k][j];
			h21[j][k] = -w*e2[j]*e1[k] - cosy*q12[k][j];
			h11[j][k] = -w*e1[j]*e1[k] - cosy*q11[k][j];
		}
	}

	delete[] e1;
	delete[] e2;
	delete[] e3;
	delete[] e4;
	free_matrix(q11, 3);
  free_matrix(q12, 3);
  free_matrix(q13, 3);
  free_matrix(q14, 3);
  free_matrix(q22, 3);
  free_matrix(q23, 3);
  free_matrix(q24, 3);
  free_matrix(q33, 3);
  free_matrix(q34, 3);
  free_matrix(q44, 3);
}

void LinearY::Hijks(int disp, int k1, int k2, int k3, int k4, C3DMatrix *h111, C3DMatrix *h112, 
                    C3DMatrix *h221, C3DMatrix *h222, C3DMatrix *h113, C3DMatrix *h123, 
                    C3DMatrix *h223, C3DMatrix *h331, C3DMatrix *h332, C3DMatrix *h333,
                    C3DMatrix *h411, C3DMatrix *h421, C3DMatrix *h422, C3DMatrix *h431,
                    C3DMatrix *h432, C3DMatrix *h433, C3DMatrix *h441, C3DMatrix *h442,
                    C3DMatrix *h443, C3DMatrix *h444)
{
  int i, j, k, l;
  double **q11, **q12, **q13, **q14, **q22, **q23, **q24, **q33, **q34, **q44;
  double *q1, *q2, *q3, *q4;
  double cosy, w, tout;

  C3DMatrix *q111 = new C3DMatrix(3, 3, 3);
  C3DMatrix *q112 = new C3DMatrix(3, 3, 3);
  C3DMatrix *q113 = new C3DMatrix(3, 3, 3);
  C3DMatrix *q114 = new C3DMatrix(3, 3, 3);
  C3DMatrix *q122 = new C3DMatrix(3, 3, 3);
  C3DMatrix *q123 = new C3DMatrix(3, 3, 3);
  C3DMatrix *q124 = new C3DMatrix(3, 3, 3);
  C3DMatrix *q133 = new C3DMatrix(3, 3, 3);
  C3DMatrix *q134 = new C3DMatrix(3, 3, 3);
  C3DMatrix *q144 = new C3DMatrix(3, 3, 3);
  C3DMatrix *q222 = new C3DMatrix(3, 3, 3);
  C3DMatrix *q223 = new C3DMatrix(3, 3, 3);
  C3DMatrix *q224 = new C3DMatrix(3, 3, 3);
  C3DMatrix *q332 = new C3DMatrix(3, 3, 3);
  C3DMatrix *q432 = new C3DMatrix(3, 3, 3);
  C3DMatrix *q442 = new C3DMatrix(3, 3, 3);
  C3DMatrix *q333 = new C3DMatrix(3, 3, 3);
  C3DMatrix *q334 = new C3DMatrix(3, 3, 3);
  C3DMatrix *q443 = new C3DMatrix(3, 3, 3);
  C3DMatrix *q444 = new C3DMatrix(3, 3, 3);
  
  q11 = init_matrix(3, 3);
  q12 = init_matrix(3, 3);
  q13 = init_matrix(3, 3);
  q14 = init_matrix(3, 3);
  q22 = init_matrix(3, 3);
  q23 = init_matrix(3, 3);
  q24 = init_matrix(3, 3);
  q33 = init_matrix(3, 3);
  q34 = init_matrix(3, 3);
  q44 = init_matrix(3, 3);

  q1 = new double[3];
  q2 = new double[3];
  q3 = new double[3];
  q4 = new double[3];

  OutOfPlane::Vect(disp, k4, k3, k2, k1, q4, q3, q2, q1, &tout);
  
  w = -(sin(tout));
  cosy = cos(tout);

  OutOfPlane::Hijs(disp, k4, k3, k2, k1, q44, q34, q24, q14, q33, q23, q13, q22, q12, q11);
  OutOfPlane::Hijks(disp, k4, k3, k2, k1, q444, q443, q334, q333, q442, q432, q332, q224, q223, q222, q144, q134, q133, q124, q123,
		    q122, q114, q113, q112, q111);

  for(i=0; i<3; i++) 
    for(j=0; j<3; j++) 
      for(k=0; k<3; k++) {
	h222->Set(i,j,k, cosy * q2[i] * q2[j] * q2[k] - cosy * q222->Get(i,j,k) - w *(q2[i] * q22[j][k] + q2[j] * q22[i][k] + 
										      q2[k] * q22[i][j]));
	h223->Set(i,j,k, cosy * q2[i] * q2[j] * q3[k] - cosy * q223->Get(i,j,k) - w *(q2[i] * q23[j][k] + q2[j] * q23[i][k] + 
										      q3[k] * q22[i][j]));
	h422->Set(i,j,k, cosy * q4[i] * q2[j] * q2[k] - cosy * q224->Get(j,k,i) - w *(q2[k] * q24[j][i] + q2[j] * q24[k][i] + 
										      q4[i] * q22[j][k]));
	h333->Set(i,j,k, cosy * q3[i] * q3[j] * q3[k] - cosy * q333->Get(i,j,k) - w *(q3[k] * q33[j][i] + q3[j] * q33[k][i] + 
										      q3[i] * q33[j][k]));
	h433->Set(i,j,k, cosy * q4[i] * q3[j] * q3[k] - cosy * q334->Get(j,k,i) - w *(q3[k] * q34[j][i] + q3[j] * q34[k][i] + 
										      q4[i] * q33[j][k]));
	h332->Set(i,j,k, cosy * q3[i] * q3[j] * q2[k] - cosy * q332->Get(i,j,k) - w *(q3[i] * q23[k][j] + q3[j] * q23[k][i] + 
										      q2[k] * q33[i][j]));
	h432->Set(i,j,k, cosy * q4[i] * q3[j] * q2[k] - cosy * q432->Get(i,j,k) - w *(q4[i] * q23[k][j] + q3[j] * q24[k][i] + 
										      q2[k] * q34[j][i]));
	h444->Set(i,j,k, cosy * q4[i] * q4[j] * q4[k] - cosy * q444->Get(i,j,k) - w *(q4[i] * q44[k][j] + q4[j] * q44[k][i] + 
										      q4[k] * q44[i][j]));
	h443->Set(i,j,k, cosy * q4[i] * q4[j] * q3[k] - cosy * q443->Get(i,j,k) - w *(q4[i] * q34[k][j] + q4[j] * q34[k][i] + 
										      q3[k] * q44[i][j]));
	h442->Set(i,j,k, cosy * q4[i] * q4[j] * q2[k] - cosy * q442->Get(i,j,k) - w *(q4[i] * q24[k][j] + q4[j] * q24[k][i] + 
										      q2[k] * q44[i][j]));
      }
  
  for(i=0; i<3; i++) 
    for(j=0; j<3; j++) 
      for(k=0; k<3; k++) {
	h221->Set(i,j,k, -(h222->Get(i,j,k) - h223->Get(i,j,k) - h422->Get(k,i,j)));
	h331->Set(i,j,k, -(h332->Get(i,j,k) - h333->Get(i,j,k) - h433->Get(k,i,j)));
	h123->Set(i,j,k, -(h332->Get(i,k,j) - h223->Get(i,j,k) - h432->Get(i,k,j)));
	h441->Set(i,j,k, -(h442->Get(i,j,k) - h443->Get(i,j,k) - h444->Get(i,j,k)));
	h431->Set(i,j,k, -(h432->Get(i,j,k) - h433->Get(i,k,j) - h443->Get(i,k,j)));
	h421->Set(i,j,k, -(h422->Get(i,j,k) - h432->Get(i,k,j) - h442->Get(i,k,j)));
      }

  for(i=0; i<3; i++) 
    for(j=0; j<3; j++) 
      for(k=0; k<3; k++) {
	h112->Set(i,j,k, -(h421->Get(i,k,j) - h123->Get(j,k,i) - h221->Get(i,k,j)));
	h113->Set(i,j,k, -(h431->Get(i,k,j) - h331->Get(i,k,j) - h123->Get(j,i,k)));
	h411->Set(i,j,k, -(h441->Get(i,j,k) - h431->Get(i,j,k) - h421->Get(i,j,k)));
      }
  
  for(i=0; i<3; i++) 
    for(j=0; j<3; j++) 
      for(k=0; k<3; k++) 
	h111->Set(i,j,k, -(h411->Get(k,i,j) - h113->Get(i,j,k) - h112->Get(i,j,k)));

  free_matrix(q11, 3);
  free_matrix(q12, 3);
  free_matrix(q13, 3);
  free_matrix(q14, 3);
  free_matrix(q22, 3);
  free_matrix(q23, 3);
  free_matrix(q24, 3);
  free_matrix(q33, 3);
  free_matrix(q34, 3);
  free_matrix(q44, 3);

  delete[] q1;
  delete[] q2;
  delete[] q3;
  delete[] q4;
  
}
