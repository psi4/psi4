/*!
    \defgroup MINTS libmints: Integral library
    \ingroup MINTS
*/
#include <stdexcept>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "extents.h"
#include "gshell.h"
#include "vector.h"
#include "matrix.h"
#include "basisset.h"

using namespace std;
using namespace psi;
using namespace boost;

namespace psi {

Extents::Extents(boost::shared_ptr<BasisSet> bs1, boost::shared_ptr<BasisSet> bs2) : 
    basis1_(bs1), basis2_(bs2)
{
}
Extents::~Extents()
{
}
boost::shared_ptr<Vector> Extents::x1()
{
    boost::shared_ptr<Vector> x(new Vector("xc", basis1_->nbf()));
    double* xp = x->pointer();   

    for (int i = 0; i < basis1_->nbf(); i++) {
        xp[i] = basis1_->shell(basis1_->function_to_shell(i))->center()[0];
    } 
 
    return x;
}
boost::shared_ptr<Vector> Extents::y1()
{
    boost::shared_ptr<Vector> y(new Vector("yc", basis1_->nbf()));
    double* yp = y->pointer();   

    for (int i = 0; i < basis1_->nbf(); i++) {
        yp[i] = basis1_->shell(basis1_->function_to_shell(i))->center()[1];
    } 
 
    return y;
}
boost::shared_ptr<Vector> Extents::z1()
{
    boost::shared_ptr<Vector> z(new Vector("zc", basis1_->nbf()));
    double* zp = z->pointer();   

    for (int i = 0; i < basis1_->nbf(); i++) {
        zp[i] = basis1_->shell(basis1_->function_to_shell(i))->center()[2];
    } 
 
    return z;
}
boost::shared_ptr<Vector> Extents::x2()
{
    boost::shared_ptr<Vector> x(new Vector("xc", basis2_->nbf()));
    double* xp = x->pointer();   

    for (int i = 0; i < basis2_->nbf(); i++) {
        xp[i] = basis2_->shell(basis2_->function_to_shell(i))->center()[0];
    } 
 
    return x;
}
boost::shared_ptr<Vector> Extents::y2()
{
    boost::shared_ptr<Vector> y(new Vector("yc", basis2_->nbf()));
    double* yp = y->pointer();   

    for (int i = 0; i < basis2_->nbf(); i++) {
        yp[i] = basis2_->shell(basis2_->function_to_shell(i))->center()[1];
    } 
 
    return y;
}
boost::shared_ptr<Vector> Extents::z2()
{
    boost::shared_ptr<Vector> z(new Vector("zc", basis2_->nbf()));
    double* zp = z->pointer();   

    for (int i = 0; i < basis2_->nbf(); i++) {
        zp[i] = basis2_->shell(basis2_->function_to_shell(i))->center()[2];
    } 
 
    return z;
}
boost::shared_ptr<Matrix> Extents::x12()
{
    boost::shared_ptr<Matrix> x(new Matrix("xc", basis1_->nbf(), basis2_->nbf()));
    double** xp = x->pointer();

    for (int P = 0; P < basis1_->nshell(); P++) {
        for (int Q = 0; Q < basis2_->nshell(); Q++) {
            boost::shared_ptr<GaussianShell> sh1 = basis1_->shell(P);    
            boost::shared_ptr<GaussianShell> sh2 = basis1_->shell(Q);
            int Pstart = sh1->function_index();    
            int Qstart = sh2->function_index();   
            int nP = sh1->nfunction();            
            int nQ = sh2->nfunction();            

            double X = 0.0;
            double X1 = sh1->center()[0];
            double X2 = sh2->center()[0];
    
            for (int K1 = 0; K1 < sh1->nprimitive(); K1++) {
                for (int K2 = 0; K2 < sh2->nprimitive(); K2++) {
                    X += sh1->coef(K1) * sh2->coef(K2) * (sh1->exp(K1) * X1 + sh2->exp(K2) * X2) 
                        / (sh1->exp(K1) + sh2->exp(K2)); 
                }
            }
 
            for (int oP = 0; oP < nP; oP++) {
                for (int oQ = 0; oQ < nQ; oQ++) {
                    xp[Pstart + oP][Qstart + oQ] = X;
                }
            }
        }
    } 

    return x;
}
boost::shared_ptr<Matrix> Extents::y12()
{
    boost::shared_ptr<Matrix> y(new Matrix("yc", basis1_->nbf(), basis2_->nbf()));
    double** yp = y->pointer();

    for (int P = 0; P < basis1_->nshell(); P++) {
        for (int Q = 0; Q < basis2_->nshell(); Q++) {
            boost::shared_ptr<GaussianShell> sh1 = basis1_->shell(P);    
            boost::shared_ptr<GaussianShell> sh2 = basis1_->shell(Q);
            int Pstart = sh1->function_index();    
            int Qstart = sh2->function_index();   
            int nP = sh1->nfunction();            
            int nQ = sh2->nfunction();            

            double X = 0.0;
            double X1 = sh1->center()[1];
            double X2 = sh2->center()[1];
    
            for (int K1 = 0; K1 < sh1->nprimitive(); K1++) {
                for (int K2 = 0; K2 < sh2->nprimitive(); K2++) {
                    X += sh1->coef(K1) * sh2->coef(K2) * (sh1->exp(K1) * X1 + sh2->exp(K2) * X2) 
                        / (sh1->exp(K1) + sh2->exp(K2)); 
                }
            }
 
            for (int oP = 0; oP < nP; oP++) {
                for (int oQ = 0; oQ < nQ; oQ++) {
                    yp[Pstart + oP][Qstart + oQ] = X;
                }
            }
        }
    } 

    return y;
}
boost::shared_ptr<Matrix> Extents::z12()
{
    boost::shared_ptr<Matrix> z(new Matrix("zc", basis1_->nbf(), basis2_->nbf()));
    double** zp = z->pointer();

    for (int P = 0; P < basis1_->nshell(); P++) {
        for (int Q = 0; Q < basis2_->nshell(); Q++) {
            boost::shared_ptr<GaussianShell> sh1 = basis1_->shell(P);    
            boost::shared_ptr<GaussianShell> sh2 = basis1_->shell(Q);
            int Pstart = sh1->function_index();    
            int Qstart = sh2->function_index();   
            int nP = sh1->nfunction();            
            int nQ = sh2->nfunction();            

            double X = 0.0;
            double X1 = sh1->center()[2];
            double X2 = sh2->center()[2];
    
            for (int K1 = 0; K1 < sh1->nprimitive(); K1++) {
                for (int K2 = 0; K2 < sh2->nprimitive(); K2++) {
                    X += sh1->coef(K1) * sh2->coef(K2) * (sh1->exp(K1) * X1 + sh2->exp(K2) * X2) 
                        / (sh1->exp(K1) + sh2->exp(K2)); 
                }
            }
 
            for (int oP = 0; oP < nP; oP++) {
                for (int oQ = 0; oQ < nQ; oQ++) {
                    zp[Pstart + oP][Qstart + oQ] = X;
                }
            }
        }
    } 

    return z;
}
boost::shared_ptr<Vector> Extents::extents1(double epsilon)
{ 
    boost::shared_ptr<Vector> e(new Vector("Extents", basis1_->nbf()));
    double* ep = e->pointer();

    double* cart = new double[basis1_->max_function_per_shell()];

    for (int P = 0; P < basis1_->nshell(); P++) {
        boost::shared_ptr<GaussianShell> sh1 = basis1_->shell(P);    
        int Pstart = sh1->function_index();   
        int nP = sh1->nfunction();            
       
        // Spherical guess (geometrically weighted)
        double numerator = 0.0;
        double denominator = 0.0;
        for (int K = 0; K < sh1->nprimitive(); K++) {
            numerator   += sh1->coef(K) * log(sh1->exp(K));
            denominator += sh1->coef(K);
        }   
        double zeta = exp(numerator / denominator);
        double R0 = sqrt(-log(epsilon) / (2.0 * zeta));

        for (int oP = 0; oP < nP; oP++) {
            int index = oP + Pstart;
            ep[index] = R0;    

            

        } 
    }

    delete[] cart;

    return e;
}
boost::shared_ptr<Vector> Extents::extents2(double epsilon)
{ 
    boost::shared_ptr<Vector> e(new Vector("Extents", basis1_->nbf()));
    return e;
}
boost::shared_ptr<Matrix> Extents::extents12(double epsilon)
{ 
    boost::shared_ptr<Matrix> e(new Matrix("Extents", basis1_->nbf(), basis2_->nbf()));
    return e;
}

}
