#ifndef QUAD_H
#define QUAD_H

#include <psi4-dec.h>

namespace psi {

class Quadrature {
    protected:
        /// Number of points in this quadrature rule
        int npoints_;
        /// Current index in the quadrature rule 
        int index_;
        /// Set of points (arbitrary domain)
        double* t_;
        /// Set of weights (Cartesian points, no spherical r^2)
        double* w_;
    public:
        /// Constructor, allocates memory
        Quadrature(int npoints);
        /// Destructor, frees memory
        virtual ~Quadrature();

        /// Prints the Quadrature rule
        virtual void print(FILE* out = outfile) = 0;

        /// Get the current quadrature weight
        double getWeight() { return w_[index_]; }

        /// Get the current quadrature point
        double getPoint() { return t_[index_]; }

        /// Move to the next point
        void nextPoint() { index_++; }

        /// Reset the quadrature
        void reset() { index_ = 0; }

        /// See if the quadrature is complete
        bool isDone() { return index_ >= npoints_; }
};

class ChebyshevIIQuadrature : public Quadrature {
    protected:

    public:
        /**!
        * Constructor, develops quadrature rule
        * \param npoints Number of points.
        * \param t0 Center point of distribution (1.0 is default)
        */
        ChebyshevIIQuadrature(int npoints, double t0 = 1.0);
        /**!
        * Destructor, does nothing
        */
        ~ChebyshevIIQuadrature() {}

        /// Prints the Quadrature rule
        void print(FILE* out = outfile);

};

}
#endif

