#ifndef MOLECULAR_MASK_H
#define MOLECULAR_MASK_H

#include <mra/mra.h>
#include <constants.h>
#include <cmath>
#include <vector>

// Distance between two points in 3D
inline double distance(const madness::coord_3d& a, const madness::coord_3d& b) {
    double x = a[0] - b[0];
    double y = a[1] - b[1];
    double z = a[2] - b[2];
    return sqrt(x*x + y*y + z*z);
}

// Basic functionality for the mask
class MolecularMaskBase {
protected:
    const double sigma;
    const std::vector<double> atomic_radii;
    const std::vector<madness::coord_3d> atomic_coords;
    const int natom;

    // signed distance function for point r relative to sphere radius R at center
    double sdf(const madness::coord_3d& r, const madness::coord_3d& center, double R) const {
        return distance(r,center) - R;
    }

    // gradient of the signed distance function
    madness::coord_3d grad_sdf(const madness::coord_3d& r, const madness::coord_3d& center) const {
        return (r - center)*(1.0/distance(r,center));
    }

    // Mask or characteristic function (argument s is the signed distance)
    double mask(double s) const {
        if (s > 6.0) return 0.0;
        else if (s < -6.0) return 1.0;
        else return 0.5*erfc(s);
    }

    // Complement of the mask or characteristic function (argument s is the signed distance)
    double cmask(double s) const {
        return mask(-s);
    }

    // Derivative of the mask w.r.t. s
    double dmask(double s) const {
        const double fac = 1.0/sqrt(madness::constants::pi);
        if (fabs(s) > 6.0) return 0.0;
        return -exp(-s*s)*fac;
    }

    // Mask or characteristic function for atom i
    double atomic_mask(const madness::coord_3d& r, unsigned int i) const {
        double s = sdf(r, atomic_coords[i], atomic_radii[i]);
        return mask(s/sigma);
    }

    // Complement of the mask or characteristic function for atom i
    // (we use this directly to avoid numerical cancellation)
    double atomic_cmask(const madness::coord_3d& r, unsigned int i) const {
        double s = sdf(r, atomic_coords[i], atomic_radii[i]);
        return cmask(s/sigma);
    }

    // Gradient of the atomic mask
    madness::coord_3d grad_atomic_mask(const madness::coord_3d& r, unsigned int i) const {
        double s = sdf(r, atomic_coords[i], atomic_radii[i]);
        return grad_sdf(r,atomic_coords[i])*(dmask(s/sigma)/sigma);
    }

    // Gradient of the molecular mask
    madness::coord_3d gradient(const madness::coord_3d& r) const {
        // precompute the atomic masks
        std::vector<double> m(natom);
        double value = 1.0;
        for (int i=0; i<natom; i++) {
            m[i] = atomic_cmask(r,i);
            value *= m[i];
        }
        
        // return 0.0 if not in the surface
        if (value<1e-12 || (1.0-value)<1e-12) return madness::coord_3d(0.0);
        
        madness::coord_3d grad(0.0);
        for (int i=0; i<natom; i++) {
            if (m[i] > 1e-12) 
                grad += grad_atomic_mask(r,i)*(value/m[i]);
        }
        return grad;
    }

public:
    MolecularMaskBase(double sigma, 
                      const std::vector<double> atomic_radii,
                      const std::vector<madness::coord_3d> atomic_coords) 
        : sigma(sigma)
        , atomic_radii(atomic_radii)
        , atomic_coords(atomic_coords)
        , natom(atomic_coords.size())
    {
        MADNESS_ASSERT(atomic_radii.size() == atomic_coords.size());
    }
};

// This functor is one inside the molecule, 1/2 on the surface, and zero
// exterior to the molecule.
class MolecularVolumeMask : private MolecularMaskBase
                          , public madness::FunctionFunctorInterface<double,3> {
public:
    MolecularVolumeMask(double sigma, 
                        const std::vector<double> atomic_radii,
                        const std::vector<madness::coord_3d> atomic_coords) 
        : MolecularMaskBase(sigma, atomic_radii, atomic_coords)
    {}

    virtual double operator()(const madness::coord_3d& r) const {
        double value = 1.0;
        for (int i=0; i<natom; i++) {
            value *= atomic_cmask(r,i);
        }
        return 1.0 - value;
    }
};

// This functor is zero inside the molecule, 1/2 on the surface, and one
// exterior to the molecule.
class MolecularVolumeComplementMask : private MolecularMaskBase
                                    , public madness::FunctionFunctorInterface<double,3> {
public:
    MolecularVolumeComplementMask(double sigma, 
                                  const std::vector<double> atomic_radii,
                                  const std::vector<madness::coord_3d> atomic_coords) 
        : MolecularMaskBase(sigma, atomic_radii, atomic_coords)
    {}
    
    virtual double operator()(const madness::coord_3d& r) const {
        double value = 1.0;
        for (int i=0; i<natom; i++) {
            value *= atomic_cmask(r,i);
        }
        return value;
    }
};

/// Switches between \a positive values \c Vint and \c Vext with special log derivative


/// Switches between \a positive values \c Vint on the interior and \c Vext on the interior.
/// It has value 
/// \f[
///  V_(r,\sigma) = \exp \left( \log V_{int} C(r,\sigma) + \log V_{ext} (1 - C(r,\sigma)  \right)
///  = V_{ext} \exp \left( \log \frac{V_{int}}{V_{ext}} C(r,\sigma) \right)
///  = V_{int} \exp \left( \log \frac{V_{ext}}{V_{int}} \overline{C}(r,\sigma) \right)
/// /f]
/// where \f$ C(r,\sigma) \f$ is the regular volume mask provided by MolecularVolumeMask,
/// and \f$ \overline{C} \f$ is its complement.
/// Its log-derivative is precisely located in the surface with value
/// \f[ 
/// \nabla \log V = \frac{\nabla V}{V} = \log \frac{V_{int}}{V_{ext}} \nabla C
/// \f] 
/// with \f$ \nabla C \f$ already being computed by MolecularVolumeMaskGrad.
/// The advantage of this is that if \f$ \| V_{ext} - V_{int} \| \f$
/// is big, the log derivative of the regular volume mask (\f$ C \f$) is
/// displaced from the surface (perhaps by multiple values of \f$
/// \sigma \f$).  This leads to slow convergence w.r.t \f$ \sigma \f$ and 
/// potential inaccuracies depending how the numerical representation is computed. 
/// The surface charge
/// in dielectric problems is controlled by \f$ \nabla \log \epsilon
/// \f$ and hence the dielectric should employ this form of the switch.
class MolecularVolumeExponentialSwitch : public madness::FunctionFunctorInterface<double,3> {
    const MolecularVolumeComplementMask cmask;
    const double Vint, Vext, fac;
public:
    MolecularVolumeExponentialSwitch(double sigma, 
                                     double Vint, 
                                     double Vext,
                                     const std::vector<double> atomic_radii,
                                     const std::vector<madness::coord_3d> atomic_coords) 
        : cmask(sigma, atomic_radii, atomic_coords)
        , Vint(Vint)
        , Vext(Vext)
        , fac(log(Vext/Vint))
    {
        if (Vint <= 0 || Vext <= 0) throw "Only works for positive values";
    }
    
    virtual double operator()(const madness::coord_3d& r) const {
        double c = cmask(r);
        if (c == 0.0) return Vint;
        else if (c == 1.0) return Vext;
        else return Vint * exp(fac * c); 
    }
};

/// Computes the reciprocal of MolecularVolumeExponentialSwitch
class MolecularVolumeExponentialSwitchReciprocal : public madness::FunctionFunctorInterface<double,3> {
    const MolecularVolumeExponentialSwitch s;
public:
    MolecularVolumeExponentialSwitchReciprocal(double sigma, 
                                               double Vint, 
                                               double Vext,
                                               const std::vector<double> atomic_radii,
                                               const std::vector<madness::coord_3d> atomic_coords) 
        : s(sigma, Vint, Vext, atomic_radii, atomic_coords)
    {}
    
    virtual double operator()(const madness::coord_3d& r) const {
        return 1.0/s(r);
    }
};


// This functor is a shell that limits to a delta function in the
// molecular surface and integrates to the molecular surface area.
class MolecularSurface : private MolecularMaskBase
                       , public madness::FunctionFunctorInterface<double,3> {
public:
    MolecularSurface(double sigma, 
                     const std::vector<double> atomic_radii,
                     const std::vector<madness::coord_3d> atomic_coords) 
        : MolecularMaskBase(sigma, atomic_radii, atomic_coords)
    {}

    virtual double operator()(const madness::coord_3d& r) const {
        madness::coord_3d grad = gradient(r);
        return sqrt(grad[0]*grad[0] + grad[1]*grad[1] + grad[2]*grad[2]);
    }
};

// Evaluates component i (0=x, 1=y, 2=z) of the gradient of the mask
class MolecularVolumeMaskGrad : private MolecularMaskBase
                              , public madness::FunctionFunctorInterface<double,3> {
    const int i;
public:
    MolecularVolumeMaskGrad(double sigma, 
                     const std::vector<double> atomic_radii,
                     const std::vector<madness::coord_3d> atomic_coords,
                     int i) 
        : MolecularMaskBase(sigma, atomic_radii, atomic_coords)
        , i(i)
    {}

    virtual double operator()(const madness::coord_3d& r) const {
        madness::coord_3d grad = gradient(r);
        return grad[i];
    }
};

/// Returns the requested component of the derivative of the log of MolecularVolumeExponentialSwitch
class MolecularVolumeExponentialSwitchLogGrad : public madness::FunctionFunctorInterface<double,3> {
    const MolecularVolumeMaskGrad g;
    const double fac;
public:
    MolecularVolumeExponentialSwitchLogGrad(double sigma, 
                                            double Vint, 
                                            double Vext,
                                            const std::vector<double> atomic_radii,
                                            const std::vector<madness::coord_3d> atomic_coords,
                                            int i) 
        : g(sigma, atomic_radii, atomic_coords, i)
        , fac(log(Vint/Vext))
    {
        if (Vint <= 0 || Vext <= 0) throw "Only works for positive values";
    }
    
    virtual double operator()(const madness::coord_3d& r) const {
        return fac * g(r);
    }
};

#endif
