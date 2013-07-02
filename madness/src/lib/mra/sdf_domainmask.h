/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680

  $Id: sdf_shape.h 1792 2010-01-26 12:58:18Z rjharrison $
*/

/*!
  \file mra/sdf_domainmask.h
  \brief Defines abstract interfaces and concrete classes signed distance
         functions and domain masks.

  Interfaces for a signed distance function (sdf) and a domain mask are
  made.  The original conception of these interfaces was for doing
  shapes and interior boundary conditions in MADNESS; however, the
  interfaces were abstracted for other applications.

  The domain mask builds on a sdf such that \f$0 <=  \mbox{mask(sdf)} <= 1\f$
  for smooth switching between domains.  That said, a mask needs a sdf.

  A general-purpose functor is given that combines a sdf and mask
  to produce MADNESS functions for various entities including
      - The domain mask (trace is the volume)
      - The gradient of mask (related to the surface)
      - The surface (trace is the surface area)
      - The normal derivative of the surface

  \ingroup mrabcint
*/

#ifndef MADNESS_MRA_SDF_DOMAINMASK_H__INCLUDED
#define MADNESS_MRA_SDF_DOMAINMASK_H__INCLUDED

#include <mra/mra.h>

namespace madness {

    /** \brief The interface for a signed distance function (sdf).

        A class implementing this interface will need to provide the sdf
        given a point in coordinate space, and also the gradient of the
        sdf at this point.

        \c NDIM is the dimensionality of the coordinate space.

        \ingroup mrabcint */
    template <std::size_t NDIM>
    class SignedDFInterface {
    public:
        /** \brief Returns the signed distance from the surface,
               - Positive is ``inside''
               - Negative is ``outside''

            \param[in] x The coordinate
            \return The signed distance */
        virtual double sdf(const Vector<double,NDIM>& x) const = 0;

        /** \brief Returns the gradient of the signed distance from the
                   surface (i.e., \c dsdf(x)/dx[i] )

            \param[in] x The coordinate
            \return The gradient */
        virtual Vector<double,NDIM> grad_sdf(const Vector<double,NDIM>& x)
            const = 0;

        virtual ~SignedDFInterface() {}
    };

    /** \brief The interface for masking functions defined by signed distance
               functions.

        This interface was initially conceived for using shapes in MADNESS
        to specify internal boundary conditions.  The signed distance function
        defines the shape, and this mask allows calculations of volumes,
        surfaces, etc.

        \ingroup mrabcint */
    class DomainMaskInterface {
    public:
        /** \brief Returns the characteristic mask function
                - 1 in interior
                - 0 in exterior
                - 1/2 on boundary

            \param[in] d The signed distance from the surface
            \return The mask or characteristic function */
        virtual double mask(double d) const = 0;

        /** \brief Returns the derivative of the characteristic mask function
                   with respect to the distance (from the surface)

            \param[in] d The signed normal distance from the surface
            \return Derivative of the characteristic mask function */
        virtual double dmask(double d) const = 0;

        /** \brief Returns the value of the normalized surface layer function

            Normalized means the volume integral of this function should
            converge to the surface area in the limit of a either an infinitely
            thin surface or zero curvature.  The function thus acts as a
            ``delta function'' located at the boundary.

            \param[in] d The signed normal distance from the surface
            \return Normalized surface layer function */
        virtual double surface(double d) const = 0;

        /** \brief Returns the derivative of the normalized surface layer
                   function

            \param[in] d The signed normal distance from the surface
            \return Derivative of the normalized surface layer function */
        virtual double dsurface(double d) const = 0;

        virtual ~DomainMaskInterface() {}
    };

    /** \brief Framework for combining a signed distance function (sdf)
               with a domain mask to produce MADNESS functions.

        This interface provides functor functionality to produce MADNESS
        functions for
           - the domain mask (given the sdf)
           - the derivative of the domain mask
           - the surface (given the sdf)
           - the normal derivative of the surface layer

        The functor defaults to the domain mask; however, member functions
        can toggle between the other options.

        \ingroup mrabcint */
    template <std::size_t NDIM>
    class DomainMaskSDFFunctor : public FunctionFunctorInterface<double,NDIM> {

    private:
        /// \brief Bury the default constructor
        DomainMaskSDFFunctor() {}

    private: // protected may be better if this becomes highly inherited
        /// The domain mask to use
        std::shared_ptr<DomainMaskInterface> mask;

        /// The signed distance function
        std::shared_ptr<SignedDFInterface<NDIM> > sdf;

        int mswitch; ///< Which masking function to use (mask, surface, etc.)

    public:
        // switch values
        static const int MASK; ///< Use the \c mask() function in \c mask
        static const int MASK_COMPLEMENT; ///< Get the complement of \c mask()
        static const int DMASK; ///< Use the \c dmask() function in \c mask
        static const int SURFACE; ///< Use the \c surface() function in \c mask
        static const int DSURFACE; ///< Use the \c dsurface() function in \c mask

        /** \brief Constructor for mask/sdf functor

            \param mask Pointer to the domain mask
            \param sdf Pointer to the signed distance function */
        DomainMaskSDFFunctor(std::shared_ptr<DomainMaskInterface> mask,
                std::shared_ptr<SignedDFInterface<NDIM> > sdf)
            : mask(mask), sdf(sdf), mswitch(MASK)
        {}

        /** \brief Constructor for mask/sdf function specifying the desired
                   function (mask, surface, etc.)

            \param mask Pointer to the domain mask
            \param sdf Pointer to the signed distance function
            \param[in] _mswitch Which function to use (MASK, DMASK, SURFACE,
                       DSURFACE, or MASK_COMPLEMENT) */
        DomainMaskSDFFunctor(std::shared_ptr<DomainMaskInterface> mask,
                std::shared_ptr<SignedDFInterface<NDIM> > sdf, int _mswitch)
            : mask(mask), sdf(sdf), mswitch(MASK) {
            if(_mswitch == MASK || _mswitch == DMASK || _mswitch == SURFACE ||
               _mswitch == DSURFACE || _mswitch == MASK_COMPLEMENT) {
                mswitch = _mswitch;
            }
            else {
                error("Unrecognized function option in DomainMaskSDFFunctor" \
                      "::DomainMaskSDFFunctor()");
            }
        }

        /** \brief Uses the functor interface to make a MADNESS function

            \param x Point to compute value
            \return Value of the desired function */
        double operator()(const Vector<double,NDIM>& x) const {
            if(mswitch == DomainMaskSDFFunctor<NDIM>::MASK)
                return mask->mask(sdf->sdf(x));
            else if(mswitch == DomainMaskSDFFunctor<NDIM>::MASK_COMPLEMENT)
                return 1.0 - mask->mask(sdf->sdf(x));
            else if(mswitch == DomainMaskSDFFunctor<NDIM>::DMASK)
                return mask->dmask(sdf->sdf(x));
            else if(mswitch == DomainMaskSDFFunctor<NDIM>::SURFACE)
                return mask->surface(sdf->sdf(x));
            else if(mswitch == DomainMaskSDFFunctor<NDIM>::DSURFACE)
                return mask->dsurface(sdf->sdf(x));
            else {
                error("Unknown function from DomainMaskInterface in " \
                      "DomainMaskSDFFunctor::operator()");
                return 0.0;
            }
        }

        /** \brief Toggles which function from DomainMaskInterface to use when
                   making the MADNESS function.

            \param _mswitch The function to use (should be MASK, DMASK,
                   SURFACE, DSURFACE, or MASK_COMPLEMENT) */
        void setMaskFunction(int _mswitch) {
            if(_mswitch == MASK || _mswitch == DMASK || _mswitch == SURFACE ||
               _mswitch == DSURFACE || _mswitch == MASK_COMPLEMENT) {
                mswitch = _mswitch;
            }
            else {
                error("Unrecognized function option in DomainMaskSDFFunctor" \
                      "::setMaskFunction()");
            }
        }

        virtual ~DomainMaskSDFFunctor() {}
    };

    template<std::size_t NDIM>
    const int DomainMaskSDFFunctor<NDIM>::MASK = 1;

    template<std::size_t NDIM>
    const int DomainMaskSDFFunctor<NDIM>::MASK_COMPLEMENT = 2;

    template<std::size_t NDIM>
    const int DomainMaskSDFFunctor<NDIM>::DMASK = 3;

    template<std::size_t NDIM>
    const int DomainMaskSDFFunctor<NDIM>::SURFACE = 4;

    template<std::size_t NDIM>
    const int DomainMaskSDFFunctor<NDIM>::DSURFACE = 5;

    /** \brief Provides the Li-Lowengrub-Ratz-Voight (LLRV) domain mask
               characteristic functions.

        \ingroup mrabcint

        See X. Li, J. Lowengrub, A. R&auml;tz, and A. Voight, ``Solving PDEs in
        Complex Geometries: A Diffuse Domain Approach,'' Commun. Math. Sci., 7,
        p81-107, 2009.

        Given a signed distance, this class implements in the domain mask
        and surface functions from the above reference.  For the domain mask,

        \f[ \varphi(d) = \frac{1}{2}\left( 1 - \tanh\left(
                         \frac{3d}{\varepsilon} \right) \right) \f]

        where \f$d\f$ is the signed distance.  The normalized surface function
        is

        \f[ B(\varphi) = \frac{36}{\varepsilon} \varphi^2 (1-\varphi)^2. \f]

        The constant \f$36/\varepsilon\f$ is chosen to fulfill

        \f[ \int_{-\infty}^\infty B(s) \, ds = 1 \f]

        This class assumes the domain mask is uniformly 0 or 1 outside
        signed distances \f$ |8 \epsilon| \f$ since the switching function
        becomes 0/1 to machine precision at these levels.  Specifically,
        for this function the parameter \f$ \epsilon \f$ is an effective
        measure of the full width of the surface layer since

        \f[ \int_{-\epsilon/2}^{\epsilon/2} B(s) \, ds \doteq 0.987 \f]

        and

        \f[ \int_{-\epsilon}^{\epsilon} B(s) \, ds \doteq 0.999963 \f] */
    class LLRVDomainMask : public DomainMaskInterface {
    private:
        LLRVDomainMask() : epsilon(0.0) {} ///< Forbidden

    protected:
        const double epsilon; ///< The width of the transition region

    public:
        /** \brief Constructor for the domain mask

            \param[in] epsilon The effective width of the surface */
        LLRVDomainMask(double epsilon)
            : epsilon(epsilon)
        {}

        /** \brief Value of characteristic function at normal distance d from
                   the surface

            \param[in] d The signed distance.  Negative is ``inside,''
                         positive is ``outside.''
            \return The domain mask */
        double mask(double d) const {
            if (d > 8.0*epsilon) {
                return 0.0; // we're safely outside
            }
            else if (d < -8.0*epsilon) {
                return 1.0; // inside
            }
            else {
                return 0.5 * (1.0 - tanh(3.0 * d / epsilon));
            }
        }

        /** \brief Derivative of characteristic function with respect to the
                   normal distance

            \param[in] d The signed distance
            \return The derivative */
        double dmask(double d) const {
            if (fabs(d) > 8.0*epsilon) {
                return 0.0; // we're safely outside or inside
            }
            else {
                double tanh3d = tanh(3.0*d/epsilon);
                return 1.5*(tanh3d*tanh3d - 1.0) / epsilon;
            }
        }

        /** \brief Value of surface function at distance d normal to surface

            \param[in] d The signed distance
            \return The value of the surface function */
        double surface(double d) const {
            double phi = mask(d);
            double phic = 1.0 - phi;
            return 36.0*phi*phi*phic*phic/epsilon;
        }

        /** \brief Value of d(surface)/ddistance

            \param[in] d The signed distance
            \return The derivative of the surface function */
        double dsurface(double d) const {
            double phi = mask(d);
            double dphi = dmask(d);
            return 72.0*phi*(1.0-phi)*dphi*(1.0 - 2.0*phi)/epsilon;
        }

        virtual ~LLRVDomainMask() {}
    };

    /** \brief Use a Gaussian for the surface function and the corresponding erf
               for the domain mask. */
    class GaussianDomainMask : public DomainMaskInterface {
    private:
        GaussianDomainMask() : epsilon(0.0) {} ///< Forbidden

    protected:
        const double epsilon; ///< The width of the transition region

    public:
        /** \brief Constructor for the domain mask

            \param[in] epsilon The effective width of the surface */
        GaussianDomainMask(double epsilon)
            : epsilon(epsilon)
        {}

        /** \brief Value of characteristic function at normal distance d from
                   the surface

            \param[in] d The signed distance.  Negative is ``inside,''
                         positive is ``outside.''
            \return The domain mask */
        double mask(double d) const {
            if (d > 8.0*epsilon) {
                return 0.0; // we're safely outside
            }
            else if (d < -8.0*epsilon) {
                return 1.0; // inside
            }
            else {
                return (1.0 - erf(d / (sqrt(2.0) * epsilon))) * 0.5;
            }
        }

        /** \brief Derivative of characteristic function with respect to the
                   normal distance

            \param[in] d The signed distance
            \return The derivative */
        double dmask(double d) const {
            if (fabs(d) > 8.0*epsilon) {
                return 0.0; // we're safely outside or inside
            }
            else {
                return -exp(-d*d*0.5/(epsilon*epsilon)) / (sqrt(2.0*constants::pi)
                    * epsilon);
            }
        }

        /** \brief Value of surface function at distance d normal to surface

            \param[in] d The signed distance
            \return The value of the surface function */
        double surface(double d) const {
            return exp(-d*d*0.5/(epsilon*epsilon)) / (sqrt(2.0*constants::pi)
                * epsilon);
        }

        /** \brief Value of d(surface)/ddistance

            \param[in] d The signed distance
            \return The derivative of the surface function */
        double dsurface(double d) const {
            return -exp(-d*d*0.5/(epsilon*epsilon)) * d / (sqrt(2.0*constants::pi)
                * epsilon*epsilon*epsilon);
        }

        virtual ~GaussianDomainMask() {}
    };

} // end of madness namespace

#endif // MADNESS_MRA_SDF_DOMAINMASK_H__INCLUDED
