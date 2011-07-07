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

  $Id: funcplot.h 2208 2011-03-04 16:33:00Z justus.c79@gmail.com $
*/
#ifndef MADNESS_MRA_FUNCPLOT_H__INCLUDED
#define MADNESS_MRA_FUNCPLOT_H__INCLUDED

/*!

  \file mra/funcplot.h
  \brief Defines/implements plotting interface for functions
  \ingroup funcplot

  @{
 */

namespace madness {
    /// Writes an OpenDX format file with a cube/slice of points on a uniform grid

    /// Collective operation but only process 0 writes the file.  By convention OpenDX
    /// files end in ".dx" but this choice is up to the user.  The binary format is
    /// more compact and vastly faster to both write and load but is not as portable.
    ///
    /// Now follow some brief tips about how to look at files inside OpenDX.
    ///
    /// To view a 1D function \c file-selector-->import-->plot-->image.
    ///
    /// To view a 2D function as a colored plane \c file-selector-->import-->autocolor-->image.
    ///
    /// To view a 2D function as a 3D surface \c file-selector-->import-->rubbersheet-->image.
    ///
    /// To view a 3D function as an isosurface \c file-selector-->import-->isosurface-->image.
    ///
    /// To select the real/imaginary/absolute value of a complex number insert a compute
    /// element after the import.
    template <typename T, std::size_t NDIM>
    void plotdx(const Function<T,NDIM>& f,
                const char* filename,
                const Tensor<double>& cell = FunctionDefaults<NDIM>::get_cell(),
                const std::vector<long>& npt = std::vector<long>(NDIM,201L),
                bool binary=true);


    /// Writes the header information of a VTK file for plotting in an external
    /// post-processing package (such as Paraview)
    //
    /// @param world World communicator
    /// @param filename String containing the filename to export to
    /// @param plotlo Vector of double values indicating the minimum coordinate to plot to in each dimension
    /// @param plothi Vector of double values indicating the maximum coordinate to plot to in each dimension
    /// @param npt Vector of long integers indicating the number of points to plot in each dimension
    /// @param binary (optional) Boolean indicating whether to print in binary

    /// The VTK routines are also designed for SERIAL data, parallel coming...
    ///
    /// This header is templated by the dimension of the data.
    ///
    /// To plot with the plotvtk_* routines:
    ///    plotvtk_begin(...)
    ///    plotvtk_data(...)
    ///    plotvtk_data(...) ...
    ///    plotvtk_end(...)
    ///
    /// NOTE: Paraview expects the structured mesh points in a particular
    /// order, which is why the LowDimIndexIterator is used...
    template<std::size_t NDIM>
    void plotvtk_begin(World &world, const char *filename,
        const Vector<double, NDIM> &plotlo, const Vector<double, NDIM> &plothi,
        const Vector<long, NDIM> &npt, bool binary = false) {

        PROFILE_FUNC;
        MADNESS_ASSERT(NDIM>=1 && NDIM<=3); // how do we plot data in more than 3-D?

        Tensor<double> cell(NDIM, 2);
        std::size_t i;
        for(i = 0; i < NDIM; ++i) {
            cell(i, 0) = plotlo[i];
            cell(i, 1) = plothi[i];
        }

        FILE *f=0;
        if(world.rank() == 0) {
            f = fopen(filename, "w");
            if(!f)
                MADNESS_EXCEPTION("plotvtk: failed to open the plot file", 0);

            fprintf(f, "<VTKFile type=\"StructuredGrid\" version=\"0.1\"" \
                " byte_order=\"LittleEndian\" compressor=\"" \
                "vtkZLibDataCompressor\">\n");
            fprintf(f, "  <StructuredGrid WholeExtent=\"");
            for(i = 0; i < NDIM; ++i)
                fprintf(f, "0 %ld ", npt[i]-1);
            for(; i < 3; ++i)
                fprintf(f, "0 0 ");
            fprintf(f, "\">\n");
            fprintf(f, "    <Piece Extent=\"");
            for(i = 0; i < NDIM; ++i)
                fprintf(f, "0 %ld ", npt[i]-1);
            for(; i < 3; ++i)
                fprintf(f, "0 0 ");
            fprintf(f, "\">\n");
            fprintf(f, "      <Points>\n");
            fprintf(f, "        <DataArray NumberOfComponents=\"3\" " \
                "type=\"Float32\" format=\"ascii\">\n");

            Vector<double, NDIM> space;
            for(i = 0; i < NDIM; ++i) {
                if(npt[i] == 1)
                    space[i] = 0.0;
                else
                    space[i] = (cell(i, 1) - cell(i, 0)) / (npt[i] - 1);
            }

            // go through the grid
            for(LowDimIndexIterator it(npt); it; ++it) {
                for(i = 0; i < NDIM; ++i)
                    fprintf(f, "%f ", plotlo[i] + it[i]*space[i]);
                for(; i < 3; ++i)
                    fprintf(f, "0.0 ");
                fprintf(f, "\n");
            }

            fprintf(f, "        </DataArray>\n");
            fprintf(f, "      </Points>\n");
            fprintf(f, "      <PointData>\n");
            fclose(f);
        }
        world.gop.fence();
    }

    /// Generic VTK data writer. Specific type instances of this function are defined for
    /// both real and complex valued functions.
    //
    /// @param function Function (real or complex) that we wish to export the data of
    /// @param fieldname A string containing the name we wish to refer to this field as in the exported data
    /// @param world World communicator
    /// @param filename String containing the filename to export to
    /// @param plotlo Vector of double values indicating the minimum coordinate to plot to in each dimension
    /// @param plothi Vector of double values indicating the maximum coordinate to plot to in each dimension
    /// @param npt Vector of long integers indicating the number of points to plot in each dimension
    /// @param binary (optional) Boolean indicating whether to print in binary

    /// This templated function won't do anything except print a warning
    /// message.  Specialized versions of this function should be used.
    template<typename T, std::size_t NDIM>
    void plotvtk_data(const T &function, const char *fieldname, World &world,
        const char *filename, const Vector<double, NDIM> &plotlo,
        const Vector<double, NDIM> &plothi, const Vector<long, NDIM> &npt,
        bool binary = false) {

        MADNESS_EXCEPTION("plotvtk only supports madness::functions", 0);
    }

    /// VTK data writer for real-valued (not complex) madness::functions.

    /// Set plot_refine=true to get a plot of the refinement levels of
    /// the given function.
    template<typename T, std::size_t NDIM>
    void plotvtk_data(const Function<T, NDIM> &function, const char *fieldname,
        World &world, const char *filename, const Vector<double, NDIM> &plotlo,
        const Vector<double, NDIM> &plothi, const Vector<long, NDIM> &npt,
        bool binary = false, bool plot_refine = false) {

        PROFILE_FUNC;
        MADNESS_ASSERT(NDIM>=1 && NDIM<=3); // no plotting high-D functions, yet...

        Tensor<double> cell(NDIM, 2);
        std::size_t i;
        for(i = 0; i < NDIM; ++i) {
            cell(i, 0) = plotlo[i];
            cell(i, 1) = plothi[i];
        }
        std::vector<long> numpt(NDIM);
        for(i = 0; i < NDIM; ++i)
            numpt[i] = npt[i];

        world.gop.barrier();

        function.verify();
        FILE *f = 0;
        if(world.rank() == 0) {
            f = fopen(filename, "a");
            if(!f)
                MADNESS_EXCEPTION("plotvtk: failed to open the plot file", 0);

            fprintf(f, "        <DataArray Name=\"%s\" format=\"ascii\" " \
                "type=\"Float32\" NumberOfComponents=\"1\">\n", fieldname);
        }

        world.gop.fence();
        Tensor<T> tmpr = function.eval_cube(cell, numpt, plot_refine);
        world.gop.fence();

        if(world.rank() == 0) {
            for(LowDimIndexIterator it(numpt); it; ++it) {
                fprintf(f, "%.6e\n", tmpr(*it));
            }
            fprintf(f, "        </DataArray>\n");
            fclose(f);
        }
        world.gop.fence();
    }

    /// VTK data writer for complex-valued madness::functions.

    /// The complex-value is written as two reals (a vector from VTK's
    /// perspective.  The first (X) component is the real part and the second
    /// (Y) component is the imaginary part.
    /// Set plot_refine=true to get a plot of the refinement levels of
    /// the given function.
    template<typename T, std::size_t NDIM>
    void plotvtk_data(const Function<std::complex<T>, NDIM> &function,
        const char *fieldname, World &world, const char *filename,
        const Vector<double, NDIM> &plotlo, const Vector<double, NDIM> &plothi,
        const Vector<long, NDIM> &npt, bool binary = false,
        bool plot_refine = false) {

        // this is the same as plotvtk_data for real functions, except the
        // real and imaginary parts are printed on the same line (needed
        // to change NumberOfComponents in the XML tag)

        PROFILE_FUNC;
        MADNESS_ASSERT(NDIM>=1 && NDIM<=3); // no plotting high-D functions, yet...

        Tensor<double> cell(NDIM, 2);
        std::size_t i;
        for(i = 0; i < NDIM; ++i) {
            cell(i, 0) = plotlo[i];
            cell(i, 1) = plothi[i];
        }
        std::vector<long> numpt(NDIM);
        for(i = 0; i < NDIM; ++i)
            numpt[i] = npt[i];

        world.gop.barrier();

        function.verify();
        FILE *f = 0;
        if(world.rank() == 0) {
            f = fopen(filename, "a");
            if(!f)
                MADNESS_EXCEPTION("plotvtk: failed to open the plot file", 0);

            fprintf(f, "        <DataArray Name=\"%s\" format=\"ascii\" " \
                "type=\"Float32\" NumberOfComponents=\"2\">\n", fieldname);
        }

        world.gop.fence();
        Tensor<std::complex<T> > tmpr = function.eval_cube(cell, numpt,
                                                           plot_refine);
        world.gop.fence();

        if(world.rank() == 0) {
            for(LowDimIndexIterator it(numpt); it; ++it) {
                fprintf(f, "%.6e %.6e\n", real(tmpr(*it)), imag(tmpr(*it)));
            }
            fprintf(f, "        </DataArray>\n");
            fclose(f);
        }
        world.gop.fence();
    }

    /// Writes the footer information of a VTK file for plotting in an external
    /// post-processing package (such as Paraview)
    //
    /// @param world World communicator
    /// @param filename Name of VTK file
    /// @param binary (Optional) Boolean indicating whether to print in binary
    template<std::size_t NDIM>
    void plotvtk_end(World &world, const char *filename, bool binary = false) {
        PROFILE_FUNC;
        MADNESS_ASSERT(NDIM>=1 && NDIM<=3);

        FILE *f = 0;
        if(world.rank() == 0) {
            f = fopen(filename, "a");
            if(!f)
                MADNESS_EXCEPTION("plotvtk: failed to open the plot file", 0);

            fprintf(f, "      </PointData>\n");
            fprintf(f, "      <CellData>\n");
            fprintf(f, "      </CellData>\n");
            fprintf(f, "    </Piece>\n");
            fprintf(f, "  </StructuredGrid>\n");
            fprintf(f, "</VTKFile>\n");
            fclose(f);
        }
        world.gop.fence();
    }

    namespace detail {
        inline unsigned short htons_x(unsigned short a) {
            return (a>>8) | (a<<8);
        }
    }

    /// Writes a Povray DF3 format file with a cube of points on a uniform grid

    /// Collective operation but only process 0 writes the file.  By convention Povray
    /// files end in ".df3" but this choice is up to the user.  The dynamic range of
    /// function values is mapped onto [0,1] and values stored in 16-bit fixed precision.
    template <typename T>
    static void plotpovray(const Function<T,3>& function,
                           const char* filename,
                           const Tensor<double>& cell = FunctionDefaults<3>::get_cell(),
                           const std::vector<long>& npt = std::vector<long>(3,201L))
    {
        using detail::htons_x;

        MADNESS_ASSERT(npt.size() == 3);
        unsigned short dims[3] = {htons_x(npt[0]),htons_x(npt[1]),htons_x(npt[2])};

        World& world = const_cast< Function<T,3>& >(function).world();
        FILE *f=0;
        if (world.rank() == 0) {
            f = fopen(filename, "w");
            if (!f) MADNESS_EXCEPTION("plotdx: failed to open the plot file", 0);
            fwrite((void*) dims, sizeof(short), 3, f);
        }
        Tensor<T> r = function.eval_cube(cell, npt);
        if (world.rank() == 0) {
            double rmax = r.max();
            double rmin = r.min();
            double rrange = rmax + rmin;
            double rmean = rrange*0.5;
            double fac = 65535.0/rrange;

            printf("plot_povray: %s: min=%.2e(0.0) mean=%.2e(0.5) max=%.2e(1.0) range=%.2e\n",
                   filename,rmin,rmean,rmax,rrange);

            std::vector<unsigned short> d(npt[0]);
            for (unsigned int i2=0; i2<npt[2]; ++i2) {
                for (unsigned int i1=0; i1<npt[1]; ++i1) {
                    for (unsigned int i0=0; i0<npt[0]; ++i0) {
                        d[i0] = (unsigned short)(htons_x((unsigned short)(fac*(r(i0,i1,i2) - rmin))));
                        //printf("%d\n",htons_x(d[i0]));
                    }
                    fwrite((void*) &d[0], sizeof(short), npt[0], f);
                }
            }

            fclose(f);
        }
    }

    static inline void plot_line_print_value(FILE* f, double_complex v) {
        fprintf(f, "    %.14e %.14e   ", real(v), imag(v));
    }

    static inline void plot_line_print_value(FILE* f, double v) {
        fprintf(f, " %.14e", v);
    }

    /// Generates ASCII file tabulating f(r) at npoints along line r=lo,...,hi

    /// The ordinate is distance from lo
    template <typename T, std::size_t NDIM>
    void plot_line(const char* filename, int npt, const Vector<double,NDIM>& lo, const Vector<double,NDIM>& hi,
                   const Function<T,NDIM>& f) {
        typedef Vector<double,NDIM> coordT;
        coordT h = (hi - lo)*(1.0/(npt-1));

        double sum = 0.0;
        for (std::size_t i=0; i<NDIM; ++i) sum += h[i]*h[i];
        sum = sqrt(sum);

        World& world = f.world();
        f.reconstruct();
        if (world.rank() == 0) {
            FILE* file = fopen(filename,"w");
            for (int i=0; i<npt; ++i) {
                coordT r = lo + h*double(i);
                fprintf(file, "%.14e ", i*sum);
                plot_line_print_value(file, f.eval(r));
                fprintf(file,"\n");
            }
            fclose(file);
        }
        world.gop.fence();
    }

    /// Generates ASCII file tabulating f(r) and g(r) at npoints along line r=lo,...,hi

    /// The ordinate is distance from lo
    template <typename T, typename U, std::size_t NDIM>
    void plot_line(const char* filename, int npt, const Vector<double,NDIM>& lo, const Vector<double,NDIM>& hi,
                   const Function<T,NDIM>& f, const Function<U,NDIM>& g) {
        typedef Vector<double,NDIM> coordT;
        coordT h = (hi - lo)*(1.0/(npt-1));

        double sum = 0.0;
        for (std::size_t i=0; i<NDIM; ++i) sum += h[i]*h[i];
        sum = sqrt(sum);

        World& world = f.world();
        f.reconstruct();
        g.reconstruct();
        if (world.rank() == 0) {
            FILE* file = fopen(filename,"w");
            for (int i=0; i<npt; ++i) {
                coordT r = lo + h*double(i);
                fprintf(file, "%.14e ", i*sum);
                plot_line_print_value(file, f.eval(r));
                plot_line_print_value(file, g.eval(r));
                fprintf(file,"\n");
            }
            fclose(file);
        }
        world.gop.fence();
    }


    /// Generates ASCII file tabulating f(r), g(r), and a(r) at npoints along line r=lo,...,hi

    /// The ordinate is distance from lo
    template <typename T, typename U, typename V, std::size_t NDIM>
    void plot_line(const char* filename, int npt, const Vector<double,NDIM>& lo, const Vector<double,NDIM>& hi,
                   const Function<T,NDIM>& f, const Function<U,NDIM>& g, const Function<V,NDIM>& a) {
        typedef Vector<double,NDIM> coordT;
        coordT h = (hi - lo)*(1.0/(npt-1));

        double sum = 0.0;
        for (std::size_t i=0; i<NDIM; ++i) sum += h[i]*h[i];
        sum = sqrt(sum);

        World& world = f.world();
        f.reconstruct();
        g.reconstruct();
        a.reconstruct();
        if (world.rank() == 0) {
            FILE* file = fopen(filename,"w");
            for (int i=0; i<npt; ++i) {
                coordT r = lo + h*double(i);
                fprintf(file, "%.14e ", i*sum);
                plot_line_print_value(file, f.eval(r));
                plot_line_print_value(file, g.eval(r));
                plot_line_print_value(file, a.eval(r));
                fprintf(file,"\n");
            }
            fclose(file);
        }
        world.gop.fence();
    }

    /// Generates ASCII file tabulating f(r), g(r), a(r), b(r) at npoints along line r=lo,...,hi

    /// The ordinate is distance from lo
    template <typename T, typename U, typename V, typename W, std::size_t NDIM>
    void plot_line(const char* filename, int npt, const Vector<double,NDIM>& lo, const Vector<double,NDIM>& hi,
                   const Function<T,NDIM>& f, const Function<U,NDIM>& g, const Function<V,NDIM>& a, const Function<W,NDIM>& b) {
        typedef Vector<double,NDIM> coordT;
        coordT h = (hi - lo)*(1.0/(npt-1));

        double sum = 0.0;
        for (std::size_t i=0; i<NDIM; ++i) sum += h[i]*h[i];
        sum = sqrt(sum);

        World& world = f.world();
        f.reconstruct();
        g.reconstruct();
        a.reconstruct();
        b.reconstruct();
        if (world.rank() == 0) {
            FILE* file = fopen(filename,"w");
            for (int i=0; i<npt; ++i) {
                coordT r = lo + h*double(i);
                fprintf(file, "%.14e ", i*sum);
                plot_line_print_value(file, f.eval(r));
                plot_line_print_value(file, g.eval(r));
                plot_line_print_value(file, a.eval(r));
                plot_line_print_value(file, b.eval(r));
                fprintf(file,"\n");
            }
            fclose(file);
        }
        world.gop.fence();
    }
}

/* @} */
#endif // MADNESS_MRA_FUNCPLOT_H__INCLUDED
