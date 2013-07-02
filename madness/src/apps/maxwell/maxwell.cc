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

  $Id$
*/
#define WORLD_INSTANTIATE_STATIC_TEMPLATES

#include <string>
#include <mra/mra.h>
#include <mra/sdf_shape_3D.h>
#include <constants.h>
#include <misc/interpolation_1d.h>
#include <linalg/gmres.h>
#include "envelopedpulse.h"
#include "frequencyhandler.h"
#include "complex_fun.h"
#include "operator-maxwell.h"

using namespace madness;

static const double thresh = 1e-4;
static const double thresh1 = thresh*0.1;
static const int k_wavelet = 6;

/// Function to setup the frequency-dependent complex dielectric constant
/// Needs the frequency, the Ag permittivity, the Ag conductivity, the
/// mask function for distinguishing Ag from vacuum (1 for Ag, 0 for vacuum),
/// and the functions to put:
///   1) epshat = I omega eps - sigma
///   2) epshatI = I omega (eps - eps0) - sigma  -> Used with the incident field
///   3) lnepshat = ln(epshat) ... easier to compute with the shape function
///
/// These three functions are assumed to be constant in a region, and to vary
/// as the shape function between the two.  Calculating lnepshat in this
/// manner, and not as a pointwise ln of epshat, produces smoother results.
inline void calcEpshat(const double omega, const double ag_eps,
	const double ag_sigma, const Function<double, 3> &ag_shape,
	Function<complexd, 3> &epshat, Function<complexd, 3> &epshatI,
	Function<complexd, 3> &lnepshat);

/// Function to solve for a scattered field...
/// Needs the frequency, the permittivity/conductivity/shape of Ag,
/// the Frequency handler for storage, the total domain mask to zero out
/// a function near the boundary, an interior mask function to damp out
/// spurious oscillations in derivatives near the boundaries, and the Green's
/// function.
void calcScatField(const double omega, const double ag_eps,
	const double ag_sigma, const Function<double, 3> &ag_shape,
	FrequencyHandler &fh, const Function<double, 3> &box_mask,
	const Function<double, 3> &grad_mask,
	const SeparatedConvolution<double, 3> &G);

//*****************************************************************************
// C++ function to solve for an electric field
int main(int argc, char **argv) {
	Vector<double, 3> pt;
	int n;

	madness::initialize(argc, argv);
	World world(MPI::COMM_WORLD);
	startup(world,argc,argv);

	// get the name of the file to use from the command-line arguments
	if(argc < 2) {
		fprintf(stderr, "Error: must provide the base name of the input file " \
			"structure.\n");
		return -1;
	}

	// get the parameters and data from the initial parameter file
	EnvelopedPulse pulse;
	string sfile(argv[1]);
	if(!pulse.read_param_file(sfile))
		return -1;

#if 0
// artificially set the frequencies for debugging
pulse.nfreq = 1;
pulse.freqs[0] = 41;
#endif

	// Function defaults
	FunctionDefaults<3>::set_k(k_wavelet);
	FunctionDefaults<3>::set_cubic_cell(pulse.sim_min, pulse.sim_max);
	FunctionDefaults<3>::set_thresh(thresh);
	FunctionDefaults<3>::set_max_refine_level(5);

	{
		Tensor<int> bc(3, 2);
		bc(_, 0) = 0; // Dirichlet in all directions
		bc(_, 1) = 0;
		FunctionDefaults<3>::set_bc(bc);
	}

	// calculate the dielectric and conductivity of the silver particle for
	// this driving wavelength (this is relative to eps0)
	double ag_eps, ag_sigma;
	{
		// use a Drude model.  This one is particularly suited for 325-375 nm.
		const double epsD = 8.926;
		const double omega_D = 0.425741;
		const double gamma_D = 7.460110e-3;

		ag_eps = omega_D*omega_D / (pulse.omegabar*pulse.omegabar +
			gamma_D*gamma_D);
		ag_sigma = ag_eps * gamma_D;
		ag_eps = epsD - ag_eps;
	}

	printf("\nSphere dielectric: %.6e au\n", ag_eps*pulse.eps0);
	printf("Sphere conductivity: %.6e au\n", ag_sigma*pulse.eps0);

	// create a World object for my local group
	World group(MPI::COMM_WORLD);
	char str[80];
	ParallelOutputArchive outarchive;
	ParallelInputArchive inarchive;

	// a domain mask for smoothing out the boundaries. --> No structures (i.e.
	// domains of interest) should be outside this box
	pt[0] = pt[1] = pt[2] = 100.0 /* nm */ / pulse.nm_per_au;
	Function<double, 3> box_mask = FunctionFactory<double,3>(group).functor(
	    std::shared_ptr< FunctionFunctorInterface<double,3> >(new SDF_Sphere<double>(
		8.0 /* nm */ / pulse.nm_per_au, thresh, 90.0 /* nm */ / pulse.nm_per_au,
		pt)));
	box_mask.truncate();

	Function<double, 3> grad_mask = FunctionFactory<double,3>(group).functor(
	    std::shared_ptr< FunctionFunctorInterface<double,3> >(new SDF_Sphere<double>(
		8.0 /* nm */ / pulse.nm_per_au, thresh, 70.0 /* nm */ / pulse.nm_per_au,
		pt)));
	grad_mask.truncate();

	// create the various masks for the spherical nanoparticle
	Function<double, 3> sphere = FunctionFactory<double,3>(group).functor(
	    std::shared_ptr< FunctionFunctorInterface<double,3> >(new SDF_Sphere<double>(
		5.0 /* nm */ / pulse.nm_per_au, thresh, 20.0 /* nm */ / pulse.nm_per_au,
		pt)));
	sphere.truncate();

	printf("\n");

	// setup the Green's function
	SeparatedConvolution<double, 3> G =
		BSHOperator<double, 3>(group, 0.0, k_wavelet, thresh, thresh1);

	// NOTE: ELECTRIC FIELDS ARE STORED IN MICRO-ATOMIC UNITS SINCE 1 au OF
	// ELECTRIC FIELD IS REALLY, REALLY LARGE

	// create an array to hold the various frequencies
	// positive/negative frequencies are conjugates of each other
	FrequencyHandler **fh = new FrequencyHandler*[pulse.nfreq];
	for(n = 0; n < pulse.nfreq; ++n)
		fh[n] = 0;

	// go through the various frequencies and do the calculations!
	for(n = 0; n < pulse.nfreq; ++n) {
		fh[n] = new FrequencyHandler(pulse.freqs[n], pulse, group);

		if(!fh[n]->read_file()) {
			fprintf(stderr, "Error reading file for frequency %d. ABORT.\n",
				pulse.freqs[n]);
			return -1;
		}

		// see if this frequency has been calculated and archived
		sprintf(str, "%s.scat.%d", pulse.basefile.c_str(), pulse.freqs[n]);
		if(ParallelInputArchive::exists(group, str)) {
			// archive is there... read the function
			printf("Reading frequency %d of %d: n = %d\n   Archive file: %s\n",
				n+1, pulse.nfreq, pulse.freqs[n], str);
			inarchive.open(group, str);
			inarchive & fh[n]->scat[0] & fh[n]->scat[1] & fh[n]->scat[2];
			inarchive.close();
		}
		else {
			printf("Calculating frequency %d of %d: n = %d\n", n+1, pulse.nfreq,
				pulse.freqs[n]);

			// perform the calculations for this frequency
			calcScatField(2.0*constants::pi*pulse.freqs[n]/pulse.T, ag_eps,
				ag_sigma, sphere, *fh[n], box_mask, grad_mask, G);

			outarchive.open(group, str);
			outarchive & fh[n]->scat[0] & fh[n]->scat[1] & fh[n]->scat[2];
			outarchive.close();
		}
	}

	// Fourier transform back, once for each time
	vecfunc scat;
	Vector<Function<double, 3>, 3> field;
	Function<double, 3> scat_real;
	complexd sinfreq;
	double t, dt = 0.01 /* fs */ / pulse.fs_per_au;
	std::shared_ptr<TimeIncident> time_inc(new TimeIncident(pulse));

	for(int ti = 0; ti < 9; ++ti) {
		t = ti * dt + 0.5 / pulse.fs_per_au;
		time_inc->setTime(t);

		// do the first frequency
		sinfreq = sin(constants::pi * pulse.freqs[0] * t / pulse.T);
		for(n = 0; n < 3; ++n) {
			scat[n] = copy(fh[0]->scat[n]);
			scat[n].scale(sinfreq);
		}

		// iterate over the other frequencies
		for(n = 1; n < pulse.nfreq; ++n) {
			sinfreq = sin(constants::pi * pulse.freqs[n] * t / pulse.T);
			for(int i = 0; i < 3; ++i) {
				scat[i].compress();
				fh[n]->scat[i].compress();
				scat[i].gaxpy(1.0, fh[n]->scat[i], sinfreq);
			}
		}

		// compute the incident pulse at this time
		field[0] = FunctionFactory<double, 3>(group).functor(time_inc);
		field[1] = copy(field[0]);
		field[2] = copy(field[0]);
		for(n = 0; n < 3; ++n) {
			field[n].compress();
			scat_real = real(scat[n]);
			scat_real.compress();
			field[n].gaxpy(pulse.e_polar[n], scat_real, 1.0);
			field[n].truncate();
			scat_real.clear();
		}

		// visualization
		char filename[100];
		sprintf(filename, "maxwell.%d.vts", ti);
		Vector<double, 3> plotlo, plothi;
		Vector<long, 3> npts;
		for(n = 0; n < 3; ++n) {
			plotlo[n] = pulse.sim_min;
			plothi[n] = pulse.sim_max;
			npts[n] = 71;
		}
		plotvtk_begin(group, filename, plotlo, plothi, npts);
		plotvtk_data(field[0], "x", group, filename, plotlo, plothi, npts);
		plotvtk_data(field[1], "y", group, filename, plotlo, plothi, npts);
		plotvtk_data(field[2], "z", group, filename, plotlo, plothi, npts);
		plotvtk_end<3>(group, filename);

		for(n = 0; n < 3; ++n)
			scat[n].clear();
	}

	// clean up
	for(n = 0; n < pulse.nfreq; ++n)
		if(fh[n])
			delete fh[n];

	delete [] fh;

	madness::finalize();

	return 0;
}

inline void calcEpshat(const double omega, const double ag_eps,
		const double ag_sigma, const Function<double, 3> &ag_shape,
		Function<complexd, 3> &epshat, Function<complexd, 3> &epshatI,
		Function<complexd, 3> &lnepshat) {

	// all eps and sigma relative to eps0
	// epshatI = I omega (ag_eps - 1) - ag_sigma in silver,
	//           0 in vacuum
	epshatI = complexd(-ag_sigma, omega*(ag_eps-1.0))*copy(ag_shape);

	// epshat = I omega ag_eps - ag_sigma in silver
	//          I omega in vacuum
	epshat = copy(epshatI);
	epshat.add_scalar(complexd(0.0, omega));

	// lnepshat = log(epshat)
	complexd vac_lnepshat(log(omega), constants::pi / 2.0);
	complexd rel_ag_lnepshat(-ag_sigma, omega*ag_eps);
	rel_ag_lnepshat = log(rel_ag_lnepshat) - vac_lnepshat;
	lnepshat = rel_ag_lnepshat*copy(ag_shape);
	lnepshat.add_scalar(vac_lnepshat);
}

void calcScatField(const double omega, const double ag_eps,
		const double ag_sigma, const Function<double, 3> &ag_shape,
		FrequencyHandler &fh, const Function<double, 3> &box_mask,
		const Function<double, 3> &grad_mask,
		const SeparatedConvolution<double, 3> &G) {

	int i;
	Function<complexd, 3> epshat, epshatI, lnepshat;

	calcEpshat(omega, ag_eps, ag_sigma, ag_shape, epshat, epshatI, lnepshat);

	// the gradient of lnepshat: multiply first by box_mask to make the function
	// 0 at the boundary, then damp out this derivative with grad_mask
	vecfunc gradlnepshat, rhs;
	for(i = 0; i < 3; ++i) {
		gradlnepshat[i] = diff(lnepshat*box_mask, i) * grad_mask;
		gradlnepshat[i].truncate();
	}

	lnepshat.clear();

	// Compute the rhs function
	Function<complexd, 3> dotp, graddotp;
	dotp = (gradlnepshat[0]*fh.inc[0] + gradlnepshat[1]*fh.inc[1] +
		gradlnepshat[2]*fh.inc[2]) * box_mask;
	for(i = 0; i < 3; ++i) {
		graddotp = diff(dotp, i) * grad_mask;
		graddotp.compress();
		rhs[i] = epshatI*fh.inc[i];
		rhs[i].compress();
		rhs[i].gaxpy(complexd(0.0, omega*fh.pulse.mu0*fh.pulse.eps0), graddotp,
			complexd(-1.0, 0.0));
		rhs[i] = apply(G, rhs[i]);
		rhs[i].scale(-1.0);
		rhs[i].truncate();
		graddotp.clear();
	}

	dotp.clear();

	// set up an initial guess for the scattered field (0)
	fh.scat[0] = complexd(0.0, 0.0) * copy(ag_shape);
	fh.scat[0].truncate();
	fh.scat[1] = copy(fh.scat[0]);
	fh.scat[2] = copy(fh.scat[0]);

	// Set up the linear operator and solve
	int maxiters = 7;
	double thresh = 1.0e-2;
	VectorOfFunctionsSpace<complexd, 3, 3> space;
	EFieldOperator efop(epshat, omega, fh.pulse.mu0, fh.pulse.eps0,
		gradlnepshat, box_mask, grad_mask, G);
	GMRES(space, efop, rhs, fh.scat, maxiters, thresh, true);
	fh.scat[0].truncate();
	fh.scat[1].truncate();
	fh.scat[2].truncate();

#if 0
		// visualization
		char filename[100];
		sprintf(filename, "maxwell.vts");
		Vector<double, 3> plotlo, plothi;
		Vector<long, 3> npts;
		for(int n = 0; n < 3; ++n) {
			plotlo[n] = 50.0 /* nm */ / fh.pulse.nm_per_au;
			plothi[n] = 150.0 /* nm */ / fh.pulse.nm_per_au;
			npts[n] = 71;
		}
		plotvtk_begin(fh.world, filename, plotlo, plothi, npts);
		plotvtk_data(fh.scat[0], "x", fh.world, filename, plotlo, plothi, npts);
		plotvtk_data(fh.scat[1], "y", fh.world, filename, plotlo, plothi, npts);
		plotvtk_data(fh.scat[2], "z", fh.world, filename, plotlo, plothi, npts);
		plotvtk_end<3>(fh.world, filename);
#endif

	epshat.clear();
	epshatI.clear();
	gradlnepshat[0].clear();
	gradlnepshat[1].clear();
	gradlnepshat[2].clear();
	rhs[0].clear();
	rhs[1].clear();
	rhs[2].clear();
}
