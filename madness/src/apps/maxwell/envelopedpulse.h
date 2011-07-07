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
  
  $Id: envelopedpulse.h 1602 2009-12-27 19:53:06Z rjharrison $
*/
#ifndef __enveloped_pulse__
#define __enveloped_pulse__

#include <string>
#include <constants.h>
#include <mra/mra.h>
#include <misc/interpolation_1d.h>

using namespace std;
using namespace madness;

typedef double_complex complexd;
typedef Function<complexd, 3> function;
typedef Vector<function, 3> vecfunc;

/// A class for storing all the pertinent data for an enveloped incident
/// pulse.
///
/// NOTE that all units are in atomic units, unless specified otherwise.
class EnvelopedPulse {
public:
	// physical constants
	const static double eps0 = 0.0795774715459477; // free-space permittivity
	const static double mu0 = 6.69176253e-4; // free-space permeability
	const static double c_vac = 1.370359991e2; // vacuum light speed

	// unit conversions
	const static double nm_per_au = 0.05291772108;
	const static double fs_per_au = 0.02418884326505;
	const static double efield_per_mau = 5.142206435e5;

	// key parameters
	double lambda; // driving wavelength
	double omegabar; // corresponding driving frequency
	double T; // Fourier period
	double stddev; // standard deviation of the wave-shape
	double offset; // offset of the wave-shape

	// box sizes / ranges
	double sim_min, sim_max; // the 3-D box sizes
	double dotp_min, dotp_max; // the dotp-product range for 1-D incident fields
	double ws_min, ws_max; // the 1-D wave-shape function

	// polarizations: these three vectors should be normalized and
	// mutually orthogonal
	Vector<double, 3> prop_dir; //< Direction of field propagation
	Vector<double, 3> e_polar; //< Initial electric field polarization
	Vector<double, 3> h_polar; //< Initial magnetization polarization

	string basefile; // the base file *.param, *.freq.n, etc.
	int nfreq; // the number of frequencies needed
	int *freqs; // the list of frequencies to be computed
	CubicInterpolationTable<double> *wave_shape;

	EnvelopedPulse() : nfreq(0), freqs(0), wave_shape(0) {
		prop_dir[0] = prop_dir[2] = e_polar[0] = e_polar[1] = h_polar[1] =
			h_polar[2] = 0.0;
		prop_dir[1] = e_polar[2] = h_polar[0] = 1.0;
	}

	~EnvelopedPulse() {
		if(freqs)
			delete [] freqs;
		if(wave_shape)
			delete wave_shape;
	}

	/// Read a *.param file to get the parameters, etc.
	/// Returns true on successful read, false otherwise
	bool read_param_file(const string &basename);
};

/// This functor recreates the incident field from an EnvelopedPulse object
/// at a given time.
class TimeIncident : public FunctionFunctorInterface<double, 3> {
protected:
	const EnvelopedPulse &pulse; //< The incident pulse information
	double t; //< The time at which to evaluate, in atomic units

public:
	TimeIncident(const EnvelopedPulse &_pulse) : pulse(_pulse), t(0.0) {}

	void setTime(double _t) {
		t = _t;
	}

	double getTime() const {
		return t;
	}

	/// This returns the magnitude of the incident electric field (still
	/// need to multiply by the unit direction vector)
	double operator() (const Vector<double, 3> &pt) const;
};

#endif
