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
  
  $Id: frequencyhandler.h 1602 2009-12-27 19:53:06Z rjharrison $
*/
#ifndef __frequency_handler__
#define __frequency_handler__

#include <string>
#include <constants.h>
#include <mra/mra.h>
#include "envelopedpulse.h"

using namespace std;
using namespace madness;

typedef double_complex complexd;
typedef Function<complexd, 3> function;
typedef Vector<function, 3> vecfunc;

/// A class for storing all the pertinent data for each frequency-specific
/// field component.
///
/// NOTE that all units are in atomic units, unless specified otherwise.
///
/// This class also extends FunctionFunctorInterface to convert the 1-D
/// incident pulse information into the 3-D incident pulse components.
/// The polarization and propagation directions are stored in the
/// EnvelopedPulse object.
class FrequencyHandler {
public:
	const int freq_index;
	Vector<Function<double, 3>, 3> inc; //< the 3-D incident pulse
	vecfunc scat; //< the 3-D scattered field
	const EnvelopedPulse &pulse; //< The incident pulse / simulation parameters
	World &world; //< The World to use when making functions

	/// Constructor: needs the frequency index, the simulation parameters,
	/// and the MPI communicator (World)
	FrequencyHandler(const int freq, const EnvelopedPulse &_pulse,
		World &_world) : freq_index(freq), pulse(_pulse), world(_world) {}

	/// Read the *.freq.n file to get the parameters, etc.
	/// The specific file name is obtained from pulse.
	/// Returns true on successful read, false otherwise
	bool read_file();
};

/// A class for turning the frequency-dependent incident field function into
/// a 3-D madness function (the dot product is used to convert the 3-D point
/// into the 1-D argument)
class FrequencyIncident : public FunctionFunctorInterface<double, 3> {
protected:
	const EnvelopedPulse &pulse; //< The incident pulse / simulation parameters
	const CubicInterpolationTable<double> &interp; //< The interpolation data

public:
	/// Constructor
	FrequencyIncident(const EnvelopedPulse &_pulse,
		const CubicInterpolationTable<double> &_interp) : pulse(_pulse),
		interp(_interp) {}

	/// The interface for turning the 1-D incident function into the 3-D
	/// vector.  This really works best if the E-field polarization is
	/// only in one of the vector directions.
	double operator() (const Vector<double, 3> &pt) const;
};

#endif
