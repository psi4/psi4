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
/// Implementation of functions in envelopedpulse.h

#include "envelopedpulse.h"
#include <vector>

bool EnvelopedPulse::read_param_file(const string &basename) {
	char filename[80];
	FILE *file;

	// open the file
	sprintf(filename, "%s.param", basename.c_str());
	file = fopen(filename, "r");
	if(file == NULL) {
		fprintf(stderr, "Unable to open file %s\n", filename);
		return false;
	}

	if(!fgets(filename, 80, file)) {
		fprintf(stderr, "Error reading file\n");
		fclose(file);
		return false;
	}

	// start reading the parameters:
	// driving wavelength, standard deviation, offset, period, 3-D box size,
	// 1-D dot product range, 1-D wave-shape range, number of frequencies
	if(fscanf(file, "%le%le%le%le%le%le%le%le%le%le%d", &lambda, &stddev,
			&offset, &T, &sim_min, &sim_max, &dotp_min, &dotp_max, &ws_min,
			&ws_max, &nfreq) < 11) {
		fprintf(stderr, "Error reading file\n");
		fclose(file);
		return false;
	}

	// calculate the equivalent driving frequency
	omegabar = 2.0 * constants::pi * c_vac / lambda;

	// read in the frequencies
	freqs = new int[nfreq];
	for(int i = 0; i < nfreq; ++i) {
		if(fscanf(file, "%d", freqs + i) < 1) {
			fprintf(stderr, "Error reading file\n");
			fclose(file);
			return false;
		}
	}

	// read in the wave shape
	int ninc;
	if(fscanf(file, "%d", &ninc) < 1) {
		fprintf(stderr, "Error reading file\n");
		fclose(file);
		return false;
	}
	vector<double> inc(ninc);
	for(int i = 0; i < ninc; ++i) {
		if(fscanf(file, "%le", &inc[i]) < 1) {
			fprintf(stderr, "Error reading file\n");
			fclose(file);
			return false;
		}
	}

	fclose(file);

	wave_shape = new CubicInterpolationTable<double>(ws_min, ws_max, ninc, inc);

	printf("\n%s\n", filename);
	printf("Driving wavelength: %.6e nm\n", lambda * nm_per_au);
	printf("Gaussian envelope std. dev.: %.6e fs / %.6e nm\n", stddev * lambda /
		(2.0 * constants::pi * c_vac) * fs_per_au, stddev * lambda /
		(2.0 * constants::pi) * nm_per_au);
	printf("Fourier period: %.6e fs\n", T * fs_per_au);
	printf("3-D box dimensions: %.6e -- %.6e nm\n", sim_min * nm_per_au,
		sim_max * nm_per_au);
	printf("Number of frequencies needed: %d\n", nfreq);
	
	basefile = basename;
	return true;
}

double TimeIncident::operator() (const Vector<double, 3> &pt) const {
	double dotp = 0.0;
	for(int i = 0; i < 3; ++i)
		dotp += pt[i]*pulse.prop_dir[i];
	dotp *= 2.0 * constants::pi / pulse.lambda;
	dotp -= pulse.omegabar*t;

	return pulse.wave_shape->operator() (dotp);
}
