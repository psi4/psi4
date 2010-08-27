/*
 * coupling.h
 *
 *  Created on: Feb 24, 2009
 *      Author: jturney
 */

#ifndef COUPLING_H_
#define COUPLING_H_

#define CLOSED  0
#define ALPHA   1
#define BETA    2
#define VIRTUAL 3

/* Coupling cofficients */
static const double f_i[] = {
		1.0,
		0.5,
		0.5,
		0.0
};

static const double alpha_ij[][4] = {
		{ 2.0, 1.0, 1.0, 0.0 },
		{ 1.0, 0.5, 0.5, 0.0 },
		{ 1.0, 0.5, 0.5, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0 },
};

static const double beta_ij[][4] = {
		{ -1.0, -0.5, -0.5, 0.0 },
		{ -0.5, -0.5,  0.5, 0.0 },
		{ -0.5,  0.5, -0.5, 0.0 },
		{  0.0,  0.0,  0.0, 0.0 }
};

#endif /* COUPLING_H_ */
