/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */
#ifndef ALGORITHMS_H_
#define ALGORITHMS_H_

namespace psi{
namespace LibParallel{

/** \brief An enumeration of all the algorithms I have coded up
 *
 *
 *  Static algorithms:
 *  -SIMPLE: MPI process i gets task i*batchsize through (i+1)*batchsize
 *  -ROUNDROBIN: Processes are given statically in a round-robin fashion
 *  -RANDOM: Processes are assigned random priorities then assigned via SIMPLE
 *
 *
 *  Dynamic algorithms:
 *  -DYNAMICRR: Tasks are assigned to masters in a round-robin fashion, then
 *      they are given to processes as requested.
 *  -DYNAMICRANDOM : Tasks are assigned random priorities then given out as
 *      requested
 *
 */

enum SchedAlgorithm{NONE,SIMPLE,STATICRR,DYNAMICRR};
enum MPIOperation{MULTIPLY,ADD,SUBTRACT,DIVIDE,MODULUS};
}}

#endif /* ALGORITHMS_H_ */