/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2016-2017 Robert A. Shaw.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/* 
 * Helpful lightweight multi-index arrays to make the code easier to write and test. 
 * These should probably be replaced at some point in the interests of speed. 
 *
 * Robert A Shaw 2016
 */

#ifndef MULTIARR_HEAD
#define MULTIARR_HEAD

#include <vector>

namespace psi {

/// Templated skeleton two index array for convenience
template<typename T>
struct TwoIndex {
	int dims[2];
	std::vector<T> data;
	T& operator()(int i, int j) { return data[i * dims[1] + j]; }
	T operator()(int i, int j) const { return data[i * dims[1] + j]; }
	void assign(int dim1, int dim2, T value) {
		dims[0] = dim1; dims[1] = dim2;
		data.resize(dim1 * dim2);
		std::fill(data.begin(), data.end(), value);
	}
	TwoIndex() { dims[0] = dims[1] = 0; }
	TwoIndex(int dim1, int dim2) {
		dims[0] = dim1; dims[1] = dim2;
		data.resize(dim1 * dim2);
	}
	TwoIndex(int dim1, int dim2, T value) { assign(dim1, dim2, value); }
	TwoIndex(const TwoIndex<T> &other) { 
		data = other.data;
		dims[0] = other.dims[0]; dims[1] = other.dims[1];
	}
};

/// Templated skeleton three index array for convenience
template<typename T>
struct ThreeIndex {
	int dims[3];
	std::vector<T> data;
	T& operator()(int i, int j, int k) { return data[i*dims[2]*dims[1] + j*dims[2] + k]; }
	T operator()(int i, int j, int k) const { return data[i*dims[2]*dims[1] + j*dims[2] + k]; }
	ThreeIndex(){ dims[0] = 0; dims[1] = 0; dims[2] = 0; }
	ThreeIndex(int dim1, int dim2, int dim3) {
		dims[0] = dim1; dims[1] = dim2; dims[2] = dim3;
		data.resize(dim1 * dim2 * dim3);
	}
	ThreeIndex(const ThreeIndex<T> &other)  { 
		data = other.data;
		for (int n = 0; n < 3; n++) dims[n] = other.dims[n]; 
	}
	void fill(T value) { std::fill(data.begin(), data.end(), value); }
};

/// Templated skeleton five index array for convenience
template<typename T>
struct FiveIndex {
	int dims[5];
	std::vector<T> data;
	T& operator()(int i, int j, int k, int l, int m) { 
		return data[m + dims[4] * (l + dims[3] * (k + dims[2] * (j + dims[1] * i)))]; 
	}
	T operator()(int i, int j, int k, int l, int m) const { 
		return data[m + dims[4] * (l + dims[3] * (k + dims[2] * (j + dims[1] * i)))];
	}
	FiveIndex() { dims[0] = dims[1] = dims[2] = dims[3] = dims[4] = 0; }
	FiveIndex(int dim1, int dim2, int dim3, int dim4, int dim5) {
		dims[0] = dim1; dims[1] = dim2; dims[2] = dim3; dims[3] = dim4; dims[4] = dim5;
		data.resize(dim1 * dim2 * dim3 * dim4 * dim5);
	}
	FiveIndex(const FiveIndex<T> &other) { 
		data = other.data;
		for (int n = 0; n < 5; n++) dims[n] = other.dims[n]; 
	}
};

/// Templated skeleton seven index array for convenience
template<typename T>
struct SevenIndex {
	int dims[7];
	std::vector<T> data;
	T& operator()(int i, int j, int k, int l, int m, int n, int p) {
		return data[p + dims[6]*(n + dims[5] * (m + dims[4] * (l + dims[3] * (k + dims[2] * (j + dims[1] * i)))))];
	}
	T operator()(int i, int j, int k, int l, int m, int n, int p) const {
		return data[p + dims[6]*(n + dims[5] * (m + dims[4] * (l + dims[3] * (k + dims[2] * (j + dims[1] * i)))))];
	}
	SevenIndex() { dims[0] = dims[1] = dims[2] = dims[3] = dims[4] = dims[5] = dims[6] = 0; }
	SevenIndex(int dim1, int dim2, int dim3, int dim4, int dim5, int dim6, int dim7) {
		dims[0] = dim1; dims[1] = dim2; dims[2] = dim3; dims[3] = dim4; dims[4] = dim5; dims[5] = dim6; dims[6] = dim7;
		data.resize(dim1 * dim2 * dim3 * dim4 * dim5 * dim6 * dim7);
	}
	SevenIndex(const SevenIndex<T> &other) {
		data = other.data;
		for (int n = 0; n < 7; n++) dims[n] = other.dims[n];
	}
};
}
#endif
