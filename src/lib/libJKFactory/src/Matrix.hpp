/*
 * JKFactory: Interface and code for highly parallel J and K
 *             builds.
 *
 *  Copyright (c) 2014 Ryan M. Richard
 *
 *  This file is part of JKFactory.
 *
 *  JKFactory is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef MATRIX_HPP_
#define MATRIX_HPP_
#include "JKFactory.hpp"
#include <boost/shared_ptr.hpp>
namespace JKFactory{

///Simple wrapper for a square matrix
class MatrixMPIData:public JKFactoryBase{
	private:
		///Routine that calculates n process per row and column
		void SplitProc();
		void DivvyMatrix(int NBasis);
	public:
		int NPRow;
		int NPCol;
		int nBlkRow;
		int nBlkCol;
		int startrow;
		int endrow;
		int startcol;
		int endcol;
		int nelements(){return (endcol-startcol+1)*(endrow-startrow+1);}
		boost::shared_ptr<double[]> MyBuffer;
		boost::shared_ptr<int[]> RowPerBlock;
		boost::shared_ptr<int[]> ColPerBlock;

		///Given i,j of the full matrix returns the element, no checks
		double& operator()(const int i, const int j){
		   return MyBuffer[(i-startrow)*(endcol-startcol+1)+(j-startcol)];
		}
		MatrixMPIData(int NBasis);
		void PrintOut()const{}
};


class Matrix:public JKFactoryBase{
	private:
		/** This pointer is the data that was given to us.
		 *  We are not responsible for it's memory management.
		 *  Do not allocate it.  Do not Delete it. This should
		 *  be the whole matrix, i.e. the gathered one, not the
		 *  distributed one.
		 */
		double* GatheredData;

		///GatherData is nbasis by nbasis
		int nbasis;

		///The data related to distributing this matrix, includes buffer
		boost::shared_ptr<MatrixMPIData> MyData;
	public:

		///Pointer to the distribution data
		boost::shared_ptr<MatrixMPIData> MPIData(){return MyData;}

		///Fills the buffer of MyData with the block and returns it
		double* GetMyBlock();

		///Switches out the pointer for GatheredData
		void Update(double* NewMat){GatheredData=NewMat;}

		void Gather();

		Matrix(double* Matrix2Store,int Dim);

		///Don't free the GatheredData
		~Matrix(){}

		///Returns the i,j-th element
		const double& operator()(const int i,const int j)const{
			return GatheredData[i*nbasis+j];
		}

		///Returns the i,j-th element
		double& operator()(const int i,const int j){
			return GatheredData[i*nbasis+j];
		}

		///Debug printing
		void PrintOut()const;
};

}




#endif /* MATRIX_HPP_ */
