/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
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

#include "psi4/libiwl/iwl.hpp"
#include "psi4/libqt/qt.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libmints/wavefunction.h"
#include "defines.h"
#include "arrays.h"
#include "dpd.h"
#include "psi4/libparallel/ParallelPrinter.h"

using namespace std;

namespace psi{ namespace occwave{

/********************************************************************************************/
/************************** SymBlockMatrix **************************************************/
/********************************************************************************************/
SymBlockMatrix::SymBlockMatrix()
{
    matrix_ = NULL;
    rowspi_ = NULL;
    colspi_ = NULL;
}

SymBlockMatrix::SymBlockMatrix(std::string name)
{
    matrix_ = NULL;
    rowspi_ = NULL;
    colspi_ = NULL;
}


SymBlockMatrix::SymBlockMatrix(int nirreps, int *ins_rowspi, int *ins_colspi)
{
    matrix_ = NULL;
    nirreps_ = nirreps;
    rowspi_ = new int[nirreps_];
    colspi_ = new int[nirreps_];
    for (int h=0; h<nirreps_; h++) {
        rowspi_[h] = ins_rowspi[h];
        colspi_[h] = ins_colspi[h];
    }
    memalloc();
}//

SymBlockMatrix::SymBlockMatrix(std::string name, int nirreps, int *ins_rowspi, int *ins_colspi)
{
    matrix_ = NULL;
    name_ = name;
    nirreps_ = nirreps;
    rowspi_ = new int[nirreps_];
    colspi_ = new int[nirreps_];
    for (int h=0; h<nirreps_; h++) {
        rowspi_[h] = ins_rowspi[h];
        colspi_[h] = ins_colspi[h];
    }
    memalloc();
}//


SymBlockMatrix::~SymBlockMatrix()
{
    release();
    if (rowspi_) delete[] rowspi_;
    if (colspi_) delete[] colspi_;
}//

int *SymBlockMatrix::rowspi()
{
   return rowspi_;
}
int *SymBlockMatrix::colspi()
{
    return colspi_;
}

SymBlockMatrix* SymBlockMatrix::generate(int nirreps, int *ins_rowspi, int *ins_colspi)
{
    return new SymBlockMatrix(nirreps, ins_rowspi, ins_colspi);
}

SymBlockMatrix* SymBlockMatrix::generate(std::string name, int nirreps, int *ins_rowspi, int *ins_colspi)
{
    return new SymBlockMatrix(name, nirreps, ins_rowspi, ins_colspi);
}

void SymBlockMatrix::init(std::string name, int nirreps, int *ins_rowspi, int *ins_colspi)
{
    if (rowspi_) delete[] rowspi_;
    if (colspi_) delete[] colspi_;
    name_ = name;
    nirreps_ = nirreps;
    rowspi_ = new int[nirreps_];
    colspi_ = new int[nirreps_];
    for (int h=0; h<nirreps_; h++) {
        rowspi_[h] = ins_rowspi[h];
        colspi_[h] = ins_colspi[h];
    }
    memalloc();
}//

void SymBlockMatrix::memalloc()
{
    if (matrix_) release();
    matrix_ = (double***)malloc(sizeof(double***) * nirreps_);
    for (int h=0; h<nirreps_; h++) {
      if (rowspi_[h] != 0 && colspi_[h] != 0) {
	matrix_[h] = block_matrix(rowspi_[h], colspi_[h]);
      }
      else matrix_[h] = NULL;
    }
}//

void SymBlockMatrix::release()
{
    if (!matrix_) return;
    for (int h=0; h<nirreps_; h++) {
      if (matrix_[h]) free_block(matrix_[h]);
    }
    matrix_ = NULL;
}//

void SymBlockMatrix::zero()
{
    size_t size;
    for (int h=0; h<nirreps_; h++) {
      size = rowspi_[h] * colspi_[h] * sizeof(double);
      if (size) {
	memset(&(matrix_[h][0][0]), 0, size);
      }
    }
}//

void SymBlockMatrix::zero_diagonal()
{
    for (int h=0; h<nirreps_; h++) {
      for (int i=0; i<MIN0(rowspi_[h], colspi_[h]); ++i) {
	matrix_[h][i][i] = 0.0;
       }
    }
}//

double SymBlockMatrix::trace()
{
    double value;
    value=0.0;
    for (int h=0; h<nirreps_; h++) {
      for (int i=0; i<MIN0(rowspi_[h], colspi_[h]); ++i) {
	value += matrix_[h][i][i];
      }
    }
    return value;
}//

SymBlockMatrix* SymBlockMatrix::transpose()
{
    SymBlockMatrix* temp;
    temp = new SymBlockMatrix(nirreps_, colspi_,rowspi_);
    temp->zero();
    for (int h=0; h<nirreps_; h++) {
      for (int i=0; i<colspi_[h]; ++i) {
	for (int j=0; j<rowspi_[h]; ++j) {
	  temp->set(h,i,j,matrix_[h][j][i]);
	}
      }
    }
    return temp;
}//


void SymBlockMatrix::copy(const SymBlockMatrix* Adum)
{
    // Make sure that matrices are in the same size
    bool same = true;
    for (int h=0; h<nirreps_; h++) {
      if (colspi_[h] != Adum->colspi_[h] || rowspi_[h] != Adum->rowspi_[h])  same = false;
    }

    if (same == false) {
        release();
        if (rowspi_) delete[] rowspi_;
        if (colspi_) delete[] colspi_;
        rowspi_ = new int[nirreps_];
        colspi_ = new int[nirreps_];
        for (int i=0; i<nirreps_; ++i) {
            rowspi_[i] = Adum->rowspi_[i];
            colspi_[i] = Adum->colspi_[i];
        }
        memalloc();
    }

    // If matrices are in the same size
    for (int h=0; h<nirreps_; h++) {
      if (rowspi_[h] != 0 && colspi_[h] != 0) {
	memcpy(&(matrix_[h][0][0]), &(Adum->matrix_[h][0][0]), rowspi_[h] * colspi_[h] * sizeof(double));
      }
    }
}//

void SymBlockMatrix::add(const SymBlockMatrix* Adum)
{
    double *lhs, *rhs;
    for (int h=0; h<nirreps_; h++) {
      size_t size = rowspi_[h] * colspi_[h];
      if (size) {
	lhs = matrix_[h][0];
        rhs = Adum->matrix_[h][0];
        for (size_t cnt=0; cnt<size; cnt++) {
	  *lhs += *rhs;
          lhs++; rhs++;
	}
      }
    }
}//

void SymBlockMatrix::add(int h, int i, int j, double value)
{
  matrix_[h][i][j]+=value;
}//

void SymBlockMatrix::subtract(const SymBlockMatrix* Adum)
{
    double *lhs, *rhs;
    for (int h=0; h<nirreps_; h++) {
      size_t size = rowspi_[h] * colspi_[h];
      if (size) {
	lhs = matrix_[h][0];
        rhs = Adum->matrix_[h][0];
        for (size_t cnt=0; cnt<size; cnt++) {
	  *lhs -= *rhs;
          lhs++; rhs++;
	}
      }
    }
}//

void SymBlockMatrix::subtract(int h, int i, int j, double value)
{
  matrix_[h][i][j]-=value;
}//

void SymBlockMatrix::scale(double a)
{
    size_t size;
    for (int h=0; h<nirreps_; h++) {
      size = rowspi_[h] * colspi_[h];
      if (size) C_DSCAL(size, a, &(matrix_[h][0][0]), 1);
    }
}//

void SymBlockMatrix::scale_row(int h, int m, double a)
{
    C_DSCAL(rowspi_[h], a, &(matrix_[h][m][0]), 1);
}//

void SymBlockMatrix::scale_column(int h, int n, double a)
{
    C_DSCAL(colspi_[h], a, &(matrix_[h][0][n]), rowspi_[h]);
}//

double SymBlockMatrix::sum_of_squares()
{
    double summ;
    summ=0.0;
    for (int h=0; h<nirreps_; h++) {
      for (int i=0; i<rowspi_[h]; ++i) {
	for (int j=0; j<colspi_[h]; ++j) {
	  summ += matrix_[h][i][j] * matrix_[h][i][j];
	}
      }
    }
    return summ;
}//

double SymBlockMatrix::rms()
{
    double summ;
    int dim;
    summ=0.0;
    dim=0;

    for (int h=0; h<nirreps_; h++) {
      if (rowspi_[h] != 0 && colspi_[h] != 0) dim += rowspi_[h] * colspi_[h];
    }

    for (int h=0; h<nirreps_; h++) {
      for (int i=0; i<rowspi_[h]; ++i) {
	for (int j=0; j<colspi_[h]; ++j) {
	  summ += matrix_[h][i][j] * matrix_[h][i][j];
	}
      }
    }
    summ=sqrt(summ)/dim;
    return summ;
}//

double SymBlockMatrix::rms(SymBlockMatrix* Atemp)
{
    double summ;
    int dim;
    summ=0.0;
    dim=0;

    for (int h=0; h<nirreps_; h++) {
      if (rowspi_[h] != 0 && colspi_[h] != 0) dim += rowspi_[h] * colspi_[h];
    }

    for (int h=0; h<nirreps_; h++) {
      for (int i=0; i<rowspi_[h]; ++i) {
	for (int j=0; j<colspi_[h]; ++j) {
	  summ += (matrix_[h][i][j] - Atemp->matrix_[h][i][j]) * (matrix_[h][i][j] - Atemp->matrix_[h][i][j]);
	}
      }
    }
    summ=sqrt(summ)/dim;
    return summ;
}//

void SymBlockMatrix::set(double value)
{
    size_t size;
    for (int h=0; h<nirreps_; h++) {
      size = rowspi_[h] * colspi_[h];
      for (size_t i=0; i<size; ++i) {
	matrix_[h][0][i] = value;
      }
    }
}//

void SymBlockMatrix::set(int h, int i, int j, double value)
{
  matrix_[h][i][j]=value;
}//

void SymBlockMatrix::set(double **Asq)
{
    int offset;
    if (Asq == NULL) return;

    offset = 0;
    for (int h=0; h<nirreps_; h++) {
      for (int i=0; i<rowspi_[h]; ++i) {
	int ii = i + offset;
        for (int j=0; j<=i; ++j) {
	  int jj = j + offset;
          matrix_[h][i][j] = Asq[ii][jj];
          matrix_[h][j][i] = Asq[jj][ii];
        }
      }
      offset += rowspi_[h];
    }
}//

void SymBlockMatrix::set(dpdbuf4 G)
{
    for(int h = 0; h < nirreps_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);
        for(int row = 0; row < G.params->rowtot[h]; ++row){
            for(int col = 0; col < G.params->coltot[h]; ++col){
                matrix_[h][row][col] = G.matrix[h][row][col];
            }
        }
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
}//

double SymBlockMatrix::get(int h, int m, int n)
{
  return matrix_[h][m][n];
}//

double *SymBlockMatrix::to_lower_triangle()
{
    int sizerow=0, sizecol=0;
    for (int h=0; h<nirreps_; h++) {
        sizerow += rowspi_[h];
        sizecol += colspi_[h];
    }
    if (sizerow != sizecol)
        return NULL;

    double *tri = new double[ioff[sizerow]];
    double **temp = to_block_matrix();
    sq_to_tri(temp, tri, sizerow);
    free_block(temp);
    return tri;
}//

double **SymBlockMatrix::to_block_matrix()
{
    int sizerow=0, sizecol=0;
    for (int h=0; h<nirreps_; h++) {
        sizerow += rowspi_[h];
        sizecol += colspi_[h];
    }

    double **temp = block_matrix(sizerow, sizecol);
    int offsetrow = 0, offsetcol=0;
    for (int h=0; h<nirreps_; h++) {
        for (int i=0; i<rowspi_[h]; ++i) {
            for (int j=0; j<colspi_[h]; ++j) {
                temp[i+offsetrow][j+offsetcol] = matrix_[h][i][j];
            }
        }
        offsetrow += rowspi_[h];
        offsetcol += colspi_[h];
    }

    return temp;
}//

void SymBlockMatrix::print(std::string OutFileRMR)
{
   std::shared_ptr<psi::PsiOutStream> printer=(OutFileRMR=="outfile"?outfile:
         std::shared_ptr<OutFile>(new OutFile(OutFileRMR,APPEND)));
   if (name_.length()) printer->Printf( "\n ## %s ##\n", name_.c_str());
    for (int h=0; h<nirreps_; h++) {
      if (rowspi_[h] != 0 && colspi_[h] != 0) {
	printer->Printf( "\n Irrep: %d\n", h+1);
	print_mat(matrix_[h], rowspi_[h], colspi_[h], OutFileRMR);
	printer->Printf( "\n");
      }
    }
}//

void SymBlockMatrix::print()
{
    if (name_.length()) outfile->Printf( "\n ## %s ##\n", name_.c_str());
    for (int h=0; h<nirreps_; h++) {
      if (rowspi_[h] != 0 && colspi_[h] != 0) {
	outfile->Printf( "\n Irrep: %d\n", h+1);
	print_mat(matrix_[h], rowspi_[h], colspi_[h], "outfile");
	outfile->Printf( "\n");
      }
    }

}//

void SymBlockMatrix::set_to_identity()
{
    size_t size;
    for (int h=0; h<nirreps_; h++) {
      size = rowspi_[h] * colspi_[h] * sizeof(double);
      if (size) {
	memset(&(matrix_[h][0][0]), 0, size);
        for (int i=0; i<MIN0(rowspi_[h], colspi_[h]); i++)
	  matrix_[h][i][i] = 1.0;
      }
    }
}//

void SymBlockMatrix::gemm(bool transa, bool transb, double alpha, const SymBlockMatrix* a, const SymBlockMatrix* b, double beta)
{
    char ta = transa ? 't' : 'n';
    char tb = transb ? 't' : 'n';
    int  m, n, k, nca, ncb, ncc;

    for (int h=0; h<nirreps_; h++) {
        m = rowspi_[h];
        n = colspi_[h];
        k = a->colspi_[h];
        nca = transa ? m : k;
        ncb = transb ? k : n;
        ncc = n;

        if (m && n && k) {
            C_DGEMM(ta, tb, m, n, k, alpha, &(a->matrix_[h][0][0]),
                    nca, &(b->matrix_[h][0][0]), ncb, beta, &(matrix_[h][0][0]),
                    ncc);
        }
    }
}//

bool SymBlockMatrix::load(PSIO* psio, int itap, const char *label, int dim)
{
    int ntri = 0.5 * dim * (dim  + 1);
    double *mybuffer=init_array(ntri);
    memset(mybuffer, 0.0, sizeof(double)*ntri);
    IWL::read_one(psio, itap, label, mybuffer, ntri, 0, 0, "outfile");

    double **Asq;
    Asq=block_matrix(dim,dim);
    memset(Asq[0], 0.0, sizeof(double)*dim*dim);
    tri_to_sq(mybuffer,Asq,dim);
    free(mybuffer);

    set(Asq);
    free_block(Asq);
    return true;
}//

bool SymBlockMatrix::load(std::shared_ptr<psi::PSIO> psio, int itap, const char *label, int dim)
{
    int ntri = 0.5 * dim * (dim  + 1);
    double *mybuffer=init_array(ntri);
    memset(mybuffer, 0.0, sizeof(double)*ntri);
    IWL::read_one(psio.get(), itap, label, mybuffer, ntri, 0, 0, "outfile");

    double **Asq;
    Asq=block_matrix(dim,dim);
    memset(Asq[0], 0.0, sizeof(double)*dim*dim);
    tri_to_sq(mybuffer,Asq,dim);
    free(mybuffer);

    set(Asq);
    free_block(Asq);
    return true;
}//

void SymBlockMatrix::diagonalize(SymBlockMatrix* eigvectors, SymBlockVector* eigvalues)
{
    for (int h=0; h<nirreps_; h++) {
      if (rowspi_[h]) {
            sq_rsp(rowspi_[h], colspi_[h], matrix_[h], eigvalues->vector_[h], 1, eigvectors->matrix_[h], 1.0e-14);
      }
    }
}//

void SymBlockMatrix::cdsyev(char jobz, char uplo, SymBlockMatrix* eigvectors, SymBlockVector* eigvalues)
{
    for (int h=0; h<nirreps_; h++) {
      if (rowspi_[h]) {
	int lwork=3*rowspi_[h];
	double **work = block_matrix(nirreps_,lwork);
	memset(work[0],0.0,sizeof(double)*nirreps_*lwork);
        C_DSYEV(jobz, uplo, rowspi_[h], &(matrix_[h][0][0]), colspi_[h], eigvalues->vector_[h], &(work[h][0]), lwork);
      }
    }
}//

void SymBlockMatrix::davidson(int n_eigval, SymBlockMatrix* eigvectors, SymBlockVector* eigvalues, double cutoff, int print)
{
    for (int h=0; h<nirreps_; h++) {
      if (rowspi_[h]) {
	david(matrix_[h], rowspi_[h], n_eigval, eigvalues->vector_[h], eigvectors->matrix_[h], cutoff, print);
      }
    }
}//

void SymBlockMatrix::cdgesv(SymBlockVector* Xvec)
{
    for (int h=0; h<nirreps_; h++) {
      if (rowspi_[h]) {
	int errcod;
	int *ipiv = init_int_array(rowspi_[h]);
	memset(ipiv,0,sizeof(int)*rowspi_[h]);
	errcod=0;
	errcod = C_DGESV(rowspi_[h], 1, &(matrix_[h][0][0]), colspi_[h], &(ipiv[0]), Xvec->vector_[h], colspi_[h]);
	delete [] ipiv;
      }
    }
}//

//void flin(double **a, double *b, int in, int im, double *det)
void SymBlockMatrix::lineq_flin(SymBlockVector* Xvec, double *det)
{
    for (int h=0; h<nirreps_; h++) {
      if (rowspi_[h]) {
	flin(matrix_[h], Xvec->vector_[h], rowspi_[h], 1, det);
      }
    }
}//

//int pople(double **A, double *x, int dimen, int num_vecs, double tolerance, std::string OutFileRMR, int print_lvl)
void SymBlockMatrix::lineq_pople(SymBlockVector* Xvec, int num_vecs, double cutoff)
{
    for (int h=0; h<nirreps_; h++) {
      if (rowspi_[h]) {
	pople(matrix_[h], Xvec->vector_[h], rowspi_[h], num_vecs, cutoff, "outfile", 0);
      }
    }
}//

void SymBlockMatrix::mgs()
{
    double rmgs1,rmgs2,sum1;

    for (int h=0; h<nirreps_; h++) {

      for (int k=0; k<rowspi_[h];k++) {// loop-1
	rmgs1=0;

	for (int i=0; i<rowspi_[h];i++) {
	  rmgs1+=matrix_[h][i][k]* matrix_[h][i][k];
	}

	rmgs1=sqrt(rmgs1);

	for (int i=0; i<rowspi_[h];i++) {
	  matrix_[h][i][k]/=rmgs1;
	}

	for (int j=(k+1); j<rowspi_[h];j++) {// loop-2
	  rmgs2=0;

	  for (int i=0; i<rowspi_[h];i++) {
	    rmgs2+=matrix_[h][i][k]* matrix_[h][i][j];
	  }

	  for (int i=0; i<rowspi_[h];i++) {
	    matrix_[h][i][j]=matrix_[h][i][j]-(rmgs2*matrix_[h][i][k]) ;
	  }
	}// end 2
      }// end 1
    }

}//

void SymBlockMatrix::gs()
{
    for (int h=0; h<nirreps_; h++) {
      if (rowspi_[h] != 0 && colspi_[h] != 0 ) {
	schmidt(matrix_[h], rowspi_[h], colspi_[h], "outfile");
      }
    }
}//


void SymBlockMatrix::write(PSIO* psio, int itap, bool saveSubBlocks)
{
    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(itap)) {
        already_open = true;
    } else {
        psio->open(itap, PSIO_OPEN_OLD);
    }

    if (saveSubBlocks) {
        for (int h=0; h<nirreps_; h++) {
           std::string str(name_);
            str += " Irrep " + to_string(h);

            // Write the sub-blocks
            if (colspi_[h] > 0 && rowspi_[h] > 0)
                psio->write_entry(itap, const_cast<char*>(name_.c_str()), (char*)matrix_[h][0], sizeof(double) * colspi_[h] * rowspi_[h]);
        }
    } else {
        double **fullblock = to_block_matrix();
        // Need to know the size
        int sizer=0, sizec=0;
        for (int h=0; h<nirreps_; h++) {
            sizer += rowspi_[h];
            sizec += colspi_[h];
        }

        // Write the full block
        if (sizer > 0 && sizec > 0)
            psio->write_entry(itap, const_cast<char*>(name_.c_str()), (char*)fullblock[0], sizeof(double) * sizer * sizec);
        free_block(fullblock);
    }

    if (!already_open) psio->close(itap, 1);     // Close and keep
}//

void SymBlockMatrix::write(std::shared_ptr<psi::PSIO> psio, int itap, bool saveSubBlocks)
{
    write(psio.get(), itap, saveSubBlocks);
}//

void SymBlockMatrix::read(std::shared_ptr<psi::PSIO> psio, int itap, bool readSubBlocks)
{
    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(itap)) {
        already_open = true;
    }
    else psio->open(itap, PSIO_OPEN_OLD);

    int sizerow=0, sizecol=0;
    for (int h=0; h<nirreps_; h++) {
        sizerow += rowspi_[h];
        sizecol += colspi_[h];
    }

    if (readSubBlocks == false) {
        double **temp = block_matrix(sizerow, sizecol);
        psio->read_entry(itap, const_cast<char*>(name_.c_str()),  (char *)temp[0],  sizeof(double) * sizerow * sizecol);
        set(temp);
        free_block(temp);
    }

    else {
        for (int h=0; h<nirreps_; h++) {
            if (colspi_[h] > 0 && rowspi_[h] > 0)
             psio->read_entry(itap, const_cast<char*>(name_.c_str()),  (char *)matrix_[h][0],  sizeof(double) * rowspi_[h] * colspi_[h]);
        }
    }

    if (already_open == false) psio->close(itap, PSIO_OPEN_OLD);

}//

void SymBlockMatrix::read(std::shared_ptr<psi::PSIO> psio, int itap, const char *label, bool readSubBlocks)
{
    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(itap)) {
        already_open = true;
    }
    else psio->open(itap, PSIO_OPEN_OLD);

    int sizerow=0, sizecol=0;
    for (int h=0; h<nirreps_; h++) {
        sizerow += rowspi_[h];
        sizecol += colspi_[h];
    }

    if (readSubBlocks == false) {
        double **temp = block_matrix(sizerow, sizecol);
        psio->read_entry(itap, label,  (char *)temp[0],  sizeof(double) * sizerow * sizecol);
        set(temp);
        free_block(temp);
    }

    else {
        for (int h=0; h<nirreps_; h++) {
            if (colspi_[h] > 0 && rowspi_[h] > 0)
             psio->read_entry(itap, label,  (char *)matrix_[h][0],  sizeof(double) * rowspi_[h] * colspi_[h]);
        }
    }

    if (already_open == false) psio->close(itap, PSIO_OPEN_OLD);

}//

void SymBlockMatrix::read_oooo(std::shared_ptr<psi::PSIO> psio, int itap, int *mosym, int *qt2pitzer, int *occ_off, int *occpi, Array3i *oo_pairidx)
{
   	IWL ERIIN(psio.get(), itap, 0.0, 1, 1);
	int ilsti,nbuf,index,fi;
	double value = 0.0;

 do
 {
        ilsti = ERIIN.last_buffer();
        nbuf = ERIIN.buffer_count();

   fi = 0;
   for (int idx=0; idx < nbuf; idx++ )
   {

        int i = ERIIN.labels()[fi];
            i = abs(i);
        int j = ERIIN.labels()[fi+1];
        int k = ERIIN.labels()[fi+2];
        int l = ERIIN.labels()[fi+3];
        value = ERIIN.values()[idx];
        fi += 4;

        int i_pitzer = qt2pitzer[i];
        int j_pitzer = qt2pitzer[j];
        int k_pitzer = qt2pitzer[k];
        int l_pitzer = qt2pitzer[l];

        int hi = mosym[i_pitzer];
        int hj = mosym[j_pitzer];
        int hk = mosym[k_pitzer];
        int hl = mosym[l_pitzer];

        int hij = hi^hj;
        int hkl = hk^hl;

        if (hij == hkl) {
            int ij = oo_pairidx->get(hij, i, j);
            int kl = oo_pairidx->get(hkl, k, l);
            matrix_[hij][ij][kl] = value;
        }

   }
        if(!ilsti)
	  ERIIN.fetch();

 } while(!ilsti);
}//

void SymBlockMatrix::read_oovv(std::shared_ptr<psi::PSIO> psio, int itap, int nocc, int *mosym, int *qt2pitzer, int *occ_off, int *vir_off, int *occpi,
                               int *virpi, Array3i *oo_pairidx, Array3i *vv_pairidx)
{
   	IWL ERIIN(psio.get(), itap, 0.0, 1, 1);
	int ilsti,nbuf,index,fi;
	double value = 0.0;

 do
 {
        ilsti = ERIIN.last_buffer();
        nbuf = ERIIN.buffer_count();

   fi = 0;
   for (int idx=0; idx < nbuf; idx++ )
   {

        int i = ERIIN.labels()[fi];
            i = abs(i);
        int j = ERIIN.labels()[fi+1];
        int a = ERIIN.labels()[fi+2];
        int b = ERIIN.labels()[fi+3];
        value = ERIIN.values()[idx];
        fi += 4;

        int i_pitzer = qt2pitzer[i];
        int j_pitzer = qt2pitzer[j];
        int a_pitzer = qt2pitzer[a];
        int b_pitzer = qt2pitzer[b];

        int hi = mosym[i_pitzer];
        int hj = mosym[j_pitzer];
        int ha = mosym[a_pitzer];
        int hb = mosym[b_pitzer];

        int hij = hi^hj;
        int hab = ha^hb;

        int A = a - nocc;
        int B = b - nocc;


        if (hij == hab) {
            int ij = oo_pairidx->get(hij, i, j);
            int ab = vv_pairidx->get(hab, A, B);
            matrix_[hij][ij][ab] = value;
        }

   }
        if(!ilsti)
	  ERIIN.fetch();

 } while(!ilsti);
}//



/********************************************************************************************/
/************************** SymBlockVector **************************************************/
/********************************************************************************************/
SymBlockVector::SymBlockVector()
{
    vector_ = NULL;
    dimvec_ = NULL;
}

SymBlockVector::SymBlockVector(string name)
{
    vector_ = NULL;
    dimvec_ = NULL;
}


SymBlockVector::SymBlockVector(int nirreps, int *ins_dimvec)
{
    vector_ = NULL;
    nirreps_ = nirreps;
    dimvec_ = new int[nirreps_];
    for (int h=0; h<nirreps_; h++) {
        dimvec_[h] = ins_dimvec[h];
    }
    memalloc();
}//

SymBlockVector::SymBlockVector(string name, int nirreps, int *ins_dimvec)
{
    vector_ = NULL;
    nirreps_ = nirreps;
    name_ = name;
    dimvec_ = new int[nirreps_];
    for (int h=0; h<nirreps_; h++) {
        dimvec_[h] = ins_dimvec[h];
    }
    memalloc();
}//


SymBlockVector::~SymBlockVector()
{
    release();
    if (dimvec_) delete[] dimvec_;
}//

int *SymBlockVector::dimvec()
{
    return dimvec_;
}//

void SymBlockVector::memalloc()
{
    if (vector_) release();
    vector_ = (double**)malloc(sizeof(double**) * nirreps_);
    for (int h=0; h<nirreps_; h++) {
      if (dimvec_[h] != 0) {
	vector_[h] = new double [dimvec_[h]];
      }
      else vector_[h] = NULL;
    }
}//

void SymBlockVector::release()
{
    if (!vector_) return;
    for (int h=0; h<nirreps_; h++) {
      if (vector_[h]) free(vector_[h]);
    }
    vector_ = NULL;
}//

void SymBlockVector::zero()
{
    size_t size;
    for (int h=0; h<nirreps_; h++) {
      size = dimvec_[h] * sizeof(double);
      if (size) {
	memset(&(vector_[h][0]), 0.0, size);
      }
    }
}//

void SymBlockVector::copy(const SymBlockVector* Adum)
{
    // Make sure that vectors are in the same size
    bool same = true;
    for (int h=0; h<nirreps_; h++) {
      if (dimvec_[h] != Adum->dimvec_[h])  same = false;
    }

    if (same == false) {
        release();
        if (dimvec_) delete[] dimvec_;
        dimvec_ = new int[nirreps_];
        for (int i=0; i<nirreps_; ++i) {
            dimvec_[i] = Adum->dimvec_[i];
        }
        memalloc();
    }

    // If vectors are in the same size
    for (int h=0; h<nirreps_; h++) {
      if (dimvec_[h] != 0) {
	memcpy(&(vector_[h][0]), &(Adum->vector_[h][0]), dimvec_[h] * sizeof(double));
      }
    }
}//

void SymBlockVector::add(const SymBlockVector* Adum)
{
    double *lhs, *rhs;
    for (int h=0; h<nirreps_; h++) {
      size_t size = dimvec_[h];
      if (size) {
	lhs = vector_[h];
        rhs = Adum->vector_[h];
        for (size_t cnt=0; cnt<size; cnt++) {
	  *lhs += *rhs;
          lhs++; rhs++;
	}
      }
    }
}//

void SymBlockVector::add(int h, int i, double value)
{
  vector_[h][i]+=value;
}//

void SymBlockVector::subtract(const SymBlockVector* Adum)
{
    double *lhs, *rhs;
    for (int h=0; h<nirreps_; h++) {
      size_t size = dimvec_[h];
      if (size) {
	lhs = vector_[h];
        rhs = Adum->vector_[h];
        for (size_t cnt=0; cnt<size; cnt++) {
	  *lhs -= *rhs;
          lhs++; rhs++;
	}
      }
    }
}//

void SymBlockVector::subtract(int h, int i, double value)
{
  vector_[h][i]-=value;
}//

void SymBlockVector::scale(double a)
{
    size_t size;
    for (int h=0; h<nirreps_; h++) {
      size = dimvec_[h];
      if (size) C_DSCAL(size, a, &(vector_[h][0]), 1);
    }
}//

double SymBlockVector::sum_of_squares()
{
    double summ;
    summ=0.0;
    for (int h=0; h<nirreps_; h++) {
	for (int j=0; j<dimvec_[h]; ++j) {
	  summ += vector_[h][j] * vector_[h][j];
	}
    }
    return summ;
}//

double SymBlockVector::rms()
{
    double summ;
    int dim;
    summ=0.0;
    dim=0;

    for (int h=0; h<nirreps_; h++) {
      if (dimvec_[h] != 0) dim += dimvec_[h];
    }

    for (int h=0; h<nirreps_; h++) {
	for (int j=0; j<dimvec_[h]; ++j) {
	  summ += vector_[h][j] * vector_[h][j];
	}
    }
    summ=sqrt(summ)/dim;
    return summ;
}//

double SymBlockVector::rms(SymBlockVector* Atemp)
{
    double summ;
    int dim;
    summ=0.0;
    dim=0;

    for (int h=0; h<nirreps_; h++) {
      if (dimvec_[h] != 0) dim += dimvec_[h];
    }

    for (int h=0; h<nirreps_; h++) {
	for (int j=0; j<dimvec_[h]; ++j) {
	  summ += (vector_[h][j] * vector_[h][j]) - (Atemp->vector_[h][j] * Atemp->vector_[h][j]);
	}
    }
    summ=sqrt(summ)/dim;
    return summ;
}//

double SymBlockVector::norm()
{
    double summ;
    summ=0.0;
    for (int h=0; h<nirreps_; h++) {
	for (int j=0; j<dimvec_[h]; ++j) {
	  summ += vector_[h][j] * vector_[h][j];
	}
    }
    summ=sqrt(summ);
    return summ;
}//

void SymBlockVector::set(double value)
{
    size_t size;
    for (int h=0; h<nirreps_; h++) {
      size = dimvec_[h];
      for (size_t i=0; i<size; ++i) {
	vector_[h][i] = value;
      }
    }
}//

void SymBlockVector::set(int h, int i, double value)
{
  vector_[h][i]=value;
}//

void SymBlockVector::set(double *Avec)
{
    int offset;
    if (Avec == NULL) return;

    offset = 0;
    for (int h=0; h<nirreps_; h++) {
      for (int i=0; i<dimvec_[h]; ++i) {
	int ii = i + offset;
        vector_[h][i] = Avec[ii];
      }
      offset += dimvec_[h];
    }
}//

double SymBlockVector::get(int h, int m)
{
  return vector_[h][m];
}//

double *SymBlockVector::to_vector()
{
    int sizecol=0;
    for (int h=0; h<nirreps_; h++) {
        sizecol += dimvec_[h];
    }

    double *temp = new double[sizecol];
    int offsetcol=0;
    for (int h=0; h<nirreps_; h++) {
      for (int j=0; j<dimvec_[h]; ++j) {
	temp[j+offsetcol] = vector_[h][j];
      }
        offsetcol += dimvec_[h];
    }

    return temp;
}//

double SymBlockVector::trace()
{
    double value;
    value=0.0;
    for (int h=0; h<nirreps_; h++) {
      if (dimvec_[h] != 0) {
	for (int j=0; j<dimvec_[h]; ++j) value += vector_[h][j];
      }
    }
    return value;
}//

void SymBlockVector::print(std::string OutFileRMR)
{
   std::shared_ptr<psi::PsiOutStream> printer=(OutFileRMR=="outfile"?outfile:
         std::shared_ptr<OutFile>(new OutFile(OutFileRMR,APPEND)));
   if (name_.length()) printer->Printf( "\n ## %s ##\n", name_.c_str());
    for (int h=0; h<nirreps_; h++) {
      if (dimvec_[h] != 0) {
	printer->Printf( "\n Irrep: %d\n", h+1);
	for (int j=0; j<dimvec_[h]; ++j) {
	  printer->Printf( "%20.14f \n", vector_[h][j]);
	}
      }
    }
}//

void SymBlockVector::print()
{
     if (name_.length()) outfile->Printf( "\n ## %s ##\n", name_.c_str());
    for (int h=0; h<nirreps_; h++) {
      if (dimvec_[h] != 0) {
	outfile->Printf( "\n Irrep: %d\n", h+1);
	for (int j=0; j<dimvec_[h]; ++j) {
	  outfile->Printf( "%20.14f \n", vector_[h][j]);
	}
      }
    }

}//

void SymBlockVector::set_to_unit()
{
    size_t size;
    for (int h=0; h<nirreps_; h++) {
      size = dimvec_[h] * sizeof(double);
      if (size) {
	memset(&(vector_[h][0]), 0.0, size);
        for (int i=0; i<dimvec_[h]; i++)
	  vector_[h][i] = 1.0;
      }
    }
}//

void SymBlockVector::gemv(bool transa, double alpha, SymBlockMatrix* A, SymBlockVector* X, double beta)
{
    char trans = transa ? 't' : 'n';

    for (int h =0; h < nirreps_; ++h) {
        C_DGEMV(trans, A->rowspi_[h], A->colspi_[h], alpha, &(A->matrix_[h][0][0]),
                A->rowspi_[h], &(X->vector_[h][0]), 1, beta, &(vector_[h][0]), 1);
    }
}//

double SymBlockVector::dot(SymBlockVector* X)
{
    double tmp = 0.0;
    for (int h=0; h<nirreps_; ++h) {
        if (dimvec_[h] != X->dimvec_[h]) {
            printf("SymBlockVector::dot: Vectors are not of the same size.\n");
            return 0.0;
        }
        for (int i=0; i<dimvec_[h]; ++i) {
            tmp += vector_[h][i] * X->vector_[h][i];
        }
    }
    return tmp;
}//

/********************************************************************************************/
/********************************************************************************************/
}} // End Namespaces
