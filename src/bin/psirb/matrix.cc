#include <ruby.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include "psirb.h"

namespace psi { namespace psirb {
	
using namespace psi;

extern "C" {
	extern void free_block(double **);
	extern double ** block_matrix(unsigned long int n, unsigned long int m);
};

//
// class Matrix
VALUE Matrix::m_rbMatrix = Qnil;

Matrix::Matrix() : m_pMatrix(NULL), m_nRows(0), m_nCols(0)
{
	
}

Matrix::~Matrix()
{
	if (m_pMatrix)
		release();
	m_pMatrix = NULL;
	m_nRows = m_nCols = 0;
}

bool Matrix::allocate(size_t rows, size_t cols)
{
	if (m_pMatrix)
		release();
		
	m_pMatrix = psi::Chkpt::matrix<double>((int)rows, (int)cols);
	m_nRows = rows;
	m_nCols = cols;
	
	return true;
}

void Matrix::release()
{
	if (m_pMatrix)
		psi::Chkpt::free(m_pMatrix);
	m_pMatrix = NULL;
	m_nRows = m_nCols = 0;
}

void Matrix::copy(Matrix *c)
{
	size_t size;
	
	release();
	allocate(c->m_nRows, m_nCols);
	
	size = m_nRows * m_nCols;
	memcpy(m_pMatrix[0], c->m_pMatrix[0], sizeof(double) * size);
}

void Matrix::create_ruby_class()
{
	// Create the Matrix class under the Psi module
	Matrix::m_rbMatrix = rb_define_class_under(Globals::g_rbPsi, "Matrix", rb_cObject);
	
	// Register the allocation function with Ruby
	rb_define_alloc_func(Matrix::m_rbMatrix, Matrix::rb_alloc);
	
	// Register the initialization function Ruby
	rb_define_method(Matrix::m_rbMatrix, "initialize", RUBYCAST(Matrix::rb_init), -1);
	rb_define_method(Matrix::m_rbMatrix, "initialize_copy", RUBYCAST(Matrix::rb_init_copy), 1);
	
	// Accessors
	rb_define_method(Matrix::m_rbMatrix, "[]", RUBYCAST(Matrix::rb_element_get), 2);
	rb_define_method(Matrix::m_rbMatrix, "[]=", RUBYCAST(Matrix::rb_element_set), 3);
	
	// Conversion routines	
	rb_define_method(Matrix::m_rbMatrix, "to_s", RUBYCAST(Matrix::rb_to_s), 0);
}

void Matrix::rb_free(void *p)
{
	Matrix *pMatrix = (Matrix*)p;
	if (pMatrix)
		delete pMatrix;
}

VALUE Matrix::rb_alloc(VALUE klass)
{
	Matrix *newMatrix = new Matrix;
	VALUE newObj;
	
	// Wrap the newly created Matrix inside a Ruby object
	newObj = Data_Wrap_Struct(klass, 0, Matrix::rb_free, newMatrix);
	// Ruby is now responsible for the object, not us.
	
	return newObj;
}

VALUE Matrix::rb_init(int argc, VALUE *argv, VALUE self)
{
	Matrix *matrix;
	int i, j;
	
	// get the object
	Data_Get_Struct(self, Matrix, matrix);
	
	// Depending on the number of arguments decide what to do.
	if (argc == 0)            // Create an empty NULL matrix
		return self;
	else if (argc == 1) {     // Check to see if the user is sending an array
		if (TYPE(argv[0]) == T_ARRAY) {
			// Yes, an array
			VALUE arr = argv[0];
			int row = RARRAY(arr)->len;
			
			// Check the length of the components of the array
			int col = 1;
			for (i=0; i<row; ++i) {
				VALUE el = RARRAY(arr)->ptr[i];
				if (TYPE(el) == T_ARRAY) {
					if (col < RARRAY(el)->len)
						col = RARRAY(el)->len;
				}
			}
			
			// We have enough information to allocate the matrix.
			matrix->allocate(row, col);
			
			// Go through and save the information to the new matrix
			for (i = 0; i<row; i++) {
				VALUE el = RARRAY(arr)->ptr[i];
				
				// Check the type: array or value
				if (TYPE(el) == T_ARRAY) {
					for (j = 0; j<RARRAY(el)->len; ++j) {
						matrix->set(i, j, NUM2DBL(RARRAY(el)->ptr[j]));
					}
				}
				else {
					matrix->set(i, 0, NUM2DBL(el));
				}
			}
		}
		else rb_raise(rb_eTypeError, "not sure what you want to do");
	}
	else if (argc == 2) // Assume that it is of the form Matrix.new row, col
		matrix->allocate(NUM2UINT(argv[0]), NUM2UINT(argv[1]));
	else
		rb_raise(rb_eTypeError, "not sure what you want to do");
		
	return self;
}

VALUE Matrix::rb_init_copy(VALUE copy, VALUE orig)
{
	Matrix *o, *c;
	
	// Do not self copy
	if (copy == orig)
		return copy;
		
	// We can only copy a Matrix to a Matrix
	if (TYPE(orig) != T_DATA ||
		RDATA(orig)->dfree != (RUBY_DATA_FUNC)Matrix::rb_free) {
		rb_raise(rb_eTypeError, "wrong argument type");	
	}
	
	Data_Get_Struct(copy, Matrix, c);
	Data_Get_Struct(orig, Matrix, o);
	
	// Deep copy
	o->copy(c);
	
	return copy;
}

VALUE Matrix::rb_element_get(VALUE self, VALUE i, VALUE j)
{
	Matrix *matrix;
	Data_Get_Struct(self, Matrix, matrix);
	
	// Check the values
	if (NUM2UINT(i) >= matrix->m_nRows || NUM2UINT(j) >= matrix->m_nCols) {
		rb_raise(rb_eRangeError, "out of range: %d [0, %d] and %d [0, %d]",
			NUM2UINT(i), matrix->m_nRows, NUM2UINT(j), matrix->m_nCols);
	}
	
	return rb_float_new(matrix->get(NUM2UINT(i), NUM2UINT(j)));
}

VALUE Matrix::rb_element_set(VALUE self, VALUE i, VALUE j, VALUE val)
{
	Matrix *matrix;
	Data_Get_Struct(self, Matrix, matrix);
	
	// Check the values
	if (NUM2UINT(i) >= matrix->m_nRows || NUM2UINT(j) >= matrix->m_nCols) {
		rb_raise(rb_eRangeError, "out of range: %d [0, %d] and %d [0, %d]",
			NUM2UINT(i), matrix->m_nRows, NUM2UINT(j), matrix->m_nCols);
	}

	matrix->set(NUM2UINT(i), NUM2UINT(j), NUM2DBL(val));
	return self;
}

#define PRINT_MATRIX_COLS 8
VALUE Matrix::rb_to_s(VALUE self)
{
	// Convert the matrix data to a std::string to be converted to be
	// converted to a Ruby string object.
	Matrix *matrix;
	std::string text;
	std::stringstream stream;
	size_t m, n, col=0, row, i, j, extra;
	double **mat;
	
	Data_Get_Struct(self, Matrix, matrix);
	m = matrix->m_nRows;
	n = matrix->m_nCols;
	mat = matrix->m_pMatrix;
	
	while (col + PRINT_MATRIX_COLS < n) {
		// Print column headers
		stream.str(std::string());
		stream << "\n ";
		for (i=1; i <= PRINT_MATRIX_COLS; ++i) {
			stream << "           " << std::setw(4) << (i+col);
		}
		stream << "\n\n";
		text += stream.str();
		
		// Now the rows of data
		stream.str(std::string());
		for (row=0; row<m; row++) {
			// Put row label
			stream << std::setw(4) << (row+1);
			for (i=0; i<PRINT_MATRIX_COLS; ++i)
				stream << " " << std::fixed << std::setw(14) << std::setprecision(7) << mat[row][col+i];
			stream << "\n";
		}
		text += stream.str();
		col += PRINT_MATRIX_COLS;
	}
	
	// Take care of the left overs
	if (col < n) {
		extra = n-col;
		stream.str(std::string());
		stream << "\n ";
		for (i=0; i < extra; ++i) {
			stream << "           " << std::setw(4) << (i+col+1);
		}
		stream << "\n\n";
		text += stream.str();
		
		// Now the rows of data
		stream.str(std::string());
		for (row=0; row<m; row++) {
			// Put row label
			stream << std::setw(4) << (row+1);
			for (i=0; i<extra; ++i)
				stream << " " << std::fixed << std::setw(14) << std::setprecision(7) << mat[row][col+i];
			stream << "\n";
		}
		text += stream.str();
	}
	text += "\n";
	
	// Convert to Ruby object
	return rb_str_new2(text.c_str());
}

}} // namespace psi::psirb
