#ifndef yeti_matrix_hpp
#define yeti_matrix_hpp

#define OLD_MATRIX 0

namespace yeti {

class MatrixGenerator;

#if OLD_MATRIX
class MatrixBlockMap;
class MatrixRow;
#else
class MatrixMap;
#endif

class Matrix;
class MatrixIndex;
class MatrixConfiguration;


typedef boost::intrusive_ptr<Matrix> MatrixPtr;
typedef boost::intrusive_ptr<const Matrix> constMatrixPtr;
typedef boost::intrusive_ptr<MatrixConfiguration> MatrixConfigurationPtr;
typedef boost::intrusive_ptr<const MatrixConfiguration> constMatrixConfigurationPtr;
#if OLD_MATRIX
typedef boost::intrusive_ptr<MatrixRow> MatrixRowPtr;
typedef boost::intrusive_ptr<const MatrixRow> constMatrixRowPtr;
typedef boost::intrusive_ptr<MatrixBlockMap> MatrixBlockMapPtr;
typedef boost::intrusive_ptr<const MatrixBlockMap> constMatrixBlockMapPtr;
#else
typedef boost::intrusive_ptr<MatrixMap> MatrixMapPtr;
#endif
typedef boost::intrusive_ptr<MatrixIndex> MatrixIndexPtr;
typedef boost::intrusive_ptr<const MatrixIndex> constMatrixIndexPtr;
typedef boost::intrusive_ptr<MatrixGenerator> MatrixGeneratorPtr;


}

#endif

