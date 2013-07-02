#ifndef yeti_matrix_hpp
#define yeti_matrix_hpp

namespace yeti {

class MatrixGenerator;
class MatrixMap;
class Matrix;
class MatrixIndex;
class MatrixConfiguration;
class MatrixMapBuilder;
class MatrixLocation;

typedef boost::intrusive_ptr<Matrix> MatrixPtr;
typedef boost::intrusive_ptr<MatrixConfiguration> MatrixConfigurationPtr;
typedef boost::intrusive_ptr<MatrixMap> MatrixMapPtr;
typedef boost::intrusive_ptr<MatrixIndex> MatrixIndexPtr;
typedef boost::intrusive_ptr<MatrixGenerator> MatrixGeneratorPtr;
typedef boost::intrusive_ptr<MatrixMapBuilder> MatrixMapBuilderPtr;
typedef boost::intrusive_ptr<MatrixLocation> MatrixLocationPtr;

}

#endif

