#ifndef yeti_diis_hpp
#define yeti_diis_hpp

namespace yeti {

class DiisParameters;
class DiisExtrapolation;

typedef boost::intrusive_ptr<DiisExtrapolation> DiisExtrapolationPtr;
typedef boost::intrusive_ptr<DiisParameters> DiisParametersPtr;

}

#endif

