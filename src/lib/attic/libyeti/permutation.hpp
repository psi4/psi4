#ifndef yeti_permuation_hpp
#define yeti_permuation_hpp

namespace yeti {

class Permutation;
class PermutationSet;
class PermutationGroup;

//typedef boost::intrusive_ptr<Permutation> PermutationPtr;
typedef boost::intrusive_ptr<PermutationSet> PermutationSetPtr;
typedef boost::intrusive_ptr<PermutationGroup> PermutationGroupPtr;

struct permutation_less : std::binary_function<Permutation*, Permutation*, bool>
{
    bool operator() (Permutation* p, Permutation* q) const;
};


}

#endif

