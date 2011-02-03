#ifndef yeti_permuation_hpp
#define yeti_permuation_hpp

namespace yeti {

class Permutation;
class PermutationSet;
class PermutationGroup;

typedef boost::intrusive_ptr<Permutation> PermutationPtr;
typedef boost::intrusive_ptr<const Permutation> constPermutationPtr;
typedef boost::intrusive_ptr<PermutationSet> PermutationSetPtr;
typedef boost::intrusive_ptr<const PermutationSet> constPermutationSetPtr;
typedef boost::intrusive_ptr<PermutationGroup> PermutationGroupPtr;
typedef boost::intrusive_ptr<const PermutationGroup> constPermutationGroupPtr;

struct permutation_less : std::binary_function<PermutationPtr, PermutationPtr, bool>
{
    bool operator() (const PermutationPtr& p, const PermutationPtr& q) const;
};


}

#endif

