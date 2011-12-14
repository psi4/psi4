#ifndef yeti_aobasis_hpp
#define yeti_aobasis_hpp

namespace yeti {

class AOBasis;
class Atom;
class Shell;
class PartitioningPolicy;

typedef boost::intrusive_ptr<AOBasis> AOBasisPtr;
typedef boost::intrusive_ptr<Atom> AtomPtr;
typedef boost::intrusive_ptr<Shell> ShellPtr;
typedef boost::intrusive_ptr<PartitioningPolicy> PartitioningPolicyPtr;

class MultiShellMap;

typedef boost::intrusive_ptr<MultiShellMap> MultiShellMapPtr;

}

#endif
