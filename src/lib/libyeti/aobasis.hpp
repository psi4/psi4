#ifndef yeti_aobasis_hpp
#define yeti_aobasis_hpp

namespace yeti {

class AOBasis;
class AtomBasis;
class ShellBasis;
class AtomShell;

typedef boost::intrusive_ptr<AOBasis> AOBasisPtr;
typedef boost::intrusive_ptr<AtomBasis> AtomBasisPtr;
typedef boost::intrusive_ptr<ShellBasis> ShellBasisPtr;
typedef boost::intrusive_ptr<AtomShell> AtomShellPtr;

}

#endif
