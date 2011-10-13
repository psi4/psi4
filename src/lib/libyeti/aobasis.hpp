#ifndef yeti_aobasis_hpp
#define yeti_aobasis_hpp

namespace yeti {

class AOBasis;
class Atom;
class Shell;

typedef boost::intrusive_ptr<AOBasis> AOBasisPtr;
typedef boost::intrusive_ptr<Atom> AtomPtr;
typedef boost::intrusive_ptr<Shell> ShellPtr;

class MultiShellMap;

typedef boost::intrusive_ptr<MultiShellMap> MultiShellMapPtr;

}

#endif
