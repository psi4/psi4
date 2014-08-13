#include <stdio.h>
#include "DADFTemplate.h"
#include "DADFAxisInfo.h"

/** This is our implementation of DistArrayTemplate.  It is kept
    private within the DistArrayDescriptorFactory.  Only
    DistArrayTemplate is exposed to the outside. 

    $Id: DADFTemplate.cxx,v 1.1 2003-08-01 00:41:53 manoj Exp $
 */

/******************************************************************************
 * Constructors and destructors
 *****************************************************************************/

/** Basic constructor.  

    Probably should not be used, since we prefer to force the user to
    provide a name for the object. 
*/
//  DADFTemplate::DADFTemplate() {

//    _name = "_UNNAMED";
//    _frozen = false;
//    _rank = -1; 
//    _volume = 1; // Allows use of *=
//    _volDefined = 0;
//  }

/** Constructor sets template name. */
DADFTemplate::DADFTemplate( const std::string name) {
  _name = name;
  _frozen = false;
  _rank = -1; 
  _volume = 1; // Allows use of *=
  _volDefined = 0;
}

/** Construct a new template as a copy of an old one.  Of course a
    name must be provided for the new one.
*/
DADFTemplate::DADFTemplate( const std::string name, DADFTemplate & original)
  : _rank( original._rank ), 
  _volume( original._volume ),
  _volDefined( original._volDefined ),
  _lowerBounds( original._lowerBounds ),
  _upperBounds( original._upperBounds ),
  _topology( original._topology ),
  _dist( original._dist )
{
  // Use name in argument instead of copied from original
  _name = name;

  // Regardless of whether original was frozen, ths one should not be
  _frozen = false;

  /** Setting up _axisInfo is a little complicated.  We expect it to
      be of size original._rank unless original hasn't had setRank()
      called on it yet.  This should take care of things.
  */
  _axisInfo.resize( original._axisInfo.size() );

  unsigned int i;
  BlockAxisInfo * bai;       // Used for casting to copy _axisInfo
  GenBlockAxisInfo * gbai;   // Used for casting to copy _axisInfo
  ImplicitAxisInfo * iai;    // Used for casting to copy _axisInfo

  // Regular distributions are per-axis
  for ( i=0; i < original._axisInfo.size(); ++i ) {
    if ( original._axisInfo[i] ) {
      switch ( (original._axisInfo[i])->getDistType() ) {
      case Collapsed:
	// Nothing to copy, so no need for a copy constructor here
	_axisInfo[i] = new CollapsedAxisInfo; 
	break;
      case Block:
	bai = dynamic_cast<BlockAxisInfo *>(original._axisInfo[i]);
	_axisInfo[i] = new BlockAxisInfo( *bai );
	break;
      case GenBlock:
	gbai = dynamic_cast<GenBlockAxisInfo *>(original._axisInfo[i]);
	_axisInfo[i] = new GenBlockAxisInfo( *gbai );
	break;
      case Implicit:
	iai = dynamic_cast<ImplicitAxisInfo *>(original._axisInfo[i]);
	_axisInfo[i] = new ImplicitAxisInfo( *iai );
	break;
      }
    } else { 
      _axisInfo[i] = 0; 
    }
  }

  /** _regionList must be treated with similar care to _axisInfo. */

  DADFRegionInfo * dri;      // Used to copy _regionList
  std::list<DADFRegionInfo *>::iterator driter;

  for ( driter=original._regionList.begin(); 
  	driter != original._regionList.end(); ++driter ) {

    dri = new DADFRegionInfo( *(*driter) );
    _regionList.push_back( dri );
  }

}

/** Construct a new template as an identical copy of an old one. */
DADFTemplate::DADFTemplate( DADFTemplate & original)
  : _name( original._name ),
    _rank( original._rank ), 
    _volume( original._volume ),
    _volDefined( original._volDefined ),
    _lowerBounds( original._lowerBounds ),
    _upperBounds( original._upperBounds ),
    _topology( original._topology ),
    _dist( original._dist ),
    _frozen( original._frozen )
{
  /** Setting up _axisInfo is a little complicated.  We expect it to
      be of size original._rank unless original hasn't had setRank()
      called on it yet.  This should take care of things.
  */
  _axisInfo.resize( original._axisInfo.size() );

  unsigned int i;
  BlockAxisInfo * bai;       // Used for casting to copy _axisInfo
  GenBlockAxisInfo * gbai;   // Used for casting to copy _axisInfo
  ImplicitAxisInfo * iai;    // Used for casting to copy _axisInfo

  // Regular distributions are per-axis
  for ( i=0; i < original._axisInfo.size(); ++i ) {
    if ( original._axisInfo[i] ) {
      switch ( (original._axisInfo[i])->getDistType() ) {
      case Collapsed:
	// Nothing to copy, so no need for a copy constructor here
	_axisInfo[i] = new CollapsedAxisInfo; 
	break;
      case Block:
	bai = dynamic_cast<BlockAxisInfo *>(original._axisInfo[i]);
	_axisInfo[i] = new BlockAxisInfo( *bai );
	break;
      case GenBlock:
	gbai = dynamic_cast<GenBlockAxisInfo *>(original._axisInfo[i]);
	_axisInfo[i] = new GenBlockAxisInfo( *gbai );
	break;
      case Implicit:
	iai = dynamic_cast<ImplicitAxisInfo *>(original._axisInfo[i]);
	_axisInfo[i] = new ImplicitAxisInfo( *iai );
	break;
      }
    } else { 
      _axisInfo[i] = 0; 
    }
  }

  /** _regionList must be treated with similar care to _axisInfo. */

  DADFRegionInfo * dri;      // Used to copy _regionList
  std::list<DADFRegionInfo *>::iterator driter;

  for ( driter=original._regionList.begin(); 
  	driter != original._regionList.end(); ++driter ) {

    dri = new DADFRegionInfo( *(*driter) );
    _regionList.push_back( dri );
  }

}

DADFTemplate::~DADFTemplate() {
  std::vector<DADFAxisInfo *>::iterator iter;
  std::list<DADFRegionInfo *>::iterator driter;

  // Clean up anything in _axisInfo
  for ( iter = _axisInfo.begin(); iter != _axisInfo.end(); ++iter ) {
    if ( *iter != 0 ) { delete *iter; }
  }

  // Clean up _regionList
  for ( driter= _regionList.begin(); driter != _regionList.end();
	++driter ) {
    delete *driter;
  }
  _regionList.resize(0);
}

/******************************************************************************
 * Define the template
 *****************************************************************************/

/** Name associated with this distribution. */
int DADFTemplate::setName(const std::string name) {
  // Insure we haven't been commit()ed
  if ( _frozen ) { return -10; }

  _name = name;
  return 0;
}

/** Set rank (number of dimensions) of distribution template. */
int DADFTemplate::setRank(const int rank) {
  int i;
  bool modifying = false; // Are we modifying the clone of an existing templ?

  // Insure we haven't been commit()ed
  if ( _frozen ) { return -10; }

  // Check for invalid rank
  if ( rank < 1 ) { return -3; }

  /** If we're modifying a cloned template, we want to do some nice
      things, like resize all the arrays of size "rank" while
      preserving as much of the original information as possible.  To
      do this, we depend on the constructor to initialize _rank to an
      invalid value */
  if ( _rank > 0 ) { modifying = true; }

  // Rank defines the size of many arrays, so allocate them now.

  _lowerBounds.resize( rank );
  _upperBounds.resize( rank );
  _topology.resize( rank );
  _dist.resize( rank );
  _axisInfo.resize( rank );
  
  // Finish initializing anything new.  Do something mostly harmless.
  for ( i = (modifying)? _rank : 0 ; i < rank; ++i ) {
    _lowerBounds[ i ] = 0;
    _upperBounds[ i ] = -1;
    _topology[ i ]    = 1;
    _dist[ i ]        = Collapsed;
    _axisInfo[ i ]    = 0; 
  }

  /** If the rank is changing, kill any explicit regions that may have
      been defined. 
  */

  if ( modifying && rank != _rank ) {
    std::list<DADFRegionInfo *>::iterator driter;

    for ( driter= _regionList.begin(); 
	  driter != _regionList.end(); ++driter ) {
      delete *driter;
    }
    _regionList.resize(0);
  }

  /** Don't forget to store the rank itself! 
      Intentionally saved for last, since above we use the fact that
      _rank holds the original length if we're modifying.
  */
  _rank = rank;

  return 0;
}

/** Set the global upper and lower bounds of the array. */
int DADFTemplate::setGlobalBounds(const int lower[], const int upper[]) {
  int i;

  // Insure we haven't been commit()ed
  if ( _frozen ) { return -10; }

  // Sanity check our state
  if ( _rank < 1 ) { return -8; }

  // Check bounds for validity and compute volume
  for ( i=0; i < _rank; ++i ) {
    if ( lower[i] > upper[i] ) { return -1; }

    _volume *= ( upper[i]-lower[i]+1 );
  }

  // Copy arguments into our data structures
  for ( i=0; i < _rank; ++i ) {
    _lowerBounds[ i ] = lower[ i ];
    _upperBounds[ i ] = upper[ i ];
  }

  return 0;
}

/** Sets process topology. */
int DADFTemplate::setProcTopology(const int topology[] ) {
  int i;

  // Insure we haven't been commit()ed
  if ( _frozen ) { return -10; }

  // Sanity check our state
  if ( _rank < 1 ) { return -8; }

  // Check bounds for validity
  for ( i=0; i < _rank; ++i ) { if ( topology[i] < 1 ) { return -6; } }

  // Copy arguments into our data structure
  for ( i=0; i < _rank; ++i ) { _topology[ i ] = topology[ i ]; }

  return 0;
}

/** Sets distribution type on each axis. */
int DADFTemplate::setDistType(const enum DistType dist[] ) {
  int i;
  //static const std::string distLabels[5] = { "Coll", "Bloc", "GenB", "Impl", "Expl" };

  // Insure we haven't been commit()ed
  if ( _frozen ) { return -10; }

  // Sanity check our state
  if ( _rank < 1 ) { return -8; }

  // Check input for validity: if one is explicit, all must be
  bool oneExpl = false;
  bool allExpl = true;
  for ( i=0; i < _rank; ++i) {
    oneExpl = ( oneExpl || dist[i] == Explicit );
    allExpl = ( allExpl && dist[i] == Explicit );
  }
  if ( oneExpl && ! allExpl ) { return -2; }

  // Copy arguments into our data structure, make sure axisInfo is consistent
  for ( i=0; i < _rank; ++i ) { 
    _dist[ i ] = dist[ i ]; 

    /** If we're modifying a cloned object, we only want to preserve
	those axisInfos that are consistent 
    */
    if ( _axisInfo[i] && (_axisInfo[i])->getDistType() != _dist[i] ) {
      delete _axisInfo[i];
      _axisInfo[i] = 0;
    }
  }

  return 0;
}

/** Set distribution parameters for an axis with a regular distributions. */
int DADFTemplate::setDistParameters(int axis, int blockSize,
			      int first) {
  // Insure we haven't been commit()ed
  if ( _frozen ) { return -10; }

  // Sanity check inputs
  if ( axis < 0 || axis >= _rank ) { return -4; }
  if ( _dist[ axis ] != Block ) { return -2; }
  if ( blockSize < 1 || blockSize > (_upperBounds[ axis ] -
				     _lowerBounds[ axis ] + 1) ) {
    return -5; 
  }
  if ( first < 0 || first >= _topology[ axis ] ) { return -6; }

  // Create axis info object

  _axisInfo[ axis ] = new BlockAxisInfo;

  BlockAxisInfo* bai = dynamic_cast<BlockAxisInfo *>(_axisInfo[axis]);

  bai->setDistParameters( blockSize, first);

  return 0;
}

/** Set distribution parameters for a GenBlock axis. */
int DADFTemplate::setGenBlock(int axis, int blockSizes[]) {
  int i;
  int total = 0;

  // Insure we haven't been commit()ed
  if ( _frozen ) { return -10; }

  // Sanity check inputs
  if ( axis < 0 || axis >= _rank ) { return -4; }
  if ( _dist[ axis ] != GenBlock ) { return -2; }

  for ( i=0; i < _topology[ axis ]; ++i ) {
    if ( blockSizes[i] < 0 ) { return -5; }
    total += blockSizes[i];
  }
  if ( total > (_upperBounds[ axis ] - _lowerBounds[ axis ] + 1) ) {
    return -5; 
  }

  _axisInfo[ axis ] = new GenBlockAxisInfo(  _topology[axis] );

  GenBlockAxisInfo* gbai = dynamic_cast<GenBlockAxisInfo *>(_axisInfo[axis]);
  gbai->setDistParameters( blockSizes );

  return 0;
}

/** Set distribution parameters for an Implicit axis. */
int DADFTemplate::setImplicitMap(int axis, int map[]) {
  int i;

  // Insure we haven't been commit()ed
  if ( _frozen ) { return -10; }

  // Sanity check inputs
  if ( axis < 0 || axis >= _rank ) { return -4; }
  if ( _dist[ axis ] != Implicit ) { return -2; }

  for ( i=0; i < (_upperBounds[ axis ] - _lowerBounds[ axis ] + 1) ; ++i ) {
    if ( map[i] < 0 || map[i] >= _topology[axis] ) { return -6; }
  }

  _axisInfo[ axis ] = new ImplicitAxisInfo( (_upperBounds[ axis ] -
					     _lowerBounds[ axis ] + 1)
					    );

  ImplicitAxisInfo* iai = dynamic_cast<ImplicitAxisInfo *>(_axisInfo[axis]);
  iai->setDistParameters( map );

  return 0;
}

/** Add a region to an Explicit distribution. */
int DADFTemplate::addExplicitRegion(int lower[], int upper[]) {
  int i;
  int regVolume = 1;

  // Insure we haven't been commit()ed
  if ( _frozen ) { return -10; }

  // Sanity check inputs
  for ( i=0; i < _rank; ++i ) {
    if ( lower[i] < _lowerBounds[i] ) { return  -1; }
    if ( upper[i] > _upperBounds[i] ) { return -1; }
    regVolume *= upper[i] - lower[i] + 1;
  }

  // Check for overlaps.  Note: strictly local to this processr for now
  std::list<DADFRegionInfo *>::iterator driter;
    
  driter= _regionList.begin();
  for ( _regionList.begin(); driter != _regionList.end(); ++driter ) {
    if ( (*driter)->compareBounds( lower, upper) != 0  ) { break; } 
  }
  // If we reached the end, there's no overlaps, otherwise there are!
  if ( driter != _regionList.end() ) { return -7; }

  // If everything is okay, add this region
  DADFRegionInfo * dri = new DADFRegionInfo( _rank );
  dri->setBounds( lower, upper );
  _regionList.push_back( dri );

  // Eventually in parallel: Check for completeness to give proper return
  _volDefined += regVolume;
  return 0;
}

/** Signal that template is completely defined. */
int DADFTemplate::commit() {
  // Insure we haven't been commit()ed already
  if ( _frozen ) { 
    return -10; 
  } else {
    _frozen = true;
    return 0;
  }

  // Perform global consistency checks

  // Should test topology if dist is explicit
  // Test if area of explicit regions == area of template
}

/******************************************************************************
 * Query the template
 *****************************************************************************/

/** Name associated with this distribution. (default value is "_UNNAMED") */
std::string DADFTemplate::getName() {
  return _name;
}

/** Get rank (number of dimensions) of distributed object. */
int DADFTemplate::getRank() {
  return _rank;
}

/** The global upper and lower bounds of the array.

    @param lower (Out) array of global lower bounds of array
    @param upper (Out) array of global upper bounds of array
*/
int DADFTemplate::getGlobalBounds(int lower[], int upper[]) {
  int i;

  // Sanity check our state
  if ( _rank < 1 ) { return -8; }

  // Copy arguments out of our data structures
  for ( i=0; i < _rank; ++i ) {
    lower[ i ] = _lowerBounds[ i ];
    upper[ i ] = _upperBounds[ i ];
  }

  return 0;
}

/** Returns process topology. */
int DADFTemplate::getProcTopology(int topology[] ) {
  int i;

  // Sanity check our state
  if ( _rank < 1 ) { return -8; }

  // Copy arguments out of our data structure
  for ( i=0; i < _rank; ++i ) { topology[ i ] = _topology[ i ]; }

  return 0;
}

/** Set distribution parameters for an axis with a regular distributions. */
int DADFTemplate::getDistParameters(int axis, int blockSize,
			      int first) {
  // Sanity check our state
  if ( _rank < 1 ) { return -8; }

  // Sanity check inputs
  if ( axis < 0 || axis >= _rank ) { return -4; }
  if ( _dist[ axis ] != Block ) { return -2; }

  // Get axis info object

  BlockAxisInfo* bai = dynamic_cast<BlockAxisInfo *>(_axisInfo[axis]);

  bai->getDistParameters( blockSize, first);

  return 0;
}

/** Set distribution parameters for a GenBlock axis. */
int DADFTemplate::getGenBlock(int axis, int blockSizes[]) {
  // Sanity check our state
  if ( _rank < 1 ) { return -8; }

  // Sanity check inputs
  if ( axis < 0 || axis >= _rank ) { return -4; }
  if ( _dist[ axis ] != GenBlock ) { return -2; }

  // Get axis info object

  GenBlockAxisInfo* gbai = dynamic_cast<GenBlockAxisInfo *>(_axisInfo[axis]);
  gbai->getDistParameters( blockSizes );

  return 0;
}

/** Set distribution parameters for an Implicit axis. */
int DADFTemplate::getImplicitMap(int axis, int map[]) {
  // Sanity check our state
  if ( _rank < 1 ) { return -8; }

  // Sanity check inputs
  if ( axis < 0 || axis >= _rank ) { return -4; }
  if ( _dist[ axis ] != Implicit ) { return -2; }

  // Get axis info object

  ImplicitAxisInfo* iai = dynamic_cast<ImplicitAxisInfo *>(_axisInfo[axis]);
  iai->getDistParameters( map );

  return 0;
}

/** Returns distribution type on each axis. */
int DADFTemplate::getDistType(enum DistType dist[] ) {
  int i;

  // Sanity check our state
  if ( _rank < 1 ) { return -8; }

  // Copy our data structure to the argument
  for ( i=0; i < _rank; ++i ) { dist[ i ] = _dist[ i ]; }

  return 0;
}

bool DADFTemplate::isDefined() {
    if ( _frozen ) {
      return true;
    } else {
      return false;
    }
}

/** Print the contents of the template (for debugging) */
void DADFTemplate::printTemplate() {
  int i;
  static const std::string distLabels[5] 
    = { "Coll", "Bloc", "GenB", "Impl", "Expl" };

  cerr << "Array distribution template `" << _name << "' rank " <<
    _rank << ((_frozen)? " (" : " (not" ) << " frozen)" << endl; 

  if ( _rank < 1 ) { return; }

  cerr << "      Bounds:             " ;
  for ( i=0; i < _rank-1 ; ++i ) {
    cerr << _lowerBounds[i] << ":" << _upperBounds[i] << ", ";
  }
  cerr << _lowerBounds[i] << ":" << _upperBounds[i] << endl;
  
  cerr << "      Process topology:   " ;
  for ( i=0; i < _rank-1 ; ++i ) {
    cerr << _topology[i] << ", ";
  }
  cerr << _topology[i] << endl;
  
  cerr << "      Distribution Types: " ;
  for ( i=0; i < _rank-1 ; ++i ) {
    cerr << distLabels[ _dist[i] ] << ", ";
  }
  cerr << distLabels[ _dist[i] ] << endl;
  
  // Print out axiswise details
  for ( i=0; i < _rank; ++i ) {
    if ( _axisInfo[ i ] ) {
      cerr << "      Axis " << i << " " << 
	distLabels[ (_axisInfo[ i ])->getDistType() ];
      
      cerr << " ";
      (_axisInfo[i])->printAxisInfo();
      cerr << endl;
	
    } else {
      cerr << "      Axis " << i << " (no info supplied)" << endl;
    }
  }

  // Print out the region list
  std::list<DADFRegionInfo *>::iterator driter;

  cerr << "      Regions registered: " << _regionList.size() << endl; 

  for ( driter= _regionList.begin(); 
  	driter != _regionList.end(); ++driter ) {
    cerr << "      Region ";
    (*driter)->printRegionInfo();
    cerr << endl;
  }
}
