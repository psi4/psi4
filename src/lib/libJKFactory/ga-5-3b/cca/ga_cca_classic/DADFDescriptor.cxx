#include "DADFDescriptor.h"
#include "DistArrayTemplate.h"

/** This is our implementation of DistArrayDescriptor.  It is kept
    private within the DistArrayDescriptorFactory.  Only
    DistArrayDescriptor is exposed to the outside. 

    $Id: DADFDescriptor.cxx,v 1.1 2003-08-01 00:41:53 manoj Exp $
 */

/******************************************************************************
 * Constructors and destructors
 *****************************************************************************/

/** Constructor sets descriptor name. */
DADFDescriptor::DADFDescriptor( const std::string name) {
  _name = name;
  _frozen = false;
  _rank = -1;
  _type = stv_Int;
  _templ = 0;
  _isExplicitDist = false;
}

/** Construct a new descriptor as a copy of an old one.
*/
DADFDescriptor::DADFDescriptor(const std::string name,  DADFDescriptor & original)
  : _rank( original._rank ), 
  _type( original._type),
  _lowerBounds( original._lowerBounds ),
  _upperBounds( original._upperBounds ),
  _topology( original._topology ),
  _procCoords( original._procCoords ),
  _isExplicitDist( original._isExplicitDist )
{
  // Use name provided by user rather than from original
  _name = name;

  // Regardless of whether original was frozen, ths one should not be
  _frozen = false;

  // Duplicate template
  DADFTemplate * originalTemplDADF = 
    dynamic_cast<DADFTemplate *>(original._templ);
  _templ = new DADFTemplate( *originalTemplDADF );

  /** _regionList must be treated with similar care to _axisInfo. */

  DADFRegionInfo * dri;      // Used to copy _regionList
  std::list<DADFRegionInfo *>::iterator driter;

  for ( driter=original._regionList.begin(); 
    	driter != original._regionList.end(); ++driter ) {

    dri = new DADFRegionInfo( *(*driter) );
    _regionList.push_back( dri );
  }

}

  DADFDescriptor::~DADFDescriptor() {
    delete _templ;

    std::list<DADFRegionInfo *>::iterator driter;
    
    // Clean up _regionList
    for ( driter= _regionList.begin(); driter != _regionList.end();
	  ++driter ) {
      delete *driter;
    }
    _regionList.resize(0);
}

/******************************************************************************
 * Define the descriptor
 *****************************************************************************/

/** Set data type. */
int DADFDescriptor::setDataType(const enum DataType type) {
  // Insure we haven't been commit()ed
  if ( _frozen ) { return -11; }

  _type = type;

  return 0;
}

/** Associate this data object with a distribution template. */
int DADFDescriptor::setTemplate(DistArrayTemplate * & templ) {
  // Insure we haven't been commit()ed
  if ( _frozen ) { return -11; }

  // If there's already one there, replace it
  if ( _templ != 0 ) { delete _templ; }

  DADFTemplate * templDADF = dynamic_cast<DADFTemplate *>(templ);
  if ( templDADF == 0 ) {
    cerr << "DADFDescriptor:setTemplate: " <<
      "Template is of the wrong type." << endl;
    return -13;
  }

  // Make a private copy of the template
  DADFTemplate * newtempl = new DADFTemplate( *templDADF );
  _templ = dynamic_cast<DistArrayTemplate *>(newtempl);

  return 0;
}

/** Sets this process's location in the process topology. */
int DADFDescriptor::setMyProcCoords(const int procCoords[] ) {
  int i;

  // Insure we haven't been commit()ed
  if ( _frozen ) { return -11; }

  // Check bounds for validity
  for ( i=0; i < _rank; ++i ) { 
    if ( procCoords[i] < 0 || procCoords[i] > _topology[i] ) { 
      return -6;
    }
  }

  // Copy arguments into our data structure
  for ( i=0; i < _rank; ++i ) { _procCoords[ i ] = procCoords[ i ]; }

  return 0;
}

/** Align object to template with identity mapping. */
int DADFDescriptor::setIdentityAlignmentMap() {
  int i, err;
  int *proc;
  int *lower, *upper; 
  DistArrayTemplate::DistType *dist;

  // Insure we haven't been commit()ed
  if ( _frozen ) { return -11; }

  // Insure we have a template
  if ( _templ == 0 ) { return -14; }

  /** At the moment, this simply means we're allowed to extract stuff
      from the template more or less with impunity. */
  _rank = _templ->getRank();

  // Resize our internal vectors
  _lowerBounds.resize( _rank );
  _upperBounds.resize( _rank );
  _topology.resize( _rank );
  _procCoords.resize( _rank );

  proc = new int[ _rank ];
  err = _templ->getProcTopology( &(*proc) );
  if ( err != 0 ) { return err; }

  for ( i=0; i < _rank; ++i ) {
    _topology[i] = proc[i];
  }
  delete [] proc;
  
  lower = new int[ _rank ];
  upper = new int[ _rank ];
  err = _templ->getGlobalBounds( &(*lower), &(*upper) );
  if ( err != 0 ) { return err; }

  for ( i=0; i < _rank; ++i ) {
    _lowerBounds[i] = lower[i];
    _upperBounds[i] = upper[i];
  }
  
  delete [] lower;
  delete [] upper;

  dist = new DistArrayTemplate::DistType[ _rank ];
  err = _templ->getDistType( &(*dist) );
  if ( err != 0 ) { return err; }

  _isExplicitDist = ( dist[0] == DistArrayTemplate::Explicit );

  delete [] dist;

  return 0;
}

/** Set pointer the local piece of the data object. */
int DADFDescriptor::setLocalDataPointer(void* data, const int
			 strides[]) {
  int i;

  // Insure we haven't been commit()ed
  if ( _frozen ) { return -11; }

  // This call is only valid for non-explicit distributions
  if ( _isExplicitDist ) { return -2; }

  /** Sanity check the strides. Until someone come up with a use case
      to demonstrate why we need non-positive strides, we'll say
      they're bad.
  */
  for (i=0; i < _rank; ++i) {
    if ( strides[i] < 1 ) { return -16; }
  }
  
  // Create a region info and set the point & stride (but not the bounds)
  DADFRegionInfo * dri = new DADFRegionInfo( _rank );
  dri->setDataLocation( data, strides );
  _regionList.push_back( dri );

  return 0;
}

/** Set pointer for a region of an explicitly distributed data object. */
int DADFDescriptor::setRegionDataPointer(const int lower[], const int
			 upper[], void* data, const int
			 strides[]) {
  int i;

  // Insure we haven't been commit()ed
  if ( _frozen ) { return -11; }

  // This call is only valid for explicit distributions
  if ( ! _isExplicitDist ) { return -2; }

  /** Sanity check the strides. Until someone come up with a use case
      to demonstrate why we need non-positive strides, we'll say
      they're bad.
  */
  for (i=0; i < _rank; ++i) {
    if ( strides[i] < 1 ) { return -16; }
  }
  
  /** This is incorrect -- it checks our region list, not the
      template's.  At the moment we don't have query functions on the
      template.  Alternatively, we could initialize our region list
      from the template's when we make the alignment map, but once
      again we really need the template queries to do that.  So for
      now, we just punt, and hope the user is careful!
  */
//    // Check that specified region matches one registered
//    std::list<DADFRegionInfo *>::iterator driter;
    
//    driter= _regionList.begin();
//    for ( _regionList.begin(); driter != _regionList.end(); ++driter ) {
//      cerr << " compare returns " << 
//        (*driter)->compareBounds( lower, upper) << endl;
//      if ( (*driter)->compareBounds( lower, upper) == 1  ) { break; } 
//    }
//    // If we reached the end, there's no match
//    if ( driter == _regionList.end() ) { return -101; }

  // Create a region info and set everything
  DADFRegionInfo * dri = new DADFRegionInfo( _rank );
  dri->setBounds( lower, upper );
  dri->setDataLocation( data, strides );
  _regionList.push_back( dri );

  return 0;
}

/** Signal that descriptor is completely defined. */
int DADFDescriptor::commit() {
  // Insure we haven't been commit()ed already
  if ( _frozen ) { 
    return -11; 
  } else {
    _frozen = true;
    return 0;
  }

  // Perform global consistency checks

  // Should test topology if dist is explicit
  // Test if area of explicit regions == area of descriptor
}

/******************************************************************************
 * Query the descriptor
 *****************************************************************************/

std::string DADFDescriptor::getName() {
  return _name;
}

bool DADFDescriptor::isDefined() {
    if ( _frozen ) {
      return true;
    } else {
      return false;
    }
}

/** Get data type. */
DistArrayDescriptor::DataType DADFDescriptor::getDataType() {
  return _type;
}

/** Return pointer to distribution template associated with this
    data object. 
*/
DistArrayTemplate * DADFDescriptor::getTemplate() {
  return _templ;
}

/** Get this process's location in the process topology. */
int DADFDescriptor::getMyProcCoords(int procCoords[] ) {
  int i;

  // Sanity check our state
  if ( _rank < 1 ) { return -8; }

  // Copy arguments out of our data structures
  for ( i=0; i < _rank; ++i ) { procCoords[ i ] = _procCoords[ i ]; }

  return 0;
}

// Kludges to check GenBlock and Explicit descriptors.
int DADFDescriptor::getNumLocalRegions() {
  //  if ( ! _frozen ) { return -100; }

  if ( _isExplicitDist ) {
    return _regionList.size();
  } else {
    return 1; 
  }
}
// Kludges to check GenBlock and Explicit descriptors.
int DADFDescriptor::getLocalRegionInfo(int region, int lower[], int upper[],
				       void * & data, int strides[]) {
  int i;

  //  if ( ! _frozen ) { return -100; }

  if ( _isExplicitDist ) {
    if ( region >= (int)_regionList.size() ) { return -15; }
  } else {
    if ( region != 0 ) { return -15; }
  }
    
  // Find the right region (assumes list not rearranged btw calls!)
  std::list<DADFRegionInfo *>::iterator driter;
  driter= _regionList.begin();     
  for ( i = 0; i < region; ++i ) { ++driter; }

  (*driter)->getBounds( lower, upper );
  (*driter)->getDataLocation(data, strides );

  return 0;
}

/** Print the contents of the descriptor (for debugging) */
void DADFDescriptor::printDescriptor() {
  int i;
  static const std::string typeLabels[12] 
    = { "stvInt", "stvFloat", "stvCplx", "stvDouble", "stvDcplx",
	"stvLong", "stvShort", "stvStr", "stvUshort", "stvUint",
	"stvUlong", "stvByte" };

  cerr << "Distributed array descriptor `" << _name << "' rank " <<
    _rank << " type `" << typeLabels[ _type ] << "'" <<
    ((_frozen)? " (" : " (not" ) << " frozen)" << endl; 

  if ( _templ != 0 ) {
    cerr << "      Associated with template: `" << _templ->getName()
	 << "'" << endl;
  }

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
  
  cerr << "      Process coordinates:   " ;
  for ( i=0; i < _rank-1 ; ++i ) {
    cerr << _procCoords[i] << ", ";
  }
  cerr << _procCoords[i] << endl;
  
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
