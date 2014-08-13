#include "gacca.h"
#include "GA_DADFArray.h"
#include "DistArrayTemplate.h"

/** This is our implementation of DistArray.  It is kept
    private within the GAServices.  Only
    DistArray is exposed to the outside. 
 */

/******************************************************************************
 * Constructors and destructors
 *****************************************************************************/

/** Constructor sets distributed array name. */
DADFArray::DADFArray( const std::string name) {
  _name = name;
  _frozen = false;
  _rank = -1;
  _type = stv_Int;
  _templ = 0;
  _isExplicitDist = false;
}

/** Construct a new distributed array as a copy of an old one.
*/
DADFArray::DADFArray(const std::string name,  DADFArray & original)
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

  _frozen = true;

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

  _handle = GA_Duplicate(original._handle, (char *)_name.c_str());
  if(!_handle) GA_Error(" DADFArray::DADFArray(): GA creation failed",0);

}

DADFArray::~DADFArray() {
    delete _templ;

    std::list<DADFRegionInfo *>::iterator driter;
    
    // Clean up _regionList
    for ( driter= _regionList.begin(); driter != _regionList.end();
	  ++driter ) {
      delete *driter;
    }
    _regionList.resize(0);

    // Destroy the Global Array
    GA_Destroy(_handle);
}

/******************************************************************************
 * Define the distributed array
 *****************************************************************************/

/** Set data type. */
int DADFArray::setDataType(const enum DataType type) {
  // Insure we haven't been commit()ed
  if ( _frozen ) { return -11; }

  _type = type;

  return 0;
}

/** Associate this data object with a distribution template. */
int DADFArray::setTemplate(DistArrayTemplate * & templ) {
  // Insure we haven't been commit()ed
  if ( _frozen ) { return -11; }

  // If there's already one there, replace it
  if ( _templ != 0 ) { delete _templ; }

  DADFTemplate * templDADF = dynamic_cast<DADFTemplate *>(templ);
  if ( templDADF == 0 ) {
    cerr << "DADFArray:setTemplate: " <<
      "Template is of the wrong type." << endl;
    return -13;
  }

  // Make a private copy of the template
  DADFTemplate * newtempl = new DADFTemplate( *templDADF );
  _templ = dynamic_cast<DistArrayTemplate *>(newtempl);

  return 0;
}

/** Sets this process's location in the process topology. */
int DADFArray::setMyProcCoords(const int procCoords[] ) {
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
  // for ( i=0; i < _rank; ++i ) { _procCoords[ i ] = procCoords[ i ]; }

  GA_Error("DADFArray::setMyProcCoords(): This method is currently not supported in GA\n", 0);
  
  return 0;
}

/** Align object to template with identity mapping. */
int DADFArray::setIdentityAlignmentMap() {
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

/** Signal that distributed array is completely defined. */
int DADFArray::commit() {
  // Insure we haven't been commit()ed already
  if ( _frozen ) { 
    return -11; 
  } else {
    _frozen = true;
    // return 0;
  }
  
  int err, me = GA_Nodeid(), type;
  int chunk[GA_MAX_DIM], dims[GA_MAX_DIM], strides[GA_MAX_DIM];
  int lower[GA_MAX_DIM], upper[GA_MAX_DIM]; 
  void *data;
  DistArrayTemplate::DistType dist[GA_MAX_DIM];

  // get the name of the array
  _name = _templ->getName();

  err = _templ->getGlobalBounds( &(*lower), &(*upper) );
  if ( err != 0 ) { return err; }
  
  err = _templ->getDistType( dist );
  if ( err != 0 ) { return err; }

  // Identify GA Data Type
  switch(_type) {
  case DistArray::stv_Int:
    type = C_INT;
    break;
  case DistArray::stv_Float:
    type = C_FLOAT;
    break;
  case DistArray::stv_Double:
    type = C_DBL;
    break;
  case DistArray::stv_Dcplx:
    type = C_DCPL;
    break;
  case DistArray::stv_Cplx:
    type = C_SCPL;
    break;
  case DistArray::stv_Long:
    type = C_LONG;
    break;
  default:
    GA_Error("DADFArray::commit(): Invalid Data Type\n", 0);
  }
  
  if(dist[0] == DistArrayTemplate::Block) {
    int first;
    for(int i=0; i<_rank; ++i) {
      _templ->getDistParameters(i, chunk[i], first);
      dims[i] = upper[i] - lower[i];// + 1;
    }
    _handle = NGA_Create(type, _rank, dims, (char *)_name.c_str(), chunk);
    if(!_handle) GA_Error(" GA creation failed",0);
    if(!me) cout << "NGA_Create called\n";
  }
  else if(dist[0] == DistArrayTemplate::GenBlock){
    int *blockSize[GA_MAX_DIM], total=0;
    for(int i=0; i<_rank; ++i) {
      blockSize[i] = new int[_topology[i]];
      _templ->getGenBlock(i, blockSize[i]);
      total += _topology[i];
      dims[i] = upper[i] - lower[i]; // + 1;
      chunk[i] = _topology[i];
    }
    int *ga_map  = new int [total];
    int offset=0;
    for(int i=0; i<_rank; i++) {
      for(int j=0; j<_topology[i]; ++j) {
        if(j!=0) ga_map[offset] = ga_map[offset-1] + blockSize[i][j-1];
        else ga_map[offset] = 0;
        offset++;
      }
    }
    _handle = NGA_Create_irreg(type, _rank, dims, (char *)_name.c_str(),
			       chunk, ga_map);
    if(!_handle) GA_Error(" GA creation failed",0);
    GA_Print_distribution(_handle);
    for(int i=0; i<_rank; i++) { delete  blockSize[i]; blockSize[i] = NULL; }
    delete[] ga_map; ga_map = NULL;
    if(!me) cout << "NGA_Create_irreg called\n";
  }
  else
    cerr << "Error: Invalid Distribution Type\n";
  
  // Create a region info and set everything
  NGA_Distribution(_handle, me, lower, upper);
  NGA_Access(_handle, lower, upper, &data, strides);
  DADFRegionInfo * dri = new DADFRegionInfo( _rank );
  dri->setBounds( lower, upper );
  dri->setDataLocation( data, strides );
  _regionList.push_back( dri );
  
  return 0;
  
  // Perform global consistency checks

  // Should test topology if dist is explicit
  // Test if area of explicit regions == area of distributed array
}

/******************************************************************************
 * Query the distributed array
 *****************************************************************************/

std::string DADFArray::getName() {
  return _name;
}

bool DADFArray::isDefined() {
    if ( _frozen ) {
      return true;
    } else {
      return false;
    }
}

/** Get data type. */
DistArray::DataType DADFArray::getDataType() {
  return _type;
}

/** Return pointer to distribution template associated with this
    data object. 
*/
DistArrayTemplate * DADFArray::getTemplate() {
  return _templ;
}

/** Get this process's location in the process topology. */
int DADFArray::getMyProcCoords(int procCoords[] ) {
  int me=GA_Nodeid();

  // Sanity check our state
  if ( _rank < 1 ) { return -8; }

  // Copy arguments out of our data structures
  // for ( int i=0; i < _rank; ++i ) { procCoords[ i ] = _procCoords[ i ]; }
  NGA_Proc_topology(_handle, me, procCoords);

  return 0;
}

// Kludges to check GenBlock and Explicit distributed arrays.
int DADFArray::getNumLocalRegions() {
  //  if ( ! _frozen ) { return -100; }

  if ( _isExplicitDist ) {
    return _regionList.size();
  } else {
    return 1; 
  }
}
// Kludges to check GenBlock and Explicit distributed arrays.
int DADFArray::getLocalRegionInfo(int region, int lower[], int upper[],
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

/** Print the contents of the distributed array (for debugging) */
void DADFArray::printArray() {
  GA_Print(_handle);
}


/** Print the distribution info of the distributed array (for debugging) */
void DADFArray::printArrayDistribution() {
  int i;
  static const std::string typeLabels[12] 
    = { "stvInt", "stvFloat", "stvCplx", "stvDouble", "stvDcplx",
	"stvLong", "stvShort", "stvStr", "stvUshort", "stvUint",
	"stvUlong", "stvByte" };

  cerr << "Distributed array `" << _name << "' rank " <<
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
