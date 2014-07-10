#ifndef DistArrayDescrFactoryPort_h_seen
#define DistArrayDescrFactoryPort_h_seen

#include <string>
#include "DistArrayDescriptor.h"
#include "DistArray.h"

namespace classic {
  namespace gov {
    namespace cca {
      
      /** An interface for creating distributed array descriptors.
	  
      $Id: DistArrayDescrFactoryPort.h,v 1.1 2003-08-01 00:41:53 manoj Exp $
      
      This interface allows users to create, clone, and destroy
      descriptors for distributed array objects.
      */
      
      class DistArrayDescrFactoryPort: public virtual ::classic::gov::cca::Port {
	
      public:
	
	/** Return an uninitialized descriptor object.
	    
	@param name (in) Name associated with this descriptor.  This
	is a "forced convenience" for the user to allow intelligible
	error messages, etc.
	
	@returns a reference to an uninitialized descriptor object.
	*/
	virtual DistArrayDescriptor * createDescriptor(std::string name) = 0;
	
	/** Return a descriptor object initialized with the contents of
	    another, but not frozen against modification.
	    
	    @param original (in) Descriptor from which to initialize new
	    descriptor. 
	    
	    @param cloneName (in) Name associated with this descriptor.
	    This is a "forced convenience" for the user to allow
	    intelligible error messages, etc.
	    
	    @returns a reference to an initialized but not frozen
	    descriptor object.
	*/
	virtual DistArrayDescriptor * cloneDescriptor(DistArrayDescriptor * original, std::string cloneName) = 0;
	
	/** Destroy an existing descriptor.
	    
	@param victim (in) Descriptor to be destroyed.
	
	@retval  0 Success
	@retval -1 Descriptor not created by this factory
	@retval -2 Descriptor has already been destroyed
	@retval -3 Concrete class of descriptor is not is not the type
	produced by this factory.
	
	@notes Should throw exceptions instead of returning an int.
	*/
	virtual int destroyDescriptor(DistArrayDescriptor * & victim) = 0;

	/** Return an uninitialized ga-dadf distributed array object. */
	virtual DistArray * createArray(std::string name) = 0;
 
	/** Return an distributed array object initialized with the contents of
	    another, but not frozen against modification. */
	virtual DistArray * cloneArray(DistArray* original,
				       std::string cloneName) = 0;
 
	/** Destroy an existing distributed array. */
	virtual int destroyArray(DistArray* & victim) = 0;
	
      }; // DistArrayDescrFactoryPort
      
    } // namespace cca
  } // namespace gov
} // namespace classic
  
#endif // DistArrayDescrFactoryPort_h_seen
