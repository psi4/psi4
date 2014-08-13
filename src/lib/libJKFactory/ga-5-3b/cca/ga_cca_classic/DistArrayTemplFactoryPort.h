#ifndef DistArrayTemplFactoryPort_h_seen
#define DistArrayTemplFactoryPort_h_seen

#include <string>
#include "DistArrayTemplate.h"

namespace classic {
  namespace gov {
    namespace cca {
      
      /** An interface for creating distributed array template
	  descriptors.
	  
	  $Id: DistArrayTemplFactoryPort.h,v 1.1 2003-08-01 00:41:54 manoj Exp $
	  
	  This interface allows users to create, clone, and destroy
	  descriptors for array distribution templates.
      */
      
      class DistArrayTemplFactoryPort: public virtual ::classic::gov::cca::Port {
	
      public:
	
	/** Return an uninitialized template object.
	    
	@param name (in) Name associated with this template.  This
	is a "forced convenience" for the user to allow intelligible
	error messages, etc.
	
	@returns a reference to an uninitialized template object.
	*/
	virtual DistArrayTemplate* createTemplate(std::string name) = 0;
	
	/** Return a template object initialized with the contents of
	    another, but not frozen against modification.
	    
	    @param original (in) Template from which to initialize new
	    template. 
	    
	    @param cloneName (in) Name associated with this template.
	    This is a "forced convenience" for the user to allow
	    intelligible error messages, etc.
	    
	    @returns a reference to an initialized but not frozen
	    template object.
	*/
	virtual DistArrayTemplate* cloneTemplate(DistArrayTemplate * original, std::string cloneName) = 0;
	
	/** Destroy an existing template.
	    
	@param victim (in) Template to be destroyed.
	
	@retval  0 Success
	@retval -1 Template not created by this factory
	@retval -2 Template has already been destroyed
	@retval -3 Concrete class of template is not is not the type
	produced by this factory.
	
	@notes Should throw exceptions instead of returning an int.
	*/
	virtual int destroyTemplate(DistArrayTemplate * & victim) = 0;
	
      }; // DistArrayTemplFactoryPort
      
    } // namespace cca
  } // namespace gov
} // namespace classic

#endif // DistArrayTemplFactoryPort_h_seen
