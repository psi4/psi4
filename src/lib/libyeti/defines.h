#ifndef _DEFINES_H
#define	_DEFINES_H

// Set either MPQC_INTERFACE or PSI4_INTERFACE to 1
#define MPQC_INTERFACE 1

#if MPQC_INTERFACE
   #define OPEN_NAMESPACE namespace tiled_tensor{
   #define CLOSE_NAMESPACE }
#elif PSI4_INTERFACE
   #define OUT0 ExEnv::out0()
   #define OUT0 ExEnv::out0()
#else
   #error "You must select an interface in defines.h"
#endif

#define harikari(msg) std::cerr << msg << endl; abort()


#endif	/* _DEFINES_H */

