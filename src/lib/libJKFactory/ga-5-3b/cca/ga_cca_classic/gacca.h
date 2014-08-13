#ifndef _GAPP_H
#define _GAPP_H

#if defined(__cplusplus) || defined(c_plusplus)

#include <stdio.h>
#include <string>
#include <map>
#include "ga.h"
#include "macdecls.h"
#include <mpi.h>


// CCAFFEINE Includes
#include <cca.h>
#include <stdPorts.h>
#include <EG.h>
#include "jc++/jc++.h"
#include "jc++/lang/jc++lang.h"
#include "jc++/util/jc++util.h"
#include "parameters/parametersStar.h"
#include "util/IO.h"

// DADF Includes
#include "DistArrayTemplFactoryPort.h"
#include "DistArrayDescrFactoryPort.h"


#define _GA_USENAMESPACE_ 1

#if _GA_USENAMESPACE_
#define _GA_STATIC_ 
#define _GA_EXTERN_ extern
#else
#define _GA_STATIC_ static
#define _GA_EXTERN_ 
#endif

#if _GA_USENAMESPACE_
namespace GA {
#else
class GA {
 public:
#endif
  class GAClassicPort;
  class GAServices;

#include "GlobalArray.h"
#include "GAClassicPort.h"
#include "GAServices.h"

  //GAServices SERVICES;

#if ! _GA_USENAMESPACE_
 private:
  GA() { }
#endif
};

#endif // _GAPP_H
#endif
