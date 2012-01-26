#ifndef GLOBALS_H
#define GLOBALS_H

#ifdef EXTERN
    #undef EXTERN
    #define EXTERN extern
#else
    #define EXTERN
#endif


namespace psi{ namespace plugin_cepa{

  //EXTERN boost::shared_ptr<PSIO> psio;
  //EXTERN struct calcinfo CalcInfo;
  //EXTERN struct params Parameters;
  //EXTERN int *ioff;

}}

#endif
