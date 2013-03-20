#include <string>

#include "package.h"

#if defined(OPTKING_PACKAGE_PSI)
  #include <libmints/writer_file_prefix.h>
#endif

#if defined(OPTKING_PACKAGE_QCHEM)
 #include "qcsys.h"
 extern void getPrefix(char*& pref);
#endif

const char* getIntcoFileName()
{
   static std::string strintco("");
   if (strintco.empty() ) {
#if defined(OPTKING_PACKAGE_QCHEM)
      char* pref=NULL;
      getPrefix(pref);
      strintco = std::string(pref) + "intco.dat";
#elif defined(OPTKING_PACKAGE_PSI)
      strintco = psi::get_writer_file_prefix() + ".intco";
#endif
   }
   return strintco.c_str();
}

const char* getOptdataFileName()
{
   static std::string stroptdata("");
   if (stroptdata.empty() ) {
#if defined(OPTKING_PACKAGE_QCHEM)
      char* pref=NULL;
      getPrefix(pref);
      stroptdata = std::string(pref) + "opt_data.1";
#elif defined(OPTKING_PACKAGE_PSI)
      // In PSI, the opt data file is file 1, with the name set by input and io_start
#endif
   }
   return stroptdata.c_str();
}

