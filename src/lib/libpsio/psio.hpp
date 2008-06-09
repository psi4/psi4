#ifndef _psi_src_lib_libpsio_psio_hpp_
#define _psi_src_lib_libpsio_psio_hpp_

#include <string>
#include <map>

#include <libpsio/config.h>

namespace psi {
  
  /**
   PSIO is an instance of libpsio library. Multiple instances of PSIO are supported.

   Each instance can be configured using filecfg_kwd().
   The following example best demonstrates how to configure a PSIO instance Lib:
   Lib->filecfg_kwd("DEFAULT","NAME",-1,"newwfn")      // all modules will set filename prefix to newwfn for all units
   Lib->filecfg_kwd("DEFAULT","NVOLUME",34,"2")        // all modules will stripe unit 34 over 2 volumes
   Lib->filecfg_kwd("CINTS","VOLUME1",-1,"/scratch1/") // module CINTS will access volume 1 of all units under /scratch
   etc.

   */
  class PSIO {
    public:
      PSIO();
      ~PSIO();

      /// return 1 if activated
      int state() {
        return state_;
      }
      /**
       set keyword kwd describing some aspect of configuration of PSIO file unit
       to value kwdval. kwdgrp specifies the keyword group (useful values are: "DEFAULT", "PSI", and the name of
       the current executable). If unit is set to -1, this keyword will set the default for all units (this keyword
       can be further overridden for some units). To specify a keyword that works for a specific unit, set unit to the
       appropriate number between 0 to PSIO_MAXUNIT.
       
       PSIO understands the following keywords: "name" (specifies the prefix for the filename,
       i.e. if name is set to "psi" then unit 35 will be named "psi.35"), "nvolume" (number of files over which
       to stripe this unit, cannot be greater than PSIO_MAXVOL), "volumeX", where X is a positive integer less than or equal to
       the value of "nvolume".
       */
      void filecfg_kwd(const char* kwdgrp, const char* kwd, int unit,
                       const char* kwdval);
      /// returns the keyword value. If not defined, returns empty string.
      const std::string& filecfg_kwd(const char* kwdgrp, const char* kwd,
                                     int unit);

      /// open unit. status can be PSIO_OPEN_OLD (if existing file is to be opened) or PSIO_OPEN_NEW if new file should be open
      void open(unsigned int unit, int status);
      /// close unit. if keep == 0, will remove the file, else keep it
      void close(unsigned int unit, int keep);
      /// sync up the object to the file on disk by closing and opening the file, if necessary
      void rehash(unsigned int unit);
      /// return 1 if unit is open
      int open_check(unsigned int unit);
      /** Reads data from within a TOC entry from a PSI file.
       **
       **  \param unit   = The PSI unit number used to identify the file to all
       **                  read and write functions.
       **  \param key    = The TOC keyword identifying the desired entry.
       **  \param buffer = The buffer to store the data as it is read.
       **  \param size   = The number of bytes to read.
       **  \param start  = The entry-relative starting page/offset of the desired data.
       **  \param end    = A pointer to the entry-relative page/offset for the next
       **                  byte after the end of the read request.
       */
      void read(unsigned int unit, char *key, char *buffer, ULI size,
                psio_address start, psio_address *end);
      /** Writes data to a TOC entry in a PSI file.
       **
       **  \param unit    = The PSI unit number used to identify the file to all read
       **                   and write functions.
       **  \param key     = The TOC keyword identifying the desired entry.
       **  \param buffer  = The buffer from which the data is written.
       **  \param size    = The number of bytes to write.
       **  \param start   = The entry-relative starting page/offset to write the data.
       **  \param end     = A pointer to the entry-relative page/offset for the next
       **                   byte after the end of the write request.
       */
      void write(unsigned int unit, char *key, char *buffer, ULI size,
                 psio_address start, psio_address *end);

      void read_entry(unsigned int unit, char *key, char *buffer, ULI size);
      void write_entry(unsigned int unit, char *key, char *buffer, ULI size);

      /** Central function for all reads and writes on a PSIO unit.
       **
       ** \params unit    = The PSI unit number.
       ** \params buffer  = The buffer containing the bytes for the read/write event.
       ** \params address = the PSIO global address for the start of the read/write.
       ** \params size    = The number of bytes to read/write.
       ** \params         = Indicates if the call is to read (0) or write (0) the input data.
       **
       ** \ingroup PSIO
       */
      void rw(unsigned int unit, char *buffer, psio_address address, ULI size,
              int wrt);
      /// Delete all TOC entries after the given key. If a blank key is given, the entire TOC will be wiped.
      void tocclean(unsigned int unit, char *key);
      /// Print the table of contents for the given unit
      void tocprint(unsigned int unit, FILE *output);
      /// Scans the TOC for a particular keyword and returns either a pointer to the entry or NULL to the caller.
      psio_tocentry* tocscan(unsigned int unit, char *key);
      ///  Write the table of contents for file number 'unit'. NB: This function should NOT call psio_error because the latter calls it!
      void tocwrite(unsigned int unit);

      /// Upon catastrophic failure, the library will exit() with this code. The default is 1, but can be overridden.
      static int _error_exit_code_;

    private:
      /// vector of units
      psio_ud *psio_unit;

      typedef std::map<std::string,std::string> KWDMap;
      /// library configuration is described by a set of keywords
      KWDMap files_keywords_;

#ifdef PSIO_STATS
      ULI *psio_readlen;
      ULI *psio_writlen;
#endif
      
      /// Library state variable
      int state_;

      /// grab the filename of unit and strdup into name.
      void get_filename(unsigned int unit, char **name);
      /// return the number of volumes over which unit will be striped
      unsigned int get_numvols(unsigned int unit);
      /// grab the path to volume of unit and strdup into path.
      void get_volpath(unsigned int unit, unsigned int volume, char **path);
      /// return the last TOC entry
      psio_tocentry* toclast(unsigned int unit);
      /// Compute the length of the TOC for a given unit using the in-core TOC list.
      unsigned int toclen(unsigned int unit);
      /** Read the length of the TOC for a given unit directly from the file.
       **
       ** \param unit = PSI unit number from which to read the toclen.
       **
       ** NB: Note that we do not exit if the read request of the toclen from
       ** the file fails. This is because the request may be to an new file
       ** for which the toclen has not yet been written.  (We allow the user
       ** to open files with status PSIO_OPEN_OLD even if they don't exist,
       ** because sometimes you can't know this in advance.)
       */
      ULI rd_toclen(unsigned int unit);
      /** Write the length of the TOC for a given unit directly to the file.
       **
       ** \param unit = PSI unit number to which to write the toclen.
       **
       ** \ingroup PSIO
       */
      void wt_toclen(unsigned int unit, ULI toclen);
      /// Read the table of contents for file number 'unit'.
      void tocread(unsigned int unit);

  };
  
  extern int psiopp_ipv1_config(PSIO *psio_obj);
  extern PSIO* _default_psio_lib_;
}

#endif /* header guard */
