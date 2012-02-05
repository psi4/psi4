#ifndef _psi_src_lib_libpsio_psio_hpp_
#define _psi_src_lib_libpsio_psio_hpp_

//#include <boost/thread/thread.hpp>
#include <string>
#include <map>
#include <set>
#include <queue>

#include <libpsio/config.h>

namespace boost {
template <class T>
class shared_ptr;
class thread;
}

namespace psi {

extern FILE *outfile;

class PSIO;
class PSIOManager;
extern boost::shared_ptr<PSIO> _default_psio_lib_;
extern boost::shared_ptr<PSIOManager> _default_psio_manager_;

/**
    PSIOManager is a class designed to be used as a static object to track all
    PSIO operations in a given PSI4 computation

    This will allow PSICLEAN to be trivially executed.
    Now supports a .psirc and interactive file placement
   */
class PSIOManager {
private:
    /// Default path for unspec'd file numbers (defaults to /tmp/)
    std::string default_path_;
    /// Specific paths for arbitrary file numbers
    std::map<int, std::string> specific_paths_;
    /// Default retained files
    std::set<int> specific_retains_;

    /// Map of files, bool denotes open or closed
    std::map<std::string, bool> files_;
    /// Set of files to retain after psiclean
    std::set<std::string> retained_files_;
public:
    /// Default constructor (does nothing)
    PSIOManager();
    /// Default destructor (does nothing)
    ~PSIOManager();

    /**
            * Mirror the current delete-able files to "psi.clean"
            **/
    void mirror_to_disk();
    /**
            * Build from "psi.clean"
            **/
    void build_from_disk();

    /**
            * Set the default path for files to be stored
            * \param path full path to scratch
            */
    void set_default_path(const std::string& path);
    /**
            * Set the path for specific file numbers
            * \param fileno PSI4 file number
            * \param path full path to file-specific scratch
            */
    void set_specific_path(int fileno, const std::string& path);
    /**
            * Set the the specific file number to be retained
            * \param fileno PSI4 file number
            * \param retain keep or not? (Allows override)
            */
    void set_specific_retention(int fileno, bool retain);
    /**
            * Get the path for a specific file number
            * \param fileno PSI4 file number
            * \return the appropriate full path
            */
    std::string get_file_path(int fileno);

    /**
      * Returns the default path.
      * \return the default path.
      */
    std::string get_default_path() { return default_path_; }

    /**
     * Write a string to a temporary file.  The scratch file is opened and closed by this function.
     * @param full_path The fill path to the scratch file
     * @param text The text to be placed in the file
     */
    void write_scratch_file(const std::string &full_path, const std::string &text);

    /**
            * Record the opening of a file
            * \param full_path filename
            */
    void open_file(const std::string & full_path, int fileno);
    /**
            * Record the opening of a file
            * \param fileno PSI4 file number
            * \param full_path filename
            * \param keep TRUE : the file is closed and retained by PSIO
                          FALSE: the file is closed and deleted by PSIO
            */
    void close_file(const std::string & full_path, int fileno, bool keep);
    /**
            * Move a file from one location to another, retaining status
            * Useful for changing namespaces
            * \param old_full_path old filename
            * \param new_full_path new filename
            */
    void move_file(const std::string & old_full_path, const std::string & new_full_path);
    /**
            * Mark a file to be retained after a psiclean operation, ie for use in
            * a later computation
            * \param full_path filename
            * \param retain keep or not? (Allows override)
            */
    void mark_file_for_retention(const std::string & full_path, bool retain);
    /**
            * Print the current status of PSI4 files
            * \param out, file to print fo
            */
    void print(FILE* out = outfile);
    /**
            * Print the current status of PSI4 files
            * \param out, file to print fo
            */
    void print_out() { print(outfile); }
    /**
            * Execute the psiclean protocol, deleting all recorded files
            * except for those currently marked for retention.
            *
            * Those files marked for retention are not deleted, and their
            * traces in the files_ set and retained_files set remain.
            * Deleted files are removed from the files_ set.
            *
            * This is useful for intermediate calls to psiclean
            */
    void psiclean();
    /**
            * Clean from disk-mirrored image after crash
            * NOT to be called during regular computation
            **/
    void crashclean();
    /// The one and (should be) only instance of PSIOManager for a PSI4 instance
    static boost::shared_ptr<PSIOManager> shared_object();
};

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
    void read(unsigned int unit, const char *key, char *buffer, ULI size,
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
    void write(unsigned int unit, const char *key, char *buffer, ULI size,
               psio_address start, psio_address *end);

    void read_entry(unsigned int unit, const char *key, char *buffer, ULI size);
    void write_entry(unsigned int unit, const char *key, char *buffer, ULI size);

    /** Zeros out a double precision array in a PSI file.
       ** Typically used before striping out a transposed array
       **  Total fill size is rows*cols*sizeof(double)
       **  Buffer memory of cols*sizeof(double) is used
       **
       **  \param unit    = The PSI unit number used to identify the file
       **  \param key     = The TOC keyword identifying the desired entry.
       **  \param rows    = The number of rows in the full array
       **  \param cols    = The number of columnss in the full array
       **
       */
    void zero_disk(unsigned int unit, const char *key, ULI rows, ULI cols);

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
    void tocclean(unsigned int unit, const char *key);
    /// Print the table of contents for the given unit
    void tocprint(unsigned int unit);
    /// Scans the TOC for a particular keyword and returns either a pointer to the entry or NULL to the caller.
    psio_tocentry* tocscan(unsigned int unit, const char *key);
    ///  Write the table of contents for file number 'unit'. NB: This function should NOT call psio_error because the latter calls it!
    void tocwrite(unsigned int unit);

    /// Upon catastrophic failure, the library will exit() with this code. The default is 1, but can be overridden.
    static int _error_exit_code_;

    /// Set the current namespace (for PREFIX.NAMESPACE.UNIT file numbering)
    static void set_default_namespace(const std::string &_ns) { default_namespace_ = _ns; }

    /// Get the default namespace (for PREFIX.NAMESPACE.UNIT file numbering)
    static std::string get_default_namespace() { return default_namespace_; }

    /// Change file FILENO from NS1 to NS2
    static void change_file_namespace(unsigned int fileno, const std::string & ns1, const std::string & ns2);

    /// Return the global shared object
    static boost::shared_ptr<PSIO> shared_object();

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

private:
    /// vector of units
    psio_ud *psio_unit;

    /// Process ID
    std::string pid_;

    /// Current default namespace (for PREFIX.NAMESPACE.UNIT numbering)
    static std::string default_namespace_;

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
    void get_filename(unsigned int unit, char **name, bool remove_namespace = false);
    /// return the number of volumes over which unit will be striped
    unsigned int get_numvols(unsigned int unit);
    /// grab the path to volume of unit and strdup into path.
    void get_volpath(unsigned int unit, unsigned int volume, char **path);
    /// return the last TOC entry
    psio_tocentry* toclast(unsigned int unit);
    /// Compute the length of the TOC for a given unit using the in-core TOC list.
    unsigned int toclen(unsigned int unit);
    /** Write the length of the TOC for a given unit directly to the file.
       **
       ** \param unit = PSI unit number to which to write the toclen.
       **
       ** \ingroup PSIO
       */
    void wt_toclen(unsigned int unit, ULI toclen);
    /// Read the table of contents for file number 'unit'.
    void tocread(unsigned int unit);

    friend class AIO_Handler;
};

}

#endif /* header guard */
