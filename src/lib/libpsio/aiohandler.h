#ifndef AIOHANDLER_H
#define AIOHANDLER_H

#include <boost/thread/thread.hpp>

namespace psi {

class AIOHandler {
private:
    /// What is the job type?
    std::queue<unsigned int> job_;
    /// Unit number argument
    std::queue<unsigned int> unit_;
    /// Entry Key (80-char) argument
    std::queue<const char*> key_;
    /// Memory buffer argument
    std::queue<char*> buffer_;
    /// Size argument
    std::queue<ULI> size_;
    /// Start address argument
    std::queue<psio_address> start_;
    /// End address pointer argument
    std::queue<psio_address*> end_;
    /// Matrix pointer for discontinuous I/O
    std::queue<double**> matrix_;
    /// Size argument for discontinuous I/O
    std::queue<ULI> row_length_;
    /// Size argument for discontinuous I/O
    std::queue<ULI> col_length_;
    /// Size argument for discontinuous I/O
    std::queue<ULI> col_skip_;
    /// PSIO object this AIO_Handler is built on
    boost::shared_ptr<PSIO> psio_;
    /// Thread this AIO_Handler is currently running on
    boost::shared_ptr<boost::thread> thread_;
    /// Lock variable
    boost::mutex *locked_;
public:
    /// AIO_Handlers are constructed around a synchronous PSIO object
    AIOHandler(boost::shared_ptr<PSIO> psio);
    /// Destructor
    ~AIOHandler();
    /// Thread object this AIO_Handler is currently running on
    boost::shared_ptr<boost::thread> get_thread();
    /// When called, synchronize will not return until all requested data has been read or written
    void synchronize();
    /// Asynchronous read, same as PSIO::read, but nonblocking
    void read(unsigned int unit, const char *key, char *buffer, ULI size,
              psio_address start, psio_address *end);
    /// Asynchronous write, same as PSIO::write, but nonblocking
    void write(unsigned int unit, const char *key, char *buffer, ULI size,
               psio_address start, psio_address *end);
    /// Asynchronous read_entry, same as PSIO::read_entry, but nonblocking
    void read_entry(unsigned int unit, const char *key, char *buffer, ULI size);
    /// Asynchronous read_entry, same as PSIO::write_entry, but nonblocking
    void write_entry(unsigned int unit, const char *key, char *buffer, ULI size);
    /// Asynchronous read for reading discontinuous disk space
    /// into a continuous chunk of memory, i.e.
    ///
    /// [***]        [----***-----]
    /// [***]        [----***-----]
    /// [***]        [----***-----]
    /// [***]  <<--  [----***-----]
    /// [***]        [----***-----]
    /// [***]        [----***-----]
    /// [***]        [----***-----]
    ///
    /// The buffer has dimensions row_length by col_length. The disk space
    /// has dimensions row_length by col_length + col_skip.
    ///
    /// These functions are not necessary for psio, but for aio they are.
    ///
    void read_discont(unsigned int unit, const char *key, double **matrix,
      ULI row_length, ULI col_length, ULI col_skip, psio_address start);
    /// Same as read_discont, but for writing
    void write_discont(unsigned int unit, const char *key, double **matrix,
      ULI row_length, ULI col_length, ULI col_skip, psio_address start);

    /// Zero disk
    /// Fills a double precision disk entry with zeros
    /// Total fill size is rows*cols*sizeof(double)
    /// Buffer memory of cols*sizeof(double) is used
    void zero_disk(unsigned int unit, const char* key, ULI rows, ULI cols);

    /// Generic function bound to thread internally
    void call_aio();
};

}

#endif // AIOHANDLER_H
