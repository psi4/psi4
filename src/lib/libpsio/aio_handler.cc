/*!
 \file
 \ingroup PSIO
 */

#include <cstdio>
#include <unistd.h>
#include <boost/shared_ptr.hpp>
#include <boost/thread/thread.hpp>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>

using namespace boost;

namespace psi {

AIO_Handler::AIO_Handler(shared_ptr<PSIO> psio)
    : psio_(psio)
{
}
AIO_Handler::~AIO_Handler() {
}
shared_ptr<boost::thread> AIO_Handler::get_thread()
{
    return thread_;
}
void AIO_Handler::synchronize()
{
    thread_->join();
}
void AIO_Handler::read(unsigned int unit, const char *key, char *buffer, ULI size, psio_address start, psio_address *end) 
{
  unit_ = unit;
  key_ = key;
  buffer_ = buffer;
  size_ = size;
  start_ = start;
  end_ = end;

  //thread start
  thread_ = shared_ptr<boost::thread>(new boost::thread(boost::bind(&AIO_Handler::call_read,this)));
}
void AIO_Handler::write(unsigned int unit, const char *key, char *buffer, ULI size, psio_address start, psio_address *end) 
{
  unit_ = unit;
  key_ = key;
  buffer_ = buffer;
  size_ = size;
  start_ = start;
  end_ = end;

  //thread start
  thread_ = shared_ptr<boost::thread>(new boost::thread(boost::bind(&AIO_Handler::call_write,this)));
}
void AIO_Handler::read_entry(unsigned int unit, const char *key, char *buffer, ULI size) 
{
  unit_ = unit;
  key_ = key;
  buffer_ = buffer;
  size_ = size;

  //thread start
  thread_ = shared_ptr<boost::thread>(new boost::thread(boost::bind(&AIO_Handler::call_read_entry,this)));
}
void AIO_Handler::write_entry(unsigned int unit, const char *key, char *buffer, ULI size) 
{
  unit_ = unit;
  key_ = key;
  buffer_ = buffer;
  size_ = size;

  //thread start
  thread_ = shared_ptr<boost::thread>(new boost::thread(boost::bind(&AIO_Handler::call_write_entry,this)));
}
void AIO_Handler::call_read()
{
    psio_->read(unit_,key_,buffer_,size_,start_,end_);
}
void AIO_Handler::call_write()
{
    psio_->write(unit_,key_,buffer_,size_,start_,end_);
}
void AIO_Handler::call_read_entry()
{
    psio_->read_entry(unit_,key_,buffer_,size_);
}
void AIO_Handler::call_write_entry()
{
    psio_->write_entry(unit_,key_,buffer_,size_);
}

} //Namespace psi

