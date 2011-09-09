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

using namespace std;
using namespace boost;

namespace psi {

AIOHandler::AIOHandler(boost::shared_ptr<PSIO> psio)
    : psio_(psio)
{
    locked_ = new boost::mutex();
}
AIOHandler::~AIOHandler() 
{
    delete locked_;
}
boost::shared_ptr<boost::thread> AIOHandler::get_thread()
{
    return thread_;
}
void AIOHandler::synchronize()
{
    unique_lock<mutex> lock(*locked_);
    lock.unlock();
    thread_->join();
}
void AIOHandler::read(unsigned int unit, const char *key, char *buffer, ULI size, psio_address start, psio_address *end) 
{
  unique_lock<mutex> lock(*locked_);

  job_.push(1);
  unit_.push(unit);
  key_.push(key);
  buffer_.push(buffer);
  size_.push(size);
  start_.push(start);
  end_.push(end);

  if (job_.size() > 1) return;

  //thread start
  thread_ = boost::shared_ptr<boost::thread>(new boost::thread(boost::bind(&AIOHandler::call_aio,this)));
}
void AIOHandler::write(unsigned int unit, const char *key, char *buffer, ULI size, psio_address start, psio_address *end) 
{
  unique_lock<mutex> lock(*locked_);

  job_.push(2);
  unit_.push(unit);
  key_.push(key);
  buffer_.push(buffer);
  size_.push(size);
  start_.push(start);
  end_.push(end);

  if (job_.size() > 1) return;

  //thread start
  thread_ = boost::shared_ptr<boost::thread>(new boost::thread(boost::bind(&AIOHandler::call_aio,this)));
}
void AIOHandler::read_entry(unsigned int unit, const char *key, char *buffer, ULI size) 
{
  unique_lock<mutex> lock(*locked_);

  job_.push(3);
  unit_.push(unit);
  key_.push(key);
  buffer_.push(buffer);
  size_.push(size);

  if (job_.size() > 1) return;

  //thread start
  thread_ = boost::shared_ptr<boost::thread>(new boost::thread(boost::bind(&AIOHandler::call_aio,this)));
}
void AIOHandler::write_entry(unsigned int unit, const char *key, char *buffer, ULI size) 
{
  unique_lock<mutex> lock(*locked_);

  job_.push(4);
  unit_.push(unit);
  key_.push(key);
  buffer_.push(buffer);
  size_.push(size);

  if (job_.size() > 1) return;

  //thread start
  thread_ = boost::shared_ptr<boost::thread>(new boost::thread(boost::bind(&AIOHandler::call_aio,this)));
}
void AIOHandler::read_discont(unsigned int unit, const char *key, 
  double **matrix, ULI row_length, ULI col_length, ULI col_skip, 
  psio_address start)
{
  unique_lock<mutex> lock(*locked_);

  job_.push(5);
  unit_.push(unit);
  key_.push(key);
  matrix_.push(matrix);
  row_length_.push(row_length);
  col_length_.push(col_length);
  col_skip_.push(col_skip);
  start_.push(start);

  if (job_.size() > 1) return;

  //thread start
  thread_ = boost::shared_ptr<boost::thread>(new boost::thread(boost::bind(&AIOHandler::call_aio,this)));
}
void AIOHandler::write_discont(unsigned int unit, const char *key, 
  double **matrix, ULI row_length, ULI col_length, ULI col_skip, 
  psio_address start)
{
  unique_lock<mutex> lock(*locked_);

  job_.push(6);
  unit_.push(unit);
  key_.push(key);
  matrix_.push(matrix);
  row_length_.push(row_length);
  col_length_.push(col_length);
  col_skip_.push(col_skip);
  start_.push(start);

  if (job_.size() > 1) return;

  //thread start
  thread_ = boost::shared_ptr<boost::thread>(new boost::thread(boost::bind(&AIOHandler::call_aio,this)));
}
void AIOHandler::zero_disk(unsigned int unit, const char *key,
    ULI rows, ULI cols)
{
  unique_lock<mutex> lock(*locked_);

  job_.push(7);
  unit_.push(unit);
  key_.push(key);
  row_length_.push(rows);
  col_length_.push(cols);

  if (job_.size() > 1) return;

  //thread start
  thread_ = boost::shared_ptr<boost::thread>(new boost::thread(boost::bind(&AIOHandler::call_aio,this)));
}
void AIOHandler::call_aio()
{
  unique_lock<mutex> lock(*locked_);

  while (job_.size() > 0) {
    int jobtype = job_.front();
    lock.unlock();

    if (jobtype == 1) { 

      lock.lock();

      unsigned int unit = unit_.front();
      const char* key = key_.front();
      char* buffer = buffer_.front();
      ULI size = size_.front();
      psio_address start = start_.front();
      psio_address* end = end_.front();

      job_.pop();
      unit_.pop();
      key_.pop();
      buffer_.pop();
      size_.pop();
      start_.pop();
      end_.pop();

      lock.unlock();

      psio_->read(unit,key,buffer,size,start,end);
    }
    else if (jobtype == 2) {

      lock.lock();

      unsigned int unit = unit_.front();
      const char* key = key_.front();
      char* buffer = buffer_.front();
      ULI size = size_.front();
      psio_address start = start_.front();
      psio_address* end = end_.front();

      job_.pop();
      unit_.pop();
      key_.pop();
      buffer_.pop();
      size_.pop();
      start_.pop();
      end_.pop();

      lock.unlock();

      psio_->write(unit,key,buffer,size,start,end);
    }
    else if (jobtype == 3) {

      lock.lock();

      unsigned int unit = unit_.front();
      const char* key = key_.front();
      char* buffer = buffer_.front();
      ULI size = size_.front();

      job_.pop();
      unit_.pop();
      key_.pop();
      buffer_.pop();
      size_.pop();

      lock.unlock();

      psio_->read_entry(unit,key,buffer,size);
    }
    else if (jobtype == 4) {

      lock.lock();

      unsigned int unit = unit_.front();
      const char* key = key_.front();
      char* buffer = buffer_.front();
      ULI size = size_.front();

      job_.pop();
      unit_.pop();
      key_.pop();
      buffer_.pop(); 
      size_.pop();

      lock.unlock();

      psio_->write_entry(unit,key,buffer,size);
    }
    else if (jobtype == 5) {

      lock.lock();

      unsigned int unit = unit_.front();
      const char* key = key_.front();
      double** matrix = matrix_.front();
      ULI row_length = row_length_.front();
      ULI col_length = col_length_.front();
      ULI col_skip = col_skip_.front();
      psio_address start = start_.front();

      job_.pop();
      unit_.pop();
      key_.pop();
      matrix_.pop();
      row_length_.pop();
      col_length_.pop();
      col_skip_.pop();
      start_.pop();

      lock.unlock();

      for (int i=0; i<row_length; i++) {
        psio_->read(unit,key,(char *) &(matrix[i][0]),
          sizeof(double)*col_length,start,&start);
        start = psio_get_address(start,sizeof(double)*col_skip);
      }
    }
    else if (jobtype == 6) {

      lock.lock();

      unsigned int unit = unit_.front();
      const char* key = key_.front();
      double** matrix = matrix_.front();
      ULI row_length = row_length_.front();
      ULI col_length = col_length_.front();
      ULI col_skip = col_skip_.front();
      psio_address start = start_.front();

      job_.pop();
      unit_.pop();
      key_.pop();
      matrix_.pop();
      row_length_.pop();
      col_length_.pop();
      col_skip_.pop();
      start_.pop();

      lock.unlock();

      for (int i=0; i<row_length; i++) {
        psio_->write(unit,key,(char *) &(matrix[i][0]),
          sizeof(double)*col_length,start,&start);
        start = psio_get_address(start,sizeof(double)*col_skip);
      }
    }
    else if (jobtype == 7) {

      lock.lock();

      unsigned int unit = unit_.front();
      const char* key = key_.front();
      ULI row_length = row_length_.front();
      ULI col_length = col_length_.front();

      job_.pop();
      unit_.pop();
      key_.pop();
      row_length_.pop();
      col_length_.pop();

      lock.unlock();

      double* buf = new double[col_length];
      memset(static_cast<void*>(buf),'\0',col_length*sizeof(double));

      psio_address next_psio = PSIO_ZERO;
      for (int i=0; i<row_length; i++) {
        psio_->write(unit,key,(char *) (buf),sizeof(double)*col_length,
          next_psio,&next_psio);
      }

      delete[] buf;
    }
    else {
      throw PsiException("Error in AIO: Unknown job type", __FILE__,__LINE__);
    }

    lock.lock();
  }
}

} //Namespace psi

