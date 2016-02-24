/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

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
#include "aiohandler.h"
#include <exception.h>

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
    // This join may be problematic in a multithreaded env.: while we wait for the
    // thread to finish, other threads may add work to the queue. We'd probably need 
    // a way to identify write jobs and check if they completed from external threads.
    boost::unique_lock<boost::mutex> lock(*locked_);
    lock.unlock();
    thread_->join();
}
void AIOHandler::read(unsigned int unit, const char *key, char *buffer, ULI size, psio_address start, psio_address *end)
{
  boost::unique_lock<boost::mutex> lock(*locked_);

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
  boost::unique_lock<boost::mutex> lock(*locked_);

  job_.push(2);
  unit_.push(unit);
  key_.push(key);
  buffer_.push(buffer);
  size_.push(size);
  start_.push(start);
  end_.push(end);

  //printf("Adding a write to the queue\n");

  // Ensures we do not already have a thread running jobs.
  //if (thread_ == NULL) {
  //    fprintf(stderr,"Thread is null\n");
  //} else if (not thread_->timed_join(boost::chrono::milliseconds(0)) ) {
  //    fprintf(stderr,"Thread is not joinable now, thus busy\n");
  //    return;
  //} else {
  //    if(job_.size() > 1) {
  //      fprintf(stderr,"%d jobs in the queue! Should NOT create thread!!\n", job_.size());
  //    }
  //    fprintf(stderr, "Thread is joined, need new thread\n");
  //}
  if (job_.size() > 1) return;

  //fprintf(stderr,"Starting a thread\n");
  //thread start
  thread_ = boost::shared_ptr<boost::thread>(new boost::thread(boost::bind(&AIOHandler::call_aio,this)));
}
void AIOHandler::read_entry(unsigned int unit, const char *key, char *buffer, ULI size)
{
  boost::unique_lock<boost::mutex> lock(*locked_);

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
  boost::unique_lock<boost::mutex> lock(*locked_);

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
  boost::unique_lock<boost::mutex> lock(*locked_);

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
  boost::unique_lock<boost::mutex> lock(*locked_);

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
  boost::unique_lock<boost::mutex> lock(*locked_);

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
  boost::unique_lock<boost::mutex> lock(*locked_);

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

      unit_.pop();
      key_.pop();
      buffer_.pop();
      size_.pop();
      start_.pop();
      end_.pop();

      //printf("Actually doing the psio_write from thread\n");
      //printf("Unit is %d, key is %s, size is %d\n", unit, key, size);
      lock.unlock();

      //printf("Printing number %d now\n", *((int *) buffer));
      //std::stringstream ss;
      //ss << thread_->get_id();
      //fprintf(stderr,"Printing from thread %s\n", ss.str().c_str());
      //fprintf(stderr,"Starting address is %lu, %lu\n",start.page, start.offset);
      psio_->write(unit,key,buffer,size,start,end);
      //fprintf(stderr,"End address is %lu, %lu\n",end->page, end->offset);
    }
    else if (jobtype == 3) {

      lock.lock();

      unsigned int unit = unit_.front();
      const char* key = key_.front();
      char* buffer = buffer_.front();
      ULI size = size_.front();

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
    // Only pop the job once the work is actually done and we are gonna leave the loop.
    // This way, job_.size() == 0 indicates there is no active thread.
    job_.pop(); 
  }
  //printf("End of function call_aio\n");
}

} //Namespace psi
