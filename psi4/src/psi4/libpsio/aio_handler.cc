/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
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

#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
#include "aiohandler.h"
#include "psi4/libpsi4util/exception.h"

#include <cstdio>
#include <unistd.h>
#include <memory>
#include <algorithm>
#include <functional>

using namespace std;

namespace psi {

AIOHandler::AIOHandler(std::shared_ptr<PSIO> psio)
    : psio_(psio)
{
    locked_ = new std::mutex();
    uniqueID_ = 0;
}
AIOHandler::~AIOHandler()
{
    synchronize();
    delete locked_;
}
void AIOHandler::synchronize()
{
    // Joining a thread twice is an error. Thus, we check whether
    // the thread is joinable before joining. Non-joinable threads
    // have either been joined, detached (no synchronization possible)
    // or have been moved from, and another thread object controls its execution.
//    std::unique_lock<std::mutex> lock(*locked_);
//    lock.unlock();
  if (thread_)
      if(thread_->joinable())
        thread_->join();
}
unsigned long int AIOHandler::read(unsigned int unit, const char *key, char *buffer, ULI size, psio_address start, psio_address *end)
{
  std::unique_lock<std::mutex> lock(*locked_);

  ++uniqueID_;
  job_.push(1);
  unit_.push(unit);
  key_.push(key);
  buffer_.push(buffer);
  size_.push(size);
  start_.push(start);
  end_.push(end);
  jobID_.push_back(uniqueID_);

  if (job_.size() > 1) return uniqueID_;

  // In C++11, destructors of threads that are still joinable call terminate()
  // Make sure current thread is joined before proceeding with next thread.
  synchronize();

  //thread start
  thread_ = std::make_shared<std::thread>(std::bind(&AIOHandler::call_aio,this));
  return uniqueID_;
}
unsigned long AIOHandler::write(unsigned int unit, const char *key, char *buffer, ULI size, psio_address start, psio_address *end)
{
  std::unique_lock<std::mutex> lock(*locked_);

  ++uniqueID_;
  job_.push(2);
  unit_.push(unit);
  key_.push(key);
  buffer_.push(buffer);
  size_.push(size);
  start_.push(start);
  end_.push(end);
  jobID_.push_back(uniqueID_);

  //printf("Adding a write to the queue\n");

  if (job_.size() > 1) return uniqueID_;

  //fprintf(stderr,"Starting a thread\n");
  //thread start

  // In C++11, destructors of threads that are still joinable call terminate()
  // Make sure current thread is joined before proceeding with next thread.
  synchronize();

  thread_ = std::make_shared<std::thread>(std::bind(&AIOHandler::call_aio,this));
  return uniqueID_;
}
unsigned long AIOHandler::read_entry(unsigned int unit, const char *key, char *buffer, ULI size)
{
  std::unique_lock<std::mutex> lock(*locked_);

  ++uniqueID_;
  job_.push(3);
  unit_.push(unit);
  key_.push(key);
  buffer_.push(buffer);
  size_.push(size);
  jobID_.push_back(uniqueID_);

  if (job_.size() > 1) return uniqueID_;

  // In C++11, destructors of threads that are still joinable call terminate()
  // Make sure current thread is joined before proceeding with next thread.
  synchronize();

  //thread start
  thread_ = std::make_shared<std::thread>(std::bind(&AIOHandler::call_aio,this));
  return uniqueID_;
}
unsigned long AIOHandler::write_entry(unsigned int unit, const char *key, char *buffer, ULI size)
{
  std::unique_lock<std::mutex> lock(*locked_);

  ++uniqueID_;
  job_.push(4);
  unit_.push(unit);
  key_.push(key);
  buffer_.push(buffer);
  size_.push(size);
  jobID_.push_back(uniqueID_);

  if (job_.size() > 1) return uniqueID_;

  // In C++11, destructors of threads that are still joinable call terminate()
  // Make sure current thread is joined before proceeding with next thread.
  synchronize();

  //thread start
  thread_ = std::make_shared<std::thread>(std::bind(&AIOHandler::call_aio,this));
  return uniqueID_;
}
unsigned long AIOHandler::read_discont(unsigned int unit, const char *key,
  double **matrix, ULI row_length, ULI col_length, ULI col_skip,
  psio_address start)
{
  std::unique_lock<std::mutex> lock(*locked_);

  ++uniqueID_;
  job_.push(5);
  unit_.push(unit);
  key_.push(key);
  matrix_.push(matrix);
  row_length_.push(row_length);
  col_length_.push(col_length);
  col_skip_.push(col_skip);
  start_.push(start);
  jobID_.push_back(uniqueID_);

  if (job_.size() > 1) return uniqueID_;

  // In C++11, destructors of threads that are still joinable call terminate()
  // Make sure current thread is joined before proceeding with next thread.
  synchronize();

  //thread start
  thread_ = std::make_shared<std::thread>(std::bind(&AIOHandler::call_aio,this));
  return uniqueID_;
}
unsigned long AIOHandler::write_discont(unsigned int unit, const char *key,
  double **matrix, ULI row_length, ULI col_length, ULI col_skip,
  psio_address start)
{
  std::unique_lock<std::mutex> lock(*locked_);

  ++uniqueID_;
  job_.push(6);
  unit_.push(unit);
  key_.push(key);
  matrix_.push(matrix);
  row_length_.push(row_length);
  col_length_.push(col_length);
  col_skip_.push(col_skip);
  start_.push(start);
  jobID_.push_back(uniqueID_);

  if (job_.size() > 1) return uniqueID_;

  // In C++11, destructors of threads that are still joinable call terminate()
  // Make sure current thread is joined before proceeding with next thread.
  synchronize();

  //thread start
  thread_ = std::make_shared<std::thread>(std::bind(&AIOHandler::call_aio,this));
  return uniqueID_;
}
unsigned long AIOHandler::zero_disk(unsigned int unit, const char *key,
    ULI rows, ULI cols)
{
  std::unique_lock<std::mutex> lock(*locked_);

  ++uniqueID_;
  job_.push(7);
  unit_.push(unit);
  key_.push(key);
  row_length_.push(rows);
  col_length_.push(cols);
  jobID_.push_back(uniqueID_);

  if (job_.size() > 1) return uniqueID_;

  // In C++11, destructors of threads that are still joinable call terminate()
  // Make sure current thread is joined before proceeding with next thread.
  synchronize();

  //thread start
  thread_ = std::make_shared<std::thread>(std::bind(&AIOHandler::call_aio,this));
  return uniqueID_;
}

unsigned long AIOHandler::write_iwl(unsigned int unit, const char *key,
              size_t nints, int lastbuf, char *labels, char *values,
              size_t labsize, size_t valsize, size_t *address) {
  std::unique_lock<std::mutex> lock(*locked_);
  ++uniqueID_;
  job_.push(8);
  unit_.push(unit);
  key_.push(key);
  buffer_.push(labels);
  buffer_.push(values);
  size_.push(labsize);
  size_.push(valsize);
  nints_.push(nints);
  lastbuf_.push(lastbuf);
  address_.push(address);
  jobID_.push_back(uniqueID_);

  if (job_.size() > 1) return uniqueID_;

  // In C++11, destructors of threads that are still joinable call terminate()
  // Make sure current thread is joined before proceeding with next thread.
  synchronize();

  //thread start
  thread_ = std::make_shared<std::thread>(std::bind(&AIOHandler::call_aio,this));
  return uniqueID_;

}

void AIOHandler::call_aio()
{
  std::unique_lock<std::mutex> lock(*locked_);

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
    else if (jobtype == 8) {

        lock.lock();

        unsigned int unit = unit_.front();
        const char* key = key_.front();
        char* labels = buffer_.front();
        buffer_.pop();
        char* values = buffer_.front();
        ULI lab_size = size_.front();
        size_.pop();
        ULI val_size = size_.front();
        int nints = nints_.front();
        int lastbuf = lastbuf_.front();
        size_t* address = address_.front();

        psio_address start = psio_get_address(PSIO_ZERO, *address);
        *address += val_size + lab_size + 2 * sizeof(int);

        unit_.pop();
        key_.pop();
        buffer_.pop();
        size_.pop();
        nints_.pop();
        lastbuf_.pop();
        address_.pop();

        lock.unlock();

        psio_->write(unit,key,(char*) &(lastbuf), sizeof(int),start,&start);
        psio_->write(unit,key,(char*) &(nints), sizeof(int), start, &start);
        psio_->write(unit,key,labels,lab_size,start,&start);
        psio_->write(unit,key,values,val_size,start,&start);

    }
    else {
      throw PsiException("Error in AIO: Unknown job type", __FILE__,__LINE__);
    }

    lock.lock();
    // Only pop the job once the work is actually done and we are gonna leave the loop.
    // This way, job_.size() == 0 indicates there is no active thread.
    job_.pop();
    // We also pop the jobID so that external threads may check that the job completed
    jobID_.pop_front();
    // Once it is popped, notify waiting threads to check again for their jobid.
    condition_.notify_all();
  }
  //printf("End of function call_aio\n");
}

void AIOHandler::wait_for_job(unsigned long jobid) {

    std::deque<unsigned long int>::iterator it;

    std::unique_lock<std::mutex> lock(*locked_);
    it = std::find(jobID_.begin(),jobID_.end(),jobid);
    bool found = it != jobID_.end();
    while(found) {
        condition_.wait(lock);
        it = std::find(jobID_.begin(),jobID_.end(),jobid);
        found = it != jobID_.end();
    }
    return;
}

} //Namespace psi
