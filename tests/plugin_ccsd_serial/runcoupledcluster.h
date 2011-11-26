#ifndef RUNCOUPLEDCLUSTER_H
#define RUNCOUPLEDCLUSTER_H

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

// psi headers
#include <libplugin/plugin.h>
#include"psi4-dec.h"
#include<boost/shared_ptr.hpp>
#include<liboptions/liboptions.h>
#include<libtrans/integraltransform.h>
#include<libtrans/mospace.h>
#include<libmints/matrix.h>
#include<libmints/vector.h>
#include<libchkpt/chkpt.h>
#include<libiwl/iwl.h>
#include <libpsio/psio.hpp>

#include"gpuhelper.h"

namespace boost {
template<class T> class shared_ptr;
}

// gpu ccsd class
namespace psi{
  void RunCoupledCluster(Options &options);
};

#endif
