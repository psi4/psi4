/*
 *  factory.cpp
 *  matrix
 *
 *  Created by Justin Turney on 4/8/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "factory.h"
#include <libparallel/parallel.h>
#include "libciomr/libciomr.h"


 using namespace psi;
 
MatrixFactory::MatrixFactory()
{
    nirreps_ = 0;
    rowspi_ = NULL;
    colspi_ = NULL;
}

MatrixFactory::MatrixFactory(const MatrixFactory& copy)
{
    nirreps_ = copy.nirreps_;
    rowspi_ = new int[nirreps_];
    colspi_ = new int[nirreps_];
    
    memcpy(rowspi_, copy.rowspi_, sizeof(int) * nirreps_);
    memcpy(colspi_, copy.colspi_, sizeof(int) * nirreps_);
}

MatrixFactory::~MatrixFactory()
{
    if (rowspi_)
        Chkpt::free(rowspi_);
    if (colspi_)
        Chkpt::free(colspi_);
}

bool MatrixFactory::init_with_chkpt(shared_ptr<PSIO> psio)
{
    shared_ptr<Chkpt> chkpt(new Chkpt(psio.get(), PSIO_OPEN_OLD));
    bool result = init_with_chkpt(chkpt);
    return result;
}

bool MatrixFactory::init_with_chkpt(shared_ptr<Chkpt> chkpt)
{
    if(Communicator::world->me() == 0) {
        nirreps_ = chkpt->rd_nirreps();
        rowspi_  = chkpt->rd_sopi();
        colspi_  = chkpt->rd_sopi();
        nso_     = chkpt->rd_nso();
    }
    if(Communicator::world->nproc() > 1) {
        Communicator::world->raw_bcast(&nirreps_, sizeof(int), 0);
        Communicator::world->raw_bcast(&nso_, sizeof(int), 0);
        if(Communicator::world->me() != 0) {
            rowspi_ = init_int_array(nirreps_);
            colspi_ = init_int_array(nirreps_);
        }
        Communicator::world->raw_bcast(&(rowspi_[0]), nirreps_*sizeof(int), 0);
        Communicator::world->raw_bcast(&(colspi_[0]), nirreps_*sizeof(int), 0);
    }

    return true;
}

bool MatrixFactory::init_with(int nirreps, int *rowspi, int *colspi)
{
    nirreps_ = nirreps;
    rowspi_ = new int[nirreps_];
    colspi_ = new int[nirreps_];
        
    nso_ = 0;
    for (int i=0; i<nirreps_; ++i) {
        rowspi_[i] = rowspi[i];
        colspi_[i] = colspi[i];
        nso_ += rowspi_[i];
    }
    
    return true;
}

