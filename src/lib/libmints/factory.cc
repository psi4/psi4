/*
 *  factory.cpp
 *  matrix
 *
 *  Created by Justin Turney on 4/8/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "factory.h"

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

bool MatrixFactory::init_with_chkpt(PSIO* psio)
{
    Chkpt* chkpt(new Chkpt(psio, PSIO_OPEN_OLD));
    bool result = init_with_chkpt(chkpt);
    delete chkpt;
    return result;
}

bool MatrixFactory::init_with_chkpt(PSIO& psio)
{
    Chkpt chkpt(psio, PSIO_OPEN_OLD);
    bool result = init_with_chkpt(chkpt);
    return result;
}

bool MatrixFactory::init_with_chkpt(Chkpt* chkpt)
{
    nirreps_ = chkpt->rd_nirreps();
    rowspi_  = chkpt->rd_sopi();
    colspi_  = chkpt->rd_sopi();
    nso_     = chkpt->rd_nso();
    
    return true;
}

bool MatrixFactory::init_with_chkpt(Chkpt& chkpt)
{
    nirreps_ = chkpt.rd_nirreps();
    rowspi_  = chkpt.rd_sopi();
    colspi_  = chkpt.rd_sopi();
    nso_     = chkpt.rd_nso();
    
    return true;
}

bool MatrixFactory::init_with(int nirreps, int *rowspi, int *colspi)
{
    nirreps_ = nirreps;
    rowspi_ = new int[nirreps_];
    colspi_ = new int[nirreps_];
        
    for (int i=0; i<nirreps_; ++i) {
        rowspi_[i] = rowspi[i];
        colspi_[i] = colspi[i];
    }
    
    return true;
}

