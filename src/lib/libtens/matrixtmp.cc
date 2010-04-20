#include "matrix.h"
#include "matrixtmp.h"

namespace psi{ namespace psimrcc{

CCMatTmp::CCMatTmp(CCMatrix* Matrix,DiskOpt disk_option):Matrix_(Matrix),disk_option_(disk_option)
{
}

CCMatTmp::~CCMatTmp()
{
  if(disk_option_ == dump)
    Matrix_->dump_to_disk();
  else if(disk_option_ == release)
    Matrix_->free_memory();
}

CCMatIrTmp::CCMatIrTmp(CCMatrix* Matrix,int irrep,DiskOpt disk_option):Matrix_(Matrix),irrep_(irrep),disk_option_(disk_option)
{
}

CCMatIrTmp::~CCMatIrTmp()
{
  if(disk_option_ == dump)
    Matrix_->dump_to_disk();
  else if(disk_option_ == release)
    Matrix_->free_memory();
}

}} /* End Namespaces */
