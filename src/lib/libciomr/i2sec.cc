/* converts integer file pointer to sector pointer */

#include "iomrparam.h"

extern "C" {

int i2sec(PSI_FPTR n)
    {
      int num;

      num = (int) ((n+4096)/4096);
      return(num);
    }

} /* extern "C" */
