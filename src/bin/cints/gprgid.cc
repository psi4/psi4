/*! \file gprgid.cc
  \ingroup (CINTS)
  \brief program id
*/
extern "C" {
  char *gprgid()
  {
    char *prgid = "cints";
    
    return(prgid);
  }
}
