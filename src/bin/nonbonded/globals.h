#ifdef EXTERN
# undef EXTERN
# define EXTERN extern
#else
# define EXTERN
#endif

namespace psi { namespace nonbonded {

extern "C" {
  EXTERN FILE *infile,*outfile;
  EXTERN char *psi_file_prefix;
}

}} // namespace psi::nonbonded

