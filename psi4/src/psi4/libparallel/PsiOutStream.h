

#ifndef _psi_src_lib_libparallel_psioutstream_h_
#define _psi_src_lib_libparallel_psioutstream_h_

#include <vector>
#include <string>
#include <iostream>

namespace psi {

class PsiOutStream {
   private:
        std::ostream *stream_;
        bool is_cout_;
        std::vector<char> buffer_;

    public:
        PsiOutStream(std::string fname = "", std::ios_base::openmode mode = std::ostream::trunc);
        ~PsiOutStream();

        void Printf(const char * fmt, ...);
        PsiOutStream& operator<<(std::string fp);
        void MakeBanner(std::string header);

};

} // End Psi namespace
#endif
