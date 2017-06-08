

#include "psi4/libpsi4util/exception.h"
#include "PsiOutStream.h"
#include <cstdio>
#include <fstream>

namespace psi {

PsiOutStream::PsiOutStream(std::string fname, std::ios_base::openmode mode) {
    if (fname == "") {
        stream_ = &std::cout;
        is_cout_ = true;
    } else {
        std::ofstream* tmpf = new std::ofstream(fname, mode);
        if (!tmpf->is_open()) {
            throw PSIEXCEPTION("PsiOutStream: Failed to open file.");
        }

        stream_ = tmpf;
        is_cout_ = false;
    }

    buffer_.resize(512);
}

PsiOutStream::~PsiOutStream() {
    if (!is_cout_) {
        delete stream_;
    }
}

void PsiOutStream::Printf(const char* format, ...) {
    // We don't know how long the fully expanded string is so lets guess our average print is about
    // a line
    va_list args;
    va_start(args, format);
    int left = vsnprintf(buffer_.data(), buffer_.size(), format, args);

    if (left < 0) {
        // Encoding error?!?
        throw PSIEXCEPTION("PsiOutStream: vsnprintf encoding error!");
    } else if (left >= buffer_.size()) {
        // Buffer was too small! Try again
        std::vector<char> tmp_buffer(left + 1);
        left = vsnprintf(tmp_buffer.data(), left + 1, format, args);
        if (left < 0) {
            throw PSIEXCEPTION("PsiOutStream: vsnprintf encoding error!");
        }
    }
    // Everything is cool

    va_end(args);
    (*stream_) << buffer_.data();
}

PsiOutStream& PsiOutStream::operator<<(std::string fp) {
    (*stream_) << fp;
    return (*this);
}

} // End Psi Namespace
