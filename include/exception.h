#ifndef _psi4_exception_h_
#define _psi4_exception_h_

#include <exception>
#include <sstream>
#include <stdexcept>

namespace psi {

#define CHARARR_SIZE 100
#define PSIEXCEPTION(message) PsiException(message, __FILE__, __LINE__)

class PsiException : public std::runtime_error {

    private:
        std::string msg_;
        const char* file_;
        int line_;

    public:
        PsiException(
            std::string message,
            const char* file,
            int line
        ) throw ();

        ~PsiException() throw ();

        const char* what() const throw ();

        const char* file() const throw();

        const char* location() const throw();

        int line() const throw();

};

class SanityCheckError : public PsiException {

    public:
        SanityCheckError(
            std::string message,
            const char* file,
            int line
        ) throw();

};

class FeatureNotImplemented : public PsiException {

    public:
        FeatureNotImplemented(
            std::string module, 
            std::string feature,
            const char* file,
            int line
        ) throw();
};

template <
    class T
>
class LimitExceeded : public PsiException {

    private:
        T maxval_;
        T errorval_;

    public:
        LimitExceeded(
            std::string resource_name,
            T maxval,
            T errorval,
            const char* file,
            int line) throw() : PsiException(resource_name, file, line)
        {
        }

        T max_value() const throw() {return maxval_;}

        T actual_value() const throw() {return errorval_;}
};

class StepSizeError : public LimitExceeded<double> {

    typedef LimitExceeded<double> ParentClass;

    public:
        StepSizeError(
            std::string msg,
            double max,
            double actual,
            const char* file,
            int line) throw();
};


class MaxIterationsExceeded : public LimitExceeded<int> {

    typedef LimitExceeded<int> ParentClass;

    public:
        MaxIterationsExceeded(
            std::string routine_name,
            int max,
            const char* file,
            int line) throw();

};

class ConvergenceError : public MaxIterationsExceeded {

    public:
        ConvergenceError(
            std::string routine_name,
            int max,
            double desired_accuracy,
            double actual_accuracy,
            const char* file,
            int line) throw();

        double desired_accuracy() const throw();

        double actual_accuracy() const throw();

    private:
        double desired_acc_;
        double actual_acc_;


};

class ResourceAllocationError : public LimitExceeded<size_t> {

    public:
        ResourceAllocationError(
            std::string resource_name,
            size_t max,
            size_t actual,
            const char* file,
            int line) throw ();
};

} //end namespace psi exception

#endif
