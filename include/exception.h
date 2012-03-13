#ifndef _psi4_exception_h_
#define _psi4_exception_h_

#include <exception>
#include <stdexcept>
#include <sstream>
#include <string.h>
#include <boost/current_function.hpp>

namespace psi {

#define CHARARR_SIZE 100
#define PSIEXCEPTION(message) PsiException(message, __FILE__, __LINE__)
#define NotImplementedException() NotImplementedException_(BOOST_CURRENT_FUNCTION, __FILE__, __LINE__)

/**
    Generic exception class for Psi4
*/
class PsiException : public std::runtime_error {

    private:
        std::string msg_;
        const char* file_;
        int line_;

    protected:
        /**
        * Override default message for exception throw in what
        * @param msg The message for what to throw
        */
        void rewrite_msg(std::string msg) throw();

    public:
        /**
        * Constructor
        * @param message The message that will be printed by exception
        * @param file The file that threw the exception (use __FILE__ macro)
        * @param line The line number that threw the exception (use __LINE__ macro)
        */
        PsiException(
            std::string message,
            const char* file,
            int line
        ) throw ();
        PsiException(const PsiException& copy) throw();
        virtual ~PsiException() throw ();

        PsiException& operator=(const PsiException& other) {
            if (this != &other) {
                msg_ = other.msg_;
                file_ = strdup(other.file_);
                line_ = other.line_;
            }
            return *this;
        }

        /**
        * Override of runtime_error's virtual what method
        * @return Description of exception
        */
        const char* what() const throw ();

        /**
        * Accessor method
        * @return File that threw the exception
        */
        const char* file() const throw();

        /**
        * Accessor method
        * @return A string description of line and file that threw exception
        */
        const char* location() const throw();

        /**
        * Accessor method
        * @return The line number that threw the exception
        */
        int line() const throw();

};

class NotImplementedException_ : public PsiException
{
public:
    NotImplementedException_(const std::string& message,
                             const char* file,
                             int line)
        : PsiException(message + " function not implemented", file, line) { }
};

/**
* Exception for sanity checks being performed, e.g. checking alignment of matrix multiplication.
*/
class SanityCheckError : public PsiException {

    public:
        /**
        * Constructor
        * @param message The message that will be printed by exception
        * @param file The file that threw the exception (use __FILE__ macro)
        * @param line The line number that threw the exception (use __LINE__ macro)
        */
        SanityCheckError(
            std::string message,
            const char* file,
            int line
        ) throw();

        virtual ~SanityCheckError() throw();
};

/**
* Exception for system errors. Those that report via errno
*/
class SystemError : public PsiException {

    public:
        /**
        * Constructor
        * @param eno The errno that will be printed by exception
        * @param file The file that threw the exception (use __FILE__ macro)
        * @param line The line number that threw the exception (use __LINE__ macro)
        */
        SystemError(
            int eno,
            const char* file,
            int line
        ) throw();

        virtual ~SystemError() throw();
};

/**
* Exception for features that are not implemented yet, but will be (maybe?)
*/
class FeatureNotImplemented : public PsiException {

    public:
        /**
        * Constructor
        * @param module The module being run
        * @param feature The feature not yet implemented
        * @param file The file that threw the exception (use __FILE__ macro)
        * @param line The line number that threw the exception (use __LINE__ macro)
        */
        FeatureNotImplemented(
            std::string module,
            std::string feature,
            const char* file,
            int line
        ) throw();

        virtual ~FeatureNotImplemented() throw();
};


/**
* A generic template class for exceptions in which a min or max value is exceed
*/
template <class T>
class LimitExceeded : public PsiException {

    private:
        T maxval_;
        T errorval_;
        std::string resource_name_;

    protected:
        /**
        * Accessor method
        * @return A string description of the limit that was exceeded
        */
        const char* description() const throw()
        {
            std::stringstream sstr;
            sstr << "value for " << resource_name_ << " exceeded.\n"
                 << "allowed: " << maxval_ << " actual: " << errorval_;
            return sstr.str().c_str();
        }

    public:
        /**
        * Constructor
        * @param resource_name The name of the value that was exceeded (e.g. memory or scf max iterations)
        * @param maxval The max (or min) value allowed
        * @param errorval The actual value obtained
        * @param f The file that threw the exception (use __FILE__ macro)
        * @param l The line number that threw the exception (use __LINE__ macro)
        */
        LimitExceeded(
            std::string resource_name,
            T maxval,
            T errorval,
            const char* f,
            int l) throw() : PsiException(resource_name, f, l), maxval_(maxval), errorval_(errorval), resource_name_(resource_name)
        {
            rewrite_msg(description());
        }

        T max_value() const throw() {return maxval_;}

        T actual_value() const throw() {return errorval_;}

        virtual ~LimitExceeded<T>() throw() {};
};

/**
* Error in a step size
*/
template <class T=double>
class StepSizeError : public LimitExceeded<T> {

    typedef LimitExceeded<T> ParentClass;

    public:
        /**
        * Constructor
        * @param resource_name The name of the value that is changed (scf, geometry, cc amps)
        * @param max The max (or min) value allowed
        * @param actual The actual value obtained
        * @param file The file that threw the exception (use __FILE__ macro)
        * @param line The line number that threw the exception (use __LINE__ macro)
        */
        StepSizeError(
            std::string resource_name,
            T max,
            T actual,
            const char* file,
            int line) throw();

        virtual ~StepSizeError() throw();
};

/**
* Maximum number of iterations exceeded
*/
template <class T=int>
class MaxIterationsExceeded : public LimitExceeded<T> {

    typedef LimitExceeded<T> ParentClass;

    public:
        /**
        * Constructor
        * @param routine_name The name of the routine that is not converging (e.g. scf, ccsd, optking disp)
        * @param max The max number of iterations
        * @param file The file that threw the exception (use __FILE__ macro)
        * @param line The line number that threw the exception (use __LINE__ macro)
        */
        MaxIterationsExceeded(
            std::string routine_name,
            T max,
            const char* file,
            int line) throw();

        virtual ~MaxIterationsExceeded() throw();

};

/**
* Convergence error for routines
*/
template <class T=int>
class ConvergenceError : public MaxIterationsExceeded<T> {
    public:
        /**
        * Constructor
        * @param routine_name The name of the routine that is not converging (e.g. scf, ccsd, optking disp)
        * @param max The max number of iterations
        * @param desired_accuracy The accuracy you want to achieve
        * @param actual_accuracy The actual accuracy achieved
        * @param file The file that threw the exception (use __FILE__ macro)
        * @param line The line number that threw the exception (use __LINE__ macro)
        */
        ConvergenceError(
            std::string routine_name,
            T max,
            double desired_accuracy,
            double actual_accuracy,
            const char* file,
            int line) throw();

        /** Accessor method
        *  @return The accuracy you wish to achieve
        */
        double desired_accuracy() const throw();

        /** Accessor method
        *  @return The accuracy you actually got
        */
        double actual_accuracy() const throw();

        virtual ~ConvergenceError() throw();

    private:
        double desired_acc_;
        double actual_acc_;
};

/**
* Convergence error for routines
*/
template<class T=size_t>
class ResourceAllocationError : public LimitExceeded<T> {

    public:
        /**
        * Constructor
        * @param resource_name The name of the resource (e.g memory)
        * @param max The max number of iterations
        * @param file The file that threw the exception (use __FILE__ macro)
        * @param line The line number that threw the exception (use __LINE__ macro)
        */
        ResourceAllocationError(
            std::string resource_name,
            T max,
            T actual,
            const char* file,
            int line) throw ();

        virtual ~ResourceAllocationError() throw();
};

/**
* Exception on input values
*/
class InputException : public PsiException {

    private:
        /**
        * Template method for writing generic input exception message
        */
        template <class T> void
        write_input_msg(std::string msg, std::string param_name, T val) throw();

    public:
        /**
        * Constructor
        * @param msg A message descring why the input parameter is invalid
        * @param param_name The parameter in the input file that needs to be changed
        * @param value The value that proved to be incorrect
        * @param file The file that threw the exception (use __FILE__ macro)
        * @param line The line number that threw the exception (use __LINE__ macro)
        */
        InputException(
            std::string msg,
            std::string param_name,
            int value,
            const char* file,
            int line
        ) throw();

        /**
        * Constructor
        * @param msg A message descring why the input parameter is invalid
        * @param param_name The parameter in the input file that needs to be changed
        * @param value The value that proved to be incorrect
        * @param file The file that threw the exception (use __FILE__ macro)
        * @param line The line number that threw the exception (use __LINE__ macro)
        */
        InputException(
            std::string msg,
            std::string param_name,
            double value,
            const char* file,
            int line
        ) throw();

        /**
        * Constructor
        * @param msg A message descring why the input parameter is invalid
        * @param param_name The parameter in the input file that needs to be changed
        * @param value The value that proved to be incorrect
        * @param file The file that threw the exception (use __FILE__ macro)
        * @param line The line number that threw the exception (use __LINE__ macro)
        */
        InputException(
            std::string msg,
            std::string param_name,
            std::string value,
            const char* file,
            int line
        ) throw();

        /**
        * Constructor
        * @param msg A message descring why the input parameter is invalid
        * @param param_name The parameter in the input file that needs to be changed
        * @param file The file that threw the exception (use __FILE__ macro)
        * @param line The line number that threw the exception (use __LINE__ macro)
        */
        InputException(
            std::string msg,
            std::string param_name,
            const char* file,
            int line
        ) throw();
};

} //end namespace psi exception

#endif
