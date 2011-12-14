
#ifndef _gigide_exception_h_
#define _gigide_exception_h_

#include <exception>
#include <stdexcept>
#include <sstream>

namespace yeti {

#define CHARARR_SIZE 100
#define except(message) throw YetiException(message, __FILE__, __LINE__)

void
exception_crash(
    const std::string& excname,
    const std::string& msg,
    const std::string& file,
    int line
);

//#define raise(exc, ...) throw exc(__VA_ARGS__, __FILE__, __LINE__)
#define yeti_throw(exc, ...) exception_crash(#exc, __VA_ARGS__, __FILE__, __LINE__)

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

/**
    Generic exception class for Psi4
*/
class YetiException : public std::runtime_error {

    private:
        std::string msg_;
        std::string file_;
        int line_;

    protected:
        /**
        * Override default message for exception throw in what
        * @param msg The message for what to throw
        */
        void rewrite_msg(const std::string& msg) throw();

    public:
        /**
        * Constructor
        * @param message The message that will be printed by exception
        * @param file The file that threw the exception (use __FILE__ macro)
        * @param line The line number that threw the exception (use __LINE__ macro)
        */
        YetiException(
            const std::string& message,
            const std::string& file,
            int line
        ) throw ();

        ~YetiException() throw ();

        /**
        * Override of runtime_error's virtual what method
        * @return Description of exception
        */
        const char* what() const throw ();

        /**
        * Accessor method
        * @return File that threw the exception
        */
        std::string file() const throw();

        /**
        * Accessor method
        * @return A string description of line and file that threw exception
        */
        std::string location() const throw();

        /**
        * Accessor method
        * @return The line number that threw the exception
        */
        int line() const throw();

};

/**
* Exception for sanity checks being performed, e.g. checking alignment of matrix multiplication.
*/
class SanityCheckError : public YetiException {

    public:
        /**
        * Constructor
        * @param message The message that will be printed by exception
        * @param file The file that threw the exception (use __FILE__ macro)
        * @param line The line number that threw the exception (use __LINE__ macro)
        */
        SanityCheckError(
            const std::string& message,
            const std::string& file,
            int line
        ) throw();

        ~SanityCheckError() throw();

};

/**
* Exception for equality checks
*/
template <class T>
class ValuesNotEqual : public SanityCheckError {

    private:
        T correctval_;
        T errorval_;
        std::string resource_name_;

    protected:
        /**
        * Accessor method
        * @return A string description of the limit that was exceeded
        */
        std::string description() const throw()
        {
            std::stringstream sstr;
            sstr << "value for " << resource_name_ << " is wrong.\n"
                 << "allowed: " << correctval_ << " actual: " << errorval_;
            return sstr.str();
        }

    public:
        /**
        * Constructor
        * @param resource_name The name of the value that was exceeded (e.g. memory or scf max iterations)
        * @param maxval The max (or min) value allowed
        * @param errorval The actual value obtained
        * @param file The file that threw the exception (use __FILE__ macro)
        * @param line The line number that threw the exception (use __LINE__ macro)
        */
        ValuesNotEqual(
            const std::string& resource_name,
            T correctval,
            T errorval,
            const char* file,
            int line) throw() : SanityCheckError(resource_name, file, line), correctval_(correctval), errorval_(errorval), resource_name_(resource_name)
        {
            rewrite_msg(description());
        }

        /** Accessor method 
        *  @return The max (or min) value allowed
        */
        T correct_value() const throw() {return correctval_;}

        /** Accessor method 
        *  @return The actual value
        */
        T actual_value() const throw() {return errorval_;}

        ~ValuesNotEqual() throw(){}
};

/**
* Exception for features that are not implemented yet, but will be (maybe?)
*/
class FeatureNotImplemented : public YetiException {

    public:
        /**
        * Constructor
        * @param module The module being run
        * @param feature The feature not yet implemented
        * @param file The file that threw the exception (use __FILE__ macro)
        * @param line The line number that threw the exception (use __LINE__ macro)
        */
        FeatureNotImplemented(
            const std::string& module, 
            const std::string& feature,
            const std::string& file,
            int line
        ) throw();

        ~FeatureNotImplemented() throw();
};

/**
* A generic template class for exceptions in which a min or max value is exceed
*/
template <class T>
class LimitExceeded : public YetiException {

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
        * @param file The file that threw the exception (use __FILE__ macro)
        * @param line The line number that threw the exception (use __LINE__ macro)
        */
        LimitExceeded(
            const std::string& resource_name,
            T maxval,
            T errorval,
            const std::string& file,
            int line) throw() : YetiException(resource_name, file, line), maxval_(maxval), errorval_(errorval), resource_name_(resource_name)
        {
            rewrite_msg(description());
        }

        /** Accessor method 
        *  @return The max (or min) value allowed
        */
        T max_value() const throw() {return maxval_;}

        /** Accessor method 
        *  @return The actual value
        */
        T actual_value() const throw() {return errorval_;}

        ~LimitExceeded() throw(){}
};

/**
* Convergence error for routines
*/
class ResourceAllocationError :
    public LimitExceeded<size_t> {

    public:
        /**
        * Constructor
        * @param resource_name The name of the resource (e.g memory)
        * @param max The maximum resource size
        * @param actual The resource size you tried to allocate
        * @param file The file that threw the exception (use __FILE__ macro)
        * @param line The line number that threw the exception (use __LINE__ macro)
        */
        ResourceAllocationError(
            const std::string& resource_name,
            size_t max,
            size_t actual,
            const std::string& file,
            int line) throw ();

        ~ResourceAllocationError() throw();
};

/**
* Exception on input values
*/
class InputException : public YetiException {

    private:
        /**
        * Template method for writing generic input exception message
        */
        template <class T> void
        write_input_msg(const std::string& msg, const std::string& param_name, T val) throw();

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
            const std::string& msg,
            const std::string& param_name,
            int value,
            const std::string& file,
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
            const std::string& msg,
            const std::string& param_name,
            double value,
            const std::string& file,
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
            const std::string& msg,
            const std::string& param_name,
            const std::string& value,
            const std::string& file,
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
            const std::string& msg,
            const std::string& param_name,
            const std::string& file,
            int line
        ) throw();
};

} //end namespace sc

#ifdef redefine_size_t
#undef size_t
#endif

#endif
