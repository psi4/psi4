#include <exception>

namespace psi {

#define CHARARR_SIZE 100

class PsiException : public std::exception {

    private:
        const char* msg_;
        const char* file_;
        int line_;

    public:
        PsiException(
            const char* message,
            const char* file,
            int line
        ) throw ()
        {
            msg_ = message;
            file_ = file;
            line_ = line;
        }

        ~PsiException() throw ()
        {
        }

        const char* what() const throw ()
        {
            return msg_;
        }

        const char* file() const throw()
        {
            return file_;
        }

        int line() const throw()
        {
            return line_;
        }


};

class SanityCheckError : public PsiException {

    public:
        SanityCheckError(
            const char* message,
            const char* file,
            int line
            ) : PsiException(message, file, line)
        {
        }
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
            const char* msg,
            T maxval,
            T errorval,
            const char* file,
            int line) : 
           PsiException(msg, file, line),
           maxval_(maxval), errorval_(errorval)
        {
        }

        T max_value(){return maxval_;}

        T actual_value(){return errorval_;}
};

class StepSizeError : public LimitExceeded<double> {

    typedef LimitExceeded<double> ParentClass;

    public:
        StepSizeError(
            const char* msg,
            double max,
            double actual,
            const char* file,
            int line)
            : ParentClass(msg, max, actual, file, line)
        {
        }

};


class MaxIterationsExceeded : public LimitExceeded<int> {

    typedef LimitExceeded<int> ParentClass;

    public:
        MaxIterationsExceeded(
            const char* msg,
            int max,
            const char* file,
            int line) :
            ParentClass(msg, max, max + 1, file, line)
        {
        }
};

class ConvergenceError : public MaxIterationsExceeded {

    public:
        ConvergenceError(
            const char* msg,
            int max,
            double desired_accuracy,
            double actual_accuracy,
            const char* file,
            int line) : 
            MaxIterationsExceeded(msg, max, file, line), desired_acc_(desired_accuracy), actual_acc_(actual_accuracy)
        {
        }

        double desired_accuracy(){return desired_acc_;}

        double actual_accuracy(){return actual_acc_;}

    private:
        double desired_acc_;
        double actual_acc_;


};

} //end namespace psi exception

