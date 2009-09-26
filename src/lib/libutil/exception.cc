#include <exception.h>

using namespace psi;
using namespace std;

PsiException::PsiException(
    string msg,
    string file,
    int line
) throw() : runtime_error(msg)
{
    msg_ = msg;
    file_ = file;
    line_ = line;
}
    
const char* 
PsiException::what() const throw()
{
    return msg_.c_str();
}

const char* 
PsiException::file() const throw()
{
    return file_.c_str();
}

int 
PsiException::line() const throw()
{
    return line_;
}

PsiException::~PsiException() throw()
{
}

SanityCheckError::SanityCheckError(
    const char* message,
    const char* file,
    int line
    ) : PsiException(message, file, line)
{
}

template <
    class T
>
LimitExceeded<T>::LimitExceeded(
    const char* msg,
    T maxval,
    T errorval,
    const char* file,
    int line) : 
   PsiException(msg, file, line),
   maxval_(maxval), errorval_(errorval)
{
}


template <
    class T
> T 
LimitExceeded<T>::max_value(){return maxval_;}

template <
    class T
> T 
LimitExceeded<T>::actual_value(){return errorval_;}

StepSizeError::StepSizeError(
    const char* msg,
    double max,
    double actual,
    const char* file,
    int line)
    : ParentClass(msg, max, actual, file, line)
{
}

MaxIterationsExceeded::MaxIterationsExceeded(
    const char* msg,
    int max,
    const char* file,
    int line) :
    ParentClass(msg, max, max + 1, file, line)
{
}

ConvergenceError::ConvergenceError(
    const char* msg,
    int max,
    double desired_accuracy,
    double actual_accuracy,
    const char* file,
    int line) : 
    MaxIterationsExceeded(msg, max, file, line), desired_acc_(desired_accuracy), actual_acc_(actual_accuracy)
{
}

double 
ConvergenceError::desired_accuracy(){return desired_acc_;}

double 
ConvergenceError::actual_accuracy(){return actual_acc_;}

