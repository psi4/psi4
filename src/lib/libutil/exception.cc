#include <exception.h>

using namespace psi;
using namespace std;

PsiException::PsiException(
    string msg,
    const char* file,
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
    return file_;
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
    string message,
    const char* file,
    int line
    ) throw() 
  : PsiException(message, file, line)
{
}

StepSizeError::StepSizeError(
    string value_name,
    double max,
    double actual,
    const char* file,
    int line) throw()
    : ParentClass(value_name, max, actual, file, line)
{
}

MaxIterationsExceeded::MaxIterationsExceeded(
    string routine_name,
    int max,
    const char* file,
    int line)  throw()
    : ParentClass(routine_name, max, max + 1, file, line)
{
}

ConvergenceError::ConvergenceError(
    string routine_name,
    int max,
    double desired_accuracy,
    double actual_accuracy,
    const char* file,
    int line) throw()
    : MaxIterationsExceeded(routine_name, max, file, line), desired_acc_(desired_accuracy), actual_acc_(actual_accuracy)
{
}

double 
ConvergenceError::desired_accuracy() const throw() {return desired_acc_;}

double 
ConvergenceError::actual_accuracy() const throw() {return actual_acc_;}

ResourceAllocationError::ResourceAllocationError(
    string resource_name,
    size_t max,
    size_t actual,
    const char* file,
    int line)  throw()
    : LimitExceeded<size_t>(resource_name, max, actual, file, line)
{
}

