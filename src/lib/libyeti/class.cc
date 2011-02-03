#include "class.h"
#include "yeti.h"

using namespace yeti;
using namespace std;

const char* TypeInfo<double>::printf_str = "%12.8f";
const char* TypeInfo<float>::printf_str = "%12.8f";
const char* TypeInfo<quad>::printf_str = "%Lf";
const char* TypeInfo<int>::printf_str = "%6d";

#ifdef CUSTOM_INTRUSIVE_PTR_FXN
void
boost::intrusive_ptr_add_ref(smartptr::Countable *obj)
{
    if (YetiRuntime::is_threaded_runtime())
    {
        cerr << "Code is not thread-safe! Non-thread safe object is being modified"
                " in threaded routine!" << endl;
        abort();
    }
    else
    {
        obj->incref();
    }
}

void
boost::intrusive_ptr_add_ref(const smartptr::Countable *obj)
{
    intrusive_ptr_add_ref(const_cast<smartptr::Countable*>(obj));
}


void
boost::intrusive_ptr_release(smartptr::Countable *obj)
{
    if (YetiRuntime::is_threaded_runtime())
    {
        cerr << "Code is not thread-safe! Non-thread safe object is being modified"
                " in threaded routine!" << endl;
        abort();
    }
    else
    {
        if (obj->decref() == 0)
            delete obj;
    }
}

void
boost::intrusive_ptr_release(const smartptr::Countable *obj)
{
    intrusive_ptr_release(const_cast<smartptr::Countable*>(obj));
}

#endif


