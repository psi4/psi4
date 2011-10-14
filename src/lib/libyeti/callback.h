#ifndef yeti_CALLBACK_H
#define yeti_CALLBACK_H

namespace yeti {

class Callback {

    public:
        virtual void callback() = 0;
};


template <class T>
class tmpl_Callback :
    public Callback
{
    private:
        typedef void (T::*void_fxn_ptr)(void);

        T* obj_;

        void_fxn_ptr fxn_;

    public:
        tmpl_Callback(
            T* obj,
            void_fxn_ptr fxn
        ) : obj_(obj),
            fxn_(fxn)
        {
        }

        void callback()
        {
            ((*obj_).*fxn_)();
        }

        void* operator new(size_t size, char* ptr)
        {
            return ptr;
        }

        void operator delete(void* ptr, size_t size)
        {
            //do nothing
        }
};

}

#endif // CALLBACK_H
