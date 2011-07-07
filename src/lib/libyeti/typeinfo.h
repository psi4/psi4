#ifndef yeti_typeinfo_h
#define yeti_typeinfo_h

namespace yeti {

/**
 * @class MatMultIterSelect
 */
template <class T>
class constPtrType {

    public:
        typedef boost::intrusive_ptr<const T> type;

};

template <class T>
class constPtrType<const T>
{
    public:
        typedef boost::intrusive_ptr<const T> type;
};

}

#endif
