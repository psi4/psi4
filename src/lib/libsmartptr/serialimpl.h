#ifndef smartptr_serialimpl_h
#define smartptr_serialimpl_h

#include "xmlarchive.h"
#include "timer.h"

/**
    @def serial_save(x)
    Macro for serializing objects.  This assumes that an archive exists in the
    local namespace named arch and that the member variable is named x_.
*/
#define serial_save(x)        smartptr::serial_call_save(arch, x##_, #x)

/**
    @def serial_save_enum(x)
    Macro for serializing enum types. These must be treated separately since the enum
    type must be written as an unsigned int. This assumes that an archive exists in the
    local namespace named arch and that the member variable is named x_.
*/
#define serial_save_enum(x) unsigned int tmp = (unsigned int) x##_; arch->setValue(tmp, #x);

/**
    @def serial_load(x)
    Macro for deserializing objects.  This assumes that an archive exists in the
    local namespace named arch and that the member variable is named x_.
*/
#define serial_load(x)        smartptr::serial_call_load(arch, x##_, #x)

/**
    @def serial_load_enum(x)
    Macro for deserializing enum types. These must be treated separately since the enum
    type must be written as an unsigned int. This assumes that an archive exists in the
    local namespace named arch and that the member variable is named x_.
*/
#define serial_load_enum(x) unsigned int tmp; arch->getValue(tmp, #x); arch->enum_cast(tmp, x##_);

/**
    @def SetRuntime(x)
    Macro for configuring the runtime type of an object.  This takes the classname
    as the parameter.  A compile time check is included to ensure proper configuration.
*/
#define SetRuntime(x) rtinfo_ = &x##_sc; x* const compile_time_check = this; 

/**
    @def SerialDeclare(x)
    Macro for configuring static runtime info associated with a serializble type. This takes the
    classname as the parameter.  This should be included at the top of the .cc file in which the
    class constructors are located.
*/
#define SerialDeclare(x) \
template<> std::string smartptr::SerialClass<x>::classname_ = ""; \
template<> smartptr::SerialRuntime::create_function* smartptr::SerialClass<x>::fxn_ = smartptr::build<x>; \
smartptr::SerialClass<x> x##_sc(#x)

namespace smartptr {

/**
    @class SerialDecide
    Type-decision class for choosing the appropriate serialization
    method.  Generally speaking, this branches serialization calls depending on
    primitive types, STL types, and Serializable derived types.
*/
template <class T>
class SerialDecide {

    public:
        /**
            Method for serializing arbitrary types. This creates
            the xml node, writes the value, and returns to original archive
            position.
            @param arch
            @param val The object to serialize
            @param tagname The XML tagname to create
        */
        static void
        serialize(
            const XMLArchivePtr& arch,
            const T& val,
            const std::string& tagname
        )
        {
            AssertNotHere<> test;
        }

        /**
            Method for serializing arbitrary types. This writes the value
            at the current XML node.
            @param arch
            @param val The object to serialize
        */
        static void
        serialize(
            const XMLArchivePtr& arch,
            const T& val
        )
        {
            AssertNotHere<> test;
        }

        /**
            Method for deserializing arbitrary types. This creates
            the xml node, reads the value, and returns to original archive
            position.
            @param arch
            @param val The object to serialize
            @param tagname The XML tagname to create
        */
        static void
        deserialize(
            const XMLArchivePtr& arch,
            T& val,
            const std::string& tagname
        )
        {
            AssertNotHere<> test;
        }

        /**
            Method for deserializing arbitrary types. This reads the value
            at the current XML node.
            @param arch
            @param val The object to serialize
        */
        static void
        deserialize(
            const XMLArchivePtr& arch,
            T& val
        )
        {
            AssertNotHere<> test;
        }

};

template <class T>
class SerialDecide<const T>
{
    public:
        /**
            Method for serializing arbitrary types. This creates
            the xml node, writes the value, and returns to original archive
            position.
            @param arch
            @param val The object to serialize
            @param tagname The XML tagname to create
        */
        static void
        serialize(
            const XMLArchivePtr& arch,
            const T& val,
            const std::string& tagname
        )
        {
            SerialDecide<T>::serialize(arch, val, tagname);
        }

        /**
            Method for serializing arbitrary types. This writes the value
            at the current XML node.
            @param arch
            @param val The object to serialize
        */
        static void
        serialize(
            const XMLArchivePtr& arch,
            const T& val
        )
        {
            SerialDecide<T>::serialize(arch, val);
        }

        /**
            Method for deserializing arbitrary types. This creates
            the xml node, reads the value, and returns to original archive
            position.
            @param arch
            @param val The object to serialize
            @param tagname The XML tagname to create
        */
        static void
        deserialize(
            const XMLArchivePtr& arch,
            const T& val,
            const std::string& tagname
        )
        {
            T& unconst = const_cast<T&>(val);
            SerialDecide<T>::deserialize(arch, val, tagname);
        }

        /**
            Method for deserializing arbitrary types. This reads the value
            at the current XML node.
            @param arch
            @param val The object to serialize
        */
        static void
        deserialize(
          const XMLArchivePtr& arch,
          const T& val
        )
        {
            T& unconst = const_cast<T&>(val);
            SerialDecide<T>::deserialize(arch, val);
        }
};

template <class T>
class SerialDecide< boost::intrusive_ptr<T> >
{
    public:
        static void
        serialize(
            const XMLArchivePtr& arch,
            const boost::intrusive_ptr<T>& val,
            const std::string& tagname
        )
        {
            arch->serialize_obj(val, tagname);
        }

        static void
        serialize(
            const XMLArchivePtr& arch,
            const boost::intrusive_ptr<T>& val
        )
        {
            arch->serialize_obj(val);
        }

        static void
        deserialize(
            const XMLArchivePtr& arch,
            boost::intrusive_ptr<T>& val,
            const std::string& tagname
        )
        {
            arch->deserialize_obj(val, tagname);
        }

        static void
        deserialize(
            const XMLArchivePtr& arch,
            boost::intrusive_ptr<T>& val
        )
        {
            arch->deserialize_obj(val);
        }

};

template <class T>
class SerialDecide< boost::intrusive_ptr<const T> >
{
    public:
        static void
        serialize(
            const XMLArchivePtr& arch,
            const boost::intrusive_ptr<const T>& val,
            const std::string& tagname
        )
        {
            boost::intrusive_ptr<T> obj(boost::const_pointer_cast<T>(val));
            SerialDecide<boost::intrusive_ptr<T> >::serialize(arch, obj, tagname);
        }

        static void
        serialize(
            const XMLArchivePtr& arch,
            const boost::intrusive_ptr<const T>& val
        )
        {
            boost::intrusive_ptr<T> obj(boost::const_pointer_cast<T>(val));
            SerialDecide<boost::intrusive_ptr<T> >::serialize(arch, obj);
        }

        static void
        deserialize(
            const XMLArchivePtr& arch,
            boost::intrusive_ptr<const T>& val,
            const std::string& tagname
        )
        {
            boost::intrusive_ptr<T> obj;
            SerialDecide<boost::intrusive_ptr<T> >::deserialize(arch, obj, tagname);
            val = obj;
        }

        static void
        deserialize(
            const XMLArchivePtr& arch,
            boost::intrusive_ptr<const T>& val
        )
        {
            boost::intrusive_ptr<T> obj;
            SerialDecide<boost::intrusive_ptr<T> >::deserialize(arch, obj);
            val = obj;
        }

};


#define vector_item_tagname "vector_item"
template <class T>
class SerialDecide< std::vector<T> >
{
    public:

        static void
        serialize(
            const XMLArchivePtr& arch,
            const std::vector<T>& vec,
            const std::string& tagname
        )
        {
            arch->appendAndStepIn(tagname);
            typename std::vector<T>::const_iterator it(vec.begin());
            for ( ; it != vec.end(); ++it)
            {
                serial_call_save(arch, *it, vector_item_tagname);
            }
            arch->stepOut();
        }

        static void
        deserialize(
            const XMLArchivePtr& arch,
            std::vector<T>& vec,
            const std::string& tagname
        )
        {
            arch->stepIn(tagname);
            arch->stepIn(vector_item_tagname);
            for ( ; arch->nonnull(); arch->nextSibling(vector_item_tagname))
            {
                T val;
                serial_call_load(arch, val);
                vec.push_back(val);
            }
            arch->stepOut(); //step out of vector_item
            arch->stepOut(); //step out of tagname
        }
};

#define map_item_tagname "map_item"
#define map_key_tagname "key"
#define map_value_tagname "value"

template <class T, class U>
class SerialDecide< std::map<T, U> >
{
    public:
        static void
        serialize(
            const XMLArchivePtr& arch,
            const std::map<T, U>& valmap,
            const std::string& tagname
        )
        {
            arch->appendAndStepIn(tagname);
            typename std::map<T,U>::const_iterator it(valmap.begin());
            for ( ; it != valmap.end(); ++it)
            {
                arch->appendAndStepIn(map_item_tagname);
                serial_call_save(arch, it->first, map_key_tagname);
                serial_call_save(arch, it->second, map_value_tagname);
                arch->stepOut();
            }
            arch->stepOut();
        }

        static void
        deserialize(
            const XMLArchivePtr& arch,
            std::map<T, U>& valmap,
            const std::string& tagname
        )
        {
            arch->stepIn(tagname);
            arch->stepIn(map_item_tagname);
            for ( ; arch->nonnull(); arch->nextSibling(map_item_tagname))
            {
                T key; U val;
                serial_call_load(arch, key, map_key_tagname);
                serial_call_load(arch, val, map_value_tagname);
                valmap[key] = val;
            }
            arch->stepOut();
            arch->stepOut();
        }
};

#define set_item_tagname "set_item"
template <class T>
class SerialDecide< Set<T> >
{
    public:
        static void
        serialize(
            const XMLArchivePtr& arch,
            const Set<T>& vec,
            const std::string& tagname
        )
        {
            arch->appendAndStepIn(tagname);
            typename Set<T>::const_iterator it(vec.begin());
            for ( ; it != vec.end(); ++it)
            {
                serial_call_save(arch, *it, set_item_tagname);
            }
            arch->stepOut();
        }

        static void
        deserialize(
            const XMLArchivePtr& arch,
            Set<T>& vec,
            const std::string& tagname
        )
        {
            arch->stepIn(tagname);
            arch->stepIn(set_item_tagname);
            for ( ; arch->nonnull(); arch->nextSibling(set_item_tagname))
            {
                T val;
                serial_call_load(arch, val);
                vec.append(val);
            }
            arch->stepOut(); //step out of set_item
            arch->stepOut(); //step out of tagname
        }
};

#define SerialDecideInstance(x) \
template <> class SerialDecide<x> { \
    public: \
        static void  \
        serialize( \
            const XMLArchivePtr& arch, \
            const x & val, \
            const std::string& tagname \
        ) \
        { \
            arch->setValue(val, tagname); \
        } \
        static void  \
        serialize( \
            const XMLArchivePtr& arch, \
            const x & val \
        ) \
        { \
            arch->setValue(val); \
        } \
        static void  \
        deserialize( \
            const XMLArchivePtr& arch, \
            x& val, \
            const std::string& tagname \
        ) \
        { \
            arch->getValue(val, tagname); \
        } \
        static void  \
        deserialize( \
            const XMLArchivePtr& arch, \
            x & val \
        ) \
        { \
            arch->getValue(val); \
        } \
}

SerialDecideInstance(int);
SerialDecideInstance(double);
SerialDecideInstance(unsigned int);
SerialDecideInstance(unsigned long);
SerialDecideInstance(bool);
SerialDecideInstance(std::string);

#define SerialDecideSubptr(x) \
template <> class SerialDecide<x> { \
    public: \
        typedef x::element_type subtype; \
        static void  \
        serialize( \
            const XMLArchivePtr& arch, \
            const x & val, \
            const std::string& tagname \
        ) \
        { \
            boost::intrusive_ptr<subtype> tmp = val.get(); \
            arch->serialize_obj<subtype>(tmp, tagname); \
        } \
        static void  \
        deserialize( \
            const XMLArchivePtr& arch, \
            x& val, \
            const std::string& tagname \
        ) \
        { \
            boost::intrusive_ptr<subtype> tmp = val.get(); \
            arch->deserialize_obj<subtype>(tmp, tagname); \
            val = tmp.get(); \
        } \
}

template <class T> class NonConst {};
template <class T> class NonConst<const T> {
    public:
        typedef T element_type;
};

#define SerialDecideConstSubptr(x) \
template <> class SerialDecide<x> { \
    public: \
        typedef NonConst<x::element_type>::element_type subtype; \
        static void  \
        serialize( \
            const XMLArchivePtr& arch, \
            const x & val, \
            const std::string& tagname \
        ) \
        { \
            boost::intrusive_ptr<subtype> tmp = const_cast<subtype*>(val.get()); \
            arch->serialize_obj<subtype>(tmp, tagname); \
        } \
        static void  \
        deserialize( \
            const XMLArchivePtr& arch, \
            x& val, \
            const std::string& tagname \
        ) \
        { \
            boost::intrusive_ptr<subtype> tmp = const_cast<subtype*>(val.get()); \
            arch->deserialize_obj<subtype>(tmp, tagname); \
            val = tmp.get(); \
        } \
}


template <class T>
void
serial_call_save(
    const XMLArchivePtr& arch,
    const T& val,
    const std::string& tagname
)
{
    SerialDecide<T>::serialize(arch, val, tagname);
}

template <class T>
void
serial_call_save(
    const XMLArchivePtr& arch,
    const T& val
)
{
    SerialDecide<T>::serialize(arch, val);
}


template <class T>
void
serial_call_load(
    const XMLArchivePtr& arch,
    T& val,
    const std::string& tagname
)
{
   SerialDecide<T>::deserialize(arch, val, tagname);
}

template <class T>
void
serial_call_load(
    const XMLArchivePtr& arch,
    T& val
)
{
   SerialDecide<T>::deserialize(arch, val);
}

/**
    @class SerialRuntime
    Encapsulates dynamic type info for polymorphic types that derive
    from Serializable.
*/
class SerialRuntime {

    public:
        /**
            Function instance that creates a polymorphic type form the
            archive deserialization constructor.
        */
        typedef Serializable*(create_function)(const XMLArchivePtr&);

        /**
            Constructor for serial runtime info.
            @param str The classname to associate with a given type
            @param fxnptr A function instance that creates the polymorphic type
        */
        SerialRuntime(const std::string& str, create_function* fxnptr);

        /**
            Function for creating a polymorphic type for an archive.  This queries
            the registry of class names to find the appropriate #create_function.
            @param arch The xml archive node to deserialize from
            @param classname The string identifying the class
        */
        static SerializablePtr getObject(
            const XMLArchivePtr& arch,
            const std::string& classname
        );

        /**
        */
        static void free();

        /**
            Method associated with a given type containing the class name.
            @return A string identifying the polymorphic type
        */
        virtual std::string classname() const = 0;

    private:
        /**
            The map registering class names with create functions.
        */
        static std::map<std::string, create_function*> *classlist_;



};

/**
    Template function for generating the correct polymoprhic type.
    This is used by SerialRuntime to generate the appropriate
    SerialRuntime::create_function.
    @param arch The archive to deserialize from
*/
template <class T>
Serializable*
build(const XMLArchivePtr& arch)
{
    return new T(arch);
}

/**
    @class SerialClass
    A template class for generating runtime info specific to a given type.
*/
template <class T> 
class SerialClass : public SerialRuntime {

    private:
        /**
            A string identifier for the class
        */
        static std::string classname_;

        /**
            The function for deserializing the class from an XML archive
        */
        static create_function* fxn_;

    public:

        /**
            Constructor
            @param name The name to associate with the class type
        */
        SerialClass(const std::string& name) :
            SerialRuntime(name, fxn_)
        {
            //this must be set here because it is a static variable
            classname_= name;
        }

        /**
            @return The name associated with a class type
        */
        std::string
        classname() const
        {
            return classname_;
        }

};

/**
    @class SerialMap
*/
class SerialMap : public Countable {

    private:
        typedef std::map<unsigned long, boost::intrusive_ptr<const Serializable> > registry_map;

        registry_map registry_;

    public:
        SerialMap();

        boost::intrusive_ptr<Serializable> get(unsigned long id) const;

        bool has(unsigned long id) const;

        void add(unsigned long id, const ConstSerializablePtr& ptr);

};

template<class Y>
std::ostream & operator<< (std::ostream & os, boost::intrusive_ptr<Y> const & p)
{
    if (p)
        p->print(os);
    else
        os << "null object" << std::endl;
    return os;
}

std::ostream& operator<< (std::ostream& os, Serializable* s);


} //end namespace

#endif
