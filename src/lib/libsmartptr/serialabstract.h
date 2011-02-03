#ifndef smartptr_serialabstract_h_
#define smartptr_serialabstract_h_

#include <string>
#include <iostream>
#include <map>

#include "ref.h"
#include "set.h"

#define SerializablePtr boost::intrusive_ptr<smartptr::Serializable>  
#define ConstSerializablePtr boost::intrusive_ptr<const smartptr::Serializable>  


#define serid_format_str "%ld"

namespace smartptr {

template <bool t = true>
struct AssertNotHere {
  enum { N = 1 - 2 * int(!t) };
      // 1 if t is true, -1 if t is false.
  static char A[N];
};

#define XMLArchivePtr boost::intrusive_ptr<smartptr::XMLArchive>
class XMLArchive;
class SerialRuntime;

/**
 @class Serializable is the base class for all classes that can be serialized to an
 XMLArchive.  Derived classes must implement two things.  First, a serialization
 constructor must exist that takes an XMLArchive as an argument. Second, the serialize method
 must be reimplemented to serialize data members.  Derived serialize methods should
 call the parent serialize method.
*/
class Serializable : public Countable {
    
    protected:
        /** An object encapsulating runtime info about the polymorphic class type */
        SerialRuntime* rtinfo_;

    private:
        /** 0 upon construction.
            The archive id for the object. When serializing objects to an archive,
            pointer types must have a unique id assigned to ensure that the data of
            pointer types only get written once to the archive.  The archive id is
            used for deserializing objects from the archive.  */
        unsigned long archive_id_;

        /** 0 upon construction.
            The runtime id for the object.  When passing objects between processes,
            objects can be declared "global" or "replicated" by assigning a runtime id.
            No data for replicated objects will be written to an archive - only the runtime
            id. When deserializing on the remote process, the runtime id is matched to the
            local pointer rather than explicitly deserialized from the archive. */
        unsigned long runtime_id_;

    public:
        /**
            @return The archive id.  See #archive_id_;
        */
        unsigned long getArchiveID() const;

        /**
            @return The runtime id.  See #runtime_id_;
        */
        unsigned long getRuntimeID() const;

        /**
            Set the archive id for the object. Except in very rare cases,
            this should only ever be called by an archive object. If the id
            has already been set or is set to 0, routine aborts.
            @param id The new archive id. Cannot be set to 0.
        */
        void setArchiveID(unsigned long id);

        /**
            Set the runtime id for the object. The user must devise a unique
            indexing scheme for replicated objects and set this explicitly.
            If the id has already been set or is set to 0, routine aborts.
            @param id The new archive id. Cannot be set to 0.
        */
        void setRuntimeID(unsigned long id);

        /**
            Reset the archive id to 0.  Except in very rare cases, this should only
            ever be called by an archive object destructor.
        */
        void resetArchiveID() const;

        /**
            Serialize the object to the xml archive.
            @param archive The XML archive to write to
        */
        virtual void serialize(const XMLArchivePtr& archive) const;

        /**
            @return The runtime info containing polymorphic type information.
        */
        const SerialRuntime* runtime_info() const;

        virtual void print(std::ostream& os = std::cout) const;

    protected:
        Serializable();

        /**
            Serialization constructor
            @param archive The xml archive to deserialize from
        */
        Serializable(const XMLArchivePtr& archive);


        virtual ~Serializable();
};

/**
    Implicitly-instantiated template function for deserializing objects
    from an archive.  This is usually called through a macro rather
    than directly.  This steps into the xml node, reads the value, and
    steps back to original archive position.
    @param arch The archive to deserialize form
    @param val The object to deserialize
    @param tagname The node tagname in the xml
*/
template <class T>
void
serial_call_load(
    const XMLArchivePtr& arch,
    T& val,
    const std::string& tagname
);

/**
    Implicitly-instantiated template function for deserializing objects
    from an archive.  This is usually called through a macro rather
    than directly.  This writes directly to the current xml node.
    @param arch The archive to deserialize from
    @param val The object to deserialize
*/
template <class T>
void
serial_call_load(
    const XMLArchivePtr& arch,
    T& val
);

/**
    Implicitly-instantiated template function for serializing objects
    to an archive.  This is usually called through a macro rather
    than directly.  This creates the xml node, writes the value, and
    steps back to original archive position.
    @param arch The archive to serialize to
    @param val The object to serialize
    @param tagname The node tagname in the xml
*/
template <class T>
void
serial_call_save(
    const XMLArchivePtr& arch,
    const T& val,
    const std::string& tagname
);

/**
    Implicitly-instantiated template function for serializing objects
    to an archive.  This is usually called through a macro rather
    than directly.  This writes directly to the current xml node.
    @param arch The archive to serialize to
    @param val The object to serialize
*/
template <class T>
void
serial_call_save(
    const XMLArchivePtr& arch,
    const T& val
);



} //end namespace

#undef heisenbug

#endif
