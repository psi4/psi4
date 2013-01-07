#ifndef smartptr_xmlarchive_h
#define smartptr_xmlarchive_h

#include "ref.h"
#include "rapidxml.hpp"
#include "serialabstract.h"

#include <list>
#include <map>
#include <vector>

#define SerialMapPtr boost::intrusive_ptr<smartptr::SerialMap>
#define XMLArchivePtr boost::intrusive_ptr<smartptr::XMLArchive>
#define BinaryMapPtr boost::intrusive_ptr<smartptr::BinaryMap>  

#define NULL_ID 0

typedef rapidxml::xml_node<> xmlnode;
typedef rapidxml::xml_document<> xmldocument;
typedef rapidxml::xml_base<> xmlbase;

namespace smartptr {

/**
    @struct binary_node_t
    Struct that encapsulates info about a block of binary data in an
    archive.
*/
struct binary_node_t {

    public:
        void* vals;
        size_t offset;
        size_t size;

};

/**
    @class BinaryMap
    Class for registering blocks of binary data in an archive.
*/
class BinaryMap : public smartptr::Countable {

    private:
        std::map<void*, binary_node_t*> ptrmap_;

        std::map<size_t, binary_node_t*> offsetmap_;

        std::vector<binary_node_t*> nodes_;

    public:
        typedef std::vector<binary_node_t*>::iterator iterator;

        /**
            Gets information about binary data in the archive based on
            the malloc'd pointer.
            @param ptr
            @return The binary node info struct This is 0
                    if pointer is not found.
        */
        binary_node_t* get(void* ptr);

        /**
            Gets information about binary data in the archive based on
            the offset in the binary archive.
            @param offset
            @return The binary node info struct. This is 0
                    if pointer is not found.
        */
        binary_node_t* get(size_t offset);

        iterator begin();

        iterator end();

        /**
            Register a block of binary data in the archive
            @param The malloc'd pointer in memory
            @param offset The offset in the binary archive
            @param The size of the data block
        */
        void
        insert(
            void* ptr,
            size_t offset,
            size_t size
        );

};

/**
    @class XMLArchive
    Class for archiving objects in an xml file
*/
class XMLArchive : public smartptr::Countable {

    public:

        typedef enum { Checkpoint, Runtime } archive_mode_t;

        typedef std::pair<SerializablePtr, xmlnode*> object_t;

        typedef std::vector<object_t>::iterator object_iter_t;

    private:
        static std::string getIDTag(unsigned long id);

        std::list<xmlnode*> nodes_;

        xmldocument* xmldoc_;

        xmlnode* node_;
        
        xmlnode* objnode_;

        size_t binary_size_;

        BinaryMapPtr binary_map_;

        std::vector<object_t> object_nodes_;

        char* buffer_;

        size_t data_size_;

        archive_mode_t mode_;

        std::istream* binstream_;

        object_iter_t
        _getObjectNode(size_t idx);

        template <class T>
        void
        _getValue(T& val, xmlbase* node);

        void
        _getValue(bool& val, xmlbase* node);

        void
        _getValue(std::string& val, xmlbase* node);

        template <class T>
        void
        _getValue(T& val, const std::string& name);

        template <class T>
        void
        _getAttribute(T& val, const std::string& attrname);

        template <class T>
        void
        _setValue(T val, xmlbase* node);

        template <class T>
        void
        _setValue(T val, const std::string& name);

        void
        _setValue(const std::string& name, xmlbase* node);

        template <class T>
        void
        _setAttribute(T val, const std::string& name);

        void
        _stepIn(xmlnode* node, const std::string& tagname);

        size_t
        _countSiblings(xmlnode* node);

        void _appendAndStepIn(xmlnode* node, const std::string& tagname);

        SerializablePtr getObject(unsigned long id);

        void setObject(const SerializablePtr& obj);

        void _setValue(void* val, xmlbase* node);

        void _setValue(int val, xmlbase* node);

        void _setValue(long val, xmlbase* node);

        void _setValue(unsigned int val, xmlbase* node);

        void _setValue(unsigned long val, xmlbase* node);

        void _setValue(double val, xmlbase* node);

        void _setValue(bool val, xmlbase* node);

        void _setBinary(void* vals, size_t size, const std::string& name);

        void _getBinary(void*& vals, size_t& size, const std::string& name);

    public:

        void appendAndStepIn(const std::string& tagname);

        void stepIn(const std::string& tagname);

        void stepTo(
            const std::string& tagname,
            const std::string& attrname,
            const std::string& attrnval
        );

        void stepOut();

        void nextSibling(const std::string& tagname);

        XMLArchive(
            const std::string& filename,
            archive_mode_t mode
        );

        XMLArchive(
            archive_mode_t mode
        );

        ~XMLArchive();

        void toFile(const std::string& filename);

        std::string toPrettyXML();

        unsigned long getNextID();

        size_t getDataSize() const;

        bool hasAttribute(const std::string& attrname);
        
        void getAttribute(int& val, const std::string& attrname);

        void getAttribute(long& val, const std::string& attrname);

        void getAttribute(size_t& val, const std::string& attrname);

        void getAttribute(double& val, const std::string& attrname);

        void getAttribute(unsigned int& val, const std::string& attrname);

        void getAttribute(bool& val, const std::string& attrname);

        void getAttribute(std::string& val, const std::string& attrname);

        void setAttribute(int val, const std::string& attrname);

        void setAttribute(long val, const std::string& attrname);

        void setAttribute(size_t val, const std::string& attrname);

        void setAttribute(double val, const std::string& attrname);

        void setAttribute(unsigned int val, const std::string& attrname);

        void setAttribute(bool val, const std::string& attrname);

        void setAttribute(const std::string& val, const std::string& attrname);

        void getValue(void*& val);

        void getValue(int& val);

        void getValue(long& val);

        void getValue(unsigned int& val);

        void getValue(unsigned long& val);

        void getValue(double& val);

        void getValue(bool& val);

        void getValue(std::string& value);

        void setValue(void* val);

        void setValue(int val);

        void setValue(long val);

        void setValue(unsigned int val);

        void setValue(unsigned long val);

        void setValue(double val);

        void setValue(bool val);

        void setValue(const std::string& value);

        void getValue(void*& val, const std::string& tagname);

        void getValue(int& val, const std::string& tagname);

        void getValue(long& val, const std::string& tagname);

        void getValue(unsigned int& val, const std::string& tagname);

        void getValue(unsigned long& val, const std::string& tagname);

        void getValue(double& val, const std::string& tagname);

        void getValue(bool& val, const std::string& tagname);

        void getValue(std::string& value, const std::string& tagname);

        void setValue(void* val, const std::string& tagname);

        void setValue(int val, const std::string& name);

        void setValue(long val, const std::string& name);

        void setValue(unsigned int val, const std::string& tagname);

        void setValue(unsigned long val, const std::string& tagname);

        void setValue(double val, const std::string& name);

        void setValue(bool val, const std::string& tagname);

        void setValue(const std::string& val, const std::string& name);


        template <class T>
        void setBinary(const T* vals, size_t size, const std::string& name)
        {
            void* ptr = reinterpret_cast<void*>(const_cast<T*>(vals));
            _setBinary(ptr, size * sizeof(double), name);
        }

        template <class T>
        void getBinary(T*& vals, size_t& size, const std::string& name)
        {
            void* ptr;
            _getBinary(ptr, size, name);
            size /= sizeof(T); //divide the size
            vals = reinterpret_cast<T*>(ptr);
        }

        template <class T>
        void getBinary(const T*& vals, size_t& size, const std::string& name)
        {
            T* ptr;
            getBinary<T>(ptr, size, name);
            vals = ptr;
        }

        /**
            Figure out whether the current xml node is valid.
            @return If the current node does not exist (is null)
        */
        bool null();

        /**
            Figure out whether the current xml node is valid.
            @return If the current node exists (is not null)
        */
        bool nonnull();

        void registerObject(const SerializablePtr& obj);

        void configureSerialize(
            SerializablePtr obj,
            unsigned long& id,
            std::string& type
        );

        template <class T>
        void
        deserialize_obj(boost::intrusive_ptr<T>& obj)
        {
            unsigned long id;
            getValue(id);
            if (id == NULL_ID)
            {
                obj = 0;
            }
            else
            {
                SerializablePtr tmp(getObject(id));
                obj = boost::static_pointer_cast<T, Serializable>(tmp);
                if (!obj)
                {
                    std::cerr << "failed to cast object in serialize" << std::endl;
                    abort();
                }
            }
        }

        template <class T>
        void
        deserialize_obj(boost::intrusive_ptr<T>& obj, const std::string& name)
        {
            stepIn(name);
            if (null())
            {
                std::cerr << "no node " << name << std::endl;
                std::cerr << toPrettyXML() << std::endl;
                abort();
            }
            deserialize_obj<T>(obj);
            stepOut();
        }

        template <class T>
        void
        serialize_obj(const boost::intrusive_ptr<T>& obj)
        {
            if (!obj)
            {
                setValue(NULL_ID);
                return;
            }

            unsigned long id;
            std::string idtype;
            configureSerialize(obj, id, idtype);
            setAttribute(idtype, "type");


            setValue(id);
            setObject(obj);
        }

        template <class T>
        void
        serialize_obj(const boost::intrusive_ptr<T>& obj, const std::string& name)
        {
            appendAndStepIn(name);
            serialize_obj<T>(obj);
            stepOut();
        }

        template <class T>
        void
        enum_cast(unsigned int tmp, T& val)
        {
            val = (T) tmp;
        }

        void
        profile(int n);



};

}

#endif

