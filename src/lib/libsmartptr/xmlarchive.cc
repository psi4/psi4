#include <fstream>
#include "rapidxml_print.hpp"

#include "xmlarchive.h"
#include "serialimpl.h"
#include "printstream.h"
#include "timer.h"

using namespace smartptr;
using namespace rapidxml;
using namespace std;

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

XMLArchive::XMLArchive(
    const std::string& filename,
    archive_mode_t mode
)
    : xmldoc_(0), node_(0), objnode_(0),
        object_nodes_(0),
        buffer_(0),
        data_size_(0), mode_(mode),
        binstream_(0),
        binary_map_(new BinaryMap),
        binary_size_(0)
{
    ifstream file(filename.data(), ios::in | ios::ate);

    if (!file.is_open())
    {
        cerr << filename << " does not exist for xml parsing" << endl;
        abort();
    }

    size_t size = file.tellg();
    buffer_ = new char[size + 1];
    file.seekg (0, ios::beg);
    file.read (buffer_, size);
    file.close();
    *(buffer_ + size) = '\0';

    std::string binfilename = filename + ".bin";
    ifstream* fbinstream = new ifstream(binfilename.data(), ios::in | ios::binary);
    if (fbinstream->is_open())
    {
        binstream_ = fbinstream;
    }



    xmldoc_ = new xml_document<>;
    xmldoc_->parse<0>(buffer_);
    node_ = xmldoc_->first_node();
    objnode_ = node_->first_node("objects");

    if (objnode_)
    {
        xmlnode* node = objnode_->first_node();
        for  ( ; node; node = node->next_sibling())
        {
            object_t objnode(0, node);
            object_nodes_.push_back(objnode);
        }
    }
}

XMLArchive::XMLArchive(
    archive_mode_t mode
)
    : xmldoc_(0),
      node_(0),
      objnode_(0),
      buffer_(0),
      mode_(mode),
      data_size_(0),
      binstream_(0),
      binary_size_(0),
      binary_map_(new BinaryMap)
{
    xmldoc_ = new xml_document<>;
    node_ = xmldoc_->allocate_node(node_element, "document");
    xmldoc_->append_node(node_);

    objnode_ = xmldoc_->allocate_node(node_element, "objects");
    node_->append_node(objnode_);
}

XMLArchive::~XMLArchive()
{
    if (buffer_)
        delete[] buffer_;

    //clear the archive id's for all the objects
    object_iter_t it(object_nodes_.begin());
    for ( ; it != object_nodes_.end(); ++it)
    {
        SerializablePtr obj(it->first);
        if (obj)
            obj->resetArchiveID();
    }

    delete binstream_;

    delete xmldoc_;
}

void
XMLArchive::nextSibling(const std::string& name)
{
    node_ = node_->next_sibling(name.data());
}

size_t
XMLArchive::_countSiblings(xmlnode* node)
{
    xmlnode* child = node;
    size_t n = 0;
    while (child)
    {
        child = child->next_sibling();
        ++n;
    }
    return n;
}

XMLArchive::object_iter_t
XMLArchive::_getObjectNode(size_t idx)
{
    return ( object_nodes_.begin() + (idx - 1) );
    //return objnodes_[idx - 1];
}

void
XMLArchive::_stepIn(xmlnode* node, const std::string& name)
{
    xml_node<>* newnode = node->first_node(name.data());
    nodes_.push_back(node_);
    node_ = newnode;
}

void
XMLArchive::stepIn(const std::string& name)
{
    _stepIn(node_, name);
}

void
XMLArchive::stepTo(
    const std::string& tagname,
    const std::string& attrname,
    const std::string& attrval
)
{
    stepIn(tagname);

    bool match = false;
    std::string check;
    while (!match && node_)
    {
        if (hasAttribute(attrname))
        {
            getAttribute(check, attrname);
            match = (check == attrval);
        }

        if (!match)
            nextSibling(tagname);
    }
}

void
XMLArchive::stepOut()
{
    node_ = nodes_.back();
    nodes_.pop_back();
}

void
XMLArchive::_appendAndStepIn(xmlnode* node, const std::string& name)
{
    const char* namestr = node->document()->allocate_string(name.data());
    xmlnode* newnode = node->document()->allocate_node(node_element, namestr);
    node->append_node(newnode);

    nodes_.push_back(node_);
    node_ = newnode;
}

void
XMLArchive::appendAndStepIn(const std::string& name)
{
    _appendAndStepIn(node_, name);
}

bool
XMLArchive::hasAttribute(const std::string& attrname)
{
    xml_attribute<>* attr = node_->first_attribute(attrname.data());
    return attr;
}

template <class T>
void
XMLArchive::_getValue(T& val, xmlbase* node)
{
    if (!node)
    {
        cerr << "getValue called on null node" << endl;
        abort();
    }

    string s(node->value());
    if (!s.size())
    {
        cerr << "no data on node" << endl;
        abort();
    }
    stringstream sstr(s);
    sstr >> val;
}

void
XMLArchive::_getValue(std::string& val, xmlbase* node)
{
    if (!node)
    {
        cerr << "getValue called on null node" << endl;
        abort();
    }

    val = node->value();
}

template <class T>
void
XMLArchive::_getValue(T& val, const std::string& name)
{
    stepIn(name);
    if (!node_)
    {
        stepOut();
        cerr << "no node " << name <<  " found" << endl;
        cerr << toPrettyXML() << endl;
        abort();
    }
    getValue(val);
    stepOut();
}

void
XMLArchive::_getValue(bool& val, xmlbase* node)
{
    std::string text;
    _getValue(text, node);
    if (text == "true")
        val = true;
    else
        val = false;
}

template <class T>
void
XMLArchive::_setValue(T val, const std::string& name)
{
    appendAndStepIn(name);
    _setValue(val, node_);
    stepOut();
}


template <class T>
void
XMLArchive::_setAttribute(T val, const std::string& attrname)
{
    const char* namestr = xmldoc_->allocate_string(attrname.data());
    xml_attribute<>* attr = xmldoc_->allocate_attribute(namestr);
    node_->append_attribute(attr);
    _setValue(val, attr);
}

template <class T>
void
XMLArchive::_getAttribute(T& val, const std::string& attrname)
{
    xml_attribute<>* attr = node_->first_attribute(attrname.data());
    _getValue(val, attr);
}

void
XMLArchive::_setValue(const std::string& val, xmlbase* node)
{
    const char* valstr = xmldoc_->allocate_string(val.data());
    node->value(valstr);
}

void
XMLArchive::_setValue(double val, xmlbase* node)
{
    data_size_ += sizeof(double);
    _setValue(stream_printf("%20.14f", val), node);
}

void
XMLArchive::_setValue(long val, xmlbase* node)
{
    data_size_ += sizeof(long);
    _setValue(stream_printf("%ld", val), node);
}

void
XMLArchive::_setValue(unsigned int val, xmlbase* node)
{
    data_size_ += sizeof(unsigned int);
    _setValue(stream_printf("%d", val), node);
}

void
XMLArchive::_setValue(unsigned long val, xmlbase* node)
{
    data_size_ += sizeof(unsigned long);
    _setValue(stream_printf("%ld", val), node);
}

void
XMLArchive::_setValue(void* val, xmlbase* node)
{
    data_size_ += sizeof(void*);
    _setValue(stream_printf("%p", val), node);
}

void
XMLArchive::_setValue(bool val, xmlbase* node)
{
    data_size_ += sizeof(int);
    if (val)
    {
        std::string truestr("true");
        _setValue(truestr, node);
    }
    else
    {
        std::string falsestr("false");
        _setValue(falsestr, node);
    }
}

void
XMLArchive::_setValue(int val, xmlbase* node)
{
    data_size_ += sizeof(int);
    _setValue(stream_printf("%d", val), node);
}

void
XMLArchive::getValue(std::string& val)
{
    val = node_->value();
}

void
XMLArchive::setValue(const std::string& val)
{
    data_size_ += val.size();
    _setValue(val, node_);
}

void
XMLArchive::setValue(double val)
{
    _setValue(val, node_);
}

void
XMLArchive::setValue(long val)
{
    _setValue(val, node_);
}

void
XMLArchive::setValue(unsigned int val)
{
    _setValue(val, node_);
}

void
XMLArchive::setValue(unsigned long val)
{
    _setValue(val, node_);
}

void
XMLArchive::setValue(void* val)
{
    _setValue(val, node_);
}

void
XMLArchive::setValue(bool val)
{
    _setValue(val, node_);
}

void
XMLArchive::setValue(int val)
{
    _setValue(val, node_);
}

void 
XMLArchive::getAttribute(int& val, const std::string& attrname)
{
    _getAttribute(val, attrname);
}

void 
XMLArchive::getAttribute(long& val, const std::string& attrname)
{
    _getAttribute(val, attrname);
}

void 
XMLArchive::getAttribute(unsigned long& val, const std::string& attrname)
{
    _getAttribute(val, attrname);
}

void 
XMLArchive::getAttribute(double& val, const std::string& attrname)
{
    _getAttribute(val, attrname);
}

void 
XMLArchive::getAttribute(unsigned int& val, const std::string& attrname)
{
    _getAttribute(val, attrname);
}

void 
XMLArchive::getAttribute(bool& val, const std::string& attrname)
{
    _getAttribute(val, attrname);
}

void 
XMLArchive::getAttribute(std::string& val, const std::string& attrname)
{
    _getAttribute(val, attrname);
}

void 
XMLArchive::setAttribute(int val, const std::string& attrname)
{
    _setAttribute(val, attrname);
}

void 
XMLArchive::setAttribute(long val, const std::string& attrname)
{
    _setAttribute(val, attrname);
}

void 
XMLArchive::setAttribute(unsigned long val, const std::string& attrname)
{
    _setAttribute(val, attrname);
}

void 
XMLArchive::setAttribute(double val, const std::string& attrname)
{
    _setAttribute(val, attrname);
}

void 
XMLArchive::setAttribute(unsigned int val, const std::string& attrname)
{
    _setAttribute(val, attrname);
}

void 
XMLArchive::setAttribute(bool val, const std::string& attrname)
{
    _setAttribute(val, attrname);
}

void 
XMLArchive::setAttribute(const std::string& val, const std::string& attrname)
{
    data_size_ += val.size();
    _setAttribute(val, attrname);
}

void
XMLArchive::getValue(void*& val)
{
    _getValue<void*>(val, node_);
}

void
XMLArchive::getValue(int& val)
{
    _getValue<int>(val,node_);
}

void
XMLArchive::getValue(long& val)
{
    _getValue<long>(val, node_);
}

void
XMLArchive::getValue(unsigned int& val)
{
    _getValue<unsigned int>(val, node_);
}

void
XMLArchive::getValue(unsigned long& val)
{
    _getValue<unsigned long>(val, node_);
}

void
XMLArchive::getValue(double& val)
{
    _getValue<double>(val, node_);
}

void
XMLArchive::getValue(bool& val)
{
    _getValue(val, node_);
}

void 
XMLArchive::getValue(int& val, const std::string& name)
{
    _getValue<int>(val, name);
}

void 
XMLArchive::getValue(long& val, const std::string& name)
{
    _getValue<long>(val, name);
}

void 
XMLArchive::getValue(unsigned int& val, const std::string& name)
{
    _getValue<unsigned int>(val, name);
}

void 
XMLArchive::getValue(unsigned long& val, const std::string& name)
{
    _getValue<unsigned long>(val, name);
}

void 
XMLArchive::getValue(double& val, const std::string& name)
{
    _getValue<double>(val, name);
}

void 
XMLArchive::getValue(bool& val, const std::string& name)
{
    _getValue<bool>(val, name);
}

void 
XMLArchive::getValue(std::string& val, const std::string& name)
{
    _getValue<std::string>(val, name);
}

void
XMLArchive::getValue(void*& val, const std::string& name)
{
    _getValue<void*>(val, name);
}

void 
XMLArchive::setValue(int val, const std::string& name)
{
    _setValue<int>(val, name);
}

void 
XMLArchive::setValue(long val, const std::string& name)
{
    _setValue<long>(val, name);
}

void 
XMLArchive::setValue(unsigned int val, const std::string& name)
{
    _setValue<double>(val, name);
}

void 
XMLArchive::setValue(unsigned long val, const std::string& name)
{
    _setValue<unsigned long>(val, name);
}

void 
XMLArchive::setValue(double val, const std::string& name)
{
    _setValue<double>(val, name);
}

void 
XMLArchive::setValue(bool val, const std::string& name)
{
    _setValue<bool>(val, name);
}

void 
XMLArchive::setValue(const std::string& val, const std::string& name)
{
    data_size_ += val.size();
    _setValue<const std::string&>(val, name);
}

void
XMLArchive::setValue(void* val, const std::string& name)
{
    _setValue<void*>(val, name);
}

void 
XMLArchive::_setBinary(void* vals, size_t size, const std::string& name)
{
    data_size_ += size;

    appendAndStepIn(name);

    binary_node_t* node = binary_map_->get(vals);
    size_t offset;
    if (node == 0) //doesn't yet exist
    {
        offset = binary_size_;
        binary_size_ += size;
        binary_map_->insert(vals, offset, size);
    }

    setValue(size, "size");
    setValue(offset, "offset");

    stepOut();
}

void 
XMLArchive::_getBinary(void*& vals, size_t& size, const std::string& name)
{
    if (!binstream_)
    {
        cerr << "no binary data exists in archive" << endl;
        abort();
    }

    stepIn(name);

    size_t offset;
    getValue(offset, "offset");
    //look up the offset in the map to see if we have already read this in
    binary_node_t* node = binary_map_->get(offset);

    if (node)
    {
        vals = node->vals;
        size = node->size;
        stepOut();
        return;
    }

    //new pointer
    getValue(size, "size");
    char* tmp = new char[size];
    binstream_->seekg(offset);
    binstream_->read(tmp, size);
    vals = reinterpret_cast<void*>(tmp);

    binary_map_->insert(tmp, offset, size);

    stepOut();
}

void
XMLArchive::toFile(const std::string& filename)
{
    ofstream file(filename.data());
    std::string s;
    print(std::back_inserter(s), *xmldoc_, 1);
    file << s;
    file.close();

    //now commit the binary
    std::string binfilename = filename + ".bin";
    ofstream binfile(binfilename.data(), ios::binary);

    BinaryMap::iterator it(binary_map_->begin());
    for ( ; it != binary_map_->end(); ++it)
    {
        binary_node_t* node = *it;
        //map value, pair 2nd
        size_t size = node->size;
        const char* ptr = reinterpret_cast<const char*>(node->vals);
        binfile.write(ptr, size);
    }
    binfile.close();
}

SerializablePtr
XMLArchive::getObject(unsigned long id)
{
    object_iter_t iter(_getObjectNode(id));
    SerializablePtr obj(iter->first);
    if (obj)
        return obj;


    xmlnode* node = iter->second;
    nodes_.push_back(node_);
    node_ = node;

    if (null())
    {
        cerr << id << " is null object" << endl;
        abort();
    }

    std::string classname; getAttribute(classname, "classname");

    obj = SerialRuntime::getObject(this, classname);

    stepOut();

    return obj;
}

void
XMLArchive::setObject(const SerializablePtr& obj)
{
    unsigned long id = obj->getArchiveID();
    if (id <= object_nodes_.size()) //this is already in archive
        return;

    std::string name = getIDTag(id);
    _appendAndStepIn(objnode_, name);

    object_t node(obj, node_);
    object_nodes_.push_back(node);

    const SerialRuntime* rtinfo = obj->runtime_info();
    setAttribute(rtinfo->classname(), "classname");
    obj->serialize(this);
    stepOut();
}

unsigned long
XMLArchive::getNextID()
{
    return object_nodes_.size() + 1;
}

std::string
XMLArchive::getIDTag(unsigned long id)
{
    return stream_printf("%ld", id);
}

bool
XMLArchive::null()
{
    bool check = node_;
    return !check;
}

bool
XMLArchive::nonnull()
{
    return node_;
}

std::string
XMLArchive::toPrettyXML()
{
    if (!node_)
        return "null";

    std::string s;
    print(std::back_inserter(s), *node_, 1);
    return s;
}

void
XMLArchive::profile(int n)
{
    tstart("next sibling");
    for (int i=0; i < n; ++i)
    {
        objnode_->next_sibling();
    }
    tstop("next sibling");

    tstart("compare");
    const char* p1 = "hello there";
    const char* p2 = "hello world";
    size_t size = 11;
    for (int i=0; i < n; ++i)
    {
        internal::compare(p1, size, p2, size, true);
    }
    tstop("compare");
}

size_t
XMLArchive::getDataSize() const
{
    return data_size_;
}

void
XMLArchive::registerObject(const SerializablePtr& obj)
{
    unsigned long id = obj->getArchiveID();
    if (id > object_nodes_.size())
    {
        cerr << "id " << id << " does not exist in archive" << endl;
        abort();
    }

    object_iter_t iter(_getObjectNode(id));
    iter->first = obj;
}

void
XMLArchive::configureSerialize(
    SerializablePtr obj,
    unsigned long &id,
    std::string &type
)
{
    if (mode_ == Runtime)
    {
        //check for a runtime id
        id = obj->getRuntimeID();
        if (id) //valid runtime object
        {
            type = "runtime";
            return;
        }
    }

    //simple archive object
    id = obj->getArchiveID();
    if (!id)
    {
        id = getNextID();
        obj->setArchiveID(id);
    }
    type = "archive";
}

binary_node_t*
BinaryMap::get(void *ptr)
{
    std::map<void*, binary_node_t*>::const_iterator it(ptrmap_.find(ptr));
    if (it == ptrmap_.end())
        return 0;

    return it->second;
}

binary_node_t*
BinaryMap::get(size_t offset)
{
    std::map<size_t, binary_node_t*>::const_iterator it(offsetmap_.find(offset));
    if (it == offsetmap_.end())
        return 0;

    return it->second;
}

void
BinaryMap::insert(void *ptr, size_t offset, size_t size)
{
     binary_node_t* node = new binary_node_t;
     node->vals = ptr;
     node->offset = offset;
     node->size = size;

     nodes_.push_back(node);
     offsetmap_[offset] = node;
     ptrmap_[ptr] = node;
}

BinaryMap::iterator
BinaryMap::begin()
{
    return nodes_.begin();
}

BinaryMap::iterator
BinaryMap::end()
{
    return nodes_.end();

}
