#include "serialabstract.h"
#include "serialimpl.h"
#include "xmlarchive.h"
#include <sstream>

using namespace smartptr;
using namespace std;

#undef heisenbug
#define heisenbug cout << "Heisenbug: " << __FILE__ << " " << __LINE__ << endl

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

map<string, SerialRuntime::create_function*> *SerialRuntime::classlist_ = 0;

Serializable::Serializable()
    : Countable(), archive_id_(0), runtime_id_(0)
{
}



Serializable::Serializable(const XMLArchivePtr& arch)
    : Countable(), archive_id_(0), runtime_id_(0)
{ 

    if (!arch)
    {
        cerr << "Serializable class received null node in constructor" << endl;
        abort();
    }

    serial_load(archive_id);
    serial_load(runtime_id);
    arch->registerObject(this);
}

void
Serializable::serialize(const XMLArchivePtr& arch) const
{
    serial_save(archive_id);
    serial_save(runtime_id);
}

void
Serializable::print(std::ostream &os) const
{
    std::cerr << "print not implemented for " << rtinfo_->classname() << endl;
    abort();
}

void
Serializable::setArchiveID(unsigned long id)
{
    if (archive_id_)
    {
        cerr << "Object already has archive id " << archive_id_
             << " and can not be reassigned. Only one archive should"
             << " link to an object at a time"
             << endl;
        abort();
    }
    archive_id_ = id;
}

void
Serializable::setRuntimeID(unsigned long id)
{
    if (runtime_id_)
    {
         cerr << "Object already has runtime id " << runtime_id_
             << " and can not be reassigned. Only one archive should"
             << " link to an object at a time"
             << endl;
        abort();
    }
    runtime_id_ = id;
}

void
Serializable::resetArchiveID() const
{
    Serializable* me = const_cast<Serializable*>(this);
    me->archive_id_ = 0;
}

unsigned long
Serializable::getArchiveID() const
{
    return archive_id_;
}

unsigned long
Serializable::getRuntimeID() const
{
    return runtime_id_;
}

Serializable::~Serializable() //delete from the registry
{
}

const SerialRuntime* 
Serializable::runtime_info() const
{
    return rtinfo_;
}

SerialRuntime::SerialRuntime(const std::string& name, create_function* fxn)
{
    if (classlist_ == 0)
        classlist_ = new map<string, create_function*>;

    classlist_->insert(pair<string,create_function*>(name,fxn));
}

boost::intrusive_ptr<Serializable>
SerialRuntime::getObject(
    const XMLArchivePtr& arch,
    const std::string& classname
)
{
    map<string, create_function*>::const_iterator it(classlist_->find(classname));
    if (it == classlist_->end())
    {
        cerr << classname << " is not a valid serializable class name" << endl;
        abort();
    }
    Serializable* obj = (*it->second)(arch);
    return obj;
}

void
SerialRuntime::free()
{
    delete classlist_;
}

SerialMap::SerialMap()
{
}

bool
SerialMap::has(unsigned long id) const
{
    registry_map::const_iterator it(registry_.find(id));
    return it != registry_.end();
}

SerializablePtr
SerialMap::get(unsigned long id) const
{
    registry_map::const_iterator it(registry_.find(id));
    if (it == registry_.end())
        return 0;

    return boost::const_pointer_cast<Serializable>(it->second);
}

void
SerialMap::add(unsigned long id, const ConstSerializablePtr& obj)
{
    registry_[id] = obj;
}

std::ostream&
smartptr::operator<< (std::ostream& os, Serializable* s)
{
    if (s)
        s->print(os);
    else
        os << "null object" << endl;
    return os;
}
