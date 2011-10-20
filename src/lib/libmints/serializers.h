#ifndef libmints_serializers_h
#define libmints_serializers_h

// Ensure we have access to everything.
#include "mints.h"
#include "coordentry.h"

// Include boost headers to help us out.
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/binary_object.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/export.hpp>

// Add our serializers to the boost::serialization namespace:
namespace boost { namespace serialization {

////////////////////////////////////////////////////////////////////////////////
//
// Psi Dimension
//
////////////////////////////////////////////////////////////////////////////////
template<class Archive>
void save(Archive & ar, const psi::Dimension & t, unsigned int /*version*/)
{
    ar & make_nvp("name", t.name());
    ar & make_nvp("rank", t.n());
    const int *data = t; // need to come up with a better scheme
    ar & make_nvp("data", make_array(data, t.n()));
}

template<class Archive>
void load(Archive & ar, psi::Dimension & t, unsigned int /*version*/)
{
    int rank;
    std::string name;

    ar & make_nvp("name", name);
    ar & make_nvp("rank", rank);
    t.init(name, rank);
    int *data = t;
    ar & make_nvp("data", make_array(data, t.n()));
}

////////////////////////////////////////////////////////////////////////////////
//
// Psi Matrix
//
////////////////////////////////////////////////////////////////////////////////
template<class Archive>
void save(Archive & ar, const psi::Matrix & t, unsigned int /*version*/)
{
    ar & make_nvp("name", t.name());
//    int temp = t.nirrep();
    ar & make_nvp("nirrep", t.nirrep());
//    temp = t.symmetry();
    ar & make_nvp("symmetry", t.symmetry());
    ar & make_nvp("rowspi", t.rowspi());
    ar & make_nvp("colspi", t.colspi());

    for (int h=0; h < t.nirrep(); ++h) {
        if (t.rowdim(h) != 0 && t.coldim(h))
            ar & make_nvp("data", make_binary_object(t.pointer(h)[0], sizeof(double)*t.rowdim(h)*t.coldim(h)));
    }
}

template<class Archive>
void load(Archive &ar, psi::Matrix & t, unsigned int /*version*/)
{
    std::string name;
    psi::Dimension rows, cols;
    int nirrep, symmetry;
    ar & make_nvp("name", name);
    ar & make_nvp("nirrep", nirrep);
    ar & make_nvp("symmetry", symmetry);
    ar & make_nvp("rowspi", rows);
    ar & make_nvp("colspi", cols);

    t.init(rows, cols, name, symmetry);
    for (int h=0; h < t.nirrep(); ++h) {
        if (t.rowdim(h) != 0 && t.coldim(h))
            ar & make_nvp("data", make_binary_object(t.pointer(h)[0], sizeof(double)*t.rowdim(h)*t.coldim(h)));
    }
}


////////////////////////////////////////////////////////////////////////////////
//
// Psi Vector3
//
////////////////////////////////////////////////////////////////////////////////
template<class Archive>
void save(Archive & ar, const psi::Vector3 & t, unsigned int /*version*/)
{
    ar & t[0] & t[1] & t[2];
}

template<class Archive>
void load(Archive & ar, psi::Vector3 & t, unsigned int /*version*/)
{
    ar & t[0] & t[1] & t[2];
}

////////////////////////////////////////////////////////////////////////////////
//
// Psi CoordEntry
//
////////////////////////////////////////////////////////////////////////////////
template<class Archive>
void save(Archive & ar, const psi::CoordEntry & t, unsigned int /*version*/)
{
    ar & make_nvp("entry_number", t.entry_number_);
    ar & make_nvp("symbol", t.symbol_);
    ar & make_nvp("label", t.label_);
    ar & make_nvp("ghosted", t.ghosted_);
    ar & make_nvp("Z", t.Z_);
    ar & make_nvp("charge", t.charge_);
    ar & make_nvp("mass", t.mass_);
    ar & make_nvp("coordinates", t.coordinates_);
}

template<class Archive>
void load(Archive & ar, psi::CoordEntry & t, unsigned int /*version*/)
{
    ar & make_nvp("entry_number", t.entry_number_);
    ar & make_nvp("symbol", t.symbol_);
    ar & make_nvp("label", t.label_);
    ar & make_nvp("ghosted", t.ghosted_);
    ar & make_nvp("Z", t.Z_);
    ar & make_nvp("charge", t.charge_);
    ar & make_nvp("mass", t.mass_);
    ar & make_nvp("coordinates", t.coordinates_);
}

////////////////////////////////////////////////////////////////////////////////
//
// Psi CartesianEntry
//
////////////////////////////////////////////////////////////////////////////////
//template<class Archive>
//void save(Archive & ar, const psi::CartesianEntry & t, unsigned int /*version*/)
//{
//    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(psi::CoordEntry);
//}

//template<class Archive>
//void load(Archive & ar, psi::CartesianEntry & t, unsigned int /*version*/)
//{
//    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(psi::CoordEntry);
//}

////////////////////////////////////////////////////////////////////////////////
//
// Psi ZMatrixEntry
//
////////////////////////////////////////////////////////////////////////////////
//template<class Archive>
//void save(Archive & ar, const psi::ZMatrixEntry & t, unsigned int /*version*/)
//{
//    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(psi::CoordEntry);
//}

//template<class Archive>
//void load(Archive & ar, psi::ZMatrixEntry & t, unsigned int /*version*/)
//{
//    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(psi::CoordEntry);
//}

////////////////////////////////////////////////////////////////////////////////
//
// Psi Molecule
//
////////////////////////////////////////////////////////////////////////////////

template<class Archive>
void save(Archive & ar, const psi::Molecule & t, unsigned int /*version*/)
{
    ar & make_nvp("name", t.name());
    ar & make_nvp("nallatom", t.nallatom());
    std::string sym = t.point_group()->symbol();
    ar & make_nvp("symmetry", sym);
    psi::Matrix full_geometry = t.full_geometry();
//    ar & make_nvp("coordinates", full_geometry);
//    for (int n=0; n<t.nallatom(); ++n) {
//        const boost::shared_ptr<psi::CoordEntry> temp = t.atom_entry(n);
//        ar & temp;
//    }
}

template<class Archive>
void load(Archive & ar, psi::Molecule & t, unsigned int /*version*/)
{
}

}} // end of namespace boost and serialization

//
// BOOST_SERIALIZATION_SPLIT_FREE
// These are needed to tell boost::serialization that our serializers
// are separate from the actual class (FREE) and separated into save
// and load functions (SPLIT)
//
// BOOST_SERIALIZATION_SHARED_PTR
// If you intent to serialize a shared_ptr object you MUST have a
// BOOST_SERIALIZATION_SHARED_PTR(class_name) listed or it will not
// work.
//
//BOOST_SERIALIZATION_SPLIT_FREE(psi::CartesianEntry)
//BOOST_SERIALIZATION_SHARED_PTR(psi::CartesianEntry)
//BOOST_CLASS_EXPORT(psi::CartesianEntry)
//BOOST_SERIALIZATION_SPLIT_FREE(psi::CoordEntry)
//BOOST_SERIALIZATION_SHARED_PTR(psi::CoordEntry)
BOOST_SERIALIZATION_SPLIT_FREE(psi::Dimension)
BOOST_SERIALIZATION_SPLIT_FREE(psi::Matrix)
BOOST_SERIALIZATION_SHARED_PTR(psi::Matrix)
//BOOST_SERIALIZATION_SPLIT_FREE(psi::BasisSet)
//BOOST_SERIALIZATION_SHARED_PTR(psi::BasisSet)
BOOST_SERIALIZATION_SPLIT_FREE(psi::Molecule)
BOOST_SERIALIZATION_SHARED_PTR(psi::Molecule)
//BOOST_SERIALIZATION_SPLIT_FREE(psi::ZMatrixEntry)
//BOOST_SERIALIZATION_SHARED_PTR(psi::ZMatrixEntry)
//BOOST_CLASS_EXPORT(psi::ZMatrixEntry)

#endif // libmints_serializers_h
