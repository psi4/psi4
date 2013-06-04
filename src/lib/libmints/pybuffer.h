/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#ifndef _psi_src_lib_libmints_pybuffer_h
#define _psi_src_lib_libmints_pybuffer_h

#include <boost/python/dict.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/list.hpp>
#include <sstream>

#define xtostring( x ) dynamic_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

namespace psi {

inline int is_big_endian(void)
{
    union {
        uint32_t i;
        char c[4];
    } bint = {0x01020304};

    return bint.c[0] == 1;
}

template <typename vectType>
inline boost::python::tuple vector_to_tuple(std::vector<vectType> v){
	boost::python::list tmp;
	for(int i = 0; i < v.size(); ++i) {
		tmp.append(v[i]);
	}
	return boost::python::tuple(tmp);
}

/*! \ingroup MINTS
 *  \class PyBuffer
 *  \brief A thin wrapper to the NumPy array interface.
 *
 *  PyBuffer instances can be used to expose the pointer of a given
 *  C++ buffer to the NumPy array interface without copying any data.
 *  For more details, see http://docs.scipy.org/doc/numpy/reference/arrays.interface.html
 */
template <class T>
class PyBuffer {

	boost::python::dict interface_dict_;

	void init_data_type_()
	{
		std::string data_str = "";
		// TODO this should be determined at compile time, or at worst only once at runtime
		if(is_big_endian())
			data_str += ">";
		else
			data_str += "<";
		data_str += data_typestr() + xtostring((int)sizeof(T));
		interface_dict_["typestr"] = data_str;
	}

	void init_(T* data, std::vector<long> shape, bool read_only)
	{

		set_shape(shape);
		read_only_ = read_only;
		set_data_ptr(data);
		init_data_type_();
	}

	void init_(T* data, boost::python::tuple shape, bool read_only)
	{
		set_shape(shape);
		read_only_ = read_only;
		set_data_ptr(data);
		init_data_type_();
	}

	inline std::string data_typestr();

	void set_data_ptr(T* dataptr){
		data_ptr_ = dataptr;
		interface_dict_["data"] = boost::python::make_tuple((long)dataptr, read_only_);
	}

	T** data_tracker_;
	T* data_ptr_;
	bool read_only_;

public:

	PyBuffer() {};

	PyBuffer(T* data, bool read_only = true) :
		interface_dict_()
	{
		data_tracker_ = &data;
		init_(data, boost::python::make_tuple(0), read_only);
	}

	PyBuffer(T* data, boost::python::tuple shape, bool read_only = true) :
		interface_dict_()
	{
		data_tracker_ = &data;
		init_(data, shape, read_only);
	}

	PyBuffer(T* data, std::vector<long> shape, bool read_only = true) :
		interface_dict_()
	{
		data_tracker_ = &data;
		init_(data, shape, read_only);
	}

	/// "Track" a changing buffer by passing a pointer to a pointer
	PyBuffer(T** data_tracker, bool read_only) :
		interface_dict_()
	{
		data_tracker_ = data_tracker;
		init_(*data_tracker_, boost::python::make_tuple(0), read_only);
	}

	~PyBuffer() {};

	boost::python::dict array_interface() {
		// First make sure our tracking pointer is up to date
		set_data_ptr(*data_tracker_);
		// Now return the array interface dict
		return interface_dict_;
	}
	void set_shape(std::vector<long> shape) {
		set_shape(vector_to_tuple(shape));
	}
	void set_shape(boost::python::tuple shape) {
		interface_dict_["shape"] = shape;
	}

};


template <> inline std::string PyBuffer<double>::data_typestr(){ return "f"; }
template <> inline std::string PyBuffer<int>::data_typestr(){ return "i"; }
template <> inline std::string PyBuffer<long>::data_typestr(){ return "i"; }
template <> inline std::string PyBuffer<long long>::data_typestr(){ return "i"; }
template <> inline std::string PyBuffer<bool>::data_typestr(){ return "b"; }

} // end namespace psi

#endif /* _psi_src_lib_libmints_pybuffer_h */
