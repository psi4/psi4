/*
 * export_libfrag.cc
 *
 *  Created on: Jun 4, 2014
 *      Author: richard
 */

#include "../lib/libfrag/LibFragDriver.h"
void export_libfrag(){
    using namespace psi::LibFrag;
	using namespace boost::python;
	class_<LibFragDriver>("LibFragDriver");
}

