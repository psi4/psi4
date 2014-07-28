/*
 * export_libfrag.cc
 *
 *  Created on: Jun 4, 2014
 *      Author: richard
 */

#include <boost/python.hpp>
#include "../lib/libfrag/libfrag.h"
#include "../lib/libmints/molecule.h"
void export_libfrag(){
	using namespace boost::python;
	using namespace LibFrag;
	class_<LibFragHelper>("LibFragHelper")
	    .def("FragHelper",&LibFragHelper::Fragment_Helper)
	    .def("EmbedHelper",&LibFragHelper::Embed_Helper)
	    .def("NMerHelper",&LibFragHelper::NMer_Helper)
	    .def("CapHelper",&LibFragHelper::Cap_Helper)
	    .def("CalcEnergy",&LibFragHelper::CalcEnergy)
	    .def("GetNNMers",&LibFragHelper::GetNNMers)
	    .def("GetGhostsNMerN",&LibFragHelper::GetGhostNMerN)
	    .def("WriteMOs",&LibFragHelper::WriteMOs)
	    .def("ReadMOs",&LibFragHelper::ReadMOs)
	    .def("GetNFrags",&LibFragHelper::GetNFrags)
	    .def("RunFrags",&LibFragHelper::RunFrags)
	    .def("Sync",&LibFragHelper::Synchronize)
	    .def("GetNMerN",&LibFragHelper::GetNMerN);
}

