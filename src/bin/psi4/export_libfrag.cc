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
    using namespace psi::LibFrag;
	using namespace boost::python;
	class_<LibFragHelper>("LibFragHelper")
	    .def("FragHelper",&LibFragHelper::Fragment_Helper)
	    .def("EmbedHelper",&LibFragHelper::Embed_Helper)
	    .def("NMerHelper",&LibFragHelper::NMer_Helper)
	    .def("CapHelper",&LibFragHelper::Cap_Helper)
	    .def("CalcEnergy",&LibFragHelper::CalcEnergy)
	    .def("GetNNMers",&LibFragHelper::GetNNMers)
	    .def("GetGhostsNMerN",&LibFragHelper::GetGhostNMerN)
	    .def("WriteMOs",&LibFragHelper::WriteMOs)
	    .def("GatherData",&LibFragHelper::GatherData)
	    .def("GetNFrags",&LibFragHelper::GetNFrags)
	    .def("RunFrags",&LibFragHelper::RunFrags)
	    .def("IsGMBE",&LibFragHelper::IsGMBE)
	    .def("Sync",&LibFragHelper::Synchronize)
	    .def("Iterate",&LibFragHelper::Iterate)
	    .def("GetNMerN",&LibFragHelper::GetNMerN);
}

