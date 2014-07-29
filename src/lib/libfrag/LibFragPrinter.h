/*
 * LibFragPrinter.h
 *
 *  Created on: May 21, 2014
 *      Author: richard
 */

#ifndef LIBFRAGPRINTER_H_
#define LIBFRAGPRINTER_H_
#include <string>
#include "../../../include/psi4-dec.h"
#include "../libparallel/parallel.h"

class LibFragPrinter{
	private:
		int me;
	public:
		LibFragPrinter(const int I=psi::WorldComm->me()):me(I){}
		void print(const std::string& message){if(me==0){psi::fprintf(psi::outfile,message.c_str());}}
};



#endif /* LIBFRAGPRINTER_H_ */
