/*
 * ParseInput.h
 *
 *  Created on: May 19, 2014
 *      Author: richard
 */

#ifndef PARSEINPUT_H_
#define PARSEINPUT_H_
#include <sstream>
#include <boost/shared_ptr.hpp>
#include "../libpsio/psio.hpp" //This is the default psiomanager definition
namespace psi{
class Molecule;
}
namespace LibFrag{
class Set;

///The InputManip object aids in manipulating the one central Psi input, into millions of fragment inputs
class InputManip{
	private:
		///Where we are writing our directories to
		std::string ScratchDir;
		//boost::shared_ptr<PSIO> OldPSIO;
		void clone(const InputManip& other){this->ScratchDir=other.ScratchDir;}
	public:
		///Default is to set ScratchDir to whatever it was
		InputManip();
		///Copies other InputManip Object
		InputManip(const InputManip& other){this->clone(other);}
		///Assignment operator
		const InputManip& operator=(const InputManip& other){if(this!=&other)this->clone(other);return *this;}

		///Switches Psi's ScratchDir to whatever you give it
		void SwitchScratch(const std::string& folder);

		///Switches Psi's ScratchDir back to it's default value
		void SwitchBack(){SwitchScratch(ScratchDir);}

		///Given a name, returns the path to that folder in the ScratchDir
		std::string FolderPath(const std::string& name)const;

		///Given a path, makes the directory with that path
		void MakeDir(const std::string& path)const;

		///Makes our new inputfile in Folder, with the fragment in Frag and atoms in Mol
		void MakeInput(const std::string& Folder,const std::string& MolName,const std::string& FileBase,
				const Set* Frag,const psi::Molecule* Mol)const;

		///Removes the molecule section and .MBE() call, rest of file is returned as FileBase and the name given to the molecule is MolName
		void ParseInput(std::stringstream &FileBase,std::stringstream &MolName)const;
};

}


#endif /* PARSEINPUT_H_ */
