/*
 * ParseInput.cc
 *
 *  Created on: May 19, 2014
 *      Author: richard
 */
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include<stdio.h>
#include<stdlib.h>
#include <fstream>
#include <iomanip>
#include <time.h>
#include "../libmints/molecule.h"
#include "../../bin/psi4/psi4.h" //This is where the inputfile is defined
#include "InputManip.h"
#include "Set.h"
typedef psi::Molecule Mol;

namespace LibFrag{

typedef Set* pSet;
typedef std::stringstream SS;
typedef std::string str;

void FromFILE2stream(SS& InputFile){
	fseek(psi::infile , 0 , SEEK_END);
	long lSize = ftell (psi::infile);
	char* buffer=new char[lSize+1];
	rewind(psi::infile);
	size_t length=fread(buffer,1,lSize,psi::infile);
	buffer[lSize]=0;
	InputFile<<buffer;
	delete [] buffer;
}

void InputManip::MakeInput(const str& Folder,const str& MolName,const str& FileBase,const Set* Frag,const Mol* AMol)const{
	std::ofstream InputFile(Folder.c_str());
	//Write the molecule line
	InputFile<<"molecule "<<MolName<<" {"<<std::endl;
	InputFile<<"units bohr"<<std::endl;
	for(int i=0;i<Frag->size();i++){
		int MyAtom=(*Frag)[i];
		InputFile<<std::setprecision(15)<<AMol->flabel(MyAtom)<<" ";
		InputFile<<std::setprecision(15)<<AMol->fx(MyAtom)<<" "<<std::setprecision(15)<<AMol->fy(MyAtom)<<" ";
		InputFile<<std::setprecision(15)<<AMol->fz(MyAtom)<<std::endl;
	}
	InputFile<<"}"<<std::endl;
	std::string temp(MolName);
	boost::trim(temp);
	//InputFile<<"psi4.set_global_option(\"basis\",\"sto-3G"<<temp<<"\")"<<std::endl;
	InputFile<<FileBase;
	InputFile.close();
}

void InputManip::ParseInput(SS& FileBase,SS& MolName)const {
	SS InputFile;
	FromFILE2stream(InputFile);
	InputFile<<std::endl;
	str line;
	bool InMolecule=false;
	while(std::getline(InputFile,line)){
		str LINE(boost::to_upper_copy(line));
		if(LINE.substr(0, LINE.find(' '))=="MOLECULE"){
			InMolecule=true;
			bool isbrack=false;
			int i=line.find('m')+8;
			while(line[i]!='{'){
				if(line[i]!=' ')MolName<<line[i++];
				else i++;
			}
		}
		//.MBE() is already case sensitive
		bool isMBE=(str::npos!=line.find(".MBE("));
		bool isRaise=(str::npos!=line.find("raise"));
		//bool isBasis=(str::npos!=line.find("basis"));
		bool isAbove=(isMBE||isRaise);
		if(!InMolecule&&!isAbove){
			FileBase<<line<<std::endl;
		}
		else if(line.find('}')!=str::npos)InMolecule=false;
	}
}

InputManip::InputManip(){
	ScratchDir=psi::_default_psio_manager_->get_default_path();
}

void InputManip::MakeDir(const str &name)const{
	boost::filesystem::path dir(name);
	boost::filesystem::create_directory(dir);
}

void InputManip::SwitchScratch(const str &folder){
	psi::_default_psio_manager_->set_default_path(folder);
}

std::string InputManip::FolderPath(const str &name)const{
	SS folder;
	folder<<ScratchDir<<name;
	return folder.str();
}
}

