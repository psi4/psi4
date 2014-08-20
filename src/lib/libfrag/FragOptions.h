/*
 * FragOptions.h
 *
 *  Created on: Jun 3, 2014
 *      Author: richard
 */

#ifndef FRAGOPTIONS_H_
#define FRAGOPTIONS_H_


#include <string>
#include <boost/shared_ptr.hpp>
#include "LibFragTypes.h"

namespace psi{
namespace LibFrag{
enum FragMethods {USER_DEFINED,BOND_BASED,DISTANCE_BASED};
enum EmbedMethods {NO_EMBED,POINT_CHARGE,ITR_POINT_CHARGE,DENSITY,
                   ITR_DENSITY};
enum CapMethods {NO_CAPS,H_REPLACE,H_SHIFTED};
enum BSSEMethods {NO_BSSE,FULL,MBCPN,VMFCN};

class GMBE;
class Fragmenter;
class BSSEer;
class Capper;
class Embedder;
class FragOptions{
	private:
        ///Sets all members to the default options
        void DefaultOptions();
        ///Actually preforms copy
		void copy(const FragOptions& other);
		std::string ToString(const FragMethods& F);
		std::string ToString(const EmbedMethods& E);
		std::string ToString(const CapMethods& C);
		std::string ToString(const BSSEMethods& B);
	public:
        ///Options
		FragMethods FMethod;
        EmbedMethods EMethod;
        CapMethods CMethod;
        BSSEMethods BMethod;
        int MBEOrder;

        ///Given a lowercase string, the following will set things up right
        void SetFMethod(const std::string& FMethodIn);
        void SetEMethod(const std::string& EMethodIn);
        void SetCMethod(const std::string& CMethodIn);
        void SetBMethod(const std::string& BMethodIn);

        ///Returns a pointer to the Fragmenter determined by this->FMethod
        boost::shared_ptr<Fragmenter> MakeFragFactory()const;
        /** \brief Returns a pointer to the BSSE factory determined by
         *         this->BMethod
         *
         *
         *   Most BSSE methods need to know the number of atoms in the
         *   supersystem so you will need to pass that to this fxn.
         *
         *   \param[in] natoms The number of atoms in the supersystem
         *
         */
        boost::shared_ptr<BSSEer> MakeBSSEFactory(const int natoms) const;
        ///Returns a capping factory based on CMethod
        boost::shared_ptr<Capper> MakeCapFactory(SharedMol& AMol)const;
        ///Returns an embedding factory fased on EMethod
        boost::shared_ptr<Embedder> MakeEmbedFactory(SharedMol& AMol)const;


        ///Nice, pretty printing of all the desired options
        void PrintOptions();
		///Constructor, calls DefaultOptions() for initialization
		FragOptions(){DefaultOptions();}
		~FragOptions(){}
		FragOptions(const FragOptions& other){this->copy(other);}
		FragOptions operator=(const FragOptions& other){
		   if(this!=&other)this->copy(other);return *this;
		}
};

}}// End namespaces


#endif /* FRAGOPTIONS_H_ */
