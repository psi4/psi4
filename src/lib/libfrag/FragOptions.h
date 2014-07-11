/*
 * FragOptions.h
 *
 *  Created on: Jun 3, 2014
 *      Author: richard
 */

#ifndef FRAGOPTIONS_H_
#define FRAGOPTIONS_H_


#include <string>
namespace LibFrag{
enum FragMethods {USER_DEFINED,BOND_BASED,DISTANCE_BASED};
enum EmbedMethods {NO_EMBED,POINT_CHARGE,ITR_POINT_CHARGE,DENSITY,
                   ITR_DENSITY};
enum CapMethods {NO_CAPS,H_REPLACE,H_SHIFTED};
enum BSSEMethods {NO_BSSE,FULL,MBCPN,VMFCN};


class FragOptions{
	private:
        ///Sets all members to the default options
        void DefaultOptions();
        ///Actually preforms copy
		void copy(const FragOptions& other);
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


		///Constructor
		FragOptions(){DefaultOptions();}
		~FragOptions(){}
		FragOptions(const FragOptions& other){this->copy(other);}
		FragOptions operator=(const FragOptions& other){
		   if(this!=&other)this->copy(other);return *this;
		}
};

}


#endif /* FRAGOPTIONS_H_ */
