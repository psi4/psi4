// This is a generated file. Do not commit to CVS.
#include <cca.h>
#include <stdPorts.h>
#include <EG.h>
#include "gacca.h"
#include "GAServices.h"

extern "C" {

classic::gov::cca::Component *create_GAServices() {
	classic::gov::cca::Component *wanker;
	GA::GAServices *component;
	component = new GA::GAServices();
	wanker = dynamic_cast<classic::gov::cca::Component *>(component);
	return wanker;
}    

char **getComponentList() {
	static char *list[2];
	list[0] = "create_GAServices GA::GAServices";
	list[1] = 0;
	return list;
}

}
static char id[]="$Id: GAServices_wrapper.cxx,v 1.1 2003-08-01 00:41:54 manoj Exp $";
