// This is a generated file. Do not commit to CVS.
#include <cca.h>
#include <stdPorts.h>
#include <EG.h>
#include "TestComponent.h"

extern "C" {

classic::gov::cca::Component *create_TestComponent() {
	classic::gov::cca::Component *wanker;
	TestComponent *component;
	component = new TestComponent();
	wanker = dynamic_cast<classic::gov::cca::Component *>(component);
	return wanker;
}    

char **getComponentList() {
	static char *list[2];
	list[0] = "create_TestComponent TestComponent";
	list[1] = 0;
	return list;
}

}
static char id[]="$Id: TestComponent_wrapper.cxx,v 1.1 2003-08-01 00:41:02 manoj Exp $";
