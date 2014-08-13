#ifndef _TESTCOMPONENT_H
#define _TESTCOMPONENT_H

class TestComponent : public virtual classic::gov::cca::Component, 
              public virtual classic::gov::cca::GoPort {
  
public: 
  TestComponent();
  ~TestComponent();

  // Component interface:
  virtual void setServices(classic::gov::cca::Services *cc);
  
  // GoPort interface:
  virtual int go();
  
private:
  classic::gov::cca::Services *svc;
  classic::gov::cca::JPrintfPort *pfp;
};


#endif // _TESTCOMPONENT_H
