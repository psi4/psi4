#include <iostream>
#include <cstdio>
#include <cmath>

#include <cca.h>
#include <stdPorts.h>
#include <EG.h>
#include "jc++/jc++.h"
#include "jc++/util/jc++util.h"
#include "parameters/parametersStar.h"
#include "util/IO.h"

#include "../gacca.h"
#include "TestComponent.h"


using std::cout;
using std::printf;
using std::sqrt;

#define N   100
#define DIM 2
#define GA_DATA_TYPE MT_F_DBL

#define CHECKERR(err) if (err < 0) cerr<<"Line "<<__LINE__<<": Error # " << err << endl;

void doWork(GA::GAClassicPort *ga_port);
void doWork_DADF(GA::GAClassicPort * ga_port, 
		 classic::gov::cca::DistArrayTemplFactoryPort * templFactory, 
		 classic::gov::cca::DistArrayDescrFactoryPort * descrFactory); 


TestComponent::TestComponent() {
  svc = 0;
}

TestComponent::~TestComponent() {
  svc = 0;
}

void 
TestComponent::setServices(classic::gov::cca::Services *cc) {

  svc = cc;
  
  // Contact the PrintfService
  classic::gov::cca::PortInfo* pinfo = cc->createPortInfo("pSvc", "gov.cca.JPrintfService", 0);
  cc->registerUsesPort(pinfo);
  pinfo = 0;
  pfp = dynamic_cast<classic::gov::cca::JPrintfPort*>(cc->getPort("pSvc"));
  CHECKDC(pfp);
  if(pfp == 0) {
    cc->addProvidesPort(this, cc->createPortInfo("DEAD=NoJPrintf", "classic::gov::cca::GoPort", 0));
    ::printf("!!! No JPrintfService available from framework.");
    return;
  }
  
  // register "GA" uses port
  pinfo = svc->createPortInfo("ga_classic_port", "GA::GAClassicPort",0);
  svc->registerUsesPort(pinfo);
  pinfo = 0;
  // register "GA" uses port
  pinfo = svc->createPortInfo("TemplateFactory",
			      "DistArrayTemplFactoryPort",0);
  svc->registerUsesPort(pinfo);
  pinfo = 0;
  // register "GA" uses port
  pinfo = svc->createPortInfo("DescriptorFactory",
                              "DistArrayDescrFactoryPort",0);
  svc->registerUsesPort(pinfo);
  pinfo = 0;

  // Provides "Go" Port
  pinfo = svc->createPortInfo("go", "classic::gov::cca::GoPort",0);
  svc->addProvidesPort(this, pinfo);
  pinfo = 0;  
}


int 
TestComponent::go() {
  classic::gov::cca::Port *port = 0, *pTemplate = 0, *pDescr = 0; 
  
  /* get the ga classic port */
  port = svc->getPort("ga_classic_port");
  if (pfp && port == 0) {
    pfp->en("TestComponent::go(): ga_clasic_port not apparently connected");
  }
  /* get the ga dadf template port */
  pTemplate = svc->getPort("TemplateFactory");
  if (pfp && pTemplate == 0) {
    pfp->en("TestComponent::go(): TemplateFactory not apparently connected");
  }
  /* get the ga dadf descriptor port */
  pDescr = svc->getPort("DescriptorFactory");
  if (pfp && pDescr == 0) {
    pfp->en("TestComponent::go(): DescriptorFactory not apparently connected");
  }

  /* type-casting */
  GA::GAClassicPort *ga_port;
  ga_port = dynamic_cast < GA::GAClassicPort *> (port);
  if(ga_port == 0) {
    if (pfp) {
      pfp->en("BSTest::go(): ga_classic_port not castable to correct type!");
    }
    return -1;
  }
  /* type-casting */
  classic::gov::cca::DistArrayTemplFactoryPort * templFactory;
  templFactory = dynamic_cast < classic::gov::cca::DistArrayTemplFactoryPort *> (pTemplate);
  if(templFactory == 0) {
    if (pfp) {
      pfp->en("BSTest::go(): TemplateFactory not castable to correct type!");
    }
    return -1;
  }
  /* type-casting */
  classic::gov::cca::DistArrayDescrFactoryPort * descrFactory;
  descrFactory = dynamic_cast < classic::gov::cca::DistArrayDescrFactoryPort *> (pTemplate);
  if(descrFactory == 0) {
    if (pfp) {
      pfp->en("BSTest::go(): DescriptorFactory not castable to correct type!");
    }
    return -1;
  }
  
 
  /****************************************/
 
  int me = ga_port->nodeid();
  int nproc = ga_port->nodes();
  int len; char proc_name[MPI_MAX_PROCESSOR_NAME];

  MPI_Get_processor_name(proc_name, &len);
  cout << proc_name << " : Rank = " << me << " : Size = " << nproc << "\n";

  cout << "\n---------------------------------------------------\n";
  cout << "           TESTING GA CLASSIC COMPONENT\n";
  cout << "---------------------------------------------------\n\n";
  doWork(ga_port);
  cout << "After doWork()\n";

  cout << "\n---------------------------------------------------\n";
  cout << "           TESTING GA DADF COMPONENT\n";
  cout << "---------------------------------------------------\n\n";
  doWork_DADF(ga_port, templFactory, descrFactory);
  cout << "After doWork_DADF()\n";
  
  /****************************************/

  cout << "\nReleasing Port(s) ...\n";
  svc->releasePort("ga_classic_port");
  svc->releasePort("TemplateFactory");
  svc->releasePort("DescriptorFactory");
  
  return 0;
  
}

void doWork(GA::GAClassicPort *ga_port) {

  int ONE=1;/* useful constants */
  int n=N, type=MT_F_DBL;
  int me=ga_port->nodeid(), nproc=ga_port->nodes();
  int i, row;
  int dims[2]={N,N};
  int lo[2], hi[2];
  
  /* Note: on all current platforms DoublePrecision == double */
  double buf[N], err, alpha, beta;

  if(me==0)printf("size = %d\n", nproc);
  
  if(me==0)printf("Creating matrix A\n");
  GA::GlobalArray *g_a = ga_port->createGA(type, 2, dims, "A", NULL);
  if(me==0)printf("OK\n");
  
  if(me==0)printf("Creating matrix B\n");
  /* create matrix B  so that it has dims and distribution of A*/
  GA::GlobalArray *g_b = ga_port->createGA(g_a, "B");
  if(me==0)printf("OK\n");
  
  g_a->zero();   /* zero the matrix */
  
  if(me==0)printf("Initializing matrix A\n");
  /* fill in matrix A with random values in range 0.. 1 */ 
  lo[1]=0; hi[1]=n-1;
  for(row=me; row<n; row+= nproc){
    /* each process works on a different row in MIMD style */
    lo[0]=hi[0]=row;   
    for(i=0; i<n; i++) buf[i]=sin((double)i + 0.1*(row+1));
    g_a->put(lo, hi, buf, &n);
  }
  
  
  if(me==0)printf("Symmetrizing matrix A\n");
  g_a->symmetrize();   /* symmetrize the matrix A = 0.5*(A+A') */
  
  
  /* check if A is symmetric */ 
  if(me==0)printf("Checking if matrix A is symmetric\n");
  g_a->transpose(g_b); /* B=A' */
  alpha=1.; beta=-1.;
  g_b->add(&alpha, g_a, &beta, g_b);  /* B= A - B */
  err= g_b->ddot(g_b);
  
  if(me==0)printf("Error=%lf\n",(double)err);
  
  if(me==0)printf("\nChecking atomic accumulate \n");
  
  g_a->zero();   /* zero the matrix */
  for(i=0; i<n; i++) buf[i]=(double)i;
  
  /* everybody accumulates to the same location/row */
  alpha = 1.0;
  row = n/2;
  lo[0]=hi[0]=row;
  lo[1]=0; hi[1]=n-1;
  g_a->acc(lo, hi, buf, &ONE, &alpha );
  ga_port->sync();
  
  if(me==0){ /* node 0 is checking the result */
    
    g_a->get(lo, hi, buf,&ONE);
    for(i=0; i<n; i++) if(buf[i] != (double)nproc*i)
      ga_port->error("failed: column=",i);
    printf("OK\n\n");
    
  }
  
  g_a->destroy();
  g_b->destroy();
  
}

void doWork_DADF(GA::GAClassicPort * ga_port, 
		 classic::gov::cca::DistArrayTemplFactoryPort * templFactory, 
		 classic::gov::cca::DistArrayDescrFactoryPort * descrFactory) {

  int i;
  int me=ga_port->nodeid();//  nproc=ga_port->nodes(); 
  int lo[DIM] = {0, 0};
  int hi[DIM] = {N, N};
  int chunk[DIM] = {2,2};
  int topology[DIM] = {2, 2};
  DistArrayTemplate * templ;
  DistArray * darr;
  
  /********* Creating Template ********/  
  if(!me) cout << "\nCreating a DADF Template:\n\n";
  templ = templFactory->createTemplate( "one" );  
  DistArrayTemplate::DistType dist[DIM] ;

  for(i=0; i<DIM; ++i) dist[i] = DistArrayTemplate::Block;
  int ierr = templ->setRank(DIM); CHECKERR(ierr);
  ierr = templ->setGlobalBounds(lo, hi); CHECKERR(ierr);
  ierr = templ->setProcTopology(topology); CHECKERR(ierr);
  ierr = templ->setDistType(dist); CHECKERR(ierr);
  for(i=0; i<DIM; ++i) {
    ierr = templ->setDistParameters(i, chunk[i], i); CHECKERR(ierr);
  }
  templ->commit();
  templ->printTemplate();
  
  /** Create Distributed Array **/
  if(!me) cout << "\n\nCreating a GA Style Distributed array:\n\n";
  darr = descrFactory->createArray("My_Array");
  ierr = darr->setDataType(DistArray::stv_Double); CHECKERR(ierr);
  ierr = darr->setTemplate(templ); CHECKERR(ierr);
  ierr = darr->setIdentityAlignmentMap(); CHECKERR(ierr);
  ierr = darr->commit(); CHECKERR(ierr);
  darr->printArrayDistribution();

  templFactory->destroyTemplate(templ);
  descrFactory->destroyArray(darr);  
}
