/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680


  $Id: world.cc 2363 2011-06-12 02:58:06Z rjharrison $
*/


#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <world/world.h>
#include <world/worldobj.h>
#include <world/worlddc.h>

#ifdef WORLD_TAU_TRACE
#include <TAU.h>
#endif

using namespace madness;
using namespace std;

void test0(World& world) {
    PROFILE_FUNC;
    const size_t n=555;
    char buf[n+2];
    buf[0] = buf[n+1] = 127;

    archive::BufferOutputArchive arout(buf+1,n-2);

    arout & 1L & 99.9;
    arout.store("hello",6);
    arout & 7L & 77.7;
    MADNESS_ASSERT(buf[0]==127 && buf[n+1]==127);

    long i;
    double a;

    archive::BufferInputArchive arin(buf+1,n-2);
    char s[8];
    s[0] = s[7] = 66;
    arin & i & a;
    MADNESS_ASSERT(i==1 && a==99.9);
    arin.load(s+1,6);
    MADNESS_ASSERT(s[0]==66 && s[1]=='h' && s[5]=='o' && s[7]==66);
    arin & i & a;
    MADNESS_ASSERT(i==7 && a==77.7);

    if (world.rank() == 0) print("test0 (serialization to/from buf) seems to be working");
}

class B {
    long b;
public:
    B(long b=0) : b(b) {};
    Void set(long value) {
        b=value;
        return None;
    };
    long get() const {
        return b;
    };
    ~B() {
        print("B destructor");
    };
    template <typename Archive> void serialize(Archive &ar) {
        ar&b;
    }
};


#include <complex>
typedef std::complex<double> double_complex;

class TestTask : public TaskInterface {
public:

    using PoolTaskInterface::run;

    void run(World& world) {
        print("Hi, I am running!");
    }
};

class TTT {
private:
    int state;
public:
    TTT() : state(0) {};

    static Void fred() {
        print("Oops-a-daisy!");
        return None;
    };

    static int mary() {
        return 99;
    };

    static int carl() {
        return 88;
    };

    static int dave(World* world) {
        return world->rank();
    };

    static int bert(int input) {
        return input+1;
    };

    static double_complex sara(double a, const double_complex& b) {
        return a*b;
    };

    static string kate(World* world, const string& msg, double d) {
        ostringstream s;
        s << "Process " << world->rank() << " says '" << msg << "' and "
        << d << " right back at you!";
        return s.str();
    };

    double jody(double a, double b, double c) {
        return a+b+c + state;
    };

    double hugh(vector< Future<int> >& a) {
        double sum = 0.0;
        for (int i=0; i<(int)a.size(); ++i) sum += a[i].get();
        return sum;
    };
};

double dumb(int a1, int a2, int a3, int a4, int a5, int a6, int a7) {
    return a1+a2+a3+a4+a5+a6+a7;
}



void test5(World& world) {
    PROFILE_FUNC;
    int nproc = world.size();
    ProcessID me = world.rank();
    ProcessID right = (me+1)%nproc;
    TaskInterface* task = new TestTask();
    task->inc();
    task->inc();
    world.taskq.add(task);
    print("added the task ... about to dec");
    task->dec();
    task->dec();
    print("fencing");
    world.gop.fence();
    print("done with fence");

    world.taskq.add(TTT::fred);
    Future<int> mary = world.taskq.add(TTT::mary);
    Future<int> carl = world.taskq.add(right,TTT::carl);
    Future<int> dave = world.taskq.add(right,TTT::dave, &world);
    Future<int> bert_input;
    Future<int> bert = world.taskq.add(TTT::bert,bert_input);
    Future<double> sara1;
    Future<double_complex> sara2;
    Future<double_complex> sara = world.taskq.add(TTT::sara,sara1,sara2);
    Future<string> kate2;
    Future<double> kate3;
    Future<string> kate = world.taskq.add(TTT::kate,&world,kate2,kate3);
    Future<string> cute = world.taskq.add(right,TTT::kate,&world,"Boo!",-42.0);
    TTT ttt;
    Future<double> jody = world.taskq.add(ttt,&TTT::jody,1.0,2.0,3.0);
    Future<double> duh = world.taskq.add(me,dumb,0,1,2,3,4,5,6);
    print("done with making futs");

    bert_input.set(7);
    sara1.set(3.0);
    sara2.set(double_complex(2.1,1.2));
    kate2.set(string("Who's your daddy?"));
    kate3.set(3.14);

    vector< Future<int> > futv = future_vector_factory<int>(7);
    Future<double> hugh = world.taskq.add(ttt,&TTT::hugh,futv);
    for (int i=0; i<7; ++i) {
        print("assigning",i,futv[i]);
        futv[i].set(i);
    }

    print("about to fence again");
    world.gop.fence();
    print("finished fence again");

    MADNESS_ASSERT(mary.probe());
    MADNESS_ASSERT(carl.probe());
    MADNESS_ASSERT(dave.probe());
    MADNESS_ASSERT(bert.probe());
    MADNESS_ASSERT(sara.probe());
    MADNESS_ASSERT(kate.probe());
    MADNESS_ASSERT(cute.probe());
    MADNESS_ASSERT(jody.probe());
    MADNESS_ASSERT(hugh.probe());
    MADNESS_ASSERT(duh.probe());

    MADNESS_ASSERT(mary.get() == 99);
    MADNESS_ASSERT(carl.get() == 88);
    MADNESS_ASSERT(dave.get() == right);
    MADNESS_ASSERT(bert.get() == 8);
    MADNESS_ASSERT(hugh.get() == 21.0);
    MADNESS_ASSERT(duh.get() == 21.0);
    print("Sara says",sara.get().real(),sara.get().imag());
    print("Kate says",kate.get());
    print("Cute says",cute.get());
    print("Jody says",jody.get());

    if (me == 0) print("test5 (tasks and futures) OK");
}

class TestBarrier : public TaskInterface {
    volatile int count;
public:

    using PoolTaskInterface::run;

    TestBarrier(const madness::TaskAttributes& attr)
        : TaskInterface(attr)
        , count(0)
    {
        print("Testing barrier with nthread", attr.get_nthread());
    }

    void run(World& world, const TaskThreadEnv& env) {
        // Using the barrier each thread takes turns to update
        // the shared counter.

        env.barrier();

        int nthread = env.nthread();
        int id = env.id();
        for (int i=0; i<100; ++i) {
            for (int p=0; p<nthread; ++p) {
                env.barrier();
                if (p == id) count += (p+1);
            }
        }
        env.barrier();
        if (id == 0)
            print("     result from sum", count, "expected", 100*nthread*(nthread+1)/2);
    }
};

class TimeBarrier : public TaskInterface {
    volatile int count;
public:

    using PoolTaskInterface::run;

    TimeBarrier(const madness::TaskAttributes& attr)
        : TaskInterface(attr)
        , count(0)
    {
        print("Timing barrier with nthread", attr.get_nthread());
    }

    void run(World& world, const TaskThreadEnv& env) {
        // Barrier a zillion times

		for (int i=0; i<1000000; ++i) {
	        env.barrier();
		}
    }
};


// test multithreaded tasks
void test_multi(World& world) {
    // Test the correctness and performance of the barrier
    for (unsigned int i=1; i<=ThreadPool::size()+1; ++i) {
        world.taskq.add(new TestBarrier(TaskAttributes::multi_threaded(i)));
        double start = cpu_time();
        world.taskq.add(new TimeBarrier(TaskAttributes::multi_threaded(i)));
        double used = cpu_time()-start;
        print("barrier took", used*10.0,"micro seconds per call");
        world.gop.fence();
    }
}


class Foo : public WorldObject<Foo> {
    int a;
public:
    Foo(World& world, int a)
            : WorldObject<Foo>(world)
            , a(a) {
        process_pending();
    }

    virtual ~Foo() { }

    int get0() {
        return a;
    }
    int get1(int a1) {
        return a+a1;
    }
    int get2(int a1, char a2) {
        return a+a1+a2;
    }
    int get3(int a1, char a2, short a3) {
        return a+a1+a2+a3;
    }
    int get4(int a1, char a2, short a3, long a4) {
        return a+a1+a2+a3+a4;
    }
    int get5(int a1, char a2, short a3, long a4, short a5) {
        return a+a1+a2+a3+a4+a5;
    }

    int get0c() const {
        return a;
    }
    int get1c(int a1) const {
        return a+a1;
    }
    int get2c(int a1, char a2) const {
        return a+a1+a2;
    }
    int get3c(int a1, char a2, short a3) const {
        return a+a1+a2+a3;
    }
    int get4c(int a1, char a2, short a3, long a4) const {
        return a+a1+a2+a3+a4;
    }
    int get5c(int a1, char a2, short a3, long a4, short a5) const {
        return a+a1+a2+a3+a4+a5;
    }

    Future<int> get0f() {
        return Future<int>(a);
    }
};

void test6(World& world) {
    PROFILE_FUNC;
    ProcessID me = world.rank();
    ProcessID nproc = world.nproc();
    Foo a(world, me*100);

    if (me == 0) {
        print(a.id());
        for (ProcessID p=0; p<nproc; ++p) {
            MADNESS_ASSERT(a.send(p,&Foo::get0).get() == p*100);
            MADNESS_ASSERT(a.task(p,&Foo::get0).get() == p*100);

            MADNESS_ASSERT(a.send(p,&Foo::get0f).get() == p*100);
            MADNESS_ASSERT(a.task(p,&Foo::get0f).get() == p*100);

            MADNESS_ASSERT(a.send(p,&Foo::get1,1).get() == p*100+1);
            MADNESS_ASSERT(a.task(p,&Foo::get1,Future<int>(1)).get() == p*100+1);

            MADNESS_ASSERT(a.send(p,&Foo::get2,1,2).get() == p*100+3);
            MADNESS_ASSERT(a.task(p,&Foo::get2,1,2).get() == p*100+3);

            MADNESS_ASSERT(a.send(p,&Foo::get3,1,2,3).get() == p*100+6);
            MADNESS_ASSERT(a.task(p,&Foo::get3,1,2,3).get() == p*100+6);

            MADNESS_ASSERT(a.send(p,&Foo::get4,1,2,3,4).get() == p*100+10);
            MADNESS_ASSERT(a.task(p,&Foo::get4,1,2,3,4).get() == p*100+10);

            MADNESS_ASSERT(a.send(p,&Foo::get5,1,2,3,4,5).get() == p*100+15);
            MADNESS_ASSERT(a.task(p,&Foo::get5,1,2,3,4,5).get() == p*100+15);

            MADNESS_ASSERT(a.send(p,&Foo::get0c).get() == p*100);
            MADNESS_ASSERT(a.task(p,&Foo::get0c).get() == p*100);

            MADNESS_ASSERT(a.send(p,&Foo::get1c,1).get() == p*100+1);
            MADNESS_ASSERT(a.task(p,&Foo::get1c,1).get() == p*100+1);

            MADNESS_ASSERT(a.send(p,&Foo::get2c,1,2).get() == p*100+3);
            MADNESS_ASSERT(a.task(p,&Foo::get2c,1,2).get() == p*100+3);

            MADNESS_ASSERT(a.send(p,&Foo::get3c,1,2,3).get() == p*100+6);
            MADNESS_ASSERT(a.task(p,&Foo::get3c,1,2,3).get() == p*100+6);

            MADNESS_ASSERT(a.send(p,&Foo::get4c,1,2,3,4).get() == p*100+10);
            MADNESS_ASSERT(a.task(p,&Foo::get4c,1,2,3,4).get() == p*100+10);

            MADNESS_ASSERT(a.send(p,&Foo::get5c,1,2,3,4,5).get() == p*100+15);
            MADNESS_ASSERT(a.task(p,&Foo::get5c,1,2,3,4,5).get() == p*100+15);
        }
    }
    world.gop.fence();


    print("test 6 (world object active message and tasks) seems to be working");
}


class TestFutureForwarding : public WorldObject<TestFutureForwarding> {
public:
    TestFutureForwarding(World& world)
            : WorldObject<TestFutureForwarding>(world) {
        this->process_pending();
    }

    virtual ~TestFutureForwarding() { }

    Future<int> test(int state) {
        if (state < 99) {
            return send(world.random_proc(), &TestFutureForwarding::test, state+1);
        }
        else {
            return Future<int>(state+1);
        }
    }
};

void test6a(World& world) {
    PROFILE_FUNC;

    if (world.size() < 2) return;

    TestFutureForwarding t(world);
    if (world.rank() == 0) {
        Future<int> fred = t.test(0);
        world.gop.fence();
        MADNESS_ASSERT(fred.get() == 100);
    }
    else {
        world.gop.fence();
    }
    if (world.rank() == 0) {
        print("If got here test6a is OK!");
    }
}


void test7(World& world) {
    PROFILE_FUNC;
    int nproc = world.size();
    ProcessID me = world.rank();
    WorldContainer<int,double> c(world);

    typedef WorldContainer<int,double>::iterator iterator;
    typedef WorldContainer<int,double>::const_iterator const_iterator;
    typedef WorldContainer<int,double>::futureT futureT;

    // Everyone inserts distinct values 0..1000 into the container,
    // fences, and then tries to read all values back

    // Note that insertion with key or accessor should be safe
    for (int i=me; i<1000; i+=nproc) c.replace(i,(double) i);
    world.gop.fence();

    for (int i=999; i>=0; --i) {
        futureT fut = c.find(i);
        iterator it = fut.get();
        MADNESS_ASSERT(it != c.end());
        double j = it->second;
        MADNESS_ASSERT(j == i);
    }
    world.gop.fence();

    // Check that unset keys return end correctly
    for (int i=10001; i<10020; ++i) {
        MADNESS_ASSERT(c.find(i).get() == c.end());
    }

    // Check that other iterators compare correctly
    MADNESS_ASSERT(c.find(10).get() == c.find(10).get());
    MADNESS_ASSERT(c.find(11).get() != c.find(12).get());
    MADNESS_ASSERT(c.end() == c.end());
    MADNESS_ASSERT(c.find(12).get() != c.end());

    // Loop thru local stuff
    for (iterator it=c.begin(); it != c.end(); ++it) {
        MADNESS_ASSERT(it->first == it->second);
    };

    // Check shallow copy and const iterator
    const WorldContainer<int,double> d(c);

    // Loop thru local stuff with a const iterator
    for (const_iterator it=d.begin(); it != d.end(); ++it) {
        MADNESS_ASSERT(it->first == it->second);
    };

    world.gop.fence();
    if (me == 0) print("test7 (world container basics) OK");
}

void test8(World& world) {
    PROFILE_FUNC;
    vector<unsigned char> v;
    archive::VectorOutputArchive arout(v);
    arout & &world;

    World* p;
    archive::VectorInputArchive arin(v);
    arin & p;
    MADNESS_ASSERT(p==&world);
    if (world.rank() == 0) print("test8 (serializing world pointer) OK");
}

Void null_func() {
    return None;
}

int val_func() {
    return 1;
}

int val1d_func(int input) {
    return input+1;
}

void test9(World& world) {
    PROFILE_FUNC;
    const int ntask = 100000;

    double used = -cpu_time();
    for (int i=0; i<ntask; ++i) world.taskq.add(null_func);
    used += cpu_time();
    print("Time to add",ntask,"null, local tasks",used,"time/task",used/ntask);

    used = -cpu_time();
    world.taskq.fence();
    used += cpu_time();
    print("Time to run",ntask,"null, local tasks",used,"time/task",used/ntask);

    vector< Future<int> > v = future_vector_factory<int>(ntask);
    used = -cpu_time();
    for (int i=0; i<ntask; ++i) v[i] = world.taskq.add(val_func);
    used += cpu_time();
    print("Time to add",ntask,"value, local tasks",used,"time/task",used/ntask);

    used = -cpu_time();
    print("AAAAAAAAAAAAAAAA0");
    std::cout.flush();
    world.taskq.fence();
    print("AAAAAAAAAAAAAAAA1");
    std::cout.flush();
    used += cpu_time();
    print("Time to run",ntask,"value, local tasks",used,"time/task",used/ntask);
    v.clear();
    print("AAAAAAAAAAAAAAAA");
    std::cout.flush();
    Future<int> input;
    Future<int> result = input;
    used = -cpu_time();
    print("AAAAAAAAAAAAAAAA2");
    std::cout.flush();
    for (int i=0; i<ntask; ++i) {
        result = world.taskq.add(val1d_func,result);
    }
    used += cpu_time();
    print("AAAAAAAAAAAAAAAA3");
    std::cout.flush();
    print("Time to make",ntask,"chain of tasks",used,"time/task",used/ntask);
    input.set(0);
    used = -cpu_time();
    print("AAAAAAAAAAAAAAAA4");
    std::cout.flush();
    world.taskq.fence();
    print("AAAAAAAAAAAAAAAA5");
    std::cout.flush();
    used += cpu_time();
    print("Time to  run",ntask,"chain of tasks",used,"time/task",used/ntask);
    MADNESS_ASSERT(result.get() == ntask);
    if (world.rank() == 0) print("test9 (time task creation and processing) OK");
}


class Mary {
private:
    mutable uint64_t val;
public:
    Mary() : val(0) {}

    Void inc() const {
        val++;
        return None;
    };
    Void add(int i) {
        val += i;
        return None;
    };
    Void fred(int i, double j) {
        val += i*(int)j;
        return None;
    };

    string cary0() {
        return string("Cary0 sends greetings");
    };

    string cary(int i) {
        ostringstream s;
        val += i;
        s << "Cary sends greetings: " << i << " " << val << endl;
        return s.str();
    };

    string alan(int i, int j) {
        ostringstream s;
        val += i*j;
        s << "Alan sends greetings: " << i << " " << j << " " << val << endl;
        return s.str();
    };

    double galahad(const string& str, int j, double z) {
        istringstream s(str);
        int i;
        s >> i;
        //val += i*j*z;
        print("Galahad",str,i,j,z,val);
        return val;
    };


    uint64_t get() const {
        return val;
    };

    bool get_me_twice(World* world, const WorldContainer<int,Mary>& d) {
        return true;
    };

    template <typename Archive>
    void serialize(const Archive& ar) {
        ar & val;
    }
};

Void pounder(const WorldContainer<int,Mary>& m, int ind) {
    for (int i=0; i<1000; ++i)
        m.send(ind, &Mary::inc);
    return None;
}

void test10(World& world) {
    PROFILE_FUNC;
    // test forwarding methods to an item
    ProcessID me = world.rank();
    int nproc = world.size();
    WorldContainer<int,Mary> m(world);
    typedef WorldContainer<int,Mary>::iterator iterator;
    //world.gop.fence();

    for (int i=0; i<nproc; ++i)
        m.send(i,&Mary::inc);
    world.gop.fence();

    for (iterator it=m.begin(); it!=m.end(); ++it) {
        print("mary",it->first,it->second.get());
        MADNESS_ASSERT(int(it->second.get()) == nproc);
    }
    world.gop.fence();

    for (int i=0; i<nproc; ++i)
        m.send(i,&Mary::add,me);
    world.gop.fence();

    for (iterator it=m.begin(); it!=m.end(); ++it) {
        print("mary",it->first,it->second.get());
        MADNESS_ASSERT(long(it->second.get()) == nproc*(nproc+1)/2);
    }
    world.gop.fence();

    for (int i=0; i<nproc; ++i)
        m.send(i,&Mary::fred,2,me);
    world.gop.fence();

    for (iterator it=m.begin(); it!=m.end(); ++it) {
        print("mary",it->first,it->second.get());
        MADNESS_ASSERT(long(it->second.get()) == nproc*(3*nproc-1)/2);
    }
    world.gop.fence();

    // Test that item methods are executed atomically by having
    // everyone pound on one item

    const int ind = 9999999;
    if (world.rank() == 0) m.replace(std::pair<int,Mary>(ind,Mary()));
    world.gop.fence();
    world.taskq.add(pounder, m, ind);
    world.taskq.add(pounder, m, ind);
    world.taskq.add(pounder, m, ind);
    world.taskq.add(pounder, m, ind);
    world.taskq.add(pounder, m, ind);
    world.taskq.add(pounder, m, ind);
    world.taskq.add(pounder, m, ind);
    world.gop.fence();
    if (world.rank() == 0)
        MADNESS_ASSERT(long(m.find(ind).get()->second.get()) == nproc * 1000 * 7);

    world.gop.fence();

    Future<double>  galahad = m.task(ProcessID(0),&Mary::galahad,string("1"),me,3.14);
    world.gop.fence();
    print("result of galahad",galahad.get());


    print("main making vector of results");
    //vector< Future<string> > results(nproc,Future<string>::default_initializer());
    vector< Future<string> > results = future_vector_factory<string>(nproc);
    vector< Future<bool> > b = future_vector_factory<bool>(nproc);
    print("main finished making vector of results");
    for (int i=0; i<nproc; ++i) {
        print("main making task",i);
        results[i] = m.task(i,&Mary::alan,3,4);
        b[i] = m.send(i,&Mary::get_me_twice,&world,m);
        print("main finished making task",i);
    }
    print("about to fence");
    world.gop.fence();

    for (int i=0; i<nproc; ++i) {
        MADNESS_ASSERT(results[i].probe());
        MADNESS_ASSERT(b[i].probe());
        print("results",i,results[i].get(),b[i].get());
    };

    world.gop.fence();

    if (me == 0) print("test10 (messaging to world container items) OK");
}


struct Key {
    typedef unsigned long ulong;
    ulong n, i, j, k;
    hashT hashval;

    Key() {};  // Empty default constructor for speed - but is therefore junk

    Key(ulong n, ulong i, ulong j, ulong k)
            : n(n), i(i), j(j), k(k), hashval(0)
    {
        madness::hash_combine(hashval, n);
        madness::hash_combine(hashval, i);
        madness::hash_combine(hashval, j);
        madness::hash_combine(hashval, k);
    }

    hashT hash() const {
        return hashval;
    }

    template <typename opT>
    void foreach_child(const opT& op) const {
        ulong n2 = n+1;
        ulong i2 = i<<1;
        ulong j2 = j<<1;
        ulong k2 = k<<1;
        for (int p=0; p<2; ++p)
            for (int q=0; q<2; ++q)
                for (int r=0; r<2; ++r)
                    op(Key(n2,i2+p,j2+q,k2+r));
    }

    template <typename Archive>
    void serialize(const Archive& ar) {
        ar & n & i & j & k & hashval;
    }

    bool operator==(const Key& b) const {
        // It's extremely probable that different keys will have a different hash
        return hashval==b.hashval && n==b.n && i==b.i && j==b.j && k==b.k;
    }
};


ostream& operator<<(ostream& s, const Key& key) {
    s << "Key(" << key.n << "," << key.i << "," << key.j << "," << key.k << "," << key.hash() << ")";
    return s;
}

struct Node {
    typedef WorldContainer<Key,Node> dcT;
    Key key;
    double value;
    bool isleaf;
    Node() : value(0.0), isleaf(true) {};
    Node(double value) : value(value), isleaf(true) {};
    Node(const Node& node) : key(node.key), value(node.value), isleaf(node.isleaf) {};

    struct do_random_insert {
        dcT& d;
        double value;
        do_random_insert(dcT& d, double value)
                : d(d), value(value) {};
        void operator()(const Key& key) const {
            d.task(key,&Node::random_insert,d, key, value);
        };
    };

    Void random_insert(const dcT& constd, const Key& keyin, double valin) {
        dcT& d = const_cast<dcT&>(constd);
        //print("inserting",keyin,valin);
        key = keyin;
        value = valin;
        isleaf = true;
        if (value > 0.25 && d.size()<4000) {
            isleaf = false;
            World& world = d.get_world();
            double ran = world.drand();
            key.foreach_child(do_random_insert(d,value*ran));
        }
        return None;
    };

    template <class Archive>
    void serialize(Archive& ar) {
        ar & key & value & isleaf;
    }

    bool is_leaf() const {
        return isleaf;
    };

    double get() const {
        return value;
    };

    Void set(double v) {
        value = v;
        return None;
    };
};

ostream& operator<<(ostream& s, const Node& node) {
    s << "Node(" << node.get() << "," << node.is_leaf() << ")" << endl;
    return s;
}


void walker1(WorldContainer<Key,Node>& d, const Key& key);

struct Walker1 {
    WorldContainer<Key,Node>& d;
    Walker1(WorldContainer<Key,Node>& d) : d(d) {};
    void operator()(const Key& key) const {
        walker1(d,key);
    };
};


void walker1(WorldContainer<Key,Node>& d, const Key& key) {
    static double counter = 0;
    WorldContainer<Key,Node>::iterator it = d.find(key).get();
    if (it != d.end()) {
        Node node = it->second;
        node.set(++counter);
        d.erase(key);
        d.replace(key,node);
        it = d.find(key).get();
        MADNESS_ASSERT(it != d.end());
        MADNESS_ASSERT(it->second.get() == counter);
        if (!node.is_leaf()) {
            key.foreach_child(Walker1(d));
        }
    }
}

void walker2(WorldContainer<Key,Node>& d, const Key& key) {
    static double counter = 1;
    WorldContainer<Key,Node>::iterator it = d.find(key).get();
    if (it != d.end()) {
        Node node = it->second;
        node.set(++counter);
        d.replace(key,node);
        it = d.find(key).get();
        MADNESS_ASSERT(it != d.end());
        if (it->second.get() != counter) {
            print("failing",it->second.get(),counter,key,d.owner(key));
        }
        MADNESS_ASSERT(it->second.get() == counter);
        if (!node.is_leaf()) {
            key.foreach_child(Walker1(d));
        }
    }
}

void test11(World& world) {
    PROFILE_FUNC;
    // Test the various flavours of erase
    ProcessID me = world.rank();
    WorldContainer<Key,Node> d(world);

    // First build an oct-tree with random depth
    world.srand();
    print("first ran#",world.drand());
    world.gop.fence();
    if (me == 0) {
        Key root = Key(0,0,0,0);
        d.task(root,&Node::random_insert,d,root,1.0);
    }
    world.gop.fence();

    print("size before erasing",d.size());
    //d.clear();
    d.erase(d.begin(),d.end());
    print("size after erasing",d.size());
    world.srand();
    print("first ran#",world.drand());
    world.gop.fence();
    // rebuild the tree in the same container
    if (me == 0) {
        Key root = Key(0,0,0,0);
        d.task(root,&Node::random_insert,d,root,1.0);
    }
    world.gop.fence();
    print("size after rebuilding",d.size());

    // Test get, erase, and re-insert of nodes with new value by node 0
    if (me == 0) {
        Key root = Key(0,0,0,0);
        walker1(d,root);
    }
    world.gop.fence();
    print("walker1 done");


    // Test get and re-insert of nodes with new value by node 0
    if (me == 0) {
        Key root = Key(0,0,0,0);
        walker2(d,root);
    }
    world.gop.fence();
    print("walker2 done");

    print("size before clearing",d.size());
    d.clear();
    print("size after clearing",d.size());
    if (me == 0) print("test11 (erasing and inserting in world containers) OK");
}


void test12(World& world) {
    PROFILE_FUNC;
    if (world.size() != 1) return;
    // Test file IO
    ProcessID me = world.rank();
    WorldContainer<int,double> d(world);

    // Everyone puts 100 distinct entries in the container
    for (int i=0; i<100; ++i) d.replace(me*100 + i, me*100+i);

    world.gop.fence();

    archive::BinaryFstreamOutputArchive out("testme.ar");
    out & d;
    out.close();

    world.gop.fence();

    archive::BinaryFstreamInputArchive in("testme.ar");
    WorldContainer<int,double> c(world);
    in & c;

    world.gop.fence();

    for (int i=0; i<100; ++i) {
        int key = me*100+i;
        MADNESS_ASSERT(c.probe(key));
        MADNESS_ASSERT(c.find(key).get()->second == key);
    }

    world.gop.fence();

    if (world.rank() == 0) print("test12 (container archive I/O) OK");
}

void test13(World& world) {
    PROFILE_FUNC;
    // Basic functionality with 1 (default) writer
    archive::ParallelOutputArchive fout(world, "fred");
    fout & 1.0 & "hello";
    fout.close();

    double v;
    char s[6];
    archive::ParallelInputArchive fin(world, "fred");
    fin & v & s;
    fin.close();
    fin.remove();

    print("This is what I read", v, s);
    world.gop.fence();


    // Store and load an archive with multiple writers

    int nio = (world.size()-1)/2 + 1;
    if (nio > 10) nio = 10;

    print("nio",nio);

    ProcessID me = world.rank();
    WorldContainer<int,double> d(world);
    // Everyone puts 100 distinct entries in the container
    for (int i=0; i<100; ++i) {
        int key = me*100+i;
        d.replace(key, double(key));
    }

    world.gop.fence();

    fout.open(world,"fred",nio);
    fout & d;
    fout.close();

    WorldContainer<int,double> c(world);
    fin.open(world,"fred");
    fin & c;

    for (int i=0; i<100; ++i) {
        int key = me*100+i;
        MADNESS_ASSERT(c.find(key).get()->second == key);
    }

    fin.close();
    archive::ParallelOutputArchive::remove(world, "fred");

    print("Test13 OK");
    world.gop.fence();
}

inline bool is_odd(int i) {
    return i & 0x1;
}

inline bool is_even(int i) {
    return !is_odd(i);
}

void work_odd(World& world) {
    test5(world);
    test6(world);
    test6a(world);
    test7(world);
    test8(world);
    test9(world);
    test10(world);
    test11(world);
    // test12(world); cannot run due to filename collision
    // test13(world);
    world.gop.fence();
}

void work_even(World& world) {
    test5(world);
    test6(world);
    test6a(world);
    test7(world);
    test8(world);
    test9(world);
    test10(world);
    test11(world);
    // test12(world); cannot run due to filename collision
    // test13(world);
    world.gop.fence();
}

void test_multi_world(World& world) {
    if (world.size() < 2) return;

    // Make two more worlds: odd processes, even processes

    // a) make list of ranks of processes in the subgroups
    //
    // Only process belonging to the subgroups participate
    // in the next steps
    //
    // b) make MPI group and hence new MPI sub-communcator
    //
    // c) make new worlds and do work

    std::cout << "\n\nREPEATING TESTS IN MULTI-WORLD\n\n" << std::endl;

    std::vector<int> odd, even;
    for (int i=0; i<world.size(); ++i) {
        if (is_odd(i))
            odd.push_back(i);
        else
            even.push_back(i);
    }

    if (world.rank() & 0x1) {   // Odd processes
        MPI::Group g_odd = world.mpi.comm().Get_group().Incl(odd.size(), &odd[0]);
        MPI::Intracomm comm_odd = world.mpi.comm().Create(g_odd);
        {
            World world_odd(comm_odd);
            work_odd(world_odd);
        }
        comm_odd.Free();

    }
    else {                      // Even processes
        MPI::Group g_even = world.mpi.comm().Get_group().Incl(even.size(),&even[0]);
        MPI::Intracomm comm_even = world.mpi.comm().Create(g_even);
        {
            World world_even(comm_even);
            work_even(world_even);
        }
        comm_even.Free();
    }

    world.gop.fence();
}


int main(int argc, char** argv) {
    initialize(argc,argv);

    World world(MPI::COMM_WORLD);

    redirectio(world);
    print("The processor frequency is",cpu_frequency());
    print("There are",world.size(),"processes and I am process",world.rank(),"with",ThreadPool::size(),"threads");

    world.args(argc,argv);

    world.gop.fence();

    try {
        PROFILE_BLOCK(main_program);

        test0(world);
        // ??????  When/why did these tests get deleted?
//       if (world.nproc() > 1) {
//         test1(world);
//         test2(world);
//         test3(world);
//       }
//       test4(world);
//       test4a(world);
        test5(world);
        test6(world);
        test6a(world);
        test7(world);
        test8(world);
        test9(world);
        test10(world);
        //test11(world);
        test12(world);
        test13(world);

        // for (int i=0; i<100; ++i) {
        //     print("REPETITION",i);
        //     test_multi_world(world);
        // }
    }
    catch (MPI::Exception e) {
        error("caught an MPI exception");
    }
    catch (madness::MadnessException e) {
        print(e);
        error("caught a MADNESS exception");
    }
    catch (const char* s) {
        print(s);
        error("caught a string exception");
    }
    catch (...) {
        error("caught unhandled exception");
    }

    print("entering final fence");
    world.gop.fence();
    print("done with final fence");

    //ThreadPool::end();
    //print_stats(world);
    finalize();
    return 0;
}
