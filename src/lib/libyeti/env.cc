#include "env.h"
#include <string.h>
#include <iostream>
#include <fstream>

using namespace yeti;
using namespace std;

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

#define ENV_INDENT_SIZE 4

static IndentTracker _env_indent_tracker;

int Env::initialized_ = 0;
int Env::me_= -1;
bool Env::out0Allocated_ = false;
bool Env::outnAllocated_ = false;
ostream *Env::out0_;
ostream *Env::outn_;
IndentTracker& Env::indent = _env_indent_tracker;

// TODO figure out a better method that will free memory automagically
std::ostream&
yeti::operator<<(std::ostream& stream, const IndentTracker& indent)
{
    indent.add_indent(stream);
    return stream;
}

IndentTracker::IndentTracker()
    : n_(0), indentsize_(ENV_INDENT_SIZE)
{
}


void
IndentTracker::add_indent(std::ostream& os) const
{
    if (!n_)
        return;

    for (int i=0; i < n_ * indentsize_; ++i)
        os << " ";
}


/**
 * Set up the output stream using a filename, which should be
 * omitted if std::cout is to be used.
 */
void
Env::init(int me, const char *outFileName){
    if(initialized_) return;
    me_ = me;
    if(strcmp("", outFileName)){
        out0_ = new ofstream(outFileName, ios_base::app | ios_base::out );
        out0Allocated_ = true;
    }else{
        out0_ = &cout;
    }
    if(me){
        outn_ = new ofstream("/dev/null");
        outnAllocated_ = true;
    }else{
        outn_ = out0_;
    }
    initialized_ = 1;
}

/**
 * Set up the output stream using an exising stream.
 */
void
Env::init(int me, ostream *stream){
    if(initialized_) return;
    me_ = me;
    outn_ = stream;
    if(me){
        out0_ = new ofstream("/dev/null");
        outnAllocated_ = true;
    }else{
        out0_ = outn_;
    }
    initialized_ = 1;
}

void
IndentTracker::operator++()
{
    ++n_;
}

void
IndentTracker::operator--()
{
    --n_;
}

/**
 * Free the memory associated with Env
 */
void
Env::free()
{
    if(initialized_){
        // We have to make sure that we allocated these before deleting
        if(out0Allocated_){
            delete out0_;
            out0_ = 0;
        }
        if(outnAllocated_){
            delete outn_;
            outn_ = 0;
        }
        initialized_ = false;
    }
}

