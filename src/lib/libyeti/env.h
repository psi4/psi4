#ifndef _tiled_tensor_Env_h
#define _tiled_tensor_Env_h

#include <iostream>
#include <fstream>


namespace yeti {

class IndentTracker {

    private:
        int indentsize_;

        int n_;

    public:
        IndentTracker();

        void
        operator++();

        void
        operator--();

        void
        add_indent(std::ostream& os) const;

};



/**
 * The Env class determines how I/O should be performed
 * based on the current node number and the desired
 * locality of the print statement
 */
class Env {

 protected:
   /// Whether the Env object has been initialized
   static int initialized_;
   /// The current node number
   static int me_;
   /// The stream object for output from node 0 only
   static std::ostream *out0_;
   /// The stream object for output from all nodes
   static std::ostream *outn_;
   /// Whether the out0 stream was allocated, or assigned
   static bool out0Allocated_;
   /// Whether the outn stream was allocated, or assigned
   static bool outnAllocated_;

 public:
   static IndentTracker& indent;

   static void init(int me, const char *outFileName = "")
{
    if(initialized_) return;
    me_ = me;
    if(strcmp("", outFileName)){
        out0_ = new std::ofstream(outFileName, std::ios_base::app | std::ios_base::out );
        out0Allocated_ = true;
    }else{
        out0_ = &std::cout;
    }
    if(me){
        outn_ = new std::ofstream("/dev/null");
        outnAllocated_ = true;
    }else{
        outn_ = out0_;
    }
    initialized_ = 1;
}

   static void init(int me, std::ostream *stream)
{
    if(initialized_) return;
    me_ = me;
    out0_ = stream;
    if(me){
        outn_ = new std::ofstream("/dev/null");
        outnAllocated_ = true;
    }else{
        outn_ = out0_;
    }
    initialized_ = 1;
}
   static void free()
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
   /// Return nonzero if Env has been initialized.
   static int initialized() { return initialized_; }
   /// Return an ostream that writes from all nodes.
   static std::ostream &outn() { return *outn_; }
   /// Return an ostream for error messages that writes from all nodes.
   static std::ostream &errn() { return *outn_; }
   /// Return an ostream that writes from node 0.
   static std::ostream &out0() { return *out0_; }
   /// Return an ostream for error messages that writes from node 0.
   static std::ostream &err0() { return *out0_; }
   
};

std::ostream&
operator<<(std::ostream& os, const IndentTracker& indent);


}

#endif // Header Guard
