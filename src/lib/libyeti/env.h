#ifndef _tiled_tensor_Env_h
#define _tiled_tensor_Env_h

#include <iostream>

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

namespace yeti {

class IndentTracker {

    private:
        int indentsize_;

        int n_;

    public:
        IndentTracker();

        IndentTracker
        operator++();

        IndentTracker
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

 public:
   static IndentTracker& indent;

   static void init(int me, const char *outFileName = "");
   static void init(int me, std::ostream *stream);
   static void free();
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

#ifdef redefine_size_t
#undef size_t
#endif

#endif // Header Guard
