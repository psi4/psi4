#ifndef _psi_src_lib_libscript_script_h_
#define _psi_src_lib_libscript_script_h_

namespace psi {
    
    /* 
     * Hopefully I can keep this abstract enough to handle multiple
     * scripting languages.
     */
    class Script {
    public:
        Script();
        virtual ~Script();
        
        /** Starts up the scripting language environment and adds
         * PSI4 specific functionality to the environment.
         */
        virtual void initialize() = 0;
        
        /** Shuts down the scripting language and frees any memory
         * that may have been allocated by the scripting language.
         */
        virtual void finalize() = 0;
    };

    class Python : public Script {
    public:
        Python();
        virtual ~Python();
        
        virtual void initialize();
        virtual void finalize();
    };
}

#endif /* end of include guard: _psi_src_lib_libscript_script_h_ */
