#include <psiconfig.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>

namespace psi {

    // Useful typedef's
    typedef PsiReturnType (*plugin_t)(Options &);
    typedef int (*read_options_t)(std::string, Options&);
    typedef void (*init_plugin_t)(const shared_ptr<Communicator>& comm, const Process::Environment& env);

    // The following lines are used in plugins
    extern "C" void init_plugin(const shared_ptr<psi::Communicator>& comm, const psi::Process::Environment& env);
    #define INIT_PLUGIN psi::init_plugin_t init_plugin_p = &psi::init_plugin;

    // The following lines are used by the PSI4 driver.
    typedef struct plugin_info_struct {
        // Uppercased plugin name.
        std::string name;

        // File handle to plugin
        void *plugin_handle;

        // Used to call the user provided plugin code
        plugin_t plugin;

        // Function pointer to the plugin's read_options function
        read_options_t read_options;

        // Function pointer to the plugin's init_plugin function found in libplugin
        init_plugin_t init_plugin;
    } plugin_info;

    /**
      Loads a plugin from disk and calls its initialization routine.

      @param plugin_path Absolute or relative path and filename name to the plugin.
      @returns A plugin_info structure containing function pointers to entry points in the plugin.
    */
    plugin_info plugin_load(std::string& fullpathname);

    /**
      Frees the plugin from memory.

      @param plugin A plugin_info struct that contains the plugin to free.
    */
    void plugin_close(const plugin_info& plugin);
}

