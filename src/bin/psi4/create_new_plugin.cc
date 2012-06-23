#include <boost/filesystem.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/regex.hpp>
#include <boost/xpressive/xpressive.hpp>
#include <boost/xpressive/regex_actions.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>

#include <libmints/basisset.h>

#include <cstdio>
#include <cstdlib>
#include <string>
#include <psi4-dec.h>
#include <psiconfig.h>

using namespace std;
using namespace psi;
using namespace boost;

namespace {
std::string make_filename(const std::string& name)
{
    // Modify the name of the basis set to generate a filename: STO-3G -> sto-3g
    string filename = name;

    string format_underscore("_"); // empty string
    // Replace all '(' with '_'
    xpressive::sregex match_format = xpressive::as_xpr("(");
    filename = regex_replace(filename, match_format, format_underscore);

    // Replace all ')' with '_'
    match_format = xpressive::as_xpr(")");
    filename = regex_replace(filename, match_format, format_underscore);

    // Replace all ',' with '_'
    match_format = xpressive::as_xpr(",");
    filename = regex_replace(filename, match_format, format_underscore);

    // Replace all '*' with 's'
    match_format = xpressive::as_xpr("*");
    string format_star("s");
    filename = regex_replace(filename, match_format, format_star);

    // Replace all '+' with 'p'
    match_format = xpressive::as_xpr("+");
    string format_plus("p");
    filename = regex_replace(filename, match_format, format_plus);

    // Replace all '-' with '_'
    match_format = xpressive::as_xpr("-");
    string format_hyphen("_");
    filename = regex_replace(filename, match_format, format_hyphen);

    return filename;
}
}

namespace psi {

/**
 * 
 */
class PluginFileManager{
  protected:
    std::string plugin_name_;
    std::vector<std::pair<std::string, std::string> > files_;
  public:
    PluginFileManager(const std::string &plugin_name):
       plugin_name_(plugin_name)
    {


    }

    /*
     * Adds a file to be copied over from the psi4/lib/plugin directory to the target
     * @param source_name: The name of the file as is appears in the psi4/lib/plugin directory
     * @param target_name: The name of the file as it will appear in the new directory.  If omitted,
     * defaults to the same name as provided for source_name.
     */
    void add_file(const std::string &source_name, const std::string &target_name = "")
    {
        if(target_name == "")
            files_.push_back(std::make_pair(source_name, source_name));
        else
            files_.push_back(std::make_pair(source_name, target_name));
    }

    void process()
    {
        // The location of the plugin templates, in the Psi4 source
        std::string psiDataDirName = Process::environment("PSIDATADIR");
        std::string psiDataDirWithPlugin = psiDataDirName + "/plugin";

        boost::filesystem::path bf_path;
        bf_path = boost::filesystem::system_complete(psiDataDirWithPlugin);
        if(!boost::filesystem::is_directory(bf_path)) {
            printf("Unable to read the PSI4 plugin folder - check the PSIDATADIR environmental variable\n"
                    "      Current value of PSIDATADIR is %s\n", psiDataDirName.c_str());
            exit(1);
        }

        // Make a faux camel-case of the name
        string Name = plugin_name_;
        Name[0] = ::toupper(Name[0]);

        // Formatted strings, to be substituted in later
        std::string format_top_srcdir(PSI_TOP_SRCDIR);
        std::string format_top_objdir(PSI_TOP_OBJDIR);
        std::string format_plugin(plugin_name_);
        std::string format_PLUGIN = boost::algorithm::to_upper_copy(plugin_name_);

        std::vector<std::pair<std::string, std::string> >::const_iterator iter;
        for(iter = files_.begin(); iter != files_.end(); ++iter){
            std::string source_name = psiDataDirWithPlugin + "/" + iter->first;
            std::string target_name   = plugin_name_ + "/" + iter->second;

            // Load in Makefile.template
            FILE* fp = fopen(source_name.c_str(), "r");
            if (fp == NULL) {
                printf("create_new_plugin: Unable to open Makefile template.\n");
                exit(1);
            }
            // Stupid way to read in entire file.
            char line[256];
            std::stringstream file;
            while(fgets(line, sizeof(line), fp))
                file << line;
            std::string filestring = file.str();
            fclose(fp);

            // Search and replace placeholders in the string
            boost::xpressive::sregex match_format = xpressive::as_xpr("@top_srcdir@");
            filestring = xpressive::regex_replace(filestring, match_format, format_top_srcdir);
            match_format = boost::xpressive::as_xpr("@top_objdir@");
            filestring = xpressive::regex_replace(filestring, match_format, format_top_objdir);
            match_format = boost::xpressive::as_xpr("@plugin@");
            filestring = xpressive::regex_replace(filestring, match_format, format_plugin);
            match_format = boost::xpressive::as_xpr("@Plugin@");
            filestring = xpressive::regex_replace(filestring, match_format, Name);
            match_format = boost::xpressive::as_xpr("@PLUGIN@");

            // Write the new file out
            fp = fopen(target_name.c_str(), "w");
            if (fp == 0) {
                boost::filesystem::remove_all(plugin_name_);
                printf("Unable to create %s\n", target_name.c_str());
                exit(1);
            }
            fputs(filestring.c_str(), fp);
            fclose(fp);

            printf("\tCreated: %s\n", iter->second.c_str());
        }


    }
};

void create_new_plugin(std::string name, const std::string& template_name)
{
    std::string template_name_lower(template_name);
    // First make it lower case
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    transform(template_name_lower.begin(), template_name_lower.end(), template_name_lower.begin(), ::tolower);

    // Start == check to make sure the plugin name is valid
    string plugin_name = make_filename(name);
    smatch results;
    regex check_name("^[A-Za-z].*");
    if (!regex_match(plugin_name, results, check_name)) {
        printf("Plugin name must begin with a letter.\n");
        exit(1);
    }
    // End == check to make sure the plugin name is valid


    if(template_name_lower.empty())
        template_name_lower = "plugin";

    // Make a directory with the name plugin_name
    if (!boost::filesystem::create_directory(plugin_name)) {
        printf("Plugin directory %s already exists.\n", plugin_name.c_str());
        exit(1);
    }
    printf("Created new plugin directory, %s, using '%s' template.\n", plugin_name.c_str(),  template_name_lower.c_str());

    // Process the files
    PluginFileManager file_manager(plugin_name);
    file_manager.add_file("/Makefile.template", "Makefile");
    file_manager.add_file("/input.dat.template", "input.dat");
    file_manager.add_file("pymodule.py.template", "pymodule.py");
    file_manager.add_file("__init__.py.template", "__init__.py");
    file_manager.add_file("inputalt.dat.template", "inputalt.dat");
    file_manager.add_file("doc.rst.template", "doc.rst");
    file_manager.add_file(template_name_lower + ".cc.template", name + ".cc");
    if(template_name_lower == "scf"){
        // The SCF file has multiple files
        file_manager.add_file("scf.scf.h.template", "scf.h");
        file_manager.add_file("scf.scf.cc.template", "scf.cc");
        file_manager.add_file("scf.cc.template", name + ".cc");
        // Overwrite the existing input file with a more appropriate one
        file_manager.add_file("scf.input.dat.template", "input.dat");
    }
    file_manager.process();

}

}
