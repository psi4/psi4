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

namespace psi {

std::string make_filename(const std::string& name)
{
    // Modify the name of the basis set to generate a filename: STO-3G -> sto-3g
    string filename = name;

    // First make it lower case
    transform(filename.begin(), filename.end(), filename.begin(), ::tolower);

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

void create_new_plugin(const std::string& name, const std::string& templ)
{
    // Start == check to make sure the plugin name is valid
    string plugin_name = make_filename(name);
    smatch results;
    regex check_name("^[A-Za-z].*");
    if (!regex_match(plugin_name, results, check_name)) {
        printf("Plugin name must begin with a letter.\n");
        exit(1);
    }
    // End == check to make sure the plugin name is valid

    FILE* fp = 0;

    string template_name = templ;
    if (template_name.empty())
        template_name = "plugin";

    std::string psiDataDirName = Process::environment("PSIDATADIR");
    std::string psiDataDirWithPlugin = psiDataDirName + "/plugin";
    std::string fileMakefile = psiDataDirWithPlugin + "/Makefile.template";
    std::string fileInput    = psiDataDirWithPlugin + "/input.dat.template";
    std::string fileSource   = psiDataDirWithPlugin + "/" + template_name + ".cc.template";

    boost::filesystem::path bf_path;
    bf_path = boost::filesystem::system_complete(psiDataDirWithPlugin);
    if(!boost::filesystem::is_directory(bf_path)) {
        printf("Unable to read the PSI4 plugin folder - check the PSIDATADIR environmental variable\n"
                "      Current value of PSIDATADIR is %s\n", psiDataDirName.c_str());
        exit(1);
    }

    // Load in Makefile.template
    fp = fopen(fileMakefile.c_str(), "r");
    if (fp == NULL) {
        printf("create_new_plugin: Unable to open Makefile template.\n");
        exit(1);
    }
    // Stupid way to read in entire file.
    char line[256];
    std::stringstream file;
    while(fgets(line, sizeof(line), fp)) {
        file << line;
    }
    fclose(fp);
    std::string makefile = file.str();

    // Load in input.dat.template
    fp = fopen(fileInput.c_str(), "r");
    if (fp == NULL) {
        printf("create_new_plugin: Unable to open input.dat template.\n");
        exit(1);
    }

    // Stupid way to read in entire file.
    std::stringstream file2;
    while(fgets(line, sizeof(line), fp)) {
        file2 << line;
    }
    fclose(fp);
    std::string input = file2.str();

    // Load in plugin.cc.template
    fp = fopen(fileSource.c_str(), "r");
    if (fp == NULL) {
        printf("create_new_plugin: Unable to open %s.cc template.\n", template_name.c_str());
        exit(1);
    }

    // Stupid way to read in entire file.
    std::stringstream file3;
    while(fgets(line, sizeof(line), fp)) {
        file3 << line;
    }
    fclose(fp);
    std::string source = file3.str();

    // Search and replace tags

    std::string format_top_srcdir(PSI_TOP_SRCDIR);
    std::string format_top_objdir(PSI_TOP_OBJDIR);
    std::string format_plugin(plugin_name);
    std::string format_PLUGIN = boost::algorithm::to_upper_copy(plugin_name);

    // Replace all '@top_srcdir@' with the top source directory
    boost::xpressive::sregex match_format = xpressive::as_xpr("@top_srcdir@");
    makefile = xpressive::regex_replace(makefile, match_format, format_top_srcdir);
    input    = xpressive::regex_replace(input,    match_format, format_top_srcdir);
    source   = xpressive::regex_replace(source,   match_format, format_top_srcdir);

    match_format = boost::xpressive::as_xpr("@top_objdir@");
    makefile = xpressive::regex_replace(makefile, match_format, format_top_objdir);
    input    = xpressive::regex_replace(input,    match_format, format_top_objdir);
    source   = xpressive::regex_replace(source,   match_format, format_top_objdir);

    match_format = boost::xpressive::as_xpr("@plugin@");
    makefile = xpressive::regex_replace(makefile, match_format, format_plugin);
    input    = xpressive::regex_replace(input,    match_format, format_plugin);
    source   = xpressive::regex_replace(source,   match_format, format_plugin);

    match_format = boost::xpressive::as_xpr("@PLUGIN@");
    makefile = xpressive::regex_replace(makefile, match_format, format_PLUGIN);
    input    = xpressive::regex_replace(input,    match_format, format_PLUGIN);
    source   = xpressive::regex_replace(source,   match_format, format_PLUGIN);

    // Make a directory with the name plugin_name
    if (!boost::filesystem::create_directory(plugin_name)) {
        printf("Plugin directory already exists.\n");
        exit(1);
    }

    fp = fopen((plugin_name + "/Makefile").c_str(), "w");
    if (fp == 0) {
        boost::filesystem::remove_all(plugin_name);
        printf("Unable to create new Makefile.\n");
        exit(1);
    }
    fputs(makefile.c_str(), fp);
    fclose(fp);

    fp = fopen((plugin_name + "/input.dat").c_str(), "w");
    if (fp == 0) {
        boost::filesystem::remove_all(plugin_name);
        printf("Unable to create new input.dat.\n");
        exit(1);
    }
    fputs(input.c_str(), fp);
    fclose(fp);

    fp = fopen((plugin_name + "/" + plugin_name + ".cc").c_str(), "w");
    if (fp == 0) {
        boost::filesystem::remove_all(plugin_name);
        printf("Unable to create new %s.cc\n", plugin_name.c_str());
        exit(1);
    }
    fputs(source.c_str(), fp);
    fclose(fp);

    printf("Created new plugin directory, %s, using '%s' template.\n", plugin_name.c_str(),  template_name.c_str());
    printf("\tcreated: Makefile\n\tcreated: input.dat\n\tcreated: %s.cc\n", plugin_name.c_str());
}

}
