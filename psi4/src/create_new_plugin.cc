/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include <cstdio>
#include <cstdlib>
#include <string>
#include <regex>
#include <sys/stat.h>

#include "psi4/psi4-dec.h"
#include "psi4/libfilesystem/path.h"
#include "psi4/libpsi4util/libpsi4util.h"

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

namespace {
std::string make_filename(const std::string &name)
{
    // Modify the name of the basis set to generate a filename: STO-3G -> sto-3g
    std::string filename = name;

    // Replace ( ) , - with _
    filename = std::regex_replace(filename, std::regex("\\(|\\)|,|\\-"), "_");

    // Replace * with s
    filename = std::regex_replace(filename, std::regex("\\*"), "s");

    // Replace + with p
    filename = std::regex_replace(filename, std::regex("\\+"), "p");

    return filename;
}
}

namespace psi {

/**
 *
 */
class PluginFileManager
{
protected:
    std::string plugin_name_;
    bool cd_into_directory_;
    std::vector<std::pair<std::string, std::string> > files_;
    std::vector<std::string> source_files_;
public:
    PluginFileManager(const std::string &plugin_name, bool cd_into_directory = true) :
            plugin_name_(plugin_name), cd_into_directory_(cd_into_directory)
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
        if (target_name == "")
            files_.push_back(std::make_pair(source_name, source_name));
        else
            files_.push_back(std::make_pair(source_name, target_name));

        std::string ext(filesystem::path(target_name).extension());
        if (ext == "h" || ext == "cc")
            source_files_.push_back(target_name);
    }

    void process()
    {
        // The location of the plugin templates, in the Psi4 source
        std::string psiDataDirName = Process::environment("PSIDATADIR");
        std::string psiDataDirWithPlugin = psiDataDirName + "/plugin";

        std::string fpath = filesystem::path(psiDataDirWithPlugin).make_absolute().str();
        struct stat sb;
        if (::stat(fpath.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode) == false) {
            printf("Unable to read the Psi4 plugin folder - check the PSIDATADIR environmental variable\n"
                           "      Current value of PSIDATADIR is %s\n", psiDataDirName.c_str());
            exit(1);
        }

        // Make a faux camel-case of the name
        std::string Name = plugin_name_;
        Name[0] = ::toupper(Name[0]);

        // Formatted strings, to be substituted in later
        std::ostringstream imploded;
        std::copy(source_files_.begin(), source_files_.end(),
                  std::ostream_iterator<std::string>(imploded, " "));
        std::string format_source_list = imploded.str();
        std::string format_plugin(plugin_name_);
        std::string format_PLUGIN = plugin_name_;
        std::transform(format_PLUGIN.begin(), format_PLUGIN.end(), format_PLUGIN.begin(), ::toupper);
        std::string format_ldflags(TOSTRING(PLUGIN_LDFLAGS));

        trim_spaces(format_source_list);

        std::vector<std::pair<std::string, std::string> >::const_iterator iter;
        for (iter = files_.begin(); iter != files_.end(); ++iter) {
            std::string source_name = psiDataDirWithPlugin + "/" + iter->first;
            std::string target_name = cd_into_directory_ ? (plugin_name_ + "/" + iter->second) : iter->second;

            // Load in Makefile.template
            FILE *fp = fopen(source_name.c_str(), "r");
            if (fp == NULL) {
                printf("create_new_plugin: Unable to open %s template.\n", source_name.c_str());
                exit(1);
            }
            // Stupid way to read in entire file.
            char line[256];
            std::stringstream file;
            while (fgets(line, sizeof(line), fp))
                file << line;
            std::string filestring = file.str();
            fclose(fp);

            filestring = std::regex_replace(filestring, std::regex("@plugin@"), format_plugin);
            filestring = std::regex_replace(filestring, std::regex("@Plugin@"), Name);
            filestring = std::regex_replace(filestring, std::regex("@PLUGIN@"), format_PLUGIN);
            filestring = std::regex_replace(filestring, std::regex("@sources@"), format_source_list);

            // Write the new file out
            fp = fopen(target_name.c_str(), "w");
            if (fp == 0) {
                // boost::filesystem::remove_all(plugin_name_);
                printf("Unable to create %s\n", target_name.c_str());
                exit(1);
            }
            fputs(filestring.c_str(), fp);
            fclose(fp);

            printf("\tCreated: %s\n", iter->second.c_str());
        }


    }
};

void create_new_plugin(std::string name, const std::string &template_name)
{
    std::string template_name_lower(template_name);
    // First make it lower case
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    transform(template_name_lower.begin(), template_name_lower.end(), template_name_lower.begin(), ::tolower);

    // Start == check to make sure the plugin name is valid
    std::string plugin_name = make_filename(name);
    if (!std::isalpha(plugin_name[0])) {
        printf("Plugin name must begin with a letter.\n");
        exit(1);
    }
    // End == check to make sure the plugin name is valid


    if (template_name_lower.empty())
        template_name_lower = "plugin";

    // Make a directory with the name plugin_name
    if (!filesystem::create_directory(plugin_name)) {
        printf("Plugin directory %s already exists.\n", plugin_name.c_str());
        exit(1);
    }
    printf("Created new plugin directory, %s, using '%s' template.\n", plugin_name.c_str(), template_name_lower.c_str());

    // Process the files
    PluginFileManager file_manager(plugin_name);
    file_manager.add_file("CMakeLists.txt.template", "CMakeLists.txt");
    file_manager.add_file("input.dat.template", "input.dat");
    file_manager.add_file("pymodule.py.template", "pymodule.py");
    file_manager.add_file("__init__.py.template", "__init__.py");
    file_manager.add_file("doc.rst.template", "doc.rst");
    file_manager.add_file(template_name_lower + ".cc.template", name + ".cc");
    if (template_name_lower == "scf") {
        // The SCF file has multiple files
        file_manager.add_file("scf.scf.h.template", "scf.h");
        file_manager.add_file("scf.scf.cc.template", "scf.cc");
        // Overwrite the existing pymodule file with a more appropriate one
        file_manager.add_file("scf.pymodule.py.template", "pymodule.py");
    }
    if (template_name_lower == "ambit") {
        file_manager.add_file("ambit.input.dat.template", "input.dat");
    }
    file_manager.process();
}

void create_new_plugin_makefile()
{
    printf("Creating new plugin Makefile in the current directory.\n");

    filesystem::path cwd = filesystem::path::getcwd();
    std::string name = make_filename(cwd.stem());
    PluginFileManager file_manager(name, false);
    file_manager.add_file("CMakeLists.txt.template", "CMakeLists.txt");
    file_manager.process();
}

}
