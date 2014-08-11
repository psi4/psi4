/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

#include "PluginMaker.h"
#include "libparallel/ParallelPrinter.h"
#include "exception.h"

namespace psi {
typedef vector<pair<string,string> > FileVector;

boost::filesystem::path PluginMaker::PluginDir(){
   string psiDataDirName=Process::environment("PSIDATADIR");
   string psiDataDirWithPlugin=psiDataDirName+"/plugin";
   return boost::filesystem::system_complete(psiDataDirWithPlugin);
}


void GMakePlugin::process() {
   // The location of the plugin templates, in the Psi4 source
   string psiDataDirName=Process::environment("PSIDATADIR");
   string psiDataDirWithPlugin=psiDataDirName+"/plugin";

   boost::filesystem::path bf_path=PluginDir();
   if (!boost::filesystem::is_directory(bf_path)) {
      string error="Unable to read the PSI4 plugin folder - check the "
            "PSIDATADIR environmental variable\n"
            "      Current value of PSIDATADIR is "+psiDataDirName;
      throw PSIEXCEPTION(error.c_str());
   }

   // Make a faux camel-case of the name
   string Name=plugin_name_;
   Name[0]=::toupper(Name[0]);

   // Formatted strings, to be substituted in later
   string format_top_srcdir(PSI_TOP_SRCDIR);
   string format_top_objdir(PSI_TOP_OBJDIR);
   string format_plugin(plugin_name_);
   string format_PLUGIN=boost::algorithm::to_upper_copy(plugin_name_);

   FileVector::const_iterator iter;
   for (iter=files_.begin(); iter!=files_.end(); ++iter) {
      string source_name=psiDataDirWithPlugin+"/"+iter->first;
      string target_name=plugin_name_+"/"+iter->second;

      // Load in Makefile.template
      FILE* fp=fopen(source_name.c_str(), "r");
      if (fp==NULL) {
         printf("create_new_plugin: Unable to open Makefile template.\n");
         exit(1);
      }
      // Stupid way to read in entire file.
      char line[256];
      std::stringstream file;
      while (fgets(line, sizeof(line), fp))
         file<<line;
      std::string filestring=file.str();
      fclose(fp);

      // Search and replace placeholders in the string
      boost::xpressive::sregex match_format=xpressive::as_xpr("@top_srcdir@");
      filestring=xpressive::regex_replace(filestring, match_format,
            format_top_srcdir);
      match_format=boost::xpressive::as_xpr("@top_objdir@");
      filestring=xpressive::regex_replace(filestring, match_format,
            format_top_objdir);
      match_format=boost::xpressive::as_xpr("@plugin@");
      filestring=xpressive::regex_replace(filestring, match_format,
            format_plugin);
      match_format=boost::xpressive::as_xpr("@Plugin@");
      filestring=xpressive::regex_replace(filestring, match_format, Name);
      match_format=boost::xpressive::as_xpr("@PLUGIN@");
      filestring=xpressive::regex_replace(filestring, match_format,
            format_PLUGIN);

      // Write the new file out
      fp=fopen(target_name.c_str(), "w");
      if (fp==0) {
         boost::filesystem::remove_all (plugin_name_);
         printf("Unable to create %s\n", target_name.c_str());
         exit(1);
      }
      fputs(filestring.c_str(), fp);
      fclose(fp);

      printf("\tCreated: %s\n", iter->second.c_str());
   }
}

void CMakePlugin::process() {
   OutFile("TestCMakeLists.txt",T);
   OutFile<<"add_subdirectory("<<plugin_name<<")"<<std::endl;

}

} //End namespace psi
