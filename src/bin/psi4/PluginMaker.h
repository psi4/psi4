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
#ifndef PLUGINMAKER_H_
#define PLUGINMAKER_H_

#include <string>
#include <vector>
#include <utility>
#include <boost/filesystem.hpp>

namespace psi{
using std::string;
using std::vector;
using std::pair;

class PluginMaker{
   protected:
      string plugin_name_;

      vector<string,string> files_;
      ///Returns the plugin directory
      boost::filesystem::path PluginDir();
   public:
      PluginMaker(const string &plugin_name):
         plugin_name_(plugin_name){}
      virtual ~PluginMaker(){}
      void add_file(const string &source,const string &target=""){
               files_.push_back(
                     std::make_pair(source,(target==""?source:target)));
      }
      virtual void process()=0;
};

class GMakePlugin: PluginMaker{
   public:
      GMakePlugin(const string &plugin_name):PluginMaker(plugin_name){}
      void process();
};
}//End namespace psi


#endif /* PLUGINMAKER_H_ */
