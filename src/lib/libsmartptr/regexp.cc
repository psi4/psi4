#include <stdio.h>
#include <string>
#include <iostream>
#include <sstream>
#include <ostream>
#include <fstream>

#include "regexp.h"

#include "boost/regex.hpp"
#include "boost/algorithm/string.hpp"

using namespace std;
using namespace smartptr;

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

bool 
smartptr::has_regexp_match(const string &regexp, const string &text, int flags)
{
    boost::smatch matches;
    
    boost::match_flag_type searchFlags = boost::match_default;
    // The default behavior is multiline, so we'd better fix that
    if(! flags & IncludeLineBreaks) searchFlags |= boost::match_single_line;
    // Make the search case-insensitive, if needed
    boost::regbase::flag_type  reFlags = boost::regex_constants::normal;
    if(flags & LowerCase  || flags & UpperCase) reFlags |= boost::regex_constants::icase;
    boost::regex re(regexp, reFlags);
    return boost::regex_search(text, matches, re, searchFlags);
}



double 
smartptr::get_regexp_double(const string &regexp, const string &text, int flags)
{
    boost::match_flag_type searchFlags = boost::match_default;
    // The default behavior is multiline, so we'd better fix that
    if(! flags & IncludeLineBreaks) searchFlags |= boost::match_single_line;
    double val = 0.0;
    boost::smatch matches;
    boost::regex re(regexp);
    if(boost::regex_search(text, matches, re, searchFlags)){
        istringstream iss(matches[1]);
        if((iss >> val).fail()){
            cout << "Unable to convert " << matches[1] << " to a double in get_regexp_double" << endl;
            abort();
        }
    }else{
        cout << text << " does not look like a double in get_regexp_double\n";
        abort();
    }
    return val;
}


string 
smartptr::get_regexp_string(const std::string &regexp, const string &text, int flags)
{
    boost::match_flag_type searchFlags = boost::match_default;
    // The default behavior is multiline, so we'd better fix that
    if(! flags & IncludeLineBreaks) searchFlags |= boost::match_single_line;
    boost::regbase::flag_type  reFlags = boost::regex_constants::normal;
    if(flags & LowerCase  || flags & UpperCase) reFlags |= boost::regex_constants::icase;
    boost::regex re(regexp, reFlags);
    string str;
    boost::smatch matches;
    if(boost::regex_search(text, matches, re, searchFlags)){
        str = matches[1];
    }else{
        cout << "Match failed in get_regexp_string" << endl;
    }
    return str;
}


double* 
smartptr::get_regexp_double_array(const string &regexp, const string &text, size_t& length, int flags)
{
    istringstream iss;
    double val;
    if(flags) cout << "WARNING: flags are ignored in get_regexp_double_array" << endl;
    boost::smatch matches;
    boost::match_flag_type searchFlags = boost::match_default;
    // The default behavior is multiline, so we'd better fix that
    if(! flags & IncludeLineBreaks) searchFlags |= boost::match_single_line;
    boost::regex re(regexp);
    string::const_iterator start = text.begin();
    string::const_iterator end   = text.end();
    vector<double> tempVals;
    while(regex_search(start, end, matches, re, searchFlags)){
        iss.clear();
        size_t nMatches = matches.size();
        for(size_t m = 1; m < nMatches; ++m){
            iss.clear();
            iss.str(matches[m]);
            if((iss >> val).fail()){
                cout << "Unable to convert " << iss.str() << " to a double in get_regexp_double_array" << endl;
                abort();
            }
            tempVals.push_back(val);
        }
        // update search position: 
        start = matches[0].second; 
        searchFlags |= boost::match_prev_avail; 
        searchFlags |= boost::match_not_bob;  
    } 
    length = tempVals.size();
    double* vals = new double[length];
    for(size_t i = 0; i < length;  ++i){
        vals[i] = tempVals[i];
    }
    return vals;
}


void
smartptr::findmatch(vector<string>& matches, const string &regexp, const string &text, int flags)
{
    boost::smatch reMatches;
    boost::match_flag_type searchFlags = boost::match_default;
    // The default behavior is multiline, so we'd better fix that
    if(! flags & IncludeLineBreaks) searchFlags |= boost::match_single_line;
    
    boost::regbase::flag_type  reFlags = boost::regex_constants::normal;
    if(flags & LowerCase  || flags & UpperCase) reFlags |= boost::regex_constants::icase;
    boost::regex re(regexp, reFlags);
    string::const_iterator start = text.begin();
    string::const_iterator end   = text.end();
    vector<double> tempVals;
    while(regex_search(start, end, reMatches, re, searchFlags)){
        size_t nMatches = reMatches.size();
        for(size_t m = 1; m < nMatches; ++m){
             matches.push_back(reMatches[m]);
             if(flags & UpperCase) boost::to_upper(matches.back());
             if(flags & LowerCase) boost::to_lower(matches.back());
             if(flags & StripWhitespace) boost::algorithm::trim(matches.back());
        }
        // update search position: 
        start = reMatches[0].second; 
        searchFlags |= boost::match_prev_avail; 
        searchFlags |= boost::match_not_bob;  
    } 
    return;
}

