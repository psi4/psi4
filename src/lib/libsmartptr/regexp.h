#ifndef smartptr_regexp_h
#define smartptr_regexp_h

#include <string>
#include <vector>

namespace smartptr {

enum RegexpFlags { IgnoreCase=1, UpperCase=2, FindAll=4, LowerCase=8, IncludeLineBreaks=16, StripWhitespace=32 };

/**
 * Extracts a double precision number from text.
 * @param regexp The regular expression containing the (captured) number.
 * @param text The text to apply the RegEx to.
 * @param flags An integer comprising bitwise options (see the PyFlags enum).
 * @return The value as a double.
 */
double get_regexp_double(const std::string &regexp, const std::string &text, int flags = 0);

/**
 * Extracts a string from text.
 * @param regexp The regular expression containing the (captured) string.
 * @param text The text to apply the RegEx to.
 * @param flags An integer comprising bitwise options (see the PyFlags enum).
 * @return The string.
 */
std::string get_regexp_string(const std::string &regexp, const std::string &text, int flags = 0);

/**
 * Extracts an array of doubles from text, according to a given RegEx.
 * Memory is allocated by the function and ownership is passed to calling routine.
 * @param regexp The regular expression containing the (captured) values.
 * @param text The text to apply the RegEx to.
 * @param flags An integer comprising bitwise options (see the PyFlags enum).
 * @return an array of all captured values.
 */
double* get_regexp_double_array(const std::string& regexp, const std::string &text, size_t& length, int flags = 0);

/**
 * Extracts strings from text, according to a given RegEx.
 * @param matches The vector of strings to hold the results.
 * @param regexp The regular expression containing the (captured) strings.
 * @param text The text to apply the RegEx to.
 * @param flags An integer comprising bitwise options (see the PyFlags enum).
 */
void
findmatch(std::vector<std::string>& matches, const std::string &regexp, const std::string &text, int flags = 0);

/**
 * Tests a string to see if it contains a certain RegEx.
 * @param regexp The regular expression to test.
 * @param text The text to apply the RegEx to.
 * @param flags An integer comprising bitwise options (see the PyFlags enum).
 * @return whether the text contains a match.
 */
bool has_regexp_match(const std::string& regexp, const std::string& text, int flags = 0);

}

#endif
