/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680


  $Id: textfsar.h 2173 2011-02-23 21:40:46Z justus.c79@gmail.com $
*/

#include <world/textfsar.h>

namespace madness {
    namespace archive {

        void TextFstreamOutputArchive::store(const char* t, long /*n*/) const {
            // store character string, escaping &, < and > along the way
            while (*t) {
                char c = *t++;
                if (c == '\\') {
                    os.put('\\');
                    os.put('\\');
                }
                else if (c == '<') {
                    os.put('\\');
                    os.put('l');
                }
                else if (c == '>') {
                    os.put('\\');
                    os.put('r');
                }
                else {
                    os.put(c);
                }
            }
            os << std::endl;
        }

        void TextFstreamOutputArchive::open(const char* filename,
                  std::ios_base::openmode mode)
        {
            os.open(filename, mode);
            os.setf(std::ios::scientific);
            os.precision(17);
            char tag[256];
            os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>" << std::endl;
            sprintf(tag,"<archive major_version=\"%d\" minor_version=\"%d\">",
                    ARCHIVE_MAJOR_VERSION, ARCHIVE_MINOR_VERSION);
            os << tag << std::endl;
            os << "<typemap>" << std::endl;
            for (int i=0; i<256; ++i) {
                sprintf(tag,"%d \"%s\"",i,archive_type_names[i]);
                store(tag,strlen(tag)); // Must use store to escape characters
            }
            os << "</typemap>" << std::endl;
        }


        void TextFstreamOutputArchive::close() {
            if (os.is_open()) {
                os << "</archive>" << std::endl;
                os.close();
            }
        }


        // Eat EOL after each entry to enable char-by-char read of strings
        void TextFstreamInputArchive::eat_eol() const {
            char eol;
            is.get(eol);
            if (eol != '\n')
                MADNESS_EXCEPTION("TextFstreamInputArchive: eat_eol: indigestion", static_cast<int>(eol));
        }


        void TextFstreamInputArchive::load(unsigned char* t, long n) const {
            for (long i=0; i<n; ++i) {
                unsigned int x;
                is >> x;
                t[i] = (unsigned char) x;
            }
            eat_eol();
        }

        void TextFstreamInputArchive::load(char* t, long n) const {
            for (long i=0; i<n; ++i) {
                char c0;
                is.get(c0);
                if (c0 == '\\') {
                    char c1;
                    is.get(c1);
                    if (c1 == '\\') {
                        t[i] = '\\';
                    }
                    else if (c1 == 'l') {
                        t[i] = '<';
                    }
                    else if (c1 == 'r') {
                        t[i] = '>';
                    }
                    else {
                      MADNESS_EXCEPTION("TextFstreamInputArchive: malformed string?",
                              static_cast<int>(c1));
                    }
                }
                else {
                    t[i] = c0;
                }
            }
            eat_eol();
        }

        void TextFstreamInputArchive::open(const char* filename, std::ios_base::openmode mode) {
            is.open(filename, mode);
            char buf[256];
            is.getline(buf,256);        // skip xml header
            is.getline(buf,256);

            char tag[256];
            sprintf(tag,"<archive major_version=\"%d\" minor_version=\"%d\">",
                    ARCHIVE_MAJOR_VERSION, ARCHIVE_MINOR_VERSION);
            if (strcmp(buf,tag)) {
                std::cout << "TextFstreamInputArchive: not an archive/bad version?" << std::endl;
                std::cout << "Found this: " << buf;
                std::cout << "Expected  : " << tag;
                MADNESS_EXCEPTION("TextFstreamInputArchive: not an archive/bad version?", 1);
            }

            // For now just skip over typemap
            for (int i=0; i<258; ++i) is.getline(buf,256);
        }

    }
}
