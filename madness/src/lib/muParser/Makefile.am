include $(top_srcdir)/config/MakeGlobal.am

lib_LIBRARIES = libMADmuparser.a

thisincludedir = $(includedir)/muParser
thisinclude_HEADERS = muParserBase.h muParserBytecode.h muParserCallback.h \
                      muParserDef.h muParserDLL.h muParserError.h \
                      muParserFixes.h muParser.h muParserInt.h \
                      muParserStack.h muParserTest.h muParserToken.h \
                      muParserTokenReader.h


libMADmuparser_a_SOURCES = muParserBase.cpp muParserBytecode.cpp \
                      muParserCallback.cpp \
                      muParser.cpp muParserError.cpp muParserTokenReader.cpp \
                      muParserBase.h muParserBytecode.h muParserCallback.h \
                      muParserDef.h muParserDLL.h muParserError.h \
                      muParserFixes.h muParser.h muParserInt.h \
                      muParserStack.h muParserTest.h muParserToken.h \
                      muParserTokenReader.h

libraries:	$(lib_LIBRARIES)
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir $@ ; done
