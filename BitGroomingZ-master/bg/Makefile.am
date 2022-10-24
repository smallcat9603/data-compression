#AM_CFLAGS = -I./include -I../zlib
#LDFLAGS=-fPIC -shared

AUTOMAKE_OPTIONS=foreign
include_HEADERS=include/defines.h include/rw.h include/conf.h include/callZlib.h\
		include/dictionary.h include/iniparser.h include/ByteToolkit.h include/bg.h

lib_LTLIBRARIES=libBG.la
libBG_la_CFLAGS=-I./include -I../zlib
libBG_la_LDFLAGS = -version-info  0:1:0
libBG_la_LIDADD=../zlib/.libs/libzlib.a
libBG_la_SOURCES=src/rw.c src/conf.c src/callZlib.c src/dictionary.c src/iniparser.c src/ByteToolkit.c src/bg.c

libBG_la_LINK= $(AM_V_CC)$(LIBTOOL) --tag=CC --mode=link $(CCLD) $(libBG_la_CFLAGS) -O3 $(libBG_la_LDFLAGS) -o $(lib_LTLIBRARIES)
