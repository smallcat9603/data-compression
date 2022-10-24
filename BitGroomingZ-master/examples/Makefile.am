AUTOMAKE_OPTIONS=foreign
AM_CFLAGS = -I../bg/include -I../zlib
AM_FCFLAGS = -I../bg/include -I../bg -I../zlib

bin_PROGRAMS=testfloat_compress testfloat_decompress

testfloat_compress_SOURCES=testfloat_compress.c
testfloat_compress_LDADD=../bg/.libs/libBG.a ../zlib/.libs/libzlib.a -lm
testfloat_decompress_SOURCES=testfloat_decompress.c
testfloat_decompress_LDADD=../bg/.libs/libBG.a ../zlib/.libs/libzlib.a -lm
