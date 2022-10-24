add_library(ZLIB SHARED 
  ./gzclose.c
  ./uncompr.c
  ./trees.c
  ./gzread.c
  ./adler32.c
  ./gzwrite.c
  ./compress.c
  ./inftrees.c
  ./crc32.c
  ./inffast.c
  ./zutil.c
  ./gzlib.c
  ./infback.c
  ./inflate.c
  ./deflate.c
  )

target_include_directories(ZLIB
  PUBLIC 
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  )

install(TARGETS ZLIB ${thirdparty_export} EXPORT ZLIBConfig
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE
  DESTINATION ${CMAKE_INSTALL_LIBDIR})

install (EXPORT ZLIBConfig DESTINATION share/ZLIB/cmake)
install(FILES zlib.h zconf.h DESTINATION include)
export(TARGETS ZLIB FILE ZLIBConfig.cmake)
