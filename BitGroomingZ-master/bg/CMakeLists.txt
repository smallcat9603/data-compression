add_library (bg
  src/bg.c  
  src/ByteToolkit.c  
  src/callZlib.c  
  src/conf.c  
  src/dictionary.c  
  src/iniparser.c
  src/rw.c
)

target_include_directories(bg 
  PUBLIC 
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/bg>
  )

target_link_libraries (bg PUBLIC ${ZLIB_dep} m)

target_compile_options(bg
	PRIVATE $<$<CONFIG:Debug>:-Wextra -Wpedantic -Wno-unused-parameter>
	)

install (TARGETS bg ${thirdparty_export} EXPORT bgConfig
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

install (EXPORT bgConfig DESTINATION share/bg/cmake)
install (DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/bg)
export(TARGETS bg FILE bg.cmake)