INCLUDE(CheckIncludeFile)
INCLUDE(CheckFunctionExists)
INCLUDE(CheckLibraryExists)
INCLUDE(CheckCCompilerFlag)

ADD_DEFINITIONS(-DHAVE_CONFIG_H)

string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_lc)

IF(cmake_build_type_lc STREQUAL "debug")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
ELSE()
    CHECK_C_COMPILER_FLAG(-Wno-stringop-truncation, COMPILER_NO_STRINGOP_TRUNCATION)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unused-variable")
    IF(COMPILER_NO_STRINGOP_TRUNCATION)
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-stringop-truncation")
    ENDIF()
ENDIF()

CHECK_INCLUDE_FILE("unistd.h" HAVE_UNISTD_H)
CHECK_INCLUDE_FILE("libgen.h" HAVE_LIBGEN_H)
CHECK_FUNCTION_EXISTS("strlcpy" HAVE_STRLCPY)
CHECK_FUNCTION_EXISTS("strlcat" HAVE_STRLCAT)
CHECK_FUNCTION_EXISTS("basename" HAVE_BASENAME)
CHECK_FUNCTION_EXISTS("dirname" HAVE_DIRNAME)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h)
INCLUDE_DIRECTORIES(BEFORE "${CMAKE_CURRENT_BINARY_DIR}")

ADD_EXECUTABLE(spagedi ${SPAGEDI_MACOSX_BUNDLE}
  main.c Autoccomp.c Autoccomp.h Xatools.c Xatools.h Autocio.c Autocio.h
  compat.c compat.h
  ${CMAKE_CURRENT_BINARY_DIR}/config.h
  ${RES_SRCS})

TARGET_LINK_LIBRARIES(spagedi m)

INSTALL(TARGETS spagedi RUNTIME DESTINATION bin ${SPAGEDI_INSTALL_EXTRA_ARGS})
