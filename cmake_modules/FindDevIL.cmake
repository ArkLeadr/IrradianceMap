MESSAGE("ENV DEVILDIR is $ENV{DEVILDIR}")

 FIND_PATH(DEVIL_INCLUDE_DIR IL/il.h
   HINTS
   $ENV{DEVILDIR}
   PATH_SUFFIXES include/IL include
   PATHS
   ~/Library/Frameworks
   /Library/Frameworks
   /usr/local/include/GL
   /usr/include/GL
   /sw # Fink
   /opt/local # DarwinPorts
   /opt/csw # Blastwave
   /opt
 )
 MESSAGE("DEVIL_INCLUDE_DIR is ${DEVIL_INCLUDE_DIR}")

 FIND_LIBRARY(DEVIL_LIBRARY
   NAMES DevIL IL
   HINTS
   $ENV{DEVILDIR}
   PATH_SUFFIXES lib
   PATHS
   /sw
   /opt/local
   /opt/csw
   /opt
 )

 MESSAGE("DEVIL_LIBRARY is ${DEVIL_LIBRARY}")

 SET(DEVIL_FOUND "NO")
 IF(DEVIL_LIBRARY)
   # For SDL2main

   SET(DEVIL_FOUND "YES")
 ENDIF(DEVIL_LIBRARY)


INCLUDE(FindPackageHandleStandardArgs)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(DevIL
                                  REQUIRED_VARS DEVIL_LIBRARY DEVIL_INCLUDE_DIR)
