 FIND_PATH(FREE_IMAGE_INCLUDE_DIR FreeImage.h
   /usr/local/include
   /usr/include
   /usr/include/Cellar/freeimage/3.16.0/include
 )

 FIND_LIBRARY(FREE_IMAGE_LIBRARY NAMES freeimage PATHS
   /usr/local/lib
   /usr/lib
   /usr/include/Cellar/freeimage/3.16.0/lib
 )
