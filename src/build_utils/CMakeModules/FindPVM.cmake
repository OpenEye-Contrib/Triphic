# attempts to find PVM using current values of PVM_ROOT and PVM_ARCH
# it will define
# PVM_FOUND it if finds everything it needs
# PVM_INCLUDE_DIR
# PVM_LIBRARIES

set( PVM_ROOT $ENV{PVM_ROOT})
set( PVM_ARCH $ENV{PVM_ARCH})

# include dir
find_path(PVM_INCLUDE_DIR pvm3.h PATH ${PVM_ROOT}/include)

message( "PVM_INCLUDE_DIR : ${PVM_INCLUDE_DIR}" )

if (NOT ${PVM_INCLUDE_DIR} MATCHES "-NOTFOUND$")
  message("PVM_INCLUDE_DIR found as ${PVM_INCLUDE_DIR}")
  set(PVM_FOUND True)
else ()
  message(FATAL_ERROR "PVM_INCLUDE_DIR not found.")
endif ()

# library
find_library( PVM_LIBRARY pvm3 ${PVM_ROOT}/lib/${PVM_ARCH} NO_DEFAULT_PATH)
if (NOT ${PVM_LIBRARY} MATCHES "-NOTFOUND$")
  message("Found pvm3 as ${PVM_LIBRARY}")
else ()
  message( FATAL_ERROR "Didn't find pvm3")
endif()
set(PVM_LIBRARIES ${PVM_LIBRARIES} ${PVM_LIBRARY})

message( "PVM_LIBRARIES : ${PVM_LIBRARIES}" )
