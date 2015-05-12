# Attempts to find OEToolkits, using the current value of $OE_DIR.
# It assumes that it's all under $OE_DIR and should cope whether OE_DIR
# points to the toolkits dir or just the top level
# It will define:
# OEToolkits_FOUND if it finds everything it needs
# OEToolkits_INCLUDE_DIRS
# OEToolkits_LIBRARIES

set(OE_DIR $ENV{OE_DIR})

#message("OE_DIR ${OE_DIR}")
#message("COMPONENTS ${OEToolkits_FIND_COMPONENTS}")

#include dir
find_path(OEToolkits_INCLUDE_DIR
  oechem.h
  PATHS ${OE_DIR}/toolkits/include ${OE_DIR}/include NO_DEFAULT_PATH)

message("OEToolkits_INCLUDE_DIR : ${OEToolkits_INCLUDE_DIR}")

if (NOT ${OEToolkits_INCLUDE_DIR} MATCHES "-NOTFOUND$")
  message("OEToolkits_INCLUDE_DIR found as ${OEToolkits_INCLUDE_DIR}")
  set(OEToolkits_FOUND True)
else ()
  message(FATAL_ERROR "OEToolkits_INCLUDE_DIR not found.")
endif ()
set(OEToolkits_INCLUDE_DIRS ${OEToolkits_INCLUDE_DIR})

# libraries
foreach(component ${OEToolkits_FIND_COMPONENTS})
  message("Looking OEToolkits for ${component}")
  find_library(OEToolkits_LIBRARY_${component}
    ${component}
    PATHS ${OE_DIR}/toolkits/lib ${OE_DIR}/lib NO_DEFAULT_PATH)
  if (NOT ${OEToolkits_LIBRARY_${component}} MATCHES "-NOTFOUND$")
    message("Found ${component} as ${OEToolkits_LIBRARY_${component}}")
  else ()
    message( FATAL_ERROR "Didn't find library ${component}")
  endif()
  set(OEToolkits_LIBRARIES ${OEToolkits_LIBRARIES} ${OEToolkits_LIBRARY_${component}})
endforeach(component)

message("OEToolkits_INCLUDE_DIRS ${OEToolkits_INCLUDE_DIRS}")
message("OEToolkits_LIBRARIES : ${OEToolkits_LIBRARIES}")
message("OEToolkits_FOUND : ${OEToolkits_FOUND}")