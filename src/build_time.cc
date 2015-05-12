//
// file build_time.cc
// David Cosgrove
// AstraZeneca
// 6th January 2005

// a weedy little file to be compiled each time that just makes a new string
// holding the build date and time from the preprocessor.  Because this date is
// only updated when the file is compiled.

#include <string>

std::string BUILD_TIME = std::string( __DATE__ ) + std::string( " at " ) +
                         std::string( __TIME__ );

