//
// file crash.H
// David Cosgrove
// AstraZeneca
// 3rd September 2008
//
// Crashes a program deliberately, so a debugger can find the exit point
// easily.

#ifndef DAC_CRASH_
#define DAC_CRASH_

#include <iostream>

namespace DACLIB {

  static void crash() {

    std::cout << "Deliberate crash from DACLIB::crash()" << std::endl;
    std::cerr << "Deliberate crash from DACLIB::crash()" << std::endl;
    int *i = 0;
    *i = 0;

  }

} // EO namespace DACLIB

#endif
