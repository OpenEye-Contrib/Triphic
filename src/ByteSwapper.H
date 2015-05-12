//
// file ByteSwapper.H
// Dave Cosgrove
// AstraZeneca
// 9th May 2001
//
// This is a little template function that takes a variable and swaps the
// bytes in it from little-endian to big-endian and vice versa.

#include <algorithm>

#ifndef DAC_BYTE_SWAPPER
#define DAC_BYTE_SWAPPER

namespace DACLIB {

template <class T> void byte_swapper( T &x ) {

  union { T i; char by[sizeof( T )]; } bs;

  bs.i = x;
  for( unsigned int i = 0 ; i < sizeof( T ) / 2 ; i++ )
    std::swap( bs.by[i] , bs.by[sizeof( T ) - i - 1] );

  x = bs.i;

}

}

#endif
