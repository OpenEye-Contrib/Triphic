//
// file MagicInts.H
// Dave Cosgrove
// AstraZeneca
// 3rd February 2009
//
// Integer consts used in Flush fingerprint files for determining endianness.

#ifndef MAGIC_INTS_DEFINED
#define MAGIC_INTS_DEFINED

namespace DAC_FINGERPRINTS {

  // FP_MAGIC_INT denotes an original fp file
  // on a bigendian machine
  static const unsigned int FP_MAGIC_INT = 'F' << 24 | '0' << 16 | '0' << 8 | '1';
  static const unsigned int FN_MAGIC_INT = 'N' << 24 | '0' << 16 | '0' << 8 | '1';
  static const unsigned int BUGGERED_FP_MAGIC_INT = '1' << 24 | '0' << 16 | '0' << 8 | 'F';
  static const unsigned int BUGGERED_FN_MAGIC_INT = '1' << 24 | '0' << 16 | '0' << 8 | 'N';
  
  // as it appears on a littleendian machine

}

#endif
