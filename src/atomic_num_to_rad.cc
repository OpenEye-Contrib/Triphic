//
// file atomic_num_to_rad.cc
// David Cosgrove
// AstraZeneca
// 25th May 2006
//
// Takes an atomic number and returns a reasonable radius for the atom.
// These numbers are the ones used in skinny, taken from the Sybyl forcefield
// I believe.

// *****************************************************************************

namespace DACLIB {

  float atomic_num_to_rad( int atomic_num ) {

    switch( atomic_num ) {
      case 1 : return 1.0;
      case 3 : return 1.1;
      case 6 : return 1.7;
      case 7 : return 1.65;
      case 8 : return 1.6;
      case 9 : return 1.47;
      case 11 : return 1.1;
      case 13 : return 2.0;
      case 14 : return 2.1;
      case 15 : return 1.8;
      case 16 : return 1.9;
      case 17 : return 1.75;
      case 19 : return 1.3;
      case 20 : return 1.5;
      case 35 : return 1.85;
      case 53 : return 1.98;
    }

    return 1.5F; // default value

  }

}
