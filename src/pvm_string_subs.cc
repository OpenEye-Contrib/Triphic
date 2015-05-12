//
// file pvm_string_subs.cc
// David Cosgrove
// AstraZeneca
// 6th August 2007
//
// This file contains stuff for passing STL strings with pvm.

#include <string>
#include <vector>
#include <pvm3.h>

#include "stddefs.H" // for make_buffer_big_enough

using namespace std;

namespace DACLIB {

  // **************************************************************************
  // pack a C++ string into a pvm buffer
  void pack_string( const string &str ) {

    int i = str.length() + 1;
    pvm_pkint( &i , 1 , 1 );
    pvm_pkstr( (char *) str.c_str() );

  }

  // **************************************************************************
  // unpack a C++ string from pvm buffer
  void unpack_string( string &str ) {

    int      i;
    pvm_upkint( &i , 1 , 1 );

    static char *buf = new char[1000];
    static int buf_len = 1000;
    DACLIB::make_buffer_big_enough( buf , i , buf_len );

    pvm_upkstr( buf );
    str = string( buf );

  }

  // **************************************************************************
  // pack a C++ string into a pvm buffer as bytes
  void pack_string_raw( const string &str ) {

    int i = str.length() + 1;
    pvm_pkint( &i , 1 , 1 );
    pvm_pkbyte( (char *) str.c_str() , i , 1 );
  
  }

  // **************************************************************************
  // unpack a C++ string from pvm buffer as bytes
  void unpack_string_raw( string &str ) {

    int      i;
    pvm_upkint( &i , 1 , 1 );

    static char *buf = new char[1000];
    static int buf_len = 1000;
    DACLIB::make_buffer_big_enough( buf , i , buf_len );

    pvm_upkbyte( buf , i , 1 );
    str = string( buf , i );

  }

  // **************************************************************************
  void pack_strings_vector( const vector<string> &strs ) {

    int num_to_send = strs.size();
    pvm_pkint( &num_to_send , 1 , 1 );
    for( int i = 0 ; i < num_to_send ; ++i )
      pack_string( strs[i] );

  }

  // **************************************************************************
  void unpack_strings_vector( vector<string> &strs ) {

    int num_to_rec;
    pvm_upkint( &num_to_rec , 1 , 1 );
    strs = vector<string>( num_to_rec );
    for( int i = 0 ; i < num_to_rec ; ++i )
      unpack_string( strs[i] );

  }

} // end of namespace
