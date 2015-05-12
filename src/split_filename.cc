//
// file split_filename.cc
// David Cosgrove
// AstraZeneca
// 7th August 2007
// 
// Splits a filename into root and extension, based on the last .

#include <iostream>
#include <string>

using namespace std;

namespace DACLIB {

  // ********************************************************************
  void split_filename( const string &filename ,
		       string &file_root , string &file_ext ) {

    size_t i = filename.rfind( '.' );
    if( i == string::npos ) {
      file_root = filename;
      file_ext = "";
    } else {
      file_root = filename.substr( 0 , i );
      file_ext = filename.substr( i + 1 );
    }

    // cout << filename << " : " << file_root << " : " << file_ext << endl;

  }

} // end of namespace DACLIB

