//
// file PharmPointPVM.cc
// David Cosgrove
// AstraZeneca
// 23rd May 2006
//
// This is the implementation for the class PharmPointPVM which derives
// from PharmPoint but has the ability to pack itself into and unpack itself
// from a PVM buffer.  It's separate from PharmPoint so that GUIs that use a
// PharmPoint don't need to link to pvm.

#include "PharmPointPVM.H"
#include "stddefs.H"

#include <pvm3.h>

using namespace std;

namespace DACLIB {
  // In pvm_string_subs.cc
  // pack a C++ string into a pvm buffer
  void pack_string( const string &str );
  // unpack a C++ string from pvm buffer
  void unpack_string( string &str );
  void pack_strings_vector( const vector<string> &strs );
  void unpack_strings_vector( vector<string> &strs );
}

// *****************************************************************************
// ****************************************************************************
// put the info into a previously initialised pvm buffer, but don't
// send it
void PharmPointPVM::pack_into_pvm_buffer() {

  int          num_to_send;

  map<string,vector<string> >::iterator p;
  num_to_send = points_defs_.size();
  pvm_pkint( &num_to_send , 1 , 1 );
  for( p = points_defs_.begin() ; p != points_defs_.end() ; p++ ) {
    DACLIB::pack_string( p->first );
    DACLIB::pack_strings_vector( p->second );
  }

  DACLIB::pack_strings_vector( h_vector_points_ );
  DACLIB::pack_strings_vector( lp_vector_points_ );
  DACLIB::pack_strings_vector( ring_normal_points_ );

  set<string>::iterator s;
  num_to_send = unique_smarts_.size();
  pvm_pkint( &num_to_send , 1 , 1 );
  for( s = unique_smarts_.begin() ; s != unique_smarts_.end() ; s++ )
    DACLIB::pack_string( *s );

  num_to_send = smarts_to_points_.size();
  pvm_pkint( &num_to_send , 1 , 1 );
  pvm_pkint( &smarts_to_points_[0] , num_to_send , 1 );

  num_to_send = int( hphobes_itmoc_ );
  pvm_pkint( &num_to_send , 1 , 1 );
  DACLIB::pack_string( itmoc_label_ );
  num_to_send = int( hphobes_itmoc_alo_ );
  pvm_pkint( &num_to_send , 1 , 1 );
  DACLIB::pack_string( itmoc_alo_label_ );

}

// ****************************************************************************
// pull it off a pvm buffer
void PharmPointPVM::unpack_from_pvm_buffer() {

  int      num_to_rec;

  pvm_upkint( &num_to_rec , 1 , 1 );
  string s1;
  vector<string> v1;
  for( int i = 0 ; i < num_to_rec ; ++i ) {
    DACLIB::unpack_string( s1 );
    DACLIB::unpack_strings_vector( v1 );
    points_defs_.insert( make_pair( s1 , v1 ) );
  }

  DACLIB::unpack_strings_vector( h_vector_points_ );
  DACLIB::unpack_strings_vector( lp_vector_points_ );
  DACLIB::unpack_strings_vector( ring_normal_points_ );

  pvm_upkint( &num_to_rec , 1 ,1 );
  for( int i = 0 ; i < num_to_rec ; ++i ) {
    DACLIB::unpack_string( s1 );
    unique_smarts_.insert( s1 );
  }

  pvm_upkint( &num_to_rec , 1 , 1 );
  smarts_to_points_ = vector<int>( num_to_rec , -1 );
  pvm_upkint( &smarts_to_points_[0] , num_to_rec , 1 );

  pvm_upkint( &num_to_rec , 1 , 1 );
  hphobes_itmoc_ = bool( num_to_rec );
  DACLIB::unpack_string( itmoc_label_ );
  pvm_upkint( &num_to_rec , 1 , 1 );
  hphobes_itmoc_alo_ = bool( num_to_rec );
  DACLIB::unpack_string( itmoc_alo_label_ );

}
