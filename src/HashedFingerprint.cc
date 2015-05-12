//
// file HashedFingerprint.cc
// Dave Cosgrove
// AstraZeneca
// 29th January 2009
//

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>

#include "ByteSwapper.H"
#include "HashedFingerprint.H"

using namespace std;

#if defined(__GNUC__) && defined(__SSSE3__)
#include "popcount_ssse3.cc"
#elif defined(__GNUC__) && defined(__SSE2__)
#include "popcount_ssse2.cc"
#endif
namespace DAC_FINGERPRINTS {

pHDC HashedFingerprint::dist_calc_ = &HashedFingerprint::tanimoto;
pHTDC HashedFingerprint::threshold_dist_calc_ = &HashedFingerprint::tanimoto;
unsigned int HashedFingerprint::num_ints_ = 0;

// **************************************************************************
HashedFingerprint::HashedFingerprint() :
  FingerprintBase() , finger_bits_( 0 ) , num_bits_set_( 0 ) {

  make_zero_fp( num_ints_ * sizeof( unsigned int ) );

}

// **************************************************************************
HashedFingerprint::HashedFingerprint( const string &name ) :
  FingerprintBase( name ) , finger_bits_( 0 ) , num_bits_set_( 0 ) {

  make_zero_fp( num_ints_ * sizeof( unsigned int ) );

}

// **************************************************************************
// this one from the contents of a string built with get_string_rep.
HashedFingerprint::HashedFingerprint( const string &name ,
                                      const string &rep ) :
  FingerprintBase( name ) , finger_bits_( 0 ) , num_bits_set_( 0 ) {

  // 10 is the number of chars needed to write
  // numeric_limits<unsigned int>::max() in ascii. At least when I built
  // this in Dec 2014. If I wanted to be really sure, I'd work it out
  // at runtime and put it in a static variable, I guess
  unsigned int new_num_ints = rep.length() / 10;
  if( rep.length() % 10 ) {
    ++new_num_ints;
  }

  if( !num_ints_ ) {
    num_ints_ = new_num_ints;
  }

  if( new_num_ints != num_ints_ ) {
    throw HashedFingerprintLengthError( num_ints_ , new_num_ints );
  }

  finger_bits_ = new unsigned int[num_ints_];

  int next_bit = 0;
  for( int i = num_ints_ - 1 ; i >= 0 ; --i , next_bit += 10 ) {
    finger_bits_[i] = boost::lexical_cast<unsigned int>( rep.substr( next_bit , 10 ) );
  }

  count_bits();

}

// **************************************************************************
HashedFingerprint::HashedFingerprint( const string &name ,
                                      unsigned int *new_ints ) :
  FingerprintBase( name ) , finger_bits_( 0 ) , num_bits_set_( 0 ) {

  if( num_ints_ > 0 ) {
    finger_bits_ = new unsigned int[num_ints_];
    copy( new_ints , new_ints + num_ints_ , finger_bits_ );
  }

  count_bits();

}

// **************************************************************************
HashedFingerprint::HashedFingerprint( const HashedFingerprint &fp ) :
  FingerprintBase() , finger_bits_( 0 ) , num_bits_set_( 0 ) {

  copy_data( fp );

}

// **************************************************************************
HashedFingerprint::~HashedFingerprint() {

  delete [] finger_bits_;

}

// **************************************************************************
HashedFingerprint &HashedFingerprint::operator=( const HashedFingerprint &fp ) {

  if( &fp == this ) {
    return *this;
  }

  copy_data( fp );

  return *this;

}

// *******************************************************************************
// bitwise operators for combining fingerprints.
HashedFingerprint HashedFingerprint::operator&( const HashedFingerprint &rhs ) const {

  HashedFingerprint ret_val( "" );
  for( unsigned int i = 0 ; i < HashedFingerprint::num_ints() ; ++i ) {
    ret_val.finger_bits_[i] = finger_bits_[i] & rhs.finger_bits_[i];
  }
  ret_val.num_bits_set_ = 0; // to force a re-count if count_bits() is called

  return ret_val;

}

// *******************************************************************************
HashedFingerprint HashedFingerprint::operator|( const HashedFingerprint &rhs ) const {

  HashedFingerprint ret_val( "" );
  for( unsigned int i = 0 ; i < HashedFingerprint::num_ints() ; ++i ) {
    ret_val.finger_bits_[i] = finger_bits_[i] | rhs.finger_bits_[i];
  }
  ret_val.num_bits_set_ = 0; // to force a re-count if count_bits() is called

  return ret_val;

}

// *******************************************************************************
HashedFingerprint &HashedFingerprint::operator&=( const HashedFingerprint &rhs ) {

  for( unsigned int i = 0 ; i < num_ints_ ; ++i ) {
    finger_bits_[i] &= rhs.finger_bits_[i];
  }
  num_bits_set_ = 0; // to force a re-count if count_bits() is called

  return *this;

}

// *******************************************************************************
HashedFingerprint &HashedFingerprint::operator|=( const HashedFingerprint &rhs ) {

  for( unsigned int i = 0 ; i < num_ints_ ; ++i ) {
    finger_bits_[i] |= rhs.finger_bits_[i];
  }
  num_bits_set_ = 0; // to force a re-count if count_bits() is called

  return *this;

}

// *******************************************************************************
// work out how many chars are required to accommodate the number of
// bits requested.
unsigned int HashedFingerprint::calc_num_ints_req( int num_bits ) const {

  unsigned int nc = num_bits / sizeof( unsigned int );
  if( num_bits % sizeof( unsigned int ) ) {
    ++nc;
  }

  return nc;

}

// *******************************************************************************
// make empty fingerprint of the required size.
void HashedFingerprint::make_zero_fp( int num_bits ) {

  unsigned int new_num_ints = calc_num_ints_req( num_bits );

  if( num_ints_ && new_num_ints != num_ints_ ) {
    throw HashedFingerprintLengthError( num_ints_ , new_num_ints );
  }

  if( !num_ints_ ) {
    num_ints_ = new_num_ints;
  }

  if( num_ints_ ) {
    if( finger_bits_ && new_num_ints != num_ints_ ) {
      delete [] finger_bits_;
      finger_bits_ = 0;
    }
    if( !finger_bits_ ) {
      finger_bits_ = new unsigned int[num_ints_];
    }
    std::fill_n( finger_bits_ , num_ints_ , 0 );
  }

}

// **************************************************************************
double HashedFingerprint::tanimoto( const HashedFingerprint &f ) const {

  if( !num_bits_set_ && !f.num_bits_set_ ) {
    return 0.0; // otherwise, we'll get a NaN.
  }

  int num_in_common = num_bits_in_common( f );

  double dist = 1.0 - ( double( num_in_common ) /
                        double( num_bits_set_ + f.num_bits_set_ - num_in_common ));

  return dist;

}

// **************************************************************************
// If the distance is predicted to be above the threshold, return 1.0
double HashedFingerprint::tanimoto( const HashedFingerprint &f ,
                                    float threshold ) const {

  float min_dist;
  if( num_bits_set_ < f.num_bits_set_ ) {
    min_dist = 1.0 - float( num_bits_set_ ) / float( f.num_bits_set_ );
  } else {
    min_dist = 1.0 - float( f.num_bits_set_ ) / float( num_bits_set_ );
  }
  if( min_dist > threshold ) {
    return 1.0;
  } else {
    return tanimoto( f );
  }

}

// **************************************************************************
double HashedFingerprint::tversky( const HashedFingerprint &f ) const {

  int num_a_not_b , num_b_not_a;

  int num_in_common = num_bits_in_common( f , num_a_not_b , num_b_not_a );

  double dist = 1.0 - ( double( num_in_common ) /
                        ( tversky_alpha_ * double( num_a_not_b ) +
                          ( 1.0 - tversky_alpha_ ) * double( num_b_not_a )
                          + double( num_in_common ) ) );

  return dist;

}

// **************************************************************************
double HashedFingerprint::tversky( const HashedFingerprint &f ,
                                   float threshold ) const {

  // don't have a calculation for the threshold, so just return the normal
  // one
  return tversky( f );

}

// ***************************************************************************
int HashedFingerprint::num_bits_in_common( const HashedFingerprint &f ) const {

  static unsigned int *common_bits = 0;
  if( !common_bits ) {
    common_bits = new unsigned int[num_ints_];
  }

  const unsigned int *these = finger_bits_ , *those = f.finger_bits_;
  unsigned int *cb = common_bits;

  transform( these , these + num_ints_ , those , cb , std::bit_and<unsigned int>() );

  return count_bits_set( common_bits , num_ints_ );

}

// ***************************************************************************
int HashedFingerprint::num_bits_in_common( const HashedFingerprint &f ,
                                           int &num_in_a_not_b ,
                                           int &num_in_b_not_a ) const {

  static unsigned int *common_bits = 0 , *in_a_not_b = 0 , *in_b_not_a = 0;
  if( !common_bits ) {
    common_bits = new unsigned int[num_ints_];
    in_a_not_b = new unsigned int[num_ints_];
    in_b_not_a = new unsigned int[num_ints_];
  }

  num_in_a_not_b = num_in_b_not_a = 0;

  unsigned int *cb = common_bits , *anb = in_a_not_b , *bna = in_b_not_a;
  const unsigned int *these = finger_bits_ , *those = f.finger_bits_;
  for( unsigned int i = 0 ; i < num_ints_ ; ++i , ++cb , ++anb , ++bna , ++these , ++those ) {
    *cb = *these & *those;
    *anb = *these & ( ~ *those );
    *bna = *those & ( ~ *these );
  }

  num_in_a_not_b = count_bits_set( in_a_not_b , num_ints_ );
  num_in_b_not_a = count_bits_set( in_b_not_a , num_ints_ );
  return count_bits_set( common_bits , num_ints_ );

}

// ***************************************************************************
int HashedFingerprint::num_set_in_this_and_not_in_2( const HashedFingerprint &fp2 ) const {

  static unsigned int *in_this_not_2 = 0;
  if( !in_this_not_2 ) {
    in_this_not_2 = new unsigned int[num_ints_];
  }

  unsigned int *itn2 = in_this_not_2;
  const unsigned int *these = finger_bits_ , *those = fp2.finger_bits_;

  for( unsigned int i = 0 ; i < num_ints_ ; ++i , ++itn2 , ++these , ++those ) {
    *itn2 = *these & ( ~ *those );
  }

  return count_bits_set( in_this_not_2 , num_ints_ );

}

// **************************************************************************
string HashedFingerprint::get_string_rep() const {

  // 10 is the number of chars needed to write
  // numeric_limits<unsigned int>::max() in ascii. At least when I built
  // this in Dec 2014. If I wanted to be really sure, I'd work it out
  // at runtime and put it in a static variable, I guess
  ostringstream oss;
  oss << "__FP";
  for( int i = num_ints_ - 1 ; i >= 0 ; --i )
    oss << setfill( '0' ) << setw( 10 ) << finger_bits_[i];

  return oss.str();

}

// *************************************************************************
void HashedFingerprint::set_similarity_calc( SIMILARITY_CALC sc ) {

  switch( sc ) {
  case TANIMOTO :
    dist_calc_ = &HashedFingerprint::tanimoto;
    threshold_dist_calc_ = &HashedFingerprint::tanimoto;
    break;
  case TVERSKY :
    dist_calc_ = &HashedFingerprint::tversky;
    threshold_dist_calc_ = &HashedFingerprint::tversky;
    break;
  }

}

// **************************************************************************
// count the number of set bits in the fingerprint
int HashedFingerprint::count_bits() const {

  if( !num_bits_set_ ) {
    num_bits_set_ = count_bits_set( finger_bits_ , num_ints_ );
  }

  return num_bits_set_;

}

// **************************************************************************
// calculate the distance between this fingerprint and the one passed in
// using dist_calc_
double HashedFingerprint::calc_distance( const FingerprintBase &f ) const {

  return f.calc_distance( *this );

}

// **************************************************************************
// calculate the distance between this fingerprint and the one passed in
// using dist_calc_
double HashedFingerprint::calc_distance( const FingerprintBase &f ,
                                         float threshold ) const {

  return f.calc_distance( *this , threshold );

}

// **************************************************************************
// calculate the distance between this fingerprint and the one passed in
// using dist_calc_
double HashedFingerprint::calc_distance( const HashedFingerprint &f ) const {

  return (this->*dist_calc_)( f );

}

// **************************************************************************
// calculate the distance between this fingerprint and the one passed
// in using threshold_dist_calc_.  If the distance is predicted to be above
// the threshold, return 1.0
double HashedFingerprint::calc_distance( const HashedFingerprint &f ,
                                         float threshold ) const {

  return (this->*threshold_dist_calc_)( f , threshold );

}

// **************************************************************************
// calculate the distance between this fingerprint and the one passed in
// using dist_calc_
double HashedFingerprint::calc_distance( const NotHashedFingerprint &f ) const {

  throw IncompatibleFingerprintError( "calc_distance" );

}

// **************************************************************************
// calculate the distance between this fingerprint and the one passed
// in using threshold_dist_calc_.  If the distance is predicted to be above
// the threshold, return 1.0
double HashedFingerprint::calc_distance( const NotHashedFingerprint &f ,
                                         float threshold ) const {

  throw IncompatibleFingerprintError( "calc_distance" );

}

// **************************************************************************
bool HashedFingerprint::binary_read( gzFile fp , bool byte_swapping ) {

  int name_len;
  gzread( fp , &name_len , sizeof( int ) );
  if( gzeof( fp ) ) {
    return false;
  }
  if( byte_swapping ) DACLIB::byte_swapper<int>( name_len );

  finger_name_.resize( name_len , ' ' );
  gzread( fp , &finger_name_[0] , name_len );
  gzread( fp , reinterpret_cast<void *>( finger_bits_ ) , num_ints_ * sizeof( unsigned int ) );

  num_bits_set_ = 0; // to force a re-count
  count_bits();

  return true;

}

// **************************************************************************
// write to a previously initialised flush output file
void HashedFingerprint::binary_write( gzFile fp ) const {

  int name_len = finger_name_.length();
  gzwrite( fp , reinterpret_cast<void *>( &name_len ) , sizeof( int ) );
  gzwrite( fp , &finger_name_[0] , name_len );
  gzwrite( fp , reinterpret_cast<void *>( finger_bits_ ) , num_ints_ * sizeof( unsigned int ) );

}

// **************************************************************************
// write to a previously initialised flush output file
void HashedFingerprint::binary_write( FILE *fp ) const {

  int name_len = finger_name_.length();
  fwrite( reinterpret_cast<void *>( &name_len ) , sizeof( int ) , 1 , fp );
  fwrite( finger_name_.c_str() , 1 , name_len , fp );
  fwrite( reinterpret_cast<void *>( finger_bits_ ) , 1 , num_ints_ * sizeof( unsigned int ), fp );

}

// **************************************************************************
// write an ascii representation
bool HashedFingerprint::ascii_read( gzFile fp , const string &sep ) {

  // in FingerprintBase
  string full_line = read_full_line( fp );
  if( gzeof( fp ) && full_line.empty() ) {
    return false;
  }

  string name , rest;
  // if sep is empty, split into 2 strings on space
  if( sep.empty() ) {
    size_t space_pos = full_line.find( ' ' );
    if( string::npos == space_pos ) {
      cerr << "Error reading fingerprint line : " << full_line << endl;
      exit( 1 );
    }
    name = full_line.substr( 0 , space_pos );
    rest = full_line.substr( space_pos + 1 );
  } else {
    // it's a bit more complicated. change sep to space, take out the 1st
    // string as name, stick the rest back together again and proceed
    // convert_sep_to_new_sep in FingerprintBase
    string new_fl = convert_sep_to_new_sep( full_line , sep , " " );
    size_t space_pos = new_fl.find( ' ' );
    if( string::npos == space_pos ) {
      cerr << "Error reading fingerprint line : " << full_line << endl;
      exit( 1 );
    }
    name = full_line.substr( 0 , space_pos );
    rest = convert_sep_to_new_sep( full_line.substr( space_pos + 1 ) , " " , "" );
  }

  // first time through, num_chars() should be zero, as we won't know at this
  // stage what we're dealing with.
  if( !num_ints() ) {
    try {
      make_zero_fp( rest.length() );
    } catch( HashedFingerprintLengthError &e ) {
      cerr << "Caught : " << e.what() << endl;
      exit( 1 );
    }
  }
  build_fp_from_bitstring( name , rest );

  return true;

}

// **************************************************************************
// write an ascii representation
void HashedFingerprint::ascii_write( gzFile fp , const string &sep ) const {

  gzprintf( fp , "%s" , finger_name_.c_str() );
  if( sep.empty() ) {
    gzprintf( fp , " " );
  }
  for( unsigned int i = 0 ; i < num_ints_ ; ++i ) {
    for( unsigned int t = 1 << (sizeof( unsigned int ) - 1) ; t ; t = t >> 1 ) {
      if( finger_bits_[i] & t ) {
        gzprintf( fp , "%s1" , sep.c_str() );
      } else {
        gzprintf( fp , "%s0" , sep.c_str() );
      }
    }
  }
  gzprintf( fp , "\n" );

}

// **************************************************************************
// write an ascii representation
void HashedFingerprint::ascii_write( FILE *fp , const string &sep ) const {

#ifdef NOTYET
  cout << "HashedFingerprint::ascii_write, num_ints_ = " << num_ints_ << endl;
#endif

  // set the bits
  int num_bits_per_int = sizeof( int ) * 8;
  static unsigned int *bit_masks = 0;
  if( !bit_masks ) {
    bit_masks = new unsigned int[num_bits_per_int];
    int b = 1;
    for( int i = 0 ; i < num_bits_per_int ; ++i ) {
      bit_masks[i] = b;
      b <<= 1;
    }
  }

  fprintf( fp , "%s" , finger_name_.c_str() );
  if( sep.empty() ) {
    fprintf( fp , " " );
  }
  for( unsigned int i = 0 ; i < num_ints_ ; ++i ) {
    for( int j = num_bits_per_int - 1 ; j >= 0 ; --j ) {
      if( finger_bits_[i] & bit_masks[j] ) {
        fprintf( fp , "%s1" , sep.c_str() );
      } else {
        fprintf( fp , "%s0" , sep.c_str() );
      }
    }
  }
  fprintf( fp , "\n" );

}

// **************************************************************************
void HashedFingerprint::copy_data( const HashedFingerprint &fp ) {

  FingerprintBase::copy_data( fp );

  if( num_ints_ && !finger_bits_ ) {
    finger_bits_ = new unsigned int[num_ints_];
  }
  copy( fp.finger_bits_ , fp.finger_bits_ + num_ints_ , finger_bits_ );

  num_bits_set_ = fp.num_bits_set_;

}

// **************************************************************************
void HashedFingerprint::build_fp_from_bitstring( const string &name ,
                                                 const string &bitstring ) {

  finger_name_ = name;
  unsigned int new_num_ints = calc_num_ints_req( bitstring.length() );
  if( num_ints_ && new_num_ints != num_ints_ ) {
    throw HashedFingerprintLengthError( new_num_ints , num_ints_ );
  }
  // obviously the bit string is the 'wrong way round' wrt least significant
  // bit

  // make sure that there are an integral number of blocks of sizeof( unsigned int )
  static string spare_bits;
  if( spare_bits.empty() ) {
    for( unsigned int i = 0 ; i < sizeof( unsigned int ) ; ++i ) {
      spare_bits += "0";
    }
  }
  int padding_needed = sizeof( unsigned int ) - bitstring.length() % sizeof( unsigned int );
  string bits_to_use = padding_needed ?
        spare_bits.substr( 0 , padding_needed ) + bitstring : bitstring;
  num_bits_set_ = 0;
  for( unsigned int i = 0 , j = bits_to_use.length() - sizeof( unsigned int ) ;
       i < num_ints_ ; ++i , j -= sizeof( unsigned int) ) {
    string next_bits = bits_to_use.substr( j , sizeof( unsigned int ) );
    finger_bits_[i] = 0;
    int k = 1 << (sizeof( unsigned int ) - 1);
    for( unsigned int l = 0 ; l < sizeof( unsigned int ) ; ++l ) {
      if( '1' == next_bits[l] ) {
        ++num_bits_set_;
        finger_bits_[i] |= k;
      }
      k >>= 1;
    }
  }

}

// ****************************************************************************
int count_bits_set( unsigned int *bits , int num_ints ) {

#if defined(__GNUC__) && defined(__SSSE3__)
  return popcount_ssse3( bits , num_ints );
#elif defined(__GNUC__) && defined(__SSE2__)
  return popcount_ssse2( bits , num_ints );
#else
  return 0;
#endif

}

// ****************************************************************************
// thrown when the program tries to change num_chars_.
HashedFingerprintLengthError::HashedFingerprintLengthError( int new_len ,
                                                            int old_len ) {

  msg_ = string( "Changed fingerprint length from " ) +
      boost::lexical_cast<string>( old_len ) + string( " to " ) +
      boost::lexical_cast<string>( new_len ) + string( "." );

}

} // end of namespace DAC_FINGERPRINTS

