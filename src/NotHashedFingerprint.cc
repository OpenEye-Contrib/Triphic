//
// file NotHashedFingerprint.cc
// Dave Cosgrove
// AstraZeneca
// 3rd February 2009
//

#include <fstream>
#include <iostream>
#include <set>
#include <sstream>

#include <boost/lexical_cast.hpp>

#include "ByteSwapper.H"
#include "NotHashedFingerprint.H"

using namespace std;

namespace DAC_FINGERPRINTS {


  pNHDC NotHashedFingerprint::dist_calc_ = &NotHashedFingerprint::tanimoto;
  pNHTDC NotHashedFingerprint::threshold_dist_calc_ =
    &NotHashedFingerprint::tanimoto;
  
  // ****************************************************************************
  NotHashedFingerprint::NotHashedFingerprint() :
  FingerprintBase() , num_frag_nums_( 0 ) , frag_nums_( 0 ) {

  }

  // ****************************************************************************
  NotHashedFingerprint::NotHashedFingerprint( const string &name ) :
    FingerprintBase( name ) , num_frag_nums_( 0 ) , frag_nums_( 0 ) {

  }

  // ****************************************************************************
  // from the product of get_string_rep.
  NotHashedFingerprint::NotHashedFingerprint( const string &name ,
					      const string &rep ) :
    FingerprintBase( name ) , num_frag_nums_( 0 ) , frag_nums_( 0 ) {

    vector<uint32_t> tmp_bits;
    istringstream iss( rep.substr( 4 ) );
    while( 1 ) {
      uint32_t next_bit;
      iss >> next_bit;
      if( iss.fail() ) {
	break;
      }
      tmp_bits.push_back( next_bit );
    }

    build_from_vector( tmp_bits );

  }

  // ****************************************************************************
  NotHashedFingerprint::NotHashedFingerprint( const string &name ,
					      const std::vector<uint32_t> in_nums ) :
    FingerprintBase( name ) , num_frag_nums_( 0 ) , frag_nums_( 0 ) {
    
    build_from_vector( in_nums );

  }

  // ****************************************************************************
  NotHashedFingerprint::NotHashedFingerprint( const NotHashedFingerprint &fp ) :
    FingerprintBase() , num_frag_nums_( 0 ) , frag_nums_( 0 ) {

    copy_data( fp );

  }

  // ****************************************************************************
  NotHashedFingerprint::~NotHashedFingerprint() {

    delete [] frag_nums_;

  }

  // ****************************************************************************
  NotHashedFingerprint &NotHashedFingerprint::operator=( const NotHashedFingerprint &fp ) {

    if( this == &fp ) {
      return *this;
    }

    copy_data( fp );

    return *this;

  }

  // *******************************************************************************
  // bitwise operators for combining fingerprints.
  NotHashedFingerprint NotHashedFingerprint::operator&( const NotHashedFingerprint &rhs ) const {

    NotHashedFingerprint ret_val( "" );

    vector<uint32_t> v1( frag_nums_ , frag_nums_ + num_frag_nums_ );
    sort( v1.begin() , v1.end() );

    vector<uint32_t> v2( rhs.frag_nums_ , rhs.frag_nums_ + rhs.num_frag_nums_ );
    sort( v2.begin() , v2.end() );

    vector<uint32_t> v1_or_2( num_frag_nums_ + rhs.num_frag_nums_ , 0 );
    vector<uint32_t>::iterator un_end =
      set_intersection( v1.begin() , v1.end() , v2.begin() ,
			v2.end() , v1_or_2.begin() );

    delete [] ret_val.frag_nums_;
    ret_val.num_frag_nums_ = distance( v1_or_2.begin() , un_end );
    ret_val.frag_nums_ = new uint32_t[ret_val.num_frag_nums_];
    copy( v1_or_2.begin() , un_end , ret_val.frag_nums_ );

    return ret_val;

  }

  // *******************************************************************************
  NotHashedFingerprint NotHashedFingerprint::operator|( const NotHashedFingerprint &rhs ) const {

    NotHashedFingerprint ret_val( "" );

    vector<uint32_t> v1( frag_nums_ , frag_nums_ + num_frag_nums_ );
    sort( v1.begin() , v1.end() );

    vector<uint32_t> v2( rhs.frag_nums_ , rhs.frag_nums_ + rhs.num_frag_nums_ );
    sort( v2.begin() , v2.end() );

    vector<uint32_t> v1_or_2( num_frag_nums_ + rhs.num_frag_nums_ , 0 );
    vector<uint32_t>::iterator un_end =
      set_union( v1.begin() , v1.end() , v2.begin() , v2.end() , v1_or_2.begin() );

    delete [] ret_val.frag_nums_;
    ret_val.num_frag_nums_ = distance( v1_or_2.begin() , un_end );
    ret_val.frag_nums_ = new uint32_t[ret_val.num_frag_nums_];
    copy( v1_or_2.begin() , un_end , ret_val.frag_nums_ );

    return ret_val;

  }

  // ****************************************************************************
  NotHashedFingerprint &NotHashedFingerprint::operator&=( const NotHashedFingerprint &rhs ) {

    vector<uint32_t> v1( frag_nums_ , frag_nums_ + num_frag_nums_ );
    sort( v1.begin() , v1.end() );

    vector<uint32_t> v2( rhs.frag_nums_ , rhs.frag_nums_ + rhs.num_frag_nums_ );
    sort( v2.begin() , v2.end() );

    vector<uint32_t> v1_or_2( num_frag_nums_ + rhs.num_frag_nums_ , 0 );
    vector<uint32_t>::iterator un_end =
      set_intersection( v1.begin() , v1.end() , v2.begin() , v2.end() ,
			v1_or_2.begin() );

    delete [] frag_nums_;
    num_frag_nums_ = distance( v1_or_2.begin() , un_end );
    frag_nums_ = new uint32_t[num_frag_nums_];
    copy( v1_or_2.begin() , un_end , frag_nums_ );

    return *this;

  }

  // ****************************************************************************
  NotHashedFingerprint &NotHashedFingerprint::operator|=( const NotHashedFingerprint &rhs ) {

    vector<uint32_t> v1( frag_nums_ , frag_nums_ + num_frag_nums_ );
    sort( v1.begin() , v1.end() );

    vector<uint32_t> v2( rhs.frag_nums_ , rhs.frag_nums_ + rhs.num_frag_nums_ );
    sort( v2.begin() , v2.end() );

    vector<uint32_t> v1_or_2( num_frag_nums_ + rhs.num_frag_nums_ , 0 );
    vector<uint32_t>::iterator un_end =
      set_union( v1.begin() , v1.end() , v2.begin() , v2.end() , v1_or_2.begin() );

    delete [] frag_nums_;
    num_frag_nums_ = distance( v1_or_2.begin() , un_end );
    frag_nums_ = new uint32_t[num_frag_nums_];
    copy( v1_or_2.begin() , un_end , frag_nums_ );

    return *this;

  }

  // ****************************************************************************
  double NotHashedFingerprint::tanimoto( const NotHashedFingerprint &fp ) const {

    int num_comm = num_bits_in_common( fp );
    return( 1.0 - ( double( num_comm ) /
		    double( num_frag_nums_ + fp.num_frag_nums_ - num_comm ) ) );

  }

  // ****************************************************************************
  double NotHashedFingerprint::tanimoto( const NotHashedFingerprint &f ,
					 float thresh ) const {
    
    float min_dist;
    if( num_frag_nums_ < f.num_frag_nums_ )
      min_dist = 1.0 - float( num_frag_nums_ ) / float( f.num_frag_nums_ );
    else
      min_dist = 1.0 - float( f.num_frag_nums_ ) / float( num_frag_nums_ );
    if( min_dist > thresh ) {
      return 1.0;
    } else
      return tanimoto( f );

  }

  // ****************************************************************************
  double NotHashedFingerprint::tversky( const NotHashedFingerprint &f ) const {

    int num_in_a_not_b , num_in_b_not_a;
    int num_in_common = num_bits_in_common( f , num_in_a_not_b , num_in_b_not_a );
    
    double dist = 1.0 - ( double( num_in_common ) /
			  ( tversky_alpha_ * double( num_in_a_not_b ) +
			    ( 1.0 - tversky_alpha_ ) * double( num_in_b_not_a )
			    + double( num_in_common ) ) ); 

    return dist;

  }

  // ****************************************************************************
  double NotHashedFingerprint::tversky( const NotHashedFingerprint &f ,
					float thresh ) const {

    return tversky( f );

  }

  // ****************************************************************************
  void NotHashedFingerprint::binary_write( gzFile fp ) const {

    int name_len = finger_name_.length();
    gzwrite( fp , reinterpret_cast<char *>( &name_len ) , sizeof( int ) );
    gzwrite( fp , &finger_name_[0] , name_len );

    gzwrite( fp , reinterpret_cast<const char *>( &num_frag_nums_ ) ,
	     sizeof( int ) );
    gzwrite( fp , reinterpret_cast<const char *>( frag_nums_ ) ,
	     num_frag_nums_ * sizeof( uint32_t ) );
  
  }

  // ****************************************************************************
  void NotHashedFingerprint::binary_write( FILE *fp ) const {

    int name_len = finger_name_.length();
    fwrite( reinterpret_cast<char *>( &name_len ) , sizeof( int ) , 1 , fp );
    fwrite( &finger_name_[0] , 1 , name_len , fp );

    fwrite( reinterpret_cast<const char *>( &num_frag_nums_ ) , sizeof( int ) ,
	    1 , fp );
    fwrite( reinterpret_cast<const char *>( frag_nums_ ) , sizeof( uint32_t ) ,
	    num_frag_nums_ , fp );
  
  }

  // ****************************************************************************
  bool NotHashedFingerprint::binary_read( gzFile fp , bool byte_swapping ) {

    int name_len;
    gzread( fp , reinterpret_cast<void *>( &name_len ) , sizeof( int ) );
    if( gzeof( fp ) ) {
      return false;
    }
    if( byte_swapping ) {
      DACLIB::byte_swapper<int>( name_len );
    }
    finger_name_.resize( name_len , ' ' );
    gzread( fp , reinterpret_cast<void *>( &finger_name_[0] ) , name_len );

    gzread( fp , reinterpret_cast<void *>( &num_frag_nums_ ) , sizeof( int ) );
    if( byte_swapping ) {
      DACLIB::byte_swapper<int>( num_frag_nums_ );
    }

    delete [] frag_nums_;
    if( num_frag_nums_ ) {
      frag_nums_ = new uint32_t[num_frag_nums_];
      gzread( fp , reinterpret_cast<void *>( frag_nums_ ) ,
	      num_frag_nums_ * sizeof( uint32_t ) );
    } else {
      frag_nums_ = 0;
    }

    return true;

  }

  // ****************************************************************************
  bool NotHashedFingerprint::ascii_read( gzFile fp , const string &sep ) {

    // in FingerprintBase.cc
    string full_line = read_full_line( fp );
    if( gzeof( fp ) && full_line.empty() ) {
      return false;
    }

    full_line = convert_sep_to_new_sep( full_line , sep , " " );
    istringstream iss( full_line );
    iss >> finger_name_;

    vector<uint32_t> fns;
    while( 1 ) {
      uint32_t next_fn;
      iss >> next_fn;
      if( iss.fail() ) {
	break;
      }
      fns.push_back( next_fn );
    }

    build_from_vector( fns );

    return true;

  }

  // ****************************************************************************
  void NotHashedFingerprint::ascii_write( gzFile fp , const string &sep ) const {

    string act_sep = sep.empty() ? " " : sep;
    gzprintf( fp , "%s" , finger_name_.c_str() );

    // cowardly skirting round of the issue of what format string to use for
    // uint32_t in printf!  There seems to be a buffer size issue here,
    // strings greater than 4K are being truncated, so write each frag_num_
    // separately.
    ostringstream oss;
    for( int i = 0 ; i < num_frag_nums_ ; ++i ) {
      oss << act_sep << frag_nums_[i];
      gzprintf( fp , "%s" , oss.str().c_str() );
      oss.str( "" );
    }
    gzprintf( fp , "\n" );

  }

  // ****************************************************************************
  void NotHashedFingerprint::ascii_write( FILE *fp , const string &sep ) const {

    string act_sep = sep.empty() ? " " : sep;
    fprintf( fp , "%s" , finger_name_.c_str() );
    
    // cowardly skirting round of the issue of what format chars to use for
    // uint32_t in printf!
    ostringstream oss;
    for( int i = 0 ; i < num_frag_nums_ ; ++i ) {
      oss << act_sep << frag_nums_[i];
    }
    fprintf( fp , "%s" , oss.str().c_str() );
    fprintf( fp , "\n" );

  }

  // ****************************************************************************
  string NotHashedFingerprint::get_string_rep() const {

    ostringstream os;
    os << "__FN";
    for( int i = 0 ; i < num_frag_nums_ - 1 ; ++i ) {
      os << frag_nums_[i] << " ";
    }
    if( num_frag_nums_ ) {
      os << frag_nums_[num_frag_nums_ - 1];
    }

    return os.str();

  }

  // ****************************************************************************
  // count the number of bits in common between the fingerprint passed in
  // and this one
  int NotHashedFingerprint::num_bits_in_common( const NotHashedFingerprint &fp ) const {

    // both sets of frag_nums are sorted, so can walk through them in sequence
    uint32_t *these = frag_nums_ , *those = fp.frag_nums_;
    uint32_t *these_stop = frag_nums_ + num_frag_nums_;
    uint32_t *those_stop = fp.frag_nums_ + fp.num_frag_nums_;
    int num_comm = 0;
    while( these != these_stop && those != those_stop ) {
      if( *these < *those ) {
	while( these != these_stop && *these < *those ) {
	  ++these;
	}
      } else {
	while( those != those_stop && *those < *these ) {
	  ++those;
	} 
      }
      if( these == these_stop || those == those_stop ) {
	break;
      }
      if( *these == *those ) {
	++num_comm;
	++these;
	++those;
      }
    }

    return num_comm;

  }

  // ****************************************************************************
  // count the number of bits in common between the fingerprint passed in
  // and this one
  int NotHashedFingerprint::num_bits_in_common( const NotHashedFingerprint &fp ,
						int &num_in_a_not_b ,
						int &num_in_b_not_a ) const {

    num_in_a_not_b = num_in_b_not_a = 0;

    // both sets of frag_nums are sorted, so can walk through them in sequence
    uint32_t *these = frag_nums_ , *those = fp.frag_nums_;
    uint32_t *these_stop = frag_nums_ + num_frag_nums_;
    uint32_t *those_stop = fp.frag_nums_ + fp.num_frag_nums_;
    int num_comm = 0;
    while( these != these_stop && those != those_stop ) {
      // move these or those ( whichever has lower value) until they're equal
      // or one hits the end
      if( *these < *those ) {
	while( *these < *those ) {
	  ++these;
	  ++those;
	  ++num_in_a_not_b;
	}
      } else {
	while( *those < *these ) {
	  ++these;
	  ++those;
	  ++num_in_b_not_a;
	}
      }
      if( those == those_stop ) {
	while( these < these_stop ) {
	  ++these;
	  ++num_in_a_not_b;
	}
	break;
      } 
      if( these == these_stop ) {
	while( those < those_stop ) {
	  ++these;
	  ++num_in_b_not_a;
	}
	break;
      }
      if( *those == *these ) {
	++num_comm;
      }
      ++these;
      ++those;
    }

    return num_comm;

  }

  // *************************************************************************
  void NotHashedFingerprint::set_similarity_calc( SIMILARITY_CALC sc ) {

    switch( sc ) {
      case TANIMOTO :
	dist_calc_ = &NotHashedFingerprint::tanimoto;
	threshold_dist_calc_ = &NotHashedFingerprint::tanimoto;
	break;
      case TVERSKY :
	dist_calc_ = &NotHashedFingerprint::tversky;
	threshold_dist_calc_ = &NotHashedFingerprint::tversky;
	break;
    }

  }

  // **************************************************************************
  // calculate the distance between this fingerprint and the one passed in
  // using dist_calc_
  double NotHashedFingerprint::calc_distance( const FingerprintBase &f ) const {

    return f.calc_distance( *this );

  }

  // **************************************************************************
  // calculate the distance between this fingerprint and the one passed in
  // using dist_calc_
  double NotHashedFingerprint::calc_distance( const FingerprintBase &f ,
					      float threshold ) const {

    return f.calc_distance( *this , threshold );

  }

  // **************************************************************************
  // calculate the distance between this fingerprint and the one passed in
  // using dist_calc_
  double NotHashedFingerprint::calc_distance( const HashedFingerprint &f ) const {
    
    throw IncompatibleFingerprintError( "calc_distance" );

  }

  // **************************************************************************
  // calculate the distance between this fingerprint and the one passed
  // in using threshold_dist_calc_.  If the distance is predicted to be above
  // the threshold, return 1.0
  double NotHashedFingerprint::calc_distance( const HashedFingerprint &f ,
					      float threshold ) const {

    throw IncompatibleFingerprintError( "calc_distance" );

  }

  // ****************************************************************************
  // calculate the distance between this fingerprint and the one passed in
  // using dist_calc_
  double NotHashedFingerprint::calc_distance( const NotHashedFingerprint &f ) const {

    return (this->*dist_calc_)( f );

  }

  // ****************************************************************************
  // calculate the distance between this fingerprint and the one passed
  // in using threshold_dist_calc_.  If the distance is predicted to be above
  // the threshold, return 1.0
  double NotHashedFingerprint::calc_distance( const NotHashedFingerprint &f ,
					      float threshold ) const {
    
    return (this->*threshold_dist_calc_)( f , threshold );
  
  }

  // ****************************************************************************
  void NotHashedFingerprint::copy_data( const NotHashedFingerprint &fp ) {

    FingerprintBase::copy_data( fp );

    delete [] frag_nums_;

    num_frag_nums_ = fp.num_frag_nums_;
    if( fp.frag_nums_ && num_frag_nums_ ) {
      frag_nums_ = new uint32_t[num_frag_nums_];
      copy( fp.frag_nums_ , fp.frag_nums_ + num_frag_nums_ , frag_nums_ );
    } else {
      frag_nums_ = 0;
    }

  }

  // ****************************************************************************
  void NotHashedFingerprint::build_from_vector( const vector<uint32_t> &in_nums ) {

    delete [] frag_nums_;

    vector<uint32_t> tmp( in_nums );
    sort( tmp.begin() , tmp.end() );
    tmp.erase( unique( tmp.begin() , tmp.end() ) , tmp.end() );

    num_frag_nums_ = tmp.size();
    if( num_frag_nums_ ) {
      frag_nums_ = new uint32_t[num_frag_nums_];
      copy( tmp.begin() , tmp.end() , frag_nums_ );
    } else {
      frag_nums_ = 0;
    }

  }

} // end of namespace DAC_FINGERPRINTS
