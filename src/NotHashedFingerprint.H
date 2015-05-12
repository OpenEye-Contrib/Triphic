//
// file NotHashedFingerprint.H
// Dave Cosgrove
// AstraZeneca
// 3rd February 2009
//
// This class, derived from FingerprintBase stores fingerprints as the original
// fragment numbers, without splitting them into a bit-string.  They're a bit
// slower to use, but shouldn't suffer from data loss by being dropped into
// a small bitstring.

#ifndef DAC_NOT_HASHED_FINGERPRINT
#define DAC_NOT_HASHED_FINGERPRINT

#include <stdint.h>
#include <vector>

#include "FingerprintBase.H"

namespace DAC_FINGERPRINTS {

// ***********************************************************************

class NotHashedFingerprint;

// similarity calc, either tversky or tanimoto, returned as a distance
// i.e. 1.0 - sim.
typedef double (NotHashedFingerprint::*pNHDC) ( const NotHashedFingerprint & ) const;
typedef double (NotHashedFingerprint::*pNHTDC)( const NotHashedFingerprint & , float ) const;

// ***********************************************************************

class NotHashedFingerprint : public FingerprintBase {

public :

  NotHashedFingerprint();
  NotHashedFingerprint( const std::string &name );
  // this one from the contents of a string built with get_string_rep.
  NotHashedFingerprint( const std::string &name , const std::string &rep );
  NotHashedFingerprint( const std::string &name ,
			const std::vector<uint32_t> in_nums );
  NotHashedFingerprint( const NotHashedFingerprint &fp  );
  ~NotHashedFingerprint();

  virtual FingerprintBase *clone() const {
    return new NotHashedFingerprint( *this );
  }

  NotHashedFingerprint &operator=( const NotHashedFingerprint &fp );
  // bitwise operators for combining fingerprints.
  NotHashedFingerprint operator&( const NotHashedFingerprint &rhs ) const;
  NotHashedFingerprint operator|( const NotHashedFingerprint &rhs ) const;
  NotHashedFingerprint &operator&=( const NotHashedFingerprint &rhs );
  NotHashedFingerprint &operator|=( const NotHashedFingerprint &rhs );

  double tanimoto( const NotHashedFingerprint &f ) const;
  double tanimoto( const NotHashedFingerprint &f , float threshold ) const;
  double tversky( const NotHashedFingerprint &f ) const;
  double tversky( const NotHashedFingerprint &f , float threshold ) const;

  // binary read and write, possibly to a compressed file
  bool binary_read( gzFile fp , bool byte_swapping );
  void binary_write( gzFile fp ) const;
  void binary_write( std::FILE *fp ) const;
  // likewise for ascii representation
  bool ascii_read( gzFile fp , const std::string &sep );
  void ascii_write( gzFile fp , const std::string &sep ) const;
  void ascii_write( std::FILE *fp , const std::string &sep ) const;

  virtual std::string get_string_rep() const;

  // frag_nums_, suitable for sending over pvm
  virtual char *data_for_pvm( int &num_bytes ) {
    num_bytes = num_frag_nums_ * sizeof( uint32_t );
    return reinterpret_cast<char *>( frag_nums_ );
  }

  static void set_similarity_calc( SIMILARITY_CALC sc );

  // count the number of bits in the fingerprint - quite easy in this case
  int count_bits() const {
    return num_frag_nums_;
  }

  // count the number of bits in common between the fingerprint passed in
  // and this one
  int num_bits_in_common( const NotHashedFingerprint &f ) const;
  int num_bits_in_common( const NotHashedFingerprint &f ,
			  int &num_in_a_not_b  , int &num_in_b_not_a ) const;

  // calculate the distance between this fingerprint and the one passed in
  // using dist_calc_
  double calc_distance( const FingerprintBase &f ) const;
  // calculate the distance between this fingerprint and the one passed
  // in using threshold_dist_calc_.  If the distance is predicted to be above
  // the threshold, return 1.0
  double calc_distance( const FingerprintBase &f , float threshold ) const;

  // calculate the distance between this fingerprint and the one passed in
  // using dist_calc_
  double calc_distance( const HashedFingerprint &f ) const;
  // calculate the distance between this fingerprint and the one passed
  // in using threshold_dist_calc_.  If the distance is predicted to be above
  // the threshold, return 1.0
  double calc_distance( const HashedFingerprint &f , float threshold ) const;
  // calculate the distance between this fingerprint and the one passed in
  // using dist_calc_
  double calc_distance( const NotHashedFingerprint &f ) const;
  // calculate the distance between this fingerprint and the one passed
  // in using threshold_dist_calc_.  If the distance is predicted to be above
  // the threshold, return 1.0
  double calc_distance( const NotHashedFingerprint &f , float threshold ) const;

  virtual void build_from_vector( const std::vector<uint32_t> &in_nums );

protected :

  int num_frag_nums_;
  uint32_t *frag_nums_;

  static pNHDC dist_calc_;
  static pNHTDC threshold_dist_calc_;

  virtual void copy_data( const NotHashedFingerprint &fp );

};

} // end of namespace DAC_FINGERPRINTS

#endif
