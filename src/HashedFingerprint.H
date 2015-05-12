//
// file HashedFingerprint.H
// Dave Cosgrove
// AstraZeneca
// 29th January 2009
//
// This is the class for hashed fingerprints, which was formally class
// Fingerprint. It derives from FingerprintBase, and stores bitstrings.

#ifndef DAC_HASHED_FINGERPRINT
#define DAC_HASHED_FINGERPRINT

#include <vector>

#include "FingerprintBase.H"
#include "MagicInts.H"

namespace DAC_FINGERPRINTS {

// ************************************************************************

class HashedFingerprint;

// similarity calc, either tversky or tanimoto, returned as a distance
// i.e. 1.0 - sim.
typedef double (HashedFingerprint::*pHDC) ( const HashedFingerprint & ) const;
typedef double (HashedFingerprint::*pHTDC)( const HashedFingerprint & , float ) const;

// *****************************************************************************

class HashedFingerprint : public FingerprintBase {

public:

  HashedFingerprint();
  HashedFingerprint( const std::string &name );
  // this one from the contents of a string built with get_string_rep.
  HashedFingerprint( const std::string &name , const std::string &rep );

  // new_ints is copied, so a memory leak will ensue if it isn't freed by
  // the caller
  HashedFingerprint( const std::string &name , unsigned int *new_ints );
  HashedFingerprint( const HashedFingerprint &fp  );
  virtual ~HashedFingerprint();
  
  virtual FingerprintBase *clone() const {
    return new HashedFingerprint( *this );
  }

  HashedFingerprint &operator=( const HashedFingerprint &fp );
  // bitwise operators for combining fingerprints - SWIG didn't
  // like the friend way of doing this.
  HashedFingerprint operator&( const HashedFingerprint &rhs ) const;
  HashedFingerprint operator|( const HashedFingerprint &rhs ) const;
  HashedFingerprint &operator&=( const HashedFingerprint &rhs );
  HashedFingerprint &operator|=( const HashedFingerprint &rhs );

  // work out how many ints are required to accommodate the number of
  // bits requested.
  unsigned int calc_num_ints_req( int num_bits ) const;

  // make empty fingerprint of the required size.
  void make_zero_fp( int num_bits );

  double tanimoto( const HashedFingerprint &f ) const;
  double tanimoto( const HashedFingerprint &f , float threshold ) const;
  double tversky( const HashedFingerprint &f ) const;
  double tversky( const HashedFingerprint &f , float threshold ) const;

  virtual std::string get_string_rep() const;

  // frag_nums_, suitable for sending over pvm
  virtual char *data_for_pvm( int &num_bytes ) {
    num_bytes = num_ints_ * sizeof( unsigned int );
    return reinterpret_cast<char *>( finger_bits_ );
  }

  int num_bits_in_common( const HashedFingerprint &fp ) const;
  int num_bits_in_common( const HashedFingerprint &fp ,
                          int &num_in_a_not_b , int &num_in_b_not_a ) const;
  int num_set_in_this_and_not_in_2( const HashedFingerprint &fp1 ) const;

  static void set_similarity_calc( SIMILARITY_CALC sc );

  static void set_num_ints( unsigned int new_val ) { num_ints_ = new_val; }
  static unsigned int num_ints() { return num_ints_; }
  int num_bits_set() const { return num_bits_set_; }
  const unsigned int *get_finger_bits() const { return finger_bits_; }

  // count the number of set bits in the fingerprint
  int count_bits() const;

  // calculate the distance between this fingerprint and the one passed in
  // using dist_calc_
  virtual double calc_distance( const FingerprintBase &f ) const;
  // calculate the distance between this fingerprint and the one passed
  // in using threshold_dist_calc_.  If the distance is predicted to be above
  // the threshold, return 1.0
  virtual double calc_distance( const FingerprintBase &f , float threshold ) const;

  // calculate the distance between this fingerprint and the one passed in
  // using dist_calc_
  virtual double calc_distance( const HashedFingerprint &f ) const;
  // calculate the distance between this fingerprint and the one passed
  // in using threshold_dist_calc_.  If the distance is predicted to be above
  // the threshold, return 1.0
  virtual double calc_distance( const HashedFingerprint &f , float threshold ) const;
  // calculate the distance between this fingerprint and the one passed in
  // using dist_calc_
  virtual double calc_distance( const NotHashedFingerprint &f ) const;
  // calculate the distance between this fingerprint and the one passed
  // in using threshold_dist_calc_.  If the distance is predicted to be above
  // the threshold, return 1.0
  virtual double calc_distance( const NotHashedFingerprint &f ,
                                float threshold ) const;

  // binary read and write, possibly to a compressed file
  bool binary_read( gzFile fp , bool byte_swapping );
  void binary_write( gzFile fp ) const;
  void binary_write( std::FILE *fp ) const;
  // likewise for ascii representation
  bool ascii_read( gzFile fp , const std::string &sep );
  void ascii_write( gzFile fp , const std::string &sep ) const;
  void ascii_write( std::FILE *fp , const std::string &sep ) const;

protected:

  static unsigned int num_ints_; /* the number of chars in the fingerprint. static
             because we can't do anything with 2 fps of
             different lengths in the same run */
  unsigned int *finger_bits_; /* the unsigned ints that hold the bits in the
          fingerprint */
  mutable int      num_bits_set_; // the number of set bits in the fingerprint

  static pHDC dist_calc_;
  static pHTDC threshold_dist_calc_;

  virtual void copy_data( const HashedFingerprint &fp );

  void build_fp_from_bitstring( const std::string &name ,
                                const std::string &bitstring );

};

int count_bits_set( unsigned int *bits , int num_ints );

// exception thrown when an attempt is made to change the fingerprint
// length except through set_num_chars().
class HashedFingerprintLengthError {
public :
  explicit HashedFingerprintLengthError( int new_len , int old_len );
  const char *what() const { return msg_.c_str(); }
private :
  std::string msg_;
};

} // end of namespace DAC_FINGERPRINTS

#endif

