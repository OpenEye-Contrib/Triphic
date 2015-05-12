//
// file FingerprintBase.H
// Dave Cosgrove
// AstraZeneca
// 29th Jan 2009
//
// This is an abstract base class for fingerprints, defining the interface
// that fingerprints must observe.
// It's not a true interface class, because it has finger_name_.  All fingerprints
// will have a name.

#ifndef DAC_FINGER_PRINT_BASE
#define DAC_FINGER_PRINT_BASE

#include <cstdio>
#include <iosfwd>
#include <string>
#include <vector>

#include <zlib.h> // compressed I/O

// used by derived classes
namespace DAC_FINGERPRINTS {

  typedef enum { TANIMOTO , TVERSKY } SIMILARITY_CALC;
  typedef enum { FLUSH_FPS , BIN_FRAG_NUMS , BITSTRINGS ,
		 FRAG_NUMS } FP_FILE_FORMAT;
  typedef enum { NO_HASH , OLD_DENSE , NEW_SPARSE } HASH_METHOD;
  typedef enum { ALFI , ECFI , FCFI , FOYFI , LIBFI } CREATION_TYPE;

  class HashedFingerprint;
  class NotHashedFingerprint;

  // ************************************************************************

  class FingerprintBase {

  public:

    // default c'tor
    FingerprintBase() : finger_name_( "Dummy" ) {}
    FingerprintBase( const std::string &name ) : finger_name_( name ) {}
    // c'tor from the product of get_string_rep
    FingerprintBase( const std::string &name , const std::string &rep );

    // d'tor
    virtual ~FingerprintBase();

    virtual FingerprintBase *clone() const = 0;

    // count the number of set bits in the fingerprint
    virtual int count_bits() const = 0;

    // set name
    void set_name( std::string new_name ) {
      finger_name_ = new_name;
    }
    // return name
    std::string get_name() const {
      return finger_name_;
    }

    // binary read and write, possibly to a compressed file
    virtual bool binary_read( gzFile fp , bool byte_swapping ) = 0;
    virtual void binary_write( gzFile fp ) const = 0;
    virtual void binary_write( std::FILE *fp ) const = 0;
    // likewise for ascii representation
    virtual bool ascii_read( gzFile fp , const std::string &sep ) = 0;
    virtual void ascii_write( gzFile fp , const std::string &sep ) const = 0;
    virtual void ascii_write( std::FILE *fp , const std::string &sep ) const = 0;
    
    virtual std::string get_string_rep() const = 0;

    // fp data suitable for sending over pvm
    virtual char *data_for_pvm( int &num_bytes ) = 0;

    static void set_tversky_alpha( double a ) { tversky_alpha_ = a; }
    static double get_tversky_alpha() { return tversky_alpha_; }

    // calculate the distance between this fingerprint and the one passed in
    // using dist_calc_
    virtual double calc_distance( const FingerprintBase &f ) const = 0;
    // calculate the distance between this fingerprint and the one passed
    // in using threshold_dist_calc_.  If the distance is predicted to be above
    // the threshold, return 1.0
    virtual double calc_distance( const FingerprintBase &f ,
				  float threshold ) const = 0;
    // calculate the distance between this fingerprint and the one passed in
    // using dist_calc_
    virtual double calc_distance( const HashedFingerprint &f ) const = 0;
    // calculate the distance between this fingerprint and the one passed
    // in using threshold_dist_calc_.  If the distance is predicted to be above
    // the threshold, return 1.0
    virtual double calc_distance( const HashedFingerprint &f ,
				  float threshold ) const = 0;
    // calculate the distance between this fingerprint and the one passed in
    // using dist_calc_
    virtual double calc_distance( const NotHashedFingerprint &f ) const = 0;
    // calculate the distance between this fingerprint and the one passed
    // in using threshold_dist_calc_.  If the distance is predicted to be above
    // the threshold, return 1.0
    virtual double calc_distance( const NotHashedFingerprint &f ,
				  float threshold ) const = 0;

  protected:
    
    std::string finger_name_;

    static double tversky_alpha_;

    virtual void copy_data( const FingerprintBase &fp ) {
      finger_name_ = fp.finger_name_;
    }

  private:

    // we won't be allowing assignments of the base.
    FingerprintBase &operator=( const FingerprintBase &b );

  };

  // open a possibly compressed fingerprint file for reading.  zlib can read
  // an uncompressed file with the same routines as a compressed one. Throws a
  // DACLIB::FileReadOpenError if it gets the mood.
  // The file will be binary, and contain either HashedFingerprints
  // NotHashedFingerprints, but will cope with an ASCII format by calling
  // the 2nd overloaded version.
  void open_fp_file_for_reading( const std::string &flush_file ,
                                 FP_FILE_FORMAT expected_format ,
                                 bool &byte_swapping , gzFile &fp );
  // open the ASCII equivalent.
  void open_fp_file_for_reading( const std::string &fp_file , gzFile &fp );

  // open a compressed file for writing
  void open_fp_file_for_writing( const std::string &flush_file ,
                                 int num_chars_in_fp ,
                                 FP_FILE_FORMAT output_format , gzFile &fp );
  // open an uncompressed file for writing
  void open_fp_file_for_writing( const std::string &flush_file ,
                                 int num_chars_in_fp ,
                                 FP_FILE_FORMAT output_format , std::FILE *&fp );

  // read fp up to next
  std::string read_full_line( gzFile fp );
  std::string read_full_line( std::FILE *fp );

  // for use in ascii_read, to convert the separator to a different one
  std::string convert_sep_to_new_sep( const std::string &ins ,
                                      const std::string &sep ,
                                      const std::string &new_sep );

  void read_fp_file( gzFile &fp , bool byteswapping ,
                     FP_FILE_FORMAT file_format ,
                     const std::string &bitstring_separator ,
                     std::vector<FingerprintBase *> &fps );

  void read_fp_file( const std::string &file , FP_FILE_FORMAT input_format ,
                     const std::string &bitstring_separator ,
                     std::vector<FingerprintBase *> &fps );

  FingerprintBase *read_next_fp_from_file( gzFile &fp , bool byteswapping ,
                                           FP_FILE_FORMAT file_format ,
                                           const std::string &bitstring_separator );
  void read_next_fps_from_file( gzFile &fp , bool byteswapping ,
                                FP_FILE_FORMAT file_format ,
                                const std::string &bitstring_separator ,
                                int chunk_size ,
                                std::vector<FingerprintBase *> &fps );

  void read_fps_from_file( gzFile &fp_file , bool byteswapping ,
                           FP_FILE_FORMAT file_format ,
                           const std::string &bitstring_separator ,
                           unsigned int first_fp , unsigned int num_fps ,
                           std::vector<FingerprintBase *> &fps );

  void decode_format_string( const std::string &format_string ,
                             FP_FILE_FORMAT &fp_file_format ,
                             bool &binary_file ,
                             std::string &bitstring_separator );

  unsigned int count_fps_in_file( const std::string &filename ,
                                  DAC_FINGERPRINTS::FP_FILE_FORMAT fp_format ,
                                  const std::string &bitstring_separator );

  void get_fp_names( const std::string &filename ,
                     DAC_FINGERPRINTS::FP_FILE_FORMAT fp_format ,
                     const std::string &bitstring_separator ,
                     std::vector<std::string> &fp_names );

  // ***********************************************************************
  class FingerprintFileError {
  public :
    explicit FingerprintFileError() : msg_( "EndOfFile" ) {}
    explicit FingerprintFileError( const std::string &file_name ,
				   FP_FILE_FORMAT expected_format ,
				   const std::string &apparent_format );
    const char *what() const { return msg_.c_str(); }
  private :
    std::string msg_;
  };

  // **********************************************************************
  class IncompatibleFingerprintError {
  public :
    explicit IncompatibleFingerprintError( const std::string &oper ) {
      msg_ = std::string( "Incompatible fingerprint types for operation " )
	+ oper + ".";
    }
    const char *what() const { return msg_.c_str(); }
  private :
    std::string msg_;
  };

  // **********************************************************************
  class FingerprintInputFormatError {
  public :
    explicit FingerprintInputFormatError( const std::string &bad_format_string ) :
      msg_( std::string( "Bad input format specifier string : " ) +
	    bad_format_string ) {}
    const char *what() {
      return msg_.c_str();
    }
  private :
    std::string msg_;
  };

  // **********************************************************************
  class FingerprintDistCalcError {
  public :
    explicit FingerprintDistCalcError( const std::string &bad_string ) :
      msg_( std::string( "Bad distance calculation specifier string : " ) +
	    bad_string ) {}
    const char *what() {
      return msg_.c_str();
    }
  private :
    std::string msg_;
  };

} // end of namespace DAC_FINGERPRINTS

#endif
