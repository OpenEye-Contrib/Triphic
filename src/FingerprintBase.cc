//
// file FingerprintBase.cc
// Dave Cosgrove
// AstraZeneca
// 23rd February 2009
//
// This file has some utility functions used by classes derived from
// FingerprintBase.

#include "ByteSwapper.H"
#include "FileExceptions.H"
#include "FingerprintBase.H"
#include "HashedFingerprint.H"
#include "NotHashedFingerprint.H"
#include "MagicInts.H"

#include <iostream>
#include <sstream>

#include <boost/regex.hpp>

using namespace boost;
using namespace std;

namespace DAC_FINGERPRINTS {

double FingerprintBase::tversky_alpha_ = 0.5F;

// **************************************************************************
FingerprintBase::~FingerprintBase() {

}

// **************************************************************************
// open a possibly compressed fingerprint file for reading.  zlib can read
// an uncompressed file with the same routines as a compressed one.
// Throws a DACLIB::FileReadOpenError if it gets the mood.
void open_fp_file_for_reading( const string &fp_file ,
                               FP_FILE_FORMAT expected_format ,
                               bool &byte_swapping , gzFile &fp ) {

  if( FRAG_NUMS == expected_format || BITSTRINGS == expected_format ) {
    open_fp_file_for_reading( fp_file , fp );
    return;
  }

  fp = gzopen( fp_file.c_str() , "rb" );
  if( !fp ) {
    throw DACLIB::FileReadOpenError( fp_file.c_str() );
  }

  // take the integer off the top to get things set up correctly.
  // The very first integer will be either MAGIC_INT
  // or BUGGERED_MAGIC_INT and indicates whether the machine reading and the
  // machine writing were both in the same bigendian/littleendian format.
  unsigned int file_type;
  gzread( fp , reinterpret_cast<void *>( &file_type ) , sizeof( int ) );
  byte_swapping = false;
  if( FP_MAGIC_INT == file_type || BUGGERED_FP_MAGIC_INT == file_type ) {
    if( expected_format != FLUSH_FPS ) {
      throw FingerprintFileError( fp_file , expected_format ,
                                  "Flush Fingerprints" );
    }
    if( BUGGERED_FP_MAGIC_INT == file_type ) {
      byte_swapping = true;
    }
    // read fingerprint size - back in the day, they used to be handled as
    // unsigned chars, but now they're unsigned ints
    int num_chars_in_fp;
    gzread( fp , reinterpret_cast<char *>( &num_chars_in_fp ) , sizeof( int ) );
    if( byte_swapping ) DACLIB::byte_swapper<int>( num_chars_in_fp );
    int num_ints_in_fp = num_chars_in_fp / sizeof( unsigned int );
    if( num_chars_in_fp % sizeof( unsigned int ) ) {
      ++num_ints_in_fp;
    }
    HashedFingerprint::set_num_ints( num_ints_in_fp );
  } else if( FN_MAGIC_INT == file_type || BUGGERED_FN_MAGIC_INT == file_type ) {
    if( expected_format != BIN_FRAG_NUMS ) {
      throw FingerprintFileError( fp_file , expected_format ,
                                  "Binary Fragment Numbers" );
    }
    if( BUGGERED_FP_MAGIC_INT == file_type ) {
      byte_swapping = true;
    }
  } else {
    cerr << "Unrecognised file format for " << fp_file << endl
         << "Possibly it's an older format no longer supported." << endl;
    exit( 1 );
  }

}

// **************************************************************************
// open the ASCII fingerprint file. It can't know what sort of fingerprints it's
// up against.
void open_fp_file_for_reading( const std::string &fp_file , gzFile &fp ) {

  fp = gzopen( fp_file.c_str() , "r" );
  if( !fp ) {
    throw DACLIB::FileReadOpenError( fp_file.c_str() );
  }
  
}

// **************************************************************************
// open a compressed fingerprint file for writing. Throws a
// DACLIB::FileReadOpenError if it gets the mood.
void open_fp_file_for_writing( const std::string &fp_file , int num_chars_in_fp ,
                               FP_FILE_FORMAT file_format , gzFile &fp ) {

  string fmt = "w";
  if( FLUSH_FPS == file_format || BIN_FRAG_NUMS == file_format ) {
    fmt += "b";
  }
  fp = gzopen( fp_file.c_str() , fmt.c_str() );

  if( !fp ) {
    throw DACLIB::FileWriteOpenError( fp_file.c_str() );
  }

  if( FLUSH_FPS == file_format ) {
    gzwrite( fp , reinterpret_cast<const void *>( &FP_MAGIC_INT ) ,
             sizeof( unsigned int ) );
    // flush file needs extra info
    gzwrite( fp , reinterpret_cast<void *>( &num_chars_in_fp ) ,
             sizeof( int ) );
  } else if( BIN_FRAG_NUMS == file_format ) {
    gzwrite( fp , reinterpret_cast<const void *>( &FN_MAGIC_INT ) ,
             sizeof( unsigned int ) );
  }

}

// **************************************************************************
// open an uncompressed fingerprint file for writing. Throws a
// DACLIB::FileReadOpenError if it gets the mood.
void open_fp_file_for_writing( const std::string &fp_file , int num_chars_in_fp ,
                               FP_FILE_FORMAT file_format , FILE *&fp ) {

  string fmt = "w";
  if( FLUSH_FPS == file_format || BIN_FRAG_NUMS == file_format ) {
    fmt += "b";
  }
  fp = fopen( fp_file.c_str() , fmt.c_str() );

  if( !fp ) {
    throw DACLIB::FileWriteOpenError( fp_file.c_str() );
  }

  if( FLUSH_FPS == file_format ) {
    fwrite( reinterpret_cast<const void *>( &FP_MAGIC_INT ) ,
            sizeof( unsigned int ) , 1 , fp );
    fwrite( reinterpret_cast<void *>( &num_chars_in_fp ) , sizeof( int ) ,
            1 , fp );
  } else if( BIN_FRAG_NUMS == file_format ) {
    fwrite( reinterpret_cast<const void *>( &FN_MAGIC_INT ) ,
            sizeof( unsigned int ) , 1 , fp );
  }

}

// **************************************************************************
// read fp up to next
string read_full_line( gzFile fp ) {

  char c;
  string full_line;
  while( 1 ) {
    c = gzgetc( fp );
    if( gzeof( fp ) || c =='\n' ) {
      break;
    }
    full_line += c;
  }

  return full_line;

}

// **************************************************************************
// read fp up to next
string read_full_line( FILE *fp ) {

  char c;
  string full_line;
  while( 1 ) {
    c = fgetc( fp );
    if( feof( fp ) || c =='\n' ) {
      break;
    }
    full_line += c;
  }

  return full_line;

}

// **************************************************************************
// for use in ascii_read, to convert the separator
string convert_sep_to_new_sep( const string &ins , const string &sep ,
                               const string &new_sep ) {

  ostringstream t;
  ostream_iterator<char,char> oi( t );
  // set it up not to use any special characters
  regex e1( sep , boost::regex::literal );
  regex_replace( oi , ins.begin() , ins.end() , e1 , new_sep ,
                 boost::match_default | boost::format_all );

  return t.str();

}

// ***************************************************************************
void read_flush_fp_file( gzFile &fp , bool byteswapping ,
                         vector<FingerprintBase *> &fps ) {

  HashedFingerprint next_fp( "DUMMY" );

  while( 1 ) {
    if( !next_fp.binary_read( fp , byteswapping ) ) {
      break;
    }
    fps.push_back( new HashedFingerprint( next_fp ) );
  }

}

// ***************************************************************************
void read_bitstrings_file( gzFile &fp , const string &bitstring_separator ,
                           vector<FingerprintBase *> &fps ) {

  HashedFingerprint next_fp( "DUMMY" );

  while( 1 ) {
    if( !next_fp.ascii_read( fp , bitstring_separator ) ) {
      break;
    }
    fps.push_back( new HashedFingerprint( next_fp ) );
  }

}

// ***************************************************************************
void read_frag_nums_file( gzFile &fp ,
                          const string &bitstring_separator ,
                          vector<FingerprintBase *> &fps ) {

  NotHashedFingerprint next_fp( "Dummy" );
  while( 1 ) {
    if( !next_fp.ascii_read( fp , bitstring_separator ) ) {
      break;
    }
    fps.push_back( new NotHashedFingerprint( next_fp ) );
  }

}

// ***************************************************************************
void read_bin_frag_nums_file( gzFile &fp , bool byteswapping ,
                              vector<FingerprintBase *> &fps ) {

  NotHashedFingerprint next_fp( "Dummy" );
  while( 1 ) {
    if( !next_fp.binary_read( fp , byteswapping ) ) {
      break;
    }
    fps.push_back( new NotHashedFingerprint( next_fp ) );
  }

}

// ***************************************************************************
void read_fp_file( gzFile &fp , bool byteswapping ,
                   FP_FILE_FORMAT file_format ,
                   const string &bitstring_separator ,
                   vector<FingerprintBase *> &fps ) {

  switch( file_format ) {
  case FLUSH_FPS :
    read_flush_fp_file( fp , byteswapping , fps );
    break;
  case BITSTRINGS :
    read_bitstrings_file( fp , bitstring_separator , fps );
    break;
  case FRAG_NUMS :
    read_frag_nums_file( fp , bitstring_separator , fps );
    break;
  case BIN_FRAG_NUMS :
    read_bin_frag_nums_file( fp , byteswapping , fps );
    break;
  }

}

// **************************************************************************
void read_fp_file( const string &file , FP_FILE_FORMAT input_format ,
                   const string &bitstring_separator ,
                   vector<FingerprintBase *> &fps ) {

  gzFile gzfp = 0;
  bool byteswapping;
  if( FLUSH_FPS == input_format || BIN_FRAG_NUMS == input_format ) {
    open_fp_file_for_reading( file , input_format , byteswapping , gzfp );
  } else {
    open_fp_file_for_reading( file , gzfp );
  }

  read_fp_file( gzfp , byteswapping , input_format , bitstring_separator ,
                fps );

}

// ***************************************************************************
// probably in a grown-up world, this would be done more elegantly using
// a better object model and virtual functions. It's by far not the
// rate-limiting step in most fingerprint calculations, though.
FingerprintBase *read_next_fp_from_file( gzFile &fp , bool byteswapping ,
                                         FP_FILE_FORMAT file_format ,
                                         const string &bitstring_separator ) {

  FingerprintBase *new_fp = 0;
  try {
    switch( file_format ) {
    case FLUSH_FPS :
      new_fp = new HashedFingerprint;
      if( !new_fp->binary_read( fp , byteswapping ) ) {
        delete new_fp;
        new_fp = 0;
      }
      break;
    case BIN_FRAG_NUMS :
      new_fp = new NotHashedFingerprint;
      if( !new_fp->binary_read( fp , byteswapping ) ) {
        delete new_fp;
        new_fp = 0;
      }
      break;
    case BITSTRINGS :
      new_fp = new HashedFingerprint;
      if( !new_fp->ascii_read( fp , bitstring_separator ) ) {
        delete new_fp;
        new_fp = 0;
      }
      break;
    case FRAG_NUMS :
      new_fp = new NotHashedFingerprint;
      if( !new_fp->ascii_read( fp , bitstring_separator ) ) {
        delete new_fp;
        new_fp = 0;
      }
      break;
    }
  } catch( HashedFingerprintLengthError &e ) {
    cerr << e.what() << endl;
    exit( 1 );
  }

  return new_fp;

}

// **************************************************************************
void read_next_fps_from_file( gzFile &fp , bool byteswapping ,
                              FP_FILE_FORMAT file_format ,
                              const std::string &bitstring_separator ,
                              int chunk_size ,
                              std::vector<FingerprintBase *> &fps ) {

  for( int i = 0 ; i < chunk_size ; ++i ) {
    FingerprintBase *next_fp = 0;
    next_fp = read_next_fp_from_file( fp , byteswapping , file_format ,
                                      bitstring_separator );
    if( !next_fp ) {
      break;
    }
    fps.push_back( next_fp );
  }

}

// **************************************************************************
void read_fps_from_file( gzFile &fp_file , bool byteswapping ,
                         FP_FILE_FORMAT file_format ,
                         const std::string &bitstring_separator ,
                         unsigned int first_fp , unsigned int num_fps ,
                         std::vector<FingerprintBase *> &fps ) {

  // spin through to first fp of interest
  for( unsigned int i = 0 ; i < first_fp ; ++i ) {
    FingerprintBase *fp = read_next_fp_from_file( fp_file , byteswapping ,
                                                  file_format ,
                                                  bitstring_separator );
    if( !fp ) {
      // bad end
      return;
    }
    delete fp;
  }

  for( unsigned int i = 0 ; i < num_fps ; ++i ) {
    FingerprintBase *fp = read_next_fp_from_file( fp_file , byteswapping ,
                                                  file_format ,
                                                  bitstring_separator );
    if( !fp ) {
      break;
    }
    fps.push_back( fp );
  }

}

// **************************************************************************
void decode_format_string( const string &format_string ,
                           FP_FILE_FORMAT &fp_file_format ,
                           bool &binary_file ,
                           string &bitstring_separator ) {

  if( !format_string.empty() && "FLUSH_FPS" == format_string ) {
    fp_file_format = FLUSH_FPS;
    binary_file = true;
  } else if( !format_string.empty() && "BITSTRINGS" == format_string ) {
    fp_file_format = BITSTRINGS;
  } else if( !format_string.empty() && "BIN_FRAG_NUMS" == format_string ) {
    fp_file_format = BIN_FRAG_NUMS;
    binary_file = true;
  } else if( !format_string.empty() && "FRAG_NUMS" == format_string ) {
    fp_file_format = FRAG_NUMS;
    if( bitstring_separator.empty() ) {
      bitstring_separator = " ";
    }
  } else {
    throw FingerprintInputFormatError( format_string );
  }

}

// **************************************************************************
unsigned int count_fps_in_file( const string &filename ,
                                DAC_FINGERPRINTS::FP_FILE_FORMAT fp_format ,
                                const string &bitstring_separator ) {

  gzFile fpfile;
  bool byteswapping = false;
  open_fp_file_for_reading( filename , fp_format , byteswapping , fpfile );

  unsigned int num_fps = 0;
  while( 1 ) {
    FingerprintBase *fp = read_next_fp_from_file( fpfile , byteswapping ,
                                                  fp_format ,
                                                  bitstring_separator );
    if( !fp ) {
      break;
    }
    delete fp;
    ++num_fps;
  }

  return num_fps;

}

// **************************************************************************
void get_fp_names( const std::string &filename ,
                   DAC_FINGERPRINTS::FP_FILE_FORMAT fp_format ,
                   const std::string &bitstring_separator ,
                   std::vector<std::string> &fp_names ) {
  gzFile fpfile;
  bool byteswapping = false;
  open_fp_file_for_reading( filename , fp_format , byteswapping , fpfile );

  while( 1 ) {
    FingerprintBase *fp = read_next_fp_from_file( fpfile , byteswapping ,
                                                  fp_format ,
                                                  bitstring_separator );
    if( !fp ) {
      break;
    }
    fp_names.push_back( fp->get_name() );
    delete fp;
  }

}

// **************************************************************************
FingerprintFileError::FingerprintFileError( const string &file_name ,
                                            FP_FILE_FORMAT expected_format ,
                                            const string &apparent_format ) {

  ostringstream oss;
  string fmt;
  if( FLUSH_FPS == expected_format ) {
    fmt = "Flush Fingerprints";
  } else if( BIN_FRAG_NUMS == expected_format ) {
    fmt = "Binary Fragment Numbers";
  } else {
    fmt = "Unknown Format";
  }

  oss << "Error for file " << file_name
      << ". Expected a " << fmt << " file," << endl
      << "but it appears to be a " << apparent_format << " file.";

  msg_ = oss.str();

}

} // end of DAC_FINGERPRINTS namespace
