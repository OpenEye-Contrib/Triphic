//
// file FileExceptions.H
// David Cosgrove
// AstraZeneca
// 1st March 2007
//
// Exceptions for use when open files for reading or writing fail.

#ifndef DAC_FILE_EXCEPTIONS
#define DAC_FILE_EXCEPTIONS

#include <string>

namespace DACLIB {

  // ***************************************************************************
  class FileReadOpenError {

  public :
    explicit FileReadOpenError( const char *filename ) {
      msg_ = std::string( "Couldn't open file " ) + std::string( filename ) +
	std::string( " for reading." );
    }
    virtual ~FileReadOpenError() {}
    virtual const char *what() {
      return msg_.c_str();
    }
  private : 
    std::string msg_;

  };

  // ***************************************************************************
  class FileWriteOpenError {

  public :
    explicit FileWriteOpenError( const char *filename ) {
      msg_ = std::string( "Couldn't open file " ) + std::string( filename ) +
	std::string( " for writing." );
    }
    virtual ~FileWriteOpenError() {}
    virtual const char *what() {
      return msg_.c_str();
    }
  private : 
    std::string msg_;

  };

} // end of namespace


#endif
