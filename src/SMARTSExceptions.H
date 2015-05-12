//
// file SMARTSExceptions.H
// David Cosgrove
// AstraZeneca
// 2nd October 2007
//
// Exceptions to use when things go wrong with a SMARTS file or SMARTS

#ifndef DAC_SMARTS_EXCEPTIONS
#define DAC_SMARTS_EXCEPTIONS

#include <stdexcept>
#include <string>

#include <boost/lexical_cast.hpp>

// **************************************************************************

namespace DACLIB {

  class SMARTSDefnError : public std::runtime_error {
  public :
    explicit SMARTSDefnError( const char *msg = "SMARTS definition error." ) :
      std::runtime_error( msg ) {}
  };

  class SMARTSFileError {

  public :
    explicit SMARTSFileError( const std::string &msg , int line_num = -1 ) {
      msg_ = std::string( "Bad SMARTS line " );
      if( -1 != line_num )
	msg_ += boost::lexical_cast<std::string>( line_num ) + std::string( " " );
      msg_ += std::string( " : \'" ) + msg + std::string( "\'" );
    }
    virtual ~SMARTSFileError() {}
    virtual const char *what() {
      return msg_.c_str();
    }
			      
  private :
    std::string msg_;

  };

  class SMARTSSubDefnError {

  public :
    explicit SMARTSSubDefnError( const std::string &label ) {
      msg_ = "SMARTS name " + label + " not defined.";
    }
    explicit SMARTSSubDefnError( const std::string &smt ,
				 const std::string &label ) {
      msg_ = "Can't expand vector bindings in SMARTS string "
	+ smt + " label " + label + ".";
    }
    explicit SMARTSSubDefnError( int line_num , const std::string &input_line ,
				 const std::string &rep_name ,
				 const std::string &this_def ,
				 const std::string &prev_def ) {
      msg_ = std::string( "ERROR reading line " ) +
	boost::lexical_cast<std::string>( line_num ) + std::string( " : " ) +
	input_line + std::string( "\nSMARTS sub-definition label " ) + rep_name +
	std::string( " already used.\n" ) +
	std::string( "This definition is " ) + this_def +
	std::string( ".\nPrevious definition is " ) + prev_def +
	std::string( ".\nPlease check your file. " ) +
	std::string( "Perhaps you've read the same file twice.\n" );
    }
    virtual ~SMARTSSubDefnError() {}
    virtual const char *what() {
      return msg_.c_str();
    }
			      
  private :
    std::string msg_;

  };

}

#endif
