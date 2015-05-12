//
// file read_site_from_file.cc
// David Cosgrove
// AstraZeneca
// 30th July 2007
//
// Reads a Sites file written, for example, by Grappel, deduces the type
// and returns an object of the correct type. At present, that will be
// either A SinglePPhoreSite or MultiPPhoreSite, but we're returning
// pointers of type BasePPhoreSite. BasePPhoreSite is an abstract base class.
// Returns a zero pointer if there's an error. Maybe it should throw an exception?

#include "SinglePPhoreSite.H"

#include <istream>

using namespace std;

// ****************************************************************************
BasePPhoreSite *read_site_from_file( istream &is ) {

  string type_name;
  is >> type_name;
  if( is.fail() || is.eof() )
    return 0;

  BasePPhoreSite *new_site = 0;

  if( string( "SinglePPhoreSite" ) == type_name ) {
    new_site = new SinglePPhoreSite;
  }

  if( new_site ) {
    new_site->read_from_stream( is );
  }
  return new_site;

}
