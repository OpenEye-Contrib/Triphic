//
// file ParseSMARTSXML.cc
// David Cosgrove
// AstraZeneca
// 1st May 2012
//

#include "FileExceptions.H"
#include "ParseSMARTSXML.H"

#include <fstream>
#include <iostream>
#include <sstream>

#include <boost/foreach.hpp>
// hints for using boost::ptree for parsing XML came from
// http://akrzemi1.wordpress.com/2011/07/13/parsing-xml-with-boost
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

using namespace std;

// ****************************************************************************
ParseSMARTSXML::ParseSMARTSXML() {

}

// ****************************************************************************
//throws a DACLIB::FileReadOpenError if it can't open the file
void ParseSMARTSXML::parse_file( const string &xml_file ,
                                 vector<pair<string,string> > &smts ,
                                 vector<pair<string,string> > &vbs ) {

  ifstream ifs( xml_file.c_str() );
  if( !ifs || !ifs.good() ) {
    throw DACLIB::FileReadOpenError( xml_file.c_str() );
  }

  parse_stream( ifs , smts , vbs );

}

// ****************************************************************************
void ParseSMARTSXML::parse_string( const string &xml ,
                                   vector<pair<string,string> > &smts ,
                                   vector<pair<string,string> > &vbs ) {

#ifdef NOTYET
  cout << "Reading XML string : " << xml << endl;
#endif

  istringstream iss( xml );
  parse_stream( iss , smts , vbs );

}

// ****************************************************************************
// this is the one that does the work
void ParseSMARTSXML::parse_stream( istream &is ,
                                   vector<pair<string,string> > &smts ,
                                   vector<pair<string,string> > &vbs ) {

  // populate tree structure pt
  using boost::property_tree::ptree;
  ptree pt;
  boost::property_tree::read_xml( is, pt );

  BOOST_FOREACH( ptree::value_type const &v, pt.get_child( "features" ) ) {
    if( v.first == "smarts") {
      string name = v.second.get<string>( "name" );
      string value = v.second.get<string>( "value" );
      // all SMARTS are used as vector bindings
      vbs.push_back( make_pair( name , value ) );
      // this is just a full SMARTS defn.
      if( string( "NOT_FOUND") == v.second.get( "vector_binding" , "NOT_FOUND" ) ) {
        smts.push_back( make_pair( name , value ) );
        // cout << "full SMARTS : " << name << " : " << value << endl;
      }
    }
  }

}

