//
// file read_smarts_file.cc
// David Cosgrove
// AstraZeneca
// 19th September 2005
//
// This file contains a function that reads a SMARTS file in PWK's format
// and splits it into full SMARTS and vector bindings. Throws
// DACLIB::SMARTSFileError or DACLIB::FileReadOpenError exceptions if
// there's an error.
// 

#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <oechem.h>

#include <boost/tuple/tuple.hpp>

#include "FileExceptions.H"
#include "SMARTSExceptions.H"

using namespace std;
using namespace OEChem;

namespace DACLIB {

  // ****************************************************************************
  void read_smarts_from_string( const char *smarts_string ,
                                vector<boost::tuple<string,string,int> > &smarts ,
                                vector<pair<int,string> > &file_lines ) {

    istringstream smarts_stream( smarts_string );

    string input_line;
    int    line_count = 0;

    while( 1 ) {

      input_line.clear();
      getline( smarts_stream , input_line , '\n' );
      // if the last line ended without '\n', eof is returned even though
      // something was read

      if( input_line.empty() &&
          ( smarts_stream.eof() || !smarts_stream.good() ) ) {
        smarts_stream.clear();
        break;
      }
      line_count++;
      if( input_line.empty() || '#' == input_line[0] ) {
        continue;
      }
      string smarts_name , smarts_def;
      int  flag1 , flag2;
      istringstream is( input_line );
      is >> smarts_name >> smarts_def >> flag1 >> flag2;
      if( is.fail() ) {
        smarts.clear();
        file_lines.clear();
        throw( DACLIB::SMARTSFileError( input_line , line_count ) );
      }

      smarts.push_back( boost::make_tuple( smarts_name , smarts_def , flag2 ) );
      file_lines.push_back( make_pair( line_count , input_line ) );

    }

  }

  // ****************************************************************************
  void read_smarts_from_string( const char *smarts_string ,
                                vector<pair<string,string> > &input_smarts ,
                                vector<pair<string,string> > &smarts_sub_defn ) {

    vector<boost::tuple<string,string,int> > smarts;
    vector<pair<int,string> > file_lines;
    read_smarts_from_string( smarts_string , smarts , file_lines );

    for( int i = 0 , is = smarts.size() ; i < is ; ++i ) {

      string smarts_name = smarts[i].get<0>();
      string smarts_def = smarts[i].get<1>();
      int flag2 = smarts[i].get<2>();

      if( 1 == flag2 ) {
        // it's a full definition, store it in the smarts vector
        input_smarts.push_back( make_pair( smarts_name , smarts_def ) );
      }

      // otherwise it's a vector binding or sub-definition - full definitions can
      // be vector bindings also
      vector<pair<string,string> >::iterator q;
      for( q = smarts_sub_defn.begin() ; q != smarts_sub_defn.end() ; ++q ) {
        if( q->first == smarts_name ) {
          throw( DACLIB::SMARTSSubDefnError( file_lines[i].first ,
                                             file_lines[i].second ,
                                             smarts_name , smarts_def ,
                                             q->second ) );
        }
      }

      smarts_sub_defn.push_back( make_pair( string( smarts_name ) ,
                                            string( smarts_def ) ) );

    }

  }

  // ****************************************************************************
  void slurp_smarts_file( const string &smarts_file ,
                          vector<char> &file_contents ) {

    ifstream ifs;
    ifs.open( smarts_file.c_str() );
    if( !ifs.good() ) {
      throw( DACLIB::FileReadOpenError( smarts_file.c_str() ) );
    }
    copy( istreambuf_iterator<char>( ifs ) , istreambuf_iterator<char>() ,
          back_inserter( file_contents ) );
    file_contents.push_back( '\0' );

  }

  // ****************************************************************************
  void read_smarts_file( const string &smarts_file ,
                         vector<pair<string,string> > &input_smarts ,
                         vector<pair<string,string> > &smarts_sub_defn ) {

    if( smarts_file.empty() ) {
      throw( DACLIB::FileReadOpenError( "No SMARTS file specified." ) );
    }

    vector<char> file_contents;
    slurp_smarts_file( smarts_file , file_contents );
    if( file_contents.empty() ) {
      return;
    }

    read_smarts_from_string( &file_contents[0] , input_smarts ,
                             smarts_sub_defn );

  }

  // ****************************************************************************
  // overloaded version, using boost tuples to return the SMARTS
  void read_smarts_file( const string &smarts_file ,
                         vector<boost::tuple<string,string,int> > &smarts ) {

    if( smarts_file.empty() ) {
      throw( DACLIB::FileReadOpenError( "No SMARTS file specified." ) );
    }

    vector<char> file_contents;
    slurp_smarts_file( smarts_file , file_contents );
    if( file_contents.empty() ) {
      return;
    }

    vector<pair<int,string> > file_lines;
    read_smarts_from_string( &file_contents[0] , smarts , file_lines );

  }

  // ****************************************************************************
  // do the expansion of any vector bindings
  void expand_smarts_defs( const vector<pair<string,string> > &input_smarts ,
                           vector<pair<string,string> > &smarts_sub_defn ,
                           vector<pair<string,string> > &exp_smarts ) {

    string exp_smts;
    for( int i = 0 , is = input_smarts.size() ; i < is ; ++i ) {
      exp_smts = input_smarts[i].second;
      if( string::npos != exp_smts.find( "$" ) &&
          !SmartsLexReplace( exp_smts , smarts_sub_defn ) ) {
        exp_smarts.clear();
        throw( DACLIB::SMARTSSubDefnError( input_smarts[i].second ,
                                           input_smarts[i].first  ) );
      }
      exp_smarts.push_back( make_pair( input_smarts[i].first , exp_smts ) );
    }

  }

} // end of namespace DACLIB
