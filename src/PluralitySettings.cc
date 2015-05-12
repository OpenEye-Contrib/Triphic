//
// file PluralitySettings.cc
// Dave Cosgrove
// 24th September 2007
//

#include "PluralitySettings.H"

#include <fstream>
#include <iostream>

#include <pvm3.h>

#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

using namespace std;
namespace po = boost::program_options;

namespace DACLIB {
  // In pvm_string_subs.cc
  // pack a C++ string into a pvm buffer
  void pack_string( const string &str );
  // unpack a C++ string from pvm buffer
  void unpack_string( string &str );
  void pack_strings_vector( const vector<string> &strs );
  void unpack_strings_vector( vector<string> &strs );
}

// ***********************************************************************
PluralitySettings::PluralitySettings() : scores_only_( false ) ,
					 comma_output_( false ) ,
					 hits_to_output_str_( "BEST_HITS_ONLY" ) ,
					 hits_to_output_( BEST_HITS_ONLY ) ,
           num_slave_procs_( 0 ) {
}

// ***********************************************************************
PluralitySettings::PluralitySettings( int argc , char **argv ) :
  scores_only_( false ) , comma_output_( false ) ,
  hits_to_output_str_( "BEST_HITS_ONLY" ) ,
  hits_to_output_( BEST_HITS_ONLY ) ,
  num_slave_procs_( 0 ) {

  po::options_description desc( "Allowed Options" );
  build_program_options( desc );

  po::variables_map vm;
  po::store( po::parse_command_line( argc , argv , desc ) , vm );
  po::notify( vm );

  if( argc < 2 || vm.count( "help" ) ) {
    cout << desc << endl;
    exit( 1 );
  }

  read_database_files_file();

  ostringstream oss;
  oss << desc;
  usage_text_ = oss.str();

}

// ***************************************************************************
bool PluralitySettings::operator!() const {

  return false;

}

// ***************************************************************************
void PluralitySettings::print_usage( ostream &os ) const {

  os << usage_text_ << endl;

}

// ***********************************************************************
bool PluralitySettings::set_hits_to_output_from_string() const {

  static vector<pair<string,HITS_TO_OUTPUT> > valid_values;
  if( valid_values.empty() ) {
    valid_values.push_back( make_pair( "BEST_HITS_ONLY" , BEST_HITS_ONLY ) );
    valid_values.push_back( make_pair( "ONE_HIT_ONLY" , ONE_HIT_ONLY ) );
    valid_values.push_back( make_pair( "ALL_HITS" , ALL_HITS ) );
  }

  for( int i = 0 , is = valid_values.size() ; i < is ; ++i ) {
    if( valid_values[i].first == hits_to_output_str_ ) {
      hits_to_output_ = valid_values[i].second;
      return true;
    }
  }

  return false;

}

// **************************************************************************
void PluralitySettings::pack_contents_into_pvm_buffer() {

  DACLIB::pack_string( pphore_file_ );
  DACLIB::pack_strings_vector( db_files_ );
  DACLIB::pack_string( smarts_file_ );
  DACLIB::pack_string( points_file_ );
  DACLIB::pack_string( protein_file_ );
  DACLIB::pack_string( subset_file_ );
  DACLIB::pack_string( not_smarts_file_ );
  DACLIB::pack_strings_vector( grid_vol_files_ );
  int i = hits_to_output_;
  pvm_pkint( &i , 1 ,1 );

}

// **************************************************************************
void PluralitySettings::unpack_contents_from_pvm_buffer() {

  DACLIB::unpack_string( pphore_file_ );
  DACLIB::unpack_strings_vector( db_files_ );
  DACLIB::unpack_string( smarts_file_ );
  DACLIB::unpack_string( points_file_ );
  DACLIB::unpack_string( protein_file_ );
  DACLIB::unpack_string( subset_file_ );
  DACLIB::unpack_string( not_smarts_file_ );
  DACLIB::unpack_strings_vector( grid_vol_files_ );

  int i;
  pvm_upkint( &i , 1 , 1 );
  hits_to_output_ = HITS_TO_OUTPUT( i );

}

// **************************************************************************
void PluralitySettings::build_program_options( po::options_description &desc ) {

  desc.add_options()
    ( "help" , "Produce this help text" )
    ( "db-file,D" , po::value<vector<string> >( &db_files_ ) ,
      "Database file name. Multiple instances allowed." )
      ( "db-files-file" , po::value<string>( &db_files_file_ ) ,
        "File containing list of database files to search." )
    ( "pphore-file,Q" , po::value<string>( &pphore_file_ ) ,
      "Pharmacophore query file name." )
    ( "smarts-file,S" , po::value<string>( &smarts_file_ ) ,
      "SMARTS file name." )
    ( "points-file,P" , po::value<string>( &points_file_ ) ,
      "Points file name." )
    ( "output-file,O" , po::value<string>( &output_file_ ) ,
      "Output file name." )
    ( "protein-file" , po::value<string>( &protein_file_ ) ,
      "File containing protein active site." )
    ( "subset-file" , po::value<string>( &subset_file_ ) ,
      "Names of molecules to be searched in the database file." )
    ( "grid-score-file,G" , po::value<vector<string> >( &grid_vol_files_ ) ,
      "File containing grid for calculating overlap volume with hit. Multiple files allowed, but only 1 per option." )
    ( "scores-only" , po::value<bool>( &scores_only_ )->zero_tokens() ,
      "Output just the scores for the hits, no molecule file." )
    ( "hits-to-output" , po::value<string>( &hits_to_output_str_ )->default_value( "BEST_HITS_ONLY" ) ,
      "Which hits to output for each molecule, best subset (default), first one found (ONE_HIT_ONLY), or all (ALL_HITS)." )
    ( "comma-output" , po::value<bool>( &comma_output_ )->zero_tokens() ,
      "Puts comma-separated list in .scores file." )
    ( "dont-give-me" , po::value<string>( &not_smarts_file_ ) ,
      "File of SMARTS that the hits must not match." )
    ( "i-already-know-about" , po::value<string>( &not_smarts_file_ ) ,
      "File of SMARTS that the hits must not match." )
    ( "slave-name" , po::value<string>( &slave_name_ ) ,
      "PVM slave name." )
    ( "num-slave-procs" , po::value<int>( &num_slave_procs_ )->default_value( 0 ) ,
      "Number of PVM slaves to launch." )
      ( "pvm-hosts-file" , po::value<string>( &pvm_hosts_file_ ) ,
        "File giving number of PVM slaves to run on each node." );

}

// **************************************************************************
void PluralitySettings::read_database_files_file() {

  if( db_files_file_.empty() ) {
    return;
  }

  ifstream iff( db_files_file_.c_str() );
  string next_file;
  while( true ) {
    iff >> next_file;
    if( !iff.good() || iff.eof() ) {
      break;
    }
    db_files_.push_back( next_file );
  }

}
