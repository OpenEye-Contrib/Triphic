//
// file PluralitySettings.cc
// David Cosgrove
// 24th September 2007
//

#include "PluralitySettings.H"

#include <fstream>
#include <iostream>

#include <mpi.h>

#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

using namespace std;
namespace po = boost::program_options;

namespace DACLIB {
// In mpi_string_subs.cc
void mpi_send_string( const string &str , int dest_rank );
void mpi_send_strings_vector( const vector<string> &strs , int dest_rank );
void mpi_rec_string( int source_rank , string &str );
void mpi_rec_strings_vector( int source_rank , vector<string> &strs );
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
void PluralitySettings::send_contents_via_mpi( int dest_slave ) {

  DACLIB::mpi_send_string( pphore_file_ , dest_slave );
  DACLIB::mpi_send_strings_vector( db_files_ , dest_slave );
  DACLIB::mpi_send_string( smarts_file_ , dest_slave );
  DACLIB::mpi_send_string( points_file_ , dest_slave );
  DACLIB::mpi_send_string( protein_file_ , dest_slave );
  DACLIB::mpi_send_string( subset_file_ , dest_slave );
  DACLIB::mpi_send_string( not_smarts_file_ , dest_slave );
  DACLIB::mpi_send_strings_vector( grid_vol_files_ , dest_slave );
  int i = hits_to_output_;
  MPI_Send( &i , 1 , MPI_INT , dest_slave , 0 , MPI_COMM_WORLD );

}

// **************************************************************************
void PluralitySettings::receive_contents_via_mpi() {

  DACLIB::mpi_rec_string( 0 , pphore_file_ );
  DACLIB::mpi_rec_strings_vector( 0 , db_files_ );
  DACLIB::mpi_rec_string( 0 , smarts_file_ );
  DACLIB::mpi_rec_string( 0 , points_file_ );
  DACLIB::mpi_rec_string( 0 , protein_file_ );
  DACLIB::mpi_rec_string( 0 , subset_file_ );
  DACLIB::mpi_rec_string( 0 , not_smarts_file_ );
  DACLIB::mpi_rec_strings_vector( 0 , grid_vol_files_ );

  int i;
  MPI_Recv( &i , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
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
