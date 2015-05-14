//
// file LoobSettings.cc
// Dave Cosgrove
// AstraZeneca
// 2nd October 2007
//

#include <algorithm>
#include <iostream>

#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "stddefs.H"
#include "LoobSettings.H"

namespace DACLIB {
  void split_filename( const std::string &filename ,
		       std::string &file_root , std::string &file_ext );
}

using namespace std;
namespace po = boost::program_options;

// ****************************************************************************
LoobSettings::LoobSettings( int argc , char **argv ) :
  pairs_( false ) ,
  triplets_( false ) ,
  quadruplets_( false ) ,
  single_conf_mols_( false ) ,
  dont_compress_fps_( false ) ,
  chiral_fps_( false ) ,
  ascii_fps_( false ) ,
  compact_fps_( false ) ,
  labels_not_bits_( false ) ,
  bit_labels_to_ascii_file_( false ) ,
  compressed_counts_cutoff_( 1 ) {
  
  po::options_description desc( "Allowed Options" );
  build_program_options( desc );

  po::variables_map vm;
  po::store( po::parse_command_line( argc , argv , desc ) , vm );
  po::notify( vm );

  if( argc < 2 || vm.count( "help" ) ) {
    cout << desc << endl;
    exit( 1 );
  }

  polish_args();

  ostringstream oss;
  oss << desc;
  usage_text_ = oss.str();

}

// ****************************************************************************
void LoobSettings::build_program_options( boost::program_options::options_description &desc ) {

  desc.add_options()
    ( "help" , "Produce this help text" )
    ( "db-file,D" , po::value<string>( &db_file_ ) ,
      "Database file name." )
    ( "smarts-file,S" , po::value<string>( &smarts_file_ ) ,
      "SMARTS file name" )
    ( "points-file,P" , po::value<string>( &points_file_ ) ,
      "Points file name" )
    ( "subset-file" , po::value<string>( &subset_file_ ) ,
      "File containing subset of db-file that is to be processed." )
    ( "pairs" , po::value<bool>( &pairs_ )->zero_tokens() ,
      "Output bits for pharmacophore pairs." )
    ( "triplets" , po::value<bool>( &triplets_ )->zero_tokens() ,
      "Output bits for pharmacophore triplets (triplets only is default)." )
    ( "quadruplets" , po::value<bool>( &quadruplets_ )->zero_tokens() ,
      "Output bits for pharmacophore quadruplets." )
    ( "single-conf-mols" , po::value<bool>( &single_conf_mols_ )->zero_tokens() ,
      "Process each conformation as a separate molecule." )
    ( "chiral-fps" , po::value<bool>( &chiral_fps_ )->zero_tokens() ,
      "For quadruplet features, whether to take account of chirality." )
    ( "bit-separator" , po::value<string>( &bit_separator_ )->default_value( "" ) ,
      "Separator between bits in ASCII files." )
    ( "distance-bound,B" , po::value<vector<float> >( &dist_bounds_ ) ,
      "Distance bin boundaries - default 4.5, 7.0, 10.0, 14.0, 19.0, 24.0." )
    ( "dont-compress-fps" , po::value<bool>( &dont_compress_fps_ )->zero_tokens() ,
      "Output all bits in the fingerprint, not just the ones that are set at least compressed-counts-cutoff times. The default is only to put out bits that are set at least once." )
    ( "bit-labels-to-ascii-file" , po::value<bool>( &bit_labels_to_ascii_file_ )->zero_tokens() ,
      "Puts the short codes for the bits on the first line at the top of the ASCII file, to easy the importing of the data into things like JMP and CART." )
    ( "compressed-counts-cutoff" , po::value<int>( &compressed_counts_cutoff_ )->default_value( 1 ) ,
      "Minimum number of compounds that set a bit for it to appear in the compressed fingerprint." )
    ( "ascii-fps-file" , po::value<string>( &ascii_fp_file_ ) ,
      "Filename for fingerprints in ASCII file." )
    ( "compact-fps-file" , po::value<string>( &compact_fp_file_ ) ,
      "Filename for compact fingerprints (just the numbers of the set bits)." )
    ( "names-fps-file" , po::value<string>( &names_fp_file_ ) ,
      "Filename for names fingerprint file (names of set bits, comma separated." )
    ( "log-file" , po::value<string>( &log_file_ ) ,
      "Name for capture of information about the run." )
    ( "bit-names-file" , po::value<string>( &bit_names_file_ ) ,
      "Name for file containing full names of all bits output, in correct order." );

}

// *******************************************************************
void LoobSettings::polish_args() {

  ascii_fps_ = !ascii_fp_file_.empty();
  compact_fps_ = !compact_fp_file_.empty();
  labels_not_bits_ = !names_fp_file_.empty();

  if( !ascii_fps_ && !compact_fps_ && !labels_not_bits_ ) {
    cerr << "You haven't specified an output file, so not doing anything." << endl;
    exit( 1 );
  }

  // the boost::program_options functions are a bit literal with strings
  if( ( bit_separator_[0] == '\\' && bit_separator_[1] == 't' ) ||
      bit_separator_ == "TAB" ) {
    bit_separator_ = '\t';
  } else if( ( bit_separator_[0] == '\\' && bit_separator_[1] == 'n' ) ||
	     bit_separator_ == "NEWLINE" ) {
    bit_separator_ = '\n';
  } else if( bit_separator_ == "SPACE" ) {
    bit_separator_ = ' ';
  } else if( bit_separator_ == "NOTHING" ) {
    bit_separator_ = "";
  }

  if( !pairs_ && !triplets_ && !quadruplets_ ) {
    triplets_ = true;
  }

  if( dist_bounds_.empty() ) {
    dist_bounds_.push_back( 4.5F );
    dist_bounds_.push_back( 7.0F );
    dist_bounds_.push_back( 10.0F );
    dist_bounds_.push_back( 14.0F );
    dist_bounds_.push_back( 19.0F );
    dist_bounds_.push_back( 24.0F );
  }

  transform( dist_bounds_.begin() , dist_bounds_.end() ,
	     inserter( sq_dist_bounds_ , sq_dist_bounds_.begin() ) ,
	     DACLIB::square<float> );

}
