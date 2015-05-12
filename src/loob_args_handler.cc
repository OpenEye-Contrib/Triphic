//
// file loob_args_handler.cc
// Dave Cosgrove
// AstraZeneca
// 2nd October 2007
//
// Uses boost program_options to process the command-line arguments and
// pull out the required stuff.

#include <iostream>
#include <limits>
#include <sstream>
#include <string>

#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "LoobSettings.H"

using namespace std;
namespace po = boost::program_options;

// ****************************************************************************
void build_program_options( boost::program_options::options_description &desc ,
			    LoobSettings &ls ) {

  desc.add_options()
    ( "help" , "Produce this help text" )
    ( "db-file,D" , po::value<string>( &ls.db_file_ ) ,
      "Database file name." )
    ( "smarts-file,S" , po::value<string>( &ls.smarts_file_ ) ,
      "SMARTS file name" )
    ( "points-file,P" , po::value<string>( &ls.points_file_ ) ,
      "Points file name" )
    ( "subset-file" , po::value<string>( &ls.subset_file_ ) ,
      "File containing subset of db-file that is to be processed." )
    ( "pairs" , po::value<bool>( &ls.pairs_ )->zero_tokens() ,
      "Output bits for pharmacophore pairs." )
    ( "triplets" , po::value<bool>( &ls.triplets_ )->zero_tokens() ,
      "Output bits for pharmacophore triplets (triplets only is default)." )
    ( "quadruplets" , po::value<bool>( &ls.quadruplets_ )->zero_tokens() ,
      "Output bits for pharmacophore quadruplets." )
    ( "single-conf-mols" , po::value<bool>( &ls.single_conf_mols_ )->zero_tokens() ,
      "Process each conformation as a separate molecule." )
    ( "chiral-fps" , po::value<bool>( &ls.chiral_fps_ )->zero_tokens() ,
      "For quadruplet features, whether to take account of chirality." )
    ( "bit-separator" , po::value<string>( &ls.bit_separator_ )->default_value( "" ) ,
      "Separator between bits in ASCII files." )
    ( "distance-bound,B" , po::value<vector<float> >( &ls.dist_bounds_ ) ,
      "Distance bin boundaries - default 4.5, 7.0, 10.0, 14.0, 19.0, 24.0." )
    ( "dont-compress-fps" , po::value<bool>( &ls.dont_compress_fps_ )->zero_tokens() ,
      "Output all bits in the fingerprint, not just the ones that are set at least compressed-counts-cutoff times. The default is only to put out bits that are set at least once." )
    ( "bit-labels-to-ascii-file" , po::value<bool>( &ls.bit_labels_to_ascii_file_ )->zero_tokens() ,
      "Puts the short codes for the bits on the first line at the top of the ASCII file, to easy the importing of the data into things like JMP and CART." )
    ( "compressed-counts-cutoff" , po::value<int>( &ls.compressed_counts_cutoff_ )->default_value( 1 ) ,
      "Minimum number of compounds that set a bit for it to appear in the compressed fingerprint." )
    ( "flush-fps-file" , po::value<string>( &ls.flush_fp_file_ ) ,
      "Filename for fingerprints in flush format." )
    ( "ascii-fps-file" , po::value<string>( &ls.ascii_fp_file_ ) ,
      "Filename for fingerprints in ASCII file." )
    // not working at the moment
    ( "compact-fps-file" , po::value<string>( &ls.compact_fp_file_ ) ,
      "Filename for compact fingerprints (just the numbers of the set bits)." )
    ( "names-fps-file" , po::value<string>( &ls.names_fp_file_ ) ,
      "Filename for names fingerprint file (names of set bits, comma separated." )
    ( "log-file" , po::value<string>( &ls.log_file_ ) ,
      "Name for capture of information about the run." )
    ( "bit-names-file" , po::value<string>( &ls.bit_names_file_ ) ,
      "Name for file containing full names of all bits output, in correct order." );

}

// ********************************************************************
void loob_args_handler( int argc , char **argv , LoobSettings &ls ) {

  po::options_description desc( "Allowed Options" );
  build_program_options( desc , ls );

  ostringstream oss;
  oss << desc;
  ls.usage_text_ = oss.str();

  po::variables_map vm;
  po::store( po::parse_command_line( argc , argv , desc ) , vm );
  po::notify( vm );

  if( 1 == argc || vm.count( "help" ) ) {
    cout << desc << endl;
    exit( 1 );
  }

}
