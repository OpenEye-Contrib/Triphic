//
// file TriphicSettings.cc
// David Cosgrove
// 2nd August 2007
//
// some functions for parsing string input

#include <fstream>
#include <iostream>
#include <iterator>

#include <mpi.h>

#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "BasePPhoreSite.H"
#include "TriphicSettings.H"

#include "oechem.h"
#include "oeomega2.h"

using namespace std;
namespace po = boost::program_options;

namespace DACLIB {
// In mpi_string_subs.cc
void mpi_send_string( const string &str , int dest_rank );
void mpi_send_strings_vector( const vector<string> &strs , int dest_rank );
void mpi_rec_string( int source_rank , string &str );
void mpi_rec_strings_vector( int source_rank , vector<string> &strs );
}

namespace DACLIB {
// in eponymous file
void split_filename( const string &filename , string &file_root ,
                     string &file_ext );
}

// **************************************************************************
TriphicSettings::TriphicSettings() :
  desc_( "Allowed Options" ) ,
  primary_score_( "RMS_AND_SIZE" ) ,
  max_hits_( 1000 ) ,
  opt_lig_rigid_( false ) ,
  opt_lig_flexi_( false ) ,
  gauss_shape_tani_( true ) ,
  grid_shape_tani_( false ) ,
  do_omega_( false ) ,
  do_flipper_( false ) ,
  do_ions_and_tauts_( false ) ,
  no_warts_( false ) ,
  all_best_hits_( false ) ,
  dont_do_sub_cliques_( false ) ,
  restart_( false ) ,
  chunk_size_( 5 ) ,
  score_method_( GtplDefs::RMS_AND_SIZE ) ,
  ring_norm_usage_( GtplDefs::IGNORE ) ,
  h_vector_usage_( GtplDefs::IGNORE ) ,
  lp_usage_( GtplDefs::IGNORE ) ,
  dist_tol_( 1.0F ) ,
  scaled_dist_tol_( -1.0F ) ,
  ring_norm_tol_( 30.0F ) ,
  h_vector_tol_( 30.0F ) ,
  lp_tol_( 30.0F ) ,
  max_rms_( 2.5F ) ,
  min_clique_size_( 3 ) ,
  min_grid_shape_tani_( -1.0F ) ,
  min_gauss_shape_tani_( -1.0F ) ,
  max_mmff_nrg_( std::numeric_limits<float>::max() ) ,
  min_surf_vol_( std::numeric_limits<float>::max() ) ,
  max_prot_clash_( std::numeric_limits<float>::max() ) ,
  min_inc_vol_( -1.0F ) ,
  vol_grid_spacing_( 0.5F ) ,
  restart_number_( 0 ) ,
  do_grid_shape_tani_cutoff_( false ) ,
  do_gauss_shape_tani_cutoff_( false ) ,
  do_mmff_nrg_cutoff_( false ) ,
  do_surf_vol_cutoff_( false ) ,
  do_prot_clash_cutoff_( false ) ,
  do_inc_vol_cutoff_( false ) ,
  no_hit_conf_number_( false ) ,
  output_query_to_hits_( false ) ,
  print_sites_and_stop_( false ) ,
  comma_output_( false ) ,
  output_scores_only_( false ) ,
  single_conf_mols_( false ) {

  build_program_options();
  ostringstream oss;
  oss << desc_;
  usage_text_ = oss.str();

}

// **************************************************************************
void TriphicSettings::parse_args( int argc , char **argv ) {
  
  po::variables_map vm;
  po::store( po::parse_command_line( argc , argv , desc_ ) , vm );
  po::notify( vm );

  if( argc < 2 ) {
    print_usage( cerr );
    exit( 1 );
  }
  if( vm.count( "help" ) ) {
    print_usage( cout );
    exit( 0 );
  }

  read_database_files_file();

}

// ***************************************************************************
bool TriphicSettings::operator!() const {

  if( db_files_.empty() && db_files_file_.empty() ) {
    cerr << "You need to specify at least one database file to search." << endl;
    return true;
  }

  if( query_file_.empty() ) {
    cerr << "You need to specify a query file to search with." << endl;
    return true;
  }

  if( points_file_.empty() ) {
    cerr << "No points file specified. Will use defaults." << endl;
  }

  if( smarts_file_.empty() ) {
    cerr << "No SMARTS file specified. Will use defaults." << endl;
  }

  if( output_file_.empty() ) {
    cerr << "You need to specify an output file." << endl;
    return true;
  }

  if( string( "ROBINS_PARETO" ) == primary_score_ &&  sites_score_file_.empty() &&
      ( ( query_file_.length() > 3 && string( ".ph4" ) == query_file_.substr( query_file_.length() - 4 ) ) ||
        ( query_file_.length() > 4 && string( ".sites" ) == query_file_.substr( query_file_.length() - 5 ) ) ) ) {
    cerr << "You can't do ROBINS_PARETO scoring without either a molecule query file or a sites score file."
         << endl;
    return true;
  }

  if( !set_score_method_from_primary_score() ) {
    cerr << primary_score_ << " is not a valid primary score value." << endl;
    return true;
  }
  // force calculation of grids for scoring
  if( score_method_ == GtplDefs::GRID_SHAPE_TANI ) {
    grid_shape_tani_ = true;
  }

  if( !set_vector_usage_from_string( ring_norm_string_ , ring_norm_usage_ ) ) {
    cerr << "Invalid option " << ring_norm_string_ << " given for --ring-normal-vector" << endl;
    print_usage( cerr );
    return true;
  }

  if( !set_vector_usage_from_string( h_vector_string_ ,	h_vector_usage_ ) ) {
    cerr << "Invalid option " << h_vector_string_ << " given for --h-bond-donor-vector" << endl;
    print_usage( cerr );
    return true;
  }

  if( !set_vector_usage_from_string( lp_string_ , lp_usage_ ) ) {
    cerr << "Invalid option " << lp_string_ << " given for --h-bond-acceptor-vector" << endl;
    print_usage( cerr );
    return true;
  }

  if( ( opt_lig_rigid_ || opt_lig_flexi_ ) && protein_file_.empty() ) {
    cerr << "You must specify a protein file in order to use --optimise-ligand*"
         << endl;
    return true;
  }

  if( opt_lig_rigid_ && opt_lig_flexi_ ) {
    cerr << "You can't have both of --optimise-ligand-rigid and" << endl
         << "--optimise-ligand-flexi." << endl;
    return true;
  }

  if( do_omega_ && !OEConfGen::OEOmegaIsLicensed( "toolkit" ) ) {
    cerr << "You have asked to use Omega, but you don't seem to have a toolkit license for it."
         << endl;
    return true;
  }

  if( do_flipper_ && !do_omega_ ) {
    cerr << "You have asked to use the flipper algorithm, but you're not using Omega." << endl
         << "That's not a fatal error, but you might want to check you getting what you wanted." << endl;
  }

  if( !sites_score_file_.empty() ) {
    string query_root , query_ext;
    DACLIB::split_filename( query_file_ , query_root , query_ext );
    if( query_ext != "sites" ) {
      cerr << "You can only specify --sites-score-sfile for a SITES query."
           << "I'm assuming the query file for a SITES query will have the"
           << "extension .sites." << endl;
      return true;
    }
  }

  set_cutoff_flags();

  return false;

}

// ***************************************************************************
void TriphicSettings::print_usage( ostream &os ) const {

  os << usage_text_ << endl;

}

// **************************************************************************
bool TriphicSettings::set_score_method_from_primary_score() const {

  static vector<pair<string,GtplDefs::SCORE_METHOD> > valid_values;
  if( valid_values.empty() ) {
    valid_values.push_back( make_pair( string( "RMS_AND_SIZE" ) ,
                                       GtplDefs::RMS_AND_SIZE ) );
    valid_values.push_back( make_pair( string( "ROBINS_PARETO" ) ,
                                       GtplDefs::ROBINS_SCORE_PARETO ) );
    valid_values.push_back( make_pair( string( "OVERALL_PARETO" ) ,
                                       GtplDefs::OVERALL_SCORE_PARETO ) );
    valid_values.push_back( make_pair( string( "GRID_SHAPE_TANI" ) ,
                                       GtplDefs::GRID_SHAPE_TANI ) );
    valid_values.push_back( make_pair( string( "GAUSS_SHAPE_TANI" ) ,
                                       GtplDefs::GAUSS_SHAPE_TANI ) );
    valid_values.push_back( make_pair( string( "SURFACE_OVLP_VOLUME" ) ,
                                       GtplDefs::SURFACE_OVLP_VOLUME ) );
    valid_values.push_back( make_pair( string( "PROTEIN_CLASH" ) ,
                                       GtplDefs::PROTEIN_CLASH ) );
    valid_values.push_back( make_pair( string( "MMFF_ENERGY" ) ,
                                       GtplDefs::MMFF_NRG ) );
  }

  for( int i = 0 , is = valid_values.size() ; i < is ; ++i ) {
    if( valid_values[i].first == primary_score_ ) {
      score_method_ = valid_values[i].second;
      return true;
    }
  }

  return false;

}

// **************************************************************************
bool TriphicSettings::set_vector_usage_from_string( const string &str ,
                                                    GtplDefs::DIRS_USAGE &du ) const {

  static vector<pair<string,GtplDefs::DIRS_USAGE> > valid_values;
  if( valid_values.empty() ) {
    valid_values.push_back( make_pair( string( "IGNORE" ) , GtplDefs::IGNORE ) );
    valid_values.push_back( make_pair( string( "SCREEN" ) , GtplDefs::SCREEN ) );
    valid_values.push_back( make_pair( string( "ALIGN" ) , GtplDefs::ALIGN ) );
  }

  for( int i = 0 , is = valid_values.size() ; i < is ; ++i ) {
    if( valid_values[i].first == str ) {
      du = valid_values[i].second;
      return true;
    }
  }

  return false;

}

// **************************************************************************
void TriphicSettings::set_cutoff_flags() const {

  if( min_grid_shape_tani_ >= 0.0F && min_grid_shape_tani_ <= 1.0F )
    do_grid_shape_tani_cutoff_ = true;
  else
    do_grid_shape_tani_cutoff_ = false;

  if( min_gauss_shape_tani_ >= 0.0F && min_gauss_shape_tani_ <= 1.0F )
    do_gauss_shape_tani_cutoff_ = true;
  else
    do_gauss_shape_tani_cutoff_ = false;

  if( max_mmff_nrg_ < numeric_limits<float>::max() )
    do_mmff_nrg_cutoff_ = true;
  else
    do_mmff_nrg_cutoff_ = false;

  if( min_surf_vol_ != 0.0F )
    do_surf_vol_cutoff_ = true;
  else
    do_surf_vol_cutoff_ = false;

  if( max_prot_clash_ < numeric_limits<float>::max() )
    do_prot_clash_cutoff_ = true;
  else
    do_prot_clash_cutoff_ = false;

  if( min_inc_vol_ > -1.0F )
    do_inc_vol_cutoff_ = true;
  else
    do_inc_vol_cutoff_ = false;

}

// **************************************************************************
void TriphicSettings::print_required_sites_and_points( ostream &os ) {

  if( !req_sites_or_.empty() ) {
    os << "Required sites OR'd list : ";
    copy( req_sites_or_.begin() , req_sites_or_.end() ,
          ostream_iterator<string>( os , " " ) );
    os << endl;
  }

  if( !req_sites_and_.empty() ) {
    os << "Required sites AND'd list : ";
    copy( req_sites_and_.begin() , req_sites_and_.end() ,
          ostream_iterator<string>( os , " " ) );
    os << endl;
  }

  if( !req_points_or_.empty() ) {
    os << "Required points OR'd list : ";
    copy( req_points_or_.begin() , req_points_or_.end() ,
          ostream_iterator<string>( os , " " ) );
    os << endl;
  }

  if( !req_points_and_.empty() ) {
    os << "Required points AND'd list : ";
    copy( req_points_and_.begin() , req_points_and_.end() ,
          ostream_iterator<string>( os , " " ) );
    os << endl;
  }

}

// **************************************************************************
// make sure the sites in req_sites_or_ and req_sites_and_ appear in query_sites
bool TriphicSettings::check_required_sites( vector<BasePPhoreSite *> &query_sites ) {

  bool ret_val1 = check_required_sites( query_sites , req_sites_or_ );
  bool ret_val2 = check_required_sites( query_sites , req_sites_and_ );

  return( ret_val1 && ret_val2 );

}

// **************************************************************************
bool TriphicSettings::check_required_sites( vector<BasePPhoreSite *> &query_sites ,
                                            const vector<string> &req_sites ) {

  bool ret_val = true;
  for( int k = 0 , ks = req_sites.size() ; k < ks ; ++k ) {
    bool found_it = false;
    for( int j = 0 , js = query_sites.size(); j < js ; ++j ) {
      if( req_sites[k] == query_sites[j]->get_full_name() ) {
        found_it = true;
        break;
      }
    }
    if( !found_it ) {
      ret_val = false;
      cerr << "Error : required site " << req_sites[k]
              << " not found in query sites." << endl;
      if( string::npos == req_sites[k].find( ':' ) ) {
        cerr << "The site name must be given in full, e.g. "
             << query_sites.front()->get_full_name() << endl;
      }
    }
  }

  return ret_val;

}

// **************************************************************************
// make sure the sites in req_points_or_ and req_points_and_ appear in query_sites
bool TriphicSettings::check_required_points( vector<BasePPhoreSite *> &query_sites ) {

  bool ret_val1 = check_required_points( query_sites , req_points_or_ );
  bool ret_val2 = check_required_points( query_sites , req_points_and_ );

  return( ret_val1 && ret_val2 );

}

// **************************************************************************
bool TriphicSettings::check_required_points( vector<BasePPhoreSite *> &query_sites ,
                                             const vector<string> &req_points ) {

  bool ret_val = true;
  for( int k = 0 , ks = req_points.size() ; k < ks ; ++k ) {
    bool found_it = false;
    for( int j = 0 , js = query_sites.size(); j < js ; ++j ) {
      if( req_points[k] == query_sites[j]->get_type_string() ) {
        found_it = true;
        break;
      }
    }
    if( !found_it ) {
      ret_val = false;
      cerr << "Error : required point " << req_points[k]
              << " not found in query sites." << endl;
    }
  }

  return ret_val;

}

// ****************************************************************************
void TriphicSettings::send_contents_via_mpi( int dest_slave ) {

  using namespace DACLIB;
  mpi_send_string( query_file_ , dest_slave );
  mpi_send_strings_vector( db_files_ , dest_slave );
  mpi_send_string( output_file_ , dest_slave );
  mpi_send_string( smarts_file_ , dest_slave );
  mpi_send_string( points_file_ , dest_slave );
  mpi_send_string( protein_file_ , dest_slave );
  mpi_send_string( subset_file_ , dest_slave );
  mpi_send_string( primary_score_ , dest_slave );
  mpi_send_string( ring_norm_string_ , dest_slave );
  mpi_send_string( h_vector_string_ , dest_slave );
  mpi_send_string( lp_string_ , dest_slave );

  MPI_Send( &max_hits_ , 1 , MPI_UNSIGNED , dest_slave , 0 , MPI_COMM_WORLD );

  int i = static_cast<int>( opt_lig_rigid_ );
  MPI_Send( &i , 1 , MPI_INT , dest_slave , 0 , MPI_COMM_WORLD );
  i = static_cast<int>( opt_lig_flexi_ );
  MPI_Send( &i , 1 , MPI_INT , dest_slave , 0 , MPI_COMM_WORLD );

  i = static_cast<int>( gauss_shape_tani_ );
  MPI_Send( &i , 1 , MPI_INT , dest_slave , 0 , MPI_COMM_WORLD );
  i = static_cast<int>( grid_shape_tani_ );
  MPI_Send( &i , 1 , MPI_INT , dest_slave , 0 , MPI_COMM_WORLD );
  i = static_cast<int>( do_omega_ );
  MPI_Send( &i , 1 , MPI_INT , dest_slave , 0 , MPI_COMM_WORLD );
  i = static_cast<int>( do_flipper_ );
  MPI_Send( &i , 1 , MPI_INT , dest_slave , 0 , MPI_COMM_WORLD );
  i = static_cast<int>( do_ions_and_tauts_ );
  MPI_Send( &i , 1 , MPI_INT , dest_slave , 0 , MPI_COMM_WORLD );
  i = static_cast<int>( no_warts_ );
  MPI_Send( &i , 1 , MPI_INT , dest_slave , 0 , MPI_COMM_WORLD );
  i = static_cast<int>( all_best_hits_ );
  MPI_Send( &i , 1 , MPI_INT , dest_slave , 0 , MPI_COMM_WORLD );
  i = static_cast<int>( dont_do_sub_cliques_ );
  MPI_Send( &i , 1 , MPI_INT , dest_slave , 0 , MPI_COMM_WORLD );
  MPI_Send( &chunk_size_ , 1 , MPI_INT , dest_slave , 0 , MPI_COMM_WORLD );

  i = static_cast<int>( score_method_ );
  MPI_Send( &i , 1 , MPI_INT , dest_slave , 0 , MPI_COMM_WORLD );

  i = static_cast<int>( ring_norm_usage_ );
  MPI_Send( &i , 1 , MPI_INT , dest_slave , 0 , MPI_COMM_WORLD );
  i = static_cast<int>( h_vector_usage_ );
  MPI_Send( &i , 1 , MPI_INT , dest_slave , 0 , MPI_COMM_WORLD );
  i = static_cast<int>( lp_usage_ );
  MPI_Send( &i , 1 , MPI_INT , dest_slave , 0 , MPI_COMM_WORLD );

  MPI_Send( &ring_norm_tol_ , 1 , MPI_FLOAT , dest_slave , 0 , MPI_COMM_WORLD );
  MPI_Send( &h_vector_tol_ , 1 , MPI_FLOAT , dest_slave , 0 , MPI_COMM_WORLD );
  MPI_Send( &lp_tol_ , 1 , MPI_FLOAT , dest_slave , 0 , MPI_COMM_WORLD );
  MPI_Send( &dist_tol_ , 1 , MPI_FLOAT , dest_slave , 0 , MPI_COMM_WORLD );
  MPI_Send( &scaled_dist_tol_ , 1 , MPI_FLOAT , dest_slave , 0 , MPI_COMM_WORLD );
  MPI_Send( &max_rms_ , 1 , MPI_FLOAT , dest_slave , 0 , MPI_COMM_WORLD );

  MPI_Send( &min_clique_size_ , 1 , MPI_INT , dest_slave , 0 , MPI_COMM_WORLD );

  MPI_Send( &min_grid_shape_tani_ , 1 , MPI_FLOAT , dest_slave , 0 , MPI_COMM_WORLD );
  MPI_Send( &min_gauss_shape_tani_ , 1 , MPI_FLOAT , dest_slave , 0 , MPI_COMM_WORLD );
  MPI_Send( &max_mmff_nrg_ , 1 , MPI_FLOAT , dest_slave , 0 , MPI_COMM_WORLD );
  MPI_Send( &min_surf_vol_ , 1 , MPI_FLOAT , dest_slave , 0 , MPI_COMM_WORLD );
  MPI_Send( &max_prot_clash_ , 1 , MPI_FLOAT , dest_slave , 0 , MPI_COMM_WORLD );
  MPI_Send( &min_inc_vol_ , 1 , MPI_FLOAT , dest_slave , 0 , MPI_COMM_WORLD );
  MPI_Send( &vol_grid_spacing_ , 1 , MPI_FLOAT , dest_slave , 0 , MPI_COMM_WORLD );
  MPI_Send( &restart_number_ , 1 , MPI_UNSIGNED , dest_slave , 0 , MPI_COMM_WORLD );

  i = static_cast<int>( do_grid_shape_tani_cutoff_ );
  MPI_Send( &i , 1 , MPI_INT , dest_slave , 0 , MPI_COMM_WORLD );
  i = static_cast<int>( do_gauss_shape_tani_cutoff_ );
  MPI_Send( &i , 1 , MPI_INT , dest_slave , 0 , MPI_COMM_WORLD );
  i = static_cast<int>( do_mmff_nrg_cutoff_ );
  MPI_Send( &i , 1 , MPI_INT , dest_slave , 0 , MPI_COMM_WORLD );
  i = static_cast<int>( do_surf_vol_cutoff_ );
  MPI_Send( &i , 1 , MPI_INT , dest_slave , 0 , MPI_COMM_WORLD );
  i = static_cast<int>( do_prot_clash_cutoff_ );
  MPI_Send( &i , 1 , MPI_INT , dest_slave , 0 , MPI_COMM_WORLD );
  i = static_cast<int>( do_inc_vol_cutoff_ );
  MPI_Send( &i , 1 , MPI_INT , dest_slave , 0 , MPI_COMM_WORLD );

  mpi_send_strings_vector( req_sites_or_ , dest_slave );
  mpi_send_strings_vector( req_sites_and_ , dest_slave );
  mpi_send_strings_vector( req_points_or_ , dest_slave );
  mpi_send_strings_vector( req_points_and_ , dest_slave );
  mpi_send_strings_vector( ignore_sites_ , dest_slave );

}

// ****************************************************************************
void TriphicSettings::receive_contents_via_mpi() {

  using namespace DACLIB;
  mpi_rec_string( 0 , query_file_ );
  mpi_rec_strings_vector( 0 , db_files_ );
  mpi_rec_string( 0 , output_file_ );
  mpi_rec_string( 0 , smarts_file_ );
  mpi_rec_string( 0 , points_file_ );
  mpi_rec_string( 0 , protein_file_ );
  mpi_rec_string( 0 , subset_file_ );
  mpi_rec_string( 0 , primary_score_ );
  mpi_rec_string( 0 , ring_norm_string_ );
  mpi_rec_string( 0 , h_vector_string_ );
  mpi_rec_string( 0 , lp_string_ );

  unsigned int ui = 0;
  MPI_Recv( &ui , 1 , MPI_UNSIGNED , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  max_hits_ = ui;

  int i = 0;
  MPI_Recv( &i , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  opt_lig_rigid_ = static_cast<bool>( i );
  MPI_Recv( &i , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  opt_lig_flexi_ = static_cast<bool>( i );

  MPI_Recv( &i , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  gauss_shape_tani_ = static_cast<bool>( i );
  MPI_Recv( &i , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  grid_shape_tani_ = static_cast<bool>( i );
  MPI_Recv( &i , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  do_omega_ = static_cast<bool>( i );
  MPI_Recv( &i , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  do_flipper_ = static_cast<bool>( i );
  MPI_Recv( &i , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  do_ions_and_tauts_ = static_cast<bool>( i );
  MPI_Recv( &i , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  no_warts_ = static_cast<bool>( i );
  MPI_Recv( &i , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  all_best_hits_ = static_cast<bool>( i );
  MPI_Recv( &i , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  dont_do_sub_cliques_ = static_cast<bool>( i );
  MPI_Recv( &chunk_size_ , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );

  MPI_Recv( &i , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  score_method_ = static_cast<GtplDefs::SCORE_METHOD>( i );

  MPI_Recv( &i , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  ring_norm_usage_ = static_cast<GtplDefs::DIRS_USAGE>( i );

  MPI_Recv( &i , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  h_vector_usage_ = static_cast<GtplDefs::DIRS_USAGE>( i );

  MPI_Recv( &i , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  lp_usage_ = static_cast<GtplDefs::DIRS_USAGE>( i );

  MPI_Recv( &ring_norm_tol_ , 1 , MPI_FLOAT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  MPI_Recv( &h_vector_tol_ , 1 , MPI_FLOAT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  MPI_Recv( &lp_tol_ , 1 , MPI_FLOAT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  MPI_Recv( &dist_tol_ , 1 , MPI_FLOAT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  MPI_Recv( &scaled_dist_tol_ , 1 , MPI_FLOAT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  MPI_Recv( &max_rms_ , 1 , MPI_FLOAT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );

  MPI_Recv( &min_clique_size_ , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );

  MPI_Recv( &min_grid_shape_tani_ , 1 , MPI_FLOAT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  MPI_Recv( &min_gauss_shape_tani_ , 1 , MPI_FLOAT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  MPI_Recv( &max_mmff_nrg_ , 1 , MPI_FLOAT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  MPI_Recv( &min_surf_vol_ , 1 , MPI_FLOAT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  MPI_Recv( &max_prot_clash_ , 1 , MPI_FLOAT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  MPI_Recv( &min_inc_vol_ , 1 , MPI_FLOAT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  MPI_Recv( &vol_grid_spacing_ , 1 , MPI_FLOAT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  MPI_Recv( &restart_number_ , 1 , MPI_UNSIGNED , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );

  MPI_Recv( &i , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  do_grid_shape_tani_cutoff_ = static_cast<bool>( i );
  MPI_Recv( &i , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  do_gauss_shape_tani_cutoff_ = static_cast<bool>( i );
  MPI_Recv( &i , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  do_mmff_nrg_cutoff_ = static_cast<bool>( i );
  MPI_Recv( &i , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  do_surf_vol_cutoff_ = static_cast<bool>( i );
  MPI_Recv( &i , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  do_prot_clash_cutoff_ = static_cast<bool>( i );
  MPI_Recv( &i , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  do_inc_vol_cutoff_ = static_cast<bool>( i );

  mpi_rec_strings_vector( 0 , req_sites_or_ );
  mpi_rec_strings_vector( 0 , req_sites_and_ );
  mpi_rec_strings_vector( 0 , req_points_or_ );
  mpi_rec_strings_vector( 0 , req_points_and_ );
  mpi_rec_strings_vector( 0 , ignore_sites_ );

}

// **************************************************************************
void TriphicSettings::add_req_sites_and( string &new_site ) {

  req_sites_and_.push_back( new_site );

}

// **************************************************************************
void TriphicSettings::build_program_options() {

  desc_.add_options()
      ( "help" , "Produce this help text" )
      ( "db-file,D" , po::value<vector<string> >( &db_files_ ) ,
        "Database file name. Multiple instances allowed" )
      ( "db-files-file" , po::value<string>( &db_files_file_ ) ,
        "File containing list of database files to search." )
      ( "query-file,Q" , po::value<string>( &query_file_ ) ,
        "Query file name." )
      ( "smarts-file,S" , po::value<string>( &smarts_file_ ) ,
        "SMARTS file name" )
      ( "points-file,P" , po::value<string>( &points_file_ ) ,
        "Points file name" )
      ( "output-file,O" , po::value<string>( &output_file_ ) ,
        "Output file name" )
      ( "max-hits,M" , po::value<unsigned int>( &max_hits_ )->default_value( 1000 ) ,
        "Maximum hits" )
      ( "protein-file" , po::value<string>( &protein_file_ ) ,
        "File containing protein active site." )
      ( "subset-file" , po::value<string>( &subset_file_ ) ,
        "File containing subset of db-file that is to be searched." )
      ( "grid-score-file,G" , po::value<vector<string> >( &grid_vol_files_ ) ,
        "File containing grid for calculating overlap volume with hit. Multiple files allowed, but only 1 per option." )
      ( "primary-score" ,
        po::value<string>( &primary_score_ )->default_value( "RMS_AND_SIZE" ) ,
        "Primary Score method - RMS_AND_SIZE|ROBINS_PARETO|OVERALL_PARETO|GAUSS_SHAPE_TANI|GRID_SHAPE_TANI|SURFACE_OVLP_VOLUME|PROTEIN_CLASH|MMFF_ENERGY" )
      ( "ring-normal-vector" , po::value<string>( &ring_norm_string_ )->default_value( "IGNORE" ) ,
        "How to deal with ring normals - IGNORE|SCREEN|ALIGN." )
      ( "ring-normal-tolerance" , po::value<float>( &ring_norm_tol_ )->default_value( 30.0F ) ,
        "Tolerance on ring normal alignment." )
      ( "h-bond-donor-vector" , po::value<string>( &h_vector_string_ )->default_value( "IGNORE" ) ,
        "How to deal with H-bond donor vectors - IGNORE|SCREEN|ALIGN." )
      ( "h-bond-donor-tolerance" , po::value<float>( &h_vector_tol_ )->default_value( 30.0F ) ,
        "Tolerance on H-bond donor vector alignment." )
      ( "h-bond-acceptor-vector" , po::value<string>( &lp_string_ )->default_value( "IGNORE" ) ,
        "How to deal with H-bond acceptor vectors - IGNORE|SCREEN|ALIGN." )
      ( "h-bond-acceptor-tolerance" , po::value<float>( &lp_tol_ )->default_value( 30.0F ) ,
        "Tolerance on H-bond acceptor vector alignment." )
      ( "optimise-ligand-rigid" , po::value<bool>( &opt_lig_rigid_ )->zero_tokens() ,
        "Rigid-body optimisation of the ligand in the protein." )
      ( "optimise-ligand-flexi" , po::value<bool>( &opt_lig_flexi_ )->zero_tokens() ,
        "Optimisation of the ligand in the protein with torsional flexibility." )
      ( "use-omega" , po::value<bool>( &do_omega_ )->zero_tokens()->default_value( false ) ,
        "Treat each database molecule as a single conformer, and generate conformers using omega. Requires licence for Omega2 toolkit." )
      ( "use-flipper" , po::value<bool>( &do_flipper_ )->zero_tokens()->default_value( false ) ,
        "If using Omega, run flipper algorithm on unspecified stereo centres." )
      ( "do-ionisation-and-tautomers" , po::value<bool>( &do_ions_and_tauts_ )->zero_tokens()->default_value( false ) ,
        "Apply the officially blessed ionisation and tautomer model to the structures." )
      ( "no-warts" , po::value<bool>( &no_warts_ )->zero_tokens()->default_value( false ) ,
        "Don't add info to molecules names about tautomers, ionisation and flipper stereochemical states." )
      ( "min-clique-size" , po::value<int>( &min_clique_size_ )->default_value( 3 ) ,
        "Minimum clique size for hits" )
      ( "dist-tol" , po::value<float>( &dist_tol_ )->default_value( 1.0F ) , "Maximum distance for matching points in cliques." )
      ( "scaled-dist-tol" , po::value<float>( &scaled_dist_tol_ )->default_value( -1.0F ) , "Scaling for tolerance of matching distance in cliques." )
      ( "max-rms" , po::value<float>( &max_rms_ )->default_value( 2.5F ) ,
        "Maximum tolerance on RMS for overlaid sites." )
      ( "min-grid-shape-tanimoto" , po::value<float>( &min_grid_shape_tani_ )->default_value( -1.0F ) ,
        "Minimum cutoff for grid shape tanimoto." )
      ( "min-gauss-shape-tanimoto" , po::value<float>( &min_gauss_shape_tani_ )->default_value( -1.0F ) ,
        "Minimum cutoff for gaussian shape tanimoto." )
      ( "max-mmff-energy" , po::value<float>( &max_mmff_nrg_ )->default_value( numeric_limits<float>::max() ) ,
        "Cutoff for MMFF energy between hit and protein." )
      ( "min-surface-volume" , po::value<float>( &min_surf_vol_ )->default_value( 0.0F ) ,
        "Cutoff for surface overlap volume." )
      ( "max-protein-clash" , po::value<float>( &max_prot_clash_ )->default_value( numeric_limits<float>::max() ) ,
        "Cutoff for volume overlap between hit and protein." )
      ( "min-included-volume" , po::value<float>( &min_inc_vol_ )->default_value( -1.0F ) ,
        "Cutoff for overlap volume between query and hit." )
      ( "volume-grid-spacing" , po::value<float>( &vol_grid_spacing_ )->default_value( 0.5F ) ,
        "Spacing for grid points in volume calculations." )
      ( "no-hit-conf-number" , po::value<bool>( &no_hit_conf_number_ )->zero_tokens() ,
        "Whether to add the conformation number to the hit molecule name." )
      ( "output-query-to-hits" , po::value<bool>( &output_query_to_hits_ )->zero_tokens() ,
        "Whether to put output the query structure to the top of the hits output file, (not Grappel state file output format)." )
      ( "print-sites-and-stop" , po::value<bool>( &print_sites_and_stop_ )->zero_tokens() ,
        "As it says on the tin." )
      ( "comma-output" , po::value<bool>( &comma_output_ )->zero_tokens() ,
        "Puts comma-separated list in .scores file." )
      ( "output-scores-only" , po::value<bool>( &output_scores_only_ )->zero_tokens() ,
        "Output scores file only, not the overlaid hit molecules." )
      ( "single-conf-mols" , po::value<bool>( &single_conf_mols_ )->zero_tokens() ,
        "Specifies that database molecules should be treated as single conformations." )
      ( "req-sites-or" , po::value<vector<string> >( &req_sites_or_ ) ,
        "Required sites - OR'd list." )
      ( "req-sites-and" , po::value<vector<string> >( &req_sites_and_ ) ,
        "Required sites - AND'd list." )
      ( "req-points-or" , po::value<vector<string> >( &req_points_or_ ) ,
        "Required points - OR'd list." )
      ( "req-points-and" , po::value<vector<string> >( &req_points_and_ ) ,
        "Required points - AND'd list." )
      ( "ignore-sites" , po::value<vector<string> >( &ignore_sites_ ) ,
        "Sites in the query that are to be ignored during the search." )
      ( "dont-give-me" , po::value<string>( &not_smarts_file_ ) ,
        "File of SMARTS that the hits must not match." )
      ( "i-already-know-about" , po::value<string>( &not_smarts_file_ ) ,
        "File of SMARTS that the hits must not match." )
      ( "sites-score-file" , po::value<string>( &sites_score_file_ ) ,
        "File containing a structure to be used for the ligand scores, when searching with a sites query file." )
      ( "gaussian-shape-tanimoto" , po::value<bool>( &gauss_shape_tani_ )->zero_tokens() ,
        "Calculate the shape tanimoto by Gaussian overlap (default=true)" )
      ( "grid-shape-tanimoto" , po::value<bool>( &grid_shape_tani_ )->zero_tokens() ,
        "Calculate the shape tanimoto by grid method (default=false)" )
      ( "all-best-hits" , po::value<bool>( &all_best_hits_ )->zero_tokens() ,
        "Keep best hits per clique rather than best hits per target molecule." )
      ( "dont-do-sub-cliques" , po::value<bool>( &dont_do_sub_cliques_ )->zero_tokens() ,
        "Don't generate and score the sub-cliques of cliques bigger than min_clique_size." )
      ( "restart" , po::value<bool>( &restart_ )->zero_tokens() ,
        "Restart from last point in restart file." )
      ( "restart-number" , po::value<unsigned int>( &restart_number_ ) ,
        "Molecule number from which to restart, over-riding what's in the restart file." )
      ( "parallel-chunk-size" , po::value<unsigned int>( &chunk_size_ ) ,
        "Number of molecules to do in one step on each slave when running in parallel (default 5)." );

}

// **************************************************************************
void TriphicSettings::read_database_files_file() {

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
