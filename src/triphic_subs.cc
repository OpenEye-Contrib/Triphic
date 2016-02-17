//
// file triphic_subs.cc
// David Cosgrove
// 31st July 2007
//
// Some of the triphic-specific functions for new triphic (v3.0).

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <oechem.h>
#include <oeomega2.h>
#include <oeszybki.h>

#include <mpi.h>

#include <boost/cast.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/construct.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/cast.hpp>
#include <boost/shared_ptr.hpp>

#include "stddefs.H"
#include "BasePPhoreSite.H"
#include "DefaultPointsDefs.H"
#include "DiamondOverlay.H"
#include "FileExceptions.H"
#include "OverlayScore.H"
#include "ParseSMARTSXML.H"
#include "PharmPoint.H"
#include "SinglePPhoreSite.H"
#include "SMARTSExceptions.H"
#include "TriphicSettings.H"
#include "VolumeGrid.H"

using namespace boost;
using namespace std;
using namespace OEChem;
using namespace OESz;
using namespace OEConfGen;

typedef shared_ptr<OEMol> pOEMol;

namespace DACLIB {
void read_smarts_file( const string &smarts_file ,
                       vector<pair<string,string> > &input_smarts ,
                       vector<pair<string,string> > &smarts_sub_defn );
void make_pphore_sites( OEMol &mol , PharmPoint &pharm_points ,
                        const vector<pair<string,string> > &input_smarts ,
                        vector<pair<string,string> > &smarts_sub_defn ,
                        vector<vector<SinglePPhoreSite *> > &pharm_sites );
VolumeGrid *prepare_mol_grid( OEMolBase *mol );
void set_environment_var_from_mpi();

// In mpi_string_subs.cc
void mpi_send_string( const string &str , int dest_rank );
void mpi_rec_string( int source_rank , string &str );

void split_filename( const string &filename , string &file_root ,
                     string &file_ext );
string create_cansmi( const OEMolBase &in_mol );
}

// find_cliques is in the eponymous file.
void find_cliques( const vector<BasePPhoreSite *> &sites1 ,
                   const vector<SinglePPhoreSite *> &sites2 ,
                   float dist_tol , float scaled_dist_tol ,
                   bool dont_do_sub_cliques ,
                   int min_clique_size , vector<vector<int> > &cliques );
// score_and_store_cliques is in the eponymous file
void score_and_store_cliques( const string &query_name , int query_conf_num ,
                              OEMolBase *query_conf ,
                              vector<BasePPhoreSite *> &query_sites ,
                              vector<SinglePPhoreSite *> &query_score_sites ,
                              vector<SinglePPhoreSite *> &target_sites ,
                              OEMol *target_mol , int target_conf_num ,
                              shared_ptr<DACLIB::VolumeGrid> &query_solid_grid ,
                              shared_ptr<DACLIB::VolumeGrid> &protein_grid ,
                              const vector<pair<string,DACLIB::VolumeGrid *> > &score_vol_grids ,
                              const vector<vector<int> > &cliques ,
                              GtplDefs::DIRS_USAGE ring_norm_usage ,
                              GtplDefs::DIRS_USAGE h_vector_usage ,
                              GtplDefs::DIRS_USAGE lp_usage ,
                              float ring_norm_tol , float h_vector_tol ,
                              float lp_tol , bool do_size_co , int size_cutoff ,
                              bool do_rms_co , float rms_cutoff ,
                              bool do_grid_shape_tani_co ,
                              float grid_shape_tani_cutoff ,
                              bool do_gauss_shape_tani_co ,
                              float gauss_shape_tani_cutoff ,
                              bool do_surf_ovlp_co , float surface_overlap_cutoff ,
                              bool do_inc_vol_co , float inc_vol_cutoff ,
                              bool do_prot_clash_co , float prot_clash_cutoff ,
                              bool do_mmff_co , float mmff_cutoff ,
                              const vector<string> &required_sites_or ,
                              const vector<string> &required_sites_and ,
                              const vector<string> &required_points_or ,
                              const vector<string> &required_points_and ,
                              shared_ptr<OESzybki> &szybki ,
                              list<OverlayScore *> &clique_overlays );
// this also in score_and_store_cliques.cc
bool add_overlay_score_to_list( OverlayScore *new_score ,
                                list<OverlayScore *> &clique_overlays );
// and so is this
void overlay_mols_and_sites( OEMol *target_mol ,
                             const vector<BasePPhoreSite *> &query_sites ,
                             const vector<SinglePPhoreSite *> &target_sites ,
                             OverlayScore *ov_score , bool use_ring_norms ,
                             bool use_h_vectors , bool use_lps ,
                             OEMolBase *&target_conf );
// and this
float optimise_overlay( shared_ptr<OESzybki> &szybki ,
                        OEMolBase *target_conf ,
                        vector<SinglePPhoreSite *> &ov_target_sites );

// in send_messages_to_tp_slaves.cc
void send_openeye_license_to_slaves( int num_slaves );
void send_cwd_to_slaves( int num_slaves );
void send_finished_messages( int num_slaves );
void receive_new_cwd();
void receive_database_details( int &db_start , int &db_step );
void send_progress_to_master( const string &progress_report );

// in step_oemolstream.cc
namespace DACLIB {
void open_databasefile( const string &db_file ,
                        bool single_conf_mols ,
                        oemolistream *&ims );
pOEMol read_nth_mol_from_oemolistream( unsigned int next_mol ,
                                       const vector<string> &db_files ,
                                       bool single_conf_mols );
}

void read_moe_ph4_file( TriphicSettings &ts ,
                        PharmPoint &pharm_points ,
                        vector<pair<string,string> > &input_smarts ,
                        vector<pair<string,string> > &smarts_sub_defn ,
                        vector<BasePPhoreSite *> &query_sites ,
                        vector<pair<string,DACLIB::VolumeGrid *> > &score_vol_grids );

// molecule standardisation and tautomer enumeration, functions in
// libtautenum.a
// void prepare_molecule( OEMolBase &mol );
OEMolBase *standardise_tautomer( OEMolBase &in_mol );
vector<OEMolBase *> enumerate_ions( OEMolBase &in_mol ,
                                    const string &prot_stand_smirks ,
                                    const string &prot_enum_smirks ,
                                    const string &prot_smirks_vbs );
vector<OEMolBase *> enumerate_tautomers( OEMolBase &in_mol ,
                                         const string &smirks_defs ,
                                         const string &smirks_vbs );

// ***************************************************************************
void read_smarts_file( const string &smarts_file ,
                       vector<pair<string,string> > &input_smarts ,
                       vector<pair<string,string> > &smarts_sub_defn ) {

  try {
    DACLIB:: read_smarts_file( smarts_file , input_smarts , smarts_sub_defn );
  } catch( DACLIB::SMARTSFileError &e ) {
    cerr << e.what() << endl;
    exit( 1 );
  } catch( DACLIB::SMARTSSubDefnError &e ) {
    cerr << e.what() << endl;
    exit( 1 );
  } catch( DACLIB::FileReadOpenError &e ) {
    cerr << e.what() << endl;
    exit( 1 );
  }

}

// ***************************************************************************
void read_query_sites( const string &query_file ,
                       vector<BasePPhoreSite *> &query_sites ) {

  // in eponymous file
  BasePPhoreSite *read_site_from_file( istream &is );

  ifstream ifs( query_file.c_str() );
  if( !ifs.good() ) {
    throw DACLIB::FileReadOpenError( query_file.c_str() );
  }

  while( 1 ) {
    query_sites.push_back( read_site_from_file( ifs ) );
    if( !query_sites.back() ) {
      query_sites.pop_back();
      break;
    }
  }

}

// ***************************************************************************
void create_query_sites( vector<pair<string,string> > &input_smarts ,
                         vector<pair<string,string> > &smarts_sub_defn ,
                         PharmPoint &pharm_points ,
                         OEMolBase *&query_mol ,
                         vector<BasePPhoreSite *> &query_sites ,
                         vector<SinglePPhoreSite *> &query_score_sites ) {

  // make_pphore_sites expects a multi-conf molecule
  vector<vector<SinglePPhoreSite *> > next_sites;
  OEMol next_mol( *query_mol );
  try {
    DACLIB::make_pphore_sites( next_mol , pharm_points , input_smarts ,
                               smarts_sub_defn , next_sites );
  } catch( DACLIB::SMARTSDefnError &e ) {
    cerr << e.what() << endl;
    exit( 1 );
  }
  // we know that there's only 1 conformation per molecule
  query_sites = vector<BasePPhoreSite*>( next_sites.front().begin() ,
                                         next_sites.front().end() );
  query_score_sites = next_sites.front();

}

// ***************************************************************************
// delete from query_sites any sites named in ignore_sites
void apply_ignore_sites( const vector<string> &ignore_sites ,
                         vector<BasePPhoreSite *> &query_sites ) {

  for( int i = 0 , is = ignore_sites.size() ; i < is ; ++i ) {
    bool found_it = false;
    for( int k = 0 , ks = query_sites.size() ; k < ks ; ++k ) {
      if( ignore_sites[i] == query_sites[k]->get_full_name() ) {
        delete query_sites[k];
        query_sites[k] = 0;
        found_it = true;
      }
    }
    if( found_it ) {
      query_sites.erase( remove( query_sites.begin() , query_sites.end() ,
                                 static_cast<BasePPhoreSite *>( 0 ) ) );
      break;
    }
    if( !found_it ) {
      cout << "Warning: ignore site " << ignore_sites[i] << " not found in query."
           << endl;
    }
  }

}

// ***************************************************************************
void read_mol_file( const string &mol_file , OEMolBase *&mol ) {

  oemolistream ims( mol_file );
  if( !ims ) {
    throw DACLIB::FileReadOpenError( mol_file.c_str() );
  }

  // by default, we'll get single-conformer molecules, which is what we
  // want for the query
  mol = OENewMolBase( OEMolBaseType::OEDefault );
  ims >> *mol;
  OEAssignAromaticFlags( *mol , OEAroModelDaylight );

}

// ***************************************************************************
void read_query_file( TriphicSettings &ts ,
                      vector<pair<string,string> > &input_smarts ,
                      vector<pair<string,string> > &smarts_sub_defn ,
                      PharmPoint &pharm_points ,
                      OEMolBase *&query_mol ,
                      vector<BasePPhoreSite *> &query_sites ,
                      vector<SinglePPhoreSite *> &query_score_sites ,
                      vector<pair<string,DACLIB::VolumeGrid *> > &score_vol_grids ) {

  string query_root , query_ext;
  DACLIB::split_filename( ts.query_file() , query_root , query_ext );
  if( query_ext == "sites" || query_ext == "ph4" ) {
    // extension .sites, it's a SITES file.
    if( query_ext == "sites" ) {
      read_query_sites( ts.query_file() , query_sites );
    } else if( query_ext == "ph4" ) {
      read_moe_ph4_file( ts , pharm_points , input_smarts , smarts_sub_defn ,
                         query_sites , score_vol_grids );
    }
    if( query_sites.empty() ) {
      cout << "No sites read." << endl;
      exit( 0 );
    } else {
      cout << "Read " << query_sites.size() << " sites for query." << endl;
    }
    if( ts.sites_score_file().empty() ) {
      // make some empty molecules
      for( int i = 0 , is = query_sites.size() ; i < is ; ++i ) {
        query_mol = OENewMolBase( OEMolBaseType::OEDefault );
      }
    } else {
      // read some structures to be used in scoring
      read_mol_file( ts.sites_score_file() , query_mol );
    }
  } else {
    read_mol_file( ts.query_file() , query_mol );
    create_query_sites( input_smarts , smarts_sub_defn , pharm_points ,
                        query_mol , query_sites , query_score_sites );
  }

  apply_ignore_sites( ts.ignore_sites() , query_sites );

}

// ********************************************************************
void read_protein_file( const string &protein_file ,
                        shared_ptr<OEMolBase> &protein ) {

  if( protein_file.empty() )
    return;

  oemolistream ims( protein_file.c_str() );
  if( !ims )
    throw DACLIB::FileReadOpenError( protein_file.c_str() );
  protein = shared_ptr<OEMolBase>( OENewMolBase( OEMolBaseType::OEDefault ) );
  ims >> *protein;

}

// ********************************************************************
void read_subset_file( const string &subset_file ,
                       vector<string> &subset_names ) {

  ifstream ifs( subset_file.c_str() );
  if( !ifs.good() )
    throw DACLIB::FileReadOpenError( subset_file.c_str() ) ;

  string tmp;
  while( 1 ) {
    ifs >> tmp;
    if( ifs.eof() || !ifs.good() )
      break;
    subset_names.push_back( tmp );
  }

  if( subset_names.empty() ) {
    cerr << "Subset file " << subset_file << " was empty, so it's a short job."
         << endl;
    exit( 1 );
  }

  sort( subset_names.begin() , subset_names.end() );

}

// ********************************************************************
void prepare_query_volume_grids( OEMolBase *query_mol ,
                                 shared_ptr<DACLIB::VolumeGrid> &query_grid ,
                                 shared_ptr<OEMolBase> &protein ,
                                 shared_ptr<DACLIB::VolumeGrid> &protein_grid ) {

  // prepare_mol_grid puts surface and core markers into the grid, so can
  // do single conf against single conf volume and surface volume calcs in
  // one go. Even if we're not sorting or filtering using GRID_SHAPE_TANI,
  // we'll want to calculate the final number of the grid_shape_tani,
  // included and total vols.
  if( query_mol && *query_mol ) {
    query_grid = shared_ptr<DACLIB::VolumeGrid>( DACLIB::prepare_mol_grid( query_mol ) );
  } else {
    query_grid = shared_ptr<DACLIB::VolumeGrid>( static_cast<DACLIB::VolumeGrid *>( 0 ) );
  }

  if( !protein ) {
    return; // no protein
  }

  protein_grid = shared_ptr<DACLIB::VolumeGrid>( DACLIB::prepare_mol_grid( protein.get() ) );

}

// ********************************************************************
void prepare_szybki_optimiser( const TriphicSettings &ts ,
                               shared_ptr<OEMolBase> &protein ,
                               shared_ptr<OESzybki> &szybki ) {

  szybki = shared_ptr<OESzybki>( new OESzybki );
  szybki->SetProtein( *protein ); // rigid by default
  szybki->SetProteinElectrostaticModel( OEProteinElectrostatics::NoElectrostatics );

  if( ts.opt_lig_rigid() ) {
    szybki->SetRunType( OERunType::SolidBodyOpt );
  } else if( ts.opt_lig_flexi() ) {
    szybki->SetRunType( OERunType::TorsionsOpt );
  } else {
    szybki->SetRunType( OERunType::SinglePoint );
  }

}

// ********************************************************************
void add_overlay_scores_to_mol( OEMolBase &mol , OverlayScore &hit_score ) {

  vector<pair<string,string> > scores;
  hit_score.put_scores_in_vector( scores );
  for( int i = 0 , is = scores.size() ; i < is ; ++i ) {
    OESetSDData( mol , scores[i].first , scores[i].second );
  }

}

// ********************************************************************
void overlay_hit_and_store( OEMol &target_mol ,
                            vector<BasePPhoreSite *> &query_sites ,
                            vector<SinglePPhoreSite *> &query_score_sites ,
                            vector<vector<SinglePPhoreSite *> > &target_sites ,
                            shared_ptr<DACLIB::VolumeGrid> &query_grid ,
                            shared_ptr<DACLIB::VolumeGrid> protein_grid ,
                            const vector<pair<string,DACLIB::VolumeGrid *> > &score_vol_grids ,
                            bool use_ring_norms , bool use_h_vectors ,
                            bool use_lps ,
                            OEMolBase *query_mol , bool gauss_shape_tani ,
                            float max_rms ,
                            OverlayScore *hit_score ) {

  // final overlaid results for hit
  OEMolBase *hit_conf = 0;

  try {
    overlay_mols_and_sites( &target_mol , query_sites ,
                            target_sites[hit_score->get_moving_conf()] ,
        hit_score , use_ring_norms , use_h_vectors , use_lps ,
        hit_conf );
  } catch( DACLIB::DiamondOverlayError &e ) {
    cout << e.what() << " : " << target_mol.GetTitle() << endl;
    return;
  }

  if( !hit_conf ) {
    return;
  }

  string hit_name( hit_score->get_moving_mol_name() + "_" +
                   lexical_cast<string>( hit_score->get_moving_conf() + 1 ) );
  hit_conf->SetTitle( hit_name );

  hit_score->calc_volume_scores( hit_conf , query_grid , protein_grid ,
                                 score_vol_grids );
  if( query_mol && gauss_shape_tani ) {
    hit_score->calc_gauss_shape_tanimoto( *query_mol , *hit_conf );
  }
  hit_score->calc_clip_score( query_sites , max_rms );
  hit_score->calc_robins_scores( query_score_sites , query_grid );

  add_overlay_scores_to_mol( *hit_conf , *hit_score );

  delete hit_conf;

}

// ********************************************************************
void write_output_to_mol_and_score_files( OEMolBase *query_mol ,
                                          vector<BasePPhoreSite *> &query_sites ,
                                          list<OverlayScore *> &hits ,
                                          const string &output_file_root ,
                                          const string &output_file_ext ,
                                          const TriphicSettings &ts ) {

  string output_molfile = output_file_root + string( "." ) + output_file_ext;

  oemolostream oms;
  if( !ts.output_scores_only() ) {
    oms.open( output_molfile.c_str() );
    if( !oms ) {
      throw DACLIB::FileWriteOpenError( output_molfile.c_str() );
    }
  }

  string output_scoresfile = output_file_root + string( ".scores" );
  ofstream ofs( output_scoresfile.c_str() );
  if( !ofs.good() ) {
    throw DACLIB::FileWriteOpenError( output_scoresfile.c_str() );
  }
  string sep = ts.comma_output() ? "," : " ";

  // headers to file, including matching points
  if( !hits.empty() ) {
    hits.front()->write_scores_headers_to_stream( ofs , sep );
  } else {
    // put some dummy ones in, so as not to produce an empty file.
    OverlayScore ovs( query_mol->GetTitle() , "Dummy" );
    ovs.write_scores_headers_to_stream( ofs , sep );
  }
  for( int j = 0 , js = query_sites.size() ; j < js ; ++j ) {
    ofs << sep << query_sites[j]->label();
  }
  ofs << sep << "Point_Names " << endl;

  if( ts.output_query_to_hits() ) {
    oms << *query_mol;
  }

  list<OverlayScore *>::iterator p , ps;
  for( p = hits.begin() , ps = hits.end() ; p != ps ; ++p ) {
    if( !ts.output_scores_only() ) {
      OEMolBase *hit_mol = (*p)->get_ov_conf();
      if( ts.no_hit_conf_number() && !ts.do_omega() ) {
        // if we've used omega for the conformations, we number won't have been
        // added in the first place.
        hit_mol->SetTitle( remove_conf_num_from_mol_name( hit_mol->GetTitle() ) );
      }
      oms << *hit_mol;
    }

    ofs << (*p)->get_moving_mol_name();
    if( !ts.no_hit_conf_number() ) {
      ofs << "_" << (*p)->get_moving_conf() + 1;
    }
    ofs << sep;
#ifdef NOTYET
    cout << (*p)->get_moving_mol_name();
    if( !ts.no_hit_conf_number() ) {
      cout << "_" << (*p)->get_moving_conf() + 1;
    }
    cout << sep;
#endif

    (*p)->write_scores_to_stream( ofs , sep );
#ifdef NOTYET
    (*p)->write_scores_to_stream( cout , sep );
    cout << endl;
#endif
    // now the query sites that were matched - 1 for a match, 0 for a miss
    vector<int> ms( query_sites.size() , 0 );
    for( int j = 0 , js = (*p)->num_sites() ; j < js ; ++j ) {
      ms[(*p)->sites()[2*j]] = 1;
    }
    for( int j = 0 , js = ms.size() ; j < js ; ++j ) {
      ofs << sep << ms[j];
    }
    ofs << sep;

    // finally, : separated list of names of query sites that matched
    for( int j = 0 , js = (*p)->num_sites() ; j < js ; ++j ) {
      ofs << query_sites[(*p)->sites()[2*j]]->label();
      if( j < js - 1 ) {
        ofs << ":";
      }
    }
    ofs << endl;
  }

}

// ********************************************************************
void make_interim_filename_root( const string &output_file ,
                                 string &int_file_root ) {

  string output_file_root , output_file_ext;
  DACLIB::split_filename( output_file , output_file_root , output_file_ext );
  int_file_root = output_file_root + string( "_INT" );

}

// ********************************************************************
void tidy_interim_files( const string &output_file ) {

  string int_file_root;
  make_interim_filename_root( output_file , int_file_root );

  string output_file_root , output_file_ext;
  DACLIB::split_filename( output_file , output_file_root , output_file_ext );
  string int_scores_file = int_file_root + ".scores";
  string int_mol_file = int_file_root + string( "." ) + output_file_ext;
  string res_file = output_file_root + string( ".RESTART" );

  using namespace boost::filesystem;
  // boost's remove checks for existence before removing, so we don't
  // need to.
  remove( int_scores_file );
  remove( int_mol_file );
  remove(res_file );

}

// ********************************************************************
void output_results( OEMolBase *query_mol ,
                     vector<BasePPhoreSite *> &query_sites ,
                     list<OverlayScore *> &hits ,
                     TriphicSettings &ts , bool interim ) {

  if( interim ) {
    cout << "Writing interim results : ";
  } else {
    cout << "Writing final results : ";
  }
  if( 1 == hits.size() ) {
    cout << "There was 1 hit for search." << endl;
  } else {
    cout << "There were " << hits.size() << " hits for search." << endl;
  }

  string output_file_root , output_file_ext;
  DACLIB::split_filename( ts.output_file() , output_file_root , output_file_ext );
  if( interim ) {
    make_interim_filename_root( ts.output_file() , output_file_root );
  }

  write_output_to_mol_and_score_files( query_mol , query_sites , hits ,
                                       output_file_root , output_file_ext , ts );

  if( !interim ) {
    // if the final results have been written, we no longer need the
    // interim ones or the restart file.
    tidy_interim_files( ts.output_file() );
  }

}

// ********************************************************************
bool not_smarts_hit( vector<shared_ptr<OESubSearch> > &not_smarts_subs ,
                     OEMol &target_mol ) {

  vector<shared_ptr<OESubSearch> >::iterator p , ps;
  for( p = not_smarts_subs.begin() , ps = not_smarts_subs.end() ;
       p != ps ; ++p ) {
    if( (*p)->SingleMatch( target_mol ) )
      return true;
  }

  return false;

}

// **************************************************************************
// save as many of the current hits as needed in the final hit list.
void store_hits( list<OverlayScore *> &this_query_hits ,
                 vector<BasePPhoreSite *> &query_sites ,
                 vector<SinglePPhoreSite *> &query_score_sites ,
                 vector<vector<SinglePPhoreSite *> > &target_sites ,
                 shared_ptr<DACLIB::VolumeGrid> &query_vol_grid ,
                 shared_ptr<DACLIB::VolumeGrid> &protein_grid ,
                 const vector<pair<string,DACLIB::VolumeGrid *> > &score_vol_grids ,
                 OEMolBase *query_mol , OEMol &target_mol ,
                 TriphicSettings &ts ,
                 list<OverlayScore *> &full_hit_list ) {

  if( this_query_hits.empty() ) {
    return;
  }

#ifdef NOTYET
  cout << query_mol->GetTitle() << endl;
  cout << "All best hits : " << ts.all_best_hits() << " this_query_hits.size() " << this_query_hits.size() << endl;
#endif


  list<OverlayScore *>::iterator p , ps;
  for( p = this_query_hits.begin() , ps = this_query_hits.end() ; p != ps ; ++p ) {
    OverlayScore *this_score = new OverlayScore( *(*p) );
    overlay_hit_and_store( target_mol , query_sites , query_score_sites ,
                           target_sites , query_vol_grid , protein_grid ,
                           score_vol_grids ,
                           ts.ring_norm_usage() == GtplDefs::ALIGN ,
                           ts.h_vector_usage() == GtplDefs::ALIGN ,
                           ts.lp_usage() == GtplDefs::ALIGN ,
                           query_mol , ts.gauss_shape_tanimoto() ,
                           ts.max_rms() , this_score );
    add_overlay_score_to_list( this_score , full_hit_list );
    if( 0 != ts.max_hits() && static_cast<int>( full_hit_list.size() ) > ts.max_hits() ) {
      delete full_hit_list.back();
      full_hit_list.pop_back();
    }
    if( !ts.all_best_hits() ) {
      break; // just want the best hit per target
    }
  }

  // hose them all
  for_each( this_query_hits.begin() , this_query_hits.end() ,
            lambda::bind( lambda::delete_ptr() ,
                          lambda::_1 ) );

}

// **************************************************************************
bool does_query_have_ext_points( const vector<BasePPhoreSite *> &query_sites ) {

  // MOE pharmacophores can have extension points defined without the
  // associated site. If the extension point has a site in close proximity,
  // it will have been merged as a direction into the base site, but if
  // that wasn't possible, we need to add them to the database molecule
  // as well. The sites will be Acc2 or Don2 if this is the case.
  for( int j = 0 , js = query_sites.size() ; j < js ; ++j ) {
    if( string( "Acc2" ) == query_sites[j]->get_type_string() ||
        string( "Don2" ) == query_sites[j]->get_type_string() ) {
      return true;
    }
  }

  return false;

}

// **************************************************************************
// take the sites and make a new set of sites at the positions of the
// virtual sites in the first set, to be used as extension sites when searching
// appropriate queries.
void make_pphore_ext_sites( PharmPoint &pharm_points ,
                            vector<vector<SinglePPhoreSite *> > &sites ,
                            vector<vector<SinglePPhoreSite *> > &ext_sites ) {

  double cds[3];
  double dir[3] = { 0.0 , 0.0 , 0.0 };
  static const string acc2_type_string( "Acc2" );
  static const string don2_type_string( "Don2" );
  static const int acc2_type_code( pharm_points.type_code_from_string( acc2_type_string ) );
  static const int don2_type_code( pharm_points.type_code_from_string( don2_type_string ) );

  const int *type_code = 0;
  const string *type_string = 0;
  for( int i = 0 , is = sites.size() ; i < is ; ++i ) {
    ext_sites.push_back( vector<SinglePPhoreSite *>() );
    for( int j = 0 , js = sites[i].size() ; j < js ; ++j ) {
      if( string( "Acc" ) == sites[i][j]->get_type_string() ||
          string( "Don" ) == sites[i][j]->get_type_string() ) {
        if( string( "Acc" ) == sites[i][j]->get_type_string() ) {
          type_code = &acc2_type_code;
          type_string = &acc2_type_string;
        } else if( string( "Don" ) == sites[i][j]->get_type_string() ) {
          type_code = &don2_type_code;
          type_string = &don2_type_string;
        }
        for( int k = 0 , ks = sites[i][j]->num_virt_sites() ; k < ks ; ++k ) {
          string label = string( "MOE_" ) + lexical_cast<string>( sites[i].size() + ext_sites.back().size() + 1 ) + string( "_" ) + *type_string;
          cds[0] = sites[i][j]->coords()[0] + sites[i][j]->virt_sites()[3*k];
          cds[1] = sites[i][j]->coords()[1] + sites[i][j]->virt_sites()[3*k+1];
          cds[2] = sites[i][j]->coords()[2] + sites[i][j]->virt_sites()[3*k+2];
          SinglePPhoreSite *ns = new SinglePPhoreSite( cds , dir , *type_code , *type_string ,
                                                       label , false ,
                                                       sites[i][j]->num_site_atoms() ,
                                                       sites[i][j]->site_atoms() ,
                                                       sites[i][j]->parent_mol_name() );
          ext_sites.back().push_back( ns );
        }
      }
    }
  }

}

// **************************************************************************
void search_target_molecule( vector<pair<string,string> > &input_smarts ,
                             vector<pair<string,string> > &smarts_sub_defn ,
                             PharmPoint &pharm_points ,
                             OEMolBase *query_mol ,
                             vector<BasePPhoreSite *> &query_sites ,
                             vector<SinglePPhoreSite *> &query_score_sites ,
                             TriphicSettings &ts ,
                             shared_ptr<DACLIB::VolumeGrid> &query_vol_grid ,
                             shared_ptr<DACLIB::VolumeGrid> &protein_grid ,
                             const vector<pair<string,DACLIB::VolumeGrid *> > &score_vol_grids ,
                             shared_ptr<OESzybki> &szybki ,
                             vector<shared_ptr<OESubSearch> > &not_smarts_subs ,
                             OEMol &target_mol ,
                             list<OverlayScore *> &hit_list ) {

#ifdef NOTYET
  cout << "Searching " << target_mol.GetTitle() << endl;
#endif

  static bool need_ext_points = false , first_time = true;
  if( first_time ) {
    need_ext_points = does_query_have_ext_points( query_sites );
    first_time = false;
  }
  
  // make sure this molecule doesn't hit something in not_smarts_subs - if it
  // does, we don't want it.
  if( not_smarts_hit( not_smarts_subs , target_mol ) ) {
    return;
  }

  vector<vector<SinglePPhoreSite *> > target_sites;
  vector<vector<SinglePPhoreSite *> > target_ext_sites;
  try {
    DACLIB::make_pphore_sites( target_mol , pharm_points , input_smarts ,
                               smarts_sub_defn , target_sites );
    if( need_ext_points ) {
      make_pphore_ext_sites( pharm_points , target_sites , target_ext_sites );
    }
  } catch( DACLIB::SMARTSDefnError &e ) {
    cerr << e.what() << endl;
    exit( 1 );
  } catch( DACLIB::DiamondOverlayError &e ) {
    cerr << "Failed to create pharmacophore sites correctly for " << target_mol.GetTitle() << endl;
    cout << "Failed to create pharmacophore sites correctly for " << target_mol.GetTitle() << endl;
    cout << e.what() << endl;
    return;
  }

  vector<vector<int> > cliques;
  vector<vector<SinglePPhoreSite *> > full_target_sites;

  list<OverlayScore *> clique_overlays;
  for( int j = 0 , js = target_sites.size() ; j < js ; ++j ) {

    full_target_sites.push_back( target_sites[j] );
    if( need_ext_points ) {
      full_target_sites.back().insert( full_target_sites.back().end() , target_ext_sites[j].begin() ,
                                       target_ext_sites[j].end() );
    }
    if( static_cast<int>( full_target_sites[j].size() ) < ts.min_clique_size() ) {
      continue; // no point doing more
    }
#ifdef NOTYET
    cout << "Number target sites for conf " << j << " of " << target_mol.GetTitle()
         << " = " << full_target_sites.back().size() << endl;
    cout << "Full target sites for conf " << j << " of " << target_mol.GetTitle() << endl;
    BOOST_FOREACH( SinglePPhoreSite *ts , full_target_sites.back() ) {
      ts->brief_report( cout );
    }
    cout << endl;
#endif

    find_cliques( query_sites , full_target_sites[j] , ts.dist_tol() ,
                  ts.scaled_dist_tol() , ts.dont_do_sub_cliques() ,
                  ts.min_clique_size() , cliques );
#ifdef NOTYET
    cout << "Number of cliques for conf " << j << " of " << target_mol.GetTitle()
         << " = " << cliques.size() << endl;
#endif
    score_and_store_cliques( query_mol->GetTitle() , 0 ,
                             query_mol , query_sites ,
                             query_score_sites , full_target_sites[j] ,
                             &target_mol , j , query_vol_grid ,
                             protein_grid , score_vol_grids ,
                             cliques ,
                             ts.ring_norm_usage() , ts.h_vector_usage() ,
                             ts.lp_usage() , ts.ring_norm_tol() ,
                             ts.h_vector_tol() , ts.lp_tol() , true ,
                             ts.min_clique_size() , true , ts.max_rms() ,
                             ts.do_grid_shape_tani_cutoff() ,
                             ts.min_grid_shape_tani() ,
                             ts.do_gauss_shape_tani_cutoff() ,
                             ts.min_gauss_shape_tani() ,
                             ts.do_surf_vol_cutoff() , ts.min_surf_vol() ,
                             ts.do_inc_vol_cutoff() , ts.min_inc_vol() ,
                             ts.do_prot_clash_cutoff() , ts.max_prot_clash() ,
                             ts.do_mmff_nrg_cutoff() , ts.max_mmff_nrg() ,
                             ts.req_sites_or() , ts.req_sites_and() ,
                             ts.req_points_or() , ts.req_points_and() ,
                             szybki , clique_overlays );
  }

  store_hits( clique_overlays , query_sites , query_score_sites ,
              full_target_sites , query_vol_grid , protein_grid ,
              score_vol_grids , query_mol ,	target_mol , ts ,
              hit_list );

  for( int i = 0 , is = target_sites.size() ; i < is ; ++i ) {
    for( int j = 0 , js = target_sites[i].size() ; j < js ; ++j ) {
      delete target_sites[i][j];
    }
  }
  for( int i = 0 , is = target_ext_sites.size() ; i < is ; ++i ) {
    for( int j = 0 , js = target_ext_sites[i].size() ; j < js ; ++j ) {
      delete target_ext_sites[i][j];
    }
  }

}

// **************************************************************************
void report_progress( unsigned int &i , unsigned int chunk_size ,
                      unsigned int restart_point ,
                      unsigned int &report_step ,
                      list<OverlayScore *> &hit_list , ostream &os ) {

  i += chunk_size;
  unsigned int i_tmp = ( i / report_step ) * report_step;

  if( i - i_tmp < chunk_size ) {
    if( 1 == i + restart_point ) {
      os << "Done 1 molecule so far with ";
    } else {
      os << "Done " << i + restart_point << " molecules so far";
      if( restart_point ) {
        cout << " (" << i << " since restart)";
      }
      cout << " with ";
    }
    os << hit_list.size() << " hit";
    if( hit_list.size() != 1 ) {
      os << "s";
    }
    if( !hit_list.empty() ) {
      os << " against " << hit_list.front()->get_fixed_mol_name()
         << "." << endl;
    } else {
      os << "." << endl;
    }
  }
  if( i >= report_step * 10 && report_step < 100000 ) {
    report_step *= 10;
  }

}

// **************************************************************************
void report_progress( unsigned int &i , unsigned int chunk_size ,
                      unsigned int &report_step ,
                      OEMolBase &target_mol ,
                      list<OverlayScore *> &hit_list ,
                      ostream &os ) {

  i += chunk_size;
  unsigned int i_tmp = ( i / report_step ) * report_step;
  if( i - i_tmp < chunk_size ) {
    if( 1 == i ) {
      os << "Done 1 molecule so far with ";
    } else {
      os << "Done " << i << " molecules (to " << target_mol.GetTitle()
         << ") so far with ";
    }
    os << hit_list.size() << " hit";
    if( hit_list.size() != 1 ) {
      os << "s";
    }
    if( !hit_list.empty() ) {
      os << " against " << hit_list.front()->get_fixed_mol_name()
         << "." << endl;
    } else {
      os << "." << endl;
    }
  }
  if( i >= report_step * 10 && report_step < 100000 ) {
    report_step *= 10;
  }

}

// ********************************************************************
void write_restart_file( TriphicSettings &ts , unsigned int mol_num ,
                         list<OverlayScore *> &hits ) {

  string output_file_root , output_file_ext;
  DACLIB::split_filename( ts.output_file() , output_file_root , output_file_ext );
  string res_file = output_file_root + string( ".RESTART" );
  ofstream ofs( res_file.c_str() );
  ofs << "<START_MOL>" << mol_num << "</START_MOL>" << endl;
  ofs << "<HIT_SET_SIZE>" << hits.size() << "</HIT_SET_SIZE>" << endl;
  for( list<OverlayScore *>::iterator p = hits.begin() ; p != hits.end() ; ++p ) {
    ofs << "<NEXT_HIT>" << endl;
    ofs << (*p)->write_contents_to_string();
    ofs << "</NEXT_HIT>" << endl;
  }

}

// ********************************************************************
// for parallel jobs, lowest_mol_done is the important one for restarts
void dump_results_so_far( unsigned int mols_done ,
                          unsigned int lowest_mol_done ,
                          unsigned int report_step ,
                          OEMolBase *query_mol ,
                          vector<BasePPhoreSite *> &query_sites ,
                          list<OverlayScore *> &hits ,
                          TriphicSettings &ts ) {

  // interim results need to be done fairly frequently.
  if( report_step > 10000 ) {
    report_step = 10000;
  }
  if( mols_done > 10 && !( mols_done % report_step ) ) {
    output_results( query_mol , query_sites , hits , ts , true );
    write_restart_file( ts , lowest_mol_done , hits );
  }

}

// ********************************************************************
void send_results_to_master( list<OverlayScore *> &hits ) {

  DACLIB::mpi_send_string( string( "Results" ) , 0 );

  unsigned int j = hits.size();
  MPI_Send( &j , 1 , MPI_UNSIGNED , 0 , 0 , MPI_COMM_WORLD );
  list<OverlayScore *>::iterator p , ps;
  for( p = hits.begin() , ps = hits.end() ; p != ps ; ++p ) {
    DACLIB::mpi_send_string( (*p)->write_contents_to_string() , 0 );
  }

}

// ********************************************************************
void do_omega_if_necessary( vector<pOEMol> &prepped_target_mols ,
                            shared_ptr<OEOmega> &omega ,
                            bool do_omega , bool do_flipper ,
                            bool do_warts ,
                            unsigned int max_versions ,
                            const string &orig_mol_title ) {

  vector<pOEMol> out_mols;

  // first do flipper if necessary for unspecified stereo centres
  if( do_flipper ) {
    for( int i = 0 , is = prepped_target_mols.size() ; i < is ; ++i ) {

      OESetDimensionFromCoords( *prepped_target_mols[i] );

      int j = 1;
      for( OESystem::OEIter<OEMolBase> stereo = OEFlipper( *prepped_target_mols[i] ) ; stereo ; ++stereo , ++j ) {
        out_mols.push_back( pOEMol( new OEMol( *stereo ) ) );
        if( do_warts ) {
          string new_name = out_mols.back()->GetTitle() + string( "_f" ) + lexical_cast<string>( j );
          out_mols.back()->SetTitle( new_name );
        }
      }
    }
  } else {
    for( int i = 0 , is = prepped_target_mols.size() ; i < is ; ++i ) {
      OESetDimensionFromCoords( *prepped_target_mols[i] );
      out_mols.push_back( pOEMol( new OEMol( *prepped_target_mols[i] ) ) );
    }
  }

  if( out_mols.size() > max_versions ) {
    cout << orig_mol_title << " exceeded maximum number of versions ("
         << out_mols.size() << " vs " << max_versions << ")." << endl;
    prepped_target_mols.clear();
    out_mols.clear();
    return;
  }

  for( int i = 0 , is = prepped_target_mols.size() ; i < is ; ++i ) {

    if( do_omega || ( prepped_target_mols[i]->NumAtoms() &&
                      3 != prepped_target_mols[i]->GetDimension() ) ) {
      if( !omega && OEOmegaIsLicensed( "toolkit" ) ) {
        omega = shared_ptr<OEOmega>( new OEOmega );
      }
      (*omega)( *(out_mols[i]) );
      OESetDimensionFromCoords( *out_mols[i] );
      if( 1 == out_mols[i]->NumConfs() && 3 != out_mols[i]->GetDimension() ) {
        cout << "Omega failed for " << prepped_target_mols[i]->GetTitle() << " ("
             << DACLIB::create_cansmi( *prepped_target_mols[i] ) << ") so skipping." << endl;
        cerr << "Omega failed for " << prepped_target_mols[i]->GetTitle() << " so skipping." << endl;
        prepped_target_mols.clear();
        return;
      }
    } else {
      out_mols.push_back( prepped_target_mols[i] );
    }
  }

  prepped_target_mols.clear();
  string last_report;
  int num_confs = 0;

  for( int i = 0 , is = out_mols.size() ; i < is ; ++i ) {
    OESetDimensionFromCoords( *out_mols[i] );
    if( 3 == out_mols[i]->GetDimension() ) {
      // if not 3, omega failed for some unknown reason
      prepped_target_mols.push_back( out_mols[i] );
      num_confs += out_mols[i]->NumConfs();
    } else {
      if( last_report != out_mols[i]->GetTitle() ) {
        cout << "Omega failed for " << out_mols[i]->GetTitle() << " for some unknown reason" << endl;
        last_report = out_mols[i]->GetTitle();
      }
    }
  }

}

// ********************************************************************
void prepare_search_molecule( OEMol &target_mol , shared_ptr<OEOmega> &omega ,
                              bool do_omega , bool do_flipper ,
                              bool do_ions_and_tauts , bool do_warts ,
                              unsigned int max_versions ,
                              vector<pOEMol> &prepped_target_mols ) {

  if( do_ions_and_tauts ) {
    // Suppress irritating warnings from OELibraryGen (and everything else, of course, but
    // OELibraryGen gives a lot of very irritating stuff)
    OESystem::OEThrow.SetLevel( OESystem::OEErrorLevel::Error );
    shared_ptr<OEMolBase> sc_target_mol( OENewMolBase( target_mol.SCMol() ,
                                                       OEMolBaseType::OEDefault ) );
    vector<OEMolBase *> tm = enumerate_tautomers( *sc_target_mol , string( "" ) ,
                                                  string( "" ) );

    vector<pOEMol> taut_mols;
    for( int i = 0 , is = tm.size() ; i < is ; ++i ) {
      taut_mols.push_back( pOEMol( new OEMol( *tm[i] ) ) );
      delete tm[i];
    }
    tm.clear();
    if( taut_mols.size() > max_versions ) {
      cout << target_mol.GetTitle() << " exceeded maximum number of versions ("
           << taut_mols.size() << " vs " << max_versions << ")." << endl;
      prepped_target_mols.clear();
      taut_mols.clear();
      // put it back like it was
      OESystem::OEThrow.SetLevel( OESystem::OEErrorLevel::Warning );
      return;
    }

    vector<pOEMol> all_tauts;
    for( int i = 0 , is = taut_mols.size() ; i < is ; ++i ) {
      if( do_warts ) {
        string new_name = taut_mols[i]->GetTitle() + string( "_t" ) + lexical_cast<string>( i + 1 );
        taut_mols[i]->SetTitle( new_name );
      }
      vector<OEMolBase *> ti = enumerate_ions( *taut_mols[i] , string( "" ) ,
                                               string( "" ) , string( "" ) );
      vector<pOEMol> these_ions;
      for( int ii = 0 , iis = ti.size() ; ii < iis ; ++ii ) {
        these_ions.push_back( pOEMol( new OEMol( *ti[ii] ) ) );
        delete ti[ii];
      }
      ti.clear();
      for( int j = 0 , js = these_ions.size() ; j < js ; ++j ) {
        if( do_warts ) {
          string new_name = these_ions[j]->GetTitle() + string( "_i" ) + lexical_cast<string>( j + 1 );
          these_ions[j]->SetTitle( new_name );
        }
        all_tauts.push_back( these_ions[j] );
      }
    }
    if( all_tauts.size() > max_versions ) {
      cout << target_mol.GetTitle() << " exceeded maximum number of versions ("
           << all_tauts.size() << " vs " << max_versions << ")." << endl;
      prepped_target_mols.clear();
      all_tauts.clear();
      // put it back like it was
      OESystem::OEThrow.SetLevel( OESystem::OEErrorLevel::Warning );
      return;
    }

    for( int i = 0 , is = all_tauts.size() ; i < is ; ++i ) {
      prepped_target_mols.push_back( all_tauts[i] );
    }
    // put it back like it was
    OESystem::OEThrow.SetLevel( OESystem::OEErrorLevel::Warning );
  } else {
    prepped_target_mols.push_back( pOEMol( new OEMol( target_mol ) ) );
  }

  do_omega_if_necessary( prepped_target_mols , omega ,
                         do_omega , do_flipper , do_warts ,
                         max_versions , target_mol.GetTitle() );

#ifdef NOTYET
  cout << target_mol.GetTitle() << " produced " << prepped_target_mols.size()
       << " molecules." << endl;
#endif

}

// ********************************************************************
// empty the individual hit lists, not the vector
void empty_hit_list( list<OverlayScore *> &hit_list ,
                     bool delete_hits ) {

  if( delete_hits ) {
    list<OverlayScore *>::iterator p , ps;
    for( p = hit_list.begin() , ps = hit_list.end() ; p != ps ; ++p ) {
      delete *p;
    }
  }
  hit_list.clear();

}

// **************************************************************************
// search the next_mol'th molecule in the databases
void search_database( unsigned int next_mol , int step_size ,
                      vector<pair<string,string> > &input_smarts ,
                      vector<pair<string,string> > &smarts_sub_defn ,
                      PharmPoint &pharm_points ,
                      OEMolBase *query_mol ,
                      vector<BasePPhoreSite *> &query_sites ,
                      vector<SinglePPhoreSite *> &query_score_sites ,
                      const vector<string> &mol_subset ,
                      TriphicSettings &ts ,
                      shared_ptr<DACLIB::VolumeGrid> &query_vol_grid ,
                      shared_ptr<DACLIB::VolumeGrid> &protein_grid ,
                      const vector<pair<string,DACLIB::VolumeGrid *> > &score_vol_grids ,
                      shared_ptr<OESzybki> &szybki , shared_ptr<OEOmega> &omega ,
                      vector<shared_ptr<OESubSearch> > &not_smarts_subs ,
                      list<OverlayScore *> &hit_list ) {

  for( int i = 0 ; i < step_size ; ++i ) {

    pOEMol tmol = DACLIB::read_nth_mol_from_oemolistream( next_mol + i , ts.db_files() ,
                                                          ts.single_conf_mols() );

    if( !tmol ) {

      // send any hits to master
      send_results_to_master( hit_list );
      empty_hit_list( hit_list , true );
      // we're done
      DACLIB::mpi_send_string( string( "Search_Finished" ) , 0 );
      return;

    }

#ifdef NOTYET
    int world_rank;
    MPI_Comm_rank( MPI_COMM_WORLD , &world_rank );
    cout << world_rank << " Now Searching " << next_mol + i << " : " << tmol->GetTitle() << " : "
         << tmol->NumConfs() << endl;
#endif

    if( mol_subset.empty() ||
        binary_search( mol_subset.begin() , mol_subset.end() , tmol->GetTitle() ) ) {

      vector<pOEMol> prepped_target_mols;
      prepare_search_molecule( *tmol , omega , ts.do_omega() ,
                               ts.do_flipper() , ts.do_ions_and_tauts() ,
                               !ts.no_warts() , ts.max_versions() ,
                               prepped_target_mols );

      for( int j = 0 , js = prepped_target_mols.size() ; j < js ; ++j ) {
        search_target_molecule( input_smarts , smarts_sub_defn , pharm_points ,
                                query_mol , query_sites , query_score_sites ,
                                ts , query_vol_grid , protein_grid ,
                                score_vol_grids , szybki , not_smarts_subs ,
                                *prepped_target_mols[j] , hit_list );
      }
    }
  }

  // send any hits to master
  send_results_to_master( hit_list );
  empty_hit_list( hit_list , true );
  // ask for next molecule
  DACLIB::mpi_send_string( string( "Send_Mol_Num" ) , 0 );

}

// ********************************************************************
void make_not_smarts_queries( const vector<pair<string,string> > &not_smarts_list ,
                              vector<pair<string,string> > &smarts_sub_defn ,
                              vector<shared_ptr<OESubSearch> > &not_smarts_subs ) {

  for( int i = 0 , is = not_smarts_list.size() ; i < is ; ++i ) {
    string exp_smarts( not_smarts_list[i].second );
    // do any vector binding expansion
    if( string::npos != exp_smarts.find( "$" ) &&
        !OESmartsLexReplace( exp_smarts , smarts_sub_defn ) ) {
      // We shouldn't have got this far without all the necessary sub-definitions
      // so just bail
      ostringstream oss;
      oss << "AWOOGA - something screwed with the SMARTS sub-definitions for "
          << endl << not_smarts_list[i].first << " : "
          << not_smarts_list[i].second << endl
          << "in make_not_smarts_queries()" << endl;
      throw( DACLIB::SMARTSDefnError( oss.str().c_str() ) );
    }
    not_smarts_subs.push_back( shared_ptr<OESubSearch>( new OESubSearch( exp_smarts.c_str() ) ) );
  }

}

// ****************************************************************************
// overlay the query onto itself to get the best possible score for a hit,
// which will be used in the Pareto rankings.
OverlayScore *make_reference_overlay( OEMolBase *query_mol ,
                                      vector<BasePPhoreSite *> &query_sites ,
                                      shared_ptr<DACLIB::VolumeGrid> &query_vol_grid ,
                                      shared_ptr<DACLIB::VolumeGrid> &protein_grid ,
                                      vector<pair<string,DACLIB::VolumeGrid *> > &score_vol_grids ) {

  vector<int> clique = vector<int>( 2 * query_sites.size() , 0 );
  vector<SinglePPhoreSite *> query_score_sites , query_sites2;
  for( int i = 0 , is = query_sites.size() ; i < is ; ++i ) {
    clique[2*i] = clique[2*i+1] = i;
    query_score_sites.push_back( polymorphic_downcast<SinglePPhoreSite *>( query_sites[i] ) );
    query_sites2.push_back( polymorphic_downcast<SinglePPhoreSite *>( query_sites[i] ) );
  }

  OEMol query_cp( *query_mol );
  shared_ptr<OESz::OESzybki> szybki( static_cast<OESz::OESzybki *>( 0 ) );
  try {
    OverlayScore *ret_val = new OverlayScore( string( query_mol->GetTitle() ) , string( query_mol->GetTitle() ) ,
                                              0 , 0 , clique , query_sites ,
                                              query_score_sites , query_sites2 ,
                                              *query_mol , query_cp ,
                                              query_vol_grid , protein_grid ,
                                              score_vol_grids , szybki ,
                                              false , false , false , false );
#ifdef NOTYET
    cout << "Scores for reference overlay" << endl;
    cout << ret_val->rms() << " : " << ret_val->num_sites() << endl;
    cout << ret_val->hphobe_score() << " : " << ret_val->hbond_score() << " : "
         << ret_val->vol_score() << endl;
#endif
    return ret_val;
  } catch( OverlayScoreError &e ) {
    cout << "Error : " << e.what() << endl;
    exit( 1 );
  }


}

// ****************************************************************************
void setup_search( TriphicSettings &ts ,
                   vector<pair<string,string> > &input_smarts ,
                   vector<pair<string,string> > &smarts_sub_defn ,
                   PharmPoint &pharm_points ,
                   OEMolBase *&query_mol ,
                   vector<BasePPhoreSite *> &query_sites ,
                   vector<SinglePPhoreSite *> &query_score_sites ,
                   shared_ptr<OEMolBase> &protein ,
                   vector<string> &mol_subset ,
                   shared_ptr<DACLIB::VolumeGrid> &query_vol_grid ,
                   shared_ptr<DACLIB::VolumeGrid> &protein_grid ,
                   vector<pair<string,DACLIB::VolumeGrid *> > &score_vol_grids ,
                   shared_ptr<OESzybki> &szybki ,
                   shared_ptr<OEOmega> &omega ,
                   vector<shared_ptr<OESubSearch> > &not_smarts_subs ,
                   bool report_sites ) {

  if( !ts.smarts_file().empty() && !ts.points_file().empty() ) {
    read_smarts_file( ts.smarts_file() , input_smarts , smarts_sub_defn );
    if( !ts.not_smarts_file().empty() ) {
      read_smarts_file( ts.not_smarts_file() , ts.not_smarts_list() ,
                        smarts_sub_defn );
    }
    pharm_points.read_points_file( ts.points_file() );
  } else {
    pharm_points.clear_data();
    pharm_points.read_points_xml_string( DACLIB::DEFAULT_POINT_DEFS );
    ParseSMARTSXML psx;
    psx.parse_string( DACLIB::DEFAULT_POINT_DEFS , input_smarts , smarts_sub_defn );
  }

  if( ts.score_method() == GtplDefs::ROBINS_SCORE_PARETO &&
      ( !pharm_points.has_point_of_name( string( "Hydrophobe" ) ) ||
        !pharm_points.has_point_of_name( string( "Acceptor" ) ) ||
        !pharm_points.has_point_of_name( string( "Donor" ) ) ) ) {
    cerr << "For a ROBINS_PARETO score, you need points of type Hydrophobe, Acceptor and Donor (case insensitive)." << endl;
    exit( 1 );
  }

  // read the query file and, if necessary, create the sites and volumes
  DACLIB::VolumeGrid::set_grid_spacing( ts.vol_grid_spacing() );
  read_query_file( ts , input_smarts , smarts_sub_defn , pharm_points ,
                   query_mol , query_sites , query_score_sites ,
                   score_vol_grids );

  if( report_sites ) {
    brief_report_sites( cout , query_sites );
    cout << endl;
    ts.print_required_sites_and_points( cout );
  }

  if( ts.print_sites_and_stop() ) {
    return;
  }

  if( query_sites.empty() ) {
    cerr << "No query sites found." << endl;
    cout << "No query sites found." << endl;
    exit( 0 );
  }

  read_protein_file( ts.protein_file() , protein );

  if( !ts.subset_file().empty() ) {
    read_subset_file( ts.subset_file() , mol_subset );
  }

  // these are arbitrary grids for volume scores, prepared outside the program
  // probably by psg or prepare_scoring_grid
  DACLIB::read_grids( ts.grid_vol_files() , score_vol_grids );

  // get the scoring volume grids. If score_vol_grids is not empty, these
  // grids must dictate the grid spacing, otherwise make it coarser than
  // the default because the volume score becomes the rate-limiting step
  if( !score_vol_grids.empty() ) {
    if( ts.vol_grid_spacing() != DACLIB::VolumeGrid::get_grid_spacing() ) {
      cerr << "Warning : volume grid spacing set at " << DACLIB::VolumeGrid::get_grid_spacing()
           << " by scoring volume files over-rides user request." << endl;
    }
  }

  prepare_query_volume_grids( query_mol , query_vol_grid , protein ,
                              protein_grid );

  if( OESzybkiIsLicensed() && protein &&
      ( ts.max_mmff_nrg() < numeric_limits<float>::max() ||
        ts.opt_lig_flexi() || ts.opt_lig_rigid() ) ) {
    prepare_szybki_optimiser( ts , protein , szybki );
  }

  if( ts.do_omega() ) {
    omega = shared_ptr<OEOmega>( new OEOmega );
  }

  // SMARTS that the hits can't match
  make_not_smarts_queries( ts.not_smarts_list() , smarts_sub_defn ,
                           not_smarts_subs );

  OverlayScore::set_score_method( ts.score_method() );
  OverlayScore *ref_ov = make_reference_overlay( query_mol , query_sites ,
                                                 query_vol_grid , protein_grid ,
                                                 score_vol_grids );
  OverlayScore::set_reference_overlay( ref_ov );

}

// ****************************************************************************
unsigned int read_restart_value( const string &full_line ,
                                 const string &tag ,
                                 const string &res_file ) {

  string start_tag = string( "<" ) + tag + string( ">" );
  string end_tag = string( "</" ) + tag + string( ">" );

  unsigned int tag_val;
  if( full_line.substr( 0 , start_tag.length() ) == start_tag ) {
    size_t n = full_line.find( end_tag );
    if( n != string::npos ) {
      tag_val = lexical_cast<unsigned int>( full_line.substr( start_tag.length() ,
                                                              n - end_tag.length() + 1 ) );
    } else {
      cerr << "Error reading restart file " << res_file << " : missing end tag in line : " << full_line << endl;
      exit( 1 );
    }
  } else {
    cerr << "Error reading restart file " << res_file << " : missing start tag in line : " << full_line << endl;
    exit( 1 );
  }

  return tag_val;

}

// ********************************************************************
void read_restart_hits( TriphicSettings &ts ,
                        unsigned int &mol_num ,
                        list<OverlayScore *> &hit_list ) {

  string output_file_root , output_file_ext;
  DACLIB::split_filename( ts.output_file() , output_file_root , output_file_ext );
  string res_file = output_file_root + string( ".RESTART" );
  ifstream ifs( res_file.c_str() );
  if( !ifs ){
    cerr << "Couldn't open restart file " << res_file << "." << endl;
    exit( 1 );
  }

  string next_line;
  getline( ifs , next_line );
  mol_num = read_restart_value( next_line , string( "START_MOL" ) , res_file );

  hit_list = list<OverlayScore *>();
  getline( ifs , next_line );
  unsigned int hit_set_size = read_restart_value( next_line , string( "HIT_SET_SIZE" ) , res_file );
  for( unsigned int j = 0 ; j < hit_set_size ; ++j ) {
    string hit_string;
    getline( ifs , next_line );
    while( true ) {
      getline( ifs , next_line );
      if( next_line == string( "</NEXT_HIT>" ) ) {
        break;
      }
      hit_string += next_line + string( "\n" );
    }
    hit_list.push_back( new OverlayScore( hit_string ) );
  }

}

// ********************************************************************
int serial_triphic_search( TriphicSettings &ts ) {

  vector<pair<string,string> > input_smarts , smarts_sub_defn;
  PharmPoint pharm_points;
  OEMolBase *query_mol;
  vector<BasePPhoreSite *> query_sites;
  vector<SinglePPhoreSite *> query_score_sites;
  shared_ptr<OEMolBase> protein;
  vector<string> mol_subset;
  shared_ptr<DACLIB::VolumeGrid> query_vol_grid;
  shared_ptr<DACLIB::VolumeGrid> protein_grid;
  vector<pair<string,DACLIB::VolumeGrid *> > score_vol_grids;
  shared_ptr<OESzybki> szybki;
  shared_ptr<OEOmega> omega;
  vector<shared_ptr<OESubSearch> > not_smarts_subs;
  oemolistream *ims = 0;

  setup_search( ts , input_smarts , smarts_sub_defn , pharm_points ,
                query_mol , query_sites , query_score_sites , protein ,
                mol_subset , query_vol_grid , protein_grid , score_vol_grids ,
                szybki , omega , not_smarts_subs , true );

  if( ts.print_sites_and_stop() ) {
    return 0;
  }

  ts.print_required_sites_and_points( cout );
  if( !ts.check_required_sites( query_sites ) ||
      !ts.check_required_points( query_sites ) ) {
    return 1;
  }
  ts.set_chunk_size( 1 );

  OEMol target_mol;
  list<OverlayScore *> hit_list;

  unsigned int start_mol_num = 0 , i = 0;
  if( ts.restart() ) {
    read_restart_hits( ts , start_mol_num , hit_list );
    if( ts.restart_number() ) {
      start_mol_num = ts.restart_number(); // over-ride what was in the file
    }
    cout << "Search restarts at molecule " << start_mol_num << "." << endl;
  }

  vector<string> db_files = ts.db_files();

  BOOST_FOREACH( string db_file , db_files ) {

    cout << "Searching database file : " << db_file << endl;

    unsigned int report_step = 1;
    ims = new oemolistream;
    DACLIB::open_databasefile( db_file , ts.single_conf_mols() , ims );

    while( *ims >> target_mol ) {
#ifdef NOTYET
      cout << "Searching " << target_mol.GetTitle() << " : " << target_mol.NumConfs() << endl;
#endif

      if( i >= start_mol_num && ( mol_subset.empty() ||
                                  binary_search( mol_subset.begin() , mol_subset.end() , target_mol.GetTitle() ) ) ) {

        vector<pOEMol> prepped_target_mols;
        prepare_search_molecule( target_mol , omega , ts.do_omega() ,
                                 ts.do_flipper() , ts.do_ions_and_tauts() ,
                                 !ts.no_warts() , ts.max_versions() ,
                                 prepped_target_mols );

        for( int j = 0 , js = prepped_target_mols.size() ; j < js ; ++j ) {
          search_target_molecule( input_smarts , smarts_sub_defn , pharm_points ,
                                  query_mol , query_sites , query_score_sites ,
                                  ts , query_vol_grid , protein_grid ,
                                  score_vol_grids , szybki , not_smarts_subs ,
                                  *prepped_target_mols[j]  , hit_list );
        }
      }
      report_progress( i , ts.chunk_size() , report_step , target_mol , hit_list , cout );
      dump_results_so_far( i , i , report_step , query_mol , query_sites ,
                           hit_list , ts );
    }
    ims->close();
    delete ims;

  }

  try {
    output_results( query_mol , query_sites , hit_list , ts , false );
  } catch( DACLIB::FileWriteOpenError &e ) {
    cerr << e.what() << endl;
  }

  // tidy stuff up, mostly to stop valgrind whinging. triphic dates from before
  // I discovered the joys of boost smart pointers.
  empty_hit_list( hit_list , true );

  return 0;

}

// ********************************************************************
void send_search_details_to_slaves( TriphicSettings &ts ,
                                    int world_size ) {

  send_openeye_license_to_slaves( world_size );
  send_cwd_to_slaves( world_size );

  for( int i = 1 ; i < world_size ; ++i ) {
    DACLIB::mpi_send_string( string( "Search_Details" ) , i );
    ts.send_contents_via_mpi( i );
  }

}

// ********************************************************************
void receive_results_from_slave( int slave_rank ,
                                 list<OverlayScore *> &slave_hit_list ) {

  unsigned int num_hits_to_rec;
  MPI_Recv( &num_hits_to_rec , 1 , MPI_UNSIGNED , slave_rank , 0 ,
            MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  for( unsigned int j = 0 ; j < num_hits_to_rec ; ++j ) {
    string hit_msg;
    DACLIB::mpi_rec_string( slave_rank , hit_msg );
    slave_hit_list.push_back( new OverlayScore( hit_msg ) );
#ifdef NOTYET
    slave_hit_list.back()->write_scores_to_stream( cout , " , " );
    cout << endl;
#endif
  }

}

// ********************************************************************
void merge_hit_lists( unsigned int max_hits ,
                      list<OverlayScore *> &sub_hit_list ,
                      list<OverlayScore *> &full_hit_list ) {

#ifdef NOTYET
  cout << "Merging " << sub_hit_list.size() << " hits" << endl;
#endif
  list<OverlayScore *>::iterator p , ps;
  for( p = sub_hit_list.begin() , ps = sub_hit_list.end() ; p != ps ; ++p ) {
    add_overlay_score_to_list( *p , full_hit_list );
    // add_overlay_score_to_list no longer maintains the list size.
    if( max_hits && full_hit_list.size() > max_hits ) {
      delete full_hit_list.back();
      full_hit_list.pop_back();
    }
  }
#ifdef NOTYET
    cout << "Number of hits now " << full_hit_list.size() << endl;
#endif

}

// ********************************************************************
void send_parallel_searches( int world_size , TriphicSettings &ts ,
                             list<OverlayScore *> &hit_list ,
                             OEMolBase *query_mol ,
                             vector<BasePPhoreSite *> &query_sites ,
                             unsigned int curr_mol ) {

  int num_finished = 0;
  int num_slaves = world_size - 1; // process 0 is the master
  unsigned int report_step = 1;
  unsigned int num_mols_done = 0;
  unsigned int restart_point = curr_mol;

  vector<unsigned int> mols_being_done( world_size , 0 );
  vector<unsigned int> last_mol_done( world_size , 0 );

  while( num_finished < num_slaves ) {

    MPI_Status status;
    MPI_Probe( MPI_ANY_SOURCE , 0 , MPI_COMM_WORLD , &status );
    string msg;
    DACLIB::mpi_rec_string( status.MPI_SOURCE , msg );

#ifdef NOTYET
    cout << "Message : " << msg << " from " << status.MPI_SOURCE << endl;
    cout << "mols_being_done : ";
    copy( mols_being_done.begin() + 1 , mols_being_done.end() , uintOut );
    cout << endl;
    cout << "last_mol_done : ";
    copy( last_mol_done.begin() + 1 , last_mol_done.end() , uintOut );
    cout << endl;
#endif

    if( string( "Results" ) == msg ) {
#ifdef NOTYET
      cout << "Results for molecule " << mols_being_done[status.MPI_SOURCE] << " from slave "
           << status.MPI_SOURCE << endl;
      cout << "mols being done : ";
      copy( mols_being_done.begin() + 1 , mols_being_done.end() , uintOut );
      cout << endl;
#endif
      list<OverlayScore *> slave_hit_list;
      receive_results_from_slave( status.MPI_SOURCE , slave_hit_list );
      merge_hit_lists( ts.max_hits() , slave_hit_list , hit_list );
      empty_hit_list( slave_hit_list , false );

#ifdef NOTYET
      cout << "YYYYYYYYYYYYYYYYYYY" << endl;
      cout << "Hits list now : " << endl;
      list<OverlayScore *> &hl = hit_lists.front();
      for( list<OverlayScore *>::iterator p = hl.begin() ; p != hl.end() ; ++p ) {
        cout << (*p)->get_moving_mol_name() << " : " << (*p)->hbond_score() << " : "
                << (*p)->hphobe_score() << " : " << (*p)->vol_score() << endl;
      }
#endif

      last_mol_done[status.MPI_SOURCE] = mols_being_done[status.MPI_SOURCE];
      unsigned int lowest_done = *min_element( last_mol_done.begin() + 1 , last_mol_done.end() );
      report_progress( num_mols_done , ts.chunk_size() , restart_point ,
                       report_step , hit_list , cout );
      dump_results_so_far( num_mols_done , lowest_done , report_step ,
                           query_mol , query_sites , hit_list , ts );
    } else if( string( "Send_Mol_Num" ) == msg ) {
      DACLIB::mpi_send_string( string( "Next_Mol_Num" ) , status.MPI_SOURCE );
      MPI_Send( &curr_mol , 1 , MPI_UNSIGNED , status.MPI_SOURCE , 0 , MPI_COMM_WORLD );
#ifdef NOTYET
      cout << "Sent " << curr_mol << " to " << status.MPI_SOURCE << endl;
#endif
      mols_being_done[status.MPI_SOURCE] = curr_mol;
      curr_mol += ts.chunk_size();
    } else if( string( "Search_Finished" ) == msg ) {
#ifdef NOTYET
      cout << "Search_Finished from " << status.MPI_SOURCE << endl;
#endif
      ++num_finished;
      cout << "Number of processes now finished : " << num_finished << " out of " << num_slaves << endl;

    }

  }

}

// ********************************************************************
void parallel_triphic_search( TriphicSettings &ts , int world_size ) {

  send_search_details_to_slaves( ts , world_size );

  // the input verification and output needs the query mols and sites,
  // so for laziness do a whole search setup. It wastes a bit of time,
  // but similar setups are being done on slaves, so it shouldn't
  // slow things down too much.
  vector<pair<string,string> > input_smarts , smarts_sub_defn;
  PharmPoint pharm_points;
  OEMolBase *query_mol;
  vector<BasePPhoreSite *> query_sites;
  vector<SinglePPhoreSite *> query_score_sites;
  shared_ptr<OEMolBase> protein;
  vector<string> mol_subset;
  shared_ptr<DACLIB::VolumeGrid> query_vol_grid;
  shared_ptr<DACLIB::VolumeGrid> protein_grid;
  vector<pair<string,DACLIB::VolumeGrid *> > score_vol_grids;
  shared_ptr<OESzybki> szybki;
  shared_ptr<OEOmega> omega;
  vector<shared_ptr<OESubSearch> > not_smarts_subs;
  setup_search( ts , input_smarts , smarts_sub_defn , pharm_points ,
                query_mol , query_sites , query_score_sites , protein ,
                mol_subset , query_vol_grid , protein_grid ,
                score_vol_grids , szybki , omega , not_smarts_subs ,
                true );

  ts.print_required_sites_and_points( cout );
  if( !ts.check_required_sites( query_sites ) ||
      !ts.check_required_points( query_sites ) ) {
    exit( 1 );
  }

  list<OverlayScore *> hit_list;

  unsigned int i = 0;
  if( ts.restart() ) {
    read_restart_hits( ts , i , hit_list );
    if( ts.restart_number() ) {
      i = ts.restart_number(); // over-ride what was in the file
    }
    cout << "Search restarts at molecule " << i << "." << endl;
  }

  // send molecules 1 at at time to slaves until all are done. To accommodate
  // searches of different lengths, and so maximise throughput of all slaves
  send_parallel_searches( world_size , ts , hit_list ,
                          query_mol , query_sites , i );

  // bung the results out
  try {
    output_results( query_mol , query_sites , hit_list , ts , false );
  } catch( DACLIB::FileWriteOpenError &e ) {
    cerr << e.what() << endl;
  }

  send_finished_messages( world_size );

  cout << "leaving parallel_triphic_search" << endl;

}

// ********************************************************************
void slave_event_loop() {

  string msg;

  vector<pair<string,string> > input_smarts , smarts_sub_defn;
  PharmPoint pharm_points;
  OEMolBase *query_mol;
  vector<BasePPhoreSite *> query_sites;
  vector<SinglePPhoreSite *> query_score_sites;
  shared_ptr<OEMolBase> protein;
  TriphicSettings ts;
  vector<string> mol_subset;
  shared_ptr<DACLIB::VolumeGrid> query_vol_grid;
  shared_ptr<DACLIB::VolumeGrid> protein_grid;
  shared_ptr<OESzybki> szybki;
  shared_ptr<OEOmega> omega;
  vector<shared_ptr<OESubSearch> > not_smarts_subs;
  list<OverlayScore *> hit_list;
  vector<pair<string,DACLIB::VolumeGrid *> > score_vol_grids;

  while( 1 ) {

    DACLIB::mpi_rec_string( 0 , msg );

#if DEBUG == 1
    int world_rank;
    MPI_Comm_rank( MPI_COMM_WORLD , &world_rank );
    cout << world_rank << " received message : " << msg << endl;
#endif

    if( string( "Finished" ) == msg ) {
      break;
    } else if( string( "Set_Environment" ) == msg ) {
      DACLIB::set_environment_var_from_mpi();
    } else if( string( "Search_Details" ) == msg ) {
      ts.receive_contents_via_mpi();
      DACLIB::mpi_send_string( string ( "Send_Mol_Num" ) , 0 );
      setup_search( ts , input_smarts , smarts_sub_defn , pharm_points ,
                    query_mol , query_sites , query_score_sites , protein ,
                    mol_subset , query_vol_grid , protein_grid ,
                    score_vol_grids , szybki , omega , not_smarts_subs ,
                    false );
    } else if( string( "Next_Mol_Num" ) == msg ) {
      unsigned int next_mol;
      MPI_Recv( &next_mol , 1 , MPI_UNSIGNED , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
      search_database( next_mol , ts.chunk_size() , input_smarts ,
                       smarts_sub_defn , pharm_points , query_mol ,
                       query_sites , query_score_sites , mol_subset , ts ,
                       query_vol_grid ,  protein_grid , score_vol_grids ,
                       szybki , omega , not_smarts_subs , hit_list );
    } else if( string( "New_CWD" ) == msg ) {
      receive_new_cwd();
    }
  }

}

