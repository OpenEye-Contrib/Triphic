//
// file triphic_subs.cc
// David Cosgrove
// AstraZeneca
// 31st July 2007
//
// Some of the triphic-specific functions for new triphic (v3.0).

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <oechem.h>
#include <oeszybki.h>

#include <pvm3.h>

#include <boost/cast.hpp>
#include <boost/foreach.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/construct.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>

#include "stddefs.H"
#include "BasePPhoreSite.H"
#include "DefaultPointsDefs.H"
#include "FileExceptions.H"
#include "OverlayScore.H"
#include "ParseSMARTSXML.H"
#include "PharmPointPVM.H"
#include "SinglePPhoreSite.H"
#include "SMARTSExceptions.H"
#include "TriphicSettings.H"
#include "VolumeGrid.H"
#include "VolumeGridPVM.H"

using namespace boost;
using namespace std;
using namespace OEChem;
using namespace OESz;

namespace DACLIB {
void read_smarts_file( const string &smarts_file ,
                       vector<pair<string,string> > &input_smarts ,
                       vector<pair<string,string> > &smarts_sub_defn );
void make_pphore_sites( OEMol &mol , PharmPoint &pharm_points ,
                        const vector<pair<string,string> > &input_smarts ,
                        vector<pair<string,string> > &smarts_sub_defn ,
                        vector<vector<SinglePPhoreSite *> > &pharm_sites );
VolumeGrid *prepare_mol_grid( OEMolBase *mol );
void send_environment_var_to_pvm_slaves( const string &var_name ,
                                         const string &var_val ,
                                         vector<int> &slave_tids );
void set_environment_var_from_pvm();
bool was_it_a_pvm_failure_message( int bufid , int &dead_tid );

// In pvm_string_subs.cc
// pack a C++ string into a pvm buffer
void pack_string( const string &str );
// unpack a C++ string from pvm buffer
void unpack_string( string &str );
void pack_strings_vector( const vector<string> &strs );
void unpack_strings_vector( vector<string> &strs );

void split_filename( const string &filename , string &file_root ,
                     string &file_ext );
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
                              vector<vector<SinglePPhoreSite *> > &query_score_sites ,
                              vector<SinglePPhoreSite *> &target_sites ,
                              OEMol *target_mol , int target_conf_num ,
                              shared_ptr<DACLIB::VolumeGrid> &query_solid_grid ,
                              shared_ptr<DACLIB::VolumeGrid> &protein_grid ,
                              const vector<pair<string,DACLIB::VolumeGrid *> > &score_vol_grids ,
                              const vector<vector<int> > &cliques ,
                              GtplDefs::SCORE_METHOD score_method ,
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
                              unsigned int max_num_hits ,
                              shared_ptr<OESzybki> &szybki ,
                              list<OverlayScore *> &clique_overlays );
// this also in score_and_store_cliques.cc
bool add_overlay_score_to_list( GtplDefs::SCORE_METHOD score_method ,
                                unsigned int max_num_hits ,
                                OverlayScore *new_score ,
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
void send_openeye_license_to_slaves( const vector<int> &slave_tids );
void send_cwd_to_slaves( const vector<int> &slave_tids );
void send_finished_messages( const vector<int> &slave_tids );
void send_database_details_to_slaves( const vector<int> &slave_tids );
void receive_new_cwd();
void receive_database_details( int &db_start , int &db_step );
void send_progress_to_master( const string &progress_report );

// in step_oemolstream due to laziness
void open_databasefile( const string &db_file ,
                        bool single_conf_mols ,
                        oemolistream *&ims );

void read_moe_ph4_file( TriphicSettings &ts ,
                        PharmPoint &pharm_points ,
                        vector<pair<string,string> > &input_smarts ,
                        vector<pair<string,string> > &smarts_sub_defn ,
                        vector<vector<BasePPhoreSite *> > &query_sites ,
                        vector<pair<string,DACLIB::VolumeGrid *> > &score_vol_grids );

// Message that a PVM process has died as suggested in the PVM manual
const int TASK_DIED = 11;

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
void create_query_sites( vector<pair<string,string> > &input_smarts ,
                         vector<pair<string,string> > &smarts_sub_defn ,
                         PharmPoint &pharm_points ,
                         vector<OEMolBase *> &query_mols ,
                         vector<vector<BasePPhoreSite *> > &query_sites ,
                         vector<vector<SinglePPhoreSite *> > &query_score_sites ) {

  for( int i = 0 , is = query_mols.size() ; i < is ; ++i ) {
    vector<vector<SinglePPhoreSite *> > next_sites;
    OEMol next_mol( *query_mols[i] );
    try {
      DACLIB::make_pphore_sites( next_mol , pharm_points , input_smarts ,
                                 smarts_sub_defn , next_sites );
    } catch( DACLIB::SMARTSDefnError &e ) {
      cerr << e.what() << endl;
      exit( 1 );
    }
    // we know that there's only 1 conformation per molecule
    query_sites.push_back( vector<BasePPhoreSite*>( next_sites.front().begin() ,
                                                    next_sites.front().end() ) );
    query_score_sites.push_back( next_sites.front() );
  }

}

// ***************************************************************************
// delete from query_sites any sites named in ignore_sites
void apply_ignore_sites( const vector<string> &ignore_sites ,
                         vector<vector<BasePPhoreSite *> > &query_sites ) {

  for( int i = 0 , is = ignore_sites.size() ; i < is ; ++i ) {
    bool found_it = false;
    for( int j = 0 , js = query_sites.size() ; j < js ; ++j ) {
      for( int k = 0 , ks = query_sites[j].size() ; k < ks ; ++k ) {
        if( ignore_sites[i] == query_sites[j][k]->get_full_name() ) {
          delete query_sites[j][k];
          query_sites[j][k] = 0;
          found_it = true;
        }
      }
      if( found_it ) {
        query_sites[j].erase( remove( query_sites[j].begin() ,
                                      query_sites[j].end() ,
                                      static_cast<BasePPhoreSite *>( 0 ) ) );
        break;
      }
    }
    if( !found_it ) {
      cout << "Warning: ignore site " << ignore_sites[i] << " not found in query."
           << endl;
    }
  }

}

// ***************************************************************************
void read_mol_file_into_vector( const string &mol_file ,
                                vector<OEMolBase *> &mols ) {

  oemolistream ims( mol_file );
  if( !ims )
    throw DACLIB::FileReadOpenError( mol_file.c_str() );
  
  // by default, we'll get single-conformer molecules, which is what we
  // want for the query
  OEMol mol;
  while( ims >> mol ) {
    OEAssignAromaticFlags( mol , OEAroModelDaylight );
    mols.push_back( OENewMolBase( mol , OEMolBaseType::OEDefault ) );
  }

  cout << "Read " << mols.size() << " molecule";
  if( mols.size() != 1 )
    cout << "s." << endl;
  else
    cout << "." << endl;

}

// ***************************************************************************
void read_query_file( TriphicSettings &ts ,
                      vector<pair<string,string> > &input_smarts ,
                      vector<pair<string,string> > &smarts_sub_defn ,
                      PharmPoint &pharm_points ,
                      vector<OEMolBase *> &query_mols ,
                      vector<vector<BasePPhoreSite *> > &query_sites ,
                      vector<vector<SinglePPhoreSite *> > &query_score_sites ,
                      vector<pair<string,DACLIB::VolumeGrid *> > &score_vol_grids ) {

  string query_root , query_ext;
  DACLIB::split_filename( ts.query_file() , query_root , query_ext );
  if( query_ext == "ph4" ) {
    read_moe_ph4_file( ts , pharm_points , input_smarts , smarts_sub_defn ,
		       query_sites , score_vol_grids );
    if( query_sites.empty() ) {
      cout << "No sites read." << endl;
      exit( 0 );
    } else {
      cout << "Read " << query_sites.back().size() << " sites for query." << endl;
    }
    if( ts.sites_score_file().empty() ) {
      // make some empty molecules
      for( int i = 0 , is = query_sites.size() ; i < is ; ++i ) {
        query_mols.push_back( OENewMolBase( OEMolBaseType::OEDefault ) );
      }
    } else {
      // read some structures to be used in scoring
      read_mol_file_into_vector( ts.sites_score_file() , query_mols );
      if( query_mols.size() != query_sites.size() ) {
        cerr << "Bad error - a different number of sets of query sites ("
             << query_sites.size() << ") than molecules for scoring ("
             << query_mols.size() << "). Can't continue." << endl;
      }
    }
  } else {
    read_mol_file_into_vector( ts.query_file() , query_mols );
    create_query_sites( input_smarts , smarts_sub_defn , pharm_points ,
                        query_mols , query_sites , query_score_sites );
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
void prepare_query_volume_grids( vector<OEMolBase *> &query_mols ,
                                 vector<shared_ptr<DACLIB::VolumeGrid> > &query_grids ,
                                 shared_ptr<OEMolBase> &protein ,
                                 shared_ptr<DACLIB::VolumeGrid> &protein_grid ) {

  // prepare_mol_grid puts surface and core markers into the grid, so can
  // do single conf against single conf volume and surface volume calcs in
  // one go. Even if we're not sorting or filtering using GRID_SHAPE_TANI,
  // we'll want to calculate the final number of the grid_shape_tani,
  // included and total vols.
  for( int i = 0 , is = query_mols.size() ; i < is ; ++i ) {
    if( *query_mols[i] ) {
      query_grids.push_back( shared_ptr<DACLIB::VolumeGrid>( DACLIB::prepare_mol_grid( query_mols[i] ) ) );
    } else {
      query_grids.push_back( shared_ptr<DACLIB::VolumeGrid>( static_cast<DACLIB::VolumeGrid *>( 0 ) ) );
    }
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
                            vector<vector<SinglePPhoreSite *> > &query_score_sites ,
                            vector<vector<SinglePPhoreSite *> > &target_sites ,
                            shared_ptr<DACLIB::VolumeGrid> &query_grid ,
                            shared_ptr<DACLIB::VolumeGrid> protein_grid ,
                            const vector<pair<string,DACLIB::VolumeGrid *> > &score_vol_grids ,
                            bool use_ring_norms , bool use_h_vectors ,
                            bool use_lps ,
                            OEMolBase *query_mol , bool gauss_shape_tani ,
                            float max_rms ,
                            shared_ptr<OESzybki> &szybki ,
                            OverlayScore *hit_score ) {

  // final overlaid results for hit
  OEMolBase *hit_conf = 0;

  overlay_mols_and_sites( &target_mol , query_sites ,
                          target_sites[hit_score->get_moving_conf()] ,
			  hit_score , use_ring_norms , use_h_vectors , use_lps ,
			  hit_conf );
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
  hit_score->calc_robins_scores( query_score_sites ,
                                 target_sites[hit_score->get_moving_conf()] ,
				 query_grid );

  add_overlay_scores_to_mol( *hit_conf , *hit_score );

  delete hit_conf;

}

// ********************************************************************
string remove_conf_num_from_mol_name( const string &mol_name ) {

  string ret( mol_name );
  size_t us_pos = ret.rfind( "_" );
  if( string::npos == us_pos ) {
    return ret;
  } else {
    return ret.substr( 0 , us_pos );
  }

}

// ********************************************************************
void write_output_to_mol_and_score_files( vector<OEMolBase *> &query_mols ,
                                          vector<vector<BasePPhoreSite *> > &query_sites ,
                                          vector<list<OverlayScore *> > &hits ,
                                          const string &output_file_root ,
                                          const string &output_file_ext ,
                                          const TriphicSettings &ts ) {

  for( unsigned int i = 0 ; i < hits.size() ; ++i ) {
    string this_root = output_file_root;
    if( hits.size() > 1 ) {
      this_root += string( "_" ) + query_mols[i]->GetTitle();
    }
    string output_molfile = this_root + string( "." ) + output_file_ext;

    oemolostream oms;
    if( !ts.output_scores_only() ) {
      oms.open( output_molfile.c_str() );
      if( !oms )
        throw DACLIB::FileWriteOpenError( output_molfile.c_str() );
    }

    string output_scoresfile = this_root + string( ".scores" );
    ofstream ofs( output_scoresfile.c_str() );
    if( !ofs.good() ) {
      throw DACLIB::FileWriteOpenError( output_scoresfile.c_str() );
    }
    string sep = ts.comma_output() ? "," : " ";

    // headers to file, including matching points
    if( !hits[i].empty() ) {
      hits[i].front()->write_scores_headers_to_stream( ofs , sep );
    } else {
      // put some dummy ones in, so as not to produce an empty file.
      OverlayScore ovs( query_mols[i]->GetTitle() , "Dummy" );
      ovs.write_scores_headers_to_stream( ofs , sep );
    }
    for( int j = 0 , js = query_sites[i].size() ; j < js ; ++j ) {
      ofs << sep << query_sites[i][j]->label();
    }
    ofs << sep << "Point_Names " << endl;

    if( ts.output_query_to_hits() ) {
      oms << *query_mols[i];
    }

    list<OverlayScore *>::iterator p , ps;
    for( p = hits[i].begin() , ps = hits[i].end() ; p != ps ; ++p ) {
      if( !ts.output_scores_only() ) {
        OEMolBase *hit_mol = (*p)->get_ov_conf();
        if( ts.no_hit_conf_number() ) {
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
      vector<int> ms( query_sites[i].size() , 0 );
      for( int j = 0 , js = (*p)->num_sites() ; j < js ; ++j ) {
        ms[(*p)->sites()[2*j]] = 1;
      }
      for( int j = 0 , js = ms.size() ; j < js ; ++j ) {
        ofs << sep << ms[j];
      }
      ofs << sep;

      // finally, : separated list of names of query sites that matched
      for( int j = 0 , js = (*p)->num_sites() ; j < js ; ++j ) {
        ofs << query_sites[i][(*p)->sites()[2*j]]->label();
        if( j < js - 1 ) {
          ofs << ":";
        }
      }
      ofs << endl;
    }
    // cout << "finished writing hit " << i << " of " << hits.size() << endl;

  }

}

// ********************************************************************
void output_results( vector<OEMolBase *> &query_mols ,
                     vector<vector<BasePPhoreSite *> > &query_sites ,
                     vector<list<OverlayScore *> > &hits ,
                     TriphicSettings &ts ) {

  for( int i = 0 , is = hits.size() ; i < is ; ++i ) {
    if( 1 == hits[i].size() )
      cout << "There was 1 hit for search " << i + 1 << "."
           << endl;
    else
      cout << "There were " << hits[i].size() << " hits for search " << i + 1
           << "." << endl;
  }

  string output_file_root , output_file_ext;
  DACLIB::split_filename( ts.output_file() , output_file_root , output_file_ext );

  write_output_to_mol_and_score_files( query_mols , query_sites , hits ,
				       output_file_root , output_file_ext , ts );

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
                 vector<vector<SinglePPhoreSite *> > &query_score_sites ,
                 vector<vector<SinglePPhoreSite *> > &target_sites ,
                 shared_ptr<DACLIB::VolumeGrid> &query_vol_grid ,
                 shared_ptr<DACLIB::VolumeGrid> &protein_grid ,
                 const vector<pair<string,DACLIB::VolumeGrid *> > &score_vol_grids ,
                 OEMolBase *query_mol , OEMol &target_mol ,
                 shared_ptr<OESzybki> &szybki ,
                 TriphicSettings &ts ,
                 list<OverlayScore *> &full_hit_list ) {

  if( this_query_hits.empty() ) {
    return;
  }

#ifdef NOTYET
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
                           ts.max_rms() , szybki , this_score );
    add_overlay_score_to_list( ts.score_method() , ts.max_hits() ,
                               this_score , full_hit_list );
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
            lambda::bind( lambda::delete_ptr() , lambda::_1 ) );

}

// **************************************************************************
vector<char> do_queries_have_ext_points( const vector<vector<BasePPhoreSite *> > &query_sites ) {

  // MOE pharmacophores can have extension points defined without the
  // associated site. If the extension point has a site in close proximity,
  // it will have been merged as a direction into the base site, but if
  // that wasn't possible, we need to add them to the database molecule
  // as well. The sites will be Acc2 or Don2 if this is the case.
  vector<char> ret_val = vector<char>( query_sites.size() , 0 );

  for( int i = 0 , is = query_sites.size() ; i < is ; ++i ) {
    for( int j = 0 , js = query_sites[i].size() ; j < js ; ++j ) {
      if( string( "Acc2" ) == query_sites[i][j]->get_type_string() ||
          string( "Don2" ) == query_sites[i][j]->get_type_string() ) {
        ret_val[i] = 1;
        break;
      }
    }
  }

  return ret_val;

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
                             vector<OEMolBase *> &query_mols ,
                             vector<vector<BasePPhoreSite *> > &query_sites ,
                             vector<vector<SinglePPhoreSite *> > &query_score_sites ,
                             shared_ptr<OEMolBase> &protein ,
                             const vector<string> &mol_subset ,
                             TriphicSettings &ts ,
                             vector<shared_ptr<DACLIB::VolumeGrid> > &query_vol_grids ,
                             shared_ptr<DACLIB::VolumeGrid> &protein_grid ,
                             const vector<pair<string,DACLIB::VolumeGrid *> > &score_vol_grids ,
                             shared_ptr<OESzybki> &szybki ,
                             vector<shared_ptr<OESubSearch> > &not_smarts_subs ,
                             OEMol &target_mol ,
                             vector<list<OverlayScore *> > &hit_lists ) {

#ifdef NOTYET
  cout << "Searching " << target_mol.GetTitle() << endl;
#endif

  if( !mol_subset.empty() &&
      !binary_search( mol_subset.begin() , mol_subset.end() , target_mol.GetTitle() ) ) {
    return;
  }

  static vector<char> need_ext_points;
  if( need_ext_points.empty() ) {
    need_ext_points = do_queries_have_ext_points( query_sites );
  }
  
  // make sure this molecule doesn't hit something in not_smarts_subs - if it
  // does, we don't want it.
  if( not_smarts_hit( not_smarts_subs , target_mol ) )
    return;

  vector<vector<SinglePPhoreSite *> > target_sites;
  vector<vector<SinglePPhoreSite *> > target_ext_sites;
  try {
    DACLIB::make_pphore_sites( target_mol , pharm_points , input_smarts ,
                               smarts_sub_defn , target_sites );
    if( count( need_ext_points.begin() , need_ext_points.end() , 1 ) ) {
      make_pphore_ext_sites( pharm_points , target_sites , target_ext_sites );
    }
  } catch( DACLIB::SMARTSDefnError &e ) {
    cerr << e.what() << endl;
    exit( 1 );
  }

  vector<vector<int> > cliques;
  vector<vector<SinglePPhoreSite *> > full_target_sites;
  for( int i = 0 , is = query_sites.size() ; i < is ; ++i ) {

    list<OverlayScore *> clique_overlays;
    for( int j = 0 , js = target_sites.size() ; j < js ; ++j ) {

      full_target_sites.push_back( target_sites[j] );
      if( need_ext_points[i] ) {
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

      find_cliques( query_sites[i] , full_target_sites[j] , ts.dist_tol() ,
                    ts.scaled_dist_tol() , ts.dont_do_sub_cliques() ,
                    ts.min_clique_size() , cliques );
      score_and_store_cliques( query_mols[i]->GetTitle() , 0 ,
                               query_mols[i] , query_sites[i] ,
                               query_score_sites , full_target_sites[j] ,
                               &target_mol , j , query_vol_grids[i] ,
                               protein_grid , score_vol_grids ,
                               cliques , ts.score_method() ,
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
                               ts.max_hits() , szybki , clique_overlays );
    }

    store_hits( clique_overlays , query_sites[i] , query_score_sites ,
                full_target_sites , query_vol_grids[i] , protein_grid ,
                score_vol_grids , query_mols[i] ,	target_mol , szybki , ts ,
                hit_lists[i] );

  }

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
void report_progress( unsigned int &i , unsigned int &report_step ,
                      ostream &os ) {

  ++i;
  if( i && !(i % report_step) ) {
    os << "Done " << i << " molecules so far." << endl;
  }

  if( i >= report_step * 10 ) {
    report_step *= 10;
  }

}

// **************************************************************************
void report_progress( unsigned int &i , unsigned int &report_step ,
                      OEMolBase &target_mol ,
                      vector<list<OverlayScore *> > &hit_lists ,
                      ostream &os ) {

  ++i;
  if( i && !(i % report_step) ) {
    os << "Done " << i << " molecules (to " << target_mol.GetTitle()
         << ") so far with ";
    for( int j = 0 , js = hit_lists.size() ; j < js ; ++j ) {
      os << hit_lists[j].size() << " hit";
      if( hit_lists[j].size() != 1 ) {
        os << "s";
      }
      if( !hit_lists[j].empty() ) {
        os << " against " << hit_lists[j].front()->get_fixed_mol_name();
      }
      if( j < int( hit_lists.size() ) - 1 && hit_lists.size() > 1 ) {
        os << " and ";
      } else {
        os << "." << endl;
      }
    }
  }
  if( i >= report_step * 10 ) {
    report_step *= 10;
  }

}

// ********************************************************************
void send_results_to_master( int master_tid , int my_tid ,
                             vector<list<OverlayScore *> > &hits ) {

  pvm_initsend( PvmDataDefault );
  pvm_pkstr( const_cast<char *>( string( "Results" ).c_str()) );
  pvm_pkint( &my_tid , 1 , 1 );
  int num_to_send = hits.size();
  pvm_pkint( &num_to_send , 1 , 1 );

  for( int i = 0 ; i < num_to_send ; ++i ) {

    int j = hits[i].size();
    pvm_pkint( &j , 1 , 1 );
    list<OverlayScore *>::iterator p , ps;
    for( p = hits[i].begin() , ps = hits[i].end() ; p != ps ; ++p ) {
      DACLIB::pack_string( (*p)->write_contents_to_string() );
    }
  }

  pvm_send( master_tid , 0 );

}

// ********************************************************************
void send_all_results_to_master( int master_tid , int my_tid ,
                                 vector<list<OverlayScore *> > &hits ) {

#ifdef NOTYET
  vector<list<OverlayScore *> > tmp_hits;
  convert_hit_lists_for_output( hits , tmp_hits );
#endif

  send_results_to_master( master_tid , my_tid , hits );

}

// ********************************************************************
// send those results to master that are better than the one the master
// has sent, which must be unpacked.
void send_subset_results_to_master( int master_tid , int my_tid ,
                                    GtplDefs::SCORE_METHOD score_method ,
                                    vector<list<OverlayScore *> > &hits ) {

  vector<OverlayScore *> scores_to_beat;
  for( int i = 0 , is = hits.size() ; i < is ; ++i ) {
    string to_beat_msg;
    DACLIB::unpack_string( to_beat_msg );
    if( string( "SEND_ALL_HITS" ) == to_beat_msg ) {
      scores_to_beat.push_back( new OverlayScore() );
    } else {
      scores_to_beat.push_back( new OverlayScore( to_beat_msg ) );
    }
  }

  vector<list<OverlayScore *> > hits_to_send;
  for( int i = 0 , is = hits.size() ; i < is ; ++i ) {
    hits_to_send.push_back( list<OverlayScore *>() );
    list<OverlayScore *>::iterator p , ps;
    for( p = hits[i].begin() , ps = hits[i].end() ; p != ps ; ++p ) {
      if( *scores_to_beat[i] < *(*p) )
        hits_to_send[i].push_back( *p );
    }
  }

#ifdef NOTYET
  vector<list<TriphicOverlayScore *> > tmp_hits;
  convert_hit_lists_for_output( hits_to_send , tmp_hits );
#endif

  send_results_to_master( master_tid , my_tid , hits_to_send );

}

// **************************************************************************
void tell_master_slave_is_finished( int master_tid , int my_tid ) {

  pvm_initsend( PvmDataDefault );
  pvm_pkstr( const_cast<char *>( string( "Search_Finished" ).c_str() ) );
  pvm_pkint( &my_tid , 1 , 1 );
  pvm_send( master_tid , 0 );

}

// **************************************************************************
// read the molecule of given sequence number from the set of database files
// returns empty shared_ptr if unsuccessful
shared_ptr<OEMol> read_nth_mol_from_oemolistream( unsigned int next_mol ,
                                                  TriphicSettings &ts ) {

  static oemolistream *db_ims = 0;
  static unsigned int curr_db_file = 0;
  static unsigned int curr_mol = 0;

  if( !db_ims ) {
    db_ims = new oemolistream;
    cout << "Searching database file " << ts.db_files()[curr_db_file] << endl;
    open_databasefile( ts.db_files()[curr_db_file] ,
                       ts.single_conf_mols() , db_ims );
  }

  static OEMol mol;
  while( curr_mol < next_mol ) {
    if( !(*db_ims >> mol) ) {
      db_ims->close();
      ++curr_db_file;
      if( curr_db_file == ts.db_files().size() ) {
        // we're done
        return shared_ptr<OEMol>();
      }
      cout << "Moving search to database file " << ts.db_files()[curr_db_file] << endl;
      open_databasefile( ts.db_files()[curr_db_file] ,
                         ts.single_conf_mols() , db_ims );
    }
    ++curr_mol;
  }

  return shared_ptr<OEMol>( new OEMol( mol ) );

}

// **************************************************************************
// search the next_mol'th molecule in the databases
void search_database( unsigned int next_mol , int master_tid , int my_tid ,
                      vector<pair<string,string> > &input_smarts ,
                      vector<pair<string,string> > &smarts_sub_defn ,
                      PharmPoint &pharm_points ,
                      vector<OEMolBase *> &query_mols ,
                      vector<vector<BasePPhoreSite *> > &query_sites ,
                      vector<vector<SinglePPhoreSite *> > &query_score_sites ,
                      shared_ptr<OEMolBase> &protein ,
                      const vector<string> &mol_subset ,
                      TriphicSettings &ts ,
                      vector<shared_ptr<DACLIB::VolumeGrid> > &query_vol_grids ,
                      shared_ptr<DACLIB::VolumeGrid> &protein_grid ,
                      const vector<pair<string,DACLIB::VolumeGrid *> > &score_vol_grids ,
                      shared_ptr<OESzybki> &szybki ,
                      vector<shared_ptr<OESubSearch> > &not_smarts_subs ,
                      vector<list<OverlayScore *> > &hit_lists ) {

  unsigned int report_step = 10;
  static unsigned int num_searched = 0;

  shared_ptr<OEMol> tmol = read_nth_mol_from_oemolistream( next_mol , ts );

  if( !tmol ) {

    // we're done
    tell_master_slave_is_finished( master_tid , my_tid );

  } else {

#ifdef NOTYET
    cout << "Searching " << tmol->GetTitle() << " : " << tmol->NumConfs() << endl;
#endif
    OEAssignAromaticFlags( *tmol , OEAroModelDaylight );
    search_target_molecule( input_smarts , smarts_sub_defn , pharm_points ,
                            query_mols , query_sites , query_score_sites ,
                            protein , mol_subset , ts , query_vol_grids ,
                            protein_grid , score_vol_grids , szybki ,
                            not_smarts_subs , *tmol , hit_lists );
    ostringstream oss;
    report_progress( num_searched , report_step , *tmol , hit_lists , oss );
    if( !oss.str().empty() ) {
      send_progress_to_master( oss.str() );
    }

    // ask for next molecule
    pvm_initsend( PvmDataDefault );
    pvm_pkstr( const_cast<char *>( string( "Send_Mol_Num" ).c_str() ) );
    pvm_pkint( &my_tid , 1 , 1 );
    pvm_send( master_tid , 0 );

  }

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

// ********************************************************************
void setup_search( TriphicSettings &ts ,
                   vector<pair<string,string> > &input_smarts ,
                   vector<pair<string,string> > &smarts_sub_defn ,
                   PharmPointPVM &pharm_points ,
                   vector<OEMolBase *> &query_mols ,
                   vector<vector<BasePPhoreSite *> > &query_sites ,
                   vector<vector<SinglePPhoreSite *> > &query_score_sites ,
                   shared_ptr<OEMolBase> &protein ,
                   vector<string> &mol_subset ,
                   vector<shared_ptr<DACLIB::VolumeGrid> > &query_vol_grids ,
                   shared_ptr<DACLIB::VolumeGrid> &protein_grid ,
                   vector<pair<string,DACLIB::VolumeGrid *> > &score_vol_grids ,
                   shared_ptr<OESzybki> &szybki ,
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

  // read the query file and, if necessary, create the sites and volumes
  DACLIB::VolumeGrid::set_grid_spacing( ts.vol_grid_spacing() );
  read_query_file( ts , input_smarts , smarts_sub_defn , pharm_points ,
                   query_mols , query_sites , query_score_sites ,
                   score_vol_grids );

  if( report_sites ) {
    for( int i = 0 , is = query_sites.size() ; i < is ; ++i ) {
      cout << endl << "Query " << i + 1 << endl;
      brief_report_sites( cout , query_sites[i] );
      cout << endl;
    }
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

  prepare_query_volume_grids( query_mols , query_vol_grids ,
                              protein , protein_grid );

  if( OESzybkiIsLicensed() && protein &&
      ts.max_mmff_nrg() < numeric_limits<float>::max() ) {
    prepare_szybki_optimiser( ts , protein , szybki );
  }

  // SMARTS that the hits can't match
  make_not_smarts_queries( ts.not_smarts_list() , smarts_sub_defn ,
                           not_smarts_subs );

  OverlayScore::set_score_method( ts.score_method() );

}

// ********************************************************************
void serial_triphic_search( TriphicSettings &ts ) {

  vector<pair<string,string> > input_smarts , smarts_sub_defn;
  PharmPointPVM pharm_points;
  vector<OEMolBase *> query_mols;
  vector<vector<BasePPhoreSite *> > query_sites;
  vector<vector<SinglePPhoreSite *> > query_score_sites;
  shared_ptr<OEMolBase> protein;
  vector<string> mol_subset;
  vector<shared_ptr<DACLIB::VolumeGrid> > query_vol_grids;
  shared_ptr<DACLIB::VolumeGrid> protein_grid;
  vector<pair<string,DACLIB::VolumeGrid *> > score_vol_grids;
  shared_ptr<OESzybki> szybki;
  vector<shared_ptr<OESubSearch> > not_smarts_subs;
  oemolistream *ims;

  setup_search( ts , input_smarts , smarts_sub_defn , pharm_points ,
                query_mols , query_sites , query_score_sites , protein ,
                mol_subset , query_vol_grids , protein_grid , score_vol_grids ,
                szybki , not_smarts_subs , true );

  if( ts.print_sites_and_stop() ) {
    exit( 0 );
  }

  ts.print_required_sites_and_points( cout );
  if( !ts.check_required_sites( query_sites ) ||
      !ts.check_required_points( query_sites ) ) {
    exit( 1 );
  }

  OEMol target_mol;
  vector<list<OverlayScore *> > hit_lists( query_sites.size() ,
                                           list<OverlayScore *>() );

  vector<string> db_files = ts.db_files();
  BOOST_FOREACH( string db_file , db_files ) {

    cout << "Searching database file : " << db_file << endl;

    unsigned int i = 0;
    unsigned int report_step = 1;
    ims = new oemolistream;
    open_databasefile( db_file , ts.single_conf_mols() , ims );

    while( *ims >> target_mol ) {
#ifdef NOTYET
      cout << "Searching " << target_mol.GetTitle() << " : " << target_mol.NumConfs() << endl;
#endif
      OEAssignAromaticFlags( target_mol , OEAroModelDaylight );
      search_target_molecule( input_smarts , smarts_sub_defn , pharm_points ,
                              query_mols , query_sites , query_score_sites ,
                              protein , mol_subset , ts , query_vol_grids ,
                              protein_grid , score_vol_grids , szybki ,
                              not_smarts_subs , target_mol , hit_lists );

      report_progress( i , report_step , target_mol , hit_lists , cout );
    }
    ims->close();

  }

  try {
    output_results( query_mols , query_sites , hit_lists , ts );
  } catch( DACLIB::FileWriteOpenError &e ) {
    cerr << e.what() << endl;
  }

  // tidy stuff up, mostly to stop valgrind whinging. triphic dates from before
  // I discovered the joys of boost smart pointers.
  for( int i = 0 , is = hit_lists.size() ; i < is ; ++i ) {
    for( list<OverlayScore *>::iterator p = hit_lists[i].begin() , ps = hit_lists[i].end() ;
         p != ps ; ++p ) {
      delete *p;
    }
  }

}

// ********************************************************************
// PVM-specific bits
// ********************************************************************
// ********************************************************************
void send_search_details_to_slaves( TriphicSettings &ts ,
                                   vector<int> &slave_tids ) {

  send_openeye_license_to_slaves( slave_tids );
  send_cwd_to_slaves( slave_tids );

  pvm_initsend( PvmDataDefault );
  pvm_pkstr( const_cast<char *>( string( "Search_Details" ).c_str() ) );
  ts.pack_contents_into_pvm_buffer();

  pvm_mcast( &slave_tids[0] , slave_tids.size() , 0 );

}

// ********************************************************************
void receive_results_from_slave( int slave_tid , vector<list<OverlayScore *> > &slave_hit_list ) {


  int num_to_rec;
  pvm_upkint( &num_to_rec , 1 , 1 );

  slave_hit_list = vector<list<OverlayScore *> >( num_to_rec ,
                                                  list<OverlayScore *>() );
  for( int i = 0 ; i < num_to_rec ; ++i ) {
    int num_hits_to_rec;
    pvm_upkint( &num_hits_to_rec , 1 , 1 );
    for( int j = 0 ; j < num_hits_to_rec ; ++j ) {
      string hit_msg;
      DACLIB::unpack_string( hit_msg );
      slave_hit_list[i].push_back( new OverlayScore( hit_msg ) );
#ifdef NOTYET
      slave_hit_list[i].back()->write_scores_to_stream( cout , " , " );
      cout << endl;
#endif
    }
  }

}

// ********************************************************************
void merge_hit_lists( unsigned int max_hits ,
                      GtplDefs::SCORE_METHOD score_method ,
                      vector<list<OverlayScore *> > &sub_hit_lists ,
                      vector<list<OverlayScore *> > &full_hit_lists ) {

  // assume that sub_hit_lists and full_hit_lists have the same number of lists
  for( int i = 0 , is = sub_hit_lists.size() ; i < is ; ++i ) {
#ifdef NOTYET
    cout << "Merging " << sub_hit_lists[i].size() << " hits" << endl;
#endif
    list<OverlayScore *>::iterator p , ps;
    for( p = sub_hit_lists[i].begin() , ps = sub_hit_lists[i].end() ; p != ps ; ++p ) {
      add_overlay_score_to_list( score_method , max_hits ,
                                 *p , full_hit_lists[i] );
      // add_overlay_score_to_list no longer maintains the list size.
      if( max_hits && full_hit_lists[i].size() > max_hits ) {
        delete full_hit_lists[i].back();
        full_hit_lists[i].pop_back();
      }
    }
#ifdef NOTYET
    cout << "Number of hits now " << full_hit_lists[i].size() << endl;
#endif
  }

}

// ********************************************************************
// all slaves should have finished by now
void get_results_from_slaves( vector<int> &slave_tids ,
                              unsigned int max_hits ,
                              GtplDefs::SCORE_METHOD score_method ,
                              vector<list<OverlayScore *> > &hit_lists ) {

  while( !slave_tids.empty() ) {

    pvm_initsend( PvmDataDefault );
    pvm_pkstr( const_cast<char *>( string( "Send_Results" ).c_str() ) );
    pvm_send( slave_tids.front() , 0 );

    int bufid = pvm_recv( slave_tids.front() , -1 );
    int done_tid;
    if( DACLIB::was_it_a_pvm_failure_message( bufid , done_tid ) ) {

      cerr << "Process " << done_tid << " has gone belly up, taking all" << endl
           << "its results with it. Carrying on, but there will be missing"
           << endl
           << "hits." << endl;
      cout << "Process " << done_tid << " has gone belly up, taking all" << endl
           << "its results with it. Carrying on, but there will be missing"
           << endl
           << "hits." << endl;
      slave_tids.erase( find( slave_tids.begin() , slave_tids.end() , done_tid ) );
      if( slave_tids.empty() ) {
        cerr << "All slaves are now dead, so that's it for now. This could" << endl
             << "be a problem with the program, it could be a problem with" << endl
             << "your database, or it could be a problem with the machines" << endl
             << "the slaves were running on, including the possibility that" << endl
             << "the OpenEye license isn't set correctly on the slave" << endl
             << "machines." << endl
             << "All results received will be written, but it'll be an incomplete set."
             << endl;
        break; // we can still write what we got.
      }

    } else {

      // if we're here, the message wasn't a failure message.
      char msg[1000]; // it'll be big enough for the message header
      pvm_upkstr( msg );
      cout << "Message : " << msg << " from " << done_tid << endl;
#ifdef NOTYET
#endif
      if( !strcmp( msg , "Results" ) ) {
        int slave_tid;
        pvm_upkint( &slave_tid , 1 , 1 );
        vector<list<OverlayScore *> > slave_hit_lists;
        receive_results_from_slave( slave_tid , slave_hit_lists );
        merge_hit_lists( max_hits , score_method , slave_hit_lists , hit_lists );
        slave_tids.erase( find( slave_tids.begin() , slave_tids.end() , slave_tid ) );
        cout << "Number of slaves still to report : " << slave_tids.size() << endl;
      }

    }
  }

}

// ********************************************************************
void send_parallel_searches( vector<int> &slave_tids ) {

  unsigned int num_finished = 0;
  unsigned int curr_mol = 0;
  unsigned int report_step = 1;

  while( num_finished < slave_tids.size() ) {
    if( slave_tids.empty() ) {
      break;
    }

    int bufid = pvm_recv( -1 , -1 );
    int done_tid = -1;
    if( DACLIB::was_it_a_pvm_failure_message( bufid , done_tid ) ) {

      slave_tids.erase( find( slave_tids.begin() , slave_tids.end() , done_tid ) );
      cerr << "Process " << done_tid << " has gone belly up, taking all" << endl
           << "its results with it. Carrying on, but there will be missing"
           << endl
           << "hits." << endl << endl
           << "Number of slaves still running : " << slave_tids.size() << "." << endl;
      cout << "Process " << done_tid << " has gone belly up, taking all" << endl
           << "its results with it. Carrying on, but there will be missing"
           << endl
           << "hits." << endl
           << "Number of slaves still running : " << slave_tids.size() << "." << endl;
      if( slave_tids.empty() ) {
        cerr << "All slaves are now dead, so that's it for now. This could" << endl
             << "be a problem with the program, it could be a problem with" << endl
             << "your database, or it could be a problem with the machines" << endl
             << "the slaves were running on, including the possibility that" << endl
             << "the OpenEye license isn't set correctly on the slave" << endl
             << "machines." << endl;
        pvm_exit();
        exit( 1 );
      }

    } else {

      // if we're here, the message wasn't a failure message, so it should be
      // a message from the slave that it needs another target mol number
      // or that it's finished.
      char msg[1000]; // it'll be big enough for the message header
      pvm_upkstr( msg );
#ifdef NOTYET
      cout << "Message : " << msg << " from " << done_tid << endl;
#endif
      if( !strcmp( "Send_Mol_Num" , msg ) ) {
        int slave_tid;
        pvm_upkint( &slave_tid , 1 , 1 );
        pvm_initsend( PvmDataDefault );
        pvm_pkstr( const_cast<char *>( string( "Next_Mol_Num" ).c_str() ) );
        pvm_pkuint( &curr_mol , 1 , 1 );
        pvm_send( slave_tid , 0 );
#ifdef NOTYET
        cout << "Sent " << curr_mol << " to " << slave_tid << endl;
#endif
        // report_progress increments curr_mol
        report_progress( curr_mol , report_step , cout );
      } else if( !strcmp( "Progress Report" , msg ) ) {
        pvm_upkint( &done_tid , 1 , 1 );
        string prog_rep;
        DACLIB::unpack_string( prog_rep );
        cout << done_tid << " : " << prog_rep;
      } else if( !strcmp( msg , "Search_Finished" ) ) {
        int slave_tid;
        pvm_upkint( &slave_tid , 1 , 1 );
#ifdef NOTYET
        cout << "Search_Finished from " << slave_tid << endl;
#endif
        ++num_finished;
      }

    }
  }

}

// ********************************************************************
void parallel_triphic_search( TriphicSettings &ts ,
                              vector<int> &slave_tids ) {

  send_search_details_to_slaves( ts , slave_tids );

  // the input verification and output needs the query mols and sites,
  // so for laziness do a whole search setup. It wastes a bit of time,
  // but similar setups are being done on slaves, so it shouldn't
  // slow things down too much.
  vector<pair<string,string> > input_smarts , smarts_sub_defn;
  PharmPointPVM pharm_points;
  vector<OEMolBase *> query_mols;
  vector<vector<BasePPhoreSite *> > query_sites;
  vector<vector<SinglePPhoreSite *> > query_score_sites;
  shared_ptr<OEMolBase> protein;
  vector<string> mol_subset;
  vector<shared_ptr<DACLIB::VolumeGrid> > query_vol_grids;
  shared_ptr<DACLIB::VolumeGrid> protein_grid;
  vector<pair<string,DACLIB::VolumeGrid *> > score_vol_grids;
  shared_ptr<OESzybki> szybki;
  vector<shared_ptr<OESubSearch> > not_smarts_subs;
  setup_search( ts , input_smarts , smarts_sub_defn , pharm_points ,
                query_mols , query_sites , query_score_sites , protein ,
                mol_subset , query_vol_grids , protein_grid ,
                score_vol_grids , szybki , not_smarts_subs , true );

  ts.print_required_sites_and_points( cout );
  if( !ts.check_required_sites( query_sites ) ||
      !ts.check_required_points( query_sites ) ) {
    exit( 1 );
  }

  // send molecules 1 at a time to slaves until all are done. To accommodate
  // searches of different lengths, and so maximise throughput of all slaves.
  send_parallel_searches( slave_tids );

  // get the results from all slaves
  vector<list<OverlayScore *> > hit_lists( query_sites.size() ,
                                           list<OverlayScore *>() );
  // get_results_from_slaves destroys slave_tids
  vector<int> slave_tids_cp( slave_tids );
  get_results_from_slaves( slave_tids , ts.max_hits() ,
                           ts.score_method() , hit_lists );

  send_finished_messages( slave_tids_cp );

  // bung the results out
  try {
    output_results( query_mols , query_sites , hit_lists , ts );
  } catch( DACLIB::FileWriteOpenError &e ) {
    cerr << e.what() << endl;
  }

}

// ********************************************************************
void slave_event_loop() {

  char     msg1[1000];

  vector<pair<string,string> > input_smarts , smarts_sub_defn;
  PharmPointPVM pharm_points;
  vector<OEMolBase *> query_mols;
  vector<vector<BasePPhoreSite *> > query_sites;
  vector<vector<SinglePPhoreSite *> > query_score_sites;
  shared_ptr<OEMolBase> protein;
  TriphicSettings ts;
  vector<string> mol_subset;
  vector<shared_ptr<DACLIB::VolumeGrid> > query_vol_grids;
  shared_ptr<DACLIB::VolumeGrid> protein_grid;
  shared_ptr<OESzybki> szybki;
  vector<shared_ptr<OESubSearch> > not_smarts_subs;
  vector<list<OverlayScore *> > hit_lists;
  vector<pair<string,DACLIB::VolumeGrid *> > score_vol_grids;

  int master_tid = pvm_parent();
  if( PvmNoParent == master_tid ) {
    // It's an error, though probably a deliberate one
    cerr << "ERROR - strange PVM behaviour for triphic." << endl;
    exit( 1 );
  }
  int my_tid = pvm_mytid();

  // we want to know if the master dies
  pvm_notify( PvmTaskExit , TASK_DIED , 1 , &master_tid );

  while( 1 ) {

    // see if there's a message from the parent
    int bufid = pvm_recv( -1 , -1 );
    int dead_tid;
    if( DACLIB::was_it_a_pvm_failure_message( bufid , dead_tid ) &&
        dead_tid == master_tid ) {
      cerr << "AWOOGA : The master is dead.  Following him on into the"
           << " afterlife." << endl;
      pvm_exit();
      exit( 1 );
    }

    pvm_upkstr( msg1 );
#if DEBUG == 1
    cout << my_tid << "  received message : " << msg1 << endl;
#endif

    if( !strcmp( "Finished" , msg1 ) ) {
      break;
    } else if( !strcmp( "Set_Environment" , msg1 ) ) {
      DACLIB::set_environment_var_from_pvm();
    } else if( !strcmp( "Search_Details" , msg1 ) ) {
      ts.unpack_contents_from_pvm_buffer();
      setup_search( ts , input_smarts , smarts_sub_defn , pharm_points ,
                    query_mols , query_sites , query_score_sites , protein ,
                    mol_subset , query_vol_grids , protein_grid ,
                    score_vol_grids , szybki , not_smarts_subs , true );
      hit_lists = vector<list<OverlayScore *> >( query_sites.size() ,
                                                 list<OverlayScore *>() );
      pvm_initsend( PvmDataDefault );
      pvm_pkstr( const_cast<char *>( string( "Send_Mol_Num" ).c_str() ) );
      pvm_pkint( &my_tid , 1 , 1 );
      pvm_send( master_tid , 0 );
    } else if( !strcmp( "Next_Mol_Num" , msg1 ) ) {
      unsigned int next_mol;
      pvm_upkuint( &next_mol , 1 , 1 );
      search_database( next_mol , master_tid , my_tid , input_smarts ,
                       smarts_sub_defn , pharm_points , query_mols ,
                       query_sites , query_score_sites , protein ,
                       mol_subset , ts , query_vol_grids , protein_grid ,
                       score_vol_grids , szybki , not_smarts_subs ,
                       hit_lists );
    } else if( !strcmp( "Send_Results" , msg1 ) ) {
      send_all_results_to_master( master_tid , my_tid , hit_lists );
    } else if( !strcmp( "New_CWD" , msg1 ) ) {
      receive_new_cwd();
    } else if( !strcmp( "Send_Results_Better_Than" , msg1 ) ) {
      send_subset_results_to_master( master_tid , my_tid ,
                                     ts.score_method() , hit_lists );
    }
  }

  pvm_exit();
  exit( 0 );

}

