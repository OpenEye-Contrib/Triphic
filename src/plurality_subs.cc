//
// file plurality_subs.cc
// Dave Cosgrove
// 31st July 2007
//
// Some of the plurality-specific functions for new plurality.

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <oechem.h>
#include <oeszybki.h>

#include <pvm3.h>

#include <boost/bind.hpp>
#include <boost/cast.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/construct.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>

#include "BasePPhoreSite.H"
#include "FileExceptions.H"
#include "PharmPointPVM.H"
#include "PluralityHit.H"
#include "PluralitySettings.H"
#include "PPhoreQuery.H"
#include "SinglePPhoreSite.H"
#include "SMARTSExceptions.H"
#include "VolumeGridPVM.H"

using namespace std;
using namespace OEChem;

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

// in send_messages_to_tp_slaves.cc
void pack_smarts_defs_into_pvm_buffer( vector<pair<string,string> > &input_smarts ,
                                       vector<pair<string,string> > &smarts_sub_defn );
void unpack_smarts_defs_from_pvm_buffer( vector<pair<string,string> > &input_smarts ,
                                         vector<pair<string,string> > &smarts_sub_defn );
void send_openeye_license_to_slaves( const vector<int> &slave_tids );
void send_cwd_to_slaves( const vector<int> &slave_tids );
void send_finished_messages( const vector<int> &slave_tids );
void send_database_details_to_slaves( const vector<int> &slave_tids );
void receive_new_cwd();
void receive_database_details( int &db_start , int &db_step );
void send_progress_to_master( const string &progress_report );

// in eponymous file
unsigned step_oemolstream( oemolistream &ims , int step );
// also in step_oemolstream due to laziness
void open_databasefile( const string &db_file ,
                        bool single_conf_mols ,
                        oemolistream &ims );

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

// ********************************************************************
void read_protein_file( const string &protein_file ,
                        boost::shared_ptr<OEMolBase> &protein ) {

  if( protein_file.empty() )
    return;

  oemolistream ims( protein_file.c_str() );
  if( !ims )
    throw DACLIB::FileReadOpenError( protein_file.c_str() );
  protein =
      boost::shared_ptr<OEMolBase>( OENewMolBase( OEMolBaseType::OEDefault ) );
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
// keep the best (by RMS) distinct hit for the conformations. 2 hits are
// indistinct if they have exactly the same sites (same type, same atoms).
// Need to do a full test on the sites as different conformations may not
// always have the same number of sites, most likely if using ITMOC
// hydrophobes
void keep_best_hits( vector<boost::shared_ptr<PluralityHit> > &hits ) {

  for( int i = 0 , is = hits.size() - 1 ; i < is ; ++i ) {
    if( !hits[i] )
      continue; // already going
    boost::shared_ptr<PluralityHit> best_i = hits[i];
    for( int j = i + 1 ; j < is + 1 ; ++j ) {
      if( !hits[j] )
        continue; // already going
      if( hits[i]->has_same_sites( *hits[j] ) ) {
        if( hits[j]->get_ov_rms() < hits[i]->get_ov_rms() ) {
          best_i = hits[j];
        }
        hits[j].reset(); // flag this for removal
      }
    }
    hits[i] = best_i;
  }

  // take out the empty ones
  hits.erase( remove_if( hits.begin() , hits.end() ,
                         boost::bind( logical_not<bool>() , _1 ) ) ,
              hits.end() );

}

// ********************************************************************
bool not_smarts_hit( vector<boost::shared_ptr<OESubSearch> > &not_smarts_subs ,
                     OEMol &target_mol ) {

  vector<boost::shared_ptr<OESubSearch> >::iterator p , ps;
  for( p = not_smarts_subs.begin() , ps = not_smarts_subs.end() ;
       p != ps ; ++p ) {
    if( (*p)->SingleMatch( target_mol ) )
      return true;
  }

  return false;

}

// ********************************************************************
// returns true if a search was done, regardless of outcome
bool search_target_mol( OEMol &target_mol , PPhoreQuery &query ,
                        vector<pair<string,string> > &input_smarts ,
                        vector<pair<string,string> > &smarts_sub_defn ,
                        PharmPoint &pharm_points ,
                        const vector<string> &mol_subset ,
                        boost::shared_ptr<DACLIB::VolumeGrid> &protein_grid ,
                        boost::shared_ptr<DACLIB::VolumeGrid> &soft_exc_vol_grid ,
                        vector<pair<string,DACLIB::VolumeGrid *> > score_vol_grids ,
                        HITS_TO_OUTPUT hits_to_output ,
                        vector<boost::shared_ptr<OESubSearch> > &not_smarts_subs ,
                        vector<boost::shared_ptr<PluralityHit> > &hits ) {

  if( !mol_subset.empty() &&
      !binary_search( mol_subset.begin() , mol_subset.end() ,
                      target_mol.GetTitle() ) )
    return false;

  // make sure this molecule doesn't hit something in not_smarts_subs - if it
  // does, we don't want it.
  if( not_smarts_hit( not_smarts_subs , target_mol ) )
    return true; // this counts as a search for these purposes

  vector<vector<SinglePPhoreSite *> > target_sites;
  try {
    DACLIB::make_pphore_sites( target_mol , pharm_points , input_smarts ,
                               smarts_sub_defn , target_sites );
  } catch( DACLIB::SMARTSDefnError &e ) {
    cerr << e.what() << endl;
    exit( 1 );
  }

  for( int i = 0 , is = target_sites.size() ; i < is ; ++i ) {
    query.search_pphore_sites( target_sites[i] , target_mol , i ,
                               protein_grid , soft_exc_vol_grid ,
                               score_vol_grids ,
                               hits_to_output == ONE_HIT_ONLY , hits );
    if( hits_to_output == ONE_HIT_ONLY && !hits.empty() )
      break;
  }

  for( int i = 0 , is = target_sites.size() ; i < is ; ++i )
    for( int j = 0 , js = target_sites[i].size() ; j < js ; ++j )
      delete target_sites[i][j];

  // keep the best (by RMS) distinct hits for the conformations.
  if( hits.size() > 1 ) {
    switch( hits_to_output ) {
    case BEST_HITS_ONLY :
      keep_best_hits( hits );
      break;
    case ONE_HIT_ONLY :
      hits.erase( hits.begin() + 1 , hits.end() );
      break;
    case ALL_HITS :
      break;
    }
  }

  return true;

}


// ********************************************************************
void initialise_mol_and_scores_outstreams( const string &output_file_root ,
                                           const string &mol_file ,
                                           bool scores_only ,
                                           bool comma_output ,
                                           const vector<string> &grid_vol_files ,
                                           PPhoreQuery &query ,
                                           oemolostream &oms ,
                                           ofstream &scores_stream ) {

  string scores_file = output_file_root + ".scores";

  if( !scores_only )
    oms.open( mol_file.c_str() );

  scores_stream.open( scores_file.c_str() );
  if( !scores_stream.good() )
    throw DACLIB::FileWriteOpenError( scores_file.c_str() );
  string sep = comma_output ? "," : " ";
  scores_stream << "Molecule_Name" << sep << "RMS" << sep << "Excluded_Vol"
                << sep << "Soft_Exc_Vols";

  for( int i = 0 , is = query.distances().size() ; i < is ; ++i )
    scores_stream << sep << "Distance_" << i + 1;

  for( int i = 0 , is = query.angles().size() ; i < is ; ++i )
    scores_stream << sep << "Angle_" << i + 1;

  for( int i = 0 , is = query.torsions().size() ; i < is ; ++i )
    scores_stream << sep << "Torsion_" << i + 1;

  for( int i = 0 , is = grid_vol_files.size() ; i < is ; ++i ) {
    boost::filesystem::path vp( grid_vol_files[i] );
    scores_stream << sep << vp.filename();
  }

  scores_stream << endl;

}

// ********************************************************************
void output_hits_to_mol_and_scores( vector<boost::shared_ptr<PluralityHit> > &hits ,
                                    bool scores_only , bool comma_output ,
                                    oemolostream &oms , ofstream &scores_stream ) {

  string sep = comma_output ? "," : " ";

  for( int i = 0 , is = hits.size() ; i < is ; ++i ) {
    if( !scores_only )
      oms << *(hits[i]->get_ov_conf());
    hits[i]->write_to_stream( scores_stream , sep );
  }

}

// ********************************************************************
void prepare_protein_grid( const string &protein_file ,
                           boost::shared_ptr<DACLIB::VolumeGrid> &protein_grid ) {

  if( !protein_file.empty() ) {
    oemolistream ims( protein_file.c_str() );
    if( !ims )
      throw DACLIB::FileReadOpenError( protein_file.c_str() );
    boost::shared_ptr<OEMolBase> protein( OENewMolBase( OEMolBaseType::OEDefault ) );
    ims >> *protein;
    if( *protein )
      protein_grid = boost::shared_ptr<DACLIB::VolumeGrid>( DACLIB::prepare_mol_grid( protein.get() ) );
  }

#ifdef NOTYET
  if( protein_grid )
    cout << "Protein volume : " << protein_grid->solid_volume()
         << " and protein surface volume : " << protein_grid->surface_volume()
         << endl;
#endif

}

// ********************************************************************
void prepare_soft_exc_vol_grid( PPhoreQuery &query ,
                                boost::shared_ptr<DACLIB::VolumeGrid> &soft_exc_vol_grid ) {

  if( query.soft_exc_vols().empty() )
    return;

  float blf[3] = { numeric_limits<float>::max() , numeric_limits<float>::max() ,
		   numeric_limits<float>::max() };
  float trb[3] = { -numeric_limits<float>::max() , -numeric_limits<float>::max() ,
                   -numeric_limits<float>::max() };

  vector<BTV>::const_iterator p , ps;
  for( p = query.soft_exc_vols().begin() , ps = query.soft_exc_vols().end() ;
       p != ps ; ++p ) {
    if( p->get<0>() - p->get<3>() < blf[0] ) blf[0] = p->get<0>() - p->get<3>();
    if( p->get<1>() - p->get<3>() < blf[1] ) blf[1] = p->get<1>() - p->get<3>();
    if( p->get<2>() - p->get<3>() < blf[2] ) blf[2] = p->get<2>() - p->get<3>();

    if( p->get<0>() + p->get<3>() > trb[0] ) trb[0] = p->get<0>() + p->get<3>();
    if( p->get<1>() + p->get<3>() > trb[1] ) trb[1] = p->get<1>() + p->get<3>();
    if( p->get<2>() + p->get<3>() > trb[2] ) trb[2] = p->get<2>() + p->get<3>();
  }

  soft_exc_vol_grid = boost::shared_ptr<DACLIB::VolumeGrid>( new DACLIB::VolumeGrid( blf , trb ) );
  soft_exc_vol_grid->drop_spheres_in( query.soft_exc_vols() );

#ifdef NOTYET
  cout << "soft_exc_vol_grid vols : " << soft_exc_vol_grid->solid_volume()
       << " and " << soft_exc_vol_grid->surface_volume() << endl;
#endif

}

// ********************************************************************
void prepare_volume_grids( const string &protein_file , PPhoreQuery &query ,
                           boost::shared_ptr<DACLIB::VolumeGrid> &protein_grid ,
                           boost::shared_ptr<DACLIB::VolumeGrid> &soft_exc_vol_grid ) {

  prepare_protein_grid( protein_file , protein_grid );

  prepare_soft_exc_vol_grid( query , soft_exc_vol_grid );

}

// ********************************************************************
void output_setup( PluralitySettings &plurality_settings ,
                   PharmPoint &pharm_points , PPhoreQuery &query ,
                   oemolostream &mol_out_stream , ofstream &out_stream ) {

  string output_file_root , output_file_ext;
  DACLIB::split_filename( plurality_settings.output_file() , output_file_root ,
                          output_file_ext );

  initialise_mol_and_scores_outstreams( output_file_root ,
					plurality_settings.output_file() ,
					plurality_settings.scores_only() ,
					plurality_settings.comma_output() ,
					plurality_settings.grid_vol_files() ,
					query , mol_out_stream , out_stream );

}

// ********************************************************************
void final_remarks( int mol_count , int conf_count , int hit_count ,
                    int hit_mol_count ) {

  cout << "Searched " << mol_count << " molecules with " << conf_count
       << " conformations giving " << hit_count
       << " hits from " << hit_mol_count << " of them." << endl;

}

// ********************************************************************
void make_not_smarts_queries( const vector<pair<string,string> > &not_smarts_list ,
                              vector<pair<string,string> > &smarts_sub_defn ,
                              vector<boost::shared_ptr<OESubSearch> > &not_smarts_subs ) {

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
    not_smarts_subs.push_back( boost::shared_ptr<OESubSearch>( new OESubSearch( exp_smarts.c_str() ) ) );
  }

}

// ********************************************************************
void setup_search( PluralitySettings &plurality_settings ,
                   vector<pair<string,string> > &input_smarts ,
                   vector<pair<string,string> > &smarts_sub_defn ,
                   PharmPoint &pharm_points ,
                   PPhoreQuery &query ,
                   boost::shared_ptr<OEMolBase> &protein ,
                   vector<string> &mol_subset ) {

  read_smarts_file( plurality_settings.smarts_file() , input_smarts ,
                    smarts_sub_defn );
  if( !plurality_settings.not_smarts_file().empty() ) {
    read_smarts_file( plurality_settings.not_smarts_file() ,
                      plurality_settings.not_smarts_list() , smarts_sub_defn );
  }

  try {
    pharm_points.read_points_file( plurality_settings.points_file() );
  } catch( DACLIB::FileReadOpenError &e ) {
    cerr << e.what() << endl;
    exit( 1 );
  }

  read_protein_file( plurality_settings.protein_file() , protein );

  if( !plurality_settings.subset_file().empty() )
    read_subset_file( plurality_settings.subset_file() , mol_subset );

  try {
    query.read_query_file( plurality_settings.pphore_file() );
  } catch( DACLIB::FileReadOpenError &e ) {
    cerr << e.what() << endl << endl;
    exit( 1 );
  } catch( PPhoreQueryFileReadError &e ) {
    cerr << e.what() << endl << endl;
    exit( 1 );
  } catch( PPhoreQueryFileError &e ) {
    cerr << e.what() << endl << endl;
    exit( 1 );
  }

}

// ********************************************************************
void report_progress( int mol_count , int conf_count ,
                      int hit_count , int hit_mol_count ,
                      ostream &os ) {

  if( ( mol_count < 100 && !( mol_count % 10 ) ) ||
      ( mol_count < 1000 && !( mol_count % 100 ) ) ||
      ( mol_count >= 1000 && mol_count < 10000 && !( mol_count % 1000 ) ) ||
      ( mol_count >= 10000 && mol_count < 100000 && !( mol_count % 10000 ) ) ||
      !( mol_count % 100000 ) ) {
    os << "Done " << mol_count << " molecules with " << conf_count
       << " conformations giving " << hit_count
       << " hits from " << hit_mol_count << " of them." << endl;
  }

}

// ********************************************************************
void serial_plurality_search( PluralitySettings &plurality_settings ) {

#ifdef NOTYET
  cout << "serial plurality search" << endl;
#endif

  vector<pair<string,string> > input_smarts , smarts_sub_defn;
  PharmPoint pharm_points;
  PPhoreQuery query;
  boost::shared_ptr<OEMolBase> protein;
  vector<string> mol_subset;
  oemolistream ims;

  setup_search( plurality_settings , input_smarts , smarts_sub_defn ,
                pharm_points , query , protein , mol_subset );

  OEMol target_mol;

  oemolostream mol_out_stream;
  ofstream out_stream;
  output_setup( plurality_settings , pharm_points , query ,
                mol_out_stream , out_stream );

  boost::shared_ptr<DACLIB::VolumeGrid> protein_grid;
  boost::shared_ptr<DACLIB::VolumeGrid> soft_exc_vol_grid;
  prepare_volume_grids( plurality_settings.protein_file() ,
                        query , protein_grid , soft_exc_vol_grid );

  // these are arbitrary grids for volume scores, prepared outside the program
  // probably by psg or prepare_scoring_grid
  vector<pair<string,DACLIB::VolumeGrid *> > score_vol_grids;
  DACLIB::read_grids( plurality_settings.grid_vol_files() , score_vol_grids );

  // SMARTS that the hits can't match
  vector<boost::shared_ptr<OESubSearch> > not_smarts_subs;
  try {
    make_not_smarts_queries( plurality_settings.not_smarts_list() ,
			     smarts_sub_defn , not_smarts_subs );
  } catch( DACLIB::SMARTSDefnError &e ) {
    cerr << e.what() << "From " << plurality_settings.not_smarts_file()
         << endl;
    exit( 1 );
  }

  int mol_count = 0 , conf_count = 0 , hit_mol_count = 0 , hit_count = 0;
  vector<string> db_files = plurality_settings.db_files();
  BOOST_FOREACH( string db_file , db_files ) {

    cout << "Database file : " << db_file << endl;
    open_databasefile( db_file , false , ims );

    while( ims >> target_mol ) {
      OEAssignAromaticFlags( target_mol , OEAroModelDaylight );
      vector<boost::shared_ptr<PluralityHit> > hits;
#ifdef NOTYET
      cout << "Searching " << target_mol.GetTitle() << endl;
#endif
      if( search_target_mol( target_mol , query , input_smarts , smarts_sub_defn ,
                             pharm_points , mol_subset , protein_grid ,
                             soft_exc_vol_grid , score_vol_grids ,
                             plurality_settings.hits_to_output() ,
                             not_smarts_subs , hits ) ) {
        hit_count += hits.size();
        if( !hits.empty() )
          ++hit_mol_count;
        ++mol_count;
        conf_count += target_mol.NumConfs();
      }
      output_hits_to_mol_and_scores( hits , plurality_settings.scores_only() ,
				     plurality_settings.comma_output() ,
				     mol_out_stream , out_stream );

      report_progress( mol_count , conf_count , hit_count , hit_mol_count , cout );
    }
    ims.close();

  }

  final_remarks( mol_count , conf_count , hit_count , hit_mol_count );

}

// ********************************************************************
void send_query_details_to_slaves( PluralitySettings &plurality_settings ,
                                   vector<int> &slave_tids ) {

  send_openeye_license_to_slaves( slave_tids );
  send_cwd_to_slaves( slave_tids );

  pvm_initsend( PvmDataDefault );
  pvm_pkstr( const_cast<char *>( string( "Query_Details" ).c_str() ) );

  plurality_settings.pack_contents_into_pvm_buffer();

  pvm_mcast( &slave_tids[0] , slave_tids.size() , 0 );

}

// ********************************************************************
// send the given slave a request for results and store the answer
void get_result_from_slave( vector<boost::shared_ptr<PluralityHit> > &hits ,
			    int &mol_count , int &conf_count ,
			    int &hit_count , int &hit_mol_count ) {

  int num_to_rec;
  pvm_upkint( &num_to_rec , 1 , 1 );
  for( int i = 0 ; i < num_to_rec ; ++i ) {
    hits.push_back( boost::shared_ptr<PluralityHit>( new PluralityHit ) );
    hits.back()->unpack_from_pvm_buffer();
  }

  int stats[4];
  pvm_upkint( stats , 4 , 1 );
  mol_count += stats[0];
  conf_count += stats[1];
  hit_count += stats[2];
  hit_mol_count += stats[3];

}

// ********************************************************************
void get_search_stats_from_slave( int &mol_count , int &conf_count ,
                                  int &hit_count , int &hit_mol_count ) {

  // message buffer already open
  int stats[4];
  pvm_upkint( stats , 4 , 1 );
  mol_count += stats[0];
  conf_count += stats[1];
  hit_count += stats[2];
  hit_mol_count += stats[3];

}

// ********************************************************************
// receives hits from slaves, and writes them to file.
void get_results_from_slaves( vector<int> &slave_tids ,
                              PluralitySettings &ps ,
                              oemolostream &mol_out_stream ,
                              ofstream &out_stream ,
                              int &mol_count , int &conf_count ,
                              int &hit_count , int &hit_mol_count ) {

#ifdef NOTYET
  cout << "get_results_from_slave : " << pvm_mytid() << endl;
#endif

  while( 1 ) {
    if( slave_tids.empty() ) {
      break;
    }

    int bufid = pvm_recv( -1 , -1 );
#ifdef NOTYET
    cout << "\nNew message : bufid : " << bufid << endl;
#endif
    int nbytes , msgtag , dead_tid;
    pvm_bufinfo( bufid , &nbytes , &msgtag , &dead_tid );
#ifdef NOTYET
    cout << "nbytes = " << nbytes << "  msgtag = " << msgtag << "  dead_tid = " << dead_tid << endl;
#endif
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
             << "machines." << endl;
        pvm_exit();
        exit( 1 );
      }

    } else {

      // if we're here, the message wasn't a failure message. It could be
      // a progress report or a load of results. done_tid isn't
      // filled in this case, so need to get it from the message
      char msg[1000]; // it'll be big enough for the message header
      pvm_upkstr( msg );
#ifdef NOTYET
      cout << "Message : " << msg << endl;
#endif
      if( !strcmp( msg , "Results" ) ) {
        int slave_tid;
        pvm_upkint( &slave_tid , 1 , 1 );
        cout << "Results from slave : " << slave_tid << endl;
        vector<boost::shared_ptr<PluralityHit> > hits;
        get_result_from_slave( hits , mol_count , conf_count , hit_count ,
			       hit_mol_count );
	output_hits_to_mol_and_scores( hits , ps.scores_only() ,
				       ps.comma_output() , mol_out_stream ,
				       out_stream );
        slave_tids.erase( find( slave_tids.begin() , slave_tids.end() , slave_tid ) );
        cout << "Number of slaves still to report : " << slave_tids.size() << endl;
      } else if( !strcmp( msg , "Search Stats" ) ) {
        get_search_stats_from_slave( mol_count , conf_count ,
				     hit_count , hit_mol_count );
      } else if( !strcmp( "Progress Report" , msg ) ) {
        pvm_upkint( &done_tid , 1 , 1 );
        string prog_rep;
        DACLIB::unpack_string( prog_rep );
        cout << done_tid << " : " << prog_rep;
      }

    }
  }

}

// ********************************************************************
void parallel_plurality_search( PluralitySettings &plurality_settings ,
                                vector<int> &slave_tids ) {

  send_query_details_to_slaves( plurality_settings , slave_tids );

  // send the start mol and step size to each slave, which will
  // also kick-start the search
  send_database_details_to_slaves( slave_tids );

  // now the slaves are running, we'll be kicking our heals for a bit
  vector<pair<string,string> > input_smarts , smarts_sub_defn;
  PharmPoint pharm_points;
  PPhoreQuery query;
  boost::shared_ptr<OEMolBase> protein;
  vector<string> mol_subset;
  setup_search( plurality_settings , input_smarts , smarts_sub_defn ,
                pharm_points , query , protein , mol_subset );

  oemolostream mol_out_stream;
  ofstream out_stream;
  output_setup( plurality_settings , pharm_points , query ,
                mol_out_stream , out_stream );

  int mol_count = 0 , conf_count = 0 , hit_count = 0 , hit_mol_count = 0;
  get_results_from_slaves( slave_tids , plurality_settings ,
                           mol_out_stream , out_stream ,
                           mol_count , conf_count , hit_count , hit_mol_count );
  send_finished_messages( slave_tids );

  final_remarks( mol_count , conf_count , hit_count , hit_mol_count );

}

// ********************************************************************
void send_results_to_master( int master_tid ,
			     int mol_count , int conf_count ,
			     int hit_count , int hit_mol_count ,
                             vector<boost::shared_ptr<PluralityHit> > &hits ) {

  pvm_initsend( PvmDataDefault );
  pvm_pkstr( const_cast<char *>( string( "Results" ).c_str() ) );
  int my_tid = pvm_mytid();
  pvm_pkint( &my_tid , 1 , 1 );

  int num_to_send = hits.size();
  cout << "sending " << num_to_send << " hits " << endl;
  pvm_pkint( &num_to_send , 1 , 1 );

  for( int i = 0 ; i < num_to_send ; ++i ) {
    hits[i]->pack_into_pvm_buffer();
  }

  int stats[4] = { mol_count , conf_count , hit_count , hit_mol_count };
  pvm_pkint( stats , 4 , 1 );
  pvm_send( master_tid , 0 );

}

// ********************************************************************
void send_search_stats_to_master( int master_tid , int mol_count ,
                                  int conf_count , int hit_count ,
                                  int hit_mol_count ) {

  pvm_initsend( PvmDataDefault );
  pvm_pkstr( const_cast<char *>( string( "Search Stats" ).c_str() ) );
  int stats[4] = { mol_count , conf_count , hit_count , hit_mol_count };
  pvm_pkint( stats , 4 , 1 );
  pvm_send( master_tid , 0 );

}

// ********************************************************************
void search_databases( int db_start , int db_step , const vector<string> &db_files ,
                       PPhoreQuery &query ,
                       vector<pair<string,string> > &input_smarts ,
                       vector<pair<string,string> > &smarts_sub_defn ,
                       PharmPoint &pharm_points ,
                       const vector<string> &mol_subset ,
                       boost::shared_ptr<DACLIB::VolumeGrid> &protein_grid ,
                       boost::shared_ptr<DACLIB::VolumeGrid> &soft_exc_vol_grid ,
                       vector<pair<string,DACLIB::VolumeGrid *> > score_vol_grids ,
                       HITS_TO_OUTPUT hits_to_output ,
                       vector<boost::shared_ptr<OESubSearch> > &not_smarts_subs ,
                       vector<boost::shared_ptr<PluralityHit> > &hits ,
                       int &mol_count , int &conf_count ,
                       int &hit_mol_count , int &hit_count ) {

  OEMol mol;
  oemolistream db_ims;

  int master_tid = pvm_parent();

  BOOST_FOREACH( string db_file , db_files ) {

    cout << "Database file : " << db_file << endl;
    open_databasefile( db_file , false , db_ims );
    if( static_cast<int>( step_oemolstream( db_ims , db_start ) ) < db_start ) {
      cout << "Problem with database " << db_file
           << " : couldn't read to " << db_start << " record." << endl;
      pvm_exit();
      exit( 1 );
    }

    while( db_ims >> mol ) {
      vector<boost::shared_ptr<PluralityHit> > these_hits;
      if( search_target_mol( mol , query , input_smarts ,
                             smarts_sub_defn , pharm_points , mol_subset ,
                             protein_grid , soft_exc_vol_grid , score_vol_grids ,
                             hits_to_output , not_smarts_subs , these_hits ) ) {
        hit_count += these_hits.size();
        if( !these_hits.empty() )
          ++hit_mol_count;
        ++mol_count;
        conf_count += mol.NumConfs();
      }
      if( !these_hits.empty() ) {
        hits.insert( hits.end() , these_hits.begin() , these_hits.end() );
      }
      // move forward step-1 molecules, unless it hits the end of file early
      if( static_cast<int>( step_oemolstream( db_ims , db_step - 1 ) ) < db_step - 1 ) {
        break;
      }
      ostringstream oss;
      report_progress( mol_count , conf_count , hit_count , hit_mol_count , oss );
      if( !oss.str().empty() ) {
        send_progress_to_master( oss.str() );
      }
      // see if the master has died, in which case stop
      // the documentation suggests that I should be able to listen just to
      // messages from the master process by passing master_tid as the 1st
      // argument of pvm_nrecv, but that doesn't seem to work (at least in
      // 3.4.6). Nothing else should be passing messages in any case.
      int bufid = pvm_nrecv( -1 , -1 );
      if( bufid > 0 ) {
        // bufid of 0 means nothing received
        int dead_tid;
        if( DACLIB::was_it_a_pvm_failure_message( bufid , dead_tid ) &&
            dead_tid == master_tid ) {
          cerr << "AWOOGA : The master is dead.  Following him on into the"
               << " afterlife." << endl;
          pvm_exit();
          exit( 1 );
        }
      }

    }

    db_ims.close();
  }

}

// ********************************************************************
void slave_event_loop() {

  char     msg1[1000];

  PluralitySettings plurality_settings;
  vector<pair<string,string> > input_smarts , smarts_sub_defn;
  PharmPointPVM pharm_points;
  PPhoreQuery query;
  boost::shared_ptr<OEMolBase> protein;
  vector<string> mol_subset;
  boost::shared_ptr<DACLIB::VolumeGrid> protein_grid;
  boost::shared_ptr<DACLIB::VolumeGrid> soft_exc_vol_grid;
  vector<boost::shared_ptr<OESubSearch> > not_smarts_subs;
  vector<boost::shared_ptr<PluralityHit> > hits;
  vector<pair<string,DACLIB::VolumeGrid *> > score_vol_grids;
  vector<string> db_files;
  HITS_TO_OUTPUT hits_to_output = BEST_HITS_ONLY;
  int db_start , db_step;

  int master_tid = pvm_parent();
  if( PvmNoParent == master_tid ) {
    // It's an error, though probably a deliberate one
    cerr << "ERROR - strange PVM behaviour for plurality." << endl;
    exit( 1 );
  }

  // we want to know if the master dies
  pvm_notify( PvmTaskExit , TASK_DIED , 1 , &master_tid );

  int mol_count = 0 , conf_count = 0 , hit_mol_count = 0 , hit_count = 0;
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
    } else if( !strcmp( "Query_Details" , msg1 ) ) {
      plurality_settings.unpack_contents_from_pvm_buffer();
      setup_search( plurality_settings , input_smarts , smarts_sub_defn ,
                    pharm_points , query , protein , mol_subset );
    } else if( !strcmp( "Database_Steps" , msg1 ) ) {
      receive_database_details( db_start , db_step );
      search_databases( db_start , db_step , plurality_settings.db_files() ,
                        query , input_smarts ,
                        smarts_sub_defn , pharm_points , mol_subset ,
                        protein_grid , soft_exc_vol_grid ,
                        score_vol_grids , hits_to_output ,
                        not_smarts_subs , hits , mol_count , conf_count ,
                        hit_mol_count , hit_count );
      send_results_to_master( master_tid , mol_count , conf_count ,
			      hit_count , hit_mol_count , hits );
    } else if( !strcmp( "New_CWD" , msg1 ) ) {
      receive_new_cwd();
    }

  }

}
