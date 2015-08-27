//
// file plurality_subs.cc
// David Cosgrove
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

#include <mpi.h>

#include <boost/bind.hpp>
//#include <boost/cast.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>

#include "BasePPhoreSite.H"
#include "FileExceptions.H"
#include "PharmPoint.H"
#include "PluralityHit.H"
#include "PluralitySettings.H"
#include "PPhoreQuery.H"
#include "SinglePPhoreSite.H"
#include "SMARTSExceptions.H"

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
void set_environment_var_from_mpi();

// In mpi_string_subs.cc
void mpi_send_string( const string &str , int dest_rank );
void mpi_send_strings_vector( const vector<string> &strs , int dest_rank );
void mpi_rec_string( int source_rank , string &str );
void mpi_rec_strings_vector( int source_rank , vector<string> &strs );

void split_filename( const string &filename , string &file_root ,
                     string &file_ext );
}

// in send_messages_to_tp_slaves.cc
void pack_smarts_defs_into_pvm_buffer( vector<pair<string,string> > &input_smarts ,
                                       vector<pair<string,string> > &smarts_sub_defn );
void unpack_smarts_defs_from_pvm_buffer( vector<pair<string,string> > &input_smarts ,
                                         vector<pair<string,string> > &smarts_sub_defn );
void send_openeye_license_to_slaves( int world_size );
void send_cwd_to_slaves( int world_size );
void send_finished_messages( int world_size );
void send_database_details_to_slaves( int world_size );
void receive_new_cwd();
void receive_database_details( int &db_start , int &db_step );

typedef enum { MOLS_AND_SCORES } OUTPUT_FORMAT;

typedef boost::shared_ptr<OEMol> pOEMol;

// in step_oemolstream.cc
namespace DACLIB {
void open_databasefile( const string &db_file ,
                        bool single_conf_mols ,
                        oemolistream *&ims );
pOEMol read_nth_mol_from_oemolistream( unsigned int next_mol ,
                                       const vector<string> &db_files ,
                                       bool single_conf_mols );
}

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
OUTPUT_FORMAT initialise_output_streams( const string &output_file_root ,
                                         const string &output_file_ext ,
                                         PluralitySettings &plurality_settings ,
                                         PharmPoint &pharm_points ,
                                         PPhoreQuery &query ,
                                         oemolostream  &mol_out_stream ,
                                         ofstream &out_stream ) {

  OUTPUT_FORMAT output_format;

  // default to mols and scores style, hoping that we got some sensible
  // file names.
  output_format = MOLS_AND_SCORES;
  initialise_mol_and_scores_outstreams( output_file_root ,
					plurality_settings.output_file() ,
					plurality_settings.scores_only() ,
					plurality_settings.comma_output() ,
					plurality_settings.grid_vol_files() ,
					query , mol_out_stream , out_stream );

  return output_format;

}

// ********************************************************************
void output_setup( PluralitySettings &plurality_settings ,
                   PharmPoint &pharm_points , PPhoreQuery &query ,
                   OUTPUT_FORMAT &output_format ,
                   oemolostream &mol_out_stream , ofstream &out_stream ) {

  string output_file_root , output_file_ext;
  DACLIB::split_filename( plurality_settings.output_file() , output_file_root ,
                          output_file_ext );

  output_format =
      initialise_output_streams( output_file_root , output_file_ext ,
                                 plurality_settings , pharm_points , query ,
                                 mol_out_stream , out_stream );

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
void report_progress( unsigned int mol_count , unsigned int conf_count ,
                      unsigned int hit_count , unsigned int hit_mol_count ,
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
  oemolistream *ims = 0;

  setup_search( plurality_settings , input_smarts , smarts_sub_defn ,
                pharm_points , query , protein , mol_subset );

  OEMol target_mol;

  oemolostream mol_out_stream;
  ofstream out_stream;
  OUTPUT_FORMAT output_format;
  output_setup( plurality_settings , pharm_points , query ,
                output_format , mol_out_stream , out_stream );

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

  unsigned int mol_count = 0 , conf_count = 0 , hit_mol_count = 0 , hit_count = 0;
  vector<string> db_files = plurality_settings.db_files();
  BOOST_FOREACH( string db_file , db_files ) {

    cout << "Database file : " << db_file << endl;
    ims = new oemolistream;
    DACLIB::open_databasefile( db_file , false , ims );

    while( *ims >> target_mol ) {
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
        if( !hits.empty() ) {
          ++hit_mol_count;
        }
        ++mol_count;
        conf_count += target_mol.NumConfs();
      }
      switch( output_format ) {
      case MOLS_AND_SCORES :
        output_hits_to_mol_and_scores( hits , plurality_settings.scores_only() ,
                                       plurality_settings.comma_output() ,
                                       mol_out_stream , out_stream );
        break;
      }

      report_progress( mol_count , conf_count , hit_count , hit_mol_count , cout );
    }
    ims->close();
    delete ims;

  }

  final_remarks( mol_count , conf_count , hit_count , hit_mol_count );

}

// ********************************************************************
void send_query_details_to_slaves( PluralitySettings &plurality_settings ,
                                   int world_size ) {

  send_openeye_license_to_slaves( world_size );
  send_cwd_to_slaves( world_size );

  for( int i = 1 ; i < world_size ; ++i ) {
    DACLIB::mpi_send_string( string( "Query_Details" ) , i );
    plurality_settings.send_contents_via_mpi( i );
  }

}

// ********************************************************************
// send the given slave a request for results and store the answer
void get_result_from_slave( int source_rank ,
                            vector<boost::shared_ptr<PluralityHit> > &hits ,
                            int &mol_count , int &conf_count ,
                            int &hit_count , int &hit_mol_count ) {

  unsigned int num_to_rec;
  MPI_Recv( &num_to_rec , 1 , MPI_UNSIGNED , source_rank , 0 ,
            MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  for( unsigned int i = 0 ; i < num_to_rec ; ++i ) {
    hits.push_back( boost::shared_ptr<PluralityHit>( new PluralityHit ) );
    hits.back()->receive_via_mpi( source_rank );
  }

  int stats[4];
  MPI_Recv( stats , 4 , MPI_INT , source_rank , 0 ,
            MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  mol_count += stats[0];
  conf_count += stats[1];
  hit_count += stats[2];
  hit_mol_count += stats[3];

}

// ********************************************************************
void get_search_stats_from_slave( int source_rank ,
                                  int &mol_count , int &conf_count ,
                                  int &hit_count , int &hit_mol_count ) {

  int stats[4];
  MPI_Recv( stats , 4 , MPI_UNSIGNED , source_rank , 0 ,
            MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  mol_count += stats[0];
  conf_count += stats[1];
  hit_count += stats[2];
  hit_mol_count += stats[3];

}

// ********************************************************************
// receives hits from slave, and writes them to file.
void get_results_from_slave( int source_rank ,
                             PluralitySettings &ps ,
                             OUTPUT_FORMAT output_format ,
                             oemolostream &mol_out_stream ,
                             ofstream &out_stream ,
                             unsigned int &mol_count , unsigned int &conf_count ,
                             unsigned int &hit_count , unsigned int &hit_mol_count ) {

  vector<boost::shared_ptr<PluralityHit> > hits;
  unsigned int num_to_rec;
  MPI_Recv( &num_to_rec , 1 , MPI_UNSIGNED , source_rank , 0 ,
            MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  for( unsigned int i = 0 ; i < num_to_rec ; ++i ) {
    hits.push_back( boost::shared_ptr<PluralityHit>( new PluralityHit ) );
    hits.back()->receive_via_mpi( source_rank );
  }

  unsigned int stats[4];
  MPI_Recv( stats , 4 , MPI_UNSIGNED , source_rank , 0 ,
            MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  mol_count += stats[0];
  conf_count += stats[1];
  hit_count += stats[2];
  hit_mol_count += stats[3];

  switch( output_format ) {
  case MOLS_AND_SCORES :
    output_hits_to_mol_and_scores( hits , ps.scores_only() ,
                                   ps.comma_output() , mol_out_stream ,
                                   out_stream );
    break;
  }

}

// ********************************************************************
void send_parallel_searches( int world_size , PluralitySettings &ps ,
                             OUTPUT_FORMAT output_format ,
                             oemolostream &mol_out_stream ,
                             ofstream &out_stream ,
                             unsigned int &mol_count , unsigned int &conf_count ,
                             unsigned int &hit_count , unsigned int &hit_mol_count ) {

  int num_finished = 0;
  int num_slaves = world_size - 1; // process 0 is the master
  unsigned int curr_mol = 0;

  while( num_finished < num_slaves ) {

    MPI_Status status;
    MPI_Probe( MPI_ANY_SOURCE , 0 , MPI_COMM_WORLD , &status );
    string msg;
    DACLIB::mpi_rec_string( status.MPI_SOURCE , msg );
#ifdef NOTYET
    cout << "Message : " << msg << " from " << status.MPI_SOURCE << endl;
#endif
    if( string( "Results" ) == msg ) {
      get_results_from_slave( status.MPI_SOURCE , ps , output_format ,
                              mol_out_stream , out_stream ,
                              mol_count , conf_count , hit_count ,
                              hit_mol_count );
      report_progress( mol_count , conf_count , hit_count , hit_mol_count , cout );
    } else if( string( "Send_Mol_Num" ) == msg ) {
      DACLIB::mpi_send_string( string( "Next_Mol_Num" ) , status.MPI_SOURCE );
      MPI_Send( &curr_mol , 1 , MPI_UNSIGNED , status.MPI_SOURCE , 0 , MPI_COMM_WORLD );
      ++curr_mol;
#ifdef NOTYET
      cout << "Sent " << curr_mol - 1 << " to " << status.MPI_SOURCE << endl;
#endif
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
void parallel_plurality_search( PluralitySettings &plurality_settings ,
                                int world_size ) {

  send_query_details_to_slaves( plurality_settings , world_size );

  // now the slaves are running, we'll be kicking our heels for a bit
  vector<pair<string,string> > input_smarts , smarts_sub_defn;
  PharmPoint pharm_points;
  PPhoreQuery query;
  boost::shared_ptr<OEMolBase> protein;
  vector<string> mol_subset;
  setup_search( plurality_settings , input_smarts , smarts_sub_defn ,
                pharm_points , query , protein , mol_subset );

  oemolostream mol_out_stream;
  ofstream out_stream;
  OUTPUT_FORMAT output_format;
  output_setup( plurality_settings , pharm_points , query ,
                output_format , mol_out_stream , out_stream );

  unsigned int mol_count = 0 , conf_count = 0 , hit_count = 0 , hit_mol_count = 0;
  send_parallel_searches( world_size , plurality_settings ,
                          output_format , mol_out_stream , out_stream ,
                          mol_count , conf_count , hit_count ,
                          hit_mol_count );

  send_finished_messages( world_size );

  final_remarks( mol_count , conf_count , hit_count , hit_mol_count );

}

// ********************************************************************
void send_results_to_master( unsigned int mol_count , unsigned int conf_count ,
                             unsigned int hit_count , unsigned int hit_mol_count ,
                             vector<boost::shared_ptr<PluralityHit> > &hits ) {

  DACLIB::mpi_send_string( string( "Results" ) , 0 );
  unsigned int num_to_send = hits.size();
  MPI_Send( &num_to_send , 1 , MPI_UNSIGNED , 0 , 0 , MPI_COMM_WORLD );

  for( unsigned int i = 0 ; i < num_to_send ; ++i ) {
    hits[i]->send_via_mpi( 0 );
  }

  unsigned int stats[4] = { mol_count , conf_count , hit_count , hit_mol_count };
  MPI_Send( stats , 4 , MPI_UNSIGNED , 0 , 0 , MPI_COMM_WORLD );

}

// ********************************************************************
bool search_database( int next_mol , PluralitySettings &ps ,
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
                      unsigned int &mol_count , unsigned int &conf_count ,
                      unsigned int &hit_mol_count , unsigned int &hit_count ) {

  pOEMol tmol = DACLIB::read_nth_mol_from_oemolistream( next_mol ,
                                                        ps.db_files() , false );

  if( !tmol ) {

    // we're done
    cout << "sending Search_Finished" << endl;
    DACLIB::mpi_send_string( string( "Search_Finished" ) , 0 );
    return false;

  } else {

#ifdef NOTYET
    cout << "Searching " << tmol->GetTitle() << endl;
#endif
    if( search_target_mol( *tmol , query , input_smarts ,
                           smarts_sub_defn , pharm_points , mol_subset ,
                           protein_grid , soft_exc_vol_grid , score_vol_grids ,
                           hits_to_output , not_smarts_subs , hits ) ) {
      hit_count += hits.size();
      if( !hits.empty() ) {
        ++hit_mol_count;
      }
      ++mol_count;
      conf_count += tmol->NumConfs();

    }

  }

  return true;

}

// ********************************************************************
void slave_event_loop() {

  string msg;

  PluralitySettings plurality_settings;
  vector<pair<string,string> > input_smarts , smarts_sub_defn;
  PharmPoint pharm_points;
  PPhoreQuery query;
  boost::shared_ptr<OEMolBase> protein;
  vector<string> mol_subset;
  boost::shared_ptr<DACLIB::VolumeGrid> protein_grid;
  boost::shared_ptr<DACLIB::VolumeGrid> soft_exc_vol_grid;
  vector<boost::shared_ptr<OESubSearch> > not_smarts_subs;
  vector<pair<string,DACLIB::VolumeGrid *> > score_vol_grids;
  vector<string> db_files;
  HITS_TO_OUTPUT hits_to_output = BEST_HITS_ONLY;

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
    } else if( string( "Query_Details" ) == msg ) {
      plurality_settings.receive_contents_via_mpi();
      setup_search( plurality_settings , input_smarts , smarts_sub_defn ,
                    pharm_points , query , protein , mol_subset );
      DACLIB::mpi_send_string( string( "Send_Mol_Num" ) , 0 );
    } else if( string( "Next_Mol_Num" ) == msg ) {
      vector<boost::shared_ptr<PluralityHit> > hits;
      unsigned int next_mol;
      MPI_Recv( &next_mol , 1 , MPI_UNSIGNED , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
      unsigned int mol_count = 0 , conf_count = 0 , hit_mol_count = 0 , hit_count = 0;
      if( search_database( next_mol , plurality_settings ,
                           query , input_smarts ,
                           smarts_sub_defn , pharm_points , mol_subset ,
                           protein_grid , soft_exc_vol_grid ,
                           score_vol_grids , hits_to_output ,
                           not_smarts_subs , hits , mol_count , conf_count ,
                           hit_mol_count , hit_count ) ) {
        send_results_to_master( mol_count , conf_count ,
                                hit_count , hit_mol_count , hits );
        hits.clear();
        DACLIB::mpi_send_string( string ( "Send_Mol_Num" ) , 0 );
      }
    } else if( string( "New_CWD" ) == msg ) {
      receive_new_cwd();
    }

  }

}
