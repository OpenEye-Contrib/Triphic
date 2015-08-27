//
// file plurality.cc
// David Cosgrove
// 13th September 2007
//
// This is main function for the program plurality.
// It is new plurality. v3.0, built to use all the stuff developed for
// grappel. That's why it's here, in the same directory tree. They have a lot
// in common.
// 

#include <iostream>
#include <string>
#include <vector>

#include <oechem.h>

#include <boost/shared_ptr.hpp>
// for the exceptions thrown by triphic_args_handler
#include <boost/program_options/options_description.hpp>

#include "BasePPhoreSite.H"
#include "FileExceptions.H"
#include "PharmPoint.H"
#include "PluralitySettings.H"
#include "PPhoreQuery.H"

#include <mpi.h>

using namespace std;
using namespace OEChem;

extern string BUILD_TIME; // in build_time.cc

// in plurality_subs.cc
void read_smarts_file( const string &smarts_file ,
                       vector<pair<string,string> > &input_smarts ,
                       vector<pair<string,string> > &smarts_sub_defn );
void read_protein_file( const string &protein_file ,
                        boost::shared_ptr<OEMolBase> &protein );
void read_subset_file( const string &subset_file ,
                       vector<string> &subset_names );
void serial_plurality_search( PluralitySettings &plurality_settings );
void parallel_plurality_search( PluralitySettings &plurality_settings ,
                                int world_size );

void slave_event_loop();

// **************************************************************************
int main( int argc , char **argv ) {

  cout << "Plurality v3.2." << endl
       << "Built " << BUILD_TIME << " using OEToolkits version "
       << OEChemGetRelease() << "." << endl << endl
       << "Copyright AstraZeneca 2009, 2010, 2011, 2012, 2014, 2015." << endl << endl;

  // sort out the  MPI environment
  MPI_Init( NULL , NULL );
  int world_rank , world_size;
  MPI_Comm_rank( MPI_COMM_WORLD , &world_rank );
  MPI_Comm_size( MPI_COMM_WORLD , &world_size );

#ifdef NOTYET
  cout << "world_rank : " << world_rank << endl
       << "world size : " << world_size << endl;
#endif

  // use 0th rank process as the master, the rest as slaves
  if( world_size > 1 && world_rank ) {
    slave_event_loop();
  } else {

    PluralitySettings plurality_settings( argc , argv );
    if( !plurality_settings ) {
      plurality_settings.print_usage( cerr );
      MPI_Finalize();
      exit( 1 );
    }

    vector<pair<string,string> > input_smarts , smarts_sub_defn;
    read_smarts_file( plurality_settings.smarts_file() , input_smarts ,
                      smarts_sub_defn );
    if( !plurality_settings.not_smarts_file().empty() )
      read_smarts_file( plurality_settings.not_smarts_file() ,
                        plurality_settings.not_smarts_list() , smarts_sub_defn );

    PharmPoint pharm_points;
    try {
      pharm_points.read_points_file( plurality_settings.points_file() );
    } catch( DACLIB::FileReadOpenError &e ) {
      cerr << e.what() << endl;
      exit( 1 );
    }

    boost::shared_ptr<OEMolBase> protein;
    read_protein_file( plurality_settings.protein_file() , protein );

    vector<string> mol_subset;
    if( !plurality_settings.subset_file().empty() )
      read_subset_file( plurality_settings.subset_file() , mol_subset );

    PPhoreQuery query;
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

    vector<string> query_pts = query.point_types();
    map<string,vector<string> > defined_pts = pharm_points.points_defs();
    bool err_flag( false );
    for( int i = 0 , is = query_pts.size() ; i < is ; ++i ) {
      if( defined_pts.end() == defined_pts.find( query_pts[i] ) ) {
        cerr << "ERROR : pharmacophore file calls for feature named " << query_pts[i]
                << " which is not defined in the points file." << endl;
        err_flag = true;
      }
    }
    if( err_flag ) {
      MPI_Finalize();
      exit( 1 );
    }

    if( 1 == world_size ) {
      try {
        serial_plurality_search( plurality_settings );
      } catch( DACLIB::FileReadOpenError &e ) {
        cerr << e.what() << endl;
        MPI_Finalize();
        exit( 1 );
      }
    } else {
      try {
        parallel_plurality_search( plurality_settings , world_size );
      } catch( DACLIB::FileReadOpenError &e ) {
        cerr << e.what() << endl;
        MPI_Finalize();
        exit( 1 );
      }
    }
  }
  MPI_Finalize();

}
