//
// file plurality.cc
// Dave Cosgrove
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

#include "BasePPhoreSite.H"
#include "DefaultPointsDefs.H"
#include "FileExceptions.H"
#include "ParseSMARTSXML.H"
#include "PharmPointPVM.H"
#include "PluralitySettings.H"
#include "PPhoreQuery.H"

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
                                vector<int> &slave_tids );

void slave_event_loop();

namespace DACLIB {
  // in eponymous file
  void launch_pvm_slaves( const string &slave_name , int num_slave_procs ,
                          vector<int> &slave_tids );
  void launch_pvm_slaves( const string &slave_name , const string &hosts_file ,
                          vector<int> &slave_tids );
}

// **************************************************************************
int main( int argc , char **argv ) {

  cout << "Plurality v3.1." << endl
       << "Built " << BUILD_TIME << " using OEToolkits version "
       << OEChemGetRelease() << "." << endl << endl
       << "Copyright AstraZeneca 2009, 2010, 2011, 2012, 2014." << endl << endl;

  // was this launched as a PVM slave?  --IM-A-SLAVE is provided by
  // my launcher. We need to intercept it if it was, and send it straight into
  // listening mode.
  if( 2 == argc && !strcmp( "--IM-A-SLAVE" , argv[1] ) ) {
    slave_event_loop();
    exit( 0 );
  }

  PluralitySettings plurality_settings( argc , argv );
  if( !plurality_settings ) {
    plurality_settings.print_usage( cerr );
    exit( 1 );
  }

  vector<int> slave_tids;
  if( !plurality_settings.pvm_hosts_file().empty() ) {
    // launch pvm slaves using details in hosts file
    DACLIB::launch_pvm_slaves( plurality_settings.slave_name() ,
                               plurality_settings.pvm_hosts_file() ,
                               slave_tids );
    if( slave_tids.empty() ) {
      cerr << "No slaves launched." << endl;
      exit( 1 );
    }
  } else if( plurality_settings.num_slave_procs() > 0 ) {
    DACLIB::launch_pvm_slaves( plurality_settings.slave_name() ,
			       plurality_settings.num_slave_procs() , slave_tids );
    if( slave_tids.empty() ) {
      cerr << "No slaves launched." << endl;
      exit( 1 );
    }
  }

  vector<pair<string,string> > input_smarts , smarts_sub_defn;
  PharmPointPVM pharm_points;
  if( !plurality_settings.smarts_file().empty() && !plurality_settings.points_file().empty() ) {
    read_smarts_file( plurality_settings.smarts_file() , input_smarts ,
		      smarts_sub_defn );

    try {
      pharm_points.read_points_file( plurality_settings.points_file() );
    } catch( DACLIB::FileReadOpenError &e ) {
      cerr << e.what() << endl;
      exit( 1 );
    }
  } else {
    pharm_points.clear_data();
    pharm_points.read_points_xml_string( DACLIB::DEFAULT_POINT_DEFS );
    ParseSMARTSXML psx;
    psx.parse_string( DACLIB::DEFAULT_POINT_DEFS , input_smarts , smarts_sub_defn );
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
    exit( 1 );
  }

  if( slave_tids.empty() ) {
    try {
      serial_plurality_search( plurality_settings );
    } catch( DACLIB::FileReadOpenError &e ) {
      cerr << e.what() << endl;
      exit( 1 );
    }
  } else {
    try {
      parallel_plurality_search( plurality_settings , slave_tids );
    } catch( DACLIB::FileReadOpenError &e ) {
      cerr << e.what() << endl;
      exit( 1 );
    } 
  }

}
