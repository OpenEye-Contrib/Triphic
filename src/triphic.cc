//
// file triphic.cc
// David Cosgrove
// AstraZeneca
// 31st July 2007
//
// This is the main function for the program triphic.
// This is new triphic, v3.0.

#include <iterator>
#include <iostream>
#include <string>
#include <vector>

#include <oechem.h>

#include <boost/shared_ptr.hpp>

#include "BasePPhoreSite.H"
#include "FileExceptions.H"
#include "PharmPointPVM.H"
#include "SMARTSExceptions.H"
#include "TriphicSettings.H"
#include "VolumeGrid.H"

using namespace std;
using namespace OEChem;

extern string BUILD_TIME; // in build_time.cc

// in triphic_subs.cc
void serial_triphic_search( TriphicSettings &triphic_settings );
void parallel_triphic_search( TriphicSettings &triphic_settings ,
                              vector<int> &slave_tids );

void slave_event_loop();

namespace DACLIB {
// in eponymous file
void launch_pvm_slaves( const string &slave_name , int num_slave_procs ,
                        vector<int> &slave_tids );
void launch_pvm_slaves( const string &slave_name , const string &hosts_file ,
                        vector<int> &slave_tids );
}

// ***********************************************************************
int main( int argc , char **argv ) {

  cout << "Triphic v3.1." << endl
       << "Built " << BUILD_TIME << " using OEToolkits version "
       << OEChemGetRelease() << "." << endl << endl
       << "Copyright AstraZeneca 2009, 2014." << endl << endl;

  // a strange quirk of the threading code in OEChem means that by default
  // it accumulates memory in a cache, unless you tell it not to. With
  // triphic running for days on millions of molecules in a typical run,
  // this can cause it to run out of memory which is bad for it and anything
  // else that's running at the time.  This forces it to use the standard
  // system-supplied memory management, which is a bit slower but essential
  // for big virtual screening runs.
  OESystem::OESetMemPoolMode( OESystem::OEMemPoolMode::System );

  // was this launched as a PVM slave?  --IM-A-SLAVE is provided by
  // my launcher. We need to intercept it if it was, and send it straight into
  // listening mode.
  if( 2 == argc && !strcmp( "--IM-A-SLAVE" , argv[1] ) ) {
    slave_event_loop();
    exit( 0 );
  }

  TriphicSettings ts( argc , argv );
  if( !ts ) {
    ts.print_usage( cerr );
    exit( 1 );
  }

  vector<int> slave_tids;
  if( !ts.pvm_hosts_file().empty() ) {
    // launch pvm slaves using details in hosts file
    DACLIB::launch_pvm_slaves( ts.slave_name() , ts.pvm_hosts_file() ,
                               slave_tids );
  } else if( ts.num_slave_procs() > 0 ) {
    // launch slaves using PVM default allocation method
    DACLIB::launch_pvm_slaves( ts.slave_name() , ts.num_slave_procs() ,
                               slave_tids );
    if( slave_tids.empty() ) {
      cerr << "No slaves launched." << endl;
      exit( 1 );
    }
  }

  if( slave_tids.empty() ) {
    try {
      serial_triphic_search( ts );
    } catch( DACLIB::FileReadOpenError &e ) {
      cerr << e.what() << endl;
      exit( 1 );
    } catch( DACLIB::SMARTSDefnError &e ) {
      cerr << e.what() << "From " << ts.not_smarts_file() << endl;
      exit( 1 );
    }
  } else {
    try {
      parallel_triphic_search( ts , slave_tids );
    } catch( DACLIB::FileReadOpenError &e ) {
      cerr << e.what() << endl;
      exit( 1 );
    }
  }

}
