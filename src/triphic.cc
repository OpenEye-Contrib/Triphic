//
// file triphic.cc
// David Cosgrove
// 31st July 2007
//
// This is the main function for the program triphic.
// This is new triphic, v3.0, built to use all the new stuff developed for
// grappel. That's why it's in the same part of the SVN repository - they share
// a lot of code.

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

#include <mpi.h>

using namespace std;
using namespace OEChem;

extern string BUILD_TIME; // in build_time.cc

bool verify_args( TriphicSettings &triphic_settings );

// in triphic_subs.cc
int serial_triphic_search( TriphicSettings &triphic_settings );
void parallel_triphic_search( TriphicSettings &triphic_settings ,
                              int num_slaves );

void slave_event_loop();

// ***********************************************************************
int main( int argc , char **argv ) {

  cout << "Triphic v3.2." << endl
       << "Built " << BUILD_TIME << " using OEToolkits version "
       << OEChemGetRelease() << "." << endl << endl
       << "Copyright AstraZeneca 2009, 2014, 2015, 2016." << endl << endl;

  // sort out the  MPI environment
  MPI_Init( NULL , NULL );
  int world_rank , world_size;
  MPI_Comm_rank( MPI_COMM_WORLD , &world_rank );
  MPI_Comm_size( MPI_COMM_WORLD , &world_size );

#ifdef NOTYET
  cout << "world_rank : " << world_rank << endl
       << "world size : " << world_size << endl;
#endif

  // a strange quirk of the threading code in OEChem means that by default
  // it accumulates memory in a cache, unless you tell it not to. With
  // triphic running for days on millions of molecules in a typical run,
  // this can cause it to run out of memory which is bad for it and anything
  // else that's running at the time.  This forces it to use the standard
  // system-supplied memory management, which is a bit slower but essential
  // for big virtual screening runs.
  OESystem::OESetMemPoolMode( OESystem::OEMemPoolMode::System );

  TriphicSettings ts;
  try {
    ts.parse_args( argc , argv );
  } catch( boost::program_options::error &e ) {
    cerr << e.what() << endl << endl;
    ts.print_usage( cerr );
    MPI_Finalize();
    exit( 1 );
  }

  if( !ts ) {
    MPI_Finalize();
    exit( 1 );
  }

  int err_code = 0;
  if( 1 == world_size ) {
    try {
      err_code = serial_triphic_search( ts );
    } catch( DACLIB::FileReadOpenError &e ) {
      cerr << e.what() << endl;
      MPI_Finalize();
      exit( 1 );
    } catch( DACLIB::SMARTSDefnError &e ) {
      cerr << e.what() << "From " << ts.not_smarts_file() << endl;
      MPI_Finalize();
      exit( 1 );
    }
  } else {
    if( world_rank ) {

      // use process of rank 0 as master
      slave_event_loop();

    } else {

      try {
        parallel_triphic_search( ts , world_size );
      } catch( DACLIB::FileReadOpenError &e ) {
        cerr << e.what() << endl;
        MPI_Finalize();
        exit( 1 );
      }
    }

  }

  MPI_Finalize();
  exit( err_code );

}
