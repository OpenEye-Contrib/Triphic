//
// file launch_pvm_slaves.cc
// David Cosgrove
// AstraZeneca
// 6th August 2007
//
// Function to launch the required number of slave processes. Expects a sensible
// PVM virtual machine to have been created ahead of time by, for example,
// the SGE job starter.

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <pvm3.h>

// Message that a PVM process has died as suggested in the PVM manual
const int TASK_DIED = 11;

using namespace std;

namespace DACLIB {

// ***********************************************************************
// decode the error from pvm_spawn and report to stderr
void report_pvm_spawn_error( int err_code ) {

  switch( err_code ) {
  case PvmBadParam :
    cerr << "Giving an invalid argument value" << endl;
    break;
  case PvmNoHost :
    cerr << "Specified host is not in the virtual machine." << endl;
    break;
  case PvmNoFile :
    cerr << "Specified executable cannot be found." << endl;
    break;
  case PvmNoMem :
    cerr << "Malloc failed. Not enough memory on host." << endl;
    break;
  case PvmSysErr :
    cerr << "pvmd not responding." << endl;
    break;
  case PvmOutOfRes :
    cerr << "Out of resources." << endl;
    break;
  default :
    cerr << "Got an error not in the v3.4 manual." << endl;
    break;
  }
}

// **************************************************************************
// launch PVM slaves using the details in the file, which should be
// the host followed by the number of copies on that host
void launch_pvm_slaves( const string &slave_name , const string &hosts_file ,
                        vector<int> &slave_tids ) {

  ifstream ifs( hosts_file.c_str() );

  if( !ifs || !ifs.good() ) {
    cerr << "Couldn't read PVM hosts file." << endl;
    exit( 1 );
  }

  vector<string> node_names;
  vector<int> num_slave_procs;
  int tot_num_procs = 0;
  while( true ) {
    string nextline;
    getline( ifs , nextline );
    if( !ifs.good() || ifs.eof() ) {
      break;
    }
    if( nextline.empty() ) {
      continue;
    }
    istringstream iss( nextline );
    string host;
    int num_procs;
    iss >> host >> num_procs;
    if( iss.fail() ) {
      continue;
    }
    node_names.push_back( host );
    num_slave_procs.push_back( num_procs );
    tot_num_procs += num_procs;
  }

  slave_tids = vector<int>( tot_num_procs , -1 );
  static char **slave_args = 0;
  if( !slave_args ) {
    slave_args = new char*[2];
    slave_args[0] = new char[13];
    strcpy( slave_args[0] , "--IM-A-SLAVE" );
    slave_args[1] = 0;
  }

  int procs_so_far = 0;
  for( int i = 0 , is = node_names.size() ; i < is ; ++i ) {
    int num_spawned = pvm_spawn( const_cast<char *>( slave_name.c_str() ) ,
                                 slave_args , PvmTaskDefault ,
                                 const_cast<char *>( node_names[i].c_str() ) ,
                                 num_slave_procs[i] , &slave_tids[0] + procs_so_far );
    if( num_spawned < num_slave_procs[i] ) {
      cerr << "Successfully spawned " << num_spawned << " copies of "
           << slave_name << " on " << node_names[i] << " but wanted "
           << num_slave_procs[i] << " so failing." << endl;
      for( int j = num_spawned ; j < num_slave_procs[j] ; ++j ) {
        report_pvm_spawn_error( slave_tids[j] );
      }
      exit( 1 );
    }

    cerr << "Successfully spawned " << num_spawned << " copies of "
         << slave_name << " on " << node_names[i] << endl;
    procs_so_far += num_spawned;
  }

  // we want to hear from any slaves that die
  pvm_notify( PvmTaskExit , TASK_DIED , slave_tids.size() , &slave_tids[0] );

}

// **************************************************************************
void launch_pvm_slaves( const string &slave_name , int num_slave_procs ,
                        vector<int> &slave_tids ) {

  static char **slave_args = 0;
  if( !slave_args ) {
    slave_args = new char*[2];
    slave_args[0] = new char[13];
    strcpy( slave_args[0] , "--IM-A-SLAVE" );
    slave_args[1] = 0;
  }

  slave_tids = vector<int>( num_slave_procs , -1 );
  int num_spawned = pvm_spawn( const_cast<char *>( slave_name.c_str() ) ,
                               slave_args , PvmTaskDefault ,
                               static_cast<char *>( 0 ) ,
                               num_slave_procs , &slave_tids[0] );
  if( num_spawned < 1 ) {
    cerr << "Successfully spawned " << num_spawned << " copies of "
         << slave_name << endl;
    exit( 1 );
  }

  for( int i = num_spawned ; i < num_slave_procs ; ++i )
    report_pvm_spawn_error( slave_tids[i] );
  slave_tids.erase( slave_tids.begin() + num_spawned , slave_tids.end() );
  cerr << "Successfully spawned " << num_spawned << " copies of "
       << slave_name << endl;

  // we want to hear from any slaves that die
  pvm_notify( PvmTaskExit , TASK_DIED , num_slave_procs , &slave_tids[0] );

}

} // end of namespace DACLIB
