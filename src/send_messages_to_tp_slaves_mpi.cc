//
// file send_messages_to_tp_slaves_mpi.cc
// David Cosgrove
// AstraZeneca
// 29th May 2015
// Various functions for sending messages to the triphic and plurality slaves
// using OpenMPI

#include <iostream>
#include <string>
#include <vector>

#include <mpi.h>

#include <sys/param.h>
#include <cstdlib>

using namespace std;

namespace DACLIB {
void send_environment_var_to_mpi_slave( const string &var_name ,
                                        const string &var_val ,
                                        int dest_rank );
// In mpi_string_subs.cc
void mpi_send_string( const string &str , int dest_rank );
void mpi_rec_string( int source_rank , string &str );
}

// ********************************************************************
// spell from t'Interweb
string getcwd(){

  char *buffer = new char[MAXPATHLEN];
  getcwd(buffer,MAXPATHLEN);
  if(buffer != NULL){
    string ret(buffer);
    delete[] buffer;
    return ret;
  } else {
    return string();
  }

}

// ********************************************************************
void send_openeye_license_to_slaves( int num_slaves ) {

  if( getenv( "OE_LICENSE" ) ) {
    string oe_lic( getenv( "OE_LICENSE" ) );
    for( int i = 1 ; i < num_slaves ; ++i ) {
      DACLIB::send_environment_var_to_mpi_slave( string( "OE_LICENSE" ) ,
                                                 oe_lic , i );
    }
  }

}

// ********************************************************************
void send_cwd_to_slaves( int num_slaves ) {

  string cwd = getcwd();
  if( cwd.length() ) {
    string msg( "New_CWD" );
    for( int i = 1 ; i < num_slaves ; ++i ) {
      DACLIB::mpi_send_string( msg , i );
      DACLIB::mpi_send_string( cwd , i );
    }
  }

}

// ********************************************************************
void send_finished_messages( int num_slaves ) {

  static const string msg( "Finished" );
  for( int i = 1 ; i < num_slaves ; ++i ) {
    DACLIB::mpi_send_string( msg , i );
  }

}

// ********************************************************************
void receive_new_cwd() {

  string new_cwd;
  DACLIB::mpi_rec_string( 0 , new_cwd );
  chdir( new_cwd.c_str() );

}

// **************************************************************************
void send_progress_to_master( const string &progress_report ) {

  const static string msg_header( "Progress Report" );
  DACLIB::mpi_send_string( msg_header , 0 );
  DACLIB::mpi_send_string( progress_report , 0 );

}

// ********************************************************************
void send_database_details_to_slaves( int world_size ) {

  static const string db_steps_msg( "Database_Steps" );
  int step_size = world_size - 1;
  for( int i = 0 ; i < step_size ; ++i ) {
    DACLIB::mpi_send_string( db_steps_msg , i + 1 );
    // send the start molecule for this slave, and the step size
    MPI_Send( &i , 1 , MPI_INT , i + 1 , 0 , MPI_COMM_WORLD );
    MPI_Send( &step_size , 1 , MPI_INT , i + 1 , 0 , MPI_COMM_WORLD );
  }

}

// ********************************************************************
void receive_database_details( int &db_start , int &db_step ) {

  MPI_Recv( &db_start , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD ,
            MPI_STATUS_IGNORE );
  MPI_Recv( &db_step , 1 , MPI_INT , 0 , 0 , MPI_COMM_WORLD ,
            MPI_STATUS_IGNORE );

}
