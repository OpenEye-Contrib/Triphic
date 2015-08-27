//
// file send_environment_var_to_mpi_slave.cc
// David Cosgrove
// AstraZeneca
// 28th May 2015

#include <string>
#include <cstdlib>

#include <mpi.h>

using namespace std;

namespace DACLIB {

// in mpi_string_subs.cc
void mpi_send_string( const string &str , int dest_rank );
void mpi_rec_string( int source_rank , string &str );

// **************************************************************************
// send a message to the given slave to set the given environment variable
void send_environment_var_to_mpi_slave( const string &var_name ,
                                        const string &var_val ,
                                        int dest_rank ) {

  string msg1( "Set_Environment" );
  mpi_send_string( msg1 , dest_rank );
  string msg = var_name + "=" + var_val;
  mpi_send_string( msg , dest_rank );

}

// **************************************************************************
// receive a message about and environment variable and set it.
void set_environment_var_from_mpi() {

  string msg;
  mpi_rec_string( 0 , msg );
  char *cmsg = (char *) malloc( msg.length() + 1 ); // putenv uses malloc/free
  strcpy( cmsg , msg.c_str() );
  cmsg[msg.length()] = '\0';
  putenv( cmsg );

}

} // end of namespace DACLIB

