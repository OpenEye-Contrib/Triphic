//
// file mpi_string_subs.cc
// David Cosgrove
// AstraZeneca
// 28th May 2015
//
// This file contains stuff for passing STL strings with MPI

#include <string>
#include <vector>

#include <mpi.h>

namespace DACLIB {

// ****************************************************************************
void mpi_send_string( const std::string &str , int dest_rank ) {

  MPI_Send( const_cast<char *>( &str[0] ) , str.length() , MPI_CHAR ,
      dest_rank , 0 , MPI_COMM_WORLD );

}

// ****************************************************************************
void mpi_send_strings_vector( const std::vector<std::string> &strs ,
                              int dest_rank ) {

  unsigned int num_to_send = strs.size();
  MPI_Send( &num_to_send , 1 , MPI_UNSIGNED , dest_rank , 0 , MPI_COMM_WORLD );
  for( unsigned int i = 0 ; i < num_to_send ; ++i ) {
    mpi_send_string( strs[i] , dest_rank );
  }

}

// ****************************************************************************
void mpi_rec_string( int source_rank , std::string &str ) {

  int      msg_len;
  MPI_Status status;
  MPI_Probe( source_rank , 0 , MPI_COMM_WORLD , &status );
  MPI_Get_count( &status , MPI_CHAR , &msg_len );

  str.resize( msg_len );

  MPI_Recv( &str[0] , msg_len , MPI_CHAR , source_rank , 0 , MPI_COMM_WORLD ,
      MPI_STATUS_IGNORE );

}

// ****************************************************************************
void mpi_rec_strings_vector( int source_rank , std::vector<std::string> &strs ) {

  unsigned int num_to_rec;
  MPI_Recv( &num_to_rec , 1 , MPI_UNSIGNED , 0 , 0 , MPI_COMM_WORLD ,
            MPI_STATUS_IGNORE );
  strs = std::vector<std::string>( num_to_rec );
  for( unsigned int i = 0 ; i < num_to_rec ; ++i ) {
    mpi_rec_string( source_rank , strs[i] );
  }

}

}
