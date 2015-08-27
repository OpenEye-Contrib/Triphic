//
// file send_pphore_sites_by_mpi.cc
// David Cosgrove
// AstraZeneca
// May 29th 2015
//

#include <sstream>
#include <string>
#include <vector>

#include <mpi.h>

#include "BasePPhoreSite.H"

using namespace std;

namespace DACLIB {
// In mpi_string_subs.cc
void mpi_send_string( const string &str , int dest_rank );
void mpi_rec_string( int source_rank , string &str );
}

// in eponymous file
BasePPhoreSite *read_site_from_file( istream &is );


// ****************************************************************************
void send_pphore_site_by_mpi( BasePPhoreSite &pphore_site ,
                              int dest_slave ) {

  ostringstream oss;
  pphore_site.write_to_stream( oss );
  DACLIB::mpi_send_string( oss.str() , dest_slave );

}

// ****************************************************************************
void send_pphore_sites_by_mpi( vector<BasePPhoreSite *> &pphore_sites ,
                               int dest_slave ) {

  unsigned int num_to_send = pphore_sites.size();
  MPI_Send( &num_to_send , 1 , MPI_UNSIGNED , dest_slave , 0 , MPI_COMM_WORLD );
  for( unsigned int i = 0 ; i < num_to_send ; ++i ) {
    send_pphore_site_by_mpi( *pphore_sites[i] , dest_slave );
  }

}

// ****************************************************************************
BasePPhoreSite *receive_pphore_site_by_mpi() {

  string str;
  DACLIB::mpi_rec_string( 0 , str );
  istringstream iss( str );
  return read_site_from_file( iss );

}

// ****************************************************************************
void receive_pphore_sites_by_mpi( vector<BasePPhoreSite *> &pphore_sites ) {

  unsigned int num_to_rec;
  MPI_Recv( &num_to_rec , 1 , MPI_UNSIGNED , 0 , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  for( unsigned int i = 0 ; i < num_to_rec ; ++i ) {
    pphore_sites.push_back( receive_pphore_site_by_mpi() );
  }

}
