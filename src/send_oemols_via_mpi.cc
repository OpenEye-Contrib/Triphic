//
// file send_oemols_via_mpi.cc
// David Cosgrove
// AstraZeneca
// 3rd June 2015
//

#include <vector>

#include <oechem.h>

#include <mpi.h>

using namespace std;
using namespace OEChem;

namespace DACLIB {

void mpi_send_string( const string &str , int dest_rank );
void mpi_rec_string( int source_rank , string &str );

// ****************************************************************************
void send_oemol_via_mpi( OEMolBase *mol , int dest_rank ) {

  unsigned int num_to_send = *mol ? 1 : 0;
  MPI_Send( &num_to_send , 1 , MPI_UNSIGNED , dest_rank , 0 , MPI_COMM_WORLD );
  if( num_to_send ) {
    oemolostream oms;
    oms.SetFormat( OEFormat::OEB );
    oms.openstring();
    oms << *mol;
    mpi_send_string( oms.GetString() , dest_rank );
  }

}

// ****************************************************************************
void rec_oemol_via_mpi( int source_rank , OEMolBase *mol ) {

  unsigned int num_to_rec;
  MPI_Recv( &num_to_rec , 1 , MPI_UNSIGNED , source_rank , 0 , MPI_COMM_WORLD ,
            MPI_STATUS_IGNORE );
  if( num_to_rec ) {
    string mol_conts;
    mpi_rec_string( source_rank , mol_conts );
    oemolistream ims;
    ims.SetFormat( OEFormat::OEB );
    ims.openstring( mol_conts );
    ims >> *mol;
  }

}

} // EO namespace DACLIB

