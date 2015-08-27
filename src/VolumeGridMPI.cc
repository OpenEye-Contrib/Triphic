//
// file VolumeGridMPI.cc
// David Cosgrove
// AstraZeneca
// 1st June 2015
//

#include "VolumeGrid.H"
#include "VolumeGridMPI.H"

#include <mpi.h>
#include <boost/scoped_ptr.hpp>

namespace DACLIB {

// ****************************************************************************
void send_by_mpi( const VolumeGrid &grid , int dest_slave ) {

  int i = grid.num_x();
  MPI_Send( &i , 1 , MPI_INT , dest_slave , 0 , MPI_COMM_WORLD );
  i = grid.num_y();
  MPI_Send( &i , 1 , MPI_INT , dest_slave , 0 , MPI_COMM_WORLD );
  i = grid.num_z();
  MPI_Send( &i , 1 , MPI_INT , dest_slave , 0 , MPI_COMM_WORLD );

  float org[3] = { grid.origin()[0] , grid.origin()[1] , grid.origin()[2] };
  MPI_Send( org , 3 , MPI_FLOAT , dest_slave , 0 , MPI_COMM_WORLD );
  float gs = VolumeGrid::get_grid_spacing();
  MPI_Send( &gs , 1 , MPI_FLOAT , dest_slave , 0 , MPI_COMM_WORLD );
  int ***gvs = grid.grid();
  MPI_Send( gvs[0][0] , grid.num_x() * grid.num_y() * grid.num_z() ,
      MPI_INT , dest_slave , 0 , MPI_COMM_WORLD );

}

// ****************************************************************************
void rec_by_mpi( int source_rank , VolumeGrid &grid ) {

  int nx , ny , nz;
  MPI_Recv( &nx , 1 , MPI_INT , source_rank , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  MPI_Recv( &ny , 1 , MPI_INT , source_rank , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  MPI_Recv( &nz , 1 , MPI_INT , source_rank , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );

  float org[3] , gs;
  MPI_Recv( org , 3 , MPI_FLOAT , source_rank , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
  MPI_Recv( &gs , 1 , MPI_FLOAT , source_rank , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );

  int npt = nx * ny * nz;
  boost::scoped_ptr<int> gp( new int[npt] );
  MPI_Recv( gp.get() , npt , MPI_INT , source_rank , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );

  VolumeGrid::set_grid_spacing( gs );
  grid = VolumeGrid( org , nx , ny , nz , gp.get() );

  for( int i = 0 ; i < npt ; ++i ) {
    grid.set_grid_value( i , *(gp.get() + i) );
  }

}

} // end of namespace DACLIB


