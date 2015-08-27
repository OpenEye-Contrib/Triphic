//
// file VolumeGridMPI.H
// David Cosgrove
// AstraZeneca
// 1st June 2015
//
// This file declares functions for sending and receiving VolumeGrids
// via MPI

#ifndef VOLUMEGRIDMPI_H
#define VOLUMEGRIDMPI_H

namespace DACLIB {

  class VolumeGrid;
  void send_by_mpi( const DACLIB::VolumeGrid &grid , int dest_slave );
  void rec_by_mpi( int source_rank , DACLIB::VolumeGrid &grid );

}

#endif // VOLUMEGRIDMPI_H
