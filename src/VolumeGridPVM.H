//
// file VolumeGridPVM.H
// David Cosgrove
// AstraZeneca
// 12th February 2009
//
// This file declares functions for putting VolumeGrids into a previously
// defined pvm buffer.

#ifndef DAC_VOLUME_GRID_PVM
#define DAC_VOLUME_GRID_PVM

namespace DACLIB {

  class VolumeGrid;
  void pack_into_pvm_buffer( const DACLIB::VolumeGrid &grid );
  void unpack_from_pvm_buffer( DACLIB::VolumeGrid &grid );

}

#endif
