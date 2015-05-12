//
// file VolumeGridPVM.cc
// David Cosgrove
// AstraZeneca
// 12th February 2009
//

#include "VolumeGrid.H"
#include "VolumeGridPVM.H"

#include <pvm3.h>
#include <boost/scoped_ptr.hpp>

namespace DACLIB {

  // **********************************************************************
  void pack_into_pvm_buffer( const VolumeGrid &grid ) {
    
    int i = grid.num_x();
    pvm_pkint( &i , 1 , 1 );
    i = grid.num_y();
    pvm_pkint( &i , 1 , 1 );
    i = grid.num_z();
    pvm_pkint( &i , 1 , 1 );

    float org[3] = { grid.origin()[0] , grid.origin()[1] , grid.origin()[2] };
    pvm_pkfloat( org , 3 , 1 );
    float gs = VolumeGrid::get_grid_spacing();
    pvm_pkfloat( &gs , 1 , 1 );
    int ***gvs = grid.grid();
    pvm_pkint( gvs[0][0] , grid.num_x() * grid.num_y() * grid.num_z() , 1 );

  }

  // **********************************************************************
  void unpack_from_pvm_buffer( VolumeGrid &grid ) {

    int nx , ny , nz;
    float org[3] , gs;

    pvm_upkint( &nx , 1 , 1 );
    pvm_upkint( &ny , 1 , 1 );
    pvm_upkint( &nz , 1 , 1 );
	  
    pvm_upkfloat( org , 3 , 1 );
    pvm_upkfloat( &gs , 1 , 1 );

    int npt = nx * ny * nz;
    boost::scoped_ptr<int> gp( new int[npt] );
    pvm_upkint( gp.get() , npt , 1 );

    VolumeGrid::set_grid_spacing( gs );
    grid = VolumeGrid( org , nx , ny , nz , gp.get() );

    for( int i = 0 ; i < npt ; ++i ) {
      grid.set_grid_value( i , *(gp.get() + i) );
    }

  }

} // end of namespace DACLIB
