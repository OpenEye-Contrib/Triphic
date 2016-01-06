//
// File AtomSphere.cc
// David Cosgrove
// AstraZeneca
// 19th February 2007

#include <iostream>

#include "gtpl_defs.H"
#include "stddefs.H"
#include "AtomSphere.H"

using namespace std;

namespace DACLIB {

float AtomSphere::surface_thickness_ = 0.5F;

// ********************************************************************************
AtomSphere::AtomSphere( int an , int sr , float par_grid_spacing , float radius ) :
  atomic_num_( an ) , step_ratio_( sr ) , radius_( radius ) ,
  grid_( 0 ) {

  float outer_rad = radius_ + surface_thickness_ / 2.0F;
  grid_spacing_ = par_grid_spacing / float( step_ratio_ );
  int num_pts = int( outer_rad / grid_spacing_ );
  if( float( num_pts ) * grid_spacing_ < outer_rad )
    num_pts++;

  points_per_side_ = 2 * num_pts + 1;

  grid_ = DACLIB::make_3d_matrix<int>( points_per_side_ , points_per_side_ ,
				       points_per_side_ );
  fill( grid_[0][0] ,
	grid_[0][0] + points_per_side_ * points_per_side_ * points_per_side_ ,
	GtplDefs::OUTSIDE );

  // fill the points in - only need to do the calculation for 1 octant, the rest
  // will be the same by symmetry.
  float rad_sq = radius_ * radius_;
  float outer_rad_sq = DACLIB::square( outer_rad );
  float inner_rad_sq = DACLIB::square( radius_ - surface_thickness_ / 2.0F );

  for( int i = 0 ; i <= num_pts ; ++i ) {
    float x = float( i ) * grid_spacing_;
    for( int j = 0 ; j <= num_pts ; ++j ) {
      float y = float( j ) * grid_spacing_;
      for( int k = num_pts ; k >= 0 ; --k ) {
	float z = float( k ) * grid_spacing_;
	float dist = x * x + y * y + z * z;
	if( dist <= inner_rad_sq ) {
	  // because it's a sphere, everything down to 0 must now be CORE
	  while( k >= 0 ) {
	    set_symmetric_grid_points( i , j , k , num_pts , GtplDefs::CORE );
	    --k;
	  }
	} else if( dist > inner_rad_sq && dist <= rad_sq ) {
	  set_symmetric_grid_points( i , j , k , num_pts , GtplDefs::INSIDE_SHELL );
	} else if( dist > rad_sq && dist <= outer_rad_sq ) {
	  set_symmetric_grid_points( i , j , k , num_pts , GtplDefs::OUTSIDE_SHELL );
	}
      }
    }
  }

}

// ********************************************************************************
AtomSphere::~AtomSphere() {

  DACLIB::destroy_3d_matrix<int>( grid_ );

}

// ********************************************************************************
int AtomSphere::grid( int x , int y , int z ) const {

  if( x < 0 || x >= points_per_side_ || y < 0 || y >= points_per_side_ ||
      z < 0 || z >= points_per_side_ )
    return 0;
  else
    return grid_[x][y][z];

}

// ********************************************************************************
  void AtomSphere::set_symmetric_grid_points( int i , int j , int k ,
					      int num_pts , int val ) {

    grid_[num_pts + i][num_pts + j][num_pts + k] = val;
    grid_[num_pts + i][num_pts + j][num_pts - k] = val;
    grid_[num_pts + i][num_pts - j][num_pts + k] = val;
    grid_[num_pts + i][num_pts - j][num_pts - k] = val;
    grid_[num_pts - i][num_pts - j][num_pts + k] = val;
    grid_[num_pts - i][num_pts - j][num_pts - k] = val;
    grid_[num_pts - i][num_pts + j][num_pts + k] = val;
    grid_[num_pts - i][num_pts + j][num_pts - k] = val;

  }

} // end of namespace
