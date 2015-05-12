//
// file AtomSphere.H
// David Cosgrove
// AstraZeneca
// 19th February 2007
//
// Makes an atom sphere centred on the origin, to the given grid spacing. The
// step ratio is how many points there are per grid spacing in the VolumeGrid
// that this will be used with.

#ifndef DAC_ATOM_SPHERE
#define DAC_ATOM_SPHERE

#include <string>

// *****************************************************************************

namespace DACLIB {

class AtomSphere {

private :

  static float surface_thickness_;

  unsigned int   atomic_num_;
  float grid_spacing_; // for this grid, not the VolumeGrid
  int   step_ratio_;
  int   points_per_side_;
  float radius_;

  int   ***grid_;

  void set_symmetric_grid_points( int i , int j , int k , int num_pts ,
				  int val );

public :

  AtomSphere( int atomic_num , int step_ratio , float par_grid_spacing ,
	      float radius );
  ~AtomSphere();

  unsigned int atomic_num() const { return atomic_num_; }
  int step_ratio() const { return step_ratio_; }
  int points_per_side() const { return points_per_side_; }
  int points_per_octant() const { return points_per_side_ / 2; }
  float grid_spacing() const { return grid_spacing_; }
  int grid( int x , int y , int z ) const;
  static float surface_thickness() { return surface_thickness_; }
  int ***grid_points() const { return grid_; }

};

} // end of namespace

#endif
