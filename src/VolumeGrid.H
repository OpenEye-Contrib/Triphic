//
// file VolumeGrid.H
// David Cosgrove
// AstraZeneca
// 12th July 2006
//
// Builds a 3D grid in which molecules can be dropped to calculate included volumes
// and the like.

#ifndef DAC_VOLUME_GRID__
#define DAC_VOLUME_GRID__

#include "gtpl_defs.H"

#include <oechem.h>

#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>

#include <iosfwd>

// ********************************************************************************

class AtomSphere;
typedef boost::tuple<float,float,float,float> BTV; // 1st 3 are coords, last is rad

namespace DACLIB {

class VolumeGrid {

protected :

  int     ***grid_;
  
  int     num_x_ , num_y_ , num_z_;
  float   origin_[3];
  mutable float solid_volume_ , surface_volume_;

  static float grid_spacing_;

  void make_empty_grid();

  void copy_data( const VolumeGrid &source_grid );

  void make_grid_enclosing_grids( const std::vector<VolumeGrid *> &grids );

  // find which cells we need to look at for the given atom coords and radius.
  void find_cells_for_atom( float at_cds[3] , float rad ,
  int &x_start , int &x_stop ,
  int &y_start , int &y_stop ,
  int &z_start , int &z_stop ) const;
  // fill the given cells
  virtual void fill_cells( int x_start , int x_stop , int y_start ,
                           int y_stop , int z_start , int z_stop ,
                           boost::shared_ptr<AtomSphere > &at , float *at_cds );

  // take the number and round it to an integral number of grid_spacing_s.
  float round_origin( float val );

  void calc_volumes() const;

  void write_to_compressed_file( const std::string &file ) const;
  void write_to_uncompressed_file( const std::string &file ) const;

public :

  VolumeGrid();
  // c'tor that takes just the bottom, left, front and top, right, back points
  // and works out the sizes from the grid spacing.
  VolumeGrid( float blf[3] , float trb[3] );
  // and one that takes an origin and number of grid points. This one assumes
  // the origin is correct i.e. it doesn't snap it to and integral number of
  // grid spacings from 0,0,0. Used, for example, to rebuild a VolumeGrid
  // that has been sent through the ether by PVM or some such.
  VolumeGrid( float org[3] , int nx , int ny , int nz , int *gps = 0 );

  // make a grid that's the same size as all the grids passed in. Each point
  // in the new grid counts the number of times the corresponding point in the
  // input grid is CORE or INSIDE_SHELL if grid_type is 1, INSIDE_SHELL or
  // OUTSIDE_SHELL if grid_type is 2.
  VolumeGrid( const std::vector<VolumeGrid *> &grids , int grid_type );
  VolumeGrid( const VolumeGrid &vg );

  virtual ~VolumeGrid();
  VolumeGrid &operator=( const VolumeGrid &vg );

  static void set_grid_spacing( float gs ) { grid_spacing_ = gs; }
  static float get_grid_spacing() { return grid_spacing_; }

  int num_x() const { return num_x_; }
  int num_y() const { return num_y_; }
  int num_z() const { return num_z_; }
  const float *origin() const { return origin_; }
  int ***grid() const { return grid_; } /* compiler won't let this be const int ***
             for reasons that I don't understand. */
  int grid( int x , int y , int z ) const;
  void set_grid_value( int x , int y , int z , int new_val );
  void set_grid_value( int point_num , int new_val );

  void zero_grid();
  void drop_molecule_in( OEChem::OEMolBase *mol ,
                         const std::vector<float> &atom_rads );
  // drop in spheres of arbitrary position and radius. All points inside
  // will be CORE
  void drop_spheres_in( const std::vector<BTV> &arb_spheres );
  float solid_volume() const;
  float surface_volume() const;

  // return solid and surface volumes in common
  void common_volume( const VolumeGrid &target , float &sol_vol ,
                      float &surf_vol ) const;

  // this is always on the solid volume
  float shape_tanimoto( const VolumeGrid &target , float &inc_vol ,
                        float &tot_vol ) const;

  // convert the grid to one in which a point flags whether the original value
  // was >= trigger_num i.e. convert it from a counts type to a binary type. All
  // points are made to type flag_num if triggered.
  void counts_to_flag( int trigger_num , int flag_num );

  int count_points_of_type( GtplDefs::VOL_MARKER vtype );

  // set all grid points to OUTSIDE except for largest contiguous piece.
  void keep_largest_contiguous_bit();
  // find a contiguous piece containing piece_pt, incrementing all members by
  // done_incr and returning the size
  int build_contiguous_piece( int piece_pt[3] , int done_incr );

  void write_to_stream( std::ostream &os ) const;
  // file will be compressed with zlib if name ends ".gz"
  void write_to_file( const std::string &file ) const;
  void read_from_stream( std::istream &is );
  // file may be compressed with zlib
  void read_from_file( const std::string &file );

  // decide if the given position is inside the vol of given type. Defaults
  // to hard sphere volume
  bool point_in_volume( float pos[3] , int vol_type = GtplDefs::CORE | GtplDefs::INSIDE_SHELL ) const;

};

void find_grid_extents( const std::vector<VolumeGrid *> &grids ,
                        float *blf , float *trb );

// work out which bits of the two grids span a common range in space, and
// return the start and stop values for the common cells
void find_grid_subsection( float org1 , float org2 , int num1 , int num2 ,
                           int &start1 , int &stop1 , int &start2 );

// read a bunch of grids into a vector, one per file
void read_grids( const std::vector<std::string> &grid_files ,
                 std::vector<std::pair<std::string,VolumeGrid *> > &vol_grids );

boost::shared_ptr<AtomSphere> get_atom_sphere( OEChem::OEMolBase *mol ,
                                               OEChem::OEAtomBase *oeatom ,
                                               float &rad , float at_cds[3] );

} // end of namespace

#endif
