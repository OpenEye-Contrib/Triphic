//
// file VolumeGrid.cc
// David Cosgrove
// AstraZeneca
// 12th July 2006
//
// Builds a 3D grid in which molecules can be dropped to calculate included volumes
// and the like.

#include "gtpl_defs.H"
#include "stddefs.H"
#include "AtomSphere.H"
#include "FileExceptions.H"
#include "VolumeGrid.H"

#include <iostream>
#include <sstream>
#include <zlib.h>

#include <boost/shared_ptr.hpp>

using namespace std;
using namespace OEChem;
using namespace OESystem;

namespace DACLIB {

float atomic_num_to_rad( int atomic_num );

// ********************************************************************************
VolumeGrid::VolumeGrid() : solid_volume_( -1.0 ) , surface_volume_( -1.0 ) {

  make_empty_grid();

}

// ********************************************************************************
VolumeGrid::VolumeGrid( float blf[3] , float trb[3] ) : solid_volume_( -1.0 ) ,
  surface_volume_( -1.0 ) {

  // want the grid origin to be an integral number of grid_spacings from (0,0,0)
  // so that 2 grids can be compared even if they are different sizes with different
  // origins, so long as the spacings are the same.
  float blfx = round_origin( blf[0] );
  float blfy = round_origin( blf[1] );
  float blfz = round_origin( blf[2] );
  float trbx = round_origin( trb[0] );
  float trby = round_origin( trb[1] );
  float trbz = round_origin( trb[2] );

  float x_ext = trbx - blfx;
  float y_ext = trby - blfy;
  float z_ext = trbz - blfz;

  origin_[0] = blfx;
  origin_[1] = blfy;
  origin_[2] = blfz;

  num_x_ = 1 + int( x_ext / grid_spacing_ );
  num_y_ = 1 + int( y_ext / grid_spacing_ );
  num_z_ = 1 + int( z_ext / grid_spacing_ );

#ifdef NOTYET
  cout << "VolumeGrid::c'tor" << endl;
  cout << "origin : "; DACLIB::vec_print( origin_ );
  cout << "making grid of size " << num_x_ << " by " << num_y_ << " by "
       << num_z_ << endl;
  cout << "blf : " << blf[0] << " , " << blf[1] << " , " << blf[2]
       << "  trb : " << trb[0] << " , " << trb[1] << " , " << trb[2] << endl;
#endif

  grid_ = DACLIB::make_3d_matrix<int>( num_x_ , num_y_ , num_z_ );

  zero_grid();

}

// ********************************************************************************
// and one that takes an origin and number of grid points. This one assumes
// the origin is correct i.e. it doesn't snap it to an integral number of
// grid spacings from 0,0,0. Used, for example, to rebuild a VolumeGrid
// that has been sent through the ether by PVM or some such.
VolumeGrid::VolumeGrid( float org[3] , int nx , int ny , int nz ,
int *gps ) :
  solid_volume_( -1.0 ) , surface_volume_( -1.0 ) {

  origin_[0] = org[0];
  origin_[1] = org[1];
  origin_[2] = org[2];
  num_x_ = nx;
  num_y_ = ny;
  num_z_ = nz;

  grid_ = DACLIB::make_3d_matrix<int>( num_x_ , num_y_ , num_z_ );
  if( gps ) {
    copy( gps , gps + num_x_ * num_y_ * num_z_ , grid_[0][0] );
  } else {
    zero_grid();
  }

}


// ********************************************************************************
// make a grid that's the same size as all the grids passed in.  Each grid point
// counts the number of times the corresponding point in the input grid is
// CORE or INSIDE_SHELL if grid_type is 1, INSIDE_SHELL or OUTSIDE_SHELL if
// grid_type is 2.
VolumeGrid::VolumeGrid( const vector<VolumeGrid *> &grids , int grid_type ) :
  solid_volume_( -1.0 ) , surface_volume_( -1.0 ) {

  if( grids.empty() ) {
    make_empty_grid();
    return;
  }

  make_grid_enclosing_grids( grids );

  zero_grid();

  for( int i = 0 , is = grids.size() ; i < is ; ++i ) {

    int start1[3] , start2[3] , stop1[3];

    find_grid_subsection( origin_[0] , grids[i]->origin_[0] , num_x_ ,
        grids[i]->num_x_ , start1[0] , stop1[0] , start2[0] );
    find_grid_subsection( origin_[1] , grids[i]->origin_[1] , num_y_ ,
        grids[i]->num_y_ , start1[1] , stop1[1] , start2[1] );
    find_grid_subsection( origin_[2] , grids[i]->origin_[2] , num_z_ ,
        grids[i]->num_z_ , start1[2] , stop1[2] , start2[2] );

    int ***tg = grids[i]->grid_;
    for( int i1 = start1[0] , i2 = start2[0] ; i1 < stop1[0] ; ++i1 , ++i2 ) {
      for( int j1 = start1[1] , j2 = start2[1] ; j1 < stop1[1] ; ++j1 , ++j2 ) {
        for( int k1 = start1[2] , k2 = start2[2] ; k1 < stop1[2] ; ++k1 , ++k2 ) {
          if( 1 == grid_type && tg[i2][j2][k2] > GtplDefs::OUTSIDE_SHELL )
            ++grid_[i1][j1][k1];
          else if( 2 == grid_type ) {
            int pt = tg[i2][j2][k2];
            if( GtplDefs::INSIDE_SHELL == pt || pt == GtplDefs::OUTSIDE_SHELL )
              ++grid_[i1][j1][k1];
          }
        }
      }
    }

  }

}

// ********************************************************************************
VolumeGrid::VolumeGrid( const VolumeGrid &vg ) {

  make_empty_grid();
  copy_data( vg );

}

// ********************************************************************************
VolumeGrid::~VolumeGrid() {

  DACLIB::destroy_3d_matrix<int>( grid_ ); // zeroes grid_ also

}

// ********************************************************************************
VolumeGrid &VolumeGrid::operator=( const VolumeGrid &vg ) {

  if( this == &vg )
    return *this;

  DACLIB::destroy_3d_matrix<int>( grid_ );

  make_empty_grid();
  copy_data( vg );

  return *this;

}

// ********************************************************************************
void VolumeGrid::make_empty_grid() {

  grid_ = 0;
  num_x_ = num_y_ = num_z_ = 0;
  origin_[0] = origin_[1] = origin_[2] = 0.0F;

}

// ********************************************************************************
void VolumeGrid::make_grid_enclosing_grids( const vector<VolumeGrid *> &grids ) {

  float trb[3];
  find_grid_extents( grids , origin_ , trb );

  float x_ext = trb[0] - origin_[0];
  float y_ext = trb[1] - origin_[1];
  float z_ext = trb[2] - origin_[2];
  num_x_ = 1 + int( x_ext / grid_spacing_ );
  num_y_ = 1 + int( y_ext / grid_spacing_ );
  num_z_ = 1 + int( z_ext / grid_spacing_ );

  grid_ = DACLIB::make_3d_matrix<int>( num_x_ , num_y_ , num_z_ );

}

// ********************************************************************************
void VolumeGrid::copy_data( const VolumeGrid &source_grid ) {

  DACLIB::destroy_3d_matrix<int>( grid_ );

  num_x_ = source_grid.num_x_;
  num_y_ = source_grid.num_y_;
  num_z_ = source_grid.num_z_;
  
  origin_[0] = source_grid.origin_[0];
  origin_[1] = source_grid.origin_[1];
  origin_[2] = source_grid.origin_[2];

  solid_volume_ = source_grid.solid_volume_;
  surface_volume_ = source_grid.surface_volume_;

  grid_ = DACLIB::make_3d_matrix<int>( num_x_ , num_y_ , num_z_ );
  copy( source_grid.grid_[0][0] ,
      source_grid.grid_[0][0] + num_x_ * num_y_ * num_z_ , grid_[0][0] );

}

// ********************************************************************************
// find which cells we need to look at for the given atom coords and radius.
void VolumeGrid::find_cells_for_atom( float at_cds[3] , float rad ,
int &x_start , int &x_stop ,
int &y_start , int &y_stop ,
int &z_start , int &z_stop ) const {

  int rad_disp = 1 + int( rint( rad / grid_spacing_ ) );

  // find the nearest grid point to the atom
  int near_x = int( rint( ( at_cds[0] - origin_[0] ) / grid_spacing_ ) );
  int near_y = int( rint( ( at_cds[1] - origin_[1] ) / grid_spacing_ ) );
  int near_z = int( rint( ( at_cds[2] - origin_[2] ) / grid_spacing_ ) );

  x_start = near_x - rad_disp;
  if( x_start < 0 ) x_start = 0;
  x_stop = near_x + rad_disp;
  if( x_stop > num_x_ ) x_stop = num_x_;

  y_start = near_y - rad_disp;
  if( y_start < 0 ) y_start = 0;
  y_stop = near_y + rad_disp;
  if( y_stop > num_y_ ) y_stop = num_y_;

  z_start = near_z - rad_disp;
  if( z_start < 0 ) z_start = 0;
  z_stop = near_z + rad_disp;
  if( z_stop > num_z_ ) z_stop = num_z_;

}

// ********************************************************************************
void VolumeGrid::fill_cells( int x_start , int x_stop , int y_start ,
                             int y_stop , int z_start , int z_stop ,
                             boost::shared_ptr<AtomSphere> &at , float *at_cds ) {

  int at_grid_step = at->step_ratio();

#ifdef NOTYET
  cout << "Dropping in atom " << at->atomic_num() << endl;
  cout << "at_grid_step = " << at_grid_step << endl;
  cout << "points_per_side : " << at->points_per_side()
       << " and per octant : " << at->points_per_octant() << endl;
#endif

  // find the nearest point in at to (xstart,ystart,zstart)
  float at_org_x = at_cds[0] - at->grid_spacing() * float( at->points_per_octant() );
  float at_org_y = at_cds[1] - at->grid_spacing() * float( at->points_per_octant() );
  float at_org_z = at_cds[2] - at->grid_spacing() * float( at->points_per_octant() );

  float grid_cds_x = origin_[0] + float( x_start ) * grid_spacing_;
  float grid_cds_y = origin_[1] + float( y_start ) * grid_spacing_;
  float grid_cds_z = origin_[2] + float( z_start ) * grid_spacing_;

  int x_disp = int( ( grid_cds_x - at_org_x ) / at->grid_spacing() );
  int y_disp = int( ( grid_cds_y - at_org_y ) / at->grid_spacing() );
  int z_disp = int( ( grid_cds_z - at_org_z ) / at->grid_spacing() );

  // make sure we don't go off the end of the AtomSphere grid. Probably
  // this could/should be done by a simple calculation but I tried that
  // and cocked it up.
  int z_stop_new = z_start;
  for( int m = z_start , ma = z_disp ; m < z_stop && ma < at->points_per_side() ; ++m , ma += at_grid_step ) {
    ++z_stop_new;
  }

#ifdef NOTYET
  if( z_stop_new != z_stop ) {
    cout << "AWOOGA : " << z_stop_new << " and " << z_stop << endl;
  }
#endif

  int ***at_grid = at->grid_points();
  for( int k = x_start , ka = x_disp ; k < x_stop ; ++k , ka += at_grid_step ) {
    if( ka < 0 || ka >= at->points_per_side() ) {
      continue;
    }
    for( int l = y_start , la = y_disp ; l < y_stop ; ++l , la += at_grid_step ) {
      if( la < 0 || la >= at->points_per_side() ) {
        continue;
      }
      int m = z_start , ma = z_disp;
      while( ma < 0 ) {
        ma += at_grid_step;
        ++m;
      }
      for( ; m < z_stop_new ; ++m , ma += at_grid_step ) {
        if( at_grid[ka][la][ma] > grid_[k][l][m] ) {
          grid_[k][l][m] = at_grid[ka][la][ma];
        }
      }
    }
  }

}

// ********************************************************************************
// take the number and round it to an integral number of grid_spacing_s.
float VolumeGrid::round_origin( float val ) {

  if( val == float( int( val / grid_spacing_ ) ) * grid_spacing_ )
    return val; // already an integral number

  int nd = int( fabs( val ) / grid_spacing_ ) + 1;
  if( val < 0.0 ) nd = -nd;

  return float( nd ) * grid_spacing_;

}

// ********************************************************************************
void VolumeGrid::calc_volumes() const {

  float cell_vol = grid_spacing_ * grid_spacing_ * grid_spacing_;

  int *pg = grid_[0][0];
  int num_sol = 0 , num_surf = 0;
  for( int i = num_x_ * num_y_ * num_z_ ; i ; --i , ++pg ) {
    // mostly, they'll be 0
    if( !*pg ) {
      continue;
    }
    if( *pg == GtplDefs::CORE || *pg == GtplDefs::INSIDE_SHELL )
      ++num_sol;
    if( *pg == GtplDefs::INSIDE_SHELL || *pg == GtplDefs::OUTSIDE_SHELL )
      ++num_surf;
  }

  solid_volume_ = float( num_sol ) * cell_vol;
  surface_volume_ = float( num_surf ) * cell_vol;

}

// ********************************************************************************
void VolumeGrid::write_to_compressed_file( const string &file ) const {

  ostringstream oss;
  write_to_stream( oss );
  string buf = oss.str();

  gzFile fp = gzopen( file.c_str() , "w" );
  gzwrite( fp , &buf[0] , buf.length() );
  gzclose( fp );

}

// ********************************************************************************
void VolumeGrid::write_to_uncompressed_file( const string &file ) const {

  ofstream ofs( file.c_str() );
  write_to_stream( ofs );

}

// ********************************************************************************
int VolumeGrid::grid( int x , int y , int z ) const {

  if( x >= 0 && x < num_x_ && y >= 0 && y < num_y_ && z >= 0 && z < num_z_ )
    return grid_[x][y][z];
  else
    return 0;

}

// ********************************************************************************
void VolumeGrid::set_grid_value( int x , int y , int z , int new_val ) {

  if( x >= 0 && x < num_x_ && y >= 0 && y < num_y_ && z >= 0 && z < num_z_ )
    grid_[x][y][z] = new_val;

  solid_volume_ = surface_volume_ = -1.0F; // previous value no longer valid

}

// ********************************************************************************
void VolumeGrid::set_grid_value( int point_num , int new_val ) {

  if( point_num >= 0 && point_num < num_x_ * num_y_ * num_z_ )
    *(grid_[0][0] + point_num) = new_val;

  solid_volume_ = surface_volume_ = -1.0F; // previous value no longer valid

}

// ********************************************************************************
void VolumeGrid::zero_grid() {

  // tests suggest that whilst using memset instead of fill is faster with
  // default compiler optimisation, with -O3 they're pretty much the same
  // speed, so leaving it with fill as it's more obvious what's going on.
  fill( grid_[0][0] , grid_[0][0] + num_x_ * num_y_ * num_z_ , 0 );

  solid_volume_ = surface_volume_ = -1.0F; // previous value no longer valid

}

// ********************************************************************************
void VolumeGrid::drop_molecule_in( OEMolBase *mol ,
                                   const vector<float> &atom_rads ) {

  int   x_start , x_stop , y_start , y_stop , z_start , z_stop;

  int j = 0;

  OEIter<OEAtomBase> atom;
  for( atom = mol->GetAtoms() ; atom ; ++atom , ++j ) {

    if( 1 == atom->GetAtomicNum() ) {
      continue;
    }

    float rad = atom_rads.empty() ? -1.0F : atom_rads[j];
    float at_cds[3];
    boost::shared_ptr<AtomSphere> at = get_atom_sphere( mol , atom , rad ,
                                                        at_cds );

    find_cells_for_atom( at_cds , rad , x_start , x_stop , y_start , y_stop ,
                         z_start , z_stop );

    fill_cells( x_start , x_stop , y_start , y_stop , z_start , z_stop ,
                at , at_cds );

  }

  // flag volumes as no longer correct
  solid_volume_ = surface_volume_ = -1.0F;

}

// ********************************************************************************
// drop in spheres of arbitrary position and radius. All points inside
// will be CORE
void VolumeGrid::drop_spheres_in( const vector<BTV> &arb_spheres ) {

  float cds[3] , rad;
  int   x_start , x_stop , y_start , y_stop , z_start , z_stop;

  vector<BTV>::const_iterator p , ps;
  for( p = arb_spheres.begin() , ps = arb_spheres.end() ; p != ps ; ++p ) {
    cds[0] = p->get<0>();
    cds[1] = p->get<1>();
    cds[2] = p->get<2>();
    rad = p->get<3>();
    find_cells_for_atom( cds , rad , x_start , x_stop , y_start , y_stop ,
                         z_start , z_stop );
    boost::shared_ptr<AtomSphere> sp( new AtomSphere( 0 , 4 , grid_spacing_ ,
                                                      rad ) );
    fill_cells( x_start , x_stop , y_start , y_stop , z_start , z_stop ,
                sp , cds );
  }

  // set everything to CORE
  for( int i = 0 , is = num_x_ * num_y_ * num_z_ ; i < is ; ++i )
    if( grid_[0][0][i] != GtplDefs::OUTSIDE )
      grid_[0][0][i] = GtplDefs::CORE;

}

// *****************************************************************************
// calculate the occupied volume;
float VolumeGrid::solid_volume() const {

  if( surface_volume_ < 0.0F )
    calc_volumes();

  return solid_volume_;

}

// *****************************************************************************
float VolumeGrid::surface_volume() const {

  if( surface_volume_ < 0.0F )
    calc_volumes();

  return surface_volume_;

}

// *****************************************************************************
// calculate the volume common to both grids.  Now allows the grids to be
// different sizes with different origins, because the meshes will coincide (that's
// the way they've been set up).
void VolumeGrid::common_volume( const VolumeGrid &target , float &sol_vol ,
                                float &surf_vol ) const {
  
  float cell_vol = grid_spacing_ * grid_spacing_ * grid_spacing_;

  int start1[3] , start2[3] , stop1[3];

  find_grid_subsection( origin_[0] , target.origin_[0] , num_x_ , target.num_x_ ,
      start1[0] , stop1[0] , start2[0] );
  find_grid_subsection( origin_[1] , target.origin_[1] , num_y_ , target.num_y_ ,
      start1[1] , stop1[1] , start2[1] );
  find_grid_subsection( origin_[2] , target.origin_[2] , num_z_ , target.num_z_ ,
      start1[2] , stop1[2] , start2[2] );

  int num_in_sol = 0 , num_in_surf = 0;
  for( int i1 = start1[0] , i2 = start2[0] ; i1 < stop1[0] ; ++i1 , ++i2 ) {
    for( int j1 = start1[1] , j2 = start2[1] ; j1 < stop1[1] ; ++j1 , ++j2 ) {
      for( int k1 = start1[2] , k2 = start2[2] ; k1 < stop1[2] ; ++k1 , ++k2 ) {
        int g1 = grid_[i1][j1][k1];
        int g2 = target.grid_[i2][j2][k2];
        if( ( g1 == GtplDefs::CORE || g1 == GtplDefs::INSIDE_SHELL ) &&
            ( g2 == GtplDefs::CORE || g2 == GtplDefs::INSIDE_SHELL ) )
          ++num_in_sol;
        if( ( g1 == GtplDefs::OUTSIDE_SHELL ||
              g1 == GtplDefs::INSIDE_SHELL ) &&
            ( g2 == GtplDefs::OUTSIDE_SHELL ||
              g2 == GtplDefs::INSIDE_SHELL ) )
          ++num_in_surf;
      }
    }
  }

  sol_vol = float( num_in_sol ) * cell_vol;
  surf_vol = float( num_in_surf ) * cell_vol;

}

// *****************************************************************************
// calculate the shape tanimoto for the two grids, based on CORE and
// INSIDE_SHELL only
float VolumeGrid::shape_tanimoto( const VolumeGrid &target , float &inc_vol ,
                                  float &tot_vol ) const {
  
  int start1[3] , start2[3] , stop1[3];

  find_grid_subsection( origin_[0] , target.origin_[0] , num_x_ , target.num_x_ ,
      start1[0] , stop1[0] , start2[0] );
  find_grid_subsection( origin_[1] , target.origin_[1] , num_y_ , target.num_y_ ,
      start1[1] , stop1[1] , start2[1] );
  find_grid_subsection( origin_[2] , target.origin_[2] , num_z_ , target.num_z_ ,
      start1[2] , stop1[2] , start2[2] );

  int num_inc = 0 , num_tot = 0;
  for( int i1 = start1[0] , i2 = start2[0] ; i1 < stop1[0] ; ++i1 , ++i2 ) {
    for( int j1 = start1[1] , j2 = start2[1] ; j1 < stop1[1] ; ++j1 , ++j2 ) {
      for( int k1 = start1[2] , k2 = start2[2] ; k1 < stop1[2] ; ++k1 , ++k2 ) {
        int g1 = grid_[i1][j1][k1];
        int g2 = target.grid_[i2][j2][k2];
        if( g1 < GtplDefs::INSIDE_SHELL || g2 < GtplDefs::INSIDE_SHELL )
          continue;
        if( g1 || g2 )
          ++num_tot;
        if( g1 && g2 )
          ++num_inc;
      }
    }
  }

  float cell_vol = grid_spacing_ * grid_spacing_ * grid_spacing_;
  inc_vol = float( num_inc ) * cell_vol;
  tot_vol = solid_volume() + target.solid_volume() - inc_vol;

  return inc_vol / tot_vol;
  
}

// ***************************************************************************
// convert the grid to one in which a point flags whether the original value
// was >= trigger_num i.e. convert it from a counts type to a binary type.
void VolumeGrid::counts_to_flag( int trigger_num , int flag_num ) {

  int *g = grid_[0][0] , *gs = g + num_x_* num_y_ * num_z_;

  for( ; g != gs ; ++g )
    *g = *g >= trigger_num ? flag_num : GtplDefs::OUTSIDE;

}

// ***************************************************************************
int VolumeGrid::count_points_of_type( GtplDefs::VOL_MARKER vtype ) {

  return count( grid_[0][0] , grid_[0][0] + num_x_ * num_y_ * num_z_ , vtype );
  
}

// ***************************************************************************
// set all grid points to OUTSIDE except for largest contiguous piece.
void VolumeGrid::keep_largest_contiguous_bit() {

  int biggest_piece_pt[3] = { 0 , 0 , 0 };
  int biggest_piece = 0;

  int done_incr = 10000;
  // find the next starting point
  for( int i = 0 ; i < num_x_ ; ++i ) {
    for( int j = 0 ; j < num_y_ ; ++j ) {
      for( int k = 0 ; k < num_z_ ; ++k ) {
        // we'll be flagging it as done by adding done_incr
        if( grid_[i][j][k] && grid_[i][j][k] < done_incr ) {
          int piece_pt[3] = { i , j , k };
          int piece_size = build_contiguous_piece( piece_pt , done_incr );
          if( piece_size > biggest_piece ) {
            biggest_piece = piece_size;
            biggest_piece_pt[0] = i;
            biggest_piece_pt[1] = j;
            biggest_piece_pt[2] = k;
          }
        }
      }
    }
  }

  // now set them all back to original values
  for( int i = 0 ; i < num_x_ ; ++i ) {
    for( int j = 0 ; j < num_y_ ; ++j ) {
      for( int k = 0 ; k < num_z_ ; ++k ) {
        if( grid_[i][j][k] ) {
          grid_[i][j][k] -= done_incr;
        }
      }
    }
  }

  // re-build the biggest piece
  build_contiguous_piece( biggest_piece_pt , done_incr );

  // now zero bits not in piece and return piece values to original
  for( int i = 0 ; i < num_x_ ; ++i ) {
    for( int j = 0 ; j < num_y_ ; ++j ) {
      for( int k = 0 ; k < num_z_ ; ++k ) {
        if( grid_[i][j][k] < done_incr  ) {
          grid_[i][j][k] = 0;
        } else {
          grid_[i][j][k] -= done_incr;
        }
      }
    }
  }

}

// *****************************************************************************
// find a contiguous piece containing piece_pt, incrementing all members by
// done_incr and returning the size
int VolumeGrid::build_contiguous_piece( int piece_pt[3] , int done_incr ) {

  vector<int> piece_cds;
  bool added_pt = true;
  piece_cds.push_back( piece_pt[0] );
  piece_cds.push_back( piece_pt[1] );
  piece_cds.push_back( piece_pt[2] );
  grid_[piece_pt[0]][piece_pt[1]][piece_pt[2]] += done_incr;

  while( added_pt ) {

    added_pt = false;
    for( unsigned int i = 0 ; i < piece_cds.size() ; i += 3 ) {
      int pt_i = piece_cds[i];
      int pt_j = piece_cds[i+1];
      int pt_k = piece_cds[i+2];
      if( pt_i > 0 && grid_[pt_i-1][pt_j][pt_k] &&
          grid_[pt_i-1][pt_j][pt_k] < done_incr ) {
        added_pt = true;
        piece_cds.push_back( pt_i-1 );
        piece_cds.push_back( pt_j );
        piece_cds.push_back( pt_k );
        grid_[pt_i-1][pt_j][pt_k] += done_incr;
      }
      if( pt_i < num_x_ - 1 && grid_[pt_i+1][pt_j][pt_k] &&
          grid_[pt_i+1][pt_j][pt_k] < done_incr ) {
        added_pt = true;
        piece_cds.push_back( pt_i+1 );
        piece_cds.push_back( pt_j );
        piece_cds.push_back( pt_k );
        grid_[pt_i+1][pt_j][pt_k] += done_incr;
      }

      if( pt_j > 0 && grid_[pt_i][pt_j-1][pt_k] &&
          grid_[pt_i][pt_j-1][pt_k] < done_incr ) {
        added_pt = true;
        piece_cds.push_back( pt_i );
        piece_cds.push_back( pt_j-1 );
        piece_cds.push_back( pt_k );
        grid_[pt_i][pt_j-1][pt_k] += done_incr;
      }
      if( pt_j < num_y_ - 1 && grid_[pt_i][pt_j+1][pt_k] &&
          grid_[pt_i][pt_j+1][pt_k] < done_incr ) {
        added_pt = true;
        piece_cds.push_back( pt_i );
        piece_cds.push_back( pt_j+1 );
        piece_cds.push_back( pt_k );
        grid_[pt_i][pt_j+1][pt_k] += done_incr;
      }

      if( pt_k > 0 && grid_[pt_i][pt_j][pt_k-1] &&
          grid_[pt_i][pt_j][pt_k-1] < done_incr ) {
        added_pt = true;
        piece_cds.push_back( pt_i );
        piece_cds.push_back( pt_j );
        piece_cds.push_back( pt_k-1 );
        grid_[pt_i][pt_j][pt_k-1] += done_incr;
      }
      if( pt_k < num_z_ - 1 && grid_[pt_i][pt_j][pt_k+1] &&
          grid_[pt_i][pt_j][pt_k+1] < done_incr ) {
        added_pt = true;
        piece_cds.push_back( pt_i );
        piece_cds.push_back( pt_j );
        piece_cds.push_back( pt_k+1 );
        grid_[pt_i][pt_j][pt_k+1] += done_incr;
      }

    }

  }

  return piece_cds.size() / 3;

}

// *****************************************************************************
void VolumeGrid::write_to_stream( ostream &os ) const {

  os << num_x_ << " " << num_y_ << " " << num_z_ << endl;
  os << origin_[0] << " " << origin_[1] << " " << origin_[2] << endl;
  os << solid_volume_ << " " << surface_volume_ << " " << grid_spacing_
     << endl;
  for( int i = 0 ; i < num_x_ ; ++i ) {
    for( int j = 0 ; j < num_y_ ; ++j ) {
      for( int k = 0 ; k < num_z_ ; ++k ) {
        os << grid_[i][j][k] << " ";
      }
      os << endl;
    }
  }

}

// *****************************************************************************
void VolumeGrid::write_to_file( const string &file ) const {

  if( file.length() > 3 &&
      file.substr( file.length() - 3 ) == string( ".gz" ) ) {
    write_to_compressed_file( file );
  } else {
    write_to_uncompressed_file( file );
  }

}

// ***************************************************************************
void VolumeGrid::read_from_stream( istream &is ) {

  if( grid_ ) {
    DACLIB::destroy_3d_matrix<int>( grid_ );
    make_empty_grid();
  }

  is >> num_x_ >> num_y_ >> num_z_;
  is >> origin_[0] >> origin_[1] >> origin_[2];
  is >> solid_volume_ >> surface_volume_ >> grid_spacing_;

  grid_ = DACLIB::make_3d_matrix<int>( num_x_ , num_y_ , num_z_ );
  for( int i = 0 ; i < num_x_ ; ++i ) {
    for( int j = 0 ; j < num_y_ ; ++j ) {
      for( int k = 0 ; k < num_z_ ; ++k ) {
        is >> grid_[i][j][k];
      }
    }
  }

}

// ***************************************************************************
void VolumeGrid::read_from_file( const string &file ) {

  gzFile fp = gzopen( file.c_str() , "r" );
  if( !fp ) {
    throw DACLIB::FileReadOpenError( file.c_str() );
  }

  ostringstream oss;
  int next_c = 0;;
  while( -1 != ( next_c = gzgetc( fp ) ) ) {
    oss << char( next_c );
  }

  istringstream iss( oss.str() );
  read_from_stream( iss );

}

// *****************************************************************************
// decide if the given position is inside the vol of given type. Defaults
// to hard sphere volume
bool VolumeGrid::point_in_volume( float pos[3] , int vol_type ) const {

  float test_rad = grid_spacing_ / 2.0;
  int x_start , x_stop , y_start , y_stop , z_start , z_stop;

  find_cells_for_atom( pos , test_rad , x_start , x_stop ,
                       y_start , y_stop , z_start , z_stop );

  int num_in = 0 , num_out = 0;
  for( int i = x_start ; i < x_stop ; ++i ) {
    for( int j = y_start ; j < y_stop ; ++j ) {
      for( int k = z_start ; k < z_stop ; ++k ) {
#ifdef NOTYET
        cout << i << " , " << j << " , " << k << " : " << grid( i , j , k )
             << " : " << vol_type << " : " << int( vol_type & grid( i , j , k ) ) << endl;
#endif
        if( vol_type & grid( i , j , k ) ) {
          ++num_in;
        } else {
          ++num_out;
        }
      }
    }
  }

#ifdef NOTYET
  if( num_in || num_out ) {
    cout << "num_in = " << num_in << "  num_out = " << num_out << endl;
  }
#endif
  return num_in && num_in >= num_out;

}

// *****************************************************************************
void find_grid_extents( const vector<VolumeGrid *> &grids ,
                        float *blf , float *trb ) {

  blf[0] = blf[1] = blf[2] = numeric_limits<float>::max();
  trb[0] = trb[1] = trb[2] = -numeric_limits<float>::max();

  for( int i = 0 , is = grids.size() ; i < is ; ++i ) {
    const float *org = grids[i]->origin();
    if( org[0] < blf[0] ) blf[0] = org[0];
    if( org[1] < blf[1] ) blf[1] = org[1];
    if( org[2] < blf[2] ) blf[2] = org[2];
    float ext = org[0] + grids[i]->num_x() * VolumeGrid::get_grid_spacing();
    if( ext > trb[0] ) trb[0] = ext;
    ext = org[1] + grids[i]->num_y() * VolumeGrid::get_grid_spacing();
    if( ext > trb[1] ) trb[1] = ext;
    ext = org[2] + grids[i]->num_z() * VolumeGrid::get_grid_spacing();
    if( ext > trb[2] ) trb[2] = ext;
  }

}

// ********************************************************************************
// work out which bits of the two grids span a common range in space, and return
// the start and stop values for the common cells
void find_grid_subsection( float org1 , float org2 , int num1 , int num2 ,
                           int &start1 , int &stop1 , int &start2 ) {

  float gs = VolumeGrid::get_grid_spacing();

  if( org1 < org2 ) {
    start1 = int( ( org2 - org1 ) / gs );
    start2 = 0;
  } else {
    start1 = 0;
    start2 = int( ( org1 - org2 ) / gs );
  }

  float end1 = org1 + num1 * gs;
  float end2 = org2 + num2 * gs;
  if( end1 < end2 )
    stop1 = num1;
  else
    stop1 = start1 + num2 - start2;

  // account for numerical hiccups
  stop1 = stop1 > num1 ? num1 : stop1;
  if( start2 + ( stop1 - start1 ) > num2 ) {
    stop1 = num2;
  }

#ifdef NOTYET
  cout << org1 << " , " << num1 << " : " << start1 << " , " << stop1 << endl;
  cout << org2 << " , " << num2 << " : " << start2 << endl;
#endif

}

// ****************************************************************************
void read_grids( const vector<string> &grid_files ,
                 vector<pair<string,DACLIB::VolumeGrid *> > &vol_grids ) {

  for( int i = 0 , is = grid_files.size() ; i < is ; ++i ) {
    vol_grids.push_back( make_pair( grid_files[i] , new DACLIB::VolumeGrid() ) );
    vol_grids[i].second->read_from_file( grid_files[i] );
  }

}

// **************************************************************************
boost::shared_ptr<AtomSphere> get_atom_sphere( OEMolBase *mol ,
                                               OEAtomBase *oeatom ,
                                               float &rad ,
                                               float at_cds[3] ) {

  static const int ATOM_SPHERE_STEP_RATIO = 4;
  static vector<boost::shared_ptr<AtomSphere> > atom_spheres;

  // make sure that, if we already have atom spheres, they are of the correct
  // grid_spacing, since the VolumeGrid grid_spacing can be changed during
  // the program's life
  if( !atom_spheres.empty() ) {
    if( VolumeGrid::get_grid_spacing() / float( ATOM_SPHERE_STEP_RATIO ) !=
        atom_spheres.front()->grid_spacing() ) {
      cout << "Clearing Atom spheres : "
           << VolumeGrid::get_grid_spacing() / float( ATOM_SPHERE_STEP_RATIO )
           << " vs "
           << atom_spheres.front()->grid_spacing() << endl;
      atom_spheres.clear();
    }
  }

  if( rad < 0.0F ) {
    rad = DACLIB::atomic_num_to_rad( oeatom->GetAtomicNum() );
  }
  mol->GetCoords( oeatom , at_cds );

  vector<boost::shared_ptr<AtomSphere> >::iterator p , ps;
  boost::shared_ptr<AtomSphere> at;
  const unsigned int atnum = oeatom->GetAtomicNum();
  for( p = atom_spheres.begin() , ps = atom_spheres.end() ; p != ps ; ++p ) {
    if( atnum == (*p)->atomic_num() ) {
      at = *p;
      break;
    }
  }
  if( !at ) {
    at = boost::shared_ptr<AtomSphere>( new AtomSphere( oeatom->GetAtomicNum() ,
                                                        ATOM_SPHERE_STEP_RATIO ,
                                                        VolumeGrid::get_grid_spacing() ,
                                                        rad ) );
    atom_spheres.push_back( at );
  }

  return at;

}

} // end of namespace DACLIB

