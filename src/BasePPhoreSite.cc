//
// file BasePPhoreSite.cc
// David Cosgrove
// AstraZeneca
// 10th October 2006
//
// Implementation of BasePPhoreSite.

#include <algorithm>
#include <numeric>

#include <boost/bind.hpp>

#include "stddefs.H"
#include "BasePPhoreSite.H"
#include "OverlayTrans.H"

namespace DACLIB {
void rotate_about_axis( const double axis[3] , const double angle ,
double vec[3] );

double angle_about_axis( const double axis[3] , const double vec1[3] ,
const double vec2[3] );
float torsion( const float *cds1 , const float *cds2 , const float *cds3 ,
               const float *cds4 , bool degs );
}

using namespace std;

// *******************************************************************************
BasePPhoreSite::BasePPhoreSite() {

  reset_data();

}

// *******************************************************************************
// tc and ts are the type code and type string respectively.
BasePPhoreSite::BasePPhoreSite( const double cds[3] , const double dir[3] , int tc ,
const string &ts , const string &lbl , bool use_dir ) {

  label_ = lbl;
  num_dirs_ = use_dir ? 1 : 0;
  selected_ = false;

  type_code_ = tc;
  type_string_ = ts;

  cds_[0] = cds[0]; cds_[1] = cds[1]; cds_[2] = cds[2];

  if( use_dir ) {
    dir_[0] = dir[0]; dir_[1] = dir[1]; dir_[2] = dir[2];
  } else {
    dir_[0] = dir_[1] = dir_[2] = 0.0;
  }

  dir_moves_ = GtplDefs::NONE;

  num_virt_sites_ = 0;
  active_virt_site_ = -1;

}

// *******************************************************************************
BasePPhoreSite::~BasePPhoreSite() {

}

// ****************************************************************************
void BasePPhoreSite::reset_data() {

  label_ = "";

  num_dirs_ = 0;
  selected_ = false;

  type_code_ = -1;
  type_string_ = "";

  cds_[0] = cds_[1] = cds_[2] = 0.0;
  dir_[0] = dir_[1] = dir_[2] = 0.0;
  dir_moves_ = GtplDefs::NONE;

  num_virt_sites_ = 0;
  active_virt_site_ = -1;

}

// ****************************************************************************
void BasePPhoreSite::copy_data( const BasePPhoreSite &c ) {

  label_ = c.label_;

  type_code_ = c.type_code_;
  type_string_ = c.type_string_;

  cds_[0] = c.cds_[0]; cds_[1] = c.cds_[1]; cds_[2] = c.cds_[2];

  num_dirs_ = c.num_dirs_;
  if( c.num_dirs_ ) {
    copy( c.dir_ , c.dir_ + num_dirs_ * 3 , dir_ );
    copy( c.dir_types_ , c.dir_types_ + num_dirs_ , dir_types_ );
  }

  selected_ = c.selected_;

  dir_moves_ = c.dir_moves_;
  twiddle_axis_[0] = c.twiddle_axis_[0];
  twiddle_axis_[1] = c.twiddle_axis_[1];
  twiddle_axis_[2] = c.twiddle_axis_[2];

  num_virt_sites_ = c.num_virt_sites_;
  copy( c.virt_sites_ , c.virt_sites_ + 3 * num_virt_sites_ , virt_sites_ );
  active_virt_site_ = c.active_virt_site_;

}

// ****************************************************************************
void BasePPhoreSite::set_coords( double *new_cds ) {

  cds_[0] = new_cds[0];
  cds_[1] = new_cds[1];
  cds_[2] = new_cds[2];

}

// ****************************************************************************
void BasePPhoreSite::set_coords( float *new_cds ) {

  cds_[0] = new_cds[0];
  cds_[1] = new_cds[1];
  cds_[2] = new_cds[2];

}

// ****************************************************************************
int BasePPhoreSite::get_num_dirs( GtplDefs::DIRS_TYPE dirs_type ) const {

  int ret_num = 0;
  for( int i = 0 ; i < num_dirs_ ; ++i ) {
    if( dirs_type == dir_types_[i] ) {
      ++ret_num;
    }
  }

  return ret_num;

}

// ****************************************************************************
const double *BasePPhoreSite::direction( int dir_num ) const {

  if( dir_num < 0 || dir_num >= num_dirs_ ) {
    return 0;
  } else {
    return dir_ + 3 * dir_num;
  }

}

// ****************************************************************************
GtplDefs::DIRS_TYPE BasePPhoreSite::direction_type( int dir_num ) const {

  if( dir_num < 0 || dir_num >= num_dirs_ ) {
    return GtplDefs::UNKNOWN;
  } else {
    return dir_types_[dir_num];
  }

}

// ****************************************************************************
void BasePPhoreSite::set_direction( const double *new_dir ,
                                    GtplDefs::DIRS_TYPE t , int dir_num  ) {

  if( -1 == dir_num && num_dirs_ < 3 ) {
    dir_[3 * num_dirs_] = new_dir[0];
    dir_[3 * num_dirs_ + 1] = new_dir[1];
    dir_[3 * num_dirs_ + 2] = new_dir[2];
    dir_types_[num_dirs_] = t;
    ++num_dirs_;
  }
  if( -1 != dir_num && dir_num < 3 ) {
    dir_[3 * dir_num] = new_dir[0];
    dir_[3 * dir_num + 1] = new_dir[1];
    dir_[3 * dir_num + 2] = new_dir[2];
    dir_types_[dir_num] = t;
  }

}

// ****************************************************************************
void BasePPhoreSite::set_direction( const float *new_dir ,
                                    GtplDefs::DIRS_TYPE t , int dir_num ) {

  if( -1 == dir_num && num_dirs_ < 3 ) {
    dir_[3 * num_dirs_] = new_dir[0];
    dir_[3 * num_dirs_ + 1] = new_dir[1];
    dir_[3 * num_dirs_ + 2] = new_dir[2];
    dir_types_[num_dirs_] = t;
    ++num_dirs_;
  }
  if( -1 != dir_num && dir_num < 3 ) {
    dir_[3 * dir_num] = new_dir[0];
    dir_[3 * dir_num + 1] = new_dir[1];
    dir_[3 * dir_num + 2] = new_dir[2];
    dir_types_[dir_num] = t;
  }

}

// ****************************************************************************
float BasePPhoreSite::square_distance( const double *cds ) const {

  return ( DACLIB::square( cds[0] - cds_[0] ) +
	   DACLIB::square( cds[1] - cds_[1] ) +
	   DACLIB::square( cds[2] - cds_[2] ) );

}

// ****************************************************************************
void BasePPhoreSite::translate( float x_trans , float y_trans , float z_trans ) {

  cds_[0] += x_trans;
  cds_[1] += y_trans;
  cds_[2] += z_trans;

}

// ****************************************************************************
void BasePPhoreSite::rotate( const float rot[3][3] ) {

  float cds[3];
  cds[0] = rot[0][0] * cds_[0] + rot[0][1] * cds_[1] + rot[0][2] * cds_[2];
  cds[1] = rot[1][0] * cds_[0] + rot[1][1] * cds_[1] + rot[1][2] * cds_[2];
  cds[2] = rot[2][0] * cds_[0] + rot[2][1] * cds_[1] + rot[2][2] * cds_[2];

  cds_[0] = cds[0];
  cds_[1] = cds[1];
  cds_[2] = cds[2];

  cds[0] = rot[0][0] * twiddle_axis_[0] + rot[0][1] * twiddle_axis_[1] + rot[0][2] * twiddle_axis_[2];
  cds[1] = rot[1][0] * twiddle_axis_[0] + rot[1][1] * twiddle_axis_[1] + rot[1][2] * twiddle_axis_[2];
  cds[2] = rot[2][0] * twiddle_axis_[0] + rot[2][1] * twiddle_axis_[1] + rot[2][2] * twiddle_axis_[2];

  twiddle_axis_[0] = cds[0];
  twiddle_axis_[1] = cds[1];
  twiddle_axis_[2] = cds[2];

  for( int i = 0 ; i < num_dirs_ ; ++i ) {
    cds[0] =
        rot[0][0] * dir_[3*i] + rot[0][1] * dir_[3*i+1] + rot[0][2] * dir_[3*i+2];
    cds[1] =
        rot[1][0] * dir_[3*i] + rot[1][1] * dir_[3*i+1] + rot[1][2] * dir_[3*i+2];
    cds[2] =
        rot[2][0] * dir_[3*i] + rot[2][1] * dir_[3*i+1] + rot[2][2] * dir_[3*i+2];

    dir_[3*i] = cds[0];
    dir_[3*i+1] = cds[1];
    dir_[3*i+2] = cds[2];
  }

  for( int i = 0 ; i < num_virt_sites_ ; ++i ) {
    cds[0] =
        rot[0][0] * virt_sites_[3*i] + rot[0][1] * virt_sites_[3*i+1] + rot[0][2] * virt_sites_[3*i+2];
    cds[1] =
        rot[1][0] * virt_sites_[3*i] + rot[1][1] * virt_sites_[3*i+1] + rot[1][2] * virt_sites_[3*i+2];
    cds[2] =
        rot[2][0] * virt_sites_[3*i] + rot[2][1] * virt_sites_[3*i+1] + rot[2][2] * virt_sites_[3*i+2];

    virt_sites_[3*i] = cds[0];
    virt_sites_[3*i+1] = cds[1];
    virt_sites_[3*i+2] = cds[2];
  }

}

// ****************************************************************************
// find the order of directions in site that has the most of the given
// type lined up with directions in this. Dir_order maps dirs_ vectors in
// site onto dirs_ vectors in this. Returns the sum of the dot-products for
// the order found.
float BasePPhoreSite::best_dir_alignments( BasePPhoreSite &site ,
                                           GtplDefs::DIRS_TYPE dirs_type ,
                                           bool twiddle_if_poss ,
                                           vector<int> &dir_order ) {
  
  // needs to be done with num_dirs_ >= site.num_dirs_
  if( site.num_dirs_ > num_dirs_ )
    return site.best_dir_alignments( *this , dirs_type , twiddle_if_poss ,
                                     dir_order );

  vector<int> combs;
  for( int i = 0 ; i < site.num_dirs_ ; ++i )
    combs.push_back( i );

  float best_dot = -1.0;
  dir_order = combs;
  do {
    // line up on the first direction if appropriate
    if( twiddle_if_poss && site.get_twiddlable() )
      site.twiddle( dir_ , combs[0] );

    float cum_dot = 0.0F;
    for( int i = 0 ; i < site.num_dirs_ ; ++i ) {
      if( dir_types_[i] == dirs_type )
        cum_dot = DACLIB::dot_product( dir_ + 3 * i , site.dir_ + 3 * combs[i] );
    }
    if( cum_dot > best_dot ) {
      best_dot = cum_dot;
      dir_order = combs;
    }
  } while ( next_permutation( combs.begin() , combs.end() ) );

  return best_dot;

}

// ****************************************************************************
// put the direction vectors into the new order.
void BasePPhoreSite::reorder_dirs( const vector<int> &new_order ) {

  double dirs2[9];
  GtplDefs::DIRS_TYPE dt[3];

  for( int i = 0 ; i < std::min( num_dirs_ , int( new_order.size() ) ) ; ++i ) {
    dirs2[3 * i] = dir_[3 * new_order[i]];
    dirs2[3 * i + 1] = dir_[3 * new_order[i] + 1];
    dirs2[3 * i + 2] = dir_[3 * new_order[i] + 2];
    dt[i] = dir_types_[new_order[i]];
  }

  for( int i = 0 ; i < get_num_dirs() ; ++i ) {
    dir_[3 * i] = dirs2[3 * i];
    dir_[3 * i + 1] = dirs2[3 * i + 1];
    dir_[3 * i + 2] = dirs2[3 * i + 2];
    dir_types_[i] = dt[i];
  }

}

// ****************************************************************************
// twiddle the given direction to point along the given vector dir
// if it's twiddlable leaving it unchanged if it isn't
void BasePPhoreSite::twiddle( const double *dir , int dir_num ) {

  if( dir_moves_ != GtplDefs::FREE || dir_num < 0 || dir_num >= num_dirs_ ) {
    return;
  }

  double in_dir[3] = { dir[0] , dir[1] , dir[2] };
  double this_dir[3] = { dir_[3 * dir_num] , dir_[3 * dir_num + 1] ,
                         dir_[3 * dir_num + 2] };

  double angle = DACLIB::angle_about_axis( twiddle_axis_ , in_dir , this_dir );
  DACLIB::rotate_about_axis( twiddle_axis_ , -angle , this_dir );
  double angle_to_use = -angle;
  // sometimes the angle comes back with the wrong sign, for reasons which
  // I don't understand and haven't delved into deeply. If it's right, the
  // rotation will put the vectors with essentially a zero angle.  If not, try
  // it the other way.
  double new_angle = DACLIB::angle_about_axis( twiddle_axis_ , in_dir , this_dir );
  if( new_angle > 1.0e-6 ) {
    double othis_dir[3] = { dir_[3 * dir_num] , dir_[3 * dir_num + 1] ,
                           dir_[3 * dir_num + 2] };
    DACLIB::rotate_about_axis( twiddle_axis_ , angle , this_dir );
    double rnew_angle = DACLIB::angle_about_axis( twiddle_axis_ , in_dir , othis_dir );
    if( rnew_angle < new_angle ) {
      angle_to_use = angle;
    }
  }

  for( int i = 0 ; i < num_dirs_ ; ++i ) {
    DACLIB::rotate_about_axis( twiddle_axis_ , angle_to_use , dir_ + 3 * i );
  }

}

// ****************************************************************************
// this one aligns all directions of given type as best as possible, and
// rearranges the directions to match those in site.
void BasePPhoreSite::twiddle( BasePPhoreSite &site ,
                              GtplDefs::DIRS_TYPE dirs ) {

  if( dir_moves_ != GtplDefs::FREE || !num_dirs_ || !site.num_dirs_ ) {
    return;
  }

  // needs to be done with num_dirs_ >= site.num_dirs_
  if( site.num_dirs_ > num_dirs_ ) {
    site.twiddle( *this , dirs );
    return;
  }

  vector<int> best_comb;
  best_dir_alignments( site , dirs , true , best_comb );

  // make best_comb the directions
  site.twiddle( dir_ , best_comb[0] );

  site.reorder_dirs( best_comb );

}

// ****************************************************************************
// flip directions by 180
void BasePPhoreSite::flip() {

  if( dir_moves_ != GtplDefs::FLIP )
    return;

  for( int i = 0 ; i < num_dirs_ ; ++i )
    DACLIB::rotate_about_axis( twiddle_axis_ , M_PI , dir_ + 3 * i );

}

// ****************************************************************************
// flip directions by 180 if the alignments between dirs_ of type dirs
// gives a better score that way.
void BasePPhoreSite::flip( BasePPhoreSite &site , GtplDefs::DIRS_TYPE dirs ) {

  if( dir_moves_ != GtplDefs::FLIP )
    return;

  vector<int> score1_comb;
  float score1 = best_dir_alignments( site , dirs , false , score1_comb );

  flip();
  vector<int> score2_comb;
  float score2 = best_dir_alignments( site , dirs , false , score2_comb );

  if( score1 > score2 ) {
    // put it back again
    site.reorder_dirs( score1_comb );
    flip();
  } else {
    site.reorder_dirs( score2_comb );
  }

}

// ****************************************************************************
void BasePPhoreSite::set_virt_sites( const double *vs, int num_vs ) {

  num_virt_sites_ = num_vs;
  copy( vs , vs + 3 * num_vs , virt_sites_ );

}

// ****************************************************************************
void BasePPhoreSite::write_to_stream( ostream &os ) const {

  os << label_ << endl << type_code_ << endl << type_string_ << endl
     << cds_[0] << " " << cds_[1] << " " << cds_[2] << endl
     << num_dirs_ << endl;
  for( int i = 0 ; i < num_dirs_ ; ++i ) {
    os << dir_[3*i] << " " << dir_[3*i+1] << " " << dir_[3*i+2]
                    << " " << dir_types_[i] << endl;
  }

  os << dir_moves_ << endl;
  if( dir_moves_ != GtplDefs::NONE ) {
    os << twiddle_axis_[0] << " " << twiddle_axis_[1] << " "
                           << twiddle_axis_[2] << endl;
  }

}

// ****************************************************************************
void BasePPhoreSite::read_from_stream( istream &is ) {

  is >> label_ >> type_code_ >> type_string_
      >> cds_[0] >> cds_[1] >> cds_[2] >> num_dirs_;
  for( int i = 0 ; i < num_dirs_ ; ++i ) {
    is >> dir_[3*i] >> dir_[3*i+1] >> dir_[3*i+2];
    int j;
    is >> j;
    dir_types_[i] = (GtplDefs::DIRS_TYPE) j;
  }

  int i;
  is >> i;
  dir_moves_ = (GtplDefs::DIR_MOVES) i;
  if( dir_moves_ != GtplDefs::NONE ) {
    is >> twiddle_axis_[0] >> twiddle_axis_[1] >> twiddle_axis_[2];
  }

}

// ****************************************************************************
void BasePPhoreSite::brief_report( ostream &os ) const {

  os << "Site " << get_full_name() << " :: ";
  os << "(" << cds_[0] << " , " << cds_[1] << " , " << cds_[2] << ")";
  for( int j = 0 , js = num_dirs_ ; j < js ; ++j ) {
    os << " (" << dir_[3 * j + 0] << " , " << dir_[3 * j + 1]
       << " , " << dir_[3 * j + 2] << ")";
  }
  os << endl;

}

// ****************************************************************************
void brief_report_sites( ostream &os , const vector<BasePPhoreSite *> &sites ) {

  for( unsigned int i = 0 ; i < sites.size() ; i++ ) {
    sites[i]->brief_report( os );
  }

}

// ****************************************************************************
// calculate the overlay transformation to move the given pairs of the second
// vector of PPhoreSites onto the first, returning the RMS of the overlay, but
// not moving the sites.
float calc_overlay_trans( const vector<BasePPhoreSite *> &sites1 , 
                          const vector<BasePPhoreSite *> &sites2 ,
                          const vector<int> &pairs ,
                          OverlayTrans &overlay_trans ,
                          bool use_ring_norm_dirs ,
                          bool use_h_vec_dirs , bool use_lp_dirs ) {

  vector<float> weights( pairs.size() / 2 , 1.0F );
  return calc_overlay_trans( sites1 , sites2 , pairs , weights , overlay_trans ,
                             use_ring_norm_dirs , use_h_vec_dirs , use_lp_dirs );

}

// ****************************************************************************
// calculate the overlay transformation to move the given pairs of the second
// vector of PPhoreSites onto the first, returning the RMS of the overlay, but
// not moving the sites.
float calc_overlay_trans( const vector<BasePPhoreSite *> &sites1 , 
                          const vector<BasePPhoreSite *> &sites2 ,
                          const vector<int> &pairs ,
                          vector<float> &weights ,
                          OverlayTrans &overlay_trans , bool use_ring_norm_dirs ,
                          bool use_h_vec_dirs , bool use_lp_dirs ) {

  vector<float> cds1 , cds2;

  for( unsigned int i = 0 , is = pairs.size() ; i < is ; i += 2 ) {
    const double *sds1 = sites1[pairs[i]]->coords();
    const double *sds2 = sites2[pairs[i+1]]->coords();
    cds1.insert( cds1.end() , sds1 , sds1 + 3 );
    cds2.insert( cds2.end() , sds2 , sds2 + 3 );
  }

  if( use_ring_norm_dirs ) {
    add_ring_norm_dir_sites( sites1 , sites2 , pairs , cds1 , cds2 );
  }
  if( use_h_vec_dirs ) {
    add_dir_sites( GtplDefs::H_VECTOR , sites1 , sites2 , pairs ,
                   cds1 , cds2 );
  }
  if( use_lp_dirs ) {
    add_dir_sites( GtplDefs::LP_VECTOR , sites1 , sites2 , pairs ,
                   cds1 , cds2 );
  }
  vector<float> cds1_cp( cds1 ) , cds2_cp( cds2 );

  // OverlayTrans c'tor moves first coords onto 2nd. It centres the coords, so
  // changes the originals, hence the need for cds1_cp and cds2_cp.
  int num_cds = cds2.size() / 3;
  // might need some more weights
  for( int i = weights.size() , is = num_cds ; i < is ; ++i ) {
    weights.push_back( 1.0F );
  }

  overlay_trans = OverlayTrans( &cds2[0] , &cds1[0] , &weights[0] , num_cds );
  // do the overlay so we can calculate the RMS.
  overlay_trans.overlay( num_cds , &cds2_cp[0] );

  // don't want to calculate RMS on the extension points as well
  int num_sites = pairs.size() / 2;
  float rms = inner_product( cds1_cp.begin() , cds1_cp.begin() + 3 * num_sites , cds2_cp.begin() ,
                             0.0F , plus<float>() ,
                             boost::bind( multiplies<float>() ,
                                          boost::bind( minus<float>() , _1 , _2 ) ,
                                          boost::bind( minus<float>() , _1 , _2 ) ) );

  rms = sqrt( rms / float( num_cds ) );

  return rms;

}

// **************************************************************************
void add_ring_norm_dir_sites( const vector<BasePPhoreSite *> &sites1 , 
                              const vector<BasePPhoreSite *> &sites2 ,
                              const vector<int> &pairs , vector<float> &cds1 ,
                              vector<float> &cds2 ) {

  float pcds1[3] , pcds2[3];

  for( int i = 0 , is = pairs.size() ; i < is ; i += 2 ) {
    BasePPhoreSite *site1 = sites1[pairs[i]];
    BasePPhoreSite *site2 = sites2[pairs[i+1]];
    if( site1->get_num_dirs() && site2->get_num_dirs() ) {
      // find the end of the normal for each site
      const double *dir1 = 0 , *dir2 = 0;
      int site1_dir_num = -1;
      for( int j = 0 ; j < 3 ; ++j ) {
        if( GtplDefs::RING_NORMAL == site1->direction_type( j ) ) {
          dir1 = site1->direction( j );
          site1_dir_num = j;
          break;
        }
      }
      if( !dir1 ) {
        continue;
      }
      int site2_dir_num = -1;
      for( int j = 0 ; j < 3 ; ++j ) {
        if( GtplDefs::RING_NORMAL == site2->direction_type( j ) ) {
          dir2 = site2->direction( j );
          site2_dir_num = j;
          break;
        }
      }
      if( !dir2 ) {
        continue;
      }
      double cdir2[3] = { dir2[0] , dir2[1] , dir2[2] };
      const double *sds1 = site1->coords();
      const double *sds2 = site2->coords();
      pcds1[0] = sds1[0] + dir1[0];
      pcds1[1] = sds1[1] + dir1[1];
      pcds1[2] = sds1[2] + dir1[2];
      // depending on arbitrary factors, ring normals might be pointing in opposite
      // directions. This is why this only works on overlaid sites (not that
      // we're checking this anywhere).
      if( DACLIB::dot_product( site1->direction( site1_dir_num ) ,
                               site2->direction( site2_dir_num ) ) <= 0.0 ) {
        cdir2[0] *= -1.0;
        cdir2[1] *= -1.0;
        cdir2[2] *= -1.0;
        site2->set_direction( cdir2 , site1->direction_type( site1_dir_num ) ,
                              site2_dir_num );
      }
      pcds2[0] = sds2[0] + cdir2[0];
      pcds2[1] = sds2[1] + cdir2[1];
      pcds2[2] = sds2[2] + cdir2[2];
      cds1.insert( cds1.end() , pcds1 , pcds1 + 3 );
      cds2.insert( cds2.end() , pcds2 , pcds2 + 3 );
    }
  }

}

// **************************************************************************
// return the directions of the given type
void get_site_dirs( GtplDefs::DIRS_TYPE dt , BasePPhoreSite &site ,
                    double dirs[9] , int &num_dirs ) {

  num_dirs = 0;
  for( int i = 0 , is = site.get_num_dirs() ; i < is ; ++i ) {
    if( dt == site.direction_type( i ) ) {
      dirs[num_dirs++] = site.direction( i )[0];
      dirs[num_dirs++] = site.direction( i )[1];
      dirs[num_dirs++] = site.direction( i )[2];
    }
  }

  num_dirs /= 3;

}

// **************************************************************************
void add_site_coords( const double *cds , const double *dir ,
                      vector<float> &site_cds ) {

  site_cds.push_back( cds[0] + dir[0] );
  site_cds.push_back( cds[1] + dir[1] );
  site_cds.push_back( cds[2] + dir[2] );

}

// **************************************************************************
void add_dir_sites( GtplDefs::DIRS_TYPE dirs_type ,
                    const vector<BasePPhoreSite *> &sites1 ,
                    const vector<BasePPhoreSite *> &sites2 ,
                    const vector<int> &pairs , vector<float> &cds1 ,
                    vector<float> &cds2 ) {

  // this is all a bit complicated, as there may be more than one vector per site
  // and the site itself may be twiddlable or flippable

  double site1_dirs[9] , site2_dirs[9];
  int nsd1 , nsd2;
  for( int i = 0 , is = pairs.size() ; i < is ; i += 2 ) {
    BasePPhoreSite *site1 = sites1[pairs[i]];
    BasePPhoreSite *site2 = sites2[pairs[i+1]];

    get_site_dirs( dirs_type , *site1 , site1_dirs , nsd1 );
    if( !nsd1 ) continue;
    get_site_dirs( dirs_type , *site2 , site2_dirs , nsd2 );
    if( !nsd2 ) continue;

    if( site1->get_twiddlable() ) {
      site1->twiddle( *site2 , dirs_type );
    }
    if( site1->get_flippable() ) {
      site1->flip( *site2 , dirs_type );
    }

    // having twiddled and flipped site1, find the best combination of
    // site1_dirs and site2_dirs of the same type and add them to site
    // coords. This is because sometimes there can be more than one
    // dir of the same type that is neither flippable nor twiddlable
    // the first example seen being aniline NH2 for DONORS
    float best_dot = -1.0;
    int best1 = -1 , best2 = -1;
    for( int i = 0 ; i < nsd1 ; ++i ) {
      if( dirs_type != site1->direction_type( i ) ) {
        continue;
      }
      for( int j = 0 ; j < nsd2; ++j ) {
        if( dirs_type != site2->direction_type( j ) ) {
          continue;
        }
        float this_dot = DACLIB::dot_product( site1_dirs + 3 * i , site2_dirs + 3 * j );
        if( this_dot > best_dot ) {
          best1 = i;
          best2 = j;
          best_dot = this_dot;
        }
      }
    }
    if( -1 != best1 && -1 != best2 ) {
      add_site_coords( site1->coords() , site1->direction( best1 ) , cds1 );
      add_site_coords( site2->coords() , site2->direction( best2 ) , cds2 );
    }
  }

}

