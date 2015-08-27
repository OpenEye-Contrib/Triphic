//
// file BasePPhoreSite.H
// David Cosgrove
// AstraZeneca
// 10th October 2006
//
// This is the interface for the abstract base class BasePPhoreSite, from
// which will be derived at least SinglePPhoreSite and MultiPPhoreSite, the latter
// being created from 1 or more of the former after, for example, clustering.

#ifndef DAC_BASE_PPHORE_SITE__
#define DAC_BASE_PPHORE_SITE__

#include <cmath>
#include <string>
#include <vector>

#include "gtpl_defs.H"

// *****************************************************************************

class OverlayTrans;

class BasePPhoreSite {

protected :

  std::string   label_;

  int      type_code_; /* the type of the pharmacophore site, indexing into the
        external list of types */
  std::string   type_string_; // the string equivalent, won't be used for searching

  double   cds_[3];
  mutable double dir_[9]; /* direction for the site if, e.g., it's a ring
           centroid. We're allowing up to 3 now, for
           H-vectors and lone pairs as well. */
  int      num_dirs_; // number of dirs to use
  GtplDefs::DIRS_TYPE dir_types_[3];

  GtplDefs::DIR_MOVES dir_moves_; // whether the direction vector can be rotated
  double    twiddle_axis_[3]; /* axis around which it may be rotated for instance
         if the vector is an OH, the C-O bond */

  // virtual sites, as defined in JCAMD 26, 451-472 (2012) Robin's pphore paper.
  // there can be up to 12 sites.  The positions here are, like the dir_ stuff,
  // relative to cds_ so to get their actual positions they need to be added
  // to cds_.
  double   virt_sites_[36];
  int      num_virt_sites_;
  // as part of the scoring of overlays, 1 of the virtual sites might be selected
  // as the key one, in that it clusters with virtual sites from other sites.
  // If there is such a site, gives its number, otherwise -1.
  int      active_virt_site_;

  bool     selected_;

  virtual void reset_data();
  virtual void copy_data( const BasePPhoreSite &c );

public :

  BasePPhoreSite();
  // tc and ts are the type code and type string respectively.
  BasePPhoreSite( const double cds[3] , const double dir[3] , int tc , const std::string &ts ,
  const std::string &lbl , bool use_dir );
  BasePPhoreSite( const BasePPhoreSite &c ) { copy_data( c ); };

  virtual ~BasePPhoreSite() = 0;

  BasePPhoreSite &operator=( const BasePPhoreSite &c ) {
    if( this == &c )
      return *this;
    copy_data( c );
    return *this;
  }
  
  // access functions
  std::string label() const { return label_; }
  virtual std::string parent_mol_name() const = 0;
  virtual std::string get_full_name() const = 0;
  int get_num_dirs() const { return num_dirs_; }
  int get_num_dirs( GtplDefs::DIRS_TYPE ) const;

  const double *coords() { return cds_; }
  void set_coords( const double *new_cds );
  void set_coords( const float *new_cds );

  const double *direction( int dir_num = 0 ) const;
  GtplDefs::DIRS_TYPE direction_type( int dir_num ) const;

  void set_direction( const double *new_dir , GtplDefs::DIRS_TYPE t ,
                      int dir_num = -1 );
  void set_direction( const float *new_dir , GtplDefs::DIRS_TYPE t ,
                      int dir_num = -1 );

  void set_selected( bool new_val ) { selected_ = new_val; }
  bool get_selected() const { return selected_; }
  int  get_type_code() const { return type_code_; }
  const std::string &get_type_string() const { return type_string_; }
  // A cheap (but not always completely reliable) way of distinguishing a
  // MultiPPhoreSite from others.
  virtual int num_base_sites() const { return 1; }
  
  float square_distance( const double *cds ) const;
  float distance( const double *cds ) const {
    return sqrt( square_distance( cds ) );
  }

  void translate( float x_trans , float y_trans , float z_trans );
  void rotate( const float rot[3][3] );

  bool get_twiddlable() const { return dir_moves_ == GtplDefs::FREE; }
  bool get_flippable() const { return dir_moves_ == GtplDefs::FLIP; }
  GtplDefs::DIR_MOVES get_dir_moves() const { return dir_moves_; }
  void set_twiddlable( bool is_twid ) {
    dir_moves_ = is_twid ? GtplDefs::FREE : GtplDefs::NONE;
  }
  void set_flippable( bool is_flip ) {
    dir_moves_ = is_flip ? GtplDefs::FLIP : GtplDefs::NONE;
  }
  void set_dir_moves( GtplDefs::DIR_MOVES dir_moves ) {
    dir_moves_ = dir_moves;
  }
  const double *twiddle_axis() const { return twiddle_axis_; }
  void set_twiddle_axis( const float axis[3] ) {
    twiddle_axis_[0] = axis[0]; twiddle_axis_[1] = axis[1];
    twiddle_axis_[2] = axis[2];
  }
  void set_twiddle_axis( const double axis[3] ) {
    twiddle_axis_[0] = axis[0]; twiddle_axis_[1] = axis[1];
    twiddle_axis_[2] = axis[2];
  }

  // find the order of directions in site that has the most of the given
  // type lined up with directions in this. Dir_order maps dirs_ vectors in
  // site onto dirs_ vectors in this. Returns the sum of the dot-products for
  // the order found.
  float best_dir_alignments( BasePPhoreSite &site ,
                             GtplDefs::DIRS_TYPE dirs_type ,
                             bool twiddle_if_poss ,
                             std::vector<int> &dir_order );

  // put the direction vectors into the new order.
  void reorder_dirs( const std::vector<int> &new_order );

  // twiddle the given direction to point along the given vector dir
  // if it's twiddlable leaving it unchanged if it isn't
  void twiddle( const double *dir , int dir_num );
  // this one aligns all directions of given type as best as possible, and
  // rearranges the directions to match those in site
  void twiddle( BasePPhoreSite &site , GtplDefs::DIRS_TYPE dirs );
  // flip directions by 180
  void flip();
  // flip directions by 180 if the alignments between dirs_ of type dirs
  // gives a better score that way.
  void flip( BasePPhoreSite &site , GtplDefs::DIRS_TYPE dirs );

  int num_virt_sites() const { return num_virt_sites_; }
  const double *virt_sites() const { return virt_sites_; }
  void set_virt_sites( const double *vs , int num_vs );
  int active_virt_site() const { return active_virt_site_; }
  void set_active_virt_site( int new_site = -1 ) { active_virt_site_ = new_site; }

  virtual void write_to_stream( std::ostream &os ) const;
  virtual void read_from_stream( std::istream &is );
  virtual void brief_report( std::ostream &os ) const;

  virtual std::string type_string_for_stream() = 0;

};

void brief_report_sites( std::ostream &os ,
                         const std::vector<BasePPhoreSite *> &sites );

// calculate the overlay transformation to move the given pairs of the second
// vector of PPhoreSites onto the first, returning the RMS of the overlay, but
// not moving the sites.
float calc_overlay_trans( const std::vector<BasePPhoreSite *> &sites1 , 
                          const std::vector<BasePPhoreSite *> &sites2 ,
                          const std::vector<int> &pairs ,
                          OverlayTrans &overlay_trans , bool use_ring_norm_dirs ,
                          bool use_h_vec_dirs , bool use_lp_dirs );
float calc_overlay_trans( const std::vector<BasePPhoreSite *> &sites1 , 
                          const std::vector<BasePPhoreSite *> &sites2 ,
                          const std::vector<int> &pairs ,
                          std::vector<float> &weights ,
                          OverlayTrans &overlay_trans , bool use_ring_norm_dirs ,
                          bool use_h_vec_dirs , bool use_lp_dirs );

void add_ring_norm_dir_sites( const std::vector<BasePPhoreSite *> &sites1 , 
                              const std::vector<BasePPhoreSite *> &sites2 ,
                              const std::vector<int> &pairs ,
                              std::vector<float> &cds1 ,
                              std::vector<float> &cds2 );
void add_dir_sites( GtplDefs::DIRS_TYPE dirs_type ,
                    const std::vector<BasePPhoreSite *> &sites1 ,
                    const std::vector<BasePPhoreSite *> &sites2 ,
                    const std::vector<int> &pairs ,
                    std::vector<float> &cds1 ,
                    std::vector<float> &cds2 );

#endif
