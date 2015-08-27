//
// file SinglePPhoreSite.cc
// David Cosgrove
// AstraZeneca
// 10th October 2006.
//
// Implementation of SinglePPhoreSite, a concrete class from BasePPhoreSite.

#include <fstream>

#include "gtpl_defs.H"
#include "stddefs.H"
#include "OverlayTrans.H"
#include "SinglePPhoreSite.H"

using namespace boost;
using namespace std;

// *******************************************************************************
SinglePPhoreSite::SinglePPhoreSite() : BasePPhoreSite() {

  reset_data();

}

// *******************************************************************************
// tc and ts are the type code and type string respectively. Throws an int
// exception if num_atoms is too large, giving the maximum number
SinglePPhoreSite::SinglePPhoreSite( const double cds[3] , const double dir[3] , int tc ,
				    const string &ts , const string &lbl ,
				    bool use_dir , int num_atoms ,
				    const int *atoms ,
				    const string &parent_mol_name ) :
  BasePPhoreSite( cds , dir , tc , ts , lbl , use_dir ) {

  if( num_atoms > MAX_SITE_ATOMS )
    throw( MAX_SITE_ATOMS );

  parent_mol_name_ = parent_mol_name;
  
  num_site_atoms_ = num_atoms;
  copy( atoms , atoms + num_atoms , site_atoms_ );
  
}

// *******************************************************************************
// this version builds the SinglePPhoreSite without the atoms and therefore won't
// throw an exception.
SinglePPhoreSite::SinglePPhoreSite( const double cds[3] , const double dir[3] , int tc ,
				    const string &ts , const string &lbl ,
				    bool use_dir ,
				    const string &parent_mol_name ) :
  BasePPhoreSite( cds , dir , tc , ts , lbl , use_dir ) {

  parent_mol_name_ = parent_mol_name;  
  num_site_atoms_ = 0;
  
}

// ****************************************************************************
void SinglePPhoreSite::reset_data() {

  BasePPhoreSite::reset_data();

  parent_mol_name_ = "";
  num_site_atoms_ = 0;

}

// ****************************************************************************
void SinglePPhoreSite::copy_data( const SinglePPhoreSite &c ) {

  BasePPhoreSite::copy_data( c );

  parent_mol_name_ = c.parent_mol_name_;
  num_site_atoms_ = c.num_site_atoms_;
  copy( c.site_atoms_ , c.site_atoms_ + MAX_SITE_ATOMS , site_atoms_ );

}

// ****************************************************************************
string SinglePPhoreSite::get_full_name() const {

  return parent_mol_name_ + ":" + label_;

}

// ****************************************************************************
void SinglePPhoreSite::write_to_stream( ostream &os ) const {

  os << "SinglePPhoreSite" << endl;

  BasePPhoreSite::write_to_stream( os );

  os << parent_mol_name_ << endl;
  os << num_site_atoms_ << " ";
  copy( site_atoms_ , site_atoms_ + num_site_atoms_ ,
	ostream_iterator<int>( os , " " ) );
  os << endl;
  
}

// ****************************************************************************
void SinglePPhoreSite::read_from_stream( istream &is ) {

  BasePPhoreSite::read_from_stream( is );

  is >> parent_mol_name_ >> num_site_atoms_;
  for( int i = 0 ; i < num_site_atoms_ ; ++i )
    is >> site_atoms_[i];

}

// ****************************************************************************
void SinglePPhoreSite::brief_report( ostream &os ) const {

  os << "Site " << get_full_name() << " ::";
  for( int l = 0 ; l < num_site_atoms() ; ++l ) {
    os << " " << site_atoms()[l] + 1;
  }
  os << " :: (" << cds_[0] << " , " << cds_[1] << " , " << cds_[2] << ")";
  for( int j = 0 , js = num_dirs_ ; j < js ; ++j ) {
    os << " (" << dir_[3 * j + 0] << " , " << dir_[3 * j + 1] << " , "
       << dir_[3 * j + 2] << ")";
  }
  os << endl;

}

// ****************************************************************************
ostream &operator<<( ostream &os , const SinglePPhoreSite &site ) {

  site.write_to_stream( os );
  return os;

}

// ****************************************************************************
istream &operator>>( istream &is , SinglePPhoreSite &site ) {

  site.read_from_stream( is );
  return is;

}

// ****************************************************************************
// calculate the overlay transformation to move the given pairs of the second
// vector of PPhoreSites onto the first, returning the RMS of the overlay, but
// not moving the sites.
float calc_overlay_trans( const vector<BasePPhoreSite *> &sites1 , 
                          const vector<SinglePPhoreSite *> &sites2 ,
                          const vector<int> &pairs ,
                          OverlayTrans &overlay_trans ,
                          bool use_ring_norm_dirs ,
                          bool use_h_vec_dirs , bool use_lp_dirs ) {

  vector<BasePPhoreSite *> base_sites2( sites2.begin() , sites2.end() );
  return calc_overlay_trans( sites1 , base_sites2 , pairs ,
                             overlay_trans , use_ring_norm_dirs ,
                             use_h_vec_dirs , use_lp_dirs );

}

// ****************************************************************************
float calc_overlay_trans( const vector<SinglePPhoreSite *> &sites1 , 
                          const vector<SinglePPhoreSite *> &sites2 ,
                          const vector<int> &pairs ,
                          OverlayTrans &overlay_trans ,
                          bool use_ring_norm_dirs ,
                          bool use_h_vec_dirs , bool use_lp_dirs ) {

  vector<BasePPhoreSite *> base_sites1( sites1.begin() , sites1.end() );
  return calc_overlay_trans( base_sites1 , sites2 , pairs , overlay_trans ,
                             use_ring_norm_dirs , use_h_vec_dirs , use_lp_dirs );

}

// **************************************************************************
float calculate_site_rms( const vector<BasePPhoreSite *> &sites1 ,
                          const vector<SinglePPhoreSite *> &sites2 ,
                          const vector<int> &pairs ) {

  float rms = 0.0F;
  for( int i = 0 , is = pairs.size() ; i < is ; i += 2 ) {
    BasePPhoreSite *site1 = sites1[pairs[i]];
    BasePPhoreSite *site2 = sites2[pairs[i+1]];
    rms += DACLIB::sq_distance( site1->coords() , site2->coords() );
  }

  return sqrt( rms / float( pairs.size() / 2 ) );

}

// **************************************************************************
float calculate_site_rms( const vector<BasePPhoreSite *> &sites1 ,
                          const vector<shared_ptr<SinglePPhoreSite> > &sites2 ,
                          const vector<int> &pairs ) {

  float rms = 0.0F;
  for( int i = 0 , is = pairs.size() ; i < is ; i += 2 ) {
    rms += DACLIB::sq_distance( sites1[pairs[i]]->coords() ,
                                sites2[pairs[i+1]]->coords() );
  }

  return sqrt( rms / float( pairs.size() / 2 ) );

}

// **************************************************************************
void overlay( SinglePPhoreSite &site , const OverlayTrans &ot ) {

  float trans1[3] , trans2[3] , rot[3][3];
  ot.get_trans1( trans1 );
  ot.get_trans2( trans2 );
  ot.get_rot( rot );

  site.translate( -trans1[0] , -trans1[1] , -trans1[2] );
  site.rotate( rot );
  site.translate( trans2[0] , trans2[1] , trans2[2] );

}

// **************************************************************************
void overlay( vector<SinglePPhoreSite *> &sites , const OverlayTrans &ot ) {

  float trans1[3] , trans2[3] , rot[3][3];
  ot.get_trans1( trans1 );
  ot.get_trans2( trans2 );
  ot.get_rot( rot );

  vector<SinglePPhoreSite *>::iterator p;
  for( p = sites.begin() ; p != sites.end() ; ++p ) {
    overlay( *(*p) , ot );
  }

}

// **************************************************************************
void overlay( vector<boost::shared_ptr<SinglePPhoreSite> > &sites , const OverlayTrans &ot ) {

  float trans1[3] , trans2[3] , rot[3][3];
  ot.get_trans1( trans1 );
  ot.get_trans2( trans2 );
  ot.get_rot( rot );

  vector<boost::shared_ptr<SinglePPhoreSite> >::iterator p;
  for( p = sites.begin() ; p != sites.end() ; ++p ) {
    overlay( *(*p) , ot );
  }

}

// **************************************************************************
// confirm that the sites named in the given clique have directions within the
// tolorance. The sites will already have been overlaid.
bool confirm_site_vectors( const vector<int> &clique ,
                           const vector<BasePPhoreSite *> &query_sites ,
                           vector<shared_ptr<SinglePPhoreSite> > &target_sites ,
                           bool check_ring_norms , float ring_norm_tol ,
                           bool check_h_vectors , float h_vector_tol ,
                           bool check_lps , float lp_tol ) {

  if( !check_ring_norms && !check_h_vectors && !check_lps ) {
    return true; // we don't care what they're doing, so the clique must be ok
  }

  for( int i = 0 , is = clique.size() ; i < is ; i += 2 ) {

    const BasePPhoreSite &qs = *query_sites[clique[i]];
    SinglePPhoreSite &ts = *target_sites[clique[i+1]];

    if( check_ring_norms && !confirm_ring_norms( qs , ts , ring_norm_tol ) ) {
      return false;
    }

    if( check_h_vectors && !confirm_dir_vectors( qs , ts , h_vector_tol ,
                                                 GtplDefs::H_VECTOR ) ) {
      return false;
    }

    if( check_lps && !confirm_dir_vectors( qs , ts , lp_tol ,
                                           GtplDefs::LP_VECTOR ) ) {
      return false;
    }

  }
  return true;

}

// **************************************************************************
bool confirm_ring_norms( const BasePPhoreSite &site1 ,
                         SinglePPhoreSite &site2 , float tol ) {

  for( int i = 0 , is = site1.get_num_dirs() ; i < is ; ++i ) {
    if( GtplDefs::RING_NORMAL != site1.direction_type( i ) ) {
      continue;
    }
    for( int j = 0 , js = site2.get_num_dirs() ; j < js ; ++j ) {
      if( GtplDefs::RING_NORMAL != site2.direction_type( j ) ) {
        continue;
      }
      float dotp = DACLIB::dot_product( site1.direction( i ) ,
                                        site2.direction( j ) );
      if( fabs( dotp ) <= tol ) {
        return false;
      }
    }
  }

  return true;

}

// **************************************************************************
bool confirm_dir_vectors( const BasePPhoreSite &site1 ,
                          SinglePPhoreSite &site2 , float tol ,
                          GtplDefs::DIRS_TYPE dirs_type ) {

  bool tried_a_dir = false;
  for( int i = 0 , is = site1.get_num_dirs() ; i < is ; ++i ) {
    if( dirs_type != site1.direction_type( i ) ) {
      continue;
    }

    for( int j = 0 , js = site2.get_num_dirs() ; j < js ; ++j ) {
      if( dirs_type != site2.direction_type( j ) ) {
        continue;
      }

      tried_a_dir = true;

      if( site2.get_twiddlable() ) {
        site2.twiddle( site1.direction( i ) , j );
      }
      float dotp = DACLIB::dot_product( site1.direction( i ) ,
                                        site2.direction( j ) );
      if( dotp > tol ) {
        return true; // got a match, can return straight away
      }
      // might be flippable
      if( site2.get_flippable() ) {
        site2.flip();
        float dotp = DACLIB::dot_product( site1.direction( i ) ,
                                          site2.direction( j ) );
        if( dotp > tol ) {
          return true; // got a match, can return straight away
        }
      }
    }
  }

  // if we tested at least on 1 vector, and we got here, then it failed. If we
  // didn't and we got here, it passed.
  return !tried_a_dir;

}


