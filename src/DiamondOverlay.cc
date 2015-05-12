//
// file DiamondOverlay.cc
// David Cosgrove
// AstraZeneca
// 23rd January 2007
//

#include <iterator>
#include <iostream>

#include "stddefs.H"
#include "DiamondOverlay.H"

#include <boost/lexical_cast.hpp>

namespace DACLIB {
  void eigen_solve_4b4f( float **a , float &value , float *vector );
}

using namespace std;

namespace DACLIB {

  // ****************************************************************************
  DiamondOverlay::DiamondOverlay( const vector<float> &fixed ,
                                  const vector<float> &moving ) :
    rot_matrix_( 0 ) , p_( 0 ) , e_zero_( -1.0 ) , e_( -1.0 ){

    if( fixed.size() != moving.size() ) {
      throw DiamondOverlayError( "fixed and moving vectors must be same size." );
    }

    fixed_ = fixed;
    moving_ = moving;

  }

  // ****************************************************************************
  DiamondOverlay::DiamondOverlay( const vector<float> &fixed ,
                                  const vector<float> &moving ,
                                  const vector<float> &weights ) :
    rot_matrix_( 0 ) , p_( 0 ) , e_zero_( -1.0 ) , e_( -1.0 ) {

    if( fixed.size() != moving.size() ) {
      throw DiamondOverlayError( "fixed and moving vectors must be same size." );
    }

    fixed_ = fixed;
    moving_ = moving;
    weights_ = weights;

  }

  // ****************************************************************************
  DiamondOverlay::DiamondOverlay( const float *fixed , const float *moving ,
                                  int num_points ) :
    rot_matrix_( 0 ) , p_( 0 ) , e_zero_( -1.0 ) , e_( -1.0 ) {

    fixed_ = vector<float>( fixed , fixed + 3 * num_points );
    moving_ = vector<float>( moving , moving + 3 * num_points );

  }

  // ****************************************************************************
  DiamondOverlay::DiamondOverlay( const float *fixed , const float *moving ,
                                  const float *weights , int num_points ) :
    rot_matrix_( 0 ) , p_( 0 ) , e_zero_( -1.0 ) , e_( -1.0 ) {

    fixed_ = vector<float>( fixed , fixed + 3 * num_points );
    moving_ = vector<float>( moving , moving + 3 * num_points );
    weights_ = vector<float>( weights , weights + num_points );

  }

  // ****************************************************************************
  // individual sets of coords.
  DiamondOverlay::DiamondOverlay( const vector<vector<float> > &fixed ,
                                  const vector<vector<float> > &moving ) :
    rot_matrix_( 0 ) , p_( 0 ) , e_zero_( -1.0 ) , e_( -1.0 ) {

    fixed_.reserve( 3 * fixed.size() );
    moving_.reserve( 3 * moving.size() );
    for( int i = 0 , is = fixed.size() ; i < is ; ++i ) {
      if( fixed[i].size() != 3 ) {
        throw DiamondOverlayError( "Bad coords in fixed coords." );
      }
      if( moving[i].size() != 3 ) {
        throw DiamondOverlayError( "Bad coords in moving coords." );
      }
      
      fixed_.insert( fixed_.end() , fixed[i].begin() , fixed[i].end() );
      moving_.insert( moving_.end() , moving[i].begin() , moving[i].end() );
    }

  }

  // ****************************************************************************
  // individual sets of coords.
  DiamondOverlay::DiamondOverlay( const vector<vector<float> > &fixed ,
                                  const vector<vector<float> > &moving ,
                                  const vector<float> &weights ) :
    rot_matrix_( 0 ) , p_( 0 ) , e_zero_( -1.0 ) , e_( -1.0 ) {

    fixed_.reserve( 3 * fixed.size() );
    moving_.reserve( 3 * moving.size() );
    for( int i = 0 , is = fixed.size() ; i < is ; ++i ) {
      if( fixed[i].size() != 3 ) {
        throw DiamondOverlayError( "Bad coords in fixed coords." );
      }
      if( moving[i].size() != 3 ) {
        throw DiamondOverlayError( "Bad coords in moving coords." );
      }
      
      fixed_.insert( fixed_.end() , fixed[i].begin() , fixed[i].end() );
      moving_.insert( moving_.end() , moving[i].begin() , moving[i].end() );
    }

    weights_ = weights;

  }

  // ****************************************************************************
  DiamondOverlay::~DiamondOverlay() {

    DACLIB::destroy_square_matrix( rot_matrix_ );
    DACLIB::destroy_square_matrix( p_ );

  }

  // ****************************************************************************
  // calculate the initial residual error
  void DiamondOverlay::calc_e_zero() {

    e_zero_ = 0.0F;
    for( int i = 0 , is = fixed_.size() ; i < is ; i += 3 ) {
      float sq_dist = sq_distance( &fixed_[i] , &moving_[i] );
      if( !weights_.empty() ) {
        sq_dist *= weights_[i/3];
      }
      e_zero_ += sq_dist;
    }

  }

  // ****************************************************************************
  // assemble the M matrix (eq. 16, D88). It is not a symmetric matrix.
  void DiamondOverlay::calc_m( float m[3][3] ) {

    for( int i = 0 ; i < 3 ; ++i ) {
      for( int j = 0 ; j < 3 ; ++j ) {
        m[i][j] = 0.0;
        for( int k = 0 , ks = fixed_.size() ; k < ks ; k += 3 ) {
          if( weights_.empty() ) {
            m[i][j] += moving_[k+i] * fixed_[k+j];
          } else {
            m[i][j] += moving_[k+i] * fixed_[k+j] * weights_[k/3];
	  }
        }
      }
    }

  }

  // ****************************************************************************
  // calculate matrix Q, equation 17, D88.
  void DiamondOverlay::calc_q( const float m[3][3] , float q[3][3] ) {

    float alpha = m[0][0] + m[1][1] + m[2][2];

    for( int i = 0 ; i < 3 ; ++i ) {
      for( int j = 0 ; j < 3 ; ++j ) {
        q[i][j] = m[i][j] + m[j][i];
      }
      q[i][i] -= 2.0 * alpha;
    }

  }

  // ****************************************************************************
  // calculate matrix Q, equation 17, D88.
  void DiamondOverlay::calc_e() {

    if( !rot_matrix_ ) {
      calc_rot_matrix();
    }

    float res[4];
    float e;

    res[0] = p_[0][0] * vector_[0] + p_[0][1] * vector_[1] +
        p_[0][2] * vector_[2] + p_[0][3] * vector_[3];

    res[1] = p_[1][0] * vector_[0] + p_[1][1] * vector_[1] +
        p_[1][2] * vector_[2] + p_[1][3] * vector_[3];

    res[2] = p_[2][0] * vector_[0] + p_[2][1] * vector_[1] +
        p_[2][2] * vector_[2] + p_[2][3] * vector_[3];

    res[3] = p_[3][0] * vector_[0] + p_[3][1] * vector_[1] +
        p_[3][2] * vector_[2] + p_[3][3] * vector_[3];

    e = vector_[0] * res[0] + vector_[1] * res[1] + vector_[2] * res[2] +
        vector_[3] * res[3];

    e_ = e_zero_ - 2.0F * e;

  }

  // ****************************************************************************
  // make sure that there are at least 3 distinct sets of coordinates
  // for at least 1 of the input coord sets.  This includes cases
  // where 2 coords are on top of each other. Throws exception if not
  void DiamondOverlay::check_input_coords() const {

    if( !(at_least_3_coords( fixed_ ) && at_least_3_coords( moving_ )) ) {
      throw DiamondOverlayError( "Not enough distinct coords for overlay" );
    }

  }

  // ****************************************************************************
  // checks given coords for problem above.
  bool DiamondOverlay::at_least_3_coords( const std::vector<float> &cds ) const {

    vector<float> eff_cds = cds;
    bool it_changed( false );

    do {
      it_changed = false;
      for( int i = 0 , is = eff_cds.size() - 3 ; i < is ; i +=3 ) {
        for( int j = i + 3 , js = is + 3 ; j < js ; j += 3 ) {
          if( sq_distance( &cds[i] , &cds[j] ) < 1.0e-3 ) {
            eff_cds.erase( eff_cds.begin() + i , eff_cds.begin() + i + 3 );
            it_changed = true;
            break;
          }
        }
        if( it_changed ) {
          break;
        }
      }

    } while( it_changed && eff_cds.size() > 6 );

    return( eff_cds.size() > 6 );

  }

  // ****************************************************************************
  void DiamondOverlay::calc_p() {

    float m[3][3] , q[3][3];

    calc_e_zero();
    calc_m( m );
    calc_q( m , q );

    // compute v, eq. 18, D88.
    float v[3];
    v[0] = m[1][2] - m[2][1];
    v[1] = m[2][0] - m[0][2];
    v[2] = m[0][1] - m[1][0];

    // now do p
    DACLIB::make_square_matrix( p_ , 4 );

    for( int i = 0 ; i < 3 ; ++i ) {
      for( int j = 0 ; j < 3 ; ++j ) {
        p_[i][j] = q[i][j];
      }
    }

    p_[0][3] = p_[3][0] = v[0];
    p_[1][3] = p_[3][1] = v[1];
    p_[2][3] = p_[3][2] = v[2];
    p_[3][3] = 0.0;

  }

  // ****************************************************************************
  // build the rotation matrix, using eq. 7, D88.
  void DiamondOverlay::calc_rot_matrix() {

    check_input_coords();
    if( !p_ ) {
      calc_p();
    }

    if( !rot_matrix_ ) {
      DACLIB::make_square_matrix( rot_matrix_ , 3 );
    }

    float value;

    // eigen_solve destroys the matrix it works on, so do it on a copy
    static float **p_tmp = 0;
    if( !p_tmp ) {
      DACLIB::make_square_matrix( p_tmp , 4 );
    }

    copy( p_[0] , p_[0] + 16 , p_tmp[0] );
    DACLIB::eigen_solve_4b4f( p_tmp , value , vector_ );

    for( int i = 0 ; i < 4 ; ++i ) {
      try{
        boost::lexical_cast<float>( vector_[i] );
      } catch( boost::bad_lexical_cast &e ) {
        throw DiamondOverlayError( "Failed to find to valid overlay." );
      }
    }

    // use the notation of eq. 1, D88.
    const float lambda = vector_[0];
    const float mu = vector_[1];
    const float nu = vector_[2];
    const float sigma = vector_[3];
    const float mu_sq = mu * mu;
    const float lambda_sq = lambda * lambda;
    const float ups_sq = nu * nu;
    const float sigma_sq = sigma * sigma;

    float *rr = rot_matrix_[0];
    rr[0] = lambda_sq - mu_sq - ups_sq + sigma_sq;
    rr[4] = -lambda_sq + mu_sq - ups_sq + sigma_sq;
    rr[8] = -lambda_sq - mu_sq + ups_sq + sigma_sq;
    rr[3] = 2.0 * ( lambda * mu + nu * sigma );
    rr[6] = 2.0 * ( lambda * nu - mu * sigma );
    rr[1] = 2.0 * ( lambda * mu - nu * sigma );
    rr[7] = 2.0 * ( mu * nu + lambda * sigma );
    rr[2] = 2.0 * ( lambda * nu + mu * sigma );
    rr[5] = 2.0 * ( mu * nu - lambda * sigma );

  }

  // ****************************************************************************
  void DiamondOverlay::get_rot_matrix( float rot[3][3] ) {

    if( !p_ || !rot_matrix_ ) {
      calc_rot_matrix();
    }

    // don't assume that rot is created with DACLIB::make_square_matrix, though
    // we know that rot_matrix__ was.
    float *rr = rot_matrix_[0];
    rot[0][0] = rr[0];
    rot[0][1] = rr[1];
    rot[0][2] = rr[2];
    rot[1][0] = rr[3];
    rot[1][1] = rr[4];
    rot[1][2] = rr[5];
    rot[2][0] = rr[6];
    rot[2][1] = rr[7];
    rot[2][2] = rr[8];

  }

  // ****************************************************************************
  void DiamondOverlay::get_p( float **p ) {

    if( !p_ ) {
      calc_p();
    }

    // don't assume that p is created with DACLIB::make_square_matrix, though
    // we know that p_ was.
    float *pp = p_[0];
    p[0][0] = pp[0];
    p[0][1] = pp[1];
    p[0][2] = pp[2];
    p[0][3] = pp[3];

    p[1][0] = pp[4];
    p[1][1] = pp[5];
    p[1][2] = pp[6];
    p[1][3] = pp[7];

    p[2][0] = pp[8];
    p[2][1] = pp[9];
    p[2][2] = pp[10];
    p[2][3] = pp[11];

    p[3][0] = pp[12];
    p[3][1] = pp[13];
    p[3][2] = pp[14];
    p[3][3] = pp[15];
    
  }

  // ****************************************************************************
  float DiamondOverlay::get_e_zero() {

    if( e_zero_ < 0.0F ) {
      calc_e_zero();
    }

    return e_zero_;

  }

  // ****************************************************************************
  float DiamondOverlay::get_e() {
    if( e_ < 0.0F ) {
      calc_e() ;
    }

    return e_;

  }

} // end of namespace DACLIB
