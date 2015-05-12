//
// file PPhoreQuery.cc
// Dave Cosgrove
// 13th September 2007
//

#include <fstream>
#include <iostream>
#include <sstream>

#include <boost/tuple/tuple_io.hpp>

#include <oechem.h>
#include <pvm3.h>

#include "stddefs.H"

#include "FileExceptions.H"
#include "PPhoreQuery.H"
#include "PluralityHit.H"
#include "SinglePPhoreSite.H"

using namespace std;
using namespace OEChem;

namespace DACLIB {
  double angle( const double *cds1 , const double *cds2 , const double *cds3 ,
		bool degs = true );
  double torsion( const double *cds1 , const double *cds2 , const double *cds3 ,
		  const double *cds4 , bool degs = true );
  void pack_strings_vector( const vector<string> &strs );
  void unpack_strings_vector( vector<string> &strs );
}

//****************************************************************************
PPhoreQuery::PPhoreQuery() : sq_dist_vals_( 0 ) , overlay_cds_( 0 ) {

}

//****************************************************************************
PPhoreQuery::~PPhoreQuery() {

  delete [] overlay_cds_;
  DACLIB::destroy_square_matrix( sq_dist_vals_ );

}

//****************************************************************************
void PPhoreQuery::read_query_file( const string &filename ) {

  filename_ = filename;
  ifstream ifs( filename_.c_str() );
  if( !ifs.good() )
    throw DACLIB::FileReadOpenError( filename_.c_str() );

  string next_line , line_type;
  while( 1 ) {
    next_line.clear();
    getline( ifs , next_line , '\n' );
    // if the last line ended without '\n', eof is returned even though
    // something was read
    
    if( next_line.empty() ) {
      if( ifs.eof() || !ifs.good() )
	break;
      else
	continue; // it's a blank line
    }
    if( '#' == next_line[0] )
      continue; // it's a comment line
    istringstream ss( next_line );
    
    ss >> line_type;
 
    if( string( "Point" ) == line_type ) {
      parse_point_line( ss );
    } else if( string( "Coords" ) == line_type ) {
      parse_coords_line( ss );
    } else if( string( "Distance" ) == line_type ) {
      parse_distance_line( ss );
    } else if( string( "Angle" ) == line_type ) {
      parse_angle_line( ss );
    } else if( string( "Torsion" ) == line_type ) {
      parse_torsion_line( ss );
    } else if( string( "hxvol" ) == line_type ) {
      parse_exc_vol_line( ss , hard_exc_vols_ );
    } else if( string( "sxvol" ) == line_type ) {
      parse_exc_vol_line( ss , soft_exc_vols_ );
    } else {
      throw PPhoreQueryFileReadError( filename_ , next_line );
    }
  }

  verify_query_details();
  assemble_overlay_coords();
  assemble_distance_constraints();

}

//****************************************************************************
void PPhoreQuery::search_pphore_sites( vector<SinglePPhoreSite *> &target_sites ,
                                       OEMol &target_mol , int conf_num ,
                                       boost::shared_ptr<DACLIB::VolumeGrid> &protein_grid ,
                                       boost::shared_ptr<DACLIB::VolumeGrid> &soft_exc_vol_grid ,
                                       const vector<pair<string,DACLIB::VolumeGrid *> > &score_vol_grids ,
                                       bool first_hit_only ,
                                       vector<boost::shared_ptr<PluralityHit> > &hits ) {

#ifdef NOTYET
  cout << "PPhoreQuery::search_pphore_sites" << endl;
#endif

  if( target_sites.size() < point_types_.size() ) {
#ifdef NOTYET
    cout << "Not enough sites : " << target_sites.size() << " vs " << point_types_.size() << endl;
#endif
    return; // can't be a solution, there aren't enough sites
  }

  // make a squared distance matrix
  float **sq_dist_mat = make_squared_distance_matrix( target_sites );
  short **m = form_m0( target_sites );

  if( !solution_possible( m , point_types_.size() , target_sites.size() ) ) {
    DACLIB::destroy_2d_matrix( m );
    DACLIB::destroy_square_matrix( sq_dist_mat );
    return;
  }

  vector<vector<int> > hit_points;
  extend_subgraph( 0 , m , point_types_.size() , target_sites.size() ,
                   const_cast<const float **>( sq_dist_mat ) , hit_points );

#ifdef NOTYET
  if( !hit_points.empty() ) {
    cout << "Number of hits : " << hit_points.size() << endl;
  }
#endif

  vector<float> hit_angles , hit_torsions;
  for( int i = 0 , is = hit_points.size() ; i < is ; ++i ) {
#ifdef NOTYET
    cout << endl << "Doing hit " << i << " : ";
    for( int jj = 0 , jjs = hit_points[i].size() ; jj < jjs ; ++jj )
      cout << hit_points[i][jj] << " ";
    cout << endl;
#endif
    if( !angles_ok( target_sites , hit_points[i] , hit_angles ) ||
        !torsions_ok( target_sites , hit_points[i] , hit_torsions ) ) {
#ifdef NOTYET
      cout << "fails angles and torsions" << endl;
#endif
      continue;
    }

#ifdef NOTYET
    cout << "angles and torsions ok, making hit" << endl;
#endif
    boost::shared_ptr<PluralityHit> hit( new PluralityHit( target_mol ,
                                                           conf_num ,
                                                           target_sites ,
                                                           hit_points[i] ,
                                                           protein_grid ,
                                                           soft_exc_vol_grid ,
                                                           hard_exc_vols_ ,
                                                           score_vol_grids ,
                                                           coords_ ,
                                                           distances_ ,
                                                           hit_angles ,
                                                           hit_torsions ) );
    // no ov_conf means no hit, for some reason such as hard_exc_vols_ were
    // breached
    if( hit->get_ov_conf() ) {
#ifdef NOTYET
      cout << "ov conf ok, so storing hit" << endl;
#endif
      hits.push_back( hit );
    } else {
#ifdef NOTYET
      cout << "ov conf not ok, hit not stored" << endl;
#endif
    }
    if( first_hit_only && !hits.empty() ) {
      break;
    }
  }

  DACLIB::destroy_2d_matrix( m );
  DACLIB::destroy_square_matrix( sq_dist_mat );

}

//****************************************************************************
void PPhoreQuery::clear_data() {

  DACLIB::destroy_square_matrix( sq_dist_vals_ );
  sq_dist_vals_ = 0;
  delete [] overlay_cds_;
  overlay_cds_ = 0;

  point_types_.clear();
  coords_.clear();
  distances_.clear();
  angles_.clear();
  torsions_.clear();
  hard_exc_vols_.clear();
  soft_exc_vols_.clear();

}

//****************************************************************************
void PPhoreQuery::parse_point_line( istringstream &ss ) {

  string point_type;
  ss >> point_type;
  if( ss.fail() )
    throw PPhoreQueryFileReadError( filename_ , ss.str() );

  point_types_.push_back( point_type );

}

//****************************************************************************
void PPhoreQuery::parse_coords_line( istringstream &ss ) {

  int i;
  float x , y , z;
  ss >> i >> x >> y >> z;

  if( ss.fail() )
    throw PPhoreQueryFileReadError( filename_ , ss.str() );

  coords_.push_back( boost::make_tuple( i , x , y , z ) );

}

//****************************************************************************
void PPhoreQuery::parse_distance_line( istringstream &ss ) {

  int i , j;
  float ld , hd;
  ss >> i >> j >> ld >> hd;
  if( ss.fail() )
    throw PPhoreQueryFileReadError( filename_ , ss.str() );

  distances_.push_back( boost::make_tuple( i , j , ld , hd ) );

}

//****************************************************************************
void PPhoreQuery::parse_angle_line( istringstream &ss ) {

  int i , j , k;
  float la , ha;
  ss >> i >> j >> k >> la >> ha;
  if( ss.fail() )
    throw PPhoreQueryFileReadError( filename_ , ss.str() );

  angles_.push_back( boost::make_tuple( i , j , k , la , ha ) );

}

//****************************************************************************
void PPhoreQuery::parse_torsion_line( istringstream &ss ) {

  int i , j , k , l;
  float lt , ht;
  ss >> i >> j >> k >> l >> lt >> ht;
  if( ss.fail() )
    throw PPhoreQueryFileReadError( filename_ , ss.str() );

  torsions_.push_back( boost::make_tuple( i , j , k , l , lt , ht ) );

}

//****************************************************************************
void PPhoreQuery::parse_exc_vol_line( istringstream &ss ,
				      vector<BTV> &exc_vols ) {

  float x , y , z , r;
  ss >> x >> y >> z >> r;
  if( ss.fail() )
    throw PPhoreQueryFileReadError( filename_ , ss.str() );

  exc_vols.push_back( boost::make_tuple( x , y , z , r ) );

}

//****************************************************************************
// make sure the query details are all consistent
void PPhoreQuery::verify_query_details() {

  if( point_types_.empty() ) {
    throw PPhoreQueryFileError( filename_ ,
				string( "No points in query." ) );
  }
  verify_coord_details();
  verify_distance_details();
  verify_angle_details();
  verify_torsion_details();

}

//****************************************************************************
void PPhoreQuery::verify_coord_details() {

  if( !coords_.empty() && coords_.size() != point_types_.size() ) {
    string msg( "Error - number of coordinates not equal to number of points." );
    throw PPhoreQueryFileReadError( msg , filename_ );
  }

  for( int i = 0 , is = coords_.size() ; i < is ; ++i ) {
    if( coords_[i].get<0>() < 1 ||
	coords_[i].get<0>() > static_cast<int>( point_types_.size() ) ) {
      ostringstream oss1;
      oss1 << "Coords " << boost::tuples::set_open( ' ' )
	   << boost::tuples::set_close( ' ' ) << coords_[i];
      ostringstream oss2;
      oss2 << filename_ << " : Point number is greater than number of points.";
      throw PPhoreQueryFileReadError( oss2.str() , oss1.str() );
    }
    --coords_[i].get<0>();
  }

}

//****************************************************************************
void PPhoreQuery::verify_distance_details() {

  for( int i = 0 , is = distances_.size() ; i < is ; ++i ) {
    if( distances_[i].get<0>() < 1 || distances_[i].get<1>() < 1 ||
	distances_[i].get<0>() > static_cast<int>( point_types_.size() ) ||
	distances_[i].get<1>() > static_cast<int>( point_types_.size() ) ) {
      ostringstream oss1;
      oss1 << "Distance " << boost::tuples::set_open( ' ' )
	   << boost::tuples::set_close( ' ' ) << distances_[i];
      ostringstream oss2;
      oss2 << filename_ << " : Point number is greater than number of points.";
      throw PPhoreQueryFileReadError( oss2.str() , oss1.str() );      
    }
    --distances_[i].get<0>();
    --distances_[i].get<1>();
  }

}

//****************************************************************************
void PPhoreQuery::verify_angle_details() {

  for( int i = 0 , is = angles_.size() ; i < is ; ++i ) {
    if( angles_[i].get<0>() < 1 || angles_[i].get<1>() < 1 ||
	angles_[i].get<2>() < 1 ||
	angles_[i].get<0>() > static_cast<int>( point_types_.size() ) ||
	angles_[i].get<1>() > static_cast<int>( point_types_.size() ) ||
	angles_[i].get<2>() > static_cast<int>( point_types_.size() ) ) {
      ostringstream oss1;
      oss1 << "Angle " << boost::tuples::set_open( ' ' )
	   << boost::tuples::set_close( ' ' ) << angles_[i];
      ostringstream oss2;
      oss2 << filename_ << " : Point number is greater than number of points.";
      throw PPhoreQueryFileReadError( oss2.str() , oss1.str() );      
    }
    --angles_[i].get<0>();
    --angles_[i].get<1>();
    --angles_[i].get<2>();
  }

}

//****************************************************************************
void PPhoreQuery::verify_torsion_details() {

  for( int i = 0 , is = torsions_.size() ; i < is ; ++i ) {
    if( torsions_[i].get<0>() < 1 || torsions_[i].get<1>() < 1 ||
	torsions_[i].get<2>() < 1 || torsions_[i].get<3>() < 1 ||
	torsions_[i].get<0>() > static_cast<int>( point_types_.size() ) ||
	torsions_[i].get<1>() > static_cast<int>( point_types_.size() ) ||
	torsions_[i].get<2>() > static_cast<int>( point_types_.size() ) ||
	torsions_[i].get<3>() > static_cast<int>( point_types_.size() ) ) {
      ostringstream oss1;
      oss1 << "Torsion " << boost::tuples::set_open( ' ' )
	   << boost::tuples::set_close( ' ' ) << torsions_[i];
      ostringstream oss2;
      oss2 << filename_ << " : Point number is greater than number of points.";
      throw PPhoreQueryFileReadError( oss2.str() , oss1.str() );      
    }
    --torsions_[i].get<0>();
    --torsions_[i].get<1>();
    --torsions_[i].get<2>();
    --torsions_[i].get<3>();
  }

}

// ****************************************************************************
void PPhoreQuery::assemble_overlay_coords() {

  if( coords_.empty() )
    return;

  overlay_cds_ = new float[3 * coords_.size()];
  for( int i = 0 , j = 0 , is = coords_.size() ; i < is ; ++i ) {
    overlay_cds_[j++] = coords_[i].get<1>();
    overlay_cds_[j++] = coords_[i].get<2>();
    overlay_cds_[j++] = coords_[i].get<3>();
  }
      
}

// ****************************************************************************
void PPhoreQuery::assemble_distance_constraints() {

  DACLIB::make_square_matrix( sq_dist_vals_ , point_types_.size() );
  sq_dist_vals_[0][0] = 0.0F;
  for( int i = 1 , is = point_types_.size() ; i < is ; ++i ) {
    sq_dist_vals_[i][i] = 0.0F;
    for( int j = i + 1 , js = point_types_.size() ; j < js ; ++j ) {
      sq_dist_vals_[i][j] = -numeric_limits<float>::max();
      sq_dist_vals_[j][i] = numeric_limits<float>::max();
    }
  }

  for( int i = 0 , is = distances_.size() ; i < is ; ++i ) {
    int j = distances_[i].get<0>();
    int k = distances_[i].get<1>();
    sq_dist_vals_[j][k] = DACLIB::square( distances_[i].get<3>() );
    sq_dist_vals_[k][j] = DACLIB::square( distances_[i].get<2>() );
  }

#ifdef NOTYET
  cout << "SQ DISTANCE CONSTRAINTS" << endl;
  for( int i = 0 , is = point_types_.size() ; i < is ; ++i ) {
    for( int j = 0 ; j < is ; ++j ) {
      cout << sq_dist_vals_[i][j] << " ";
    }
    cout << endl;
  }
#endif

}

// ****************************************************************************
float **PPhoreQuery::make_squared_distance_matrix( vector<SinglePPhoreSite *> &target_sites ) {

  float **ret_mat;
  const int ns = target_sites.size();
  DACLIB::make_square_matrix( ret_mat , ns );

  for( int i = 0 ; i < ns - 1 ; ++i ) {
    ret_mat[i][i] = 0.0F;
    for( int j = i + 1 ; j < ns ; ++j ) {
      ret_mat[i][j] = ret_mat[j][i] =
          target_sites[i]->square_distance( target_sites[j]->coords() );
    }
  }
  ret_mat[ns - 1][ns - 1] = 0.0F;

#ifdef NOTYET
  cout << "SQ DISTANCE MATRIX" << endl;
  for( int ii = 0 ; ii < ns ; ++ii ) {
    for( int jj = 0 ; jj < ns ; ++jj ) {
      cout << ret_mat[ii][jj] << " ";
    }
    cout << endl;
  }
#endif

  return ret_mat;
}

// ****************************************************************************
// initial correspondence matrix m0. m0[i][j] is 1 if point_types_[i] is the
// same as target_sites[j]->type_string().
short **PPhoreQuery::form_m0( vector<SinglePPhoreSite *> &target_sites ) {

  const int npt = point_types_.size();
  const int nts = target_sites.size();
  short **m0 = DACLIB::make_2d_matrix<short>( npt , nts );
  fill_n( m0[0] , npt * nts , short( 0 ) );

  for( int i = 0 ; i < nts ; ++i ) {
    for( int j = 0 ; j < npt ; ++j ) {
      if( point_types_[j] == target_sites[i]->get_type_string() ) {
	m0[j][i] = 1;
      }
    }
  }

#ifdef NOTYET
  cout << "m0" << endl;
  for( int i = 0 ; i < npt ; i++ ) {
    for( int j = 0 ; j < nts ; j++ )
      cout << m0[i][j] << " ";
    cout << endl;
  }
#endif

  return m0;

}

// ****************************************************************************
// decide whether m contains a possible solution. If a row is all 0, then there
// isn't. If there are more columns all zero than
// num_rows - num_cols there also isn't, as there aren't enough
// target points to match all the query points. Targets are in the columns (there
// what we trying to match this PPhoreQuery with)
// 
bool PPhoreQuery::solution_possible( short **m , int num_rows ,
				     int num_cols ) {

  for( int i = 0 ; i < num_rows ; ++i ) {
    bool all_zero = true;
    for( int j = 0 ; j < num_cols ; ++j ) {
      if( m[i][j] ) {
	all_zero = false;
	break;
      }
    }
    if( all_zero )
      return false;
  }

  int zero_col_count = 0;
  for( int i = 0 ; i < num_cols ; ++i ) {
    bool zero_col = true;
    for( int j = 0 ; j < num_rows ; ++j ) {
      if( m[j][i] ) {
	zero_col = false;
	break;
      }
    }
    if( zero_col )
      ++zero_col_count;
    if( zero_col_count > num_cols - num_rows )
      return false;
  }

  return true;

}

// *****************************************************************************
// use an Ullmann-type procedure to find the hits
void PPhoreQuery::extend_subgraph( int row_num , short **m , const int num_rows ,
				   const int num_cols , const float **sq_dist_mat ,
				   vector<vector<int> > &hit_points ) {

  short **this_m = DACLIB::make_2d_matrix<short>( num_rows , num_cols );

  // at this level of recursion, we're investigating which query points in
  // columns (target sites) can match row row_num (point_types_[row_num]) in
  // the query pphore.

  for( int i = 0 ; i < num_cols ; ++i ) {
    if( m[row_num][i] ) {
      // this site in target could match point_types_[row_num]. Take a copy
      // of the current m and refine it to take account of this possible
      // match
      copy( m[0] , m[0] + num_rows * num_cols , this_m[0] );
      if( !refine_m( num_cols , sq_dist_mat , row_num , i , this_m ) )
	continue;

      if( check_for_solution( this_m , num_rows , num_cols ) ) {
	store_solution( this_m , num_rows , num_cols , hit_points );
	continue;
      }
      if( solution_possible( this_m , num_rows , num_cols ) &&
	  row_num < num_rows - 1 )
	extend_subgraph( row_num + 1 , this_m , num_rows , num_cols ,
			 sq_dist_mat , hit_points );
    }
  }

  DACLIB::destroy_2d_matrix( this_m );

}

// ***************************************************************************
// Ullmann's refinement process :
// refine the matrix until all possible incompatabilities have been removed
// or no solution can result.
bool PPhoreQuery::refine_m( const int num_cols , const float **sq_dist_mat ,
			    const int row_num , const int col_num , short **m ) {

  const int num_rows = point_types_.size();

  // cross everything but m[row_num][col_num] out of row row_num, col col_num
  fill( m[row_num] , m[row_num] + num_cols , 0 );
  for( int i = 0 ; i < num_rows ; ++i )
    m[i][col_num] = 0;
  m[row_num][col_num] = 1;

#ifdef NOTYET
  cout << "Refining for " << row_num << " , " << col_num << endl;
  for( int i = 0 ; i < num_rows ; ++i ) {
    for( int j = 0 ; j < num_cols ; ++j )
      cout << m[i][j] << " ";
    cout << endl;
  }
#endif

  bool sol_poss = true;

  while( 1 ) {
#ifdef NOTYET
    cout << "NEXT ITER" << endl;
    cout << "Refining" << endl;
    for( int i = 0 ; i < num_rows ; ++i ) {
      for( int j = 0 ; j < num_cols ; ++j )
	cout << m[i][j] << " ";
      cout << endl;
    }
#endif
    bool change = false;
    for( int i = 0 ; i < num_rows ; ++i ) {
      for( int j = 0 ; j < num_cols ; ++j ) {
	if( m[i][j] ) {
#ifdef NOTYET
	  cout << "checking q " << i << " matches t " << j << endl;
#endif
	  // j in the target could correspond to i in the query. Check this by
	  // ensuring that for all query contraints involving i, there is at
	  // least one match in the target. If there isn't, then i and j can't
	  // match, so set m[i][j] to 0.
	  bool early_k_break = false;
	  for( int k = 0 ; k < num_rows ; ++k ) {
	    if( k == i )
	      continue;
	    float upper_c , lower_c;
	    if( i < k ) {
	      upper_c = sq_dist_vals_[i][k];
	      lower_c = sq_dist_vals_[k][i];
	    } else {
	      upper_c = sq_dist_vals_[k][i];
	      lower_c = sq_dist_vals_[i][k];
	    }
#ifdef NOTYET
	    cout << i << " -> " << k << " contraints : "
		 << lower_c << " -> " << upper_c << endl;
#endif
	    bool early_l_break = false;
	    for( int l = 0 ; l < num_cols ; l++ ) {
	      // we know that i and j are the same type because m[i][j]
	      // wouldn't be set if they weren't. Check that k and l are also
	      // also a possible match (set in m) and that the constraint is
	      // met
	      if( m[k][l] ) {
#ifdef NOTYET
		cout << j << " --> " << l << " dist : " << sq_dist_mat[j][l] << endl;
#endif
		// check the distance between j and l is within range for i,k
		if( sq_dist_mat[j][l] > lower_c && sq_dist_mat[j][l] < upper_c ) {
#ifdef NOTYET
		  cout << "query " << i << " could match target " << j << endl;
#endif
		  early_l_break = true;
		  break;
		}
	      }
	    }
	    if( !early_l_break ) {
#ifdef NOTYET
	      cout << "no j->l matches " << i << " -> " << k << endl;
#endif
	      early_k_break = true;
	      break;
	    }
	  }
	  // not all the i constraints can be met by a j, so j can't match i, and
	  // we know this because we broke out of k loop early
	  if( early_k_break ) {
	    m[i][j] = 0;
	    change = true;
	  }
	}
      }
    }
#ifdef NOTYET
    cout << "Becomes" << endl;
    for( int i = 0 ; i < num_rows ; ++i ) {
      for( int j = 0 ; j < num_cols ; ++j )
	cout << m[i][j] << " ";
      cout << endl;
    }
#endif
    sol_poss = solution_possible( m , num_rows , num_cols );
    if( !change || !sol_poss )
      break;

  }

#ifdef NOTYET
  cout << "Done Refining for " << row_num << " , " << col_num << endl;
  for( int i = 0 ; i < num_rows ; ++i ) {
    for( int j = 0 ; j < num_cols ; ++j )
      cout << m[i][j] << " ";
    cout << endl;
  }
#endif

  return sol_poss;

}

// ***************************************************************************
// it's a solution if all of the rows have 1 1.
bool PPhoreQuery::check_for_solution( short **m , const int num_rows ,
				      const int num_cols ) {

#ifdef NOTYET
  cout << "Checking solution : " << endl;
  for( int i = 0 ; i < num_rows ; ++i ) {
    for( int j = 0 ; j < num_cols ; ++j )
      cout << m[i][j] << " ";
    cout << endl;
  }
#endif

  for( int i = 0 ; i < num_rows ; ++i ) {
    bool found_a_one = false;
    for( int j = 0 ; j < num_cols ; ++j ) {
      if( m[i][j] ) {
	if( found_a_one )
	  return false;
	else
	  found_a_one = true;
      }
    }
  }

#ifdef NOTYET
  cout << "HIT" << endl;
#endif

  return true;

}

// ***************************************************************************
void PPhoreQuery::store_solution( short **m , const int num_rows ,
				  const int num_cols ,
				  vector<vector<int> > &hit_points ) {

#ifdef NOTYET
  cout << "Storing solution : " << endl;
  for( int i = 0 ; i < num_rows ; ++i ) {
    for( int j = 0 ; j < num_cols ; ++j )
      cout << m[i][j] << " ";
    cout << endl;
  }
#endif

  hit_points.push_back( vector<int>( num_rows , -1 ) );
  for( int i = 0 ; i < num_rows ; ++i ) {
    for( int j = 0 ; j < num_cols ; ++j ) {
      if( m[i][j] ) {
	hit_points.back()[i] = j;
	break;
      }
    }
  }

#ifdef NOTYET
  cout << "HIT :: ";
  for( int i = 0 , is = hit_points.back().size() ; i < is ; ++i )
    cout << hit_points.back()[i] << " ";
  cout << endl;
#endif

  // if this is a repeat of one already found dump it.
  for( int i = 0 , is = hit_points.size() - 1 ; i < is ; ++i ) {
    if( hit_points[i] == hit_points.back() ) {
      hit_points.pop_back();
      break;
    }
  }

}

// ***************************************************************************
// check that the hit identified by ullman satisfies the angle contraints
bool PPhoreQuery::angles_ok( vector<SinglePPhoreSite *> &target_sites ,
			     const vector<int> &hit_points ,
			     vector<float> &hit_angles ) {

  hit_angles.clear();
  for( int i = 0 , is = angles_.size() ; i < is ; ++i ) {
    if( angles_[i].get<0>() != angles_[i].get<1>() &&
	angles_[i].get<1>() != angles_[i].get<2>() ) {
      SinglePPhoreSite *site1 = target_sites[hit_points[angles_[i].get<0>()]];
      SinglePPhoreSite *site2 = target_sites[hit_points[angles_[i].get<1>()]];
      SinglePPhoreSite *site3 = target_sites[hit_points[angles_[i].get<2>()]];
      float ang = DACLIB::angle( site1->coords() , site2->coords() ,
				 site3->coords() , true );
      if( ang < angles_[i].get<3>() || ang > angles_[i].get<4>() ) {
	hit_angles.clear();
	return false;
      }
      hit_angles.push_back( ang );
    } else {
      float ang;
      if( !angles_with_dirs_ok( angles_[i] , target_sites , hit_points ,
				ang ) ) {
	hit_angles.clear();
	return false;
      }
      hit_angles.push_back( ang );
    }
  }

  return true;

}

// ***************************************************************************
// check the angle with one of the end sites defined by the direction off the
// central site
bool PPhoreQuery::angles_with_dirs_ok( BTA &angle ,
				       vector<SinglePPhoreSite *> &target_sites ,
				       const vector<int> &hit_points ,
				       float &hit_ang ) {

  int p1 = hit_points[angle.get<0>()];
  int p2 = hit_points[angle.get<1>()];
  int p3 = hit_points[angle.get<2>()];

  SinglePPhoreSite *site1 = target_sites[p1];
  SinglePPhoreSite *site2 = target_sites[p2];
  SinglePPhoreSite *site3 = target_sites[p3];

  if( !site2->get_num_dirs() ) {
    return false; // but something's mucked up
  }

  if( p2 == p3 )
    swap( site1 , site3 );

  int num_dirs;
  double dirs[9] , cds[9];

  // direction end goes into cds[0-2], 2 sets of coords into cds[3-8]
  cds[3] = site2->coords()[0];
  cds[4] = site2->coords()[1];
  cds[5] = site2->coords()[2];
  cds[6] = site3->coords()[0];
  cds[7] = site3->coords()[1];
  cds[8] = site3->coords()[2];

  if( GtplDefs::RING_NORMAL == site2->direction_type( 0 ) ) {
    dirs[0] = site2->direction( 0 )[0];
    dirs[1] = site2->direction( 0 )[1];
    dirs[2] = site2->direction( 0 )[2];
    dirs[3] = -dirs[0];
    dirs[4] = -dirs[1];
    dirs[5] = -dirs[2];
    num_dirs = 2;
  } else {
    num_dirs = site2->get_num_dirs();
    for( int i = 0 , j = 0 ; i < num_dirs ; ++i ) {
      dirs[j++] = site2->direction( i )[0];
      dirs[j++] = site2->direction( i )[1];
      dirs[j++] = site2->direction( i )[2];
    }
  }

  for( int i = 0 ; i < num_dirs ; ++i ) {
    cds[0] = cds[3] + dirs[3 * i + 0];
    cds[1] = cds[4] + dirs[3 * i + 1];
    cds[2] = cds[5] + dirs[3 * i + 2];
    float ang = DACLIB::angle( cds , cds + 3 , cds + 6 , true );
    if( ang >= angle.get<3>() && ang <= angle.get<4>() ) {
      hit_ang = ang;
      return true;
    }
  }

  return false;

}

// ***************************************************************************
// check that the hit identified by ullman satisfies the torsion contraints
bool PPhoreQuery::torsions_ok( vector<SinglePPhoreSite *> &target_sites ,
			       const vector<int> &hit_points ,
			       vector<float> &hit_torsions ) {

  for( int i = 0 , is = torsions_.size() ; i < is ; ++i ) {
    if( torsions_[i].get<0>() == torsions_[i].get<1>() ||
	torsions_[i].get<2>() == torsions_[i].get<3>() ) {
      float tors;
      if( !torsions_with_dirs_ok( torsions_[i] , target_sites , hit_points ,
				  tors ) ) {
	hit_torsions.clear();
	return false;
      }
      hit_torsions.push_back( tors );
    } else {
      SinglePPhoreSite *site1 = target_sites[hit_points[torsions_[i].get<0>()]];
      SinglePPhoreSite *site2 = target_sites[hit_points[torsions_[i].get<1>()]];
      SinglePPhoreSite *site3 = target_sites[hit_points[torsions_[i].get<2>()]];
      SinglePPhoreSite *site4 = target_sites[hit_points[torsions_[i].get<3>()]];
      float tors = DACLIB::torsion( site1->coords() , site2->coords() ,
				    site3->coords() , site4->coords() , true );
      if( tors < torsions_[i].get<4>() || tors > torsions_[i].get<5>() ) {
	hit_torsions.clear();
	return false;
      } else {
	hit_torsions.push_back( tors );
      }
    }
  }

  return true;

}

// ***************************************************************************
bool PPhoreQuery::torsions_with_dirs_ok( BTT &torsion ,
					 vector<SinglePPhoreSite *> &target_sites ,
					 const vector<int> &hit_points ,
					 float &hit_tors ) {

  int p1 = hit_points[torsion.get<0>()];
  int p2 = hit_points[torsion.get<1>()];
  int p3 = hit_points[torsion.get<2>()];
  int p4 = hit_points[torsion.get<3>()];

  SinglePPhoreSite *site1 = target_sites[p1];
  SinglePPhoreSite *site2 = target_sites[p2];
  SinglePPhoreSite *site3 = target_sites[p3];
  SinglePPhoreSite *site4 = target_sites[p4];

  int num_dirs1 , num_dirs2;
  double dirs1[9] , dirs2[9] , cds[12];

  // points 2 and 3 always good
  cds[3] = site2->coords()[0];
  cds[4] = site2->coords()[1];
  cds[5] = site2->coords()[2];
  cds[6] = site3->coords()[0];
  cds[7] = site3->coords()[1];
  cds[8] = site3->coords()[2];

  if( p1 == p2 ) {
    if( GtplDefs::RING_NORMAL == site2->direction_type( 0 ) ) {
      dirs1[0] = site2->direction( 0 )[0];
      dirs1[1] = site2->direction( 0 )[1];
      dirs1[2] = site2->direction( 0 )[2];
      dirs1[3] = -dirs1[0];
      dirs1[4] = -dirs1[1];
      dirs1[5] = -dirs1[2];
      num_dirs1 = 2;
    } else {
      num_dirs1 = site2->get_num_dirs();
      for( int i = 0 , j = 0 ; i < num_dirs1 ; ++i ) {
	dirs1[j++] = site2->direction( i )[0];
	dirs1[j++] = site2->direction( i )[1];
	dirs1[j++] = site2->direction( i )[2];
      }
    }
  } else {
    // fake the directions
    num_dirs1 = 1;
    dirs1[0] = cds[3] - site1->coords()[0];
    dirs1[1] = cds[4] - site1->coords()[1];
    dirs1[2] = cds[5] - site1->coords()[2];
  }

  // same for p3 and p4
  if( p4 == p3 ) {
    if( GtplDefs::RING_NORMAL == site3->direction_type( 0 ) ) {
      dirs2[0] = site3->direction( 0 )[0];
      dirs2[1] = site3->direction( 0 )[1];
      dirs2[2] = site3->direction( 0 )[2];
      dirs2[3] = -dirs2[0];
      dirs2[4] = -dirs2[1];
      dirs2[5] = -dirs2[2];
      num_dirs2 = 2;
    } else {
      num_dirs2 = site3->get_num_dirs();
      for( int i = 0 , j = 0 ; i < num_dirs2 ; ++i ) {
	dirs2[j++] = site3->direction( i )[0];
	dirs2[j++] = site3->direction( i )[1];
	dirs2[j++] = site3->direction( i )[2];
      }
    }
  } else {
    // fake the directions
    num_dirs2 = 1;
    dirs2[0] = cds[6] - site4->coords()[0];
    dirs2[1] = cds[7] - site4->coords()[1];
    dirs2[2] = cds[8] - site4->coords()[2];
  }

  for( int i = 0 ; i < num_dirs1 ; ++i ) {
    cds[0] = cds[3] + dirs1[3 * i + 0];
    cds[1] = cds[4] + dirs1[3 * i + 1];
    cds[2] = cds[5] + dirs1[3 * i + 2];
    for( int j = 0 ; j < num_dirs2 ; ++j ) {
      cds[9] = cds[6] + dirs2[3 * j + 0];
      cds[10] = cds[7] + dirs2[3 * j + 1];
      cds[11] = cds[8] + dirs2[3 * j + 2];
      float tors = DACLIB::torsion( cds , cds + 3 , cds + 6 , cds + 9 , true );
      if( tors >= torsion.get<4>() && tors <= torsion.get<5>() ) {
	hit_tors = tors;
	return true;
      }
    }
  }

  return false;

}

// ***************************************************************************
// put the details into an already initialised PVM buffer
void PPhoreQuery::pack_into_pvm_buffer() {

  DACLIB::pack_strings_vector( point_types_ );

  int ints_to_go[5];
  float flts_to_go[5];

  int num_to_send = coords_.size();
  pvm_pkint( &num_to_send , 1 , 1 );
  for( int i = 0 ; i < num_to_send ; ++i ) {
    ints_to_go[0] = coords_[i].get<0>() + 1; // count from 1, like it does in    
    flts_to_go[0] = coords_[i].get<1>();     // input file. Makes unpacking easier.
    flts_to_go[1] = coords_[i].get<2>();
    flts_to_go[2] = coords_[i].get<3>();
    pvm_pkint( ints_to_go , 1 , 1 );
    pvm_pkfloat( flts_to_go , 3 , 1 );
  }

  num_to_send = distances_.size();
  pvm_pkint( &num_to_send, 1 , 1 );
  for( int i = 0 ; i < num_to_send ; ++i ) {
    ints_to_go[0] = distances_[i].get<0>() + 1;
    ints_to_go[1] = distances_[i].get<1>() + 1;
    flts_to_go[0] = distances_[i].get<2>();
    flts_to_go[1] = distances_[i].get<3>();
    pvm_pkint( ints_to_go , 2 , 1 );
    pvm_pkfloat( flts_to_go , 2 , 1 );
  }

  num_to_send = angles_.size();
  pvm_pkint( &num_to_send, 1 , 1 );
  for( int i = 0 ; i < num_to_send ; ++i ) {
    ints_to_go[0] = angles_[i].get<0>() + 1;
    ints_to_go[1] = angles_[i].get<1>() + 1;
    ints_to_go[2] = angles_[i].get<2>() + 1;
    flts_to_go[0] = angles_[i].get<3>();
    flts_to_go[1] = angles_[i].get<4>();
    pvm_pkint( ints_to_go , 3 , 1 );
    pvm_pkfloat( flts_to_go , 2 , 1 );
  }

  num_to_send = torsions_.size();
  pvm_pkint( &num_to_send, 1 , 1 );
  for( int i = 0 ; i < num_to_send ; ++i ) {
    ints_to_go[0] = torsions_[i].get<0>() + 1;
    ints_to_go[1] = torsions_[i].get<1>() + 1;
    ints_to_go[2] = torsions_[i].get<2>() + 1;
    ints_to_go[3] = torsions_[i].get<3>() + 1;
    flts_to_go[0] = torsions_[i].get<4>();
    flts_to_go[1] = torsions_[i].get<5>();
    pvm_pkint( ints_to_go , 4 , 1 );
    pvm_pkfloat( flts_to_go , 2 , 1 );
  }

  num_to_send = hard_exc_vols_.size();
  pvm_pkint( &num_to_send, 1 , 1 );
  for( int i = 0 ; i < num_to_send ; ++i ) {
    flts_to_go[0] = hard_exc_vols_[i].get<0>();
    flts_to_go[1] = hard_exc_vols_[i].get<1>();
    flts_to_go[2] = hard_exc_vols_[i].get<2>();
    flts_to_go[3] = hard_exc_vols_[i].get<3>();
    pvm_pkfloat( flts_to_go , 4 , 1 );
  }

  num_to_send = soft_exc_vols_.size();
  pvm_pkint( &num_to_send, 1 , 1 );
  for( int i = 0 ; i < num_to_send ; ++i ) {
    flts_to_go[0] = soft_exc_vols_[i].get<0>();
    flts_to_go[1] = soft_exc_vols_[i].get<1>();
    flts_to_go[2] = soft_exc_vols_[i].get<2>();
    flts_to_go[3] = soft_exc_vols_[i].get<3>();
    pvm_pkfloat( flts_to_go , 4 , 1 );
  }

}

// ***************************************************************************
void PPhoreQuery::unpack_from_pvm_buffer() {

  clear_data();
  DACLIB::unpack_strings_vector( point_types_ );

  int ints_in[5];
  float flts_in[5];
  int num_to_rec;

  pvm_upkint( &num_to_rec , 1 , 1 );
  for( int i = 0 ; i < num_to_rec ; ++i ) {
    pvm_upkint( ints_in , 1 , 1 );
    pvm_upkfloat( flts_in , 3 , 1 );
    coords_.push_back( boost::make_tuple( ints_in[0] , flts_in[0] ,
					  flts_in[1] , flts_in[2] ) );
  }

  pvm_upkint( &num_to_rec , 1 , 1 );
  for( int i = 0 ; i < num_to_rec ; ++i ) {
    pvm_upkint( ints_in , 2 , 1 );
    pvm_upkfloat( flts_in , 2 , 1 );
    distances_.push_back( boost::make_tuple( ints_in[0] , ints_in[1] ,
					     flts_in[0] , flts_in[1] ) );
  }

  pvm_upkint( &num_to_rec , 1 , 1 );
  for( int i = 0 ; i < num_to_rec ; ++i ) {
    pvm_upkint( ints_in , 3 , 1 );
    pvm_upkfloat( flts_in , 2 , 1 );
    angles_.push_back( boost::make_tuple( ints_in[0] , ints_in[1] ,
					  ints_in[2] ,
					  flts_in[0] , flts_in[1] ) );
  }

  pvm_upkint( &num_to_rec , 1 , 1 );
  for( int i = 0 ; i < num_to_rec ; ++i ) {
    pvm_upkint( ints_in , 4 , 1 );
    pvm_upkfloat( flts_in , 2 , 1 );
    torsions_.push_back( boost::make_tuple( ints_in[0] , ints_in[1] ,
					    ints_in[2] , ints_in[3] ,
					    flts_in[0] , flts_in[1] ) );
  }

  pvm_upkint( &num_to_rec , 1 , 1 );
  for( int i = 0 ; i < num_to_rec ; ++i ) {
    pvm_upkfloat( flts_in , 4 , 1 );
    hard_exc_vols_.push_back( boost::make_tuple( flts_in[0] , flts_in[1] ,
						 flts_in[2] , flts_in[3] ) );
  }

  pvm_upkint( &num_to_rec , 1 , 1 );
  for( int i = 0 ; i < num_to_rec ; ++i ) {
    pvm_upkfloat( flts_in , 4 , 1 );
    soft_exc_vols_.push_back( boost::make_tuple( flts_in[0] , flts_in[1] ,
						 flts_in[2] , flts_in[3] ) );
  }

  verify_query_details();
  assemble_overlay_coords();
  assemble_distance_constraints();

}
