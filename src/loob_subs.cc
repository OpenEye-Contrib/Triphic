//
// file loob_subs.cc
// David Cosgrove
// AstraZeneca
// 2nd October 2007
//
// Where all the hard work is done for loob. Most of the guts pulled out of the
// previous version pretty much verbatim.

#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <string>
#include <vector>

#include <oechem.h>

#include <boost/bind.hpp>
#include <boost/regex.hpp>
#include <boost/scoped_ptr.hpp>

#include <cstdio> // for uncompressed FP I/O
#include <zlib.h> // for compressed FP I/O

#include "stddefs.H"
#include "FileExceptions.H"
#include "LoobSettings.H"
#include "PharmPoint.H"
#include "SinglePPhoreSite.H"
#include "SMARTSExceptions.H"

using namespace std;
using namespace OEChem;

namespace DACLIB {
  // in eponymous file - throws DACLIB::FileReadOpenError or
  // DACLIB::SMARTSFileError
  void read_smarts_file( const string &smarts_file ,
			 vector<pair<string,string> > &input_smarts ,
			 vector<pair<string,string> > &smarts_sub_defn );
  // in eponymous file. Throws DACLIB::SMARTSDefnError if the mood takes it.
  void make_pphore_sites( OEMol &mol , PharmPoint &pharm_points ,
			  const vector<pair<string,string> > &input_smarts ,
			  vector<pair<string,string> > &smarts_sub_defn ,
			  vector<vector<SinglePPhoreSite *> > &pharm_sites );
  void split_filename( const string &filename ,
		       string &file_root , string &file_ext );
}

const static float MIN_DIST = 1.0F; /* minimum distance between two points i.e.
				       shortest length that a side can have. So
				       that points e.g. acceptor/donor that can
				       be both on the same atom don't form two
				       corners of a triangle. */

// ***************************************************************************
void read_smarts_file( const string &smarts_file ,
		       vector<pair<string,string> > &input_smarts ,
		       vector<pair<string,string> > &smarts_sub_defn ) {

  try {
    DACLIB:: read_smarts_file( smarts_file , input_smarts , smarts_sub_defn );
  } catch( DACLIB::SMARTSFileError &e ) {
    cerr << e.what() << endl;
    exit( 1 );
  } catch( DACLIB::SMARTSSubDefnError &e ) {
    cerr << e.what() << endl;
    exit( 1 );
  } catch( DACLIB::FileReadOpenError &e ) {
    cerr << e.what() << endl;
    exit( 1 );
  }

}

// ********************************************************************
void read_subset_file( const string &subset_file ,
		       vector<string> &subset_names ) {

  ifstream ifs( subset_file.c_str() );
  if( !ifs.good() ) {
    throw DACLIB::FileReadOpenError( subset_file.c_str() ) ;
  }

  string tmp;
  while( 1 ) {
    ifs >> tmp;
    if( ifs.eof() || !ifs.good() ) {
      break;
    }
    subset_names.push_back( tmp );
  }

  if( subset_names.empty() ) {
    cerr << "Subset file " << subset_file << " was empty, so it's a short job."
	 << endl;
    exit( 1 );
  }

  sort( subset_names.begin() , subset_names.end() );

}

// **************************************************************************
// check the triangle passed in is valid. Returns true if so. d1, d2 and d3 are
// the indices into dist_bounds of the three side lengths
bool valid_triangle( int d1 , int d2 , int d3 ,
		     const vector<float> &dist_bounds ) {

  int sides[3] = { d1 , d2 , d3 };

  // find the long and short sides, making sure they're different
  int long_side = 0 , mid_side = 1 , short_side = 2;
  if( sides[1] > sides[long_side] ) long_side = 1;
  if( sides[2] > sides[long_side] ) long_side = 2;
  if( sides[0] < sides[short_side] ) short_side = 0;
  if( sides[1] < sides[short_side] ) short_side = 1;

  // and the mid side
  if( 0 == long_side && 1 == short_side ) mid_side = 2;
  else if( 1 == long_side && 0 == short_side ) mid_side = 2;
  else if( 0 == long_side && 2 == short_side ) mid_side = 1;
  else if( 2 == long_side && 0 == short_side ) mid_side = 1;
  else if( 1 == long_side && 2 == short_side ) mid_side = 0;
  else if( 2 == long_side && 1 == short_side ) mid_side = 0;
  
  // if the lower bound of the long side is more than the sum of the
  // upper bounds of the other two sides, then the triangle is impossible
  // and so is the polygon
  float long_dist =
    0 == sides[long_side] ? MIN_DIST : dist_bounds[sides[long_side] - 1];
  float short_dist = dist_bounds[sides[short_side]] +
    dist_bounds[sides[mid_side]];

  return bool( long_dist <= short_dist );

}

// ***************************************************************************
void make_2point_fp_map( map<string,int> &fp_map , vector<string> &bit_decoding ,
			 int num_point_types , LoobSettings &loob_settings ) {

  string this_label( "   " );
  int bit_count = fp_map.size();
  const int num_dists = loob_settings.dist_bounds_.size();

  for( int f1 = 0 ; f1 < num_point_types ; ++f1 ) {
    for( int f2 = 0 ; f2 < num_point_types ; ++f2 ) {
      this_label[0] = 'a' + f1;
      this_label[1] = 'a' + f2;
      if( this_label[0] > this_label[1] ) {
	swap( this_label[0] , this_label[1] );
      }
      for( int d1 = 0 ; d1 < num_dists ; ++d1 ) {
	this_label[2] = 'a' + d1;
	if( fp_map.insert( make_pair( this_label , bit_count ) ).second ) {
	  bit_decoding.push_back( this_label );
	  ++bit_count;
	}
      }
    }
  }

}

// ***************************************************************************
// make all possible 3 point codes, not including invalid triangles (by the
// triangle rule) and redundant distance combinations (by the encoding rules)
void make_3point_fp_map( map<string,int> &fp_map , vector<string> &bit_decoding ,
			 int num_point_types , LoobSettings &loob_settings ) {

  string   this_label( "      " );
  int bit_count = fp_map.size();
  const int num_dists = loob_settings.dist_bounds_.size();

  for( int d3 = 0 ; d3 < num_dists ; ++d3 ) {
    this_label[2] = 'a' + d3;
    for( int d2 = d3 ; d2 < num_dists ; ++d2 ) {
      this_label[1] = 'a' + d2;
      for( int d1 = d2 ; d1 < num_dists ; ++d1 ) {
	if( !valid_triangle( d1 , d2 , d3 , loob_settings.dist_bounds_ ) ) {
	  continue;
	}
	this_label[0] = 'a' + d1;
	for( int f1 = 0 ; f1 < num_point_types ; ++f1 ) {
	  this_label[3] = 'a' + f1;
	  for( int f2 = 0 ; f2 < num_point_types ; ++f2 ) {
	    this_label[4] = 'a' + f2;
	    for( int f3 = 0 ; f3 < num_point_types ; ++f3 ) {
	      this_label[5] = 'a' + f3;
	      if( fp_map.insert( make_pair( this_label , bit_count ) ).second ) {
		bit_decoding.push_back( this_label );
		++bit_count;
	      }
	    }
	  }
	}
      }
    }	      
  }

}

// ***************************************************************************
void make_4point_fp_map( map<string,int> &fp_map , vector<string> &bit_decoding ,
			 int num_point_types , LoobSettings &loob_settings ) {

  string   this_label( "          " );
  int bit_count = fp_map.size();
  const int num_dists = loob_settings.dist_bounds_.size();

  for( int d3 = 0 ; d3 < num_dists ; ++d3 ) {
    for( int d2 = d3 ; d2 < num_dists ; ++d2 ) {
      for( int d1 = d2 ; d1 < num_dists ; ++d1 ) {
	if( !valid_triangle( d1 , d2 , d3 , loob_settings.dist_bounds_ ) ) {
	  continue;
	}
	// d1 is the longest side, so d4 to d6 can't exceed it
	for( int d4 = 0 ; d4 <= d1 ; ++d4 ) {
	  for( int d5 = 0 ; d5 <= d1 ; ++d5 ) {
	    if( !valid_triangle( d1 , d4 , d5 , loob_settings.dist_bounds_ ) ) {
	      continue;
	    }
	    for( int d6 = 0 ; d6 <= d1 ; ++d6 ) {
	      if( !valid_triangle( d2 , d4 , d6 , loob_settings.dist_bounds_ ) ) {
		continue;
	      }
	      if( !valid_triangle( d3 , d5 , d6 , loob_settings.dist_bounds_ ) ) {
		continue;
	      }
	      for( int f1 = 0 ; f1 < num_point_types ; ++f1 ) {
		for( int f2 = 0 ; f2 < num_point_types ; ++f2 ) {
		  for( int f3 = 0 ; f3 < num_point_types ; ++f3 ) {
		    for( int f4 = 0 ; f4 < num_point_types ; ++f4 ) {
		      this_label[9] = 'a' + f4;
		      this_label[8] = 'a' + f3;
		      this_label[7] = 'a' + f2;
		      this_label[6] = 'a' + f1;
		      this_label[5] = 'a' + d6;
		      this_label[4] = 'a' + d5;
		      this_label[3] = 'a' + d4;
		      this_label[2] = 'a' + d3;
		      this_label[1] = 'a' + d2;
		      this_label[0] = 'a' + d1;
		      if( fp_map.insert( make_pair( this_label ,
						    bit_count ) ).second ) {
			bit_decoding.push_back( this_label );
			++bit_count;
		      }
		      if( loob_settings.chiral_fps_ ) {
			this_label[8] = 'a' + f1;
			this_label[6] = 'a' + f3;
			this_label[5] = 'a' + d4;
			this_label[3] = 'a' + d6;
			this_label[2] = 'a' + d1;
			this_label[0] = 'a' + d3;
			if( fp_map.insert( make_pair( this_label ,
						      bit_count ) ).second ) {
			  bit_decoding.push_back( this_label );
			  ++bit_count;
			}
		      }
		    }
		  }
		}
	      }		      
	    }
	  }
	}
      }
    }	      
  }

}

// ***************************************************************************
void make_fp_map( map<string,int> &fp_map , vector<string> &bit_decoding ,
		  int num_point_types , LoobSettings &loob_settings ) {

  if( loob_settings.pairs_ ) {
    make_2point_fp_map( fp_map , bit_decoding , num_point_types ,
			loob_settings );
  }
  if( loob_settings.triplets_ ) {
    make_3point_fp_map( fp_map , bit_decoding , num_point_types ,
			loob_settings );
  }
  if( loob_settings.quadruplets_ ) {
    make_4point_fp_map( fp_map , bit_decoding , num_point_types ,
			loob_settings );
  }

  cout << num_point_types << " point types over "
       << loob_settings.dist_bounds_.size()
       << " distances, giving fingerprint size : "
       << fp_map.size() << " bits." << endl;

}

// ***************************************************************************
int dist_bin( float sq_dist , const vector<float> &sq_dist_bounds ) {

  if( sq_dist > sq_dist_bounds.back() || sq_dist < MIN_DIST ) {
    return sq_dist_bounds.size();
  }

  for( int i = sq_dist_bounds.size() - 1 ; i >= 0 ; --i ) {
    if( sq_dist > sq_dist_bounds[i] ) {
      return i + 1;
    }
  }
    
  return 0;

}

// ****************************************************************************
// see if the 4 points passed in form a negatively chiral quadruplet according
// to the rules of Abrahamian et al.
bool negative_chiral( int f1 , int f2 , int f3 , int f4 ,
		      const vector<SinglePPhoreSite *> &sites ) {

  double f2_f1[3] , f2_f3[3], f2_f4[3] , x[3];

  DACLIB::join_vector( sites[f2]->coords() , sites[f1]->coords() , f2_f1 );
  DACLIB::join_vector( sites[f2]->coords() , sites[f3]->coords() , f2_f3 );
  DACLIB::join_vector( sites[f2]->coords() , sites[f4]->coords() , f2_f4 );

  DACLIB::cross_product( f2_f1 , f2_f3 , x );
  return( DACLIB::dot_product( x , f2_f4 ) <= 0.0 );

}

// ***************************************************************************
void make_sq_dist_matrix( vector<SinglePPhoreSite *> &sites ,
			  float **&sq_dist_mat ) {

  DACLIB::make_square_matrix( sq_dist_mat , sites.size() );
  int is = sites.size() - 1;
  for( int i = 0 ; i < is ; ++i ) {
    sq_dist_mat[i][i] = 0.0F;
    for( int j = i + 1 , js = is + 1 ; j < js ; ++j ) {
      sq_dist_mat[i][j] = sq_dist_mat[j][i] =
	sites[i]->square_distance( sites[j]->coords() );
    }
  }

  sq_dist_mat[is][is] = 0.0F;

}

// **************************************************************************
// take the triangle points given, distances and ends, and produce the order
// according to Abrahamian et al.
void order_triangle_points( int d[6][2] , int e[6][2] , int &f1 , int &f2 ,
			    int &f3 , int &d1 , int &d2 , int &d3 ) {

  int l , m , r[6];
  r[0] = 0; r[1] = 1; r[2] = 2; r[3] = 3; r[4] = 4; r[5] = 5;
  // rank the distances, using the second entry as a tie-breaker
  for( l = 0 ; l < 5 ; l++ ) {
    for( m = l + 1 ; m < 6 ; m++ ) {
      if( d[r[l]][0] < d[r[m]][0] ||
	  ( d[r[l]][0] == d[r[m]][0] && d[r[l]][1] > d[r[m]][1] ) ) {
	swap( r[l] , r[m] );
      }
    }
  }
  // f2 is the element of e[r[0]] and e[r[5]] that is shared
  if( e[r[0]][0] == e[r[5]][1] ) {
    f2 = e[r[0]][0];
  } else if( e[r[0]][0] == e[r[5]][0] ) {
    f2 = e[r[0]][0];
  } else if( e[r[0]][1] == e[r[5]][1] ) {
    f2 = e[r[0]][1];
  } else if( e[r[0]][1] == e[r[5]][0] ) {
    f2 = e[r[0]][1];
  }
  
  // f1 is at the other end of the longer side from f2
  for( l = 0 ; l < 6 ; l++ ) {
    if( e[r[l]][0] == f2 ) {
      f1 = e[r[l]][1];
      d1 = d[r[l]][0];
      break;
    }
  }
  // f3 is at the shorter
  for( l = 5 ; l >= 0 ; l-- ) {
    if( e[r[l]][0] == f2 ) {
      f3 = e[r[l]][1];
      d3 = d[r[l]][0];
      break;
    }
  }
  // d2 is the distance between f1 and f3
  for( l = 0 ; l < 6 ; l++ ) {
    if( e[l][0] == f1 && e[l][1] == f3 ) {
      d2 = d[l][0];
      break;
    }
  }

}

// ****************************************************************************
// rank the 4 triangles passed in using the rules of Abrahmian et al. for
// deciding which is the base triangle of the quadruplet.
void rank_triangles( int r[4] , int f1[4] , int f2[4] , int f3[4] ,
		     int fd1[4] , int fd2[4] , int fd3[4] ,
		     const vector<SinglePPhoreSite *> &sites ) {

  int      i , j;

  r[0] = 0; r[1] = 1; r[2] = 2; r[3] = 3;

  // the rules are:
  // first, order by longest side, which is now fd1
  // second, by shortest side fd3
  // third, by third side fd2
  // 4th, 5th and 6th by the points_list values of f1, f2 and f3
  // if that doesn't sort things out, they are equivalent, so it doesn't matter.
  for( i = 0 ; i < 3 ; i++ ) {
    for( j = i + 1 ; j < 4 ; j++ ) {
      if( fd1[r[i]] < fd1[r[j]] ) {
	swap( r[i] , r[j] );
      } else if( fd1[r[i]] == fd1[r[j]] ) {
	if( fd3[r[i]] > fd3[r[j]] ) {
	  swap( r[i] , r[j] );
	} else if( fd3[r[i]] == fd3[r[j]] ) {
	  if( fd2[r[i]] < fd2[r[j]] ) {
	    swap( r[i] , r[j] );
	  } else if( fd2[r[i]] == fd2[r[j]] ) {
	    if( sites[f1[r[i]]]->get_type_code() >
		sites[f1[r[j]]]->get_type_code() ) {
	      swap( r[i] , r[j] );
	    } else if( sites[f1[r[i]]]->get_type_code() ==
		       sites[f1[r[j]]]->get_type_code() ) {
	      if( sites[f2[r[i]]]->get_type_code() >
		  sites[f2[r[j]]]->get_type_code() ) {
		swap( r[i] , r[j] );
	      } else if( sites[f2[r[i]]]->get_type_code() ==
			 sites[f2[r[j]]]->get_type_code() ) {
		if( sites[f3[r[i]]]->get_type_code() >
		    sites[f3[r[j]]]->get_type_code() )
		  swap( r[i] , r[j] );
	      }
	    }
	  }
	}
      }
    }

  }

}

// ***************************************************************************
// make an ordered triangle from the i, j and k points
bool build_triangle( int d[6][2] , int e[6][2] , int &f1 , int &f2 , int &f3 ,
		     int &d1 , int &d2 , int &d3 , int i , int j , int k ,
		     float **sq_dist_mat , const vector<float> &sq_dist_bounds ,
		     const vector<SinglePPhoreSite *> &sites ) {

  const int num_dists = sq_dist_bounds.size();

  d[0][0] = dist_bin( sq_dist_mat[i][j] , sq_dist_bounds );
  if( num_dists == d[0][0] ) return false; // beyond range
  d[1][0] = dist_bin( sq_dist_mat[i][k] , sq_dist_bounds );
  if( num_dists == d[1][0] ) return false; // beyond range
  d[2][0] = d[0][0];
  d[3][0] = dist_bin( sq_dist_mat[j][k] , sq_dist_bounds );
  if( num_dists == d[3][0] ) return false; // beyond range
  d[4][0] = d[1][0];
  d[5][0] = d[3][0];

  d[0][1] = d[5][1] = sites[j]->get_type_code();
  d[1][1] = d[3][1] = sites[k]->get_type_code();
  d[2][1] = d[4][1] = sites[i]->get_type_code();

  e[0][0] = i; e[0][1] = j;
  e[1][0] = i; e[1][1] = k;
  e[2][0] = j; e[2][1] = i;
  e[3][0] = j; e[3][1] = k;
  e[4][0] = k; e[4][1] = i;
  e[5][0] = k; e[5][1] = j;
  
  order_triangle_points( d , e , f1 , f2 , f3 , d1 , d2 , d3 );

  return true;

}

// ***************************************************************************
void make_pair_fingerprint( vector<SinglePPhoreSite *> &sites ,
			    LoobSettings &loob_settings ,
			    float **&sq_dist_mat ,
			    set<string> &these_bits ) {

  const int num_dists = loob_settings.sq_dist_bounds_.size();

  if( sites.size() < 2 ) {
    return;
  }

  // we'll always be making it for pairs
  make_sq_dist_matrix( sites , sq_dist_mat );

  string bit_label( "   " );
  for( int i = 0 , is = sites.size() - 1 ; i < is ; ++i ) {
    char bl1 = 'a' + sites[i]->get_type_code();
    for( int j = i + 1 , js = sites.size() ; j < js ; ++j ) {
      char bl2 = 'a' + sites[j]->get_type_code();
      int db = dist_bin( sq_dist_mat[i][j] , loob_settings.sq_dist_bounds_ );
      if( db == num_dists ) {
	continue; // beyond range
      }
      if( bl1 > bl2 ) {
	bit_label[0] = bl2;
	bit_label[1] = bl1;
	bit_label[2] = 'a' + db;
      } else {
	bit_label[0] = bl1;
	bit_label[1] = bl2;
	bit_label[2] = 'a' + db;
      }
      these_bits.insert( bit_label );
    }
  }

}

// ***************************************************************************
void make_triplet_fingerprint( vector<SinglePPhoreSite *> &sites ,
			       LoobSettings &loob_settings ,
			       float **&sq_dist_mat ,
			       set<string> &these_bits ) {

  if( sites.size() < 3 ) {
    return;
  }

  if( !sq_dist_mat ) {
    // may have been made for pairs 
    make_sq_dist_matrix( sites , sq_dist_mat );
  }

  string bit_label( "      " );
  int      d[6][2] , e[6][2] , f1 , f2 , f3 , d1 , d2 , d3;

  for( int i = 0 , is = sites.size() - 2 ; i < is ; ++i ) {

    for( int j = i + 1 , js = sites.size() - 1 ; j < js ; ++j ) {

      for( int k = j + 1 , ks = sites.size() ; k < ks ; ++k ) {

	if( !build_triangle( d , e , f1 , f2 , f3 , d1 , d2 , d3 , i , j , k ,
			     sq_dist_mat , loob_settings.sq_dist_bounds_ ,
			     sites ) ) {
	  continue; // one of the distances was out of the bin ranges
	}

	// f[123] come back as the ordered i,j and k. Need to get the point type
	// for the bit label
	int ff1 = sites[f1]->get_type_code();
	int ff2 = sites[f2]->get_type_code();
	int ff3 = sites[f3]->get_type_code();

	bit_label[0] = 'a' + d1;
	bit_label[1] = 'a' + d2;
	bit_label[2] = 'a' + d3;
	bit_label[3] = 'a' + ff1;
	bit_label[4] = 'a' + ff2;
	bit_label[5] = 'a' + ff3;
	these_bits.insert( bit_label );

      }
      
    }

  }

}

// ***************************************************************************
void make_quadruplet_fingerprint( vector<SinglePPhoreSite *> &sites ,
				  LoobSettings &loob_settings ,
				  float **&sq_dist_mat ,
				  set<string> &these_bits ) {

  if( sites.size() < 4 ) {
    return;
  }

  if( !sq_dist_mat ) {
    // may have been made for pairs or triplets
    make_sq_dist_matrix( sites , sq_dist_mat );
  }

  int      d[4][6][2] , e[4][6][2] , f1[4] , f2[4] , f3[4] , fd1[4] , fd2[4] ,
           fd3[4] , f4 = -1 , d4 = -1 , d5 = -1 , d6 = -1;

  string bit_label( "          " );

  vector<float> &sq_dist_bounds = loob_settings.sq_dist_bounds_;

  for( int i = 0 , is = sites.size() - 3 ; i < is ; ++i ) {

    for( int j = i + 1 , js = sites.size() - 2 ; j < js ; ++j ) {

      for( int k = j + 1 , ks = sites.size() -1 ; k < ks ; ++k ) {

	for( int l = k + 1 , ls = sites.size() ; l < ls ; ++l ) {

	  // make the 4 triangles for this quadruplet, ordered by normal rules
	  // they will be from points i,j,k  i,j,l  i,k,l and j,k,l in that order
	  // the f's will contain the ordered point numbers (i,j,k or l) and the
	  // fd's will have the distance bins
	  if( !build_triangle( d[0] , e[0] , f1[0] , f2[0] , f3[0] , fd1[0] ,
			       fd2[0] , fd3[0] , i , j , k ,
			       sq_dist_mat , sq_dist_bounds , sites ) ) {
	    continue; // one of the distances was out of bounds
	  }
	  if( !build_triangle( d[1] , e[1] , f1[1] , f2[1] , f3[1] , fd1[1] ,
			       fd2[1] , fd3[1] , i , j , l ,
			       sq_dist_mat , sq_dist_bounds , sites ) ) {
	    continue; // one of the distances was out of bounds
	  }
	  if( !build_triangle( d[2] , e[2] , f1[2] , f2[2] , f3[2] , fd1[2] ,
			       fd2[2] , fd3[2] , i , k , l ,
			       sq_dist_mat , sq_dist_bounds , sites ) ) {
	    continue; // one of the distances was out of bounds
	  }
	  if( !build_triangle( d[3] , e[3] , f1[3] , f2[3] , f3[3] , fd1[3] ,
			       fd2[3] , fd3[3] , j , k , l ,
			       sq_dist_mat , sq_dist_bounds , sites ) ) {
	    continue; // one of the distances was out of bounds
	  }
	  // rank the triangles by the rules of Abrahamian et al., to find the
	  // base triangle
	  int r[4];
	  rank_triangles( r , f1 , f2 , f3 , fd1 , fd2 , fd3 , sites );

	  // decide what f4 is
	  switch( r[0] ) {
	    case 0 : f4 = l; break;
	    case 1 : f4 = k; break;
	    case 2 : f4 = j; break;
	    case 3 : f4 = i; break;
	  }
	  // and make d4, d5 and d6 which are to f4 respectively from f1, f2, f3
	  d4 = dist_bin( sq_dist_mat[f1[r[0]]][f4] ,
			 sq_dist_bounds ); // we know distance is ok
	  d5 = dist_bin( sq_dist_mat[f2[r[0]]][f4] ,
			 sq_dist_bounds ); // we know distance is ok
	  d6 = dist_bin( sq_dist_mat[f3[r[0]]][f4] ,
			 sq_dist_bounds ); // we know distance is ok
	  bit_label[0] = 'a' + fd1[r[0]];
	  bit_label[1] = 'a' + fd2[r[0]];
	  bit_label[2] = 'a' + fd3[r[0]];
	  bit_label[3] = 'a' + d4;
	  bit_label[4] = 'a' + d5;
	  bit_label[5] = 'a' + d6;
	  bit_label[6] = 'a' + sites[f1[r[0]]]->get_type_code();
	  bit_label[7] = 'a' + sites[f2[r[0]]]->get_type_code();
	  bit_label[8] = 'a' + sites[f3[r[0]]]->get_type_code();
	  bit_label[9] = 'a' + sites[f4]->get_type_code();
	  if( loob_settings.chiral_fps_ &&
	      negative_chiral( f1[r[0]] , f2[r[0]] , f3[r[0]] ,
			       f4 , sites ) ) {
	    swap( bit_label[0] , bit_label[2] );
	    swap( bit_label[3] , bit_label[5] );
	    swap( bit_label[6] , bit_label[8] );
	  }
	  these_bits.insert( bit_label );
	}

      }

    }

  }

}

// ***************************************************************************
// make the pharmacophore fingerprint from the sites
void make_pharm_fingerprint( vector<vector<SinglePPhoreSite *> > &sites ,
			     LoobSettings &loob_settings ,
			     set<string> &these_bits ) {

  for( int i = 0 , is = sites.size() ; i < is ; ++i ) {

    float **sq_dist_mat = 0;

    if( loob_settings.pairs_ ) {
      make_pair_fingerprint( sites[i] , loob_settings , sq_dist_mat ,
			     these_bits );
    }
    if( loob_settings.triplets_ ) {
      make_triplet_fingerprint( sites[i] , loob_settings , sq_dist_mat ,
				these_bits );
    }
    if( loob_settings.quadruplets_ ) {
      make_quadruplet_fingerprint( sites[i] , loob_settings ,
				   sq_dist_mat , these_bits );
    }

    DACLIB::destroy_square_matrix( sq_dist_mat );

  }

}

// ***************************************************************************
// convert the bit names into a sorted vector of set bits
void bits_strings_to_ints( const set<string> &these_bits ,
			   const map<string,int> &fp_map ,
			   vector<int> &these_set_bits ) {

  set<string>::const_iterator p , ps;
  map<string,int>::const_iterator q;

  for( p = these_bits.begin() , ps = these_bits.end() ; p != ps ; ++p ) {
    q = fp_map.find( *p );
    if( q == fp_map.end() ) {
      cerr << "AWOOGA - " << *p << " not in fingerprint map. This is such a bad"
	   << endl << "thing, I'm taking my ball home forthwith." << endl;
      exit( 1 );
    }
    these_set_bits.push_back( q->second );
  }

  sort( these_set_bits.begin() , these_set_bits.end() );

}

// ***************************************************************************
void process_molecule( vector<pair<string,string> > &input_smarts ,
		       vector<pair<string,string> > &smarts_sub_defn ,
		       PharmPoint &pharm_points ,
		       const vector<string> &mol_subset ,
		       LoobSettings &loob_settings , OEMol &mol ,
		       const map<string,int> &fp_map ,
		       vector<int> &these_set_bits ) {

  if( !mol_subset.empty() &&
      !binary_search( mol_subset.begin() , mol_subset.end() ,
		      mol.GetTitle() ) ) {
    return;
  }

  vector<vector<SinglePPhoreSite *> > sites;
  try {
    DACLIB::make_pphore_sites( mol , pharm_points , input_smarts ,
			       smarts_sub_defn , sites );
  } catch( DACLIB::SMARTSDefnError &e ) {
    cerr << e.what() << endl;
    exit( 1 );
  }

  set<string> these_bits;
  make_pharm_fingerprint( sites , loob_settings , these_bits );

  // convert the bit names into a sorted vector of set bits
  bits_strings_to_ints( these_bits , fp_map , these_set_bits );

  for( int i = 0 , is = sites.size() ; i < is ; ++i ) {
    for( int j = 0 , js = sites[i].size() ; j < js ; ++j ) {
      delete sites[i][j];
    }
  }

}

// *************************************************************************
void prepare_compressed_fp( int counts_cutoff , const vector<int> &fp_counts ,
			    const vector<int> &full_to_compressed ,
			    vector<int> &bits ) {

  vector<int>::iterator p , ps;
  for( p = bits.begin() , ps = bits.end() ; p < ps ; ++p ) {
    if( fp_counts[*p] < counts_cutoff ) {
      *p = -1; // flag for removal
    } else {
      *p = full_to_compressed[*p];
    }
  }

  bits.erase( remove( bits.begin() , bits.end() , -1 ) ,
	      bits.end() );

}

// *************************************************************************
// convert set_bits, containing the numbers of the bits to be set,
// to a vector or 0 or 1 showing which bits are set. Assumes that fp_bits
// is already big enough, but might have been used already
void bits_to_fp( const vector<int> &set_bits , vector<int> &fp_bits ) {

  fill( fp_bits.begin() , fp_bits.end() , 0 );

  vector<int>::const_iterator p , ps;
  for( p = set_bits.begin() , ps = set_bits.end() ; p < ps ; ++p ) {
    fp_bits[*p] = 1;
  }
  
}

// *************************************************************************
// if we're doing compress fps, make a conversion from the full bit number
// to the bit number in the compressed fingerprint
void map_full_to_compressed( const vector<int> &fp_counts , int counts_cutoff ,
			     vector<int> &full_to_compressed ) {

  int i = 0 , j = 0;
  vector<int>::const_iterator p , ps;
  full_to_compressed = vector<int>( fp_counts.size() , -1 );
  for( p = fp_counts.begin() , ps = fp_counts.end() ; p != ps ; ++p , ++i ) {
    if( *p >= counts_cutoff ) {
      full_to_compressed[i] = j++;
    }
  }

}
			  
// *************************************************************************
// make the composite fp, with a bit set for every position where
// fp_counts exceeds compressed_counts_cutoff_. Only sensible for uncompressed
// fps, as for compressed ones it should be all bits set.
void make_composite_fp( LoobSettings &loob_settings ,
			const vector<int> &fp_counts ,
			vector<int> &composite_fp_bits ) {

  if( !loob_settings.dont_compress_fps_ ) {
    int num_bits_to_set =
      count_if( fp_counts.begin() , fp_counts.end() ,
		boost::bind( greater<int>() , _1 ,
			     loob_settings.compressed_counts_cutoff_ - 1 ) );
    composite_fp_bits = vector<int>( num_bits_to_set , 1 );
    return;
  }

  composite_fp_bits = vector<int>( fp_counts.size() , 0 );
  vector<int>::const_iterator p , ps;
  vector<int>::iterator q;
  for( p = fp_counts.begin() , ps = fp_counts.end() ,
	 q = composite_fp_bits.begin() ; p != ps ; ++p , ++q ) {
    if( *p >= loob_settings.compressed_counts_cutoff_ ) {
      *q = 1;
    }
  }

}

// *************************************************************************
void write_ascii_fp( LoobSettings &loob_settings ,
		     const vector<pair<string,vector<int> > > &bits_set ,
		     const vector<int> &fp_counts ,
		     const vector<int> &full_to_compressed ,
		     const vector<int> &composite_fp_bits ,
		     const vector<string> &bit_names ) {

  ofstream ofs( loob_settings.ascii_fp_file_.c_str() );
  if( !ofs ){
    throw DACLIB::FileWriteOpenError( loob_settings.ascii_fp_file_.c_str() );
  }

  int num_bits_to_set = !loob_settings.dont_compress_fps_ ?
    count_if( fp_counts.begin() , fp_counts.end() ,
	      boost::bind( greater<int>() , _1 ,
			   loob_settings.compressed_counts_cutoff_ - 1 ) ) :
    fp_counts.size();

  if( loob_settings.bit_labels_to_ascii_file_ ) {
    string sep = loob_settings.bit_separator_;
    if( sep.empty() ) {
      sep = " ";
    }

    vector<int>::const_iterator p , ps;
    int i = 0;
    ofs << "MOLECULE";
    for( p = composite_fp_bits.begin() , ps = composite_fp_bits.end() ;
	 p != ps ; ++p , ++i ) {
      if( *p || loob_settings.dont_compress_fps_ ) {
	ofs << sep << bit_names[i];
      }
    }
    ofs << endl;
  }

  vector<pair<string,vector<int> > >::const_iterator p , ps;
  vector<int> fp_bits( num_bits_to_set , 0 );
  for( p = bits_set.begin() , ps = bits_set.end() ; p != ps ; ++p ) {
    ofs << p->first;
    if( loob_settings.bit_separator_.empty() ) {
      ofs << " ";
    }
    // copy p->second as we can't want to change the original
    vector<int> these_bits( p->second );
    // if we're doing compressed fps, cut these_bits down
    if( !loob_settings.dont_compress_fps_ ) {
      prepare_compressed_fp( loob_settings.compressed_counts_cutoff_ , fp_counts ,
			     full_to_compressed , these_bits );
    }
    // convert these_bits, with the bit numbers, to the fp representation
    // (0s and 1s)
    bits_to_fp( these_bits , fp_bits );
    vector<int>::const_iterator q , qs;
    for( q = fp_bits.begin() , qs = fp_bits.end() ; q != qs ; ++q ) {
      ofs << loob_settings.bit_separator_.c_str() << *q;
    }
    ofs << endl;

  }

  ofs << "COMPOSITE_FINGERPRINT";
  if( loob_settings.bit_separator_.empty() ) {
    ofs << " ";
  }
  vector<int>::const_iterator q , qs;
  for( q = composite_fp_bits.begin() , qs = composite_fp_bits.end() ; q != qs ; ++q ) {
    ofs << loob_settings.bit_separator_.c_str() << *q;
  }
  ofs << endl;

}

// *************************************************************************
// write out compact fingerprints - just the bit numbers of the set bits.
void write_compact_fp( LoobSettings &loob_settings ,
		       const vector<pair<string,vector<int> > > &bits_set ,
		       const vector<int> &full_to_compressed ,
		       const vector<int> &composite_fp_bits ) {

  cout << "Writing compact fps to " << loob_settings.compact_fp_file_ << endl;

  ofstream ofs( loob_settings.compact_fp_file_.c_str() );
  if( !ofs ){
    throw DACLIB::FileWriteOpenError( loob_settings.compact_fp_file_.c_str() );
  }

  vector<pair<string,vector<int> > >::const_iterator p , ps;
  vector<int>::const_iterator q , qs;
  for( p = bits_set.begin() , ps = bits_set.end() ; p != ps ; ++p ) {
    ofs << p->first;
    for( q = p->second.begin() , qs = p->second.end() ; q != qs ; ++q ) {
      if( !loob_settings.dont_compress_fps_ )
	ofs << "," << full_to_compressed[*q];
      else
	ofs << "," << *q;
    }
    ofs << endl;
  }

  // if this is a compressed fingerprint, all bits in the composite fingerprint
  // will be set, which makes this very easy
  ofs << "COMPOSITE_FINGERPRINT";
  for( int i = 0 , is = composite_fp_bits.size() ; i < is ; ++i ) {
    if( composite_fp_bits[i] ) {
      ofs << "," << i;
    }
  }
  ofs << endl;

}

// *************************************************************************
void write_labels_not_bits_file( LoobSettings &loob_settings ,
				 const vector<pair<string,vector<int> > > &bits_set ,
				 const vector<string> &bit_names ) {

  ofstream ofs( loob_settings.names_fp_file_.c_str() );
  if( !ofs ) {
    throw DACLIB::FileWriteOpenError( loob_settings.names_fp_file_.c_str() );
  }

  vector<int> bits_used( bit_names.size() , 0 );

  for( int i = 0 , is = bits_set.size() ; i < is ; ++i ) {
    ofs << bits_set[i].first;
    for( int j = 0 , js = bits_set[i].second.size() ; j < js ; ++j ) {
      ofs << "," << bit_names[bits_set[i].second[j]];
      bits_used[bits_set[i].second[j]] = 1;
    }
    ofs << endl;
  }

  ofs << "COMPOSITE_FINGERPRINT";
  for( int i = 0 , is = bits_used.size() ; i < is ; ++i ) {
    if( bits_used[i] ) {
      ofs << "," << bit_names[i];
    }
  }
  ofs << endl;

}

// *************************************************************************
// assume fp_counts is set to the right size and zeroed
void create_fp_counts( vector<pair<string,vector<int> > > &bits_set ,
		       vector<int> &fp_counts ) {

  vector<pair<string,vector<int> > >::iterator p , ps;
  vector<int>::iterator q , qs;

  for( p = bits_set.begin() , ps = bits_set.end() ; p != ps ; ++p ) {
    for( q = p->second.begin() , qs = p->second.end() ; q != qs ; ++q ) {
      fp_counts[*q]++;
    }
  }

}

// ****************************************************************************
string decode_2point_bit_code( const string &short_name ,
			       const map<string,vector<string> > &points_defs ,
			       const vector<float> &dist_bounds ) {

  static vector<string> points_names;
  if( points_names.empty() ) {
    map<string,vector<string> >::const_iterator q;
    for( q = points_defs.begin() ; q != points_defs.end() ; q++ ) {
      points_names.push_back( q->first );
    }
  }

  string label = short_name;
  ostringstream os;

  os << points_names[label[0]-'a'] << "_" << points_names[label[1]-'a'] << ":";
  int d = label[2] - 'a';
  if( 0 == d ) {
    os << "0";
  } else{
    os << dist_bounds[d-1];
  }
  os << "-" << dist_bounds[d];
  
  return os.str();

}

// ****************************************************************************
string decode_3point_bit_code( const string &short_name ,
			       const map<string,vector<string> > &points_defs ,
			       const vector<float> &dist_bounds ) {

  static vector<string> points_names;
  if( points_names.empty() ) {
    map<string,vector<string> >::const_iterator q;
    for( q = points_defs.begin() ; q != points_defs.end() ; q++ ) {
      points_names.push_back( q->first );
    }
  }

  // The 3 point codes are 3 distances and 3 point defs in a particular order.
  // They will be put out as f1->f2 dist, f1->f3 dist, f2->f3 dist.  It will
  // then be relatively straightforward to convert it to a Plurality query,
  // for instance.

  string label = short_name;
  ostringstream os;

  os << points_names[label[3]-'a'] << "_" << points_names[label[4]-'a'] << ":";
  int d = label[0] - 'a';
  if( 0 == d )
    os << "0";
  else
    os << dist_bounds[d-1];
  os << "-" << dist_bounds[d] << "::";
  
  os << points_names[label[3]-'a'] << "_" << points_names[label[5]-'a'] << ":";
  d = label[1] - 'a';
  if( 0 == d )
    os << "0";
  else
    os << dist_bounds[d-1];
  os << "-" << dist_bounds[d] << "::";
  
  os << points_names[label[4]-'a'] << "_" << points_names[label[5]-'a'] << ":";
  d = label[2] - 'a';
  if( 0 == d )
    os << "0";
  else
    os << dist_bounds[d-1];
  os << "-" << dist_bounds[d];
  
  return( os.str() );

}

// ****************************************************************************
string decode_4point_bit_code( const string &short_name ,
			       const map<string,vector<string> > &points_defs ,
			       const vector<float> &dist_bounds ) {

  // the 4 point codes are 6 distances and 4 features, again in a pre-determined
  // order: 1->2, 1->3, 2->3, 1->4, 2->4, 3->4
  static vector<string> points_names;
  if( points_names.empty() ) {
    map<string,vector<string> >::const_iterator q;
    for( q = points_defs.begin() ; q != points_defs.end() ; q++ ) {
      points_names.push_back( q->first );
    }
  }

  string label = short_name;
  ostringstream os;

  os << points_names[label[6]-'a'] << "_" << points_names[label[7]-'a'] << ":";
  int d = label[0] - 'a';
  if( 0 == d )
    os << "0";
  else
    os << dist_bounds[d-1];
  os << "-" << dist_bounds[d] << "::";

  os << points_names[label[6]-'a'] << "_" << points_names[label[8]-'a'] << ":";
  d = label[1] - 'a';
  if( 0 == d )
    os << "0";
  else
    os << dist_bounds[d-1];
  os << "-" << dist_bounds[d] << "::";

  os << points_names[label[7]-'a'] << "_" << points_names[label[8]-'a'] << ":";
  d = label[2] - 'a';
  if( 0 == d )
    os << "0";
  else
    os << dist_bounds[d-1];
  os << "-" << dist_bounds[d] << "::";

  os << points_names[label[6]-'a'] << "_" << points_names[label[9]-'a'] << ":";
  d = label[3] - 'a';
  if( 0 == d )
    os << "0";
  else
    os << dist_bounds[d-1];
  os << "-" << dist_bounds[d] << "::";

  os << points_names[label[7]-'a'] << "_" << points_names[label[9]-'a'] << ":";
  d = label[4] - 'a';
  if( 0 == d )
    os << "0";
  else
    os << dist_bounds[d-1];
  os << "-" << dist_bounds[d] << "::";

  os << points_names[label[8]-'a'] << "_" << points_names[label[9]-'a'] << ":";
  d = label[5] - 'a';
  if( 0 == d )
    os << "0";
  else
    os << dist_bounds[d-1];
  os << "-" << dist_bounds[d];
  
  return os.str();
  
}

// *************************************************************************
string decode_bit_name( const string &short_name ,
			const map<string,vector<string> > &points_defs ,
			const vector<float> &dist_bounds ) {

  if( 3 == short_name.length() ) {
    return decode_2point_bit_code( short_name , points_defs , dist_bounds );
  } else if( 6 == short_name.length() ) {
    return decode_3point_bit_code( short_name , points_defs , dist_bounds );
  } else if( 10 == short_name.length() ) {
    return decode_4point_bit_code( short_name , points_defs , dist_bounds );
  }

  return short_name;

}

// *************************************************************************
void report_fingerprint_counts( ostream &os , LoobSettings &loob_settings ,
				int mols_done ,
				const map<string,vector<string> > &points_defs ,
				const vector<string> &bit_names ,
				const vector<int> &fp_counts ) {

  vector<int>::const_iterator p , ps;
  vector<string>::const_iterator q;
  vector<pair<string,int> > output_set;
  for( p = fp_counts.begin() , ps = fp_counts.end() , q = bit_names.begin() ;
       p != ps ; ++p , ++q ) {
    if( *p >= loob_settings.compressed_counts_cutoff_ ) {
      output_set.push_back( make_pair( decode_bit_name( *q , points_defs , loob_settings.dist_bounds_ ) , *p ) );
    }
  }

  stable_sort( output_set.begin() , output_set.end() ,
	       boost::bind( greater<int>() ,
			    boost::bind( &pair<string,int>::second , _1 ) ,
			    boost::bind( &pair<string,int>::second , _2 ) ) );
  
  vector<pair<string,int> >::iterator r , rs;
  for( r = output_set.begin() , rs = output_set.end() ; r != rs ; ++r ){
    os << r->first
       << " in " << r->second << " of " << mols_done << " molecules." << endl;
  }

}

// *************************************************************************
void write_fingerprint_bit_names( ostream &os ,
				  LoobSettings &loob_settings ,
				  const vector<int> &fp_counts ,
				  const map<string,vector<string> > &points_defs ,
				  const vector<string> &bit_names ) {

  vector<string>::const_iterator p , ps;
  vector<int>::const_iterator q = fp_counts.begin();
  for( p = bit_names.begin() , ps = bit_names.end() ; p != ps ; ++p , ++q) {
    if( loob_settings.dont_compress_fps_ ||
	( !loob_settings.dont_compress_fps_ &&
	  *q >= loob_settings.compressed_counts_cutoff_ ) ) {
      os << decode_bit_name( *p , points_defs , loob_settings.dist_bounds_ )
	 << endl;
    }
  }

}
				  
// *************************************************************************
void fp_map_to_bit_names( const map<string,int> &fp_map ,
			  vector<string> &bit_names ) {

  bit_names = vector<string>( fp_map.size() );
  map<string,int>::const_iterator p , ps;
  for( p = fp_map.begin() , ps = fp_map.end() ; p != ps ; ++p ) {
    bit_names[p->second] = p->first;
  }

}

// *************************************************************************
void write_bit_names_file( LoobSettings &loob_settings ,
			   const vector<string> &bit_names ,
			   const vector<int> &fp_counts ,
			   const map<string,vector<string> > &points_defs ) {

  ofstream ofs( loob_settings.bit_names_file_.c_str() );
  if( !ofs.good() ) {
    throw DACLIB::FileWriteOpenError( loob_settings.bit_names_file_.c_str() );
  }

  write_fingerprint_bit_names( ofs , loob_settings , fp_counts ,
			       points_defs , bit_names );

}

// *************************************************************************
void write_log_file( LoobSettings &loob_settings , int mols_done ,
		     const map<string,vector<string> > &points_defs ,
		     const vector<string> &bit_names ,
		     const vector<int> &fp_counts ) {

  if( loob_settings.log_file_.empty() ) {
    return;
  }

  ofstream logstream( loob_settings.log_file_.c_str() );
  if( !logstream.good() ) {
    throw( DACLIB::FileWriteOpenError( loob_settings.log_file_.c_str() ) );
  }

  report_fingerprint_counts( logstream , loob_settings , mols_done ,
			     points_defs , bit_names , fp_counts );
  write_fingerprint_bit_names( logstream , loob_settings , fp_counts ,
			       points_defs , bit_names );
  
}

// *************************************************************************
void output_fingerprints( LoobSettings &loob_settings , int mols_done ,
			  vector<pair<string,vector<int> > > &bits_set ,
			  const map<string,int> &fp_map ,
			  const map<string,vector<string> > &points_defs ) {

  vector<int> fp_counts( fp_map.size() , 0 );
  create_fp_counts( bits_set , fp_counts );

  // if we're doing compressed fps, make a conversion from the full bit number
  // to the bit number in the compressed fingerprint
  vector<int> full_to_compressed;
  if( !loob_settings.dont_compress_fps_ ) {
    map_full_to_compressed( fp_counts , loob_settings.compressed_counts_cutoff_ ,
			    full_to_compressed );
  }
  vector<int> composite_fp_bits;
  make_composite_fp( loob_settings , fp_counts , composite_fp_bits );

  vector<string> bit_names;
  fp_map_to_bit_names( fp_map , bit_names );

  if( loob_settings.ascii_fps_ ) {
    try {
      write_ascii_fp( loob_settings , bits_set , fp_counts ,
		      full_to_compressed , composite_fp_bits , bit_names );
    } catch( DACLIB::FileWriteOpenError &e ) {
      cerr << "Error opening file " << e.what() << " for writing ASCII fingerprints." << endl;
    }
  }
  if( loob_settings.compact_fps_ ) {
    try {
      write_compact_fp( loob_settings , bits_set ,  full_to_compressed , composite_fp_bits );
    } catch( DACLIB::FileWriteOpenError &e ) {
      cerr << "Error opening file " << e.what() << " for writing compact fingerprints." << endl;
    }
  }
  if( loob_settings.labels_not_bits_ ) {
    try {
      write_labels_not_bits_file( loob_settings , bits_set , bit_names );
    } catch( DACLIB::FileWriteOpenError &e ) {
      cerr << "Error opening file " << e.what() << " for writing labels file." << endl;
    }

  }
  if( !loob_settings.bit_names_file_.empty() ) {
    try {
      write_bit_names_file( loob_settings , bit_names , fp_counts , points_defs );
    } catch( DACLIB::FileWriteOpenError &e ) {
      cerr << "Error opening file " << e.what() << " for writing bit names file." << endl;
    }

  }

  report_fingerprint_counts( cout , loob_settings , mols_done ,
			     points_defs , bit_names , fp_counts );

  write_log_file( loob_settings , mols_done , points_defs , bit_names ,
		  fp_counts );

}
