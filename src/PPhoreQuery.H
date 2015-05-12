//
// file PPhoreQuery.H
// Dave Cosgrove
// 13th September 2007
//
// This class holds the details of a pharmacophore query to be used by plurality.

#ifndef DAC_PPHORE_QUERY
#define DAC_PPHORE_QUERY

#include "VolumeGrid.H"

#include <iosfwd>
#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>

// **************************************************************************

class PPhoreQueryFileReadError {

public :

  explicit PPhoreQueryFileReadError( const std::string &filename ,
				     const std::string &file_line ) {
    msg_ = std::string( "Error reading line \n" ) + file_line +
      std::string( "\nof file " ) + filename;
  }
  virtual ~PPhoreQueryFileReadError() {}

  virtual const char *what() {
    return msg_.c_str();
  }

private :

  std::string msg_;

};

// ******************************************************************************

class PPhoreQueryFileError {

public :

  explicit PPhoreQueryFileError( const std::string &filename ,
				 const std::string &file_problem ) {
    msg_ = std::string( "Error reading file : " ) + filename + std::string( " :\n" ) +
      file_problem;

  }
  virtual ~PPhoreQueryFileError() {}

  virtual const char *what() {
    return msg_.c_str();
  }

private :

  std::string msg_;

};

// ******************************************************************************

typedef boost::tuple<int,float,float,float> BTC;
typedef boost::tuple<int,int,float,float> BTD;
typedef boost::tuple<int,int,int,float,float> BTA;
typedef boost::tuple<int,int,int,int,float,float> BTT;
typedef boost::tuple<float,float,float,float> BTV; // 1st 3 are coords, last is rad

namespace OEChem {
  class OEMol;
}

class PluralityHit;
class SinglePPhoreSite;

class PPhoreQuery {

public :

  PPhoreQuery();
  ~PPhoreQuery();

  // throws a FileReadOpenError exception if there's a problem
  void read_query_file( const std::string &filename );

  void search_pphore_sites( std::vector<SinglePPhoreSite *> &target_sites ,
			    OEChem::OEMol &target_mol , int conf_num ,
			    boost::shared_ptr<DACLIB::VolumeGrid> &protein_grid ,
			    boost::shared_ptr<DACLIB::VolumeGrid> &soft_exc_vol_grid ,
			    const std::vector<std::pair<std::string,DACLIB::VolumeGrid *> > &score_vol_grids ,
			    bool first_hit_only ,
			    std::vector<boost::shared_ptr<PluralityHit> > &hits );

  const std::vector<std::string> &point_types() const { return point_types_; }
  const std::vector<BTC> &coords() const { return coords_; }
  const std::vector<BTD> &distances() const { return distances_; }
  const std::vector<BTA> &angles() const { return angles_; }
  const std::vector<BTT> &torsions() const { return torsions_; }
  const std::vector<BTV> &hard_exc_vols() const { return hard_exc_vols_; }
  const std::vector<BTV> &soft_exc_vols() const { return soft_exc_vols_; }

  void pack_into_pvm_buffer();
  void unpack_from_pvm_buffer();

private :

  std::string filename_; // where the query came from
  std::vector<std::string> point_types_;
  std::vector<BTC> coords_;
  std::vector<BTD> distances_;
  std::vector<BTA> angles_;
  std::vector<BTT> torsions_;
  std::vector<BTV> hard_exc_vols_; // excluded volumes - hard can't be hit
  std::vector<BTV> soft_exc_vols_; // soft measures overlay

  float **sq_dist_vals_; /* square matrix of distance constraints, upper triangle
			    is upper limit, lower triangle is lower limit */
  float *overlay_cds_; // extracted from coords_ or 0 if there aren't any

  void clear_data();

  void parse_point_line( std::istringstream &ss );
  void parse_coords_line( std::istringstream &ss );
  void parse_distance_line( std::istringstream &ss );
  void parse_angle_line( std::istringstream &ss );
  void parse_torsion_line( std::istringstream &ss );
  void parse_exc_vol_line( std::istringstream &ss , std::vector<BTV> &exc_vols );

  // make sure the query details are all consistent
  void verify_query_details();
  void verify_coord_details();
  void verify_distance_details();
  void verify_angle_details();
  void verify_torsion_details();

  void assemble_overlay_coords();
  void assemble_distance_constraints();

  //
  // *********************************************************
  // functions used for searching target_sites to see if there is a match
  // to this query
  float **make_squared_distance_matrix( std::vector<SinglePPhoreSite *> &target_sites );
  // initial correspondence matrix m0
  short **form_m0( std::vector<SinglePPhoreSite *> &target_sites );
  bool solution_possible( short **m , int num_query_points ,
			  int num_target_points );
  // use an Ullmann-type procedure to find the hits
  void extend_subgraph( int row_num , short **m , const int num_rows ,
			const int num_cols , const float **sq_dist_mat ,
			std::vector<std::vector<int> > &hit_points );
  // Ullmann's refinement process
  bool refine_m( const int num_cols , const float **sq_dist_mat ,
		 const int row_num , const int col_num , short **m );
  bool check_for_solution( short **m , const int num_rows , const int num_cols );
  void store_solution( short **m , const int num_rows , const int num_cols ,
		       std::vector<std::vector<int> > &hit_points );

  // check that the hit identified by ullman satisfies the angle contraints
  bool angles_ok( std::vector<SinglePPhoreSite *> &target_sites ,
		  const std::vector<int> &hit_points ,
		  std::vector<float> &hit_angles );
  bool angles_with_dirs_ok( BTA &angle ,
			    std::vector<SinglePPhoreSite *> &target_sites ,
			    const std::vector<int> &hit_points ,
			    float &hit_ang );
  // check that the hit identified by ullman satisfies the torsion contraints
  bool torsions_ok( std::vector<SinglePPhoreSite *> &target_sites ,
		    const std::vector<int> &hit_points ,
		    std::vector<float> &torsions );
  bool torsions_with_dirs_ok( BTT &torsion ,
			      std::vector<SinglePPhoreSite *> &target_sites ,
			      const std::vector<int> &hit_points ,
			      float &hit_tors );

};

#endif
