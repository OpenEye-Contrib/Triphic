
#include "stddefs.H"

// **************************************************************************
// calculate the angle between the 3 sets of coords given.

namespace DACLIB {

float angle( const float *cds1 , const float *cds2 , const float *cds3 ,
	     bool degs = true ) {

  float vec1[3] , vec2[3];
  
  join_vector( cds2 , cds1 , vec1 );
  join_vector( cds2 , cds3 , vec2 );
  float ang = atan2( sin_angle( vec1 , DACLIB::length( vec1 ) ,
				vec2 , DACLIB::length( vec2 ) ) ,
		     cos_angle( vec1 , DACLIB::length( vec1 ) ,
				vec2 , DACLIB::length( vec2 ) ) );
  if( degs )
    ang = 180.0 * ang / M_PI;

  return ang;

}

double angle( const double *cds1 , const double *cds2 , const double *cds3 ,
	      bool degs = true ) {

  double vec1[3] , vec2[3];
  
  join_vector( cds2 , cds1 , vec1 );
  join_vector( cds2 , cds3 , vec2 );
  double ang = atan2( sin_angle( vec1 , DACLIB::length( vec1 ) ,
				 vec2 , DACLIB::length( vec2 ) ) ,
		      cos_angle( vec1 , DACLIB::length( vec1 ) ,
				 vec2 , DACLIB::length( vec2 ) ) );
  if( degs )
    ang = 180.0 * ang / M_PI;

  return ang;

}

}
