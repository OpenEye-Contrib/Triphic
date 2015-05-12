
#include "stddefs.H"

namespace DACLIB {

// **************************************************************************
// calculate the torsion angle between the 4 sets of coords given.
// Returns angle between -180 and 180, with clockwise positive when looking
// from 1 to 2 (counting from 0, of course)
// This code is cribbed unashamedly from some old Viking stuff used in BasilI.
float torsion( const float *cds1 , const float *cds2 , const float *cds3 ,
               const float *cds4 , bool degs = true ) {

  float    ab[3] , bc[3] , cd[3] , abc[3] , abc_len , zax[3] , zax_len ,
      bcd[3] , cos_tors , sin_tors , tors;

  DACLIB::join_vector( cds1 , cds2 , cd );
  DACLIB::join_vector( cds2 , cds3 , bc );
  DACLIB::join_vector( cds3 , cds4 , ab );

  // take the cross of ab and bc
  DACLIB::cross_product( ab , bc , abc );
  abc_len = DACLIB::length( abc );
  if( abc_len >= 1.0e-6 ) {
    // take cross of BC and abc, call it the z axis
    DACLIB::cross_product( abc , bc , zax );
    zax_len = DACLIB::length( zax );
    if( zax_len >= 1.0e-6 ) {
      // take cross between bc and cd to give line perp to plane BCD.
      // Angle between this line and line perp to plane ABC is torsion.
      DACLIB::cross_product( bc , cd , bcd );
      cos_tors = DACLIB::dot_product( abc , bcd ) / abc_len;
      sin_tors = DACLIB::dot_product( zax , bcd ) / zax_len;
      tors = atan2( sin_tors , cos_tors );
    } else
      tors = 0.0; // B--C--D are co-linear
  } else
    tors = 0.0; // A--B--C are co-linear

  if( degs )
    return( tors * 180.0 / M_PI );
  else
    return tors;

}

double torsion( const double *cds1 , const double *cds2 , const double *cds3 ,
                const double *cds4 , bool degs = true ) {

  double    ab[3] , bc[3] , cd[3] , abc[3] , abc_len , zax[3] , zax_len ,
      bcd[3] , cos_tors , sin_tors , tors;

  DACLIB::join_vector( cds1 , cds2 , cd );
  DACLIB::join_vector( cds2 , cds3 , bc );
  DACLIB::join_vector( cds3 , cds4 , ab );

  // take the cross of ab and bc
  DACLIB::cross_product( ab , bc , abc );
  abc_len = DACLIB::length( abc );
  if( abc_len >= 1.0e-6 ) {
    // take cross of BC and abc, call it the z axis
    DACLIB::cross_product( abc , bc , zax );
    zax_len = DACLIB::length( zax );
    if( zax_len >= 1.0e-6 ) {
      // take cross between bc and cd to give line perp to plane BCD.
      // Angle between this line and line perp to plane ABC is torsion.
      DACLIB::cross_product( bc , cd , bcd );
      cos_tors = DACLIB::dot_product( abc , bcd ) / abc_len;
      sin_tors = DACLIB::dot_product( zax , bcd ) / zax_len;
      tors = atan2( sin_tors , cos_tors );
    } else
      tors = 0.0; // B--C--D are co-linear
  } else
    tors = 0.0; // A--B--C are co-linear

  if( degs )
    return( tors * 180.0 / M_PI );
  else
    return tors;

}

}
