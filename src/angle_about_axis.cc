//
// file angle_about_axis.cc
// David Cosgrove
// AstraZeneca
// 27th November 2003
//
// This file contains the function angle_about_axis which measures the angle
// between two vectors projected onto a plane perpendicular to an arbitrary
// axis.  Used, for instance, to find out what angle a O-H bond needs to be
// rotated about a C-O bond so that the H atom is on top of another H atom.

#include <cmath>

#include "stddefs.H"

namespace DACLIB {

// ****************************************************************************
double angle_about_axis( const float axis[3] , const float vec1[3] ,
const float vec2[3] ) {

  // find the component of vec1 and vec2 perpendicular to a normalised axis
  float axis_len = DACLIB::length( axis );
  float axis_norm[3];
  axis_norm[0] = axis[0] / axis_len;
  axis_norm[1] = axis[1] / axis_len;
  axis_norm[2] = axis[2] / axis_len;

  float vec1_perp[3] , vec2_perp[3] , vec1_proj , vec2_proj;
  vec1_proj = DACLIB::dot_product( vec1 , axis_norm );
  vec1_perp[0] = vec1[0] - vec1_proj * axis_norm[0];
  vec1_perp[1] = vec1[1] - vec1_proj * axis_norm[1];
  vec1_perp[2] = vec1[2] - vec1_proj * axis_norm[2];
  
  vec2_proj = DACLIB::dot_product( vec2 , axis_norm );
  vec2_perp[0] = vec2[0] - vec2_proj * axis_norm[0];
  vec2_perp[1] = vec2[1] - vec2_proj * axis_norm[1];
  vec2_perp[2] = vec2[2] - vec2_proj * axis_norm[2];

  // now find and return the angle between the two perpendicular vectors
  return DACLIB::angle( vec1_perp , DACLIB::length( vec1_perp ) ,
                        vec2_perp , DACLIB::length( vec2_perp ) );

}


// ****************************************************************************
double angle_about_axis( const double axis[3] , const double vec1[3] ,
const double vec2[3] ) {

  // find the component of vec1 and vec2 perpendicular to a normalised axis
  double axis_len = DACLIB::length( axis );
  double axis_norm[3];
  axis_norm[0] = axis[0] / axis_len;
  axis_norm[1] = axis[1] / axis_len;
  axis_norm[2] = axis[2] / axis_len;

  double vec1_perp[3] , vec2_perp[3] , vec1_proj , vec2_proj;
  vec1_proj = DACLIB::dot_product( vec1 , axis_norm );
  vec1_perp[0] = vec1[0] - vec1_proj * axis_norm[0];
  vec1_perp[1] = vec1[1] - vec1_proj * axis_norm[1];
  vec1_perp[2] = vec1[2] - vec1_proj * axis_norm[2];
  
  vec2_proj = DACLIB::dot_product( vec2 , axis_norm );
  vec2_perp[0] = vec2[0] - vec2_proj * axis_norm[0];
  vec2_perp[1] = vec2[1] - vec2_proj * axis_norm[1];
  vec2_perp[2] = vec2[2] - vec2_proj * axis_norm[2];

  // now find and return the angle between the two perpendicular vectors
  return DACLIB::angle( vec1_perp , DACLIB::length( vec1_perp ) ,
                        vec2_perp , DACLIB::length( vec2_perp ) );

}

} // end of namespace
