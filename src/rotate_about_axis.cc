//
// file rotate_about_axis.cc
// David Cosgrove
// AstraZeneca
// 27th November 2003
//
// This file contains the function rotate_about_axis which rotates a vector
// vec about an arbitrary axis by an angle expressed in radians.
// cribbed from 'Mathematics for 3D Game Programming and Computer Graphics' by
// Eric Lengyel

#include <cmath>

#include "stddefs.H"

namespace DACLIB {

// *****************************************************************************
void rotate_about_axis( const float axis[3] , float angle , float vec[3] ) {

  float ans[3] , norm_ax[3];
  float c = cos( angle );
  float s = sin( angle );
  float omc = 1.0 - c;

  norm_ax[0] = axis[0]; norm_ax[1] = axis[1]; norm_ax[2] = axis[2];
  float norm_len = DACLIB::length( norm_ax );
  norm_ax[0] /= norm_len; norm_ax[1] /= norm_len; norm_ax[2] /= norm_len;

  ans[0] = ( c + norm_ax[0] * norm_ax[0] * omc ) * vec[0] +
      ( norm_ax[0] * norm_ax[1] * omc - norm_ax[2] * s ) * vec[1] +
      ( norm_ax[0] * norm_ax[2] * omc + norm_ax[1] * s ) * vec[2];

  ans[1] = ( norm_ax[0] * norm_ax[1] * omc + norm_ax[2] * s ) * vec[0] +
      ( c + norm_ax[1] * norm_ax[1] * omc ) * vec[1] +
      ( norm_ax[1] * norm_ax[2] * omc - norm_ax[0] * s ) * vec[2];

  ans[2] = ( norm_ax[0] * norm_ax[2] * omc - norm_ax[1] * s ) * vec[0] +
      ( norm_ax[1] * norm_ax[2] * omc + norm_ax[0] * s ) * vec[1] +
      ( c + norm_ax[2] * norm_ax[2] * omc ) * vec[2];

  vec[0] = ans[0]; vec[1] = ans[1]; vec[2] = ans[2];

}

// *****************************************************************************
void rotate_about_axis( const double axis[3] , double angle , double vec[3] ) {

  double ans[3] , norm_ax[3];
  double c = cos( angle );
  double s = sin( angle );
  double omc = 1.0 - c;

  norm_ax[0] = axis[0]; norm_ax[1] = axis[1]; norm_ax[2] = axis[2];
  double norm_len = DACLIB::length( norm_ax );
  norm_ax[0] /= norm_len; norm_ax[1] /= norm_len; norm_ax[2] /= norm_len;

  ans[0] = ( c + norm_ax[0] * norm_ax[0] * omc ) * vec[0] +
      ( norm_ax[0] * norm_ax[1] * omc - norm_ax[2] * s ) * vec[1] +
      ( norm_ax[0] * norm_ax[2] * omc + norm_ax[1] * s ) * vec[2];

  ans[1] = ( norm_ax[0] * norm_ax[1] * omc + norm_ax[2] * s ) * vec[0] +
      ( c + norm_ax[1] * norm_ax[1] * omc ) * vec[1] +
      ( norm_ax[1] * norm_ax[2] * omc - norm_ax[0] * s ) * vec[2];

  ans[2] = ( norm_ax[0] * norm_ax[2] * omc - norm_ax[1] * s ) * vec[0] +
      ( norm_ax[1] * norm_ax[2] * omc + norm_ax[0] * s ) * vec[1] +
      ( c + norm_ax[2] * norm_ax[2] * omc ) * vec[2];

  vec[0] = ans[0]; vec[1] = ans[1]; vec[2] = ans[2];

}

} // end of namespace

