//
// file OverlayTrans.H
// Dave Cosgrove
// AstraZeneca UK
// 23rd March 2000
//
// This is the interface (and probably the implementation) file for the class
// OverlayTrans.  It holds two translations and a rotation, and is used for
// overlaying GraphicalObjects, most likely molecules.  The first translation
// is applied to the object, then the rotation, then the second translation.

#ifndef DAC_OVERLAY_TRANS
#define DAC_OVERLAY_TRANS

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

#include "stddefs.H"

// **************************************************************************

class OverlayTrans {

protected:

  float      trans1_[3] , trans2_[3];

  float      rot_[3][3];

public:

  OverlayTrans() {
    trans1_[0] = trans1_[1] = trans1_[2] = 0.0;
    trans2_[0] = trans2_[1] = trans2_[2] = 0.0;
    rot_[0][0] = rot_[0][1] = rot_[0][2] = 0.0;
    rot_[1][0] = rot_[1][1] = rot_[1][2] = 0.0;
    rot_[2][0] = rot_[2][1] = rot_[2][2] = 0.0;
  }

  // translate by t1, rotate and translate by t2
  OverlayTrans( const float t1[3] , const float t2[3] , const float r[3][3] ) {
    trans1_[0] = t1[0]; trans1_[1] = t1[1]; trans1_[2] = t1[2];
    trans2_[0] = t2[0]; trans2_[1] = t2[1]; trans2_[2] = t2[2];
    rot_[0][0] = r[0][0]; rot_[0][1] = r[0][1]; rot_[0][2] = r[0][2];
    rot_[1][0] = r[1][0]; rot_[1][1] = r[1][1]; rot_[1][2] = r[1][2];
    rot_[2][0] = r[2][0]; rot_[2][1] = r[2][1]; rot_[2][2] = r[2][2];
  }

  // the transformation that moves moving onto fixed
  OverlayTrans( float *moving , float *fixed , int natoms );
  // the transformation that moves cds1 onto cds2 with weights
  OverlayTrans( float *moving , float *fixed , float *weights , int natoms );

  ~OverlayTrans() {}

  void overlay( int num_points , float *cds ) const;

  friend std::ostream &operator<<( std::ostream &s , const OverlayTrans &ot );
  friend std::istream &operator>>( std::istream &s , OverlayTrans &ot );

  void get_trans1( float t1[3] ) const {
    t1[0] = trans1_[0]; t1[1] = trans1_[1]; t1[2] = trans1_[2]; 
  }
  void get_trans2( float t2[3] ) const {
    t2[0] = trans2_[0]; t2[1] = trans2_[1]; t2[2] = trans2_[2]; 
  }
  void get_rot( float r[3][3] ) const {
    r[0][0] = rot_[0][0]; r[0][1] = rot_[0][1]; r[0][2] = rot_[0][2];
    r[1][0] = rot_[1][0]; r[1][1] = rot_[1][1]; r[1][2] = rot_[1][2];
    r[2][0] = rot_[2][0]; r[2][1] = rot_[2][1]; r[2][2] = rot_[2][2];
  }

  // write to binary stream
  void binary_write( std::ostream &os ) const;
  // read from binary stream
  void binary_read( std::istream &os );

};

#endif
