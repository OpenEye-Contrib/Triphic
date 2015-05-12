//
// file OverlayTrans.cc
// David Cosgrove
// AstraZeneca UK
// 12th February 2003
//
// This is the implementation file for the class
// OverlayTrans.  It holds two translations and a rotation, and is used for
// overlaying GraphicalObjects, most likely molecules.  The first translation
// is applied to the object, then the rotation, then the second translation.

#include <iterator>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "stddefs.H"
#include "DiamondOverlay.H"
#include "OverlayTrans.H"

// in eponymous file
namespace DACLIB {

  void centre_coords( float *coords , int num_points , float avge[3] );

}

using namespace std;

// **************************************************************************
// the transformation that moves moving onto fixed
OverlayTrans::OverlayTrans( float *moving , float *fixed , int natoms ) {

  // bring both sets of coords to origin
  DACLIB::centre_coords( moving , natoms , trans1_ );
  DACLIB::centre_coords( fixed , natoms , trans2_ );

  try {
    DACLIB::DiamondOverlay dov( fixed , moving , natoms );
    dov.get_rot_matrix( rot_ );
  } catch( DACLIB::DiamondOverlayError &e ) {
    throw e;
  }

}

// **************************************************************************
// the transformation that moves cds1 onto cds2 with weights
OverlayTrans::OverlayTrans( float *moving , float *fixed , float *weights ,
			    int natoms ) {

  // bring both sets of coords to origin
  DACLIB::centre_coords( moving , natoms , trans1_ );
  DACLIB::centre_coords( fixed , natoms , trans2_ );

  DACLIB::DiamondOverlay dov( fixed , moving , weights , natoms );
  dov.get_rot_matrix( rot_ );

}

// **************************************************************************
void OverlayTrans::overlay( int num_points , float *cds ) const {

  int      i;

  float    *t;

  for( i = num_points , t = cds ; i ; i-- ) {
    *t -= trans1_[0]; t++;
    *t -= trans1_[1]; t++;
    *t -= trans1_[2]; t++;
  }

  for( i = num_points , t = cds ; i ; i-- ) {
    float cds[3];
    cds[0] = rot_[0][0] * t[0] + rot_[0][1] * t[1] + rot_[0][2] * t[2];
    cds[1] = rot_[1][0] * t[0] + rot_[1][1] * t[1] + rot_[1][2] * t[2];
    cds[2] = rot_[2][0] * t[0] + rot_[2][1] * t[1] + rot_[2][2] * t[2];
    *t = cds[0] + trans2_[0]; t++;
    *t = cds[1] + trans2_[1]; t++;
    *t = cds[2] + trans2_[2]; t++;
  }

}

// **************************************************************************
ostream &operator<<( ostream &s , const OverlayTrans &ot ) {

  s << setw( 10 ) << setprecision( 10 )
    << ot.trans1_[0] << " " << ot.trans1_[1] << " " << ot.trans1_[2] << " "
    << ot.trans2_[0] << " " << ot.trans2_[1] << " " << ot.trans2_[2] << " "
    << ot.rot_[0][0] << " " << ot.rot_[0][1] << " " << ot.rot_[0][2] << " "
    << ot.rot_[1][0] << " " << ot.rot_[1][1] << " " << ot.rot_[1][2] << " "
    << ot.rot_[2][0] << " " << ot.rot_[2][1] << " " << ot.rot_[2][2];

  return s;

}

// **************************************************************************
istream &operator>>( istream &s , OverlayTrans &ot ) {

  s >> ot.trans1_[0] >> ot.trans1_[1] >> ot.trans1_[2]
    >> ot.trans2_[0] >> ot.trans2_[1] >> ot.trans2_[2]
    >> ot.rot_[0][0] >> ot.rot_[0][1] >> ot.rot_[0][2]
    >> ot.rot_[1][0] >> ot.rot_[1][1] >> ot.rot_[1][2]
    >> ot.rot_[2][0] >> ot.rot_[2][1] >> ot.rot_[2][2];
  
  return s;

}

// **************************************************************************
// write to binary stream
void OverlayTrans::binary_write( ostream &os ) const {

  os.write( (char *) trans1_ , 3 * sizeof( float ) );
  os.write( (char *) trans2_ , 3 * sizeof( float ) );
  os.write( (char *) rot_[0] , 3 * sizeof( float ) );
  os.write( (char *) rot_[1] , 3 * sizeof( float ) );
  os.write( (char *) rot_[2] , 3 * sizeof( float ) );

}

// **************************************************************************
// read from binary stream
void OverlayTrans::binary_read( istream &os ) {

  os.read( (char *) trans1_ , 3 * sizeof( float ) );
  os.read( (char *) trans2_ , 3 * sizeof( float ) );
  os.read( (char *) rot_[0] , 3 * sizeof( float ) );
  os.read( (char *) rot_[1] , 3 * sizeof( float ) );
  os.read( (char *) rot_[2] , 3 * sizeof( float ) );

}

