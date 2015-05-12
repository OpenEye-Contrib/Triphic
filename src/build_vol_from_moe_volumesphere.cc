//
// file build_vol_from_moe_volumesphere.cc
// David Cosgrove
// AstraZeneca
// 1st April 2014
//
// Reads a volumesphere record from a MOE ph4 file and creates a VolumeGrid
// object

#include "VolumeGrid.H"

#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using namespace boost;
using namespace std;

// ****************************************************************************
DACLIB::VolumeGrid *build_vol_from_moe_volumesphere( string &vol_rec ) {

  // first line of incoming string (up to \n) is a header line
  size_t i = 0;
  while( '\n' != vol_rec[i] ) {
    ++i;
  }
  vol_rec = vol_rec.substr( i + 1 );

  vector<BTV> spheres; // defined in VolumeGrid.H

  // records are in blocks of 4 - x, y, z and radius.
  istringstream iss( vol_rec );
  float x_max = -numeric_limits<float>::max();
  float y_max = -numeric_limits<float>::max();
  float z_max = -numeric_limits<float>::max();
  float x_min = numeric_limits<float>::max();
  float y_min = numeric_limits<float>::max();
  float z_min = numeric_limits<float>::max();

  while( 1 ) {
    float x , y , z , rad;
    iss >> x >> y >> z >> rad;
    if( iss.fail() ) {
      break;
    }
    spheres.push_back( BTV( x , y , z , rad ) );
    if( x - rad < x_min ) {
      x_min = x - rad;
    }
    if( y - rad < y_min ) {
      y_min = y - rad;
    }
    if( z - rad < z_min ) {
      z_min = z - rad;
    }
    if( x + rad > x_max ) {
      x_max = x + rad;
    }
    if( y + rad > y_max ) {
      y_max = y + rad;
    }
    if( z + rad > z_max ) {
      z_max = z + rad;
    }
  }

  float blf[3] = { x_min , y_min , z_min };
  float trb[3] = { x_max , y_max , z_max };

  DACLIB::VolumeGrid *retval = new DACLIB::VolumeGrid( blf , trb );
  retval->drop_spheres_in( spheres );

  return retval;

}
