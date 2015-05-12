/*
  file centre_coords.c
  David Cosgrove
  AstraZeneca
  19th October 1998

  This function takes a set of double coords, packed into a 1D array,
  assumed to be 3D coords (ie 3 doubles per coord) and moves them
  so that their avge coordinate is at the origin

  */

namespace DACLIB {

  // ***************************************************************************
  void centre_coords( double *coords , int num_points , double avge[3] ) {

    int      i;

    avge[0] = avge[1] = avge[2] = 0.0;

    double *cds = coords;
    for( i = 0 ; i < num_points ; i++ ) {
      avge[0] += *cds; cds++;
      avge[1] += *cds; cds++;
      avge[2] += *cds; cds++;
    }

    avge[0] /= (double) num_points;
    avge[1] /= (double) num_points;
    avge[2] /= (double) num_points;

    cds = coords;
    for( i = 0 ; i < num_points ; i++ ) {
      *cds -= avge[0]; cds++;
      *cds -= avge[1]; cds++;
      *cds -= avge[2]; cds++;
    }

  }

  // ***************************************************************************
  void centre_coords( float *coords , int num_points , float avge[3] ) {

    int      i;

    avge[0] = avge[1] = avge[2] = 0.0;

    float *cds = coords;
    for( i = 0 ; i < num_points ; i++ ) {
      avge[0] += *cds; cds++;
      avge[1] += *cds; cds++;
      avge[2] += *cds; cds++;
    }

    avge[0] /= (float) num_points;
    avge[1] /= (float) num_points;
    avge[2] /= (float) num_points;

    cds = coords;
    for( i = 0 ; i < num_points ; i++ ) {
      *cds -= avge[0]; cds++;
      *cds -= avge[1]; cds++;
      *cds -= avge[2]; cds++;
    }

  }

  // ***************************************************************************
  void centre_coords( double *coords , int num_points , float avge[3] ) {

    int      i;

    avge[0] = avge[1] = avge[2] = 0.0;

    double *cds = coords;
    for( i = 0 ; i < num_points ; i++ ) {
      avge[0] += *cds; cds++;
      avge[1] += *cds; cds++;
      avge[2] += *cds; cds++;
    }

    avge[0] /= (float) num_points;
    avge[1] /= (float) num_points;
    avge[2] /= (float) num_points;

    cds = coords;
    for( i = 0 ; i < num_points ; i++ ) {
      *cds -= avge[0]; cds++;
      *cds -= avge[1]; cds++;
      *cds -= avge[2]; cds++;
    }

  }

} // end of namespace
