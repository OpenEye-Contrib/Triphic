/*
  set of functions that will take a real, symmetric matrix and compute
  all eigen vectors and eigenvalues.  Returns the eigenvectors sorted in
  ascending order of eigenvalue.
  Specialisation for for 4x4 float matrices, returning only eigenvector for
  highest eigenvalue
  Assumes that ppr_a and ppr_v have been created using
  DACLIB::make_square_matrix such that ppr_a[0] points to 16 contiguous floats
  pointed to by ppr_a[0-3].

*/

#include <iostream>
#include <cstdlib>
#include <cmath>

#include "stddefs.H"

namespace DACLIB {

#define ROTATE( a,i,j,k,l ) r_g=a[i][j];\
   r_h=a[k][l];a[i][j]=r_g-r_s*(r_h+r_g*r_tau);\
   a[k][l]=r_h+r_s*(r_g-r_h*r_tau);

  /*  function vf_jacobi

  This function takes a real, symmetric matrix ppr_a[0-(n-1)][0-(n-1)] and
  computes all its eigenvalues and eigenvectors.  The above-diagonal elements
  of a are destroyed in the process.  The eigenvalues are returned in 
  array pr_d, and the eigenvectors as the columns of ppr_v[0-(n-1)][0-(n-1)].
  nrot returns the number of rotations used.
  It uses the method of Jacobi rotations, and is copied with minor 
  alterations from 'Numerical Recipes in C', 2nd edition, p. 467 
  
  */

  /******************************************************************************/

void vf_jacobi_4b4f( float **ppr_a , float *pr_d , float **ppr_v ,
                     int &nrot ) {

  int      i , ip , iq , j;

  float   r_c , r_g , r_h , r_s , r_sm , r_t , r_tau , r_theta , r_tresh;

  float pr_b[4] , pr_z[4] = { 0.0F , 0.0F , 0.0F , 0.0F };

  float *pa = ppr_a[0];
  float *pv = ppr_v[0];

  /* initialise ppr_v to the identity matrix */
  pv[1] = pv[2] = pv[3] = pv[4] = pv[6] = pv[7] = pv[8] = pv[9]
      = pv[11] = pv[12] = pv[13] = pv[14] = 0.0F;
  pv[0] = pv[5] = pv[10] = pv[15] = 1.0F;

  /* initialise pr_b and pr_d to the diagonal of ppr_a and zero ppr_z */
  pr_b[0] = pr_d[0] = pa[0];
  pr_b[1] = pr_d[1] = pa[5];
  pr_b[2] = pr_d[2] = pa[10];
  pr_b[3] = pr_d[3] = pa[15];

  nrot = 0;

  for( i = 1 ; i <= 50 ; i++ ) {

    /* sum the off-diagonal elements */
    r_sm = pa[1] + pa[2] + pa[3] + pa[6] + pa[7] + pa[11];

    /* if r_sm has reached 0.0, the problem has converged */
    if( 0.0 == r_sm ) {
      return;
    }

    /* set r_tresh higher for 1st three sweeps */
    if( i < 4 ) {
      r_tresh = 0.2 * r_sm / ( 16.0F );
    } else {
      r_tresh = 0.0;
    }

    for( ip = 0 ; ip < 3 ; ++ip ) {
      for( iq = ip + 1 ; iq < 4 ; ++iq ) {
        r_g = 100.0 * fabs( ppr_a[ip][iq] );
        /*  after 4 sweeps, skip the rotation if the off-diagonal element
        is small */
        if( i > 4 &&  fabs( pr_d[ip] ) + r_g == fabs( pr_d[ip] ) &&
            fabs( pr_d[iq] ) + r_g == fabs( pr_d[iq] ) ) {
          ppr_a[ip][iq] = 0.0;
        } else if( fabs( ppr_a[ip][iq] ) > r_tresh ) {
          r_h = pr_d[iq] - pr_d[ip];
          if( fabs( r_h ) + r_g == fabs( r_h ) ) {
            r_t = ( ppr_a[ip][iq] ) / r_h;
          } else {
            r_theta = 0.5 * r_h / ( ppr_a[ip][iq] );
            r_t = 1.0 /( fabs( r_theta ) + sqrt( 1.0 + r_theta * r_theta ) );
            if( r_theta < 0.0 )
              r_t = - r_t;
          }
          r_c = 1.0 / sqrt( 1.0 + r_t * r_t );
          r_s = r_t * r_c;
          r_tau = r_s / ( 1.0 + r_c );
          r_h = r_t * ppr_a[ip][iq];
          pr_z[ip] -= r_h;
          pr_z[iq] += r_h;
          pr_d[ip] -= r_h;
          pr_d[iq] += r_h;
          ppr_a[ip][iq] = 0.0;

          if( 1 == ip ) {
            ROTATE( ppr_a , 0 , 1 , 0 , iq );
          } else if( 2 == ip ) {
            ROTATE( ppr_a , 0 , 2 , 0 , iq );
            ROTATE( ppr_a , 1 , 2 , 1 , iq );
          }

          for( j = ip+1 ; j < iq ; j++ ) {
            ROTATE( ppr_a , ip , j , j , iq );
          }

          if( 1 == iq ) {
            ROTATE( ppr_a , ip , 2 , iq , 2 );
            ROTATE( ppr_a , ip , 3 , iq , 3 );
          } else if( 2 == iq ) {
            ROTATE( ppr_a , ip , 3 , iq , 3 );
          }

          ROTATE( ppr_v , 0 , ip , 0 , iq );
          ROTATE( ppr_v , 1 , ip , 1 , iq );
          ROTATE( ppr_v , 2 , ip , 2 , iq );
          ROTATE( ppr_v , 3 , ip , 3 , iq );
          nrot++;
        }
      }
    }
    pr_b[0] += pr_z[0];
    pr_b[1] += pr_z[1];
    pr_b[2] += pr_z[2];
    pr_b[3] += pr_z[3];
    pr_d[0] = pr_b[0];
    pr_d[1] = pr_b[1];
    pr_d[2] = pr_b[2];
    pr_d[3] = pr_b[3];
    pr_z[0] = pr_z[1] = pr_z[2] = pr_z[3] = 0.0F;

  }
  printf( "Error : too many iterations in routine vf_jacobi\n" );

}

  /****************************************************************************/
  /* function eigen_solve - the main driver for the eigen solver

  Takes an array a, of size n by n and returns the sorted eigenvectors
  in vectors.  The top half of a is destroyed.  a must be real and
  symmetric
  
  */
  void eigen_solve_4b4f( float **a , float &values , float *vectors ) {

    int     nrot;

    float   evalues[4];
    static float **evectors = 0;
    if( !evectors ) {
      DACLIB::make_square_matrix( evectors , 4 );
    }

    vf_jacobi_4b4f( a , evalues , evectors , nrot );

    // extract the vector with highest eigenvalue
    int largest_eigval = 0;
    if( evalues[1] > evalues[largest_eigval] ) {
      largest_eigval = 1;
    } 
    if( evalues[2] > evalues[largest_eigval] ) {
      largest_eigval = 2;
    }
    if( evalues[3] > evalues[largest_eigval] ) {
      largest_eigval = 3;
    }

    values = evalues[largest_eigval];
    vectors[0] = evectors[0][largest_eigval];
    vectors[1] = evectors[1][largest_eigval];
    vectors[2] = evectors[2][largest_eigval];
    vectors[3] = evectors[3][largest_eigval];
    
  }

} // end of namespace DACLIB
