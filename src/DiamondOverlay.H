//
// file DiamondOverlay.H
// David Cosgrove
// AstraZeneca
// 23rd January 2007
//
// This class uses Diamond's algorithm to find the rotation matrix that optimally
// superimposes one set of coordinates on another. It assumes that the coordinates
// are already centred on the origin.
// R. Diamond, Acta. Cryst., A44, p211-216 (1988). 'A Note on the Rotational
// Superposition Problem'. (D88)
// This is just a C++ class-based version of the C code we've been using for years
// that was a translation of the FORTRAN implementation JAG and I coded up in 1994.

#ifndef DAC_DIAMOND_OVERLAY__
#define DAC_DIAMOND_OVERLAY__

#include <stdexcept>
#include <vector>

namespace DACLIB {

  class DiamondOverlayError : public std::runtime_error {
  public :
    explicit DiamondOverlayError( const char *msg = "DiamondOverlay error" ) :
      runtime_error( msg ) {}
  };

  class DiamondOverlay {

  private :

    std::vector<float> fixed_ , moving_ , weights_;
    float **rot_matrix_;
    float **p_;
    float e_zero_ , e_ , vector_[4]; // error before and after transformation

    void calc_e_zero();
    void calc_m( float m[3][3] );
    void calc_q( const float m[3][3] , float q[3][3] );
    void calc_e();

    // make sure that there are at least 3 distinct sets of coordinates
    // for at least 1 of the input coord sets.  This includes cases
    // where 2 coords are on top of each other. Throws exception if not
    void check_input_coords() const;
    // checks given coords for problem above.
    bool at_least_3_coords( const std::vector<float> &cds ) const;

  public :

    // coords, packed with x,y,z in sequence for each point
    DiamondOverlay( const std::vector<float> &fixed ,
                    const std::vector<float> &moving );
    DiamondOverlay( const std::vector<float> &fixed ,
                    const std::vector<float> &moving ,
                    const std::vector<float> &weights );
    DiamondOverlay( const float *fixed , const float *moving , int num_points );
    DiamondOverlay( const float *fixed , const float *moving ,
                    const float *weights , int num_points );
    // individual sets of coords.
    DiamondOverlay( const std::vector<std::vector<float> > &fixed ,
                    const std::vector<std::vector<float> > &moving  );
    DiamondOverlay( const std::vector<std::vector<float> > &fixed ,
                    const std::vector<std::vector<float> > &moving ,
                    const std::vector<float> &weights );

    ~DiamondOverlay();

    void calc_p();
    void calc_rot_matrix();
    void get_rot_matrix( float rot[3][3] );
    void get_p( float **p );
    float get_e_zero();
    float get_e();

  };

}

#endif
