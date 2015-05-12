//
// file ShapeTanimoto.H
// David Cosgrove
// AstraZeneca
// 15th September 2009
//

#include "stddefs.H"
#include "ShapeTanimoto.H"

#include <oechem.h>

using namespace std;
using namespace OEChem;
using namespace OESystem;

namespace DACLIB {

  float atomic_num_to_rad( int atomic_num );

  // ************************************************************************
  ShapeTanimoto::ShapeTanimoto(OEMolBase &mol1 , OEMolBase &mol2) :
    mol1_( mol1 ) , mol2_( mol2 ) , p_( 2.70 ) , lambda_( 1.5514 ) ,
    shape_tani_( 0.0 ) {

    calculate_shape_tani();

  }

  // ************************************************************************
  void ShapeTanimoto::calculate_shape_tani() {

    assign_alphas( mol1_ , mol1_atom_alphas_ );
    assign_alphas( mol2_ , mol2_atom_alphas_ );

    mol1_vs_mol1_ = calc_gaussian_overlap( mol1_ , mol1_ ,
					   mol1_atom_alphas_ ,
					   mol1_atom_alphas_ );
    mol2_vs_mol2_ = calc_gaussian_overlap( mol2_ , mol2_ ,
					   mol2_atom_alphas_ ,
					   mol2_atom_alphas_ );
    mol1_vs_mol2_ = calc_gaussian_overlap( mol1_ , mol2_ ,
					   mol1_atom_alphas_ ,
					   mol2_atom_alphas_ );

    shape_tani_ = mol1_vs_mol2_ / ( mol1_vs_mol1_ + mol2_vs_mol2_ - mol1_vs_mol2_ );

  }

  // ************************************************************************
  void ShapeTanimoto::assign_alphas( OEMolBase &mol ,
				     vector<double> &alphas ) const {

    double kappa = M_PI / 1.3401; // pi / ( lambda ^ 2/3 )
    alphas = vector<double>( mol.GetMaxAtomIdx() , 0.0 );
    for( OEIter<OEAtomBase> atom = mol.GetAtoms( OEIsHeavy() ) ; atom ; ++atom ) {
      double rad = DACLIB::atomic_num_to_rad( atom->GetAtomicNum() );
      alphas[atom->GetIdx()] = kappa / ( rad * rad );
    }

  }

  // ************************************************************************
  double ShapeTanimoto::calc_gaussian_overlap( OEMolBase &mol1 ,
					       OEMolBase &mol2 ,
					       const vector<double> &alpha1s ,
					       const vector<double> &alpha2s ) const {

    double ret_val = 0.0;
    double p2 = p_ * p_;

    OEIter<OEAtomBase> at1 , at2;
    for( at1 = mol1.GetAtoms( OEIsHeavy() ) ; at1 ; ++at1 ) {
      float at1_cds[3];
      mol1.GetCoords( at1 , at1_cds );
      for( at2 = mol2.GetAtoms( OEIsHeavy() ) ; at2 ; ++at2 ) {
	float at2_cds[3];
	mol2.GetCoords( at2 , at2_cds );
	double delta12 = alpha1s[at1->GetIdx()] + alpha2s[at2->GetIdx()];
	double alpha12 = alpha1s[at1->GetIdx()] * alpha2s[at2->GetIdx()];
	double rij2[3];
	rij2[0] = DACLIB::square( at1_cds[0] - at2_cds[0] );
	rij2[1] = DACLIB::square( at1_cds[1] - at2_cds[1] );
	rij2[2] = DACLIB::square( at1_cds[2] - at2_cds[2] );
	double rij = rij2[0] + rij2[1] + rij2[2];
	double k12 = exp( -( alpha12 * rij ) / delta12 );
	double ob = M_PI / delta12;
	ret_val += p2 * k12 * sqrt( ob * ob * ob );
      }
    }
    return ret_val;

  }

} // EO namespace DACLIB
