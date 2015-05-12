//
// file ShapeTanimoto.H
// David Cosgrove
// AstraZeneca
// 15th September 2009
//
// Computes the approximate shape tanimoto coefficient for the 2 OEMolBases
// passed in.  Uses 1st order intersections only, so not a reliable way
// of calculating volumes, but good enough for the tanimoto according to
// JAG. It's how ROCS does it, or so he told me.

#ifndef DAC_SHAPE_TANIMOTO
#define DAC_SHAPE_TANIMOTO

#include <vector>

// ******************************************************************************
namespace OEChem {
  class OEMolBase;
}

namespace DACLIB {

  class ShapeTanimoto {

  public :

    ShapeTanimoto( OEChem::OEMolBase &mol1 , OEChem::OEMolBase &mol2 );
    double shape_tanimoto() { return shape_tani_; }

  private :

    OEChem::OEMolBase &mol1_ , &mol2_;

    const double p_;
    const double lambda_;
    std::vector<double> mol1_atom_alphas_;
    std::vector<double> mol2_atom_alphas_;

    double mol1_vs_mol1_ , mol2_vs_mol2_ , mol1_vs_mol2_;
    double shape_tani_;

    void calculate_shape_tani();
    void assign_alphas( OEChem::OEMolBase &mol ,
			std::vector<double> &alphas ) const;
    double calc_gaussian_overlap( OEChem::OEMolBase &mol1 ,
				  OEChem::OEMolBase &mol2 ,
				  const std::vector<double> &alpha1s ,
				  const std::vector<double> &alpha2s ) const;

  };

} // EO namespace DACLIB

#endif
