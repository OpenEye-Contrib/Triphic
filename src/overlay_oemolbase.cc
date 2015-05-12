//
// file overlay_oemolbase.cc
// David Cosgrove
// 18th September 2007
//
// Takes an OEMolBase and applies the given OverlayTrans to it.

#include "OverlayTrans.H"

#include <oechem.h>

using namespace OEChem;
using namespace OESystem;

// **************************************************************************
void overlay_oemolbase( OEMolBase &mol , const OverlayTrans &ot ) {

  float trans1[3] , trans2[3] , rot[3][3];
  ot.get_trans1( trans1 );
  ot.get_trans2( trans2 );
  ot.get_rot( rot );

  OEIter<OEAtomBase> atom;
  float cds[3];
  for( atom = mol.GetAtoms() ; atom ; ++atom ) {
    mol.GetCoords( atom , cds );
    DACLIB::translate( cds , -trans1[0] , -trans1[1] , -trans1[2] );
    DACLIB::rotate( cds , rot );
    DACLIB::translate( cds , trans2[0] , trans2[1] , trans2[2] );
    mol.SetCoords( atom , cds );
  }

}

