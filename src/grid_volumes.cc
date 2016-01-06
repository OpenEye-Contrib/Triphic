//
// file grid_volumes.cc
// David Cosgrove
// AstraZeneca
// 23rd June 2006
//
// Takes a set of OEMolBases and calculates a load of volume stats for them.
// Uses a grid-based method because it's easy. Gaussian methods don't work so
// well for multiple overlaid molecules.

#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

#include <oechem.h>

#include <boost/shared_ptr.hpp>
#include "gtpl_defs.H"
#include "stddefs.H"
#include "AtomSphere.H"
#include "VolumeGrid.H"

using namespace boost;
using namespace std;
using namespace OEChem;
using namespace OESystem;

namespace DACLIB {

float atomic_num_to_rad( int atomic_num );

float VolumeGrid::grid_spacing_ = 0.25F;

// ***************************************************************************
void find_molecule_extremes( OEMolBase *mol , float *blf , float *trb ,
                             vector<float> &atom_rads ) {

  blf[0] = blf[1] = blf[2] = numeric_limits<float>::max();
  trb[0] = trb[1] = trb[2] = -numeric_limits<float>::max();

  float cds[3];
  OEIter<OEAtomBase> atom;

  float shell_b2 = AtomSphere::surface_thickness() / 2.0F;
  for( atom = mol->GetAtoms() ; atom ; ++atom ) {
    float rad = DACLIB::atomic_num_to_rad( atom->GetAtomicNum() );
    atom_rads.push_back( rad );
    rad += shell_b2;
    if( 1 == atom->GetAtomicNum() )
      continue;
    mol->GetCoords( atom , cds );
    if( cds[0] - rad < blf[0] ) blf[0] = cds[0] - rad;
    if( cds[1] - rad < blf[1] ) blf[1] = cds[1] - rad;
    if( cds[2] - rad < blf[2] ) blf[2] = cds[2] - rad;
    if( cds[0] + rad > trb[0] ) trb[0] = cds[0] + rad;
    if( cds[1] + rad > trb[1] ) trb[1] = cds[1] + rad;
    if( cds[2] + rad > trb[2] ) trb[2] = cds[2] + rad;
  }

}

// ***************************************************************************
void find_molecule_extremes( vector<OEMolBase *> &mols , float *blf ,
                             float *trb , vector<vector<float> > &atom_rads ) {

  blf[0 ] = blf[1] = blf[2] = numeric_limits<float>::max();
  trb[0 ] = trb[1] = trb[2] = -numeric_limits<float>::max();

  float tblf[3] , ttrb[3];
  OEIter<OEAtomBase> atom;

  for( int i = 0 , is = mols.size() ; i < is ; ++i ) {
    atom_rads.push_back( vector<float>() );
    find_molecule_extremes( mols[i] , tblf , ttrb , atom_rads.back() );
    if( tblf[0] < blf[0] ) blf[0] = tblf[0];
    if( tblf[1] < blf[1] ) blf[1] = tblf[1];
    if( tblf[2] < blf[2] ) blf[2] = tblf[2];
    if( ttrb[0] > trb[0] ) trb[0] = ttrb[0];
    if( ttrb[1] > trb[1] ) trb[1] = ttrb[1];
    if( ttrb[2] > trb[2] ) trb[2] = ttrb[2];

  }

}

// ********************************************************************
// make a VolumeGrid from the given conformation of the given molecule
VolumeGrid *prepare_mol_grid( OEMolBase *mol ) {

  vector<float> atom_rads;
  float blf[3] , trb[3];

  find_molecule_extremes( mol , blf , trb , atom_rads );

  VolumeGrid *grid = 0;
  grid = new VolumeGrid( blf , trb );

  grid->drop_molecule_in( mol , atom_rads );

  return grid;
  
}

// ********************************************************************
// make a VolumeGrid from the given conformation of the given molecule
VolumeGrid *prepare_mol_grid( OEMol *mol , int conf_num ) {

  OEIter<OEConfBase> conf = mol->GetConfs();
  for( int i = 0 ; i < conf_num ; ++i ) {
    ++conf;
  }

  return prepare_mol_grid( conf );
  
}

// *********************************************************************
// prepare a counts grid containing multiple molecules that will be used for
// more than one shape tanimoto calculation.
void prepare_multi_mol_grids( vector<OEMolBase *> &mols ,
                              boost::shared_ptr<VolumeGrid> &solid_grid ) {

  vector<VolumeGrid *> grids;
  for( int i = 0 , is = mols.size() ; i < is ; ++i )
    grids.push_back( prepare_mol_grid( mols[i] ) );

  solid_grid.reset( new VolumeGrid( grids , 1 ) );

  for( int i = 0 , is = mols.size() ; i < is ; ++i ) {
    delete grids[i];
  }

  solid_grid->counts_to_flag( 1 , GtplDefs::CORE );

}

} // end of namespace
