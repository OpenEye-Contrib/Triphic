//
// file create_cansmi.cc
// David Cosgrove
// AstraZeneca
// 23rd December 2011
//
// Convenience function to return a canonical SMILES from an OEMol.

#include <oechem.h>

#include <string>

namespace DACLIB {

std::string create_cansmi( const OEChem::OEMolBase &in_mol ) {

  std::string smi;
  // default for OECreateSmiString is RGroups, Canonical and AtomMaps, but I don't think we'll
  // ever want the atom maps and they confuse matters.
  OEChem::OECreateSmiString( smi , in_mol , OEChem::OESMILESFlag::RGroups | OEChem::OESMILESFlag::Canonical | OEChem::OESMILESFlag::AtomStereo | OEChem::OESMILESFlag::BondStereo );
  return smi;

}
}
