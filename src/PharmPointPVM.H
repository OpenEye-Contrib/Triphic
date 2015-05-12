//
// file PharmPointPVM.H
// David Cosgrove
// AstraZeneca
// 23rd May 2006
//
// This is the interface for the class PharmPointPVM which derives
// from PharmPoint but has the ability to pack itself into and unpack itself
// from a PVM buffer.  It's separate from PharmPoint so that GUIs that use a
// PharmPoint don't need to link to pvm.

#ifndef DAC_PHARM_POINT_PVM__
#define DAC_PHARM_POINT_PVM__

#include "PharmPoint.H"

// ****************************************************************************
class PharmPointPVM : public PharmPoint {

public :

  PharmPointPVM() : PharmPoint() {}
  virtual ~PharmPointPVM() {}

  // put the info into a previously initialised pvm buffer, but don't send it
  virtual void pack_into_pvm_buffer();
  // pull it off a pvm buffer
  virtual void unpack_from_pvm_buffer();

};

#endif

