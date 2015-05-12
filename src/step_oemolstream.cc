//
// file step_oemolstream.cc
// David Cosgrove
// AstraZeneca
// 16th October 2014
//
// This function takes an oemolistream and moves it forward by
// the given number of molecule records. It returns step if it
// succeeds, the number of molecules it moved forwards otherwise,
// for example if it hits end of file

#include "FileExceptions.H"

#include <oechem.h>

#include <iostream>

// ****************************************************************************
unsigned int step_oemolstream( OEChem::oemolistream &ims , int step ) {

#ifdef NOTYET
  std::cout << "Stepping " << step << " steps" << std::endl;
#endif
  OEChem::OEMol mol;
  for( int i = 0 ; i < step ; ++i ) {
    if( !(ims >> mol) ) {
      return i;
    }
#ifdef NOTYET
    std::cout << i << ": Stepped past : " << mol.GetTitle() << std::endl;
#endif
  }

  return step;

}

// ********************************************************************
void open_databasefile( const std::string &db_file ,
                        bool single_conf_mols ,
                        OEChem::oemolistream &ims ) {

  if( !ims.open( db_file.c_str() ) ) {
    throw DACLIB::FileReadOpenError( db_file.c_str() );
  }

  if( single_conf_mols ) {
    // just to make it clear what test is being used - this is the default for
    // oemolistream so is not strictly necessary
    ims.SetConfTest( OEChem::OEDefaultConfTest() );
  } else {
    ims.SetConfTest( OEChem::OEAbsCanonicalConfTest() );
  }

}

// ********************************************************************
void open_databasefile( const std::string &db_file ,
                        bool single_conf_mols ,
                        OEChem::oemolistream *&ims ) {

  if( !ims->open( db_file.c_str() ) ) {
    throw DACLIB::FileReadOpenError( db_file.c_str() );
  }

  if( single_conf_mols ) {
    // just to make it clear what test is being used - this is the default for
    // oemolistream so is not strictly necessary
    ims->SetConfTest( OEChem::OEDefaultConfTest() );
  } else {
    ims->SetConfTest( OEChem::OEAbsCanonicalConfTest() );
  }

}
