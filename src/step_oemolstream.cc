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
#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

typedef boost::shared_ptr<OEChem::OEMol> pOEMol;

using namespace std;

namespace DACLIB {

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
void open_databasefile( const string &db_file ,
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

// **************************************************************************
// read the molecule of given sequence number from the set of database files
// returns empty shared_ptr if unsuccessful
pOEMol read_nth_mol_from_oemolistream( unsigned int next_mol ,
                                       const vector<string> &db_files ,
                                       bool single_conf_mols ) {

  using namespace OEChem;
  static oemolistream *db_ims = 0;
  static unsigned int curr_db_file = 0;
  static unsigned int curr_mol = 0;

  if( !db_ims ) {
    db_ims = new oemolistream;
    cout << "Searching database file " << db_files[curr_db_file] << endl;
    open_databasefile( db_files[curr_db_file] ,
                       single_conf_mols , db_ims );
  }

  static OEMol mol;
  while( curr_mol <= next_mol ) {
    while( true ) {
      if( !(*db_ims >> mol) ) {
        db_ims->close();
        ++curr_db_file;
        if( curr_db_file == db_files.size() ) {
          // we're done
          return pOEMol();
        }
        cout << "Moving search to database file " << db_files[curr_db_file] << endl;
        open_databasefile( db_files[curr_db_file] ,
                           single_conf_mols , db_ims );
      } else {
        break;
      }
    }
    ++curr_mol;
  }

  return pOEMol( new OEMol( mol ) );

}

} // EO namespace DACLIB
