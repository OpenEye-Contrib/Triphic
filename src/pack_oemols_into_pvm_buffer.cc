//
// file pack_oemols_into_pvm_buffer.cc
// David Cosgrove
// AstraZeneca
// 6th August 2007

#include <vector>

#include <oechem.h>

#include <pvm3.h>

using namespace std;
using namespace OEChem;

namespace DACLIB {

  void pack_string_raw( const string &str );
  void unpack_string_raw( string &str );

  // *************************************************************************
  void pack_oemol_into_pvm_buffer( OEMolBase *mol ) {

    int i = 1;
    if( !*mol ) {
      i = 0;
    }
    pvm_pkint( &i , 1 , 1 );
    if( i ) {
      oemolostream oms;
      oms.SetFormat( OEFormat::OEB );
      oms.openstring();
      oms << *mol;
      pack_string_raw( oms.GetString() );
    }

  }

  // *************************************************************************
  void pack_oemols_into_pvm_buffer( vector<OEMolBase *> &mols ) {

    int num_to_send = mols.size();
    pvm_pkint( &num_to_send , 1 , 1 );
    for( int i = 0 ; i < num_to_send ; ++i ) {
      pack_oemol_into_pvm_buffer( mols[i] );
    }

  }

  // *************************************************************************
  void unpack_oemol_from_pvm_buffer( OEMolBase *mol ) {

    int i;
    pvm_upkint( &i , 1 , 1 );
    // if i is 0, there was no molecule sent, so stick with empty one already
    // initialised
    if( i ) {
      string msg;
      unpack_string_raw( msg );

      oemolistream ims;
      ims.SetFormat( OEFormat::OEB );
      ims.openstring( msg );
      ims >> *mol;
    }

  }

  // *************************************************************************
  void unpack_oemols_from_pvm_buffer( vector<OEMolBase *> &mols ) {

    int num_to_rec;
    pvm_upkint( &num_to_rec , 1 , 1 );
    for( int i = 0 ; i < num_to_rec ; ++i ) {
      mols.push_back( OENewMolBase( OEMolBaseType::OEDefault ) );
      unpack_oemol_from_pvm_buffer( mols[i] );
    }

  }

} // end of namespace DACLIB
