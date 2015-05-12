//
// file get_given_oeconf.cc
// David Cosgrove
// AstraZeneca
// 31st July 2007
//
// Takes an OEMol and a conformation number and returns a pointer to a NEW
// OEMolBase containing the given conformation.

#include <string>
#include <oechem.h>

#include <boost/lexical_cast.hpp>

using namespace std;
using namespace OEChem;
using namespace OESystem;

// ****************************************************************************
OEMolBase *get_given_oeconf( OEMol &mol , int conf_num ,
                             bool add_conf_num_to_name ) {

  int i = 0;
  OEIter<OEConfBase> conf;
  for( conf = mol.GetConfs() ; conf && i < conf_num ; ++i ) {
    ++conf;
  }
  if( !conf ) {
    return static_cast<OEMolBase *>( 0 );
  }

  OEMolBase *ret_conf = OENewMolBase( *conf , OEMolBaseType::OEDefault );

  if( add_conf_num_to_name ) {
    string conf_name = ret_conf->GetTitle();
    if( string::npos == conf_name.find( "_Conf" ) ) {
      conf_name += "_Conf" + boost::lexical_cast<string>( conf_num );
    }
    ret_conf->SetTitle( conf_name.c_str() );
  }

  return ret_conf;

}

