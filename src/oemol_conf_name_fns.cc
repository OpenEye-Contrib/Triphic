//
// file oemol_conf_name_fns.cc
// David Cosgrove
// AstraZeneca
// May 17th 2007
//
// This file contains a couple of functions for dealing with conformer names
// in oemols, as used by Triphic inter alia.

#include <string>
#include <boost/lexical_cast.hpp>

using namespace std;

// ****************************************************************************
// take the _Conf<nn> off the end of the molecule name
string root_mol_name( const string &full_mol_name ) {

  size_t pos = full_mol_name.rfind( "_Conf" );
  if( string::npos == pos ||
      ( full_mol_name.length() &&
	!isdigit( full_mol_name[full_mol_name.length()-1] ) ) ) {
    return full_mol_name;
  }

  return full_mol_name.substr( 0 , pos );

}

// ****************************************************************************
// extract the name and conf number from the name - return -1 for the conf number
// if there is no _Conf<nn>
void extract_mol_name_and_conf( const string &full_mol_name ,
				string &mol_name , int &conf_num ) {

  conf_num = -1;
  mol_name = "";
  size_t pos = full_mol_name.rfind( "_Conf" );
  if( string::npos == pos ) {
    mol_name = full_mol_name;
    return;
  }

  mol_name = full_mol_name.substr( 0 , pos );
  conf_num = boost::lexical_cast<int>( full_mol_name.substr( pos + 5 ) );

}

