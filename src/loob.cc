//
// file loob.cc
// David Cosgrove
// AstraZeneca
// 2nd October 2007
//
// Main function for program loob.
// This is new loob, v 4.0, that uses the same gear as triphic and
// plurality. Functionally, should be identical to previous incarnations.

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <oechem.h>

#include <boost/shared_ptr.hpp>

#include "BasePPhoreSite.H"
#include "DefaultPointsDefs.H"
#include "FileExceptions.H"
#include "ParseSMARTSXML.H"
#include "PharmPoint.H"
#include "LoobSettings.H"

using namespace std;
using namespace OEChem;

extern string BUILD_TIME; // in build_time.cc

// in loob_subs.cc
void read_smarts_file( const string &smarts_file ,
		       vector<pair<string,string> > &input_smarts ,
		       vector<pair<string,string> > &smarts_sub_defn );
void read_subset_file( const string &subset_file ,
		       vector<string> &subset_names );
void process_molecule( vector<pair<string,string> > &input_smarts ,
		       vector<pair<string,string> > &smarts_sub_defn ,
		       PharmPoint &pharm_points ,
		       const vector<string> &mol_subset ,
		       LoobSettings &loob_settings , OEMol &mol ,
		       const map<string,int> &fp_map ,
		       vector<int> &these_set_bits );
void make_fp_map( map<string,int> &fp_map , vector<string> &bit_decoding ,
		  int num_point_types , LoobSettings &loob_settings );

void output_fingerprints( LoobSettings &loob_settings , int mols_done ,
			  vector<pair<string,vector<int> > > &bits_set ,
			  const map<string,int> &fp_map ,
			  const map<string,vector<string> > &points_defs );

// **************************************************************************
int main( int argc , char **argv ) {

  cout << "Loob v4.0." << endl
       << "Built " << BUILD_TIME << " using OEToolkits version "
       << OEChemGetRelease() << "." << endl;

  LoobSettings loob_settings( argc , argv );

  PharmPoint pharm_points;
  vector<pair<string,string> > input_smarts , smarts_sub_defn;
  if( !loob_settings.smarts_file_.empty() && !loob_settings.points_file_.empty() ) {
    read_smarts_file( loob_settings.smarts_file_ , input_smarts ,
		      smarts_sub_defn );

    try {
      pharm_points.read_points_file( loob_settings.points_file_ );
    } catch( DACLIB::FileReadOpenError &e ) {
      cerr << e.what() << endl;
      exit( 1 );
    }
  } else {
    pharm_points.clear_data();
    pharm_points.read_points_xml_string( DACLIB::DEFAULT_POINT_DEFS );
    ParseSMARTSXML psx;
    psx.parse_string( DACLIB::DEFAULT_POINT_DEFS , input_smarts , smarts_sub_defn );
  }

  vector<string> mol_subset;
  if( !loob_settings.subset_file_.empty() ) {
    read_subset_file( loob_settings.subset_file_ , mol_subset );
  }

  oemolistream ims( loob_settings.db_file_.c_str() );
  if( !ims ) {
    throw DACLIB::FileReadOpenError( loob_settings.db_file_.c_str() );
  }

  if( loob_settings.single_conf_mols_ ) {
    // just to make it clear what test is being used - this is the default for
    // oemolistream so is not strictly necessary
    ims.SetConfTest( OEDefaultConfTest() );
  } else {
    ims.SetConfTest( OEAbsCanonicalConfTest() );
  }

  cout << "Distances :";
  for( int ii = 0 , iis = loob_settings.dist_bounds_.size() ; ii < iis ; ++ii ) {
    cout << " " << loob_settings.dist_bounds_[ii];
  }
  cout << endl;

  // make the mappings from bit label to bit number
  map<string,int> fp_map;
  vector<string> bit_decoding;
  make_fp_map( fp_map , bit_decoding , pharm_points.points_defs().size() ,
	       loob_settings );

  int mols_done = 0;
  OEMol mol;
  vector<pair<string,vector<int> > > set_bits;
  while( ims >> mol ) {
    vector<int> these_set_bits;
    OEAssignAromaticFlags( mol , OEAroModelDaylight );
    process_molecule( input_smarts , smarts_sub_defn , pharm_points ,
		      mol_subset , loob_settings , mol , fp_map ,
		      these_set_bits );
    set_bits.push_back( make_pair( mol.GetTitle() , these_set_bits ) );
    ++mols_done;
    if( !( mols_done % 10000 ) ) {
      cout << "Processed " << mols_done << " molecules." << endl;
    }
  }

  cout << "Finished processing " << mols_done << " molecules." << endl;

  output_fingerprints( loob_settings , mols_done , set_bits , fp_map ,
		       pharm_points.points_defs() );

}
