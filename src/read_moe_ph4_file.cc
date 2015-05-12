//
// file read_moe_ph4_file.cc
// David Cosgrove
// AstraZeneca
// 1st April 2015
//
// Functions for reading pharmacophore files written by CCG's MOE program
// into a format usable by triphic and friends.
// The format has been deduced by examining the entrails of a number of such
// files rather than from any documentation that might exist, so it may not
// be entirely correct and is certainly incomplete.
// It reads features and volumes, but not lots of the fancy combination flags
// for the features.

#include "BasePPhoreSite.H"
#include "FileExceptions.H"
#include "PharmPoint.H"
#include "TriphicSettings.H"
#include "VolumeGrid.H"

#include <fstream>
#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

using namespace boost;
using namespace std;

// in eponymous files
// moe_features will contribute 1 vector of sites to query_sites. moe_features,
// pharm_points, input_smarts and smarts_sub_defn will be altered by function.
void build_sites_from_moe_features( string &moe_features ,
                                    PharmPoint &pharm_points ,
                                    vector<pair<string,string> > &input_smarts ,
                                    vector<pair<string,string> > &smarts_sub_defn ,
                                    vector<vector<BasePPhoreSite *> > &query_sites );
DACLIB::VolumeGrid *build_vol_from_moe_volumesphere( string &vol_rec );

// ***************************************************************************
string read_next_moe_record( const string &file_contents ,
                              size_t &fs ) {

  string next_rec;

  if( '#' == file_contents[fs] ) {
    ++fs; // we might be at the start of the record, want to go on to end of next one
  }
  // everything up to next #
  while( fs != file_contents.size() && '#' != file_contents[fs] ) {
    if( '\r' == file_contents[fs] ) {
      ; // bloody DOS users!
    } else {
      next_rec += file_contents[fs];
    }
    ++fs;
  }

  return next_rec;

}

// ***************************************************************************
void decode_selected_moe_sites( TriphicSettings &ts ,
                                vector<BasePPhoreSite *> &sites ) {

  for( int i = 0 , is = sites.size() ; i < is ; ++i ) {
    string site_name = sites[i]->get_full_name();
    if( sites[i]->get_selected() ) {
      // sites that are selected are marked essential in the feature record,
      // so need to be added to ts.required_sites_and()
      ts.add_req_sites_and( site_name );
    }
    if( site_name.substr( 0 , 22 ) == string( "MOE_PH4_Query:MOE_Aro_" ) &&
        sites[i]->get_num_dirs( GtplDefs::RING_NORMAL ) ) {
      ts.set_ring_norm_usage( GtplDefs::ALIGN );
    }
  }

}

// ***************************************************************************
void read_moe_ph4_file( TriphicSettings &ts ,
                        PharmPoint &pharm_points ,
                        vector<pair<string,string> > &input_smarts ,
                        vector<pair<string,string> > &smarts_sub_defn ,
                        vector<vector<BasePPhoreSite *> > &query_sites ,
                        vector<pair<string,DACLIB::VolumeGrid *> > &score_vol_grids ) {

#ifdef NOTYET
  cout << "Reading MOE PH4 file " << ts.query_file() << endl;
#endif

  ifstream ifs( ts.query_file().c_str() );
  if( !ifs.good() ) {
    throw DACLIB::FileReadOpenError( ts.query_file().c_str() );
  }

  // slurp the whole file up.
  string file_contents;
  char c;
  while( ifs.get( c ) ) {
    file_contents += c;
  }

  // process any feature records
  size_t fs = 1;
  string next_rec;
  while( 1 ) {
    next_rec = read_next_moe_record( file_contents , fs );
    if( next_rec.empty() ) {
      break;
    }
    if( string( "feature" ) == next_rec.substr( 0 , 7 ) ) {
#ifdef NOTYET
      cout << "MOE feature : " << next_rec << endl;
#endif
      // string feat_rec = read_next_moe_record( file_contents , fs );
      build_sites_from_moe_features( next_rec , pharm_points , input_smarts ,
                                     smarts_sub_defn , query_sites );
      // some sites might be selected for various reasons. Translate this.
      decode_selected_moe_sites( ts , query_sites.back() );
    } else if( string( "volumesphere" ) == next_rec.substr( 0 , 12 ) ) {
#ifdef NOTYET
      cout << "MOE volumesphere" << endl;
#endif
      // string vol_rec = read_next_moe_record( file_contents , fs );
      DACLIB::VolumeGrid *moe_vol = build_vol_from_moe_volumesphere( next_rec );
      string vol_label = string( "MOE_Volume_" ) + lexical_cast<string>( score_vol_grids.size() + 1 );
      score_vol_grids.push_back( make_pair( vol_label , moe_vol ) );
    }

  }

}

