//
// file build_sites_from_moe_features.cc
// David Cosgrove
// AstraZeneca
// 18th March 2015
//
// This file takes a long string containing a feature record read from a MOE
// pharmacophore file and turns it into a vector of GTPL sites for use as a
// query. The string will contribute 1 vector of sites.
// Sites will be selected if the MOE file defines them as essential

#include "stddefs.H"
#include "MOEPointsDefs.H"
#include "ParseSMARTSXML.H"
#include "PharmPoint.H"
#include "SinglePPhoreSite.H"

#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>

using namespace boost;
using namespace std;

// ****************************************************************************
string moe_name_to_point_name( const string &name ,
                               const vector<pair<string,string> > &moe_names_to_points_names ) {

  vector<pair<string,string> >::const_iterator p = find_if( moe_names_to_points_names.begin() ,
                                                            moe_names_to_points_names.end() ,
                                                            bind( equal_to<string>() ,
                                                                  name ,
                                                                  bind( &pair<string,string>::first , _1 ) ) );

  if( p == moe_names_to_points_names.end() ) {
    cerr << "Error : can't translate MOE PH4 type " << name << " to triphic equivalent." << endl;
    cout << "Error : can't translate MOE PH4 type " << name << " to triphic equivalent." << endl;
    exit( 1 );
  } else {
    return p->second;
  }

}

// ****************************************************************************
string point_name_to_moe_name( const string &name ,
                               const vector<pair<string,string> > &moe_names_to_points_names ) {

  vector<pair<string,string> >::const_iterator p = find_if( moe_names_to_points_names.begin() ,
                                                            moe_names_to_points_names.end() ,
                                                            bind( equal_to<string>() ,
                                                                  name ,
                                                                  bind( &pair<string,string>::second , _1 ) ) );

  if( p == moe_names_to_points_names.end() ) {
    cerr << "Error : can't translate triphic type " << name << " to MOE PH4 equivalent." << endl;
    cout << "Error : can't translate triphic type " << name << " to MOE PH4 equivalent." << endl;
    exit( 1 );
  } else {
    return p->second;
  }

}

// ****************************************************************************
void add_moe_site( double cds[3] , const string &type_string ,
                   PharmPoint &pharm_points ,
                   const vector<pair<string,string> > moe_names_to_points_names ,
                   bool sel_site ,
                   vector<BasePPhoreSite *> &sites ) {

  string triphic_name = moe_name_to_point_name( type_string , moe_names_to_points_names );

  int type_code = pharm_points.type_code_from_string( triphic_name );
  double dir[3] = { 0.0 , 0.0 , 0.0 };
  string label = string( "MOE_" ) + type_string + string( "_" ) + lexical_cast<string>( sites.size() + 1 );

  SinglePPhoreSite *site = new SinglePPhoreSite( cds , dir , type_code , triphic_name ,
                                                 label , false , "MOE_PH4_Query" );
  site->set_selected( sel_site );
  sites.push_back( site );

}

// ****************************************************************************
void map_moe_names_to_points_names( vector<pair<string,string> > &moe_names_to_points_names ) {

  // pharm_points should have a good mapping of the MOE features types as
  // they are set up in MOEPointsDefs.H, but there are some wrinkles.
  moe_names_to_points_names.push_back( make_pair( string( "Acc" ) , string( "Acc" ) ) );
  moe_names_to_points_names.push_back( make_pair( string( "Don" ) , string( "Don" ) ) );
  moe_names_to_points_names.push_back( make_pair( string( "Ani" ) , string( "Ani" ) ) );
  moe_names_to_points_names.push_back( make_pair( string( "Cat" ) , string( "Cat" ) ) );
  moe_names_to_points_names.push_back( make_pair( string( "Hyd" ) , string( "Hyd" ) ) );
  moe_names_to_points_names.push_back( make_pair( string( "Aro" ) , string( "Aro" ) ) );
  // some common combinations
  moe_names_to_points_names.push_back( make_pair( string( "Don!Cat" ) , string( "NeuDon" ) ) );
  moe_names_to_points_names.push_back( make_pair( string( "Acc!Ani" ) , string( "NeuAcc" ) ) );
  moe_names_to_points_names.push_back( make_pair( string( "Don|Cat" ) , string( "AllDon" ) ) );
  moe_names_to_points_names.push_back( make_pair( string( "Acc|Ani" ) , string( "AllAcc" ) ) );
  moe_names_to_points_names.push_back( make_pair( string( "Cat|Don" ) , string( "AllDon" ) ) );
  moe_names_to_points_names.push_back( make_pair( string( "Ani|Acc" ) , string( "AllAcc" ) ) );
  moe_names_to_points_names.push_back( make_pair( string( "Hyd|Aro" ) , string( "AllHyd" ) ) );
  moe_names_to_points_names.push_back( make_pair( string( "Aro|Hyd" ) , string( "AllHyd" ) ) );
  moe_names_to_points_names.push_back( make_pair( string( "Acc&Ani" ) , string( "NegAcc" ) ) );
  moe_names_to_points_names.push_back( make_pair( string( "Cat&Don" ) , string( "PosDon" ) ) );
  moe_names_to_points_names.push_back( make_pair( string( "Ani&Acc" ) , string( "NegAcc" ) ) );
  moe_names_to_points_names.push_back( make_pair( string( "Don&Cat" ) , string( "PosDon" ) ) );

  // Acc2 and Don2 are are acceptor and donor extension points and not supported for now.
  moe_names_to_points_names.push_back( make_pair( string( "Acc2" ) , string( "Acc2" ) ) );
  moe_names_to_points_names.push_back( make_pair( string( "Don2" ) , string( "Don2" ) ) );
  // normals for aromatic rings
  moe_names_to_points_names.push_back( make_pair( string( "PiN" ) , string( "PiN" ) ) );

}

// ****************************************************************************
// process all features into sites, irrespective of type
void make_raw_moe_sites( string &moe_features ,
                         PharmPoint &pharm_points ,
                         const vector<pair<string,string> > &moe_names_to_points_names ,
                         vector<BasePPhoreSite *> &new_sites ) {

  trim( moe_features );
  list<string> splits;
  split( splits , moe_features , is_any_of( " \n\t" ) );

#ifdef NOTYET
  cout << "Got " << splits.size() << " records, for " << splits.size() / 9 << " features." << endl;
#endif
  if( splits.size() % 9 ) {
    cerr << "ERROR : invalid feature record in MOE pharmacophore file." << endl;
    cerr << moe_features << endl;
    exit( 1 );
  }
  double cds[3];
  int essential = 0;
  string type_string;
  bool sel_site; // sites will be selected if essential

  while( !splits.empty() ) {
    // each feature is described by 9 fields
#ifdef NOTYET
    cout << "Next feature : ";
    list<string>::iterator p = splits.begin();
    for( int ii = 0 ; ii < 9 ; ++ii , ++p ) {
      cout << *p << " ";
    }
    cout << endl;
#endif

    type_string = splits.front(); splits.pop_front();
    splits.pop_front(); // the colour, which we're not interested in
    cds[0] = lexical_cast<double>( splits.front() ); splits.pop_front();
    cds[1] = lexical_cast<double>( splits.front() ); splits.pop_front();
    cds[2] = lexical_cast<double>( splits.front() ); splits.pop_front();
    splits.pop_front(); // the tolerance or radius, which we're not interested in
    essential = lexical_cast<int>( splits.front() ); splits.pop_front(); // 1 is yes, 0 no
    sel_site = bool( essential == 1 );
    splits.pop_front(); splits.pop_front(); // 2 internally used codes
    size_t pp = type_string.find( "|" );
    if( string::npos != pp ) {
      // an OR'd definition, make two of them
      add_moe_site( cds , type_string.substr( 0 , pp ) , pharm_points ,
                    moe_names_to_points_names , sel_site , new_sites );
      add_moe_site( cds , type_string.substr( pp + 1 ) , pharm_points ,
                    moe_names_to_points_names , sel_site , new_sites );
    } else {
      add_moe_site( cds , type_string , pharm_points ,
                    moe_names_to_points_names , sel_site , new_sites );
    }
  }

#ifdef NOTYET
  cout << "Made " << new_sites.size() << " raw MOE sites" << endl;
#endif

}

// ****************************************************************************
void find_nearest_site( BasePPhoreSite *query_site ,
                        const vector<BasePPhoreSite *> sites ,
                        const string &site_type ,
                        int &nearest_site_num , double &site_dist ) {

  nearest_site_num = -1;
  site_dist = numeric_limits<double>::max();
  for( int j = 0 , js = sites.size() ; j < js ; ++j ) {
    if( !sites[j] || query_site == sites[j] || site_type != sites[j]->get_type_string() ) {
      continue;
    }
    double dist = DACLIB::sq_distance( query_site->coords() , sites[j]->coords() );
    if( dist < site_dist ) {
      nearest_site_num = j;
      site_dist = dist;
    }
  }

#ifdef NOTYET
  if( -1 != nearest_site_num ) {
    cout << "Nearest site to " << query_site->label() << " is "
         << sites[nearest_site_num]->label() << " at " << site_dist << endl;
  }
#endif

}

// ****************************************************************************
void find_nearest_sites_within_dist( BasePPhoreSite *query_site ,
                                     const vector<BasePPhoreSite *> &sites ,
                                     const string &site_type ,
                                     double sq_cutoff_dist ,
                                     vector<int> &near_sites ) {

  near_sites.clear();
  for( int j = 0 , js = sites.size() ; j < js ; ++j ) {
    if( !sites[j] || query_site == sites[j] || site_type != sites[j]->get_type_string() ) {
      continue;
    }
    double sq_dist = DACLIB::sq_distance( query_site->coords() , sites[j]->coords() );
    if( sq_dist < sq_cutoff_dist ) {
      near_sites.push_back( j );
    }
  }

}

// ****************************************************************************
// take the coords of dir_site_num and use to define the direction for site_num
// and delete dir_site_num, setting value in sites to 0.
void combine_moe_site_and_dir( vector<BasePPhoreSite *> &sites ,
                               int site_num , int dir_site_num ) {

#ifdef NOTYET
  cout << "Combining " << sites[site_num]->label()
       << " " << sites[site_num]->get_type_string() << " with "
       << sites[dir_site_num]->label() << " " << sites[dir_site_num]->get_type_string()
       << " as dir" << endl;
#endif

  double dir_cds[3] = { sites[dir_site_num]->coords()[0] ,
                        sites[dir_site_num]->coords()[1] ,
                        sites[dir_site_num]->coords()[2]};

  dir_cds[0] -= sites[site_num]->coords()[0];
  dir_cds[1] -= sites[site_num]->coords()[1];
  dir_cds[2] -= sites[site_num]->coords()[2];
  DACLIB::normalise( dir_cds );

  GtplDefs::DIRS_TYPE t = GtplDefs::UNKNOWN;
  if( string( "Don" ) == sites[site_num]->get_type_string() ) {
    t = GtplDefs::H_VECTOR;
  } else if( string( "Acc" ) == sites[site_num]->get_type_string() ) {
    t = GtplDefs::LP_VECTOR;
  } else if( string( "Aro" ) == sites[site_num]->get_type_string() ) {
    t = GtplDefs::RING_NORMAL;
  }

  sites[site_num]->set_direction( dir_cds , t );
  delete sites[dir_site_num];
  sites[dir_site_num] = 0;

}

// ****************************************************************************
void refine_moe_sites( vector<BasePPhoreSite *> &sites ,
                       const vector<pair<string,string> > &moe_names_to_points_names ) {

  for( int i = 0 , is = sites.size() ; i < is ; ++i ) {
    if( !sites[i] ) {
      continue;
    }
    string ts = point_name_to_moe_name( sites[i]->get_type_string() ,
                                        moe_names_to_points_names );
    if( ts[ts.length() - 1] == '2' ) {
      ts = ts.substr( 0 , ts.length() - 1 );
      // It's an extension point. Find nearest corresponding site. If there is one,
      // combine the two. Otherwise, leave as is and the search will have to create
      // extension points on the database molecules.
      vector<int> near_sites;
      // 3.16 A away is plenty for an extension point
      find_nearest_sites_within_dist( sites[i] , sites , ts , 10.0F , near_sites );
      for( int j = 0 , js = near_sites.size() ; j < js ; ++j ) {
        combine_moe_site_and_dir( sites , near_sites[j] , i );
      }
    }
    if( string( "PiN" ) == ts ) {
      // It's a ring normal. Find the ring.
      int nearest_site_num;
      double site_dist;
      string pn = moe_name_to_point_name( string( "Aro" ) , moe_names_to_points_names );
      find_nearest_site( sites[i] , sites , pn , nearest_site_num , site_dist );
      if( -1 != nearest_site_num && fabs( site_dist - 2.1 * 2.1 ) ) {
        combine_moe_site_and_dir( sites , nearest_site_num , i );
      } else {
        cerr << "Warning : " << sites[i]->label()
             << " defines a ring normal without corresponding base point, so it will be ignored." << endl;
        cout << "Warning : " << sites[i]->label()
             << " defines a ring normal without corresponding base point, so it will be ignored." << endl;
        delete sites[i];
        sites[i] = 0;
      }
    }
  }

  sites.erase( remove( sites.begin() , sites.end() ,
                       static_cast<BasePPhoreSite *>( 0 ) ) , sites.end() );

}

// ****************************************************************************
void build_sites_from_moe_features( string &moe_features ,
                                    PharmPoint &pharm_points ,
                                    vector<pair<string,string> > &input_smarts ,
                                    vector<pair<string,string> > &smarts_sub_defn ,
                                    vector<BasePPhoreSite *> &query_sites ) {

#ifdef NOTYET
  cout << "build_sites_from_moe_features : " << moe_features << endl;
#endif
  if( !pharm_points.points_defs().empty() ) {
    cout << "Warning : MOE pharmacophore files use their own canned pharmacophore feature definitions." << endl
         << "The ones previously loaded will be over-ridden." << endl;
    pharm_points.clear_data();
  }
  pharm_points.read_points_xml_string( DACLIB::MOE_POINT_DEFS );
  input_smarts.clear();
  smarts_sub_defn.clear();
  ParseSMARTSXML psx;
  psx.parse_string( DACLIB::MOE_POINT_DEFS , input_smarts , smarts_sub_defn );

  vector<pair<string,string> > moe_names_to_points_names;
  map_moe_names_to_points_names( moe_names_to_points_names );

  // first line of incoming string (up to \n) is a header line
  size_t i = 0;
  while( '\n' != moe_features[i] ) {
    ++i;
  }
  moe_features = moe_features.substr( i + 1 );
  // process the features into sites, irrespective of type string
  make_raw_moe_sites( moe_features , pharm_points , moe_names_to_points_names ,
                      query_sites );
  // various refinements because some of the sites are extension points
  // which may need to be made directions on the parent site
  refine_moe_sites( query_sites , moe_names_to_points_names );

#ifdef NOTYET
  cout << "Made " << query_sites.size() << " MOE query sites" << endl;
#endif

}

