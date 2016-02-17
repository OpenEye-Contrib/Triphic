//
// file PharmPoint.cc
// David Cosgrove
// AstraZeneca
// 12th November 2002
//
// This is the implementation file for the class PharmPoint.  It holds
// definitions for pharmacophore points used by various programs including Triphic,
// Loob and Plurality.

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

#include <oechem.h>

#include "stddefs.H"
#include "FileExceptions.H"
#include "PharmPoint.H"
#include "SMARTSExceptions.H"

#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
// hints for using boost::ptree for parsing XML came from
// http://akrzemi1.wordpress.com/2011/07/13/parsing-xml-with-boost
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

using namespace std;
using namespace OEChem;
using namespace OESystem;

// ******************************************************************************
// protected member functions
// ******************************************************************************
void PharmPoint::set_defaults() {

  hphobes_itmoc_ = hphobes_itmoc_alo_ = false;

}

// ******************************************************************************
// public member functions
// ******************************************************************************
bool PharmPoint::has_point_of_name( const string &point_name ) const {

  for( map<string,vector<string> >::const_iterator p = points_defs_.begin() ,
       ps = points_defs_.end() ; p != ps ; ++p ) {
    if( boost::iequals( point_name , p->first ) ) {
      return true;
    }
  }
  return false;

}

// ******************************************************************************
// read the points definitions from file
void PharmPoint::read_points_file( const string &file_name ) {

  ifstream points_stream( file_name.c_str() );

  if( !points_stream || !points_stream.good() ) {
    throw DACLIB::FileReadOpenError( file_name.c_str() );
  }

  string input_line;
  while( 1 ) {
    input_line = string( "" );
    getline( points_stream , input_line );
    if( input_line.empty() &&
        ( points_stream.eof() || !points_stream.good() ) ) {
      points_stream.clear();
      break;
    }
    if( input_line.empty() || '#' == input_line[0] )
      continue;
    istringstream is( input_line );
    string point_name , def_name;
    is >> point_name;
    if( point_name == "H_VECTOR_POINTS" ) {
      while( 1 ) {
        is >> def_name;
        if( is.fail() )
          break;
        h_vector_points_.push_back( def_name );
      }
    } else if( point_name == "LP_VECTOR_POINTS" ) {
      while( 1 ) {
        is >> def_name;
        if( is.fail() )
          break;
        lp_vector_points_.push_back( def_name );
      }
    } else if( point_name == "RING_NORMAL_POINTS" ) {
      while( 1 ) {
        is >> def_name;
        if( is.fail() )
          break;
        ring_normal_points_.push_back( def_name );
      }
    } else {
      vector<string> these_defs;
      while( 1 ) {
        is >> def_name;
        if( is.fail() )
          break;
        if( def_name == "ITMOC" ) {
          hphobes_itmoc_ = true;
          itmoc_label_ = point_name;
        } else if( def_name == "ITMOC_ALO" ) {
          hphobes_itmoc_alo_ = true;
          itmoc_alo_label_ = point_name;
        } else {
          these_defs.push_back( def_name ); // don't want ITMOC in as a SMARTS
        }
      }
      points_defs_.insert( make_pair( point_name , these_defs ) );
    }
  }

  if( hphobes_itmoc_ && hphobes_itmoc_alo_ ) {
    hphobes_itmoc_alo_ = false; // it's redundant, so ignore it
  }

}

// ************************************************************************
void PharmPoint::read_points_xml_file( const std::string &file_name ) {

  ifstream ifs( file_name.c_str() );
  if( !ifs || !ifs.good() ) {
    throw DACLIB::FileReadOpenError( file_name.c_str() );
  }

  read_points_xml_stream( ifs );

}

// ************************************************************************
void PharmPoint::read_points_xml_string( const std::string xml ) {

  istringstream iss( xml );
  read_points_xml_stream( iss );

}

// ************************************************************************
void PharmPoint::read_points_xml_stream( std::istream &is ) {

  // populate tree structure pt
  using boost::property_tree::ptree;
  ptree pt;
  boost::property_tree::read_xml( is, pt );

  BOOST_FOREACH( ptree::value_type const &v, pt.get_child( "features" ) ) {
    if( v.first == "feature" ) {
      string point_name;
      vector<string> smt_names;
      // there can be multiple 'smarts' nodes
      BOOST_FOREACH( ptree::value_type const &vf , v.second ) {
#ifdef NOTYET
        // leaving this in to remind me how to do it
        cout << "vf.first : " << vf.first << endl;
        cout << "vf.second.data() : " << vf.second.data() << endl;
#endif
        if( vf.first == "name" ) {
          point_name = vf.second.data();
        } else if( vf.first == "smarts" ) {
          smt_names.push_back( vf.second.data() );
        } else if( vf.first == "direction" ) {
          if( vf.second.data() == "normal" ) {
            ring_normal_points_.push_back( point_name );
          } else if( vf.second.data() == "h_vector" ) {
            h_vector_points_.push_back( point_name );
          } else if( vf.second.data() == "lp_vector" ) {
            lp_vector_points_.push_back( point_name );
          }
        } else if( vf.first == "algorithmic" ) {
          if( vf.second.data() == "aliphatic_only" ) {
            itmoc_alo_label_ = point_name;
            hphobes_itmoc_alo_ = true;
          } else if( vf.second.data() == "original" ||
                     vf.second.data() == "full" ) {
            itmoc_label_ = point_name;
            hphobes_itmoc_ = true;
          }
        }
      }
      points_defs_.insert( make_pair( point_name , smt_names ) );
    }
  }

  if( hphobes_itmoc_ && hphobes_itmoc_alo_ ) {
    hphobes_itmoc_alo_ = false; // it's redundant, so ignore it
  }

#ifdef NOTYET
  report_points_defined( cout );
#endif

}

// ************************************************************************
void PharmPoint::report_points_defined( ostream &os ) const {

  os << "Points defined : " << endl;
  map<string,vector<string> >::const_iterator p;
  vector<string>::const_iterator r;
  for( p = points_defs_.begin() ; p != points_defs_.end() ; p++ ) {
    os << p->first << " : ";
    for( r = p->second.begin() ; r != p->second.end() ; r++ )
      os << *r << " ";
    os << endl;
  }

  if( !ring_normal_points_.empty() ) {
    os << "Ring normals will be checked for sites of type : ";
    copy( ring_normal_points_.begin() , ring_normal_points_.end() ,
          ostream_iterator<string>( os , " " ) );
    os << endl;
  }

  if( !h_vector_points_.empty() ) {
    os << "Hydrogen vectors will be generated for sites of type : ";
    for( r = h_vector_points_.begin() ; r != h_vector_points_.end() ; r++ )
      os << *r << " ";
    os << endl;
  }

  if( !lp_vector_points_.empty() ) {
    os << "Lone-Pair vectors will be generated for sites of type : ";
    for( r = lp_vector_points_.begin() ; r != lp_vector_points_.end() ; r++ )
      os << *r << " ";
    os << endl;
  }

  if( hphobes_itmoc_ || hphobes_itmoc_alo_ )
    os << "Hydrophobes will be generated by the Catalyst method (JCICS, 34,"
       << " 1297-1308 (1994))" << endl;

}

// ************************************************************************
void PharmPoint::clear_data() { 

  points_defs_.clear();
  h_vector_points_.clear();
  lp_vector_points_.clear();
  ring_normal_points_.clear();
  unique_smarts_.clear();
  smarts_to_points_.clear();

  set_defaults();

}

// ****************************************************************************
// this one returns the labels of any that are missing
void PharmPoint::check_points_smarts( const vector<pair<string,string> > &input_smarts ) {

  map<string,vector<string> >::const_iterator p;
  vector<string>::const_iterator q;
  vector<pair<string,string> >::const_iterator r;
  for( p = points_defs_.begin() ; p != points_defs_.end() ; p++ ) {
    for( q = p->second.begin() ; q != p->second.end() ; q++ ) {
      unique_smarts_.insert( *q );
      for( r = input_smarts.begin() ; r != input_smarts.end() ; r++ ) {
        if( r->first == *q )
          break;
      }
      if( r == input_smarts.end() ) {
        ostringstream oss;
        oss << "SMARTS definition " << *q << " not found for point definition "
            << p->first << ".";
        throw DACLIB::SMARTSDefnError( oss.str().c_str() );
      }
    }
  }

  // build the map of which point each SMARTS is in
  smarts_to_points_ = vector<int>( unique_smarts_.size() , -1 );
  set<string>::iterator s;
  int i , j;
  for( s = unique_smarts_.begin() , j = 0 ; s != unique_smarts_.end() ; ++s , ++j ) {
    for( p = points_defs_.begin() , i = 0 ; p != points_defs_.end() ; ++p , ++i ) {
      for( q = p->second.begin() ; q != p->second.end() ; ++q ) {
        if( *s == *q ) {
          smarts_to_points_[j] = i;
          break;
        }
      }
    }
  }

  //#if DEBUG == 1
  cout << "Unique SMARTS : ";
  copy( unique_smarts().begin() , unique_smarts().end() ,
        ostream_iterator<string>( cout , " " ) );
  cout << endl;
  //#endif

}

// ****************************************************************************
// take the string, a points type name, and return the corresponding integer
// code - -1 if not found.
int PharmPoint::type_code_from_string( const string &type_string ) const {

  map<string,vector<string> >::const_iterator p;
  p = points_defs_.find( type_string );
  if( p == points_defs_.end() )
    return -1;
  else
    return distance( points_defs_.begin() , p );
  
}

// ****************************************************************************
// make a map of OESubSearch objects corresponding to the input SMARTS,
// keyed on the SMARTS name.  The name is the first string in the pair,
// the SMARTS definition (with all vector bindings expanded) in the second.
void build_oesubsearches( PharmPoint &pharm_points ,
                          const vector<pair<string,string> > &smarts_defs ,
                          map<string,OESubSearch *> &subs ) {

  map<string,vector<string> > &points_defs = pharm_points.points_defs();
  map<string,vector<string> >::iterator p , ps;
  vector<string>::iterator q;
  map<string,OESubSearch *>::iterator r;
  OESubSearch *next_subs;

  for( p = points_defs.begin() , ps= points_defs.end() ; p != ps ; ++p ) {
    if( p->second.empty() )
      continue; // point defined by key word (e.g. ITMOC, ITMOC_ALO) not SMARTS.

    for( q = p->second.begin() ; q != p->second.end() ; ++q ) {
      for( int i = 0 , is = smarts_defs.size() ; i < is ; ++i ) {
        if( *q == smarts_defs[i].first ) {
          r = subs.find( *q );
          if( r == subs.end() ) {
            // allow OESubSearch to rearrange for efficiency
            next_subs = new OESubSearch( smarts_defs[i].second.c_str() , true );
            subs.insert( make_pair( *q , next_subs ) );
            break;
          }
        }
      }
    }
  }

}
